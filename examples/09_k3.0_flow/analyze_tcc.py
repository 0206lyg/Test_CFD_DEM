#!/usr/bin/env python3
"""Exclusive TCC particle populations and contacted-particle force/torque residuals.

RUN from the case containing particleData/:
    python3 analyze_tcc.py

Standard library only; no companion scripts. All nonnegative numeric .bin
snapshot times are processed in numerical order. Two outputs, one row per time:
  postProcessing/tccCounts.dat
    time 4A 5A 6Z 6A 7A unclassified total frac_4A ... frac_unclassified
  postProcessing/tccResiduals.dat
    time upper_mean_etaF_4A ... upper_mean_etaT_4A ... upper_N_F_4A ...
         upper_N_T_4A ... pore_mean_etaF_4A ... pore_mean_etaT_4A ...
         pore_N_F_4A ... pore_N_T_4A ...
    Every block has the six groups in the count-file order.

PARTICLE CLASSIFICATION
Discover motifs using the unchanged symmetry-filtered TCC-family predicates
below, then assign each particle once with priority:
    7A > 6A > 6Z > 5A > 4A > unclassified.
The five motif counts plus unclassified equal the upper-reservoir particle
count; fractions sum to one when that count is nonzero. These are PARTICLE
counts, replacing the old overlapping CLUSTER-count output. The public
classify_clusters(points, cutoff) still returns overlapping member sets, so
nested motifs remain available for constructing higher motifs and for callers.
6A-over-6Z is the agreed labeling priority, not a tetrahedral hierarchy.

REGIONS AND CUTOFF (all coordinates are m)
* Upper reservoir: x,z=[-0.020,0.020], y=[0.025,0.075]. Motif classification
  uses this entire region; all member centers must be within it. Static
  particles are included. No periodic wrapping is applied.
* Local pore prism: x,z=[-0.005,0.005], y=[0.025,0.035], center selection.
  This is the 10 mm layer immediately above the pore; it is not a first-hit
  ray selection. Its TCC labels are inherited from the upper classification.
  This preset interprets the requested near-pore layer in the case's SI
  coordinates: literal y=0.25--0.35 m lies outside this reservoir. Edit
  PORE_BOUNDS_M or pass --pore-bounds to specify another local region.
* Both ends of each bound are included, with BOUND_EPS tolerance.
* Bonds require r_ij < fixed cutoff. BOND_CUTOFF_M/--bond-cutoff takes
  precedence; otherwise cutoff=BOND_CUTOFF_FACTOR*reference diameter.
  The default factor is 1.41, matching the preceding two-snapshot analysis
  (the original GitHub script defaulted to 1.295). The reference diameter
  is the median (6V/pi)^(1/3) in the first populated upper snapshot, fixed
  for every subsequent time; --diameter can supply it explicitly.
  For spheres this equals the geometric diameter; for STL it is only a
  volume-equivalent length, not a contact distance. Mixed sizes use one
  pooled cutoff. No automatic g(r) first-minimum estimation is performed.

CONTACTED-ONLY RESIDUAL STATISTICS
TCC classification is independent of actual contact. Select particles with
at least one evaluated pp or wall contact having |Fn+Ft|>CONTACT_FORCE_EPS_N.
This threshold selects particles only. Every recorded contact, including
subthreshold contacts and partners outside either region or motif, contributes:
  Fc = Fn + Ft; tau_c = (contact_point - owner_center) cross Fc
  etaF = |FH + sum Fc| / (|FH| + sum |Fc|)
  etaT = |TH + sum tau_c| / (|TH| + sum |tau_c|)
Group means are arithmetic means of individual ratios. N_F/N_T are the
actual corresponding sample counts. Empty groups and undefined ratios yield
nan, not zero. No velocity criterion is imposed. For sphere pairs the
solver's contact-plane point is reconstructed from both radii and centers;
wall moments use the stored point and full normal+tangential force.
The bin lacks non-sphere pp contact points and a guaranteed STL torque
reference, so affected torque ratios are nan; force statistics remain valid.
Missing legacy force data likewise affect residuals, never TCC populations.
The stored loads omit separate gravity/adhesion; interpreting the ratios as
full mechanical balance assumes those omitted loads are zero, as in the
reference flow case. Do not reinterpret these as a complete residual in a
case with additional nonzero loads.

The .bin header, index and consumed arrays are checked; corrupt/partial
archives abort before replacing either previous result. Stored original files
and STL geometry are never extracted, and input files are not modified.

DEFINITIONS
4A: six bonds of K4 and all 12 vertex bond angles strictly between 50 and 70 deg.
5A: two qualifying 4A tetrahedra sharing exactly one triangular face. Their
    apices are on opposite sides of that face and are not bonded. Thus the
    five-particle core has nine internal bonds.
6Z: three qualifying tetrahedra around one common bond u-v, with the other
    four vertices forming a chordless path a-b-c-d. The tetrahedra are
    (u,v,a,b), (u,v,b,c), (u,v,c,d). Adjacent tetrahedra occupy opposite sides
    of their shared faces. Six vertices, twelve internal bonds.
6A: a chordless square with two unbonded poles, each bonded to all four ring
    vertices. Ring edge angles: 90 +/- 12 deg; ring centroid rotation angles:
    90 +/- 10 deg; adjacent oriented local normals differ by <10 deg. At the
    six-particle centroid, the 12 angles for ring-edge and ring-pole pairs
    must also lie strictly between 80 and 100 deg.
7A: a chordless pentagon with two poles, each bonded to all five ring vertices
    AND to each other, as required explicitly in the paper's SI. Ring edge
    angles: 108 +/- 12 deg; ring centroid rotation angles: 72 +/- 10 deg;
    adjacent oriented local normals differ by <10 deg. All bounds are strict.

The 4A and ring/angle predicates come from Tsurusawa & Tanaka's paper/SI.
The face-sharing 5A and common-edge 6Z constructions complete the unspecified
2-/3-tetrahedra rules using their diagrams and conventional TCC motifs. The
induced-core nonbond constraints and opposite-face condition are explicit
choices here. In particular, excluding the 6Z terminal a-d bond and the 6A
pole-pole bond is stricter than an unconstrained spindle-overlap search.
This is a documented, symmetry-filtered TCC-family implementation; it does
not claim to reproduce the authors' unavailable program byte for byte.


Sources for the unchanged geometry predicates:
https://doi.org/10.1038/s41567-023-02063-x
https://arxiv.org/abs/1307.5517
Based on Test_CFD_DEM/examples/09_k3.0_flow/analyze_tcc.py
GitHub blob b8bacdd5a1b6f982da7bfe1370bf51ce3ae85bcd.
"""

# -------------------- EDITABLE DEFAULTS --------------------
CASE_DIRECTORY = "."
# xmin, xmax, ymin, ymax, zmin, zmax; every coordinate is in m.
BOUNDS_M = (-0.020, 0.020, 0.025, 0.075, -0.020, 0.020)
PORE_BOUNDS_M = (-0.005, 0.005, 0.025, 0.035, -0.005, 0.005)
BOUND_EPS = 1.0e-12
BOND_CUTOFF_M = None           # Explicit distance overrides factor/diameter.
BOND_CUTOFF_FACTOR = 1.41      # Agreed two-snapshot analysis setting.
REFERENCE_DIAMETER_M = None    # Infer once from first populated upper snapshot.
CONTACT_FORCE_EPS_N = 1.0e-15  # Particle eligibility only, not a force-term cutoff.
OUTPUT_FILE = "postProcessing/tccCounts.dat"
RESIDUAL_OUTPUT_FILE = "postProcessing/tccResiduals.dat"
# -----------------------------------------------------------

from array import array
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import sys


_TIME_STEM = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")


def time_snapshots(directory):
    """Return every existing nonnegative numeric time, sorted numerically.

    Original filename stems remain the time labels. Non-time files are ignored;
    duplicate numeric times (e.g. 1.bin and 1.0.bin) are rejected.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError("particleData directory does not exist: {}".format(directory))
    selected = {}
    for path in directory.glob("*.bin"):
        if not path.is_file() or not _TIME_STEM.fullmatch(path.stem):
            continue
        try:
            time = Decimal(path.stem)
        except InvalidOperation:
            continue
        if not time.is_finite() or time < 0:
            continue
        if time in selected:
            raise ValueError("duplicate numeric time: {} and {}".format(
                selected[time].name, path.name))
        selected[time] = path
    if not selected:
        raise ValueError("no nonnegative numeric-time .bin files in {}".format(directory))
    return [(selected[time].stem, selected[time]) for time in sorted(selected)]


_ANALYSIS_HEADER = struct.Struct("<16sI")
_ANALYSIS_FOOTER = struct.Struct("<16sQQ32s")
_ANALYSIS_DTYPES = {"<f8": ("d", 8), "<i8": ("q", 8),
                    "<u8": ("Q", 8), "u1": ("B", 1)}


def read_analysis_snapshot(path, contact_force_threshold=1.0e-15):
    """Return geometry and per-particle residuals without filtering particles.

    ``contacted`` means at least one recorded particle or wall contact with
    norm(Fn + Ft) greater than contact_force_threshold in N. This threshold
    affects selection only: every recorded contact contributes to residuals.
    Force and torque denominators are sums of individual vector magnitudes.
    Unknown force data, zero denominators and unreconstructable non-spherical
    particle-contact torques produce NaN, never an invented zero residual.
    STL owner torques are unavailable because their stored vertex-mean center
    is not guaranteed to equal the solver's torque-reference center.
    """
    path = Path(path)
    if not math.isfinite(contact_force_threshold) or contact_force_threshold < 0:
        raise ValueError("contact force threshold must be finite and nonnegative")
    try:
        with path.open("rb") as stream:
            size = os.fstat(stream.fileno()).st_size
            if size < _ANALYSIS_HEADER.size + _ANALYSIS_FOOTER.size:
                raise ValueError("truncated binary archive")
            if _ANALYSIS_HEADER.unpack(stream.read(_ANALYSIS_HEADER.size)) != (
                    b"PDBIN-HEADER-v1!", 1):
                raise ValueError("unsupported particleData header/version")
            stream.seek(-_ANALYSIS_FOOTER.size, os.SEEK_END)
            magic, offset, length, digest = _ANALYSIS_FOOTER.unpack(
                stream.read(_ANALYSIS_FOOTER.size))
            if (magic != b"PDBIN-FOOTER-v1!" or offset < _ANALYSIS_HEADER.size
                    or offset + length != size - _ANALYSIS_FOOTER.size):
                raise ValueError("invalid index footer")
            stream.seek(offset)
            payload = stream.read(length)
            if len(payload) != length or hashlib.sha256(payload).digest() != digest:
                raise ValueError("index SHA256 mismatch")
            index = json.loads(payload)
            if (not isinstance(index, dict)
                    or index.get("format") != "openHFDIB-DEM particleData binary"
                    or type(index.get("schema_version")) is not int
                    or index.get("schema_version") != 1):
                raise ValueError("unsupported particleData schema")
            if index.get("time_name") != path.stem:
                raise ValueError("index time_name does not match filename")
            if index.get("length_unit") != "m" or index.get("time_unit") != "s":
                raise ValueError("expected length_unit=m and time_unit=s")
            if index.get("byte_order") != "little" or index.get("compression") != "none":
                raise ValueError("unsupported byte order or compression")
            datasets = index.get("datasets")
            if not isinstance(datasets, dict):
                raise ValueError("missing datasets index")
            ranges = []
            for name, item in datasets.items():
                if not isinstance(item, dict):
                    raise ValueError("invalid dataset metadata: {}".format(name))
                shape = item.get("shape")
                start, count = item.get("offset"), item.get("size")
                dtype, checksum = item.get("dtype"), item.get("sha256")
                if (not isinstance(shape, list)
                        or any(type(v) is not int or v < 0 for v in shape)
                        or type(start) is not int or type(count) is not int
                        or start < _ANALYSIS_HEADER.size or count < 0
                        or start + count > offset or not isinstance(checksum, str)
                        or re.fullmatch(r"[0-9a-f]{64}", checksum) is None):
                    raise ValueError("invalid dataset metadata: {}".format(name))
                if dtype in _ANALYSIS_DTYPES:
                    if math.prod(shape) * _ANALYSIS_DTYPES[dtype][1] != count:
                        raise ValueError("invalid dataset byte count: {}".format(name))
                elif dtype != "json-strings":
                    raise ValueError("unsupported dataset dtype: {}".format(name))
                if count:
                    ranges.append((start, start + count))
            ranges.sort()
            if any(a[1] > b[0] for a, b in zip(ranges, ranges[1:])):
                raise ValueError("overlapping binary datasets")

            def read_array(name, dtype, expected_shape):
                item = datasets.get(name)
                if (not isinstance(item, dict) or item.get("dtype") != dtype
                        or item.get("shape") != expected_shape):
                    raise ValueError("{} must have dtype {} and shape {}".format(
                        name, dtype, expected_shape))
                stream.seek(item["offset"])
                raw = stream.read(item["size"])
                if (len(raw) != item["size"]
                        or hashlib.sha256(raw).hexdigest() != item["sha256"]):
                    raise ValueError("{} SHA256 mismatch".format(name))
                code, width = _ANALYSIS_DTYPES[dtype]
                values = array(code)
                if values.itemsize != width:
                    raise ValueError("unsupported native array item size")
                values.frombytes(raw)
                if sys.byteorder != "little" and width > 1:
                    values.byteswap()
                return list(values)

            position = datasets.get("particles/position", {})
            shape = position.get("shape")
            if not isinstance(shape, list) or len(shape) != 2 or shape[1] != 3:
                raise ValueError("positions must have shape [N, 3]")
            n = shape[0]
            body_ids = read_array("particles/body_id", "<i8", [n])
            flat_points = read_array("particles/position", "<f8", [n, 3])
            volumes = read_array("particles/volume", "<f8", [n])
            if len(set(body_ids)) != n:
                raise ValueError("duplicate body IDs")
            if not all(math.isfinite(v) for v in flat_points):
                raise ValueError("nonfinite particle position")
            if not all(math.isfinite(v) and v > 0 for v in volumes):
                raise ValueError("particle volumes must be finite and positive")
            points = [tuple(flat_points[i:i + 3]) for i in range(0, 3 * n, 3)]
            warnings = []
            metadata = {"warnings": warnings, "force_unavailable_count": n,
                        "torque_geometry_unavailable_count": 0,
                        "zero_force_denominator_count": 0,
                        "zero_torque_denominator_count": 0,
                        "contact_force_threshold_N": contact_force_threshold}
            result = {"body_ids": body_ids, "points": points, "volumes": volumes,
                      "force_available": [False] * n, "contacted": [False] * n,
                      "eta_force": [float("nan")] * n,
                      "eta_torque": [float("nan")] * n,
                      "active_contact_count": [0] * n,
                      "force_schema": False, "metadata": metadata}
            expected = {
                "particles/fluid_force": ("<f8", [n, 3]),
                "particles/fluid_torque": ("<f8", [n, 3]),
                "particles/force_data_available": ("u1", [n]),
                "contacts/particle_offsets": ("<u8", [n + 1]),
                "contacts/wall_offsets": ("<u8", [n + 1]),
                "contacts/particle_partner_id": ("<i8", None),
                "contacts/particle_normal_force": ("<f8", None),
                "contacts/particle_tangential_force": ("<f8", None),
                "contacts/wall_contact_point": ("<f8", None),
                "contacts/wall_normal_force": ("<f8", None),
                "contacts/wall_tangential_force": ("<f8", None),
            }
            if not any(name in datasets for name in expected):
                warnings.append("archive has no stored force/contact schema")
                return result
            if not all(name in datasets for name in expected):
                raise ValueError("incomplete optional force/contact schema")
            force_metadata = index.get("force_output")
            if (not isinstance(force_metadata, dict)
                    or force_metadata.get("version") != 1
                    or force_metadata.get("force_unit") != "N"
                    or force_metadata.get("torque_unit") != "N m"
                    or force_metadata.get("contact_point_unit") != "m"
                    or force_metadata.get("vector_convention") != "force on owning particle"):
                raise ValueError("unsupported force metadata, units or vector convention")
            metadata["force_output"] = force_metadata
            particle_shape = datasets["contacts/particle_partner_id"]["shape"]
            wall_shape = datasets["contacts/wall_contact_point"]["shape"]
            if (len(particle_shape) != 1 or len(wall_shape) != 2
                    or wall_shape[1] != 3):
                raise ValueError("invalid contact array shapes")
            nparticle, nwall = particle_shape[0], wall_shape[0]
            for name in tuple(expected):
                dtype, expected_shape = expected[name]
                if expected_shape is None:
                    count = nparticle if name.startswith("contacts/particle_") else nwall
                    expected_shape = [count] if name.endswith("partner_id") else [count, 3]
                expected[name] = (dtype, expected_shape)
            data = {name: read_array(name, dtype, shape)
                    for name, (dtype, shape) in expected.items()}
            geometry = read_array("particles/geometry_type", "u1", [n])
            radii = read_array("particles/radius", "<f8", [n])
            if any(g not in (0, 1) for g in geometry):
                raise ValueError("unsupported particle geometry type")
            if any(g == 0 and (not math.isfinite(r) or r <= 0)
                   for g, r in zip(geometry, radii)):
                raise ValueError("sphere radii must be finite and positive")
            available = data["particles/force_data_available"]
            if any(v not in (0, 1) for v in available):
                raise ValueError("force availability must contain only 0 or 1")
            po, wo = data["contacts/particle_offsets"], data["contacts/wall_offsets"]
            for offsets, count in ((po, nparticle), (wo, nwall)):
                if (offsets[0] != 0 or offsets[-1] != count
                        or any(a > b for a, b in zip(offsets, offsets[1:]))):
                    raise ValueError("invalid contact offsets")
                if any(not available[i] and offsets[i] != offsets[i + 1] for i in range(n)):
                    raise ValueError("unavailable force output has contact records")
            for name, values in data.items():
                if name.startswith("contacts/") and expected[name][0] == "<f8":
                    if not all(math.isfinite(v) for v in values):
                        raise ValueError("nonfinite contact data: {}".format(name))
            for name in ("particles/fluid_force", "particles/fluid_torque"):
                if any(available[i] and not all(math.isfinite(v) for v in data[name][3*i:3*i+3])
                       for i in range(n)):
                    raise ValueError("available particle has nonfinite {}".format(name))
            lookup = {body_id: i for i, body_id in enumerate(body_ids)}
            partner_ids = data["contacts/particle_partner_id"]
            if any(body_id not in lookup for body_id in partner_ids):
                raise ValueError("contact partner ID is absent from snapshot")
            partner_rows = [lookup[body_id] for body_id in partner_ids]
            result["force_schema"] = True
            result["force_available"] = [bool(v) for v in available]
            metadata["force_unavailable_count"] = available.count(0)
            if metadata["force_unavailable_count"]:
                warnings.append("{} particles have unavailable stored force data".format(
                    metadata["force_unavailable_count"]))

            def vector(name, row):
                return data[name][3 * row:3 * row + 3]

            def contact_force(kind, row):
                normal = vector("contacts/{}_normal_force".format(kind), row)
                tangent = vector("contacts/{}_tangential_force".format(kind), row)
                return tuple(a + b for a, b in zip(normal, tangent))

            def cross(a, b):
                return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
                        a[0]*b[1] - a[1]*b[0])

            def residual(vectors):
                denominator = math.fsum(math.hypot(*v) for v in vectors)
                if denominator == 0:
                    return float("nan")
                total = tuple(math.fsum(v[axis] for v in vectors) for axis in range(3))
                value = math.hypot(*total) / denominator
                if not math.isfinite(value):
                    raise ValueError("nonfinite residual from stored loads")
                return value

            for i in range(n):
                if not available[i]:
                    continue
                forces = [vector("particles/fluid_force", i)]
                torques = [vector("particles/fluid_torque", i)]
                torque_geometry_available = geometry[i] == 0
                active = 0
                for row in range(po[i], po[i + 1]):
                    j = partner_rows[row]
                    if i == j:
                        raise ValueError("particle cannot contact itself")
                    force = contact_force("particle", row)
                    forces.append(force)
                    active += math.hypot(*force) > contact_force_threshold
                    if geometry[i] != 0 or geometry[j] != 0:
                        torque_geometry_available = False
                        continue
                    delta = tuple(points[j][axis] - points[i][axis] for axis in range(3))
                    distance_sq = math.fsum(v * v for v in delta)
                    if distance_sq <= 0 or not math.isfinite(distance_sq):
                        raise ValueError("coincident or invalid sphere contact centers")
                    alpha = (distance_sq + radii[i]**2 - radii[j]**2) / (2 * distance_sq)
                    # Match the solver's spherical contact plane, also for unequal radii.
                    contact_point = tuple(points[i][axis] + alpha * delta[axis] for axis in range(3))
                    lever = tuple(contact_point[axis] - points[i][axis] for axis in range(3))
                    torques.append(cross(lever, force))
                for row in range(wo[i], wo[i + 1]):
                    force = contact_force("wall", row)
                    forces.append(force)
                    active += math.hypot(*force) > contact_force_threshold
                    point = vector("contacts/wall_contact_point", row)
                    lever = tuple(point[axis] - points[i][axis] for axis in range(3))
                    # The normal-force moment is retained at the actual stored point.
                    torques.append(cross(lever, force))
                result["active_contact_count"][i] = active
                result["contacted"][i] = active > 0
                eta_force = residual(forces)
                result["eta_force"][i] = eta_force
                if math.isnan(eta_force):
                    metadata["zero_force_denominator_count"] += 1
                if torque_geometry_available:
                    eta_torque = residual(torques)
                    result["eta_torque"][i] = eta_torque
                    if math.isnan(eta_torque):
                        metadata["zero_torque_denominator_count"] += 1
                else:
                    metadata["torque_geometry_unavailable_count"] += 1
            for key, reason in (
                    ("torque_geometry_unavailable_count", "STL torque reference or non-spherical contact point unavailable"),
                    ("zero_force_denominator_count", "zero force normalisation denominators"),
                    ("zero_torque_denominator_count", "zero torque normalisation denominators")):
                if metadata[key]:
                    warnings.append("{} particles: {}".format(metadata[key], reason))
            return result
    except (OSError, ValueError, KeyError, TypeError, OverflowError, struct.error) as exc:
        raise ValueError("{}: {}".format(path, exc)) from exc


from collections import defaultdict
from itertools import combinations, product
import math


MOTIFS = ("4A", "5A", "6Z", "6A", "7A")


def vector(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def angle_within(a, b, low, high):
    """Strict angular bounds in degrees; a zero-length vector never qualifies."""
    norm_sq = dot(a, a) * dot(b, b)
    if norm_sq <= 0.0:
        return False
    cosine = dot(a, b) / math.sqrt(norm_sq)
    return math.cos(math.radians(high)) < cosine < math.cos(math.radians(low))


def bond_graph(points, cutoff):
    """Nonperiodic, strict center-distance cutoff; each pair is visited once."""
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("bond cutoff must be finite and positive")
    adjacency = [set() for _ in points]
    cells = {}
    cutoff_sq = cutoff * cutoff
    offsets = tuple(product((-1, 0, 1), repeat=3))
    for i, p in enumerate(points):
        if not all(math.isfinite(x) for x in p):
            raise ValueError("nonfinite particle position")
        key = tuple(math.floor(x / cutoff) for x in p)
        for offset in offsets:
            other = tuple(key[a] + offset[a] for a in range(3))
            for j in cells.get(other, ()):
                d = vector(p, points[j])
                if dot(d, d) < cutoff_sq:
                    adjacency[i].add(j)
                    adjacency[j].add(i)
        cells.setdefault(key, []).append(i)
    return adjacency


def regular_tetra(tetra, points):
    for center in tetra:
        arms = [vector(points[j], points[center]) for j in tetra if j != center]
        if any(not angle_within(a, b, 50.0, 70.0) for a, b in combinations(arms, 2)):
            return False
    return True


def opposite_sides(face, first, second, points):
    a, b, c = (points[i] for i in face)
    normal = cross(vector(b, a), vector(c, a))
    side1 = dot(normal, vector(points[first], a))
    side2 = dot(normal, vector(points[second], a))
    return (side1 < 0.0 < side2) or (side2 < 0.0 < side1)


def regular_ring(ring, points):
    """Paper SI edge angles, centroid angles and adjacent local plane normals."""
    count = len(ring)
    center = tuple(sum(points[i][a] for i in ring) / count for a in range(3))
    edge_target = 108.0 if count == 5 else 90.0
    rotation_target = 72.0 if count == 5 else 90.0
    normals = []
    for k, i in enumerate(ring):
        previous, following = ring[k - 1], ring[(k + 1) % count]
        left = vector(points[previous], points[i])
        right = vector(points[following], points[i])
        if not angle_within(left, right, edge_target - 12.0, edge_target + 12.0):
            return False
        if not angle_within(vector(points[i], center), vector(points[following], center),
                            rotation_target - 10.0, rotation_target + 10.0):
            return False
        normal = cross(left, right)
        norm = math.sqrt(dot(normal, normal))
        if norm == 0.0:
            return False
        normals.append(tuple(x / norm for x in normal))
    threshold = math.cos(math.radians(10.0))
    return all(dot(normals[k - 1], normal) > threshold for k, normal in enumerate(normals))


def octahedron_angles(ring, poles, points):
    members = tuple(ring) + tuple(poles)
    center = tuple(sum(points[i][a] for i in members) / 6.0 for a in range(3))
    vectors = {i: vector(points[i], center) for i in members}
    pairs = [(ring[k], ring[(k + 1) % 4]) for k in range(4)]
    pairs.extend(product(ring, poles))
    return all(angle_within(vectors[a], vectors[b], 80.0, 100.0) for a, b in pairs)


def induced_cycles(vertices, size, adjacency):
    """Enumerate each chordless 4/5 cycle once, with a canonical orientation."""
    allowed = set(vertices)
    for start in sorted(allowed):
        def extend(path):
            if len(path) == size:
                if path[1] < path[-1]:
                    yield tuple(path)
                return
            closing = len(path) == size - 1
            for nxt in sorted(adjacency[path[-1]] & allowed):
                if nxt <= start or nxt in path:
                    continue
                if closing and start not in adjacency[nxt]:
                    continue
                if any(nxt in adjacency[v] and not (closing and v == start) for v in path[:-1]):
                    continue
                yield from extend(path + [nxt])
        yield from extend([start])


def classify_clusters(points, cutoff):
    """Return all distinct member-index tuples, independently for each motif.

    External bonds do not invalidate a motif. Nested motifs are retained.
    5A and 6Z use the explicitly chosen induced-core completion documented
    in the module docstring, plus regular tetrahedra on opposite face sides.
    """
    result = {name: set() for name in MOTIFS}
    if len(points) < 4:
        return result
    adjacency = bond_graph(points, cutoff)

    # Ordered clique enumeration finds each candidate tetra exactly once.
    for a, neighbors in enumerate(adjacency):
        for b in sorted(v for v in neighbors if v > a):
            shared = neighbors & adjacency[b]
            for c in sorted(v for v in shared if v > b):
                for d in sorted(v for v in shared & adjacency[c] if v > c):
                    tetra = (a, b, c, d)
                    if regular_tetra(tetra, points):
                        result["4A"].add(tetra)

    faces = defaultdict(list)
    fans = defaultdict(lambda: defaultdict(set))
    for tetra in sorted(result["4A"]):
        for face in combinations(tetra, 3):
            cap = next(i for i in tetra if i not in face)
            faces[face].append(cap)
        for axis in combinations(tetra, 2):
            a, b = (i for i in tetra if i not in axis)
            fans[axis][a].add(b)
            fans[axis][b].add(a)

    # Five particles, two tetrahedra and nine internal bonds (K5 minus one).
    for face, caps in faces.items():
        for a, b in combinations(caps, 2):
            if b not in adjacency[a] and opposite_sides(face, a, b, points):
                result["5A"].add(tuple(sorted(face + (a, b))))

    # Six particles: common bonded axis joined to an induced four-node path.
    # Each path edge has already passed the full tetrahedral angle test.
    for axis, fan in fans.items():
        for a in sorted(fan):
            for b in sorted(fan[a]):
                for c in sorted(fan[b] - {a}):
                    if c in adjacency[a] or not opposite_sides(axis + (b,), a, c, points):
                        continue
                    for d in sorted(fan[c] - {a, b}):
                        if d <= a or d in adjacency[a] or d in adjacency[b]:
                            continue
                        if opposite_sides(axis + (c,), b, d, points):
                            result["6Z"].add(tuple(sorted(axis + (a, b, c, d))))

    # Candidate poles share ring vertices. Enumerate two-hop pole candidates
    # locally instead of all N(N-1)/2 pairs or all N-choose-7 particle subsets.
    ring_geometry = {}
    for u in range(len(points)):
        if len(adjacency[u]) < 4:
            continue
        candidates = set()
        for neighbor in adjacency[u]:
            candidates.update(adjacency[neighbor])
        for v in sorted(x for x in candidates if x > u):
            bonded_poles = v in adjacency[u]
            size = 5 if bonded_poles else 4
            common = adjacency[u] & adjacency[v]
            if len(common) < size:
                continue
            motif = "7A" if bonded_poles else "6A"
            for ring in induced_cycles(common, size, adjacency):
                members = tuple(sorted(ring + (u, v)))
                if members in result[motif]:
                    continue
                if ring not in ring_geometry:
                    ring_geometry[ring] = regular_ring(ring, points)
                if not ring_geometry[ring]:
                    continue
                if motif == "6A" and not octahedron_angles(ring, (u, v), points):
                    continue
                result[motif].add(members)
    return result




def read_snapshot(path, bounds, bound_eps):
    """Compatibility API: selected IDs, centers, volumes, and archive count."""
    bounds = validated_bounds(bounds, "bounds")
    if not math.isfinite(bound_eps) or bound_eps < 0:
        raise ValueError("boundary tolerance must be finite and nonnegative")
    snapshot = read_analysis_snapshot(path)
    selected = [i for i, p in enumerate(snapshot["points"])
                if inside_bounds(p, bounds, bound_eps)]
    return ([snapshot["body_ids"][i] for i in selected],
            [snapshot["points"][i] for i in selected],
            [snapshot["volumes"][i] for i in selected], len(snapshot["points"]))


import argparse
import statistics
import tempfile


GROUPS = MOTIFS + ("unclassified",)


def inside_bounds(point, bounds, eps=BOUND_EPS):
    return all(bounds[2*a] - eps <= point[a] <= bounds[2*a+1] + eps
               for a in range(3))


def exclusive_particle_labels(count, clusters):
    """One label per particle; later MOTIFS entries have higher priority.

    Keep overlapping cluster discovery intact: 5A/6Z construction needs 4A.
    This resolves particle membership, not ownership of entire clusters.
    """
    labels = ["unclassified"] * count
    for motif in MOTIFS:
        for cluster in clusters[motif]:
            for i in cluster:
                labels[i] = motif
    return labels


def contacted_group_statistics(snapshot, selected, labels, local_bounds=None):
    """Equal-weight means of individual ratios, without relabeling the ROI.

    Forces have already been summed using ALL owner contacts. The contact
    threshold only selects particles; it does not discard small force terms.
    N_F/N_T are each mean's actual denominator, which may differ if torque
    geometry is unavailable or its normalization denominator is zero.
    """
    values = {g: {"F": [], "T": []} for g in GROUPS}
    for i, group in zip(selected, labels):
        if local_bounds is not None and not inside_bounds(snapshot["points"][i], local_bounds):
            continue
        if not snapshot["force_available"][i] or not snapshot["contacted"][i]:
            continue
        for field, key in (("F", "eta_force"), ("T", "eta_torque")):
            value = snapshot[key][i]
            if math.isfinite(value):
                values[group][field].append(value)
    result = {}
    for g in GROUPS:
        result[g] = {}
        for field in ("F", "T"):
            samples = values[g][field]
            result[g]["N_" + field] = len(samples)
            result[g]["mean_eta" + field] = (
                math.fsum(samples)/len(samples) if samples else float("nan"))
    return result


def common_header(stream, metadata):
    stream.write("# Exclusive particle groups: 7A > 6A > 6Z > 5A > 4A > unclassified\n")
    stream.write("# TCC classification uses all centers in upper_bounds_m, then labels are held fixed for pore statistics.\n")
    stream.write("# Nonperiodic; all motif members must lie within upper_bounds_m; static particles included.\n")
    for name in ("upper_bounds_m", "pore_bounds_m"):
        stream.write("# {} xmin xmax ymin ymax zmin zmax: {}\n".format(
            name, " ".join(format(x, ".17g") for x in metadata[name])))
    stream.write("# inclusive_boundary_tolerance_m {:.17g}\n".format(BOUND_EPS))
    stream.write("# bond_cutoff_m {}; source {}\n".format(
        metadata["cutoff"], metadata["cutoff_source"]))
    stream.write("# reference_diameter_m {}; reference_time_s {}; bond_factor {}\n".format(
        metadata["diameter"], metadata["reference_time"], metadata["bond_factor"]))
    stream.write("# Single fixed cutoff for all times; r_ij < cutoff; no contact requirement for TCC bonds.\n")
    stream.write("# Time in seconds; lengths in m; frac is a fraction, not a percentage.\n")


def write_count_table(stream, rows, metadata):
    common_header(stream, metadata)
    stream.write("# Counts below are PARTICLES, not the old overlapping cluster counts.\n")
    stream.write("# sum(4A,5A,6Z,6A,7A,unclassified)=total (upper reservoir); frac_g=N_g/total.\n")
    stream.write("# Empty upper reservoir: total=0, all counts=0, fractions=nan.\n")
    columns = ["time"] + list(GROUPS) + ["total"] + ["frac_" + g for g in GROUPS]
    stream.write("# " + "\t".join(columns) + "\n")
    for row in rows:
        total = row["total"]
        fields = [row["time"]] + [str(row["counts"][g]) for g in GROUPS] + [str(total)]
        fields += [format(row["counts"][g]/total, ".17g") if total else "nan" for g in GROUPS]
        stream.write("\t".join(fields) + "\n")


def write_residual_table(stream, rows, metadata):
    common_header(stream, metadata)
    stream.write("# Contacted-only sample: >=1 recorded particle OR wall contact with |Fn+Ft| > {:.17g} N.\n".format(metadata["contact_force_eps"]))
    stream.write("# Contact may be outside the motif or either spatial region; all contacts enter each particle's residual.\n")
    stream.write("# Fc=Fn+Ft; tau_c=(contact_point-owner_center) cross Fc, including full wall normal-force moment.\n")
    stream.write("# etaF=|FH+sum(Fc)|/(|FH|+sum(|Fc|)); etaT=|TH+sum(tau_c)|/(|TH|+sum(|tau_c|)).\n")
    stream.write("# mean_etaF/mean_etaT: arithmetic means of individual dimensionless ratios.\n")
    stream.write("# N_F/N_T: actual numbers of particles used for the corresponding mean, not total TCC-group populations.\n")
    stream.write("# Empty sample / zero denominator / unavailable force data: excluded from corresponding N; mean=nan if N=0.\n")
    stream.write("# STL torque references and non-sphere particle contact points are unavailable in v1: affected etaT=nan, etaF remains usable.\n")
    stream.write("# Recorded loads only: no separate gravity or adhesion. Total-balance interpretation assumes these omitted loads are zero.\n")
    stream.write("# pore is a center-selected rectangular volume, not a first-hit surface selection.\n")
    columns = ["time"]
    for region in ("upper", "pore"):
        for field in ("mean_etaF", "mean_etaT", "N_F", "N_T"):
            columns += [region + "_" + field + "_" + g for g in GROUPS]
    stream.write("# " + "\t".join(columns) + "\n")
    for row in rows:
        fields = [row["time"]]
        for region in ("upper", "pore"):
            for field in ("mean_etaF", "mean_etaT", "N_F", "N_T"):
                for g in GROUPS:
                    value = row[region][g][field]
                    fields.append(str(value) if field.startswith("N_") else format(value, ".17g"))
        stream.write("\t".join(fields) + "\n")


def save_tables(count_path, residual_path, rows, metadata):
    """Stage both complete files after every snapshot passes, then replace.

    A corrupt later input leaves both previous results untouched. Each file
    replacement is atomic; simultaneous multi-file atomicity is not claimed.
    """
    staged = []
    try:
        for path, writer in ((count_path, write_count_table), (residual_path, write_residual_table)):
            path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                             dir=path.parent, prefix=path.name+".", suffix=".tmp",
                                             delete=False) as stream:
                temporary = Path(stream.name)
                staged.append((temporary, path))
                writer(stream, rows, metadata)
                stream.flush()
                os.fsync(stream.fileno())
        for temporary, path in staged:
            os.replace(temporary, path)
    finally:
        for temporary, _ in staged:
            if temporary.exists():
                temporary.unlink()


def validated_bounds(bounds, name):
    if (len(bounds) != 6 or not all(math.isfinite(x) for x in bounds)
            or any(bounds[2*a] >= bounds[2*a+1] for a in range(3))):
        raise ValueError(name + " must contain finite xmin<xmax, ymin<ymax, zmin<zmax")
    return tuple(bounds)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Exclusive upper-reservoir TCC particle counts/fractions and contacted-particle residuals at ALL .bin times.",
        epilog="Run python3 analyze_tcc.py from the case; no third-party packages or companion scripts required.")
    parser.add_argument("--case", type=Path, default=Path(CASE_DIRECTORY))
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_FILE), help="particle count/fraction file, relative to case")
    parser.add_argument("--residual-output", type=Path, default=Path(RESIDUAL_OUTPUT_FILE), help="contacted-only residual file, relative to case")
    parser.add_argument("--bond-cutoff", type=float, default=BOND_CUTOFF_M, help="fixed center-distance cutoff in m; overrides factor/diameter")
    parser.add_argument("--bond-factor", type=float, default=BOND_CUTOFF_FACTOR, help="cutoff/reference-diameter ratio (default: %(default)s)")
    parser.add_argument("--diameter", type=float, default=REFERENCE_DIAMETER_M, help="reference diameter in m; otherwise first populated upper snapshot median")
    parser.add_argument("--bounds", type=float, nargs=6, default=BOUNDS_M, metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"), help="whole upper-reservoir classification/analysis bounds in m")
    parser.add_argument("--pore-bounds", type=float, nargs=6, default=PORE_BOUNDS_M, metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"), help="local pore prism bounds in m; must lie within --bounds")
    parser.add_argument("--contact-force-eps", type=float, default=CONTACT_FORCE_EPS_N, help="strict minimum |Fn+Ft| in N for sample selection only")
    args = parser.parse_args(argv)
    try:
        for name, value in (("bond cutoff", args.bond_cutoff), ("bond factor", args.bond_factor), ("reference diameter", args.diameter)):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(name + " must be finite and positive")
        if not math.isfinite(BOUND_EPS) or BOUND_EPS < 0:
            raise ValueError("BOUND_EPS must be finite and nonnegative")
        if not math.isfinite(args.contact_force_eps) or args.contact_force_eps < 0:
            raise ValueError("contact-force-eps must be finite and nonnegative")
        upper_bounds = validated_bounds(args.bounds, "upper bounds")
        pore_bounds = validated_bounds(args.pore_bounds, "pore bounds")
        if any(pore_bounds[2*a] < upper_bounds[2*a] - BOUND_EPS or
               pore_bounds[2*a+1] > upper_bounds[2*a+1] + BOUND_EPS for a in range(3)):
            raise ValueError("pore bounds must lie within upper bounds; all coordinates are in m (default pore y=0.025--0.035)")
        case = args.case.resolve()
        data = case / "particleData"
        snapshots = time_snapshots(data)
        output = (args.output if args.output.is_absolute() else case / args.output).resolve()
        residual_output = (args.residual_output if args.residual_output.is_absolute() else case / args.residual_output).resolve()
        if output == residual_output:
            raise ValueError("count and residual output files must differ")
        for path in (output, residual_output):
            if path == data.resolve() or data.resolve() in path.parents:
                raise ValueError("outputs must be outside particleData to preserve input archives")
            if path == Path(__file__).resolve():
                raise ValueError("output cannot overwrite this script")
        cutoff, diameter = args.bond_cutoff, args.diameter
        source = "fixed_explicit_distance" if cutoff is not None else "factor_times_reference_diameter"
        reference_time = "explicit" if diameter is not None else "none"
        if cutoff is not None:
            diameter = None
            reference_time = "not_used"
        elif diameter is not None:
            cutoff = args.bond_factor * diameter
        if cutoff is not None and (not math.isfinite(cutoff) or cutoff <= 0):
            raise ValueError("computed bond cutoff must be finite and positive")
        print("Analyzing {} saved times; exclusive particle groups; upper bounds {} m; pore bounds {} m.".format(
            len(snapshots), upper_bounds, pore_bounds), flush=True)
        rows = []
        for time, path in snapshots:
            snapshot = read_analysis_snapshot(path, contact_force_threshold=args.contact_force_eps)
            selected = [i for i, p in enumerate(snapshot["points"]) if inside_bounds(p, upper_bounds)]
            points = [snapshot["points"][i] for i in selected]
            if cutoff is None and points:
                diameter = statistics.median((6.0*snapshot["volumes"][i]/math.pi)**(1.0/3.0) for i in selected)
                cutoff = args.bond_factor * diameter
                if not math.isfinite(cutoff) or cutoff <= 0:
                    raise ValueError("computed bond cutoff must be finite and positive")
                reference_time = time
            clusters = classify_clusters(points, cutoff) if cutoff is not None else {g: set() for g in MOTIFS}
            labels = exclusive_particle_labels(len(selected), clusters)
            counts = {g: labels.count(g) for g in GROUPS}
            if sum(counts.values()) != len(selected):
                raise ValueError("internal error: particle-group counts do not sum to upper particle count")
            row = {"time": time, "total": len(selected), "counts": counts,
                   "upper": contacted_group_statistics(snapshot, selected, labels),
                   "pore": contacted_group_statistics(snapshot, selected, labels, pore_bounds)}
            rows.append(row)
            n_pore = sum(inside_bounds(p, pore_bounds) for p in points)
            n_force = sum(row["upper"][g]["N_F"] for g in GROUPS)
            n_torque = sum(row["upper"][g]["N_T"] for g in GROUPS)
            print("t={} s: upper={}/{}, pore={}; contacted valid F/T={}/{}; {}".format(
                time, len(selected), len(snapshot["points"]), n_pore, n_force, n_torque,
                " ".join("{}={}".format(g, counts[g]) for g in GROUPS)), flush=True)
            n_missing = sum(not snapshot["force_available"][i] for i in selected)
            if n_missing:
                print("  {} upper particles have no evaluated force data; omitted from residual means only.".format(n_missing), flush=True)
            n_contact = sum(snapshot["force_available"][i] and snapshot["contacted"][i] for i in selected)
            if n_torque < n_contact:
                print("  {} contacted upper particles lack a defined torque ratio (zero denominator or missing contact geometry).".format(n_contact-n_torque), flush=True)
        metadata = dict(upper_bounds_m=upper_bounds, pore_bounds_m=pore_bounds,
                        cutoff=cutoff, cutoff_source=source if cutoff is not None else "empty_region",
                        diameter=diameter, reference_time=reference_time, bond_factor=args.bond_factor,
                        contact_force_eps=args.contact_force_eps)
        save_tables(output, residual_output, rows, metadata)
        print("Wrote {} and {} ({} times; fixed cutoff {} m).".format(output, residual_output, len(rows), cutoff), flush=True)
        return 0
    except (OSError, ValueError, TypeError, OverflowError, struct.error) as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
