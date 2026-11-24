#!/usr/bin/env python3
"""Standard TCC particle populations and contacted-particle force/torque residuals.

Run from the simulation case: python3 analyze_standard_tcc.py
Or supply --case /path/to/case. Python standard library only.

Reads every numeric-time particleData/*.bin, then writes one row per time to
postProcessing/tccStandardCounts.dat and postProcessing/tccStandardResiduals.dat.
Columns, units, sample selection, normalization, and undefined-value handling
match analyze_tcc.py (Test_CFD_DEM source blob 2ee0c09c6b93e7b2a0b35c5bd4d9cd52f8cf037d).

The single-species, nonperiodic modified Voronoi graph uses fc=0.82 by default.
The maximum candidate bond length is 1.41 times the median volume-equivalent
diameter in the first nonempty upper region. This fixed scale applies at all
times to spheres and STL centers alike; --bond-cutoff overrides it.
The five standard TCC motifs are identified by their rings and spindles
following Malins et al. (J. Chem. Phys. 139, 234506) and royallgroup/TCC
v1.0.4, revision d865986214254dfdd90038016990759246ea1ca7.
A spindle bonds to every particle in a ring. A three-member ring with one
spindle gives 4A; rings of three, four or five with two spindles give 5A,
6A or 7A. Two overlapping 5A form 6Z under the standard spindle/ring rules.
There is no separate face-angle or regularity filter.

Classify all centers inside --bounds, assign each particle once by priority
7A > 6A > 6Z > 5A > 4A > unclassified, then use those labels in --pore-bounds.
Counts/fractions describe particles. Force and torque means include only
particles with an evaluated contact above --contact-force-eps. Every contact
enters that particle's residual, including contacts outside either region.
N_F and N_T are independent sample counts. Undefined or unrecorded values
are excluded from the corresponding mean; empty samples report nan.
As in analyze_tcc.py, the archived loads omit separate gravity/adhesion.
STL torque references and non-spherical pair contact points are unavailable
in v1, so affected torque residuals are undefined.
"""

# -------------------- EDITABLE DEFAULTS --------------------
CASE_DIRECTORY = "."
# xmin, xmax, ymin, ymax, zmin, zmax; every coordinate is in m.
BOUNDS_M = (-0.020, 0.020, 0.025, 0.075, -0.020, 0.020)
PORE_BOUNDS_M = (-0.005, 0.005, 0.025, 0.035, -0.005, 0.005)
BOUND_EPS = 1.0e-12
BOND_CUTOFF_M = None           # Explicit distance overrides factor/diameter.
BOND_CUTOFF_FACTOR = 1.41      # Maximum candidate distance in reference diameters.
VORONOI_FC = 0.82              # Modified Voronoi screening parameter.
REFERENCE_DIAMETER_M = None    # Infer once from first populated upper snapshot.
CONTACT_FORCE_EPS_N = 1.0e-15  # Particle eligibility only, not a force-term cutoff.
OUTPUT_FILE = "postProcessing/tccStandardCounts.dat"
RESIDUAL_OUTPUT_FILE = "postProcessing/tccStandardResiduals.dat"
# -----------------------------------------------------------

import argparse
from array import array
from collections import defaultdict
from decimal import Decimal, InvalidOperation
import hashlib
from itertools import combinations, product
import json
import math
import os
from pathlib import Path
import re
import statistics
import struct
import sys
import tempfile


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


MOTIFS = ("4A", "5A", "6Z", "6A", "7A")


def vector(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def _candidate_graph(points, cutoff):
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


def bond_graph(points, cutoff, fc=VORONOI_FC):
    """Build the nonperiodic, single-cutoff modified Voronoi graph (Malins et al.)."""
    if not math.isfinite(fc) or not 0.0 < fc <= 1.0:
        raise ValueError("modified Voronoi fc must be finite and in (0, 1]")
    candidates = _candidate_graph(points, cutoff)
    adjacency = [set() for _ in points]
    for i, neighbors in enumerate(candidates):
        # Increasing squared distance, then original particle order for ties.
        ordered = []
        for j in neighbors:
            delta = vector(points[i], points[j])
            ordered.append((dot(delta, delta), j, delta))
        ordered.sort()

        for rank, (_, j, ij) in enumerate(ordered):
            # Even a screened k can screen a farther j. The ratio comparison
            # uses the 1e-6 tolerance of the reference detector, independent of units.
            for nearer in range(rank):
                ik2, k, ik = ordered[nearer]
                jk = vector(points[j], points[k])
                denominator = ik2 + dot(jk, jk)
                # Evaluate Eq. 6 in the same order as the reference code,
                # including the strict comparison near the numerical tolerance.
                numerator = dot(ij, ik) - dot(ij, jk)
                if denominator and numerator / denominator - fc > 1e-6:
                    break
            else:
                # TCC adds the reverse edge if either endpoint retains its bond.
                adjacency[i].add(j)
                adjacency[j].add(i)
    return adjacency


def classify_clusters(points, cutoff, fc=VORONOI_FC):
    """Return unique sorted member-index tuples for the five standard TCC motifs."""
    adjacency = bond_graph(points, cutoff, fc)
    result = {name: set() for name in MOTIFS}

    # A triangle with exactly one common neighbor is sp3b (4A). Exactly
    # two gives sp3c (5A), independently of any 4A classification.
    triangles_by_edge = defaultdict(list)
    for a, neighbors in enumerate(adjacency):
        for b in neighbors:
            if b <= a:
                continue
            shared = neighbors & adjacency[b]
            for c in shared:
                if c <= b:
                    continue
                spindles = shared & adjacency[c]
                if len(spindles) == 1:
                    result["4A"].add(tuple(sorted((a, b, c, *spindles))))
                elif len(spindles) == 2:
                    ring = (a, b, c)
                    result["5A"].add(tuple(sorted((*ring, *spindles))))
                    for u, v, tip in ((a, b, c), (a, c, b), (b, c, a)):
                        triangles_by_edge[u, v].append((tip, spindles))

    # A 6Z consists of two sp3c rings sharing an edge. Each remaining
    # ring vertex is a spindle of the other ring, and these vertices bond.
    for common_edge, rings in triangles_by_edge.items():
        for (a, caps_a), (b, caps_b) in combinations(rings, 2):
            if b not in caps_a or a not in caps_b or b not in adjacency[a]:
                continue
            if caps_a.isdisjoint(caps_b):
                result["6Z"].add(tuple(sorted((*common_edge, *caps_a, *caps_b))))

    # For 4- and 5-member rings, chordlessness enforces shortest paths.
    # Search from candidate spindle pairs, then check every possible
    # spindle in the full graph. Spindle-spindle bonds are unrestricted.
    for u, neighbors in enumerate(adjacency):
        if len(neighbors) < 4:
            continue
        two_hop = set()
        for neighbor in neighbors:
            two_hop.update(adjacency[neighbor])
        for v in two_hop:
            if v <= u:
                continue
            common = neighbors & adjacency[v]
            for size, motif in ((4, "6A"), (5, "7A")):
                if len(common) < size:
                    continue
                for ring in induced_cycles(common, size, adjacency):
                    spindles = set.intersection(*(adjacency[n] for n in ring))
                    if len(spindles) == 2:
                        result[motif].add(tuple(sorted((*ring, u, v))))
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


GROUPS = MOTIFS + ("unclassified",)


def inside_bounds(point, bounds, eps=BOUND_EPS):
    return all(bounds[2*a] - eps <= point[a] <= bounds[2*a+1] + eps
               for a in range(3))


def exclusive_particle_labels(count, clusters):
    """One label per particle; later MOTIFS entries have higher priority.

    Motif discovery can yield overlapping particle sets.
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
    stream.write("# Standard TCC: modified Voronoi neighbors, fc={:.17g}; no motif angle filters.\n".format(metadata["fc"]))
    stream.write("# Candidate neighbors: r_ij < fixed cutoff at every time; no contact requirement for structural bonds.\n")
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
        description="Standard TCC: exclusive upper-reservoir particle counts/fractions and contacted-particle residuals at ALL .bin times.",
        epilog="Run python3 analyze_standard_tcc.py from the case; no third-party packages or companion scripts required.")
    parser.add_argument("--case", type=Path, default=Path(CASE_DIRECTORY))
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_FILE), help="particle count/fraction file, relative to case")
    parser.add_argument("--residual-output", type=Path, default=Path(RESIDUAL_OUTPUT_FILE), help="contacted-only residual file, relative to case")
    parser.add_argument("--fc", type=float, default=VORONOI_FC, help="modified Voronoi screening parameter in (0, 1] (default: %(default)s)")
    parser.add_argument("--bond-cutoff", type=float, default=BOND_CUTOFF_M, help="maximum Voronoi candidate distance in m; overrides factor/diameter")
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
        if not math.isfinite(args.fc) or not 0.0 < args.fc <= 1.0:
            raise ValueError("modified Voronoi fc must be finite and in (0, 1]")
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
        print("Standard TCC: modified Voronoi fc={}, nonperiodic.".format(args.fc), flush=True)
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
            clusters = classify_clusters(points, cutoff, args.fc) if cutoff is not None else {g: set() for g in MOTIFS}
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
        metadata = dict(fc=args.fc, upper_bounds_m=upper_bounds, pore_bounds_m=pore_bounds,
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
