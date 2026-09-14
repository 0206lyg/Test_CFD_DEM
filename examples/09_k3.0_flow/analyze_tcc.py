#!/usr/bin/env python3
"""Count symmetry-filtered TCC-family motifs in the upper reservoir.

RUN, from the case containing particleData/:
    python3 analyze_tcc.py

No third-party packages, companion scripts or command-line options are needed.
Every existing nonnegative numeric time is processed, including 0.15, 0.3, ... .
Output: postProcessing/tccCounts.dat
    # time 4A 5A 6Z 6A 7A
    0      ...
    0.15   ...

DEFAULTS (editable immediately below this docstring)
* Region and boundary tolerance match upper_reservoir_solid_fraction.py:
  x,z = [-0.020,0.020] m; y = [0.025,0.075] m, both ends included with
  1e-12 m tolerance. The whole upper reservoir is analyzed, without y-slicing.
  All members of a counted cluster must have their centers in this region.
  No periodic boundary wrapping is used. All stored particle types, including
  static particles, are included when their positions satisfy this selection.
* Bonds satisfy r_ij < fixed BOND_CUTOFF_M. If BOND_CUTOFF_M is None, the
  FIRST PASS uses the paper's factor 1.295 times a reference diameter.
  The reference is the median volume-equivalent diameter (6 V/pi)^(1/3) in
  the first populated upper-reservoir snapshot; it then stays fixed across
  ALL times. For spheres this equals their geometric diameter. For STL
  particles it is a volume-equivalent scale, not a surface-contact distance.
  This is a provisional starting cutoff, NOT a minimum measured from YOUR
  g(r). Set BOND_CUTOFF_M to your measured first minimum when available.
  No automatic first-minimum selection or cutoff change across time occurs.
* The reference diameter can instead be fixed with REFERENCE_DIAMETER_M.
  Mixed sizes still use one fixed pooled center-distance cutoff in this version.
* Centers are read from particles/position. The current compressor stores
  sphere centers and the arithmetic mean of unique STL vertices for STL
  particles (the latter is not necessarily the volume centroid).

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

COUNTING
Each distinct member set counts once within each type, regardless of how
many ring/axis decompositions find it. External bonds do not invalidate it.
Nested and overlapping clusters count independently. Thus an ideal 7A yields
4A=5, 5A=5, 6Z=5, 6A=0, 7A=1. Counts are numbers of CLUSTERS, not unique
particles, particle fractions, or exclusive highest-order visualization labels.
classify_clusters(points, cutoff) returns the full member-index sets so that
future membership export can map them directly to the stored body IDs.

The v1 .bin reader checks the index and consumed arrays, with no extraction
of original files or STL geometry. Inputs are unchanged. An incomplete or
corrupt selected archive aborts processing and preserves the previous output.

Sources:
Tsurusawa & Tanaka, Nature Physics (2023), Methods and SI sections I-III:
https://doi.org/10.1038/s41567-023-02063-x
Malins et al., J. Chem. Phys. 139, 234506 (2013), Sec. IV and Table II:
https://arxiv.org/abs/1307.5517
"""


# -------------------- EDITABLE DEFAULTS --------------------
CASE_DIRECTORY = "."
BOUNDS_M = (-0.020, 0.020, 0.025, 0.075, -0.020, 0.020)
BOUND_EPS = 1.0e-12
BOND_CUTOFF_M = None            # Set to your measured g(r) first minimum, m.
BOND_CUTOFF_FACTOR = 1.295      # Provisional paper value; used only if above is None.
REFERENCE_DIAMETER_M = None     # None: infer once from first populated snapshot.
OUTPUT_FILE = "postProcessing/tccCounts.dat"
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


_BIN_HEADER = struct.Struct("<16sI")
_BIN_FOOTER = struct.Struct("<16sQQ32s")
_BIN_HEADER_MAGIC = b"PDBIN-HEADER-v1!"
_BIN_FOOTER_MAGIC = b"PDBIN-FOOTER-v1!"
_BIN_FORMAT = "openHFDIB-DEM particleData binary"
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


def read_snapshot(path, bounds, bound_eps):
    """Read IDs, centers and volumes; return selected lists and archive count.

    Bounds are xmin,xmax,ymin,ymax,zmin,zmax in m. All three axes include both
    endpoints with bound_eps tolerance, matching the upper-reservoir script.
    The archive's row order and body IDs are retained. Position values are
    absolute coordinates, not coordinates relative to the observation window.
    Only the JSON index and these three datasets are read and checksummed;
    STL geometry and stored original files are not read or extracted.
    """
    path = Path(path)
    try:
        if (len(bounds) != 6 or not all(math.isfinite(x) for x in bounds)
                or any(bounds[2 * a] >= bounds[2 * a + 1] for a in range(3))
                or not math.isfinite(bound_eps) or bound_eps < 0):
            raise ValueError("invalid bounds or boundary tolerance")
        with path.open("rb") as stream:
            size = os.fstat(stream.fileno()).st_size
            if size < _BIN_HEADER.size + _BIN_FOOTER.size:
                raise ValueError("truncated binary archive")
            if _BIN_HEADER.unpack(stream.read(_BIN_HEADER.size)) != (_BIN_HEADER_MAGIC, 1):
                raise ValueError("unsupported particleData header/version")
            stream.seek(-_BIN_FOOTER.size, os.SEEK_END)
            magic, offset, length, digest = _BIN_FOOTER.unpack(stream.read(_BIN_FOOTER.size))
            if (magic != _BIN_FOOTER_MAGIC or offset < _BIN_HEADER.size
                    or offset + length != size - _BIN_FOOTER.size):
                raise ValueError("invalid index footer")
            stream.seek(offset)
            raw_index = stream.read(length)
            if len(raw_index) != length or hashlib.sha256(raw_index).digest() != digest:
                raise ValueError("index SHA256 mismatch")
            index = json.loads(raw_index)
            if not isinstance(index, dict):
                raise ValueError("index must be a JSON object")
            if (index.get("format") != _BIN_FORMAT
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

            position = datasets.get("particles/position")
            shape = position.get("shape") if isinstance(position, dict) else None
            if (not isinstance(shape, list) or len(shape) != 2 or shape[1] != 3
                    or any(type(x) is not int or x < 0 for x in shape)):
                raise ValueError("positions must have shape [N, 3]")
            total_count = shape[0]
            ranges = []

            def read_array(name, dtype, expected_shape, typecode):
                item = datasets.get(name)
                if not isinstance(item, dict):
                    raise ValueError("missing {} dataset".format(name))
                actual_shape = item.get("shape")
                if (item.get("dtype") != dtype or actual_shape != expected_shape
                        or any(type(x) is not int for x in actual_shape)):
                    raise ValueError("{} must have dtype {} and shape {}".format(
                        name, dtype, expected_shape))
                elements = math.prod(expected_shape)
                start, count = item.get("offset"), item.get("size")
                if (type(start) is not int or type(count) is not int
                        or start < _BIN_HEADER.size or count != 8 * elements
                        or start + count > offset):
                    raise ValueError("invalid {} dataset offset/size".format(name))
                if count and any(start < end and begin < start + count
                                 for begin, end in ranges):
                    raise ValueError("overlapping consumed datasets")
                if count:
                    ranges.append((start, start + count))
                stream.seek(start)
                payload = stream.read(count)
                if len(payload) != count or hashlib.sha256(payload).hexdigest() != item.get("sha256"):
                    raise ValueError("{} dataset SHA256 mismatch".format(name))
                values = array(typecode)
                if values.itemsize != 8:
                    raise ValueError("this Python platform does not provide 8-byte numeric arrays")
                values.frombytes(payload)
                if sys.byteorder != "little":
                    values.byteswap()
                return values

            body_ids = read_array("particles/body_id", "<i8", [total_count], "q")
            positions = read_array("particles/position", "<f8", [total_count, 3], "d")
            volumes = read_array("particles/volume", "<f8", [total_count], "d")

        if len(set(body_ids)) != total_count:
            raise ValueError("body_id dataset contains duplicate IDs")
        if not all(math.isfinite(x) for x in positions):
            raise ValueError("position dataset contains nonfinite coordinates")
        if not all(math.isfinite(x) and x > 0 for x in volumes):
            raise ValueError("volume dataset contains nonpositive or nonfinite values")
        kept_ids, kept_points, kept_volumes = [], [], []
        for i, body_id in enumerate(body_ids):
            point = tuple(positions[3 * i:3 * i + 3])
            if all(bounds[2 * a] - bound_eps <= point[a] <= bounds[2 * a + 1] + bound_eps
                   for a in range(3)):
                kept_ids.append(body_id)
                kept_points.append(point)
                kept_volumes.append(volumes[i])
        return kept_ids, kept_points, kept_volumes, total_count
    except (OSError, ValueError, TypeError, struct.error, OverflowError) as exc:
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


import argparse
import statistics
import tempfile


def save_counts(path, rows, cutoff, cutoff_source, diameter, bounds, reference_time):
    """Replace the complete result only after every selected input succeeds."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                         dir=path.parent, prefix=path.name + ".", suffix=".tmp",
                                         delete=False) as stream:
            temporary = Path(stream.name)
            stream.write("# Symmetry-filtered TCC-family cluster counts; time in seconds\n")
            stream.write("# region upper reservoir; all members must be inside; nonperiodic\n")
            stream.write("# bounds_m xmin xmax ymin ymax zmin zmax: " + " ".join(format(x, ".17g") for x in bounds) + "\n")
            stream.write("# inclusive_boundary_tolerance_m {:.17g}\n".format(BOUND_EPS))
            stream.write("# bond_cutoff_m {}\n".format("not_needed_empty_region" if cutoff is None else format(cutoff, ".17g")))
            stream.write("# cutoff_source {}\n".format(cutoff_source))
            if diameter is not None:
                stream.write("# reference_diameter_m {:.17g}; reference_time_s {}\n".format(diameter, reference_time))
            stream.write("# single fixed cutoff for all times; no first-minimum estimation in this script\n")
            stream.write("# 4A all_12_vertex_angles=(50,70)_deg\n")
            stream.write("# 5A two_4A_shared_face; opposite_apices; 9_internal_bonds\n")
            stream.write("# 6Z three_4A_common_bond; induced_4_vertex_path; 12_internal_bonds\n")
            stream.write("# 6A induced_square; unbonded_poles; SI_ring_and_12_center_angle_tests\n")
            stream.write("# 7A induced_pentagon; bonded_poles; SI_ring_angle_and_planarity_tests\n")
            stream.write("# unique_member_sets_per_type; overlaps_and_nested_motifs_included\n")
            stream.write("# selected_time_s " + " ".join(row[0] for row in rows) + "\n")
            stream.write("# upper_particle_count " + " ".join(str(row[2]) for row in rows) + "\n")
            stream.write("# source_particle_count " + " ".join(str(row[3]) for row in rows) + "\n")
            stream.write("# time " + " ".join(MOTIFS) + "\n")
            for time, counts, _, _ in rows:
                stream.write(time + "\t" + "\t".join(str(counts[motif]) for motif in MOTIFS) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Count upper-reservoir TCC-family motifs at ALL .bin times; no options required.",
        epilog="Default invocation: python3 analyze_tcc.py. All settings can also be edited at the top of this file.")
    parser.add_argument("--case", type=Path, default=Path(CASE_DIRECTORY),
                        help="case containing particleData (default: current directory)")
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_FILE),
                        help="output path, relative to the case unless absolute")
    parser.add_argument("--bond-cutoff", type=float, default=BOND_CUTOFF_M,
                        help="fixed center-distance cutoff, m; overrides factor/diameter")
    parser.add_argument("--bond-factor", type=float, default=BOND_CUTOFF_FACTOR,
                        help="provisional cutoff/reference-diameter ratio (default: 1.295)")
    parser.add_argument("--diameter", type=float, default=REFERENCE_DIAMETER_M,
                        help="reference diameter, m (default: first populated snapshot median)")
    parser.add_argument("--bounds", type=float, nargs=6, default=BOUNDS_M,
                        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
                        help="override the upper-reservoir center-selection bounds, m")
    args = parser.parse_args(argv)
    try:
        for name, value in (("bond cutoff", args.bond_cutoff), ("bond factor", args.bond_factor),
                            ("reference diameter", args.diameter)):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError("{} must be finite and positive".format(name))
        if (not all(math.isfinite(v) for v in args.bounds)
                or any(args.bounds[2*a] >= args.bounds[2*a+1] for a in range(3))):
            raise ValueError("bounds must be finite, with each maximum greater than its minimum")
        case = args.case.resolve()
        data = case / "particleData"
        snapshots = time_snapshots(data)
        output = (args.output if args.output.is_absolute() else case / args.output).resolve()
        if output == data.resolve() or data.resolve() in output.parents:
            raise ValueError("output must be outside particleData to preserve the input archives")
        cutoff, diameter = args.bond_cutoff, args.diameter
        source = "fixed_explicit_distance" if cutoff is not None else "provisional_paper_factor_times_reference_diameter"
        reference_time = "explicit" if diameter is not None else "none"
        if cutoff is not None:
            diameter = None
        elif diameter is not None:
            cutoff = args.bond_factor * diameter
        if cutoff is not None and (not math.isfinite(cutoff) or cutoff <= 0):
            raise ValueError("computed bond cutoff must be finite and positive")
        print("Analyzing all {} saved times in upper reservoir.".format(len(snapshots)), flush=True)
        if cutoff is not None:
            print("Fixed bond cutoff: {:.12g} m ({})".format(cutoff, source), flush=True)
        rows = []
        for time, path in snapshots:
            body_ids, points, volumes, total = read_snapshot(path, args.bounds, BOUND_EPS)
            if cutoff is None and points:
                diameter = statistics.median((6.0 * volume / math.pi)**(1.0/3.0) for volume in volumes)
                cutoff = args.bond_factor * diameter
                if not math.isfinite(cutoff) or cutoff <= 0:
                    raise ValueError("computed bond cutoff must be finite and positive")
                reference_time = time
                print("Provisional fixed cutoff: {:.12g} m = {:.8g} x {:.12g} m; reference t={} s".format(
                    cutoff, args.bond_factor, diameter, time), flush=True)
            print("t={} s: {}/{} particles; classifying...".format(time, len(points), total), flush=True)
            clusters = classify_clusters(points, cutoff) if cutoff is not None else {name: set() for name in MOTIFS}
            counts = {name: len(clusters[name]) for name in MOTIFS}
            rows.append((time, counts, len(points), total))
            print("  " + "  ".join("{}={}".format(name, counts[name]) for name in MOTIFS), flush=True)
        if cutoff is None:
            source = "not_needed_empty_region"
        save_counts(output, rows, cutoff, source, diameter, args.bounds, reference_time)
        print("Wrote {} ({} times; counts of distinct clusters)".format(output, len(rows)), flush=True)
        return 0
    except (OSError, ValueError, TypeError, OverflowError) as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
