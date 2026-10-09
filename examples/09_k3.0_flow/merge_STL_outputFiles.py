#!/usr/bin/env python3
"""Export particleData/<time>.bin to ParaView with bodyId and TCCType.

RUN in the case directory:
    python3 merge_STL_outputFiles.py
Optional explicit center-distance cutoff, in metres:
    python3 merge_STL_outputFiles.py --bond-cutoff 0.0039
Repair an existing state without recalculating TCC or rewriting data:
    python3 merge_STL_outputFiles.py --state-only

Standalone Python 3.9+; no NumPy, VTK, analyze_tcc.py or compressor import is
required. Input is binary-v1 from compress_bodiesInfo.py. Original archives
are read only; neither extraction nor another simulation is needed.

SCOPE AND CLASSIFICATION
* ALL stored particles at every nonnegative numeric time are analyzed and
  visualized, including static particles. There is NO spatial/region filter.
* The geometric/topological predicates are copied unchanged from the supplied
  analyze_tcc.py: symmetry-filtered 4A, 5A, 6Z, 6A and 7A. No periodic wrapping.
  4A requires all 12 bond angles in (50,70) deg; 5A uses opposite face-sharing
  tetrahedra; 6Z uses three tetrahedra around a common bond and an induced
  four-vertex path. 6A and 7A use the same ring/planarity tests; the 6A poles
  are unbonded and the 7A poles are bonded. This retains that script's explicit
  motif definitions, rather than changing them to a different TCC variant.
* TCCType is one exclusive VISUALIZATION label per particle:
      0=None, 1=4A, 2=5A, 3=6Z, 4=6A, 5=7A
  Highest priority wins: 4A < 5A < 6Z < 6A < 7A. All memberships are determined
  before resolving overlaps; only the final label is exported.
* Default cutoff: 1.295 times the median volume-equivalent diameter in the
  first nonempty FULL snapshot, then fixed for the entire series. As in
  analyze_tcc.py, this is a provisional value, NOT an estimated g(r) minimum.
  Use --bond-cutoff (or BOND_CUTOFF_M below) for your measured minimum.
* Stored centers are used exactly: sphere centers, and unique-vertex mean
  positions for STL particles. For STL this is not necessarily a volume
  centroid. A single center-distance cutoff is used for all shapes/sizes.
* bodyId is the stored int64 body ID, never a row/triangle index. No ID
  reconstruction, reuse detection or correction across restarts is attempted.

OUTPUT in STLMerged/ (one merged file per time and represented geometry type):
    STL_Results0001.vtp, ...   actual STL surfaces; Cell Data: bodyId, TCCType
    STL_Results.pvd            surfaces with the original physical time values
    Sphere_Results0001.vtp, ...  centers; Point Data: bodyId, TCCType, radius
    Sphere_Results.pvd         compact spheres; radius is needed for sizing
    Particle_Results.pvsm      ready-to-load state showing all geometry types
The STL or sphere series is omitted if that type never occurs. Empty steps
are retained. Mixed STL/sphere cases use both series with one common TCC
analysis, including motifs spanning the two types. Exactly ONE state file,
Particle_Results.pvsm, is written for pure-sphere, pure-STL and mixed cases.
The PVD/VTP files are datasets referenced by this state, not alternate states.
--state-only rebuilds the state from existing PVD files and keeps its camera.
It does not read particleData, rerun classification, or rewrite any PVD/VTP.
Redundant Sphere_Results.pvsm/STL_Results.pvsm copies from the earlier version
are removed only when their contents exactly match the old combined state.

In ParaView, File > Load State > Particle_Results.pvsm shows TCC categories
and correctly sized sphere glyphs. Change Color By to bodyId or TCCType.
When moving output to another computer, use Load State's data-directory
option to locate the accompanying PVD/VTP files. Surface PVD can also be
opened directly. The state fixes TCC colors/labels for all timesteps.

VTP uses little-endian appended binary arrays (readable by ParaView/VTK).
STL triangles retain their coordinates and winding. Sphere surface meshes
are not duplicated on disk. The numerical arrays and consumed geometry are
SHA-256 checked. Reruns REBUILD the complete series in a staging directory;
completed files replace previous outputs only after all inputs succeed.
Old merged .stl files are not used or deleted. No tccCounts.dat is modified.

Algorithm sources retained from analyze_tcc.py:
    https://doi.org/10.1038/s41567-023-02063-x
    https://arxiv.org/abs/1307.5517
"""
# Original STL merge script created by OStudenik.

CASE_DIRECTORY = "."
OUTPUT_DIRECTORY = "STLMerged"
BOND_CUTOFF_M = None
BOND_CUTOFF_FACTOR = 1.295
REFERENCE_DIAMETER_M = None

import argparse
from array import array
from collections import defaultdict
from decimal import Decimal, InvalidOperation
import hashlib
from itertools import chain, combinations, islice, product, repeat
import json
import math
import os
from pathlib import Path
import re
import statistics
import struct
import sys
import tempfile
import xml.etree.ElementTree as ET


class SnapshotError(ValueError):
    """An input archive or consumed array is incomplete or malformed."""


class BinarySnapshot:
    """Read checked numerical arrays without extracting original files.

    Triangle coordinates are streamed directly into VTP, so their memory
    requirement does not grow with the total number of surface triangles.
    """
    HEADER = struct.Struct("<16sI")
    FOOTER = struct.Struct("<16sQQ32s")
    DTYPES = {"<i8": ("q", 8), "<u8": ("Q", 8),
              "<f8": ("d", 8), "u1": ("B", 1)}

    def __init__(self, path):
        self.path = Path(path)
        self.ranges = {}
        self.geometry_bounds = None

    def __enter__(self):
        self.stream = self.path.open("rb")
        try:
            self._read_index()
        except BaseException:
            self.stream.close()
            raise
        return self

    def __exit__(self, *args):
        self.stream.close()

    def fail(self, message):
        raise SnapshotError("{}: {}".format(self.path, message))

    def _read_index(self):
        file_size = self.stream.seek(0, os.SEEK_END)
        if file_size < self.HEADER.size + self.FOOTER.size:
            self.fail("truncated archive")
        self.stream.seek(0)
        if self.HEADER.unpack(self.stream.read(self.HEADER.size)) != (b"PDBIN-HEADER-v1!", 1):
            self.fail("unsupported header/version")
        self.stream.seek(-self.FOOTER.size, os.SEEK_END)
        magic, offset, length, digest = self.FOOTER.unpack(self.stream.read(self.FOOTER.size))
        if (magic != b"PDBIN-FOOTER-v1!" or offset < self.HEADER.size
                or offset + length != file_size - self.FOOTER.size):
            self.fail("invalid index footer")
        self.stream.seek(offset)
        payload = self.stream.read(length)
        if len(payload) != length or hashlib.sha256(payload).digest() != digest:
            self.fail("index SHA-256 mismatch")
        index = json.loads(payload)
        if (not isinstance(index, dict)
                or index.get("format") != "openHFDIB-DEM particleData binary"
                or type(index.get("schema_version")) is not int
                or index.get("schema_version") != 1
                or index.get("time_name") != self.path.stem
                or index.get("byte_order") != "little"
                or index.get("compression") != "none"
                or index.get("length_unit") != "m"
                or index.get("time_unit") != "s"):
            self.fail("unsupported schema, units or time name")
        if not isinstance(index.get("datasets"), dict):
            self.fail("missing datasets index")
        self.index_offset = offset
        self.datasets = index["datasets"]

    def describe(self, name, dtype, shape):
        item = self.datasets.get(name)
        if not isinstance(item, dict):
            self.fail("missing " + name)
        actual = item.get("shape")
        if (item.get("dtype") != dtype or not isinstance(actual, list)
                or any(type(n) is not int or n < 0 for n in actual)
                or actual != list(shape)):
            self.fail("invalid dtype/shape for " + name)
        typecode, width = self.DTYPES[dtype]
        if array(typecode).itemsize != width:
            self.fail("unsupported native numeric array size")
        start, size = item.get("offset"), item.get("size")
        if (type(start) is not int or type(size) is not int
                or start < self.HEADER.size or size != width * math.prod(shape)
                or start + size > self.index_offset
                or not isinstance(item.get("sha256"), str)
                or re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) is None):
            self.fail("invalid offset/size/checksum for " + name)
        if size:
            for other, (begin, end) in self.ranges.items():
                if other != name and start < end and begin < start + size:
                    self.fail("overlapping datasets: {} and {}".format(name, other))
            self.ranges[name] = (start, start + size)
        return item

    def read_array(self, name, dtype, shape):
        item = self.describe(name, dtype, shape)
        self.stream.seek(item["offset"])
        payload = self.stream.read(item["size"])
        if (len(payload) != item["size"]
                or hashlib.sha256(payload).hexdigest() != item["sha256"]):
            self.fail("SHA-256 mismatch for " + name)
        values = array(self.DTYPES[dtype][0])
        values.frombytes(payload)
        if sys.byteorder != "little":
            values.byteswap()
        return values

    def read_particles(self):
        item = self.datasets.get("particles/body_id", {})
        shape = item.get("shape") if isinstance(item, dict) else None
        if (not isinstance(shape, list) or len(shape) != 1
                or type(shape[0]) is not int or shape[0] < 0):
            self.fail("body_id must have shape [N]")
        n = shape[0]
        ids = self.read_array("particles/body_id", "<i8", [n])
        xyz = self.read_array("particles/position", "<f8", [n, 3])
        volumes = self.read_array("particles/volume", "<f8", [n])
        kinds = self.read_array("particles/geometry_type", "u1", [n])
        radii = self.read_array("particles/radius", "<f8", [n])
        offsets = self.read_array("geometry/triangle_offsets", "<u8", [n + 1])
        if len(set(ids)) != n:
            self.fail("duplicate body IDs")
        if not all(math.isfinite(x) for x in xyz):
            self.fail("nonfinite particle positions")
        if not all(math.isfinite(v) and v > 0 for v in volumes):
            self.fail("nonpositive/nonfinite particle volumes")
        if offsets[0] != 0 or any(a > b for a, b in zip(offsets, offsets[1:])):
            self.fail("invalid triangle offsets")
        for i, kind in enumerate(kinds):
            triangle_count = offsets[i + 1] - offsets[i]
            if kind == 0:
                if triangle_count or not math.isfinite(radii[i]) or radii[i] <= 0:
                    self.fail("invalid sphere geometry for bodyId {}".format(ids[i]))
            elif kind == 1:
                if triangle_count == 0:
                    self.fail("empty STL geometry for bodyId {}".format(ids[i]))
            else:
                self.fail("unsupported geometry type {}".format(kind))
        self.describe("geometry/triangles", "<f8", [offsets[-1], 3, 3])
        points = list(zip(xyz[0::3], xyz[1::3], xyz[2::3]))
        return dict(ids=ids, points=points, volumes=volumes, kinds=kinds,
                    radii=radii, offsets=offsets)

    def copy_triangles(self, output, count):
        item = self.describe("geometry/triangles", "<f8", [count, 3, 3])
        self.stream.seek(item["offset"])
        remaining, digest = item["size"], hashlib.sha256()
        while remaining:
            # Each chunk contains whole triangles and hence whole XYZ triples.
            payload = self.stream.read(min(remaining, 8192 * 9 * 8))
            if not payload or len(payload) % 72:
                self.fail("truncated triangle data")
            digest.update(payload)
            values = array("d")
            values.frombytes(payload)
            if sys.byteorder != "little":
                values.byteswap()
            if not all(math.isfinite(x) for x in values):
                self.fail("nonfinite triangle coordinates")
            axes = [values[a::3] for a in range(3)]
            bounds = [min(axis) for axis in axes] + [max(axis) for axis in axes]
            self.geometry_bounds = merge_bounds(self.geometry_bounds, bounds)
            output.write(payload)
            remaining -= len(payload)
        if digest.hexdigest() != item["sha256"]:
            self.fail("SHA-256 mismatch for geometry/triangles")


def time_snapshots(directory):
    """Keep original stems, sort numerically, reject ambiguous duplicate times."""
    if not directory.is_dir():
        raise ValueError("particleData directory does not exist: {}".format(directory))
    selected = {}
    pattern = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
    for path in directory.glob("*.bin"):
        if not path.is_file() or not pattern.fullmatch(path.stem):
            continue
        try:
            value = Decimal(path.stem)
        except InvalidOperation:
            continue
        if not value.is_finite() or value < 0:
            continue
        if not math.isfinite(float(value)):
            raise ValueError("time cannot be represented in ParaView: " + path.stem)
        if value in selected:
            raise ValueError("duplicate numeric time: {} and {}".format(
                selected[value].name, path.name))
        selected[value] = path
    if not selected:
        raise ValueError("no nonnegative numeric-time .bin files in {}".format(directory))
    return [(selected[t].stem, selected[t]) for t in sorted(selected)]


# Predicates copied unchanged from the supplied analyze_tcc.py.
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


def particle_labels(count, clusters):
    """Apply priority only after every overlapping/nested motif is identified."""
    labels = array("B", [0]) * count
    for rank, motif in enumerate(MOTIFS, start=1):
        for members in clusters[motif]:
            for index in members:
                labels[index] = rank
    return labels


def merge_bounds(current, addition):
    if addition is None:
        return current
    if current is None:
        return list(addition)
    return ([min(current[a], addition[a]) for a in range(3)]
            + [max(current[a + 3], addition[a + 3]) for a in range(3)])


def sphere_bounds(particles):
    bounds = None
    for i, kind in enumerate(particles["kinds"]):
        if kind == 0:
            center, radius = particles["points"][i], particles["radii"][i]
            bounds = merge_bounds(bounds, [x - radius for x in center]
                                  + [x + radius for x in center])
    return bounds


def write_numbers(output, values, typecode):
    """Write little-endian numeric data with bounded temporary memory."""
    iterator = iter(values)
    while True:
        chunk = array(typecode, islice(iterator, 65536))
        if not chunk:
            break
        if sys.byteorder != "little":
            chunk.byteswap()
        output.write(chunk.tobytes())


def numeric_array(section, name, vtk_type, values, count, components=1):
    typecode, width = {"Int64": ("q", 8), "Float64": ("d", 8),
                       "UInt8": ("B", 1)}[vtk_type]
    if array(typecode).itemsize != width:
        raise ValueError("unsupported native numeric array size")
    return dict(section=section, name=name, vtk_type=vtk_type,
                components=components, size=count * components * width,
                write=lambda output: write_numbers(output, values, typecode))


def write_vtp(path, number_points, number_verts, number_polys, arrays, time, cutoff):
    root = ET.Element("VTKFile", type="PolyData", version="1.0",
                      byte_order="LittleEndian", header_type="UInt64")
    poly = ET.SubElement(root, "PolyData")
    metadata = ET.SubElement(poly, "FieldData")
    for name, value in (("TimeValue", float(time)), ("TCCBondCutoff_m", cutoff or 0.0)):
        field = ET.SubElement(metadata, "DataArray", type="Float64", Name=name,
                              NumberOfTuples="1", format="ascii")
        field.text = format(value, ".17g")
    piece = ET.SubElement(poly, "Piece", NumberOfPoints=str(number_points),
                          NumberOfVerts=str(number_verts), NumberOfLines="0",
                          NumberOfStrips="0", NumberOfPolys=str(number_polys))
    sections = {}
    for name in ("PointData", "CellData", "Points", "Verts", "Polys"):
        sections[name] = ET.SubElement(piece, name)
    offset = 0
    for spec in arrays:
        attributes = dict(type=spec["vtk_type"], NumberOfComponents=str(spec["components"]),
                          format="appended", offset=str(offset))
        if spec["name"] is not None:
            attributes["Name"] = spec["name"]
        ET.SubElement(sections[spec["section"]], "DataArray", attributes)
        if spec["name"] == "TCCType":
            sections[spec["section"]].set("Scalars", "TCCType")
        offset += 8 + spec["size"]
    xml = ET.tostring(root, encoding="utf-8")
    with path.open("wb") as output:
        output.write(b'<?xml version="1.0"?>\n')
        output.write(b'<!-- TCCType: 0=None; 1=4A; 2=5A; 3=6Z; 4=6A; 5=7A -->\n')
        output.write(xml[:-len(b"</VTKFile>")])
        output.write(b'\n<AppendedData encoding="raw">\n_')
        for spec in arrays:
            output.write(struct.pack("<Q", spec["size"]))
            begin = output.tell()
            spec["write"](output)
            if output.tell() - begin != spec["size"]:
                raise ValueError("incorrect VTP array length: {}".format(spec["name"]))
        output.write(b'\n</AppendedData>\n</VTKFile>\n')


def write_spheres(path, particles, labels, time, cutoff):
    selected = [i for i, kind in enumerate(particles["kinds"]) if kind == 0]
    count = len(selected)
    arrays = [
        numeric_array("PointData", "bodyId", "Int64",
                      (particles["ids"][i] for i in selected), count),
        numeric_array("PointData", "TCCType", "UInt8",
                      (labels[i] for i in selected), count),
        numeric_array("PointData", "radius", "Float64",
                      (particles["radii"][i] for i in selected), count),
        numeric_array("Points", None, "Float64",
                      chain.from_iterable(particles["points"][i] for i in selected), count, 3),
        numeric_array("Verts", "connectivity", "Int64", range(count), count),
        numeric_array("Verts", "offsets", "Int64", range(1, count + 1), count),
    ]
    write_vtp(path, count, count, 0, arrays, time, cutoff)


def write_surfaces(path, particles, labels, snapshot, time, cutoff):
    offsets = particles["offsets"]
    count = offsets[-1]

    def per_triangle(values):
        for i, value in enumerate(values):
            yield from repeat(value, offsets[i + 1] - offsets[i])

    arrays = [
        numeric_array("CellData", "bodyId", "Int64", per_triangle(particles["ids"]), count),
        numeric_array("CellData", "TCCType", "UInt8", per_triangle(labels), count),
        dict(section="Points", name=None, vtk_type="Float64", components=3,
             size=count * 72, write=lambda output: snapshot.copy_triangles(output, count)),
        numeric_array("Polys", "connectivity", "Int64", range(3 * count), 3 * count),
        numeric_array("Polys", "offsets", "Int64", range(3, 3 * count + 1, 3), count),
    ]
    write_vtp(path, 3 * count, 0, count, arrays, time, cutoff)


def write_empty(path, family, time, cutoff):
    particles = dict(ids=[], points=[], volumes=[], kinds=[], radii=[], offsets=[0])
    if family == "Sphere":
        write_spheres(path, particles, [], time, cutoff)
    else:
        class EmptyGeometry:
            def copy_triangles(self, output, count):
                pass
        write_surfaces(path, particles, [], EmptyGeometry(), time, cutoff)


def write_pvd(path, family, snapshots):
    root = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    collection = ET.SubElement(root, "Collection")
    for index, (time, _) in enumerate(snapshots, start=1):
        ET.SubElement(collection, "DataSet", timestep=time, group="", part="0",
                      file="{}_Results{:04d}.vtp".format(family, index))
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def write_state(path, families, final_directory, bounds, times, camera=None):
    """Keep the existing compact sphere-glyph setup; add surfaces and TCC colors."""
    root = ET.Element("GenericParaViewApplication")
    manager = ET.SubElement(root, "ServerManagerState", version="5.13.1")

    def proxy(group, kind, identifier, servers="21"):
        return ET.SubElement(manager, "Proxy", group=group, type=kind,
                             id=str(identifier), servers=servers)

    def prop(parent, name, values, refs=False):
        node = ET.SubElement(parent, "Property", name=name,
                             id=parent.get("id") + "." + name,
                             number_of_elements=str(len(values)))
        for i, value in enumerate(values):
            if refs:
                attributes = {"value": str(value)}
                if name in ("Input", "GlyphType"):
                    attributes["output_port"] = "0"
                ET.SubElement(node, "Proxy", **attributes)
            else:
                ET.SubElement(node, "Element", index=str(i), value=str(value))
        return node

    source_ids = [200 + i for i in range(len(families))]
    display_ids = [500 + i for i in range(len(families))]
    keeper = proxy("misc", "TimeKeeper", 100, "16")
    prop(keeper, "TimeSources", source_ids, True)
    prop(keeper, "Views", [400], True)
    scene = proxy("animation", "AnimationScene", 110, "16")
    for name, values, refs in (("Cues", [120], True), ("PlayMode", [2], False),
                               ("TimeKeeper", [100], True), ("ViewModules", [400], True),
                               ("StartTime", [float(times[0])], False),
                               ("EndTime", [float(times[-1])], False),
                               ("AnimationTime", [float(times[0])], False)):
        prop(scene, name, values, refs)
    cue = proxy("animation", "TimeAnimationCue", 120, "16")
    prop(cue, "AnimatedPropertyName", ["Time"])
    prop(cue, "AnimatedProxy", [100], True)
    prop(cue, "Enabled", [1])
    prop(cue, "UseAnimationTime", [1])

    for family, identifier in zip(families, source_ids):
        reader = proxy("sources", "PVDReader", identifier, "1")
        prop(reader, "FileName", [str(final_directory / (family + "_Results.pvd"))])
        if family == "Sphere":
            prop(reader, "PointArrayStatus", ["bodyId", 1, "TCCType", 1, "radius", 1])
        else:
            prop(reader, "CellArrayStatus", ["bodyId", 1, "TCCType", 1])
    if "Sphere" in families:
        sphere = proxy("sources", "SphereSource", 300)
        for name, values in (("Center", [0, 0, 0]), ("Radius", [1]),
                             ("PhiResolution", [32]), ("ThetaResolution", [32])):
            prop(sphere, name, values)

    focal = [0, 0, 0] if bounds is None else [0.5 * (bounds[a] + bounds[a + 3]) for a in range(3)]
    span = [1, 1, 1] if bounds is None else [bounds[a + 3] - bounds[a] for a in range(3)]
    scale = 1.1 * max(0.5 * span[1], 0.5 * span[0] / (4.0 / 3.0), 0.05 * max(span), 1e-12)
    distance = max(4 * max(span), 4 * scale, 1e-9)
    view = proxy("views", "RenderView", 400)
    for name, values, refs in (("Representations", display_ids + [800], True),
                               ("ViewSize", [800, 600], False),
                               ("CenterOfRotation", focal, False),
                               ("CameraFocalPoint", focal, False),
                               ("CameraPosition", [focal[0], focal[1], focal[2] + distance], False),
                               ("CameraViewUp", [0, 1, 0], False),
                               ("CameraParallelProjection", [1], False),
                               ("CameraParallelScale", [scale], False)):
        prop(view, name, values, refs)
    if camera:
        for name, values in camera.items():
            previous = view.find("Property[@name='{}']".format(name))
            if previous is not None:
                view.remove(previous)
            prop(view, name, values)
    for family, source_id, display_id in zip(families, source_ids, display_ids):
        display = proxy("representations", "GeometryRepresentation", display_id)
        prop(display, "Input", [source_id], True)
        prop(display, "Representation", ["3D Glyphs" if family == "Sphere" else "Surface"])
        prop(display, "Visibility", [1])
        prop(display, "ColorArrayName", ["", "", "", 0 if family == "Sphere" else 1, "TCCType"])
        prop(display, "LookupTable", [700], True)
        prop(display, "ScalarOpacityFunction", [710], True)
        prop(display, "TransferFunction2D", [720], True)
        prop(display, "MapScalars", [1])
        prop(display, "InterpolateScalarsBeforeMapping", [0])
        if family == "Sphere":
            glyph = prop(display, "GlyphType", [300], True)
            ET.SubElement(glyph, "Domain", name="input_type", id=str(display_id) + ".GlyphType.input_type")
            domain = ET.SubElement(glyph, "Domain", name="proxy_list", id=str(display_id) + ".GlyphType.proxy_list")
            ET.SubElement(domain, "Proxy", value="300")
            for name, values in (("Masking", [0]), ("Orient", [0]), ("ScaleFactor", [1]),
                                 ("ScaleMode", [1]), ("Scaling", [1]),
                                 ("SelectScaleArray", ["radius"]),
                                 ("SetScaleArray", ["", "", "", 0, "radius"])):
                prop(display, name, values)
    table = proxy("lookup_tables", "PVLookupTable", 700)
    # XML uses internal property names, not Python/UI aliases:
    # InterpretValuesAsCategories in paraview.simple is IndexedLookup here.
    prop(table, "IndexedLookup", [1])
    prop(table, "AnnotationsInitialized", [1])
    prop(table, "AutomaticRescaleRangeMode", [-1])
    prop(table, "ScalarRangeInitialized", [1])
    prop(table, "RGBPoints", [0, 0.65, 0.65, 0.65, 5, 0.85, 0.20, 0.20])
    prop(table, "ScalarOpacityFunction", [710], True)
    prop(table, "TransferFunction2D", [720], True)
    prop(table, "Annotations", list(chain.from_iterable(
        (str(i), name) for i, name in enumerate(("None",) + MOTIFS))))
    prop(table, "IndexedColors", [0.65, 0.65, 0.65, 0.18, 0.50, 0.80,
                                   0.20, 0.70, 0.45, 0.90, 0.55, 0.10,
                                   0.60, 0.40, 0.75, 0.85, 0.20, 0.20])
    prop(table, "IndexedOpacities", [1] * 6)
    # ParaView's color editor expects the companion transfer functions even
    # though these opaque surfaces do not use opacity or 2D color mapping.
    opacity = proxy("piecewise_functions", "PiecewiseFunction", 710)
    prop(opacity, "Points", [0, 0, 0.5, 0, 5, 1, 0.5, 0])
    prop(opacity, "ScalarRangeInitialized", [1])
    transfer2d = proxy("transfer_2d_functions", "TransferFunction2D", 720)
    prop(transfer2d, "Range", [0, 5, 0, 1])
    legend = proxy("representations", "ScalarBarWidgetRepresentation", 800)
    prop(legend, "LookupTable", [700], True)
    prop(legend, "Title", ["TCCType"])
    prop(legend, "ComponentTitle", [""])
    prop(legend, "Visibility", [1])
    layout = proxy("misc", "ViewLayout", 600, "20")
    ET.SubElement(ET.SubElement(layout, "Layout", number_of_elements="1"),
                  "Item", direction="0", fraction="0.5", view="400")
    collections = {
        "animation": [(110, "AnimationScene1"), (120, "TimeAnimationCue1")],
        "layouts": [(600, "Layout #1")],
        "representations": [(i, "GeometryRepresentation" + str(i)) for i in display_ids],
        "scalar_bars": [(800, "TCCTypeScalarBar")],
        "sources": [(i, family + "_Results.pvd") for family, i in zip(families, source_ids)],
        "lookup_tables": [(700, "TCCType.PVLookupTable")],
        "piecewise_functions": [(710, "TCCType.PiecewiseFunction")],
        "transfer_2d_functions": [(720, "TCCType.TransferFunction2D")],
        "timekeeper": [(100, "TimeKeeper1")], "views": [(400, "RenderView1")],
    }
    if "Sphere" in families:
        sphere_display = display_ids[families.index("Sphere")]
        collections["pq_helper_proxies." + str(sphere_display)] = [(300, "GlyphType")]
    for name, items in collections.items():
        collection = ET.SubElement(manager, "ProxyCollection", name=name)
        for identifier, title in items:
            ET.SubElement(collection, "Item", id=str(identifier), name=title)
    ET.SubElement(manager, "CustomProxyDefinitions")
    ET.SubElement(manager, "Links")
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def duplicate_states(output):
    """Remember only byte-identical aliases emitted by the previous version."""
    primary = output / "Particle_Results.pvsm"
    if not primary.is_file():
        return []
    old = primary.read_bytes()
    return [(path, old) for path in (output / "Sphere_Results.pvsm", output / "STL_Results.pvsm")
            if path.is_file() and path.read_bytes() == old]


def remove_duplicate_states(duplicates):
    for path, expected in duplicates:
        if path.is_file() and path.read_bytes() == expected:
            path.unlink()


def repair_state(output):
    """Repair visualization metadata without touching existing particle data."""
    primary = output / "Particle_Results.pvsm"
    if not primary.is_file():
        raise ValueError("--state-only requires an existing " + str(primary))
    old = ET.parse(primary)
    families = []
    for reader in old.findall(".//Proxy[@type='PVDReader']"):
        element = reader.find("Property[@name='FileName']/Element")
        filename = Path(element.get("value", "")).name if element is not None else ""
        family = next((f for f in ("STL", "Sphere") if filename == f + "_Results.pvd"), None)
        if family is None or family in families:
            raise ValueError("unsupported or duplicate reader in existing state: " + filename)
        families.append(family)
    if not families:
        raise ValueError("existing state has no particle PVD readers")
    times = None
    for family in families:
        collection = ET.parse(output / (family + "_Results.pvd"))
        datasets = collection.findall("./Collection/DataSet")
        current = [item.get("timestep", "") for item in datasets]
        if not current or any(not math.isfinite(float(t)) or float(t) < 0 for t in current):
            raise ValueError("invalid times in " + family + "_Results.pvd")
        if any(Decimal(a) >= Decimal(b) for a, b in zip(current, current[1:])):
            raise ValueError("PVD times must be unique and increasing")
        if times is not None and list(map(Decimal, current)) != list(map(Decimal, times)):
            raise ValueError("STL and sphere PVD times do not match")
        times = current
        for item in datasets:
            if not item.get("file") or not (output / item.get("file")).is_file():
                raise ValueError("PVD references a missing dataset: " + str(item.get("file")))
    camera = {}
    view = old.find(".//Proxy[@type='RenderView']")
    if view is not None:
        for name, count in (("CenterOfRotation", 3), ("CameraFocalPoint", 3),
                            ("CameraPosition", 3), ("CameraViewUp", 3),
                            ("CameraParallelProjection", 1), ("CameraParallelScale", 1),
                            ("CameraViewAngle", 1)):
            elements = view.findall("Property[@name='{}']/Element".format(name))
            values = [e.get("value", "") for e in elements]
            if values:
                if len(values) != count or any(not math.isfinite(float(v)) for v in values):
                    raise ValueError("invalid saved camera property: " + name)
                camera[name] = values
    duplicates = duplicate_states(output)
    with tempfile.TemporaryDirectory(prefix=".particle-state-", dir=output) as staging_name:
        state = Path(staging_name) / primary.name
        write_state(state, families, output, None, times, camera=camera)
        os.replace(state, primary)
    remove_duplicate_states(duplicates)
    print("Repaired state only; existing PVD/VTP files and TCC results retained.", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", type=Path, default=Path(CASE_DIRECTORY))
    parser.add_argument("--output-dir", type=Path, default=Path(OUTPUT_DIRECTORY))
    parser.add_argument("--state-only", action="store_true",
                        help="repair existing Particle_Results.pvsm without recomputing TCC or data")
    parser.add_argument("--bond-cutoff", type=float, default=BOND_CUTOFF_M,
                        help="fixed center-distance cutoff in m (overrides factor/diameter)")
    parser.add_argument("--bond-factor", type=float, default=BOND_CUTOFF_FACTOR)
    parser.add_argument("--diameter", type=float, default=REFERENCE_DIAMETER_M,
                        help="reference volume-equivalent diameter in m")
    args = parser.parse_args(argv)
    try:
        for name, value in (("bond cutoff", args.bond_cutoff), ("bond factor", args.bond_factor),
                            ("reference diameter", args.diameter)):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(name + " must be finite and positive")
        case = args.case.resolve()
        data = case / "particleData"
        output = (args.output_dir if args.output_dir.is_absolute() else case / args.output_dir).resolve()
        if output == data.resolve() or data.resolve() in output.parents:
            raise ValueError("output directory must be outside particleData")
        if args.state_only:
            repair_state(output)
            print("Load in ParaView: {}".format(output / "Particle_Results.pvsm"), flush=True)
            return 0
        snapshots = time_snapshots(data)
        cutoff = args.bond_cutoff
        if cutoff is None and args.diameter is not None:
            cutoff = args.bond_factor * args.diameter
        if cutoff is not None and (not math.isfinite(cutoff) or cutoff <= 0):
            raise ValueError("computed cutoff must be finite and positive")
        output.mkdir(parents=True, exist_ok=True)
        duplicates = duplicate_states(output)
        print("Analyzing ALL stored particles at {} times; no spatial filter.".format(len(snapshots)), flush=True)
        if cutoff is not None:
            print("Fixed bond cutoff: {:.17g} m".format(cutoff), flush=True)
        bounds, present = None, set()
        with tempfile.TemporaryDirectory(prefix=".particle-vtp-", dir=output) as staging_name:
            staging = Path(staging_name)
            for output_index, (time, path) in enumerate(snapshots, start=1):
                with BinarySnapshot(path) as snapshot:
                    particles = snapshot.read_particles()
                    count = len(particles["ids"])
                    if cutoff is None and count:
                        diameter = statistics.median((6.0 * v / math.pi)**(1.0/3.0) for v in particles["volumes"])
                        cutoff = args.bond_factor * diameter
                        if not math.isfinite(cutoff) or cutoff <= 0:
                            raise ValueError("computed cutoff must be finite and positive")
                        print("Provisional fixed cutoff: {:.17g} m = {} x {:.17g} m; reference t={} s".format(
                            cutoff, args.bond_factor, diameter, time), flush=True)
                    print("t={} s: {} particles; classifying...".format(time, count), flush=True)
                    clusters = classify_clusters(particles["points"], cutoff) if count else {m: set() for m in MOTIFS}
                    labels = particle_labels(count, clusters)
                    print("  clusters: " + "  ".join("{}={}".format(m, len(clusters[m])) for m in MOTIFS), flush=True)
                    if any(kind == 1 for kind in particles["kinds"]):
                        present.add("STL")
                        write_surfaces(staging / "STL_Results{:04d}.vtp".format(output_index),
                                       particles, labels, snapshot, time, cutoff)
                        bounds = merge_bounds(bounds, snapshot.geometry_bounds)
                    else:
                        # Also verify the empty geometry dataset's checksum.
                        with open(os.devnull, "wb") as sink:
                            snapshot.copy_triangles(sink, 0)
                    if any(kind == 0 for kind in particles["kinds"]):
                        present.add("Sphere")
                        write_spheres(staging / "Sphere_Results{:04d}.vtp".format(output_index),
                                      particles, labels, time, cutoff)
                        bounds = merge_bounds(bounds, sphere_bounds(particles))
            families = [family for family in ("STL", "Sphere") if family in present] or ["Sphere"]
            for family in families:
                for index, (time, _) in enumerate(snapshots, start=1):
                    path = staging / "{}_Results{:04d}.vtp".format(family, index)
                    if not path.exists():
                        write_empty(path, family, time, cutoff)
                write_pvd(staging / (family + "_Results.pvd"), family, snapshots)
            state = staging / "Particle_Results.pvsm"
            write_state(state, families, output, bounds, [time for time, _ in snapshots])
            # Publish only after all input checks, classification and exports succeed.
            # Data files first; collection/state files last. Never append to old STL.
            for path in sorted(staging.iterdir(), key=lambda p: (p.suffix != ".vtp", p.name)):
                os.replace(path, output / path.name)
        remove_duplicate_states(duplicates)
        print("Load in ParaView: {}".format(output / "Particle_Results.pvsm"), flush=True)
        print("Color By: bodyId or TCCType (0=None, 1=4A, 2=5A, 3=6Z, 4=6A, 5=7A).", flush=True)
        return 0
    except (OSError, ValueError, TypeError, KeyError, OverflowError, struct.error, ET.ParseError) as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
