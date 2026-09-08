#!/usr/bin/env python3
"""Reversible bodiesInfo <-> particleData conversion using ONLY Python's stdlib.

NO installation, NumPy, h5py, module load, compression, or extraction utility.
Run after the simulation has finished, inside the case directory:
    python3 compress_bodiesInfo.py
    python3 compress_bodiesInfo.py decompress
    python3 compress_bodiesInfo.py --keep-original
    python3 compress_bodiesInfo.py --case /path/to/case

The default command writes particleData/<time>.bin, verifies its checksums,
then removes bodiesInfo/<time> (and bodiesInfo itself if empty). Decompress
restores exact original bytes, names, directories, modes and modification
times; .bin files are kept. Different existing contents are never overwritten.
A small OS lock file in the case prevents simultaneous conversions. Verified
archives survive interrupted source deletion, which the next run completes.

Each .bin is an UNCOMPRESSED, indexed binary data file (format version 1):
    fixed header | raw original bytes and numerical arrays | JSON index | footer
The JSON index stores array/file byte offsets, types, shapes and SHA-256 hashes.
Its footer stores the index offset, length and SHA-256. All numerical payloads
are little-endian contiguous arrays. No pickle or executable serialization.
Array data can be read directly, without loading/extracting original files.

Direct access from other Python scripts (this file is also the reader module):
    from compress_bodiesInfo import ParticleDataFile
    with ParticleDataFile("particleData/0.15.bin") as data:
        ids = data.read_array("particles/body_id")
        xyz = data.read_array("particles/position")  # flat x0,y0,z0,x1,y1,z1,...
        volume = data.read_array("particles/volume")
        centers = zip(xyz[0::3], xyz[1::3], xyz[2::3])
        mesh = data.read_triangles(0)  # flat 9 numbers/triangle; empty for sphere

Schema (SI units):
    particles/body_id                int64 [N], sorted by body ID
    particles/body_name              strings [N], read_strings(name)
    particles/geometry_type          uint8 [N]: 0=sphere, 1=STL
    particles/position               float64 [N,3]
    particles/radius                 float64 [N], NaN for STL
    particles/volume                 float64 [N]
    particles/velocity               float64 [N,3]
    particles/omega                  float64 [N], angular speed
    particles/axis                   float64 [N,3], angular-velocity axis
    particles/static                 uint8 [N]
    particles/time_steps_in_contact  int64 [N]
    geometry/triangles               float64 [M,3,3]
    geometry/triangle_offsets        uint64 [N+1]
    original/data                    byte [B], exact original file contents

Particle i has triangles[offsets[i]:offsets[i+1]]. The binary index lists
original filenames and byte ranges, modes, timestamps and empty directories.
Use read_original(name) for a particular original file, or verify() to check
all payloads. Normal array reads do not scan unrelated geometry or originals.

STL centers use the unique-vertex arithmetic mean, matching the original
upper_reservoir_solid_fraction.py, not a volume centroid. Mesh volume is the
absolute sum of signed tetrahedron volumes, requiring a closed, consistently
oriented STL; no geometry repair is performed. ASCII and binary STL work.
There is no equal-volume assumption. Only one particle mesh is parsed at a time.
The raw bytes and analysis arrays are both stored intentionally: reversibility
and efficient repeated analysis are prioritized over disk-space savings.

Requires Python 3.9+ only. This .bin version does not read the previous HDF5
format. The previous missing-dependency error occurred before any conversion.
"""
from __future__ import annotations

import argparse
from array import array
from contextlib import contextmanager
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import struct
import sys
import tempfile

class ParticleParseError(ValueError):
    """An input particle record or STL cannot be interpreted safely."""


def _foam_tokens(data, name):
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ParticleParseError(f"{name}: body info is not UTF-8 text") from exc
    # Preserve quoted strings, including comment-like characters inside them.
    pattern = re.compile(r'"(?:\\.|[^"\\])*"|/\*[\s\S]*?\*/|//[^\r\n]*')
    text = pattern.sub(lambda m: m[0] if m[0].startswith('"') else " ", text)
    if "/*" in re.sub(r'"(?:\\.|[^"\\])*"', "", text):
        raise ParticleParseError(f"{name}: unterminated OpenFOAM comment")
    token_pattern = re.compile(r'"(?:\\.|[^"\\])*"|[{}();]|[^\s{}();"]+')
    tokens, end = [], 0
    for match in token_pattern.finditer(text):
        if text[end:match.start()].strip():
            raise ParticleParseError(f"{name}: malformed quoted string")
        tokens.append(match[0])
        end = match.end()
    if text[end:].strip():
        raise ParticleParseError(f"{name}: malformed quoted string")
    return tokens


def _foam_dictionary(tokens, name, start=0, nested=False):
    """Parse dictionary entries without interpreting unrelated field values."""
    entries, index = {}, start
    while index < len(tokens):
        key = tokens[index]
        index += 1
        if key == "}" and nested:
            return entries, index
        if key in "{}();" or index >= len(tokens):
            raise ParticleParseError(f"{name}: malformed dictionary entry {key!r}")
        if key in entries:
            raise ParticleParseError(f"{name}: duplicate dictionary key {key!r}")
        if tokens[index] == "{":
            value, index = _foam_dictionary(tokens, name, index + 1, True)
            if index < len(tokens) and tokens[index] == ";":
                index += 1
        else:
            value, depth = [], 0
            while index < len(tokens):
                token = tokens[index]
                index += 1
                if token == ";" and depth == 0:
                    break
                if token in ("{", "}", ";"):
                    raise ParticleParseError(f"{name}: malformed value for {key}")
                if token == "(":
                    depth += 1
                elif token == ")":
                    depth -= 1
                    if depth < 0:
                        raise ParticleParseError(f"{name}: unbalanced vector for {key}")
                value.append(token)
            else:
                raise ParticleParseError(f"{name}: missing semicolon for {key}")
            if not value or depth:
                raise ParticleParseError(f"{name}: incomplete value for {key}")
        entries[key] = value
    if nested:
        raise ParticleParseError(f"{name}: unclosed dictionary block")
    return entries, index


def _info_scalar(entries, key, name):
    value = entries.get(key)
    if not isinstance(value, list) or len(value) != 1:
        raise ParticleParseError(f"{name}: missing or invalid {key}")
    return value[0]


def _info_float(entries, key, name):
    text = _info_scalar(entries, key, name)
    try:
        value = float(text)
    except ValueError as exc:
        raise ParticleParseError(f"{name}: invalid number for {key}") from exc
    if not math.isfinite(value):
        raise ParticleParseError(f"{name}: non-finite {key}")
    return value


def _info_int(entries, key, name):
    text = _info_scalar(entries, key, name)
    if re.fullmatch(r"[+-]?\d+", text) is None:
        raise ParticleParseError(f"{name}: invalid integer for {key}")
    value = int(text)
    if not -(1 << 63) <= value < (1 << 63):
        raise ParticleParseError(f"{name}: {key} is outside signed 64-bit range")
    return value


def _info_vector(entries, key, name):
    value = entries.get(key)
    if (not isinstance(value, list) or len(value) != 5
            or value[0] != "(" or value[-1] != ")"):
        raise ParticleParseError(f"{name}: missing or invalid three-component {key}")
    try:
        result = tuple(float(part) for part in value[1:4])
    except ValueError as exc:
        raise ParticleParseError(f"{name}: invalid vector for {key}") from exc
    if not all(math.isfinite(part) for part in result):
        raise ParticleParseError(f"{name}: non-finite vector for {key}")
    return result


def parse_particle_info(data: bytes, name: str) -> dict:
    """Parse the numeric particle fields; original bytes are preserved elsewhere."""
    entries, _ = _foam_dictionary(_foam_tokens(data, name), name)
    body_id = _info_int(entries, "bodyId", name)
    match = re.fullmatch(r"body(-?\d+)\.info", Path(name).name)
    if match is None or int(match[1]) != body_id:
        raise ParticleParseError(f"{name}: bodyId does not match filename")
    body_name = _info_scalar(entries, "bodyName", name)
    if body_name.startswith('"'):
        body_name = re.sub(r'\\(["\\])', r'\1', body_name[1:-1])
    if not body_name:
        raise ParticleParseError(f"{name}: empty bodyName")
    static_text = _info_scalar(entries, "static", name).lower()
    bools = {"0": 0, "false": 0, "no": 0, "off": 0,
             "1": 1, "true": 1, "yes": 1, "on": 1}
    if static_text not in bools:
        raise ParticleParseError(f"{name}: invalid static flag {static_text!r}")
    result = {
        "body_id": body_id,
        "body_name": body_name,
        "velocity": _info_vector(entries, "Vel", name),
        "omega": _info_float(entries, "omega", name),
        "axis": _info_vector(entries, "Axis", name),
        "static": bools[static_text],
        "contact_steps": _info_int(entries, "timeStepsInContWStatic", name),
        "sphere_center": None,
        "radius": None,
    }
    if "sphere" in entries:
        sphere = entries["sphere"]
        if not isinstance(sphere, dict):
            raise ParticleParseError(f"{name}: invalid sphere dictionary")
        result["sphere_center"] = _info_vector(sphere, "position", name)
        radius = _info_float(sphere, "radius", name)
        if radius <= 0:
            raise ParticleParseError(f"{name}: sphere radius must be positive")
        result["radius"] = radius
    return result


"""Standard-library STL helpers; the host defines ParticleParseError."""

from array import array
import math
import re
import struct


def parse_stl(data: bytes, name: str):
    """Return flat float64 triangles, unique-vertex mean, and volume.

    The first result is array('d'), ordered triangle/vertex/coordinate.
    The center matches the existing upper-reservoir script's unique-vertex
    mean, rather than the volume centroid. A closed, consistently wound
    surface is expected; the parser does not repair geometry.
    """
    triangles = array("d")
    unique = set()

    def add_vertex(vertex):
        if not all(math.isfinite(value) for value in vertex):
            raise ParticleParseError(f"{name}: non-finite STL vertex")
        triangles.extend(vertex)
        unique.add(vertex)

    count = None
    if len(data) >= 84:
        declared = struct.unpack_from("<I", data, 80)[0]
        if 84 + 50 * declared == len(data):
            count = declared
    if count is not None:
        if count == 0:
            raise ParticleParseError(f"{name}: binary STL has no triangles")
        # Exact record length identifies binary STL, even when its header
        # begins with the ASCII word 'solid'. Normals are not geometry.
        vertex_record = struct.Struct("<3f")
        for offset in range(84, len(data), 50):
            for vertex_offset in (offset + 12, offset + 24, offset + 36):
                add_vertex(vertex_record.unpack_from(data, vertex_offset))
    else:
        try:
            text = data.decode("ascii")
        except UnicodeDecodeError as exc:
            raise ParticleParseError(f"{name}: invalid or truncated STL") from exc
        if re.search(r"^\s*endsolid\b", text, re.MULTILINE | re.IGNORECASE) is None:
            raise ParticleParseError(f"{name}: ASCII STL is missing endsolid")
        vertices = re.finditer(r"^\s*vertex\s+([^\r\n]*)", text,
                               re.MULTILINE | re.IGNORECASE)
        for match in vertices:
            try:
                vertex = tuple(float(part) for part in match.group(1).split())
                if len(vertex) != 3:
                    raise ValueError("vertex must have three components")
            except ValueError as exc:
                raise ParticleParseError(f"{name}: invalid ASCII STL vertex") from exc
            add_vertex(vertex)
        if not triangles or len(triangles) % 9:
            raise ParticleParseError(f"{name}: incomplete ASCII STL triangles")

    try:
        center = tuple(math.fsum(vertex[axis] for vertex in unique) / len(unique)
                       for axis in range(3))

        def signed_six_volumes():
            # Shift near the mesh before integrating to avoid cancellation
            # from large absolute simulation coordinates.
            cx, cy, cz = center
            for offset in range(0, len(triangles), 9):
                ax = triangles[offset] - cx
                ay = triangles[offset + 1] - cy
                az = triangles[offset + 2] - cz
                bx = triangles[offset + 3] - cx
                by = triangles[offset + 4] - cy
                bz = triangles[offset + 5] - cz
                dx = triangles[offset + 6] - cx
                dy = triangles[offset + 7] - cy
                dz = triangles[offset + 8] - cz
                yield (ax * (by * dz - bz * dy)
                       + ay * (bz * dx - bx * dz)
                       + az * (bx * dy - by * dx))

        volume = abs(math.fsum(signed_six_volumes())) / 6.0
    except (ValueError, OverflowError) as exc:
        raise ParticleParseError(f"{name}: STL has invalid enclosed volume") from exc
    if (not all(math.isfinite(value) for value in center)
            or not math.isfinite(volume) or volume <= 0):
        raise ParticleParseError(f"{name}: STL has invalid or zero enclosed volume")
    return triangles, center, volume


def stl_triangle_count(data: bytes, name: str) -> int:
    """Count triangles for allocation; parse_stl validates all coordinates later."""
    if len(data) >= 84:
        count = struct.unpack_from("<I", data, 80)[0]
        if 84 + 50 * count == len(data):
            if not count:
                raise ParticleParseError(f"{name}: binary STL has no triangles")
            return count
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ParticleParseError(f"{name}: invalid or truncated STL") from exc
    if re.search(r"^\s*endsolid\b", text, re.MULTILINE | re.IGNORECASE) is None:
        raise ParticleParseError(f"{name}: ASCII STL is missing endsolid")
    count = sum(1 for _ in re.finditer(r"^\s*vertex\s+([^\r\n]*)", text,
                                     re.MULTILINE | re.IGNORECASE))
    if not count or count % 3:
        raise ParticleParseError(f"{name}: incomplete ASCII STL triangles")
    return count // 3



# A portable, uncompressed, indexed binary container. No pickle or extensions.
FORMAT = "openHFDIB-DEM particleData binary"
SCHEMA_VERSION = 1
HEADER = struct.Struct("<16sI")
FOOTER = struct.Struct("<16sQQ32s")
HEADER_MAGIC = b"PDBIN-HEADER-v1!"
FOOTER_MAGIC = b"PDBIN-FOOTER-v1!"
DTYPES = {"<f8": ("d", 8), "<i8": ("q", 8), "<u8": ("Q", 8), "u1": ("B", 1)}
BLOCK_BYTES = 8 * 1024 * 1024
TIME_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")
BODY_RE = re.compile(r"body(-?\d+)\.info\Z")


class ConversionError(RuntimeError):
    pass


def little_endian_buffer(values, dtype):
    code, size = DTYPES[dtype]
    result = values if isinstance(values, array) and values.typecode == code else array(code, values)
    if result.itemsize != size:
        raise ConversionError(f"Unsupported native numeric type size: {code}")
    if sys.byteorder != "little" and size > 1:
        result = array(code, result)
        result.byteswap()
    return memoryview(result).cast("B")


class BinaryWriter:
    def __init__(self, stream, time_name):
        self.stream = stream
        self.end = HEADER.size
        stream.write(HEADER.pack(HEADER_MAGIC, SCHEMA_VERSION))
        self.index = {"format": FORMAT, "schema_version": SCHEMA_VERSION,
                      "time_name": time_name, "length_unit": "m", "time_unit": "s",
                      "compression": "none", "byte_order": "little",
                      "created_utc": datetime.now(timezone.utc).isoformat(),
                      "stl_center": "arithmetic mean of unique STL vertices",
                      "stl_volume": "absolute signed tetrahedron volume sum",
                      "geometry_type_codes": {"sphere": 0, "STL": 1},
                      "datasets": {}, "original": {}}

    def write_dataset(self, name, dtype, shape, byte_count, chunks):
        padding = (-self.end) % 8
        self.stream.seek(self.end)
        self.stream.write(b"\0" * padding)
        self.end += padding
        start = self.end
        digest = hashlib.sha256()
        for chunk in chunks:
            # A generator may read earlier original bytes using this stream.
            self.stream.seek(self.end)
            self.stream.write(chunk)
            digest.update(chunk)
            self.end += len(chunk)
        if self.end - start != byte_count:
            raise ConversionError(f"Incorrect byte count for dataset {name}")
        self.index["datasets"][name] = {"offset": start, "size": byte_count,
                                        "dtype": dtype, "shape": list(shape),
                                        "sha256": digest.hexdigest()}

    def write_array(self, name, values, dtype, shape=None):
        buffer = little_endian_buffer(values, dtype)
        size = DTYPES[dtype][1]
        if shape is None:
            shape = (len(buffer) // size,)
        if math.prod(shape) * size != len(buffer):
            raise ConversionError(f"Invalid array shape for {name}")
        self.write_dataset(name, dtype, shape, len(buffer),
                           (buffer[p:p + BLOCK_BYTES] for p in range(0, len(buffer), BLOCK_BYTES)))

    def write_strings(self, name, values):
        payload = json.dumps(values, ensure_ascii=True, separators=(",", ":")).encode("ascii")
        self.write_dataset(name, "json-strings", (len(values),), len(payload), (payload,))

    def read_original(self, entry):
        offset = self.index["datasets"]["original/data"]["offset"] + entry["start"]
        self.stream.seek(offset)
        data = self.stream.read(entry["size"])
        if len(data) != entry["size"]:
            raise ConversionError(f"Truncated original file: {entry['path']}")
        return data

    def finish(self):
        payload = json.dumps(self.index, ensure_ascii=True, allow_nan=False,
                             separators=(",", ":"), sort_keys=True).encode("ascii")
        self.stream.seek(self.end)
        self.stream.write(payload)
        self.stream.write(FOOTER.pack(FOOTER_MAGIC, self.end, len(payload),
                                      hashlib.sha256(payload).digest()))
        self.stream.flush()


class ParticleDataFile:
    """Read numerical fields and exact original bytes without extracting files.

    read_array(name, start=0, stop=None) returns a flat stdlib array.array;
    start/stop are flat element indices. datasets[name]['shape'] gives its shape.
    read_strings(name) returns a list of strings. verify() reads and checks all
    payloads; normal field reads need not scan unrelated geometry/original data.
    """
    def __init__(self, path):
        self.path = Path(path)
        self.stream = self.path.open("rb")
        try:
            self._read_index()
        except BaseException:
            self.stream.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.stream.close()

    def _read_index(self):
        size = self.stream.seek(0, os.SEEK_END)
        if size < HEADER.size + FOOTER.size:
            raise ConversionError(f"Truncated particle-data file: {self.path}")
        self.stream.seek(0)
        if HEADER.unpack(self.stream.read(HEADER.size)) != (HEADER_MAGIC, SCHEMA_VERSION):
            raise ConversionError(f"Not a supported particle-data .bin file: {self.path}")
        self.stream.seek(-FOOTER.size, os.SEEK_END)
        magic, offset, length, digest = FOOTER.unpack(self.stream.read(FOOTER.size))
        if magic != FOOTER_MAGIC or offset < HEADER.size or offset + length != size - FOOTER.size:
            raise ConversionError(f"Invalid or truncated particle-data index: {self.path}")
        self.stream.seek(offset)
        payload = self.stream.read(length)
        if hashlib.sha256(payload).digest() != digest:
            raise ConversionError(f"Index checksum mismatch: {self.path}")
        self.index = json.loads(payload)
        if (self.index.get("format") != FORMAT
                or self.index.get("schema_version") != SCHEMA_VERSION
                or self.index.get("compression") != "none"):
            raise ConversionError(f"Unsupported particle-data schema: {self.path}")
        self.time_name = self.index["time_name"]
        if time_value(self.time_name) is None:
            raise ConversionError("Invalid stored time name")
        self.datasets = self.index["datasets"]
        intervals = []
        for name, item in self.datasets.items():
            shape = item["shape"]
            if (not isinstance(shape, list) or any(type(v) is not int or v < 0 for v in shape)
                    or type(item["offset"]) is not int or type(item["size"]) is not int
                    or item["offset"] < HEADER.size or item["size"] < 0
                    or item["offset"] + item["size"] > offset
                    or re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) is None):
                raise ConversionError(f"Invalid dataset metadata: {name}")
            if item["dtype"] in DTYPES:
                if math.prod(shape) * DTYPES[item["dtype"]][1] != item["size"]:
                    raise ConversionError(f"Invalid dataset size: {name}")
            elif item["dtype"] != "json-strings":
                raise ConversionError(f"Unsupported numeric type: {name}")
            if item["size"]:
                intervals.append((item["offset"], item["offset"] + item["size"]))
        intervals.sort()
        if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
            raise ConversionError("Overlapping binary datasets")
        self._prepare_original_metadata()
        n = self.datasets["particles/body_id"]["shape"]
        if len(n) != 1:
            raise ConversionError("Invalid particle ID array shape")
        for name in ("body_name", "geometry_type", "radius", "volume", "omega", "static",
                     "time_steps_in_contact"):
            if self.datasets[f"particles/{name}"]["shape"] != n:
                raise ConversionError(f"Incorrect particle array shape: {name}")
        for name in ("position", "velocity", "axis"):
            if self.datasets[f"particles/{name}"]["shape"] != [n[0], 3]:
                raise ConversionError(f"Incorrect vector array shape: {name}")
        if (self.datasets["geometry/triangle_offsets"]["shape"] != [n[0] + 1]
                or self.datasets["geometry/triangles"]["shape"][1:] != [3, 3]):
            raise ConversionError("Incorrect geometry array shape")

    def _prepare_original_metadata(self):
        original = self.index["original"]
        files, directories = original["files"], original["directories"]
        self.original_files = {entry["path"]: entry for entry in files}
        names = list(self.original_files)
        dirs = [entry["path"] for entry in directories]
        if (len(names) != len(files) or len(set(dirs)) != len(dirs)
                or set(names) & set(dirs) or "" not in dirs):
            raise ConversionError("Invalid or duplicate original paths")
        for name in names:
            checked_relative(name)
        for name in dirs:
            checked_relative(name, directory=True)
        directory_set = set(dirs)
        for name in names + [d for d in dirs if d]:
            parent = name.rsplit("/", 1)[0] if "/" in name else ""
            if parent not in directory_set:
                raise ConversionError(f"Missing stored parent for {name}")
        end, offsets = 0, [0]
        for entry in files:
            if (type(entry["start"]) is not int or entry["start"] != end
                    or type(entry["size"]) is not int or entry["size"] < 0
                    or re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) is None):
                raise ConversionError("Invalid stored file offset or checksum")
            end += entry["size"]
            offsets.append(end)
        if end != self.datasets["original/data"]["size"]:
            raise ConversionError("Original-file sizes do not match raw data")
        for entry in files + directories:
            if type(entry["mode"]) is not int or type(entry["mtime_ns"]) is not int:
                raise ConversionError("Invalid stored filesystem metadata")
        self.metadata = {"names": names, "directories": dirs, "offsets": offsets,
                         "hashes": [e["sha256"] for e in files],
                         "modes": [e["mode"] for e in files],
                         "mtime_ns": [e["mtime_ns"] for e in files],
                         "directory_modes": [e["mode"] for e in directories],
                         "directory_mtime_ns": [e["mtime_ns"] for e in directories]}

    def chunks(self, offset, size):
        while size:
            self.stream.seek(offset)
            data = self.stream.read(min(size, BLOCK_BYTES))
            if not data:
                raise ConversionError(f"Truncated payload in {self.path}")
            yield data
            offset += len(data)
            size -= len(data)

    def verify(self, expected_time=None):
        if expected_time is not None and self.time_name != expected_time:
            raise ConversionError(f"Time-name mismatch in {self.path}")
        for name, item in self.datasets.items():
            digest = hashlib.sha256()
            for chunk in self.chunks(item["offset"], item["size"]):
                digest.update(chunk)
            if digest.hexdigest() != item["sha256"]:
                raise ConversionError(f"Dataset checksum mismatch: {name} in {self.path}")
        return True

    def read_array(self, name, start=0, stop=None):
        item = self.datasets[name]
        if item["dtype"] not in DTYPES:
            raise ConversionError(f"Not a numeric array: {name}")
        code, size = DTYPES[item["dtype"]]
        count = item["size"] // size
        stop = count if stop is None else stop
        if not 0 <= start <= stop <= count:
            raise IndexError(f"Array range outside {name}")
        result = array(code)
        if result.itemsize != size:
            raise ConversionError(f"Unsupported native numeric type size: {code}")
        for chunk in self.chunks(item["offset"] + start * size, (stop - start) * size):
            result.frombytes(chunk)
        if sys.byteorder != "little" and size > 1:
            result.byteswap()
        return result

    def read_strings(self, name):
        item = self.datasets[name]
        if item["dtype"] != "json-strings":
            raise ConversionError(f"Not a string array: {name}")
        result = json.loads(b"".join(self.chunks(item["offset"], item["size"])))
        if (not isinstance(result, list) or len(result) != item["shape"][0]
                or any(not isinstance(v, str) for v in result)):
            raise ConversionError(f"Invalid string array: {name}")
        return result

    def original_chunks(self, name):
        entry = self.original_files[name]
        start = self.datasets["original/data"]["offset"] + entry["start"]
        return self.chunks(start, entry["size"])

    def read_original(self, name):
        return b"".join(self.original_chunks(name))

    def read_triangles(self, particle_index):
        offsets = self.read_array("geometry/triangle_offsets", particle_index, particle_index + 2)
        return self.read_array("geometry/triangles", 9 * offsets[0], 9 * offsets[1])


def write_original(writer, source, items):
    files, directories = items
    entries = []
    total = sum(item["size"] for item in files)
    def chunks():
        start = 0
        for item in files:
            digest, size = hashlib.sha256(), 0
            with (source / item["path"]).open("rb") as stream:
                for chunk in iter(lambda: stream.read(BLOCK_BYTES), b""):
                    digest.update(chunk)
                    size += len(chunk)
                    if size > item["size"]:
                        raise ConversionError(f"Source changed while reading: {item['path']}")
                    yield chunk
            if size != item["size"]:
                raise ConversionError(f"Source changed while reading: {item['path']}")
            entries.append({"path": item["path"], "start": start, "size": size,
                            "mode": item["mode"], "mtime_ns": item["mtime_ns"],
                            "sha256": digest.hexdigest()})
            start += size
    writer.write_dataset("original/data", "u1", (total,), total, chunks())
    writer.index["original"] = {"files": entries,
                               "directories": [{k: r[k] for k in ("path", "mode", "mtime_ns")}
                                               for r in directories]}


def write_particles(writer):
    entries = {entry["path"]: entry for entry in writer.index["original"]["files"]}
    names = list(entries)
    def raw(name):
        if name not in entries:
            raise ConversionError(f"Missing companion STL: {name}")
        return writer.read_original(entries[name])
    malformed = [name for name in names if "/" not in name and name.startswith("body")
                 and name.endswith(".info") and not BODY_RE.fullmatch(name)]
    if malformed:
        raise ConversionError(f"Unexpected particle-info filename: {malformed[0]}")
    info_names = sorted((name for name in names if BODY_RE.fullmatch(name)),
                        key=lambda name: int(BODY_RE.fullmatch(name)[1]))
    records = [parse_particle_info(raw(name), name) for name in info_names]
    ids = [r["body_id"] for r in records]
    if len(set(ids)) != len(ids):
        raise ConversionError("Duplicate body IDs in snapshot")
    expected = {f"stlFiles/{r['body_id']}.stl" for r in records if r["sphere_center"] is None}
    actual = {name for name in names if name.startswith("stlFiles/")
              and name.count("/") == 1 and name.endswith(".stl")}
    if actual != expected:
        missing, extra = sorted(expected - actual), sorted(actual - expected)
        detail = f"missing {missing[0]}" if missing else f"unmatched {extra[0]}"
        raise ConversionError(f"Particle-info/STL mismatch: {detail}")
    offsets = [0]
    for record in records:
        count = 0
        if record["sphere_center"] is None:
            name = f"stlFiles/{record['body_id']}.stl"
            count = stl_triangle_count(raw(name), name)
        offsets.append(offsets[-1] + count)
    positions, volumes, radii, types = array("d"), array("d"), array("d"), array("B")
    def geometry_chunks():
        for index, record in enumerate(records):
            if record["sphere_center"] is not None:
                types.append(0)
                positions.extend(record["sphere_center"])
                radii.append(record["radius"])
                volume = 4.0 * math.pi * record["radius"] ** 3 / 3.0
                if not math.isfinite(volume) or volume <= 0:
                    raise ConversionError(f"Invalid sphere volume: body{record['body_id']}.info")
                volumes.append(volume)
            else:
                types.append(1)
                radii.append(float("nan"))
                name = f"stlFiles/{record['body_id']}.stl"
                mesh, center, volume = parse_stl(raw(name), name)
                positions.extend(center)
                volumes.append(volume)
                if len(mesh) != 9 * (offsets[index + 1] - offsets[index]):
                    raise ConversionError(f"STL triangle count mismatch: {name}")
                buffer = little_endian_buffer(mesh, "<f8")
                for pos in range(0, len(buffer), BLOCK_BYTES):
                    yield buffer[pos:pos + BLOCK_BYTES]
                del buffer, mesh
    writer.write_dataset("geometry/triangles", "<f8", (offsets[-1], 3, 3),
                          offsets[-1] * 72, geometry_chunks())
    writer.write_array("geometry/triangle_offsets", offsets, "<u8")
    n = len(records)
    writer.write_array("particles/body_id", ids, "<i8")
    writer.write_strings("particles/body_name", [r["body_name"] for r in records])
    writer.write_array("particles/geometry_type", types, "u1")
    writer.write_array("particles/position", positions, "<f8", (n, 3))
    writer.write_array("particles/radius", radii, "<f8")
    writer.write_array("particles/volume", volumes, "<f8")
    for name in ("velocity", "axis"):
        writer.write_array(f"particles/{name}", (v for r in records for v in r[name]), "<f8", (n, 3))
    writer.write_array("particles/omega", [r["omega"] for r in records], "<f8")
    writer.write_array("particles/static", [r["static"] for r in records], "u1")
    writer.write_array("particles/time_steps_in_contact", [r["contact_steps"] for r in records], "<i8")
    return n


def finish_pending_deletions(bodies, data):
    if not bodies.exists():
        return
    for pending in sorted(bodies.glob(".packed-*")):
        name = pending.name[len(".packed-"):]
        if time_value(name) is None:
            continue
        archive = data / (name + ".bin")
        if pending.is_symlink() or archive.is_symlink() or not archive.is_file():
            raise ConversionError(f"Cannot verify interrupted deletion: {pending}")
        with ParticleDataFile(archive) as reader:
            reader.verify(name)
            compare_directory(pending, reader.metadata, subset=True)
        shutil.rmtree(pending)
        sync_directory(bodies)
        print(f"Finished interrupted cleanup: {name}", flush=True)


def pack_time(source, archive, keep_original):
    name = source.name
    if archive.is_symlink():
        raise ConversionError(f"Refusing symbolic link: {archive}")
    if archive.exists():
        with ParticleDataFile(archive) as reader:
            reader.verify(name)
            before = compare_directory(source, reader.metadata)
            count = reader.datasets["particles/body_id"]["shape"][0]
        if not keep_original:
            sync_file(archive)
            remove_verified_source(source, before)
        return count, "reused verified archive"
    before = inventory(source)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{name}.", suffix=".bin.tmp", dir=archive.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w+b") as stream:
            writer = BinaryWriter(stream, name)
            write_original(writer, source, before)
            count = write_particles(writer)
            writer.finish()
        with ParticleDataFile(temporary) as reader:
            reader.verify(name)
        if inventory_signature(before) != inventory_signature(inventory(source)):
            raise ConversionError(f"Source changed during conversion; kept: {source}")
        sync_file(temporary)
        os.link(temporary, archive)  # Atomic publication without overwrite.
        temporary.unlink()
        sync_directory(archive.parent)
        if not keep_original:
            remove_verified_source(source, before)
        return count, "written and verified"
    finally:
        if temporary.exists():
            temporary.unlink()


def unpack_time(archive, destination):
    name = archive.stem
    with ParticleDataFile(archive) as reader:
        reader.verify(name)
        meta = reader.metadata
        if destination.exists() or destination.is_symlink():
            compare_directory(destination, meta)
            return "already restored; identical contents"
        stage = Path(tempfile.mkdtemp(prefix=f".restore-{name}-", dir=destination.parent))
        try:
            for directory in meta["directories"]:
                (stage / directory).mkdir(parents=True, exist_ok=True)
            for index, relative in enumerate(meta["names"]):
                target = stage / relative
                digest = hashlib.sha256()
                with target.open("xb") as stream:
                    for chunk in reader.original_chunks(relative):
                        stream.write(chunk)
                        digest.update(chunk)
                if digest.hexdigest() != meta["hashes"][index]:
                    raise ConversionError(f"Checksum mismatch while restoring: {relative}")
                timestamp = int(meta["mtime_ns"][index])
                os.utime(target, ns=(timestamp, timestamp))
                os.chmod(target, int(meta["modes"][index]) & 0o7777)
            compare_directory(stage, meta)
            for index in sorted(range(len(meta["directories"])),
                                key=lambda i: meta["directories"][i].count("/"), reverse=True):
                target = stage / meta["directories"][index]
                timestamp = int(meta["directory_mtime_ns"][index])
                os.utime(target, ns=(timestamp, timestamp))
                os.chmod(target, int(meta["directory_modes"][index]) & 0o7777)
            if destination.exists() or destination.is_symlink():
                raise ConversionError(f"Destination appeared during restore: {destination}")
            stage.rename(destination)
            sync_directory(destination.parent)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    return "restored and verified"


def time_value(name):
    if not TIME_RE.fullmatch(name):
        return None
    value = Decimal(name)
    return value if value.is_finite() else None


def time_entries(root, archives=False):
    found = []
    for path in root.iterdir():
        if archives:
            if path.suffix != ".bin":
                continue
            name = path.stem
        else:
            name = path.name
        value = time_value(name)
        if value is None:
            continue
        if path.is_symlink():
            raise ConversionError(f"Refusing symbolic link: {path}")
        expected = path.is_file() if archives else path.is_dir()
        if not expected:
            raise ConversionError(f"Expected {'file' if archives else 'directory'}: {path}")
        found.append((value, name, path))
    found.sort()
    for previous, current in zip(found, found[1:]):
        if previous[0] == current[0]:
            raise ConversionError(f"Ambiguous duplicate times: {previous[1]} and {current[1]}")
    return [(name, path) for _, name, path in found]


def checked_relative(name, directory=False):
    if directory and name == "":
        return name
    parts = name.split("/")
    if (not name or "\\" in name or "\x00" in name or ":" in name
            or any(p in ("", ".", "..") for p in parts)):
        raise ConversionError(f"Unsafe stored path: {name!r}")
    return name


def inventory(root):
    """List every ordinary file and directory, including empty directories."""
    files, directories = [], []

    def visit(path, relative):
        status = path.lstat()
        record = {"path": relative, "size": status.st_size,
                  "mode": stat.S_IMODE(status.st_mode), "mtime_ns": status.st_mtime_ns,
                  "signature": (status.st_dev, status.st_ino, status.st_size,
                                status.st_mtime_ns, status.st_ctime_ns, status.st_mode)}
        if stat.S_ISDIR(status.st_mode):
            checked_relative(relative, directory=True)
            directories.append(record)
            for child in sorted(path.iterdir(), key=lambda p: p.name):
                visit(child, f"{relative}/{child.name}" if relative else child.name)
        elif stat.S_ISREG(status.st_mode):
            checked_relative(relative)
            files.append(record)
        else:
            raise ConversionError(f"Only regular files/directories are supported: {path}")

    visit(root, "")
    return files, directories


def inventory_signature(items):
    files, directories = items
    return tuple((r["path"], r["signature"]) for r in files + directories)


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(BLOCK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare_directory(path, meta, subset=False):
    before = inventory(path)
    files, directories = before
    stored = {name: index for index, name in enumerate(meta["names"])}
    file_names = {r["path"] for r in files}
    dir_names = {r["path"] for r in directories}
    if subset:
        matches = file_names <= set(stored) and dir_names <= set(meta["directories"])
    else:
        matches = file_names == set(stored) and dir_names == set(meta["directories"])
    if not matches:
        raise ConversionError(f"Existing directory differs from archive; kept untouched: {path}")
    for item in files:
        index = stored[item["path"]]
        size = int(meta["offsets"][index + 1]) - int(meta["offsets"][index])
        if (item["size"] != size
                or file_sha256(path / item["path"]) != meta["hashes"][index]):
            raise ConversionError(f"Existing file differs from archive; kept untouched: {path / item['path']}")
    if inventory_signature(before) != inventory_signature(inventory(path)):
        raise ConversionError(f"Source changed during verification: {path}")
    return before


def sync_file(path):
    with path.open("rb") as stream:
        os.fsync(stream.fileno())


def sync_directory(path):
    if os.name == "posix":
        fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(fd)
        finally:
            os.close(fd)


@contextmanager
def conversion_lock(case):
    """OS lock automatically releases even if the converter is interrupted."""
    path = case / ".compress_bodiesInfo.lock"
    if path.is_symlink():
        raise ConversionError(f"Refusing symbolic link: {path}")
    with path.open("a+b") as stream:
        stream.seek(0, os.SEEK_END)
        if stream.tell() == 0:
            stream.write(b"0")
            stream.flush()
        stream.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise ConversionError("Another conversion is already running in this case") from exc
        try:
            yield
        finally:
            if os.name == "nt":
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def remove_verified_source(source, expected_inventory):
    if inventory_signature(expected_inventory) != inventory_signature(inventory(source)):
        raise ConversionError(f"Source changed before deletion; kept: {source}")
    pending = source.parent / (".packed-" + source.name)
    if pending.exists():
        raise ConversionError(f"Pending deletion already exists: {pending}")
    source.rename(pending)
    sync_directory(source.parent)
    # The verified archive was published and fsynced before reaching this point.
    shutil.rmtree(pending)
    sync_directory(source.parent)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Pack bodiesInfo into one UNCOMPRESSED .bin file per time, or restore it exactly.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n  python3 compress_bodiesInfo.py\n"
               "  python3 compress_bodiesInfo.py decompress\n"
               "  python3 compress_bodiesInfo.py --keep-original\n"
               "  python3 compress_bodiesInfo.py --case /path/to/case\n\n"
               "Compress deletes each source time directory ONLY after verification.\n"
               "Decompress keeps particleData/*.bin. Run after the simulation has finished.")
    parser.add_argument("action", nargs="?", choices=("compress", "decompress"), default="compress")
    parser.add_argument("--case", default=".", help="case directory (default: current directory)")
    parser.add_argument("--keep-original", action="store_true", help="compress without deleting bodiesInfo")
    args = parser.parse_args(argv)
    if args.action == "decompress" and args.keep_original:
        parser.error("--keep-original is a compress option; decompress always keeps the binary files")
    case = Path(args.case).expanduser().resolve()
    if not case.is_dir():
        parser.error(f"Case directory not found: {case}")
    bodies, data = case / "bodiesInfo", case / "particleData"
    try:
        for path in (bodies, data):
            if path.is_symlink() or (path.exists() and not path.is_dir()):
                raise ConversionError(f"Expected an ordinary directory: {path}")
        with conversion_lock(case):
            finish_pending_deletions(bodies, data)
            if args.action == "compress":
                if not bodies.is_dir():
                    if data.is_dir() and time_entries(data, archives=True):
                        print("No bodiesInfo directory; particleData is already present. Nothing to pack.")
                        return 0
                    raise ConversionError(f"bodiesInfo directory not found: {bodies}")
                entries = time_entries(bodies)
                if not entries:
                    if (not any(bodies.iterdir()) and data.is_dir()
                            and time_entries(data, archives=True)):
                        if not args.keep_original:
                            bodies.rmdir()
                        print("All snapshots are already packed. Nothing more to pack.")
                        return 0
                    raise ConversionError(f"No numeric time directories in {bodies}")
                data.mkdir(exist_ok=True)
                for index, (name, source) in enumerate(entries, 1):
                    count, status = pack_time(source, data / (name + ".bin"), args.keep_original)
                    print(f"[{index}/{len(entries)}] {name}.bin: {count} particles; {status}; "
                          f"source {'kept' if args.keep_original else 'removed'}", flush=True)
                if not args.keep_original and not any(bodies.iterdir()):
                    bodies.rmdir()
                print(f"Done: {len(entries)} time snapshots in {data}")
            else:
                if not data.is_dir():
                    raise ConversionError(f"particleData directory not found: {data}")
                entries = time_entries(data, archives=True)
                if not entries:
                    raise ConversionError(f"No numeric .bin files in {data}")
                bodies.mkdir(exist_ok=True)
                for index, (name, archive) in enumerate(entries, 1):
                    status = unpack_time(archive, bodies / name)
                    print(f"[{index}/{len(entries)}] {name}: {status}", flush=True)
                print(f"Done: {len(entries)} time snapshots in {bodies}; binary files kept.")
        return 0
    except KeyboardInterrupt:
        print("\nInterrupted. Completed binary files are retained; rerun to continue.", file=sys.stderr)
        return 130
    except (OSError, ValueError, KeyError, TypeError, OverflowError, ConversionError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
