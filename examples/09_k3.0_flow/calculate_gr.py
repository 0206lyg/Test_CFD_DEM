#!/usr/bin/env python3
"""Calculate a 3-D radial distribution function from particleData/*.bin.

Python 3.8+ standard library only. This is a standalone reader for the v1
binary format written by compress_bodiesInfo.py; no decompression is needed.
Input archives are never changed. Only existing, nonnegative INTEGER times
are analyzed: a 0.15 s output interval gives 0, 3, 6, ... s, not 1 or 2 s.
No interpolation, nearest-time substitution, or averaging over time is done.

Run in the case directory, with no arguments:
    python3 calculate_gr.py

The defaults use the same upper reservoir as upper_reservoir_solid_fraction.py:
x,z = [-0.020, 0.020] m; y = [0.025, 0.075] m; volume = 8.0e-5 m^3;
no periodic axes. Both endpoints are included with the same 1e-12 m tolerance.
The default r_max is 0.020 m and dr is 0.0001 m (200 distance bins).
Settings can be edited in the DEFAULT_* constants below or overridden by CLI.

Optional --bounds specifies a rectangular observation window for particle
centers in the order xmin xmax ymin ymax zmin zmax. Nonperiodic axes retain
low - BOUND_EPS <= position <= high + BOUND_EPS.
Periodic bounds must span one FULL physical period; coordinates are wrapped
on those axes. Do not use the bounding box of a channel containing excluded
solid regions as though its whole volume were available to particle centers.

The output defaults to postProcessing/radialDistribution.dat under --case:
    # time_s 0 3 6
    # particle_count 1200 1190 1170
    # columns r_m g_0s g_3s g_6s
    <shell midpoint in m> <g at 0 s> <g at 3 s> <g at 6 s>
All times share the same shell edges. g is dimensionless; r is a physical
center-to-center distance, NOT r/diameter or particle-surface separation.
The final shell can be narrower than --dr. N < 2 yields nan, not zero.

Estimator, counting each unordered pair once:
    g[k] = 2 V sum(weight_ij) / (N (N-1) shell_volume[k])
    shell_volume[k] = 4*pi/3 * (edge[k+1]**3 - edge[k]**3)
    weight_ij = product(L[a] / (L[a] - abs(delta[a]))) over open axes
Periodic axes use minimum-image distances and unit edge weights.
This translation edge correction compensates for missing pairs outside a
rectangular observation window. It does not remove physical density gradients,
wall layering, or inaccessible solids. Independent uniform points have an
expected g=1, including the finite-N correction. No smoothing is applied.
The r limit is half the shortest box side, keeping complete spherical shells
and bounded edge weights. Cell lists avoid allocating an N-by-N distance array.

Position convention is exactly particles/position in the archive. For spheres
this is the center from body.info; the current compressor uses the arithmetic
mean of unique STL vertices for STL particles, not the volume centroid. All
stored particles whose positions pass the window selection are included,
including static particles. Particle sizes and orientations are not needed
for this pooled center-position RDF.

Method background (no dependency on the referenced software):
https://rdrr.io/cran/spatstat.explore/man/pcf3est.html
https://rdrr.io/cran/spatstat.explore/man/edge.Trans.html
"""

import argparse
from array import array
from bisect import bisect_right
from decimal import Decimal, InvalidOperation
import hashlib
from itertools import combinations, product
import json
import math
import os
from pathlib import Path
import struct
import sys
import tempfile


# User settings: running `python3 calculate_gr.py` uses these defaults.
DEFAULT_BOUNDS = (-0.020, 0.020, 0.025, 0.075, -0.020, 0.020)
DEFAULT_PERIODIC = "none"
DEFAULT_R_MAX = None  # None: half the shortest box side (0.020 m here).
DEFAULT_DR = None  # None: r_max / DEFAULT_RADIAL_BINS (0.0001 m here).
DEFAULT_RADIAL_BINS = 200
DEFAULT_OUTPUT = "postProcessing/radialDistribution.dat"
BOUND_EPS = 1.0e-12  # Same inclusive-boundary tolerance as upper reservoir script.


HEADER = struct.Struct("<16sI")
FOOTER = struct.Struct("<16sQQ32s")
HEADER_MAGIC = b"PDBIN-HEADER-v1!"
FOOTER_MAGIC = b"PDBIN-FOOTER-v1!"
FORMAT = "openHFDIB-DEM particleData binary"
NEIGHBOR_OFFSETS = tuple(product((-1, 0, 1), repeat=3))


def read_positions(path):
    """Read and checksum only the index and positions, not STL/original data."""
    try:
        with path.open("rb") as stream:
            size = os.fstat(stream.fileno()).st_size
            if size < HEADER.size + FOOTER.size:
                raise ValueError("truncated binary archive")
            if HEADER.unpack(stream.read(HEADER.size)) != (HEADER_MAGIC, 1):
                raise ValueError("unsupported particleData header/version")
            stream.seek(-FOOTER.size, os.SEEK_END)
            magic, offset, length, digest = FOOTER.unpack(stream.read(FOOTER.size))
            if (magic != FOOTER_MAGIC or offset < HEADER.size
                    or offset + length != size - FOOTER.size):
                raise ValueError("invalid index footer")
            stream.seek(offset)
            index_bytes = stream.read(length)
            if hashlib.sha256(index_bytes).digest() != digest:
                raise ValueError("index SHA256 mismatch")
            index = json.loads(index_bytes)
            if not isinstance(index, dict):
                raise ValueError("index must be a JSON object")
            if index.get("format") != FORMAT or index.get("schema_version") != 1:
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
            item = datasets.get("particles/position")
            if not isinstance(item, dict):
                raise ValueError("missing particles/position dataset")
            shape = item.get("shape")
            if (item.get("dtype") != "<f8" or not isinstance(shape, list)
                    or len(shape) != 2 or shape[1] != 3
                    or any(type(v) is not int or v < 0 for v in shape)):
                raise ValueError("positions must have dtype <f8 and shape [N, 3]")
            n = shape[0]
            body_ids = datasets.get("particles/body_id")
            if not isinstance(body_ids, dict) or body_ids.get("shape") != [n]:
                raise ValueError("body_id and position particle counts differ")
            start, count = item.get("offset"), item.get("size")
            if (type(start) is not int or type(count) is not int
                    or start < HEADER.size or count != 24 * n
                    or start + count > offset):
                raise ValueError("invalid position dataset offset/size")
            stream.seek(start)
            payload = stream.read(count)
            if len(payload) != count or hashlib.sha256(payload).hexdigest() != item.get("sha256"):
                raise ValueError("position dataset SHA256 mismatch")
            values = array("d")
            if values.itemsize != 8:
                raise ValueError("this Python platform does not provide 8-byte doubles")
            values.frombytes(payload)
            if sys.byteorder != "little":
                values.byteswap()
            if not all(math.isfinite(v) for v in values):
                raise ValueError("position dataset contains nonfinite coordinates")
            return values
    except (OSError, ValueError, TypeError, struct.error, OverflowError) as exc:
        raise ValueError("{}: {}".format(path, exc)) from exc


def integer_snapshots(directory):
    """Select by exact decimal filename time, before opening any archive."""
    if not directory.is_dir():
        raise ValueError("particleData directory does not exist: {}".format(directory))
    selected = {}
    skipped = 0
    for path in directory.glob("*.bin"):
        if not path.is_file():
            continue
        try:
            time = Decimal(path.stem)
        except InvalidOperation:
            skipped += 1
            continue
        if not time.is_finite() or time < 0 or time != time.to_integral_value():
            skipped += 1
            continue
        if time in selected:
            raise ValueError("duplicate integer time: {} and {}".format(selected[time].name, path.name))
        selected[time] = path
    if not selected:
        raise ValueError("no existing nonnegative integer-time .bin files in {}".format(directory))
    return [(str(int(time)), selected[time]) for time in sorted(selected)], skipped


def select_points(values, lower, upper, lengths, periodic):
    """Return selected positions relative to the lower bounds, plus total N."""
    points = []
    for i in range(0, len(values), 3):
        point = values[i:i + 3]
        if any(not periodic[a] and not lower[a] - BOUND_EPS <= point[a] <= upper[a] + BOUND_EPS
               for a in range(3)):
            continue
        relative = tuple((point[a] - lower[a]) % lengths[a] if periodic[a]
                         else point[a] - lower[a] for a in range(3))
        points.append(relative)
    return points, len(values) // 3


def shell_edges(r_max, dr):
    ratio = r_max / dr
    if not math.isfinite(ratio) or ratio > 1000000:
        raise ValueError("--dr is too small: use at most 1,000,000 distance bins")
    nearest = round(ratio)
    if nearest > 0 and math.isclose(ratio, nearest, rel_tol=1e-12):
        bins = nearest
    else:
        bins = max(1, math.ceil(ratio))
    edges = [i * dr for i in range(bins)] + [r_max]
    if any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("distance bin edges are not strictly increasing")
    return edges


def radial_distribution(points, lengths, periodic, edges):
    """Translation-corrected histogram using cells of width >= r_max."""
    n = len(points)
    if n < 2:
        return [math.nan] * (len(edges) - 1), 0
    r_max = edges[-1]
    r_max_sq = r_max * r_max
    ncell = tuple(max(1, int(length / r_max)) for length in lengths)
    widths = tuple(lengths[a] / ncell[a] for a in range(3))
    cells = {}
    for point in points:
        key = tuple(max(0, min(ncell[a] - 1, int(point[a] / widths[a]))) for a in range(3))
        cells.setdefault(key, []).append(point)

    counts = [0.0] * (len(edges) - 1)
    pair_count = 0
    lx, ly, lz = lengths
    px, py, pz = periodic
    for key, local in cells.items():
        neighbors = set()
        for delta in NEIGHBOR_OFFSETS:
            other = tuple((key[a] + delta[a]) % ncell[a] if periodic[a]
                          else key[a] + delta[a] for a in range(3))
            if other >= key and other in cells:
                neighbors.add(other)
        # Deduplicate wrapped cells (particularly when there are only 2 cells).
        for other in sorted(neighbors):
            pairs = combinations(local, 2) if other == key else product(local, cells[other])
            for p, q in pairs:
                dx, dy, dz = abs(p[0] - q[0]), abs(p[1] - q[1]), abs(p[2] - q[2])
                if px:
                    dx = min(dx, lx - dx)
                if py:
                    dy = min(dy, ly - dy)
                if pz:
                    dz = min(dz, lz - dz)
                distance_sq = dx * dx + dy * dy + dz * dz
                if distance_sq >= r_max_sq:
                    continue
                k = bisect_right(edges, math.sqrt(distance_sq)) - 1
                if k >= len(counts):  # sqrt can round up to r_max at the limit
                    continue
                weight = 1.0
                if not px:
                    weight *= lx / (lx - dx)
                if not py:
                    weight *= ly / (ly - dy)
                if not pz:
                    weight *= lz / (lz - dz)
                counts[k] += weight
                pair_count += 1
    factor = 2.0 * (lx * ly * lz) / (n * (n - 1))
    result = [factor * count / ((4.0 * math.pi / 3.0) * (hi**3 - lo**3))
              for count, lo, hi in zip(counts, edges, edges[1:])]
    return result, pair_count


def write_output(path, snapshots, counts, totals, pairs, columns, bounds, periodic, edges, dr):
    """Write complete data atomically; failures leave an existing result intact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                         dir=path.parent, prefix=path.name + ".", suffix=".tmp",
                                         delete=False) as stream:
            temporary = Path(stream.name)
            stream.write("# Radial distribution from particleData v1; g dimensionless, r in m\n")
            stream.write("# Selection: existing nonnegative integer times only; no interpolation\n")
            stream.write("# bounds_m xmin xmax ymin ymax zmin zmax: " + " ".join(format(v, ".17g") for v in bounds) + "\n")
            stream.write("# periodic_axes " + periodic + "\n")
            stream.write("# Open-axis selection: low-eps <= position <= high+eps; eps_m={:.17g}; periodic coordinates wrapped\n".format(BOUND_EPS))
            stream.write("# positions: archive particles/position; all stored particle types included\n")
            stream.write("# edge_correction translation on open axes; minimum image on periodic axes\n")
            stream.write("# g=2*V*sum_pair_weights/(N*(N-1)*shell_volume); N<2 gives nan\n")
            stream.write("# dr_m {:.17g}; r_max_m {:.17g}; final shell may be narrower\n".format(dr, edges[-1]))
            stream.write("# shell_edges_m " + " ".join(format(v, ".17g") for v in edges) + "\n")
            stream.write("# time_s " + " ".join(time for time, _ in snapshots) + "\n")
            stream.write("# source_files " + " ".join(source.name for _, source in snapshots) + "\n")
            stream.write("# particle_count " + " ".join(map(str, counts)) + "\n")
            stream.write("# source_particle_count " + " ".join(map(str, totals)) + "\n")
            stream.write("# unordered_pairs_below_r_max " + " ".join(map(str, pairs)) + "\n")
            stream.write("# columns r_m " + " ".join("g_{}s".format(time) for time, _ in snapshots) + "\n")
            for k, (lo, hi) in enumerate(zip(edges, edges[1:])):
                row = [(lo + hi) / 2.0] + [column[k] for column in columns]
                stream.write("\t".join(format(value, ".12g") for value in row) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0],
                                     epilog="See the script docstring for the estimator and full usage notes.")
    parser.add_argument("--case", type=Path, default=Path("."), help="case containing particleData (default: current directory)")
    parser.add_argument("--bounds", type=float, nargs=6, default=DEFAULT_BOUNDS,
                        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
                        help="center observation window, m (default: upper reservoir, -0.020 0.020 0.025 0.075 -0.020 0.020)")
    parser.add_argument("--periodic", default=DEFAULT_PERIODIC, choices=("none", "x", "y", "z", "xy", "xz", "yz", "xyz"),
                        help="periodic axes (default: none); periodic bounds must span their full period")
    parser.add_argument("--r-max", type=float, default=DEFAULT_R_MAX, help="maximum center distance, m (default: half the shortest box side)")
    parser.add_argument("--dr", type=float, default=DEFAULT_DR, help="radial bin width, m (default: r-max/200)")
    parser.add_argument("--output", type=Path, default=Path(DEFAULT_OUTPUT),
                        help="output file; relative paths are under --case")
    args = parser.parse_args(argv)
    try:
        lower, upper = tuple(args.bounds[::2]), tuple(args.bounds[1::2])
        lengths = tuple(hi - lo for lo, hi in zip(lower, upper))
        if (not all(math.isfinite(v) for v in args.bounds)
                or not all(math.isfinite(v) and v > 0 for v in lengths)):
            raise ValueError("--bounds must be finite and each maximum must exceed its minimum")
        volume = lengths[0] * lengths[1] * lengths[2]
        if not math.isfinite(volume) or volume <= 0:
            raise ValueError("box volume is outside the supported floating-point range")
        limit = 0.5 * min(lengths)
        r_max = limit if args.r_max is None else args.r_max
        if not math.isfinite(r_max) or r_max <= 0 or r_max > limit:
            raise ValueError("--r-max must be positive and at most half the shortest box side ({:.12g} m)".format(limit))
        dr = r_max / DEFAULT_RADIAL_BINS if args.dr is None else args.dr
        if not math.isfinite(dr) or dr <= 0:
            raise ValueError("--dr must be finite and positive")
        edges = shell_edges(r_max, dr)
        if any(not math.isfinite(b**3 - a**3) or b**3 - a**3 <= 0 for a, b in zip(edges, edges[1:])):
            raise ValueError("shell volumes are outside the supported floating-point range")
        periodic = tuple(axis in args.periodic for axis in "xyz")
        case = args.case.resolve()
        data_dir = case / "particleData"
        snapshots, skipped = integer_snapshots(data_dir)
        output = (args.output if args.output.is_absolute() else case / args.output).resolve()
        if output == data_dir.resolve() or data_dir.resolve() in output.parents:
            raise ValueError("--output must be outside particleData to preserve input archives")
        print("Selected times (s): " + ", ".join(time for time, _ in snapshots), flush=True)
        print("Skipped {} noninteger/non-time .bin files; r_max={:.8g} m, dr={:.8g} m".format(skipped, r_max, dr), flush=True)
        columns, counts, totals, pair_counts = [], [], [], []
        for time, source in snapshots:
            values = read_positions(source)
            points, total = select_points(values, lower, upper, lengths, periodic)
            del values
            print("t={} s: {}/{} particles in window; calculating...".format(time, len(points), total), flush=True)
            gr, pairs = radial_distribution(points, lengths, periodic, edges)
            columns.append(gr)
            counts.append(len(points))
            totals.append(total)
            pair_counts.append(pairs)
            print("  {} pairs below r_max{}".format(pairs, "; g(r)=nan because N<2" if len(points) < 2 else ""), flush=True)
        write_output(output, snapshots, counts, totals, pair_counts, columns,
                     args.bounds, args.periodic, edges, dr)
        print("Wrote {} ({} r rows x {} time columns)".format(output, len(edges) - 1, len(snapshots)), flush=True)
        return 0
    except (OSError, ValueError, OverflowError) as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
