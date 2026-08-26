#!/usr/bin/env python3
"""Center-assigned solid volume fraction in the 08/09 upper reservoir.

The physical geometry is fixed to the current examples/08 and examples/09:

    x, z = -0.020 ... 0.020 m
    y    =  0.025 ... 0.075 m

The reservoir is divided uniformly along y.  A particle contributes its full
volume to the single bin containing its center.  Sphere centers/radii are read
from body*.info.  For STL bodies, the center is the arithmetic mean of the
unique STL vertices, matching openHFDIB-DEM's stlBased::getCoM().

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import statistics
import struct
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence


X_MIN, X_MAX = -0.020, 0.020
Y_MIN, Y_MAX = 0.025, 0.075
Z_MIN, Z_MAX = -0.020, 0.020
DEFAULT_BINS = 5
BOUND_EPS = 1.0e-12
RESERVOIR_VOLUME_M3 = 8.0e-5

FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
BODY_FILE_RE = re.compile(r"^body(-?\d+)\.info$")
BODY_ID_RE = re.compile(r"^\s*bodyId\s+(-?\d+)\s*;", re.MULTILINE)
BODY_NAME_RE = re.compile(r"^\s*bodyName\s+\S+\s*;", re.MULTILINE)
VEL_RE = re.compile(r"^\s*Vel\s+\([^)]*\)\s*;", re.MULTILINE)
OMEGA_RE = re.compile(r"^\s*omega\s+" + FLOAT_PATTERN + r"\s*;", re.MULTILINE)
AXIS_RE = re.compile(r"^\s*Axis\s+\([^)]*\)\s*;", re.MULTILINE)
STATIC_RE = re.compile(r"^\s*static\s+\S+\s*;", re.MULTILINE)
CONTACT_STEPS_RE = re.compile(
    r"^\s*timeStepsInContWStatic\s+-?\d+\s*;", re.MULTILINE
)
SPHERE_BLOCK_RE = re.compile(r"\bsphere\s*\{(.*?)\}", re.DOTALL)
POSITION_RE = re.compile(r"\bposition\s*\(([^)]*)\)\s*;")
RADIUS_RE = re.compile(r"\bradius\s+(" + FLOAT_PATTERN + r")\s*;")
ASCII_VERTEX_RE = re.compile(
    r"^\s*vertex\s+(" + FLOAT_PATTERN + r")\s+(" + FLOAT_PATTERN
    + r")\s+(" + FLOAT_PATTERN + r")\s*$",
    re.IGNORECASE | re.MULTILINE,
)
ASCII_END_RE = re.compile(r"^\s*endsolid\b", re.IGNORECASE | re.MULTILINE)

REQUIRED_INFO_PATTERNS = (
    BODY_NAME_RE,
    VEL_RE,
    OMEGA_RE,
    AXIS_RE,
    STATIC_RE,
    CONTACT_STEPS_RE,
)


class SnapshotError(RuntimeError):
    """A bodiesInfo time directory is incomplete or malformed."""


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError("value must be a finite positive number")
    return number


def nonnegative_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise argparse.ArgumentTypeError("value must be a finite non-negative number")
    return number


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return number


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute center-assigned solid volume fractions in the fixed "
            "examples/08 and examples/09 upper reservoir."
        )
    )
    parser.add_argument(
        "case",
        nargs="?",
        default=".",
        help="OpenFOAM case directory (default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="postProcessing/upperReservoirSolidFraction.dat",
        help="output path, relative to CASE unless absolute",
    )
    parser.add_argument(
        "--bins",
        type=positive_int,
        default=DEFAULT_BINS,
        help="number of equal y intervals (default: 5)",
    )
    volume_group = parser.add_mutually_exclusive_group()
    volume_group.add_argument(
        "--particle-volume-m3",
        type=positive_float,
        help="use this monodisperse particle volume instead of auto-detection",
    )
    volume_group.add_argument(
        "--equivalent-diameter-mm",
        type=positive_float,
        help="set particle volume to pi*d_eq^3/6 using d_eq in mm",
    )
    volume_group.add_argument(
        "--k",
        type=positive_float,
        help=(
            "set d_eq = 10 mm/k and particle volume = pi*d_eq^3/6, "
            "matching the current 08/09 submission scripts"
        ),
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help=(
            "conservatively skip the newest bodiesInfo time directory; use "
            "while the solver is running"
        ),
    )
    parser.add_argument(
        "--settle-seconds",
        type=nonnegative_float,
        default=2.0,
        help=(
            "skip a snapshot whose files changed more recently than this "
            "wall-clock interval (default: 2)"
        ),
    )
    parser.add_argument(
        "--volume-tolerance",
        type=positive_float,
        default=1.0e-5,
        help=(
            "maximum relative auto-detected volume spread allowed for a "
            "monodisperse run (default: 1e-5)"
        ),
    )
    return parser.parse_args(argv)


def numeric_time_directories(root: Path) -> list[tuple[float, Path]]:
    result: list[tuple[float, Path]] = []
    try:
        entries = list(root.iterdir())
    except OSError as exc:
        raise SystemExit(f"ERROR: cannot read {root}: {exc}") from exc

    for entry in entries:
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if math.isfinite(value):
            result.append((value, entry))
    result.sort(key=lambda item: (item[0], item[1].name))
    return result


def snapshot_signature(directory: Path) -> tuple[tuple[str, int, int], ...]:
    """Return names, sizes and mtimes for files relevant to one snapshot."""
    paths = [directory]
    stl_directory = directory / "stlFiles"
    if stl_directory.exists():
        paths.append(stl_directory)
    paths.extend(sorted(directory.glob("body*.info"), key=lambda path: path.name))
    if stl_directory.is_dir():
        paths.extend(sorted(stl_directory.glob("*.stl"), key=lambda path: path.name))

    signature: list[tuple[str, int, int]] = []
    for path in paths:
        try:
            status = path.stat()
        except OSError as exc:
            raise SnapshotError(f"cannot stat {path}: {exc}") from exc
        signature.append(
            (str(path.relative_to(directory)), status.st_size, status.st_mtime_ns)
        )
    return tuple(signature)


def require_settled(
    signature: Sequence[tuple[str, int, int]], settle_seconds: float
) -> None:
    if not signature or settle_seconds == 0.0:
        return
    newest_mtime_ns = max(item[2] for item in signature)
    age_seconds = (time.time_ns() - newest_mtime_ns) / 1.0e9
    if age_seconds < settle_seconds:
        raise SnapshotError(
            f"snapshot may still be writing (newest item is {age_seconds:.3g} s old)"
        )


def vector_from_text(value: str, source: Path) -> tuple[float, float, float]:
    parts = value.split()
    if len(parts) != 3:
        raise SnapshotError(f"expected a three-component vector in {source}")
    try:
        vector = tuple(float(part) for part in parts)
    except ValueError as exc:
        raise SnapshotError(f"invalid vector in {source}") from exc
    if not all(math.isfinite(component) for component in vector):
        raise SnapshotError(f"non-finite vector in {source}")
    return vector  # type: ignore[return-value]


def read_body_info(path: Path) -> tuple[int, tuple[float, float, float] | None, float | None]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise SnapshotError(f"cannot read {path}: {exc}") from exc

    body_id_match = BODY_ID_RE.search(text)
    if body_id_match is None:
        raise SnapshotError(f"missing bodyId in {path}")
    for pattern in REQUIRED_INFO_PATTERNS:
        if pattern.search(text) is None:
            raise SnapshotError(f"incomplete body record in {path}")

    body_id = int(body_id_match.group(1))
    file_match = BODY_FILE_RE.match(path.name)
    if file_match is None or int(file_match.group(1)) != body_id:
        raise SnapshotError(f"bodyId does not match filename in {path}")

    sphere_match = SPHERE_BLOCK_RE.search(text)
    if sphere_match is None:
        return body_id, None, None

    sphere_text = sphere_match.group(1)
    position_match = POSITION_RE.search(sphere_text)
    radius_match = RADIUS_RE.search(sphere_text)
    if position_match is None or radius_match is None:
        raise SnapshotError(f"incomplete sphere block in {path}")

    center = vector_from_text(position_match.group(1), path)
    radius = float(radius_match.group(1))
    if not math.isfinite(radius) or radius <= 0.0:
        raise SnapshotError(f"invalid sphere radius in {path}")
    return body_id, center, (4.0 * math.pi * radius**3) / 3.0


def binary_stl_triangles(data: bytes, source: Path) -> list[tuple[tuple[float, float, float], ...]] | None:
    if len(data) < 84:
        return None
    triangle_count = struct.unpack_from("<I", data, 80)[0]
    expected_size = 84 + 50 * triangle_count
    if expected_size != len(data):
        return None
    if triangle_count == 0:
        raise SnapshotError(f"binary STL has no triangles: {source}")

    record = struct.Struct("<12fH")
    triangles: list[tuple[tuple[float, float, float], ...]] = []
    offset = 84
    for _ in range(triangle_count):
        values = record.unpack_from(data, offset)
        triangles.append(
            (
                (values[3], values[4], values[5]),
                (values[6], values[7], values[8]),
                (values[9], values[10], values[11]),
            )
        )
        offset += record.size
    return triangles


def ascii_stl_triangles(data: bytes, source: Path) -> list[tuple[tuple[float, float, float], ...]]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise SnapshotError(f"unrecognized STL encoding: {source}") from exc
    if ASCII_END_RE.search(text) is None:
        raise SnapshotError(f"incomplete ASCII STL (missing endsolid): {source}")

    vertices = [
        (float(match.group(1)), float(match.group(2)), float(match.group(3)))
        for match in ASCII_VERTEX_RE.finditer(text)
    ]
    if not vertices or len(vertices) % 3 != 0:
        raise SnapshotError(f"invalid ASCII STL triangle data: {source}")
    if not all(math.isfinite(value) for vertex in vertices for value in vertex):
        raise SnapshotError(f"non-finite STL vertex in {source}")
    return [tuple(vertices[index : index + 3]) for index in range(0, len(vertices), 3)]


def determinant(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def read_stl_center_and_volume(path: Path) -> tuple[tuple[float, float, float], float]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise SnapshotError(f"cannot read {path}: {exc}") from exc

    triangles = binary_stl_triangles(data, path)
    if triangles is None:
        triangles = ascii_stl_triangles(data, path)

    unique_vertices = {vertex for triangle in triangles for vertex in triangle}
    if not unique_vertices:
        raise SnapshotError(f"STL has no vertices: {path}")
    count = len(unique_vertices)
    center = tuple(
        math.fsum(vertex[axis] for vertex in unique_vertices) / count
        for axis in range(3)
    )

    # Shifting the tetrahedron reference to the vertex-average center reduces
    # cancellation from the absolute particle position.
    signed_six_volume = 0.0
    for triangle in triangles:
        shifted = tuple(
            tuple(vertex[axis] - center[axis] for axis in range(3))
            for vertex in triangle
        )
        signed_six_volume += determinant(shifted[0], shifted[1], shifted[2])
    volume = abs(signed_six_volume) / 6.0
    if not math.isfinite(volume) or volume <= 0.0:
        raise SnapshotError(f"STL is not a closed, consistently oriented solid: {path}")
    return center, volume


def particle_bin(center: tuple[float, float, float], bins: int) -> int | None:
    x_value, y_value, z_value = center
    if not (X_MIN - BOUND_EPS <= x_value <= X_MAX + BOUND_EPS):
        return None
    if not (Y_MIN - BOUND_EPS <= y_value <= Y_MAX + BOUND_EPS):
        return None
    if not (Z_MIN - BOUND_EPS <= z_value <= Z_MAX + BOUND_EPS):
        return None

    y_value = min(max(y_value, Y_MIN), Y_MAX)
    width = (Y_MAX - Y_MIN) / bins
    # The small positive shift makes decimal values that are mathematically on
    # an internal edge (for example 0.030 m) enter the upper interval despite
    # binary floating-point representation noise.
    index = math.floor((y_value - Y_MIN + BOUND_EPS) / width)
    return min(index, bins - 1)


def parse_snapshot(
    directory: Path, bins: int, settle_seconds: float
) -> tuple[list[int], list[float]]:
    before = snapshot_signature(directory)
    require_settled(before, settle_seconds)

    stl_directory = directory / "stlFiles"
    if not stl_directory.is_dir():
        raise SnapshotError(f"missing stlFiles directory in {directory}")

    info_files: list[tuple[int, Path]] = []
    for path in directory.glob("body*.info"):
        match = BODY_FILE_RE.match(path.name)
        if match is None:
            raise SnapshotError(f"unexpected body-info filename: {path}")
        info_files.append((int(match.group(1)), path))
    info_files.sort(key=lambda item: item[0])

    bin_counts = [0] * bins
    detected_volumes: list[float] = []
    expected_stl_names: set[str] = set()
    for filename_id, info_path in info_files:
        body_id, sphere_center, sphere_volume = read_body_info(info_path)
        if body_id != filename_id:
            raise SnapshotError(f"bodyId mismatch in {info_path}")

        if sphere_center is not None and sphere_volume is not None:
            center, volume = sphere_center, sphere_volume
        else:
            stl_name = f"{body_id}.stl"
            expected_stl_names.add(stl_name)
            stl_path = stl_directory / stl_name
            if not stl_path.is_file():
                raise SnapshotError(f"missing companion STL for {info_path}")
            center, volume = read_stl_center_and_volume(stl_path)

        detected_volumes.append(volume)
        index = particle_bin(center, bins)
        if index is not None:
            bin_counts[index] += 1

    actual_stl_names = {path.name for path in stl_directory.glob("*.stl")}
    if actual_stl_names != expected_stl_names:
        raise SnapshotError(
            f"body-info/STL set mismatch in {directory}; snapshot may be incomplete"
        )

    after = snapshot_signature(directory)
    if before != after:
        raise SnapshotError(f"snapshot changed while it was being read: {directory}")
    return bin_counts, detected_volumes


def choose_particle_volume(
    arguments: argparse.Namespace, detected: Sequence[float]
) -> tuple[float, str]:
    if arguments.particle_volume_m3 is not None:
        return arguments.particle_volume_m3, "command_line_volume"
    if arguments.equivalent_diameter_mm is not None:
        diameter = arguments.equivalent_diameter_mm * 1.0e-3
        return math.pi * diameter**3 / 6.0, "command_line_equivalent_diameter"
    if arguments.k is not None:
        diameter = 0.010 / arguments.k
        return math.pi * diameter**3 / 6.0, "command_line_k"
    if not detected:
        return 0.0, "no_particles_detected"

    particle_volume = statistics.fmean(detected)
    volume_min = min(detected)
    volume_max = max(detected)
    relative_spread = (volume_max - volume_min) / particle_volume
    if relative_spread > arguments.volume_tolerance:
        raise SystemExit(
            "ERROR: auto-detected particle volumes are not monodisperse: "
            f"min={volume_min:.16g}, max={volume_max:.16g}, "
            f"relative spread={relative_spread:.6g}. Supply the intended "
            "monodisperse volume explicitly or analyze each size separately."
        )
    return particle_volume, "auto_detected_from_bodiesInfo"


def output_path(case_directory: Path, requested: str) -> Path:
    path = Path(requested).expanduser()
    return path if path.is_absolute() else case_directory / path


def bin_edges(bins: int) -> list[float]:
    width = (Y_MAX - Y_MIN) / bins
    return [Y_MIN + index * width for index in range(bins + 1)]


def write_output(
    path: Path,
    rows: Iterable[tuple[float, Sequence[int]]],
    bins: int,
    particle_volume: float,
    volume_source: str,
) -> None:
    edges = bin_edges(bins)
    reservoir_volume = RESERVOIR_VOLUME_M3
    bin_volume = reservoir_volume / bins
    column_names = ["time_s", "phi_upper"] + [
        f"phi_y{1000.0 * edges[index]:g}_{1000.0 * edges[index + 1]:g}mm"
        for index in range(bins)
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as output:
            output.write("# Center-assigned solid volume fraction\n")
            output.write(
                "# reservoir_bounds_m "
                f"x=[{X_MIN:.16g},{X_MAX:.16g}] "
                f"y=[{Y_MIN:.16g},{Y_MAX:.16g}] "
                f"z=[{Z_MIN:.16g},{Z_MAX:.16g}]\n"
            )
            output.write(
                f"# bins={bins} bin_volume_m3={bin_volume:.12g} "
                f"reservoir_volume_m3={reservoir_volume:.12g}\n"
            )
            output.write("# bin_order=throat_to_top\n")
            output.write(
                "# bin_assignment=lower_edge_inclusive_upper_edge_exclusive; "
                "final_upper_edge_inclusive\n"
            )
            output.write(
                f"# particle_volume_m3={particle_volume:.16g} "
                f"source={volume_source}\n"
            )
            output.write("# " + " ".join(column_names) + "\n")

            for time_value, counts in rows:
                fractions = [count * particle_volume / bin_volume for count in counts]
                whole_fraction = sum(counts) * particle_volume / reservoir_volume
                values = [time_value, whole_fraction, *fractions]
                output.write(" ".join(f"{value:.16g}" for value in values) + "\n")
        os.replace(temporary, path)
    except Exception:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_args(argv)
    case_directory = Path(arguments.case).expanduser().resolve()
    bodies_root = case_directory / "bodiesInfo"
    if not bodies_root.is_dir():
        print(f"ERROR: bodiesInfo directory not found: {bodies_root}", file=sys.stderr)
        return 2

    time_directories = numeric_time_directories(bodies_root)
    if not time_directories:
        print(f"ERROR: no numeric time directories found in {bodies_root}", file=sys.stderr)
        return 2
    for previous, current in zip(time_directories, time_directories[1:]):
        if previous[0] == current[0]:
            print(
                "ERROR: duplicate numeric time value in bodiesInfo: "
                f"{previous[1].name!r} and {current[1].name!r}",
                file=sys.stderr,
            )
            return 2
    if arguments.live:
        skipped_live = time_directories.pop()
        print(
            f"Skipping newest snapshot {skipped_live[1]} in --live mode.",
            file=sys.stderr,
        )

    rows: list[tuple[float, list[int]]] = []
    detected_volumes: list[float] = []
    skipped = 0
    for time_value, directory in time_directories:
        try:
            counts, volumes = parse_snapshot(
                directory, arguments.bins, arguments.settle_seconds
            )
        except SnapshotError as exc:
            skipped += 1
            print(f"WARNING: skipping {directory}: {exc}", file=sys.stderr)
            continue
        rows.append((time_value, counts))
        detected_volumes.extend(volumes)

    if not rows:
        print(
            "ERROR: no complete snapshots were available. If the simulation has "
            "finished, wait briefly or use --settle-seconds 0.",
            file=sys.stderr,
        )
        return 2

    particle_volume, volume_source = choose_particle_volume(
        arguments, detected_volumes
    )
    destination = output_path(case_directory, arguments.output)
    try:
        write_output(
            destination,
            rows,
            arguments.bins,
            particle_volume,
            volume_source,
        )
    except OSError as exc:
        print(f"ERROR: cannot write {destination}: {exc}", file=sys.stderr)
        return 2

    edges = bin_edges(arguments.bins)
    print(f"Wrote {len(rows)} time rows to {destination}")
    print(
        f"Upper reservoir: x,z = [{X_MIN:g}, {X_MAX:g}] m; "
        f"y = [{Y_MIN:g}, {Y_MAX:g}] m"
    )
    print(
        f"Intervals: {arguments.bins} along y, "
        f"{1000.0 * (edges[1] - edges[0]):g} mm each"
    )
    print(f"Particle volume: {particle_volume:.16g} m^3 ({volume_source})")
    if skipped:
        print(f"Skipped {skipped} incomplete or unsettled snapshot(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
