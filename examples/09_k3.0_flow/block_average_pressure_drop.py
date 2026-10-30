#!/usr/bin/env python3
"""Create 10-point non-overlapping block averages of pressure-drop data.

Default input directory:
    CASE/postProcessing/pressureDrop

Default output file:
    CASE/postProcessing/pressureDrop_blockAverage.dat

Only the Python standard library is required. Time spacing may be nonuniform.
Each output row contains:
    column 1: arithmetic mean of the times in the block
    column 2: arithmetic mean of the pressure values in the block

All numeric time subdirectories are read in numerical order. A directory named
T supplies only rows with T <= time < the next directory's start time. The last
directory supplies all its rows at or after T. For example, if directories 0
and 15 exist, rows at or after time 15 in directory 0 are discarded, even if
the first sample in directory 15 is later than 15. No old samples are used to
fill gaps in a restarted run.

After merging, consecutive rows are grouped in blocks of 10, including across
restart boundaries. A final incomplete block is retained and averaged over
its available rows. At least 10 merged rows are required; individual restart
files may contain fewer. --input FILE retains the single-file workflow.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Sequence


BLOCK_SIZE = 10
DEFAULT_INPUT_DIRECTORY = Path("postProcessing/pressureDrop")
DEFAULT_OUTPUT_FILE = Path("postProcessing/pressureDrop_blockAverage.dat")
PREFERRED_INPUT_NAMES = ("surfaceFieldValue.dat", "pressureDrop.dat")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Average consecutive non-overlapping groups of 10 pressure-drop "
            "samples after merging all numeric restart directories in "
            "postProcessing/pressureDrop. Each restart replaces earlier-run "
            "data at and after its directory time."
        )
    )
    parser.add_argument(
        "case",
        nargs="?",
        default=".",
        help="OpenFOAM case directory (default: current directory)",
    )
    parser.add_argument(
        "--input",
        dest="input_file",
        help=(
            "read only this data file (without restart merging), relative to "
            "CASE unless absolute; by default all numeric restart directories "
            "in postProcessing/pressureDrop are merged"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=str(DEFAULT_OUTPUT_FILE),
        help=(
            "output file, relative to CASE unless absolute "
            f"(default: {DEFAULT_OUTPUT_FILE})"
        ),
    )
    return parser.parse_args(argv)


def resolve_from_case(case_directory: Path, requested: str | Path) -> Path:
    path = Path(requested).expanduser()
    return path if path.is_absolute() else case_directory / path


def detect_directory_file(input_directory: Path) -> Path:
    for name in PREFERRED_INPUT_NAMES:
        candidate = input_directory / name
        if candidate.is_file():
            return candidate

    candidates = sorted(
        path
        for path in input_directory.iterdir()
        if path.is_file() and not path.name.startswith(".")
    )
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise SystemExit(f"ERROR: no data file found in {input_directory}")

    listing = "\n  ".join(str(path) for path in candidates)
    raise SystemExit(
        "ERROR: more than one possible input file was found; specify the exact "
        f"file with --input:\n  {listing}"
    )


def detect_input_sources(
    case_directory: Path, requested: str | None,
) -> list[tuple[float | None, Path]]:
    """Return (restart time, data file); None denotes an explicit single file."""
    if requested is not None:
        path = resolve_from_case(case_directory, requested)
        if not path.is_file():
            raise SystemExit(f"ERROR: input file not found: {path}")
        return [(None, path)]

    input_directory = case_directory / DEFAULT_INPUT_DIRECTORY
    if not input_directory.is_dir():
        raise SystemExit(f"ERROR: input directory not found: {input_directory}")

    directories: list[tuple[float, Path]] = []
    for path in input_directory.iterdir():
        if not path.is_dir():
            continue
        try:
            start_time = float(path.name)
        except ValueError:
            continue
        if math.isfinite(start_time):
            directories.append((start_time, path))

    directories.sort(key=lambda item: (item[0], item[1].name))
    if not directories:
        raise SystemExit(
            f"ERROR: no numeric time directories found in {input_directory}"
        )

    for previous, current in zip(directories, directories[1:]):
        if previous[0] == current[0]:
            raise SystemExit(
                "ERROR: multiple directories represent the same restart time: "
                f"{previous[1]} and {current[1]}"
            )

    # Missing or ambiguous restart files are errors: skipping one could retain
    # invalid data from an earlier run beyond that restart point.
    return [
        (start_time, detect_directory_file(path))
        for start_time, path in directories
    ]


def read_two_column_data(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise SystemExit(f"ERROR: cannot read {path}: {exc}") from exc

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 2:
            raise SystemExit(
                f"ERROR: expected at least two columns in {path}:{line_number}"
            )

        try:
            time_value = float(fields[0])
            pressure_value = float(fields[1])
        except ValueError as exc:
            raise SystemExit(
                f"ERROR: non-numeric value in the first two columns of "
                f"{path}:{line_number}: {raw_line!r}"
            ) from exc

        if not (math.isfinite(time_value) and math.isfinite(pressure_value)):
            raise SystemExit(
                f"ERROR: non-finite value in {path}:{line_number}: {raw_line!r}"
            )
        rows.append((time_value, pressure_value))

    if not rows:
        raise SystemExit(f"ERROR: no numeric data found in {path}")

    # The intervals may vary. Only the chronological order is checked.
    for index in range(1, len(rows)):
        previous_time = rows[index - 1][0]
        current_time = rows[index][0]
        if current_time <= previous_time:
            raise SystemExit(
                f"ERROR: time values must be strictly increasing in {path}; found "
                f"{previous_time:.16g} followed by {current_time:.16g}"
            )

    return rows


def read_merged_data(
    sources: Sequence[tuple[float | None, Path]],
) -> tuple[list[tuple[float, float]], list[tuple[int, int]]]:
    """Merge disjoint restart intervals before block averaging.

    Counts are (rows read, rows retained) for each source in the same order.
    """
    merged_rows: list[tuple[float, float]] = []
    counts: list[tuple[int, int]] = []

    for index, (start_time, path) in enumerate(sources):
        rows = read_two_column_data(path)
        if start_time is None:
            retained = rows
        else:
            end_time = (
                sources[index + 1][0] if index + 1 < len(sources) else math.inf
            )
            assert end_time is not None
            retained = [row for row in rows if start_time <= row[0] < end_time]
        merged_rows.extend(retained)
        counts.append((len(rows), len(retained)))

    if len(merged_rows) < BLOCK_SIZE:
        raise SystemExit(
            f"ERROR: only {len(merged_rows)} numeric rows remain after input "
            f"selection; at least {BLOCK_SIZE} are required"
        )
    return merged_rows, counts


def calculate_block_averages(
    rows: Sequence[tuple[float, float]],
) -> tuple[list[tuple[float, float]], int]:
    output_rows: list[tuple[float, float]] = []

    for start in range(0, len(rows), BLOCK_SIZE):
        block = rows[start : start + BLOCK_SIZE]
        sample_count = len(block)
        mean_time = math.fsum(row[0] for row in block) / sample_count
        mean_pressure = math.fsum(row[1] for row in block) / sample_count
        output_rows.append((mean_time, mean_pressure))

    final_block_size = len(rows) % BLOCK_SIZE or BLOCK_SIZE
    return output_rows, final_block_size


def write_output(
    path: Path,
    sources: Sequence[tuple[float | None, Path]],
    rows: Sequence[tuple[float, float]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")

    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as output:
            output.write("# 10-point non-overlapping block average\n")
            if sources[0][0] is not None:
                output.write(
                    "# restart_merge: each directory supplies "
                    "[start_time, next_start_time); latest has no upper bound\n"
                )
            for start_time, source in sources:
                output.write(f"# source {source}\n")
                if start_time is not None:
                    output.write(f"# restart_start={start_time:.16g}\n")
            output.write(
                f"# block_size={BLOCK_SIZE} time_spacing=arbitrary "
                "final_incomplete_block=retained\n"
            )
            output.write("# Time blockAverage\n")
            for time_value, pressure_average in rows:
                output.write(f"{time_value:.16g} {pressure_average:.16g}\n")
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
    if not case_directory.is_dir():
        print(f"ERROR: case directory not found: {case_directory}", file=sys.stderr)
        return 2

    try:
        sources = detect_input_sources(case_directory, arguments.input_file)
        raw_rows, counts = read_merged_data(sources)
        averaged_rows, final_block_size = calculate_block_averages(raw_rows)
        output_file = resolve_from_case(case_directory, arguments.output)
        write_output(output_file, sources, averaged_rows)
    except SystemExit as exc:
        message = str(exc)
        if message:
            print(message, file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"ERROR: file operation failed: {exc}", file=sys.stderr)
        return 2

    for (_, source), (read_count, retained_count) in zip(sources, counts):
        print(f"Used {retained_count} of {read_count} rows from {source}")
    print(f"Selected {len(raw_rows)} rows from {len(sources)} input file(s)")
    print(f"Wrote {len(averaged_rows)} block-average rows to {output_file}")
    if final_block_size < BLOCK_SIZE:
        print(
            f"Retained the final incomplete block containing "
            f"{final_block_size} row(s)."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
