#!/usr/bin/env python3
"""Create 10-point non-overlapping block averages of pressure-drop data.

Default input directory:
    CASE/postProcessing/pressureDrop/0

Default output file:
    CASE/postProcessing/pressureDrop_blockAverage.dat

Only the Python standard library is required. Time spacing may be nonuniform.
Each output row contains:
    column 1: arithmetic mean of the times in the block
    column 2: arithmetic mean of the pressure values in the block

Consecutive input rows are grouped in blocks of 10. A final incomplete block
is retained and averaged over its available rows.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Sequence


BLOCK_SIZE = 10
DEFAULT_INPUT_DIRECTORY = Path("postProcessing/pressureDrop/0")
DEFAULT_OUTPUT_FILE = Path("postProcessing/pressureDrop_blockAverage.dat")
PREFERRED_INPUT_NAMES = ("surfaceFieldValue.dat", "pressureDrop.dat")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Average consecutive non-overlapping groups of 10 pressure-drop "
            "samples from postProcessing/pressureDrop/0."
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
            "input data file, relative to CASE unless absolute; by default the "
            "script detects the data file in postProcessing/pressureDrop/0"
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


def detect_input_file(case_directory: Path, requested: str | None) -> Path:
    if requested is not None:
        path = resolve_from_case(case_directory, requested)
        if not path.is_file():
            raise SystemExit(f"ERROR: input file not found: {path}")
        return path

    input_directory = case_directory / DEFAULT_INPUT_DIRECTORY
    if not input_directory.is_dir():
        raise SystemExit(f"ERROR: input directory not found: {input_directory}")

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

    if len(rows) < BLOCK_SIZE:
        raise SystemExit(
            f"ERROR: {path} contains only {len(rows)} numeric rows; "
            f"at least {BLOCK_SIZE} are required"
        )

    # The intervals may vary. Only the chronological order is checked.
    for index in range(1, len(rows)):
        previous_time = rows[index - 1][0]
        current_time = rows[index][0]
        if current_time <= previous_time:
            raise SystemExit(
                "ERROR: time values must be strictly increasing; found "
                f"{previous_time:.16g} followed by {current_time:.16g}"
            )

    return rows


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
    source: Path,
    rows: Sequence[tuple[float, float]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")

    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as output:
            output.write("# 10-point non-overlapping block average\n")
            output.write(f"# source {source}\n")
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
        input_file = detect_input_file(case_directory, arguments.input_file)
        raw_rows = read_two_column_data(input_file)
        averaged_rows, final_block_size = calculate_block_averages(raw_rows)
        output_file = resolve_from_case(case_directory, arguments.output)
        write_output(output_file, input_file, averaged_rows)
    except SystemExit as exc:
        message = str(exc)
        if message:
            print(message, file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"ERROR: file operation failed: {exc}", file=sys.stderr)
        return 2

    print(f"Read {len(raw_rows)} rows from {input_file}")
    print(f"Wrote {len(averaged_rows)} block-average rows to {output_file}")
    if final_block_size < BLOCK_SIZE:
        print(
            f"Retained the final incomplete block containing "
            f"{final_block_size} row(s)."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
