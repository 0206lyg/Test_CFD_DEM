#!/usr/bin/env python3

import re
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, Dict, List, Optional, TextIO


INPUT_FILE = Path("log.pimpleLYJHFDIBFoam")
OUTPUT_FILE = Path("disappeared_body_report.txt")

HISTORY_LENGTH = 100


TIME_RE = re.compile(
    r"^\s*Time\s*=\s*(\S+)"
)

BODY_START_RE = re.compile(
    r"^\s*--\s*body\s*:\s*(\d+)\s+"
    r"current center of mass position:\s*(\(.*\))"
)

CELL_SIZE_RE = re.compile(
    r"^\s*Body characteristic cell size:\s*(\S+)"
)

MASS_RATIO_RE = re.compile(
    r"^\s*--\s*body\s+(\d+)\s+current M/M0:\s*(\S+)"
)

CELL_COUNT_RE = re.compile(
    r"^\s*--\s*body\s+(\d+)\s+current cellCount:\s*(\S+)"
)


@dataclass
class BodyRecord:
    time: str
    position: Optional[str] = None
    cell_size: Optional[str] = None
    mass_ratio: Optional[str] = None
    cell_count: Optional[str] = None

    def is_complete(self) -> bool:
        return (
            self.position is not None
            and self.cell_size is not None
            and self.mass_ratio is not None
            and self.cell_count is not None
        )


def format_value(value: Optional[str]) -> str:
    return value if value is not None else "[not found]"


def write_disappearance_report(
    output: TextIO,
    body_id: int,
    previous_time: str,
    current_time: str,
    history: Deque[BodyRecord],
) -> None:
    output.write("\n")
    output.write("=" * 90 + "\n")
    output.write(f"body {body_id} disappeared!!\n")
    output.write(f"last present time : {previous_time}\n")
    output.write(f"first missing time: {current_time}\n")
    output.write(
        f"last {min(len(history), HISTORY_LENGTH)} records "
        f"before disappearance, oldest to newest:\n"
    )
    output.write("-" * 90 + "\n")

    if not history:
        output.write("[no previous history found]\n")
        output.write("=" * 90 + "\n")
        return

    for index, record in enumerate(history, start=1):
        output.write(
            f"[{index:03d}] "
            f"Time={record.time}  "
            f"M/M0={format_value(record.mass_ratio)}  "
            f"cellCount={format_value(record.cell_count)}  "
            f"cellSize={format_value(record.cell_size)}\n"
        )
        output.write(
            f"      center={format_value(record.position)}\n"
        )

    output.write("=" * 90 + "\n")


def detect_disappeared_bodies() -> None:
    if not INPUT_FILE.is_file():
        raise SystemExit(
            f"Error: input file does not exist: {INPUT_FILE}"
        )

    previous_time: Optional[str] = None
    previous_body_ids: Optional[set[int]] = None

    current_time: Optional[str] = None
    current_records: Dict[int, BodyRecord] = {}
    current_body_id: Optional[int] = None

    histories: Dict[int, Deque[BodyRecord]] = {}

    reported_bodies: set[int] = set()
    disappearance_count = 0
    time_step_count = 0

    def process_completed_time_step(output: TextIO) -> None:
        nonlocal previous_time
        nonlocal previous_body_ids
        nonlocal disappearance_count
        nonlocal time_step_count

        if current_time is None:
            return

        current_body_ids = set(current_records.keys())

        # 현재 time step에서 사라진 body를 먼저 확인한다.
        # histories에는 아직 현재 step은 들어 있지 않으므로,
        # 사라지기 직전까지의 이력만 들어 있다.
        if (
            previous_time is not None
            and previous_body_ids is not None
        ):
            disappeared_ids = sorted(
                previous_body_ids - current_body_ids
            )

            for body_id in disappeared_ids:
                if body_id in reported_bodies:
                    continue

                write_disappearance_report(
                    output=output,
                    body_id=body_id,
                    previous_time=previous_time,
                    current_time=current_time,
                    history=histories.get(
                        body_id,
                        deque(maxlen=HISTORY_LENGTH),
                    ),
                )

                reported_bodies.add(body_id)
                disappearance_count += 1

        # 현재 time step에서 발견된 body 기록을 history에 추가한다.
        for body_id, record in current_records.items():
            if body_id not in histories:
                histories[body_id] = deque(
                    maxlen=HISTORY_LENGTH
                )

            histories[body_id].append(record)

        previous_time = current_time
        previous_body_ids = current_body_ids
        time_step_count += 1

    with INPUT_FILE.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as log_file, OUTPUT_FILE.open(
        "w",
        encoding="utf-8",
    ) as output:

        output.write(f"Input log file: {INPUT_FILE}\n")
        output.write(
            "Step boundary: OpenFOAM 'Time = ...' line\n"
        )
        output.write(
            "Detection: body present at previous time but absent "
            "at current time\n"
        )
        output.write(
            f"History length: previous {HISTORY_LENGTH} time steps\n"
        )

        for line in log_file:
            time_match = TIME_RE.match(line)

            if time_match:
                process_completed_time_step(output)

                current_time = time_match.group(1)
                current_records = {}
                current_body_id = None
                continue

            if current_time is None:
                continue

            body_start_match = BODY_START_RE.match(line)

            if body_start_match:
                body_id = int(body_start_match.group(1))
                position = body_start_match.group(2)

                current_body_id = body_id

                # 동일한 Time 안에서 같은 body가 여러 번 계산되면
                # 가장 마지막 geometry 계산 결과로 교체한다.
                current_records[body_id] = BodyRecord(
                    time=current_time,
                    position=position,
                )
                continue

            cell_size_match = CELL_SIZE_RE.match(line)

            if (
                cell_size_match
                and current_body_id is not None
                and current_body_id in current_records
            ):
                current_records[current_body_id].cell_size = (
                    cell_size_match.group(1)
                )
                continue

            mass_ratio_match = MASS_RATIO_RE.match(line)

            if mass_ratio_match:
                body_id = int(mass_ratio_match.group(1))

                if body_id in current_records:
                    current_records[body_id].mass_ratio = (
                        mass_ratio_match.group(2)
                    )

                continue

            cell_count_match = CELL_COUNT_RE.match(line)

            if cell_count_match:
                body_id = int(cell_count_match.group(1))

                if body_id in current_records:
                    current_records[body_id].cell_count = (
                        cell_count_match.group(2)
                    )

                continue

        # 파일의 마지막 time step 처리
        process_completed_time_step(output)

        output.write("\n")
        output.write("-" * 90 + "\n")
        output.write(
            f"Total processed time steps: {time_step_count}\n"
        )
        output.write(
            f"Total unique disappeared bodies: "
            f"{disappearance_count}\n"
        )

        if disappearance_count == 0:
            output.write(
                "No disappeared bodies were detected.\n"
            )

    print(f"Analysis completed: {INPUT_FILE}")
    print(f"Report saved to: {OUTPUT_FILE}")
    print(
        f"Detected unique disappeared bodies: "
        f"{disappearance_count}"
    )


if __name__ == "__main__":
    detect_disappeared_bodies()