#!/usr/bin/env python3

import os
import sys
import numpy as np


CASE_DIR = os.path.dirname(os.path.abspath(__file__))
POST_DIR = os.path.join(CASE_DIR, "postProcessing")
EXIT_FILE = os.path.join(POST_DIR, "particleExit.dat")
PRESSURE_DIR = os.path.join(POST_DIR, "pressureDrop")
PRESSURE_FILE_NAME = "surfaceFieldValue.dat"


def read_threshold():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 analyzeAvalanches.py DELTA_T_C")

    try:
        threshold = float(sys.argv[1])
    except ValueError:
        raise SystemExit("DELTA_T_C must be greater than zero")

    if not np.isfinite(threshold) or threshold <= 0:
        raise SystemExit("DELTA_T_C must be greater than zero")

    return threshold


def read_exit_times():
    if not os.path.isfile(EXIT_FILE):
        return np.empty(0, dtype=float)

    times = []

    with open(EXIT_FILE, "r", encoding="utf-8") as exit_file:
        for line in exit_file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split()

            if len(columns) < 5:
                raise ValueError(f"Invalid particle-exit row: {line}")

            times.append(float(columns[0]))

    return np.sort(np.asarray(times, dtype=float))


def read_observation_end():
    times = []

    if os.path.isdir(PRESSURE_DIR):
        for root, _, files in os.walk(PRESSURE_DIR):
            if PRESSURE_FILE_NAME not in files:
                continue

            path = os.path.join(root, PRESSURE_FILE_NAME)

            with open(path, "r", encoding="utf-8") as pressure_file:
                for line in pressure_file:
                    line = line.strip()

                    if not line or line.startswith("#"):
                        continue

                    columns = line.split()

                    if len(columns) >= 2:
                        try:
                            time_value = float(columns[0])
                            pressure_value = float(columns[1])
                        except ValueError:
                            continue

                        if np.isfinite(time_value) and np.isfinite(pressure_value):
                            times.append(time_value)

    if not times:
        raise ValueError(
            "No completed pressure records were found under "
            f"{PRESSURE_DIR}"
        )

    return max(times)


def make_avalanche_table(exit_times, threshold, observation_end):
    if not exit_times.size:
        return np.empty((0, 7), dtype=float)

    gaps = np.diff(exit_times)
    tolerance = max(1e-12, threshold * 1e-12)
    starts = np.r_[0, np.flatnonzero(gaps > threshold + tolerance) + 1]
    ends = np.r_[starts[1:] - 1, exit_times.size - 1]
    rows = []

    for index, (first, last) in enumerate(zip(starts, ends)):
        start_time = exit_times[first]
        end_time = exit_times[last]

        if index == 0:
            preceding_arrest = np.nan
        else:
            preceding_arrest = start_time - exit_times[ends[index - 1]]

        if index + 1 < starts.size:
            following_arrest = exit_times[starts[index + 1]] - end_time
            censored = 0
        else:
            following_arrest = max(0.0, observation_end - end_time)
            censored = 1

        rows.append(
            (
                index + 1,
                start_time,
                end_time,
                last - first + 1,
                preceding_arrest,
                following_arrest,
                censored,
            )
        )

    return np.asarray(rows, dtype=float)


def main():
    threshold = read_threshold()
    observation_end = read_observation_end()
    exit_times = read_exit_times()
    exit_times = exit_times[exit_times <= observation_end]

    inter_exit = np.column_stack((exit_times[1:], np.diff(exit_times)))
    avalanches = make_avalanche_table(
        exit_times,
        threshold,
        observation_end,
    )

    np.savetxt(
        os.path.join(POST_DIR, "interExitTimes.txt"),
        inter_exit,
        fmt="%.12g",
        header="exit_time inter_exit_time",
    )
    np.savetxt(
        os.path.join(POST_DIR, "avalanches.txt"),
        avalanches,
        fmt=("%d", "%.12g", "%.12g", "%d", "%.12g", "%.12g", "%d"),
        header=(
            "avalanche_id start_time end_time size "
            "preceding_arrest_time following_arrest_time "
            "following_arrest_censored\n"
            f"delta_t_c {threshold:.12g}\n"
            f"observation_end {observation_end:.12g}"
        ),
    )

    print(f"Observation end : {observation_end:.12g}")
    print(f"Completed exits : {exit_times.size}")
    print(f"Avalanches      : {avalanches.shape[0]}")


if __name__ == "__main__":
    main()
