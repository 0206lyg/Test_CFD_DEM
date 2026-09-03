#!/usr/bin/env python3

import os
import numpy as np


CASE_DIR = os.path.dirname(os.path.abspath(__file__))
CONTROL_DICT = os.path.join(CASE_DIR, "system", "controlDict")
EXIT_FILE = os.path.join(CASE_DIR, "postProcessing", "particleExit.dat")
OUTPUT_FILE = os.path.join(
    CASE_DIR,
    "postProcessing",
    "particleExit_cumulative.txt",
)
INTERVAL_OUTPUT_FILE = os.path.join(
    CASE_DIR,
    "postProcessing",
    "particleExit_intervals_ccdf.dat",
)


def remove_comments(text):
    while "/*" in text:
        start = text.find("/*")
        end = text.find("*/", start + 2)

        if end < 0:
            text = text[:start]
            break

        text = text[:start] + text[end + 2:]

    return "\n".join(line.split("//", 1)[0] for line in text.splitlines())


def read_control_value(entry):
    with open(CONTROL_DICT, "r", encoding="utf-8") as control_file:
        text = remove_comments(control_file.read())

    for line in text.splitlines():
        words = line.replace(";", " ").split()

        if len(words) >= 2 and words[0] == entry:
            return float(words[1])

    raise ValueError(f"Cannot find '{entry}' in {CONTROL_DICT}")


def read_exit_times():
    if not os.path.isfile(EXIT_FILE):
        return np.empty(0, dtype=float)

    exit_times = []

    with open(EXIT_FILE, "r", encoding="utf-8") as exit_file:
        for line_number, line in enumerate(exit_file, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split()

            if len(columns) < 5:
                raise ValueError(
                    f"Invalid row {line_number} in {EXIT_FILE}: {line}"
                )

            exit_times.append(float(columns[0]))

    return np.asarray(exit_times, dtype=float)


def calculate_interval_statistics(exit_times):
    if exit_times.size == 0:
        return np.nan, np.nan, np.nan, np.empty((0, 2), dtype=float)

    if not np.all(np.isfinite(exit_times)):
        raise ValueError("Particle exit times must all be finite")

    particle_exit_start_time = exit_times[0]
    particle_exit_end_time = exit_times[-1]

    if exit_times.size < 2:
        return (
            particle_exit_start_time,
            particle_exit_end_time,
            np.nan,
            np.empty((0, 2), dtype=float),
        )

    sorted_intervals = np.sort(np.diff(exit_times))
    number_of_intervals = sorted_intervals.size
    interval_sum = np.sum(sorted_intervals)

    if interval_sum > 0.0:
        ranks = np.arange(1, number_of_intervals + 1, dtype=float)
        gini_coefficient = np.sum(
            (2.0 * ranks - number_of_intervals - 1.0) * sorted_intervals
        ) / (number_of_intervals * interval_sum)
        gini_coefficient = float(np.clip(gini_coefficient, 0.0, 1.0))
    else:
        gini_coefficient = np.nan

    time_scale = max(1.0, float(np.max(np.abs(exit_times))))
    equality_tolerance = 32.0 * np.finfo(float).eps * time_scale
    group_starts = np.concatenate(
        (
            np.array([True]),
            np.diff(sorted_intervals) > equality_tolerance,
        )
    )
    group_indices = np.flatnonzero(group_starts)
    interval_counts = np.diff(
        np.append(group_indices, number_of_intervals)
    )
    unique_intervals = (
        np.add.reduceat(sorted_intervals, group_indices) / interval_counts
    )
    ccdf = np.cumsum(interval_counts[::-1])[::-1] / number_of_intervals
    ccdf_output = np.column_stack((unique_intervals, ccdf))

    return (
        particle_exit_start_time,
        particle_exit_end_time,
        gini_coefficient,
        ccdf_output,
    )


def main():
    write_interval = read_control_value("writeInterval")
    start_time = read_control_value("startTime")
    sampling_interval = write_interval / 10.0

    if sampling_interval <= 0.0:
        raise ValueError("Invalid writeInterval in controlDict")

    exit_times = np.sort(read_exit_times())
    tolerance = max(1e-12, sampling_interval * 1e-9)

    if exit_times.size:
        last_exit_time = exit_times[-1]
        number_of_steps = max(
            0,
            int(
                np.ceil(
                    (last_exit_time - start_time - tolerance)
                    / sampling_interval
                )
            ),
        )
        sample_times = start_time + sampling_interval * np.arange(
            number_of_steps + 1
        )
        cumulative_counts = np.searchsorted(
            exit_times,
            sample_times + tolerance,
            side="right",
        )
        output = np.column_stack((sample_times, cumulative_counts))
    else:
        last_exit_time = None
        cumulative_counts = np.empty(0, dtype=int)
        output = np.empty((0, 2), dtype=float)

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    np.savetxt(
        OUTPUT_FILE,
        output,
        fmt=("%.12g", "%d"),
        header=(
            "time cumulative_discharged_particles\n"
            f"sampling_interval {sampling_interval:.12g}"
        ),
    )

    (
        particle_exit_start_time,
        particle_exit_end_time,
        gini_coefficient,
        ccdf_output,
    ) = calculate_interval_statistics(exit_times)

    np.savetxt(
        INTERVAL_OUTPUT_FILE,
        ccdf_output,
        fmt=("%.12g", "%.12g"),
        header=(
            f"particle_exit_start_time {particle_exit_start_time:.12g}\n"
            f"particle_exit_end_time {particle_exit_end_time:.12g}\n"
            f"Gini_coefficient {gini_coefficient:.12g}"
        ),
    )

    print(f"Input       : {EXIT_FILE}")
    print(f"Time spacing: {sampling_interval:.12g}")
    print(f"Last exit   : {last_exit_time}")
    print(f"Total exits : {exit_times.size}")
    print(f"Output      : {OUTPUT_FILE}")
    print(f"CCDF output : {INTERVAL_OUTPUT_FILE}")


if __name__ == "__main__":
    main()
