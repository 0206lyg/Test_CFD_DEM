#!/usr/bin/env python3
"""Regression for particle-particle tangential velocity projection."""

from __future__ import annotations

import math
import random
import re
from pathlib import Path
from typing import Tuple


Vector = Tuple[float, float, float]


def dot(a: Vector, b: Vector) -> float:
    return sum(a[index] * b[index] for index in range(3))


def subtract(a: Vector, b: Vector) -> Vector:
    return tuple(a[index] - b[index] for index in range(3))  # type: ignore[return-value]


def scale(value: float, vector: Vector) -> Vector:
    return tuple(value * component for component in vector)  # type: ignore[return-value]


def unit(vector: Vector) -> Vector:
    length = math.sqrt(dot(vector, vector))
    return scale(1.0 / length, vector)


def main() -> int:
    repo_dir = Path(__file__).resolve().parents[2]
    source_path = repo_dir / "src" / "HFDIBDEM" / "contactModels" / "prtSubContactInfo.C"
    source = source_path.read_text()
    compact = re.sub(r"\s+", "", source)

    correct_expression = (
        "prtCntVars_.contactNormal_*"
        "(relVeli&prtCntVars_.contactNormal_)"
    )
    if correct_expression not in compact:
        raise AssertionError("Correct particle normal projection is absent from source")
    if "(relVeli)*(relVeli&prtCntVars_.contactNormal_)" in compact:
        raise AssertionError("Defective relVeli-scaled projection remains in source")

    generator = random.Random(1207)
    old_nonorthogonal = 0
    for _ in range(20000):
        normal = unit(
            (
                2.0 * generator.random() - 1.0,
                2.0 * generator.random() - 1.0,
                2.0 * generator.random() - 1.0,
            )
        )
        relative_velocity = (
            4.0 * generator.random() - 2.0,
            4.0 * generator.random() - 2.0,
            4.0 * generator.random() - 2.0,
        )
        normal_projection = scale(dot(relative_velocity, normal), normal)
        tangential_velocity = subtract(relative_velocity, normal_projection)
        if abs(dot(tangential_velocity, normal)) > 2.0e-12:
            raise AssertionError("Correct projection produced non-tangential velocity")

        defective_projection = scale(dot(relative_velocity, normal), relative_velocity)
        defective_tangent = subtract(relative_velocity, defective_projection)
        if abs(dot(defective_tangent, normal)) > 1.0e-8:
            old_nonorthogonal += 1

    if old_nonorthogonal < 19000:
        raise AssertionError("Regression inputs do not expose the defective projection")

    print(
        "Particle-particle projection validated: 20000 random Vt vectors are "
        "orthogonal to the contact normal; the former expression fails the "
        "same invariant."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
