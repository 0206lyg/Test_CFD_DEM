#!/usr/bin/env python3
"""Validate the moving dodecagonal piston and polygon support followers."""

from __future__ import annotations

import math
from pathlib import Path

from generate_dodecagon_geometry import (
    LARGE_APOTHEM,
    N_SIDES,
    Y_CONTRACTION,
    Y_TOP,
    Wall,
    build_walls,
    dot,
    polygon_area_3d,
    vector_sub,
)


TOLERANCE = 2.0e-12


def moved_large_wall(wall: Wall, piston_y: float) -> Wall:
    if piston_y <= Y_CONTRACTION + TOLERANCE:
        raise ValueError("piston crossed the fixed large-wall cap")

    moved_vertices = tuple(
        (x, piston_y if abs(y - Y_TOP) <= TOLERANCE else y, z)
        for x, y, z in wall.vertices
    )
    return Wall(
        wall.name,
        wall.normal,
        wall.plane_point,
        moved_vertices,
        wall.kind,
    )


def validate_two_cap_extrusion(wall: Wall) -> None:
    y_values = [vertex[1] for vertex in wall.vertices]
    min_count = sum(abs(y - Y_CONTRACTION) <= TOLERANCE for y in y_values)
    max_count = sum(abs(y - Y_TOP) <= TOLERANCE for y in y_values)
    if min_count < 2 or max_count < 2 or min_count + max_count != len(y_values):
        raise AssertionError(f"{wall.name} is not a safe two-cap extrusion")
    if abs(wall.normal[1]) > TOLERANCE:
        raise AssertionError(f"{wall.name} normal is not perpendicular to piston motion")


def main() -> int:
    repo_dir = Path(__file__).resolve().parents[2]
    support_source = (
        repo_dir / "src" / "HFDIBDEM" / "openHFDIBDEM.C"
    ).read_text()
    if "supportLimit is not defined for polygonal finite wall" in support_source:
        raise AssertionError("Legacy polygon supportLimit rejection remains in source")
    for safety_message in (
        "must be a two-cap extrusion",
        "needs at least two vertices on each support cap",
        "crossed the fixed cap of follower",
    ):
        if safety_message not in support_source:
            raise AssertionError(f"Missing source safety check: {safety_message}")

    walls = build_walls()
    large_walls = [wall for wall in walls if wall.kind == "largeSide"]
    piston = next(wall for wall in walls if wall.kind == "piston")

    if len(large_walls) != N_SIDES:
        raise AssertionError("Expected 12 piston-limited large walls")
    for wall in large_walls:
        validate_two_cap_extrusion(wall)

    original_xz = {
        wall.name: tuple((vertex[0], vertex[2]) for vertex in wall.vertices)
        for wall in large_walls
    }
    perimeter = 2.0 * N_SIDES * LARGE_APOTHEM * math.tan(math.pi / N_SIDES)

    for piston_y in (Y_TOP, 0.12, 0.09, 0.06):
        moved_piston = tuple((x, piston_y, z) for x, _, z in piston.vertices)
        if any(abs(vertex[1] - piston_y) > TOLERANCE for vertex in moved_piston):
            raise AssertionError("Piston vertices did not translate together")

        moved_walls = [moved_large_wall(wall, piston_y) for wall in large_walls]
        total_side_area = 0.0
        for wall in moved_walls:
            if tuple((vertex[0], vertex[2]) for vertex in wall.vertices) != original_xz[wall.name]:
                raise AssertionError(f"{wall.name} changed tangential geometry")
            if sum(abs(vertex[1] - piston_y) <= TOLERANCE for vertex in wall.vertices) != 2:
                raise AssertionError(f"{wall.name} moving cap did not follow piston")
            if sum(abs(vertex[1] - Y_CONTRACTION) <= TOLERANCE for vertex in wall.vertices) != 2:
                raise AssertionError(f"{wall.name} fixed cap moved")
            for vertex in wall.vertices:
                if abs(dot(vector_sub(vertex, wall.plane_point), wall.normal)) > TOLERANCE:
                    raise AssertionError(f"{wall.name} became non-coplanar")
            total_side_area += polygon_area_3d(wall.vertices, wall.normal)

        expected_area = perimeter * (piston_y - Y_CONTRACTION)
        if abs(total_side_area - expected_area) > 2.0e-13:
            raise AssertionError(
                f"Large-wall area mismatch at y={piston_y}: "
                f"{total_side_area} vs {expected_area}"
            )

    try:
        moved_large_wall(large_walls[0], Y_CONTRACTION)
    except ValueError:
        pass
    else:
        raise AssertionError("Support crossing safety check did not reject zero span")

    print(
        "Moving-wall geometry validated: 12 safe two-cap polygon followers "
        "track the dodecagonal piston from y=0.15 to y=0.06."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
