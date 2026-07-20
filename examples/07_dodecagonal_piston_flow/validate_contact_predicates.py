#!/usr/bin/env python3
"""Verify the optimized finite-wall AABB predicates against the old form.

This is a deterministic, OpenFOAM-independent regression check.  For every
generated DEM polygon (including all twelve contraction trapezoids), it
compares the former eight-corner tests with the center/support-radius form now
used in finiteWallGeometry.
"""

from __future__ import annotations

import math
import random
from typing import Iterable, Sequence, Tuple

from generate_dodecagon_geometry import (
    Point3,
    Wall,
    build_walls,
    cross,
    dot,
    magnitude,
    vector_scale,
    vector_sub,
)


Box = Tuple[Point3, Point3]


def add(a: Point3, b: Point3) -> Point3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def box_corners(box: Box) -> Iterable[Point3]:
    lo, hi = box
    for i in range(2):
        for j in range(2):
            for k in range(2):
                yield (
                    hi[0] if i else lo[0],
                    hi[1] if j else lo[1],
                    hi[2] if k else lo[2],
                )


def box_center_half_span(box: Box) -> Tuple[Point3, Point3]:
    lo, hi = box
    center = tuple(0.5 * (lo[i] + hi[i]) for i in range(3))
    half_span = tuple(0.5 * (hi[i] - lo[i]) for i in range(3))
    return center, half_span  # type: ignore[return-value]


def support_radius(direction: Point3, half_span: Point3) -> float:
    return sum(abs(direction[i]) * half_span[i] for i in range(3))


class PolygonPredicates:
    def __init__(self, wall: Wall):
        self.wall = wall
        self.normal = vector_scale(1.0 / magnitude(wall.normal), wall.normal)

        first_edge = vector_sub(wall.vertices[1], wall.vertices[0])
        first_edge = vector_sub(
            first_edge,
            vector_scale(dot(first_edge, self.normal), self.normal),
        )
        self.tangent1 = vector_scale(1.0 / magnitude(first_edge), first_edge)
        self.tangent2 = cross(self.normal, self.tangent1)
        self.tangent2 = vector_scale(1.0 / magnitude(self.tangent2), self.tangent2)

        local = [self.local_coordinates(vertex) for vertex in wall.vertices]
        twice_area = sum(
            local[i][0] * local[(i + 1) % len(local)][1]
            - local[(i + 1) % len(local)][0] * local[i][1]
            for i in range(len(local))
        )
        if twice_area < 0:
            local.reverse()
        self.local_vertices = local

        self.local_edge_normals = []
        self.global_edge_normals = []
        self.edge_offsets = []
        for i, a in enumerate(local):
            b = local[(i + 1) % len(local)]
            edge_u = b[0] - a[0]
            edge_v = b[1] - a[1]
            inverse_length = 1.0 / math.sqrt(edge_u * edge_u + edge_v * edge_v)
            local_normal = (-edge_v * inverse_length, edge_u * inverse_length, 0.0)
            global_normal = add(
                vector_scale(local_normal[0], self.tangent1),
                vector_scale(local_normal[1], self.tangent2),
            )
            self.local_edge_normals.append(local_normal)
            self.global_edge_normals.append(global_normal)
            self.edge_offsets.append(-dot(local_normal, a))

        spans = tuple(
            max(vertex[i] for vertex in wall.vertices)
            - min(vertex[i] for vertex in wall.vertices)
            for i in range(3)
        )
        patch_scale = max(magnitude(spans), 1.0e-3)
        self.tolerance = max(1.0e-12, 1.0e-10 * patch_scale)

    def signed_distance(self, point: Point3) -> float:
        return dot(vector_sub(point, self.wall.plane_point), self.normal)

    def local_coordinates(self, point: Point3) -> Point3:
        relative = vector_sub(point, self.wall.plane_point)
        return (dot(relative, self.tangent1), dot(relative, self.tangent2), 0.0)

    def old_classify(self, box: Box) -> str:
        corners = list(box_corners(box))
        distances = [self.signed_distance(point) for point in corners]
        if not any(value >= -self.tolerance for value in distances):
            return "outside"

        all_support = True
        local_corners = [self.local_coordinates(point) for point in corners]
        for i, a in enumerate(self.local_vertices):
            b = self.local_vertices[(i + 1) % len(self.local_vertices)]
            edge_u = b[0] - a[0]
            edge_v = b[1] - a[1]
            edge_length = math.sqrt(edge_u * edge_u + edge_v * edge_v)
            values = [
                (
                    edge_u * (point[1] - a[1])
                    - edge_v * (point[0] - a[0])
                )
                / edge_length
                for point in local_corners
            ]
            if not any(value >= -self.tolerance for value in values):
                return "outside"
            all_support = all_support and all(
                value >= -self.tolerance for value in values
            )

        if all(value >= -self.tolerance for value in distances) and all_support:
            return "inside"
        return "mixed"

    def new_classify(self, box: Box) -> str:
        center, half_span = box_center_half_span(box)
        center_distance = self.signed_distance(center)
        plane_radius = support_radius(self.normal, half_span)
        min_distance = center_distance - plane_radius
        max_distance = center_distance + plane_radius
        if max_distance < -self.tolerance:
            return "outside"

        all_support = True
        relative_center = vector_sub(center, self.wall.plane_point)
        for edge_normal, offset in zip(
            self.global_edge_normals, self.edge_offsets
        ):
            center_value = dot(edge_normal, relative_center) + offset
            radius = support_radius(edge_normal, half_span)
            if center_value + radius < -self.tolerance:
                return "outside"
            all_support = (
                all_support
                and center_value - radius >= -self.tolerance
            )

        if min_distance >= -self.tolerance and all_support:
            return "inside"
        return "mixed"

    def old_plane_intersects(self, box: Box) -> bool:
        values = [self.signed_distance(point) for point in box_corners(box)]
        return min(values) <= self.tolerance and max(values) >= -self.tolerance

    def new_plane_intersects(self, box: Box) -> bool:
        center, half_span = box_center_half_span(box)
        center_value = self.signed_distance(center)
        radius = support_radius(self.normal, half_span)
        return (
            center_value - radius <= self.tolerance
            and center_value + radius >= -self.tolerance
        )

    def old_plane_penetrates(self, box: Box) -> bool:
        return any(self.signed_distance(point) >= 0.0 for point in box_corners(box))

    def new_plane_penetrates(self, box: Box) -> bool:
        center, half_span = box_center_half_span(box)
        return (
            self.signed_distance(center)
            + support_radius(self.normal, half_span)
            >= 0.0
        )


def random_boxes(wall: Wall, seed: int, count: int) -> Iterable[Box]:
    rng = random.Random(seed)
    for _ in range(count):
        center = (
            rng.uniform(-0.065, 0.065),
            rng.uniform(-0.015, 0.165),
            rng.uniform(-0.065, 0.065),
        )
        half_span = tuple(10 ** rng.uniform(-6.0, -1.65) for _ in range(3))
        yield (
            tuple(center[i] - half_span[i] for i in range(3)),
            tuple(center[i] + half_span[i] for i in range(3)),
        )  # type: ignore[misc]

    # Concentrate additional boxes around this polygon's plane and vertices.
    for _ in range(count // 4):
        anchor = wall.vertices[rng.randrange(len(wall.vertices))]
        center = tuple(anchor[i] + rng.uniform(-2.0e-3, 2.0e-3) for i in range(3))
        half_span = tuple(10 ** rng.uniform(-6.0, -2.3) for _ in range(3))
        yield (
            tuple(center[i] - half_span[i] for i in range(3)),
            tuple(center[i] + half_span[i] for i in range(3)),
        )  # type: ignore[misc]


def main() -> int:
    walls = build_walls()
    boxes_per_wall = 4000
    checks = 0

    for wall_i, wall in enumerate(walls):
        predicates = PolygonPredicates(wall)
        for box in random_boxes(wall, 90210 + wall_i, boxes_per_wall):
            old_classification = predicates.old_classify(box)
            new_classification = predicates.new_classify(box)
            if old_classification != new_classification:
                raise AssertionError(
                    f"classification mismatch for {wall.name}: "
                    f"old={old_classification}, new={new_classification}, box={box}"
                )
            if predicates.old_plane_intersects(box) != predicates.new_plane_intersects(box):
                raise AssertionError(f"plane-intersection mismatch for {wall.name}: {box}")
            if predicates.old_plane_penetrates(box) != predicates.new_plane_penetrates(box):
                raise AssertionError(f"plane-penetration mismatch for {wall.name}: {box}")
            checks += 3

    print(
        f"Finite-wall predicates equivalent: {checks} comparisons across "
        f"{len(walls)} polygons, including 12 shoulder trapezoids."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
