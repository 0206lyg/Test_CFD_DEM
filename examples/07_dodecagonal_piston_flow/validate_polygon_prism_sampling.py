#!/usr/bin/env python3
"""Deterministic regression for the polygon-prism center sampler.

This mirrors the small geometry kernel in polygonPrismGeometry.C.  It checks
area-weighted uniformity and, more importantly, verifies that every vertex of
the randomly rotated/scaled case STL remains inside the complete prism.
"""

from __future__ import annotations

import math
import random
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from generate_dodecagon_geometry import (
    LARGE_APOTHEM,
    N_SIDES,
    Y_CONTRACTION,
    Y_TOP,
    regular_polygon_vertices,
)


Point2 = Tuple[float, float]
Point3 = Tuple[float, float, float]
EPS = 2.0e-12


def signed_area(points: Sequence[Point2]) -> float:
    return 0.5 * sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    )


def normalized_polygon(points: Sequence[Point2]) -> List[Point2]:
    polygon = list(points)
    if signed_area(polygon) < 0.0:
        polygon.reverse()
    return polygon


def half_planes(polygon: Sequence[Point2]) -> List[Tuple[float, float, float]]:
    """Return unit inward half-planes nx*x + nz*z + offset >= 0."""

    planes: List[Tuple[float, float, float]] = []
    for index, a in enumerate(polygon):
        b = polygon[(index + 1) % len(polygon)]
        edge_x = b[0] - a[0]
        edge_z = b[1] - a[1]
        length = math.hypot(edge_x, edge_z)
        nx, nz = -edge_z / length, edge_x / length
        planes.append((nx, nz, -(nx * a[0] + nz * a[1])))
    return planes


def clip_half_plane(
    polygon: Sequence[Point2], nx: float, nz: float, offset: float
) -> List[Point2]:
    if not polygon:
        return []

    result: List[Point2] = []
    previous = polygon[-1]
    previous_value = nx * previous[0] + nz * previous[1] + offset
    previous_inside = previous_value >= -EPS

    for current in polygon:
        current_value = nx * current[0] + nz * current[1] + offset
        current_inside = current_value >= -EPS

        if current_inside != previous_inside:
            denominator = previous_value - current_value
            fraction = previous_value / denominator
            result.append(
                (
                    previous[0] + fraction * (current[0] - previous[0]),
                    previous[1] + fraction * (current[1] - previous[1]),
                )
            )
        if current_inside:
            result.append(current)

        previous = current
        previous_value = current_value
        previous_inside = current_inside

    return result


def feasible_center_polygon(
    domain: Sequence[Point2],
    relative_body_points: Sequence[Point3],
    isotropic_radius: float = 0.0,
) -> List[Point2]:
    center_polygon = normalized_polygon(domain)
    for nx, nz, offset in half_planes(center_polygon):
        minimum_support = min(
            (nx * point[0] + nz * point[2] for point in relative_body_points),
            default=0.0,
        )
        center_polygon = clip_half_plane(
            center_polygon,
            nx,
            nz,
            offset + minimum_support - isotropic_radius,
        )
    return center_polygon


def triangle_areas(polygon: Sequence[Point2]) -> List[float]:
    a = polygon[0]
    return [
        0.5
        * abs(
            (polygon[index + 1][0] - a[0])
            * (polygon[index + 2][1] - a[1])
            - (polygon[index + 1][1] - a[1])
            * (polygon[index + 2][0] - a[0])
        )
        for index in range(len(polygon) - 2)
    ]


def sample_polygon(
    polygon: Sequence[Point2], selector: float, barycentric_1: float, barycentric_2: float
) -> Tuple[Point2, int]:
    areas = triangle_areas(polygon)
    target = min(max(selector, 0.0), 1.0) * sum(areas)
    cumulative = 0.0
    selected = -1
    for index, area in enumerate(areas):
        if area <= EPS * EPS:
            continue
        cumulative += area
        selected = index
        if target < cumulative:
            break
    if selected < 0:
        raise AssertionError("No positive-area fan triangle")

    root = math.sqrt(min(max(barycentric_1, 0.0), 1.0))
    edge = min(max(barycentric_2, 0.0), 1.0)
    a = polygon[0]
    b = polygon[selected + 1]
    c = polygon[selected + 2]
    return (
        (
            (1.0 - root) * a[0]
            + root * (1.0 - edge) * b[0]
            + root * edge * c[0],
            (1.0 - root) * a[1]
            + root * (1.0 - edge) * b[1]
            + root * edge * c[1],
        ),
        selected,
    )


def axial_center_range(
    relative_body_points: Sequence[Point3], isotropic_radius: float = 0.0
) -> Tuple[float, float]:
    minimum_support = min((point[1] for point in relative_body_points), default=0.0)
    maximum_support = max((point[1] for point in relative_body_points), default=0.0)
    return (
        Y_CONTRACTION - minimum_support + isotropic_radius,
        Y_TOP - maximum_support - isotropic_radius,
    )


def sample_prism_center(
    domain: Sequence[Point2],
    relative_body_points: Sequence[Point3],
    values: Sequence[float],
    isotropic_radius: float = 0.0,
) -> Tuple[Point3, int]:
    polygon = feasible_center_polygon(domain, relative_body_points, isotropic_radius)
    if len(polygon) < 3:
        raise AssertionError("Body has no feasible cross-sectional center region")
    axial_min, axial_max = axial_center_range(relative_body_points, isotropic_radius)
    if axial_max < axial_min:
        raise AssertionError("Body has no feasible axial center range")
    (x, z), triangle = sample_polygon(polygon, values[0], values[1], values[2])
    y = axial_min + values[3] * (axial_max - axial_min)
    return (x, y, z), triangle


def assert_body_inside(
    domain: Sequence[Point2], center: Point3, relative_body_points: Sequence[Point3]
) -> None:
    planes = half_planes(normalized_polygon(domain))
    for point in relative_body_points:
        translated = (
            center[0] + point[0],
            center[1] + point[1],
            center[2] + point[2],
        )
        if translated[1] < Y_CONTRACTION - EPS or translated[1] > Y_TOP + EPS:
            raise AssertionError(f"Body point escaped axial interval: {translated}")
        for nx, nz, offset in planes:
            if nx * translated[0] + nz * translated[2] + offset < -EPS:
                raise AssertionError(f"Body point escaped polygon: {translated}")


def read_unique_stl_vertices(path: Path) -> List[Point3]:
    pattern = re.compile(
        r"^\s*vertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$"
    )
    vertices: List[Point3] = []
    seen = set()
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            point = tuple(float(value) for value in match.groups())
            if point not in seen:
                seen.add(point)
                vertices.append(point)  # type: ignore[arg-type]
    if not vertices:
        raise AssertionError(f"No ASCII STL vertices found in {path}")
    return vertices


def rotate(point: Point3, angle: float, axis: Point3) -> Point3:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    dot_product = sum(point[index] * axis[index] for index in range(3))
    cross_product = (
        axis[1] * point[2] - axis[2] * point[1],
        axis[2] * point[0] - axis[0] * point[2],
        axis[0] * point[1] - axis[1] * point[0],
    )
    return tuple(
        cosine * point[index]
        + sine * cross_product[index]
        + (1.0 - cosine) * dot_product * axis[index]
        for index in range(3)
    )  # type: ignore[return-value]


def centered(points: Sequence[Point3]) -> List[Point3]:
    center = tuple(sum(point[index] for point in points) / len(points) for index in range(3))
    return [
        tuple(point[index] - center[index] for index in range(3))  # type: ignore[misc]
        for point in points
    ]


def random_unit_axis(generator: random.Random) -> Point3:
    # Mirrors the existing onceScatter positive-octant random-axis behavior.
    axis = (generator.random(), generator.random(), generator.random())
    length = math.sqrt(sum(value * value for value in axis))
    return tuple(value / length for value in axis)  # type: ignore[return-value]


def chi_square(observed: Sequence[int], expected: Sequence[float]) -> float:
    return sum((actual - target) ** 2 / target for actual, target in zip(observed, expected))


def validate_uniform_point_sampling(domain: Sequence[Point2]) -> None:
    generator = random.Random(1206)
    polygon = feasible_center_polygon(domain, [])
    areas = triangle_areas(polygon)
    total_area = sum(areas)
    sample_count = 120_000
    selected_counts = [0] * len(areas)
    sector_counts = [0] * N_SIDES
    sum_x = sum_y = sum_z = sum_xz = 0.0

    for _ in range(sample_count):
        values = [generator.random() for _ in range(4)]
        center, selected = sample_prism_center(domain, [], values)
        selected_counts[selected] += 1
        angle = (math.atan2(center[2], center[0]) + 2.0 * math.pi) % (2.0 * math.pi)
        sector_counts[int(angle / (2.0 * math.pi / N_SIDES)) % N_SIDES] += 1
        sum_x += center[0]
        sum_y += center[1]
        sum_z += center[2]
        sum_xz += center[0] * center[2]

    expected_triangles = [sample_count * area / total_area for area in areas]
    triangle_statistic = chi_square(selected_counts, expected_triangles)
    sector_statistic = chi_square(sector_counts, [sample_count / N_SIDES] * N_SIDES)
    means = (sum_x / sample_count, sum_y / sample_count, sum_z / sample_count)

    if triangle_statistic > 35.0:
        raise AssertionError(f"Fan-triangle weighting failed: chi2={triangle_statistic}")
    if sector_statistic > 40.0:
        raise AssertionError(f"Rotational sector symmetry failed: chi2={sector_statistic}")
    if abs(means[0]) > 3.0e-4 or abs(means[2]) > 3.0e-4:
        raise AssertionError(f"Cross-sectional sample mean is biased: {means}")
    if abs(means[1] - 0.5 * (Y_CONTRACTION + Y_TOP)) > 3.0e-4:
        raise AssertionError(f"Axial sample mean is biased: {means[1]}")
    if abs(sum_xz / sample_count) > 1.5e-5:
        raise AssertionError("Cross-sectional sample covariance is unexpectedly biased")


def validate_rotated_case_stl(domain: Sequence[Point2], case_dir: Path) -> None:
    source_points = centered(
        read_unique_stl_vertices(case_dir / "constant" / "triSurface" / "PH_SC_reservoir.stl")
    )
    generator = random.Random(6021)

    for _ in range(600):
        angle = (2.0 * generator.random() - 1.0) * math.pi
        axis = random_unit_axis(generator)
        relative_points = [
            tuple(1.02 * value for value in rotate(point, angle, axis))  # type: ignore[misc]
            for point in source_points
        ]
        values = [generator.random() for _ in range(4)]
        center, _ = sample_prism_center(domain, relative_points, values)
        assert_body_inside(domain, center, relative_points)


def validate_generic_convex_domains() -> None:
    domains: Iterable[Sequence[Point2]] = (
        ((-0.03, -0.02), (0.04, -0.015), (0.005, 0.045)),
        ((-0.04, -0.03), (0.035, -0.025), (0.025, 0.03), (-0.02, 0.04)),
        tuple(
            (radius * math.cos(angle), radius * math.sin(angle))
            for radius, angle in (
                (0.035, 0.0),
                (0.041, 0.9),
                (0.037, 1.9),
                (0.043, 3.0),
                (0.038, 4.15),
                (0.040, 5.2),
            )
        ),
    )
    body = [(-0.0015, -0.002, -0.001), (0.001, 0.0025, -0.0015), (0.0, 0.0, 0.002)]
    generator = random.Random(206)
    for domain in domains:
        # Reverse one valid input to exercise winding independence.
        for ordered_domain in (domain, tuple(reversed(domain))):
            for _ in range(100):
                center, _ = sample_prism_center(
                    ordered_domain, body, [generator.random() for _ in range(4)]
                )
                assert_body_inside(ordered_domain, center, body)


def validate_sphere_support(domain: Sequence[Point2]) -> None:
    generator = random.Random(12)
    radius = 0.002
    planes = half_planes(normalized_polygon(domain))
    for _ in range(250):
        center, _ = sample_prism_center(
            domain, [(0.0, 0.0, 0.0)], [generator.random() for _ in range(4)], radius
        )
        if center[1] - radius < Y_CONTRACTION - EPS or center[1] + radius > Y_TOP + EPS:
            raise AssertionError("Analytic sphere escaped axial interval")
        if any(nx * center[0] + nz * center[2] + offset < radius - EPS for nx, nz, offset in planes):
            raise AssertionError("Analytic sphere escaped polygon")


def main() -> int:
    case_dir = Path(__file__).resolve().parent
    domain = regular_polygon_vertices(LARGE_APOTHEM)
    expected_area = N_SIDES * LARGE_APOTHEM**2 * math.tan(math.pi / N_SIDES)
    expected_volume = expected_area * (Y_TOP - Y_CONTRACTION)
    if abs(abs(signed_area(domain)) - expected_area) > 1.0e-15:
        raise AssertionError("Dodecagon area mismatch")
    if abs(expected_volume - 5.144624494677556e-4) > 1.0e-16:
        raise AssertionError("Dodecagon prism volume mismatch")

    validate_uniform_point_sampling(domain)
    validate_rotated_case_stl(domain, case_dir)
    validate_generic_convex_domains()
    validate_sphere_support(domain)

    print(
        "polygonPrism sampler validated: area-weighted uniform centers, "
        "600 rotated/scaled case-STL placements, sphere support, and generic "
        "convex polygons."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
