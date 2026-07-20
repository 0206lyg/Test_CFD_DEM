#!/usr/bin/env python3
"""Generate and validate the moving-piston dodecagonal contraction geometry.

The CFD mesh, dynamic-mesh zones, and all DEM finite-wall polygons are produced
from the same geometry constants.  The y=0 face remains a CFD outlet but has
no DEM wall; a moving dodecagonal piston closes the initial y=0.15 face.

The cross-section is a logically Cartesian 80 x 80 grid mapped exactly onto a
regular dodecagon.  The central 20 x 20 grid is the downstream channel.  Every
2-D cell is emitted as one block with one cell in each cross-sectional block
direction; the blocks are extruded uniformly in y.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


Point2 = Tuple[float, float]
Point3 = Tuple[float, float, float]
Face = Tuple[int, int, int, int]


N_SIDES = 12
SMALL_APOTHEM = 0.01
LARGE_APOTHEM = 0.04
Y_BOTTOM = 0.00
Y_CONTRACTION = 0.05
Y_TOP = 0.15

NY_DOWNSTREAM = 50
NY_UPSTREAM = 100

# The breakpoints put dodecagon vertices exactly on logical grid lines.
# The interval counts retain 20 cells across the small channel and 80 across
# the large chamber, matching the resolution of example 05.
TAN_15 = math.tan(math.pi / N_SIDES)
SCALE_RATIO = SMALL_APOTHEM / LARGE_APOTHEM
BASE_LOGICAL_COORDS = (
    -1.0,
    -TAN_15,
    -SCALE_RATIO,
    -SCALE_RATIO * TAN_15,
    0.0,
    SCALE_RATIO * TAN_15,
    SCALE_RATIO,
    TAN_15,
    1.0,
)
BASE_INTERVAL_CELLS = (29, 1, 7, 3, 3, 7, 1, 29)


def snap(value: float, tolerance: float = 5.0e-15) -> float:
    return 0.0 if abs(value) < tolerance else value


def fmt(value: float) -> str:
    value = snap(value)
    return "0" if value == 0.0 else f"{value:.15g}"


def fmt_point(point: Sequence[float]) -> str:
    return "(" + " ".join(fmt(value) for value in point) + ")"


def vector_add(a: Point3, b: Point3) -> Point3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vector_sub(a: Point3, b: Point3) -> Point3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vector_scale(value: float, a: Point3) -> Point3:
    return (value * a[0], value * a[1], value * a[2])


def dot(a: Point3, b: Point3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Point3, b: Point3) -> Point3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def magnitude(a: Point3) -> float:
    return math.sqrt(dot(a, a))


def logical_coordinates() -> Tuple[List[float], Dict[float, int]]:
    coordinates: List[float] = [BASE_LOGICAL_COORDS[0]]

    for interval, cell_count in enumerate(BASE_INTERVAL_CELLS):
        lo = BASE_LOGICAL_COORDS[interval]
        hi = BASE_LOGICAL_COORDS[interval + 1]
        for cell_i in range(1, cell_count + 1):
            coordinates.append(lo + (hi - lo) * cell_i / cell_count)

    if len(coordinates) != 81:
        raise RuntimeError(f"Expected 81 logical points, got {len(coordinates)}")

    landmark_indices: Dict[float, int] = {}
    for landmark in BASE_LOGICAL_COORDS:
        index = min(range(len(coordinates)), key=lambda i: abs(coordinates[i] - landmark))
        if abs(coordinates[index] - landmark) > 1.0e-13:
            raise RuntimeError(f"Missing logical landmark {landmark}")
        coordinates[index] = landmark
        landmark_indices[landmark] = index

    return coordinates, landmark_indices


FACET_NORMALS_2D: Tuple[Point2, ...] = tuple(
    (snap(math.cos(k * math.pi / 6.0)), snap(math.sin(k * math.pi / 6.0)))
    for k in range(N_SIDES)
)


def map_square_to_dodecagon(u: float, v: float) -> Point2:
    """Homogeneous map whose square contours become similar dodecagons."""

    rho = max(abs(u), abs(v))
    if rho == 0.0:
        return (0.0, 0.0)

    max_projection = max(nx * u + nz * v for nx, nz in FACET_NORMALS_2D)
    scale = rho * LARGE_APOTHEM / max_projection
    return (snap(scale * u), snap(scale * v))


def signed_area_2d(points: Sequence[Point2]) -> float:
    return 0.5 * sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    )


def polygon_area_3d(points: Sequence[Point3], normal: Point3) -> float:
    area_vector = (0.0, 0.0, 0.0)
    for i, point in enumerate(points):
        area_vector = vector_add(area_vector, cross(point, points[(i + 1) % len(points)]))
    return abs(dot(vector_scale(0.5, area_vector), normal))


def face_area_3d(points: Sequence[Point3]) -> float:
    area_vector = (0.0, 0.0, 0.0)
    for i, point in enumerate(points):
        area_vector = vector_add(area_vector, cross(point, points[(i + 1) % len(points)]))
    return 0.5 * magnitude(area_vector)


@dataclass(frozen=True)
class Block:
    vertices: Tuple[int, int, int, int, int, int, int, int]
    cells: Tuple[int, int, int]
    zone: str

    def faces(self) -> Tuple[Face, Face, Face, Face, Face, Face]:
        v = self.vertices
        return (
            (v[0], v[3], v[2], v[1]),
            (v[4], v[5], v[6], v[7]),
            (v[0], v[1], v[5], v[4]),
            (v[1], v[2], v[6], v[5]),
            (v[2], v[3], v[7], v[6]),
            (v[3], v[0], v[4], v[7]),
        )


@dataclass(frozen=True)
class Wall:
    name: str
    normal: Point3
    plane_point: Point3
    vertices: Tuple[Point3, ...]
    kind: str


@dataclass
class MeshData:
    vertices: List[Point3]
    blocks: List[Block]
    patches: Dict[str, List[Face]]
    total_cells: int


def build_mesh() -> MeshData:
    coords, landmarks = logical_coordinates()
    small_lo = landmarks[-SCALE_RATIO]
    small_hi = landmarks[SCALE_RATIO]

    vertices: List[Point3] = []
    layer_ids: Dict[str, Dict[Tuple[int, int], int]] = {}

    def add_layer(name: str, y: float, index_range: Iterable[int]) -> None:
        ids: Dict[Tuple[int, int], int] = {}
        indices = tuple(index_range)
        for i in indices:
            for j in indices:
                x, z = map_square_to_dodecagon(coords[i], coords[j])
                ids[(i, j)] = len(vertices)
                vertices.append((x, y, z))
        layer_ids[name] = ids

    add_layer("bottom", Y_BOTTOM, range(small_lo, small_hi + 1))
    add_layer("contraction", Y_CONTRACTION, range(len(coords)))
    add_layer("top", Y_TOP, range(len(coords)))

    blocks: List[Block] = []
    patches: Dict[str, List[Face]] = {
        "bottomOutlet": [],
        "pistonWall": [],
        "contractionWall": [],
        "smallWall": [],
        "largeWall": [],
    }

    def add_block(
        lower: str,
        upper: str,
        i: int,
        j: int,
        ny: int,
        region: str,
    ) -> None:
        lo = layer_ids[lower]
        hi = layer_ids[upper]

        # Clockwise in x-z when viewed from +y, followed by +y extrusion.
        block = Block(
            (
                lo[(i, j)],
                lo[(i, j + 1)],
                lo[(i + 1, j + 1)],
                lo[(i + 1, j)],
                hi[(i, j)],
                hi[(i, j + 1)],
                hi[(i + 1, j + 1)],
                hi[(i + 1, j)],
            ),
            (1, 1, ny),
            "largeReservoir" if region == "upstream" else "",
        )
        blocks.append(block)
        lower_face, upper_face, left_face, top_face, right_face, bottom_face = block.faces()

        if region == "downstream":
            patches["bottomOutlet"].append(lower_face)
            if i == small_lo:
                patches["smallWall"].append(left_face)
            if i == small_hi - 1:
                patches["smallWall"].append(right_face)
            if j == small_lo:
                patches["smallWall"].append(bottom_face)
            if j == small_hi - 1:
                patches["smallWall"].append(top_face)
        else:
            patches["pistonWall"].append(upper_face)
            inside_opening = (
                small_lo <= i < small_hi and small_lo <= j < small_hi
            )
            if not inside_opening:
                patches["contractionWall"].append(lower_face)
            if i == 0:
                patches["largeWall"].append(left_face)
            if i == len(coords) - 2:
                patches["largeWall"].append(right_face)
            if j == 0:
                patches["largeWall"].append(bottom_face)
            if j == len(coords) - 2:
                patches["largeWall"].append(top_face)

    for i in range(small_lo, small_hi):
        for j in range(small_lo, small_hi):
            add_block("bottom", "contraction", i, j, NY_DOWNSTREAM, "downstream")

    for i in range(len(coords) - 1):
        for j in range(len(coords) - 1):
            add_block("contraction", "top", i, j, NY_UPSTREAM, "upstream")

    total_cells = sum(block.cells[0] * block.cells[1] * block.cells[2] for block in blocks)
    mesh = MeshData(vertices, blocks, patches, total_cells)
    validate_mesh(mesh)
    return mesh


def validate_mesh(mesh: MeshData) -> None:
    if len(mesh.vertices) != 13563:
        raise RuntimeError(f"Unexpected vertex count: {len(mesh.vertices)}")
    if len(mesh.blocks) != 6800:
        raise RuntimeError(f"Unexpected block count: {len(mesh.blocks)}")
    if mesh.total_cells != 660000:
        raise RuntimeError(f"Unexpected cell count: {mesh.total_cells}")

    zone_blocks = [block for block in mesh.blocks if block.zone]
    if len(zone_blocks) != 6400:
        raise RuntimeError(f"Expected 6400 largeReservoir blocks, got {len(zone_blocks)}")
    if any(block.zone != "largeReservoir" for block in zone_blocks):
        raise RuntimeError("Unexpected dynamic-mesh cellZone name")
    zone_cells = sum(
        block.cells[0] * block.cells[1] * block.cells[2]
        for block in zone_blocks
    )
    if zone_cells != 640000:
        raise RuntimeError(f"Expected 640000 largeReservoir cells, got {zone_cells}")

    face_owners: Counter[Tuple[int, ...]] = Counter()
    for block in mesh.blocks:
        lower_quad = [
            (mesh.vertices[block.vertices[k]][0], mesh.vertices[block.vertices[k]][2])
            for k in range(4)
        ]
        if signed_area_2d(lower_quad) >= -1.0e-18:
            raise RuntimeError(f"Non-positive block orientation: {block.vertices}")
        for face in block.faces():
            face_owners[tuple(sorted(face))] += 1

    if any(owner_count not in (1, 2) for owner_count in face_owners.values()):
        raise RuntimeError("A mesh face has an invalid owner count")

    external_faces = {face for face, count in face_owners.items() if count == 1}
    listed_faces: List[Tuple[int, ...]] = []
    for patch_faces in mesh.patches.values():
        listed_faces.extend(tuple(sorted(face)) for face in patch_faces)

    if len(listed_faces) != len(set(listed_faces)):
        raise RuntimeError("A boundary face is listed more than once")
    if set(listed_faces) != external_faces:
        missing = external_faces.difference(listed_faces)
        extra = set(listed_faces).difference(external_faces)
        raise RuntimeError(f"Boundary mismatch: missing={len(missing)}, extra={len(extra)}")

    expected_patch_faces = {
        "bottomOutlet": 400,
        "pistonWall": 6400,
        "contractionWall": 6000,
        "smallWall": 80,
        "largeWall": 320,
    }
    actual_patch_faces = {name: len(faces) for name, faces in mesh.patches.items()}
    if actual_patch_faces != expected_patch_faces:
        raise RuntimeError(
            f"Unexpected patch face counts: {actual_patch_faces}, expected {expected_patch_faces}"
        )

    actual_patch_areas = {
        name: sum(
            face_area_3d([mesh.vertices[vertex_i] for vertex_i in face])
            for face in faces
        )
        for name, faces in mesh.patches.items()
    }
    polygon_area = lambda apothem: (
        N_SIDES * apothem**2 * math.tan(math.pi / N_SIDES)
    )
    perimeter = lambda apothem: (
        2.0 * N_SIDES * apothem * math.tan(math.pi / N_SIDES)
    )
    expected_patch_areas = {
        "bottomOutlet": polygon_area(SMALL_APOTHEM),
        "pistonWall": polygon_area(LARGE_APOTHEM),
        "contractionWall": polygon_area(LARGE_APOTHEM) - polygon_area(SMALL_APOTHEM),
        "smallWall": perimeter(SMALL_APOTHEM) * (Y_CONTRACTION - Y_BOTTOM),
        "largeWall": perimeter(LARGE_APOTHEM) * (Y_TOP - Y_CONTRACTION),
    }
    for patch_name, expected_area in expected_patch_areas.items():
        if abs(actual_patch_areas[patch_name] - expected_area) > 1.0e-11:
            raise RuntimeError(
                f"Patch area mismatch for {patch_name}: "
                f"{actual_patch_areas[patch_name]} vs {expected_area}"
            )


def regular_polygon_vertices(apothem: float) -> Tuple[Point2, ...]:
    radius = apothem / math.cos(math.pi / N_SIDES)
    return tuple(
        (
            snap(radius * math.cos(math.pi / N_SIDES + k * 2.0 * math.pi / N_SIDES)),
            snap(radius * math.sin(math.pi / N_SIDES + k * 2.0 * math.pi / N_SIDES)),
        )
        for k in range(N_SIDES)
    )


def facet_endpoints(apothem: float, facet_i: int) -> Tuple[Point2, Point2]:
    theta = facet_i * 2.0 * math.pi / N_SIDES
    radius = apothem / math.cos(math.pi / N_SIDES)
    minus = (
        snap(radius * math.cos(theta - math.pi / N_SIDES)),
        snap(radius * math.sin(theta - math.pi / N_SIDES)),
    )
    plus = (
        snap(radius * math.cos(theta + math.pi / N_SIDES)),
        snap(radius * math.sin(theta + math.pi / N_SIDES)),
    )
    return minus, plus


def build_walls() -> List[Wall]:
    walls: List[Wall] = []

    for region, apothem, y_lo, y_hi in (
        ("smallSide", SMALL_APOTHEM, Y_BOTTOM, Y_CONTRACTION),
        ("largeSide", LARGE_APOTHEM, Y_CONTRACTION, Y_TOP),
    ):
        for facet_i, (nx, nz) in enumerate(FACET_NORMALS_2D):
            minus, plus = facet_endpoints(apothem, facet_i)
            walls.append(
                Wall(
                    f"{region}_{facet_i:02d}",
                    (nx, 0.0, nz),
                    (apothem * nx, 0.0, apothem * nz),
                    (
                        (minus[0], y_lo, minus[1]),
                        (plus[0], y_lo, plus[1]),
                        (plus[0], y_hi, plus[1]),
                        (minus[0], y_hi, minus[1]),
                    ),
                    region,
                )
            )

    for facet_i in range(N_SIDES):
        inner_minus, inner_plus = facet_endpoints(SMALL_APOTHEM, facet_i)
        outer_minus, outer_plus = facet_endpoints(LARGE_APOTHEM, facet_i)
        walls.append(
            Wall(
                f"shoulder_{facet_i:02d}",
                (0.0, -1.0, 0.0),
                (0.0, Y_CONTRACTION, 0.0),
                (
                    (inner_minus[0], Y_CONTRACTION, inner_minus[1]),
                    (outer_minus[0], Y_CONTRACTION, outer_minus[1]),
                    (outer_plus[0], Y_CONTRACTION, outer_plus[1]),
                    (inner_plus[0], Y_CONTRACTION, inner_plus[1]),
                ),
                "shoulder",
            )
        )

    piston_vertices = tuple(
        (x, Y_TOP, z) for x, z in regular_polygon_vertices(LARGE_APOTHEM)
    )
    walls.append(
        Wall(
            "piston",
            (0.0, 1.0, 0.0),
            (0.0, Y_TOP, 0.0),
            piston_vertices,
            "piston",
        )
    )

    validate_walls(walls)
    return walls


def validate_walls(walls: Sequence[Wall]) -> None:
    if len(walls) != 37:
        raise RuntimeError(f"Expected 37 DEM walls, got {len(walls)}")

    fast_path_count = 0
    shoulder_area = 0.0
    small_side_area = 0.0
    large_side_area = 0.0
    piston_area = 0.0

    for wall in walls:
        normal_length = magnitude(wall.normal)
        if abs(normal_length - 1.0) > 1.0e-12:
            raise RuntimeError(f"Non-unit normal for {wall.name}: {normal_length}")

        for vertex in wall.vertices:
            plane_distance = dot(vector_sub(vertex, wall.plane_point), wall.normal)
            if abs(plane_distance) > 1.0e-12:
                raise RuntimeError(
                    f"Non-coplanar vertex in {wall.name}: distance={plane_distance}"
                )

        # Check convexity using a local 2-D basis derived from the wall normal.
        tangent1 = vector_sub(wall.vertices[1], wall.vertices[0])
        tangent1 = vector_sub(tangent1, vector_scale(dot(tangent1, wall.normal), wall.normal))
        tangent1 = vector_scale(1.0 / magnitude(tangent1), tangent1)
        tangent2 = cross(wall.normal, tangent1)
        local = [
            (
                dot(vector_sub(vertex, wall.plane_point), tangent1),
                dot(vector_sub(vertex, wall.plane_point), tangent2),
            )
            for vertex in wall.vertices
        ]
        orientation = 1.0 if signed_area_2d(local) > 0.0 else -1.0
        for i in range(len(local)):
            a = local[i]
            b = local[(i + 1) % len(local)]
            c = local[(i + 2) % len(local)]
            turn = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
            if orientation * turn < -1.0e-13:
                raise RuntimeError(f"Non-convex polygon: {wall.name}")

        nonzero_normal_components = sum(abs(value) > 1.0e-12 for value in wall.normal)
        if wall.kind in ("smallSide", "largeSide") and nonzero_normal_components == 1:
            fast_path_count += 1

        if wall.kind == "shoulder":
            shoulder_area += polygon_area_3d(wall.vertices, wall.normal)
        elif wall.kind == "smallSide":
            small_side_area += polygon_area_3d(wall.vertices, wall.normal)
        elif wall.kind == "largeSide":
            large_side_area += polygon_area_3d(wall.vertices, wall.normal)
        elif wall.kind == "piston":
            piston_area += polygon_area_3d(wall.vertices, wall.normal)

    if fast_path_count != 8:
        raise RuntimeError(f"Expected 8 axis-aligned side fast paths, got {fast_path_count}")

    expected_shoulder_area = N_SIDES * (
        LARGE_APOTHEM**2 - SMALL_APOTHEM**2
    ) * math.tan(math.pi / N_SIDES)
    if abs(shoulder_area - expected_shoulder_area) > 1.0e-13:
        raise RuntimeError(
            f"Shoulder area mismatch: {shoulder_area} vs {expected_shoulder_area}"
        )

    perimeter = lambda apothem: (
        2.0 * N_SIDES * apothem * math.tan(math.pi / N_SIDES)
    )
    expected_areas = {
        "small sides": perimeter(SMALL_APOTHEM) * (Y_CONTRACTION - Y_BOTTOM),
        "large sides": perimeter(LARGE_APOTHEM) * (Y_TOP - Y_CONTRACTION),
        "piston": N_SIDES * LARGE_APOTHEM**2 * math.tan(math.pi / N_SIDES),
    }
    actual_areas = {
        "small sides": small_side_area,
        "large sides": large_side_area,
        "piston": piston_area,
    }
    for group, expected_area in expected_areas.items():
        if abs(actual_areas[group] - expected_area) > 1.0e-13:
            raise RuntimeError(
                f"DEM area mismatch for {group}: "
                f"{actual_areas[group]} vs {expected_area}"
            )


def render_block_mesh_dict(mesh: MeshData) -> str:
    lines: List[str] = [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "| Generated regular-dodecagon moving-piston contraction mesh                |",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        "    object      blockMeshDict;",
        "}",
        "",
        "convertToMeters 1;",
        "",
        "// Generated by generate_dodecagon_geometry.py. Do not edit by hand.",
        f"// Dh_small={fmt(2*SMALL_APOTHEM)}, Dh_large={fmt(2*LARGE_APOTHEM)}",
        f"// vertices={len(mesh.vertices)}, blocks={len(mesh.blocks)}, cells={mesh.total_cells}",
        "",
        "vertices",
        "(",
    ]

    for vertex_i, vertex in enumerate(mesh.vertices):
        lines.append(f"    {fmt_point(vertex)} // {vertex_i}")

    lines.extend((");", "", "blocks", "("))
    for block in mesh.blocks:
        vertex_text = " ".join(str(vertex) for vertex in block.vertices)
        cell_text = " ".join(str(count) for count in block.cells)
        zone_text = f" {block.zone}" if block.zone else ""
        lines.append(
            f"    hex ({vertex_text}){zone_text} ({cell_text}) "
            "simpleGrading (1 1 1)"
        )

    lines.extend((");", "", "edges", "(", ");", "", "boundary", "("))

    patch_types = {
        "bottomOutlet": "patch",
        "pistonWall": "wall",
        "contractionWall": "wall",
        "smallWall": "wall",
        "largeWall": "wall",
    }
    for patch_name in (
        "bottomOutlet",
        "pistonWall",
        "contractionWall",
        "smallWall",
        "largeWall",
    ):
        lines.extend(
            (
                f"    {patch_name}",
                "    {",
                f"        type {patch_types[patch_name]};",
                "        faces",
                "        (",
            )
        )
        for face in mesh.patches[patch_name]:
            lines.append("            (" + " ".join(str(vertex) for vertex in face) + ")")
        lines.extend(("        );", "    }", ""))

    lines.extend((");", "", "mergePatchPairs", "(", ");", "", "// ************************************************************************* //", ""))
    return "\n".join(lines)


def render_collision_patches(walls: Sequence[Wall]) -> str:
    lines = [
        "// Generated by generate_dodecagon_geometry.py. Do not edit by hand.",
        "// 12 small sides + 12 large sides + 12 shoulders + 1 piston = 37.",
        "// There is intentionally no DEM wall at the CFD bottomOutlet.",
        "",
    ]

    for wall in walls:
        lines.extend(
            (
                f"{wall.name}",
                "{",
                "    material    wall;",
                f"    nVec        {fmt_point(wall.normal)};",
                f"    planePoint  {fmt_point(wall.plane_point)};",
                "",
                "    finitePatch true;",
                "    vertices",
                "    (",
            )
        )
        for vertex in wall.vertices:
            lines.append(f"        {fmt_point(vertex)}")
        lines.append("    );")

        if wall.kind == "largeSide":
            lines.extend(
                (
                    "",
                    "    supportLimit",
                    "    {",
                    "        type  followWall;",
                    "        wall  piston;",
                    "        bound max;",
                    "    }",
                )
            )
        elif wall.kind == "piston":
            lines.extend(
                (
                    "",
                    "    motion",
                    "    {",
                    "        type  meshPatchTranslation;",
                    "        patch pistonWall;",
                    "    }",
                )
            )

        lines.extend(("}", ""))

    return "\n".join(lines)


def render_addition_domain() -> str:
    """Render the upstream prism used by onceScatter.

    Keeping these vertices in the geometry generator makes the CFD boundary,
    DEM walls, and particle-addition domain share one source of truth.
    """

    lines = [
        "// Generated by generate_dodecagon_geometry.py. Do not edit by hand.",
        "// Full large dodecagon from y=0.05 m to y=0.15 m.",
        "addDomain polygonPrism;",
        "polygonPrismCoeffs",
        "{",
        "    axis       (0 1 0);",
        f"    minAxial  {fmt(Y_CONTRACTION)};",
        f"    maxAxial  {fmt(Y_TOP)};",
        "    vertices",
        "    (",
    ]

    for x, z in regular_polygon_vertices(LARGE_APOTHEM):
        lines.append(f"        {fmt_point((x, Y_CONTRACTION, z))}")

    lines.extend(("    );", "}", ""))
    return "\n".join(lines)


def render_topo_set_dict() -> str:
    margin = 1.0e-6
    lines = [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "| Generated dynamic-mesh face zones                                        |",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        "    object      topoSetDict;",
        "}",
        "",
        "// Generated by generate_dodecagon_geometry.py. Do not edit by hand.",
        "actions",
        "(",
        "    {",
        "        name    reservoirTopFaces;",
        "        type    faceSet;",
        "        action  new;",
        "        source  patchToFace;",
        "        sourceInfo",
        "        {",
        "            name pistonWall;",
        "        }",
        "    }",
        "",
        "    {",
        "        name    reservoirTop;",
        "        type    faceZoneSet;",
        "        action  new;",
        "        source  setToFaceZone;",
        "        sourceInfo",
        "        {",
        "            set reservoirTopFaces;",
        "        }",
        "    }",
        "",
        "    {",
        "        name    reservoirBottomFaces;",
        "        type    faceSet;",
        "        action  new;",
        "        source  boxToFace;",
        "        sourceInfo",
        "        {",
        "            box",
        f"                {fmt_point((-LARGE_APOTHEM-margin, Y_CONTRACTION-margin, -LARGE_APOTHEM-margin))}",
        f"                {fmt_point(( LARGE_APOTHEM+margin, Y_CONTRACTION+margin,  LARGE_APOTHEM+margin))};",
        "        }",
        "    }",
        "",
        "    {",
        "        name    reservoirBottom;",
        "        type    faceZoneSet;",
        "        action  new;",
        "        source  setToFaceZone;",
        "        sourceInfo",
        "        {",
        "            set reservoirBottomFaces;",
        "        }",
        "    }",
        ");",
        "",
        "// ************************************************************************* //",
        "",
    ]
    return "\n".join(lines)


def write_if_changed(path: Path, content: str) -> bool:
    if path.exists() and path.read_text() == content:
        return False
    path.write_text(content)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate geometry and verify that generated files are current",
    )
    args = parser.parse_args()

    case_dir = Path(__file__).resolve().parent
    block_mesh_path = case_dir / "system" / "blockMeshDict"
    collision_path = case_dir / "constant" / "dodecagonCollisionPatches"
    addition_domain_path = case_dir / "constant" / "dodecagonAdditionDomain"
    topo_set_path = case_dir / "system" / "topoSetDict"

    mesh = build_mesh()
    walls = build_walls()
    generated = {
        block_mesh_path: render_block_mesh_dict(mesh),
        collision_path: render_collision_patches(walls),
        addition_domain_path: render_addition_domain(),
        topo_set_path: render_topo_set_dict(),
    }

    if args.check:
        stale = [path for path, content in generated.items() if not path.exists() or path.read_text() != content]
        if stale:
            for path in stale:
                print(f"STALE: {path.relative_to(case_dir)}")
            return 1
        action = "validated"
    else:
        changed = [path for path, content in generated.items() if write_if_changed(path, content)]
        action = "generated" if changed else "already current"

    print(
        f"Dodecagon geometry {action}: "
        f"{len(mesh.blocks)} blocks, {mesh.total_cells} cells, "
        f"{len(walls)} DEM walls (12 polygon followers + moving piston; "
        "8 side walls on the legacy fast path)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
