#!/usr/bin/env python3
"""Run static consistency checks for the dodecagonal piston case."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import List, Set

from generate_dodecagon_geometry import (
    LARGE_APOTHEM,
    N_SIDES,
    Y_CONTRACTION,
    Y_TOP,
    build_mesh,
    build_walls,
    regular_polygon_vertices,
)


CASE_DIR = Path(__file__).resolve().parent


def without_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*?$", "", text, flags=re.MULTILINE)


def check_balanced(path: Path) -> None:
    text = without_comments(path.read_text())
    pairs = {"}": "{", ")": "("}
    stack: List[str] = []
    for char in text:
        if char in "{(":
            stack.append(char)
        elif char in "})":
            if not stack or stack.pop() != pairs[char]:
                raise AssertionError(f"Unbalanced delimiters in {path}")
    if stack:
        raise AssertionError(f"Unbalanced delimiters in {path}")


def boundary_patch_names(path: Path) -> Set[str]:
    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*|[{}]", without_comments(path.read_text()))
    try:
        index = tokens.index("boundaryField")
        index = tokens.index("{", index + 1)
    except ValueError as error:
        raise AssertionError(f"No boundaryField dictionary in {path}") from error

    patches: Set[str] = set()
    depth = 1
    index += 1
    while index < len(tokens) and depth:
        token = tokens[index]
        if token == "{":
            depth += 1
        elif token == "}":
            depth -= 1
        elif depth == 1 and index + 1 < len(tokens) and tokens[index + 1] == "{":
            patches.add(token)
        index += 1
    return patches


def main() -> int:
    mesh = build_mesh()
    walls = build_walls()
    expected_cfd_patches = set(mesh.patches)

    dictionary_paths = [
        CASE_DIR / "constant" / "HFDIBDEMDict",
        CASE_DIR / "constant" / "dodecagonAdditionDomain",
        CASE_DIR / "constant" / "dodecagonCollisionPatches",
        CASE_DIR / "constant" / "dynamicMeshDict",
        CASE_DIR / "system" / "blockMeshDict",
        CASE_DIR / "system" / "controlDict",
        CASE_DIR / "system" / "topoSetDict",
    ]
    dictionary_paths.extend(
        CASE_DIR / "0.org" / name
        for name in ("U", "p", "f", "lambda", "pointDisplacement")
    )
    for path in dictionary_paths:
        check_balanced(path)

    for field_name in ("U", "p", "f", "lambda", "pointDisplacement"):
        path = CASE_DIR / "0.org" / field_name
        actual = boundary_patch_names(path)
        if actual != expected_cfd_patches:
            raise AssertionError(
                f"{field_name} patches {sorted(actual)} != mesh patches "
                f"{sorted(expected_cfd_patches)}"
            )

    hfdib_text = (CASE_DIR / "constant" / "HFDIBDEMDict").read_text()
    if '#include "dodecagonCollisionPatches"' not in hfdib_text:
        raise AssertionError("HFDIBDEMDict does not include generated collision walls")
    if '#include "dodecagonAdditionDomain"' not in hfdib_text:
        raise AssertionError("HFDIBDEMDict does not include generated addition domain")
    if "addDomain boundBox;" in without_comments(hfdib_text):
        raise AssertionError("Legacy boundBox addition domain remains active")
    if not re.search(r"\bfieldValue\s+0\.3\s*;", without_comments(hfdib_text)):
        raise AssertionError("The intended fieldValue 0.3 target changed")
    if not re.search(r"\bnSolidsInDomain\s+2000\s*;", without_comments(hfdib_text)):
        raise AssertionError("The preserved example-04 2000-particle cap changed")
    if not re.search(r"\brandomSeed\s+1207\s*;", without_comments(hfdib_text)):
        raise AssertionError("The deterministic case-07 random seed changed")

    addition_text = without_comments(
        (CASE_DIR / "constant" / "dodecagonAdditionDomain").read_text()
    )
    if not re.search(r"\baddDomain\s+polygonPrism\s*;", addition_text):
        raise AssertionError("Generated addition domain is not polygonPrism")
    if not re.search(r"\baxis\s*\(0\s+1\s+0\)\s*;", addition_text):
        raise AssertionError("Unexpected polygonPrism axis")

    axial_min = re.search(r"\bminAxial\s+([-+0-9.eE]+)\s*;", addition_text)
    axial_max = re.search(r"\bmaxAxial\s+([-+0-9.eE]+)\s*;", addition_text)
    if not axial_min or not axial_max:
        raise AssertionError("Missing polygonPrism axial bounds")
    if abs(float(axial_min.group(1)) - Y_CONTRACTION) > 1.0e-14:
        raise AssertionError("Unexpected polygonPrism minAxial")
    if abs(float(axial_max.group(1)) - Y_TOP) > 1.0e-14:
        raise AssertionError("Unexpected polygonPrism maxAxial")

    vertex_block = re.search(r"\bvertices\s*\((.*?)\)\s*;", addition_text, re.DOTALL)
    if not vertex_block:
        raise AssertionError("Missing polygonPrism vertices")
    parsed_vertices = [
        tuple(float(value) for value in match)
        for match in re.findall(
            r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)",
            vertex_block.group(1),
        )
    ]
    expected_vertices = [
        (x, Y_CONTRACTION, z)
        for x, z in regular_polygon_vertices(LARGE_APOTHEM)
    ]
    if len(parsed_vertices) != N_SIDES:
        raise AssertionError(f"Expected {N_SIDES} addition vertices")
    for actual, expected in zip(parsed_vertices, expected_vertices):
        if max(abs(a - e) for a, e in zip(actual, expected)) > 1.0e-13:
            raise AssertionError("Addition-domain vertices drifted from mesh geometry")

    addition_volume = (
        N_SIDES
        * LARGE_APOTHEM**2
        * math.tan(math.pi / N_SIDES)
        * (Y_TOP - Y_CONTRACTION)
    )
    if abs(addition_volume - 5.144624494677556e-4) > 1.0e-16:
        raise AssertionError("Unexpected exact addition-domain volume")

    collision_text = (CASE_DIR / "constant" / "dodecagonCollisionPatches").read_text()
    generated_names = set(re.findall(r"^([A-Za-z][A-Za-z0-9_]*)\n\{", collision_text, re.MULTILINE))
    expected_names = {wall.name for wall in walls}
    if generated_names != expected_names:
        raise AssertionError("Generated DEM wall names do not match geometry model")
    if collision_text.count("finitePatch true;") != 37:
        raise AssertionError("Expected exactly 37 finite DEM collision walls")
    if collision_text.count("type  followWall;") != 12:
        raise AssertionError("Expected 12 piston-limited large side walls")
    if collision_text.count("wall  piston;") != 12:
        raise AssertionError("Every large side wall must follow piston")
    if collision_text.count("type  meshPatchTranslation;") != 1:
        raise AssertionError("Expected one moving DEM wall")
    if "patch pistonWall;" not in collision_text:
        raise AssertionError("DEM piston does not follow CFD pistonWall")
    if re.search(r"^bottom(?:Wall|Outlet)\s*\n\{", collision_text, re.MULTILINE):
        raise AssertionError("A DEM collision wall remains at the outlet")

    legacy_names = {
        "topOpen",
        "bottomWall",
        "xMinWall",
        "xMaxWall",
        "zMinWall",
        "zMaxWall",
    }
    for field_name in ("U", "p", "f", "lambda", "pointDisplacement"):
        if legacy_names.intersection(
            boundary_patch_names(CASE_DIR / "0.org" / field_name)
        ):
            raise AssertionError(f"A rectangular legacy patch remains in {field_name}")

    block_text = without_comments((CASE_DIR / "system" / "blockMeshDict").read_text())
    if len(re.findall(r"\bhex\s*\([^\n]+\)\s+largeReservoir\s*\(", block_text)) != 6400:
        raise AssertionError("Expected 6400 blocks in largeReservoir cellZone")

    topo_text = without_comments((CASE_DIR / "system" / "topoSetDict").read_text())
    for required_name in (
        "reservoirTopFaces",
        "reservoirTop",
        "reservoirBottomFaces",
        "reservoirBottom",
    ):
        if not re.search(rf"\bname\s+{required_name}\s*;", topo_text):
            raise AssertionError(f"Missing dynamic-mesh set/zone {required_name}")
    if not re.search(r"\bname\s+pistonWall\s*;", topo_text):
        raise AssertionError("reservoirTop is not sourced from pistonWall")

    dynamic_text = without_comments(
        (CASE_DIR / "constant" / "dynamicMeshDict").read_text()
    )
    for required_name in ("largeReservoir", "reservoirTop", "reservoirBottom"):
        if not re.search(rf"\b{required_name}\b", dynamic_text):
            raise AssertionError(f"dynamicMeshDict misses {required_name}")

    control_text = without_comments((CASE_DIR / "system" / "controlDict").read_text())
    if not re.search(r"\bpatches\s*\(\s*pistonWall\s*\)\s*;", control_text):
        raise AssertionError("Fluid force function does not use pistonWall")

    allrun_text = (CASE_DIR / "Allrun").read_text()
    required_commands = (
        "runApplication blockMesh",
        "runApplication topoSet",
        "runApplication checkMesh -allGeometry -allTopology",
    )
    if any(command not in allrun_text for command in required_commands):
        raise AssertionError("Allrun misses a required mesh/zone validation command")

    motion_text = (CASE_DIR / "constant" / "reservoirTopDisplacement.dat").read_text()
    if "(0  (0  0     0))" not in motion_text or "(90 (0 -0.09  0))" not in motion_text:
        raise AssertionError("Example-04 piston displacement history changed")

    print(
        "Piston case consistent: 5 CFD patches, 37 DEM walls with no outlet "
        "collision wall, 640000 dynamic-zone cells, 12 polygon followers, "
        "and a moving dodecagonal piston."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
