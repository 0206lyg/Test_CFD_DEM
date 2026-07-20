# Dodecagonal piston flow

This case preserves `examples/04_piston_flow` and combines its moving-piston
setup with the regular-dodecagon contraction and `polygonPrism` particle
initialization introduced in case 06.

## Geometry and boundaries

The small and large regular dodecagons use apothems `0.01 m` and `0.04 m`.
For a regular polygon, `Dh = 4A/P = 2a`, so these values preserve the original
hydraulic diameters of `0.02 m` and `0.08 m`.

The CFD mesh has five aggregated patches:

- `smallWall`: 12 small-channel side facets
- `largeWall`: 12 large-reservoir side facets
- `contractionWall`: the 12-facet annular shoulder
- `bottomOutlet`: the open lower face
- `pistonWall`: the moving 12-sided upper face

DEM uses 37 exact convex finite walls: 12 small sides, 12 large sides, 12
shoulder trapezoids, and the moving dodecagonal piston. There is intentionally
no DEM collision wall at `bottomOutlet`.

`generate_dodecagon_geometry.py` is the shared source for
`system/blockMeshDict`, `system/topoSetDict`,
`constant/dodecagonCollisionPatches`, and
`constant/dodecagonAdditionDomain`. The mesh has 6800 block definitions and
660000 run-time cells. Its 6400 large-reservoir blocks, or 640000 cells, are
placed in the `largeReservoir` dynamic-mesh cell zone.

## Piston motion and finite DEM walls

The `pistonWall` displacement history is preserved from case 04: it moves from
`y=0.15 m` to `y=0.06 m`, while the contraction shoulder remains at
`y=0.05 m`. The `reservoirTop` and `reservoirBottom` face zones are generated
by `topoSet` and drive `displacementLayeredMotion`.

The 12 large DEM side polygons each use:

```foam
supportLimit
{
    type  followWall;
    wall  piston;
    bound max;
}
```

Polygon support is accepted only when the moving leader is axis-aligned, the
follower surface is parallel to that motion, every follower vertex lies on
one of two end caps, and each cap has at least two vertices. The solver also
stops before a moving cap crosses the fixed cap. For a valid wall, the existing
finite-wall bounds/vertices are updated and its cached polygon geometry is
rebuilt. Contact predicates, clipping, contact area/volume integration, and
wall-force formulas are unchanged by this support extension.

## Particle initialization

`onceScatter` samples directly in the complete large dodecagonal prism from
`y=0.05 m` to `y=0.15 m`. Randomly rotated and scaled arbitrary STL particles
are accepted only when their full support fits inside the prism. The case
retains the example-04 cap `nSolidsInDomain 2000`, the target
`fieldValue 0.3`, and uses deterministic `randomSeed 1207`.

## Particle-particle projection correction

`prtSubContactInfo.C` now computes the normal component of relative velocity
as `n*(relativeVelocity & n)`. Consequently the tangential component is
orthogonal to the contact normal. This is an intentional particle-particle
force-model correction and is independent of polygon-wall support.

## Generate, validate, and run

From this directory:

```sh
./generate_dodecagon_geometry.py
./Allcheck
./Allrun
```

`Allcheck` always runs deterministic geometry, dictionary, finite-contact,
polygon-sampling, moving-wall, and velocity-projection regressions. If
OpenFOAM mesh tools are on `PATH`, it also runs `blockMesh`, `topoSet`, and
`checkMesh -allGeometry -allTopology`. `Allrun` requires those mesh checks to
pass before starting the solver.

After applying the source overlay, rebuild the library and solver in the
target OpenFOAM environment:

```sh
cd src/HFDIBDEM
wmake libso

cd ../../applications/solvers/incompressible/pimpleLYJHFDIBFoam
wmake
```

The dodecagonal piston area is `0.00514462494677556 m2`, versus `0.0064 m2`
for the former square piston. Hydraulic diameter is preserved, but raw piston
force and swept volume should therefore not be compared as if the face areas
were identical.

The existing independent-subcontact behavior at coplanar shoulder seams is
unchanged. Acceptance testing should include particle impacts at a shoulder
facet center and directly over a radial seam, as well as piston-side contacts.

## Source files changed relative to overlay v2

- `src/HFDIBDEM/openHFDIBDEM.C`: safe polygonal `supportLimit` validation and
  moving-cap crossing guard
- `src/HFDIBDEM/contactModels/prtSubContactInfo.C`: correct normal-velocity
  projection for particle-particle tangential contact

The other source files included in the cumulative v3 overlay are unchanged
from v2.
