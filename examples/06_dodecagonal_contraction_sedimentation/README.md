# Dodecagonal contraction sedimentation

This case preserves `examples/05_contraction_sedimentation` and replaces its
rectangular cross-section with a regular 12-sided approximation of a circular
4:1 contraction.

## Geometry

For a regular polygon with apothem `a`,

`Dh = 4 A / P = 2 a`.

The selected apothems therefore preserve the original hydraulic diameters:

| Region | y range (m) | Apothem (m) | Hydraulic diameter (m) |
|---|---:|---:|---:|
| Small channel | 0.00 to 0.05 | 0.01 | 0.02 |
| Large chamber | 0.05 to 0.15 | 0.04 | 0.08 |

The CFD boundary uses five aggregated patches:

- `smallWall`: all 12 downstream side facets
- `largeWall`: all 12 upstream side facets
- `contractionWall`: the complete annular shoulder
- `bottomWall`: closed downstream end
- `topOpen`: open upstream end

DEM contact uses 37 exact convex finite polygons: 12 small side rectangles,
12 large side rectangles, 12 shoulder trapezoids, and one 12-vertex bottom
polygon. The count is 37, not 38, because `topOpen` has no collision wall.

`generate_dodecagon_geometry.py` is the single geometry source for both
`system/blockMeshDict` and `constant/dodecagonCollisionPatches`. The logical
80 by 80 large-section grid maps onto the dodecagon, and its central 20 by 20
grid is shared conformally with the small channel. The resulting mesh has
6800 block definitions and 660000 run-time cells, the same cell count as case
05. One mapped cross-sectional cell per block avoids hidden interpolation or
non-conformal faces at the polygon corners and contraction opening.

## Particle initialization

`onceScatter` now supports an ordered, convex `polygonPrism` addition domain.
This case uses the complete large dodecagon from `y=0.05 m` to `y=0.15 m`.
`constant/dodecagonAdditionDomain` is generated from the same vertices as the
CFD mesh and DEM walls, then included by `HFDIBDEMDict`.

The input form is:

```foam
addDomain polygonPrism;
polygonPrismCoeffs
{
    axis       (0 1 0);
    minAxial  0.05;
    maxAxial  0.15;
    vertices  ( /* ordered coplanar convex vertices */ );
}
```

`minAxial` and `maxAxial` are coordinates along the normalized `axis`, i.e.
`point & axis`. Vertices may use either winding but must be contiguous,
coplanar, and convex. This makes the feature reusable for triangles,
trapezoids, and general convex polygonal cross-sections, including an axis not
aligned with a Cartesian direction.

For each attempt, `onceScatter` first applies the existing random rotation and
scale. It then erodes every polygon half-plane by the actual directional
support of the enlarged body and erodes the axial interval by its axial
extents. A center is drawn uniformly from the remaining prism using
area-weighted fan triangulation and square-root barycentric sampling. Therefore
all surface vertices of an arbitrary STL body, not only its center or AABB,
are inside the domain before the existing particle-particle/wall rejection is
run. Spheres use their analytic radius; a geometry type that exposes no
surface points uses a conservative bounding-box fallback.

The `nSolidsInDomain 1000` value is an intentional validation-stage cap, while
`fieldValue 0.3` remains the later production target. With the exact large
prism volume (`5.144624494677556e-4 m3`) and the current STL volume
(`about 2.0279e-8 m3`), 1000 non-overlapping particles correspond to only about
3.94% nominal geometric volume. An ideal no-overlap estimate for 30% is about
7611 particles; the actual required cap should be determined from the
discretized `lambda` fraction and packing/rejection behavior.

## Generate, validate, and run

From this case directory:

```sh
./generate_dodecagon_geometry.py
./Allcheck
./Allrun
```

`Allrun` regenerates the geometry, runs `blockMesh`, requires the full
`checkMesh -allGeometry -allTopology` check to pass, then runs the solver from
`system/controlDict`.

`Allcheck` always runs the deterministic static checks. When OpenFOAM commands
are on `PATH`, it also runs `blockMesh` and `checkMesh`. The predicate
regression compares the old eight-AABB-corner formulation with the optimized
center/support-radius formulation for every wall, including all 12 shoulder
trapezoids. The polygon-prism regression additionally checks 120000 point
samples for area weighting, 600 randomly rotated/scaled placements of the case
STL for whole-body containment, analytic sphere clearance, reversed winding,
and several other convex polygon shapes.

## Source changes in this package

The new `polygonPrismGeometry` kernel validates and caches the prism basis and
convex edge half-planes. `addModelOnceScatter` uses it only when
`addDomain polygonPrism` is selected; the existing `cellZone` and `boundBox`
placement paths are retained. On parallel runs, the new path draws random
rotation, scale, and position values on the master rank and synchronizes them,
so all ranks construct the same candidate. Addition-zone cells are found once
by a linear cell-center scan instead of invoking the legacy bounding-box graph
walk.

The general finite-polygon path now caches each exact polygon-edge half-space
and evaluates AABB extrema analytically. It also reuses an exact wall-plane
distance range in the finite contact-area path. These changes remove repeated
corner transforms, square roots, and a duplicate eight-corner loop.

The legacy axis-aligned rectangle path is unchanged; eight of the 24 side
facets automatically use it. The remaining 29 finite polygons use the general
path. Polygon clipping, fine-cell enumeration, particle volume queries,
64-point mixed-cell quadrature, contact area/volume accumulation, and force
models are unchanged.

### Existing coplanar-seam limitation

The 12 shoulder trapezoids tile the contraction plane without finite-area
overlap. Their shared radial edges therefore do not create a geometric gap or
double-counted surface (apart from the existing predicate tolerance at a
measure-zero edge). However, the current OpenHFDIB-DEM finite-wall architecture
creates an independent sub-contact for each collision patch. A particle whose
contact footprint crosses a radial shoulder seam can consequently have two
coplanar sub-contacts instead of one monolithic annular contact. Depending on
the force/history model, that is not guaranteed to be numerically identical to
a single wall with a hole.

This package does not attempt a riskier coplanar-patch force aggregation. For
acceptance, compare a single-particle drop at a shoulder-facet center with a
drop centered on a radial shoulder seam, in addition to side-facet and
dodecagon-corner contacts. The geometry predicate regression covers both
locations, while the target OpenFOAM run is needed to assess force continuity.

After applying the source files, rebuild both the library and solver in the
same OpenFOAM environment used for the project:

```sh
cd src/HFDIBDEM
wmake libso

cd ../../applications/solvers/incompressible/pimpleLYJHFDIBFoam
wmake
```

The included Python checks do not replace compilation and `checkMesh` on the
target OpenFOAM installation.
