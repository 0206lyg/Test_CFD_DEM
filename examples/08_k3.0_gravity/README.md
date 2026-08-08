# Pure-DEM icosahedron clogging in a rectangular contraction

This case isolates gravity-driven DEM contact and clogging from all fluid
equations. It runs with `HFDIBDEMFoam` and the same `openHFDIBDEM` contact
library used by the coupled solvers.

## Geometry

- Axis-aligned square reservoir: 40 mm x 40 mm, y = 25--75 mm
- Axis-aligned square throat: 10 mm x 10 mm, y = 0--25 mm
- Linear contraction ratio: 4:1
- Original `PH_SC.stl`, without scaling
- Particle volume-equivalent diameter: approximately 3.383 mm
- Opening ratio: D/d_v = 10/3.383 = 2.96 (approximately 3)
- Uniform Eulerian mesh spacing: 1 mm
- Total mesh size: 82,500 cells

The lower `bottomOutlet` and upper `topOpen` mesh patches are deliberately
absent from `DEM/collisionPatches`. Only the 12 side/shoulder patches collide
with particles.

## Runtime inlet

The existing `repeatRandomPosition` model maintains a local particle volume
fraction of `fieldValue 0.20` in the top inlet slab
`(-0.020, 0.060, -0.020)`--`(0.020, 0.075, 0.020)` m. This 15 mm high slab
contains approximately 237 unscaled icosahedra at a geometrical solid fraction
of 0.20. It is filled at time zero and replenished as particles fall out of the
slab.

Particles are monodisperse (`noScaling`) and use the new
`uniformRandomRotation` mode. It samples a uniform three-dimensional
orientation while retaining the model's original clock-seeded random-number
generator. The existing `distribution` addModel and `distributionDict` format
are unchanged and remain available for future polydisperse cases.

If the inlet becomes too crowded to accept another non-overlapping particle,
the existing addition loop abandons the current replenishment pass after 50
consecutive rejected attempts. The solver then advances normally and may try
again at the next time step; there is no runtime active-particle cap in this
case.

## Build and run

Build `src/HFDIBDEM` first, then rebuild
`applications/solvers/pureDEM/HFDIBDEMFoam`. Its `Make/options` links against
`liblyjHFDIBDEM`.

From the repository root:

```sh
(cd src/HFDIBDEM && wclean && wmake libso) && \
(cd applications/solvers/pureDEM/HFDIBDEMFoam && wclean && wmake)

ldd "$FOAM_USER_APPBIN/HFDIBDEMFoam" | grep liblyjHFDIBDEM
```

From this case directory:

```sh
./Allrun
```

The configured DEM substep is `stepDEM*deltaT = 1e-5 s`. Treat this as a
performance candidate rather than an already converged contact timestep. The
default 0.3 s run is an initial stability test; compare a short interval with a
`1e-6 s` reference before interpreting discharge or clogging results.
