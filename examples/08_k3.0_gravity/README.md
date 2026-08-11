# Parameterized pure-DEM clogging in a rectangular contraction

This case isolates gravity-driven DEM contact and clogging from all fluid
equations. It runs with `HFDIBDEMFoam` and the same `openHFDIBDEM` contact
library used by the coupled solvers. Particle shape, gravitational acceleration,
and the pore-to-particle size ratio are selected by the supplied workstation
submission script.

## Geometry

- Axis-aligned square reservoir: 40 mm x 40 mm, y = 25--75 mm
- Axis-aligned square throat: 10 mm x 10 mm, y = 0--25 mm
- Linear contraction ratio: 4:1
- Available shapes: `icosahedral.stl`, `cylinder.stl`, `cube.stl`, and an
  implicit center-radius sphere
- Common unscaled volume: approximately 20.278942696 mm3
- Common unscaled volume-equivalent diameter: 3.383365334148 mm
- Default `PH_SC.stl`: the icosahedral template
- Default scale: 0.98521235637496851
- Default volume-equivalent diameter after scaling: 3.333333333333 mm
- Default pore-to-particle ratio: k = 10/3.333333333333 = 3.00
- Uniform Eulerian mesh spacing: 1 mm
- Total mesh size: 82,500 cells

The three template meshes are closed, convex, and equal-volume to numerical
precision. The implicit sphere is normalized through its volume-equivalent
diameter, so all four shapes use the same definition of `k`. The `cylinder`
template is a regular octagonal prism with a bounding length-to-diameter ratio
of 3.

## Particle-size convention

The size ratio is defined using the square throat width and the
volume-equivalent particle diameter:

```text
k = D_pore / d_eq
d_eq,target = D_pore / k
scale = D_pore / (k d_eq,0)
R_sphere,target = D_pore / (2k)
```

Here, `D_pore = 0.010 m` and the STL reference diameter is
`d_eq,0 = 0.003383365334148 m`. For an STL particle, the insertion model uses
`randomScaling` with identical `minScale` and `maxScale`; the submission script
overwrites both limits with the requested deterministic scale. For an implicit
sphere, the script writes the target radius directly and selects `noScaling`,
preventing the radius from being scaled a second time. Both routes are
monodisperse and satisfy the same requested `k` exactly.

Finite-wall contact for the implicit sphere uses the same virtual-mesh controls
as the STL path. In particular, increasing `virtualMesh.level` or decreasing
`virtualMesh.charCellSize` refines finite-wall overlap integration. Analytic
sphere--sphere and infinite-wall contacts remain on their existing paths.

The lower `bottomOutlet` and upper `topOpen` mesh patches are deliberately
absent from `DEM/collisionPatches`. Only the 12 side/shoulder patches collide
with particles.

## Initial reservoir loading and runtime inlet

Two logical addition models use the same selected geometry. For an STL shape,
both use `constant/triSurface/PH_SC.stl`; `stlBaseName PH_SC` decouples the
addition-model names from the STL basename, so duplicate or symbolic-link STL
files are unnecessary. For a sphere, both use their local `sphere` dictionary
and the STL entry is ignored.

At time zero, `PH_SC_prefill` uses `onceScatter` to fill the complete upper
reservoir,
`(-0.020, 0.025, -0.020)`--`(0.020, 0.075, 0.020)` m, to `fieldValue 0.20`.
It never replenishes that region after time advances. `PH_SC_inlet` monitors
the upper 10 mm slab,
`(-0.020, 0.065, -0.020)`--`(0.020, 0.075, 0.020)` m, and replenishes it during
the run with `repeatRandomPosition`. The top slab is therefore a replenishment
zone nested inside the prefilled reservoir, not a second disjoint initial layer.
Neither source initializes the throat.

For the default exact `k = 3.00`, the geometric expectation is approximately
825 particles in the complete upper reservoir. The actual count is controlled
by the Eulerian `lambda` field and can differ slightly with shape, orientation,
and the final accepted particle. The `nSolidsInDomain 10000` value is a
non-binding constructor-time ceiling over the intended, mesh-resolved `k`
range; it is not a requested particle count or a runtime replenishment limit.

STL particles use `uniformRandomRotation`. The patched `onceScatter` model uses
Shoemake quaternion sampling, matching the existing uniform SO(3) proposal in
`repeatRandomPosition`. Overlap rejection can still bias the set of accepted
orientations in a crowded region; the orientation proposals themselves are
uniform. The sphere branch selects `noRotation`, since rotating a center-radius
sphere does not alter its geometry. `updateTorque true` is retained for its
frictional spin dynamics.

All particles start from rest (`velocity (0 0 0)` and `startSynced false`). A
coupled-flow comparison should use the same particle-velocity conditioning;
matching only the initial concentration does not make the startup states equal.

During time-zero initialization, 1000 consecutive rejected placements stop the
current model's constructor pass. During runtime inlet replenishment, 50
consecutive rejected placements stop only the current time-step pass; the
solver advances and may try the inlet again later.

## Workstation submission

Edit only the user-input section near the top of
`run_08_gravity_parametric.sh`. For the three STL shapes, the script copies the
selected template once to the run-local `constant/triSurface/PH_SC.stl` and
writes the same requested scale to both logical addition models. For `sphere`,
it selects `bodyGeom sphere`, writes `D_pore/(2k)` to both radius entries, and
selects `noScaling` and `noRotation`. It also writes gravity and synchronizes
`numberOfSubdomains` with `SLURM_NTASKS`. Run directories use this format:

```text
<particle>_k<k>_g<g>_<SLURM_JOB_ID>
```

For example, `sphere_k2.5_g2_14572831`.

## Build and run

The dictionary and Slurm-script changes in this control overlay do not require
compilation. Implicit-sphere contact with finite walls does require the
previous finite-wall sphere source overlay to be applied and the shared library
to be rebuilt once:

```sh
(cd src/HFDIBDEM && wmake libso)

ldd "$FOAM_USER_APPBIN/HFDIBDEMFoam" | grep liblyjHFDIBDEM
```

The public solver interface and ABI are unchanged, so an already-built
`HFDIBDEMFoam` does not need to be recompiled. Build the solver normally only
if it does not yet exist in `FOAM_USER_APPBIN`.

From this case directory:

```sh
./Allrun
```

Direct `./Allrun` uses the icosahedral `PH_SC.stl`, `k = 3.00`, and
`g = (0 -2 0) m/s2`. It performs the one-time complete-reservoir prefill and
then maintains the top inlet slab. The submission script overrides the shape,
size, and gravity in the run-local copy without changing the source case.

The configured DEM substep is `stepDEM*deltaT = 1e-5 s`. Treat this as a
performance candidate rather than an already converged contact timestep. The
default 1.0 s run should be compared over a short interval with a `1e-6 s`
reference before interpreting discharge or clogging results.
