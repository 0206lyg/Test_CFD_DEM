# Parameterized pure-DEM clogging in a rectangular contraction

This case isolates gravity-driven DEM contact and clogging from all fluid
equations. It runs with `HFDIBDEMFoam` and the same `openHFDIBDEM` contact
library used by the coupled solvers. Particle shape, gravitational acceleration,
and the pore-to-particle size ratio can be selected by the external workstation
submission script without further library changes after the one-time patch in
this overlay is built.

## Geometry

- Axis-aligned square reservoir: 40 mm x 40 mm, y = 25--75 mm
- Axis-aligned square throat: 10 mm x 10 mm, y = 0--25 mm
- Linear contraction ratio: 4:1
- Available templates: `icosahedral.stl`, `cylinder.stl`, and `cube.stl`
- Common unscaled volume: approximately 20.278942696 mm3
- Common unscaled volume-equivalent diameter: 3.383365334148 mm
- Default `PH_SC.stl`: the icosahedral template
- Default scale: 0.98521235637496851
- Default volume-equivalent diameter after scaling: 3.333333333333 mm
- Default pore-to-particle ratio: k = 10/3.333333333333 = 3.00
- Uniform Eulerian mesh spacing: 1 mm
- Total mesh size: 82,500 cells

The three template meshes are closed, convex, and equal-volume to numerical
precision. The `cylinder` template is a regular octagonal prism with a bounding
length-to-diameter ratio of 3.

## Particle-size convention

The size ratio is defined using the square throat width and the
volume-equivalent particle diameter:

```text
k = D_pore / d_eq
d_eq,target = D_pore / k
scale = D_pore / (k d_eq,0)
```

Here, `D_pore = 0.010 m` and `d_eq,0 = 0.003383365334148 m`. The insertion
model uses `randomScaling` with identical `minScale` and `maxScale`, which
makes the scale deterministic and keeps all particles monodisperse. The
workstation submission script overwrites both limits with the requested value.

The lower `bottomOutlet` and upper `topOpen` mesh patches are deliberately
absent from `DEM/collisionPatches`. Only the 12 side/shoulder patches collide
with particles.

## Initial reservoir loading and runtime inlet

Two logical addition models use the same `constant/triSurface/PH_SC.stl` file.
The optional `stlBaseName PH_SC` entry decouples the addition-model name from
the STL basename; no duplicate or symbolic-link STL files are required.

At time zero, `PH_SC_prefill` uses `onceScatter` to fill the lower 40 mm of the
upper reservoir,
`(-0.020, 0.025, -0.020)`--`(0.020, 0.065, 0.020)` m, to `fieldValue 0.20`.
It never replenishes that region after time advances. `PH_SC_inlet` controls
the disjoint upper 10 mm slab,
`(-0.020, 0.065, -0.020)`--`(0.020, 0.075, 0.020)` m, and replenishes it during
the run with `repeatRandomPosition`. Together they target mean `lambda = 0.20`
over the two subregions of the 40 mm x 40 mm x 50 mm upper reservoir without
initializing the throat.

Because each candidate particle must fit inside its own insertion box, the
artificial split at `y = 0.065 m` can produce a transient particle-depleted
band around that plane. It is not a wall and closes after particles start
moving, but the initial vertical concentration profile should be checked and
the same loading protocol should be used in paired gravity and flow cases.

For the default exact `k = 3.00`, the geometric expectation is approximately
660 particles in the lower prefill and 165 in the top slab, or 825 total.
The actual count is controlled by the Eulerian `lambda` field and can differ
slightly with shape, orientation, and the final accepted particle. The
`nSolidsInDomain 100000` value is chosen as a non-binding constructor-time
ceiling over the intended, mesh-resolved `k` range; it is not a requested
particle count or a runtime replenishment limit.

Both models are monodisperse and use `uniformRandomRotation`. The patched
`onceScatter` model uses Shoemake quaternion sampling, matching the existing
uniform SO(3) proposal in `repeatRandomPosition`. Overlap rejection can still
bias the set of accepted orientations in a crowded region; the orientation
proposals themselves are uniform.

All particles start from rest (`velocity (0 0 0)` and `startSynced false`). A
coupled-flow comparison should use the same particle-velocity conditioning;
matching only the initial concentration does not make the startup states equal.

During time-zero initialization, 1000 consecutive rejected placements stop the
current model's constructor pass. During runtime inlet replenishment, 50
consecutive rejected placements stop only the current time-step pass; the
solver advances and may try the inlet again later.

## Workstation submission

The external `run_08_gravity_parametric.sh` script keeps machine-specific
paths outside the repository. It copies the selected template once to the
run-local `constant/triSurface/PH_SC.stl`, writes the requested gravity, and
writes the same requested scale to both logical addition models. It uses the
following run-directory format:

```text
<particle>_k<k>_g<g>_<SLURM_JOB_ID>
```

For example, `cube_k3.0_g20_14572831`.

## Build and run

This update changes `addModelOnceScatter.C` and the library-side
`initializeAddModels.H`. Rebuild the shared library after applying the overlay:

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
`g = (0 -2 0) m/s2`. It performs the one-time lower-reservoir prefill and then
maintains the top inlet slab. The external submission script overrides the
shape, size, and gravity in the run-local copy without changing the source
case.

The configured DEM substep is `stepDEM*deltaT = 1e-5 s`. Treat this as a
performance candidate rather than an already converged contact timestep. The
default 1.0 s run should be compared over a short interval with a `1e-6 s`
reference before interpreting discharge or clogging results.
