# Velocity-controlled CFD-DEM clogging at k = 3.0

This case extends `08_k3.0_gravity` to a fully coupled, gravity-free water
flow.  It keeps the same particle and contraction geometry and prescribes a
downward reservoir inlet velocity of 5 mm/s.

## Case definition

- Upstream reservoir: 40 mm x 40 mm, y = 25--75 mm
- Downstream throat: 10 mm x 10 mm, y = 0--25 mm
- Default sphere radius: 1.66667 mm
- Default sphere volume: 1.93927e-8 m3
- Opening ratio: k = 10/3.33334 = 3.00
- Eulerian mesh spacing: 1 mm (82,500 cells)
- Water: rho = 1000 kg/m3, nu = 1e-6 m2/s
- Gravity: zero, so this case isolates liquid driving
- Inlet velocity: `(0 -0.005 0)` m/s
- Fixed outlet gauge pressure: zero
- Initial and inlet target solid volume fraction: 0.30

The inlet area is 1.6e-3 m2, so the imposed bulk flow is 8e-6 m3/s
(8 mL/s).  The clean-flow mean throat velocity is 0.08 m/s.  With no mean
particle slip and no clogging, the nominal particle flux is

```text
Ndot = phi A U / Vp = 123.76 particles/s
     = 7,425 particles/min = 445,529 particles/hour.
```

This is a reference flux, not a prescribed particle count.  The measured exit
rate can differ because of slip, accumulation and clogging.

## Initial fill and inlet replenishment

Two existing body-addition models are used; the placement algorithms remain
inside those models and do not introduce a separate placement class:

1. `PH_SC_prefill` uses `onceScatter + geometricVolumeBased` in the complete
   upstream reservoir, y = 25--75 mm.  It converts phi = 0.30 to the nearest
   integer count and constructs a complete, full-size, continuous-random pose
   plan before registering any immersed body.  For the supplied sphere this is
   1,238 particles and a realizable fraction of 0.300101.
2. `PH_SC_inlet` uses `repeatRandomPosition + fieldBased` to maintain
   `lambda = 0.30` in the upper 10 mm control slab, y = 65--75 mm, after t = 0.
   New particles are sampled at their final size and at continuous positions.

`nSolidsInDomain 2000` is only a constructor-time safety ceiling.  It is not a
request to insert 2,000 bodies.  Initial filling stops only after its exact
nearest-count plan has been registered; runtime inlet replenishment is not
limited by this constructor ceiling.

The prefill planner uses broad-phase bins, AABB rejection and an exact sphere or
convex-body overlap test.  If direct random placement stagnates, it removes a
small conflict neighbourhood only in the in-memory pose plan and atomically
rebuilds that neighbourhood.  Failed rebuilds roll back without touching body
IDs, contact history, Verlet lists or solver fields.  There is no lattice
fallback, temporary particle shrinking, particle growth or underfilled exit.

The inlet rebuilds its broad-phase index from the moving bodies, checks a
bounded batch of full-size candidates, and retains the framework contact test
as the final wall/mixed-geometry guard.  Exhausting the per-time-step proposal
budget defers insertion to the next time step; it never reduces particle size.
`nonContainingFaces (yMin)` treats the slab's lower face as an internal control
surface: the particle centre remains at or above y = 65 mm while its support
may extend below it.  The other faces contain the complete rotated particle.

Both sources use `uniformRandomRotation` and `startSynced true`, so orientations
are isotropic and each accepted particle starts with its local fluid velocity.
Explicit dictionary seeds make placement reproducible and MPI-consistent.

## Time integration and startup

- Maximum CFD step: 1e-3 s (`maxCo = 0.7`)
- Maximum DEM substep: `stepDEM*deltaT = 1e-5 s`
- Field write interval: 0.05 s
- Pressure-drop sample interval: every 10 CFD steps (nominally 0.01 s)
- Default end time: 30 s

`Allrun` first runs `potentialFoam` on the empty static geometry.  This gives a
divergence-free clean-flow initial field before the particles are constructed,
instead of starting the 16:1 area-ratio contraction from a uniform velocity.
OpenFOAM's automatic step adjustment monitors the fluid Courant number, not a
particle translation/rotation Courant number; monitor both in higher-velocity
runs.  The DEM step should still be compared with a 1e-6 s short reference
before quantitative production use.

## Recorded data

The existing particle-exit logger writes

```text
postProcessing/particleExit.dat
# time bodyId x y z
```

An event occurs when less than 1% of a particle's reference mass remains in the
domain.  It is therefore suitable for cumulative discharge and avalanche
statistics, but it is not an exact fixed-plane crossing time.

Because `bottomOutlet` has fixed gauge pressure `p = 0`, the area-average
pressure on `topOpen` is exactly the inlet-to-outlet pressure drop.  The
`pressureDrop` function object converts kinematic pressure to Pa and writes
only two numeric columns:

```text
postProcessing/pressureDrop/<runStart>/surfaceFieldValue.dat
time deltaP_Pa
```

A perfectly sealed throat is singular under incompressible fixed-flow control:
the pressure demand can grow until the nonlinear solve fails.  Treat sustained
pressure growth together with a censored no-exit interval as the permanent-clog
signal; the last converged pressure timestamp is the valid observation end.

## Cumulative discharge and avalanches

After a run (or an interruption), generate cumulative discharge with:

```sh
python3 binParticleExit.py
```

It writes `postProcessing/particleExit_cumulative.txt` with two columns,
`time cumulative_discharged_particles`.  This is the unchanged script from
`08_k3.0_gravity`, so its endpoint remains the bin containing the last exit.

Avalanche size alone cannot distinguish a short pause from a long-lived clog.
For exit times `t_i`, define `dt_i = t_(i+1) - t_i`; a gap larger than a chosen
`dt_c` separates two avalanches.  Analyze both the avalanche size and the
adjacent arrest duration:

```sh
python3 analyzeAvalanches.py DT_C_SECONDS
```

The script writes:

- `postProcessing/interExitTimes.txt`: `exit_time inter_exit_time`
- `postProcessing/avalanches.txt`: avalanche ID, start/end time, size,
  preceding arrest, following arrest, and a right-censor flag

The final no-exit interval is right-censored: it is the observed lower bound on
the arrest time, not proof of a permanent clog.  As an initial scale only,
`d/U_throat = 0.003333/0.08 = 0.0417 s`; production analysis should determine
`dt_c` consistently for every shape, for example from the lower edge of the
power-law tail of the inter-exit-time survival distribution.  Never tune the
threshold separately to favor one particle shape.

## Build and run

The overlay changes the existing `onceScatter` and `repeatRandomPosition`
headers and implementations, plus this case dictionary and README.  It adds no
new source file and does not change `Make/files`.  Rebuild the shared library
from the repository root; the solver does not need to be rebuilt:

```sh
(cd src/HFDIBDEM && wmake libso)
```

Then run from this case directory:

```sh
./Allrun
python3 binParticleExit.py
python3 analyzeAvalanches.py 0.0417
```

Use `./Allclean` before an independent realization.  The supplied analysis
workflow assumes one clean realization, matching the batch script that removes
old time directories, `bodiesInfo` and `postProcessing` before every run.
