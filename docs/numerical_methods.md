# Cleaner convergence methods for EXHALE steady-state winds

**Status.** This is the options survey written before a steady solver existed.
**Option 3 (PTC -> JFNK) was chosen and implemented** and is now the
`Solver: Newton` key: matrix-free Newton-Krylov with a banded preconditioner,
pseudo-transient continuation, frozen WENO weights and a non-monotone line
search (`src/modules/time_step/steady_newton.f90`; design record
`docs/steady_solver_design.md`, `docs/newton_scaling_and_base_wall.md`). It does
*not* replace the two-stage marching described below: it warms up with it and
takes over at `du < newton_du_switch`. Option 1 (local time-stepping) and
Option 2 (steady BVP relaxation) were not implemented in the main solver;
Option 2's algorithm does exist in the tree as the Wind-AE initial-condition
generator (`src/modules/wind_ae/`, `IC mode: windae`). The rest of this note is
kept as the reasoning that led to that choice.

## Motivation

EXHALE reaches its steady state by **explicit time-marching** (RK + HLLC/ROE
finite volume) until the momentum non-uniformity `du = dMdot/Mdot` (the spatial
spread of the mass flux `rho*v*r^2` over the wind region `[j_min:N]`) drops below
`du_th`. In practice this requires the somewhat awkward **two-stage workflow**:
run **PLM** from general initial conditions until `du` falls to about 0.5-1,
stop, then **restart** with `Load IC` + **WENO3** and run to full convergence
(`du < 1e-3`).

The clunkiness has two distinct roots:

1. **Explicit time-marching to steady state is inherently slow.** Reaching the
   steady state takes of order (domain sound-crossing time)/(CFL `dt`) steps.
   This is worst for the **weak, subsonic winds** (e.g. HD189733b, HD209458b),
   which can need ~50k-100k steps; strong/fast winds (WASP-121b) converge ~25x
   faster.
2. **PLM is diffusive**, so it stalls around `du ~ 0.02` and one must hand off
   to the higher-order **WENO3** to reach `du < 1e-3`: the source of the
   two-stage restart.

Below are cleaner alternatives, in increasing order of implementation effort.

## Option 1: Local (cell-by-cell) time-stepping  *(cheapest; stays inside EXHALE)*

For a *steady* state, time-accuracy is irrelevant, so give each cell its own
maximum-CFL `dt` instead of the global minimum `dt`. The slow-relaxing cells
then advance much faster, accelerating convergence to steady state, precisely
where it is needed (the weak winds), and largely removing the need for staged
restarts.

- **Pros:** minimal code change (a few lines in `eval_dt` and the update loop);
  keeps the existing reconstruction and flux routines; directly attacks the
  slow-weak-wind problem. Optionally combine with implicit residual smoothing or
  FAS multigrid for further acceleration.
- **Cons:** still only first-order (linear) convergence to the steady state, a
  speedup, not the quadratic convergence of a Newton method.

## Option 2: Steady transonic-wind BVP relaxation  *(physics-native; what wind-AE does)*

Drop time-marching entirely. Discretize the **steady** wind ODEs and Newton-
relax them (e.g. Numerical Recipes `solvde`), enforcing the **sonic-point
regularity (critical-point) condition** so the solution passes smoothly through
Mach 1. This yields quadratic convergence, `du -> 0` cleanly, and a single method
with no reconstruction switch.

- **Pros:** fast (quadratic near the solution); gives the exact steady state to
  machine precision; conceptually clean; this is the standard wind-solver
  approach (wind-AE).
- **Cons:** needs a **near-transonic initial guess** and careful critical-point
  handling; less forgiving when the thermal structure is stiff or multivalued
  (strong metal-cooling fronts can make the steady solution hard to bracket).

## Option 3: Pseudo-transient continuation -> Newton-Krylov (JFNK)  *(gold standard)*

Solve the steady residual `F(U) = 0` directly with **matrix-free Newton-Krylov
(JFNK / GMRES)**, globalized by **pseudo-transient continuation (PTC)**: start
with a small pseudo-`dt` (robust, behaves like the current time-marcher) and
ramp `dt -> infinity`, which turns the iteration into pure Newton (quadratic).

- **Pros:** robust startup *and* fast finish in one method; converges to a
  machine-precision steady state; no PLM/WENO3 hand-off; the modern default for
  steady-state CFD.
- **Cons:** requires implementing the residual evaluation and a simple
  preconditioner (and the PTC `dt` schedule).

## Recommendation

| Want | Use |
| :-- | :-- |
| A quick, low-risk win inside EXHALE | **Option 1 (local time-stepping)** (+ residual smoothing) |
| The cleanest "always converges well" solver | **Option 3 (PTC -> JFNK)** |
| The most physics-tailored / fastest | **Option 2 (BVP / Newton relaxation)** |

The fundamental tradeoff is **robustness-of-startup** (time-marching, PTC) vs.
**speed-near-the-solution** (Newton, BVP). Pseudo-transient continuation is
attractive because it provides both: time-marching-like robustness far from the
solution that continuously morphs into Newton's quadratic convergence as the
residual falls.

A pragmatic path: start with **local time-stepping** (immediate speedup, keeps
the code), and if a fully clean single-method solver is desired later, move to
**PTC -> JFNK**.
