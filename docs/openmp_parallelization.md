# OpenMP parallelization of the ionization sweep

_Working note (2026-06-15). Motivation: profiling showed EXHALE used ~60 threads
but only ~15-18% of their aggregate capacity (htop), i.e. poor parallel
efficiency._

## Diagnosis

Before this change the only `!$omp parallel do` in the code was the
photoheating / photoionization-rate energy-integral loop in
`radiation/util_ion_eq.f90`. Everything else — the hydro RK stages, the
per-cell ionization-equilibrium Newton solves, cooling, charge exchange,
temperature, advection — ran serial.

Measured thread scaling (no-metals HD209458b cold IC, steps/s):

| threads | 1 | 2 | 4 | 8 | 16 | 32 | 60 |
|---|---|---|---|---|---|---|---|
| baseline | 42.5 | 50.2 | 54.8 | 57.8 | **59.4** | 59.4 | 56.3 |

Only ~1.40x from 1->16 threads, nothing beyond 16, and **60 threads is slower
than 16** (fork/join overhead on a ~500-cell grid). Amdahl back-out => ~70%
serial. A lightweight phase profiler (env `EXHALE_PROFILE=1`, prints the
`ioniz_eq` wall-time fraction) showed `ioniz_eq` is **31.7%** of a no-metals
step but **60.9%** of a full-physics (He 2^3S + metals) step — the per-cell
Newton solves dominate the heavy runs.

## What was changed

### (2) Thread-count default (`init.f90`)
`omp_get_max_threads()` capped at 60 grabbed every core for ~no gain. Now: an
explicit `OMP_NUM_THREADS` is honored as-is; otherwise the default is
`min(cores, 16)` (the measured knee). Frees cores and is no slower.

### (1) Parallel ionization-equilibrium cell sweep
The two per-cell solver loops in `radiation/ionization_equilibrium.f90` are now
`!$omp parallel do` over cells. The sweep is embarrassingly parallel: at
`count > 0` each cell's Newton warm-start is its **own** previous-step value, so
cells are independent and the backward loop order is irrelevant. `count == 0`
runs serial (the `if(count > 0)` clause) because its first-step warm-start reads
the just-solved neighbour cell.

Thread-safety required making the per-cell scratch and module state per-thread:

- **Global NL scratch** `sys_x, sys_sol, wa, info` (`init/parameters.f90`) made
  `!$omp threadprivate`. Serial regions (setup, post-processing) resolve to the
  master's copy = original behavior; each thread lazily allocates its own inside
  `ioniz_eq`.
- **Metal coefficients** `met_*` (`System_HeH_metals.f90`, also `use`d by
  `System_HeH_TR_metals`) made `threadprivate` — they are rebuilt per cell by
  `set_metal_coeffs`. Lazy per-thread allocation (the existing `if(.not.
  allocated)` guard runs per thread).
- **Charge exchange** `cx_kc` (per-cell rates) and `cx_metal_base` (the 4<->5
  row toggle) made `threadprivate`; `cx_set_cell` lazily allocates `cx_kc` per
  thread; `cx_metal_base` is broadcast with `copyin`. The setup-once
  `cx_act / cx_nact` stay shared (read-only during the sweep).
- **Newton usage counters** `nt_calls / nt_fallback` (`newton_solver.f90`)
  updated with `!$omp atomic`. The lazy `nt_init` block does not race because it
  is set on the serial `count == 0` step. The MINPACK routines (hybrd1, fdjac1,
  qrfac, ...) are pure (arg-based, no save state) and need no change.

The subroutine-local scratch (`params`, `usednt`, `i0`, `top`, `im`, the `meg_*`
metal-coefficient temporaries) is in the `private` clause; `schedule(dynamic,8)`
balances the uneven per-cell solve cost.

## Verification (correctness first)

Because each cell solve is independent and deterministic, a correct
parallelization must be **bit-identical** to serial. Verified with a
deterministic step cap (env `EXHALE_MAXSTEPS=N`): 3000-step HD209458b runs at
OMP_NUM_THREADS = 1 vs 8 vs 16 produced **byte-identical** `Hydro_ioniz.txt` and
`Ion_species.txt`, for **both** metals-off and full-physics (He 2^3S + metals).
Bit-identical per step => the full converged trajectory is identical, so the
existing production results are unaffected (no need to re-run them).

## Result

| case | 1 thr | 16 thr | speedup |
|---|---|---|---|
| **full-physics** (He 2^3S + metals) | 27.1 | 56.2 | **2.07x** |
| no-metals | 41.8 | 60.2 | 1.44x |

The win is concentrated in the heavy full-physics runs (where `ioniz_eq` is
~61% of a step): they roughly **halve** in wall-time. No-metals is little
changed (its ionization solve is cheap; its serial bottleneck is the hydro /
radiation-column / temperature work, still serial — a candidate for a later
pass). 16 threads remains the practical knee.

## Diagnostics added (env-gated, zero cost when off)

- `EXHALE_PROFILE=1` — prints the `ioniz_eq` wall-time fraction every 500 steps.
- `EXHALE_MAXSTEPS=N` — stop after N steps and write output (serial-vs-parallel
  bit-comparison; also handy for fixed-length benchmarks).
