# OpenMP parallelization of the ionization sweep

_Working note, opened 2026-06-15, last measured 2026-09-07. Motivation:
profiling showed EXHALE used ~60 threads but only ~15-18% of their aggregate
capacity (htop), i.e. poor parallel efficiency._

**Every number in sections "Diagnosis", "Verification" and "Result" below was
measured on 2026-06-15 and describes the code of that date.** Two parallel
regions have been added since, and they change both the speedup and the
thread-count identity: the coupled source step, and the one heating assembly
(the photoionization heating and rate integrals are now assembled once per
step instead of once per solver pass). The measured state of the code as of
2026-09-07 is in "Current measured state" at the end; read that section for
what the code does now, and the dated sections for what it did then.

## Diagnosis

Before this change the only `!$omp parallel do` in the code was the
photoheating / photoionization-rate energy-integral loop in
`radiation/util_ion_eq.f90`. Everything else (the hydro RK stages, the
cell-by-cell ionization-equilibrium Newton solves, cooling, charge exchange,
temperature, advection) ran serial.

Measured thread scaling (no-metals HD209458b cold IC, steps/s):

| threads | 1 | 2 | 4 | 8 | 16 | 32 | 60 |
|---|---|---|---|---|---|---|---|
| baseline | 42.5 | 50.2 | 54.8 | 57.8 | **59.4** | 59.4 | 56.3 |

Only ~1.40x from 1->16 threads, nothing beyond 16, and **60 threads is slower
than 16** (fork/join overhead on a ~500-cell grid). Amdahl back-out => ~70%
serial. A lightweight phase profiler (env `EXHALE_PROFILE=1`, prints the
`ioniz_eq` wall-time fraction) showed `ioniz_eq` is **31.7%** of a no-metals
step but **60.9%** of a full-physics (He 2^3S + metals) step: the cell-by-cell
Newton solves dominate the heavy runs.

## What was changed

### (2) Thread-count default (`init.f90`)
`omp_get_max_threads()` capped at 60 grabbed every core for ~no gain. Now: an
explicit `OMP_NUM_THREADS` is honored as-is; otherwise the default is
`min(cores, 16)` (the measured knee). Frees cores and is no slower.

### (1) Parallel ionization-equilibrium cell sweep
The two solver loops over cells in `radiation/ionization_equilibrium.f90` are now
`!$omp parallel do` over cells. The sweep is embarrassingly parallel: at
`count > 0` each cell's Newton warm-start is its **own** previous-step value, so
cells are independent and the backward loop order is irrelevant. `count == 0`
runs serial (the `if(count > 0)` clause) because its first-step warm-start reads
the just-solved neighbor cell.

Thread-safety required making the scratch for each cell and the module state thread-local:

- **Global NL scratch** `sys_x, sys_sol, wa, info` (`init/parameters.f90`) made
  `!$omp threadprivate`. Serial regions (setup, post-processing) resolve to the
  master's copy = original behavior; each thread lazily allocates its own inside
  `ioniz_eq`.
- **Metal coefficients** `met_*` (`System_HeH_metals.f90`, also `use`d by
  `System_HeH_TR_metals`) made `threadprivate`: they are rebuilt per cell by
  `set_metal_coeffs`. Lazy allocation for each thread (the existing `if(.not.
  allocated)` guard runs per thread).
- **Charge exchange** `cx_kc` (rates in each cell) and `cx_metal_base` (the 4<->5
  row toggle) made `threadprivate`; `cx_set_cell` lazily allocates `cx_kc` per
  thread; `cx_metal_base` is broadcast with `copyin`. The setup-once
  `cx_act / cx_nact` stay shared (read-only during the sweep).
- **Newton usage counters** `nt_calls / nt_fallback` (`newton_solver.f90`)
  updated with `!$omp atomic`. The lazy `nt_init` block does not race because it
  is set on the serial `count == 0` step. The MINPACK routines (hybrd1, fdjac1,
  qrfac, ...) are pure (arg-based, no save state) and need no change.

The subroutine-local scratch (`params`, `usednt`, `i0`, `top`, `im`, the `meg_*`
metal-coefficient temporaries) is in the `private` clause; `schedule(dynamic,8)`
balances the uneven cost of solving each cell.

## Verification (correctness first)

_(2026-06-15 measurement; superseded, see "Current measured state".)_

Because each cell solve is independent and deterministic, a correct
parallelization must be **bit-identical** to serial. Verified with a
deterministic step cap (env `EXHALE_MAXSTEPS=N`): 3000-step HD209458b runs at
OMP_NUM_THREADS = 1 vs 8 vs 16 produced **byte-identical** `Hydro_ioniz.txt` and
`Ion_species.txt`, for **both** metals-off and full-physics (He 2^3S + metals).
Bit-identical per step => the full converged trajectory is identical, so the
existing production results are unaffected (no need to re-run them).

## Result

_(2026-06-15 measurement, steps/s; superseded, see "Current measured state".)_

| case | 1 thr | 16 thr | speedup |
|---|---|---|---|
| **full-physics** (He 2^3S + metals) | 27.1 | 56.2 | **2.07x** |
| no-metals | 41.8 | 60.2 | 1.44x |

The win is concentrated in the heavy full-physics runs (where `ioniz_eq` is
~61% of a step): they roughly **halve** in wall-time. No-metals is little
changed (its ionization solve is cheap; its serial bottleneck is the hydro /
radiation-column / temperature work, still serial: a candidate for a later
pass). 16 threads remains the practical knee.

## Diagnostics added (env-gated, zero cost when off)

- `EXHALE_PROFILE=1`: prints the `ioniz_eq` wall-time fraction every 500 steps.
- `EXHALE_MAXSTEPS=N`: stop after N steps and write output (serial-vs-parallel
  bit-comparison; also handy for fixed-length benchmarks).

## Current measured state (2026-09-07)

MEASURED by the COST3 item on the regression case `wasp_full` (He 2^3S +
metals), 300 steps:

| threads | wall time | speedup |
|---|---|---|
| 1 | 141.4 s | 1.0x |
| 16 | 17.25 s | **8.2x** |

Four times the 2.07x of 2026-06-15, because the work that was still serial
then has since been put in parallel regions of its own: the coupled source
step, and the one heating assembly.

**The thread count is no longer bit-neutral.** The same case at the two thread
counts differs by 6.6e-14 relative, i.e. at the last bits of a double, not at
any digit a result is read to. The mechanism is arithmetic grouping, not a
race: the parallel regions of `chemical_rate_coefficients` and `eval_cool`
(`radiation/util_ion_eq.f90`) partition the cells into `nblk = max(1,
min(2*nthr, ncell))` blocks with `nthr = omp_get_max_threads()`, so the block
boundaries themselves move with the thread count. Which of the routines called
inside a block turns a moved boundary into a moved last bit has not been
localized.

Two consequences for how the code is checked:

- A regression comparison across thread counts is a tolerance comparison, not
  a bitwise one. `backup/regression/run_check.sh` runs every case at
  `OMP_NUM_THREADS=1` for exactly this reason: one thread count is
  reproducible bit for bit, and that is what the goldens are snapshots of.
- A quantitative result quoted from a run should name the thread count it was
  produced at, as it already names the compiler.

## Thread-count determinism restored (2026-09-07, item THREAD-DET)

The 6.6e-14 relative difference between 1 and 16 threads recorded above was
localized and removed. Cause: `Cool_coeff.o` calls glibc libmvec's two-lane
`exp`, `log` and `pow` in the vectorized rate loops; the vector and scalar
variants differ in the last bit; `chemical_rate_coefficients` and `eval_cool`
cut the grid into `2 x omp_get_max_threads()` blocks, so the block starts,
hence the lane pairing, moved with the thread count. Both sweeps now use a
fixed even block length (`xuv_rate_block = 32`, `util_ion_eq.f90`). MEASURED on
`wasp_full` 300 steps: 1 and 16 threads byte-identical in every output file,
and the single-thread result unchanged from before the fix. The remaining
statement of this document, that a thread count cannot move a result, holds
again for these two sweeps; any new parallel region over cells must follow
the same rule (fixed block length, not a block count from the thread count).
