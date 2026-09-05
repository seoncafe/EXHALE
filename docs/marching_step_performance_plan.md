# Where a marching step spends its time, and the plan to make it scale (2026-09-04)

**Judgement first. The hydrodynamic update is not the bottleneck, and
parallelizing it would gain 3 percent. The step is limited by the cooling
evaluation `eval_cool`, called four times per step, and 80 percent of that
evaluation is spent on metal ions that this run does not contain.** The plan
below is therefore (0) stop computing absent species, (1) make `eval_cool`
run in parallel over cells, (2) remove one of its four calls, and only then
(3) touch the hydro stages.

## 1. Measurement

Hot-Uranus molecular run (`Molecular chemistry: True`, carrier transport on,
metals off: every metal column identically zero), 300 marching steps,
`EXHALE_PROFILE=1`, the phase timers of `EXHALE_main.f90` plus the five
sub-block timers added to `eval_cool` (`util_ion_eq.f90`). Times are wall
seconds over the 300 steps.

| phase | 1 thread | 16 threads | scales? |
|---|---|---|---|
| ionization sweep `ioniz_eq` (includes one `eval_cool` call) | 13.2 | 2.5 | yes, 5.2x |
| semi-implicit energy update (three `eval_cool` calls) | 4.5 | 4.3 | **no** |
| carrier transport (block-tridiagonal Newton) | 0.7 | 0.7 | no |
| hydro: three SSP-RK3 stages (reconstruction, Riemann, positivity, BC) | 0.6 | 0.6 | no |
| everything else | 0.1 | 0.1 | -- |
| **whole step** | **19.0** | **8.2** | 2.3x |

At 16 threads the step is 27 ms, of which 14 ms is the energy update, 8 ms
the sweep, 2 ms the carrier Newton and 2 ms the hydro. Amdahl on the present
code: the serial 5.7 s caps the speed-up at 3.3x however many threads are
added -- which is why 8, 16, 32 and 48 threads all give 2000-2400 steps per
minute and each core reads 20 percent busy.

`eval_cool` itself, summed over all four calls per step (all callers):

| sub-block | 1 thread | share |
|---|---|---|
| metal recombination / ionization / cooling table interpolation (`n_mion` ions) | 3.71 s | **61 %** |
| fine-structure line transfer + CNO cooling | 1.09 s | 18 % |
| H, He recombination / ionization / excitation coefficient fits | 1.03 s | 17 % |
| free-free Gaunt factors | 0.16 s | 3 % |
| H3+ and molecular infrared | 0.03 s | 1 % |
| **total** | **6.03 s** (20 ms per step, 5 ms per call) | |

**In a run with no metals, 79 percent of the cooling evaluation is
interpolating tables for ions whose density is exactly zero in every cell.**

## 2. The plan, in the order of return on effort

### Step 0 -- skip the metal blocks when there are no metals (serial, no threads)

Gate the metal table loops, the fine-structure transfer and the CNO cooling
in `eval_cool` (and the matching rate blocks in the ionization sweep) on the
metals actually being present (`thereis_metals`, or `maxval(nm) > 0`), with
the metal outputs (`rec_m`, `aion_m`, `c_metal`) set to zero so no caller reads
an undefined array. For a metals-off run the skipped terms are products with
zero and the result is bitwise identical; for a metals-on run nothing changes.

Expected: `eval_cool` 20 -> 4 ms per step; the 16-thread step 27 -> 11 ms,
i.e. **2.4x faster marching on every H/He run, with no parallel code at
all.** This is the first thing to do and the only one that is free.

### Step 1 -- `eval_cool` in parallel over cells

Every quantity `eval_cool` computes is a function of the cell's own
(T, densities); there is no coupling between cells and no reduction, so the
per-cell arithmetic is identical in any thread order and the result is
bitwise reproducible -- the same property the ionization sweep already relies
on. Two ways to do it; the first is recommended:

* **(1a) Block decomposition at the top.** Turn `eval_cool` into a driver
  that splits `1-Ng:N+Ng` into contiguous blocks and, inside one
  `!$omp parallel do`, calls a kernel `eval_cool_block(j_lo, j_hi, ...)`
  containing the present body written over `j_lo:j_hi`. The ~30 coefficient
  routines it calls are explicit-shape array routines over the whole grid;
  they either become elemental / take an index range, or the kernel calls
  them on its block through contiguous array sections of the same declared
  extent. This is the change with the best fork/join ratio: one parallel
  region per call.
* **(1b) Parallelize each coefficient routine's own loop.** About 30 small
  `!$omp parallel do` regions per call; simpler edits, worse overhead
  (~30 x 5 us against 5 ms of work is still acceptable at 16 threads, but
  it does not improve with more threads).

Caveat measured on this code before (`Update_EXHALE.md` section 145 and the
Codex review item 4.1): changing explicit-shape dummies to assumed-shape
moved the `-O3` arithmetic and broke byte-identity. The kernel must keep
explicit-shape arrays (pass `j_lo, j_hi` and the full arrays, loop over the
range), so the compiled arithmetic per cell is unchanged.

Expected at 16 threads: the three energy-update calls 14 -> ~2 ms and the
sweep's own call 5 -> ~0.6 ms; with step 0 already in, the step goes
11 -> ~6 ms. **Combined with step 0, about 4.5x over today's marching.**

### Step 2 -- three `eval_cool` calls per step instead of four: REFUTED (2026-09-05)

The premise was wrong. The sweep's one `eval_cool` (`ioniz_eq`, line 849)
runs BEFORE the cell sweep, on the entry composition, because the sweep needs
the recombination and ionization coefficients it returns; the sweep then
overwrites `nhi`, `nhii`, the molecular arrays and `rho` itself
(`calc_rho`), and the energy update builds its species from that post-sweep
`(rho, f_sp)` and evaluates at a `T_old` that has been round-tripped through
`p` against the post-sweep particle count. The two evaluations are of
different states by construction, and passing the sweep's `cool` in would
have replaced the new composition's cooling by the old one's -- a lagged term,
not a removed duplicate. What is genuinely wasted is smaller: in the
semi-implicit energy path (the default) the cooling SUM the sweep's call
returns is never read (`solve_energy_semi_implicit` overwrites `cool` with
its own `cool_trial`), while the coefficients that same call produces are
needed. Skipping only the sum would save under 2 per cent of a step (the
whole energy phase is 4.8 per cent, the sweep 84.5 per cent, measured
2026-09-05 at 1 thread).

What the reading exposed instead was a defect, fixed the same day
(`Update_EXHALE.md` section 170): the heating and cooling `ioniz_eq` returned
were those of the composition BEFORE its sweep. With the assembly moved after
the sweep, the sweep's returned `cool` IS the cooling of the state the energy
update starts from, and the energy update's first `eval_cool` was removed on
that basis -- the count is four evaluations per step as before (coefficients
before the sweep, the sum after it, the derivative and the final one in the
energy update), the energy phase fell 31 per cent, and the whole step rose
5 per cent for the second traversal. Step 2 as a speed-up is closed.

### Step 3 -- the hydro stages (last, small)

Reconstruction, Riemann fluxes and the positivity limiter are face- and
cell-local and parallelize the same way (no reductions except the CFL
minimum in `eval_dt`, which stays serial or uses a reduction that is then
fixed-order). At 2 ms per step out of 27 the gain is 3 percent today and
~30 percent only after steps 0-2 have removed everything else. Do it then, or
not at all.

### What stays serial

The carrier Newton (block Thomas, 2 ms), `Apply_BC`, and the CFL minimum. The
Jacobian assembly inside the carrier Newton (per-cell finite differences of
`carrier_source`) is cell-local and could join step 1's pattern if it ever
matters.

## 3. Verification that goes with each step

* Step 0: `make check` must be 11 of 11 byte-identical (five of the matrix
  cases are metals-off; the metals-on cases are the control that nothing
  physical moved). `run_fcheck.sh` for the zero-initialised outputs.
* Step 1: byte-identity serial vs parallel at 1, 8, 16 threads on the same
  case (the ionization sweep's existing test, applied to `eval_cool`), plus
  `make check`.
* Step 2: refuted by reading the code before any assertion was run (see
  above); nothing to verify.
* Every step: the `EXHALE_PROFILE=1` phase table before and after, on the
  same 300-step case, so the gain is a measured number and not a projection.

## 4. What this does not address

The number of steps. A marched hot-Uranus layer needs 1e6-1e7 s of model
time to relax thermally at ~2 s per step (`Update_EXHALE.md` section 163.5);
a 4-5x faster step turns days into a day, and no per-step speed-up turns it
into an hour. That needs an implicit treatment of the slow rows or a
steady solver that reaches the layer, which is a different plan.
