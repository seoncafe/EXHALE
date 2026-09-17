# L24: the `atomic_elem_newton` reload at two entry pseudo-times

Item L24 of `docs/PLAN_20260913_lhs_stationary.md`, planned in
`docs/PLAN_20260916_rev3.md` section 4 (review row R20). Measured 2026-09-16
with the delivered tree binary `EXHALE.x`, md5
`c2e9c9990b9f14f1be8cd77abca68945` (verified before the runs). No source was
changed and nothing was built for this item.

## 1. What is being measured and why

`backup/regression/atomic_elem_newton` was pinned on 2026-09-09 for a
standstill of a different kind: the element-row coupled solve made no Krylov
direction at all, `||R||` stood at 1.892 to every printed digit for eleven
iterations and the run ended `info = 2`. On the delivered binary that
obstruction is gone: the partial run of 2026-09-16
(`backup/regression/atomic_elem_newton/run.log.20260916_partial`, the cold
start, 8 threads) descends 1.892 to 9.186e-02 with trust ratios near one and
then stops moving at iterations 60 and 61 with the full Krylov budget spent
(`gm` = 40) and the pseudo-time at 9.60e-05.

That shape is the pseudo-time ramp collapse item L4h is about
(`docs/lhs1140b_stationary_L4h_20260915.md`), and L4h's guard was measured on
the LHS 1140 b cases and on this reload, but the ramp's ENTRY was never varied
on this fixture. L14 (`docs/lhs1140b_stationary_L14_20260914.md` section 2.2
and its section 4.6 recipe) enters mapped seeds at `EXHALE_PTC_DTAU0=1e8` for
exactly this reason: where the dtau0 = 1 solve cannot leave the ramp floor,
the 1e8 entry takes the same state to the solve that works. L24 asks whether
that entry is what this reload needs.

The pseudo-time enters the step and not the residual, so the fixed point
`F(Y) = 0` is the same at either entry (L4h section 4.2); what the entry
changes is which states the iteration walks through and how far the ramp may
fall. READ from `steady_newton.f90` lines 15873-15881 and 16993: the floor of
every cut is `min(dt_CFL of the entry state, dtau0)`, so the two runs share
the same floor, and the ceiling is `1e14*dtau0`, which is the only bound the
entry moves.

## 2. Recipe

Both runs are of the reload, on scratch copies; the fixture directory itself
was not written to. From `backup/regression/atomic_elem_newton`:

- `input.inp` and `metals.inp` copied unchanged except for `Load IC? False` ->
  `Load IC? True` and the added line `Restart intent: stationary`;
- `IC/Hydro_ioniz_IC.txt` and `IC/Ion_species_IC.txt` copied into `output/` of
  the copy.

`Coupled carrier solve` is left at the fixture's own `True`. That is the
configuration the 2026-09-16 partial run was taken in and the one whose stall
this item re-measures; the partitioned recipe of the fixture README, which
sets that key to False, is a different measurement and is not this one.

Environment of both runs:

```
OMP_NUM_THREADS=1        the 8-thread path is not reproducible before item L15
EXHALE_JUDGED_ROWS=1     the linear rows divided by the certification scales (L4e)
EXHALE_JFNK_MAXIT=200    the cap of the run's one stationary solve
EXHALE_OUTER_PASSES=1    one outer pass
```

and, for run (B) only, `EXHALE_PTC_DTAU0=1e8`. Run (A) takes the default of
`Restart intent: stationary`, dtau0 = 1.0, READ from the log line
`(EXHALE_main) stationary restart: dtau0 = 1.00E+00` and its 1.00E+08 twin.

The state as loaded is the same in both runs and REFUSES on 11 entries of the
inventory (MEASURED, identical in the two logs): hydrodynamic mass 1.936E-01
at cell 285 against its own rounding anchor 3.1E-12, momentum 5.611E-02 at
cell 368, energy 1.674E+00 at cell 181, and the eight elemental transport rows
at cell 291 (He/H 5.840E-01, C 4.704E-01, O 1.468E-01, N 3.112E-01, Mg
5.034E-01, Ca 3.573E-01, Na 4.869E-01, Fe 2.954E-01).

## 3. Run (B), the 1e8 entry: it does not descend at all

MEASURED. The solve stops on the stagnation detector at iteration 64 of the
200 asked, "neither the judged distance nor the merit improved in 40
iterations, monotone search included", at `||R||` 1.903E+00, judged distance
6.209E+10 (its own best), `info = 2`, 2727.52 s for the pass. The lowest
`||R||` any iterate of the run carried is 1.634E+00 at iteration 9, against
1.674E+00 for the state as loaded: nothing was bought.

| it | `\|\|R\|\|` | `\|\|Fs\|\|2` | lam | dtau | gm | worst r | worst row | TR ratio | `\|\|s\|\|` | step |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.666E+00 | 1.68E+01 | 1.00E+00 | 2.00E+08 | 4 | 1.069 | energy of cell 181 | 1.059 | 1.410E+00 | accepted |
| 2 | 1.649E+00 | 1.68E+01 | 1.00E+00 | 4.00E+08 | 5 | 1.069 | energy of cell 181 | 0.994 | 2.821E+00 | accepted |
| 3 | 1.648E+00 | 1.69E+01 | 1.00E+00 | 8.00E+08 | 5 | 1.069 | energy of cell 181 | 1.022 | 3.578E-01 | accepted |
| 10 | 1.774E+00 | 1.69E+01 | 1.00E+00 | 1.02E+11 | 10 | 1.069 | energy of cell 181 | 0.999 | 1.430E+00 | accepted |
| 20 | 1.804E+00 | 1.70E+01 | 1.00E+00 | 1.05E+14 | 9 | 1.263 | mass of cell 290 | 1.000 | 5.278E-01 | accepted |
| 30 | 1.682E+00 | 1.71E+01 | 1.00E+00 | 2.62E+13 | 10 | 1.079 | mass of cell 191 | -20.959 | 7.200E-03 | rejected |
| 40 | 1.682E+00 | 1.70E+01 | 1.00E+00 | 1.60E+09 | 10 | 1.079 | mass of cell 191 | 1.000 | 5.493E-08 | accepted |
| 50 | 1.903E+00 | 1.69E+01 | 1.00E+00 | 3.20E+09 | 10 | 1.263 | mass of cell 290 | 1.000 | 1.099E-07 | accepted |
| 60 | 1.903E+00 | 1.68E+01 | 1.00E+00 | 8.00E+08 | 10 | 1.263 | mass of cell 290 | -24.774 | 2.197E-07 | rejected |
| 64 | 1.903E+00 | 1.68E+01 | 1.00E+00 | 1.60E+09 | 10 | 1.263 | mass of cell 290 | 1.000 | 5.493E-08 | accepted |

The events: the first rejection is at iteration 26 and from there the
rejections come in runs, ten of them "the reduction ratio is below eta" and
eleven "the ray difference is below its cancellation floor", the second being
the finite-difference probe of the Jacobian-vector product falling into
rounding because the radius is too small to form a difference. Steps accepted
44, rejected 21, model refused 18; the final trust radius is 2.747E-08 and the
accepted steps of the last twenty-five iterations have `||s||` of 5E-08 to
2E-07, which move nothing.

The state handed back refuses on all 11 entries, the same 11 the state as
loaded refuses on: mass 2.132E-01 at cell 285 against that cell's anchor
3.4E-12 (the largest mass measure is 1.903E+00 at cell 290, r = 1.2629),
momentum 5.540E-02 at cell 368, energy 1.601E+00 at cell 181, and the eight
elemental rows at cell 290 (He/H 5.849E-01, C 4.743E-01, O 1.454E-01, N
3.094E-01, Mg 5.085E-01, Ca 3.615E-01, Na 4.921E-01, Fe 3.001E-01). Flux
spread 3.526E+00.

What the trace says about the linear solve and the residual at this entry: all
65 Krylov cycles reached the requested tolerance, in 4 to 10 products of the
40 available, and the forcing term stayed at its 1.00E-01 ceiling because
`||F||` never fell. 5229 residual samples were admitted and NONE refused, in
any class (negative species unknown, element fraction above one, shared
element budget, non-finite sweep, non-finite row, no chemical root); no trial
or probe state was refused by an admissibility screen; there was no
out-of-domain closure activation and no cell without a chemical root. 23
trials were written onto the faces of the species box.

## 4. Run (A), the default entry: four and a half decades, then the element row

MEASURED. The solve descends `||R||` 1.674E+00 -> 1.795E-05 (its lowest, at
iteration 52) and stops on the stagnation detector at iteration 90 of the 200
asked, "neither the judged distance nor the merit improved in 40 iterations",
`||R||` 4.892E-05, judged distance 7.182E+03, `info = 2`, 5925.34 s for the
pass (its first 45 minutes shared the machine with run (B), so the seconds are
not a claim about cost).

| it | `\|\|R\|\|` | `\|\|Fs\|\|2` | lam | dtau | gm | worst r | worst row | TR ratio | `\|\|s\|\|` | step |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.678E+00 | 1.68E+01 | 1.00E+00 | 2.00E+00 | 3 | 1.069 | energy of cell 181 | 1.023 | 7.429E-02 | accepted |
| 2 | 1.688E+00 | 1.67E+01 | 1.00E+00 | 4.00E+00 | 3 | 1.069 | energy of cell 181 | 1.001 | 1.486E-01 | accepted |
| 3 | 1.860E+00 | 1.66E+01 | 1.00E+00 | 8.00E+00 | 2 | 1.069 | energy of cell 181 | 1.001 | 2.972E-01 | accepted |
| 10 | 1.022E+00 | 3.58E+00 | 1.00E+00 | 1.28E+02 | 4 | 1.001 | energy of cell 3 | 0.938 | 2.377E+00 | accepted |
| 20 | 3.939E-02 | 1.62E-01 | 1.00E+00 | 1.64E+04 | 1 | 1.006 | energy of cell 29 | 0.994 | 2.780E+00 | accepted |
| 30 | 7.105E-05 | 2.50E-03 | 1.00E+00 | 4.10E+03 | 40 | 1.447 | element He of cell 335 | 0.956 | 1.650E-03 | accepted |
| 40 | 2.699E-05 | 2.49E-03 | 1.00E+00 | 1.02E+03 | 40 | 1.447 | element He of cell 335 | -1.107 | 6.437E-04 | rejected |
| 50 | 2.082E-05 | 2.49E-03 | 1.00E+00 | 2.05E+03 | 40 | 1.447 | element He of cell 335 | 0.984 | 8.046E-05 | accepted |
| 60 | 4.230E-05 | 2.48E-03 | 1.00E+00 | 2.62E+05 | 40 | 1.447 | element He of cell 335 | 0.530 | 1.609E-04 | accepted |
| 70 | 4.894E-05 | 2.49E-03 | 1.00E+00 | 6.55E+04 | 40 | 1.447 | element He of cell 335 | 1.228 | 2.011E-05 | rejected |
| 80 | 4.892E-05 | 2.49E-03 | 1.00E+00 | 5.00E-01 | 40 | 1.447 | element He of cell 335 | -232.665 | 1.535E-10 | rejected |
| 90 | 1.334E-05 | 2.48E-03 | 1.00E+00 | 1.56E-02 | 40 | 1.447 | element He of cell 335 | 1.001 | 5.995E-13 | accepted |

The shape: `||R||` falls four and a half decades in the first thirty
iterations with trust ratios near one and the worst row in the base layer
(energy of cells 3 to 29), and at iteration 29 the worst row moves to the
helium element row of cell 335 (r = 1.447, a wind cell) and stays there for
the remaining sixty-two iterations while `||R||` moves only between
1.795E-05 and 4.9E-05. Steps accepted 58, rejected 33, model refused 19; the
rejections end on "the reduction ratio is below eta" 26 times, "the model and
the true slope differ in size" 6 times and "the ray difference is below its
cancellation floor" once; the final trust radius is 2.997E-13 and the last
accepted steps have `||s||` of 6E-13.

The state handed back refuses on 4 of the 11 entries the loaded state refused
on: hydrodynamic mass 2.155E-08 at cell 214 (r = 1.10501) against the fixed
3.0E-12, a judged distance of 7.182E+03 which is the whole of the run's judged
distance; momentum 2.995E-08 at cell 500 against 1.0E-08; energy 4.892E-05 at
cell 24 against 1.0E-06; and the He/H elemental transport row 4.031E-04 at
cell 335 against 1.0E-05. The seven other elemental rows are now inside their
tolerance, and the flux spread is 4.562E-08 against the 2.0E-05 gate. So the
default entry takes the reload from eleven refusing entries to four and from a
flux spread of 3.5E+00 to 4.6E-08, and does not certify it.

The linear solve is where this stall sits: 62 of the 91 Krylov cycles were
ENDED BY THE SUBSPACE, "the subspace was exhausted with a finite approximate
step", reaching a relative residual of 9.1E-01 to 9.99E-01 against the
1.00E-01 asked, and every one of those is an iteration in the plateau. 29
reached the requested tolerance, and they are the descent. Admissibility and
the residual branches are silent exactly as in run (B): 9112 residual samples
admitted and none refused in any class, no trial or probe refused by an
admissibility screen, no trial written onto a face of the species box, no
out-of-domain closure activation, no cell without a chemical root.

## 5. L4h's rule (4): it does not fire in either run

MEASURED: neither log carries the line "no full step in 8 consecutive
iterations; return ... to the best iterate" or its stopping twin, and in
neither run does the pseudo-time enter the band the rule acts in
(`dtau <= 10 dt_CFL`). Run (A) ranges 2.00E+00 to 8.39E+06 and ends
oscillating between 7.81E-03 and 3.12E-02; run (B) ranges 2.00E+08 to
6.71E+15 and ends at 1.60E+09.

That is a difference from the cold start this item was opened on. READ from
`backup/regression/atomic_elem_newton/run.log.20260916_partial`, whose first
stationary solve also entered at a pseudo-time of 1.0 (`dtau` 2.00E+00 at
iteration 1) but from the MARCHED state at `||R||` 2.419E-01: there the
pseudo-time does reach its explicit-stable floor, 9.60E-05, and rule (4) fires
once, "no full step in 8 consecutive iterations; return 1 to the best iterate
with dtau= 3.91E-04, ||R||= 9.186E-02", which is the 3.91E-04 that shows at
iteration 61.

So the ramp collapse rev3 section 4 cites belongs to the cold start's own
state and not to the reload: the same fixture, the same entry pseudo-time and
the same binary, reloaded from `IC/`, descends two and a half decades further
(4.892E-05 against 9.186E-02) and never reaches the band where rule (4) can
act.

## 6. The candidates rev3 names, read off the two traces

**The entry pseudo-time. EXCLUDED as the repair, and the direction is the
answer.** The two runs differ only in it, and 1e8 is the worse of the two by
four and a half decades of `||R||` and by seven of the eleven refusing
entries. L14's recipe, which enters mapped seeds at 1e8 where the dtau0 = 1
solve cannot leave the ramp floor, does not transfer to this reload for the
plain reason that this reload's dtau0 = 1 solve never reaches the ramp floor
(section 5). At 1e8 the shift `1/dtau` is 1e-8 of the band and the step is the
unshifted Newton step from a state four and a half decades away from the root,
which the trust region then refuses down to a radius of 2.7E-08.

**The preconditioner and the linear solve. SUPPORTED, in run (A) and only
there.** The 62 cycles that exhaust the 40-vector subspace at 0.91 to 0.999 of
their right-hand side are exactly the iterations of run (A)'s plateau, and the
29 that reach the 1.00E-01 forcing term are exactly its descent. In run (B)
all 65 cycles reach that tolerance in 4 to 10 products, so the linear
obstruction is a property of the states run (A) reaches and not of the
configuration: it appears when the solve has descended into the 1e-5 decade
and the binding row has become an element row. This is the obstruction items
N26c and N27 recorded ("the binding helium row of cell 246 is held by the
linear solve, Krylov 40 of 40 at 0.99"), two decades lower in `||R||` and at a
different cell.

**Admissibility. EXCLUDED in both runs.** 9112 and 5229 residual samples
admitted, none refused in any of the six classes, and no trial or probe state
refused by an admissibility screen. No step of either run was lost to the
describable set.

**A state-dependent residual branch. EXCLUDED in both runs.** No out-of-domain
closure activation, no cell without a chemical root, no non-finite sweep or
row, and the eliminated-species closure row is inside its 1.0E-06 in both
handed-back states (6.692E-12 at cell 231 in (A), 4.395E-12 at cell 248 in
(B)). No branch message appears at all.

**The coupling, and which row holds. In run (A) it is the element/hydrodynamic
pair at two different radii.** The worst row of the linear model from
iteration 29 on is the helium element row of cell 335 (r = 1.447), and the
certification refuses the He/H elemental transport row there at 4.031E-04
against 1.0E-05, while the judged distance 7.182E+03 is entirely the
hydrodynamic mass row of cell 214 (r = 1.105) at 2.155E-08 against the fixed
3.0E-12. The two are not the same cell, so the solve is being asked to hold a
wind element row and a base-layer mass row at once, and the linear solve that
would move both is the one exhausting its subspace. In run (B) the coupling
never gets a chance: the mass row of cell 290 holds from iteration 20 on at
1.9E+00.

## 7. What stays open

- **No entry defect is proved and none is refuted for the general case.** What
  is measured is this fixture: the 200 iterations rev3 asked for were not
  spent, both solves ended earlier on the stagnation detector (90 and 64), so
  the reading is of two stagnations and not of two caps.
- **Run (A) is not certified and the reason is not diagnosed here.** The
  linear solve exhausting its subspace is where the plateau sits; whether the
  band preconditioner, the 40-vector subspace or the element row's scaling is
  what puts it there is the question items N25, N26c and N27 left open, and
  this measurement adds only that it is still open two decades lower.
- **One outer pass only.** Both runs were given `EXHALE_OUTER_PASSES=1`, so
  the elemental relaxation ran once. The twelve-pass partitioned ladder of
  2026-09-11 CERTIFIES this fixture (fixture README, items P16 and P17), and
  nothing here contradicts that: it is a different route and a different
  budget. What this item measured is the coupled single pass at two entries.
- **The 2026-09-16 cold-start stall is a different state.** Section 5 shows
  the reload does not reproduce it. Whether the cold start's own collapse is
  removed by any entry was not measured, because a cold start is not a reload
  and `EXHALE_PTC_DTAU0` does not reach the marching hand-off.
- The fixture `backup/regression/atomic_elem_newton` is UNCHANGED apart from
  one appended paragraph; its `IC/`, `output/`, `run.log` and
  `run.log.20260916_partial` are the records items N7, N26b and L15 quote.

## 8. Where the runs are

Both run directories are kept under the item's scratch directory,
`<scratchpad>/L24/A` (the default entry) and `<scratchpad>/L24/B` (1e8), each
with its `input.inp`, `metals.inp`, `run.log` and `output/`. The variant
fixture `backup/regression/atomic_elem_newton_dtau1e8/` carries only a
`README.md` and a `run.sh` that reproduces run (B) on a copy; it holds no
initial condition and no output of its own.
