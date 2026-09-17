# L4d: the two step-control changes of L4b, implemented and measured

Item L4d of `docs/PLAN_20260913_lhs_stationary.md`, with the corrections
required by `docs/PLAN_20260913_lhs_stationary_review.md` section 1. Every
number below is MEASURED on this tree unless it is marked READ or DERIVED.
The file changed is `src/modules/time_step/steady_newton.f90` and nothing
else.

## 1. Verdict

**The stagnation detector was the refusal L4b measured, and repairing it
lifts that refusal.** On the reproduction the entry text stops the
stationary solve as STAGNATED after 40 iterations and hands back the state
it was given; with the counter reading the merit as well as the judged
distance, the same solve from the same state reaches `info = 0` at iteration
98, `||R||` 1.121e-1 -> 2.799e-8, with every hydrodynamic row inside its
tolerance. On the full outer loop of the isolating control case of L4b,
`atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`, passes 2 to 6 -- where the entry
text was `info = 2` in every configuration L4b tried and the loop REFUSED at
pass 4 -- are now `info = 0` with `omega` never halved, and the elemental
transport row reaches its 1e-5 tolerance at pass 7 (6.37e-6). This is not a
certified solution: what the loop then stands on is the hydrodynamic energy
row of the OUTERMOST cell, a different entry from the one L4b diagnosed and
one this item did not investigate (section 10).

**The forcing term is NOT that repair, and L4b's reasoning for it is refuted
by measurement -- but it is not idle either.** L4b inferred that the fixed
Krylov tolerance `1.0d-1` let the linear solve leave ten percent of the
scaled residual in the continuity rows. It does not: with the banded
preconditioner this route carries, the cycle reaches 7.9e-4 to 1.8e-3 of its
right-hand side in ONE product, two decades inside the tolerance it is asked
for, at every iteration of the STALLED pass. On that pass the forcing term
stays at its ceiling, because the merit barely moves, and it changes
nothing. Where the merit DOES fall it acts: on outer pass 1 of the same case
it leaves the ceiling at 11 of 17 cycles, down to 1.83e-3, and the pass ends
at `||R||` 2.720e-8 against the control's 1.692e-7, a factor 6.2, in the
same 18 outer iterations. So the forcing term is a tighter linear solve
where a tighter one is affordable, and it is not the answer to the
continuity excursion; the row-resolved measurement of section 5 says why
no scalar tolerance can be.

**What the continuity excursion actually is.** MEASURED, row by row, at the
first Newton step of the reproduction: the linear residual the cycle leaves,
divided by the same row of `F`, is 2.6e+05 in the mass row, 2.0e+05 in the
momentum row and 0.44 in the energy row. The 2-norm the cycle minimizes is
the energy rows; the mass and momentum rows of the entry state stand at
5.0e-9 and 4.2e-13, ten decades below, so a linear residual that is
negligible in the norm is five decades larger than the rows it is left in.
The step then raises the mass row to 1.17e-2 and the momentum row to
3.1e-3. This is not "the Krylov tolerance was not reached" and not "the
nonlinear model failed": it is a scalar tolerance on a norm that does not
see the rows the certification reads.

## 2. What changed in the source

Both changes are in `solve_steady_jfnk`.

**(1) The stagnation counter reads both functionals the solve carries.** The
ledger that chooses the returned state still ranks iterates on the judged
distance alone, and the certification is untouched. What changed is the
counter under it: it is reset when the judged distance improves, as before,
OR when the MERIT improves. The merit it reads is not the one the line
search descends on. Three properties were required by the review and all
three are in the text:

- **a fixed scale.** The line search's merit is `|| F/D ||_2` with `D`
  rebuilt from every iterate, so two of its values are two different norms.
  MEASURED on the reproduction: the entry state reports 3.18e-2 at the top
  of the solve and the same state reports 1.97e-2 after a rebuild of the
  scales, a factor 1.6 with the state untouched. `Dfix` is the column scale
  of the state the solve was entered at, held for the life of the solve, and
  the counter reads `|| F/Dfix ||_2`.
- **a window of the counter's own length.** The reference is the
  fixed-scale merit of `n_stall_best` iterations ago. Comparing against the
  running minimum of the window instead makes the test a fall per
  ITERATION, because a descending merit's minimum IS its previous value, and
  the length of the window then never enters; that was measured and
  discarded (section 3).
- **a margin, not any strict decrease.** The current value must stand below
  that reference by more than `merit_fall_min`, default 0.1, overridable for
  measurement by `EXHALE_MERIT_FALL_MIN`.

**(2) A forcing term in place of the literal `1.0d-1`.** The relative
tolerance of the two `pgmres` calls of the three-unknown branch is now
Eisenstat & Walker (1996), SIAM J. Sci. Comput. 17, 16 (doi
10.1137/0917003, ADS `1996SJSC...17...16E`; bibliographic record verified
through ADS, the article text not obtained), their choice 2,
`eta_k = gamma (||F_k||/||F_{k-1}||)^alpha` with `gamma = 0.9` and
`alpha = (1+sqrt(5))/2`, floored at 1e-4 and capped at the 1e-1 the branch
carried before. The safeguard `eta_k <- max(eta_k, gamma eta_{k-1}^alpha)`
is applied WITHOUT the paper's `> 0.1` gate, and that is the one departure
from the published form: under a ceiling of 1e-1 that quantity is at most
`0.9 * 0.1^1.618 = 0.0217` (DERIVED, and the review's own derivation), so
the gate would make the safeguard unreachable at every iterate. The comment
in the source says this.

**(3) A row-resolved diagnostic, off by default** (`EXHALE_LINEAR_ROWS=1`):
the unpreconditioned linear residual of the PTC system row by row against
the same rows of `F`, the three hydrodynamic rows entering and leaving the
step, and the judged slots of the state and of the trial. It costs one
Jacobian action an iteration and nothing when it is not armed. It is what
section 5 is measured with.

## 3. The window and the margin, and what anchors them

**The residual evaluation carries no run-to-run uncertainty on this case.**
MEASURED: three runs of the reproduction at 8 threads and one each at 4 and
16 threads, each stopped after two Newton iterations, wrote
`output/Hydro_ioniz.txt` identical in every column to every digit (maximum
relative difference 0.0 over all 504 x 7 entries, all five pairs). So the
margin cannot be anchored on scatter between runs; the only noise floor of
the quantity is the rounding of the assembly itself, near `sqrt(3N)` times
the double-precision epsilon, 4e-15 at N = 500 (DERIVED). The margin is
therefore not a noise threshold at all, and the two regimes below are what
sets it.

**A first form of the window was measured and discarded.** Taking the
reference to be the MINIMUM over the last `n_stall_best` fixed-scale merits
makes the test a fall per iteration whenever the merit descends, because the
minimum of a descending sequence is its last entry. MEASURED: with that
form at a margin of 0.1 the reproduction STAGNATES again at iteration 39 and
returns the entry state, because its merit falls about 1.8 percent an
iteration and never 10 percent in one. The delivered form compares with the
entry `n_stall_best` iterations back, which is the window the counter's own
length asks about.

**The margin separates two measured regimes of this route.**

| regime | fall of the fixed-scale merit over a 20-iteration window | where |
|---|---|---|
| a solve that goes on to certify | 25 to 35 percent | the reproduction, iterations 1 to 98 |
| a grind four decades from certification | about 6 percent | outer pass 7 of the same case's loop, iterations 136 to 215, `||R||` 2.102e-3 -> 1.879e-3 |

The grind is what the first margin tried, 1e-3, let run: MEASURED, outer
pass 7 reached iteration 250 at `||R||` 1.810e-3 still falling 4 percent per
25 iterations, which extrapolates past the 3000-iteration cap. 0.1 stands
between the two and is what the default carries. The reproduction reaches
`info = 0` at iteration 98 and `||R||` 2.799e-08 at that default, the same
numbers the first form gave.

## 4. RED and GREEN

The reproduction is the isolating control of L4b section 5 (e),
`LHS1140b/models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`: the certified
scalar family's own configuration with the C, N, O reservoirs added, no
profile, He/H = 2.13. The state loaded is the one that case's campaign run
handed back, which is what every outer pass after the first is given
(L4b section 5 f measures the same structure). `EXHALE_OUTER_PASSES=1`
isolates one hydrodynamic solve. 8 threads, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_JFNK_MAXIT=300`.

Entry state: `||R||` 1.121e-1, judged distance 1.121e+05 held by the energy
row at cell 245 (r = 1.3278); mass 5.002e-9, momentum 4.183e-13, both
inside.

| build | outcome | iterations | `\|\|R\|\|` returned | mass | momentum | energy |
|---|---|---|---|---|---|---|
| entry text (RED) | **STAGNATED**, entry state returned | 40 of 300 | 1.121e-1 | 5.00e-9 | 4.18e-13 | 1.121e-1 |
| change (1) only | `info = 0` | 98 | 3.143e-8 | 9.31e-10 | 1.60e-14 | 3.14e-8 |
| changes (1) and (2), as delivered | `info = 0` | 98 | 2.799e-8 | 1.23e-9 | 4.65e-14 | 2.80e-8 |

Under the entry text the judged distance never improves on its entry value
(1.121e+05 at entry, 1.112e+09 at iteration 1 and never below the entry
again), the counter fires at 20, the solve restarts from the best iterate
and fires again at 40 -- while the merit falls from 3.18e-2 to 1.68e-2 and
is still falling. That is the RED.

**The full outer loop.** The same case from its own seed, default outer pass
budget, 8 threads, `EXHALE_PTC_DTAU0=1.0`.

Both builds were run here; the control's ladder reproduces the campaign's
own run of this case (`LHS1140b/models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13/run.log`)
in every printed digit of every pass, which is what says the control build is
the campaign's binary.

| pass | entry text: hydro | element row (cell) | energy | omega | as delivered: hydro | element row (cell) | energy | omega |
|---|---|---|---|---|---|---|---|---|
| 1 | `info = 0` | 8.82e-4 (250) | 1.69e-7 | 0.500 | `info = 0` | 8.82e-4 (250) | 2.72e-8 | 0.500 |
| 2 | **`info = 2`** STAGNATED | 4.73e-4 (260) | **8.64e-2** | 0.250 | **`info = 0`** | 2.86e-4 (253) | 3.05e-8 | 0.500 |
| 3 | **`info = 2`** STAGNATED | 3.70e-4 (264) | **1.05e-1** | 0.125 | **`info = 0`** | 1.06e-4 (309) | 3.79e-8 | 0.500 |
| 4 | **`info = 2`** STAGNATED, then REFUSED | 3.31e-4 (265) | **1.12e-1** | 0.125 | **`info = 0`** | 4.14e-5 (313) | 3.36e-8 | 0.500 |

The entry text ends at pass 4, "REFUSED -- the joint distance of the state
has not fallen in 3 consecutive passes". As delivered the energy row stays
at 3e-8 instead of climbing to 1e-1, `omega` is never halved, and the loop
goes on:

| pass | 5 | 6 | 7 |
|---|---|---|---|
| hydro | `info = 0` | `info = 0` | `info = 2` |
| element row (cell) | 1.96e-5 (217) | 1.15e-5 (217) | **6.37e-6 (283), inside its 1e-5** |
| energy | 3.48e-8 | 3.73e-8 | 1.03e-3 |

**The elemental side of this case reaches its tolerance at outer pass 7**,
where the entry text refused at pass 4 with that row at 3.31e-4 and the
energy row at 1.12e-1. What still refuses at pass 7 is the hydrodynamic
energy row at 1.03e-3, at cell 500, the outer boundary cell -- a different
entry from the one L4b diagnosed, and one the review's section 6 (the
subsonic outer boundary) is about. The loop was still running at pass 8
when this item closed.

That is the coupling L4b section 7 describes -- the progress control halving
the element step because the hydrodynamic energy row held the joint distance
-- and it is released.

**How far the margin may move without changing the outcome** is measured in
section 6.4.

## 5. What the row-resolved measurement says, and what it asks for next

MEASURED with `EXHALE_LINEAR_ROWS=1` on the reproduction, 8 threads. The
columns are mass, momentum, energy.

| iteration | rows entering | linear residual left, over the same row of `F` | rows leaving | cycle reached / asked |
|---|---|---|---|---|
| 1 | 5.002e-9, 4.183e-13, 1.121e-1 | **2.608e+05, 1.996e+05, 4.447e-1** | 1.169e-2, 3.119e-3, 6.486e-2 | 7.900e-4 of 1.0e-1, 1 product of 40 |
| 2 | 1.169e-2, 3.119e-3, 6.486e-2 | 7.506e-3, 4.086e-2, 3.677e-2 | 1.057e-2, 6.226e-3, 6.566e-2 | 1.425e-3 of 1.0e-1, 1 product |
| 3 | 1.057e-2, 6.226e-3, 6.566e-2 | 3.875e-3, 3.570e-1, 2.697e-2 | 1.017e-2, 8.994e-3, 6.760e-2 | 1.453e-3 of 1.0e-1, 1 product |
| 5 | 1.002e-2, 1.087e-2, 6.723e-2 | 2.648e-3, 8.572e-1, 1.676e-2 | 9.861e-3, 1.156e-2, 6.309e-2 | 1.591e-3 of 1.0e-1, 1 product |

Read it in that order. The first step is entered with the mass and momentum
rows at their floors, the cycle satisfies the tolerance it is given with
room to spare, and the linear residual it leaves is five decades larger than
the mass and momentum rows themselves. From the second step on those rows
are no longer at their floors and the same linear residual is a few parts in
a thousand of them; the momentum row's relative linear residual is of order
unity throughout, and it does not matter because that row is inside its
tolerance whatever the step does to it.

So the discriminator the review asked for is answered: the obstruction is
neither an unreached Krylov tolerance -- the cycle reaches two decades
inside what it is asked for, in one product of the forty available -- nor a
failure of the nonlinear model, whose trial the line search accepts at
`lam = 1` at every one of these iterations. It is that the quantity the
linear solve controls, one scalar 2-norm of the scaled residual, is the
energy rows, and the rows the certification reads are ten decades below it.

**An existing option that is not the row scaling this needs.** The tree
already carries `EXHALE_MODEL_ROW_EQUIL=1`, which scales the rows of the
linear model to unit infinity norm of the band
(`unit_infinity_norm_row_scaling_of_the_band`). MEASURED on the same first
step: the linear residual left in the mass row is 2.610e+05 of that row
against 2.608e+05 without it, and the step still raises the mass row to
1.171e-2 against 1.171e-2. It equilibrates the band, not the residual.

**What to investigate, NOT implemented here** (the item's scope is the two
step-control changes):

- a row scaling of the linear system by the size each row is JUDGED at --
  the cell's own continuity tolerance, the anchored rounding floor, rather
  than by the band's infinity norm -- so that the cycle's 2-norm is a norm
  in which the continuity rows are visible;
- or a continuity-constrained step: solve the Newton system with the mass
  rows imposed as a constraint rather than as three of five hundred terms of
  a norm, which is the natural statement of "this state's continuity is
  already at its floor and the step must not leave it".

Tightening the scalar tolerance further is measured not to be a route: the
cycle already delivers 7.9e-4 where 1e-1 is asked, and the mass row's own
linear residual at that point is still 2.6e+05 of itself.

## 6. The control fixtures

All at 8 threads. The control build and the measured build differ only in
`steady_newton.f90`. Both fixtures were run twice on the measured side, with
the first form of the window (section 3) and with the delivered one, and
every number below is the same in both: the element reload reports 843 outer
iterations, 9 solves ending "no descent direction exists", 0 STAGNATED and
the identical element ladder in each, and the carrier reload reports the
identical row values pass for pass. The tables are the delivered build's.

### 6.1 `backup/regression/atomic_elem_newton`, the element reload

HD 209458 b, atomic with solar trace metals and binary element diffusion,
reloaded from `IC/` with `Load IC? True`, `Coupled carrier solve: False`,
`Restart intent: stationary`; `EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80
EXHALE_DIFF_OMEGA=0.5`, 8 threads (the README's own recipe).

**CERTIFIED at outer pass 12 in both builds**, "every active equation of
this state is within its own tolerance", which is what the README records
for this fixture since item P16.

The elemental side is IDENTICAL pass for pass, row, value and cell:

| pass | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| worst gated species row, both builds | 2.61e-3 | 3.81e-3 | 1.76e-3 | 9.31e-4 | 4.61e-4 | 2.31e-4 | 1.19e-4 | 5.76e-5 | 3.11e-5 | 1.41e-5 | 1.04e-5 | 7.73e-6 |

at cell 268 (O) through pass 10 and cell 289 (He/H) at 11 and 12 in both,
with `omega` 0.500 at pass 1 and 0.250 from pass 2 in both. The refusing
entries of the final certification are the same three rows to four digits
(He/H 1.371e-5 and 1.038e-5, O 1.411e-5). So the change stays inside the
hydrodynamic solve and does not touch the relaxation.

What moved is inside that solve: the hydrodynamic `info` of passes 1 and 11
goes from 1 and 2 to 0 and 0, solves ending "no descent direction exists"
fall from 11 to 9, and neither build reports a STAGNATED solve at all. The
hydrodynamic rows scatter at the 1e-10 to 1e-9 level (mass 2.09e-9 -> 1.16e-9
at pass 1, 9.26e-10 -> 1.26e-9 at pass 12), which is the rounding floor this
fixture's README already records as carrying no signal.

### 6.2 `backup/regression/carrier_model_a_newton`, the carrier reload

`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40 OMP_NUM_THREADS=8 ./run.sh`,
the README's recipe.

**Identical verdict and identical rows.** Twelve passes, every hydrodynamic
solve `info = 0` in both builds; the worst gated species row of the last
three passes is 1.23e-3 at cell 304, 1.68e-3 at cell 248 and 7.42e-4 at cell
323 in both, to three digits; both end NOT CERTIFIED on the same two
entries, carrier balance H+ 7.425e-4 at cell 323 and carrier balance H2
6.996e-4 at cell 363, to four digits. The hydrodynamic rows move at 1e-11 to
1e-12, the last-bit level of this fixture.

### 6.3 The regression cases that reach this code path

`REGRESSION_EXE=<binary> REGRESSION_GOLDEN_DIR=<the tree's golden>
./run_check.sh check wasp_full mol_base_handoff wasp_full_newton`, run on a
copy of the harness and of the three case directories under the scratch
directory, never in `backup/regression/` itself. With `REGRESSION_EXE` set
the script skips both the rebuild and the `make -q` up-to-date test (its own
lines 344 to 350), so nothing in the delivered tree is built or written; that
is what makes it safe to run while a campaign holds `EXHALE.x`.

Of the three, `wasp_full` carries no `Solver:` key at all, so it never
enters `solve_steady_jfnk`. MEASURED rather than argued: run at
`EXHALE_MAXSTEPS=400`, one thread, with the control build and with the
delivered build, it writes `Hydro_ioniz.txt` and `Ion_species.txt` BITWISE
IDENTICAL between the two and reports no Newton iteration at all. The
uncapped case marches for well over 13 000 steps and was stopped once the
capped pair had settled it. `mol_base_handoff` and `wasp_full_newton` both
run `Solver: Newton`.

## 7. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
D=$EX/LHS1140b/models/atomic_scalarCNO_gj1132_kzz1e9/HeH2.13
W=$EX/scratchpad/L4d/run                 # any scratch directory

mkdir -p $W/output
\cp -f $D/input.inp $D/base.inp $W/
sed -i "s|^Spectrum file: .*|Spectrum file: $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt|" \
    $W/input.inp

# THE REPRODUCTION: the state every outer pass after the first is given
\cp -f $D/output/Hydro_ioniz.txt $W/output/Hydro_ioniz_IC.txt
\cp -f $D/output/Ion_species.txt $W/output/Ion_species_IC.txt
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 \
         EXHALE_JFNK_MAXIT=300 <binary>
#   entry text:  STAGNATED at 40, entry state returned
#   as delivered: info = 0 at 98

# the row-resolved measurement of section 5 (twelve iterations are enough)
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 \
         EXHALE_JFNK_MAXIT=12 EXHALE_LINEAR_ROWS=1 <binary>

# THE FULL OUTER LOOP from this case's own seed
\cp -f $D/output/Hydro_ioniz_IC.txt $D/output/Ion_species_IC.txt $W/output/
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 <binary>

# the reproducibility of the residual assembly (section 3): same state,
# 4, 8 and 16 threads, two Newton iterations, compare output/Hydro_ioniz.txt
```

The fixtures are their own READMEs' recipes, run on scratch copies:
`backup/regression/atomic_elem_newton` with `Load IC? True`,
`Coupled carrier solve: False`, `Restart intent: stationary` and
`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80 EXHALE_DIFF_OMEGA=0.5`;
`backup/regression/carrier_model_a_newton` with
`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40` through its `run.sh`.

## 8. Scope, and what was not measured

- **What was verified is the path the change touches**: the stationary solve
  of the three-unknown (partitioned) route, on the reproduction, on the full
  outer loop of that case, on the two stationary fixtures, and on the
  regression cases that carry `Solver: Newton`. The trust-region branch
  (`use_tr`, the coupled route with species rows) shares the stagnation
  counter and was NOT run here; the forcing term does not reach it, its
  Krylov tolerance being `tr_gm_rtol` and unchanged.
- **The control and the measured binary** were both built from the tree with
  only `steady_newton.f90` differing, and both before another worker's edits
  to `EXHALE_main.f90`, `load_IC.f90` and `write_output.f90` reached the tree
  at 16:22 on 2026-09-13. A binary built from the delivered source after
  that time would carry those edits, so the final check was made from an
  isolated copy of the tree (`git archive HEAD` plus this file only): it
  reproduces the measured binary's ladder on the reproduction to every
  printed digit over its first 40 iterations.
- **The LHS 1140 b closure rungs of L4b were not re-run.** The reproduction
  measured here is L4b's own isolating control, which fails at the same
  outer pass, in the same row and with the same ledger signature as the
  `HeH9` rung (L4b section 5 e); the rung itself needs the Photochem column
  at every iteration and is a campaign-scale run.
- **No golden was refreshed and none was compared against by this item's
  own decision**: the regression verdicts of section 6.3 are read from the
  harness against the tree's `golden/`.
- Wall times are contended: the machine ran a campaign of this tree
  throughout. Compare named outcomes and row values, not seconds.

## 9. The risk the change carries, and what it measured

A counter that is also reset by the merit can in principle keep a solve
running on micro-improvements until the outer cap. MEASURED, it did the
opposite on both fixtures: the total number of outer iterations over the
twelve passes falls from 909 to 843 on the element reload and from 133 to
131 on the carrier reload, and neither build reports a single STAGNATED
solve on either. The two limits that do not depend on the counter are
unchanged and still hold: `n_no_descent_max = 12` consecutive iterations
with no acceptable step, and the outer cap itself.

| case | control build | as delivered |
|---|---|---|
| `wasp_full_newton` | `info = 0`, `||R||` 4.213e-9, 13 outer iterations, log Mdot 13.30; `Hydro_ioniz.txt` and `Ion_species.txt` **bitwise identical to `golden/`** | `info = 0`, `||R||` 7.860e-9, 12 outer iterations, log Mdot **13.30**; WITHIN the 1e-3 gate at a maximum relative difference of **1.3e-8** (`Hydro_ioniz.txt`) and **2.5e-8** (`Ion_species.txt`) |
| `mol_base_handoff` | bitwise identical to `golden/` | **bitwise identical to `golden/`**; the case is a fixed-step snapshot (`maxsteps`) that ends at `count = 12000`, `du = 7.9073e-1`, log Mdot 10.57 in both and never reaches the stationary solve |

The control build reproducing both goldens bitwise is what says it is the
build the goldens were taken from. The movement of `wasp_full_newton` is
2.5e-8 relative, five decades below the 0.1 percent the standing rule calls
identical, so **no golden is touched and none needs to be.**

### 6.4 How far the margin may move

The reproduction at `EXHALE_MERIT_FALL_MIN` = 0.05, 0.10 (the default) and
0.20, a factor of four apart: `info = 0` at iteration 98, `||R||` 2.799e-08,
flux spread 1.706e-12, the same in all three to every printed digit. The
threshold is inert over that range on this reproduction, which is what says
it was not chosen to recover a count; what it does choose is the grind of
section 3, which the earlier per-iteration form at 1e-3 let run past 250
iterations.

## 10. What this item did not settle

The rung's outer pass 8 was still running when the item closed and was
stopped there, at outer iteration 279, `||R||` moving between 2.23e-4 and
4.30e-4 over the last twenty iterations, with the binding row the
hydrodynamic energy row of cell 500 -- the outermost cell, at the subsonic outer
boundary, and not the band of C I cooling L4b diagnosed. Two statements
follow and only the first is measured here:

- the refusal L4b measured (the stagnation detector stopping every pass
  after the first) is gone, and the elemental side of this case reaches its
  tolerance;
- what the loop now stands on is a row at the outer boundary. The review's
  section 6 asks for an acceptance item there (the isothermal hydrostatic
  continuation into the outer ghosts is written for an outflow with no
  incoming characteristic, and this wind is subsonic to 30 R_p, where `v-c`
  enters the domain). Nothing in this item measured that boundary, and no
  claim is made about it.

The grind of section 3 is also not fully cured, only bounded: outer pass 8
descends slowly enough to keep the counter resetting (its merit falls 12 to
45 percent over a 20-iteration window against the 10 percent margin) while
`||R||` oscillates rather than falls. A counter that reads the
merit cannot separate "slow" from "stuck" better than the margin lets it;
separating them properly is the row-scaling question of section 5, not a
threshold.

## 11. Noticed and fixed, outside this item's question

`steady_newton.f90` carried the word "pipeline" in eight comments and in the
name of one diagnostic routine. The word is on the forbidden list of
`docs/worker_rules.md` and of the user's standing rules, and a routine named
for a software structure rather than for what it computes is against the
naming rule of the project file. The routine is renamed
`residual_terms_in_evaluation_order`, which is what it returns -- every
quantity the row of one cell is built from, in the order the evaluation
builds them -- and the comments now say "the evaluation" or "the terms".
The rename is confined to this file; no other source refers to it. MEASURED:
a build of the delivered source reproduces the measured build's ladder on
the reproduction over twelve iterations to every printed digit.

Two occurrences of the same word remain outside this item's file and are
NOT fixed here: `src/modules/functions/binary_element_diffusion.f90` line
861 (already reported by L5b) and two comments in
`src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90`.
