# LHS 1140 b, item L7e: the H2 carrier balance at r = 1.6 R_p, term by term

Item L7e of `docs/PLAN_20260913_lhs_stationary.md`, carried out on the live
tree while `make check` and a golden refresh were running in it. Nothing in
`build/`, `EXHALE.x`, `backup/regression/` or any catalogue case directory of
`LHS1140b/models/` was written by this item: every binary below is a private
build in its own object directory, every run is under
`LHS1140b/models/.L7e/`, and the regression matrix was run on a copy of
`backup/regression/` in the scratch directory.

| build | what it is | md5 |
|---|---|---|
| `EXHALE_L7e_ctl.x` | the tree as delivered to this item (the control) | `a11038c245050d4f11852af829e13fc7` |
| `EXHALE_L7e_diag.x` | the control plus the row dump of section 3 alone, nothing else | `d910f1330eb06d7bb3fec43d43c7c4e3` |
| `EXHALE_L7e_s9.x` | the dump, the two closure corrections of section 5 and the reports of section 4; the build the section 6 record was made with | `6b753c3dc7047a8ddedcff813e4095d9` |
| `EXHALE_L7e.x` | the delivered source text: the above plus section 11 (the progress-control correction and the three refusal reports; the e-fold clause of 11.1 is NOT in it) | `95e06b34bfc83421b5b5796c85ced2cf` |

Every number is MEASURED (produced by a run of one of those builds, or by
arithmetic stated here on a file such a run wrote) unless it is marked READ
(from a source file or a document). Two comment blocks were corrected after
`EXHALE_L7e.x` was built; the final source text was rebuilt into its own
object directory and gives the same md5, so the delivered text is the text
every measurement below was made with.

## 1. Verdict

**The measure of exactly 1.000 is neither a masked chemistry, nor a defective
row scale, nor a floor. The H2 row at cell 280 carries one term that nothing
in the row can balance: the chemical PRODUCTION of H2, 3.598e-01 cm^-3 s^-1,
which stands 5.2 decades above the total H2 loss (2.154e-06) and 7.8 decades
above the transport divergence (5.5e-09). The residual IS that production
term, so the measure is 1.000 on the scale the code uses and 0.99999 on a
scale that sums the magnitudes of every term instead — the scale is not the
question.** All of the production is one reaction, R9, `H2+ + H -> H2 + H+`:
`k9 n(H2+) n(H I)` reproduces the dumped production to four digits at every
cell examined. The H2+ that feeds it is not made from H2 — at these cells
n(H2+) stands 5e3 times above n(H2) — but from the atomic gas, through
He(2^3S) + H -> HeH+ + e (the associative branch, 1 - f_penning of the
He(2^3S)+H ionization) followed by HeH+ + H -> H2+ + He (R19). The row's own
root is n(H2) = 1.45e+02 cm^-3, a factor 1.7e+05 above the 8.69e-04 the state
carries.

**The carrier relaxation does not reach that root because on this route it
never takes a single transport step.** `relax_photochemical_composition`
refuses every trial, down to the shortest admissible one, and hands back its
entry composition; the refusal is a property of the ENTRY state and not of
the trial, so it recurs identically on every outer pass and the alternation's
fixed point is the seed. Two independent gates do it, both measured here:

1. **The chemistry closure read a count that describes the root search and
   not the state it hands back.** `equilibrate_chemistry_at_fixed_conserved_state`
   refused whenever the sweep's ledger reported `n_offsimplex > 0`. That
   counter counts cells at which no starting point stayed inside the element
   simplex — and BOTH branches that raise it end in an admissible
   composition, the molecular one by projecting the closest root onto the
   element budget and the atomic one by handing back the uncoupled ionization
   balance, each rechecking its reaction residual. A root that sits ON A FACE
   of the simplex raises it by construction, and a fully dissociated H2 in a
   6000 K wind cell is exactly such a root. MEASURED: 12 cells, all of them
   molecular clamps, worst element-budget excursion 6.54e-04.

2. **The closure's tolerance sits below the equilibrium sweep's own noise.**
   Traced cycle by cycle, the closure contracts for two or three cycles and
   then stops: 7.6e-05 -> 1.1e-05 -> 2.8e-06, and from there a band 1.7e-06
   to 3.4e-06 that does not narrow over thirty further cycles. The cell that
   carries it is not the molecular layer but the far wind at r = 10 to
   11 R_p. Against `chem_cycle_tol` = 1e-6 and `chem_cycles_max` = 5 every
   call spent its budget.

With both corrected the carrier moves for the first time on this route
(section 6).

## 2. What was measured, and on which state

The state is the one item L7d left: the `local` partition of
`molecular_scalar_gj1132_kzz1e9/HeH2.13` after the partitioned stationary
route, reloaded as an initial condition in `LHS1140b/models/.L7e/dump_local`
(`output/Hydro_ioniz.txt` and `Ion_species.txt` of `.L7d/t_local` copied to
`*_IC.txt`). Run with `EXHALE_L7e_diag.x`, one outer pass, one JFNK
iteration, so that the certification measures essentially that state; it
reproduces L7d's refusal exactly: `carrier balance H2: gated row measure
1.000E+00 above 1.0E-05 at cell 280 (a wind cell)`.

## 3. The row dump

The code carried no cell-by-cell record of a carrier row's terms (item L7c
section 5.2 says so). One was added, default off:

```
EXHALE_CARRIER_ROW_TERMS=1   ->   output/carrier_row_terms.txt
```

written by the certification's own isolated evaluation, so the terms belong
to the state that received the verdict. One line per (cell, solved carrier):

```
cell r[R_p] T[K] carrier n_c x_c diffusive advective production loss
photo_loss net_source residual floor scale_net measure_net
scale_terms measure_terms
```

`scale_net` is the scale the code uses, |transport| + |production - loss| +
floor; `scale_terms` is |transport| + production + loss + floor, the sum of
the magnitudes of the terms the row contains. Both are printed so that the
difference can be read rather than argued.

MEASURED on the state of section 2 (all rates cm^-3 s^-1):

| cell | r [R_p] | T [K] | n(H2) | production | loss | photo-loss | transport | residual | measure_net | measure_terms |
|---|---|---|---|---|---|---|---|---|---|---|
| 194 | 1.1340 | 3911 | 5.579e+01 | 8.761e+00 | 1.152e-01 | 1.641e-05 | -6.49e-03 | -8.652e+00 | 1.000000 | 0.974059 |
| 218 | 1.2042 | 5231 | 2.290e-01 | 4.834e+00 | 8.890e-04 | 7.271e-08 | -1.18e-05 | -4.833e+00 | 1.000000 | 0.999632 |
| 250 | 1.3578 | 6258 | 4.652e-03 | 1.405e+00 | 1.997e-05 | 1.608e-09 | -7.42e-08 | -1.405e+00 | 1.000000 | 0.999972 |
| 280 | 1.6044 | 6032 | 8.692e-04 | 3.598e-01 | 2.154e-06 | 3.176e-10 | -3.20e-09 | -3.598e-01 | 1.000000 | 0.999988 |

(transport is the diffusive plus the advective divergence; the floor of the
scale is 5.2e-12 at cell 194 and 6.2e-14 at cell 280, so it decides nothing.)

### 3.1 The four hypotheses the item named

- **(a) the chemistry is masked off in that cell** — REFUTED. Production and
  loss are both finite and nonzero at every cell of the table, and the loss
  carries the temperature dependence a thermal dissociation should: the loss
  rate per H2 molecule is 2.06e-03 s^-1 at 3911 K, 3.88e-03 at 5231 K,
  4.29e-03 at 6258 K, 2.48e-03 at 6032 K.
- **(b) the row scale counts only the transport term** — REFUTED. The scale
  is |transport| + |net source| + floor and the net source is the term that
  dominates it; replacing it by the sum of the magnitudes of every term moves
  the measure from 1.000000 to 0.999988 at cell 280 and to 0.974059 at cell
  194. **The scale is therefore left alone, and that is a judgment and not an
  omission**: |transport| + |production - loss| is what the equation the row
  states — transport equals net source — is made of, while a scale built on
  production + loss would call a cell converged when the fast chemistry is
  nearly balanced and the H2 CONTENT is arbitrary, which is the pathology
  `input_read.f90` already refuses a coupled solve for (the measured 0.77 of
  any seed perturbation surviving every sweep).
- **(c) the row is at a floor** — REFUTED. The absolute floor of the scale is
  6.2e-14 at cell 280 against terms of 3.6e-01, twelve decades below, and the
  carrier is not marked absent anywhere in the table.
- **(d) something else** — CONFIRMED, and it is section 1.

### 3.2 Where the production comes from

`k9 = 6.4e-10` cm^3 s^-1 (R9, `H2+ + H -> H2 + H+`, READ from
`mol_rates.f90`), and the H2+ and H I of the same state (READ from
`output/Ion_species.txt`):

| cell | n(H2+) | n(H I) | k9 n(H2+) n(H I) | production dumped |
|---|---|---|---|---|
| 194 | 4.2402e+00 | 3.2282e+09 | 8.760e+00 | 8.761e+00 |
| 218 | 6.4860e+00 | 1.1646e+09 | 4.834e+00 | 4.834e+00 |
| 250 | 6.1857e+00 | 3.5492e+08 | 1.405e+00 | 1.405e+00 |
| 280 | 4.4353e+00 | 1.2675e+08 | 3.598e-01 | 3.598e-01 |

so R9 is the whole production and the other three channels of the H2 row
(`k6 n_e n(H3+)`, `k11 n(H3+) n(H I)`, `k15 n(H I)^2`) contribute nothing
measurable. n(H2+) is 5.1e+03 times n(H2) at cell 280, so the H2+ is not
supplied by H2; the atomic channels that supply it are
He(2^3S) + H -> HeH+ + e and HeH+ + H -> H2+ + He (R19, `k19 = 9.1e-10`),
with n(HeH+) = 5.84 cm^-3 and n(He 2^3S) = 45.6 cm^-3 at that cell.

**So the row is not asking for something impossible: it is asking for the H2
that this network's own non-thermal formation chain makes.** The thermal pair
in the row (R15 three-body association and R12 dissociation, R12 built as
R15 divided by the equilibrium constant, so the pair satisfies detailed
balance by construction) would give the thermochemical H2 the `local` seed
carries; the ion chain adds a source the thermochemistry does not have, and
the row's root sits five decades above it.

## 4. Why the relaxation never moves: the measurement

`relax_photochemical_composition` is the carrier's half of the outer
alternation (`EXHALE_main.f90`); `photochemical_transport_step` is called
nowhere else on the stationary route, and `carrier_in_newton` is F for this
case (READ from `EXHALE_resolved.out`), so it is the only thing that can move
H2. Every outer pass of L7d's run reports

```
carrier relaxation ended on the chemistry of the shortest admissible trial
did not close; drift 0.00E+00 in 0 transport steps
```

and the code named no reason. A report was added — one line, printed only
when the relaxation ends that way — giving the closure's reason and the
ledger counts behind it, and a trace of the closure at each cycle under the
existing `EXHALE_CARRIER_DEBUG=1`. With them:

```
the closure refused it: a cell left the element simplex; cells nonfinite 0,
off the element simplex 12, molecular clamps 12, worst violation 6.54E-04,
last |dT|/T 1.80+308
```

and after the first correction, on the same state:

```
the closure refused it: cycle budget spent; cells nonfinite 0, off the
element simplex 0, molecular clamps 0, worst violation 4.11E-04,
last |dT|/T 3.68E-05
```

The cycle trace of the second (`EXHALE_CARRIER_DEBUG=1`, a probe build with
the budget raised to 40 for the measurement only):

| cycle | 1 | 2 | 3 | 4 | 5 | 6 | 7 | ... | 22 |
|---|---|---|---|---|---|---|---|---|---|
| max abs(dT)/T | 2.16e-04 | 8.83e-05 | 2.65e-05 | 4.92e-05 | 2.37e-05 | 3.59e-06 | 2.15e-06 | 1.7e-06 to 3.4e-06 | 7.92e-07 |

with the cell that carries the maximum at r = 7 to 11 R_p throughout. Five
further closures of the same run reach the band at cycle 3 and then need 8,
22, 9 and more than 31 cycles to dip below 1e-6; one did not dip within 31.

## 5. The corrections

### 5.1 The closure reads the state it hands back, not the root search

`src/modules/lower_atmosphere/diffusive_photochemistry.f90`,
`equilibrate_chemistry_at_fixed_conserved_state`: the refusal on
`ledger%n_offsimplex > 0` is removed, with the reasoning and the measurement
above it in the code. What still refuses is unchanged: a cell whose accepted
state broke a balance row out of the reals (`n_nonfinite`), a composition
that is not a number, a thermal state that is not a gas, and a temperature
that has not stopped moving. The reason code `chem_closure_off_simplex` can
no longer be returned and is retired at its declaration (the value is left
unused rather than reassigned, so the other reasons keep their numbers).

### 5.2 The budget and the tolerance are set from the measurement

Same file, the two parameters, with the trace of section 4 recorded at the
declaration:

```
chem_cycles_max   5     -> 12      (twice the 6 cycles the slowest closure of
                                    the measurement needed to reach the band)
chem_cycle_tol    1.0e-6 -> 1.0e-5 (a factor of three above the 3.4e-6 the
                                    band reaches, which is the equilibrium
                                    sweep's own acceptance noise in the far
                                    wind and not a state still moving)
```

### 5.3 The row dump, and the production/loss split it needs

- `src/modules/nonlinear_system_solver/System_HeH_mol.f90`: `mol_heh_rows`
  gains optional `p_Hp, l_Hp, p_H2, l_H2, l_H2_phot` and
  `oxygen_carrier_rows` optional `p_OH, l_OH, p_H2O, l_H2O, p_H2_oxy,
  l_H2_oxy`. The rows themselves are now assembled FROM those names, in the
  same order and grouping they were written in, so the residual is unchanged
  to the last bit and the row and its terms cannot drift apart.
- `src/modules/lower_atmosphere/diffusive_photochemistry.f90`:
  `carrier_source` gains optional `sprod, sloss, sphot`; `carrier_residual`
  records the terms of every row when the switch is on;
  `carrier_row_terms_on` and `carrier_row_terms_write` are the switch and the
  writer.
- `src/modules/time_step/certification.f90`: `carrier_rows_of_state` writes
  the file after its own evaluation.
- `src/EXHALE_main.f90`: the closure-refusal report of section 4.
- `src/tests/carrier_retry/carrier_retry.f90`: a comment that stated the
  retired refusal as the expected verdict, corrected.

**Bitwise check of the split.** `EXHALE_L7e_diag.x` (the split and the dump,
nothing else) and `EXHALE_L7e_ctl.x` were run on the same reloaded state with
the same environment: `output/Ion_species.txt` is byte-identical and
`output/Hydro_ioniz.txt` differs in one line, the provenance timestamp.

## 6. The case

`molecular_scalar_gj1132_kzz1e9/HeH2.13`, run by `run_case.sh` on a copy in
`LHS1140b/models/.L7e/`, seeded from the frozen atomic state
`.L7d/src_atomic_HeH2.13` with `SEED_X2=local`, `EXHALE_BIN=EXHALE_L7e.x`,
8 threads. (The copy's `Spectrum file:` was made absolute: `run_case.sh`
resolves a case relative to `models/<group>/<case>` and the extra `.L7e`
level breaks the relative path. Nothing else in the copy differs from the
catalogue case.)

### 6.1 What the outer loop does now

MEASURED, the outer passes as `run.log` records them (`hydro info`, the
worst gated species row, and what the carrier half of the alternation did):

| pass | hydro info | worst gated row | at cell | carrier relaxation ended on | drift | transport steps |
|---|---|---|---|---|---|---|
| 1 | 0 | 1.000 | 279 | the shortest admissible trial was not covered | 0 | **0** |
| 2 | 0 | 1.000 | 280 | the movement bound | 0 | **0** |
| 3 | 0 | 1.000 | 280 | the movement bound | 0 | **0** |
| 4 | 2 | 1.000 | 280 | the movement bound | 9.87e-03 | **6** |
| 5 | 0 | 9.97e-01 | 217 | the movement bound | 9.81e-03 | **10** |
| 6 | 0 | 9.94e-01 | 217 | the shortest admissible trial was not covered | 4.14e-03 | **16** |
| 7 | 0 | 9.93e-01 | 217 | the shortest admissible trial was not covered | 2.44e-03 | **9** |

**The carrier transports for the first time on this route.** Against every
pass of L7d and of every run before the correction -- 0 steps, drift
identically zero, the state a fixed point of the alternation because nothing
moved it -- the carrier now takes 6 to 16 transport steps a pass and the
gated row leaves 1.000 for the first time. The hydrodynamic rows are
unaffected (mass 7.2e-07, momentum 5.4e-11, energy 6.4e-06 throughout, the
numbers L7d reports), and the count of refusing certification entries falls
from 5 at entry to 1 by pass 7: the H2 carrier row alone.

**The case is NOT CERTIFIED**, and the row that refuses it is still the H2
carrier balance. That is expected from section 1 and is not a surprise: the
row's root is 1.7e+05 times the carried value, and what now sets the pace is
the third gate, the movement bound.

### 6.2 The gate that now holds it, named and not fixed

`relax_photochemical_composition` refuses a trial whose returned state moves
the carrier by more than `trust` of the LARGEST H2 mixing ratio of the entry
state, one number over the whole grid
(`carrier_composition_displacement`), with `trust` starting at 1e-2 and
floored at 1e-3 (`EXHALE_main.f90`, anchored on the hot-Uranus hand-off
measurement of `docs/p50_carrier_wind_alternation.md`). The largest H2
mixing ratio of this column is the BASE value, 0.19, so the bound admits an
absolute change of 1.9e-03 anywhere -- while the wind cells that carry the
refusing row need a change of order 1e-06 in the same units and can never be
what exceeds it. Passes 2 to 5 end on that bound, and the H2 front does not
move at all (x2 = 0.5 at r = 1.0583 and x2 = 1e-2 at r = 1.0722 at every
pass from 4 to 7), so the trial length is being set by a part of the column
that is already at its fixed point.

Two further observations, both MEASURED and neither acted on here:

- The outer loop's response to a pass that made no progress is to HALVE the
  movement bound (`trust_pass = max(0.5 trust_pass, 1e-3)`). Where the
  carrier took no step at all that cannot help, and it tightens the only
  thing preventing one; the response would need to tell "the carrier moved
  and the alternation did not contract" from "the carrier could not move".
- Passes 1, 6 and 7 end instead on `the shortest admissible trial was not
  covered`, the transport operator's own substep budget. L7d's run reports
  the same thing from the marching side (28 uncovered intervals, worst row
  cell 411 at |res|/physical terms 1.36e-05).

**Both are solver-control questions about where a bound is measured, not
physics, and they are left for a decision rather than changed here.**

### 6.3 The run this section reports, and why it was stopped

The table above is passes 1 to 7 as
`LHS1140b/models/.L7e/molecular_scalar_gj1132_kzz1e9/HeH2.13/run.log` has
them. **The run was stopped by hand during outer pass 8 of the cap of 20**:
its trajectory was already flat -- the hydrodynamic rows unmoved to three
digits from pass 2 and the carrier advancing under the movement bound at a
rate the next twelve passes could not change -- and it was holding eight
threads against the measurement of section 12, which is the one that
decides anything. The directory holds the log up to that point and no
written state, since the run never reached its post-processing pass.

## 7. Tests and regression

### 7.1 The suites the change touches

Each built into its own object directory, so the tree's `build/` was never
written; the module tree they compile against is `build_L7e`, the delivered
one.

| suite | why it is in scope | result |
|---|---|---|
| `carrier_retry` | the closure's reasons and the relaxation's retry controller | PASSED, every assertion |
| `carrier_reference_scales` | `carrier_source`'s signature and the row scales | PASSED |
| `carrier_constraint_attribution` | `carrier_source` again | PASSED |
| `carrier_returned_state_acceptance` | the returned-state measure of a carrier row | PASSED |
| `certification` | the carrier rows of a state and the isolated evaluation the dump is written from | PASSED, every assertion |
| `molecular_seed` | the seed the case is built from, run against `EXHALE_L7e.x` | PASSED |
| `steady_species_rows` | the elemental and carrier rows of the stationary system | 195 PASS, 0 FAIL |
| `physics_probe` | it calls `mol_heh_rows` directly (`co_helium_ion_sink.f90`), which the production/loss split touches | 1491 PASS, 0 FAIL |

`carrier_retry` carries a comment that named the retired refusal as the
expected verdict of one of its passes; it was corrected with the change
(section 5.3) and the suite's assertions are unchanged.

### 7.2 Regression

`REGRESSION_EXE=EXHALE_L7e.x` on `mol_carrier`, `mol_base_handoff` and
`mol_diffusion`, run on a copy of `backup/regression/` because the tree's own
`make check` held the lock:

```
==> REGRESSION PASS (identical, or within 1e-3 relative)
```

with the largest relative movement of any column 1.272e-07 on
`Hydro_ioniz.txt`, 1.814e-06 on `Ion_species.txt` and 2.1e-05 to 2.5e-05 on
the `*_adv.txt` profiles, in all three cases. **All of that is below 0.1
percent and by the standing rule counts as identical; no golden is refreshed
by this item.**

**The movement is not this item's.** The same three cases were run again
with `EXHALE_L7e_ctl.x`, the control, for a binary-against-binary comparison
rather than a comparison against a golden the tree's other uncommitted work
has already moved. `mol_carrier`, the case that carries the transported H2:
`Hydro_ioniz.txt`, `Ion_species.txt` and both `*_adv.txt` profiles are
**byte-identical** between the two binaries (the provenance timestamp
excepted), and every `PASS` line of the two matrix runs -- all twelve, for
all three cases -- carries the same numbers to the last digit printed. So
the 1e-07 to 2e-05 above is the distance between the tree and its goldens,
not the distance this item added.

None of the three cases can reach the code this item changed: each is a
relaxation snapshot pinned at 12000 marching steps (`maxsteps`), their
`run.log` shows the JFNK hand-off announced and never entered, and
`equilibrate_chemistry_at_fixed_conserved_state` has exactly one caller,
`relax_photochemical_composition`, which only the stationary route enters.

## 8. How to reproduce

`EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00`.

```bash
cd $EX
make -j8 OBJDIR=build_L7e_ctl EXE=EXHALE_L7e_ctl.x    # before any edit
make -j8 OBJDIR=build_L7e     EXE=EXHALE_L7e.x        # after

# the row dump, on the state item L7d left
cd $EX/LHS1140b/models
mkdir -p .L7e/dump_local/output
\cp -f .L7d/t_local/input.inp .L7d/t_local/base.inp .L7e/dump_local/
\cp -f .L7d/t_local/output/Hydro_ioniz.txt .L7e/dump_local/output/Hydro_ioniz_IC.txt
\cp -f .L7d/t_local/output/Ion_species.txt .L7e/dump_local/output/Ion_species_IC.txt
( cd .L7e/dump_local && OMP_NUM_THREADS=8 EXHALE_CARRIER_ROW_TERMS=1 \
     EXHALE_OUTER_PASSES=1 EXHALE_JFNK_MAXIT=1 $EX/EXHALE_L7e_diag.x > run.log 2>&1 )
#   -> .L7e/dump_local/output/carrier_row_terms.txt

# the closure trace of section 4 (the same directory, with the trace on)
#   OMP_NUM_THREADS=8 EXHALE_CARRIER_DEBUG=1 EXHALE_OUTER_PASSES=2 \
#     EXHALE_JFNK_MAXIT=1 $EX/EXHALE_L7e.x
# the 40-cycle probe is the same build with chem_cycles_max raised to 40 in
# diffusive_photochemistry.f90, built to its own object directory and reverted.

# the case
cd $EX/LHS1140b/models
OMP_NUM_THREADS=8 EXHALE_BIN=$EX/EXHALE_L7e.x \
  ATOMIC_SEED=$PWD/.L7d/src_atomic_HeH2.13 SEED_X2=local \
  ./run_case.sh .L7e/molecular_scalar_gj1132_kzz1e9/HeH2.13

# the suites (their own object directories, so the tree's build/ is untouched)
EXHALE_TEST_OBJDIR=<scratch>/<name> EXHALE_OBJDIR=$EX/build_L7e \
   src/tests/<name>/run.sh          # carrier_retry, carrier_reference_scales,
                                    # carrier_constraint_attribution,
                                    # carrier_returned_state_acceptance
EXHALE_TEST_OUT=<scratch>/certification EXHALE_OBJDIR=$EX/build_L7e \
   src/tests/certification/run.sh
EXHALE_TEST_OUT=<scratch>/molseed EXHALE_EXE=$EX/EXHALE_L7e.x \
   src/tests/molecular_seed/run.sh
EXHALE_TEST_OUT=<scratch>/ssr EXHALE_OBJDIR=$EX/build_L7e \
   src/tests/steady_species_rows/run.sh

# the regression, on a copy of backup/regression/ because the tree's own
# make check was running and the harness refuses two runs in one tree
REGRESSION_EXE=$EX/EXHALE_L7e.x <scratch>/reg/backup/regression/run_check.sh \
   check mol_carrier mol_base_handoff mol_diffusion
```

## 9. Scope, and what was not measured

- The dump and the production/loss split are inert unless
  `EXHALE_CARRIER_ROW_TERMS=1` is set, and the split was shown byte-identical
  through the code (section 5.3).
- The two closure corrections change the stationary route with a transported
  carrier and nothing else: `equilibrate_chemistry_at_fixed_conserved_state`
  is called only from `relax_photochemical_composition`, which the marching
  path never enters.
- **What this item did NOT judge**: whether the H2 the ion chain makes is the
  right amount. The chain is He(2^3S) + H -> HeH+ + e at
  `1 - f_penning_HeI23S` of the He(2^3S)+H ionization, then R19 and R9; each
  reaction is in the network with its published source, but the item measured
  the row and not the branching ratios behind it, and n(H2) = 1.45e+02 cm^-3
  at 6000 K is a prediction of that chain that nothing here checks against a
  published model. It is the obvious next question and it is a physics
  question, not a solver one.
- The `handoff` partition was not run: L7d's four independent measurements
  that it returns its entry state stand, and `SEED_X2` defaults to `local`.

## 10. The chain that makes H2 in a 6000 K wind, and what it means

The row of section 3 is not asking for something the network cannot supply.
Its production is the last step of a chain that starts in the atomic gas and
never passes through molecular hydrogen. MEASURED at cell 280 (r = 1.6044
R_p, T = 6032 K) on the state of section 2, with the densities READ from
`output/Ion_species.txt` and the rate coefficients READ from their source
files:

| step | reaction | rate coefficient | densities | rate [cm^-3 s^-1] |
|---|---|---|---|---|
| 1 | He(2^3S) + H -> HeH+ + e | (1 - f_penning) Q31, f_penning = 0.9, Q31 = 1.337e-09 (Taylor et al. 2025 fit, `ioniz_HeI23S_H`) | n(He 2^3S) = 45.56, n(H I) = 1.2675e+08 | **7.72e-01** |
| 2 | HeH+ + H -> H2+ + He | k19 = 9.1e-10 (R19) | n(HeH+) = 5.840, n(H I) | **6.74e-01** |
| 3 | H2+ + H -> H2 + H+ | k9 = 6.4e-10 (R9) | n(H2+) = 4.435, n(H I) | **3.60e-01** |
| — | H2 destruction, all channels | 2.478e-03 s^-1, dominated by R12 (thermal dissociation, built as the three-body association R15 divided by the equilibrium constant, so the pair satisfies detailed balance) | n(H2) = 8.692e-04 | 2.15e-06 |

The chain closes quantitatively: step 2 is 0.87 of step 1 (the rest of the
HeH+ goes to its own recombination), step 3 is 0.53 of step 2 (the rest of
the H2+ goes to dissociative recombination, `k5 n_e` at n_e ~ 1.1e+07 is 9.5
times `k9 n(H I)`), and step 3 reproduces the dumped production to four
digits. So **about half of every associative ionization of H by metastable
helium ends as a hydrogen molecule**, and the H2 balance of that cell is
"the He(2^3S) reservoir makes it, the thermal dissociation destroys it":
n(H2) = 3.598e-01/2.478e-03 = **1.45e+02 cm^-3**, a mixing ratio of 1.1e-06
by hydrogen nucleus.

### Judgment

**H2 can stand at its chemical root well above the molecular front, and that
root is set locally and not by transport.** At r = 1.6 R_p the transport
divergence of the carrier is 5.5e-09 cm^-3 s^-1 against chemical terms of
3.6e-01, so the H2 of that cell is made and destroyed where it stands: it is
a photochemical steady state of the hot gas, not molecular hydrogen advected
out of the layer below and slowly destroyed. Its level, 1.1e-06 by nucleus,
is about five decades above the thermochemical (detailed-balance) value the
`local` seed carries at the same (p, T), because the ion chain is a source
the thermochemistry does not have; the same is true at 1.13, 1.20 and 1.36
R_p, where the row demands 76, 5.4e+03 and 7.0e+04 times what is carried.
The practical consequence is that the H2 content of a hot wind is tied to the
metastable helium population -- the same population the He 10830 line
measures -- and not to how much H2 survives the front. **What this does NOT
establish**: the 145 cm^-3 is the root of the carrier row with H2+, HeH+ and
the ion stages frozen at the sweep's values, so it is a root of the
alternation's carrier half and not of the fully coupled system, which will
move when the sweep re-solves those ions at the higher H2; the branching
that carries the chain (the 10 percent associative branch of He(2^3S) + H,
Garcia Munoz 2025) and the two rate coefficients R19 and R9 are taken from
their sources and were not re-derived here; and no published model of a hot
wind has been compared against this H2 level. Those are the next questions,
and they are physics questions.

## 11. The movement bound and the progress control (approved follow-up)

Section 6.2 named two solver-control questions and left them for a decision.
Both were taken up; one of the two proposals was tried, MEASURED and
WITHDRAWN, and the memo records that because the measurement is the reason.

### 11.1 The e-fold clause: tried, measured, removed

The proposal was to admit every cell a change of at least a factor e of its
own value on top of the absolute bound, so that a wind cell could not be
frozen by a bound written on the base's abundance: allowed change =
max(trust x_ref, an e-fold of the local value). It was implemented and
measured, and the measurement refutes the premise it rested on.

**The e-fold is looser exactly where x_j is LARGE.** The absolute bound
admits an absolute change of trust x_ref = 1.9e-03 anywhere, and an e-fold
of the local value exceeds that wherever x_j > trust x_ref/(e - 1), about
1.1e-03 of the column maximum. So the clause does nothing for the wind cells
it was meant to free -- their whole required change, of order 1e-06 in
absolute units, was already far inside the absolute bound and was never what
the bound held -- and it frees the base and the front instead, which is the
opposite of what a bound on the wind's linear response is for. MEASURED on
the `src/tests/carrier_retry` column: with the clause in, the pass returns
the same advance at trust = 1e-2, 1e-3 and 1e-4 (drift 6.31969e-01 in 6 kept
steps at all three) and 3.0200e-01 in 5 steps at trust = 0, so **the setting
was inert over three decades**; with it out, the same column gives
9.3814e-03 in 2 steps at 1e-2 and 6.7222e-04 in 1 step at 1e-3, which is the
bound behaving as a bound.

The clause was removed. The bound is the absolute one it always was, and
`carrier_worst_composition_change` (new) returns the cell and carrier that
attain the largest change beside it, so that a pass refused on the bound
names the cell it was refused on. `src/tests/carrier_retry`'s movement-bound
rows are back to the contract they stated, and the reasoning above is at the
declaration of that routine so the experiment is not repeated.

### 11.2 A pass that kept no step does not tighten the bound (kept)

The outer loop's response to a pass that made no progress was to halve the
movement bound. That is the prescription for a carrier that moved and did
not help; for a carrier that kept no step at all the bound is the only thing
standing between the relaxation and a step, and halving it tightens exactly
that. `EXHALE_main.f90` now carries the previous pass's kept-step count and
leaves the bound where it is when that count is zero, saying so in the pass
report. The pass is still counted as one without a fall, so the loop still
ends on `outer_no_fall_max` consecutive such passes.

### 11.3 What a refused pass now says

Three reports were added, all printed only when the relaxation ends on the
reason they describe, so that the gate can be identified from the run log
instead of inferred:

- the chemistry closure's reason and the ledger counts behind it (section 4);
- the cell, the carrier, the absolute change and the entry value that
  attained the movement bound;
- the cell, the carrier and the row imbalance on which the transport
  operator could not cover the interval (`carrier_exhausted_record`).

Under `EXHALE_CARRIER_DEBUG=1` the relaxation also reports, once per pass,
the kept step count, the drift, and for three probe radii in the wind
(r = 1.20, 1.36 and 1.60 R_p) the H2 fraction the cell entered the pass
with, the one it leaves with, their ratio, and the physical interval the
kept steps covered there -- the last being what says whether a wind cell's
advance is limited by the interval the relaxation takes or by the bound.

## 12. The three gates, and the time scales they act on

The relaxation can end for four reasons and three of them have been seen on
this case. What each one means, stated from the code and from the state:

| ending | what it means | where it was seen |
|---|---|---|
| the chemistry did not close | the sweep after the trial did not reach a temperature fixed point, or left a cell non-finite | every pass before section 5 |
| the movement bound | the largest absolute carrier change of the returned state exceeds `trust` x (the largest H2 mixing ratio of the entry state) | section 6.1, passes 2 to 5 |
| the interval was not covered | `carrier_transport_interval` could not accept a substep even after halving the requested interval `carrier_retry_max = 8` times, i.e. at 1/256 of it; the entry state is restored and `frac_done` says how much was covered | section 6.1, passes 1, 6, 7 |
| the step budget / the fixed point | the pass ran out of trials, or the carriers stopped moving | — |

**The interval the relaxation asks for carries no chemical time.** It is
`dt_code(j) = min(dr_j^2/(D_i + K_zz), dr_j/|v_j|)`, the cell's diffusive or
advective crossing time, whichever is shorter. MEASURED on the state of
section 2 (`dr_j` from the `r` column, `v` from `Hydro_ioniz.txt`, the H2
loss rate per molecule from the row dump of section 3):

| cell | r [R_p] | v [cm/s] | dr [cm] | dr/\|v\| [s] | chemical time 1/(loss rate) [s] | ratio |
|---|---|---|---|---|---|---|
| 194 | 1.1340 | 7.57e+01 | 2.66e+06 | 3.51e+04 | 4.84e+02 | 73 |
| 218 | 1.2042 | 2.05e+02 | 4.04e+06 | 1.97e+04 | 2.58e+02 | 76 |
| 250 | 1.3578 | 6.04e+02 | 7.06e+06 | 1.17e+04 | 2.33e+02 | 50 |
| 280 | 1.6044 | 1.43e+03 | 1.19e+07 | 8.32e+03 | 4.04e+02 | 21 |

(the diffusive time is longer than the advective one at these cells --
dr^2/K_zz = 1.4e+05 s at cell 280 with K_zz = 1e9 -- so `dt_code` is the
advective crossing time there.) So every trial asks the transport operator
to integrate 21 to 76 chemical times of the H2 row in one interval, and the
operator may subdivide it 8 times, down to 1/256, which at cell 280 is 32 s
and is below the chemical time. The substep floor is therefore not itself
the obstruction, and which of the three gates sets the pace is a question
the run answers rather than the arithmetic.

## 13. The third gate, measured: one cell at r = 7 R_p refuses the column

`LHS1140b/models/.L7e/meas/HeH2.13`, the case on `EXHALE_L7e.x` with
`EXHALE_CARRIER_DEBUG=1`. Outer pass 1 (763 s, hydro `info = 0`, mass
7.20e-07, momentum 5.35e-11, energy 6.42e-06):

```
(carrier relaxation) kept 0 step(s); drift 0.000E+00
  probe cell 217 r 1.2007: x(H2) 3.08084E-11 -> 3.08084E-11, ratio 1.000, kept interval 0 s
  probe cell 250 r 1.3578: x(H2) 2.05927E-12 -> 2.05927E-12, ratio 1.000, kept interval 0 s
  probe cell 280 r 1.6044: x(H2) 1.27813E-12 -> 1.27813E-12, ratio 1.000, kept interval 0 s
carrier relaxation ended on the shortest admissible trial was not covered
  the interval was refused on carrier H2 at cell 412, r 7.0456:
  |res| over its physical terms 6.10E-05, physical terms 1.17E+02
```

Outer pass 2 (185 s, hydro `info = 0`) refuses on the OTHER gate and at the
neighbouring cell:

```
(carrier relaxation) kept 0 step(s); drift 0.000E+00
  probe cell 217 r 1.2007: x(H2) 2.85858E-11 -> 2.85858E-11, ratio 1.000, kept interval 0 s
  probe cell 250 r 1.3578: x(H2) 1.91142E-12 -> 1.91142E-12, ratio 1.000, kept interval 0 s
  probe cell 280 r 1.6044: x(H2) 1.18653E-12 -> 1.18653E-12, ratio 1.000, kept interval 0 s
carrier relaxation ended on the movement bound, with the last step inside it
  the bound was refused on carrier H2 at cell 411, r 6.9412:
  absolute change 1.53E-03, from an entry value of 6.58E-02
```

so the movement bound IS reached -- by a cell at 6.94 R_p moving 2.3 percent
of its own value, which is the operator correctly taking the hump down --
and not by the base, not by the front, and not by the cells the
certification refuses on. (The probe values move between the two passes only
because the sweep and the hydrodynamic solve move n_tot; their ratio across
each relaxation is 1.000 to every digit printed.)

**So it is neither of the two the question was framed between.** It is not
the movement bound -- the bound is never reached, because no trial is ever
kept -- and it is not the pass budget: the wind cells the certification
refuses on do not move by one part in 1e+16, because the relaxation keeps no
step at all. **The gate is the transport operator's returned-state
acceptance, `carrier_accept_tol = newton_floor = 1.0e-08` applied as a
maximum over every cell, and the cell that refuses it stands at r = 7.05
R_p**, five radii beyond anything the certification is about.

### 13.1 What cell 412 carries, and why no substep can solve it

MEASURED on the same state (`output/Ion_species.txt` and the row dump of
section 3):

| cell | r [R_p] | T [K] | n(H2) | x(H2) | n(H I) | n(H II) | 2n(H2)/n_H |
|---|---|---|---|---|---|---|---|
| 411 | 6.941 | 1413 | 7.71e+04 | 1.43e-01 | 2.55e+01 | 1.96e+05 | 0.440 |
| 412 | 7.046 | 1397 | 7.21e+04 | 1.42e-01 | 2.70e+01 | 1.85e+05 | 0.438 |
| 413 | 7.152 | 1382 | 6.75e+04 | 1.42e-01 | 3.15e+01 | 1.75e+05 | 0.436 |

**Forty-four percent of the hydrogen nuclei at 7 R_p are in H2, while atomic
H stands four decades BELOW H2 and more than half the hydrogen is ionized.**
A gas cannot hold that: H2 is made from atomic H, and the row says so --
production 2.08e-04 against a loss of 7.47e-02 cm^-3 s^-1, a factor 359, of
which 2.32e-02 is photodissociation and photoionization by the stellar
field. The root of that row is n(H2) = 2.0e+02 cm^-3; the state carries
3.6e+02 times it.

And it is not one cell. The same dump over the whole column:

| cell | r [R_p] | T [K] | x(H2) | production | loss | root n(H2) | carried/root |
|---|---|---|---|---|---|---|---|
| 300 | 1.857 | 5364 | 4.98e-12 | 1.34e-01 | 9.32e-07 | 9.73e+01 | 7.0e-06 |
| 340 | 2.723 | 3638 | 1.42e-10 | 1.27e-02 | 9.59e-07 | 4.50e+01 | 7.6e-05 |
| 380 | 4.460 | 2145 | 5.68e-07 | 5.20e-04 | 4.77e-05 | 1.99e+01 | 9.2e-02 |
| 400 | 5.904 | 1568 | 1.16e-03 | 7.74e-05 | 6.19e-03 | 1.60e+01 | **8.0e+01** |
| 412 | 7.046 | 1397 | 1.42e-01 | 2.08e-04 | 7.47e-02 | 2.01e+02 | **3.6e+02** |
| 430 | 9.274 | 1155 | 1.63e-01 | 6.70e-05 | 1.86e-02 | 1.06e+02 | **2.8e+02** |
| 450 | 12.73 | 944 | 2.10e-01 | 3.13e-05 | 5.40e-03 | 7.85e+01 | **1.7e+02** |
| 470 | 17.62 | 733 | 1.72e-01 | 8.79e-06 | 1.53e-03 | 2.58e+01 | **1.7e+02** |
| 490 | 24.55 | 581 | 1.29e-01 | 1.88e-06 | 5.07e-04 | 5.61e+00 | **2.7e+02** |

**The state carries a molecular hump over the whole outer wind**, from about
5.9 R_p to the outer boundary at 29 R_p: 13 to 21 percent of the gas in H2,
80 to 360 times the root of its own row, in a gas at 580 to 1570 K that the
star is irradiating and that is half ionized. Inside 4.5 R_p the same
column is the opposite -- the H2 is 1e-6 to 1e-1 of its root, which is
section 10 -- so the carrier is wrong in both directions and by many decades
at both ends.

### 13.2 Where the hump comes from

`SEED_X2 = local` gives every cell the thermochemical equilibrium
q_H2(p, T) of `q_h2_equilibrium`. That fit is a statement about a gas whose
H2 is set by the balance of the three-body association against the thermal
dissociation, and it is the right statement at the base. In the outer wind
it is not: the gas is at 580 to 1570 K only because it expanded, it is half
ionized, and the stellar field dissociates H2 there -- the fit knows none of
that, and it returns a large H2 fraction because the temperature alone is
low. The seed therefore asserts an equilibrium that does not hold over
roughly a fifth of the column.

### 13.3 Judgment

**Not a budget, not a bound: the relaxation is blocked by a seeded
composition that is not a possible state of the gas, through an acceptance
that is a maximum over cells.** Two consequences follow, and they are
separate:

1. **The acceptance is doing its job.** A cell 360 times from the root of
   its own row, with atomic hydrogen four decades below the molecule it is
   supposed to make, is exactly the state a returned-state acceptance should
   refuse. Loosening `carrier_accept_tol`, or exempting the outer wind from
   it, would hide the seed's error rather than remove it; the code's own
   principle (an absent carrier is REPORTED and does not gate, but a carrier
   holding 14 percent of the gas is not absent) refuses that reading.
2. **The seed is the thing to fix.** The `local` extension must not impose
   the thermochemical H2 where the thermochemistry is not the balance that
   holds. The quantity that does hold there is the one the row already
   computes: the cell's own production over its own loss rate, which in the
   hump is 80 to 360 times smaller and at the base is the fit itself.

**The proposal is therefore to give the molecular seed the network's own
balance wherever it is the smaller of the two**, with the fit kept where it
is the statement that holds, and to record the crossover radius each case
lands on. It touches which state every molecular case starts from, so it is
proposed here and not executed.

## 14. The seed rewritten, and the question it opens at the base

Approved and implemented: `SEED_X2 = local` now means **the smaller of the
thermochemical fit q_H2(p, T) and the root of that cell's own H2 carrier
row** (production over loss rate, with the photodissociation, the
photoionization and the ion channels in it). The root is taken from
`carrier_h2_chemical_root`, which calls `carrier_source` -- the routine that
assembles the row -- so there is no second statement of the chemistry. The
molecular ions that carry the production do not exist in an atomic state, so
the sequence is: apply the fit, take ONE equilibrium sweep, read the root,
throw that sweep's composition away, and apply the revised partition to the
ATOMIC state, so nothing the sweep decided about the ionization reaches the
file. The crossover radius and the number of cells the root governs are on
the run log and in a `# molecular_partition_local:` line of both state
files.

### 14.1 What the two statements actually say

MEASURED on `.L7d/src_atomic_HeH2.13` with He/H = 2.13
(`EXHALE_CARRIER_DEBUG=1` prints this table):

| cell | r [R_p] | T [K] | x2 fit | x2 root | adopted |
|---|---|---|---|---|---|
| 1 | 1.0002 | 418 | **1.000** | 8.915e-02 | 8.915e-02 |
| 25 | 1.0048 | 632 | **1.000** | 1.374e-01 | 1.374e-01 |
| 100 | 1.0250 | 1230 | **1.000** | 4.773e-02 | 4.773e-02 |
| 150 | 1.0615 | 1899 | 1.727e-01 | 6.139e-02 | 6.139e-02 |
| 200 | 1.1489 | 4030 | 8.350e-09 | 2.353e-06 | 8.350e-09 |
| 250 | 1.3578 | 5811 | 2.555e-11 | 2.130e-06 | 2.555e-11 |
| 300 | 1.8571 | 4873 | 2.056e-11 | 4.056e-06 | 2.056e-11 |
| 375 | 4.1715 | 2152 | 4.588e-07 | 3.043e-05 | 4.588e-07 |
| 400 | 5.9043 | 1515 | 4.060e-03 | 8.127e-05 | 8.127e-05 |
| 425 | 8.5835 | 1057 | **1.000** | 5.890e-03 | 5.890e-03 |
| 450 | 12.726 | 748 | **1.000** | 7.303e-03 | 7.303e-03 |
| 500 | 29.031 | 459 | **1.000** | 1.493e-03 | 1.493e-03 |

**The fit is saturated at x2 = 1 -- "every hydrogen nucleus is in H2" -- over
BOTH ends of the column**, cells 1 to about 125 and 425 to 500, and carries
information only in the 1.06 to 5.9 R_p band between them. The old `local`
seed therefore was not "each cell's chemical equilibrium" at either end; it
was the element-ratio ceiling, capped by whatever neutral hydrogen the
atomic state had (L7d reports 235 cells so capped, this seed 0).

**The outer hump is removed**, which is what the item measured and what the
change was made for: at cell 412 the seeded x2 falls from 4.4e-01 to
5.2e-03, a factor 85, and the whole 5.9 to 29 R_p region now carries the
root of its own row.

### 14.2 The base layer: a finding, and not a rule

Applied everywhere, the min also revises the cells nearest the base, and
there it disagrees with the run's own input: `base.inp` states
`q_H2_base = 0.190`, which for He/H = 2.13 is x2 = 0.99998, and the row's own
root says 4.8e-02 to 1.4e-01 through the cells up to r = 1.0533 R_p (cells 1
to 142, where the thermochemical fit falls under `q_H2_base` and thermal
dissociation takes over).

**Carrying the handoff through that layer was tried and withdrawn.** The
reason it is wrong is about the SOLUTION and not about the seed: the handoff
partition is imposed during a run on the GHOST cells alone
(`ionization_equilibrium`, the `base_h2_composition_imposed` block, `j <= 0`),
so cells 1 to 142 are free and the fixed point of the solve puts their x2 at
the network's own root, not at the handoff's. A seed that imposes the
handoff there is not seeding the fixed point; it adds a transient the
relaxation must then undo -- a 50 percent change of the particle count of
those cells -- and the movement bound takes that at one percent a pass,
holding the whole column while it runs. MEASURED: on that seed the bound was
refused at cell 143 on every pass of `.L7e/fixed3` and `.L7e/fixed4` alike.

**What is left is a finding.** Inside that layer the handoff stands **7.24 to
21.00 times** the root of the wind network's own H2 row: the lower-atmosphere
model and the wind network do not agree about a layer they both describe, and
one of the two is wrong. The seed reports the ratio in its log and in the
`# molecular_partition_local:` line of both state files, and plan item L7f
asks which chemistry is right -- the answer decides something this item does
not, namely whether the base boundary condition should keep imposing the
handoff on the ghosts alone or over the layer.

## 15. The seed WAS the block: the carrier moves

MEASURED, outer pass 1 of `.L7e/fixed2/HeH2.13` -- the seed with the outer
hump removed -- against every pass of every run before it, which kept ZERO
steps and moved the wind cells by nothing:

```
(carrier relaxation) kept 3 step(s); drift 9.966E-03
  probe cell 217 r 1.2007: x(H2) 3.27136E-11 -> 1.81034E-08, ratio 5.534E+02, kept interval 6.021E+01 s
  probe cell 250 r 1.3578: x(H2) 2.20821E-12 -> 1.74155E-08, ratio 7.887E+03, kept interval 2.838E+01 s
  probe cell 280 r 1.6044: x(H2) 1.36205E-12 -> 1.33330E-08, ratio 9.789E+03, kept interval 2.171E+01 s
(EXHALE_main) outer pass 1: hydro info = 0, mass 3.37E-07, momentum 7.62E-11, energy 1.41E-06
carrier relaxation ended on the movement bound; drift 9.97E-03 in 3 transport steps
  the bound was refused on carrier H2 at cell 3, r 1.0006: absolute change 7.45E-05,
  from an entry value of 6.82E-03
```

**The wind cells rise by 2.7 to 4.0 decades in one pass**, toward the root
section 10 measured for them, in a kept interval of 22 to 60 s against their
own chemical time of about 400 s -- so the advance is set by the chemistry
and by how much of it one trial may cover, not by a starvation of interval.
The hydrodynamic rows are unaffected (mass 3.4e-07, energy 1.4e-06, better
than the 7.2e-07 and 6.4e-06 of every earlier run). **So the answer to the
question the item was opened on is: the carrier relaxation did not reach its
fixed point because the seed asserted a molecular hump the operator could not
integrate, and with that gone it moves.**

Outer pass 2 gives the rate the question was framed in terms of:

```
(carrier relaxation) kept 4 step(s); drift 9.863E-03
  probe cell 217: x(H2) 1.73740E-08 -> 3.54071E-08, ratio 2.038, kept interval 6.552E+01 s
  probe cell 250: x(H2) 1.68065E-08 -> 3.39889E-08, ratio 2.022, kept interval 3.108E+01 s
  probe cell 280: x(H2) 1.28388E-08 -> 2.63108E-08, ratio 2.049, kept interval 2.376E+01 s
outer pass 2: hydro info = 0, worst gated species row 9.97E-01 at cell 218
```

with the gated row leaving 1.000 for the first time (9.97e-01), and pass 3
adds x(H2) 2.579e-08 -> 4.002e-08 at cell 280, a ratio of 1.55 and the row at
9.93e-01.

**THE RATE IS A CONSTANT ABSOLUTE INCREMENT AND NOT A CONSTANT RATIO**, and
the distinction decides whether the cap is enough. At cell 280 the three
passes give 1.28e-08 -> 2.63e-08 -> 4.00e-08: the ratios fall, 2.05 then
1.52, while the increments are 1.35e-08 and 1.37e-08, the same number. That
is the movement bound: it admits a fixed absolute change of `trust` x_ref per
pass, and the wind cells ride whatever trial length the bound leaves. At
1.37e-08 a pass, reaching the 1.1e-06 root of section 10 would take about
seventy more passes against a cap of twenty -- unless the cell that saturates
the bound stops moving first.

**And that cell is the base.** The bound is refused at cell 3, r = 1.0006, on
every pass, at an absolute change of 7.4e-05 from an entry value of 6.8e-03
-- the allowance `trust` x_ref with x_ref the column's largest entry H2
fraction. The wind's own change, 1.4e-08, is four decades below the
allowance and is never what the bound measures; it is a passenger on the trial
length the base's relaxation leaves. (This run's base is the ROOT's and is
therefore relaxing away from its seed; the delivered seed carries the handoff
there, which is what `fixed3` measures.) (This run's base is the ROOT's, because it is
the control of the base decision of section 14.2; the delivered seed carries
the handoff there.)

## 16. The chemistry is already implicit; the limiter is the bound

A backward-Euler treatment of the chemical source was proposed for the
carrier substep. **The code already does it, and more than the proposal
asks.** READ, `carrier_residual`: every row is

    nrho (fc - fc_old)/dt  +  transport(fc)  -  src(fc nrho)   =   0

with `src` evaluated at the TRIAL `fc`, which is backward Euler with an
implicit source; `solve_carriers` Newtons on it, builds the diagonal block
as `bb = ... - (d src/d n) nrho` by finite differences of `carrier_source`,
solves the block-tridiagonal system with `block_thomas` and line-searches
the step. The proposed `x_new = (x_old + dt P)/(1 + dt L)` is the
linearized, single-species special case of that, and dropping the coupling
between the carriers and the transport terms would make it weaker, not
stronger.

MEASURED, and this is the decisive number: outer pass 1 of `.L7e/fixed3`
lifts x(H2) at cell 280 from 1.30233e-12 to 1.88355e-09, **a factor 1446, in
a kept interval of 3.505 s**, against a local chemical time of about 400 s
(cell 250: a factor 1085 in 3.807 s; cell 217: 80.8 in 4.958 s). An explicit
source cannot move a species by more than about dt/tau, here 0.9 percent, in
that interval. The substep log says the same from the other side: `newton 2
resid 1.26E-13`, two iterations to 1e-13 per substep, and the interval's own
verdict is `the residual reached the absolute floor`.

**So the physical interval of a pass is short because the MOVEMENT BOUND
cuts the trial, not because the integrator cannot cross a chemical time.**
The bound is refused, on every pass, at a cell of the base or the front and
never at the cells the certification is about:

| run | base of the seed | bound refused at | absolute change | entry value |
|---|---|---|---|---|
| `.L7e/fixed2` | the root's | cell 3, r = 1.0006 | 7.4e-05 | 6.8e-03 |
| `.L7e/fixed3` | the handoff's | cell 143, r = 1.0543 | 5.97e-04 | 4.73e-03 |

and in `fixed2` the consequence is arithmetic: a constant absolute increment
of 1.37e-08 per pass at cell 280, which is the share of `trust` x_ref the
wind gets once the saturating cell has taken its own.

### 16.1 What would actually lift it, and it is not the integrator

The bound exists because the wind is held fixed while the carriers relax, so
that one pass may move the composition only as far as the wind's response to
it stays linear. **What the wind responds to is the particle count and the
mean molecular mass of the cell**, not the carrier's own fraction: an H2 at
1e-08 of the gas changes neither, an H2 at tens of percent changes both. A
bound written on the quantity the wind feels -- the relative change of
n_tot + n_e of each cell, equivalently of 1/mu -- is therefore the physical
statement of the bound's own purpose, and it is automatically loose in the
wind cells and tight at the base and the front, with no threshold beyond the
`trust` already anchored. That is the repair the e-fold clause of section
11.1 was reaching for and got wrong by writing it on the carrier fraction
instead.

**Proposed and not implemented**, because it is the same solver-control
surface the e-fold attempt was withdrawn from, and because `fixed3` may
still certify within its pass cap, in which case nothing needs changing.

## 17. The movement bound, rewritten on what the wind feels

Approved and implemented. The bound is now the relative change of each
cell's PARTICLE COUNT, `n_tot + n_e`, against the composition the pass was
entered with, at most `trust` -- the same threshold, on the quantity the
fixed wind responds to. `carrier_particle_count_change` (new) measures it and
carries the reasoning and the two measurements that retired the old one;
`EXHALE_CARRIER_BOUND_FRACTION=1` restores the carrier-fraction measure.

**Why this quantity.** The bound exists because the wind is held while the
carriers relax, so a pass may move the composition only as far as the wind's
response stays linear. What the wind responds to is the particle count and
the mean molecular mass -- they are what set the pressure at the conserved
thermal energy -- not how much a carrier moved relative to another cell's
abundance. An H2 at 1e-08 of the gas moves neither; an H2 at tens of percent
moves both. The bound is therefore loose exactly where the carrier cannot
affect the wind and tight exactly where it can, with no threshold beyond
`trust`.

**And it answers to its setting**, which the e-fold clause of section 11.1
did not. MEASURED on the `src/tests/carrier_retry` column:

| `trust` | particle-count change | carrier drift | kept steps |
|---|---|---|---|
| 1e-2 | 9.9238e-03 | 1.6874e-02 | 4 |
| 1e-3 | 7.7835e-04 | 1.3440e-03 | 1 |

both falling by about 12.7 when the setting falls by 10, against the e-fold
clause's inertness over three decades. At `trust` = 1e-2 the same column now
admits a carrier drift of 1.69e-02 where the old measure admitted 9.38e-03,
because the carrier is free to move where the wind cannot feel it. The
suite's movement-bound rows are restated on the new quantity and pass (143
assertions).

## 18. The control's record, and why it was stopped

`.L7e/fixed2` -- the seed whose base layer carries the row's own root, the
control of the base decision of section 14.2 -- was stopped by hand after
outer pass 8 of 20. Its record, which is what it was run for:

| pass | x(H2) at cell 280 | gated row | refusing entries | hydro mass / energy |
|---|---|---|---|---|
| 1 | 1.333e-08 (from 1.362e-12) | 1.000 | 2 | 3.37e-07 / 1.41e-06 |
| 2 | 2.631e-08 | 9.97e-01 | 2 | 3.35e-07 / 1.41e-06 |
| 3 | 4.002e-08 | 9.93e-01 | 2 | 3.34e-07 / 1.41e-06 |
| 4 | 5.524e-08 | 9.89e-01 | 2 | 3.35e-07 / 1.40e-06 |
| 5 | 7.174e-08 | 9.82e-01 | 2 | 3.37e-07 / 1.40e-06 |
| 6 | 9.012e-08 | 9.74e-01 | 2 | 3.35e-07 / 1.40e-06 |
| 7 | 1.108e-07 | 9.61e-01 | **1** | 3.32e-07 / 1.39e-06 |
| 8 | 1.346e-07 | 9.42e-01 | **1** | 3.31e-07 / 1.39e-06 |

so the row falls monotonically, the hydrodynamic rows do not move, and the
count of entries refusing the certification goes 5 at entry to 2 to 1 -- the
carrier row alone. The increments at cell 280 are 1.35, 1.37, 1.52, 1.65,
1.84, 2.08 and 2.38e-08: **the bound's allowance per pass, rising as the
cell that saturates it relaxes**. It was stopped because its trajectory is
established and the machine was carrying another session's work at a load of
200; the run that decides anything is `.L7e/fixed4`, on the new bound.

## 19. The new bound against the old, on one seed and one pass

`.L7e/fixed3` and `.L7e/fixed4` are the same case, the same delivered seed
and the same binary text; they differ in the movement bound alone. Outer
pass 1 of each:

| | old bound (carrier fraction) | new bound (particle count) |
|---|---|---|
| x(H2) at cell 250 | 2.099e-12 -> 2.278e-09, x1085 | 2.099e-12 -> **1.452e-08, x6920** |
| x(H2) at cell 280 | 1.302e-12 -> 1.884e-09, x1446 | 1.302e-12 -> **1.244e-08, x9549** |
| kept interval at cell 280 | 3.505 s | **24.17 s** |
| kept transport steps | 2 | 2 |
| carrier drift of the pass | 9.81e-03 | **5.86e-02** |
| refused at | cell 143, carrier change 5.97e-04 of an entry 4.73e-03 | cell 143, particle-count change 1.01e-02 of an entry 4.41e-03 |
| hydrodynamic mass / energy | 7.38e-07 / 6.54e-06 | 7.38e-07 / 6.54e-06 |

**The same two kept steps carry the wind cells 6.6 times further, over a
physical interval 6.9 times longer**, and the pass moves the carriers by
5.86e-02 against 9.81e-03 -- six times as much -- while the hydrodynamic rows
are identical to the last digit printed, which is the change touching the
carrier half of the alternation and nothing else. The bound is still what
ends the pass and it is still cell 143, the joint between the seed's handoff
base layer and the root above it; but it now refuses on a cell whose particle
count really does move by one percent, which is a statement the wind can feel,
and not on a carrier fraction measured against another cell's abundance.

## 20. What the new bound does over eleven passes

`.L7e/fixed4`, the delivered seed on the particle-count bound. MEASURED:

| pass | x(H2) at cell 280 | ratio | kept interval at 280 | worst gated row | at cell |
|---|---|---|---|---|---|
| 1 | 1.244e-08 | x9549 | 24.2 s | 1.000 | 279 |
| 2 | 2.491e-08 | 2.14 | 29.1 s | 9.95e-01 | 217 |
| 3 | 3.906e-08 | 1.63 | 36.8 s | 9.89e-01 | 217 |
| 4 | 5.692e-08 | 1.49 | 53.0 s | 9.80e-01 | 217 |
| 5 | 8.013e-08 | 1.42 | 80.0 s | 9.63e-01 | 217 |
| 6 | 1.129e-07 | 1.42 | 149 s | 9.31e-01 | 217 |
| 7 | 1.429e-07 | 1.27 | 220 s | 8.67e-01 | **366** |
| 8 | 1.677e-07 | 1.18 | 370 s | 8.03e-01 | 367 |
| 9 | 1.849e-07 | 1.11 | 529 s | 7.21e-01 | 369 |
| 10 | 1.929e-07 | 1.045 | 811 s | 6.31e-01 | 370 |
| 11 | **1.944e-07** | **1.010** | **1208 s** | **5.31e-01** | 371 |

Two things happen together and they are the same thing. **The cell the
certification used to refuse on has converged**: the ratio at cell 280 falls
to 1.010 and its kept interval reaches 1208 s, three of its own chemical
times, so the bound is no longer what stops it there -- the cell has reached
the root of its own row. **And the worst row moves outward**, from cell 217
to 366 at pass 7 and then one cell a pass to 371: the relaxation is sweeping
a front through the column, leaving converged cells behind it. The gated row
falls 9.31e-01, 8.67, 8.03, 7.21, 6.31, 5.31e-01, by a step that grows
every pass.

Against the same case on the old bound at the same passes, `.L7e/fixed3`
stands at cell 280 x(H2) = 1.99e-08 with its worst row still at cell 217 and
9.91e-01: a factor 9.8 behind in the carrier and not moving its front at
all.

The hydrodynamic rows are unmoved by any of it (mass 7.1e-07 to 7.4e-07,
energy 6.4e-06), except at passes 8 and 9 where they fall to 5.8e-08 and
5.4e-07 -- an order of magnitude -- and come back, which is the composition
passing through consistency with the wind.

## 21. The carrier row is closed; what refuses is the base cell's mass row

`.L7e/fixed5/HeH2.13` -- the delivered seed (min of the fit and the root in
every cell) on the particle-count bound, 40 outer passes, the molecular
`EXHALE_OUTER_PASSES` raised to the atomic cases' 40. Its certification:

| entry | measure | tolerance | verdict |
|---|---|---|---|
| carrier balance H2 | **1.93e-08** | 1.0e-05 | within |
| elemental transport He/H | 8.5e-09 | 1.0e-05 | within |
| level balance He 2^3S, eliminated-species closure | within | | within |
| hydrodynamic momentum | 2.848e-11 | 1.0e-08 | within |
| hydrodynamic energy | 6.937e-07 | 1.0e-06 | within |
| **hydrodynamic mass, cell 1** | **1.555e-07** | **3.7e-09** | **ABOVE, 41.9x** |

**The item's own target is met.** The H2 carrier row, which opened this memo
at 1.000 and was the reason L7e exists, stands at 1.93e-08, five hundred
times inside its tolerance. What refuses the state is one hydrodynamic row
in one cell.

### 21.1 It is stalled, not descending

MEASURED, the cell-1 mass row over the forty passes:

| pass | 1 | 5 | 9 | 13 | 17 | 21 | 25 | 29 | 33 | 37 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mass row | 2.86e-07 | 1.56e-07 | 1.55e-07 | 1.55e-07 | 1.55e-07 | 1.56e-07 | 1.57e-07 | 1.57e-07 | 1.57e-07 | 1.57e-07 | 1.56e-07 |
| momentum | 5.23e-11 | 2.84e-11 | 2.85e-11 | 2.85e-11 | 2.85e-11 | 2.85e-11 | 2.84e-11 | 2.84e-11 | 2.84e-11 | 2.85e-11 | 2.85e-11 |
| energy | 1.27e-06 | 6.97e-07 | 6.96e-07 | 6.95e-07 | 6.95e-07 | 6.95e-07 | 6.91e-07 | 6.92e-07 | 6.92e-07 | 6.94e-07 | 6.94e-07 |

It falls once, over the first four passes, and then **sits at 1.55 to
1.57e-07 for thirty-five passes**. So a continuation is not what it needs:
the alternation is at its fixed point and this row is where the fixed point
stands.

### 21.2 What the row is made of, and why it cannot close

The mass row of cell 1 is the difference of two face mass fluxes, the base
face and the 1-2 face. The written state (`output/`):

| | r [R_p] | rho [mH/cm3] | v [cm/s] | x2 = 2n(H2)/n_H |
|---|---|---|---|---|
| ghost | 0.99981 | 2.938e+13 | -0.8174 | **0.99998** |
| ghost | 1.00000 | 2.903e+13 | -0.8269 | **0.99998** |
| cell 1 | 1.00019 | **2.161e+13** | -0.1374 | **0.3315** |
| cell 2 | 1.00039 | 1.787e+13 | +0.1572 | 0.3317 |

**The two sides of the base face are not the same gas.** The boundary
condition imposes the lower-atmosphere handoff on the GHOSTS -- x2 =
0.99998, the `q_H2_base` = 0.19011 of `base.inp` -- while the first physical
cell, which the boundary does not impose and the solve is free to move,
stands at x2 = 0.3315, the network's own root. At He/H = 2.13 that is 2.63
against 2.964 particles per hydrogen nucleus, so the mean molecular mass
differs by **12.8 percent across one face**, and the density steps 26 percent
from 2.903e+13 to 2.161e+13 while the velocity steps from -0.827 to -0.137.
The inflow the boundary specifies and the column that receives it therefore
carry different mass fluxes at the same face, and no relaxation of the
composition can make them agree while the two sides state different
chemistries.

### 21.3 The atomic case says the same thing from the other side

`atomic_scalar_gj1132_kzz1e9/HeH2.13`, the same planet, the same He/H, the
state this seed was built from, CERTIFIED: its mass row reads **1.487e-09 at
cell 4**, and the cell that decides its verdict is 178 at 1.351e-11 against
2.2e-11 -- within. Its base ghosts step from rho = 5.217e+13 to 4.997e+13 at
cell 1, **4.2 percent**, with no composition discontinuity at all because the
gas is atomic on both sides.

So the molecular case's cell-1 mass row is **105 times** the atomic case's
worst mass row anywhere, and the difference is exactly the composition step
the molecular base boundary condition creates.

### 21.4 Judgment

**The L7f disagreement is not an abstract one: it is this row.** The base
boundary condition imposes the lower-atmosphere model's H2 on the ghost
cells; the wind network's own chemistry puts the first physical cell at a
third of it; and the mass row of the face between them is the difference.
Until one of the two chemistries is settled -- which is L7f -- the molecular
column cannot close its base cell, and the choice that follows is the one
L7f already names: whether the base boundary condition should keep imposing
the handoff on the ghosts alone, over the first cells as well, or not at all.
Nothing in the carrier operator, the seed or the movement bound can remove
it, and this item's work on those three is what made it the only thing left.

## 22. Where L7e stands

**Closed, on its own question.** The H2 carrier balance that opened the item
at exactly 1.000 stands at 1.93e-08 against a tolerance of 1.0e-05 on the
state `.L7e/fixed5/HeH2.13` writes. Four things were wrong and each was
found by measurement:

1. the chemistry closure refused every carrier trial on a count that
   describes the root search and not the state (section 5.1);
2. its tolerance sat below the equilibrium sweep's own noise (5.2);
3. the molecular seed asserted a thermochemical equilibrium over the whole
   outer wind, 80 to 360 times the root of the row it seeds (13, 14);
4. the movement bound was written on the carrier fraction against the
   column maximum, so it was refused on base and front cells and never on
   the cells the certification is about (16, 17).

Two proposals were tried and withdrawn on measurement rather than argument
-- the e-fold clause on the bound (11.1) and the handoff through the base
layer (14.2) -- and one was refuted before it was written, the implicit
chemistry the code already had (16).

**Open, and named:** the base cell's mass row, which is the lower-atmosphere
handoff and the wind network disagreeing about the same layer across one
face (21), and which plan item L7f asks the physics of.

Tests on the delivered binary: `carrier_retry` 143, `certification` 84,
`steady_species_rows` 195, `molecular_seed` 24, and
`carrier_reference_scales`, `carrier_constraint_attribution`,
`carrier_returned_state_acceptance` with no failures. Regression
`mol_carrier`, `mol_base_handoff`, `mol_diffusion`: **REGRESSION PASS**,
worst relative movement 2.548e-05, the same numbers the control binary gave,
so no golden is refreshed.

**Not started**: the eight remaining molecular cases, which were gated on a
certification this state does not reach.

### 22.1 The two controls, finished

Both ran to their end and both are on disc; neither is running any longer.

| run | seed | bound | ended | refusing entries |
|---|---|---|---|---|
| `.L7e/fixed3` | handoff through the base layer | carrier fraction | its cap, then the continuation | **4**: mass 1.025e-06 at cell 1, energy 8.411e-06 at cell 2, **carrier balance H2 4.108e-02 at cell 500**, He/H 3.535e-05 at cell 281 |
| `.L7e/fixed4` | handoff through the base layer | particle count | pass 34, the outer loop refusing it ("the joint distance has not fallen in 3 consecutive passes") | **5**: mass 1.997e+00 at cell 139, momentum 7.906e-02 at cell 129, energy 1.101e+00 at cell 162, **carrier balance H2 4.081e-01 at cell 253** |
| `.L7e/fixed5` | min in every cell | particle count | 40 passes | **1**: mass 1.555e-07 at cell 1 |

So the seed decision and the bound decision are both visible in the outcome
and not only in the rate: with the handoff carried through the base layer the
carrier row ends at 4.1e-02 or 4.1e-01 and three or four hydrodynamic rows
refuse with it, and `fixed4` is refused by the outer loop itself for not
approaching a joint fixed point. The delivered pair -- the min everywhere and
the particle-count bound -- is the only one of the three whose carrier row
closes, and it leaves exactly one entry, the base cell's mass row of section
21.

---

## 23. Review of 2026-09-15, items 5 to 7: the boundary derivatives and what actually holds the two wellmixed columns

### 23.1 The carrier matrix was missing two boundary derivatives and carried one wrong factor

The review asked for the first-order terms of the carrier matrix to be
measured, not argued. A new test suite,
`src/tests/carrier_boundary_jacobian/`, forms a central difference of
`carrier_advective_divergence` in one cell's carrier density and compares it
against the coefficient the assembly writes into `aa`, `bb`, `cc`. The test
state is built so that the mass fraction is exactly constant
(`fc(j,ic) = yconst*msum(j)/carrier_mass_amu(ic)`, outer ghosts flattened to
`msum(N)`), which makes the limited reconstruction exact and the measured
derivative the analytic one to round-off.

Three defects were found and fixed in `diffusive_photochemistry.f90`:

| # | What was missing | Where |
|---|---|---|
| 1 | the outward face's derivative at the last cell: with `Frho(N) < 0` the donor is the outer ghost and nothing was written, so `bb(N)` carried no outward term | carrier Jacobian, advective block |
| 2 | the inward face's derivative at the first cell when the base composition is NOT imposed: the base face's donor is then cell 1 itself and its derivative belongs on the diagonal | same |
| 3 | a spurious division by `carrier_mass_amu(ic)`: the carrier mass cancels between the face flux and the cell's own count, so the coded coefficient was a factor `m_c` too small (`2.0088` measured for H2 before the fix) | same |

All ten assertions of the new suite pass after the fix. The measured
coefficient is `advj(j)*msum(j)/msum(donor)`.

This is a statement about the Jacobian, and only that. It is **not** claimed
to be the cause of any refusal: the two columns below were measured after
the fix and still refuse.

### 23.2 The row decomposition of the two wellmixed columns

`EXHALE_CARRIER_ROW_TERMS=1` now writes, beside the row's own terms and the
H2 chemistry reaction by reaction, the TWO FACES of every row separately —
the diffusive and the advective flux of each face with its sign, and the
mass flux of each face — so that a reader can see which face carries what.
The record is written inside the certification's own isolated evaluation, so
the file and the verdict describe one state.

Both states were re-read from their written profiles and re-measured with
`EXHALE_L7f.x` (md5 `defab12df01629792c99d0d5b6d47228`); the reload
reproduces the verdict exactly (`2.462E-02` at cell 500 and `2.477E-02` at
cell 306, the same numbers the campaign runs refused on), so the
decomposition below is of the refused state and not of a neighbour of it.

The row is `transport - (production - loss)` and the certification scale is
the sum of the row's own terms with the faces counted separately,
`|F_dif,in| + |F_dif,out| + |F_adv,in| + |F_adv,out| + |net source| + floor`.
Every rate is cm^-3 s^-1.

**`molecular_scalar_gj1132_wellmixed/HeH0.083`** (24 outer passes, refused at
cell 500):

| | cell 500 (worst row) | cell 499 (neighbour) | cell 231 (worst drift) |
|---|---|---|---|
| r [R_p] | 29.031 | 28.549 | 1.2565 |
| T [K] | 285.3 | 286.5 | 1976.6 |
| x(H2) | 7.807e-03 | 8.094e-03 | 2.872e-01 |
| residual | +6.747e-04 | +6.864e-04 | -3.908e+02 |
| row scale | 2.741e-02 | 3.161e-02 | 5.910e+04 |
| **measure** | **2.462e-02** | 2.171e-02 | 6.612e-03 |
| diffusive, inner face | +9.786e-04 | +2.189e-03 | -7.849e+03 |
| diffusive, outer face | 0 (domain edge) | -1.021e-03 | +7.526e+03 |
| advective, inner face | -1.337e-02 | -1.444e-02 | -2.190e+04 |
| advective, outer face | +1.304e-02 | +1.394e-02 | +2.152e+04 |
| mass flux, inner / outer face | +1.418e-10 / +1.372e-10 | +1.467e-10 / +1.418e-10 | +7.471e-08 / +7.418e-08 |
| production | 7.637e-08 | 8.395e-08 | 1.248e+02 |
| loss | 1.965e-05 | 2.106e-05 | 4.396e+02 |
| largest loss channels | LW photodissociation 1.07e-05, H2 photoionization 8.90e-06 | 1.15e-05, 9.54e-06 | H+ charge transfer 1.42e+02, H2+ + H2 1.27e+02, photoionization 8.76e+01, He+ + H2 7.20e+01 |

At cell 500 the chemistry is **3 % of the residual** (1.96e-05 against
6.75e-04). The two advective face fluxes are each 1.3e-02 and cancel to
-3.23e-04; the diffusive inflow 9.79e-04 has no outward diffusive face
because cell 500 is the domain edge, and the advective divergence removes
only a third of it. The mass flux is outward at both faces, so the outer
boundary is a clean outflow and the advective donor is the interior.

**`molecular_scalar_gj1132_wellmixed/HeH0.55`** (40 outer passes, refused at
cell 306):

| | cell 305 | cell 306 (worst row) | cell 307 | cell 225 (worst drift) |
|---|---|---|---|---|
| r [R_p] | 1.9353 | 1.9517 | 1.9685 | 1.2309 |
| T [K] | 4493.7 | 4461.5 | 4429.0 | 2747.0 |
| x(H2) | 2.967e-06 | 2.664e-06 | 2.434e-06 | 8.880e-02 |
| residual | -7.282e-03 | -5.492e-03 | -4.146e-03 | -2.419e+02 |
| row scale | 2.945e-01 | 2.218e-01 | 1.678e-01 | 4.755e+04 |
| **measure** | 2.473e-02 | **2.477e-02** | 2.471e-02 | 5.088e-03 |
| diffusive, inner / outer face | -1.293e-01 / +9.841e-02 | -9.508e-02 / +7.219e-02 | -6.974e-02 / +5.272e-02 | -1.821e+04 / +1.774e+04 |
| advective, inner / outer face | -2.156e-02 / +1.919e-02 | -1.854e-02 / +1.680e-02 | -1.623e-02 / +1.495e-02 | -5.687e+03 / +5.460e+03 |
| mass flux, inner / outer face | +4.046e-08 / +3.978e-08 | +3.978e-08 / +3.911e-08 | +3.911e-08 / +3.845e-08 | +9.950e-08 / +9.885e-08 |
| production (H2+ + H) | 2.527e-02 | 2.386e-02 | 2.253e-02 | 3.120e+02 |
| loss | 5.127e-02 | 4.301e-02 | 3.667e-02 | 7.694e+02 |
| largest loss channel | He+ + H2 4.95e-02 (97 %) | 4.16e-02 (97 %) | 3.55e-02 (97 %) | H+ charge transfer 4.26e+02 |
| photodissociation + photoionization | 4.57e-05 | 3.95e-05 | 3.47e-05 | 2.03e+01 |

At cell 306 the balance is the opposite of cell 500's: the diffusive supply
(-2.29e-02) against the chemical destruction (-1.91e-02), with advection
1.7e-03, 7 %. Radiation is 0.09 % of the loss there; the destruction is the
He+ + H2 channel and the production the H2+ + H channel.

### 23.3 What holds both columns is the movement bound, not the rows

The question the review put — advective flux difference, diffusive face
switch, chemistry, or under-relaxation of the coupled update — is answered by
the pass trajectory and the relaxation's own ending, not by the terms alone.

Every carrier relaxation of both runs, in every pass, ends the same way:

```
carrier relaxation ended on the movement bound, with the last step
inside it; drift 1.67E-02 in 13 transport steps      (HeH0.083, early)
                               ... drift 1.68E-03 in 10 transport steps  (late)
```

Not one pass of either run ended on its own residual. The bound is halved
whenever the joint distance fails to fall — 5.0e-03, 2.5e-03, 1.25e-03,
then the floor 1.0e-03 — and the composition movement of a pass falls with
it in exact proportion (1.67e-02 -> 8.41e-03 -> 4.21e-03 -> 2.10e-03 ->
1.68e-03). The worst row measure does not fall with it:

| pass | HeH0.083 worst row | cell | HeH0.55 worst row | cell |
|---|---|---|---|---|
| 1 | 4.79e-02 | 342 | 1.27e-02 | 233 |
| 10 | 3.16e-02 | 348 | 3.93e-02 | 258 |
| 15 | 2.48e-02 | 500 | 3.08e-02 | 259 |
| 20 | 2.42e-02 | 500 | 2.70e-02 | 269 |
| 24 | 2.46e-02 | 500 | 2.29e-02 | 274 |
| 30 | — | — | 2.27e-02 | 293 |
| 40 | — | — | 2.48e-02 | 306 |

Two separate readings, and they are different states of the same block:

* **HeH0.083 is stalled.** The refusing cell reached the domain edge at pass
  14 and has stayed there; the measure has been flat between 2.41e-02 and
  2.46e-02 for ten passes, rising in the last four. The H2 front
  (x2 = 0.5) has not moved more than one cell since pass 19 (r = 1.2477 to
  1.2521) and x2 never falls below 1e-2 anywhere in the column. More passes
  at this bound do not descend.
* **HeH0.55 has not arrived.** The refusing cell migrates outward pass by
  pass (233 -> 306) at the same rate as the x2 = 1e-2 radius of the front
  (1.4816 -> 1.5075 R_p over the last five passes, about one cell a pass),
  and the measure has risen since pass 27. The state is still following a
  front that is moving, at a bound that admits about one cell of front
  travel per pass.

The cell that carries the drift and the cell that carries the refusing row
are **not the same cell and are not in the same regime**: the drift is at
the H2 front (cell 231 at 1.2565 R_p, x2 = 0.29 for HeH0.083; cell 225 at
1.2309 R_p, x2 = 0.089 for HeH0.55), whose row measure is 6.6e-03 and
5.1e-03 — a factor 4 to 5 INSIDE the refusing one — while the refusing rows
sit at 29 R_p and at 1.95 R_p in gas whose H2 is 2.7e-06 of the particles.
The bound is one number over the column, written on the particle-count
change, and it is the front that attains it; the cells whose rows refuse
move by orders of magnitude less than the allowance and are passengers on
whatever step the front leaves.

**Judgment.** The refusal of both columns is neither the chemistry nor a
face-flux switch nor a defect of the outer boundary. It is the carrier
relaxation being stopped by the composition movement bound in every pass
while the bound is cut geometrically to its floor — the same shape of defect
item L7e found in the OLD bound and fixed once, now visible in the new one:
the bound is a single scalar over a column that contains both a front that
must move slowly and a far wind that could converge at once. HeH0.083 is
stalled at that floor and will not descend with more passes of the same
kind; HeH0.55 is still travelling and its residual is not a statement about
its equations yet.

This is recorded as a measurement, not a fix. Nothing was changed in the
physics path for it, the row decomposition is a diagnostic behind
`EXHALE_CARRIER_DEBUG` and `EXHALE_CARRIER_ROW_TERMS` and is off by default,
and no cell was excluded from certification and no local equilibrium was
imposed anywhere.

### 23.4 What the three wellmixed cases are recorded as

| case | record |
|---|---|
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | **not solved**: one scalar movement bound over a column holding a slow H2 front and a far wind (the decomposition above) |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | the same |
| `molecular_scalar_gj1132_wellmixed/HeH2.13` | **not solved**: no certified molecular wellmixed state to seed it from, the two that would have been being the rows above |

A regional or coupled relaxation is the obvious thing to try and is NOT
tried here. It is opened as plan item L22 as a proposal only: what it would
have to answer before anything is written is over what region the fixed
wind's linear response is one statement, whether the right object is a bound
per region, a bound weighted by each cell's own contribution to that
response, or a relaxation that advances the front and the wind together so
that no bound of this kind is needed, and how any of the three is kept from
becoming a knob that certifies a state by loosening what holds it.
