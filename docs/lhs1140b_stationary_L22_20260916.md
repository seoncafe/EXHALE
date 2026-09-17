# LHS 1140 b stationary series, item L22: the composition movement bound, the progress measures and the relaxation contract

Item L22 of `docs/PLAN_20260916_rev3.md` section 3. This memo carries the
measurements of that item, one section per step. Every number is labeled
MEASURED (this memo ran it) or READ (taken from a source file, a run log
written earlier, or another document).

## Step 1: the progress measures, named and decided

Written 2026-09-16 (KST). The code is `src/EXHALE_main.f90` (the outer
iteration of `steady_wind_with_element_diffusion`) and
`src/modules/functions/binary_element_diffusion.f90` (the element
relaxation's returned measures and the new residual norm). Private build
`EXHALE_L22a.x`; control binary the tree's `EXHALE.x`
(md5 `c2e9c9990b9f14f1be8cd77abca68945`). Every run below is
`OMP_NUM_THREADS=1`.

### 1. What the four quantities are, and why they are four

| quantity | what it measures | what can move it |
|---|---|---|
| `carrier_residual_returned` | the residual of the carrier operator at the composition the pass hands back, `\|res\|` over the sum of that row's own physical terms (`carrier_steady_residual`) | the state alone |
| `element_residual_returned` | the residual of the elemental transport operator at the same composition, `\|res\|` over the sum of that row's own terms (`element_transport_residual_norm`, new) | the state alone |
| `element_map_distance`, `*_displacement_pass` | the distance to the endpoint of the element relaxation's inner map, and what each half of the pass actually applied | the movement bound `trust_pass` and the damping factor `comp_omega`, both of which the same loop shortens |
| `prog_worst` | the coupled residual: the largest certification entry over its own tolerance | the state alone |

The first two are what one pass of the alternation removes, the fourth is
what the acceptance is taken on, and the third is what a control shortens and
therefore decides nothing. Reading a displacement as progress reports "the
composition is closing" exactly when the loop has stopped letting it move,
which is the defect item R7 names and which the carrier half was already
free of (R23).

The element half was NOT free of it in the way the carrier half is. Its
measure was `element_composition_distance(Ynow, Ypass, Yres)`, formed before
`omega` damps the applied update, so it is undamped (R24) -- but it is the
distance to the endpoint of a FINITE inner relaxation, which exits on
`he_relax_maxstep` as readily as on its own movement measure, so it is not
the residual of the state handed back. It is now the elemental transport
residual, evaluated on the composition the pass returns after the damping,
the projection and the state refreshes, in the same shape the carrier half
is evaluated in. Both halves are dimensionless ratios of a row's residual to
the sum of that row's own physical terms; neither scale carries an iteration
length. The absolute residual and the absolute scale are printed beside each
ratio, so a movement of the scale cannot be read as a movement of the
balance.

`element_transport_residual_norm` reduces the rows `element_transport_residual`
returns by the reduction `certification_row_measure` forms for the same rows,
so the number the progress control reads and the elemental entries of the
certification of the same state are one measure. It is measured over the
whole column, as the carrier residual is; the certification's own GATE is at
`r >= 1.2 R_p` and is a separate statement about acceptance.

### 2. Unavailable is an outcome

`crc_max .eq. crc_max` rejects a NaN and nothing else: an infinity compares
equal to itself and passed it (R44), and nothing asked for a sign. The test
is now `residual_norm_is_admissible(x)`, which is
`x .eq. x .and. abs(x) .le. huge(x) .and. x .ge. 0`, applied to both halves.
A half whose measure is not admissible, or whose frozen background is not
ready (`bg_ready` false), or which had nothing to advance, is reported
UNAVAILABLE in the log; the joint progress rule is not applied to it and no
other quantity is read in its place. A pass in which neither half could be
measured records nothing in `comp_drift_last` / `comp_drift_before`, for the
same reason a pass that took no update records nothing: it is not evidence
either way about the contraction.

The control keys restore the superseded measures one half at a time, for a
comparison without a rebuild: `EXHALE_CARRIER_DRIFT_IS_DISPLACEMENT=1` (the
carrier displacement the pass kept) and `EXHALE_ELEMENT_DRIFT_IS_MAP_DISTANCE=1`
(the endpoint distance of the element relaxation's inner map). Both are in
the environment-key table of `README_HOWTO.md`.

### 3. The carrier relaxation's ending, in the outcome table

The outer loop reads `carrier_outcome` into the table of the contract
(`docs/PLAN_20260916_rev3.md` section 3 step 3) and prints the classification
beside the ending:

| inner ending | classification | what the outer iteration does |
|---|---|---|
| `carrier_relax_fixed_point` | the carrier block closed at this background | continue to the coupled check |
| the excursion bound, the step budget, a refused interval, a refused closure or the mask safety stop, with transport steps kept | an admissible PARTIAL update, not stationary | continue to the hydrodynamic solve; the bounding cell is named in the log as before |
| any of those with no transport step kept | NO ADMISSIBLE ADVANCE | the exact rejection condition is already printed; the movement bound is left where it is (shortening it is the wrong response to a pass that kept no step) and the pass counts as one without a fall |
| a returned composition that is not finite | unusable | the outer iteration ends with `outer_state_not_finite` |
| `carrier_relax_nothing_to_advance` | nothing to advance | the carrier half is UNAVAILABLE |

A nonzero physical row still cannot pass the FINAL certification; that is a
statement about the state and not about a pass.

### 4. The tests

`src/tests/element_operator` (RED before the change: the driver does not
compile against the tree objects, `Symbol 'element_transport_residual_norm'
... not found in module 'binary_element_diffusion'` and `Keyword argument
'displacement' ... is not in the procedure`; MEASURED). GREEN after, on the
synthetic element column of that suite:

| row | measured |
|---|---|
| `an_undamped_pass_displacement_equals_its_map_distance` | 0.00000E+00 of 1e-14 |
| `the_element_residual_of_this_state_is_measured` | 1 |
| `element_residual_norm_is_the_certification_row_measure` | 0.00000E+00 of 1e-14 (norm 1.60445E-13, certification row measure 1.60445E-13, at cell 40) |
| `the_two_row_scale_floors_are_one_number` | 0.00000E+00 of 0 |
| `the_relaxed_state_still_carries_an_element_row` | 1.60445E-13 > 0 |
| `a_damped_pass_displacement_is_omega_times_its_map_distance` | 4.44089E-16 of 1e-10 (map distance 5.00000E-01, displacement 2.50000E-01 at omega = 1/2) |

`src/tests/carrier_retry` carries the guard rows, with the superseded guard
transcribed beside the new one on the same values, which is the RED and the
GREEN in one run (MEASURED): `superseded_guard_accepts_an_infinite_residual`
= 1 and `superseded_guard_accepts_a_negative_residual` = 1, against
`an_infinite_residual_norm_is_not_admissible` = 0,
`a_negative_infinite_residual_norm_is_not_admissible` = 0,
`a_nan_residual_norm_is_not_admissible` = 0,
`a_negative_residual_norm_is_not_admissible` = 0, and the three that must
still be admitted, `a_zero_residual_norm_is_admissible` = 1,
`a_finite_positive_residual_norm_is_admissible` = 1,
`the_largest_finite_residual_norm_is_admissible` = 1.

Every suite whose driver links `binary_element_diffusion` passes with the
private build (MEASURED, exit 0 each): `element_operator`, `carrier_retry`,
`certification` (84 rows), `steady_species_rows` (195), `species_face_flux`,
`ionization_stage_flux` (15), `physics_probe` (1515).

### 5. The impact, measured

A shared acceptance interface changes -- the progress decision of the outer
iteration -- so the gate is a scoped measurement and not byte identity.

**The certified molecular case.**
`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13` reloaded on scratch
copies of its own `input.inp` and `base.inp` (only the spectrum path made
absolute) and its own written state as the IC pair, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=3`, single-threaded, the tree binary against
`EXHALE_L22a.x`. Both runs CERTIFY at outer pass 3 and the written state is
byte-identical except the `# provenance:` run timestamp (MEASURED:
`Ion_species.txt` and `Lyman_Werner.txt` identical, `Hydro_ioniz.txt`
differing in the one timestamp line). Every printed number of every pass
agrees -- `hydro info`, the three hydrodynamic rows, the worst gated species
row and its cell, the relaxation step counts and the endings.

What the new lines say, and what the superseded measure said on the same
passes (MEASURED):

| pass | element residual (new) | element map distance (old) | carrier residual (new default, and the tree's) | carrier displacement kept | coupled residual | decided by |
|---|---|---|---|---|---|---|
| 1 | 1.16E-05 = 5.81E-22 / 5.02E-17 at cell 2 | 5.56E-07 | 2.67E-09 = 1.69E-06 / 6.34E+02 | 1.31E-07 | 1.03E+00 | the coupled residual fell |
| 2 | 5.84E-06 = 2.93E-22 / 5.02E-17 at cell 2 | 2.81E-07 | 1.58E-09 = 1.00E-06 / 6.34E+02 | 5.70E-08 | 1.19E+00 | the composition residual, not yet measured twice, and the state's distance is not carried by a hydrodynamic row |
| 3 | not measured: the pass took no composition update | | | | 6.14E-01 | not asked: the pass did not reach the progress control |

The two element measures differ by a factor 21 in size and fall by the same
factor per pass on this case (1.986 for the residual, 1.979 for the map
distance, MEASURED), so the progress decision is the same at every pass and
the accepted states agree. The case does not separate the two measures; it
states that replacing one by the other did not move a certified solution.

**The control keys reproduce the superseded measure exactly.** The same
reload with `EXHALE_CARRIER_DRIFT_IS_DISPLACEMENT=1` and
`EXHALE_ELEMENT_DRIFT_IS_MAP_DISTANCE=1` reads element 5.56E-07 and 2.81E-07
and carrier 1.31E-07 and 5.70E-08 at passes 1 and 2, which are the tree
binary's own printed numbers to every digit (MEASURED), and its written
state is the tree's but for the timestamp.

**The marching snapshot.** `backup/regression/mol_carrier` (12000 steps,
`EXHALE_MAXSTEPS` from the case's `maxsteps`), run on scratch copies with
both binaries (MEASURED): `Ion_species.txt` and `Ion_species_adv.txt` are
byte-identical, and `Hydro_ioniz.txt` and `Hydro_ioniz_adv.txt` differ in
one line each, the `# provenance:` run timestamp. Every data line of every
file agrees. The case never enters the stationary outer iteration (no
`outer pass` line in either run log, MEASURED), so no line of this item is
on its path; the run is the standing check that nothing else moved.

**The pre-existing failing row is not moved.**
`src/tests/grid_and_gates/output_state_consistency.sh` fails
`outer_iteration_ending_is_the_stagnation_one` with `measured=pass_budget`
against `reference=no_progress` on BOTH binaries (MEASURED, the tree's
`EXHALE.x` and `EXHALE_L22a.x`); its other three rows pass identically
(`heat` 1.678011e-16, `cool` 2.135820e-14). The row is not this item's and
was not touched.

## Step 2: the causal comparisons

Written 2026-09-16 (KST). MEASURED on the three frozen states named in the
plan, reloaded and re-solved single-threaded with a private build. Nothing
in the physics path was changed for it: the comparison is made by an
environment-gated diagnostic that changes which cells choose the length of a
carrier trial, is off unless `EXHALE_L22_MASK` is set, and excludes no cell
from the certification.

**Verdict.** The refusing cell does not descend, in any of the three
interval-selection sets and in either mode, on any of the three states. In
`veto` mode the comparison is exactly the control, trial for trial and byte
for byte in the row dump, which is what that mode is for. In `waive` mode
the trial interval does grow, the control of the bound moves to the edge of
the mask rather than disappearing, and the refusing row is WORSE after the
next coupled pass than in the control. One run changed qualitatively
(`wellmixed/HeH0.55`, far-wind set, `waive`: the carrier relaxation reached
its own fixed point, carrier residual 1.05e-09 against 2.62e-01), and the
coupled solve that had to follow it did not converge.

### 1. The question, and what a result of it can say

The plan's question is whether the cell whose carrier row refuses the
certification descends when the cells that attain the composition movement
bound no longer choose the trial multiplier. The three states are those of
`docs/lhs1140b_stationary_L7e_20260915.md` section 23 (READ): in both
well-mixed states every carrier relaxation of every pass ended on the
movement bound and none on its own residual, and the cell that attains the
bound is not the cell that refuses.

### 2. The diagnostic, and what still holds a trial

`EXHALE_L22_MASK=<file>` names a list of physical cells, one integer per
line, that are OMITTED FROM THE INTERVAL SELECTION.
`EXHALE_L22_MASK_MODE` says what the omission means.

| mode | what the omitted cells do |
|---|---|
| `veto` | the particle-count bound is still evaluated on them and still refuses the trial, so the feasible set, and with it the accepted interval, is the one of the unmasked run; only the naming of the controlling cell changes. This is the null comparison, and it is reported as such. |
| `waive` | the bound is not evaluated on them, so the interval is chosen by the selection set alone. The run says in its log that it is a diagnostic experiment. |

What holds a trial in `waive` mode, all of it already in the code path and
none of it waived: the transport step's own interval coverage
(`carrier_interval_covered`), the positivity and the element and charge
feasibility of the chemistry closure
(`equilibrate_chemistry_at_fixed_conserved_state`), the finiteness and gas
conditions of `thermal_state_admissible` (finite p, T, heat, cool and eta
with p > 0 and T > 0), the element census of the pass
(`element_census_verify`), and a safety stop added for this experiment: a
trial in which any omitted cell's particle count `n_tot + n_e` changes by
more than a factor two ends the relaxation at the last accepted state with
the outcome `carrier_relax_mask_safety_stop`. An empty selection set is not
a comparison and ends the relaxation the same way with its own message.

The mask file is read once and frozen for the whole run, which is the
condition the comparison is read under. The code is in
`src/modules/lower_atmosphere/diffusive_photochemistry.f90`
(`l22_interval_mask_setup`, `l22_particle_count_change_by_set`, and the
bound block of `relax_photochemical_composition`).

### 3. The masks, and how they were measured

The mask is not an abundance threshold. It is read from the rejected trials
of the FIRST pass of the control run of each case, which the trial-by-trial
log now names (`EXHALE_CARRIER_DEBUG=1`, default off).

The width used to extend the bound-controlling cells to their neighbors is
the FRONT'S WIDTH measured as the cells over which x2 = 2 n(H2)/n_H falls
from 0.9 of its base-cell value to 0.1 of it, on the frozen state (this is
the first of the two definitions the plan offers; the second, the cells that
attain the bound in any trial of the pass, is reported beside it because it
is what the mask is centered on). The narrow mask is then the interval
spanned by the bound-controlling cells, extended by floor(width/2) cells on
each side and clipped to the physical grid.

| case | bound-controlling cells of the rejected trials of pass 1 (MEASURED) | cells attaining the bound in any trial of pass 1 (MEASURED) | front width, 0.9 to 0.1 of the base x2 (MEASURED) | narrow mask | mask radii |
|---|---|---|---|---|---|
| `wellmixed/HeH0.083` | 224 | 223, 224 | 131 cells (167 to 297) | cells 159 to 289, 131 cells | 1.0722 to 1.7073 R_p |
| `wellmixed/HeH0.55` | 219, 220, 221 | 219, 220, 221, 223 to 227 | 97 cells (149 to 245) | cells 171 to 269, 99 cells | 1.0893 to 1.4987 R_p |
| `kzz1e9/HeH0.083` | 51, 53, 54, 55 | 49 to 55, and 461 to 500 | 300 cells (198 to 497) | cells 1 to 205, 205 cells | 1.0002 to 1.1626 R_p |

The third interval-selection set is the FAR WIND: the cells at r >= 5 R_p,
which on this grid is cells 389 to 500 (r = 5.0484 R_p at cell 389,
MEASURED); the mask file omits cells 1 to 388. The three cases share one
grid, so one file serves all three. This is the "cells above 5 R_p" option
of the plan and not the Kzz case's own bound cell, and it is stated here
because the two are different sets.

Worth recording before the tables: in `kzz1e9/HeH0.083` the bound of the
first pass of this reload is attained at cells 49 to 55 (r = 1.014 to
1.018 R_p), the near-base molecular layer, and not at the far wind cell 466
(16.5 R_p) that the campaign run of 2026-09-15 recorded at its pass 40
(READ, `docs/lhs1140b_stationary_L7e_20260915.md` section 23). Which cell
carries the bound is a property of the pass and not of the case.

### 4. What was run

Each case was reloaded from its own last written state with `Load IC? True`
and the case's own `input.inp` and `base.inp` copied line for line (only the
spectrum path was made absolute for the scratch directory), with
`EXHALE_PTC_DTAU0=1.0e8` as its `REPRODUCE.md` states for the continuation
(READ), `EXHALE_OUTER_PASSES=3`, `OMP_NUM_THREADS=1`,
`EXHALE_CARRIER_DEBUG=1` and `EXHALE_CARRIER_ROW_TERMS=1`. The horizon of
the experiment is TWO outer passes: the first gives the carrier residual
immediately after the carrier update, the second gives the certification and
the carrier residual after a complete coupled outer update, which is what
the plan asks for. Each run was stopped once its second pass had written its
diagnostics, so the final state file was not written and the pressure column
of the written state is therefore not part of this measurement; T, the
composition, the transport terms, the sources and the row residuals over the
whole column come from the `carrier_row_terms.txt` dump of the certification,
before (the frozen state) and after (the certification at the head of the
second pass).

### 5. The comparisons

All numbers below are MEASURED in this item's runs. The measure of the
refusing row is the one the certification uses, the residual over the sum of
the row's own terms with the two faces counted separately, and the
certification's own verdict line is quoted where the pass reached one.

#### `molecular_scalar_gj1132_wellmixed/HeH0.083`, refusing cell 500 (r = 29.0312 R_p)

Frozen state, the row the certification refuses: residual 6.747E-04, faces scale measure 2.462E-02, x2 7.807E-03, T 285.3 K.

| set / mode | omitted | interval-selection cells | accepted trials, grow min..max | rejections (reason x count) | cell that set or vetoed the bound | kept steps | drift | relaxation ending |
|---|---|---|---|---|---|---|---|---|
| full grid (control), pass 1 | 0 | 500 | 14, 1.485E-03..2.563E+01 | the movement bound x20 | 224 | 14 | 0.01693 | the movement bound, with the last step inside it |
| full grid (control), pass 2 | 0 | 500 | 16, 3.341E-03..2.563E+01 | the movement bound x21 | 226,227 | 16 | 0.01697 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 1 | 131 | 369 | 14, 1.485E-03..2.563E+01 | the movement bound x20 | 224 | 14 | 0.01693 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 2 | 131 | 369 | 16, 3.341E-03..2.563E+01 | the movement bound x21 | 226,227 | 16 | 0.01697 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 1 | 131 | 369 | 16, 3.341E-03..3.844E+01 | the movement bound x21 | 290 | 16 | 0.03529 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 2 | 131 | 369 | 17, 1.253E-03..3.844E+01 | the movement bound x21 | 290 | 17 | 0.03137 | the movement bound, with the last step inside it |
| far wind, veto, pass 1 | 388 | 112 | 14, 1.485E-03..2.563E+01 | the movement bound x20 | 224 | 14 | 0.01693 | the movement bound, with the last step inside it |
| far wind, veto, pass 2 | 388 | 112 | 16, 3.341E-03..2.563E+01 | the movement bound x21 | 226,227 | 16 | 0.01697 | the movement bound, with the last step inside it |
| far wind, waive, pass 1 | 388 | 112 | 26, 1.505E-03..2.919E+02 | the movement bound x27 | 389 | 26 | 0.2132 | the movement bound, with the last step inside it |
| far wind, waive, pass 2 | 388 | 112 | 24, 1.338E-03..6.568E+02 | the movement bound x26 | 389 | 24 | 0.179 | the movement bound, with the last step inside it |

| set / mode | refusing cell row measure, pass 1 entry | after the carrier update (carrier steady residual, volume / worst / cell) | refusing cell row measure after the coupled pass | after that pass's carrier update | hydrodynamic rows of pass 2 (mass / momentum / energy) | H2 front x2=0.5 [R_p], pass 1 -> 2 | safety stops |
|---|---|---|---|---|---|---|---|
| full grid (control) | 2.462E-02 | 7.640E-03 / 1.910E-01 / 1 | 3.345E-02 at cell 500 | 7.440E-03 / 1.830E-01 / 1 | 8.340E-09 / 6.790E-13 / 2.690E-08 | 1.2565 -> 1.2657 | 0 |
| narrow mask, veto | 2.462E-02 | 7.640E-03 / 1.910E-01 / 1 | 3.345E-02 at cell 500 | 7.440E-03 / 1.830E-01 / 1 | 8.340E-09 / 6.790E-13 / 2.690E-08 | 1.2565 -> 1.2657 | 0 |
| narrow mask, waive | 2.462E-02 | 7.350E-03 / 1.750E-01 / 1 | 4.537E-02 at cell 500 | 6.940E-03 / 1.610E-01 / 1 | 2.030E-08 / 2.230E-13 / 1.180E-08 | 1.2704 -> 1.2850 | 0 |
| far wind, veto | 2.462E-02 | 7.640E-03 / 1.910E-01 / 1 | 3.345E-02 at cell 500 | 7.440E-03 / 1.830E-01 / 1 | 8.340E-09 / 6.790E-13 / 2.690E-08 | 1.2565 -> 1.2657 | 0 |
| far wind, waive | 2.462E-02 | 5.240E-03 / 6.070E-02 / 1 | 2.362E-01 at cell 500 | 4.840E-03 / 4.150E-02 / 1 | 2.150E-08 / 1.200E-13 / 1.590E-08 | 1.3905 -> 1.6151 | 0 |

Column-wide change of the row dump against the frozen state (largest relative change over the 500 physical cells; the dump is the certification of the second pass where that pass was reached):

| set / mode | T | x2 | diffusive | advective | production | loss | net source | residual |
|---|---|---|---|---|---|---|---|---|
| full grid (control) | 1.61E-02 | 1.69E-01 | 6.72E+01 | 4.29E+00 | 1.32E-01 | 1.02E-01 | 8.31E-01 | 4.54E+01 |
| narrow mask, veto | 1.61E-02 | 1.69E-01 | 6.72E+01 | 4.29E+00 | 1.32E-01 | 1.02E-01 | 8.31E-01 | 4.54E+01 |
| narrow mask, waive | 3.38E-02 | 3.70E-01 | 1.30E+02 | 9.58E+00 | 2.56E-01 | 1.99E-01 | 1.67E+00 | 9.37E+01 |
| far wind, veto | 1.61E-02 | 1.69E-01 | 6.72E+01 | 4.29E+00 | 1.32E-01 | 1.02E-01 | 8.31E-01 | 4.54E+01 |
| far wind, waive | 2.04E-01 | 4.06E+00 | 4.04E+02 | 7.69E+01 | 8.35E-01 | 7.31E-01 | 7.36E+00 | 4.39E+02 |

#### `molecular_scalar_gj1132_wellmixed/HeH0.55`, refusing cell 306 (r = 1.9517 R_p)

Frozen state, the row the certification refuses: residual -5.492E-03, faces scale measure -2.477E-02, x2 2.664E-06, T 4461.5 K.

| set / mode | omitted | interval-selection cells | accepted trials, grow min..max | rejections (reason x count) | cell that set or vetoed the bound | kept steps | drift | relaxation ending |
|---|---|---|---|---|---|---|---|---|
| full grid (control), pass 1 | 0 | 500 | 17, 9.766E-04..5.767E+01 | the movement bound x22 | 219,220,221 | 17 | 0.02904 | the movement bound, with the last step inside it |
| full grid (control), pass 2 | 0 | 500 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 226,227 | 18 | 0.02898 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 1 | 99 | 401 | 17, 9.766E-04..5.767E+01 | the movement bound x22 | 219,220,221 | 17 | 0.02904 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 2 | 99 | 401 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 226,227 | 18 | 0.02898 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 1 | 99 | 401 | 18, 9.766E-04..8.650E+01 | the movement bound x22 | 170 | 18 | 0.04328 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 2 | 99 | 401 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 278 | 18 | 0.03058 | the movement bound, with the last step inside it |
| far wind, veto, pass 1 | 388 | 112 | 17, 9.766E-04..5.767E+01 | the movement bound x22 | 219,220,221 | 17 | 0.02904 | the movement bound, with the last step inside it |
| far wind, veto, pass 2 | 388 | 112 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 226,227 | 18 | 0.02898 | the movement bound, with the last step inside it |
| far wind, waive, pass 1 | 388 | 112 | 70, 1.000E+00..1.414E+12 | none | (none refused) | 70 | 0.1599 | the carriers stopped moving at the fixed wind |

| set / mode | refusing cell row measure, pass 1 entry | after the carrier update (carrier steady residual, volume / worst / cell) | refusing cell row measure after the coupled pass | after that pass's carrier update | hydrodynamic rows of pass 2 (mass / momentum / energy) | H2 front x2=0.5 [R_p], pass 1 -> 2 | safety stops |
|---|---|---|---|---|---|---|---|
| full grid (control) | 2.477E-02 | 7.010E-03 / 2.620E-01 / 1 | 2.981E-02 at cell 297 | 6.460E-03 / 2.420E-01 / 1 | 1.150E-08 / 4.230E-13 / 2.280E-08 | 1.1626 -> 1.1684 | 0 |
| narrow mask, veto | 2.477E-02 | 7.010E-03 / 2.620E-01 / 1 | 2.981E-02 at cell 297 | 6.460E-03 / 2.420E-01 / 1 | 1.150E-08 / 4.230E-13 / 2.280E-08 | 1.1626 -> 1.1684 | 0 |
| narrow mask, waive | 2.477E-02 | 6.650E-03 / 2.430E-01 / 1 | 3.425E-02 at cell 300 | 6.140E-03 / 2.290E-01 / 1 | 1.850E-08 / 2.440E-13 / 1.200E-08 | 1.1654 -> 1.1714 | 0 |
| far wind, veto | 2.477E-02 | 7.010E-03 / 2.620E-01 / 1 | 2.981E-02 at cell 297 | 6.460E-03 / 2.420E-01 / 1 | 1.150E-08 / 4.230E-13 / 2.280E-08 | 1.1626 -> 1.1684 | 0 |
| far wind, waive | 2.477E-02 | 7.320E-10 / 1.050E-09 / 309 | not reached | not reached | not reached | 1.1542 | 0 |

Column-wide change of the row dump against the frozen state (largest relative change over the 500 physical cells; the dump is the certification of the second pass where that pass was reached):

| set / mode | T | x2 | diffusive | advective | production | loss | net source | residual |
|---|---|---|---|---|---|---|---|---|
| full grid (control) | 4.00E-02 | 3.41E+00 | 2.93E+01 | 3.89E+01 | 1.67E+00 | 2.34E+00 | 4.48E+00 | 4.11E+01 |
| narrow mask, veto | 4.00E-02 | 3.41E+00 | 2.93E+01 | 3.89E+01 | 1.67E+00 | 2.34E+00 | 4.48E+00 | 4.11E+01 |
| narrow mask, waive | 5.73E-02 | 6.43E+00 | 4.05E+01 | 7.19E+01 | 2.64E+00 | 3.91E+00 | 8.01E+00 | 1.15E+02 |
| far wind, veto | 4.00E-02 | 3.41E+00 | 2.93E+01 | 3.89E+01 | 1.67E+00 | 2.34E+00 | 4.48E+00 | 4.11E+01 |
| far wind, waive | 0.00E+00 | 0.00E+00 | 7.81E-08 | 2.58E-03 | 0.00E+00 | 0.00E+00 | 0.00E+00 | 9.56E-08 |

#### `molecular_scalar_gj1132_kzz1e9/HeH0.083`, refusing cell 500 (r = 29.0312 R_p)

Frozen state, the row the certification refuses: residual 4.638E-04, faces scale measure 2.972E-03, x2 3.271E-02, T 216.5 K.

| set / mode | omitted | interval-selection cells | accepted trials, grow min..max | rejections (reason x count) | cell that set or vetoed the bound | kept steps | drift | relaxation ending |
|---|---|---|---|---|---|---|---|---|
| full grid (control), pass 1 | 0 | 500 | 34, 9.643E-03..7.482E+03 | the movement bound x31 | 51,53,54,55 | 34 | 0.02234 | the movement bound, with the last step inside it |
| full grid (control), pass 2 | 0 | 500 | 26, 1.505E-03..4.379E+02 | the movement bound x27 | 455,456,457 | 26 | 0.02623 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 1 | 205 | 295 | 34, 9.643E-03..7.482E+03 | the movement bound x31 | 51,53,54,55 | 34 | 0.02234 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 2 | 205 | 295 | 26, 1.505E-03..4.379E+02 | the movement bound x27 | 455,456,457 | 26 | 0.02623 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 1 | 205 | 295 | 34, 1.205E-03..7.482E+03 | the movement bound x31 | 206 | 34 | 0.02277 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 2 | 205 | 295 | 22, 9.514E-03..4.379E+02 | the movement bound x24 | 455,456,457 | 22 | 0.0262 | the movement bound, with the last step inside it |
| far wind, veto, pass 1 | 388 | 112 | 34, 9.643E-03..7.482E+03 | the movement bound x31 | 51,53,54,55 | 34 | 0.02234 | the movement bound, with the last step inside it |
| far wind, veto, pass 2 | 388 | 112 | 26, 1.505E-03..4.379E+02 | the movement bound x27 | 455,456,457 | 26 | 0.02623 | the movement bound, with the last step inside it |
| far wind, waive, pass 1 | 388 | 112 | 34, 2.411E-03..1.122E+04 | the movement bound x31 | 407,415,416,417,418 | 34 | 0.03551 | the movement bound, with the last step inside it |
| far wind, waive, pass 2 | 388 | 112 | 21, 9.766E-04..1.946E+02 | the movement bound x24 | 457,458 | 21 | 0.02575 | the movement bound, with the last step inside it |

| set / mode | refusing cell row measure, pass 1 entry | after the carrier update (carrier steady residual, volume / worst / cell) | refusing cell row measure after the coupled pass | after that pass's carrier update | hydrodynamic rows of pass 2 (mass / momentum / energy) | H2 front x2=0.5 [R_p], pass 1 -> 2 | safety stops |
|---|---|---|---|---|---|---|---|
| full grid (control) | 2.972E-03 | 2.130E-03 / 1.210E-02 / 1 | 1.391E-02 at cell 500 | 2.390E-03 / 1.310E-02 / 2 | 1.660E-01 / 1.210E-03 / 4.270E-02 | 1.3705 -> 1.3641 | 0 |
| narrow mask, veto | 2.972E-03 | 2.130E-03 / 1.210E-02 / 1 | 1.391E-02 at cell 500 | 2.390E-03 / 1.310E-02 / 2 | 1.660E-01 / 1.210E-03 / 4.270E-02 | 1.3705 -> 1.3641 | 0 |
| narrow mask, waive | 2.972E-03 | 2.110E-03 / 1.200E-02 / 1 | 1.542E-02 at cell 500 | 2.110E-03 / 1.450E-02 / 3 | 3.780E-01 / 2.870E-04 / 9.940E-01 | 1.3770 -> 1.3641 | 0 |
| far wind, veto | 2.972E-03 | 2.130E-03 / 1.210E-02 / 1 | 1.391E-02 at cell 500 | 2.390E-03 / 1.310E-02 / 2 | 1.660E-01 / 1.210E-03 / 4.270E-02 | 1.3705 -> 1.3641 | 0 |
| far wind, waive | 2.972E-03 | 1.860E-03 / 1.070E-02 / 1 | 2.174E-02 at cell 500 | 2.300E-03 / 1.370E-02 / 1 | 1.270E-08 / 3.980E-13 / 7.890E-09 | 1.3905 -> 1.3973 | 0 |

Column-wide change of the row dump against the frozen state (largest relative change over the 500 physical cells; the dump is the certification of the second pass where that pass was reached):

| set / mode | T | x2 | diffusive | advective | production | loss | net source | residual |
|---|---|---|---|---|---|---|---|---|
| full grid (control) | 2.81E-01 | 3.56E-02 | 1.23E+01 | 2.44E+04 | 6.27E-01 | 1.42E-01 | 1.04E+02 | 4.58E+02 |
| narrow mask, veto | 2.81E-01 | 3.56E-02 | 1.23E+01 | 2.44E+04 | 6.27E-01 | 1.42E-01 | 1.04E+02 | 4.58E+02 |
| narrow mask, waive | 2.13E-01 | 3.61E-02 | 1.41E+01 | 1.65E+04 | 5.48E-01 | 1.33E-01 | 1.12E+02 | 4.89E+02 |
| far wind, veto | 2.81E-01 | 3.56E-02 | 1.23E+01 | 2.44E+04 | 6.27E-01 | 1.42E-01 | 1.04E+02 | 4.58E+02 |
| far wind, waive | 3.04E-02 | 5.68E-02 | 2.29E+01 | 2.20E+00 | 1.99E-01 | 1.03E-01 | 1.67E+02 | 7.36E+02 |

#### 5.1 What the three cases say

**The `veto` mode is the null comparison, and it is exactly null.** In all
three cases the `veto` runs reproduce the control trial for trial: the same
number of accepted trials, the same `grow` at every one of them, the same
bound measure, the same kept steps, the same drift, the same certification
verdict, and a `carrier_row_terms.txt` dump byte-identical to the control's
(MEASURED, `cmp` of the control's dump against both `veto` dumps in all
three cases: six comparisons, six byte-identical). That is the designed behavior: with the bound still
evaluated on the omitted cells the feasible set is unchanged, so the
experiment is unchanged, and the mask changes only which cell the log names.

**In `waive` mode the control of the bound moves to the edge of the mask and
does not disappear.** MEASURED:

| case | set | cell that set the bound in the control, pass 1 | cell that set it with the bound waived | mask boundary |
|---|---|---|---|---|
| `wellmixed/HeH0.083` | narrow | 224 | 290 | mask ends at 289 |
| `wellmixed/HeH0.083` | far wind | 224 | 389 | selection set begins at 389 |
| `wellmixed/HeH0.55` | narrow | 219, 220, 221 | 170 | mask begins at 171 |
| `wellmixed/HeH0.55` | far wind | 219, 220, 221 | no trial refused | selection set begins at 389 |
| `kzz1e9/HeH0.083` | narrow | 51, 53, 54, 55 | 206 | mask ends at 205 |
| `kzz1e9/HeH0.083` | far wind | 51, 53, 54, 55 | 407, 415 to 418 | selection set begins at 389 |

In five of the six the cell that takes over is the first cell of the
selection set or its immediate neighbor, so the bound is not carried by an
isolated cell but by a band whose edge the mask merely relocates. The
accepted interval grows by a factor 1.5 (`HeH0.083` narrow, `grow` 25.6 to
38.4) to 11 (`HeH0.083` far wind, 25.6 to 292; `HeH0.55` far wind, 57.7 to
1.4e+12), the pass keeps one to three more steps, and the composition drift
of the pass rises in the same proportion.

**The refusing cell does not descend, in any set and in either mode.** Over
the two-pass horizon its row measure RISES in every run, and it rises MORE
where the bound was waived (MEASURED, entry of pass 1 to the certification
at the head of pass 2):

| case | full grid (control) | narrow mask, waive | far wind, waive |
|---|---|---|---|
| `wellmixed/HeH0.083`, cell 500 | 2.462e-02 -> 3.345e-02 | 2.462e-02 -> 4.537e-02 | 2.462e-02 -> 2.362e-01 |
| `wellmixed/HeH0.55`, cell 306 | 2.477e-02 -> 2.981e-02 (worst row moved to cell 297) | 2.477e-02 -> 3.425e-02 (cell 300) | pass 2 not reached, see below |
| `kzz1e9/HeH0.083`, cell 500 | 2.972e-03 -> 1.391e-02 | 2.972e-03 -> 1.542e-02 | 2.972e-03 -> 2.174e-02 |

The carrier residual measured immediately after the carrier update does fall
slightly with the longer interval (`HeH0.083`, volume-weighted 7.64e-03
control against 7.35e-03 narrow-waive and 5.24e-03 far-waive; worst cell
1.91e-01 against 1.75e-01 and 6.07e-02), and it is still the base cell 1
that carries it. After the coupled update of the next pass the same ordering
holds and the certified row is worse, so the small gain of the carrier half
does not survive the wind's response to it.

**The one qualitative change, and why it is not an answer.** On
`wellmixed/HeH0.55` with the far-wind selection set in `waive` mode the
carrier relaxation ended, for the first time on this route, on its OWN
criterion rather than on the bound: 70 kept steps, `grow` from 1.0 to
1.4e+12, not one trial refused, ending "the carriers stopped moving at the
fixed wind", and the carrier steady residual after that update
7.32e-10 volume-weighted with a worst cell of 1.05e-09 at cell 309, against
7.01e-03 and 2.62e-01 in the control. The state it handed back is one the
coupled solve could not follow: the pass-2 hydrodynamic solve was still at
`||R||` = 1.118 after 205 JFNK iterations with the line-search damping at
9.5e-07, the pseudo-time step at 7.3e-05 and all 500 mass rows outside their
tolerance, having plateaued near `||R||` = 1.1 to 1.4 from iteration 100 on,
so the run was stopped there and the second pass has no certification. The
carrier block closing at a FIXED wind is not the same statement as the
coupled system closing, and this run is the clearest instance of the
difference.

**The `kzz1e9` case, where the control's own hydrodynamic solve is the
weaker half.** Its pass-2 solve returns `info = 2` in the control and in
both `veto` runs and in the narrow `waive` run, with pass-2 hydrodynamic
rows of 1.66e-01 (mass), 1.21e-03 (momentum), 4.27e-02 (energy); the
far-wind `waive` run is the only one of the five whose pass-2 solve returns
`info = 0`, with rows 1.27e-08, 3.98e-13 and 7.89e-09. The refusing row
nonetheless rose further there than in the control (2.174e-02 against
1.391e-02), so the better hydrodynamic solve did not carry the carrier row
with it.

**The safety stop and the empty-set outcome were not exercised.** No trial
in any of the twelve masked runs moved an omitted cell's particle count by
more than a factor two, so the safety stop fired zero times, and no
selection set was empty (the smallest was 112 cells).

### 6. Byte identity with the gate off

`backup/regression/mol_carrier` was run on a scratch copy, single-threaded
(`OMP_NUM_THREADS=1`, `EXHALE_MAXSTEPS=12000` from the case's `maxsteps`),
once with the tree binary `EXHALE.x` (md5 `c2e9c9990b9f14f1be8cd77abca68945`)
and once with `EXHALE_L22.x` (md5 `f23da8d84fafc2f0587b610add56e714`), with
`EXHALE_L22_MASK` unset. MEASURED:

| file | result |
|---|---|
| `output/Ion_species.txt` | byte-identical (md5 `eb197f9cad02ac854eb68900a4bee425` both) |
| `output/Hydro_ioniz.txt` | identical in every data line; the two files differ only in the `# provenance:` header, which records the run's own wall-clock start time (`19:00:54` against `19:00:51`) |

No golden was refreshed and no golden was read for this.

### 7. What this measurement does not settle

The plan's wording is adopted as it stands: no improvement was observed for
the tested masks and intervals; the named competing explanations (outer
boundary carrier treatment, the deferred reconstruction terms of the
Jacobian) are tested before a production design is chosen.

Where a number DID move (the far-wind selection set in `waive` mode on
`wellmixed/HeH0.55`, where the carrier relaxation reached its own fixed point
for the first time on this route), the reading is the one the plan requires:
a positive result may come from a region leaving the small-response regime
rather than from the strategy, so the bound and its control stay the leading
hypothesis, not a proven cause, until an intervention separates them. In
that particular run the front's cells were free of the bound and the
composition moved by a factor to which the fixed wind cannot be expected to
respond linearly, which is exactly the regime question and not an answer to
it.

No production design is chosen here, and step 3 of the plan is not entered.

### 8. Noticed beside the measurement, not acted on

The `carrier_row_terms.txt` dump carries two columns named `scale_terms` and
`measure_terms`, and they are NOT the scale and the measure the
certification gates on. The dump's scale is the sum of the magnitudes of the
NET divergences, |diffusive| + |advective| + |net source|, while the
certification's scale counts the two faces separately,
|F_dif,in| + |F_dif,out| + |F_adv,in| + |F_adv,out| + |net source|. MEASURED
on `wellmixed/HeH0.083` cell 500 of the frozen state: the dump's columns read
1.3216e-03 and 5.105e-01, while the faces give 2.7419e-02 and 2.4619e-02,
the second of which is the number the certification prints. A reader who
compares the dump's `measure_terms` against the certification verdict gets a
different number for the same row. The dump's face lines are written beside
every row, so the certification's scale is recoverable from the file; this
is a naming question and not a defect of either quantity, and it is reported
and not changed.

### 9. Where the record is

`docs/audit_20260905/L22_step2_20260916/`: `IDENTITIES.md` (the md5 of every
frozen state file reloaded, of the measuring binary and of the edited
source), the four mask files with their definitions in their headers,
`TABLES.md` (the tables above as the analysis produced them), `logs/` (the
run log of all fifteen runs and of the three evaluations of the frozen
states) and `row_terms/` (the `carrier_row_terms.txt` of each).

## Step 2b: the two named competing explanations, tested

Written 2026-09-16 (KST). Item L22 of `docs/PLAN_20260916_rev3.md` section 3
step 2 requires the two named competing explanations to be tested before any
production design is chosen. This section carries that test. The record is
`docs/audit_20260905/L22_step2b_20260916/`. Private build `EXHALE_L22b.x`
(md5 `1152b33c11612518dc78a54b11afe5a1`); the tree binary `EXHALE.x`
(md5 `c2e9c9990b9f14f1be8cd77abca68945`) is untouched. Every run is
`OMP_NUM_THREADS=1`. Every number is labeled MEASURED or READ.

**Verdict, one line each.**

| explanation | verdict | the numbers it rests on |
|---|---|---|
| (A) the outer-boundary carrier treatment | SUPPORTED as a description of the cell-500 row and EXCLUDED as its cause or its remedy | the outer face carries no diffusive, eddy or settling flux at all, so the row at cell 500 is a one-sided influx balanced against advection alone; but all three interventions make that row WORSE on the frozen state, by 1.12x to 20.6x on `wellmixed/HeH0.083` and 2.7x to 182x on `kzz1e9/HeH0.083`, and none of them moves the `HeH0.55` row at cell 306 by a single digit |
| (B) the deferred reconstruction terms of the Jacobian | MEASURED and EXCLUDED as the cause | the assembled action reproduces a central difference of the full residual to a mean relative error of 1.5e-03 over the direction's support, the entries the tridiagonal band cannot hold reach 9.4e-02 of a row's band, and the relaxation's Newton already ends every substep at a row residual of 1e-13 to 1e-17 with every interval covered in one attempt, so the Jacobian is not what limits a pass |
| (C) the coupled failure after the far-wind waive | MEASURED, and it is a displacement of the state the wind reads | the composition that relaxation hands back moves the pressure by 5.0e-02 and the temperature by 1.9e-02 at the BASE against 8.9e-03 and 5.0e-03 in the control, and the hydrodynamic family that plateaued is the MASS row, all 500 cells outside their tolerance with the worst at cell 123, r = 1.038, in that same inner layer |

No production design is chosen here.

### 1. What the outer boundary of the carrier column actually is

READ from `src/modules/lower_atmosphere/diffusive_photochemistry.f90`, and
this is the first thing the item had to establish, because the plan's wording
(a "copied ghost" at the outer diffusive face) does not describe the code:

1. `carrier_face_coefficients` fills the face arrays over `j = 1, N-1` only.
   `Agrd` and `Bdrf` are dimensioned `0:N` and initialized to zero, so the
   face at `r_{N+1/2}` and the face at `r_{1/2}` carry NO molecular diffusion,
   NO eddy diffusion and NO settling drift. The carrier column is a CLOSED
   BOX for that flux at both ends.
2. `carrier_residual` forms the face fluxes over the same `j = 1, N-1`, so
   `Jf(0)` and `Jf(N)` stay at their initialized zero. There is no ghost in
   the diffusive term to copy or to extrapolate: the term is not formed.
3. The outer ghost enters the row through the ADVECTIVE face alone, as the
   reconstruction `species_face_fraction` takes of the carrier mass fraction
   `Y = m_c f_c / msum` over the cells `N-1, N, N+1`.
4. And the two evaluations of that row do not use the same ghost.
   `carrier_state` fills `fc` over the whole array out of `f_sp`, whose outer
   ghost cells the ionization sweep solves (`ionization_equilibrium.f90`
   loops `1-Ng, N+Ng`), so in the CERTIFICATION the outer ghost is DATA.
   Inside `solve_carriers` every line-search trial refills
   `ftry(N+1:N+Ng,ic) = ftry(N,ic)`, a COPY. MEASURED consequence, in the new
   suite: with a copied ghost the forward difference entering the MC limiter
   of PLM is exactly zero, the limited slope of cell N is exactly zero and
   the reconstruction returns the donor cell average
   (`the_copied_outer_ghost_leaves_no_slope_on_the_outflow_face` = 0.0), while
   a ghost that is not a copy leaves a slope of 3.62e-02 of the cell value
   (`a_ghost_that_is_not_a_copy_leaves_a_slope_on_the_outflow_face`). These
   cases run `Reconstruction scheme: PLM` (READ from their `input.inp`).

So the certification measures an outer advective face that the relaxation's
own trials do not use. That is a defect of its own and it is reported here;
it is NOT fixed by this item.

### 2. The refusing rows, term by term, on the frozen states

MEASURED, from the `carrier_row_terms.txt` of the certification of each frozen
state (`EXHALE_CARRIER_ROW_TERMS=1`). The two faces of each flux are separate
columns of that dump, so the row is recovered exactly:
`residual = (F_dif,in + F_dif,out) + (F_adv,in + F_adv,out) - net source`.
All rates are volumetric, cm^-3 s^-1.

| case, cell | r [R_p] | T [K] | x_c | F_dif,in | F_dif,out | F_adv,in | F_adv,out | net source | residual | measure |
|---|---|---|---|---|---|---|---|---|---|---|
| `wellmixed/HeH0.083`, 500 | 29.031 | 285.3 | 7.807e-03 | +9.7857e-04 | **+0.0** | -1.3366e-02 | +1.3042e-02 | -1.9574e-05 | +6.7473e-04 | 2.4619e-02 |
| `kzz1e9/HeH0.083`, 500 | 29.031 | 216.5 | 3.271e-02 | +3.4783e-03 | **+0.0** | -7.7813e-02 | +7.4743e-02 | -5.5024e-05 | +4.6382e-04 | 2.9715e-03 |
| `wellmixed/HeH0.55`, 306 | 1.952 | 4461.5 | 2.664e-06 | -9.5084e-02 | +7.2185e-02 | -1.8544e-02 | +1.6801e-02 | -1.9149e-02 | -5.4924e-03 | 2.4767e-02 |

The two cell-500 rows have an outer diffusive flux of EXACTLY zero, and in
both the diffusive influx is the largest single term of the row: 1.45 times
the residual on `HeH0.083` and 7.5 times it on the Kzz case, with the
advective divergence taking up the rest and the chemistry contributing 2.9 and
11.9 percent of the residual. Cell 306 of `HeH0.55` is an interior cell with
both diffusive faces live and nearly cancelling, and its residual is the
residue of that cancellation against the He+ + H2 destruction; the boundary
has nothing to do with it.

How much of the refusal is the LAST CELL, cell by cell (MEASURED, the
certification measure over the gated column `r >= 1.2 R_p`):

| case | worst gated row | the same excluding cell N | ratio | the cells beside it |
|---|---|---|---|---|
| `wellmixed/HeH0.083` | 2.4619e-02 at 500 | 2.1713e-02 at 499 | 1.13 | a band at cells 339-344 reads 2.01e-02 to 2.06e-02, so the whole outer column sits at that level |
| `kzz1e9/HeH0.083` | 2.9715e-03 at 500 | 3.5188e-04 at 499 | 8.44 | 493-499 fall smoothly 2.27e-04 to 3.52e-04, so cell 500 IS a single-cell spike |
| `wellmixed/HeH0.55` | 2.4767e-02 at 306 | the same row | 1.00 | 302-309 all read 2.42e-02 to 2.48e-02, a band and not a cell |

### 3. The interventions (A1), (A2), (A3), measured on the frozen state

Each is one environment key, off by default, and the operator with all of them
off is the one the tree integrates (section 7 below). The three frozen states
were reloaded exactly as step 2 reloaded them, and every intervention was
measured on a state whose carrier densities agree with the control's to every
digit printed (largest relative difference of n(H2) over the 500 cells 0.0,
MEASURED), so what follows compares operators and not states. The control
column reproduces step 2's numbers for the same rows to every digit printed.

| case, cell | control | (A1) copy, X_ghost = X_N | (A1) extrap, X_ghost = 2 X_N - X_{N-1} | (A3) donor-cell outflow face |
|---|---|---|---|---|
| `wellmixed/HeH0.083`, 500 | 2.4619e-02 | 5.0629e-01 | 5.3605e-02 | 2.7614e-02 |
| `kzz1e9/HeH0.083`, 500 | 2.9715e-03 | 5.4027e-01 | 6.8795e-02 | 8.0880e-03 |
| `wellmixed/HeH0.55`, 306 | 2.4767e-02 | 2.4767e-02 | not run | 2.4767e-02 |

**(A1) states the outer face and makes the row far worse, and the physics says
why.** With the copied ghost the gradient part of the flux vanishes and the
settling drift is what crosses the face; H2 is heavier than the mean particle
of an ionized hydrogen wind, so `G = (m_c - mbar) g/(kT)` is positive and the
drift is INWARD. MEASURED at cell 500 of `wellmixed/HeH0.083`, that inward
flux is -2.9471e-02 against an interior influx of +9.7857e-04, thirty times
the row's largest term, and the row measure goes to 5.06e-01. The linear
continuation is milder but still adds an outward flux of +8.3938e-04, nearly
as large as the influx it is supposed to pass on, and doubles the row measure.
On the Kzz case, where the eddy coefficient makes `A_grd` larger, the same two
give 182x and 23x. A stated outer diffusive face therefore is NOT a remedy;
what keeps the cell-500 row as small as it is, is that the face carries
nothing.

**(A2), excluding cell N from the carrier gate, is measured and it is not
enough, and the hydrodynamic rows exclude nothing.** READ from
`certification.f90`: `hydro_row_entry` fills `rr(j)` and `ss(j)` over
`do j = 1, N` and hands them to `certification_row_measure(N, ...)`, whose own
loop is `do j = 1, nc` with `nc = N`; the hydrodynamic mass, momentum and
energy rows therefore include cell N like any other cell, and the carrier
entry does the same. Dropping cell N from the carrier gate would leave
2.1713e-02, 3.5188e-04 and 2.4767e-02, which are 2170, 35 and 2477 times the
1.0e-05 wind tolerance: the refusal stands in all three.

**(A3) is a no-op wherever the ghost is a copy, and a real change where it is
data.** MEASURED in the new suite: with the outer ghost copied from cell N the
reconstruction already returns the donor cell average, so the intervention
removes nothing (0.0 exactly). The changes in the table above are therefore
entirely in the CERTIFICATION, whose ghost is data: forcing the donor average
there raises the cell-500 row from 2.4619e-02 to 2.7614e-02 (`HeH0.083`, 12
percent) and from 2.9715e-03 to 8.0880e-03 (Kzz, 172 percent). Read the other
way round: the row the relaxation's own trials balance at cell 500 is 1.12 and
2.7 times the row the certification refuses it on.

### 4. Two outer passes

The horizon step 2 used is no longer cheap on this tree. MEASURED on the same
reload of `wellmixed/HeH0.55`: the loaded state now reports an energy row of
1.0004 at cell 246 under the current restart contract ("the work state: the
conserved density, momentum and energy of the loaded state are held and p, T
follow the refreshed composition"), where the binary of step 2 read
3.4597e-08 for the same file. The stationary restart solve that follows is
therefore a real solve and not a formality, and the two-outer-pass horizon
costs hours rather than minutes. What was reached within the item, on
`wellmixed/HeH0.083`, is the certification at the head of outer pass 1, after
that restart solve and before any carrier relaxation (MEASURED: the carrier
densities of the four runs are identical there, largest relative difference
0.0, so the carrier operator had not yet acted):

| | frozen state | at the head of outer pass 1 |
|---|---|---|
| control | 2.4619e-02 at cell 500 | 2.6503e-02 at cell 500 |
| (A1) copy | 5.0629e-01 | 5.2811e-01 |
| (A1) extrap | 5.3605e-02 | 5.9571e-02 |
| (A3) donor cell | 2.7614e-02 | 2.9486e-02 |

The refusing cell does not move inward under any of them and the ordering is
unchanged. The full two-pass comparison of step 2 was NOT completed here and
is not claimed.

### 5. (B) The deferred reconstruction terms of the Jacobian

**Exactly which terms are lagged** (READ, `solve_carriers`): the advective
entries of the block-tridiagonal matrix are the FIRST-ORDER donor-cell
linearization of the face-flux divergence,
`d adv(j)/d f_c(don) = advj(j) msum(j)/msum(don)` on the diagonal or on the
neighbour according to the sign of the face mass flux, with the boundary
cases of the copied outer ghost and of the base face. The limited slopes of
`species_face_fraction` -- the PLM or WENO3 reconstruction and the scaling
that holds the face fraction on [0,1] -- are not differentiated. Nothing else
is lagged: the diffusive and drift entries are the exact derivatives
`carrier_face_flux` returns, and the chemistry block is a forward difference
of `carrier_source` itself.

**The action, measured.** `EXHALE_L22B_JAC_ACTION=<file>` with
`EXHALE_L22B_JAC_CELLS=300,312` writes, at the first Jacobian assembly of the
run, the action of the assembled matrix on a direction that is each cell's own
carrier fraction over cells 300-312 and zero elsewhere, beside a central
difference of the full residual along the same direction at a step of 1e-06,
with the outer ghosts filled by the rule the line search uses. On the frozen
`wellmixed/HeH0.55` state (MEASURED):

| | mean relative error over cells 299-313 | worst in that window | worst over the column |
|---|---|---|---|
| the assembly as it stands | 1.542e-03 | 6.709e-03 at cell 299 | 1.000 at cell 314 |
| with the reconstruction terms, banded | 3.463e-03 | 5.306e-03 at cell 302 | 1.000 at cell 314 |

Row by row at the edge of the direction's support, where the reconstruction
enters most (MEASURED): cells 299, 300, 301 read 6.71e-03, 1.88e-03, 4.88e-03
with the first-order entries and 9.6e-11, 5.8e-11, 2.7e-10 with the banded
derivative of the operator itself, so the banded entries ARE the operator's
own derivative there. In the interior of the support the comparison does not
improve, and the reason is the row the mean cannot show: **cell 314, where the
assembled action is exactly zero and the operator's is 1.679e-03.** The
direction stops at cell 312 and the reconstruction stencil of the face at
`r_{313+1/2}` reaches it, so the entry the action needs is
`d res(314)/d f_c(312)`, two cells off the diagonal, which a block-tridiagonal
matrix has no place for. Measured over the whole column, the largest advective
entry outside the band is 9.435e-02 of the band of its own row.

**And it is not what limits a pass.** MEASURED in the control relaxation of
outer pass 1 of the same state: every transport step ends with
"the residual reached the absolute floor" at a row residual between 1e-13 and
1e-17 after 2 to 5 Newton iterations, every interval is covered in one attempt
(the "interval was covered in N attempts" line is never printed), and all 22
refused trials of the pass are refused on the MOVEMENT BOUND. A Jacobian whose
action is right to 1.5e-03 is already far better than the iteration needs.

**And the pass with the reconstruction terms is the control pass, to every
digit.** MEASURED, outer pass 1 of the same reload with
`EXHALE_L22B_JAC_RECON=1` against the control: 17 kept transport steps in
both, 22 refused trials in both, a kept displacement of 2.904e-02 in both, the
same ending (the movement bound, with the last step inside it), the same
carrier steady residual afterwards (volume-weighted 7.01e-03, worst cell
2.62e-01 at cell 1), and the same probe ratios at r = 1.2007, 1.3578 and
1.6044 R_p (1.079, 1.275, 2.240) with the same kept intervals. Cell 306 does
not descend, and the Newton takes no fewer rejected trials, because none of
its trials was ever rejected on the Newton.

### 6. (C) The coupled failure after a converged relaxation

MEASURED with `EXHALE_L22B_DISPLACEMENT=1`, which reports what the composition
a pass hands back did to the primitive state the wind reads. Both runs
reproduce step 2's relaxation exactly (17 kept steps and a kept displacement
of 2.904e-02 in the control; 70 kept steps and 1.599e-01 with the far-wind set
waived), so the displacement below is of the same two compositions step 2
handed to the hydrodynamic solve:

| | control | far-wind set waived |
|---|---|---|
| kept transport steps | 17 | 70 |
| carrier displacement of the pass | 2.904e-02 | 1.599e-01 |
| max \|dp\|/p | 8.862e-03 at cell 230, r 1.2521 | 5.016e-02 at cell 10, r 1.0019 |
| max \|dT\|/T | 4.971e-03 at cell 275, r 1.5539 | 1.870e-02 at cell 1, r 1.0002 |
| max \|dmbar\|/mbar | 9.804e-03 at cell 217, r 1.2007 | 6.448e-02 at cell 3, r 1.0006 |
| max \|d(n_tot+n_e)\|/(n_tot+n_e) | 1.000e-02 at cell 221, r 1.2153 | 6.893e-02 at cell 3, r 1.0006 |

The cell that carries the displacement moves from the front (r = 1.20 to 1.55)
to the BASE (r = 1.0006 to 1.0019) and the displacement grows by factors 3.8 to
6.9. The last row of the table is the movement bound itself: in the control it
is attained exactly, 1.000e-02 against `trust` = 1.0e-02, and in the waived run
it stands at 6.9 times that value at cell 3.

**And that names what the far-wind experiment actually did.** The mask file of
that experiment omits cells 1 to 388, so the set whose particle count no longer
held a trial was not the far wind at all: it was everything BELOW 5 R_p, the
base and the front included. What the waive bought was a base free to move by
6.9 percent in one pass, and the coupled solve that had to follow it plateaued
with its worst mass row at cell 123, r = 1.038, in exactly that layer. The
carrier residual of 1.05e-09 that run reached is a statement about the carrier
block at a FIXED wind and says nothing about a wind that has to carry a base
which moved by seven percent.

READ from the step-2 log of that run
(`docs/audit_20260905/L22_step2_20260916/logs/wm055_far_waive.log`), the
hydrodynamic family that plateaued is the MASS row: from iteration 100 to the
205 at which the run was stopped, all 500 mass rows stand outside their
tolerance with the worst at **cell 123, r = 1.038**, the momentum row has 483
to 492 cells outside (worst at 496) and the energy row all 500 (worst at 59),
with the residual norm oscillating between 1.100 and 1.162 and the line-search
factor alternating 1.0 and 9.54e-07 at a pseudo-time step of 1.46e-04 and
7.32e-05. That is a limit cycle and not a slow descent.

The same log records what the pass did to the composition there (READ): the
base H2 partition fell from 0.934 to 0.784 over the 70 steps, and x(H2) at
r = 1.2007, 1.3578 and 1.6044 R_p was multiplied by 1.103, 1.855 and 8.833.

So the measurement rev3 step 3 (iii) asked for is this: the carrier block
closing at a fixed wind returned a composition whose pressure and temperature
in the inner layer are 5.0e-02 and 1.9e-02 away from the ones the
hydrodynamic rows were last solved at, and it is in that inner layer that the
coupled solve then failed to move. The measurement does not by itself say a
coupled block is required -- a shorter alternation with the same bound on the
inner layer would move the same way -- and this section does not choose one.

### 7. Byte identity with the gates off, and why the comparison has no control

`backup/regression/mol_carrier` was run on scratch copies, single-threaded
(`OMP_NUM_THREADS=1`, `EXHALE_MAXSTEPS=12000` from the case's `maxsteps`),
once with the tree binary `EXHALE.x` and once with `EXHALE_L22b.x`, with none
of the `EXHALE_L22B_*` keys set. **The four output files DIFFER** (MEASURED),
and the difference is NOT this item's: it cannot be attributed to this item's
file at all, and the comparison has no control.

What the comparison measured (MEASURED, largest relative movement over the 504
written rows, column by column): `Hydro_ioniz.txt` 1.98e-02 to 7.50e-02 across
every column it carries, `Ion_species.txt` 1.86e-02 to 5.85e-02 on the ATOMIC
columns H I, H II, He I, He II, He III and 1.26e-02 to 2.87e-02 on the
molecular carriers H2, H2+, H3+ and HeH+. And `Hydro_ioniz_adv.txt` differs in
its SCHEMA: the tree binary writes `# coupling:` and `# provenance:` headers
where the private build writes `# derived_from:`.

Three things say this is the tree moving under other workers and not this
item:

1. That header is written by `src/modules/files_IO/write_output.f90` (READ,
   the `# derived_from:` line), a file this item never opened and which
   another worker changed at 21:44 today. A private build compiles the LIVE
   tree, so it carries every such change.
2. The movement is LARGEST in the hydrodynamic columns and the atomic species
   and SMALLEST in the molecular carriers. Nothing in this item's file
   computes a hydrodynamic column or an atomic ionization stage.
3. This item's file cannot move a marching run at all with the keys off. Every
   line it adds to a path the marching loop takes is either inside a branch
   that is false (`l22b_outer_diff .ne. l22b_outer_diff_none`,
   `l22b_outer_adv_upwind`, `l22b_jac_recon`, `l22b_jac_action`), or a loop
   bound that takes its previous value (`jlast = N - 1`), or a declaration
   (`dadv`), or a call to `l22b_setup`, which reads the environment once and
   writes no state a residual reads. The relaxation, where the three
   unconditional array copies of the displacement report live, is not entered
   by a marching case (step 2 measured that `mol_carrier` never prints an
   `outer pass` line).

**What would settle it and was not available**: the gate the worker rules ask
for is a control binary built from the same live tree with this item's file at
its entry text, and this item has no copy of that text (the tree's own
`EXHALE.x` was built at 10:43, before several other workers' current edits).
The statement this section can make is therefore the code-path one above and
the measurement beside it, and not a byte comparison. No golden was refreshed
and no golden was read.

### 8. The tests

`src/tests/carrier_outer_boundary/` is new. Its driver states the outer
boundary of the carrier column and the deferred reconstruction terms as
measurements, on a synthetic twelve-cell column with the molecular chemistry
and the carrier transport on, and `run.sh` runs it once per setting of the
keys, which are read once per process. All rows pass (MEASURED):

| row | setting | measured |
|---|---|---|
| `the_outer_diffusive_face_carries_no_gradient_flux` | none | 0.0 of 0.0 |
| `the_outer_diffusive_face_carries_no_settling_drift` | none | 0.0 of 0.0 |
| `the_face_inside_it_carries_a_gradient_flux` | none | 1.490836E+09 > 0 |
| `the_outer_diffusive_face_is_formed` | copy, extrap | 1.261965E+09 > 0 |
| `the_copied_ghost_carries_cell_N_mixing_ratio` | copy | 6.030000E-01 of 6.030000E-01, 1e-14 |
| `the_extrapolated_ghost_continues_the_mixing_ratio` | extrap | 6.467819E-01 of 6.467819E-01, 1e-14 |
| `the_copied_outer_ghost_leaves_no_slope_on_the_outflow_face` | any | 0.0 within 1e-13 |
| `a_ghost_that_is_not_a_copy_leaves_a_slope_on_the_outflow_face` | not upwind | 3.620727E-02, at least 1e-3 |
| `the_donor_cell_outer_face_carries_no_slope` | upwind | 0.0 within 1e-13 |
| `the_first_order_entry_is_the_derivative_on_a_constant_profile` | none | 1.009449E+05 of 1.009449E+05, 1e-6 |
| `the_first_order_entry_misses_the_reconstruction_on_a_curved_profile` | none | 3.250776E-01, at least 1e-3 |
| `the_banded_derivative_is_the_operators_own` | none | 7.618035E+04 of 7.618035E+04, 1e-5 |

Every suite whose driver links `diffusive_photochemistry` was run with the
private build, single-threaded, into this item's own object directories
(MEASURED, exit 0 each, no FAIL line): `carrier_boundary_jacobian` 10,
`carrier_outer_boundary` 28, `carrier_retry` 152, `carrier_reference_scales`
14, `carrier_constraint_attribution` 10, `carrier_returned_state_acceptance`
36, `certification` 84, `steady_species_rows` 195, `species_face_flux` 1,
`element_operator` 34, `molecular_seed` 25.

Two rows of the new suite were written wrong and were corrected by
measurement and not by tolerance, and both corrections are results: a row
asserting that the first-order advective entry is the derivative on a LINEAR
mass fraction FAILED at 1.009449E+05 against 7.618035E+04, because the MC
limiter of PLM returns the donor value only on a CONSTANT one; and a row
asserting that the donor-cell intervention changes the outflow face FAILED at
exactly 0.0, which is the finding that a copied outer ghost already reduces
that face to the donor average.

The suite is a new one and its driver names three module entities this item
introduces (`l22b_setup`, `l22b_outer_ghost_state`,
`l22b_advective_band_derivative`), so it cannot be compiled against the entry
text; there is no RED run of it to quote and none is claimed. What the rows do
is state the operator's present behavior, including the two facts the plan's
wording did not have: that the outer diffusive face is not formed at all, and
that a copied outer ghost already reduces the outflow face to the donor cell
average.

### 9. Where the record is

`docs/audit_20260905/L22_step2b_20260916/`: `IDENTITIES.md` (the md5 of every
frozen state file reloaded, of the measuring binary and of the source, and the
table of keys), `TABLES.md` (the tables above as the analysis produced them),
`row_terms/` (the certification row dump of each frozen state under each
setting), `jacobian/` (the two Jacobian-action files) and `logs/`.

## Step 2c: one outer ghost rule, and the outer boundary of the diffusive and settling fluxes

Written 2026-09-17 (KST). Item L22 of `docs/PLAN_20260916_rev3.md` section 3
step 2c. What step 2b left was a defect of its own, reported there and not
fixed: the certification and the fixed-wind relaxation evaluated the carrier
row of one state through two different outer ghosts. This step states ONE
rule, applies it in one place, asserts it in a test, and measures what it
moves. Every number below is labeled MEASURED or READ, and every run is
`OMP_NUM_THREADS=1`.

### 1. The rule, and why it is the physical one

**(a) The ghost composition of the carrier faces is the interior continued,
X_ghost = X_N, in both directions of the face mass flux.** The outer face of
the domain is an outflow face of a transported scalar: its characteristics
leave the domain, so the boundary states no composition of its own. The
equilibrium composition the ionization sweep solves in the ghost cell is the
composition a cell of that density and temperature would relax to, not the
composition the wind carried out of cell N, and a transported quantity may not
take it as a boundary value.

Where the face mass flux reverses, `Frho(N) < 0` (item L7e sections 20-22,
R4), the rule does not change, and the reason is READ from the hydrodynamic
closure at the same face: `free_outflow_ghost` in `src/modules/states/Apply_BC.f90`
continues cell N's own state outward, the isothermal hydrostatic density at
cell N's temperature with the velocity copied, and states no reservoir at all.
Its own comment says so: an inflow at that face "would make the boundary an
INFLOW, carrying characteristics inward that only a reservoir may state, and
this closure has none". So the material an inflowing outer face returns is the
material the domain released, and the carrier rows are closed with the ghost
the hydrodynamic rows are closed with.

**(b) No diffusive, eddy or settling flux crosses the outer face**, which is
option (i) of the two the step named.

The justification is that the outermost cells are not collisional, so a
Chapman-Enskog transport coefficient there is not a fluid quantity. READ from
`docs/collisional_validity.md`: on LHS 1140 b the bulk Knudsen number is
0.13-0.30 at 8-9.5 `R_p`, 0.66-1.11 at 20 `R_p` and 1.4-2.4 at the domain top
of 30 `R_p`, with the exobase at 18.8-25.4 `R_p`. The composition still leaves
through that face by advection, which is a statement about the bulk motion and
not about the collisions.

Option (ii), forming the outer face from the ghost the way an interior face is
formed, is not defensible, and the measurement says why. With the stated
outflow ghost the gradient part of the flux vanishes identically, so only the
settling drift would cross the face; H2 is heavier than the mean particle of an
ionized hydrogen wind, `G = (m_c - mbar) g/(kT) > 0`, and the drift is INWARD.
The domain would gain carriers through its outer boundary out of nothing, at a
rate set by a coefficient that is undefined there. READ from step 2b section 3,
on the frozen `wellmixed/HeH0.083` state at cell 500, r = 29.0 `R_p`: that
inward flux is -2.9471e-02 cm^-3 s^-1 against an interior influx of
+9.7857e-04, thirty times the largest term of the row, and the row measure of
the cell goes from 2.4619e-02 to 5.0629e-01; the linear continuation of the
ghost gives 5.3605e-02, and on `kzz1e9/HeH0.083`, where the eddy coefficient
is larger, the two give 5.4027e-01 and 6.8795e-02 against a control of
2.9715e-03. The two variants that measured this are removed by this step, so
the delivered code carries one rule and no key that decides it.

**(c) The inner face carries the same statement, and it was already coded
this way.** `carrier_face_coefficients` forms no diffusive or settling
coefficient at `r_{1/2}` either, and what the layer below hands over is a
COMPOSITION on the inflowing base face, which `carrier_face_mass_fraction`
puts there: the handoff partition where a handoff states one, the base cell's
own where none does. A diffusive flux through the same face would state the
same handoff a second time. The two ends of the column are therefore one rule
and one closed box, stated at `carrier_face_coefficients`.

### 2. Where the rule lives now

Both ghost rules of the carrier column are applied in
`carrier_face_mass_fraction`, which every evaluation of the advective term
goes through: the stationary rows the certification measures, the trials of
the fixed-wind relaxation, and the derivative probes. The outer ghost of the
carrier state array itself is `carrier_outflow_ghost`, one routine that the
element budget limiter and the two derivative probes call in place of the
three inline copies they carried. The certification row and the relaxation row
are therefore the same function of the interior state by construction, and no
face of the column reads the ghost the ionization sweep leaves in the species
vector.

Removed with the rule: the two step 2b keys the unified rule makes
meaningless, `EXHALE_L22B_OUTER_DIFF` and `EXHALE_L22B_OUTER_ADV`, with
`l22b_outer_ghost_state` and the Jacobian branch that carried the experimental
outer diffusive entry. Kept, because they remain measurements and not
operators: `EXHALE_L22B_JAC_RECON`, `EXHALE_L22B_JAC_ACTION`,
`EXHALE_L22B_JAC_CELLS`, `EXHALE_L22B_DISPLACEMENT`, now documented in
`README_HOWTO.md` beside `EXHALE_CARRIER_ROW_TERMS`.

### 3. The test, RED and GREEN

`src/tests/carrier_outer_boundary/` states the rule as measurements. The row
that decides the step is
`one_outer_ghost_rule_gives_one_outflow_face`: the outward advective face term
of cell N is formed once from a column whose outer ghost is DATA, the kind of
ghost the stationary evaluation used to read out of the species vector, and
once from the same column with the ghost the relaxation's line search writes,
and the two are compared at tolerance ZERO.

MEASURED, against the entry text of `diffusive_photochemistry.f90` (a driver
compiled against it, with the ghost fill written in the driver, since the new
module routine does not exist there): **FAIL at 2.1724e-02**, the face mass
fraction of a column whose donor cell average is 3.0e-01, so the two
evaluations of one row put face values 7.2 percent apart. MEASURED with the
delivered code: **PASS at exactly 0.0**, with
`the_outflow_face_is_the_donor_cell_average` at 0.0 within 1e-13 beside it.

All sixteen rows of the suite pass (MEASURED): the four rows that state the
two end faces carry neither gradient flux nor settling drift, the row that
states the face inside them does (1.4908e+09), the two rows above, and the
Jacobian rows of step 2b unchanged (1.0094e+05 against 1.0094e+05 on the
constant profile, 3.2508e-01 of deferred reconstruction on the curved one,
7.6180e+04 against 7.6180e+04 for the banded derivative).

`src/tests/steady_species_rows/` carries the row that compares the advective
term of a carrier row with the operator the Runge-Kutta stages apply. It
FAILED at 6.0175e-01 with the delivered code and PASSES at 7.5340e-14
(MEASURED), and the correction is a result and not a tolerance: that test
already wrote the INNER ghost of the stage state by hand, "which is what makes
the two comparable", and now writes the outer ghost the same way. What it
compares is the two forms of the term on one state under one boundary rule.

### 4. What it moves

**The control is a second build of ONE snapshot of the tree**, because the
tree moved under this item while it ran: `System_HeH_mol.f90`,
`ionization_equilibrium.f90`, `util_ion_eq.f90`, `steady_newton.f90`,
`write_output.f90` and `input_read.f90` were all written by other workers
between 00:19 and 00:22 on 2026-09-17 (MEASURED from their modification
times), and a first control built at 23:52 gave differences of up to a factor
of 1.4 in the H2 PRODUCTION column at cells 1-3 that no boundary rule can
produce. The two binaries below are built from one copy of the tree that
differs in this item's file alone.

| item | value |
|---|---|
| measured build | `EXHALE_A.x`, md5 `ae7050e059b102c8a0ec8436572229d7` |
| control build, the same snapshot at the entry text of this item's file | `EXHALE_B.x`, md5 `62d046cc76a9cf140841f5514cb4a353` |
| this item's file before the edit | md5 `4f9687fc1c64d045454ee5c71c77ce3e` |

**The three frozen states: bitwise identical, and the reason is a
measurement.** MEASURED, `Restart intent: stationary evaluate` with
`EXHALE_CARRIER_ROW_TERMS=1`: `output/carrier_row_terms.txt` is byte-identical
between the two builds on `wellmixed/HeH0.083`, on `kzz1e9/HeH0.083` and on
`kzz1e9/HeH2.13`, and so is every certification line.

| case | worst gated carrier row, control | measured | refusing cell |
|---|---|---|---|
| `wellmixed/HeH0.083` | 2.462e-02 | 2.462e-02 | 500, unchanged |
| `kzz1e9/HeH0.083` | 2.972e-03 | 2.972e-03 | 500, unchanged |
| `kzz1e9/HeH2.13` | 7.817e-01 | 7.817e-01 | 243, unchanged |

The refusing cell does not move inward, and the row at cells 495-500 of each
state is unchanged to every digit printed.

**Why the rows do not move, MEASURED and not assumed.** A temporary probe in
`carrier_steady_residual` (built privately, removed before delivery) printed
the carrier state of cell N and of its two outer ghosts on the reload of
`wellmixed/HeH0.083`: `f_c(N) = f_c(N+1) = f_c(N+2) = 6.309600E-03` and
`Frho(N) = +1.371665E-10`. The outer ghost of the species vector is ALREADY
cell N's own value on these states, and the face is an outflow face, so the
two evaluations happened to agree and the unification changes no number on
them. The defect step 2b reported is therefore real as a property of the CODE
and not observable on these states today; what this step removes is the
dependence of a certified carrier row on an upstream accident, and the test of
section 3 is what holds it.

**The molecular regression cases.** MEASURED on scratch copies, single
threaded, `EXHALE_MAXSTEPS` from each case's `maxsteps`, the two builds of the
one snapshot: `mol_carrier`, `mol_base_handoff` and `mol_diffusion` each ran 12000 steps with both builds, and every output file of every case is byte-identical except the `# provenance:` line, which carries the wall-clock time of the run: `Hydro_ioniz.txt`, `Ion_species.txt`, their `_adv` twins, `Cooling_breakdown.txt`, `Heating_breakdown.txt` and `IC_dump.txt` all show zero differing lines with that header excluded (MEASURED). The goldens of these cases therefore do not move, and none was refreshed or read.

**The suites.** Every suite whose driver links `diffusive_photochemistry` was
run against the delivered tree, single-threaded, into this item's own object
directories (MEASURED, exit 0 each, no FAIL line): `carrier_boundary_jacobian`
10, `carrier_outer_boundary` 16, `carrier_retry` 152,
`carrier_reference_scales` 14, `carrier_constraint_attribution` 10,
`carrier_returned_state_acceptance` 36, `certification` 84,
`steady_species_rows` 195, `species_face_flux` 1, `element_operator` 34,
`molecular_seed` 25.

### 5. What this step found and did not fix

`carrier_mass_fractions` in `src/modules/functions/binary_element_diffusion.f90`
is the stage path's spelling of `carrier_face_mass_fraction`: it states the
INNER ghost rule of the carrier column and states no outer one, so the
Runge-Kutta stages read whatever the species vector holds in the outer ghosts.
That is the same defect this step removed on the stationary side, in a file
this item may not edit. It is not observable today for the reason section 4
gives, the species vector's outer ghost being the outflow copy on the states
measured, but the rule belongs to the equation and not to one evaluation of
it. Reported, not fixed.

No golden was refreshed and no golden was read.

---

## Step 3, increments I1 and I2

Written 2026-09-17 (KST). The design is
`docs/lhs1140b_stationary_L22_step3_design_20260916.md`, approved by the user
on 2026-09-17 with all nine decisions D1 to D9 as recommended; this section
records the two increments that decision opened, I1 (the coupling
measurement, a diagnostic that changes no solve) and I2 (the three-valued key
and the handover). I3 to I5 are gated on I1 and are not run here.

The source identities the measurements were made against (MEASURED,
`md5sum`): `src/modules/time_step/steady_newton.f90`
`2f0ca95722aaa7ed3eae7f20469c66fb` after this increment,
`src/modules/lower_atmosphere/diffusive_photochemistry.f90`
`f04ed90a1f1b34bdc619ede52793df6f`, which is the tree WITH step 2c landed.
The control binary is a private build of the tree as this increment found it.

### 1. I1, the coupling of the block, measured

`EXHALE_COUPLED_JAC_ACTION=<file>` with `EXHALE_COUPLED_JAC_CELLS=lo,hi` and
`EXHALE_COUPLED_JAC_DIR=mass|momentum|energy|carrier` writes, on the state
`solve_steady_jfnk` is entered at, row by row for every unknown of the block:
the action the solve uses (`jacobian_action_of_direction`, hence `jv_product`,
a forward difference of the full residual at the probe length the route sets),
a CENTRAL difference of the same residual at that length and at a tenth of it,
and the row of the banded preconditioner `build_banded_jac_full` on the same
direction. It then stops the run, having taken no step. The direction is the
named unknown's own column scale (`cell_state_scales`) in the cells named and
zero elsewhere. With the key unset the routine returns at its first statement.

The key also registers the block's rows for the measurement
(`coupled_block_jacobian_action_requested`), which is what lets it be made on
a state the alternation wrote: `carrier_newton` is a restart token a
`Restart option change:` line may not name, so a state written under
`Coupled carrier solve: False` cannot be reloaded under `True` at all
(MEASURED, the reload is refused with "The state in the file is a state of
another equation set"). That is also why `On stall` leaves `carrier_in_newton`
false.

**MEASURED on the frozen `molecular_scalar_gj1132_wellmixed/HeH0.55` state**
(its `input.inp`, its `base.inp` and its `output/`, spectrum path absolute,
`Restart intent: stationary`, single-threaded, `nvar_jac` = 4, `nspec_row` = 1,
`kl_jac` = `ku_jac` = 11, `ncolor_jac` = 23, N = 500), direction on cells 300
to 312:

| | carrier direction | energy direction |
|---|---|---|
| probe length the route uses | 1.8313e-06 | 2.0955e-02 |
| fraction to the boundary of the species box | 1.3796e+01 | no component points out of the box |
| species components with no room at any step length | 0 | 0 |
| the action was sampled | yes | yes |
| unresolved preconditioner colors | 0 | 0 |
| mean relative error, action against the central difference, over the rows the direction moves | 3.692e-04 (40 rows) | 3.792e-03 (21 rows) |
| worst of those | 2.513e-03 | 1.259e-02 |
| the same over EVERY row of the support | 1.365e-01, worst 1.004 | 7.150e-03, worst 7.685e-02 |

**The two means differ because of the rows the direction does not move.** The
mass row of a cell is a statement about the total density at a fixed conserved
state and does not move with the composition at all: at cell 306 the two sides
of the carrier direction's mass row are 1.1e-13 and 4.6e-16 against a momentum
row of 2.4e-10 in the same cell, so the quotient is a ratio of two round-off
numbers and reads 1 whatever the action is. A row counts as moved when its
response reaches 1e-3 of the largest response in the support; both numbers are
written and the gate is on the moved rows.

**The verdict I1 was opened for: the difference quotient is usable.** The
action reproduces a central difference of the same map to 3.7e-04 and 3.8e-03
in the mean, against the 1.542e-03 that step 2b measured for the carrier
operator's own assembled action and accepted. Neither of the two ways the
coupled route was set aside before appears on this state: no species component
is blocked, the fraction-to-the-boundary length is 13.8 times the direction
(the probe uses at most 0.9 of the room to a face), and every color of the
preconditioner resolved.

**The cross terms are both live, and of a stated size.** The smallest energy
row over the carrier direction's cells is 4.472e-09 in code units (the
composition reaching the energy row through the pressure and the temperature
at the conserved state, and through the heating and the cooling of the sweep
the evaluation ran); the smallest carrier row over the energy direction's
cells is 1.721e-11 (the hydrodynamic state reaching the carrier row through
the advective face flux, the temperature of every rate, the density of the
three-body terms, the diffusivities and the photolysis).

**The band holds the two-cell reconstruction entry.** `d res_c(j+2)/d f_c(j)`,
the entry the carrier relaxation's block-tridiagonal matrix has no place for
and at which step 2b measured an assembled action of exactly zero against an
operator action of 1.679e-03, stands in the block at a flat distance of
`2*nvar_jac` = 8 inside `kl_jac` = 11. All thirteen of those entries over the
direction's cells are nonzero, at 1.0012e-12 to 1.2725e-12, falling smoothly
outward. The block needs no special treatment for the reconstruction stencil.

**The species rows are on the certification's scale.** The carrier row of
`eval_residual`, divided by `tscale_code = R0/(v0*n0)` = 1.7452e-10, is the
row `carrier_steady_residual` returns on the state that evaluation produced,
bitwise, at all thirteen cells.

**The outer ghost rule is the same in both evaluations.** The design expected
this row to FAIL while two rules were in force. It PASSES: at cells 499 and
500 the two evaluations agree bitwise (6.433052715436995e-05 and
6.245054851935448e-05). That is step 2c, which landed on 2026-09-17 in the
tree these measurements were made on; the row is the standing statement that
it stays landed, and no tolerance was widened to obtain it.

The suite is `src/tests/coupled_block_jacobian/` (`run.sh` + `README.md`),
eight rows, GREEN on the increment's binary and not stateable against the
entry text, which has no such key: the control binary ignores it, writes no
file and goes on solving (MEASURED; the suite bounds each invocation with
`EXHALE_JAC_TIMEOUT`, so a binary that does not stop is a failed row rather
than a suite that never returns).

### 2. I2, the three-valued key and the handover

`Coupled carrier solve:` takes `False` (the alternation, the default and
unchanged), `True` (the block from the first pass, unchanged) and `On stall`
(new). A word the key has no meaning for now STOPS the run: the entry-text
parser read word 4 alone and tested it against `True`, so `On stall` fell
through to `False` silently, which leaves an input file asking for the block
and a run that never enters it (MEASURED on the control binary: the run starts
and `EXHALE_resolved.out` reads `carrier_in_newton F` with no other record of
the line). `On stall` reaches the same block one pass later, so the
eliminated-H2 refusal of item B5i holds for it as well.

**The handover.** In `steady_wind_with_element_diffusion` the route of a pass
is a local, `block_now`, which is `carrier_in_newton` at entry. At the end of
a pass the state is handed to the block when all of the following hold, which
is the condition of section 6.1 of the design:

- the outer loop's own ending for this pass is `outer_no_progress`, that is,
  the pass did not count as progress under the joint rule of L14 and that has
  now happened on `outer_no_fall_max` = 3 consecutive passes, which is the
  ending the loop REFUSES the state on today; AND
- the carrier relaxation of this pass ended on the composition movement bound
  (`carrier_relax_movement_bound`); AND
- so did the relaxations of the two passes before it (three consecutive,
  counted in `n_bound_endings`, reset by any other ending and by a pass in
  which no relaxation runs).

The ending is then a handover instead of a refusal: `block_now` becomes true,
`outer_ending` returns to `outer_running`, the no-fall counter is cleared and
the remaining passes register the carrier and element rows in the Newton.
The log prints the reason, the two counts, the movement bound that is no
longer read and the entry that refuses the state. The acceptance is untouched:
the certification of the refreshed state, entry by entry, against the same
tolerances, exactly as before.

**With the key absent or `False` every path is the alternation, bit for bit.**
MEASURED, control binary against the increment's, single-threaded, on scratch
copies: `backup/regression/lower_profile` (`EXHALE_MAXSTEPS` 12000 from its
`maxsteps`), `backup/regression/wasp_full_newton`, and the certified
`molecular_scalar_gj1132_kzz1e9/HeH2.13` reloaded at `EXHALE_OUTER_PASSES=2`:
the control is a build of the SAME live tree with this increment's four files, and only those, put back to their entry text (the tree moved under this item while it ran: another worker landed L22 step 2c in `diffusive_photochemistry.f90` between the first control build and the measured one, which by itself moved the marching `mol_carrier` trajectory in the seventh digit, so the first control was no control at all). Against the reverted-source control, single-threaded: `backup/regression/mol_carrier` (marching, 1500 steps, the case the design names for the I1 key with the key unset) reproduces EVERY DATA LINE of all six output files, the four that carry a provenance or source header differing in that header alone (the control is built out of tree and stamps `git=unknown tree=clean`); the certified `molecular_scalar_gj1132_kzz1e9/HeH2.13` reload on the stationary alternation route prints an IDENTICAL log, line for line and digit for digit, over 152 lines including the whole JFNK trace, the one differing line being the same provenance statement, and its output files are bit-identical; `backup/regression/lower_profile` (marching, 1500 steps, then its stationary route) reproduces every data line of all ten output files, the four provenance-carrying headers again the only difference, and its whole 1716-line log differs in the two provenance lines alone; `backup/regression/wasp_full_newton` and the certified reload were stopped before their own ends, and over the 841 and 152 lines they had printed the two logs differ in 0 and 2 lines, the 2 being that same provenance statement, so their whole printed trajectory, JFNK iterations included, is digit for digit the control's

`src/tests/grid_and_gates/coupled_carrier_h2_row.sh` grew three stages and
five rows (the third value accepted and echoed as itself, an unknown value
refused with the three values named, `On stall` on an eliminated-H2 molecular
configuration refused, and the handover sentence present in the binary).
MEASURED: four of the five RED against the control binary, all ten rows GREEN
against the increment's.

### 3. The exploratory `On stall` solve

Run at 8 threads with a pass cap of 5 on the frozen
`wellmixed/HeH0.55` state (D8 permits one exploratory solve at 8 threads).
The handover cannot fire before pass 4 by construction, and after 25 minutes
the run stood at iteration 38 of the hydrodynamic solve of PASS 1
(`||R||` 1.636, merit 4.39e-03, line-search factor 0.25, `dtau` 1.02). It was
stopped there. **So no handover was observed, and none was expected in the
time available; whether the criterion fires on this state is I4's
measurement and is not claimed here.**

### 4. What this step found and did not fix

`load_IC.f90` records the restart option token `carrier_newton` from
`carrier_in_newton` alone, which `On stall` leaves false. A state written
AFTER a handover is therefore written with `carrier_newton: F` although the
block produced it. The route that produced a state belongs in `cert_reason`
under the L23 metadata contract, which is that item's to carry; the token
itself is a layout statement about which equation set the file holds, and
after a handover the equation set of the written state is the block's. Named
here, in a file this increment may not edit.

## Step 3, increments I3 to I5

Written 2026-09-17 (KST). The design is
`docs/lhs1140b_stationary_L22_step3_design_20260916.md`; I1 and I2 are the
section above. Increments I3 (T2, the block on a state that already has an
answer), I4 (T3, the three refusing cases under the handover) and I5 (the
record) follow here, with the one code change they needed first.

The identities the measurements were made against are in
`docs/audit_20260905/L22_step3_20260916/identities.txt`: the private binary
`EXHALE_L22e.x` and a control built out of tree from the same snapshot with
this increment's two source files, and only those, at their entry text.

### 1. The one code change: a route is not an equation set

I1 found, and I2 reported in a file it could not edit, that the reload of an
alternation state under `Coupled carrier solve: True` was REFUSED, with "The
state in the file is a state of another equation set", so the block could not
be entered on a reload at all; and that a state written after an `On stall`
handover recorded `carrier_newton=F` although the block produced it. Both are
one mistake about what the `carrier_newton` token of a state file's
`# options` line says. The whole argument and what changed are
`docs/lhs1140b_stationary_L23_20260916.md` section 6c; in short:

- `carrier_newton` adds no row and no column. The five tokens it stood beside
  (`metals`, `mol`, `oxychem`, `carrier`, `iontrans`) each decide how many
  unknowns a state has; this one decides only whether the transported balances
  are unknowns of the Newton vector or are relaxed at a held wind. The
  balances are the same balances and the certification evaluates the same rows
  against the same tolerances, so under the L23 metadata contract the route is
  metadata OF a state and a difference in it is admissible with nothing named.
- `opt_is_route` is added beside `opt_changes_layout` in `load_IC.f90`, true
  for this token alone; the refusal stands unchanged for the tokens that DO
  change the equation set, `iontrans` among them.
- The load prints `the restart changes the ROUTE and not the equations` and
  writes a `# route_change` line into both halves of the state it produces,
  inherited by the rungs after it.
- The written token is now `carrier_in_newton .or. carrier_rows_entered_newton`,
  the second set at the handover in `EXHALE_main.f90`, so a state states the
  route that produced it.

**The test.** Five rows in `src/tests/grid_and_gates/restart_option_change.sh`
on `backup/regression/carrier_model_a_newton`, a molecular case with the H2
carrier transported and `Coupled carrier solve: False`, hence an alternation
state; every run takes `Restart intent: stationary evaluate`, so what is under
test is the load. MEASURED: four rows RED against the control binary and all
twenty rows of the suite GREEN against the private one
(`restart_option_change_control.log` and `_increment.log` in the record). The
fifth row, the same state reloaded under `On stall`, passes on both binaries
and is stated as such: `On stall` leaves the token at `F` until a handover
fires, so before one there is nothing to admit.

**What the reload of a certified state now measures on entry.** The states
this step starts from were written by the campaign and are NOT roots of
today's equations, which is the impact section of
`docs/lhs1140b_stationary_L7g_model_20260916.md`. The entry evaluation states
it in numbers (MEASURED, the same certification the gate uses, on the state as
loaded):

| case, the state reloaded | hydrodynamic mass row | hydrodynamic energy row | worst gated carrier row |
|---|---|---|---|
| `kzz1e9/HeH2.13`, the CERTIFIED state | within | 4.267e-01 of 1.0e-06 at cell 203 | 7.817e-01 of 1.0e-05 at cell 243 |
| `wellmixed/HeH0.083` | 5.791e-08 of 1.8e-08 at cell 1 | 1.525e-01 of 1.0e-06 at cell 282 | 2.462e-02 of 1.0e-05 at cell 500 |
| `wellmixed/HeH0.55` | 6.640e-08 of 2.2e-08 at cell 1 | 1.000e+00 of 1.0e-06 at cell 246 | 2.364e-02 of 1.0e-05 at cell 306 |
| `kzz1e9/HeH0.083` | 3.766e-08 of 1.6e-08 at cell 1 | 3.305e-01 of 1.0e-06 at cell 411 | 2.972e-03 of 1.0e-05 at cell 500 |

The energy row is five to six decades outside its tolerance on every one of
them, so what these solves are handed is not a state near a root whose carrier
row alone refuses it. **That changes what T2 and T3 can say** and it is stated
before the results: the block is asked to solve a problem the design did not
pose it, and a case that does not stall on the composition movement bound
cannot reach the handover at all.

### 2. I3 (T2): the block against the alternation on the certified case

`molecular_scalar_gj1132_kzz1e9/HeH2.13`, its `REPRODUCE.md` recipe
(`Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=40`), the certified state copied to a scratch directory
(`LHS1140b/models/.L22/`, not the catalog), both runs started together on
`lart3` at `OMP_NUM_THREADS=8`, one with `Coupled carrier solve: True` (the
block) and one with `False` (the alternation, the control on today's physics).
The reload under `True` is possible only because of the change of section 1.

**Verdict: the alternation converges on this state and the block does not.**
MEASURED over the same 2 h 45 min of wall clock:

| | the alternation (`False`) | the block (`True`) |
|---|---|---|
| what the route did | 6 outer passes | 7 trust-region iterations of one solve |
| worst gated carrier row, pass by pass | 7.82e-01, 5.30e-01, 1.49e-01, 2.28e-02, 2.49e-03, 4.72e-04 (cells 243, 221, 221, 222, 242, 221) | one iterate; the row stands at 1.317 times its tolerance window at cell 243 |
| `\|\|R\|\|` | 4.27e-01 at entry, 1.15e-03 inside pass 6 | 4.267e-01 at entry, 4.248e-01 at iteration 7 |
| hydrodynamic `info` | 2 at pass 1, then 0 at passes 2 to 6 | the solve has not returned |
| Krylov cycles that exhausted their subspace | 0 of 202 printed iterations | 7 of 7 |
| species unknowns on a bound | none printed | 0 |

**The alternation CERTIFIED at outer pass 12** (`info` = 0, "every active
equation of this state is within its own tolerance"), its worst gated carrier
row reaching 8.18e-06 of 1.0e-05 at cell 221 after 1.42e-05 at pass 11 and
2.16e-05 at pass 10. The block, from the same state and over a longer wall
clock, stood at trust-region iteration 14 with `||R||` = 4.241e-01 against the
4.267e-01 it started from, 0.6 per cent, and **every one of its thirteen linear
solves ended with the subspace exhausted**, 40 products of 40 with the relative
residual at 2.597e-01 against the 1.00e-01 asked for. That is the L24 failure
mode, which the design's risk table names and asks T3 to report; it appears
here, on the case the design expected to be the easy one. The species box is
not what holds it: 0 unknowns on a bound, and I1 measured the
fraction-to-the-boundary length at 13.8 times the direction.

**The cost, per the design's section 8.1.** At `nvar_jac` = 4 and N = 500 the
block pays `ncolor_jac` + 1 = 24 full residual evaluations for the
preconditioner and one per Krylov product, and it took 40 products at every
iteration, so an iteration cost about 64 full residual evaluations, each one an
`Apply_BC`, an `excited_H_update`, an `ioniz_eq` sweep over 500 cells and a
`carrier_steady_residual`. The alternation's hydrodynamic solve carries three
unknowns per cell and took 1 or 2 products per iteration on this state. The
wall clocks are quoted from a host carrying other work at a load average of
109 to 128 on 72 cores, so they are not a clean cost measure and the counted
evaluations are; the two routes ran side by side under the same load.

**The certified state beside the archived one** (review row R37: reported
beside the gate and not as the gate, the archived state being no oracle of a
changed operator). The post-processing pass and `EXHALE_transit.py` were run
on the state the alternation certified, by the case's own `REPRODUCE.md`
recipe, with the WINERED HIRES-Y kernel:

| | this item, certified at pass 12 | the archived catalog state |
|---|---|---|
| base T [K] at r = 1.000193 R_p | 778.79 | 808.32 |
| x2 at cell 1 | 3.3536e-01 | 3.5530e-01 |
| outermost radius with x2 >= 1e-2 (the H2 front) | 1.1597 R_p (cell 204) | 1.1626 R_p (cell 205) |
| log10 Mdot [g/s] | 7.9024 | 7.9091 |
| He I 10830 red-pair equivalent width [%A] over a 1.6 A window | 1.5599 | 1.5754 |
| red-pair depth [%] | 5.5712 | 5.6282 |
| FWHM [A] | 0.2628 | 0.2628 |

The archived state's numbers are reproduced here from its own files by the
same tool and agree with the ones its `REPRODUCE.md` quotes (1.5758 %A,
5.6282 %, 0.2628 A), so the two columns are one measurement of two states and
not two measurements. **The re-solve on today's physics moves the observable
by about one per cent**: the equivalent width by -1.0 per cent, the mass-loss
rate by -1.5 per cent, the front by one cell, and the base temperature by
-3.7 per cent, which is the largest of them and is where the L7g corrections
land. The block produced no state to put in a third column.

Both halves of the comparison carry `adv_input_certified F` in the transit
header: the post-processing pass rewrites the state through `Do only PP` and
that write states no stationary claim, which is a property of that route and
not of the state the solve certified (the certified state is kept beside it as
`output/*_IC.txt`).

### 3. I4 (T3): the three refusing cases under `On stall`

Each from the state the campaign left, its own recipe with
`EXHALE_PTC_DTAU0=1.0e8` (the continuation value its `REPRODUCE.md` states),
`EXHALE_OUTER_PASSES=40`, the row-kind ledger on, on `lart3`.

**No handover fired in any of the three, and in two of them the reason is that
the case is not stalling from this state.** MEASURED, worst gated carrier row
pass by pass:

| case | passes reached | worst gated carrier row, pass by pass | hydrodynamic `info` |
|---|---|---|---|
| `wellmixed/HeH0.083` | 2 | 2.65e-02 (cell 500), 3.85e-02 (cell 500) | 0, 2 |
| `wellmixed/HeH0.55` | 1 | 5.92e-02 (cell 315) | 0 |
| `kzz1e9/HeH0.083` | 5 | 2.97e-03, 2.52e-02, 7.98e-03, 3.35e-03, 4.10e-03 (cells 500, 500, 500, 499, 500) | 2, 0, 0, 0, 0 |

These three ran first on `lart3`, which carried other work at a load average of
109 to 128 on 72 cores, and were restarted from the same states on the idle
`lart4` at 04:30 KST on the user's instruction that a new long solve goes
there; their `lart3` logs are kept beside them as `run_lart3.log`. The pass
trajectory above is the same on both hosts pass for pass, which is what L15
stage 2 says it should be; the wall clock per pass fell by a factor of about
8.5 (1161 s to 136 s on the first pass of `wellmixed/HeH0.083`, 1821 s to
213 s on `kzz1e9/HeH0.083`).

The handover asks for `outer_no_progress` (the joint rule of L14 failing on
three consecutive passes) AND three consecutive carrier relaxations ending on
the composition movement bound. `kzz1e9/HeH0.083` has run five passes without
either count reaching three, and its carrier row is a decade below where the
campaign's forty passes left it. The two `wellmixed` cases are earlier than
the earliest pass at which a handover can fire by construction, which is
pass 4. So **the criterion was not exercised** by these runs in the time
available: that is a statement about the runs and not about the criterion, and
it is not a claim that the criterion does not fire.

No named risk of the design's section 8.3 was met in these three: 0 species
unknowns on a bound, 0 Krylov cycles that exhausted their subspace, no
positivity or box event and no chemical-root refusal in any of the three
printed trajectories. The one risk that WAS met is in I3, and it is named
there.

### 4. I5: the record, the fixtures, the suites

The record is `docs/audit_20260905/L22_step3_20260916/` (its `README.md` lists
what each file holds): the identities, the RED and GREEN suite logs, the suite
table, the fixture comparison and the five run logs.

**The fixtures, with the new key value absent.** `carrier_elem_newton`,
`carrier_model_a_newton` and `atomic_elem_newton`, each run on a scratch copy
with the private binary and with the control, single-threaded, both started
together and both stopped after 25 minutes; what is compared is the whole
printed trajectory over the lines both had printed, the method I2 used.
MEASURED: `carrier_elem_newton` 1036 lines compared, 0 differing;
`atomic_elem_newton` 4782 lines, 0 differing; `carrier_model_a_newton` 3310
lines, 40 differing, every one of them an outer-pass line differing in its
trailing wall-clock seconds field alone, and 0 differing with that field cut.
**Every number the two binaries printed is the same number.** No golden was
refreshed and no golden is proposed for refresh.

**The suites of the design's section 7**, run with the private objects
`build_L22e` and the private binary: `krylov_and_dogleg` 334,
`attempted_step` 70, `steady_species_rows` 195, `certification` 84,
`carrier_retry` 152, `carrier_returned_state_acceptance` 36,
`carrier_reference_scales` 14, `carrier_constraint_attribution` 10,
`carrier_boundary_jacobian` 10, `carrier_outer_boundary` 16,
`element_operator` 34, `molecular_seed` 25, all PASS with 0 FAIL; the whole
`grid_and_gates` suite 250 PASS. `steady_completion_flag` and
`steady_selfconsistent_residual` are log readers; given the log of the
certified `i3_alt` solve through `EXHALE_STEADY_RUNLOGS` they give 3 PASS and
2 PASS with one FAIL, `handback_matches_the_accepted_iterate`, 2.83103e-08
against a reference of 1 within a factor of 2. Its input is a log and not a
binary, so the same log gives the same verdict whichever build produced it and
the row says nothing about this increment; whether that number is a property
of a solve that converged this far or a defect was not established here.

### 5. What this step found and did not fix

**One `grid_and_gates` row fails on the control binary as well and is not this
increment's.** `output_state_consistency.sh`,
`outer_iteration_ending_is_the_stagnation_one`, expects the outer loop to end
its bounded run on `no_progress` and measures `pass_budget`. It runs with
`Coupled carrier solve: False`, so no route branch and no handover is reached.
MEASURED on the control binary, which carries this increment's two files at
their entry text: the same FAIL, the same measured value. It is the ending of
the outer loop having moved under the row, and the file that states it is not
one this item may edit.

**The two questions I3 and I4 were opened to answer are not answered by these
runs**, and the reason is a property of the states and not of the runs: on
today's equations none of the four is a state whose carrier row alone refuses
it, the energy row being five to six decades out on every one of them. What
these measurements do settle is narrower and is stated as such: the block,
entered from such a state, is held by its linear solve (7 of 7 subspaces
exhausted) while the alternation from the same state falls three decades in
six passes. That supports decision D3, the alternation first and the block on
stall, and it argues against `True` as a route to be entered from far away.
