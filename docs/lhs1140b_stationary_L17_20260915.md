# L17: why the 0.02-XUV hydrodynamic solve of the LHS 1140 b wind stalls

Item L17 of `docs/PLAN_20260913_lhs_stationary.md`, opened by the user on
2026-09-15 16:50 on the evidence L4h section 5.3 handed it. Every number below
is MEASURED on this tree at ONE thread unless it is marked READ. The file
changed is `src/modules/time_step/steady_newton.f90` and nothing else.

## 1. Verdict

**The refusal is neither the finite-difference Jacobian action nor the
certification anchor. It is the Newton model of a wind that is subsonic to
Mach 1e-8 over its whole lower atmosphere.**

- **The anchor is not the cause, and it is generous, not tight.** At the cell
  that binds, the continuity row stands at 1.6189e-04 of its own flux against
  a tolerance of 3.0e-07, and the row's MEASURED response to one ulp of every
  cell density there is 7.65e-14. The row is therefore 5.3e+02 times its
  anchored tolerance and 2.1e+09 ulps of its own arithmetic away from where
  the certification asks it to be: a resolvable imbalance, not a rounding
  floor. The base velocity at 0.02 is the base velocity of the CERTIFIED 0.03
  state to three digits (-1.28048e-01 against -1.27875e-01 cm/s, Mach 8.595e-07
  against 8.589e-07), so nothing about the base collapsed between the two rungs.
- **The finite-difference action is rounding-limited at the arc the Krylov
  cycle uses, and that is not the cause either.** The additivity defect of the
  first two basis directions is 1.510e-02 and scales as the arc to the power
  -0.976, which is rounding (N27's measurement reproduced). Multiplying the
  probe arc by a thousand takes the defect to 1.243e-03 and the exponent to
  +0.750, i.e. twelve times a better operator -- and the solve then stalls in
  the same place with the same rows.
- **What fails is the linear model and the step built on it.** The Krylov cycle
  reaches its forcing term and leaves the rows that bind at 1.7 and 21 times
  the values they entered the step with; asked for ten times more (200-vector
  subspace, tolerance reached in 50 products) it leaves them at 3.0 and 24.
  The scaled trust region, made available on this route by this item, predicts
  a three per cent fall of the squared merit at its first radius and measures
  the merit nearly DOUBLING (ratio -26.8), cuts the radius twelve times to
  2.3e-12, accepts no step at all, and hands back the state it was given.
- **Five step controls, twenty outer iterations each, leave the state where it
  was** (section 5). The merit of this state is 153 times the merit of the
  certified 0.03 state beside it and no variant moves it by 30 per cent.

**The case is NOT solved.** The chain to 0.015 and 0.01 was therefore not run
and the three catalogue cases at 0.01 XUV are NOT certified.

**What the tree gained** is the step control the hydrodynamic half of the
partitioned route did not have (section 6), default off and byte-identical with
the option off. It is delivered because it is what produced the measurement in
the third bullet, not because it repairs the case: it does not.

## 2. The fixture

`LHS1140b/models/.L17/` throughout; `LHS1140b/models/.L4h/` and `.L14/` were
read and not written to.

- **The stalled state**: `LHS1140b/models/.L14/x002_HeH2.13/output/{Hydro_ioniz,
  Ion_species}.txt`, the state item L14's own 0.02 run handed back
  (`certified=F`, `cert_reason=no_stationary_claim`). It is the L4h stall in a
  cheaper form: the same rung, the same rows, the same cell, and its elemental
  row is already inside its tolerance (7.60e-08 gated against 1.0e-05), so on
  this state ONLY the two hydrodynamic rows refuse and the Picard alternation
  is not in the way. The L4h run of the same rung ran to outer pass 7 and was
  stopped in its eighth (2026-09-15 17:41; it wrote no completion line and what
  ended it is NOT established here, so its 40 passes were never taken). Over
  the seven it did take, its hydrodynamic rows are flat from pass 3 while its
  elemental row falls by a decade and a half:

  | pass | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  |---|---|---|---|---|---|---|---|
  | mass | 2.41e-04 | 2.73e-04 | 2.11e-04 | 1.87e-04 | 1.98e-04 | 1.93e-04 | 1.93e-04 |
  | energy | 6.95e-04 | 8.62e-04 | 6.77e-04 | 4.47e-04 | 5.81e-04 | 5.91e-04 | 5.99e-04 |
  | elemental, gated | 2.95e-02 | 2.14e-02 | 1.44e-02 | 9.06e-03 | 5.40e-03 | 3.06e-03 | 1.81e-03 |

  the same division this memo is about. `.L4h/x002_HeH9.7` was stopped in its
  third pass, its mass row 4.30e-05 then 4.56e-05 and its energy row 4.64e-05
  then 3.86e-05.
- **The certified state beside it**: `.L14/x003_HeH2.13/output/*` (0.03 XUV,
  `certified=T`).
- **Control binary** `EXHALE_L17_ctl.x`, md5 `9a944cdf573541a47a50ec9599703c07`
  -- the same md5 as `EXHALE_L4h.x`, i.e. the tree exactly as L17 found it.
  **Measured binary** `EXHALE_L17.x`, md5 `922dec0fa183c2611d1bfb33c218a6cd`.
  Both built as `make OBJDIR=build_L17 EXE=...`; the tree's `EXHALE.x` and
  `build/` were not touched.
- Every comparison at `OMP_NUM_THREADS=1` (L4h section 5.0: this solve is not
  reproducible at eight).

## 3. The rows of the two states, and the anchor

`Restart intent: stationary evaluate` with `EXHALE_MASS_FLOOR_SCAN=1`.

| at cell 1 (r = 1.00019) | 0.02, stalled | 0.03, certified |
|---|---|---|
| continuity row `\|R_1\|/s_1` | **1.6189e-04** | 2.4467e-09 |
| the row's step under one ulp of every cell density | 7.6543e-14 | 1.4537e-08 |
| `mass_row_rounding_floor` (the estimate) | 3.0263e-08 | 2.0028e-08 |
| `cert_tol_mass_at` = 10 x that floor | 3.0e-07 | 2.0e-07 |
| Mach | 8.595e-07 | 8.589e-07 |
| cell-centred v [cm/s] | -1.28048e-01 | -1.27875e-01 |

The rows of the whole state as loaded:

| row | 0.02, stalled | 0.03, certified |
|---|---|---|
| mass | 1.619e-04 at cell 1, tol 3.0e-07, **ABOVE** | 3.444e-08 at cell 20, within |
| momentum | 8.045e-12 at cell 492, tol 1e-08, within | 4.495e-13 at cell 19, within |
| energy | 5.216e-04 at cell 1, tol 1e-06, **ABOVE** | 1.623e-06 at cell 20, ABOVE |
| elemental transport He/H | 7.602e-08 gated, tol 1e-05, within | 6.703e-06 gated, within |
| `\|\|R\|\|` / merit `\|\|Fs\|\|2` | 5.216e-04 / 1.58e-08 | 1.623e-06 / 1.03e-10 |

Two things follow at once. The anchor stands **six decades above** the
measured rounding of the cell it anchors, so it is not what holds the 0.02
row; and the 0.02 state is 153 times the certified state in the merit the
solve descends, so the solve is not at a root of anything either.

The energy row's terms at the cell that binds, `write_residual_breakdown`:
`R = -1.7186e-09`, flux divergence `3.264262e-06`, heat `3.294767e-06`,
cool `2.878645e-08` -- the row is a cancellation of its largest term by
2.6e+02, which is what puts its rounding where it is.

## 4. The residual is a state function, and the operator is rounding-limited

- `EXHALE_RESID_DETERMINISM=1` on the stalled state: **"the residual is a state
  function (bitwise identical)"**, 0 of 1500 entries differ after three trials,
  after a discarded trust-region trial and after a Jacobian-vector probe. The
  stall is not evaluation noise between calls.
- `EXHALE_JV_ADDITIVITY=1` at outer iteration 1:

| directions | probe arc | relative additivity defect | defect ~ arc^p |
|---|---|---|---|
| first two of the Arnoldi basis | 8.004e-12 | **1.510e-02** | p = -0.976 (rounding) |
| two off the basis | 3.720e-07 | 6.495e-03 | p = +0.943 (curvature) |
| first two of the basis, `EXHALE_JV_PROBE_ARC=1.0e3` | 8.075e-09 | **1.243e-03** | p = +0.750 |

  Nine tenths of the squared defect sits in 33 of 500 cells between cell 1 and
  cell 69 in the first row of the table and in cells 1 and 2 in the second; the
  second difference of the residual along the short arc has ratios 1.06 and
  1.33 where a residual smooth on that arc gives 4 and 4, i.e. the short arc is
  at its rounding floor. This is N27's finding on this state. **It is not what
  holds the case**: the third row of the table is a twelve times better
  operator and its solve is the entry text's (section 5).

## 5. The linear solve, the step, and five step controls

`EXHALE_LINEAR_ROWS=1` separates what the cycle left from what the step did.
At outer iteration 1 of the stalled state:

| | forcing term 1e-1, what the cycle reached | linear residual left, over the row of F (mass, momentum, energy) | rows after the step |
|---|---|---|---|
| entry text, 40-vector subspace | 4.013e-01 in 40 of 40, subspace exhausted | **1.680, 1.500e+01, 2.126e+01** | 1.619e-04, 8.042e-12, 5.208e-04 at `lam` 7.63e-06 |
| `EXHALE_GM_M=200` | 9.203e-02 in 50 products, tolerance reached | **2.963, 1.554e+01, 2.359e+01** | unchanged at `lam` 1.22e-04 |

A cycle that reduces the 2-norm it minimizes by a factor eleven leaves the
three rows the certification reads at three to twenty-four times their own
values. The norm the cycle minimizes -- the state column scales, since this
route's `Drow` is `D` -- does not contain the rows that bind.

Twenty outer iterations from the same state, one thread, `EXHALE_PTC_DTAU0=1.0e8`,
`EXHALE_OUTER_PASSES=1`:

| step control | `\|\|R\|\|` at iteration 20 | merit at 20 | state handed back | how it ended |
|---|---|---|---|---|
| entry text | 4.045e-04 | 1.43e-08 | 5.216e-04 | the iteration cap |
| `EXHALE_JUDGED_ROWS=1` | 1.765e-04 | 1.46e-08 | 5.318e-04 | cap; `dtau` collapsed to 4e-02 |
| `EXHALE_JUDGED_ROWS=1 EXHALE_KRYLOV_TOL_ABS=1` | 5.940e-04 | 1.08e-08 | 5.318e-04 | cap |
| `EXHALE_JV_PROBE_ARC=1.0e3` | 3.704e-04 | 1.34e-08 | 5.228e-04 | cap; `dtau` collapsed to 4e-01 |
| `EXHALE_TRUST_REGION=1` (new) | 5.216e-04 | 1.58e-08 | 5.216e-04 | **no step accepted at all**, aborted at 12 |

The state handed back is the state each was given, to within the first
iteration's own evaluation. For comparison, the certified 0.03 state re-entered
the same way is at `||R||` 1.623e-06 and merit 1.03e-10, and ONE Newton
iteration at the same `dtau0` takes it to 9.271e-07 and 2.49e-11 with the
acceptance gate met.

**What the trust region says, and it is the sharpest statement in this memo.**
Its first radius is measured on the Krylov step, 3.873e-05 against a step of
4.130e-04; at that radius the banded model predicts a reduction of the squared
merit of `pred` 4.177e-18 -- three per cent of the 1.25e-16 the state carries --
and the measured change is `actual` -1.120e-16, the merit nearly doubling,
`ratio` -26.804. The region is cut twelve times, to 2.308e-12, and every trial
is refused; below `||s||` ~ 1e-09 the measured change stops shrinking and
changes sign from one radius to the next, which is its own resolution. **The
banded model of this state is wrong by a factor 27 at a tenth of the Newton
step.**

## 6. What the physics of this column is, and where the model loses it

The wind is subsonic everywhere, and by eight decades at its base. Mach from
the state:

| cell | 1 | 2 | 3 | 4 | 5 | 6 | 50 | 200 | 300 | 400 | 500 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| r [R_p] | 1.0002 | 1.0004 | 1.0006 | 1.0008 | 1.0010 | 1.0012 | 1.0097 | 1.1489 | 1.8571 | 5.9043 | 29.03 |
| Mach | 8.59e-07 | 1.02e-06 | 1.40e-07 | 8.86e-08 | 3.63e-08 | 2.92e-09 | 2.73e-08 | 3.17e-06 | 5.30e-05 | 1.46e-03 | 9.78e-03 |
| sign of v | - | - | + | - | + | - | + | + | + | + | + |

Two consequences, both measured above and both stated here because the next
item needs them.

1. **The rows the certification reads and the norm the linear cycle minimizes
   are different functionals on this column, and by the reciprocal Mach
   number.** The continuity row is judged against the flux it differences
   (`mass_flux_row_scale`), the merit and the linear model divide it by the
   density (`cell_state_scales`), and those stand 1e+06 to 1e+08 apart here.
   `judged_row_scaling_on` repairs the second for the LINEAR MODEL only, and
   section 5 measures that this is not enough on this state; the merit form of
   the same repair is what section 7 delivers, and the trust region that needs
   it rejects every step, so the repair is available and is not the cure.
2. **The HLLC flux selects its branch on the SIGN of the contact speed**
   (`Num_Fluxes.f90`: `SL < 0 and S_star >= 0` takes the left star flux,
   `S_star < 0 and SR >= 0` the right), and `S_star ~ v` sits at 1e-09 to
   1e-02 of the wave span over this whole domain and changes sign from one
   base cell to the next. The mass component of the star flux is
   `rho* S_star` with `rho*_L` and `rho*_R` differing by the density jump
   across the face, 4 per cent per cell at this base: the flux is continuous
   at `S_star = 0` and its derivative is not, and every base cell sits on that
   kink. The measured signature is in section 4: nine tenths of the second
   difference of the residual along a probe direction in cells 1 and 2, with
   ratios 1.06 and 1.33 where a smooth residual gives 4 and 4.

**A C1 blend of the two star fluxes was implemented, measured and REVERTED.**
Over `|S_star| < b (S_R - S_L)` with a cubic Hermite weight, at `b = 1e-3`:
the band covers **409 of 500 cells, out to 6.74 R_p**, and the residual of the
same state moves from (1.619e-04, 8.045e-12, 5.216e-04) to
(3.309e-03, 1.284e-07, 1.370e-03). That is not a regularization of this
discretization, it is a different one. The reason a narrower band does not
help is in the algebra: inside the band the flux is perturbed by
`(1 - w)(rho*_R - rho*_L)/rho*_L`, about half the density jump across the
face, **independently of b** -- about two per cent of the face mass flux, which
in a stationary state IS the mass-loss rate. Only the NUMBER of faces shrinks
with b. `src/modules/flux/Num_Fluxes.f90`, `src/modules/init/parameters.f90`
and `src/modules/files_IO/input_read.f90` were restored to the tree's state
and the delivered binary was rebuilt from them.

## 7. What changed in the source

All of it in `src/modules/time_step/steady_newton.f90`, all of it default off,
and the delivered binary with no option set is **bit-identical in behaviour to
the control** (section 8).

1. **`merit_reads_the_model_row_scaling`** (`EXHALE_JUDGED_MERIT=1`) and the
   single reader **`merit_is_on_the_model_row_scaling()`**. Wherever the system
   carries a species row the merit is already read on `Drow`, the row scaling
   of the linear model, and `cell_row_scales` makes that the certification's
   own, so descending the merit is descending the gate. On the three-unknown
   route the merit is `|| F/D ||_2` on the state column scales while the gate
   is `max_k,j |F| / (tol_k s_k)`: two functionals, and section 5 measures how
   far apart they are here. The new switch reads the merit on `Drow` there as
   well. Every `|| F/. ||_2` in the module now asks one question, through that
   function: the initial merit, the true merit of each iterate, `Dfix` of the
   stagnation window, the line-search trial, the Levenberg-Marquardt gradient
   and its trials, and the damped Gauss-Newton escape. With
   `judged_row_scaling_on` off, `Drow` IS `D` on this route and every one of
   those is the number it was.
2. **`trust_region_on_the_hydrodynamic_rows`** (`EXHALE_TRUST_REGION=1`). The
   scaled trust region was the step control of the coupled route only:
   `use_tr = (nspec_row > 0)`, and the environment variable could only turn it
   OFF, so the hydrodynamic half of the partitioned route had no step control
   but the pseudo-transient backtracking line search -- whose failure,
   "no descent direction exists for the banded model at this state", is the one
   the trust region was MADE the default for on the other route (the
   measurement is at the declaration of `use_tr`). `=1` now makes it the step
   control here too. **It carries the merit with it**: the region's model
   residual is `F/Drow` and its acceptance compares the reduction it predicts
   with the reduction of the merit, so the two must be one functional; the
   switch therefore sets `merit_reads_the_model_row_scaling` as well.
3. `use_tr` is now assigned AFTER `read_species_unknown_space_controls`, which
   is where the switches it depends on are read. It was assigned before, which
   with (2) would have read the previous solve's value.

Delivered because it produced the measurement of section 5 -- the banded
model's misprediction -- which the route could not state before. **It does not
solve this case** and it is not proposed as a default.

## 8. Verification

- **Control against measured, option off, one thread.** Same state, same
  environment (`EXHALE_PTC_DTAU0=1.0e8 EXHALE_OUTER_PASSES=1
  EXHALE_JFNK_MAXIT=12`), `EXHALE_L17_ctl.x` against `EXHALE_L17.x`: every
  `(JFNK)` line identical to every digit over the twelve iterations compared
  (`diff` of the two logs' JFNK lines: empty).
- **Suites**, run from the item's own object directory `build_L17/` with
  `EXHALE_SPECIES_EXE` pointed at the measured binary:

| suite | PASS | FAIL |
|---|---|---|
| `krylov_and_dogleg` | 334 | 0 |
| `attempted_step` | 70 | 0 |
| `certification` | 84 | 0 |
| `steady_species_rows` | 199 | 0 |

- The regression matrix was not run: no case of it carries this route's new
  switches, both default off, and the control/measured identity above is the
  statement that the operator did not move.

## 9. A separate finding, reported and not acted on

**The certified 0.03 state does not reproduce its certification through its own
file.** Written at outer pass 19 of `.L14/x003_HeH2.13` with the mass row at
3.764e-08 (cell 7) and the energy row at 8.255e-07 (cell 7), both within, it
is re-entered here with `Restart intent: stationary evaluate` and measures the
mass row at 3.444e-08 (cell 20, within) and the energy row at **1.623e-06
(cell 20), above its 1.0e-06**; the run exits 2 and the state is refused.

It is not the digits and it is not the composition: the file carries 17
significant figures and the loaded composition reproduces the sweep's own root
to 1.782e-14. What the round trip rebuilds is the radiative terms, and the
energy row at this base cancels its largest term by 2.6e+02, so a part in 1e+03
of the heating is the whole of the row. **Every rung of the 0.03 -> 0.02 ->
0.015 -> 0.01 ladder is seeded from such a file**, so the state each rung is
handed is a state whose energy row is already at one to two times the
tolerance the rung is asked to reach.

## 10. Reproduce

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd $EX && make OBJDIR=build_L17 EXE=EXHALE_L17.x          # md5 922dec0f...
M=$EX/LHS1140b/models;  L=$M/.L17

# the fixture
mkdir -p $L/diag_x002/output
cp $M/.L14/x002_HeH2.13/input.inp                $L/diag_x002/
cp $M/.L14/x002_HeH2.13/output/Hydro_ioniz.txt   $L/diag_x002/output/Hydro_ioniz_IC.txt
cp $M/.L14/x002_HeH2.13/output/Ion_species.txt   $L/diag_x002/output/Ion_species_IC.txt

# section 3, the rows and the rounding floor (add "evaluate" to the intent)
sed -i 's/^Restart intent:.*/Restart intent: stationary evaluate/' input.inp
OMP_NUM_THREADS=1 EXHALE_MASS_FLOOR_SCAN=1 $EX/EXHALE_L17.x

# section 4, the operator; section 5, the linear solve and the step
OMP_NUM_THREADS=1 EXHALE_PTC_DTAU0=1.0e8 EXHALE_OUTER_PASSES=1 \
  EXHALE_JFNK_MAXIT=20 EXHALE_LINEAR_ROWS=1 EXHALE_JV_ADDITIVITY=1 \
  $EX/EXHALE_L17.x
OMP_NUM_THREADS=1 EXHALE_RESID_DETERMINISM=1 $EX/EXHALE_L17.x

# section 5, the trust region on this route
OMP_NUM_THREADS=1 EXHALE_PTC_DTAU0=1.0e8 EXHALE_OUTER_PASSES=1 \
  EXHALE_JFNK_MAXIT=40 EXHALE_TRUST_REGION=1 $EX/EXHALE_L17.x
```

## 11. What the next item has to decide

Stated as questions the measurements above pose, not as a plan.

1. The linear cycle minimizes a 2-norm that does not contain the rows the state
   is judged by, and making it contain them (`EXHALE_JUDGED_ROWS=1`, with and
   without an absolute target) does not move this state. The remaining freedom
   is the COLUMN scaling and the preconditioner, neither of which has been
   measured on a column this subsonic.
2. The banded model mispredicts the merit by a factor 27 at a tenth of the
   Newton step. The band omits the non-local radiative coupling, and at this
   base the heating of a cell is set by the column above it; that is a dense
   row, and the energy row it feeds cancels by 2.6e+02.
3. The contact branch of the HLLC flux is not differentiable and the whole
   column sits on it. A blend of the two star fluxes is not the instrument
   (section 6); what would be is a contact whose star density is not chosen by
   a sign -- a question about the DISCRETIZATION of a Mach 1e-8 wind, to be
   settled before it is coded.
4. Whether a rung should be seeded from a file at all, given section 9.
