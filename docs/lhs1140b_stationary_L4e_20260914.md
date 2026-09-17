# L4e: the linear solve and the line search of the partitioned route judged by the rows the certification judges

Item L4e of `docs/PLAN_20260913_lhs_stationary.md`, the proposals P1, P2 and
P4 that `docs/lhs1140b_stationary_L4d_20260913.md` section 5 and
`docs/lhs1140b_stationary_L4f_20260913.md` section 8 left open. Every number
is MEASURED on this tree unless it is marked READ. The file changed is
`src/modules/time_step/steady_newton.f90` and nothing else.

## 1. Verdict

**P1 is implemented and measured, and it is NOT the repair: the continuity
excursion is the Newton direction of these states, not the tolerance or the
row scaling of the linear solve.** With the rows of the linear system and of
the band that preconditions it divided by the size each row is judged at,
`tol_k(j) s_k(j)`, the linear residual left in the continuity row of the
LHS 1140 b C/N/O pass-2 state falls from 2.608e+05 of that row to 4.334e+03,
and with the cycle additionally asked for a residual INSIDE the tolerances it
falls to 1.584e+02 at a relative residual of 8.179e-10 reached in four
products. The step that direction produces still raises the continuity row
from 5.002e-09 to 1.093e-02, within seven percent of the 1.169e-02 the entry
text's loose cycle produced, and the twelve-iteration ladder is the entry
text's to three digits in every row. An essentially exact solve of the same
system gives essentially the same step, so what moves the continuity row is
the direction and not the residual left beside it.

**What the obstruction IS, measured: the pseudo-transient start that
`Restart intent: stationary` fixes at `dtau0 = 1.0`.** The two states L4d and
L4f measured as insoluble are solved, by the ENTRY TEXT, in five and six
Newton iterations at `dtau0 = 1e8`:

| state | `dtau0` | entry text | with the judged rows |
|---|---|---|---|
| pass-2 entry (L4d's reproduction) | 1.0 | STAGNATED at 40, `info = 0` at 98 after L4d's detector (READ, L4d section 1) | ladder identical to the entry text's, 12 iterations |
| pass-2 entry | 1e8 | `info = 0` at iteration 6, `\|\|R\|\|` 1.121e-01 -> 3.674e-08 | `info = 0` at iteration 5, 2.718e-08 |
| pass-8 entry (L4f's reproduction) | 1.0 | `info = 1` at 60, THE ENTRY STATE RETURNED, 58 non-monotone accepts | `info = 2` at 6, the same state returned |
| pass-8 entry | 1e8 | `info = 0` at iteration 5, 1.032e-03 -> 6.972e-08 | `info = 0` at iteration 4, 2.718e-08 |

`lam = 1` at every iteration of all four `dtau0 = 1e8` solves and not one
non-monotone accept. The grind L4f took for a property of the outer boundary
cell, and the continuity excursion L4d took for a property of the linear
tolerance, are both the pseudo-transient shift: at `dtau0 = 1` the system the
step solves is `(J + I/dtau) dY = -F` with the shift comparable to the
operator, and the sequence it generates is a relaxation, not a Newton
iteration.

**Both cases the item ends on then certify, with the ENTRY TEXT binary and
that one variable** (section 7): `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`,
resumed at `dtau0 = 1e8` from the state its first seven passes produce,
is ACCEPTED at outer pass 1 in 22.8 s -- "every active equation of this state
is within its own tolerance", written `certified=T` -- where the same state
at `dtau0 = 1` spends 783.6 s on a pass that ends `info = 2` with the energy
row exactly where it started; and `atomic_scalar_gj1132_kzz1e9/HeH4.0`,
resumed from the state its REFUSED four-pass run ends on, is ACCEPTED at
outer pass 3 in 93 s with `omega` never halved. The first is a certified
LHS 1140 b solution carrying the C, N and O reservoirs, which no item of this
plan had reached.

**P2 is not implemented.** The brief's condition for it was that P1 be
insufficient, which it is; but the measurement above says the line search is
not where this is decided either -- at `dtau0 = 1e8` the SAME line search
accepts `lam = 1` at every iteration and the solve certifies. Changing the
acceptance functional would be a change to the step control of a route whose
obstruction has now been located elsewhere, and it is not made here.

**P4 is implemented**: the count of cells outside each hydrodynamic row's
tolerance, with the cell furthest outside, is printed beside the judged-rows
line at every outer iteration.

## 2. What changed in the source

All in `src/modules/time_step/steady_newton.f90`. The file stands at +589/-108
against `43bc28c`, of which L4d's change is +276/-36 (READ, L4d section 2);
this item's increment is the six blocks below.

1. **`judged_size_of_the_hydrodynamic_rows(u, Dj)`**, new. `Dj(k,j) =
   tol_k(j) max(s_k(j), cert_scale_floor)` with `s_k = residual_row_scale`
   and `tol_k` the certification's own: `cert_tol_mass_at(j,u)` cell by cell,
   `cert_tol_momentum`, `cert_tol_energy`. `|F|/Dj` is then the distance
   `hydrodynamic_distance_from_certification_by_cell` takes the maximum of.
2. **`row_scaling_of_the_linear_model(Y, D, Dr, u)`**, new, the one statement
   of what the rows of the Newton system are divided by: the certification
   row scales where a species row is carried (`cell_row_scales`, unchanged),
   the judged sizes on the three-unknown route, and `D` where
   `judged_row_scaling_on` is off. Both call sites of `cell_row_scales` in
   `solve_steady_jfnk` now go through it.
3. **The three-unknown branch of the band build** divides its rows by `Drow`
   instead of by `D` and carries the pseudo-transient term as
   `idtau*(D(jc)/Drow(jc))`, which is the same diagonal `E = D/Drow` on both
   sides of the system, so the step it solves for is unchanged and only the
   norm the truncated cycle minimizes moves. With the option off `D/Drow` is
   one and the branch is the entry text, bitwise (section 6.1).
4. **The band that is factorized is equilibrated on that route too**
   (`equilibrate_the_preconditioner`, already used on the four-unknown route
   and already on by default) when the judged scaling is on. The reason is
   that partial pivoting compares entries ACROSS rows, and rows carrying
   1/tol_k(j) stand twelve decades apart, so the pivot sequence would be
   chosen by the scaling and not by the operator. MEASURED on the pass-2
   state it changes nothing there: the cycle reaches 2.849e-03 in two
   products with it and without it, and leaves 4.334e+03 against 4.379e+03
   of the mass row (`EXHALE_PRECON_EQUIL=0`). It is kept because the
   argument is about the factorization and not about this state, and it
   costs one pass over the band.
5. **`pgmres` and the four routines that form the same product** lost the
   three-unknown copy of the action; `w = w/Drow + idtau*(D/Drow)*z` is now
   the one expression, and it is the old one wherever `Drow` holds `D`.
   `levenberg_marquardt_descent` forms the gradient of the merit IT
   descends in the coordinates `ab` is stored in, `A^T((Drow/D) F/D)`, since
   the model's rows are no longer the merit's.
6. **P4**: `certified_row_measures` gains two optional outputs, the count of
   cells outside each hydrodynamic row's tolerance and the cell furthest
   outside; `eval_residual` carries them through its existing optional
   arguments; and `solve_steady_jfnk` prints one line beside the judged-rows
   line at every outer iteration:

   ```
   (JFNK) cells outside the tolerance of their row, of 500: mass 0 (worst cell 1), momentum 0 (1), energy 500 (245)
   ```

   The counts come from the quantities that routine has already formed, so
   nothing is evaluated twice, and the continuity row's verdict is
   `mass_row_cell_verdict`'s, the one expression of that rule. The first
   form of this print took its own sweep of `residual_row_scale` and
   `mass_row_rounding_floor`, and that MOVED the carrier reload at the
   1e-11 level (section 6.3); a diagnostic that changes the answer is not a
   diagnostic.

Two options, both read once per solve by
`read_species_unknown_space_controls`, **both off by default, so the
delivered binary is the entry text bitwise on both routes** (section 6):

- `judged_row_scaling_on`, `EXHALE_JUDGED_ROWS=1` turns it on. Off because
  it does not repair what it was proposed for and costs one outer pass on
  the element reload (sections 4 and 6.2).
- `krylov_stops_inside_the_tolerances`, `EXHALE_KRYLOV_TOL_ABS=1` turns it
  on, and it acts only with the judged rows: the cycle is then asked for
  `eta = 1/||F/Drow||_2`, a linear residual whose 2-norm is below one, which
  in the judged rows means inside every row's tolerance in every cell. Off
  because 1/||b|| is 1e-10 on these states, below the accuracy of a
  finite-difference Jacobian action, and because it buys nothing (section 4).

## 3. RED: the three reproductions, on the entry text

The control build is `43bc28c` plus the working tree's `steady_newton.f90`,
which is L4d's change and is not in HEAD; the measured build differs from it
in that file alone. Both were built in an isolated copy of the tree, because
seven source files were being edited by other workers while this item ran
(section 7). 8 threads, `EXHALE_PTC_DTAU0=1.0`.

### (a) The pass-8 entry state of `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`, 60 iterations

`info = 1`, THE ENTRY STATE RETURNED, `||R||` 1.032e-03 in and out, 58
non-monotone accepts, 35 accepted steps that raised the judged rows above
their five-iterate window, best judged distance 1.032e+03. Row by row, the
linear residual left over the same row of `F`:

| iteration | mass | momentum | energy |
|---|---|---|---|
| 1 | 1.139e+03 | 4.443e+03 | 3.613e-01 |
| 60 | 2.034e-01 | 5.130e+02 | 5.632e-01 |

Every one of those is the number L4f section 4 records, to every printed
digit. The P4 print on the same state: mass 0 cells outside (worst cell 1),
momentum 0 (1), energy 499 of 500 (worst cell 157) -- which is the count L4f
had to measure by hand. After the first step it reads mass 500 (worst cell
258), momentum 269 (496), energy 498 (174): the step that leaves the judged
rows where they were has moved the continuity row outside its tolerance in
every cell of the column, which no line of the entry text says.

### (b) The pass-2 entry state, the L4d reproduction

| iteration | rows entering (mass, momentum, energy) | linear residual left, over the same row of `F` | rows leaving |
|---|---|---|---|
| 1 | 5.002e-09, 4.183e-13, 1.121e-01 | **2.608e+05, 1.996e+05, 4.447e-01** | 1.169e-02, 3.119e-03, 6.486e-02 |

and the judged distance goes 1.121e+05 -> 1.112e+09 across that step: the
first Newton step raises the continuity row by seven decades and the state's
distance from certification by four. Again L4d section 5 to every digit.

### (c) `atomic_scalar_gj1132_kzz1e9/HeH4.0` from its own certified seed (L11)

Outer pass 1 `hydro info = 0` (the seed is in the basin), worst gated species
row 8.93e-05 of 1.0e-05 at cell 292, mass 5.63e-10, momentum 1.12e-14, energy
1.18e-08, 250.7 s. Outer pass 2's hydrodynamic solve then grinds: run with
no iteration cap it stood at iteration 306, `||R||` 4.033e-03 falling by
parts in a thousand, worst row the energy row of cell 500, judged distance
5.51e+05, after 1 h 40 min, and was stopped. Capped at 120 iterations the
whole loop runs to its end and reproduces L11:

| pass | hydro `info` | worst gated species row | energy row | `omega` |
|---|---|---|---|---|
| 1 | 0 | 8.93e-05 at cell 292 | 1.18e-08 | 0.500 |
| 2 | 1 | 4.80e-05 at cell 294 | 5.36e-03 | 0.250 |
| 3 | 1 | 3.68e-05 at cell 295 | 6.73e-03 | 0.125 |
| 4 | 1 | 3.25e-05 at cell 295 | 7.25e-03 | 0.125 |

and pass 4 is REFUSED, "the joint distance of the state has not fallen in 3
consecutive passes". L11 records 3.3e-05, 7.2e-03 at cell 500 and `omega`
0.125 at pass 4 (READ), which is this run. So the HeH4.0 stall is the same
entry: a hydrodynamic solve that cannot reach its tolerance at `dtau0 = 1`,
and an element relaxation that then halves its step because the joint
distance does not fall.

## 4. GREEN: what the judged rows and the absolute target do

Same three states, same seeds, measured build.

### The linear residual, first Newton step

| state | quantity | entry text | judged rows | judged rows + absolute target |
|---|---|---|---|---|
| pass-2 | linear residual left in the mass row, over that row of `F` | 2.608e+05 | 4.334e+03 | 1.584e+02 |
| pass-2 | the same, momentum | 1.996e+05 | 3.348e+03 | 3.677e+01 |
| pass-2 | the same, energy | 4.447e-01 | 7.467e-03 | 2.453e-04 |
| pass-2 | what the cycle reached, in how many products | 7.900e-04, 1 | 2.849e-03, 2 | 8.179e-10, 4 |
| pass-8 | linear residual left in the mass row | 1.139e+03 | 1.601e+01 | 4.538e-01 |
| pass-8 | what the cycle reached, in how many products | 1.048e-03, 1 | 2.242e-03, 2 | 1.969e-06, 3 |

So the row scaling does what it was asked for: two decades off the continuity
row's linear residual, four with the absolute target, and the cycle's relative
tolerance now means the same thing in rows whose tolerances stand ten decades
apart.

### And what the step does with it

| state | quantity | entry text | judged rows | judged rows + absolute target |
|---|---|---|---|---|
| pass-2 | mass row after the first step | 1.169e-02 | 1.095e-02 | 1.093e-02 |
| pass-2 | mass row after 12 steps | 8.930e-03 | 8.934e-03 | 8.932e-03 |
| pass-2 | energy row after 12 steps | 5.112e-02 | 5.111e-02 | 5.111e-02 |
| pass-8 | mass row after the first step | 4.383e-05 | 4.090e-05 | 4.086e-05 |
| pass-8 | outcome of the 60-iteration solve | `info = 1` at 60, entry state | `info = 2` at 6, entry state | `info = 2` at 6, entry state |

**That is the refutation.** At a linear residual of 8.2e-10 of the judged
right-hand side -- four products, four decades inside anything the entry text
ever reached -- the step raises the continuity row to 1.093e-02 where the
entry text's 7.9e-04 raised it to 1.169e-02. The twelve-iteration ladder of
the three builds agrees to three digits in every row of every iteration. The
excursion is first order in a direction that is nearly the same direction in
all three, and no tolerance on the linear residual reaches it.

## 5. What the obstruction is: the pseudo-transient start

The discriminating measurement, and it is one line of environment: the same
two states, the same binaries, `EXHALE_PTC_DTAU0=1.0e8` in place of the 1.0
the stationary restart fixes.

| state | build | iterations to `info = 0` | `\|\|R\|\|` reached | `lam` | non-monotone accepts |
|---|---|---|---|---|---|
| pass-2 entry | entry text | 6 | 3.674e-08 | 1.0 throughout | 0 |
| pass-2 entry | judged rows | 5 | 2.718e-08 | 1.0 throughout | 0 |
| pass-8 entry | entry text | 5 | 6.972e-08 | 1.0 throughout | 0 |
| pass-8 entry | judged rows | 4 | 2.718e-08 | 1.0 throughout | 0 |

Every one of the four ends on the loop-top stop, "the acceptance gate is met
and every one of the 3 certified rows of the system this solve carries is
within its own tolerance". The state L4f measured as a sixty-iteration grind
that returns what it was given is four Newton iterations from a certified
hydrodynamic solution; the state L4d measured as needing 98 iterations after
a repaired stagnation detector is six.

The continuity excursion is still there at `dtau0 = 1e8` -- the first step
raises the mass row to 3.29e-02, three times higher than at `dtau0 = 1` --
and it is then removed by the iteration that follows, which is what a Newton
step is entitled to do and what a relaxation at `dtau0 = 1` cannot. This is
also why no tolerance on the linear residual could have repaired it: the
excursion is not an error in the step, it is the step.

`dtau0` is not changed here and its default is not this item's to move: the
1.0 that `Restart intent: stationary` fixes was itself the repair of a CFL
start that failed (READ, `backup/regression/atomic_elem_newton/README.md`:
at the CFL start the first hydrodynamic solve of that reload stopped on its
stagnation detector at `||R||` 1.490, where at 1.0 it hands back 1.167e-08).
What is measured here is that on the LHS 1140 b partitioned route 1.0 is far
below what these states can take, and that the ladder both memos diagnosed is
its consequence.

## 6. Controls

### 6.1 The delivered defaults are the entry text, bitwise

With both options at their delivered defaults, on the pass-2 reproduction
(three-unknown route) and on the carrier reload (four-unknown route),
measured build against control build: `output/Hydro_ioniz.txt` and `output/Ion_species.txt`
IDENTICAL byte for byte apart from the provenance line's timestamp, and every
printed iteration of the log the same to every digit. The control build is
bit-reproducible run to run (checked, both cases twice), so that comparison
means what it says. Getting there required one parenthesis and one rewrite of
P4 (section 2, item 6). The parenthesis:
`(idtau*D)/D` rounds twice and moved the state at 1e-11; `idtau*(D/D)` is
`idtau` exactly, and the reason is at the code site.

### 6.2 `backup/regression/atomic_elem_newton`, the element reload

HD 209458 b, atomic with solar trace metals and binary element diffusion,
reloaded from `IC/` with `Load IC? True`, `Coupled carrier solve: False`,
`Restart intent: stationary`; `EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80
EXHALE_DIFF_OMEGA=0.5`, 8 threads, the README's own recipe. The measured
build here is the judged row scaling turned ON, which is what this control
was for.

**CERTIFIED in both**, "every active equation of this state is within its own
tolerance": the control at outer pass 11, the measured build at outer pass 12.

| pass | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| worst gated species row, control | 2.91e-3 | 4.21e-3 | 1.77e-3 | 9.42e-4 | 4.66e-4 | 2.42e-4 | 1.23e-4 | 6.30e-5 | 3.22e-5 | 1.64e-5 | 8.42e-6 | -- |
| the same, judged rows | 2.61e-3 | 2.78e-4 | 3.96e-3 | 1.63e-3 | 8.53e-4 | 4.20e-4 | 2.11e-4 | 1.08e-4 | 5.25e-5 | 2.77e-5 | 1.30e-5 | 7.34e-6 |

The relaxation halves the row per pass in both. What the judged rows cost is
the pass-2 hydrodynamic solve: it ends "no descent direction exists for the
banded model at this state" at `||R||` 6.310e-01 after 14 iterations, against
the control's `info = 1` at 8.40e-10, the element relaxation then halves
`omega` to 0.250 for the rest of the solve where the control holds 0.500, and
the certification arrives one pass later. That one pass is the measured price
of the option and the reason it is delivered off (section 2).

### 6.3 `backup/regression/carrier_model_a_newton`, the carrier reload

`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40 OMP_NUM_THREADS=8 ./run.sh`.
This is the four-unknown route, which the change does not touch by
construction: **`output/Hydro_ioniz.txt` and `output/Ion_species.txt` are
BITWISE IDENTICAL** to the control build's, provenance line aside, with the
delivered defaults. Both end NOT CERTIFIED on the same two entries, carrier
balance H2 6.996e-04 at cell 363 and carrier balance H+ 7.425e-04 at cell 323.

That bitwise result is also what caught a defect in the first form of P4.
The count was taken by a second sweep of `residual_row_scale` and
`mass_row_rounding_floor` over the column at every outer iteration, and those
go through the caloric equation of state; the extra calls moved this fixture's
hydrodynamic rows at the 1e-11 level (the control build is bit-reproducible
run to run, checked, so the difference was the print). The counts now come
from `certified_row_measures`, which has the row measures and the rounding
floors in hand already, and nothing is evaluated twice.

### 6.4 The regression cases that reach this code path

With both options at their delivered defaults the binary is the entry text
BITWISE on both routes -- the pass-2 reproduction above (three-unknown) and
the carrier reload (four-unknown), `Hydro_ioniz.txt` and `Ion_species.txt`
identical byte for byte apart from the provenance timestamp. `wasp_full_newton`
and `mol_base_handoff` enter `solve_steady_jfnk` through those same two
branches, so no golden can move and none was refreshed.

The harness was nevertheless run on the CONTROL build, on scratch copies of
the two case directories, `REGRESSION_EXE` pointing at it and
`REGRESSION_GOLDEN_DIR` at the tree's `golden/`, to record where the tree
stands independently of this item:

| case | file | verdict |
|---|---|---|
| `wasp_full_newton` | `Hydro_ioniz.txt` | PASS within 1e-3, max rel 1.321e-08 |
| `wasp_full_newton` | `Ion_species.txt` | PASS within 1e-3, max rel 2.512e-08 |
| `wasp_full_newton` | `Hydro_ioniz_adv.txt` | **FAIL**, max rel 1.000e+00 at row 259 col 10, 0 against 4.001e-14 |
| `wasp_full_newton` | `Ion_species_adv.txt` | PASS within 1e-3, max rel 2.512e-08 |
| `mol_base_handoff` | all four files | PASS within 1e-3, max rel 8.703e-09 to 7.231e-08 |

That one FAIL is a zero compared with 4e-14 -- a relative measure of one by
construction on a quantity at the rounding floor -- and it is the CONTROL
build, so it is the state of the tree before this item and not something this
item did. It is reported here and not acted on.

### 6.5 The test suites

Run from the isolated tree against its private object directory, on the
measured build and on the control build:

| suite | control | measured |
|---|---|---|
| `krylov_and_dogleg` | all rows passed | all rows passed |
| `steady_species_rows` | FAIL 0 | FAIL 0 |
| `certification` | every assertion passed | every assertion passed |
| `attempted_step` | every assertion passed | every assertion passed |
| `carrier_retry` | every assertion passed | every assertion passed |
| `steady_selfconsistent_residual` | needs a run log (`EXHALE_STEADY_RUNLOGS`), not run | the same |

`steady_species_rows` and `steady_selfconsistent_residual` are the two suites
whose drivers link `steady_newton.o`; the other four were run because the
brief names them. The nine `inventory_reference_*` rows of
`steady_species_rows` fail identically on BOTH builds when the isolated tree
has no `backup/` (it is not in the git remote); with that directory in place
they pass on both.

## 7. The two cases the item ends on, and they certify

The recipe is the ENTRY TEXT binary and one environment variable: the raw
seed is given to the loop at the pseudo-time start the restart contract fixes
(`dtau0 = 1`), and the state that produces is given to the loop again at
`EXHALE_PTC_DTAU0=1.0e8`.

### `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`, the case of L4b, L4d and L4f

Resumed from the pass-8 entry state (the state its first seven passes
produce, which is L4f's reproduction), 8 threads:

```
outer pass 1: hydro info=0, worst gated species row 6.20E-06 of 1.0E-05
              at cell 217 (elemental transport He/H partition),
              mass 4.57E-09, momentum 6.91E-14, energy 7.06E-08,
              omega 0.500, 22.78 s
outer pass 1: ACCEPTED -- every active equation of this state is within
              its own tolerance.
```

`CERTIFIED: every active equation was evaluated and is within its tolerance`,
certified in the wind at `r >= 1.20`, written with `certified=T
cert_reason=certified_in_wind`, flux gate accepted 1.7637e-12. **This is a
certified LHS 1140 b solution with the C, N and O reservoirs, which no item
of this plan had reached.**

The like-for-like control is the same state, the same binary and `dtau0 = 1`:
its first pass takes 783.6 s and ends `hydro info = 2` with the elemental row
already inside (6.37e-06 of 1e-05) and the ENERGY row at 1.03e-03, which is
the number it entered at, and the loop goes on to a second pass doing the
same. 22.78 s and accepted against 783.6 s and refused, on one variable.

### `atomic_scalar_gj1132_kzz1e9/HeH4.0`, the L11 ladder step

Resumed from the state its REFUSED four-pass run at `dtau0 = 1` ends on:

| pass | hydro `info` | worst gated species row | mass | energy | `omega` | s |
|---|---|---|---|---|---|---|
| 1 | 0 | 3.13e-05 at cell 295 | 8.63e-10 | 2.09e-08 | 0.500 | 38.6 |
| 2 | 0 | 1.64e-05 at cell 297 | 9.01e-10 | 2.05e-08 | 0.500 | 32.6 |
| 3 | 0 | 8.61e-06 at cell 298 | 1.22e-09 | 1.92e-08 | 0.500 | 21.6 |

**ACCEPTED at outer pass 3**, 93 s in all, `omega` never halved. At
`dtau0 = 1` the same case halves `omega` twice and is REFUSED at pass 4 for
no joint progress (section 3 (c)).

### What does NOT work, measured

Giving the RAW seed to the loop at `dtau0 = 1e8` fails: the scalarCNO case's
first pass ends `hydro info = 2` with the mass row at 1.92 and the energy row
at 1.90, and the pseudo-time floor is then `dtau0` itself, so the cut
`dtau <- max(dtau/4, dtau0)` cannot recover. The 1.0 the restart contract
fixes is right where the seed is far from the root, and too small once it is
near it.

**The mechanism, and it is one line.** The pseudo-time is raised by
`dtau <- min(dtau*max(lam,0.1)*(f2/f2_try), 1e14*dtau0)` at every accepted
step, `f2` being the line-search merit `||F/D||_2`. On these states the merit
is nearly flat -- L4f measured a factor 8 over sixty iterations, and 58 of
those steps were non-monotone accepts -- so the ratio is about one, `dtau`
never leaves `dtau0`, and the solve spends its whole life as a relaxation at
`dtau = 1`. That it never leaves `dtau0` is directly visible: the damped
Gauss-Newton escape is entered only when `dtau <= dtau0`, and L4f counted 46
of them in 60 iterations. So the pseudo-time ramp is tied to a functional
that does not fall on exactly the states where the ramp is needed. That is
where the next item belongs, and it is NOT the acceptance test P2 proposed:
the same line search takes `lam = 1` at every iteration once `dtau` is large.

## 8. Reproduction

The measurements were made from an ISOLATED copy of the tree, `git archive
HEAD` (`43bc28c`) plus the working tree's `steady_newton.f90`, for the reason
L4f gives: other workers were editing `ion_residual_core.f90`, the `System_*`
family and five more files while this item ran, and a binary built from the
live tree would carry unfinished physics into every number above. The control
build is that tree with L4d's `steady_newton.f90`; the measured build is the
same tree with this item's.

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
T=$EX/scratchpad/L4e/tree
mkdir -p $T && cd $EX && git archive HEAD | tar -x -C $T
\cp -f $EX/src/modules/time_step/steady_newton.f90 $T/src/modules/time_step/
cd $T && PATH=/usr/bin:$PATH make OBJDIR=build_L4e_v2 EXE=EXHALE_L4e_v2.x -j12

# (a) the pass-8 entry state, 60 iterations, row-resolved
W=$EX/scratchpad/L4e/red_A          # input.inp, base.inp and output/*_IC.txt
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 \
         EXHALE_JFNK_MAXIT=60 EXHALE_LINEAR_ROWS=1 <binary>

# (b) the pass-2 entry state, twelve iterations, row-resolved
W=$EX/scratchpad/L4e/red_B
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 \
         EXHALE_JFNK_MAXIT=12 EXHALE_LINEAR_ROWS=1 <binary>

# (c) HeH4.0 from its own seed, four passes capped at 120 iterations
W=$EX/scratchpad/L4e/red_C
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=4 \
         EXHALE_JFNK_MAXIT=120 <binary>

#   EXHALE_JUDGED_ROWS=1      the judged row scaling (P1)
#   EXHALE_KRYLOV_TOL_ABS=1   and the absolute Krylov target beside it
#   neither                   the delivered default, bitwise the entry text

# THE MEASUREMENT OF SECTION 5 AND 7: the same states, one variable
cd $W && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 ... <binary>

# the fixtures, on copies of the case directories
#   atomic_elem_newton: Load IC? True, Coupled carrier solve: False,
#   Restart intent: stationary, IC/*_IC.txt into output/, then
#   EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80 EXHALE_DIFF_OMEGA=0.5
#   OMP_NUM_THREADS=8 ./run.sh <binary>
#   carrier_model_a_newton: EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40
#   OMP_NUM_THREADS=8 ./run.sh <binary>
# the suites, against the private object directory
#   EXHALE_OBJDIR=$T/build_L4e_v2 EXHALE_TEST_OUT=$T/build/tests/<suite> \
#     $T/src/tests/<suite>/run.sh
```

> Note added 2026-09-14 15:20: the private builds of this memo were made with `PATH=/usr/bin:$PATH make`, which selects the Ubuntu gfortran 13.1 and its LAPACK, not the conda-forge gfortran 16.2 + OpenBLAS the tree binary `EXHALE.x` is built with (bare `make`). The control and the measured build of the memo share one toolchain, so every comparison here stands; reproduce with bare `make` (`OBJDIR=... EXE=...` as written) to obtain the tree's toolchain.

## 9. Scope, and what was not measured

- **What was verified is the path the change touches**: the stationary solve
  of the three-unknown (partitioned) route on the three reproductions, the
  four-unknown route on the carrier reload (bitwise), the element reload with
  the option on, and the six test suites. The trust-region branch shares no
  line of the change: its Krylov tolerance is `tr_gm_rtol` and its row scaling
  `cell_row_scales`, both untouched.
- **No golden was refreshed and none needed to be**: with the delivered
  defaults the binary is the entry text bitwise on both routes.
- **The regression matrix was not run.** `wasp_full_newton` and
  `mol_base_handoff` were started on scratch copies with `REGRESSION_EXE`
  pointing at the control build and had not finished when this was written;
  they cannot move under the delivered defaults, for the reason section 6.4
  gives.
- **`EXHALE_PTC_DTAU0=1.0e8` is an experiment, not a proposal.** 1e8 is one
  decade grid point, chosen large enough to be a Newton iteration and not
  tuned; the measurement is that the states are two to five Newton iterations
  from a certified solve, not that 1e8 is the number to use. Nothing about
  `dtau0` was changed in the source.
- The two certified states of section 7 were reached from states the earlier
  passes produced, so they are certifications of THOSE cases, not of the
  runner recipe: a campaign case starts from its seed and would need the
  two-stage recipe written into `models/run_case.sh`, which is L4's file and
  not this item's.
- Wall times are contended: other workers held the machine throughout.

## 10. Noticed outside this item, reported and not acted on

- **The pseudo-time ramp is driven by the line-search merit**
  (`dtau <- dtau*max(lam,0.1)*(f2/f2_try)`), a functional that is flat on
  exactly the states where the ramp is needed (section 7). The judged
  distance is computed at every iteration already and is what the ledger
  ranks by; a ramp read from it, or from the fall of the worst judged row,
  is the natural repair and is one line. It is a change to the step control
  of every solve in the code and is not made without a decision.
- **The pseudo-time floor is `dtau0`**, so `dtau <- max(dtau/4, dtau0)`
  cannot cut below the start. That is why raising `dtau0` cannot be the
  recipe by itself: a start large enough for the near-root state leaves the
  cut with nowhere to go when a far state needs it (section 7).
- `steady_species_rows` reads `backup/regression/golden/`, which is not in
  the git remote, and reports nine FAIL rows in a tree that has no
  `backup/`. It is the same on any build; a missing reference reported as a
  failed assertion is a suite that cannot be run on a clean clone.
