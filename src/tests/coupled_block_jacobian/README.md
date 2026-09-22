# coupled_block_jacobian

T1 of `docs/lhs1140b_stationary_L22_step3_design_20260916.md` section 7: the
coupling of the coupled block, measured against a central difference of the
same map on a frozen state.

## What the block is

With `Coupled carrier solve` on, the registry adds one unknown and one row per
cell for every transported balance beside the three conserved hydrodynamic
ones, `nvar_jac = 3 + nspec_row`. The Krylov direction is built from
`jacobian_action_of_direction`, a directional finite difference of the FULL
residual, so no coupling derivative is written down anywhere: the cross terms
are in the action to first-order accuracy whether or not anyone has derived
them. What has to be measured is that the difference quotient is a usable
one, and that is this suite.

## How it is run

`run.sh` copies the frozen state
`LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55` (its `input.inp`,
its `base.inp` and its `output/`, with the spectrum path made absolute and
`Restart intent: stationary`) into the working directory four times and runs
the binary under `EXHALE_COUPLED_JAC_ACTION=<file>` with
`EXHALE_COUPLED_JAC_CELLS=300,312` and `EXHALE_COUPLED_JAC_DIR` set to
`carrier` and to `energy`, each direction at the probe arc the route sets
and again at `EXHALE_JV_PROBE_ARC=10` times that arc. The key registers the
block's rows for the measurement, writes the file and stops the run without
a step, so no invocation takes a solve and none writes a state.
`EXHALE_JV_PROBE_ARC` multiplies `probe_length_of_the_jacobian_action`, so
the second arc costs one invocation per direction (18 s MEASURED) and no
code of its own.

Every row but the power is read at the PRODUCTION arc. The probe scalar of
each file is written into it (`# probe scalar`, the scalar `eps1` the
direction is multiplied by; the displacement is `eps1*||v|| =
sqrt(epsilon)(1 + ||Y||)`, which is why the two directions print different
`eps1` for the same displacement, INSPECTED at
`steady_newton.f90:7529` and `:11613-11644`), and `run.sh` prints it beside
every statistic it reads, so that a drift can be attributed to the arc or to
the curvature.

Environment: `EXHALE_EXE` the binary, `EXHALE_TEST_OUT` the working
directory, `EXHALE_JAC_CASE` and `EXHALE_JAC_SED` the state and the spectrum.

## The rows

| row | what it states |
|---|---|
| `the_action_of_a_carrier_direction_moves_the_energy_row` | the composition reaches the energy row, through the pressure and the temperature at the conserved state (the molecular caloric equation of state moves the heat capacity with the composition) and through the heating and the cooling of the sweep the evaluation ran |
| `the_action_of_a_hydrodynamic_direction_moves_the_carrier_row` | the hydrodynamic state reaches the carrier row, through the advective face flux, the temperature of every rate, the density of the three-body terms, the diffusion and settling coefficients and the photolysis |
| `the_distance_to_the_central_difference_scales_as_the_first_power_of_the_arc_carrier` and `_energy` | the mean distance between the action the solve uses and a central difference of the same residual, over the rows the direction MOVES, moves with the FIRST power of the probe arc: measured at the production arc and at ten times it, `\|d ln e / d ln h\|` is 1 within 0.30 |
| `the_band_holds_the_two_cell_reconstruction_entry` | `d res_c(j+2)/d f_c(j)`, which the carrier relaxation's block-tridiagonal matrix has no place for, is inside the block's band: at a flat distance of `2*nvar_jac` = 8 within `kl_jac` = 11 |
| `the_species_rows_scale_is_the_certifications` | the carrier row of `eval_residual` is the row `carrier_steady_residual` returns converted by the code time scale `tscale_code = R0/(v0*n0)`, bitwise. Both sides are read in the units of the residual vector, the certification's row multiplied by `tscale_code` and never the residual divided back: for a double `x` and a scale `t`, `(x*t)/t` is not in general `x` again, so a bitwise statement made across that round trip fails on the rounding of the conversion alone |
| `the_outer_ghost_rule_is_the_same_in_both_evaluations` | the same comparison at cells N-1 and N |

## What is gated, and why it is a power and not a number

The statistic the file carries,

```
e(h) = mean over the moved rows of  |Jv(h) - cd1(h)| / max(|Jv(h)|, |cd1(h)|)
```

compares a FORWARD difference and a CENTRAL difference of one residual AT THE
SAME ARC `h`, so to leading order it is `(h/2)|F''|/|F'|` plus a rounding
term of order `eps|F|/h`: **the truncation of the forward difference, and not
a property of the assembly**. Its value at one arc therefore states that arc
and little else, and a gate on that value is a gate on the arc. MEASURED on
this fixture, these cells and this statistic, energy direction
(`docs/lhs1140b_coupled_action_20260921.md` section 10.2):

| arc multiplier | 1e-4 | 1e-3 | 1e-2 | 1e-1 | 1 | 10 |
|---|---|---|---|---|---|---|
| `e` | 1.2390e-06 | 1.2388e-05 | 1.2396e-04 | 1.2381e-03 | 1.2379e-02 | 1.1897e-01 |

an exact first power over four decades with no floor, while the central
difference it is measured against holds to 0.19 per cent across the same
range. The carrier direction is the other half of the same V: 2.9313e-01,
1.8570e-02, 2.3836e-03, 4.6314e-04, 5.0222e-05 at multipliers 1e-3 to 10,
falling as `1/h` on the rounding branch, its worst moved row at multiplier
1e-3 standing at exactly 2.0, the largest value this relative measure can
take, three decades below its production arc.

So the gate is on the POWER, which is free of the arc: `e` moves with the
first power of `h`, upward on the truncation branch and downward on the
rounding branch, and **a constant error of the assembly is exactly what
cannot do that**. `e` is taken at the production arc and at ten times it, and
`|log10(e(10h)/e(h))|` is gated at 1 within **0.30**.

The tolerance is stated from the ladder above, before the suite was run with
it: it admits every decade of that ladder, the widest departure being 0.289
at the carrier decade 1e-1 to 1 where the V turns, so the curvature it admits
is a factor 2.0 in the ratio over one decade; and it refuses a constant error
by 0.70. MEASURED, by adding a constant `C` to both statistics of a direction
in the reader: the carrier row refuses at `C = 1.0e-04` (power 0.574) and
passes at `5.0e-05` (0.709), the energy row refuses at `C = 3.0e-02` (0.546)
and passes at `1.4e-02` (0.703); the thresholds `(e(10h)+C)/(e(h)+C) =
10^0.70` are 5.3e-05 and 1.4e-02. That is the size of assembly error each
direction's row has teeth against, and it is bounded by how far that
direction's own `e` sits above its rounding floor, not by anything chosen
here.

**What the gate assumes, and when it would refuse for a physical reason.**
The two arcs must lie on one branch of the V. They do on this fixture: the
energy direction is on the truncation branch from 1e-4 to 10 and the carrier
one on the rounding branch from 1e-3 to 10, each with its turn more than a
decade away. A direction whose optimum fell between `h` and `10h` would read
a power near zero and refuse although its action were exact, which is why the
two probe scalars are printed and why the ladder above is kept here: the
first thing to read after a refusal is whether the V has moved under the
production arc.

`e` is read over the rows the direction MOVES, a row counting as moved when
its response reaches 1e-3 of the largest response in the support. A row whose
two sides both stand at the rounding floor carries no derivative to compare:
the mass row of a cell is a statement about the total density at a fixed
conserved state and does not move with the composition at all, so under a
carrier direction its quotient is a ratio of two round-off numbers and reads
1 whatever the action is. The mean over EVERY row of the support is printed
beside the gated one and is not the gate.

## The verdict at the time of writing (2026-09-17, MEASURED)

Carrier direction on cells 300-312: mean relative error over the rows the
direction moves 3.692e-04, worst 2.513e-03 over 40 rows; over every row of
the support 1.365e-01, worst 1.004 at the mass row of cell 306, whose two
sides are 1.1e-13 and 4.6e-16 against a momentum row of 2.4e-10 in the same
cell. Energy direction on the same cells: 3.792e-03 and 1.259e-02 over 21
moved rows; 7.150e-03 and 7.685e-02 over every row of the support. All
thirteen two-cell reconstruction entries nonzero, at about 1.0e-12. Thirteen
of thirteen species-row scale comparisons bitwise equal, and both outer ghost
rows bitwise equal.

The outer ghost row was expected to FAIL while two ghost rules were in force
(L22 step 2b section 1). It passes on the tree of 2026-09-17, whose
`diffusive_photochemistry.f90` is `f04ed90a1f1b34bdc619ede52793df6f`.

The numbers above are of that tree and were not re-measured on it, and **no
arc is recorded with them**, which is why the drift to the numbers below
cannot be attributed to the arc or to the curvature. Recording the arc, here
and in the printed rows, begins with the measurement below.

## The verdict on the tree of 2026-09-21 (MEASURED)

gfortran, `EXHALE.x` md5 `89ec67149aa452eea3deeb19052f83d1`, the same case
and cells.

| direction | probe scalar | `e`, moved rows | worst | moved rows | `e`, every row of the support | worst |
|---|---|---|---|---|---|---|
| carrier | 1.8041485e-06 | 4.6313904e-04 | 8.955e-03 | 30 | 2.8150786e-02 | 4.495e-01 |
| carrier | 1.8041485e-05 | 5.0222275e-05 | | 30 | | |
| energy | 6.8675961e-02 | 1.2378866e-02 | 4.103e-02 | 21 | 2.1526659e-02 | 1.200e-01 |
| energy | 6.8675961e-01 | 1.1897242e-01 | | 21 | | |

so the power is **0.9648** (carrier, on the rounding branch) and **0.9828**
(energy, on the truncation branch), 0.035 and 0.017 from 1 against a
tolerance of 0.30, and both rows PASS. The energy statistic at the
production arc, 1.238e-02, is what the earlier gate of 1.0e-2 refused; it is
printed and not gated, because at arc 1e-1 the same quantity reads 1.2381e-03
and at arc 10 it reads 1.190e-01 without anything in the tree having changed.
The question that gate was written to ask is answered by the power together
with the central difference itself, which is flat to 4.8e-03 over three
decades of the arc, from 1e-3 to 1, with the production arc inside that
interval (`docs/lhs1140b_coupled_action_20260921.md` section 10.2).

All thirteen two-cell reconstruction entries are nonzero, and the thirteen
species-row scale comparisons and the two outer ghost rows are bitwise
equal. Two of the thirteen scale comparisons read unequal until the
comparison was moved into the units of the residual vector: their two sides
differed in the last digit the file prints, which is the rounding of the
conversion and not a difference between the two operators.

## RED before the increment

The entry text has no `EXHALE_COUPLED_JAC_ACTION`: no file is written and no
row of this suite can be stated against it.

RED of the power row, MEASURED 2026-09-21 by doctoring the long-arc file in a
scratch copy and reading it with the same reader. With the long-arc statistic
set to the production one, which is what a constant error of the assembly
produces, the power reads 0.0000 and both directions FAIL. With the long
arc's forward difference replaced by the central difference of the same arc,
so that the statistic is exactly zero, the power is `unread` and both
directions FAIL rather than reaching the logarithm.

---

# `coupled_block_linear_system.f90`: the linear system of one frozen iterate

Added for item L32 of `docs/PLAN_20260917.md`. A program, not a key of the
production binary: it links the production objects and calls the same
`eval_residual` and `jacobian_action_of_direction` the solve calls, so that
several samples of ONE operator can be taken between one initialization and
one exit. That is necessary because the radiation and chemistry caches the
residual assembly keeps outside `steady_newton` have no accessor a second
process could be restored through, which is what `l15_dump_state` says in
its own comment.

Run it with `run_linear_system.sh`. `EXHALE_OBJDIR` names the build whose
objects are linked, `EXHALE_TEST_OBJDIR` where this program's objects go,
`EXHALE_L32_CASE` the run directory of the state (`input.inp`, `base.inp`,
`output/` with the state). With `EXHALE_L32_CASE` unset it builds and stops.

What it measures, all at the state the run directory holds:

- repeated residual and directional-action evaluations, and the same again
  after an unrelated state has been evaluated through the caches in
  between. The reference is the control of `src/tests/residual_determinism`,
  `1.6e-12` in row-scale units;
- a perturbation ladder for the action at four arcs, the maximum and the
  RMS of the error of the forward difference against a central difference
  over the same arc, reported separately for the mass, momentum, energy,
  elemental and carrier rows, with the cell of each maximum and the room
  the species box leaves along the direction. Three directions: the scaled
  right-hand side, the carrier unknowns, the energy unknowns;
- homogeneity, `A(3v)` against `3 A(v)`, and additivity,
  `A(v1+v2) - A(v1) - A(v2)`;
- with `EXHALE_L32_ASSEMBLE=1` and a small grid, the Jacobian assembled
  column by column from the same action, a direct solve by `dgesv`, and the
  two residuals of that step: against the assembled matrix, which says
  whether the direct solve solved what it was given, and against one fresh
  product of the matrix-free action, which says whether the assembled
  matrix is that action at all. `EXHALE_L32_IDTAU` sets the pseudo-time
  shift; the default 0 is the shift the trust-region Krylov leg uses
  (`steady_newton.f90` line 15328).

The rows it asserts are that the residual and the action are one map
(against the noise floor) and that the action is homogeneous; the ladder
and the additivity defect are printed to be read, not gated, because no
tolerance for them was stated before they were measured. The verdict of the
first measurement is `docs/lhs1140b_stationary_L32_20260917.md`.
