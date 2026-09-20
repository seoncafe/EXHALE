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
`Restart intent: stationary`) into the working directory twice and runs the
binary under `EXHALE_COUPLED_JAC_ACTION=<file>` with
`EXHALE_COUPLED_JAC_CELLS=300,312` and `EXHALE_COUPLED_JAC_DIR` set to
`carrier` and to `energy`. The key registers the block's rows for the
measurement, writes the file and stops the run without a step, so neither
invocation takes a solve and neither writes a state.

Environment: `EXHALE_EXE` the binary, `EXHALE_TEST_OUT` the working
directory, `EXHALE_JAC_CASE` and `EXHALE_JAC_SED` the state and the spectrum.

## The rows

| row | what it states |
|---|---|
| `the_action_of_a_carrier_direction_moves_the_energy_row` | the composition reaches the energy row, through the pressure and the temperature at the conserved state (the molecular caloric equation of state moves the heat capacity with the composition) and through the heating and the cooling of the sweep the evaluation ran |
| `the_action_of_a_hydrodynamic_direction_moves_the_carrier_row` | the hydrodynamic state reaches the carrier row, through the advective face flux, the temperature of every rate, the density of the three-body terms, the diffusion and settling coefficients and the photolysis |
| `the_assembled_action_reproduces_the_central_difference_carrier` and `_energy` | the action the solve uses against a central difference of the same residual, over the rows the direction MOVES, at or below 1.0e-2 |
| `the_band_holds_the_two_cell_reconstruction_entry` | `d res_c(j+2)/d f_c(j)`, which the carrier relaxation's block-tridiagonal matrix has no place for, is inside the block's band: at a flat distance of `2*nvar_jac` = 8 within `kl_jac` = 11 |
| `the_species_rows_scale_is_the_certifications` | the carrier row of `eval_residual`, divided by the code time scale `tscale_code = R0/(v0*n0)` it was converted by, is the row `carrier_steady_residual` returns, bitwise |
| `the_outer_ghost_rule_is_the_same_in_both_evaluations` | the same comparison at cells N-1 and N |

## The tolerance, and which rows it is read over

1.0e-2, stated before the measurement. It is the order of the 1.542e-03 that
L22 step 2b section 5 measured for the carrier operator's own assembled
action and accepted, with headroom, and the criterion of section 4.1 of the
design is read against it.

It is read over the rows the direction MOVES, a row counting as moved when
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

## RED before the increment

The entry text has no `EXHALE_COUPLED_JAC_ACTION`: no file is written and no
row of this suite can be stated against it.

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
