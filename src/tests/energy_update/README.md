# Tests of the implicit energy source step

What is under test: the residual-controlled temperature solve of
`src/modules/time_step/energy_semi_implicit.f90` (step B2 of
`docs/PLAN_20260906_rev2.md`, replacing item C3 of
`docs/d0_governing_system_20260906.md`).

The equation of one cell, at frozen composition and frozen heating, is

```text
R(T) = u(T) - u(T_old) - a [ heat - cool(T) ] = 0 ,   a = dt/(n_tot + n_e) ,
```

with `u` the internal energy per `(n_tot + n_e)` from the caloric equation of
state. The solve must return a temperature whose residual passes
`|R|/scale <= energy_res_tol` with
`scale = |E_old| + dt (|heat| + |cool|)`, must write back the cooling
evaluated at that same temperature, and must report a named failure rather
than a clamped or unconverged state.

Run it:

```bash
src/tests/energy_update/run.sh
# with a private build:
EXHALE_OBJDIR=$PWD/build_x EXHALE_TEST_OUT=/tmp/eu src/tests/energy_update/run.sh
```

The driver links the production objects and drives the solver through its
three iteration routines (`energy_balance_init`, `energy_balance_update`,
`energy_balance_finish`), supplying analytic cooling laws in place of
`eval_cool`. That is what gives each case a closed form or an independently
integrable root; the production wrapper differs only in where the cooling
comes from and in what it does with the answer.

A failed update stops nothing: the routine reports its status and assembles
no state, and the caller refuses the attempted step on it. The failure cases
below therefore run in the driver exactly as they run in production.

| assertion | what it establishes |
|---|---|
| `linear_cooling_status`, `linear_cooling_closed_form` | with `cool = c (T - T_eq)` and no heating the implicit step is linear in `T`; the solve reproduces its closed-form root |
| `linear_cooling_returned_cool_at_returned_T` | the cooling written back is `cool(T_returned)`, not the cooling of a previous iterate |
| `linear_cooling_residual` | the returned residual passes the acceptance tolerance |
| `falling_branch_root` | on a cooling that DECREASES with temperature, as metal line cooling does above its peak, the solve reaches the root of an independent bisection of the same equation |
| `falling_branch_cool_at_returned_T`, `falling_branch_residual` | the same two conditions on that branch |
| `falling_branch_two_iteration_residual_red` | the previous update (two iterations, cooling derivative frozen and entering as `1 + c abs(dC/dT)`, cooling refreshed once, no residual test) leaves a residual ABOVE the acceptance tolerance on the same problem. This is the red measurement the step exists for, kept in the suite as a standing reference |
| `frozen_heating_residual_with_frozen_heat` | the heating is frozen by construction; the balance the solve closes is the one with the frozen value, and it closes to tolerance. The diagnostic line beside it measures how far that state is from a balance whose heating followed `T`, which is the coupled source step of B3c, not this solve |
| `no_bracket_status`, `no_bracket_reason` | a heating so large that the balance is not reached below the ceiling of the physical range returns `ENERGY_UPDATE_NO_BRACKET` with reason `ENERGY_REASON_CEILING` |
| `floor_status`, `floor_reason`, `floor_last_iterate_is_the_floor` | a balance that asks for a temperature below the floor of the equation-of-state validity returns `ENERGY_UPDATE_FLOOR`, and no state is built from the floor |
| `floor_residual_energy`, `floor_missing_source` | the reported energy deficit, and the volumetric rate a specified external reservoir would have to supply for the floor to be a solution, match the analytic values. Specifying such a reservoir is B1's business, not this step's |
| `floor_counters_owned_by_the_wrapper` | the core does not touch `n_energy_floor_hits`; the production wrapper counts floor FAILURES there |

---

# Tests of the Crank-Nicolson transport stage

What is under test: `viscous_conduction_step` of
`src/modules/time_step/viscous_conduction.f90` (step B2b of
`docs/PLAN_20260906_rev2.md`), the second implicit temperature update of the
marching loop and operation 11 of the enumerated step of
`docs/b1_target_system_20260906.md` section 7.

With the conduction key on and the viscosity key off the equation of one
cell is

```text
C_j (T*_j - T_j)/dt = (1/2) [ Q_j(T) + Q_j(T*) ] ,
Q_j = (1/dV_j) [ F_{j+1/2} (T_{j+1} - T_j) - F_{j-1/2} (T_j - T_{j-1}) ] ,
```

with `F` the face area times the interpolated conductivity over the centre
spacing, the base ghost held fixed (a Dirichlet anchor) and zero diffusive
flux through the outer face. The stage must return a state only when every
solved temperature is finite and above the floor of the range in which the
equation of state and the transport coefficients are defined; a cell below
that floor is a failure with a named status, not a clamp.

Run it: the same `run.sh` builds and runs this driver after the energy-update
one. It sets `conduction_stop_on_failure = .false.` so that the failure
statuses can be exercised; production runs leave it `.true.` and stop.

The driver assembles the same equation from the flux form and solves it with
its own elimination, taking the conductivity from the module's public
`thermal_conductivity`, so the reference is of the discretization and the
solve rather than of the coefficient fit.

| assertion | what it establishes |
|---|---|
| `above_floor_status`, `above_floor_temperature_vs_reference` | a column whose solution stays above the floor returns `CONDUCTION_OK` and the temperature of the independent Crank-Nicolson solve |
| `above_floor_no_floor_hit`, `above_floor_min_temperature_over_floor` | nothing is counted, and every returned temperature is strictly above the floor, which is what makes this the admissible branch |
| `above_floor_energy_row` | the energy row returned is the one the pressure row it was rebuilt from implies |
| `floor_status`, `floor_reason`, `floor_cell` | a column drained against a cold base anchor, so that the operator asks the base cell for a temperature below the floor, returns `CONDUCTION_FLOOR` with reason `FLOOR_REACHED` and names the first failing cell |
| `floor_reported_solution` | the temperature reported is the unclamped value the operator solved for, matching the independent solve |
| `floor_energy_deficit`, `floor_energy_deficit_is_positive` | the reported energy is `c_v (T_floor - T_solved)` per volume in erg/cm3, the energy a floor would have injected with no source term behind it, and it matches the analytic value |
| `floor_state_unchanged_W`, `floor_state_unchanged_u` | no state is assembled from the floor: the caller still holds, bit for bit, what it passed in |
| `floor_counted_as_an_attempt`, `floor_failing_cell_count` | the floor counters keep working, now as counters of refused attempts, one per cell the operator asked to go below the floor |
| `floor_clamp_would_have_moved_the_state_red` | the standing red measurement: the clamp this step removed would have returned a temperature that differs from the solved one, which is the unbudgeted accepted correction the step exists to remove |
| `nonfinite_status`, `nonfinite_reason`, `nonfinite_cell`, `nonfinite_state_unchanged_*` | a non-finite temperature on entry returns `CONDUCTION_NONFINITE` with the cell, and no state |
| `closed_column_operator_conserves` | the conduction operator is in flux form, so on a column with no flux through either end the total thermal energy it moves is zero to round-off whatever the interior profile is |
| `stage_conserves_up_to_the_base_face_flux` | the stage as boundary-conditioned in a run is NOT closed at the base: the ghost is a Dirichlet anchor, and the total thermal energy change equals the heat crossing that one face at the Crank-Nicolson average of the two times, to round-off. That exchange with the anchored lower atmosphere is what the term is for |
