# D2b: the lower boundary as a contact, upwinded after the acoustic matching

Item D2b of `docs/PLAN_20260918_rev2.md`, on the user's decision of
2026-09-18 (`docs/PLAN_20260918_rev2.md`, "Decisions made", row D2a): choice
**A** at the pressure-balanced stationary contact, with the rest of the
Codex D2a text of `docs/DECISION_D2a_D8_codex.md` sections 3.1 to 3.5 kept.
Every number below is labelled MEASURED (run here) or READ (from a source
file or a memo).

---

## 1. Verdict, in three sentences

The lower boundary is now a **contact between the lower atmosphere and the
first interior cell, closed in three steps in that order**: the two sides
are matched acoustically by the interior's outgoing relation (C-) closed
with the reservoir's pressure; the direction of the contact is read off the
matched state; the entropy and composition of the face are then upwinded on
that direction, with no intermediate value and no transition width. **At
inflow and at a pressure-balanced contact at rest the reservoir owns the
thermodynamic state of the level** (choice A), the physical assumption being
that the lower atmosphere is a heat bath on the time scales of a stationary
solution, so a column at rest above it takes the level's temperature; only a
**reverse** flow carries the interior's own entropy and composition out, and
the reservoir then states the single condition its one entering
characteristic allows, the pressure. **Each transport operator carries its
own base condition and none of them follows the contact's direction**:
thermal conduction holds the level's own temperature in both directions,
element diffusion holds the reservoir's composition as a Dirichlet value and
now reports the diffusive flux that drives, and the molecular carrier
transport carries no diffusive flux across the base face.

The boundary model identity is `..._contact_upwind_v2` in place of
`..._smoothstep_v1`; the reservoir prescription stays version 1.

---

## 2. The measurement that shaped the implementation

The decision says the direction is decided "from the RESULT of that matching
(the sign of the face velocity), not from the interior cell velocity".
Taken as a bare sign test on the matched face velocity `v_b`, that rule puts
**every certified LHS 1140 b state on the reversal branch**, which is the
opposite of what choice A is for. MEASURED with the control build on the
three catalog states, through `base_closure_candidates_probe`:

| state | `v_b` [cm/s] | `M_i` | `M_wind` | face flux / window flux | `rho_rev/rho_res` |
|---|---:|---:|---:|---:|---:|
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | **-7.254** | -5.605e-06 | +2.992e-07 | +0.9999844 | 0.9269 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | **-2.098** | -3.071e-05 | +5.757e-07 | +0.9996434 | 0.2531 |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | **-3.415** | -2.891e-05 | +4.789e-07 | +0.9995324 | 0.4047 |

All three carry the wind's own outward mass flux through the base face to
better than five parts in 1e4, and all three have a negative matched face
velocity. The reason is arithmetic, not physical: `v_b = v_i + (p_res -
p_i)/(rho_i c_i)`, and a converged state leaves the interior's continued
face pressure about 1e-4 above the reservoir's, which at a base Mach number
of 3e-7 is a velocity **two orders of magnitude above the flow's own** and
of the opposite sign. The matched face velocity is a resolved statement of
the direction only where the flow is itself resolved; at this operating
point it is dominated by the boundary's own pressure residual.

So the direction is read from **the mass flux the level carries**, which is
what the contact's speed is, and from the quantity that resolves it:

- the **wind window** `r >= r_flux`, where the cell-centred product IS the
  flux, when it carries a flux at all (`|M_wind| >= 1e-12`) and carries ONE
  flux (`d_window <= 2e-3`). Both constants and their calibration are
  unchanged and are READ from the module: the certified stationary states
  sit at `d_window` 4.06e-05 to 1.82e-04 and the states that are not winds
  at 1.39e-02 and above, a factor 73 gap with the threshold in its middle.
  They are now thresholds on the window's STANDING rather than weights in a
  blend;
- the **matched face velocity** where the window has no standing: a cold
  start, a column at rest, a state that is not a wind.

The first interior cell's own velocity is not a candidate at all, which is
the part of the decision that is implemented literally: READ from item L21,
the cell-centred product reads -2.00 F_wind at cell 1 of the certified
atomic state while every face of the grid carries +1.00000 F_wind.

**Zero to the matching's own precision is zero, and the reservoir owns it.**
At a pressure-balanced contact `p_b - p_i` is a cancellation, so the sign of
`v_b` there is the sign of an arithmetic remainder. MEASURED on a discrete
hydrostatic column built on the reservoir's own isentrope with the molecular
equation of state: `v_b = -1.218e-16` of the code velocity unit, which a
bare sign test reads as a reverse flow and hands the level to the interior,
the one state choice A exists to answer. The implementation therefore
compares `v_b` against a floor of a few rounding units of the terms the
subtraction is made of,

    v_floor = 4 eps_mach [ (|p_i| + |p_b|)/(rho_i c_i) + |v_i| ] ,

which is about 1e-15 of the sound speed. It is not a tuned width: it scales
with the state, it vanishes with the arithmetic precision, and it cannot
move a branch on any state whose face Mach number exceeds 1e-14.

---

## 3. The contact law as implemented

`src/modules/states/base_boundary.f90`,
`characteristic_base_face_state`, in the order of the code:

1. the interior state at the face, continued from cell 1 along cell 1's own
   hydrostatic isentrope (unchanged);
2. the reservoir state at the face, continued from `r_base_level` along its
   own hydrostatic isentrope (unchanged);
3. `p_b = p_res` (the incoming-invariant weight stays 0, unchanged);
4. **`v_b = v_i + (p_b - p_i)/(rho_i c_i)`, moved ahead of the branch.** The
   face velocity is a result of the matching and does not depend on the face
   density, so upwinding the contact on it is not circular;
5. the direction: `reverse` if the window has standing and `M_wind < 0`,
   else if `v_b < -v_floor`; `w_rev` is 1 on a reverse flow and 0 otherwise;
6. `rho_b = (1 - w_rev) rho_res + w_rev rho_rev` with `w_rev` in {0,1}, so
   the face carries one of the two isentropes and never a mixture;
7. the supersonic-outflow transition at `M_i = -1` (unchanged) and the face
   Mach cap (unchanged).

The physical assumption is stated at the head of the module, under "WHO OWNS
THE LEVEL WHEN NOTHING FLOWS THROUGH IT", with its validity range: a
subsonic base at a level lying inside the radiatively controlled lower
atmosphere (1 microbar for the hot Uranus of Koskinen et al. 2022, the
photochemical column's matching level for LHS 1140 b). It says nothing about
a supersonic base and nothing about a level above the region where radiation
controls the temperature.

### The regularization decision

**The reversal branch has no width, and the smoothstep is no longer a
physical rule.** `characteristic_branch_weight` keeps exactly one use, the
supersonic-outflow transition at `M_i = -1`, where the count of outgoing
characteristics changes from two to three; that width is declared at
`base_face_mach_blend` as a numerical regularization of a change of the
characteristic count, with its default 1e-8 and the statement that no state
of this code has been measured inside it. A convergence test in the width
therefore has nothing to converge: the row
`contact_face_state_independent_of_the_width_1e6` /`_1e10` asserts that the
face state is **bit-identical** when the width is moved two decades in each
direction, which it is on the delivered build and is not on the entry text
(1.50e-06 and 1.80e-06 relative). `EXHALE_BASE_BRANCH_ON_CELL1` and the
module flag behind it are removed: the discriminant they restored is the one
the decision rules out.

---

## 4. Diff by file

| file | what changed |
|---|---|
| `src/modules/states/base_boundary.f90` | the contact law of section 3; the order of the three steps and choice A stated at the head of the module with the heat-bath assumption and its validity range; the branch is discrete, `w_rev` in {0,1}, with the matching's rounding floor at zero; the wind window's two constants restated as thresholds on the window's STANDING (values and calibration unchanged); `base_face_mach_blend` restated as the supersonic-outflow regularization; `base_branch_on_wind_flux` and `EXHALE_BASE_BRANCH_ON_CELL1` removed; `base_reservoir_temperature_at(r)` added, the level's own temperature for an operator that exchanges energy with the lower atmosphere; the model identity moved to `characteristic_face_ps_reservoir_C_minus_contact_upwind_v2`; the diagnostic print restated |
| `src/modules/time_step/viscous_conduction.f90` | `conduction_base_level_T` added: the base Dirichlet temperature of the conduction operator, the level's own temperature, read at the three places that used the advective ghost (the conductivity at the base ghost, the explicit source row, the implicit step's base row); `conduction_base_heat_flux`, the budget entry of that condition, set from the base face's own term of the row it differences |
| `src/modules/functions/binary_element_diffusion.f90` | the Dirichlet reservoir base condition stated at the operator as a condition that drives a flux whichever way the gas moves; the face flux loop is no longer conditional on the optional output, and `element_base_diffusive_flux` records face 0 as this boundary condition's budget entry |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | the carrier base condition stated in one block: zero diffusive flux at face 0, an advective term upwinded on the FACE MASS FLUX and not on the contact's direction |
| `src/tests/boundary_state/contact_law_at_the_base.f90` | new, 22 assertions, the Codex 3.5 list |
| `src/tests/boundary_state/run.sh`, `README.md` | the `contact_law` row, its object directory and its statement |
| `src/tests/grid_and_gates/base_branch_discriminant.f90` | rewritten for the new law on a hydrostatic fixture (11 assertions) |
| `src/tests/grid_and_gates/base_boundary_continuity_probe.f90` | header only: what its derivative rows now measure |
| `docs/input_schema.md` | `T_base` is the temperature OF THE LEVEL, and what the two boundary-model identities mean |
| `LHS1140b/MODELS.md` section 4 | "What the base level is": the contact law, the transport conditions, the model identity |
| `src/modules/states/Apply_BC.f90`, `boundary_state_trace.f90`, `certification.f90` | not edited: the contract D5b-2 wrote is unchanged by this item and no gate was added |

---

## 5. RED and GREEN

Both builds are private (`make OBJDIR=build_D2b EXE=EXHALE_D2b.x` and a
control build of the entry text), single threaded, and both are deleted at
the end of the item.

### `src/tests/boundary_state/run.sh contact_law`, 22 assertions

MEASURED. On a build of the entry text three rows do not compile at all: the
transport operators had no base condition of their own, so
`conduction_base_level_T`, `conduction_base_heat_flux` and
`base_reservoir_temperature_at` do not exist and the driver's section d is
removed to build it. Of the nineteen rows that remain, seven FAIL:

| row | entry text | delivered |
|---|---|---|
| `contact_stationary_face_is_the_reservoir` | 2.500e-01 | 0.0 |
| `contact_stationary_takes_no_intermediate_trace` | 5.000e-01 | 0.0 |
| `contact_reverse_flow_carries_the_interior_trace` | 5.000e-01 | 1.0 |
| `contact_inflow_carries_the_reservoir_trace` | 5.000e-01 | 0.0 |
| `contact_inflow_face_is_the_reservoir` | 3.319e-06 | 0.0 |
| `contact_face_state_independent_of_the_width_1e6` | 1.501e-06 | 0.0 |
| `contact_face_state_independent_of_the_width_1e10` | 1.802e-06 | 0.0 |
| the other twelve | PASS | PASS |
| the three transport rows | does not compile | PASS |

The reported, not asserted, quantity of section d: the conductive heat flux
through the base face with the gas at rest and the interior standing at
twice the level's temperature, MEASURED at **-5.6474708e-05** in code units,
positive inward, i.e. out of the domain. A zero bulk velocity does not make
it zero, which is why it is measured.

### `src/tests/grid_and_gates/run.sh base_branch`, 11 assertions

MEASURED, on a hydrostatic fixture built from the face along the
reservoir's own isentrope at twice the level's temperature. Entry text: 5
FAIL. Delivered: 11 PASS.

| row | entry text | delivered |
|---|---|---|
| `base_branch_local_face_wins_where_the_window_is_not_one_wind` | 0.0 | 1.0 |
| `base_branch_reverse_face_is_the_interior_trace` | 1.000012 | 0.0 |
| `base_branch_at_rest_the_reservoir_owns_the_level` | 5.000e-01 | 0.0 |
| `base_branch_at_rest_face_is_the_reservoir` | 2.500e-01 | 0.0 |
| `base_branch_sweep_takes_only_the_two_isentropes` | 1 distinct value | 2 distinct values |
| the other six | PASS | PASS |

The discriminating row is the second group: the first interior cell reads
`M_i = +5.475e-07`, an inflow, while the matched face velocity reads
`v_b = -9.940e-06`, a reverse flow, and the face density of the two branches
differs by a factor 2.

### The suites that link the modules this item changed

MEASURED with the delivered build, single threaded:

| suite or row | result |
|---|---|
| `src/tests/boundary_state/run.sh` (all three rows) | 35 PASS, 0 FAIL, "every test program passed". The D5b-1 and D5b-2 rows run on `carrier_model_a_newton`, which has `Conduction: True`, so they also cover the conduction boundary condition this item states |
| `src/tests/grid_and_gates/run.sh base_level base_cell_width momentum_row output_state grid_width grid_window` | 39 PASS, 0 FAIL, "every test program passed" |
| `src/tests/grid_and_gates` `hydrostatic_residual` (the four-resolution ladder and the verdict pass) | 60 PASS, 0 FAIL, including `well_balanced_momentum_row_at_rounding` at 5.45e-14 and the mass and energy rows at 1e-16 for all four reconstruction and flux pairs |
| `src/tests/grid_and_gates` `base_continuity` on the certified LHS 1140 b atomic state | 35 PASS, 0 FAIL |

### `src/tests/grid_and_gates/run.sh base_continuity`

MEASURED on the certified LHS 1140 b atomic state: 35 PASS, 0 FAIL on the
delivered build. Every derivative row now reads exactly 0.0, which is the
stronger form of what it was written to measure: a state whose base sits
decades inside one branch cannot be moved by a perturbation of the window.

---

## 6. Movement, MEASURED single threaded

No golden was refreshed and no regression or catalog directory was written
to; every run is on a scratch copy.

### The four regression cases

Data rows compared number by number, not by `diff`, so that a formatting
difference could not be read as a movement.

| case | steps | `Hydro_ioniz.txt` | `Ion_species.txt` | steady-state log10 Mdot |
|---|---|---|---|---|
| `wasp_full` | 400 | **identical, every number** | **identical, every number** | not quoted (marching snapshot) |
| `mol_base_handoff` | 12000 (pinned) | largest relative change 1.185e-05, the velocity at cell 203 (856.891 against 856.901 cm/s) | 2.403e-06, the H2 column at cell 27 | 10.57 both |
| `mol_carrier` | 12000 (pinned) | 5.226e-06, the velocity at cell 196 (699.141 against 699.145 cm/s) | 5.595e-06, the H2 column at cell 27 | 10.57 both |
| `lower_profile` | 12000 (pinned) | 3.312e-04, the cooling rate at cell 495 near the outer edge; density 1.41e-04, velocity 8.89e-05, pressure 2.99e-04, temperature 1.58e-04, all at cell 495 | 2.568e-04, column 11 at cell 494 | 8.82 both |

Per column on the hot-Uranus gate (`mol_base_handoff`, largest relative
change over the 504 rows): radius 0, density 4.44e-07, velocity 1.19e-05,
pressure 4.91e-07, temperature 4.58e-07, heating 4.29e-07, cooling 1.22e-06.
`mol_carrier`: 0, 8.36e-07, 5.23e-06, 1.40e-06, 8.25e-07, 1.04e-06,
2.87e-06. Every column is below 1e-4, and the two mass-loss rates agree to
the digits the run prints.

**Which branch moved, MEASURED on the written state of `mol_base_handoff`**
with the two builds through `base_closure_candidates_probe`: the interior
cell reads `M_i = +1.2324e-04`, an inflow, while the matched face velocity
reads `v_b = -230.05 cm/s`, a reverse flow, and the wind window is silent
there (`d_window` = 1.83e-01, READ from the module's own table). The entry
text follows the cell and takes the reservoir (`w_rev` = 0, face density
1.4204e+13 mH/cm3); the delivered build follows the face and takes the
interior's trace (`w_rev` = 1, 1.4179e+13), a change of 1.74e-03 in the face
density. That is the decided change acting: the direction is read from the
result of the matching and no longer from the first interior cell. The
movement it leaves after 12000 steps is 1.2e-05 because the two isentropes
are only 0.17 per cent apart on this state.

The same mechanism, MEASURED on the written state of `lower_profile`: the
interior cell reads `M_i = +8.367e-05`, the matched face velocity `v_b =
-808.97 cm/s`, the window is silent (`d_window` = 7.87e-01, READ), and the
branch goes from the reservoir (`w_rev` = 0, face density 6.1536e+12) to the
interior's trace (`w_rev` = 1, 6.1176e+12), a change of 5.85e-03. Its
largest movement after 12000 steps is 3.3e-04, at the outer edge of the
column, below the 1e-3 relative tolerance of the regression and below the
0.1 per cent the standing rule of 2026-09-05 counts as identical; it is
reported as a physics change nonetheless, because it is one.

`wasp_full` does not move because its base and its window select the same
branch under both models over its 400 steps.

Apart from the data rows, the written products differ in the header: the
timestamp and the `# boundary_model` line, which now reads
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2`.

### The one case with thermal conduction on

`backup/regression/carrier_model_a_newton` has `Conduction: True` and a
molecular base, so it is the only fixture on which the conduction boundary
condition of section 4 can move anything. MEASURED on a scratch copy, the
marching stage capped at 30 steps, control against delivered:

| product | differing numbers | largest relative change |
|---|---:|---|
| `Hydro_ioniz.txt` | 3021 of 3528 | 2.596e-05, the velocity of cell 3 (4.73057 against 4.73044) |
| `Ion_species.txt` | 18090 | 5.212e-06, column 38 of cell 3 |

That is the conduction operator reading the level's own temperature where it
used to read the advective ghost's `p/n_part`.

**CORRECTED 2026-09-20, item P6c
([lhs1140b_p6b_p6c_20260920.md](lhs1140b_p6b_p6c_20260920.md)).** This
paragraph and section 7 below attributed that distance to "the 4.7 per cent
between the prescribed and the solved ghost particle count". That is not the
distance this operator changed by, and 4.7 per cent is not a property of this
fixture. MEASURED on this state at the radius where both temperatures are
defined, the ghost's `p/n_part` and the level's own temperature stand at
9.1e-07, which is the size the movement above reports. The 4.66 per cent of
D5a section 1 is the ratio between a count prescribed at the base LEVEL and a
count solved in the lower GHOST CELL one cell below it, on an LHS 1140 b
grid whose base cells are 7.8e-02 pressure scale heights; on the hot-Uranus
grid of this fixture the same ratio is 5.5e-03. The evaluate route of the same fixture is byte-identical in
both products, as it must be: it takes no step. No golden was refreshed;
this case is not in the default regression matrix.

### The flowing catalog states of D5a, evaluate route

MEASURED with `Restart intent: stationary evaluate` and
`EXHALE_MASS_FLOOR_SCAN=1` on scratch copies:

| state | data rows | cell-1 continuity row | verdict |
|---|---|---|---|
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | identical, every number | 1.729e-11 against 3.8e-11 at cell 162, the largest measure at cell 1 within its 7.9e-09 | CERTIFIED, both builds |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | identical, every number | 1.558e-08 against 7.8e-09 at cell 1 | NOT CERTIFIED, both builds, same entry and same number |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | identical, every number | 8.240e-08 against 7.5e-09 at cell 1 | NOT CERTIFIED, both builds, same entry and same number |

The two molecular refusals are D5b-2's consequence and not this item's: READ
from `docs/lhs1140b_stationary_D5b2_20260918.md`, those states were already
refused on the cell-1 continuity row when the ghost composition became
solved, and their certificates are stale under boundary model v1. This item
does not move them by a single digit.

### L35's rest fixture, and the states beside it

MEASURED with `base_closure_candidates_probe` on the two builds. The
perturbation sets every cell's momentum to zero at fixed mass density and
fixed thermal energy, so the wind window carries no flux and the matched
face velocity decides.

| state | perturbation | entry text face rho | delivered face rho | `w_rev` before | after |
|---|---|---:|---:|---:|---:|
| `atomic_scalar_.../HeH2.13` | none | 9.505840969499e+13 | identical | 0 | 0 |
| `atomic_scalar_.../HeH2.13` | rest | 9.158406986462e+13 | 8.810973003425e+13 | 1/2 | 1 |
| `molecular_scalar_.../HeH2.13` | none | 1.126231901978e+14 | identical | 0 | 0 |
| `molecular_scalar_.../HeH2.13` | rest | 7.056540303112e+13 | 1.126231901978e+14 | 1/2 | 0 |
| `molecular_scalar_.../HeH9.7` | none | 1.212003870518e+14 | identical | 0 | 0 |
| `molecular_scalar_.../HeH9.7` | rest | 8.512719328908e+13 | 1.212003870518e+14 | 1/2 | 0 |
| `carrier_model_a_newton/IC` | none | 1.420370160866e+13 | identical | 0 | 0 |
| `carrier_model_a_newton/IC` | rest | 1.424295106443e+13 | 1.420370160866e+13 | 1/2 | 0 |

On every state that flows the face state is bit-identical, including on
`carrier_model_a_newton`, whose window has no standing (`d_window` =
1.51e-02) and where the matched face velocity therefore decides: it reads
+36.68 cm/s, an inflow, the same branch the interior cell's own velocity
selected. The `rest` states are where the model changed, and the change is
what choice A asks for: the midpoint of the two isentropes is gone. Three of
the four take the reservoir; the fourth, the atomic fiducial at rest, takes
the interior's trace, because deleting the momentum of a wind does not make
a pressure-balanced contact and the matched face velocity there is -6.667
cm/s, a resolved reverse flow 1e16 times its own rounding floor.

**L35 measured the two closures bit-identical on flowing inflow states, and
that still holds**, for the reason L35 gave: on a state whose base carries
one coherent wind the window and the face flux select the same branch, and
under this model both are "the reservoir owns the level", so the face state
is the reservoir's either way. What changed is only the states where neither
resolves an outward flow, and there the answer is now the reservoir instead
of the average.

---

### The element-diffusion budget entry

`element_base_diffusive_flux` now records the diffusive helium mass flux
through the base face at every step, driven by the Dirichlet reservoir
composition whatever the bulk velocity there is. It is a module quantity and
is not yet printed: `element_flux_profile.txt` is written over the faces
1 to N-1 and its window statistics are a gated quantity of the
`lower_profile` case, so putting face 0 into that file is a change of a
gated product and is left as a follow-up rather than made inside this item.
The point the entry exists to make is nevertheless MEASURED in that case's
own product: at the first face above the base, r = 1.0002932 R_p, the
diffusive flux reads -5.5736040e-13 g cm^-2 s^-1 against an advective
-4.1860913e-08, i.e. it is not zero where the gas is nearly at rest.

---

## 6b. The direction rule accepted, and the check added (2026-09-19)

The user accepted the direction rule of section 2 on 2026-09-19 with one
addition, made by the advisor: where the wind window decides the direction,
the certification report compares it with the base face mass flux the
Riemann solve assembled for the same state, and states the outcome in the
face-flux block, one of four:

- read from the matched face velocity (the window has no standing);
- read from the wind window, the base face flux at its rounding (at or below
  1e-12 of the window mean) and carrying no direction;
- read from the wind window, which agrees with the base face flux;
- DISAGREEMENT, the window and the base face flux carrying opposite signs,
  which is the case the decision text guarded against and the rule does not
  cover.

The judgment is `base_contact_direction_agreement` in
`src/modules/time_step/certification.f90`; it gates nothing. Six rows of
`src/tests/certification/` give each outcome on chosen numbers (109/0, the
103 before plus these six). MEASURED on the three catalog states of section
2 through the evaluate route: all three read "from the wind window, which
agrees with the base face flux (w_rev = 0.0)", and the data rows of the
atomic state are identical before and after the check (it prints and
changes nothing).

## 7. Noticed outside this item, reported and not fixed

**`write_setup_report.f90` labels the width with a statement that is no
longer true.** Line 597 prints `entropy-branch blend window, face Mach =`
followed by `base_face_mach_blend`. That constant no longer sets an
entropy-branch window; it regularizes the supersonic-outflow transition at
`M_i = -1`. The file is not this item's, so the label stands; the fix is one
string.

**The loader already states what a model change means.**
`load_IC.f90`'s `report_boundary_model_of_the_restart` prints, for a state
written under another boundary model, that the state is loaded and that "its
residual is not the residual it was written with", which is the statement
the item asked for; it does not use the word certificate and it does not
refuse the seed, which is correct. `load_IC.f90` is not this item's file, so
nothing was changed there.

**The ghost's temperature and the level's temperature are two numbers.** The
conduction operator now reads the level's own temperature, while the ghost
cell the advection writes carries `p/n_part` at the ghost's SOLVED particle
count. **CORRECTED 2026-09-20, item P6c:** the two are 9.1e-07 apart on this
state at the ghost cell that sits at the level's own radius (MEASURED,
[lhs1140b_p6b_p6c_20260920.md](lhs1140b_p6b_p6c_20260920.md) section 6.1),
and not 4.7 per cent; the per cent belongs to two counts read one cell apart
on a different grid. The distance is the model's own, and this item does not
reconcile it; it only makes clear which of the two each operator reads.

---

## 8. How to repeat it

The two builds, the run directories and the scripts are under
`<scratchpad>/D2b/` (`states/`, `runs/`, `reg/`, `probe/`, `bs/`, `gg/`,
`mkstate.sh`, `evalrun.sh`, `regcase.sh`, `cmp.py`).

- the suites: `EXHALE_OBJDIR=<objdir> EXHALE_EXE=<binary>
  EXHALE_TEST_OUT=<scratch> src/tests/boundary_state/run.sh contact_law` and
  the same for `src/tests/grid_and_gates/run.sh base_branch`;
- one catalog reading: copy `input.inp`, `base.inp` and the `_IC` pair into a
  scratch directory, set `Restart intent: stationary evaluate`, run with
  `EXHALE_MASS_FLOOR_SCAN=1`;
- the face state of a state and of a perturbation of it: build
  `src/tests/grid_and_gates/base_closure_candidates_probe.f90` against the
  object directory and run it in the state's directory with the state name
  and one of `none`, `rest`, `weak`, `reversed`, `window`.
