# Review and remediation plan for the EXHALE open defects

Date: 2026-09-03

## 1. Scope and evidence

This review compares `docs/open_defects_20260903.md` with the implementation in
the working tree. It follows the routines that select the molecular closure,
advance the carrier equations, construct the hydrodynamic residual, apply the
lower boundary, solve the steady system, serialize restart state, and evaluate
the affected diagnostics.

No production calculation was rerun for this review. Numerical values quoted by
the defect catalog remain historical measurements from the documents named
there. Conclusions marked **verified in code** are based on the current call
paths and expressions. Conclusions marked **reported measurement** have not been
remeasured here.

The working tree contains extensive uncommitted development, including new
physics modules. Therefore, this review describes the tree as inspected, not the
last commit (`6d07d48`). Source line numbers should be treated as identifiers for
this snapshot and will move as the code changes.

## 2. Executive judgment

The first release-quality objective should not be to tune the present molecular
answer. It should be to define and enforce a single discrete steady problem that
contains carrier transport and whose accepted root is also a fixed point of the
production update. The current default molecular calculation omits an H2
transport term that is reported to exceed the local net chemical term by factors
of 44--117 through part of the H2 front. Physical correctness therefore requires
transport; the inability of the coupled steady solver to converge does not make
the local-equilibrium closure acceptable.

The recommended order is:

1. Add invariant and residual diagnostics that can identify the exact failing
   equation and cell without changing the solution.
2. Repair the oxygen-carrier conservation defect and create a real restart
   round-trip test. These are prerequisites for trusting coupled-solver results.
3. Make the steady residual composition- and boundary-consistent, then compare
   it directly with one complete production update.
4. Replace the lower boundary with a face-based characteristic condition and
   verify the base energy and mass budgets.
5. Finish the coupled carrier/wind solve with a genuine trust-region method.
6. Only after a coupled solution passes the physical gates, make carrier
   transport the molecular default and refresh affected reference results.
7. Treat Lyman--Werner geometry and H2 photoionization branching as independent
   physics improvements after the steady-state foundation is reliable.

Items 1--6 of the catalog are not independent. Items 3, 5, and 6 are different
observations of the base-layer state, and item 2 blocks the physically correct
resolution of item 1. Solving them as unrelated patches would make it difficult
to establish which equation the final state satisfies.

## 3. Corrections to the defect catalog

Several catalog statements should be corrected before it is used as an
implementation specification.

### 3.1 `Heating_breakdown.txt` is already repaired

Item 6 says that the diagnostic under-reports heating by 49 percent and that its
patch is not in the repository. This is stale. In the current
`write_heat_breakdown_eq` implementation, `heat_tot` includes photoheating,
photoelectron heating, Lyman-alpha de-excitation, helium-recombination coupling,
both Penning channels, Lyman--Werner heating, FUV photolysis, and molecular
reaction heat. `docs/p23_thermal_budget.md` section 9.1 also marks the issue
fixed and reports agreement with the solver heating column to `1.6e-6` on the
state used there. This diagnostic issue must be removed from the open part of
item 6. The base energy-closure issue itself remains open.

The catalog also labels this diagnostic as P52, but `TO_BE_DONE.md` and
`docs/Update_EXHALE_stage1.md` use P52 for convergence-stop reporting and the unified
steady gate. That P52 is closed. A unique defect identifier is needed if the
heating diagnostic is discussed historically.

### 3.2 The lower boundary is not over-specified merely because it prescribes
two thermodynamic quantities

For outward subsonic flow through the inner boundary, two characteristics enter
the computational domain and one leaves it. Two reservoir conditions and one
interior compatibility relation are therefore expected. In the default branch,
`BC_component_constrho` prescribes density and pressure, while velocity is
obtained from the first interior cell through the one-way valve. The count of
conditions is not, by itself, excessive.

The actual implementation concerns are stronger and more precise:

- The boundary is imposed component by component rather than through the
  characteristic invariants.
- `Rec_BC` calls `BC_component_constrho` for reconstructed face states, but the
  velocity formula evaluates `r(index)` at a ghost-cell center. `RK_rhs` then
  treats the result as a state at `r_edg`. This verifies the location mismatch
  described by the catalog.
- The `Base velocity: massflux` branch prescribes density, velocity, and
  pressure. For subsonic inflow this branch does prescribe all three primitive
  variables and should be checked separately for characteristic consistency.
- Prescribing both density and pressure is a legitimate reservoir model, but it
  need not match the hydrostatic or thermal state selected by the resolved
  atmosphere. The measured pressure agreement and temperature disagreement
  diagnose incompatibility with that model, not an incorrect number of boundary
  conditions.

The refinement table in item 4 labels a large negative cell-centered quantity
as “base-face flux,” while the same section says the Riemann face flux is
positive and of order the wind flux. These are different quantities. The table
heading should say `cell-centered rho*v*r^2 in cell 1` unless the source data
show otherwise.

### 3.3 The steady residual is not the complete production-update map

`assemble_residual` and the marching loop share `Reconstruct` and `RK_rhs` for
the Euler part. They do not, by that fact alone, define identical fixed points.
The production step also includes an operator-split carrier advance, an
ionization sweep, a semi-implicit energy update, repeated boundary fills,
optional viscous/conductive evolution, and an optional Shapiro filter. The
steady residual eliminates some of this physics locally and represents other
parts as source terms. Item 3's conclusion that the discrepancy must be a small
nonzero residual in the same operator is therefore too strong.

The correct discriminator is the defect of the full update map,

```text
G_dt(q) = [Phi_dt(q) - q] / dt,
```

evaluated at the accepted steady state for decreasing `dt`. Agreement between
`G_dt` and `assemble_residual` by row and cell would establish that the two
routes solve the same discrete equations. Disagreement would locate a splitting,
boundary-ordering, or state-refresh defect.

### 3.4 The direction of the Lyman--Werner geometry error is not established

The code clearly applies a trapping factor tabulated from slab calculations to
a spherical atmosphere. That is a valid open modeling issue. However, the
catalog's claim that the code necessarily traps too much is not established by
the implementation or by a two-geometry calculation. Angular path lengths,
spherical dilution, inward optical depth, and frequency redistribution all enter
the escape probability. Until slab and spherical transfer are compared, the
sign and magnitude should be reported as unknown.

### 3.5 Restart transport is a fixed implementation with an open compatibility
policy, not an open numerical defect

The current tree writes and reads the coupling header, equilibrates a loaded
composition, and recomputes temperature after the post-solve composition sweep.
These call paths verify the three repairs in item 11. Old files without a header
remain intentionally ambiguous. Move item 11 to a migration/compatibility
section and retain only the stale `heat` and `cool` output columns as a small
diagnostic-consistency issue.

## 4. Additional implementation concern: the residual gate can hide a narrow
cellwise imbalance

This concern is more concrete than the hypothesis in item 8.

With `resid_vol = .true.` (the default), `relnorm_over_cells` computes

```text
sum_j |R_j| V_j / sum_j S_j V_j,
```

where `S_j` is the row scale. This is a scale-weighted average of the cellwise
ratios, not a requirement that every cell satisfy `|R_j|/S_j < tol`. An error in
one or a few cells can pass when the rest of a large region is quiet. Splitting
the grid at `j_min` prevents the outer wind from diluting the whole lower region,
but it does not prevent the hundreds of cells within the lower region from
diluting cell 1 or a narrow front.

The alternative branch is not a cellwise infinity norm either. It computes
`max_j |R_j| / max_j S_j`; the numerator and denominator can come from different
cells. The usual local relative infinity norm is
`max_j (|R_j|/S_j)`. The comments at `residual_norms` describe the latter more
closely than the code implements it.

This may explain how a state with a reported cell-1 energy imbalance near unity
can have a small aggregate norm, but that causal link has not been measured.
Before changing acceptance, write the following for each row and for the carrier
equations:

- `max_j |R_j|/S_j`, its cell index, and radius;
- the current ratio of weighted sums;
- the numerator and denominator contributions accumulated over `r < 1.03`,
  `1.03 <= r < 1.10`, `1.10 <= r < 1.20`, and `r >= 1.20`;
- the signed dimensional residual and every term in the worst cell.

If the accepted hot-Uranus state has the reported local imbalances, acceptance
should require both an integrated norm and a cellwise norm. The integrated norm
measures a column budget; the cellwise norm prevents a narrow physical layer
from being averaged away. Their tolerances should be calibrated separately.

## 5. Review of each catalog item

| Catalog item | Code judgment | Recommended disposition |
|---|---|---|
| 1. H2 local equilibrium | **Confirmed.** `carrier_transport` defaults to false, except that an unstated value is assigned from `thereis_oxychem`. `ioniz_eq` holds H2 fixed only when transport owns it; otherwise the cell system solves a local H2 balance. | Critical physics defect for molecular results. Decouple the option from oxygen, add a closure-validity diagnostic immediately, and make transport mandatory after the coupled solve is qualified. |
| 2. Coupled carrier/wind solve | **Confirmed implementation and open convergence failure; historical outcome not rerun.** The fourth unknown, carrier residual, JFNK path, damping, line search, and stagnation return are present. | Add a scaled trust region and verify the Jacobian/model before changing defaults. Do not loosen the carrier gate to obtain a nominal success. |
| 3. Accepted root moves under marching | **Reported measurement; explanation incomplete.** The Euler spatial routines are shared, but the full operators are not identical. | Measure the complete update-map defect. Repair ordering or splitting discrepancies before attributing the motion to a physical acoustic mode. |
| 4. Base velocity and ghost state | **Location mismatch confirmed; interpretation partly rejected.** The ghost-center radius is reused in a face state. Two reservoir conditions are appropriate for subsonic inflow. | Replace component-wise face overwrites with a characteristic face solve. Treat the cell-centered velocity and Riemann face flux as distinct diagnostics. |
| 5. Grid nonconvergence | **Reported measurement.** It cannot yet be assigned uniquely to item 3. Boundary truncation, operator switching, and sampling of a transient are all possible. | Repeat only after items 3 and 4 are repaired. Compare the same discrete steady definition on each grid, not states stopped at different oscillation phases. |
| 6. Base energy closure | **Physical requirement confirmed; cause open.** `assemble_residual` contains the expected energy row, but acceptance can average a narrow imbalance. The cited heating-output defect is already fixed. | Add local and shell-integrated energy gates. Compare solver terms with the full update map and then repair the term or boundary causing the mismatch. |
| 7. Marching `du` and steady flux gate | **Different definitions confirmed, but not necessarily one defect.** They serve termination, handoff, and acceptance roles. | Do not force one threshold onto all roles. Rename and document them by role; require the steady gates for a publishable state. Recalibrate handoff separately. |
| 8. Residual weighting | **Open, with a verified formulation concern.** Both current norm branches differ from a cellwise maximum of relative residuals. | Add dual integrated/local gates and contribution diagnostics before changing weights or split radii. |
| 9. Lyman--Werner geometry | **Model mismatch confirmed; sign and size unverified.** | Add an explicit geometry provenance flag and uncertainty. Implement a two-geometry transfer comparison before replacing the table. |
| 10. Boundary fill before ionization | **Ordering confirmed; effect unmeasured.** `eval_residual` fills conservative ghosts before the ionization sweep and then assembles the residual after the composition changes. | Make residual evaluation a deterministic composition fixed point and fill the boundary from that same state. Test repeated `F(Y)` evaluations and Jacobian directional derivatives. |
| 11. Restart state | **Main repairs confirmed in code.** Legacy files remain ambiguous and `heat/cool` can lag the corrected output state. | Add a restart schema version and an executable round-trip test; either recompute all diagnostic columns from the written state or label them stale. |
| 12.1--12.3 H2 photochannels | **Approximations documented in code.** Channel stoichiometry is incomplete for double ionization and neutral dissociation. | Introduce explicit cross sections for mutually exclusive final states and derive every species/electron source from their stoichiometric vectors. |
| 12.4 round-trip case | **Confirmed.** The hydro IC is empty, the species IC is absent, and the case is outside the default matrix. | Rebuild the test fixture from a deliberately short deterministic run, not from an expensive converged product. Test write/read/write identity and one-step continuity. |
| 12.5 stale arm case | **Confirmed configuration conflict; intent cannot be inferred from code.** | Quarantine it from active regression and label the historical result invalid. Choose atomic or molecular semantics only from the original campaign intent. |
| 12.6 base-H2 fit coverage | **Confirmed by the catalog's path analysis; not rerun.** | Add a small setup/equilibrium test that omits the handoff and asserts the fitted base H2 fraction and resulting ghost composition. |
| 12.7 PLM-to-WENO3 transition | **Open path verified; diagnosis historical.** The switch changes the discrete operator in one step. | Treat PLM and WENO3 as separate discretizations. Use PLM only to generate an initial state, then apply continuation in the residual or flux blend before solving the WENO3 system. |
| 12.8 oxygen-carrier double count | **Critical reported conservation failure; current write-back code contains the relevant repartition.** Exact offending refresh not isolated here. | Add element-budget assertions around every ionization sweep, carrier write-back, and outer steady pass. Fix this before any coupled A2 result is accepted. |
| 12.9 isothermal base state | **Reported measurement and closely coupled to item 4.** | Resolve through the characteristic boundary work; do not add another empirical ghost-temperature mode first. |
| 12.10 cold molecular layer | **Symptom, not an independent defect yet.** | Reassess only after transported H2 and base energy closure are obtained. |
| 12.11 stale scientific products | **Confirmed as a provenance issue from the catalog, not a code defect.** | Mark products with executable/source revision and physics-option metadata. Do not update scientific figures until the critical physics path is accepted. |

## 6. Proposed implementation program

### Phase A: establish non-negotiable invariants

This phase should change diagnostics and tests before it changes physical
solutions.

#### A1. Centralize elemental accounting

Add one routine that computes H, He, C, O, and metal nuclei from `rho` and
`f_sp`, including H2, molecular ions, HeH+, OH, H2O, and CO with their exact
stoichiometric multiplicities. Use it in debug assertions immediately before
and after:

- `carrier_write_back`;
- `ioniz_eq` in the steady outer loop;
- `relax_photochemical_composition`;
- restart equilibration;
- acceptance and output.

The tolerance should be based on roundoff and the nonlinear cell-solver
tolerance, not on a percentage of abundance. On failure, report the element,
cell, radius, before/after totals, and each contributing species. This should
isolate item 12.8 without relying on post-processing.

#### A2. Implement dual residual reporting

Retain the current integrated norm as a budget measure, but add
`max_j |R_j|/S_j` for each row and region. Report dimensional terms at each
worst cell. Initially make the new values diagnostic only. After calibration on
accepted atomic cases and the molecular hot Uranus, make both integrated and
local conditions part of acceptance.

Do not use `max(|R|)/max(S)` as a replacement: it is not a local relative norm.

#### A3. Create focused deterministic tests

Add tests for:

- elemental closure through one carrier write-back and one ionization sweep;
- repeated residual evaluation, requiring `F(Y)` to be independent of the
  preceding trial state;
- an H2 photoevent ledger, requiring conservation of H nuclei and charge for
  every final-state channel;
- restart write/read/write preservation of state and coupling metadata;
- the chemical-equilibrium base-H2 branch without a handoff.

These tests should be small module or short-run tests. They should not regenerate
the production benchmark outputs.

### Phase B: make the residual and production map consistent

#### B1. Remove the ghost/composition lag

Refactor residual evaluation into a state function with this logical order:

1. unpack physical unknowns;
2. solve the composition/temperature fixed point for those physical cells;
3. construct boundary face and ghost states from the converged composition;
4. evaluate columns and radiation consistently with that boundary state;
5. assemble hydrodynamic and carrier residuals from the same state.

The exact order may require one additional composition/radiation sweep because
outer ghost populations enter column integrals. The acceptance condition is
not “one more `Apply_BC` call”; it is that repeating `eval_residual` at the same
`Y` produces the same residual regardless of the prior trial history.

Avoid converting every interior cell from conservative to primitive and back
merely to fill ghosts. Preserve the physical interior bytes and write only the
ghost indices. This removes the non-idempotent round trip identified by item 10
without making last-bit perturbations part of the boundary operation.

#### B2. Add the complete update-map diagnostic

Expose a diagnostic that applies one production step without stop logic or file
output and returns `G_dt`. Compare it with the steady residual at the same state
for at least three decreasing time steps. Break down differences after the Euler,
chemistry, energy, carrier, diffusion, boundary, and filter stages.

This test decides item 3:

- If `G_dt` approaches the steady residual, tighten the local residual gate and
  solve the actual root more accurately.
- If it does not, align the differing stage before changing solver tolerances.
- If the Shapiro filter is enabled, either include its steady operator in the
  residual or prohibit it in steady-root validation. A filtered marching map
  and an unfiltered residual cannot be expected to share a fixed point.

### Phase C: redesign the lower boundary at the face

Implement the lower condition at `r_edg(0)` rather than constructing a face
state by calling a ghost-cell routine.

For subsonic inflow:

1. take the outgoing acoustic invariant from the first interior reconstructed
   state;
2. prescribe two independent reservoir quantities, preferably entropy or
   temperature plus a pressure/density datum with an explicit physical meaning;
3. solve for the boundary-face primitive state;
4. construct ghost-cell averages that reproduce that face state to the order of
   the reconstruction.

For supersonic inflow, prescribe all incoming characteristics. For outflow or
flow reversal, select the appropriate characteristic count explicitly rather
than passing the interior velocity through a softplus expression.

Validation should include:

- a stationary hydrostatic atmosphere with gravity and the same source
  discretization as `RK_rhs`;
- an outgoing small-amplitude acoustic pulse, measuring the reflected amplitude;
- a smooth manufactured inflow, demonstrating the expected grid order at the
  boundary face;
- the hot-Uranus reload-and-march experiment, comparing face flux, cell-centered
  flux, and the complete update-map defect;
- three base grids after the boundary is stable.

This is a public numerical-model change. The chosen reservoir pair and behavior
under flow reversal should be agreed before implementation because they change
the mathematical boundary-value problem.

### Phase D: finish transported molecular steady states

#### D1. Decouple transport from oxygen

Replace

```fortran
if (.not. carrier_transport_stated) carrier_transport = thereis_oxychem
```

with one molecular rule. During development, the safe rule is:

- atomic run: carrier transport is irrelevant;
- molecular run with transport explicitly false: continue only with a prominent
  invalid-closure warning and report the omitted transport/chemistry ratio;
- production molecular result: require carrier transport and a carrier residual
  gate.

After Phase D2 passes, make molecular transport the default. Oxygen chemistry
must add OH/H2O/CO carriers but must not decide whether H2 is transported.

#### D2. Add a real trust region to the coupled solver

The current Levenberg--Marquardt routine searches for a descent direction after
an ascending Newton/PTC direction, but it is not a trust-region method that
controls the accepted step by model agreement. Add a scaled trust-region method
to the four-unknown branch:

- solve the constrained linear model with a dogleg or truncated Krylov step;
- compute predicted reduction from the same scaled model used by GMRES;
- compute actual reduction from the full residual;
- update the radius using the actual/predicted reduction ratio;
- enforce carrier headroom and thermodynamic admissibility within the step;
- reject a model whose directional derivative disagrees with a finite-difference
  ray test.

The radius should be initialized from the measured model-valid arc, not from a
fixed raw norm. Report radius, predicted reduction, actual reduction, ratio,
carrier worst cell, and all physical gates at every accepted iteration.

Success means all of the following, not merely `info = 0`:

- hydrodynamic integrated and local residual gates pass;
- carrier integrated and local residual gates pass;
- the mass-flux gate passes on its stated wind window;
- elemental budgets close;
- the complete production update leaves the state fixed to the calibrated
  tolerance;
- the H2 front no longer drifts under an additional transport relaxation.

#### D3. Replace the local-equilibrium molecular references

Only after D2 succeeds should molecular reference outputs and derived science
products be regenerated. The acceptance comparison is physical and numerical:
front location, H3+ column, temperature, mass-loss rate, energy closure, grid
convergence, and restart stability. Byte identity with the old local-equilibrium
results is not an acceptance criterion for this physics correction.

### Phase E: complete independent molecular microphysics

#### E1. H2 photoionization final states

Represent the absorption cross section as a sum of explicit, mutually exclusive
channels:

```text
H2 + photon -> H2+ + e-
H2 + photon -> H + H+ + e-
H2 + photon -> H+ + H+ + 2e-
H2 + photon -> H + H
```

Store or calculate one cross section for each channel. Construct opacity from
their sum and construct species, electron, and heat source terms from a
stoichiometric table. This prevents double ionization and neutral dissociation
from being hidden inside `f_di`. Above 124 eV, either obtain a defensible model
or expose the constant branching as an uncertainty option; do not silently
present it as measured.

#### E2. Lyman--Werner geometry

First preserve the current slab table with explicit provenance and a
`slab_surrogate` model name. Then calculate a matched slab/spherical comparison
with identical molecular data. The replacement table must recover the present
slab result in the slab limit and converge in angle, frequency, and radial
resolution. Until that exists, report the geometry uncertainty rather than
assigning a correction sign.

### Phase F: convergence workflow and reconstruction

Keep three concepts separate:

- marching stop: an economical indication that time evolution has slowed;
- nonlinear-solver handoff: a criterion chosen for solver robustness;
- scientific acceptance: satisfaction of all physical residual and invariant
  gates.

They may use related observables but should have different names and calibrated
tolerances. A run stopped only by `du` must not be labeled a steady solution.

For PLM-to-WENO3 cases, do not switch the discrete operator discontinuously and
expect the existing state to remain a solution. Use PLM to form an initial guess,
then solve a continuation family

```text
R_lambda = (1 - lambda) R_PLM + lambda R_WENO3,
```

with step control in `lambda`, ending at `lambda = 1`. The final gates must be
evaluated with the pure WENO3 operator. This also separates a WENO3 instability
from a transient caused only by changing operators.

## 7. Validation matrix and stopping criteria

Validation should follow the changed path rather than rerunning every stored
calculation after each patch.

| Change | Focused validation | Broader validation required before acceptance |
|---|---|---|
| Element accounting | Carrier write-back, ionization sweep, A2 oxygen case | Molecular + metals + oxygen cases because composition is shared |
| Residual norms | Saved atomic and molecular states; verify worst-cell arithmetic by an independent script | All steady-solver regression cases because acceptance changes |
| Residual ordering | Repeated `F(Y)`, directional derivative, one-step map | Atomic, molecular, excited-H, metals, and conduction/viscosity steady cases |
| Characteristic boundary | Hydrostatic state, acoustic pulse, smooth inflow | Planet matrix because every solution sees the lower boundary |
| Trust region | Coupled hot Uranus from two initial states and two grids | Atomic JFNK cases to establish no unintended solver regression |
| H2 channels | Cross-section continuity and stoichiometric conservation | Molecular cases spanning soft-EUV and hard-X-ray spectra |
| Restart schema | Short deterministic write/read/write and one-step restart | At least one staged-secondary-ionization and one transported molecular case |
| PLM/WENO continuation | Fixed-state operator difference and continuation trace | A2 molecular wind plus one atomic two-stage case |

For every physics-changing run, record source revision or working-tree hash,
input hash, reconstruction, boundary model, carrier model, restart schema,
residual definitions, grid, and all stop flags. Stored outputs without these
fields should not be used as current regression references.

## 8. Definition of completion

The molecular defect group is complete only when one transported molecular case
and one oxygen-chemistry case demonstrate all of the following:

1. elemental nuclei and charge are conserved through every operator to the
   stated nonlinear tolerance;
2. the H2 continuity equation, including advection and diffusion, passes both
   integrated and local residual gates;
3. mass, momentum, and energy pass integrated and local gates throughout the
   domain, with the boundary cells reported separately;
4. the accepted state is a fixed point of the complete production update;
5. the lower boundary is subsonic where the subsonic model is used and shows a
   convergent face state under refinement;
6. front position, H3+ column, temperature minimum, and mass-loss rate converge
   under at least three grids or have a quantified observed order and error;
7. a restart preserves coupling state and returns to the same fixed point;
8. the result is independent of whether the initial guess came from PLM
   marching, a previous WENO3 state, or a compatible restart within the basin
   expected for the physical solution.

Until these conditions are met, the code should clearly label default-off
carrier transport and local-equilibrium molecular outputs as experimental or
closure-limited. This is preferable to treating solver nonconvergence as
permission to publish a solution known not to satisfy the transported H2
continuity equation.
