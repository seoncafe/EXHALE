# EXHALE development plan (2026-09-05)

Written after reading `docs/physics_numerics_audit_20260905.md` (the external
audit of revision `35d9dd5`) and `docs/code_status_20260905.md` (the in-house
state of the code), and after re-inspecting the source of the same revision.
This plan orders the work the two documents call for, resolves where they
disagree, and attaches to every item the file it touches, the acceptance test
that closes it, and the decision it requires from the user.

It does not implement anything. Nothing below was run beyond the checks of
section 1.

## 1. What was re-verified today in the source

Every defect the audit lists as "confirmed" was re-read at HEAD `35d9dd5`
(working tree clean except for the audit directory). All of them are present
as described.

| Audit item | Location (verified) | Present |
|---|---|---|
| N1 fixed-temperature composition projection | `src/EXHALE_main.f90` 1168-1212: `comp_T_from_p` -> carriers -> `ioniz_eq` -> `comp_p_from_T(T, n_tot_new)` -> `W_to_U`, then the energy step | yes |
| N2 LLF speed omits `v - a` | `src/modules/flux/Num_Fluxes.f90` 287, `a1 = max(abs(vL+aL),abs(vR+aR))`; also reached by `positivity_limited_fluxes` (`RK_rhs.f90` 164) | yes |
| N3 ROE two-rarefaction `sqrt` and shock density ratio not multiplied by upstream density | `src/modules/flux/speed_estimate_ROE.f90` 76-77, 101-104 | yes |
| N4 H+ carrier uses the oxygen reference scale | `src/modules/lower_atmosphere/diffusive_photochemistry.f90` 1680-1686 (residual floor), 1830-1838 (Jacobian step); `ic_Hp = 5` at 377 | yes |
| N5 two-iteration energy update returns the cooling of the previous iterate | `src/modules/time_step/energy_semi_implicit.f90` 211-295 | yes |
| N6 carrier residual recorded before `limit_to_element_budget` | `diffusive_photochemistry.f90` 1889-1893 | yes |
| N7 self-consistency loop stops on a 1e-3 change of the norm | `src/modules/time_step/steady_newton.f90` 2636-2655 | yes |
| P1 51.4 eV charged for H2 double ionization, no fragment-energy recipient | `src/modules/radiation/util_ion_eq.f90` 937-968; `parameters.f90` 722 | yes |
| P2/P3 H3+ density clamp and piecewise polynomial joins | `src/modules/lower_atmosphere/h3p_cooling.f90` 100-108, 134-135 | yes |
| 5.1 transported H+ refused with Newton, coupled carrier, or atomic gas | `src/modules/files_IO/input_read.f90` 1845-1881; `steady_newton.f90` packs `nvar_jac = 3` or `4` (one H2 row) only | yes |

Two facts about the repository state that shape the ordering:

- All twelve regression goldens differ from the outputs now sitting in the case
  directories (`cmp` on `Hydro_ioniz.txt` and `Ion_species.txt`, every case).
  The goldens are the pre-section-170 snapshots; the case directories hold the
  post-170 outputs. Refreshing them is the user's action and is a precondition
  for using the harness as a control during the work below.
- The steady residual (`eval_residual`, `steady_newton.f90` 428-668) follows the
  same T-from-p, sweep, p-from-T sequence as the marching loop. At a converged
  steady state the sweep returns the seed composition, so the projection of N1
  contributes nothing there; N1 is a defect of the evolution map and of every
  state in which the composition is transported or still relaxing. This matters
  for how its effect is measured (section 3, item A4).

## 2. Where the two documents disagree, and the resolution

| Topic | `code_status` | audit | this plan |
|---|---|---|---|
| First milestone | refresh goldens, then H+ Newton row | conservation and demonstrated errors first; goldens are not a scientific milestone | Goldens are refreshed first as bookkeeping (user action, minutes), because the harness is the control for every later change. The H+ Newton row moves to Phase 3, after the update map and energy convention are fixed, since a steady system built on the present composition/energy update would have to be rebuilt afterward. |
| Koskinen heating discrepancy | property of their code; EXHALE keeps h nu - I | not attributable until EXHALE's own energy issues (N1, P1, P5) are closed and the comparison is rebuilt term by term | Follow the audit's order (section 9 of the audit) in Phase 6. Until then the 2.5-3x factor is a bracket, not a finding about either code. |
| H3+ cooling | swap in the Koskinen 2009 detailed-balance table | the table swap does not fix the missing low-density limit or the joins | Do both, in that order: continuity and the n(H2) -> 0 limit first (Phase 2), the alternative table as a keyed option afterward (Phase 6). |
| Base limit cycle | task 6, local time step or revisit the face condition | do not attribute it to the LLF defect until measured | Phase 4, after N2 is fixed and its effect on the base cells is measured with the existing `mol_sec_ion` case. |
| Atomic-line ionization transport | "decide later" | H+ should work in atomic gas; molecules should not be required | Decision for the user in Phase 3 (section 8). The plan designs the species rows so the atomic case is not excluded by construction. |

Everything else in the two documents is consistent, and the audit's six stages
are adopted as the phases below with the in-house task list folded into them.

## 3. Phase 1: correct the demonstrated errors and close the energy ledger

Goal: no reachable code path evaluates a wrong expression, and the marching
update conserves thermal plus chemical energy in a closed cell. Each item is
independent of the others unless a dependency is stated, so they can be
delegated in parallel.

### A1. LLF signal speed (N2)

- Change `Num_Fluxes.f90` 287 to `a1 = max(abs(vL)+aL, abs(vR)+aR)`.
- Callers reached: the `Numerical flux: LLF` branch and
  `positivity_limited_fluxes` in `RK_rhs.f90`. The default HLLC path is not
  touched.
- Test: a mirrored left/right Riemann pair (v -> -v must mirror the flux), an
  inward supersonic pair, and the two states of the audit driver (expected
  mass flux -4.5, not -3.5). Add them to the audit probe (section 7).
- Expected movement: goldens do not use LLF and the positivity repair is not
  triggered in the matrix (to be confirmed by a counter, not assumed); if it is
  never triggered the goldens stay bitwise. Report the counter.

### A2. ROE star-state estimate (N3)

- `speed_estimate_ROE.f90` 76: replace `sqrt(...)` by `(...)**(1/z)`, `z =
  (gam-1)/(2 gam)`. Lines 101-104: multiply each compression ratio by its
  upstream density.
- Test: scale both states' rho and p by 10 at fixed v and p/rho; star sound
  speeds must not change. Exact Riemann solutions for a shock and a rarefaction
  at gamma = 5/3.
- Add to the routine header: valid for constant gamma only; the molecular
  (caloric EOS) path should select HLLC, and `input_read` should refuse
  `Numerical flux: ROE` when the caloric EOS is active until a general-EOS
  Roe solver exists. This is a restriction, not a fix, and is stated as one.

### A3. Species reference scale in the carrier solve (N4, N6)

- Put the element reference of every carrier in metadata next to
  `ic_H2 ... ic_Hp` (`diffusive_photochemistry.f90` 376-377): H2 and H+
  reference hydrogen, OH and H2O oxygen, CO min(O, C). Use that one table in
  both the Jacobian step (1830-1838) and the residual floor (1680-1686).
- Move `pct_newton_resid` after `limit_to_element_budget` and evaluate the
  residual of the returned state; return a convergence status from
  `solve_carriers` and have the caller reject or subdivide the step when it is
  not met (today reaching the iteration cap is silent).
- Test: the source Jacobian in metal-free gas at H+ = 0, 1e-12 n_H, 1e-3 n_H;
  directional derivatives against independently scaled perturbations. A case
  that activates the CO ceiling must either satisfy the post-limiter residual
  or fail the step.
- Expected movement: `mol_carrier`, `mol_sec_ion`, the hot-Uranus H+ runs move;
  atomic cases bitwise.

### A4. Energy-conserving composition update (N1)

This is the largest item of the phase and the one the later phases depend on.

- Design (the audit's "thermal-energy-only" convention, chosen here because it
  keeps the restart file and the output columns meaning what they mean today):
  after the carrier transport and the ionization sweep, hold the thermal
  energy density fixed and invert the EOS for the new temperature, instead of
  holding T fixed and recomputing p. Concretely, in `EXHALE_main.f90` 1168-1212
  the `comp_p_from_T(T, n_tot_new, ne_new, p)` step becomes
  `T = T(u_th, composition_new)` followed by `p = p(T, composition_new)`. The
  caloric EOS already provides the inversion (`caloric_eos.f90`,
  `internal_energy_per_particle` and its inverse are what
  `energy_semi_implicit` uses).
- The reaction energies then have to enter exactly once as explicit sources:
  photoelectron heat (already), collisional reaction heat (P53, already), the
  H2 photoevent ledger of A5, and the binding energy of every composition
  change the sweep makes at fixed T (ionization by collisions, recombination,
  dissociation/association). Today the last group is partly represented by
  the projection and partly by explicit terms; the derivation must list every
  term and where it now lives before the code changes. Write that derivation
  into `docs/` first (one page), then implement.
- Do the same in `eval_residual` (`steady_newton.f90` 428-668) so the steady
  residual and the marching map remain one operator. At a converged state the
  two conventions agree; away from it they do not, and the Newton finish must
  converge to the fixed point of the corrected map.
- Measure before replacing: add a diagnostic that integrates the present
  projection energy over the column and over a run (the audit asks for this).
  Report it for one atomic case and one molecular case.
- Test: a closed homogeneous reacting cell (no transport, no radiation) with
  several initial ionization and H2 fractions conserves thermal plus chemical
  energy to round-off under time-step halving; a prescribed photoionization
  event deposits h nu minus the stored chemical energy and the escaping
  radiation.
- Expected movement: atomic converged states essentially unchanged (section 1,
  second fact); measure Mdot on the four Newton-finished planets and report.
  Molecular fixed-step snapshots move; all molecular goldens are refreshed
  once at the end of Phase 1.

### A5. H2 photoevent energy ledger (P1)

- In `util_ion_eq.f90` 937-968 separate, for the double channel: the
  absorption threshold (51.4 eV, unchanged; it is the cross-section physics),
  the chemical product energy (2 I(H) + D0(H2) = 31.675 eV), the electron
  kinetic energy, and the fragment kinetic energy (19.725 eV in the vertical
  picture). Apply the electron-degradation shares to the electron energy only;
  deposit the fragment energy as heat with a stated validity condition
  (collisional gas; record it separately where the mean free path exceeds the
  scale height).
- Replace the reuse of the single-dissociative degradation row by a row
  evaluated at the mean electron energy, or quantify the error range in the
  header if the row is kept for now.
- Test: each of the four channels closes energy, nuclei and charge at the
  threshold, at an intermediate energy and in the high-energy limit; extend
  `src/tests/e1_h2/e1_h2_channel_check.f90`, which today checks nuclei and
  charge only.

### A6. Missing reaction-heat terms (P5, the associative He(2^3S)/HeH+ term)

- Add the admitted omitted term in `molecular_reaction_heat.f90`.
- For the oxygen network: either extend `molecular_chemical_heating` to the
  oxygen reactions with the O(1D) energy handled together with its elimination
  (`ionization_equilibrium.f90` around 2094-2104), or make `input_read` mark
  `Oxygen chemistry` as "energy budget incomplete" in the setup report and in
  `EXHALE_resolved.out`. Implementing the ledger is preferred; restriction is
  the fallback if the derivation is not finished within the phase.
- Test: homogeneous H/O/He cell with photolysis off conserves energy; with
  photolysis on, incident = escaped + stored + thermal.

### A7. Energy update returns the state it converged to (N5)

- `energy_semi_implicit.f90`: iterate to a stated tolerance on the energy
  residual, re-evaluate cooling at the final T, return that cooling, and flag
  a failed cell so the caller can subdivide the step. Keep the sign-safe
  damping (|dC/dT|) as the iteration matrix, but test the residual, not the
  iteration count.
- Molecular gas: the heating array carries temperature-dependent collisional
  reaction terms; include their derivative in dF/dT or state in the header
  that they are frozen during the T solve and by how much that lags.
- Test: stiff-cooling cell, a cell on the falling cooling branch, a
  dissociating cell; dt, dt/2, dt/4 convergence.
- Expected movement: every golden moves (the second cooling evaluation changes
  the returned `cool` in all cases). Measure and report; this is the reason
  the goldens are refreshed once at the end of the phase, not item by item.

### Phase 1 closure

- Run `make check` after each item for the cases the item reaches (A1/A2
  none expected; A3 the carrier cases; A4-A7 all), report the movement, refresh
  the goldens once at the end and report the refresh.
- Run `run_fcheck.sh` once.
- Record the phase in `docs/Update_EXHALE_stage1.md` as one dated section per item.

## 4. Phase 2: one thermodynamic model and a converged source step

Goal: chemistry, EOS and cooling share one H2 state model; cooling functions
have correct limits; the source step converges instead of iterating twice.

### B1. Common H2 internal-state model (P4)

- Choose the level list of `molecular_infrared_data.f90` (the EOS side) as the
  single source, because it is the level-resolved data set and the one the
  infrared emission already uses. Derive from it: Z_rv(T), u_rv(T) = k T^2
  d ln Z/dT, c_v, and the H2 formation equilibrium constant used by
  `mol_rates::keq_H_H_to_H2`. Retire `q_rovib_H2` in `mol_rates.f90`.
- Ortho/para: state the policy in the module header (equilibrium populations
  at the gas temperature unless a conversion timescale argument is written
  down). Formation energies stay separate from thermal excitation.
- The oxygen reverse rates use Shomate data for H2; replace the H2 entry with
  the same Z_rv so a cycle that crosses modules satisfies detailed balance.
- Test: the identity u = k T^2 d ln Z/dT to 1e-6 relative on 100-10000 K;
  c_v > 0; EOS round trips; K_eq identical between the H/He and the oxygen
  modules; a closed-network cycle (H2 + O -> OH + H, OH + H2 -> H2O + H, ...)
  returns to its start with zero net energy.
- Expected movement: molecular snapshots move at the percent level (the
  table of the audit shows 0.6-8% in u_rv); atomic bitwise.

### B2. H3+ cooling continuity and low-density limit (P2, P3)

- Replace the four polynomial branches by one smooth positive function fitted
  to the published emission with the joins removed; document the maximum
  deviation from Miller 2013 Table 5 in the header.
- Below n(H2) = 1e6 cm^-3 continue the non-LTE factor so that the cooling
  scales as n(H3+) n(H2) (collisional limit) instead of a constant. The first
  version can be a power-law continuation from the lowest two table columns;
  record out-of-range evaluations in a counter reported in the setup summary.
- Test: relative jump < 1e-3 at 300, 800, 1800 K; cooling -> 0 as n(H2) -> 0
  at fixed T and n(H3+).
- Expected movement: the Koskinen gate's outer temperature (H3+ dominates the
  cooling there); `mol_*` goldens.

### B3. Coupled species-energy source step

- A7 leaves a scalar T solve with the composition frozen. Where reactions are
  stiff (dissociation front, molecular base) integrate composition and thermal
  energy together in one implicit step with the same source Jacobian A3 makes
  species-aware. Design it as the source operator both the marching loop and
  the residual call, so Phase 3 does not have to redo it.
- Test: an isolated cell against an independently integrated reference
  network (e.g. a stiff ODE solve of the same rates in Python); temporal order
  measured, not assumed.

### B4. Rate validity metadata

- Every rate routine in `mol_rates.f90`, `oxygen_rates.f90`,
  `water_photolysis.f90`, `charge_exchange.f90` states its temperature range
  and collider; evaluations outside the range increment a counter printed at
  the end of the run. No behavior change otherwise.

## 5. Phase 3: the steady system with the transported species

Goal: a Newton-finished molecular wind whose transported H+ (and later He+)
survives the finish; the same residual for marching and steady solve.

### C1. H+ continuity row as a Newton unknown (in-house task 2)

- Extend the unknown packing of `steady_newton.f90` (today `nvar_jac` = 3 or
  4) to a variable number of species rows, with the H+ row built from the same
  advection operator and the same `mol_heh_rows` source the carrier operator
  uses. Band widths `kl_jac`/`ku_jac` follow the row count.
- Remove the `input_read` refusals (1845-1881) one at a time as each
  combination becomes valid; the message for each remaining refusal must say
  which phase lifts it.
- Test: (a) Koskinen gate Newton-finished with x(H+) equal to the marching
  transported value where P r/v < 1 and to the local root where P r/v >> 1
  (limit test on `examples/18_oxygen_chemistry`); (b) the layer closed to
  ||R|| 1e-5 with base mass flux equal to the wind's; (c) atomic cases bitwise
  with the key off. These are the in-house acceptance criteria and are kept.

### C2. He+ and the molecular ions: transport or eliminate by relaxation time

- Compute the chemical Jacobian's relaxation times per cell (null modes
  excluded) and the Damkoehler number Da = (L/|v|)/tau against the advection
  time; print them as a diagnostic column first.
- Transport He+ as a carrier (the `ic_Hp` template of section 169) where Da
  is small in the target runs; keep H2+, H3+, HeH+ in local equilibrium only
  where the diagnostic shows Da >> 1 throughout the domain, and say so in the
  setup report.
- Test: He+ on the Koskinen gate within the plotted points; element census
  closed to 1e-10; key-off bitwise.

### C3. Composition convergence gate (N7)

- The self-consistency loop (`steady_newton.f90` 2636-2655) must test the
  change of every species fraction, T, the optical depths and the diffuse
  radiation rates against their tolerances, not the 1e-3 change of the norm
  alone. Report the maximum local imbalance as well as the norm.
- Test: reload a converged state, rebuild rates and columns, and require the
  same species, T, heating and cooling within tolerance (the "final-state
  residual" test of the audit's matrix). `make residual_determinism` remains
  the separate determinism check.

### C4. Nonlocal radiation in the residual

- Attenuation couples cells; either its response enters the Jacobian action
  (JFNK already applies the full residual, so the action is exact; the
  preconditioner is banded) or the radiation iteration is driven to a stated
  tolerance inside `eval_residual`. State which, and measure the Newton
  convergence rate on `wasp_full_newton` before and after.
- This is also where in-house task 7 (the base energy row holding
  `wasp_full_newton` at ||R|| 5e-4) is attacked, after A4/A7 have changed the
  energy update. Acceptance: 1e-5 on the four planets; Mdot movement reported.

## 6. Phase 4: boundaries and discrete transport consistency

### D1. Species advection from the hydrodynamic face mass flux (audit 5.4)

- Build the carrier advective flux from the face Riemann mass flux and the
  reconstructed species mass fraction, with the same spherical volume
  `(r_right^3 - r_left^3)/3` the hydro uses, in place of the cell-centered
  derivative with `r^2 dr`. Diffusion enters as an explicitly balanced extra
  flux with the sum of mass-weighted diffusive fluxes zero.
- The base species inflow then reads the actual face flux rather than the
  wind-region average of rho v r^2.
- Test: frozen passive scalar in a smooth wind (speed, order, integrated
  conservation); column element budgets close including boundary fluxes.

### D2. Base limit cycle (in-house task 6, audit N2 caveat)

- Only after A1, measure `mol_sec_ion` at the default CFL. If the odd-even
  pattern persists, apply local time-step control in the base cells (CFL below
  0.15-0.3 there, per section 161) or revisit the characteristic face
  condition of section 152. Acceptance: base velocity free of the pattern;
  goldens refreshed once.
- Test: a zero-irradiation hydrostatic column stays stationary to truncation
  error over 1e4 steps at three resolutions.

### D3. Boundary constraint count

- Document, in `base_boundary.f90`, which quantities are imposed at the base
  and show that the count equals the number of incoming characteristics plus
  the species flux conditions; today rho or p, an isothermal ghost, and the
  velocity from the mass flux are imposed, and each transported carrier adds a
  Dirichlet inflow. Adjust where the count exceeds the incoming set.

## 7. Cross-cutting: tests as gates

- Add a `make test` target that builds and runs, with nonzero exit on any
  failed criterion: `element_census_tests`, `diffusion_tests`, the E1 channel
  check, the A2 kinetics checks, `residual_determinism`, and a production copy
  of the audit probe (`docs/audit_20260905/audit_probe.f90` moved to
  `src/tests/physics_probe/` and extended item by item in Phase 1). Today
  several of these print balances without enforcing them; each acquires a
  tolerance and a verdict.
- The probe should compile from the Makefile's dependency generator, not from
  the hand-ordered list in `run_checks.sh`.
- Regression policy is unchanged: a refactor is bitwise; a physics fix moves
  goldens and the movement is reported; movement below 0.1% from flags or
  constant unification is identical and never refreshes a golden.
- Naming: new routines are named for the physics they compute
  (`energy_conserving_composition_update`, `h2_rovibrational_partition`,
  `h3p_collisional_limit`), not for their role in the call sequence.

## 8. Decisions that belong to the user

1. Golden refresh for section 170 before Phase 1 starts (twelve cases;
   `run_check.sh golden ... && run_check.sh check`).
2. Energy convention beyond Phase 1: keep the thermal-energy formulation
   (chosen here, no file-format change), or move to total energy including
   formation energies (audit 6.3; changes the meaning of the stored energy
   column and the restart file). The plan proceeds on the thermal form.
3. Whether the atomic line gets in-hydro ionization transport (Phase 3, C1
   lifts the "needs molecular chemistry" refusal or keeps it).
4. Oxygen chemistry when the energy ledger is incomplete: implement (A6
   preferred) or restrict the option with a visible label.
5. Ask Koskinen et al. for the mean solar spectrum file and the heating per
   ionization (in-house task 4). Phase 6 reconstructs the comparison either
   way, but the bracket closes only with their numbers.
6. Instruction-gated items stay gated: paper and poster re-convergence, the
   LHS 1140 b write-up, the public switch of the repository.

## 9. Phase 5 and Phase 6: extension and external validation

Phase 5 (after Phases 1-4): a kinetic CO formation/destruction network with
photodissociation and shielding in place of the equilibrium ceiling (P6);
rovibrational heat capacities for H2O, CO, OH and the molecular ions when they
are not trace species; a net-exchange infrared form that is zero by
construction at T_gas = T_rad, W = 1 (P7); band optical depth. Each with the
equilibrium, quenched, thin and thick limit tests of the audit's matrix.

Phase 6: the Koskinen 2022 Model A comparison rebuilt in the audit's order:
geometry and lower boundary; absorbed photon count and energy at the same
spectrum; heating and cooling terms at prescribed n_s(r), T(r); rate
coefficients at those states; prescribed-wind transport; then the coupled
wind, Newton-finished (C1), with spatial and source-step convergence shown.
The Koskinen 2009 H3+ detailed-balance table enters here as a keyed
alternative (in-house task 5) once B2 has fixed the function's limits. Only
then do the published profiles become acceptance targets, and the remaining
differences are listed with their sizes as `docs/koskinen2022_model_a_comparison.tex`
does now.

## 10. Order, dependencies, and size

```
Goldens (user) -> Phase 1 {A1, A2, A3, A5, A6, A7 independent; A4 last}
               -> Phase 2 {B1, B2, B4 independent; B3 after A4, A7}
               -> Phase 3 {C1 after A3, A4, B3; C3, C4 independent; C2 after C1}
               -> Phase 4 {D1 after C1; D2 after A1; D3 independent}
               -> Phase 5, Phase 6
```

Rough sizes (one "session" is a working day of the advisor/worker pattern):
A1, A2, A6, B2, B4, C3: small, one session or less each. A3, A5, A7, D2, D3:
one to two sessions. A4, B1, B3, C4, D1: two to four sessions each, with a
written derivation before code. C1, C2: the largest, four to eight sessions,
because the unknown packing, band structure, preconditioner and refusals all
change. Phases 5 and 6 are not sized here.

Housekeeping from `code_status` task 8 (`interp_ic.py` under `src/utils/`,
the callerless `Cool_coeff` wrappers and `_func` names, the H I / He I literal
constants, regenerating the manual and the `.tex` changelog from the `.md`, a
spectrum file with the 1-5 A part) is done in the phase whose work touches
each file, not as a separate pass.

## 11. What this plan does not claim

- No defect's effect on a planet's mass-loss rate is asserted; each item's
  expected movement is a hypothesis to be measured with the harness.
- The base limit cycle's cause is not attributed to any one item.
- The Koskinen heating factor is not attributed to either code.
- Existence or stability of a molecular steady state with transported H+ is
  what Phase 3 tests, not what it assumes.
