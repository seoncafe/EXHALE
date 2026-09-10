# EXHALE development plan, revision 1 (2026-09-05)

Supersedes `docs/development_plan_20260905.md` (rev 0) after the review in
`docs/development_plan_20260905_review.md`. Inputs, in order of authority for
this document: the source at HEAD `35d9dd5`, the review, the audit
(`docs/physics_numerics_audit_20260905.md`), and `docs/code_status_20260905.md`.
The review's factual claims about the source were re-checked before being
adopted (section 1). Rev 0 is kept unchanged as the record of what the review
addressed.

Nothing here is implemented. No atmosphere was run.

## 1. What the review corrected, and what was verified today

| Review item | Claim | Checked in source today | Disposition |
|---|---|---|---|
| R2 | `eval_residual` does not reconstruct the interior energy after the sweep; it derives T from the given U, sweeps, refreshes the particle count, rewrites the ghosts, and assembles the residual with the interior U unchanged | `steady_newton.f90` 555-606: `comp_T_from_p` -> `ioniz_eq` -> `get_species_densities` -> `Apply_BC` -> `assemble_residual(u, n_tot+ne, heat, cool, R)`. No `comp_p_from_T`, no `W_to_U`. Rev 0 section 1 was wrong. | Adopted. Section 3.1 below restates what each path holds fixed. |
| R3 | A two-column power law through the Table 6 entries at 1000 K has exponent 0.46, not the collisional exponent 1 | `h3p_cooling.f90` 67-70: `sLogN = 6, 8, 10, 12, 14`; the entries are the Miller Table 6 values the review quotes. The lowest tabulated columns are not in the linear asymptote. | Adopted. B2 rewritten (section 4.4). |
| R6 | `mol_sec_ion` and `mol_carrier` do not run transported H+ | `backup/regression/mol_sec_ion/input.inp`: molecular chemistry + `Secondary_ionization: Immediate`, no carrier transport. `mol_carrier`: carrier transport, no `Ionization transport`. The only inputs with `Ionization transport` are the three Koskinen benchmark directories. | Adopted. A dedicated H+ case set is Phase 0 work (section 3.3). |
| R7 | The JFNK action is a finite difference, not exact | `steady_newton.f90` 1370-1386: `Jv = (Fp - F0)/eps`, `eps = sqrt(epsilon)*(1+|Y|)/|v|`, with step halving on inadmissible states | Adopted. C4 rewritten (section 4.6). |
| R10 | `make check` passes no case list; goldens carry provenance of a dirty `6d07d48` tree | `Makefile` 294-300 calls `run_check.sh check` without arguments; the script accepts `check [case...]`. Golden headers: `git=6d07d48afd41 tree=dirty`. | Adopted. Baseline provenance rule in section 3.2. |
| R1, R4, R5 | A4 as written is not a conserving reacting update; A7's residual and derivative must describe one equation; D1/D3 must precede C1; the CO ceiling conflicts with step rejection | Design arguments, not source facts; consistent with the code read for rev 0 | Adopted. Sections 4.1, 4.3, 4.5, 4.6. |
| R8, R9, R11, R12 | B1 normalization and spin policy; A5 fixes a deficit, not a distribution; closures dropped from the phases; parallelism and sizes overstated | Design arguments | Adopted. Sections 4.3, 4.4, 6, 7. |

The review's second calculation (all twelve goldens differ from the case
outputs in numerical content) agrees with the `cmp` result of rev 0.

## 2. Principle of the revised order

The rev 0 order put a larger steady system (C1) before the transport and
boundary equations it must contain (D1, D3), and treated the energy update (A4)
as a mechanical substitution. The revised order follows the equations
downstream work consumes: derive first, then implement the bounded algebraic
fixes that need no derivation, then the thermodynamics and the coupled source
update, then the spatial operator, then the steady solver on that operator,
then validation and extension.

```
Phase 0  evidence, provenance, gates, H+ test configurations
Phase 1  bounded corrections needing no shared derivation      (A1, A2, A3)
Phase 2  the common physical system, on paper                   (derivation D0)
Phase 3  thermodynamics and the coupled source update           (A4+B3, A5, A6, A7, B1, B2, B4)
Phase 4  the spatial operator and its validity conditions       (D1, D2, D3, mixture energy transport)
Phase 5  the steady solver on that operator                     (C1, C2, C3, C4)
Phase 6  coupled validation and extension                       (Koskinen, lower/upper matching, CO, IR, observables)
```

Only Phase 0 and Phase 1 contain tasks that are independent of each other.
From Phase 2 on, the tasks share interfaces (the species list, the reaction
ledger, the energy convention, the face flux) and are implemented after the
derivation fixes those interfaces, not in parallel before it.

## 3. Phase 0: evidence, provenance and gates

### 3.1 What each code path holds fixed (corrects rev 0 section 1)

- Marching step (`EXHALE_main.f90` 1168-1212): from the post-hydro U, T is
  derived from p and the current composition; carriers are transported and
  the sweep is solved at that T; the particle count is refreshed; p is
  rebuilt at the unchanged T (`comp_p_from_T`) and U is rebuilt from it
  (`W_to_U`); then the energy step applies heat minus cool. The rebuild at
  fixed T is the projection of audit N1.
- Steady residual (`steady_newton.f90` 555-606): from the unknowns U, T is
  derived from p and the seed composition; the sweep is solved at that T; the
  particle count is refreshed; the ghosts are rewritten with the new
  composition; the hydrodynamic residual is assembled with the interior U
  unchanged. There is no interior energy rebuild. The temperature the fluxes
  see is therefore p/(n_tot+n_e) of the post-sweep composition, and the
  composition used for the rates is the seed's.
- Final closure loop (2636-2655): the post-sweep composition is fed back as
  the seed until the residual norm moves by less than 1e-3.

Consequences: (a) A4 has nothing to replace mechanically on the steady side;
the steady-side work is self-consistent evaluation of T, composition,
radiation and the energy residual at the given unknowns (Phase 5, C3). (b)
Whether the corrected marching map and the corrected residual share a fixed
point, and whether previously converged atomic solutions are preserved, are
hypotheses to be tested on the same physical equations, not consequences of
the present code. Rev 0's inference to the contrary is withdrawn.

### 3.2 Baseline provenance

- The goldens of record carry `git=6d07d48afd41 tree=dirty`. Before any
  edit, the user refreshes the twelve goldens from the post-170 case outputs
  (`run_check.sh golden <cases> && run_check.sh check`), and the old goldens
  are moved to `backup/regression/golden_pre170_<date>/` with a README that
  states the revision, toolchain, enabled physics and stopping reason of
  each case. A reference is associated with its input checksum (already in
  the header), source identity, toolchain and purpose; an output directory is
  not a reference by virtue of existing.
- Refreshing is a control for later edits, not a scientific milestone, and it
  does not block source inspection, unit tests or the bounded fixes of
  Phase 1.
- Vocabulary in every report: byte-identical; numerically equal within
  round-off (below 0.1%, from flags or constant unification); within a
  declared tolerance; moved. None of the last three is identity, and a small
  movement is never an argument for accepting a physical defect.
- `make check` runs the whole matrix. A bounded run uses the script directly:
  `backup/regression/run_check.sh check <case...>`. When the test tooling is
  touched (3.4), the Makefile target gains a `CASES=` variable.

### 3.3 Missing test configurations

- Transported H+ has no regression case. Add three to the matrix, all
  H/He-only (no oxygen), all fixed-step snapshots with `maxsteps`:
  `hp_zero_seed` (initial H+ = 0 everywhere, `Ionization transport: True`),
  `hp_trace_seed` (H+ = 1e-12 n_H), and `hp_front` (a mixed H2/H+ front from
  the hot-Uranus gate). They are the cases A3 and C1 are measured on.
- A zero-irradiation hydrostatic column (`hydrostatic_column`) for D2/D3.
- The audit driver `docs/audit_20260905/audit_probe.f90` moves to
  `src/tests/physics_probe/` and is built through the Makefile's dependency
  generator (the hand-ordered list in `run_checks.sh` is retired). Every
  printed balance acquires a tolerance and a verdict; the program exits
  nonzero on any failure.

### 3.4 Gates

- `make test`: `element_census_tests`, `diffusion_tests`, the E1 channel
  check, the A2 kinetics checks, `residual_determinism`, the physics probe.
  Nonzero exit on any failed criterion.
- A test that prints a number and does not compare it is a diagnostic, and is
  labeled so in its output.

## 4. Phases 1-5: from the bounded corrections to the steady solver

### 4.1 Phase 1: bounded corrections (A1, A2, A3)

These need no shared derivation and can be delegated in parallel.

A1, LLF signal speed (`Num_Fluxes.f90` 287): `a1 = max(abs(vL)+aL,
abs(vR)+aR)`. Reached by the LLF option and by `positivity_limited_fluxes`
(`RK_rhs.f90` 164). Tests: mirrored pair, inward supersonic pair, the audit
states (mass flux -4.5). Add a counter of positivity-repair invocations to
the run summary; the goldens' movement is predicted from that counter, not
assumed.

A2, ROE star states (`speed_estimate_ROE.f90` 76, 101-104): exponent 1/z on
the two-rarefaction bracket; compression ratios multiplied by the upstream
densities. Tests: normalization invariance (rho, p scaled by 10 at fixed v
and p/rho), exact constant-gamma shock and rarefaction. `input_read` refuses
`Numerical flux: ROE` with the caloric EOS active, with a message naming
this section; the header says the routine is valid for constant gamma.

A3, carrier reference scales and returned-state diagnostics
(`diffusive_photochemistry.f90` 376-377, 1680-1686, 1830-1838, 1889-1893):
one element-reference table for every carrier (H2, H+: hydrogen; OH, H2O:
oxygen; CO: min(O, C)) used by the Jacobian step and the residual floor;
`pct_newton_resid` evaluated on the state after `limit_to_element_budget`;
`solve_carriers` returns a status.

What the caller does with a failed status, and what "residual of the returned
state" means when the CO ceiling has acted, are settled here and not left to
step subdivision (review R5): a kinetic residual with `src(ic_CO) = 0` cannot
be satisfied by a state the ceiling has moved, and halving the step does not
change that. Phase 1 therefore defines the restricted configuration in which
CO is accepted: the ceiling is a declared algebraic constraint, the residual
reported for CO is the constraint residual (distance to the ceiling where it
binds, the transport residual elsewhere), and the setup report and
`EXHALE_resolved.out` label `Oxygen chemistry` as "CO: reservoir with
equilibrium ceiling, no kinetics". A failed status in any other carrier
rejects the step. The kinetic CO network remains Phase 6 work.

Tests: the source Jacobian at H+ = 0, 1e-12 n_H, 1e-3 n_H in metal-free gas,
directional derivatives against independently scaled perturbations, calling
the actual source assembly; the three `hp_*` cases of 3.3; a case that binds
the CO ceiling must report the constraint residual and pass or fail the step
by the stated rule.

### 4.2 Phase 2: the common physical system on paper (D0)

One document, `docs/governing_system_<date>.md`, written before the code of
Phases 3-5, fixing:

1. Species set and metadata: mass, elements, charge, formation energy
   epsilon_s^0, internal-state model, transport class. He(2^3S) is a level of
   He I, not an independent nucleus.
2. Reaction representation: stoichiometric matrix nu, progress rates R_alpha,
   S_s = sum_alpha nu R_alpha, with A nu = 0 and z^T nu = 0 checked per
   reaction. The same object yields species sources, reaction energy,
   derivatives and range diagnostics. The present hand-maintained lists
   (`molecular_reaction_heat`, the `eval_cool` consumers, the carrier
   `carrier_source`) are mapped onto it and every duplicate recipient of one
   energy term is listed.
3. Energy convention: thermal energy as the evolved variable (no change to
   the stored energy column or the restart file), with the exact
   accounting identity for a homogeneous cell,
   E_th^{n+1} - E_th^n = - sum_s epsilon_s^0 (n_s^{n+1} - n_s^n) + Q_nonchem,
   where Q_nonchem is the integrated radiative and external exchange. The
   three operations that change a cell's state are treated separately: a
   change of thermodynamic coordinates at fixed physical state (no energy
   change), a chemical reaction (the identity above), spatial transport of
   species and their enthalpy (fluxes). The total-energy convention (audit
   6.3) is recorded as the alternative and the user's decision on it stands
   (section 8).
4. Photon ledger per event: threshold, chemical product energy, electron
   energy and its degradation, fragment energy, internal excitation, escaping
   radiation. Exact sum first; the adopted distribution and its validity
   range second.
5. Species continuity in the face-flux form the hydro uses: advective flux
   from the Riemann mass flux and a reconstructed mass fraction, diffusive
   flux J_s with sum_s m_s J_s = 0 and enthalpy transport sum_s h_s J_s,
   spherical volume (r_R^3 - r_L^3)/3, base condition as a species flux
   through the actual face. Constraint count at the base equals incoming
   characteristics plus species flux conditions.
6. Steady residual as the zero of the same operator, with the independent,
   frozen and eliminated variables named, the radiation closure named, and
   the composition-convergence gate defined on species, T, optical depths
   and diffuse rates.
7. Validity diagnostics to be tracked at runtime: Knudsen number,
   ion-neutral and electron-heavy thermal equilibration times, chemical
   relaxation times from the constrained Jacobian, infrared critical
   densities, optical depths.

Acceptance of D0 is a review by the user; nothing in Phases 3-5 starts
before it.

### 4.3 Phase 3: A4 + B3, the coupled species-energy source update

Rev 0's A4 (hold thermal energy fixed, invert the EOS, add sources) is
withdrawn as the definition of the fix. It removes the particle-count
injection but does not make the reacting update conserve energy: the
composition an equilibrium sweep returns at the old T is not equilibrium at
the inverted T, and adding "binding energy of every change" on top of the
existing collisional-ionization cooling and molecular reaction heat would
count terms twice. B3 is the fix; A4 is its first increment.

Design, following D0 item 3:

- Source step for a cell: given U after the hydro stage and the transported
  composition, solve composition and thermal energy together,
  F_s = n_s^{n+1} - n_s^{*} - dt S_s(n^{n+1}, T^{n+1}) = 0 for the species
  that are integrated, algebraic equilibrium for the species that D0's
  relaxation criterion eliminates, and
  F_E = u(T^{n+1}, n^{n+1}) - u^{*} + sum_s epsilon_s^0 (n_s^{n+1} - n_s^{*})
  - dt [H_rad(T^{n+1}, n^{n+1}) - C(T^{n+1}, n^{n+1})] = 0,
  where H_rad contains only the photon-ledger heating (photoelectron, fragment,
  photolysis) and C the radiative losses. Collisional reaction energy is not
  a separate source: it is the epsilon_s^0 term. The present `P53` module and
  the ionization-potential terms inside `eval_cool` are consumed by this
  identity and removed as independent sources; D0 item 2 lists them.
- The same residual and the same iteration matrix (the derivative of F_E
  includes dH_rad/dT, dC/dT and the composition derivatives through the
  chemical Jacobian; the |dC/dT| damping of today can remain as a trial
  matrix with a residual test and a step-size fallback). Root selection: the
  root continuous with T^{*} under step halving; a cell that fails the
  tolerance is subdivided, and the count is reported (review R4).
- A7 is not a separate task: the two-iteration energy solve in
  `energy_semi_implicit.f90` is replaced by this step. Where the caller wants
  the cheap path (atomic gas, composition eliminated), the same routine
  reduces to a scalar T solve with H and C re-evaluated at the returned T
  and the residual tested.
- Increments: (i) the scalar T solve with residual test and returned-state
  cooling (this alone is what rev 0 called A7 and touches every case); (ii)
  the epsilon_s^0 accounting with the existing equilibrium sweep, i.e. the
  sweep at T^{*} followed by the energy identity and a re-sweep at T^{n+1}
  to a stated tolerance; (iii) the fully coupled step for the integrated
  species. Each increment has its own gate below; none claims the
  closed-cell gate before (iii).
- Measure before replacing: integrate the present projection energy over a
  run for one atomic and one molecular case and report it.

Gates:

- Closed reacting cell: radiation off means photoionization off and, for the
  material-energy test, radiative recombination and spontaneous emission
  either off or their escaping energy included in the measured total (review
  R1). Total thermal plus chemical energy constant to the nonlinear
  tolerance; equilibrium reached; error decreasing under dt halving
  (trajectory convergence is reported as a separate number).
- Photon event: h nu = heat + stored + escaped for each channel, at the
  threshold, at an intermediate energy, in the high-energy limit (this is
  A5's gate; A5 is implemented inside this phase).
- Homogeneous H/O/He cell with photolysis off conserves energy (A6's gate);
  with photolysis on, incident = escaped + stored + thermal.
- The four Newton-finished atomic planets: Mdot and the T, v, n_s profiles
  before and after, plotted; movement reported, not predicted.

A5 (H2 photoevent ledger, `util_ion_eq.f90` 937-968) inside this phase:
impose the exact event-energy sum; the 19.725 eV is the deficit of the
represented products, not a measured fragment energy for every event. Adopt a
stated partition (vertical picture: fragments carry the deficit; electrons
share h nu - 51.4 eV) and its validity range; the degradation share
evaluated at the mean electron energy is an approximation of a nonlinear
response and is labeled so, with the alternative (a dedicated absorber row)
recorded. Fragment thermalization uses the energy-dependent stopping length,
not the thermal mean free path.

A6 (missing reaction-heat terms): the He(2^3S)/HeH+ associative term enters
the ledger; the oxygen network's reaction energies and the O(1D) elimination
enter it together. If the oxygen ledger is not complete when the phase
closes, the option is labeled incomplete in the setup report and is not used
for oxygen-rich production runs; the label is disclosure, not closure
(review R11).

### 4.4 Phase 3, thermodynamics: B1, B2, B4

B1, one H2 thermodynamic model. The level list of `molecular_infrared_data`
is the candidate single source, subject to a completeness check over
100-10000 K (a line-derived list is not automatically a complete state sum).
From it derive Z_rv, u_rv = k T^2 d ln Z/dT, c_v, and the standard-state
functions H^0(T), S^0(T), G^0(T) with translation, standard pressure,
formation reference, statistical weights and the code's mass convention
made explicit; only then can the oxygen module's Shomate H2 entry be
replaced. Ortho/para: state the conversion timescale in the target gas and
the condition under which equilibrium populations hold; the policy is
whichever the argument supports, with the frozen ratio as the alternative.
Tests: the identity u = k T^2 d ln Z/dT; c_v > 0; K_eq identical across
modules; a closed H/O cycle returns zero net energy. The effect on a
molecular snapshot is unmeasured until run; the audit's 0.6-8% in u_rv is
not a prediction of it.

B2, H3+ cooling. Two separate gates. Physical: below the tabulated density
the non-LTE factor follows s proportional to n(H2), the exponent 1 imposed,
its normalization from a stated low-density rate model (or a matching
approximation documented as such), with a smooth transition to the table;
the table's lowest two columns are not in that asymptote (exponent 0.46 at
1000 K) and are not used to infer it. Numerical: continuity of the emission
function by construction (one smooth fit through the published data with the
joins removed), derivative continuity where the nonlinear solver needs it,
and a separately reported fit accuracy against Miller 2013 Table 5. A
relative jump below 1e-3 is not the definition of continuity and is dropped
as a criterion.

B4, rate validity metadata: temperature range and collider in every rate
routine; out-of-range counters in the run summary. No behavior change.

### 4.5 Phase 4: the spatial operator (D1, D2, D3, mixture energy transport)

Implemented from D0 items 5 and 7, before the steady system (review R5).

D1: carrier advective flux from the face mass flux and a reconstructed mass
fraction, common spherical volume, diffusion as an explicitly balanced extra
flux with sum_s m_s J_s = 0 and the enthalpy flux sum_s h_s J_s in the
energy equation; base species inflow read from the actual face flux. Tests:
frozen passive scalar (speed, order, integrated conservation); column
element budgets including boundary fluxes; the enthalpy flux closes the
column energy budget in a diffusion-only test.

D3: constraint count at the base documented in `base_boundary.f90` against
the incoming characteristics and the species flux conditions; adjusted where
the count exceeds the incoming set.

D2: base limit cycle, measured on `mol_sec_ion` and `hydrostatic_column`
after A1 and D1/D3; then local time-step control in the base cells or the
face condition revisited. Gate: no odd-even pattern; the hydrostatic column
stationary to truncation error at three resolutions.

Validity: the Knudsen and equilibration diagnostics of D0 item 7 are computed
and written per cell; the code does not yet change its closure when they
fail, and every run summary says whether they were exceeded and where.

### 4.6 Phase 5: the steady solver on that operator (C1-C4)

C1: species rows as Newton unknowns, variable in number, built from the
Phase 4 operator and the Phase 3 source step; the `input_read` refusals
(1845-1881) are lifted one at a time and each remaining one names the phase
that lifts it. Gate: a Newton-converged solution compared with a
time-converged reference of the same equations on the same grid (not with an
old fixed-step snapshot); the `hp_*` cases; atomic cases with the keys off
byte-identical.

C2: elimination by relaxation time from the constrained chemical Jacobian,
against advection, diffusion, temperature-change and irradiation-change
times; the P r/v criterion is retired as the decision rule. He+ transported
where the criterion requires; the digitized Koskinen He+ points are external
validation in Phase 6, not this gate.

C3: composition convergence gate on species, T, optical depths and diffuse
rates; the reload test (rebuild rates and columns from a converged state,
require the same state within tolerance). Where the atomic branch cannot
equilibrate on the flow time, transport off is a stated approximation, not
an equivalent option (review R11); the decision on enabling it is the
user's (section 8).

C4: the Jacobian action is a finite difference with a step, rounding and
inner-solve error; verify the directional derivative over a range of eps,
name the independent, frozen and eliminated variables, verify the inner
radiation/composition closure separately, and state what the banded
preconditioner omits. The base energy row of `wasp_full_newton` (||R||
5e-4 against 1e-5) is attacked here, after Phase 3 has changed the energy
residual.

## 5. Phase 6: coupled validation and extension

- Koskinen 2022 Model A rebuilt in the audit's order: geometry and lower
  boundary; absorbed photon count and energy at the same spectrum; heating
  and cooling terms at prescribed n_s(r), T(r); rate coefficients at those
  states; prescribed-wind transport; the coupled wind, Newton-finished, with
  spatial and source-step convergence. The Koskinen 2009 H3+ table as a keyed
  alternative after B2. The heating discrepancy is a bracket until then.
- Lower/upper matching: species and element fluxes and energy across the
  interface; a matching-pressure sweep showing an overlap region in which
  the upper solution is insensitive to the interface.
- Post-processing: molecular composition, n_e and T of the `_adv` profiles
  from the same non-equilibrium state as the hydro, or the transit
  calculation restricted to a domain where the neglected molecules are shown
  irrelevant.
- Photon conservation across cells and bands: monochromatic absorption, a
  single optical-depth jump, competing H/H2/He absorbers, the LW/oxygen
  overlap, refinement at thresholds; recombination-photon recycling counted
  once.
- CO kinetics with photodissociation and shielding, molecular heat
  capacities for H2O, CO, OH and the molecular ions where they are not trace
  species, a net-exchange infrared form zero by construction at T_gas =
  T_rad, W = 1, and band optical depth. Each with equilibrium, quenched, thin
  and thick tests.
- Multidimensional geometry, magnetic confinement and stellar-wind
  interaction are separate models and stay outside this plan.

## 6. Ownership of every audit closure

Review R11 asked that no closure disappear. Owner phase and gate for each:

| Closure | Phase | Gate |
|---|---|---|
| Reaction invariants A nu = 0, z^T nu = 0 | 2 (D0), 3 | one check for each reaction in `make test` |
| Thermodynamic identity, K_eq across modules | 3 (B1) | identity to 1e-6; cycle test |
| Closed reacting cell | 3 (A4+B3) | energy constant to tolerance |
| Photon event ledger | 3 (A5) | sum exact per channel |
| Finite-volume absorption, band overlap | 6 | absorbed = in - out |
| Infrared equilibrium cavity | 6 | zero by construction |
| H3+ low-density limit and continuity | 3 (B2) | exponent 1; continuous |
| Frozen passive scalar, species conservation | 4 (D1) | order and budget |
| Charged/neutral diffusion, zero net diffusive mass, enthalpy flux | 4 (D1) | budget closes |
| Momentum/thermal coupling validity | 4 | diagnostics written; exceedance reported |
| Riemann and mirrored flow | 1 (A1, A2) | bounds, conservation |
| Hydrostatic column | 4 (D2) | stationary to truncation |
| Source-step refinement | 3 | measured order |
| Spatial refinement | 6 | three grids, science quantities |
| Final-state residual | 5 (C3) | reload reproduces the state |
| Lower/upper interface | 6 | flux continuity, overlap region |
| Molecular post-processing consistency | 6 | same state or restricted domain |
| Molecular benchmark | 6 | Koskinen, term by term |

## 7. Parallelism and effort

Independent, delegable now: A1, A2, A3, B4, the H+ test configurations, the
probe conversion, the baseline provenance directory, the projection-energy
diagnostic. Everything from D0 onward is sequential at the interface level:
D0, then the Phase 3 increments in order, then Phase 4, then Phase 5. Inside
a phase, work on disjoint files may run in parallel once D0 has fixed the
interfaces.

No calendar estimate is given. Effort is committed per increment after its
derivation exists, and each increment ends at its gate. A6's oxygen ledger
and B2's low-density model in particular are not bounded until their
derivations are written. A schedule is never an argument for keeping
incomplete physics.

## 8. Decisions for the user

1. Golden refresh for section 170 and the archival of the old goldens with
   provenance (3.2). Bounded fixes and tests proceed regardless.
2. Energy convention: thermal energy as the evolved variable (this plan) or
   total energy including formation energies (audit 6.3). D0 is written for
   the thermal form and records the alternative.
3. D0 review and acceptance before Phase 3 starts.
4. Oxygen chemistry when its ledger is incomplete at the end of Phase 3:
   restricted with a visible label, or held until complete.
5. In-hydro ionization transport for the atomic line: where the relaxation
   criterion of C2 shows it is needed, transport off is an approximation;
   whether to enable it is the user's call in Phase 5.
6. Request to Koskinen et al. for the spectrum file and the heating per
   ionization.
7. Instruction-gated items stay gated: paper and poster re-convergence, the
   LHS 1140 b write-up, the public switch of the repository.

## 9. What this plan does not claim

- No defect's effect on Mdot, on the base cycle or on any planet's profile is
  asserted; each is measured at its gate.
- Preservation of the converged atomic solutions under Phases 3-5 is a
  hypothesis (3.1), tested on the four planets with plotted profiles.
- The Koskinen heating factor is attributed to neither code.
- Existence and stability of a molecular steady state with transported H+
  is what Phase 5 tests.
- The review's remaining design arguments (R1, R4, R5, R8, R9, R12) were
  adopted on their reasoning; their consequences for the code become facts
  only when the gates above are run.
