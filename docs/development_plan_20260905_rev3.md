# EXHALE development plan, revision 3 (2026-09-05)

Supersedes `docs/development_plan_20260905_rev2.md` after the third review,
`docs/development_plan_20260905_review3.md`. Inputs, in order of authority:
the source at HEAD `35d9dd5`, the published papers in `references/`, the
three reviews, the audit (`docs/physics_numerics_audit_20260905.md`), and
`docs/code_status_20260905.md`. Rev 0 to rev 2 are kept as the record of
what each review addressed.

Scope of today's verification (corrects the over-broad opening claim of
rev 2): every **source** claim of the third review was re-read at HEAD and
its diagnostics were re-run (`docs/audit_20260905/run_review3_checks.sh`);
of the **paper** claims, the Miller Table 4 and Table 6 entries, the Yan page
number and the Koskinen equations B9 and B10 were checked in the local PDFs.
The Chung, Dalgarno, Glassgold and Schulik and Booth readings reported by
the earlier reviews were not repeated today; they are cited as the reviews
report them and marked for confirmation when D0 is written.

Nothing here is implemented. No atmosphere was run.

## 1. What the third review established, and what was verified today

| Item | Claim | Verified today | Disposition |
|---|---|---|---|
| T1 | `Num_flux` uses the ROE star values only for the entropy correction; NaN star speeds leave it inactive and a finite flux results; the production ROE flux is not positivity preserving for a vacuum-producing expansion | `Num_Fluxes.f90` 170-175: `speed_estimate_ROE(...,v_star,aL_star,aR_star,...)`, then `l1R = v_star - aL_star`, `l3L = v_star + aR_star`, nothing else. Re-run: `(1,-1,1)/(1,1,1)` gives NaN star speeds and flux `(0, 0.5858, 0)`; `(1,-10,1)/(1,10,1)` gives finite star speeds, flux `(0, 41.84, 0)`, and one isolated Euler update at Courant 0.113 leaves thermal energy density `-2.93` | Adopted. A2 rewritten with a caller contract (section 4.1). |
| T2 | The `retry_step` loop covers the hydro RK stages only; the carrier step runs later; `carrier_write_back` follows `solve_carriers` without a status test | `EXHALE_main.f90` 1039-1145 (`retry_step: do ... enddo retry_step`), carrier call at 1189, `ioniz_eq` 1198, energy 1228; `diffusive_photochemistry.f90` 616 `solve_carriers`, 621 `carrier_write_back` unconditionally | Adopted. A3 split into scale/status (bounded) and the retry controller (section 4.1). |
| T3 | `pct_co_ceiling` is reset on every limiter call and incremented only above `limit_report = 1e-10`; the resolved configuration is rewritten at the end of the run | `diffusive_photochemistry.f90` 2022 (`pct_co_ceiling = 0`), 2075 (`if (over .gt. limit_report) ...`); `write_resolved_config` called at `EXHALE_main.f90` 765, 774 and 1897 (end) | Adopted. Cumulative activation record and whole-atmosphere restriction (section 4.1). **Closed 2026-09-06 (item CEILING-DEL):** the ceiling, `pct_co_ceiling` and the cumulative counters are deleted; the destruction model's domain record replaced them in the same reporting paths. |
| T4 | Table 6's last-density emission at 5000 K is 4.7587e-18, 0.998 of the scaled LTE and 1.395 of the unscaled; Table 4 begins at 500 K while the code evaluates down to 30 K | Miller PDF Table 6 row 5000 K, column 1e20 m^-3: `0.47587e-17`; Table 4 first row 500 K; `h3p_cooling.f90` clamps at 30 K | Adopted. B2's emission choice and temperature coverage (section 4.4). |
| T5 | "Escaping" radiation needs a local ownership rule; a two-cell emission/reabsorption example creates 10 eV if only domain-escaping emission is subtracted; C1 needs the independent species space named | Arithmetic re-run (`review3_ledger_checks.py`): correct material sum 0, wrong rule 10 eV | Adopted. D0 items 3 and 6 (section 4.2). |
| T6 | The audit sources are to stay in `docs/audit_20260905/`; `photoevent_energy_check.py` evaluates a superseded equation and cannot become an assertion; the E1 program has six sections, not two | `docs/audit_20260905/README.md` lists nine retained files and asks that they remain; `src/tests/e1_h2/e1_h2_channel_check.f90` sections: channel sum, stoichiometry, joins, `off` reproduces the old split, rate-weighted ledger, sigma continuity | Adopted. Section 3.3 rewritten. |
| T7 | `frac_H2_double_of_protons` is used only inside `h2_photo_channels`; consumers call `h2_channel_cross_sections`; channel arrays are stored by `set_energy_vectors.f90` and integrated in `util_ion_eq.f90` | `grep`: the only external callers of `h2_channel_cross_sections` are `set_energy_vectors.f90` 221 and the E1 test; no external caller of `frac_H2_double_of_protons` | Adopted. A8 gains the `q_D` admissibility rule, the rename and the integration check (section 4.1). |
| T8 | `n_faces_flux_positivity_limited` increments only when `positivity_limited_fluxes` returns `repaired`; Yan's sentence is on p. 1048; Koskinen B9-B10 show the electric term and the zero diffusive mass condition, not a zero-current equation | `RK_rhs.f90` 318-320: `if (repaired) n_faces_... = ... + n_repl`. Yan PDF: page header `1048 YAN, SADEGHPOUR, & DALGARNO` precedes the "vertical threshold of 51.4 eV" line. Koskinen PDF: B9 is the species momentum/diffusion equation, B10 is `sum rho_s w_s = 0` | Adopted. Corrections in sections 4.1, 4.2 and 4.3. |

The review's confirmation that the rev 2 energy identity gives `Delta u_th =
6.4016 eV` and `Delta u_th + Delta E_chem = 20 eV` for H photoionization at
20 eV was reproduced.

## 2. Principle of the order (unchanged)

```
Phase 0  evidence, provenance, gates, test construction
Phase 1  bounded corrections needing no shared derivation      (A1, A2, A3-scale, A8)
Phase 2  the common physical system, on paper                   (derivation D0)
Phase 3  thermodynamics and the coupled source update           (A4+B3, A5, A6, B1, B2, B4)
Phase 4  the spatial operator and its validity conditions       (D1, D2, D3, mixture and charged transport)
Phase 5  the steady solver on that operator                     (C1, C2, C3, C4)
Phase 6  coupled validation and extension                       (Koskinen, matching, CO kinetics, IR, observables)
```

What moved since rev 2: the rollback/retry contract is attached to the first
implementation that introduces failure recovery (A3's second increment or
A4+B3), not to a later D0 paragraph; the H3+ data selection becomes a named
prerequisite of B2; CO activation history and the restricted interpretation
of affected runs accompany any use of the ceiling; A2 acquires a caller
contract.

## 3. Phase 0: evidence, provenance, gates, test construction

### 3.1 What each code path holds fixed

Unchanged from rev 2 section 3.1: the marching step rebuilds p at fixed T
after the sweep (audit N1); the steady residual keeps the interior U and
rewrites the ghosts; inside `assemble_residual`, `Reconstruct` and the
viscous branch's `U_to_W` recompute pressure from U with the post-sweep
caloric composition, which differs from the seed-derived pressure only under
the caloric EOS; the closure loop stops on a 1e-3 change of the norm. The
seed-dependent quantities are the optical depths and rate coefficients; the
local chemistry evaluates trial populations; the returned heating and
cooling are of the returned composition.

### 3.2 Baseline provenance and vocabulary

Timing (user decision of 2026-09-05, replacing the rev 2 order): the golden
of record is refreshed ONCE, at the end of Phase 1, after every bounded
correction has landed; a refresh before them would be redone after them.
Until then the control for byte-identity is a scratch baseline
(`backup/regression/baseline_post170_<date>/`, a copy of the post-170 case
outputs) read through a new `REGRESSION_GOLDEN_DIR` variable of
`run_check.sh`; `golden/` itself is not touched. At the refresh: archive
`golden/` to `golden_pre170_<date>/` and verify the copy before
`run_check.sh golden` overwrites it; mark unknown provenance unknown;
promotion requires recorded identity (input checksum, binary md5,
`tree=clean` at the intended revision) and purpose. Vocabulary is
byte-identical / within the declared tolerance (1e-3 by default is a
tolerance, not round-off) / moved; a bounded run is `run_check.sh check
<case...>`.

### 3.3 Test construction (T6 corrected)

- **Retention.** The diagnostic sources stay in `docs/audit_20260905/` as
  the record of what each review measured at `35d9dd5`, with their README.
  Maintained tests are **derived** from them into `src/tests/physics_probe/`,
  each stating its origin and how its acceptance differs from the probe. The
  probes are not moved and are not converted mechanically: a probe that
  evaluates a superseded equation as a counterexample
  (`photoevent_energy_check.py`) is not an assertion of the production code;
  the maintained test evaluates the accepted equation or the production
  routine.
- **Acceptance rules by reference type.** An exact analytical reference
  (Riemann solution, energy identity) takes a relative tolerance near
  round-off with an explicit NaN and sign check; an approximate estimator
  (the ROE star state) is compared at the accuracy of its approximation; a
  published measurement takes its stated uncertainty (Chung 4-5%; Miller
  Table 4 fit errors); a solver residual takes its declared tolerance; an
  invariant (elements, charge, channel sum) takes round-off. A tolerance
  intended for a fit never covers an invariant.
- **E1** (`src/tests/e1_h2/e1_h2_channel_check.f90`) already has six sections:
  channel sum, stoichiometry per event, continuity across the joins, `off`
  reproducing the pre-E1 two-way split, a rate-weighted ledger over a
  power-law spectrum, and continuity of `sigma_H2`. What it lacks is the
  detector-observable inversion; that section is added (below), and none of
  the six is described as absent.
- **Chung detector inversion**: with `sigma_i`, the signal ratio `r` and
  `q_D = N_D/(N_S + N_D)`: `sigma_M = sigma_i/(1+r)`, `sigma_S = (1-q_D) r
  sigma_i/(1+r)`, `sigma_D = q_D r sigma_i/(1+r)`; reconstruct `(S+D)/M = r`
  from the returned channels; multiplicities from the stoichiometric table.
- **ROE**: (a) nonvacuum equal-pressure expansion `(1,-1,1)/(1,1,1)`: finite,
  positive star pressure and densities, `c_*` within the estimator's
  tolerance of 0.9577; (b) vacuum-producing expansion `(1,-10,1)/(1,10,1)`:
  the **final flux** from `Num_flux` and one conservative update must leave
  positive density and thermal energy, with and without the production
  repair path (`positivity_limited_fluxes`) engaged, and the two results
  reported separately; (c) stationary contact, asymmetric shock and
  rarefaction, normalization invariance. A gate demanding strictly positive
  star density everywhere is wrong for the vacuum case; the vacuum test
  checks zero pressure and density at the edges instead.
- **Transported H+** (`hp_zero_seed`, `hp_trace_seed`, `hp_front`): built
  from verified restart files because the carrier constraint applies only
  when `bg_ready` or a restart is loaded; state asserted immediately before
  `carrier_source` and the Jacobian assembly; derivative test and front test
  separate.
- **Hydrostatic column** (mechanical): fixed composition, cooling, chemistry
  and conduction off, thermodynamically stationary state with compatible
  boundaries; discrete pressure-gravity residual and spurious velocity
  measured. Radiative-chemical equilibrium is a separate test.
- **Two-cell radiation exchange** (T5): one cell emits, the other absorbs
  fully; the material sum over the two cells is zero; a rule that subtracts
  only domain-escaping emission fails this test by the photon energy.
- `make test` aggregates `element_census_tests`, `diffusion_tests`, E1, the
  A2 kinetics checks, `residual_determinism` and the derived probe tests;
  nonzero exit on failure.

## 4. Phases 1-5

### 4.1 Phase 1: bounded corrections (A1, A2, A3-scale, A8)

**A1, LLF signal speed** (`Num_Fluxes.f90` 287): **Landed 2026-09-05 (batch 2b, `Update_EXHALE_stage2.md` section 5).** `a1 = max(abs(vL)+aL,
abs(vR)+aR)`. Reporting (T8): the existing `n_faces_flux_positivity_limited`
counts **repaired interfaces in successful repair calls** (it increments only
when `repaired` is true); it is not a call count, and a repaired attempt can
still be rejected by a later RK stage. A1 adds three separately named
quantities: invocations of `positivity_limited_fluxes`, repaired interfaces,
and repaired interfaces in accepted steps, plus a count of direct LLF use;
each is printed even when zero. None predicts the size of an atmospheric
change.

**A2, ROE** (`speed_estimate_ROE.f90`, `Num_Fluxes.f90` 138-215). The
exponent `1/z`, the upstream densities on the shock ratios, and
admissibility-based branch selection remain. The caller contract is now
part of the task (T1):

- `speed_estimate_ROE` returns a status (`nonvacuum`, `vacuum`,
  `inadmissible`) with the star quantities; in the vacuum case it returns
  the two gas-edge speeds `v_L + 2c_L/(gamma-1)` and `v_R - 2c_R/(gamma-1)`
  under names that say so, not as "sound speeds", and zero star pressure
  and density.
- `Num_flux` uses the status: for `nonvacuum` the entropy correction as
  today; for `vacuum` a flux that is positivity preserving in this limit
  (HLLE with the edge speeds as the wave bounds is the candidate; the choice
  is made in the task and justified by test (b) of 3.3); for
  `inadmissible` the two-rarefaction estimate is retried and, failing that,
  the vacuum branch is used.
- Because the routine's public signature changes, the concrete interface is
  presented to the user before implementation (section 8, item 3).
- `input_read` refuses `Numerical flux: ROE` with the caloric EOS; the
  header states the constant-gamma validity.

Gate: the three ROE tests of 3.3, with the final flux and the conservative
update as the tested objects, not the auxiliary routine alone.

**Landed 2026-09-05 (batch 2a, `Update_EXHALE_stage2.md` section 4).** Implemented as approved in
`docs/a2_roe_interface.md`, with four accepted deviations: a positivity guard
on the input states at the top of the estimator, the `Numerical flux: ROE`
refusal placed in the consistency block after the optional keys, the old
`p_min` floor removed as superseded by that guard, and the HLLE counter
printed for ROE runs only. Thirty ROE assertions of
`src/tests/physics_probe/riemann_wave_speeds.f90` are green.

**A3, carrier reference scales and returned-state diagnostics**, split (T2):
**Landed 2026-09-05 (batch 2b, `Update_EXHALE_stage2.md` section 5; increment (a) only).** 

- Increment (a), bounded and delegable now: one element-reference table for
  every carrier used by the Jacobian step and the residual floor;
  `pct_newton_resid` evaluated on the state after `limit_to_element_budget`;
  `solve_carriers` returns a status; a failed status in an unconstrained
  carrier **terminates the run with a message naming this section**. This
  increment does not claim rejection or retry.
- Increment (b), the retry controller, is the first implementation that
  introduces failure recovery and therefore carries the rollback contract of
  D0 item 8 with it. Two designs are possible, and the choice is made before
  coding: (1) a carrier-substep retry inside `photochemical_transport_step`
  with a frozen background, a saved entry state, no `carrier_write_back`
  from a failed trial, and accepted substeps covering the original interval;
  (2) an outer attempted-step controller that saves and restores every
  operator already executed after the hydro stages (element diffusion, the
  transport backgrounds, rate and column caches) and re-enters at the
  hydro stage. The present `retry_step` loop covers the RK stages only and
  cannot be assumed to cover chemistry. Test: force a solver failure after
  mutating a trial carrier state; assert that rejected species never reach
  `carrier_write_back`; compare the accepted result with an explicitly
  subdivided calculation; check the elapsed physical interval and the
  restored background.

**CO ceiling while it exists** (T3). **Closed 2026-09-06.** D0 item 9 was
supplied (item B3b-CO: He+ + CO at the UMIST RATE22 4068 rate and shielded
photodissociation on the Lyman-Werner beam, with the domain measured cell by
cell), and the ceiling was then deleted (item CEILING-DEL), so nothing in this
paragraph describes the current tree. It is kept as the record of the
requirement. Until D0 item 9 supplies a justified
model, a run in which the ceiling acts is a **diagnostic experiment for the
whole atmosphere**, not a run with unreliable C/O columns only: the ceiling
changes atomic C and O available for cooling, rescales the other carriers
through the limiter, and changes chemical energy without a source, so
temperature, flow and mass loss depend on it. The record must be cumulative
over accepted evolution: an activation indicator, the amount of CO removed
(integrated over cells and time), and the cells and times affected; rejected
trials counted separately once retries exist. The present `pct_co_ceiling`
is reset on every call and counts only corrections above 1e-10 relative,
so a zero at the end proves nothing; the condition "the ceiling never acted"
is tested on the cumulative record, with the recorded magnitude available to
explain round-off-level corrections without redefining the condition. The
setup report states that the configuration permits an unvalidated ceiling;
actual activation is written to the final status and to
`EXHALE_resolved.out` through the end-of-run `write_resolved_config` call
(`EXHALE_main.f90` 1897), which needs the cumulative state supplied to it.
The report of such a run names which quantities, if any, are independently
supported.

**A8, Chung detector inversion** (`h2_photo_channels.f90`), extended (T7):
**Landed 2026-09-05 (batch 2b, `Update_EXHALE_stage2.md` section 5).** 

- Header and inversion rewritten to the detector observable `(N_S +
  N_D)/N_M`; the "overcount" note removed; Chung's "about 20% at 80 eV" read
  as the event fraction `q_D`.
- `yan_rho`: from Yan's `rho_D = N_D/(N_M + N_S)` and the detector ratio
  `r`, `q_D = rho_D/(1 + rho_D) * (1 + r)/r`. The result must satisfy `0 <=
  q_D <= 1`; an incompatible pair of input models produces a diagnostic and
  a stated model restriction, **not** a silent clip that redefines the
  measured observable.
- `frac_H2_double_of_protons` is renamed for what it returns (an event
  fraction among dissociative detector events); its only use is inside the
  module, so the rename is internal. External consumers call
  `h2_channel_cross_sections` (`set_energy_vectors.f90` 221 and the E1 test).
- Integration check: at a fixed state and selected photon energies, follow
  the channel arrays stored by `set_energy_vectors.f90` into the ionization,
  secondary-ionization and heating integrals of `util_ion_eq.f90`: one H2
  destroyed per ionizing event, two protons and two electrons per D event,
  opacity equal to the channel sum, and the consumers reading the corrected
  arrays. A8 changes integrated heating through changed probabilities even
  though the energy assignment per channel is A5's; the check does not claim
  A5's defects are corrected.

Gate: the detector test of 3.3; E1's six existing sections still pass; the
integration check. The local channel fractions are already measured (0.0225
against 0.0441 for D at 80 eV); the effect on an atmosphere is unmeasured
until run.

Independent and delegable in parallel: A1, A2 (after the interface decision),
A3 increment (a), A8, the test construction of 3.3, the baseline archive of
3.2, the projection-energy diagnostic.

### 4.2 Phase 2: the common physical system on paper (D0)

Items as in rev 2 (species and metadata; one ledger for reactions and
photoevents; energy convention; thermal enthalpy and chemical flux; charged
transport closure; species continuity in face-flux form; boundary constraint
count; time integration and stationary residual; CO constraint model;
validity diagnostics), with these amendments:

3\. **Energy convention, local ownership of radiation** (T5). "Escaping"
means leaving the **local material system of the cell**, not the
atmosphere. Emission is a loss of the emitting cell; absorption is a gain of
the absorbing cell; for an explicit radiation field its energy change and
boundary flux enter the domain budget; for an eliminated diffuse field or a
local recycling approximation the corresponding net material exchange is
derived and the same photon is not reintroduced as an independent source.
The two-cell test of 3.3 is the gate.

5\. **Charged transport closure** (T8 corrected). Koskinen et al. (2022)
Appendix B, eq. B9, is the species diffusion equation with the electric term;
eq. B10 is the zero net diffusive mass condition `sum rho_s w_s = 0`. They
do not state a zero-current equation. Zero net current `sum z_s e J_s = 0`
with the ambipolar field is a candidate closure that D0 **derives** for the
EXHALE model, or D0 justifies the common-velocity approximation over an
explicit domain; neither is attributed to the paper.

6\. **Independent species variables for C1** (T5). D0 names the independent
unknowns and the algebraic constraints (element totals, charge neutrality,
total mass) so that species continuity, mass continuity and the constraints
do not form a redundant stationary system; the constrained chemical
Jacobian is written on that independent space before any global Newton row
is built.

8\. **Rollback contract** attaches to A3 increment (b) or to A4+B3,
whichever lands first (section 4.1).

9\. **CO constraint model**: unchanged alternatives (one-sided instantaneous
destruction with sign, products, energy and a timescale/domain argument, or
exclusion), plus the cumulative activation record of 4.1 as the observable
that decides whether a run fell in the excluded class. An equilibrium
constant alone does not give the destruction timescale; the domain argument
is required for the model to be more than a diagnostic.

D0 acceptance is the checklist of the third review's section 5: material and
radiation inventory with a defined local boundary; independent variables and
constraints with exact invariants and consistent forward/reverse
thermodynamics; thermal plus chemical finite-volume equations summing to the
material equation with the same face fluxes; charged closure derived or the
approximation justified over a domain; CO validity argument or exclusion;
internal populations consistent across EOS, chemistry, degradation and
emission within an accepted range; boundary count for the full independent
system with physical time, rollback and stationary residual specified
separately; a verification matrix with units, tolerances and references for
every closure.

### 4.3 Phase 3, sources: A4+B3, A5, A6

Unchanged from rev 2 in structure: the coupled species-energy source step in
the external-exchange convention; increments (i) scalar T solve, (ii)
equilibrium composition solved with the energy constraint (may claim the
closed-cell gate for its equilibrium model), (iii) finite-rate integration
of the species not eliminated; gates for isolated photoionization,
recombination, photon-driven cycle, closed reacting cell, elementary
one-cell photon balance, source-step order, and the four atomic planets.

Amendments:

- The closed-cell and photon gates use the local ownership rule of D0 item
  3; the two-cell exchange test runs beside the one-cell balance.
- A5: the Yan vertical-threshold sentence is on p. 1048, section 4 (rev 2
  said 1049). The rest of A5 stands: exact event sum for the double channel
  with `Delta E_chem = 31.675 eV`; stated partition and its range; the
  47-48 eV appearance region investigated before 51.4 eV is called a
  measured onset; excited-fragment recipients for the neutral and
  dissociative channels assigned once; stopping lengths for fragments and
  for primary and secondary electrons.
- A6 unchanged.

### 4.4 Phase 3, thermodynamics: B1, B2, B4

**B1** unchanged from rev 2 (single H2 state model with completeness check,
standard-state functions, ortho/para policy by argument, population
consistency across EOS, chemistry, degradation and emission, trace-species
runtime condition).

**B2, H3+ cooling**, with the data selection as a named prerequisite (T4):

- **Emission model, one compatible choice.** Miller et al. (2013) Table 4
  gives unscaled and scaled LTE emission (3.4119e-18 and 4.7688e-18 W
  molecule^-1 sr^-1 at 5000 K); Table 6 gives non-LTE emission and factors
  whose high-density emission at 5000 K is 4.7587e-18, i.e. 0.998 of the
  scaled and 1.395 of the unscaled value. Therefore: either (a) the scaled
  LTE emission with Table 6 departures, reproducing Table 6's absolute
  emission within a stated data-consistency tolerance for rounded entries;
  or (b) the unscaled LTE emission multiplied by the Table 6 factors,
  described and validated as a different absolute-emission model. Today's
  code is (b) without saying so. The choice is recorded in the module header
  and in D0 item 10's validity list. Rounded table entries are not forced to
  exact identities.
- **Temperature coverage.** Table 4 starts at 500 K; Table 6 adds a 300 K
  anchor; the present routine is evaluated down to 30 K and has its largest
  join defect at 300 K. B2 states the supported range, the data or
  state/transition calculation supporting each interval, and the physical
  extrapolation or restriction outside it; no validated low-temperature
  extension is advertised without data or a derivation.
- **Continuity and density**, unchanged: continuity by construction with
  derivative continuity where the solver needs it; fit accuracy against the
  sampled data with the inter-sample uncertainty reported separately;
  low-density asymptote with exponent 1 and stated normalization;
  high-density end using the table value at 1e14 (0.9985 at 5000 K) with a
  smooth approach to 1; radiative and formation pumping as separate ledger
  recipients not forced to zero with n(H2). The local value of the present
  jump (0.9985 to 1) is measured; its effect on an atmosphere is not.
- The Koskinen et al. (2009) table as a keyed alternative remains Phase 6.

**B4** unchanged.

### 4.5 Phase 4: the spatial operator

Unchanged from rev 2 (D1 face-flux species advection with thermal enthalpy
and chemical fluxes and the chosen charged closure; D3 constraint count for
the full independent system; D2 discrete hydrostatic residual first, then the
base cycle; validity diagnostics per cell), with the diffusion-only column
test closing mass, elements, the specified current constraint and thermal
plus chemical energy.

### 4.6 Phase 5: the steady solver

Unchanged from rev 2 (C1 species rows on the independent space of D0 item 6,
unsplit stationary residual, refusals lifted one at a time, comparison with
a time-converged reference and an option-off isolation check against the
immediately preceding build; C2 relaxation-time elimination; C3 composition
convergence gate and reload test; C4 finite-difference Jacobian action
verified over `eps`, inner closure verified deterministically, preconditioner
omissions stated, `wasp_full_newton` base energy row).

## 5. Phase 6

Unchanged from rev 2: Koskinen term by term (the local paper suffices to
begin the process comparison; the spectrum file limits exact reproduction
only); lower/upper matching; `_adv` consistency; band photon conservation;
CO kinetics replacing the ceiling; molecular heat capacities; net-exchange
infrared form; band optical depth.

## 6. Ownership of every closure

As in rev 2 section 6, with these rows added or changed:

| Closure | Phase | Gate |
|---|---|---|
| ROE caller contract and vacuum flux | 1 (A2) | final flux and one conservative update admissible, with and without repair |
| Positivity-repair reporting: invocations, repaired interfaces, accepted repairs, direct LLF use | 1 (A1) | four quantities printed, zero included |
| Carrier failure recovery and rollback | 1 (A3 b) or 3 (A4+B3), whichever first | rejected species never written back; subdivided reference matched; interval and background restored |
| CO activation record | 1 (accompanies any ceiling use) | cumulative indicator, removed amount, cells and times, in final status and `EXHALE_resolved.out` |
| Two-cell radiation exchange | 3 | material sum zero |
| Independent species space for the stationary system | 2 (D0), 5 (C1) | no redundant rows; constrained Jacobian on that space |
| H3+ emission model choice and temperature coverage | 3 (B2 prerequisite) | stated model, stated range with data per interval |
| `q_D` admissibility across models | 1 (A8) | `0 <= q_D <= 1` or a diagnostic, never a silent clip |

## 7. Parallelism and effort

Delegable now, in parallel: A1, A3 increment (a), A8, the test construction
of 3.3, the baseline archive of 3.2, the projection-energy diagnostic, B4;
A2 after the interface decision of section 8. A3 increment (b) and the CO
model are not in this list. Everything from D0 onward is sequential at the
interface level. No calendar estimate is given.

## 8. Decisions for the user

1. Golden refresh once at the end of Phase 1 (3.2); the scratch baseline is
   the control until then; unknown provenance marked unknown. (Decided
   2026-09-05.)
2. Energy convention: `U(3)` kinetic plus thermal with the external-exchange
   identity and the local ownership rule of D0 item 3 (this plan), or a
   formation-inclusive total energy.
3. A2 interface: `speed_estimate_ROE` gains a status and returns edge speeds
   in the vacuum case; `Num_flux` gains a vacuum branch (candidate: HLLE with
   edge-speed bounds). The concrete signature is presented before coding.
   **Decided 2026-09-05: the default above. The signature and caller
   contract of `docs/a2_roe_interface.md` were approved by the user the same
   day; implementation proceeds as written there.**
4. A3 retry design: carrier-substep retry inside the transport step, or an
   outer attempted-step controller re-entering at the hydro stage.
5. Charged transport closure (D0 item 5): derived zero-current closure, or
   common-velocity approximation over an explicit domain.
6. CO constraint (D0 item 9): derive and validate the one-sided model with a
   timescale/domain argument, or exclude ceiling-active runs; in either case
   ceiling-active runs before that decision are diagnostic experiments.
7. B2 emission model: scaled LTE with Table 6 departures, or unscaled LTE
   with the factors as a separate model; supported temperature range.
8. D0 review and acceptance before Phase 3.
9. Oxygen ledger incomplete at the end of Phase 3: label and exclude, or
   hold.
10. Atomic-line in-hydro ionization transport where C2's criterion requires
    it: the restricted experiment is selectable and is reported as an
    approximation.
11. Request to Koskinen et al. for the spectrum file and the heating per
    ionization.
12. Instruction-gated items stay gated: paper and poster re-convergence, the
    LHS 1140 b write-up, the public switch of the repository.

## 9. What this plan does not claim

- Local values already measured by the retained probes: the ROE NaN and the
  negative thermal energy of one isolated update, the channel fractions at
  80 eV, the H3+ factor jump at 1e14, the Table 6 to Table 4 ratios. What
  remains unmeasured is the effect of each on an atmospheric solution, the
  base cycle and Mdot; those are measured at their gates.
- The isolated ROE update does not show that a full EXHALE step accepts a
  negative-energy state; the production loop has repair and step reduction
  that the probe deliberately bypasses. Test (b) of 3.3 measures both.
- Preservation of the converged atomic solutions under Phases 3-5 is a
  hypothesis.
- The Koskinen heating factor is attributed to neither code.
- Existence and stability of a molecular steady state with transported H+
  is what Phase 5 tests.
- Acceptance of a review finding in this plan does not repair the source;
  the confirmed defects (LLF speed, ROE branches and flux, carrier reference
  scales, Chung inversion, H3+ boundaries, fixed-T projection, H2 double
  channel energy, two-iteration energy solve) remain in HEAD `35d9dd5`.

## 10. Independent source audit of 2026-09-05: defects the three reviews did not find

Added after the rev 3 text above, at the user's request. Five auditors (Opus
workers) each read one module group of `src/` completely, with a brief that
listed every defect already recorded in sections 1-9 and in the three reviews
so that only new items would be reported; the brief, the five reports and the
verification probe are kept in `docs/audit_20260905/source_audit/`. Every P0
and P1 item below was then re-read in the source by the advisor, and the three
with the largest reach were reproduced independently (`threshold_bin_and_
sublyman_probe.py`; the heat-column comparison of 10.2 item 5 on the shipped
`mol_base_handoff` outputs). P2 items were spot-checked at their quoted lines.
No source was changed; no atmosphere was run.

### 10.1 New P0 defects (wrong answer in a default path)

**Landed 2026-09-05 (batch 2c, `Update_EXHALE_stage2.md` section 6):** items 1 (`dr_j`), 2 (thresholds on bin edges, extended to every active threshold and the H2 channel edges, and on 2026-09-06 to the loaded-SED grid as well, D0 item C5) and 3 (decision 13: Planck type in 2a, the Balmer continuum on the run's type in 2c).

1. **The cell width `dr_j(j)` is the width of cell j+1.**
   `src/modules/init/define_grid.f90` 152 and 185:
   `dr_j(2-Ng:N+Ng-1) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)`. With
   `r_edg(j) = r_{j+1/2}` (the stated convention, and the one `RK_rhs` uses:
   `rp = r_edg(j)`, `rm = r_edg(j-1)`) this is `r_edg(j+1) - r_edg(j)`, the
   width of the next cell. On the Mixed grid `dr_j` is exact in the uniform
   base and too large by the stretch factor (0.7% on `wasp_full`, 1.2% on
   `mol_base_handoff`, 1.5% on an r_max = 10 grid) everywhere else. Consumers
   that change the answer: the rectangle-rule column densities
   (`utilities.f90` 461-529, every photoionization optical depth), the PLM
   gravity source `Source.f90` 34 (the hydrostatic balance itself), the WENO3
   pressure gradient `RK_rhs.f90` 96 (the weight of pressure plus gravity
   against the advective flux divergence), `eval_dt`, the escape-probability
   path lengths, the residual volume weights, and `dr_cm` in the output.
   Inherited verbatim from ATES v2.0. Fix: shift the left-hand section by one
   (`dr_j(3-Ng:N+Ng) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)`). Every golden
   moves. Gate: `dr_j(j) = r_edg(j) - r_edg(j-1)` asserted for every cell in
   the derived grid test; the cell width used in `RK_rhs` and the one used in
   the column routine are the same object.

2. **The photon bin centered on a threshold charges the cross section over
   the half of the bin that lies below the threshold.**
   `set_energy_vectors.f90` 78-121: the grid points sit on 13.6, 24.6 and 54.4
   eV and the widths are central differences, so the 13.6 eV bin is 0.40-0.43
   eV wide when the sub-Lyman grid is present (20 logarithmic points from 4.8
   or 5.14 eV, last spacing 0.65-0.69 eV against 0.16 eV above) while only
   0.08 eV of it is physical. Reproduced today for the default power-law
   spectrum (PLind = -1, [13.6, 123.98, 1240] eV): P_HI code/exact = 1.0951
   (floor 4.80 eV), 1.0891 (floor 5.139 eV), 1.0003 (floor 13.6 eV); P_HeI
   1.0157 and P_HeII 1.0298 in all three (24.6 and 54.4 eV are grid points
   whatever the floor; measured by the Phase 0 test
   `src/tests/grid_and_gates/photon_grid_quadrature.f90` against the
   production `set_energy_vectors`); the H I heating is untouched (0.9998)
   because the photoelectron share vanishes at the threshold. The same test
   records that the X-ray block ends at `e_mid (e_top/e_mid)^((num_X-1)/num_X)`
   = 1184.19 eV, not at `e_top` = 1240 eV over which `J_inc` normalizes; the
   omitted band carries 1.9e-10 of P_HI, so it is a grid-definition defect
   without a rate consequence, to be closed with the bin edges. The bias persists through
   the ionized wind and dies only where the 13.6 eV bin is itself optically
   thick (code/exact 1.072 at N_HI = 1e17, 1.002 at 1e18). Reaches every case
   with `Include He23S? True` or a low-IP metal in `metals.inp`, i.e. the whole
   regression matrix. Fix: put thresholds on bin edges, or weight the
   straddling bin by its physical fraction. Gate: a quadrature test of P_HI,
   P_HeI, P_HeII against a fine reference integral, tolerance 1e-3, with and
   without the sub-Lyman grid.

3. **The 4.8-13.6 eV field is the EUV power law extrapolated two decades
   below its normalization band; it is the whole radiation the He 2^3S
   metastable and the low-IP metals see.** `J_inc.f90` 45-51 normalizes the
   power law on [e_low, e_mid]; `set_energy_vectors.f90` 78-83 and 258-261
   evaluate it down to `e_sub_low`. The same input states `Stellar Teff` and
   `Stellar radius`, and `excited_hydrogen.f90` 132-133 already integrates
   `pi B_nu(T_eff) (R_star/a)^2` over 3.4-13.6 eV for the Balmer continuum,
   so one run describes one band two ways. Reproduced today for `wasp_full`
   (T_eff 6459 K, R 1.458 Rsun, a 0.02544 AU, LEUV 30.42): the photospheric
   field exceeds the extrapolated power law by 1640x at 4.8 eV, 464x at 6 eV,
   40x at 8 eV, 2.7x at 10 eV, and falls below it above 12 eV; the 2^3S
   photoionization rate over 4.78-13.6 eV with the code's own cross section
   and Rate/2 dilution is 0.091 s^-1 (code) against 45.9 s^-1 (Planck), a
   factor 506. The Planck field itself over-estimates a real F star's NUV
   (line blanketing, Balmer jump), so the true correction is smaller than
   500 but is not small. The integrated grid flux is also 1.47x the nominal
   J_XUV the setup report prints. Reaches: every power-law run with the
   metastable on: the `HD209458b/`, `HD189733b/` and `WASP-121b/` planet
   folders, `examples/tutorial`, `benchmarks/hd209`, and all regression cases.
   The He 10830 results of those three planets were produced with this field
   (WASP-52b used a loaded SED and is not affected, provided its file reaches
   4.8 eV; not checked). The extrapolation is recorded in
   `docs/lower_atmosphere_coupling.md` 626-631 for the Lyman-Werner band; its
   consequence for He 2^3S is not recorded anywhere. Disposition (decision
   13, section 10.5): one spectrum type builds every band; for a power-law run
   the sub-13.6 eV field IS the power law by definition, a `Planck` type is
   added for the photospheric alternative, and the Balmer continuum follows
   the run's type too. Gate: one spectrum type per run in every band,
   including the Balmer continuum, and the setup report stating the
   sub-13.6 eV source and the integrated grid flux.

### 10.2 New P1 defects

**Landed 2026-09-05 (batch 2c, `Update_EXHALE_stage2.md` section 6):** items 1 (momentum kick), 2 (one state per output file; the heat column / breakdown mismatch itself stays for Phase 3), 3 (restart grid guard), 4 (profile base level), 5 (signed `du`), 6 (C I Voronov row), 9 and 11 (helium mass, Jupiter radius, one definition each with the thresholds and the other constants). Items 7, 8, 10, 12 landed in 2a (section 4).

1. **Momentum kick.** `EXHALE_main.f90` 1333: `u(2,:) = u(2,:) + 1.0e-16` is
   applied to the state, not to a copy, before the `dtu` diagnostic; `u_old =
   u` at 999 then carries it, so every cell and ghost gains 1e-16 code units
   of momentum density per step with no mass or energy change, outside the
   operator split and invisible to the steady residual. Measured against the
   shipped goldens: 5e-14 of the smallest momentum on `wasp_full`, 5e-11 on
   `mol_base_handoff`, 8.5e-7 on `lower_profile`. Inherited from ATES. Fix:
   floor the denominator of the diagnostic. Every golden moves at round-off
   level.

2. **The final output claims the secondary-ionization coupling was on while
   its heat and cool columns were produced with it off.** `EXHALE_main.f90`
   1784-1787 forces `sec_ion_active = .true.` after the loop without a
   re-sweep; `write_coupling_state_header` (`utilities.f90` 183) then writes
   `sec_ion=T`, and the breakdown files and the post-process run with the
   coupling armed. Measured on the shipped `mol_base_handoff` outputs: the
   `heat` column of `Hydro_ioniz.txt` and `heat_total` of
   `Heating_breakdown.txt` differ by a median of 72% and up to 82% (ratio
   0.26 in the layer, 1.0 in the wind). Reaches every run that leaves the loop
   without the staged flip (all `mol_*` and `lower_profile` snapshots), and a
   restart from such a file adopts a coupling the state was never relaxed
   under. Fix: one state per file (re-sweep once, or record the coupling the
   state was relaxed under).

3. **A restart overwrites the cell centers with the file's radii and keeps the
   run's faces, widths and window indices.** `load_IC.f90` 751 (`r(j) =
   vals(1)` in `scatter_row`) and 214; `r_edg`, `dr_j`, `j_min`, `j_flux` are
   not rebuilt. The only guard is the row count. Fix: leave `r` alone, or
   refuse a mismatch.

4. **With a lower-atmosphere profile the base level is stated twice and
   never compared.** `input_read.f90` 1652 skips the base-level consistency
   check when a profile is in use, `apply_lower_atmosphere_profile` (2326)
   adopts `p_match_bar` as `p_base_bar`, and `n0` still comes from the density
   key. The shipped `lower_profile` case ran with the profile at 1.00e-6 bar
   (`run.log` line 10) and the base at 2.0035e-5 bar (`EXHALE_setup.out` line
   38), a factor 20, so its T, r, q_H2 and elemental ratios were taken from a
   level the wind does not start at. Fix: `n0` from `p_match_bar`, or refuse
   the density key beside a profile.

5. **The `du` stop is the spread of |rho v r^2|, not of rho v r^2.**
   `EXHALE_main.f90` 1309-1316 takes the absolute value before the max and
   min and normalizes by the minimum; `flux_spread` twenty lines below
   (1368-1369) is the signed spread over the mean on the same window. On the
   `lower_profile` golden, whose window changes sign twice, the coded form
   reads 227x the signed spread (conservative there); a window with +F and -F
   halves would read zero. Fix: one functional, the signed one, with the
   arming and the hand-off reading it.

6. **The C I electron-impact ionization row of Voronov (1997) is transcribed
   with P and X shifted.** `Cool_coeff.f90` 3004-3007: `P = 0.193, X = 0.25`
   with the `(1 + P sqrt(U))` factor; the published row is P = 0, X = 0.193
   (as `p-winds/p_winds/carbon.py` 161 has it). The rate is 2.56x too high at
   2000 K, 1.69x at 1e4 K, 1.29x at 5e4 K. Every other Voronov row in the file
   is correct. Reaches every metals-on case (`wasp_full`, `wasp_he23off`,
   `mol_metals`, `lower_profile`).

7. **Opacity model P multiplies the column densities only.** `opa_pf` enters
   `calc_column_dens*` and nothing else, so the beam loses f(p) times more
   photons than ionize anything or heat the gas. Opt-in (`opacity.inp`,
   model P); no regression case ships one.
   **Landed 2026-09-05 (batch 2a, `Update_EXHALE_stage2.md` section 4).** The factor now multiplies the shared flux
   weight of every absorption integrand, so the beam's loss equals the local
   absorption; a cell with f = 10 closes energy and H I ionizations to 3e-8.

8. **The Lyman-alpha photolysis band carries the stellar flux with no atomic
   hydrogen attenuation.** `util_ion_eq.f90` 179-196 and 355-374,
   `water_photolysis.f90` 403-411: band B2 is the H I resonance line, its
   optical depth is `sigma_H2O N(H2O) + sigma_OH N(OH)` only, and `tr_out` is
   forced to 1 for every band but LW, while `lya_rt.f90` 213-226 already
   computes the penetrating stellar Lyman-alpha beam through the H I column.
   The module header (251-252) says atomic H does not absorb these bands,
   which is true of the continuum and false of the line. Opt-in (`Oxygen
   chemistry: True` with a Lyman-alpha flux; `examples/18_oxygen_chemistry`,
   where B2 is 6% of the unattenuated photolysis; unbounded for a
   Lyman-alpha-dominated star). No regression case has oxygen chemistry.
   **Landed 2026-09-05 (batch 2a, `Update_EXHALE_stage2.md` section 4).** Band B2 carries
   `lya_stellar_beam_transmission(T, tau)`, the expression `lya_rt.f90` already
   held inline. The `oxygen_chemistry` case moved as expected and its scratch
   baseline was re-snapshotted.

9. **The Wind-AE bridge uses a second Jupiter radius.** `wae_exhale_input.f90`
   14 has `RJ = 7.1492d9`; `parameters.f90` 691 has `6.9911d9`; both multiply
   the same `Planet radius [R_J]` line, so the Wind-AE seed is built for a
   planet 2.26% larger than the one EXHALE relaxes. MJ and MSUN differ by
   0.011% and 0.029%. Opt-in (`IC mode: windae`, `wind_ae_ic.x`).

10. **The same bridge converts He/H to mass fractions with 4 m_H and hands
    the solver CODATA masses** (`wae_exhale_input.f90` 71-77), so Wind-AE
    solves at 1.00717x the requested He/H. Same reach as item 9.
    **Landed 2026-09-05 (batch 2a, `Update_EXHALE_stage2.md` section 4).** The mass fractions
    are built from the same `m_He/m_H` the solver is handed, and the
    reconstructed He/H equals the input to 1e-12. The `RJ`/`MJ`/`MSUN`
    half of the pair is decision 14 and is not in this batch.

11. **Helium weighs 4.0 hydrogen atoms in the species table while the mass
    unit is the hydrogen atom.** `species_table.f90` 104-107 with
    `parameters.f90` 675 (`mu = 1.67353284d-24`, the H atom): the correct
    ratio is 3.9715, so helium is 0.72% too heavy and the mass per H nucleus
    at He/H = 0.0793 is 0.17% high, propagating into `rho_bc`, `v0`, `p0`,
    `q0` and Mdot. The metals in the same table carry their true atomic
    weights. The comment "He mass 4.0, not 4.0026" names the wrong correct
    value (4.0026 is in u, the table is in m_H). Every golden moves by ~0.17%,
    above the 0.1% rule.

12. **The constrained-equilibrium continuation ignores the transported-proton
    constraint.** `constrained_chemical_equilibrium.f90` 807-863 reads
    `x_h2_fixed` and `x_ox_fixed` and never `x_hp_fixed` (read only by the
    seed at 1464 and by the dump I/O), so with `Ionization transport: True`
    it solves the local H+ root and returns a candidate the acceptance test
    (`ion_system_HeH_mol_metals`, which does pin row 1) then rejects; each
    such cell spends up to 40 `hybrd` solves before falling to class 4.
    Opt-in. **Landed 2026-09-05 (batch 2a, `Update_EXHALE_stage2.md` section 4).** `hp_is_fixed` pins
    `species_fixed(is_HII)` like H2 and the oxygen carriers. The three `hp_*`
    cases never enter the continuation, so no run exercises it yet.

### 10.3 New P2 items (hygiene, documentation, dead code, test coverage)

| Item | Location | Note |
|---|---|---|
| Conduction stage floors T at 0.01 T0 silently | `viscous_conduction.f90` 527-532 | the explicit energy update counts the same floor; this one does not, and a floored cell is not a zero of the steady residual. Landed 2026-09-06 (B2b: a floored cell is a failure of the stage, counted by `conduction_floor_cell_hits`; `Update_EXHALE_stage2.md` section 7). |
| `minloc(..., dim=1)` used as a declared subscript | `set_IC.f90` 116-118 | index off by Ng; out of bounds near the top |
| `dF(:,1-Ng)`, `S(:,1-Ng)`, `WL(:,1-Ng)` read before written | `RK_rhs.f90` 92, `Apply_BC.f90` 231 | overwritten before use; traps under `-finit-real=snan` |
| `refresh_row_terms` drops `Smom`/`Sene` and accumulates the energy scale with `max` | `steady_residual.f90` 298-303 | systematic in the self-consistency loop |
| `it_sc` printed one too high on non-convergence | `steady_newton.f90` 2640-2649 | DO variable after normal completion |
| `j_min` has no lower clamp; the search starts in the ghosts | `define_grid.f90` 193-197 | `Escape radius <= 1` puts ghosts in every window |
| `okres` written twice, never read in the ray test | `steady_newton.f90` 1848-1860 | inadmissible probes enter `fslope` |
| `build_banded_jac` keeps the flat `max(|Y|,1)` step | `steady_newton.f90` 744 | diagnostic caller only |
| Reversal branch of the base boundary wider than its comment; ghost stratification on the base composition there | `base_boundary.f90` 346-349, 446 | `w_rev = 1` for all `M_i <= -1e-6` |
| Mixed grid and `Rec_BC` assume Ng = 2 without saying so | `define_grid.f90` 77-83, `Apply_BC.f90` 263 | |
| HLLC returns the upwind face pressure where LLF/ROE return the mean | `Num_Fluxes.f90` 117-132 | SUSPECTED effect on the WENO3 momentum operator; `p*` is available |
| `report_unmet_steady_gates` never called | `EXHALE_main.f90` 2606 | dead |
| Mdot index `N - 20` unguarded for N < 21; `log10` of `v <= 0` unguarded | `EXHALE_main.f90` 1943-1946 | |
| `get_word` returns an unallocated result for a missing word; `=` and TAB accepted as label separators but not as word separators | `input_read.f90` 2452-2551, 2361-2380 | both compilers return an empty string today |
| `close(unit=1)` on a unit never opened | `input_read.f90` 1168 | dead |
| Post-process `xion` built once per pass, reused after the advection solve rewrote the densities; the comment says the opposite | `post_process_adv.f90` 346-355, 820 | |
| The regression matrix compares no `_adv` output | `run_check.sh` 153 | `post_process_adv.f90` has no golden |
| A failed `make` inside `run_check.sh` is swallowed by `&&` under `set -e`; the summary greps every case directory | `run_check.sh` 264, 293 | the `make -q` guard still stops a stale binary |
| `EXHALE_transit.py` metals guard compares a column count that is never below 26; `max_rows=1` warns on numpy 1.26 | `EXHALE_transit.py` 872-876 | metals-off branch is dead |
| Triaxial branch keeps the superseded multiplet Balmer treatment and mixes two annulus radii | `EXHALE_transit.py` 1190-1253 | branch is dead by default |
| Rotational Doppler shift applied with the wrong sign | `EXHALE_transit.py` 774 | cancels by the azimuthal symmetry of the spherical model |
| `exhale_io.read_input` reads `input.inp` where `EXHALE_resolved.out` holds the resolved radius; docstrings wrong about `LX` and the row header | `examples/exhale_io.py` 346-441 | latent for handoff runs |
| `collisional_validity.py`: electron-ion coupling carries an extra ion-fraction factor; sound speed uses 5/3 where the run used the caloric gamma; cited line numbers stale | `src/utils/collisional_validity.py` 156-163, 290, 588-592 | |
| FUV and infrared column ledgers sum the ghost cells | `write_output.f90` 653, 727 |. Landed 2026-09-06 (LEDGER-FUV, section 7). |
| Transit output headers do not state the wavelength frame; the line list mixes air and vacuum | `exhale_transit_lib.py` 56-99 | He 10830 air/vacuum differs by 2.98 A |
| `eval_cool` charges the He I ground-state collisional coefficients to the total He I, triplet included, then adds the 2^3S terms | `util_ion_eq.f90` 1516, 1568 | <= 2.2e-4 of the He I channel on `wasp_full`. Landed 2026-09-05 (batch 2b, item 2b-HEI) |
| He II recombination cooling at alpha_B while the triplet path recombines at rec_11S + rec_23S (and the coupling replaces alpha_B) | `ionization_equilibrium.f90` 935-938, `Cool_coeff.f90` 3409 | 23-35% low on that channel. Landed 2026-09-05 (batch 2b, item 2b-HEI) |
| Chemical-equilibrium retry seed applies (1 - x_H2) twice | `ionization_equilibrium.f90` 3143, 3156 | seed only |
| Pure-hydrogen secondary-ionization path gives helium's share to hydrogen (+14%) | `electron_energy_degradation.f90` 663-670 | design choice, size not stated |
| The photon field of a cell is the field at its inner face (whole-cell column) | `utilities.f90` 461-483 | ~1.5% at the front on `wasp_full`; ATES heritage. Landed 2026-09-06 (B4-6/6b/6c cell-mean attenuation, section 7). |
| Two spectra describe 3.4-13.6 eV in one run; the Balmer continuum is unattenuated | `excited_hydrogen.f90` 132-133 | the source (Christie 2013) makes the same approximation for n = 2 columns only |
| Dead parameters `f_sing_HeI`, `Ee_sing_HeI`; `h2_channel_threshold` unused; `h2_channel_stoichiometry` header claims every source term derives from it while `System_HeH_mol.f90` writes them by hand; O I cutoff exponent quoted wrong in a header; `iscool` list comment names 8 of 12 ions; `gbar_ff` names the wrong Rydberg constant; legacy `use_2lev_cool` omits Fe I; `n_set` under-counts in `metals_input_read` | various, `group_B_report.md` "Dead and stale items" | |
| `equilibrium_constant_conc` overflows to Inf below ~10.4 K and the oxygen seed becomes NaN | `oxygen_rates.f90` 770-791, 929-959 | reachability not shown |
| `k_lw` evaluated at the inner-face column while the H2O/OH rates of the same beam use the cell mean | `util_ion_eq.f90` 337-346 | one-signed low at the H2 front; moves `mol_lyman_werner`. Landed 2026-09-05 (batch 2b, item 2b-KLW) |
| R5 and R16 dissociative recombinations deposit 10.9 and 11.8 eV as heat although H(n = 2) at 10.2 eV is an open exit | `molecular_reaction_heat.f90` 224-250 | upper bound not stated |
| `oxygen_rates.f90` header says it is not in the build; it is (`Makefile` 133) | `oxygen_rates.f90` 5-9 | |
| "FUV shell average" key described in two comments does not exist | `parameters.f90` 1208-1214, `lyman_werner.f90` 153-157 | |
| `Tcode` threaded through four carrier routines and never used; the temperature used is `bg_cell%T_K` | `diffusive_photochemistry.f90` 570-2289 | trap for a caller that passes a T the sweep has not seen |
| `adv_corr` allocated before the `weno_mode = 2` early return | `diffusive_photochemistry.f90` 1063-1080 | one reordering from an uninitialized read |
| Dead/duplicated constants and an unread `wp_ready` in `water_photolysis.f90`; public module variable `i` in `molecular_infrared_data` | 287-388; 24-53 | |
| `carrier_mass_amu` hardcodes `bsp_mass(7/11/12/2/13)` under a comment saying it reads the table | `diffusive_photochemistry.f90` 1238-1254 | |
| Harmonic-mean slope does not cancel the donor-cell truncation exactly on a quadratic, as claimed | `diffusive_photochemistry.f90` 230-235 | second order either way |
| FUV band edges leave 1-A gaps at 1201, 1230, 1450 A and no 1202-1230 A continuum | `oxygen_rates.f90` 382-386 | user-facing flux definition |
| Fluorescent-trapping ratio in the self-shielding table is 1.17-1.32 at the optically thin end instead of 1 | `h2_self_shielding_table.f90` 93-104, generator | SUSPECTED: slab-versus-sphere geometry of the CLOUDY runs |
| `element_census` computes the charge density and never checks it; mass closure never fails; metals reservoir test ignores `he_metal_diffusion` | `element_census.f90` 249-289, 375-390 | assert runs only |
| Threadprivate initializer claim ("reaches the master thread only") is false for gfortran 16; `System_HeH_metals` header says the sweep is serial | `ion_cell_state.f90` 64-67, `System_HeH_mol.f90` 170-174, `System_HeH_metals.f90` 16-19 | documentation |
| `pi` truncated to 11 digits against the block's own policy | `parameters.f90` 665 | 1.3e-11 |
| H I / He I ionization energies duplicated as single-precision literals in `T_equation.f90` 118-120 as well as `eval_cool` | | extends a known item |
| `sigma_tab` column comment lists 15 of 17 ions (Fe I, Fe II missing) | `parameters.f90` 1062-1066 | |
| One `tol` passed to `hybrd1` (xtol) and to `newton_dense` (a residual 2-norm in cm^-3 s^-1) | `newton_solver.f90` 36-137, `ionization_equilibrium.f90` 613 | measured: the residual branch is unreachable at base density |
| `calc_mmw` counts nuclei as particles (wrong in molecular gas; atomic-only caller today) | `utilities.f90` 539-574 | |
| `he_metal_diffusion` and `he_alphaT` never exercised by `diffusion_tests`; `element_census_tests` do not cover the verify/reservoir routines | `src/tests/` | coverage |
| The four A2/E1 test programs assert nothing (all `write`); no oxygen-chemistry case in the matrix | `src/tests/a2_*`, `e1_h2`, `run_check.sh` 187 | the whole oxygen path is outside the regression |
| `e_th_HeTR = 4.80` eV against the cross-section threshold 4.78 eV in `sigma_HeI23S` and the level's 4.768 eV | `parameters.f90` 725, `cross_sec.f90` 96 | grid floor and cross-section onset disagree by 0.02 eV |
| `WASP-52b/input.inp` line 12 points at `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE/WASP-52b/eps_eri_sed_fxuv1p0.txt`, a directory that no longer exists since the 2026-07-14 relocation; the planet folder's own run cannot start (`benchmarks/wasp52` uses a relative path and works) | `WASP-52b/input.inp` 12 | Fixed 2026-09-06 (C3, section 7): the path is now relative to the folder. |
| The top-level Koskinen benchmark (`benchmarks/koskinen2022_model_a/input.inp` 5 and its README) runs at 0.048 au while Model A and the SED file beside it are at 0.05 au; the three `matched_*` runs use 0.0500 | `benchmarks/koskinen2022_model_a/input.inp` 5 | geometry at 0.048, radiation at 0.05 (8.5% in flux) in that older run |
| `HD209458b/input.inp` 6 has 0.0480 au against the literature 0.047 (Salz 2016 Table 2) / 0.04707 (p-winds), 4% in flux | `HD209458b/input.inp` 6 | settle at the 2d switch. Settled 2026-09-06: 0.04707 AU, the Oklopcic and Hirata (2018) value (C3, section 7). |
| `LHS1140b/sed/README.md` last section describes integral-matched products with `A_eff` headers; the shipped file header and `make_sed.py` 58-63, 95-98 carry `A = 0.59` (catalog X-ray ratio) | `LHS1140b/sed/README.md` | stale paragraph |
| `inputdata/scaled_solar_hd189.sed`: referenced nowhere, no build script, entered in a bulk commit (586520e); its 1-912 A flux matches `HD189733b/input.inp`'s LX+LEUV to 0.03% but its hardness F(1-100)/F(1-912) = 0.089 contradicts the 0.355 those luminosities imply; parent spectrum unidentified; ends at 1899.5 A | `inputdata/scaled_solar_hd189.sed` | orphan; deletion is a user decision (`inputdata/sed/README.md` section 5) |
| `e_th_HI = 13.6` eV serves both as the photon-grid threshold and, through the heating share, as the chemical energy of an H+ + e pair; CODATA is 13.598434599 eV | `parameters.f90` 708 | 0.003 eV; a decision for A5 when the ledger separates thresholds from chemical energies |
| The LLF header claims the Perthame and Shu positivity property for a viscosity the code does not implement (`max(abs(vL+aL),abs(vR+aR))` instead of `abs(v)+a`) | `Num_Fluxes.f90` 246-258 vs 287 | documentation defect attached to A1. Landed 2026-09-05 (batch 2b, A1). |

### 10.4 What changes in the plan

- **Phase 1 grows by the bounded corrections** 10.1 item 1 (`dr_j`), 10.1
  item 2 (threshold bins), 10.2 items 1, 3, 5, 6, 9, 10, 11, 12 and the P2
  index/guard items; each is algebraic with a local test and none needs D0.
  10.2 item 2 (one state per output file) and item 4 (one base level) are
  Phase 1 as well, because they are consistency defects of the record, not
  physics choices. Together with the known items these move every golden
  (items 1, 2 of 10.1 and items 1, 11 of 10.2 are in every case), so the
  Phase 1 golden refresh of section 3.2 happens once, after all of them, with
  each item's movement measured separately first on one atomic and one
  molecular case.
- **10.1 item 3 (the sub-Lyman field) becomes, under decision 13, a
  spectrum-type consistency item of Phase 1**: a `Spectrum type: Planck`
  option, the Balmer continuum following the run's type, and the setup-report
  lines; no band is filled from a type the input did not select. The He 10830
  results of HD 209458 b, HD 189733 b and WASP-121 b were computed under the
  power-law type and therefore with a power-law NUV; that is recorded as the
  approximation they carry, and a Planck or table rerun is part of the
  instruction-gated tail of section 8, not an autonomous action.
- **Phase 0 gains tests**: grid-width identity; photon quadrature against a
  fine reference with and without the sub-Lyman grid; restart onto a
  different `Grid type` refused; `du` functional against a sign-changing
  window; header/state consistency of every output file (the `# coupling`
  line against a re-sweep); base level stated once; the `_adv` files enter
  the regression comparison; an oxygen-chemistry case enters the matrix; the
  four A2/E1 programs acquire verdicts (already in 3.3).
- **Phase 3 (A5, A6, B2) gains** the R5/R16 product-state question, the
  He II recombination-cooling coefficient of the triplet path, the He I
  singlet density in `eval_cool`, and the Lyman-alpha photolysis attenuation
  (through the existing `lya_rt` beam).
- **Phase 4 (D1)** gains the cell-mean photon field (the inner-face column
  of `calc_column_dens`) and the `k_lw` cell mean, as one discretization for
  every absorber of one beam.
- **Phase 5 (C1)** must honor `x_hp_fixed` in the constrained continuation
  before the H+ row is added, or the continuation is dead for exactly the
  cells the row exists for.

### 10.5 Additional decisions for the user

13. **Sub-Lyman field** (10.1 item 3). **Decided 2026-09-05 (final
    statement, replacing the two earlier ones): one spectrum type builds every
    band.** `Spectrum type: Power-law` means the power law on the whole grid,
    the XUV and the part below 13.6 eV alike (the present construction, now
    the stated meaning of the key rather than an accident); a new `Spectrum
    type: Planck` means `pi B_nu(T_eff) (R_star/a)^2` on the whole grid from
    `Stellar Teff` and `Stellar radius`; `Spectrum type: Load` means the SED
    table on the whole grid, which must cover the grid floor. No band is
    built from a different type than the rest, so the Balmer continuum of
    `excited_hydrogen.f90`, which today integrates Planck in every run, is
    changed to follow the run's spectrum type as well. Consequence for 10.1
    item 3: it is no longer a defect the code repairs by inserting a
    photospheric field; it is the meaning of the power-law type, and the
    setup report states it (the sub-13.6 eV field is the same power law, and
    the integrated grid flux against the nominal J_XUV), so that a user who
    wants a photospheric NUV selects Planck or a table. The threshold-bin
    quadrature (10.1 item 2) stands independently.
    **Landed 2026-09-05 (batch 2a, item 2a-PLANCK, `Update_EXHALE_stage2.md`
    section 4):** the `Planck` type, which requires `Stellar Teff` and
    `Stellar radius` and fills the whole photon grid from
    `planck_stellar_flux_eV`, and the three setup-report lines (the type, the
    source of the band below 13.6 eV, the integrated grid flux against the
    nominal J_XUV) that every run now prints. **Remaining:** the Balmer
    continuum of `excited_hydrogen.f90` still integrates Planck whatever the
    type and is to follow the run's type in batch 2c; and `J_XUV` (the total
    the Wind-AE bridge is handed and the setup report calls the XUV flux),
    `F_inc` of the Ly-alpha pumping estimate, the five FUV band fluxes and the
    monochromatic `F_XUV(1)` are still built from `LX`/`LEUV` or from input
    keys rather than from the spectrum type, so the decision is enforced for
    the photon grid only.

14. **Mass and radius constants** (10.2 items 9, 11). **Decided
    2026-09-05:** `m_He/m_H = 3.9715` in the species table (the user's note
    reads 3.9175, taken as a transposition of the physical value 3.97152;
    every golden moves ~0.17%), and **`RJ = 7.1492e9` cm** (IAU 2015 nominal
    equatorial radius) for the code and the Wind-AE bridge alike. The change
    of `RJ` rescales `R0` of every run by +2.26%, so every golden, every
    planet-folder result and every benchmark moves; the same change makes
    `MJ` and `MSUN` the IAU 2015 nominal values the bridge already carries
    (1.8982e30, 1.98842e33 g; both below the 0.1% rule). **Follow-through
    requested by the user 2026-09-05 (item 2c-CONST):** the physical masses
    themselves live in `parameters.f90` (`m_He` = 4.002602 u, `m_e`, `m_p`,
    `amu`, CODATA 2018, with the derived `m_He_over_m_H = m_He/mu`), and
    `species_table.f90` derives `bsp_mass` from them instead of carrying the
    literal 3.9715; the metal atomic weights are converted with `amu/mu`.
    The Wind-AE modules keep their copies for the standalone build, and the
    Phase 1 test asserts they equal the definitions.
15. **`du` functional** (10.2 item 5). **Decided 2026-09-05: the signed
    spread over the mean** that `flux_spread` already computes; the stop, the
    JFNK arming and the secondary-ionization flip fire at different steps
    than today, and the movement is recorded.

16. **SED of every planet example** (decided 2026-09-05). Every planet
    calculation of record (the planet folders `HD209458b/`, `HD189733b/`,
    `WASP-52b/`, `WASP-121b/`, `LHS1140b/`; every `benchmarks/` case; the
    `examples/` runs that represent a planet) uses `Spectrum type: Load` with
    the SED the reference literature used, or assumed, for that planet. One
    file per planet under `inputdata/sed/`, with a README stating the paper,
    the proxy star, the scaling to the planet's orbit, the wavelength
    coverage and the unit conversion. The power-law type remains for the
    tutorial and for the technical regression snapshots, which are not planet
    results. Present state, checked today: `WASP-52b/` and `benchmarks/wasp52`
    already load an eps Eri proxy (5.9-3700 A, covers the 4.8 eV floor);
    `benchmarks/koskinen2022_model_a` loads the Ribas 2005 solar file
    (5-1995 A, which stops at 6.2 eV, so `sed_read` warns and the metastable
    would see no field below it); `HD209458b/`, `HD189733b/`, `WASP-121b/`
    and their benchmarks run the power law; `inputdata/scaled_solar_hd189.sed`
    (0.5-1900 A) exists and is referenced nowhere. Work: (a) an inventory,
    per planet, of the SED the reference paper used or assumed (Huang et al.
    2023 for WASP-121 b; the He 10830 papers for HD 209458 b and HD 189733 b;
    Koskinen et al. 2022 for the hot Uranus), read in the PDFs, with the file
    obtained or built and its coverage extended to the grid floor where the
    metastable is on; (b) the input files switched, as a configuration change
    recorded in each folder's README (Phase 1, 2d). The reruns themselves
    stay in the instruction-gated tail unless instructed otherwise; the
    switch is what makes a later rerun comparable to the literature.
    **Inventory done 2026-09-05 (Phase 0, W4; `inputdata/sed/README.md`
    with verbatim paper quotes, build scripts and a coverage test, all
    PASS):** BUILT `wasp52b_epseri_yan2022.txt` (MUSCLES v22 eps Eri const-res
    SED scaled by (3.212 pc/0.0272 au)^2, reproduces the existing file to
    round-off; 5.9-3700 A), `hot_uranus_solar_koskinen2022.txt` (the
    benchmark's Ribas 2005 band table 5-1995 A, byte-identical, continued to
    2995 A with the Gueymard composite the table was cut from above 1180 A;
    a stand-in for SOLAR2000 + Woods and Rottman with the right band fluxes,
    stated in its header), `hd209458b_solar_whi2008.txt` (Oklopcic and Hirata
    2018: WHI 2008 solar-minimum spectrum scaled to 0.04707 au; reproduces
    p-winds' file exactly), `lhs1140b_gj1132_cherubim2026.txt` (the existing
    GJ 1132 proxy, A = 0.59). NEEDS DATA: WASP-121 b (Huang 2023: solar SED
    rescaled to F_XUV = 1.6e6 below 1700 A plus an LLmodels F5V spectrum
    above; the NUV half is exactly the metastable band, so no partial file
    was written) and HD 189733 b (Salz 2016: CHIANTI plasma + 9.4 A Gaussian
    Ly-alpha + Woods and Rottman EUV shape + blackbody, normalized to their
    Table 3; the solar table is not in the workspace; Bourrier 2020's MOVES
    III spectrum in `VULCAN/atm/stellar_flux/` is a published alternative
    whose normalization convention is unresolved). Request to the user: the
    two missing spectra, or permission to build the HD 189733 b one from
    CHIANTI plus a WHI solar shape as a documented stand-in.

17. **SED that does not reach the triplet band** (decided 2026-09-05, final
    wording). When a loaded SED ends above the grid floor an active absorber
    needs (He 2^3S at 4.8 eV, i.e. 2583 A; or a low-IP metal at its
    threshold), the run STOPS. There is no key to continue. The stop message
    names the file, the missing band in eV and A, and the absorbers that
    need it, and offers the two remedies: state `Include He23S? False` (or
    remove the metal from `metals.inp`) in the input, or supply an SED that
    covers down to that wavelength. The present warning-only branch of
    `sed_read.f90` 113-117 is replaced by this stop. Checked today: the three
    Koskinen benchmark inputs run with `Include He23S? False`, so their
    1995 A file is not caught; the WASP-52 b file reaches 3700 A. Phase 1
    (2a): a reader behavior change that moves no golden, since every
    regression case is power-law.
    **Landed 2026-09-05 (batch 2a, item 2a-SED, `Update_EXHALE_stage2.md`
    section 4):** the run stops with `error stop 1` and no key to continue,
    and the message names the file, the band not covered in eV and A, each
    absorber whose threshold lies below the lowest photon in the file (He 2^3S
    and every active neutral metal, K I at 4.341 eV included) and the two
    remedies. Six assertions in
    `src/tests/grid_and_gates/sed_coverage_stop.sh`; no configuration of
    record in the tree trips it. **Remaining:** nothing of this decision. The
    symmetric case at the top of the grid, a table whose highest row is below
    `e_top_read`, is still silent and is outside the wording of this decision
    (noticed by the item, not fixed).

### 10.6 Scope of this audit

Read completely by the auditors: every `.f90` under `src/modules/` except
`wind_ae/` (only its EXHALE bridge was read), `src/EXHALE_main.f90`, the six
test programs, `EXHALE_transit.py`, `exhale_transit_lib.py`,
`examples/exhale_io.py`, `src/utils/collisional_validity.py`, `Makefile`,
`run_check.sh`. Not checked: the MINPACK internals, the digits of the large
tables (Fe II cooling, Gaunt factors, Samson and Haddad, Chung, Dalgarno,
CHIANTI polynomials, the self-shielding cubes beyond their headers' own
claims), the Verner 1996 metal coefficients and the Mao and Kaastra
coefficients against their published tables (done 2026-09-07, item
REF-METALS: 132 photoionization and 220 recombination parameters match,
one Mg II threshold gate and one Ca I recombination row corrected), the
Shomate coefficients against NIST, the offline table generators and CLOUDY
decks, `roche_recon.py`, `EXHALE_plots.py`, `eta_approx.py`, the VULCAN and
photochem drivers, and the effect of any item on a converged atmosphere,
which no auditor ran. The clean checks each auditor recorded (Roche
potential, gravitational work flux, viscous stress, caloric EOS identities,
characteristic base condition, Verner form, Badnell and Voronov rows other
than C I, Christie 2013 rates, Koskinen Table 1 R1-R23, the H2 equilibrium
constant, the lower column against Model A, the H2 line data, analytic
Jacobians against finite differences, OpenMP privacy, input defaults against
`docs/input_schema.md`, restart field coverage, profile reader units) are
listed in the five reports so that the negative results are on record too.
