# Third review of the EXHALE development plan

Review date: 2026-09-05 (KST).

Reviewed plan: [development_plan_20260905_rev2.md](development_plan_20260905_rev2.md).
Previous assessment: [development_plan_20260905_review2.md](development_plan_20260905_review2.md).
Source inspected: `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`.

## 1. Assessment

Revision 2 is a reasonable development roadmap and a suitable basis for writing D0, the governing-system document. The main physical contradictions identified in the previous review have been corrected in the plan. In particular, the external-energy convention now includes chemical storage correctly, the Chung detector observable is correctly distinguished from particle multiplicity, and thermal enthalpy transport is derived from the chosen energy variable. These issues should be considered addressed at the planning level rather than repeatedly reopened without new evidence.

The plan still needs several concrete amendments before individual implementation tasks are accepted. The most consequential concern the interface between the proposed ROE vacuum treatment and the actual flux routine, the scope of A3's failed-step retry, and the meaning of a run in which the CO ceiling has been active. H3+ cooling also needs a precise choice of reference emission and temperature coverage before B2 can be called a complete specification.

These remaining issues do not justify discarding the phase order or postponing all work. A1, construction of the missing tests, and the A8 detector inversion can proceed with targeted validation. The A3 scale correction can proceed independently of a larger retry implementation. D0 should now resolve the shared equations and validity conditions, using the acceptance checklist in Section 5 of this review.

Production source remains unchanged. This review adds the requested document and preserves new diagnostic programs in `docs/audit_20260905/`. No atmosphere, full regression matrix, or existing reference output was regenerated.

## 2. What was checked

The full revision 2 document was read. Source inspection followed the proposed corrections into their callers and downstream consumers:

- `speed_estimate_ROE` into `Num_flux`, including its entropy correction and the final numerical flux.
- The hydrodynamic `retry_step` loop, the later carrier call, `solve_carriers`, and `carrier_write_back`.
- CO-ceiling counters, their reset and reporting thresholds, and setup/final output paths.
- The H2 channel routine, its double-fraction definition, energy-grid initialization, photoionization/heating consumers, and the existing E1 diagnostics.
- Existing test build targets and the retained audit programs.

The local published PDFs were used for the relevant source checks: Chung et al. (1993), Yan et al. (1998), Miller et al. (2013), Koskinen et al. (2022), and Schulik and Booth (2023). The checks used their experimental methods, equations, tables, and numerical-method descriptions, rather than the plan's account of those papers. The ADS identities were already verified during the second review and were not queried again solely to repeat unchanged metadata.

A new program called the actual production ROE flux for two symmetric expansions. It then applied that flux to a simple homogeneous-Euler cell update. This follows a previously demonstrated auxiliary-routine defect into its consumer; it is not a repeat of an untouched planetary regression. Separate arithmetic checks evaluated revision 2's energy identity and the compatibility of the published H3+ normalization choices.

The Fortran diagnostic used GNU Fortran 16.2.0 with `-O0 -g -fcheck=all -fbacktrace -fopenmp`. The source and runner are [roe_flux_vacuum_probe.f90](audit_20260905/roe_flux_vacuum_probe.f90) and [run_review3_checks.sh](audit_20260905/run_review3_checks.sh). The arithmetic checks are in [review3_ledger_checks.py](audit_20260905/review3_ledger_checks.py). Their scope and limitations are stated below.

## 3. Previous findings: what is now resolved in the plan

| Previous finding | Revision 2 assessment |
|---|---|
| Q1: Chung detector inversion | Correctly adopted as A8. The formulas now invert an event signal and apply particle multiplicity afterward. |
| Q2: double subtraction of photoionization energy | Corrected. The source-step identity uses full external exchange with the chemical increment. |
| Q3: excited fragments and vertical threshold | Appropriately adopted. Unknown state branching is identified as a modeling question, and the vertical threshold is no longer described as the measured appearance threshold. |
| Q4: CO constraint | The complementarity equations, products, energy, and exclusion alternative are now assigned to D0. The temporary Phase 1 treatment still needs amendment under T3. |
| Q5: thermal enthalpy and charged transport | The energy-convention derivation is correct. The electric/current closure is now an explicit D0 decision. |
| Q6: ROE rarefaction branch | The cause and needed branch test are adopted. The proposed vacuum return remains incomplete at the caller interface; see T1. |
| Q7: H3+ cooling surface | Both density boundaries and the emission choice are included. The selected constraints must be made mutually compatible; see T4. |
| Q8: time integration and stationary residual | Correctly distinguishes backward Euler, splitting, physical time, pseudo-time, and the unsplit stationary operator. |
| Q9: steady-path state distinctions | The description now distinguishes seed-derived pressure from the pressure reconstructed using the new caloric composition. The option-off reference is also identified. |
| Q10: test construction | Controlled H+ initialization and a mechanically stationary hydrostatic column are appropriate. Tests should assert the state at the actual production call, as specified. |
| Q11: provenance and counters | Archiving before overwriting and identifying `1e-3` as a tolerance are correct. The existing positivity counter's exact meaning still needs precision; see T8. |
| Q12: internal populations and early conservation gates | Population consistency, the trace-molecule limitation, and the elementary cell photon balance now have owners. These are suitable D0/Phase 3 requirements. |

For the corrected energy identity, the newly executed arithmetic gives

\[
\Delta u_{\rm th}=20-13.598434599=6.401565401\ {\rm eV},
\qquad \Delta u_{\rm th}+\Delta E_{\rm chem}=20\ {\rm eV}.
\]

This confirms that the previous double subtraction is removed from the written formulation. It does not claim that a coupled implementation already exists or conserves energy.

The first two reviews identified implementation defects in the unchanged source. Their adoption into this plan does not repair those defects; revision 2 correctly states that distinction.

## 4. Remaining amendments

### T1. A2 must define a vacuum interface that its caller can actually use

Priority: high before implementing A2. Evidence: source inspection and direct execution.

The vacuum condition in revision 2,

\[
v_R-v_L\geq\frac{2(c_L+c_R)}{\gamma-1},
\]

is the appropriate constant-gamma separating-flow condition. However, the proposed instruction to return vacuum speeds does not specify how they enter the actual solver.

`speed_estimate_ROE` returns one `u_star` and two quantities named star sound speeds. `Num_flux` uses them only to form the entropy-correction estimates `u_star-cL_star` and `u_star+cR_star`. It independently constructs Roe averages, wave strengths, and eigenvectors before evaluating the final flux. Returning finite values from the auxiliary routine therefore does not by itself make the final Roe flux admissible near a vacuum.

In a nonvacuum two-rarefaction state, the two sides share a contact velocity. In a separated vacuum, the gas edges instead have two distinct velocities:

\[
v_{\rm edge,L}=v_L+\frac{2c_L}{\gamma-1},\qquad
v_{\rm edge,R}=v_R-\frac{2c_R}{\gamma-1}.
\]

At those vacuum edges, density, pressure, and sound speed vanish. A single contact velocity with two physically zero star sound speeds cannot represent both edges when a finite vacuum gap exists. Encoding arbitrary nonzero numbers as sound speeds could supply the two arithmetic edge estimates, but then the outputs no longer mean what their names and positivity tests claim.

The direct production checks give:

| Input states, gamma=5/3 | Auxiliary sound speeds | Actual `Num_flux` result `(mass, momentum, energy)` |
|---|---|---|
| `(1,-1,1)` and `(1,+1,1)` | Both NaN | `(0, 0.585786437626905, 0)` |
| `(1,-10,1)` and `(1,+10,1)` | Both finite | `(0, 41.8392021690038, 0)` |

The first row adds useful precision to the previous review: NaN auxiliary estimates do not necessarily propagate into a NaN flux. In this build, comparisons involving them leave the entropy correction inactive. The exact interface pressure is `0.224614296367211`, but the difference between that value and an approximate Roe flux is not, by itself, a failure criterion for every nonvacuum Riemann problem.

The second row has exact vacuum edges at `-6.12701665379258` and `+6.12701665379258`. The interface lies inside the vacuum and its exact physical flux is zero. To test the actual numerical flux's admissibility, the probe advances the left cell with its uniform-state physical flux on the outer face, the measured Roe flux on the inner face, and `dt/dx=0.01`. The initial state is `(rho, momentum, energy)=(1,-10,51.5)`. The Courant number is `0.112909944487358`, and the returned state is

```text
rho                     0.900000000000000
momentum               -9.40839202169004
kinetic + thermal      46.2500000000000
thermal internal       -2.92657801877819
```

This is a local demonstration that the existing Roe flux is not positivity preserving for this case. It is not a demonstration that the full EXHALE run accepts the state: the production marching loop has later flux repair and step reduction, which this isolated test intentionally does not invoke. It is also not an execution of the unimplemented revision 2 algorithm.

Required amendment: specify the vacuum signal/return state and its treatment in `Num_flux`, or select a justified flux for this physical limit. Test the final flux and an admissible conservative update, including the actual production repair path when that path is part of the design. Preserve separate checks for nonvacuum star positivity and vacuum states with zero pressure/density. A gate requiring strictly positive star density everywhere contradicts the explicit vacuum case.

The numerical interface can be settled within A2; it does not require redesigning the entire molecular solver. If the implementation changes a public routine signature, present that concrete interface for the user review required by the workspace instructions.

### T2. A3's failed-step status needs a defined retry scope before Phase 2

Priority: high for the retry portion of A3. Evidence: source call order.

The scale correction and recomputation of the residual on the returned state remain reasonable bounded tasks. Returning a status is also useful. However, revision 2 additionally requires rejection of a failed unconstrained carrier step, while the full restoration rule is described later in D0 item 8.

The current `retry_step` block in `EXHALE_main.f90` ends immediately after the hydrodynamic RK stages. The subsequent sequence is composition/temperature extraction, optional element diffusion, `photochemical_transport_step`, ionization equilibrium, and the energy update. Inside the carrier routine, photolysis and transport backgrounds are prepared, `solve_carriers` is called, and `carrier_write_back` currently follows without a status test.

Consequently, adding a return flag alone does not define which already-completed operations are undone, how much physical time remains, or whether mutated background data survive the rejection. A retry cannot be inserted into the existing hydro-only loop by assuming that it already covers chemistry.

Two bounded choices are possible:

1. Implement a carrier-substep retry with a clearly defined frozen background, saved entry state, no write-back from a failed trial, and accepted substeps that cover the original interval.
2. Propagate a failure to an outer attempted-step controller that saves and restores every already-executed affected operator, including element diffusion and rate/column state where applicable.

The implementation can also begin with scale/residual reporting and a clear terminating failure while the selected retry controller is implemented. That increment must be described accurately rather than claiming the rejection/retry behavior is already complete.

Required test: force a solver failure after mutating a trial carrier state, verify that rejected species never reach `carrier_write_back`, and compare the accepted result with the corresponding explicitly subdivided calculation. Check the elapsed physical interval and restored background as well as the returned residual. This test belongs with the Phase 1 retry change, not solely with a later energy integrator.

### T3. CO activation must be recorded across the run, and its impact is not confined to C/O outputs

**Closed on 2026-09-06, and not by a record.** The review's own preferred
remedy was taken: a physically selected destruction mechanism with
photodissociation and shielding (item B3b-CO), after which the thermodynamic
ceiling was deleted (item CEILING-DEL) together with `pct_co_ceiling` and the
cumulative counters. There is no ceiling-active run to accept or exclude, and
what a run reports instead is the destruction model's domain. The review text
below is kept as written.

Priority: high for acceptance of any ceiling-active run. Evidence: source counter semantics and coupled energy/composition paths.

Revision 2 proposes to identify a run in which `pct_co_ceiling > 0` and exclude its outputs from quantitative C/O results. That is not sufficient to implement the stated restriction.

The actual limiter resets `pct_co_ceiling=0` at every call. It increments the counter only when the relative correction exceeds `limit_report=1e-10`, although the CO value is changed for any positive excess. Thus a zero returned count does not establish that the ceiling was never used, or even that no value was changed in that call. A later inactive call also erases the previous activation count.

Record at least an activation indicator over accepted evolution, the amount of CO removed, and the time/cells affected. Keep rejected-trial counts separate when retries are introduced. If the physical condition is “the ceiling never acts,” use that actual condition; a thresholded diagnostic count is not equivalent. A recorded correction magnitude can explain rounding or numerical tolerance without silently changing the acceptance definition.

The output restriction should also cover the coupled atmosphere. Changing the CO reservoir changes atomic C and O available for cooling, can change other carriers through the limiter, and changes chemical energy unless a consistent source is provided. Temperature, flow, mass loss, and non-C/O observables can therefore depend on this incomplete closure. Restricting only quantitative C/O reporting leaves a route to claiming a validated bulk wind from an unclosed energy/composition model.

Until the D0 constraint model is justified, treat a ceiling-active atmosphere as a diagnostic experiment and state which independently supported quantities, if any, can still be used. The sign of the complementarity sink and conservation of its products are necessary, but the one-sided equilibrium ceiling still needs a physical timescale/domain argument. An equilibrium constant alone does not determine how quickly CO can be destroyed.

The setup report can state that the selected configuration permits an unvalidated ceiling. Actual activation belongs in runtime/final status and persisted output metadata. The current code rewrites resolved configuration at the end as well as initially, so there is an available final reporting path; a cumulative state still has to be supplied to it.

### T4. B2 must choose compatible absolute emission, departure factors, and temperature coverage

Priority: medium to high before B2 implementation. Evidence: the local published Miller tables and explicit arithmetic.

Revision 2 correctly separates temperature continuity, density continuation, and emission normalization. However, its remaining choice between unscaled and scaled LTE emission cannot be combined arbitrarily with a requirement to reproduce Table 6's absolute emission and approach a departure factor of one.

At 5000 K, the published values used in the plan are:

```text
Table 4 unscaled LTE                  3.4119e-18
Table 4 scaled LTE                    4.7688e-18
Table 6 emission at the last density  4.7587e-18
Units: W molecule^-1 sr^-1
```

The directly calculated ratios are

\[
E_6/E_{\rm LTE,unscaled}=1.394736070811,
\qquad E_6/E_{\rm LTE,scaled}=0.997882066767.
\]

An unscaled LTE denominator cannot reproduce that absolute Table 6 emission with a subthermal factor approaching one. Multiplying the unscaled LTE model by the published dimensionless factors is a different absolute-emission model; it should be described and validated as such. If the goal is the published scaled emission, use the corresponding LTE definition and compatible departures. The independently tabulated emission values and rounded factors also need a stated data-consistency tolerance; do not force every rounded entry to be an exact identity. [Miller et al. 2013, Tables 4-6 and the cooling/non-LTE methods](../../references/Miller_2013_JPCA_117_9770.pdf).

There is a separate lower-temperature dependency. Table 4 begins at 500 K, while the existing LTE routine is used down to 30 K and has its largest reported join problem at 300 K. Table 6 supplies an additional 300 K anchor, but sparse anchors do not determine the complete curve below them. A smooth positive fit through Table 4 alone can satisfy continuity while remaining physically unsupported over part of the advertised range.

Required amendment: specify the temperature range to be supported, the emission data or state/transition calculation that supports each interval, and the physical extrapolation or restriction outside it. Do not advertise a validated low-temperature extension until those data or a derivation exist. Keep mathematical continuity, agreement with sampled data, and uncertainty between samples as separate gates.

These are B2 input/model decisions, not reasons to defer A1 or A8. The revised high-density continuation and collisional low-density exponent can be retained.

### T5. Define external radiative exchange locally before using the new energy identity

Priority: medium; a D0 specification requirement. Status: an ambiguity to resolve, not a confirmed new implementation error.

The corrected energy convention is sound. The word “escaping” in Section 4.3 nevertheless needs an explicit boundary: leaving the local material system is different from leaving the entire atmosphere.

For example, let one cell emit a 10 eV photon that is fully absorbed in another cell. The material exchange is -10 eV at the emitter and +10 eV at the absorber; the domain sum is zero. If only photons escaping the outer simulation boundary are subtracted at emission while all internal absorptions are added, this example creates 10 eV of material energy.

D0 should define emission, absorption, any explicit radiation-energy storage, and any eliminated diffuse field using the same spatial ownership rule. For an explicit radiation field, include its energy change and boundary flux in the domain budget. For an eliminated field or local recycling approximation, derive the corresponding net material exchange and ensure the same photon is not reintroduced as an independent source.

The plan already names reabsorbed photons and once-only recipients, so this is a clarification of the implementation contract rather than a reason to reject its energy equation. Add a two-cell emission/reabsorption test beside the elementary cell absorption gate. It is small enough to detect an ownership error before the full band-overlap validation in Phase 6.

The same D0 document should explicitly choose the independent species variables or constraints for C1. Species continuity, total mass continuity, element conservation, and charge constraints must not produce a redundant stationary system. Naming a constrained chemical Jacobian is appropriate; specifying its independent space before building the global Newton rows completes that requirement.

### T6. Preserve the requested audit sources when creating maintained tests

Priority: required for the user's archival request; medium for test interpretation.

Section 3.3 says the audit programs will move to `src/tests/physics_probe/`. The user explicitly requested that all code created for these checks be retained in `docs/audit_20260905/`. Keep those reference sources available there. Maintained regression or unit tests may be derived in `src/tests/`, with their origin and altered acceptance conditions recorded, while the historical probes remain available for comparison.

Not every printed diagnostic can be converted into a passing production assertion by attaching a tolerance. For example, `photoevent_energy_check.py` deliberately evaluates the superseded revision 1 equation as a counterexample. Asserting that its incorrect thermal increment agrees with the correct one would create a permanently failing test unrelated to a repair of the production source. The maintained test must evaluate the accepted current equation or the actual source routine instead.

Likewise, an exact analytical reference, an approximate Riemann estimate, a published measurement with uncertainty, a solver residual, and an invariant require different acceptance rules. Use dimensionally meaningful absolute and relative tolerances, explicit NaN checks, and the right reference for each quantity. A failed physical invariant must not be hidden by a tolerance intended for an approximate emission fit.

The current E1 program also does more than channel-sum and stoichiometric checks: it prints join behavior, compares the option-off construction with an older split, and evaluates a spectrum-weighted ledger. What it lacks is the detector-observable inversion test. Amend that description and add the missing test without misidentifying existing coverage.

### T7. A8 is bounded, but its integration gate must follow channel arrays into physical sources

Priority: medium. Evidence: source consumers of the channel routine.

The revised detector inversion is correct. For a measured detector ratio r and a double-to-single event ratio rho_D, the corresponding double fraction among dissociative detector events is

\[
q_D=\frac{\rho_D}{1+\rho_D}\frac{1+r}{r}.
\]

It must satisfy `0 <= q_D <= 1`; incompatible input models should produce a diagnostic/failure or a stated model restriction rather than silently redefining the measured observable by clipping. This follows the nonnegative-S gate already proposed in A8.

Rename `frac_H2_double_of_protons` when its return value becomes an event fraction. The existing name describes a different physical quantity. Source search found its production use within `h2_photo_channels`; the inspected higher-level consumers call `h2_channel_cross_sections`, so the internal rename has a narrow scope.

The routine test should be followed by a small energy-grid/source-assembly check. `set_energy_vectors.f90` stores the channel arrays, and `util_ion_eq.f90` integrates them into ionization, secondary-ionization, and heating terms. Verify that an ionizing event still destroys one H2, a D event produces two protons/electrons, the opacity remains the channel sum, and source consumers use the corrected arrays. A8 changes integrated heating through changed probabilities even if the energy assignment of an individual channel is deferred to A5.

This integration check can use a fixed state and selected photon energies. It does not require the full planetary regression suite and must not claim that A5's existing energy defects are already corrected.

### T8. Smaller factual and reporting corrections

The existing `n_faces_flux_positivity_limited` increments only if `positivity_limited_fluxes` returns `repaired=.true.`. It counts repaired interfaces in successful repair calls, not every invocation. A later RK-stage failure can still lead to rejection of an attempt containing an earlier successful repair. If A1's purpose is to distinguish execution, successful repair, and accepted evolution, report separate quantities or name the existing counter precisely. A zero count is not proof that the routine was never called.

Revision 2 attributes Yan's vertical-threshold sentence to p. 1049. In the local published PDF it is on printed p. 1048, in Section 4. The 51.4 eV value itself is correctly quoted. The previous review cited the section without this incorrect page number. [Yan et al. 1998, p. 1048](../../references/Yan_1998_ApJ_496_1044.pdf).

Koskinen's Eqs. (B9)-(B10) explicitly show the electric-field term and the separate zero diffusive mass condition. They do not by themselves state the zero-current equation written in D0. Zero current is a defensible candidate closure for an appropriately specified model, but derive it as the selected EXHALE constraint rather than claiming that those two equations alone establish it. [Koskinen et al. 2022, Appendix B](../../references/Koskinen_2022_ApJ_929_52.pdf).

The opening sentence saying every paper claim was rechecked conflicts with the later explicit statement that several paper readings were not repeated. Retain the useful disclosure and scope the opening verification claim to what was actually checked.

Finally, the plan says the A8 channel-probability change and B2 high-density change are unmeasured until run. Their local values have already been measured by the retained probes. Their effect on an atmospheric solution is unmeasured. Use those two statements separately.

## 5. What D0 should now deliver

D0 should be the next shared design artifact. Its acceptance can be made concrete without demanding a complete implementation first:

| Required D0 result | Acceptance question |
|---|---|
| Material and radiation energy inventory | Does each event close chemical, thermal, excitation, fast-particle, and radiation exchange with a defined local boundary? |
| Species and reaction representation | Are independent variables and algebraic constraints specified, with exact element/charge invariants and consistent forward/reverse thermodynamics? |
| Thermal versus chemical spatial flux | Does summing the thermal and chemical finite-volume equations recover the intended material-energy equation using the same species face fluxes? |
| Charged transport | Is the electric/current closure derived for the chosen model, or is the common-velocity approximation justified over an explicit domain? |
| CO treatment | Is there a physical validity argument, a source with the correct sign and products, and a matching energy term, or an exclusion of the affected calculation? |
| Internal populations | Do EOS, chemistry, degradation, and emission assumptions describe a consistent state within the accepted density/temperature/composition range? |
| Boundary and time discretization | Are conditions counted for the full independent system, and are physical time integration, rollback, and stationary residuals specified separately? |
| Verification matrix | Does every accepted closure have a local test with units, tolerances, and a physical reference, followed by the appropriate coupled test? |

The existing phase order should be retained. The main sequencing change is to attach the retry/rollback contract to the first implementation that introduces failure recovery, rather than relying on a later D0 paragraph to complete Phase 1 behavior. H3+ data selection is a named prerequisite of B2. CO activation history and the restricted interpretation of affected runs accompany any temporary use of that model.

## 6. Verification outcome and retained artifacts

The new diagnostic compiled and ran successfully. It reproduced a finite full flux despite NaN auxiliary estimates in the moderate expansion, and negative internal energy in the deliberately isolated Euler update for the stronger vacuum-producing expansion. These results support T1's requirement to test the full consumer and state update. They do not claim that a complete EXHALE atmosphere fails after its existing repair logic.

The arithmetic checks confirmed revision 2's corrected H photoionization identity and the incompatibility of reproducing the cited Table 6 absolute emission with the unscaled LTE denominator and a subthermal factor approaching one. The two-cell calculation illustrates the radiative-boundary distinction; it does not diagnose a new executed radiation-transfer bug.

To reproduce the added checks from `EXHALE_v1.00`:

```bash
bash docs/audit_20260905/run_review3_checks.sh
```

The [audit README](audit_20260905/README.md) lists the retained programs from all three reviews. Build products are isolated under `/tmp`; the runner prints the location. Exit status reports successful execution of these diagnostic programs, not acceptance of the physical defects they expose.

No claim is made here about newly measured mass loss, atmospheric profile changes, base oscillation amplitudes, or convergence of a molecular steady wind. Those measurements remain with the corresponding implementation and coupled-validation gates. The reviewed plan and earlier review documents were preserved.
