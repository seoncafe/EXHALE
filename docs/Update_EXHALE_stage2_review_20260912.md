# Review of recent stage-2 implementation changes

Date: September 12, 2026  
Review target: the current working tree, not the update log alone  
HEAD: `a7c18a620f1b` (full revision retained in the audit archive)

## 1. Executive assessment

The recent work makes substantive improvements: the stationary element relaxation now uses an explicit failure result, the carrier transport trial has an enforced movement test, and the Newton loop checks certified rows before stopping. Those improvements do not establish that every caller, returned thermodynamic state, or certification path is correct.

This review found two high-priority problems in the recent stationary composition path: rejected search trials can permanently prevent subsequent certification, and the fixed-pressure molecular update is not consistent with simultaneously retaining the conserved energy of the wind. It also reproduced a temperature lag after the newly added chemical refresh. An important limitation of the new tests is that their synthetic molecular column carries approximately 31.49% more mass than its prescribed density.

The strongest results below come from calls to freshly compiled production routines. The 500-cell molecular measurement uses the actual saved regression fixture and is mass-closed; it is not dependent on the defective synthetic column. Production source was not modified during this review.

| ID | Finding | Evidence | Priority |
| --- | --- | --- | --- |
| R1 | Rejected stationary transport trials contaminate the carrier-history certification flag | Reproduced on a mass-closed column; rejection changes no composition | High |
| R2 | Fixed-pressure carrier updates and unchanged conserved energy describe different molecular states | Reproduced with the 500-cell molecular fixture and production EOS maps | High: resolve the state contract |
| R3 | The chemical refresh returns a lagged temperature and does not establish a joint chemical/thermal fixed point | Source inspection and direct returned-state measurements | Medium; address with R2 |
| R4 | Element-step safety checks are conditional on requesting `status`; the marching caller does not request it | Complete call-path inspection; no failing physical trajectory reproduced here | High validation and failure-handling gap |
| R5 | The mass-tolerance ceiling permits a mass row that the stated unresolved-row policy should not certify | Direct execution of the production cell verdict | Medium; extreme low-flux regime |
| R6 | The molecular retry test column is not mass-closed, and its new consistency assertions are too weak | Direct census: relative mass discrepancy `0.3149430466`; all existing assertions still pass | High for validation quality |
| R7 | The outer failure summary pairs a maximum mass row with another cell's tolerance | Source inspection of the new cell-dependent tolerance and its consumer | Low: incorrect diagnostic attribution |

The findings do not demonstrate that a particular published atmospheric result is wrong, that every current run fails, or that partitioning is intrinsically unsuitable. They identify specific implementation and verification contracts that remain incomplete.

## 2. Scope, source verification, and execution

The principal scope is Section 8 of [Update_EXHALE_stage2.md](Update_EXHALE_stage2.md): P1–P22, especially P1–P3, P6, P14–P17, and P21. The terminology changes were checked where they affect those paths. The preceding log was used for context, not treated as evidence that a routine currently behaves as recorded.

The inspected paths include:

- `steady_wind_with_element_diffusion` and its composition updates in `src/EXHALE_main.f90`;
- element nonlinear solves, projections, optional acceptance checks, and relaxation in `binary_element_diffusion.f90`;
- carrier interval integration, relaxation, chemical refresh, and history flags in `diffusive_photochemistry.f90`;
- the temperature input and composition output of `ioniz_eq`;
- species-density reconstruction and the caloric EOS maps;
- Newton stopping, returned-state certification, mass-row scales, and cell-dependent tolerances;
- the interface data returned by the kind-generic hydrodynamic rows;
- the relevant test drivers and their build scripts.

A fresh private production build used GNU Fortran 16.2.0 with the Makefile's `-O3 -fopenmp` settings and OpenBLAS under `/opt/miniconda3/lib`. Its objects are in `/tmp/exhale_stage2_review_20260912.XJbhQc/build`. Custom diagnostic drivers used `-O0 -g -fcheck=all -fopenmp`; those checks do not retroactively instrument the optimized production objects. The existing `carrier_retry` script separately rebuilt its source dependency set with its own checking flags.

Directly executed checks:

| Check | Observed result | What it establishes |
| --- | --- | --- |
| Fresh production build | Successful | The reviewed source links through the canonical build |
| Existing `carrier_retry` suite | 127 PASS, 0 FAIL | Existing assertions pass, including the recently added checks |
| Existing `element_operator` suite | 24 PASS, 0 FAIL | The tested stationary element acceptance and boundary cases pass |
| Returned-state probe derived from `carrier_retry` | Existing driver assertions pass; additional audit measurements expose R3, R5, and R6 | Passing current assertions does not cover these contracts |
| Fresh stationary-retry history probe | History flag changes from true to false with zero accepted steps and zero composition change | R1 is reproducible independently of an invalid initial mass normalization |
| Production-module molecular fixture probe | One hydro call and one bounded carrier update complete; total approximately 9.04 s | R2 and R3 also occur on a mass-closed saved atmosphere |

The molecular fixture probe explicitly selects the hydro-only registry, then calls the current carrier relaxation interface. It does not invoke the production main program's entire outer controller. Its copied input still contains the fixture's original coupled-solve key; the standalone driver's registry calls determine the experiment. The source inspection establishes how the production main program passes the same primitive and conserved quantities.

No full regression suite, long convergence continuation, or unrelated data generation was run. The earlier log's pass counts, timings, and certified atmospheric results are not presented here as freshly repeated measurements.

## 3. R1: a rejected stationary search trial permanently blocks certification

### Implementation and reproduction

In [diffusive_photochemistry.f90](../src/modules/lower_atmosphere/diffusive_photochemistry.f90), `relax_photochemical_composition` calls `photochemical_transport_step` around line 5002. The wrapper does more than return interval success. On an exhausted interval in initialization mode, it sets `carrier_history_certifiable_flag = .false.` at line 1847.

The outer relaxation then restores the trial composition and shortens the step. This is a rejected nonlinear-search trial, not a physical interval the run adopted without integrating it. Nevertheless, the history flag remains false. In [certification.f90](../src/modules/time_step/certification.f90), lines 1204–1242, certification reads that flag and refuses the state independently of its current row residuals. There is no context distinction at that read between stationary search history and a physical trajectory.

The [mass-closed history probe](audit_20260905/stage2_review_20260912/stationary_retry_history_final.log) initializes the flag normally and injects interval rejection through the existing test hook:

```text
entry_mass_closure = 2.2204460492503131e-16
history_before = T
history_after = F
composition_change = 0
kept_steps = 0
outcome = 3  [carrier_relax_interval_refused]
```

This is a controlled rejection test, not a claim that this rejection frequency occurs in an ordinary atmosphere. It proves that the rejection semantics are incorrect even when the entry mass normalization is valid.

### Consequence and correction

A later state can satisfy every stationary equation and still be permanently disqualified because a discarded trial failed. The error makes certification too restrictive, not falsely permissive. It also makes the new retry contract depend on the wrapper's unrelated run-mode policy. Under physical mode that wrapper can stop the process before the outer relaxation gets a chance to shorten its trial, unless the stop has been suppressed.

The stationary relaxation should use the interval integrator's returned status directly, as the physical-step controller already does for its own retry decisions. Keep attempt statistics, but separate them from violations in an adopted physical history. Do not fix this by globally clearing genuine trajectory failures or by ignoring the history flag everywhere.

Required tests: rejected stationary trials leave the state and stationary eligibility unchanged; a successful retry remains eligible; an actually skipped or invalid adopted physical interval still marks the physical history. Also test the same stationary trial independently of the global run-mode setting.

## 4. R2: fixed pressure is not fixed conserved energy with molecular heat capacity

### Physical and implementation issue

With the molecular caloric EOS, pressure is a function of both thermal energy and composition. Holding `(rho, v, p)` fixed while changing molecular content generally changes the thermal energy implied by that state. Holding `(rho, rho*v, E)` fixed instead changes pressure and temperature as the composition changes. These are different conditional problems.

This is directly implemented in [caloric_eos.f90](../src/modules/states/caloric_eos.f90), lines 659–700: `pressure_from_energy_density` and `energy_density_from_pressure` use the composition-dependent molecular heat capacity. [UW_conversions.f90](../src/modules/functions/UW_conversions.f90), lines 84–109, uses the former map in `U_to_W`. The composition data those maps read are refreshed by `get_species_densities`, [composition.f90](../src/modules/functions/composition.f90), line 135.

The new carrier routine explicitly holds pressure fixed. However, its production caller around [EXHALE_main.f90](../src/EXHALE_main.f90), lines 6460–6466, retains `u(3,:)` and only refreshes particle densities and `T` from the old `p`. It does not convert that fixed-pressure candidate into a consistent new conserved energy. The subsequent hydro solve consumes the retained conserved variables, whose pressure is now different under the updated composition.

### Direct measurement on a valid atmospheric fixture

The [500-cell probe log](audit_20260905/stage2_review_20260912/current_probe_run.oY88nt/molecular_wind_state_probe.log) records the state after a hydro-only production solve and one bounded carrier relaxation with `trust = 0.01`:

| Quantity | Direct measurement |
| --- | ---: |
| Carrier steps retained | 5 |
| Returned carrier drift | `9.9955477323e-3` |
| Relative mass closure after the carrier update | `8.3323877256e-15` |
| Maximum relative difference between held pressure and pressure reconstructed from unchanged thermal energy at the new composition | `3.2993113229e-3` |
| Maximum relative change of thermal energy required to retain the held pressure at the new composition | `3.3249715434e-3` |

These are local maximum relative differences, not global energy-budget errors. The density and species normalization are valid in this reproduction. The result therefore cannot be dismissed as an artifact of the synthetic test problem in R6.

### What this establishes, and what it does not

The problem is the inconsistent state contract, not the size of the discrepancy. A fixed-pressure species predictor can be an intentional nonlinear iteration. It must then be identified as such, and the conversion into the state consumed by the full residual must be explicit and consistent. It is not the same as solving the species block at fixed conserved hydrodynamic variables.

This experiment does not prove that the final stationary root is wrong or that a certified result has violated a global energy budget. It shows that the current claims about a fixed wind and the actual conserved-variable handoff are not interchangeable. In particular, P21's successful conditional fixed point does not establish a fixed-energy species solution or a joint stationary atmosphere.

Choose and implement one contract:

1. For a block that holds conserved hydrodynamic variables fixed, recompute the caloric composition, pressure, and temperature from the unchanged `u` as chemistry and transported composition change.
2. For a block that intentionally holds primitive pressure fixed, form a consistent candidate `u(3,:)` from the updated composition before it enters the next full-state evaluation. Treat this as a nonlinear predictor, not as an unaccounted physical energy source or a physical-time step.

The choice changes the solver's conditional problem and potentially its interface; it should be agreed before implementation. Test the EOS round trip at every handoff and compare the two conditional formulations from the same valid atmosphere before attributing front movement solely to physical feedback.

## 5. R3: the new chemical refresh does not return the temperature of its output composition

`equilibrate_chemistry_at_fixed_pressure`, [diffusive_photochemistry.f90](../src/modules/lower_atmosphere/diffusive_photochemistry.f90), lines 4804–4830, does the following:

```text
get_species_densities(entry composition)
comp_T_from_p(held pressure, entry particle counts)
ioniz_eq(that temperature, composition in/out)
return
```

The temperature argument of `ioniz_eq` is `intent(in)` at [ionization_equilibrium.f90](../src/modules/radiation/ionization_equilibrium.f90), line 736. The chemical sweep changes `f_sp`, but cannot repair `T` for the changed particle count. The refresh routine performs no density/temperature reconstruction afterward.

The directly measured discrepancy between returned `T` and `p/(n_tot+n_e)` is `5.5997819434e-9` in the mass-closed 500-cell case. The synthetic test gives `2.5674706456e-7`, but that number belongs to the invalid-mass fixture described in R6 and is not used as the physical-atmosphere measurement.

The main caller subsequently reconstructs `T`, which repairs the temperature column at that location. It does not make the internal carrier steps use an exactly consistent refreshed background retroactively. `carrier_state` explicitly reads `bg_cell%T_K` and `%ntot`; the chemical sweep stores its cell background before completing the local solve. Thus a fresh sweep is an improvement over a background frozen for the entire relaxation, but is not automatically a closed thermochemical state.

Two additional checks are missing from the new relaxation contract:

- The refresh does not request the available `sweep_ledger` result or reject a failed/nonfinite chemical trial before another transport step consumes it.
- The fixed-point ending tests movement of solved carrier columns at a sufficiently large step, not the eliminated chemical balances, all thermodynamic variables, or the refreshed stationary carrier residual itself. The full outer certification still matters.

First resolve R2's held-variable choice. Then close the chemical and EOS reconstruction on that same choice, use the chemical solve's validity result, and test a defined residual for the returned state. A post-sweep temperature assignment alone does not make rate coefficients evaluated at the previous temperature exact. If a lagged sweep is intentionally retained, state its approximation and control its error by what the downstream residual reads; do not describe it as an exact joint fixed point solely because the carrier displacement is small.

## 6. R4: the element safety repair is not applied to the marching caller

The P1 repair is real but caller-dependent. In [binary_element_diffusion.f90](../src/modules/functions/binary_element_diffusion.f90):

- Line 650 sets `judged = present(status)`.
- Around line 735, failure of `solve_mass_fraction` restores the entry state only if `judged` is true.
- Around lines 855–871, finite-state, helium-bound, and mass-closure acceptance checks also run only if `judged` is true.

The stationary relaxation supplies `status=st`. The marching call at [EXHALE_main.f90](../src/EXHALE_main.f90), line 2297, does not. An unsuccessful inner solve therefore does not force rejection at that boundary on the marching path.

This is not resolved merely by the existence of the outer physical-step controller. Its acceptance routine, `certification_evaluate_physical_step`, receives the carrier interval result but no element-step result. It checks finite values and other conditions, but cannot directly reject the unreported element nonlinear failure. The elemental invariant check is explicitly skipped when element diffusion or molecular carrier transport is active; the replacement global transport-budget check is not implemented there. The recorded `he_fraction_*` diagnostics are saved by the checkpoint code, not used as an element-solve acceptance result in that path.

No failing physical trajectory was reproduced in this review. The confirmed finding is the unguarded call and the absent propagation of its solve result, not a measured error in every physical run. The 24 passing element tests exercise the checked interface and do not close this gap.

Acceptance of a numerical update should not change because a caller omitted an optional diagnostic argument. Make validity enforcement unconditional, propagate the outcome to the physical-step controller, and reject or retry before chemistry consumes an invalid element update. The existing physical-step rollback framework is the appropriate consumer. Test finite but unconverged inner solves, bound excursions, and mass-closure failures through that caller, not just a direct call with `status` present.

## 7. R5: the mass-row ceiling lacks the stated unresolved classification

The adaptive mass tolerance is an intentional user decision recorded in the log. This review does not propose silently reverting that decision. The issue is its boundary behavior.

The comment beside `cert_tol_mass_ceiling` in [certification.f90](../src/modules/time_step/certification.f90), lines 221–229, says that a state whose estimated floor reaches this regime cannot be judged, rather than passing. The implementation of `mass_row_cell_verdict`, lines 793–814, instead clips the tolerance at 1 and checks only `q/tol < 1`.

Direct execution of the production routine gives:

```text
q = 0.5, floor_q = 1
tol = 1, distance = 0.5, within = T, anchored = T
```

Therefore the mass entry can pass with a normalized imbalance of one half when the estimated arithmetic floor is too large to resolve its balance. This is a cell-verdict counterexample, not a complete atmospheric state that passed every certification entry. No ordinary fixture in this review was shown to reach that regime.

Return an explicit unresolved result when the floor is too large to support the requested physical statement. Preserve a distinction between a row that is resolved and within tolerance, a row dominated by estimated arithmetic uncertainty, and an unresolved row. Do not call the last category certified merely because a clipped number is below 1.

The existing ceiling test uses `q = 2`, which naturally fails against a tolerance of 1; it does not test the problematic `0 < q < 1` interval at the ceiling. Add that case, nonfinite estimates, and extremely small physical fluxes. Separately, the empirical rounding estimate and margin should be validated across reconstruction, flux, molecular EOS, and near-zero-flow regimes before being treated as a universal bound. A perturbation response measured on several fixtures is useful evidence, not a proof of such a bound.

## 8. R6: the recent molecular consistency tests do not start from a physical mass normalization

`build_molecular_hydrogen_column` in [carrier_retry.f90](../src/tests/carrier_retry/carrier_retry.f90), around lines 1326–1338, assigns H I approximately 0.20, H2 approximately 0.40, and He I approximately 0.0793 directly to fractions whose layout is a number count relative to the mass-density normalization. The hydrogen entries already account for approximately unit mass; adding helium without renormalizing makes the total exceed the prescribed density.

Using production `element_nuclei_and_charge` on this unchanged fixture gives:

```text
max |rho_species / (rho * n0) - 1| = 0.31494304662695316
```

The [measurement](audit_20260905/stage2_review_20260912/carrier_state_contract_final.log) is not inferred from the input literals alone. The fixture calculates reconstructed density in its existing test code, but its conservation assertions compare element counts before and after; they do not require the entry to reconstruct the prescribed mass density.

The existing suite still passes 127 assertions. Many checkpoint and rejection tests remain useful because those properties can be tested independently of a realistic atmosphere. However, the newly added pressure, chemistry, and background-consistency statements need a mass-closed fixture before their quantitative results can support physical correctness.

Two more weaknesses in those assertions matter:

- Around lines 1100–1131, background consistency is tested by being closer to the returned composition than to the initial one. That is an improvement test, not an absolute consistency or residual test.
- Around lines 1153–1161, an absolute maximum change in species entries is compared with `sqrt(epsilon)`, justified by the nonlinear solver's `xtol`. A nonlinear solver's relative step criterion is not, by itself, a bound on absolute changes in every output species or on the reaction residual. The test needs an independently justified scale and a residual check.

Use a single mass-closed column constructor, assert density reconstruction before advancing, then test the actual EOS and chemical residuals afterward. Include both H2-only and transported-ion configurations. The 500-cell R2 measurement already supplies an independent mass-closed counterexample to the state-contract claim, so correcting this fixture would not remove that finding.

## 9. R7: the outer mass-failure message combines values from different cells

After P16, `mass_row_verdict` preserves `row_max` and `jworst` for the largest normalized mass row, but stores `tol`, `row_at_bind`, `jbind`, and `dist_bind` for the cell furthest outside its own tolerance. These cells need not coincide.

The outer summary in [EXHALE_main.f90](../src/EXHALE_main.f90), around lines 6315–6340, ranks a nongated entry using `row_max / tol` and identifies `jworst`. For mass this pairs the maximum row with the binding cell's tolerance. The actual certification verdict remains the one produced by `mass_row_verdict`; the bug is in the reported worst refusal and its severity.

Use `dist_bind`, `row_at_bind`, and `jbind` for mass failure attribution, or provide a common accessor that returns the binding measure, tolerance, distance, and cell for every entry. Retain `row_max` as a separate descriptive maximum. A unit case where the maximum-row cell differs from the binding cell is sufficient to test this mapping; no atmosphere needs to be regenerated for a reporting-only correction.

## 10. Improvements that remain justified after the recent repairs

### Corrected paths should remain corrected

The checked stationary element interface restores failed trials; the current carrier movement test rejects excessive transport trials before retaining them; and the JFNK loop around line 15356 now reads certified rows before stopping. The inspected kind-generic wrappers return and store their own face departures, and reject intermediate reconstruction continuation. The old D1–D3 and stale-face-data findings should not simply be repeated as if those changes did not exist.

Their remaining limits are different: missing callers in R4, the new history side effect in R1, the thermodynamic contract in R2–R3, and untested edge classifications in R5.

### Control the full accepted state, not only the species progress indicator

The current outer progress controller follows the worst gated species row. That is useful for the front, but is not a measure of progress of the complete system when the species rows already pass and a hydro or chemical condition still fails. Use the maximum properly normalized distance of all required entries, with an explicit policy for unavailable or nonfinite entries. The mass distance must use its binding cell, not a maximum divided by an unrelated tolerance.

Helium-only `omega` remains explicit in the source; it is not a uniform damping factor for all transported elements. If the controller is intended to shorten the whole composition change, define a conservative update of the independent transported quantities rather than assuming the current helium blend does so.

### Bound and accept the post-refresh candidate

P21 tests the movement bound before the chemical refresh. The inspected chemical path intentionally holds the transported partition, and the measured default molecular update remained inside the bound. No post-refresh bound violation was reproduced here. Nevertheless, the public contract describes the returned state, so its bound, finite values, elemental budgets, and chemical validity should be checked after all transformations that produce that state. A rejected refreshed candidate needs restoration of the corresponding background state as well as `f_sp`.

### Distinguish a conditional fixed point from a certified atmosphere

P21's recorded unbounded fixed-pressure experiment and its long production continuation were not repeated here. The new measurements do not overturn those recorded numbers; they limit their interpretation. Before assigning a molecular front shift to shielding or hydrodynamic feedback, compare consistent fixed-energy and fixed-pressure formulations and require the same refreshed full residual at the final endpoint. Do not equate a small carrier displacement or a successful process exit with that endpoint.

## 11. Recommended order and acceptance criteria

1. Separate stationary retry attempts from physical-history failures, preserving real trajectory marks. Require the R1 reproduction to retain eligibility without adopting a failed step.
2. Agree on the held hydrodynamic variables for molecular composition updates. Enforce EOS-consistent conversion and chemical validity at every handoff; rerun the valid 500-cell R2 probe on the selected contract.
3. Propagate element solve failure into every affected caller, especially the physical-step controller. Do not allow optional reporting to select whether validity is enforced.
4. Replace the invalid molecular test normalization and strengthen its thermodynamic and chemical assertions. Keep the existing checkpoint tests, but do not use their pass count as a substitute for the new physical checks.
5. Define unresolved mass-row behavior at the ceiling and test the missing edge cases without undoing the user's ordinary rounding-anchor decision.
6. Correct binding-cell diagnostics and the outer progress metric. Only then spend a larger convergence budget comparing solver approaches.

No production changes, API changes, tolerance changes, or fixture edits were made in this review. The code changes suggested above require their own implementation and focused validation.

## 12. Artifacts and limitations

All new diagnostic source and run products are retained in [stage2_review_20260912](audit_20260905/stage2_review_20260912/README.md), under the audit location previously requested by the user. The archive includes the source diff relative to HEAD and SHA-256 hashes of the reviewed Fortran and include files. The source diff was unchanged between its capture and the final source check.

This is not an exhaustive review of every reaction rate or every change since September 5. No new literature-fit validation, global elemental-flux budget, global energy budget, spatial convergence study, oxygen/carbon molecular continuation, or failing physical trajectory was performed. The new numerical findings concern the directly exercised contracts described above. Existing source changes and simulation products were preserved.
