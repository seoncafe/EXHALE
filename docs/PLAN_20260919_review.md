# Critical review of PLAN_20260919

Date: 2026-09-19 (KST)

Reviewed plan: [PLAN_20260919.md](PLAN_20260919.md), MD5 `46d65e76f3cbcf528828d406d94c78ad`.

## 1. Overall judgment

**The plan identifies important unfinished work, but its causal claims, grouping of defects, and execution order need revision before implementation.** In particular, approve investigation of the molecular boundary closure, not the conclusion that a tighter ionization tolerance alone fixes it. Move termination classification ahead of further campaigns. Treat the failed Jacobian action test as an unresolved numerical question, not automatically as a harness defect.

The following parts are sound:

- Preserve existing generations and do not replace unresolved states with successful-looking labels.
- Investigate the molecular restart refusal before launching another expensive molecular campaign.
- Bound the proposed Newton investigation and distinguish diagnosis from authorization to change solver behavior.
- Retain the accepted D2b boundary decision and the pinned grid inputs while measuring these issues.
- Report failed tests rather than silently adjusting their expectations.

The principal corrections are summarized below.

| Finding | Priority | Review judgment |
| --- | --- | --- |
| R1. P1 does not isolate the cause; its recommended fix omits coupled temperature and state-refresh effects | High | Combine P1 with the relevant P6 closure work; require a controlled attribution experiment |
| R2. The current ghost diagnostics omit molecular electrons | Medium, relevant to R1 | Confirmed diagnostic inconsistency; use the existing full electron count |
| R3. P4 can classify an unfinished inner solve as a completed solve | High | Reproduced; repair before campaigns and update the parallel parser in the publisher |
| R4. A full certification report silently overwrites its last entry | High consequence, latent trigger | Reproduced through the public production API; not merely a test-harness issue |
| R5. The coupled Jacobian suite has two qualitatively different failures | High for affected numerical work | Large action discrepancy needs diagnosis; last-bit scale inequality needs a separate numerical contract |
| R6. Excluding all ghost information is not a sufficient CFL design | Medium | Minimize over evolved control volumes, but retain boundary-face wave-speed constraints |
| R7. Two P6 implementation claims are outdated or overgeneralized | Medium | Solution-state transit support and early line-search termination already exist |
| R8. P5's output-directory claim is not reproduced on the current scripts | Medium | Seventeen shell targets read the override; omitted explicit forwarding does not erase an inherited environment variable |
| R9. P4's missing provenance is not confined to old imported evaluations | Medium | A newly solved molecular generation also has no recorded parent |
| R10. P7 and the catalog counts need explicit scope | Medium | A single molecular refusal does not establish failure of every molecular case; count a named case set |

This review recommends changes but does not implement them. No production source, input, catalog state, index, or golden output was changed.

## 2. Verification performed

### 2.1 Binary and working-tree scope

Direct checks found:

- `EXHALE.x` MD5: `2c3b0acc9983aec03bb4f844294fed18`.
- `LHS1140b/models/EXHALE_2c3b0acc.x` has the same MD5.
- `make -q` returned 0. No production rebuild was performed.
- The tree already contained extensive modified and untracked files. Those changes were preserved.

Hashes of the plan and the principal inspected source files are retained in the audit logs. This identifies the reviewed text without pretending that the repository's existing commit alone describes the working tree.

### 2.2 Focused execution results

| Check | Direct result | What it establishes |
| --- | --- | --- |
| Re-evaluate molecular `HeH9.7` generation `g0004_20260919T103629Z_9b8394a3` with the current binary | Exit 2; cell-1 mass measure `8.335283e-8`, reported threshold `7.5e-9` | P1's refusal is reproducible |
| Same evaluation: molecular partition closure | Reported zero residual in at most three passes | Partition closure alone does not remove the refusal |
| Same evaluation: ionization acceptance report | 504 converged roots; largest reported accepted reaction residual `3.65e-17` | The actual accepted states are much more accurate in this residual measure than the allowed `1e-6` band |
| Extracted production `classify_ending`, only an inner `JFNK done info=0` line | `0\|solved` | The missing outer verdict is not protected |
| Extracted production `classify_ending`, outer `info=0` but no state files | `0\|solved` | The success branch bypasses state-presence checks |
| Public production certification entry writer called with `rep%n=32` | Count stays 32; old entry name remains; its measure changes from 42 to 0 | Silent overwrite is real |
| Current coupled Jacobian test | Five PASS, two FAIL | The two stated failures remain on the recorded binary |
| Production transit selector with `solution` | Selects the uncorrected hydro/species pair | The transit interface is not restricted to `_adv` files |
| Shell environment inheritance and child-script inspection | Override inherited; all 17 shell targets read `EXHALE_TEST_OUT` | The stated current harness failure is not established |

The molecular evaluation and Jacobian probes used fresh temporary directories, `OMP_NUM_THREADS=1`, and `OPENBLAS_NUM_THREADS=1`. Their original case directories were not run in place. The capacity probe linked existing production objects rather than a rewritten approximation of `add_entry`.

No full regression suite, catalog re-solve, golden refresh, transit spectrum calculation, or multithreaded numerical comparison was performed. These were unnecessary for the bounded findings above. In particular, this review has not solved P1's root cause or the five P2 states.

### 2.3 Retained audit material

- [Review driver](audit_20260905/plan_20260919_review_probe.py)
- [Production-linked certification capacity probe](audit_20260905/plan_20260919_cert_capacity.f90)
- [Main probe results](audit_20260905/plan_20260919_review_probe.log)
- [Molecular evaluation log](audit_20260905/plan_20260919_evaluate.log)
- [Boundary trace](audit_20260905/plan_20260919_boundary_trace.txt)
- [Coupled Jacobian results](audit_20260905/plan_20260919_coupled_jacobian.log)
- [Capacity probe results](audit_20260905/plan_20260919_capacity_probe.log)
- [State comparison and catalog selection results](audit_20260905/plan_20260919_state_comparison.log)

The detailed Jacobian samples and staged numerical products remain in `/tmp/exhale_plan_20260919_review_r85j12_o/`; the reusable test code and principal evidence are retained under `docs/audit_20260905/`.

## 3. R1: P1 is real, but P1a is not yet an established remedy

### 3.1 What was reproduced

The saved `g0004` certificate reports a maximum mass-row measure of `3.873e-9` and an accepted hydrodynamic mass entry. A fresh evaluation reports `8.335283e-8` at cell 1 and refuses the state. The old certificate was read, not regenerated by repeating the nine-pass solve.

The fresh evaluation also reports:

```text
ghost H2 partition closed to 0.00000E+00 in at most 3 passes
ioniz-eq acceptance: 504 converged root(s) (max res 3.65E-17)
the loaded composition against the sweep's own root:
  max |d(n_tot+n_e)|/(n_tot+n_e) = 2.527E-12
```

The last two observations matter. `ieq_res_tol = 1e-6` is an allowed acceptance threshold; it is not the residual actually left by this evaluation. The measured maximum is already near floating-point resolution in the normalized reaction measure. This does not prove that every relevant composition component is accurately determined, because conditioning and rate consistency also matter. It does mean that the explanation “the ghost stopped at `1e-6`, so tighten that threshold” is not demonstrated by the current case.

The production sweep distinguishes three controls:

1. `ieq_res_tol`: acceptance of the returned reaction residual.
2. `ieq_inner_tol`: termination of the nonlinear cell solver; the default derives from `sqrt(dpmpar(1))` and the existing diagnostic override is `EXHALE_IEQ_TOL`.
3. `base_ghost_closure_tol`: the imposed molecular partition equation, currently `1e-10`.

These are different conditions. Lowering an acceptance threshold does not necessarily cause the nonlinear solver to take a more accurate step. Increasing an iteration budget does not help if an unchanged stopping test exits first. See [ionization_equilibrium.f90](../src/modules/radiation/ionization_equilibrium.f90), lines 574–591, 1107–1115, and 2696–2740.

### 3.2 The evaluate route is not a byte-preserving measurement of every species

The current path in [EXHALE_main.f90](../src/EXHALE_main.f90), lines 5558–5700, does the following:

1. Load conserved variables and composition.
2. Build a boundary and derive pressure and temperature.
3. Run `ioniz_eq`.
4. Reconstruct pressure and temperature at fixed conserved energy using the refreshed composition.
5. Rebuild the boundary and evaluate the residual.

Thus the evaluation can change eliminated species throughout the column, even though it takes no hydrodynamic step. The code explicitly reports that its work state is one sweep from the loaded state and is not a closed fixed point.

Measured differences between the stored `g0004` state and the fresh evaluation output over physical cells were:

| Quantity | Maximum relative difference |
| --- | ---: |
| Density | 0 |
| Velocity | `2.2194536800371782e-16` |
| Pressure | `4.4942272485243727e-16` |
| Temperature | `2.5276731448484412e-12` |

These are not evidence that an interior error causes the refusal. They establish that the proposed attribution experiment must control more than the two ghost abundance rows.

### 3.3 The missing temperature leg belongs inside P1

The ghost iteration in `ionization_equilibrium.f90`, lines 1865–1930 and 2696–2740, closes the H2 partition and ionization at the cell temperature and radiation context supplied to that loop. Additional partition iterations do not independently close the loop involving the refreshed boundary pressure, particle count, temperature, temperature-dependent rates, and molecular caloric state.

For a compact mathematical description, let `z_g` denote the ghost thermochemical state and write the boundary closure as

\[
G(z_g;U_{\mathrm{phys}},Y_{\mathrm{phys}},J_{\nu},\mathcal B)=0,
\]

where `B` denotes the prescribed reservoir and the accepted boundary model. A residual evaluation then uses

\[
R=R(U_{\mathrm{phys}},Y_{\mathrm{phys}},z_g).
\]

Solving only one subequation of `G` more accurately need not make `z_g` independent of its entry state. The local temperature/rate context can still differ between solve and reload. A deterministic algorithm for a partially closed system is also not evidence that the fully coupled physical boundary equations are satisfied.

**Recommendation:** move P6's temperature/composition closure into P1's diagnosis and acceptance criteria. Do not declare P1 closed and refresh molecular reference products before this issue is resolved.

Also narrow the plan's statement that P1a “moves no interior cell.” A ghost-only update can leave the stored interior array untouched during that call, but it changes the boundary flux seen by the interior equations. A subsequent solve can therefore converge to a different interior state. That consequence must be measured rather than excluded by the edit's location.

### 3.4 Replace the two-run attribution with a controlled sequence

The proposed reuse of the last solved ghost is a useful diagnostic. Agreement in that one experiment would establish sensitivity to the supplied ghost/history, not that the ionization acceptance threshold is the sole cause.

Use the same frozen physical column, reservoir, spectrum, operator settings, executable, and thread settings for the following tests:

| Experiment | Variable changed | Required observation |
| --- | --- | --- |
| Exact replay | None, including rate and closure context | Repeatability of the signed residual and face fluxes |
| Ghost-seed family | Initial ghost composition only | Final ghost state, ionization residual, partition residual, temperature, and mass row |
| Inner-accuracy ladder | Cell-solver stopping accuracy, with acceptance reported separately | Whether the boundary state and mass row converge as accuracy increases |
| Temperature closure | Iteration of thermodynamics and temperature-dependent rates | Whether the restart difference remains after the full boundary closure converges |
| Interior-refresh isolation | Hold or refresh interior eliminated species under a declared diagnostic policy | Separate interior refresh from ghost effects |
| Reload ladder | Written state reloaded and evaluated repeatedly | Distinguish a first refresh, a fixed point, a cycle, and accumulated drift |

At each stage record `p`, `rho`, `T`, particle count, electron count, molecular fractions, charge and elemental budgets, local reaction residuals, the two face mass fluxes, their signed difference, the row scale, and the rounding estimate. Store the equation/context identity with the measurement.

The ordinary four ghost-file variants remain worthwhile, but they mainly prove that ignored file rows are not inputs. They do not substitute for a seed-independence test of the internal ghost solve.

### 3.5 P1a should be an accuracy contract, not “solve to rounding”

Prefer a fully specified boundary solve with:

- Simultaneous acceptance of the thermodynamic, partition, ionization, and inventory conditions at the returned state.
- A declared physical branch and handling of multiple admissible roots, if encountered.
- A measured effect of the remaining closure error on the hydrodynamic boundary residual.
- Explicit failure when this accuracy cannot be achieved; no silent use of an accepted non-root ghost.
- No modification of prescribed reservoir pressure to force agreement with an absolute particle density measured at a different location.

The current molecular ghost exit checks partition closure, while ionization acceptance can in general use the non-root classes allowed elsewhere in the sweep. P1 should explicitly require an admissible chemical root for the boundary it certifies, not just a small imposed-partition residual.

For a locally nonsingular closure, the relevant sensitivity is approximately

\[
\delta R\simeq R_{z_g}\delta z_g,
\qquad
\delta z_g\simeq-G_{z_g}^{-1}G.
\]

This explains why multiplying a scalar reaction tolerance by a gain measured for a particular abundance perturbation is not automatically a valid error bound. The directions, normalizations, and conditioning differ.

### 3.6 P1b is not the next default action

Do not loosen the cell-1 gate to fit the observed restart spread. A tolerance may eventually account for a justified arithmetic/evaluation uncertainty, but only after the operator is consistent and the uncertainty is measured on the relevant states.

Keep three statements separate:

- Discrete mass balance is resolved and satisfies its criterion.
- The state is physically plausible but this residual cannot be resolved to the requested accuracy.
- The residual is resolved and exceeds the criterion.

The second statement should be reported as unresolved, not converted into certification by increasing a gate until it passes. Numerical uncertainty and physical conservation accuracy are related but not interchangeable.

## 4. R2: A concrete molecular diagnostic defect omitted by the plan

At [ionization_equilibrium.f90](../src/modules/radiation/ionization_equilibrium.f90), lines 3405–3423, `ghost_ne_solved` includes H+, He+, He++, and selected metal ions, but omits H2+, H3+, and HeH+. The nearby comment calls their contribution negligible.

That is not the electron-count policy of the current implementation. [utilities.f90](../src/modules/functions/utilities.f90), `calc_ne`, lines 255–308, includes all three molecular ions. The sweep calls that routine with `nmol_eq` at line 3296, so the consistent electron array is already available before the diagnostic is built.

On the fresh HeH9.7 evaluation, direct summation of the output species gives the following molecular share of the electron density:

- Ghost cell -1: `0.010223612332093884`, approximately 1.02%.
- Ghost cell 0: `0.007763590082186174`, approximately 0.78%.

The physical judgment comes first: a diagnostic labeled as the solved electron count should include the charges included by the EOS, regardless of the effect size. The measurement shows that the omission is also not at rounding level relative to the electron count.

**Recommended correction:** use the existing consistent `ne` value, with its units explicitly converted, or call the shared counting routine. Do not maintain a second incomplete formula. Test a molecular-ion-rich composition as well as the neutral limit.

This is a confirmed diagnostic defect, not evidence that the production EOS itself omits these electrons. The inspected consumers of `base_ghost_electron_count` are diagnostic. It nevertheless compromises the measurements proposed to isolate P1.

## 5. P6's alleged 4.7% inconsistency needs a same-location comparison

The plan inherits a comparison between the prescribed reservoir count and a solved ghost count. Those are not automatically the same physical quantity: the ghost can lie below the input level, its density is hydrostatically continued, and the numerical representation uses cell averages.

The current boundary constructs `T_face` using `base_reservoir_nhat` and integrates hydrostatic pressure and density in `base_ghost_averages` ([base_boundary.f90](../src/modules/states/base_boundary.f90), lines 1009–1036). The caloric map uses the installed mixture. That shared use must be checked for consistency, but comparing two absolute particle densities alone does not measure an EOS defect.

For the fresh HeH9.7 evaluation:

| Ghost cell | Radius | Output temperature | `p / (n_particles k_B T) - 1` |
| --- | ---: | ---: | ---: |
| -1 | `0.9998066812221619` | `233.10036275288516 K` | `-2.220446049250313e-16` |
| 0 | `1.0` | `226.02805557716187 K` | `-4.440892098500626e-16` |

The last column uses the full heavy-particle and electron census, excluding a second count of He 2^3S, which is already contained in He I. It confirms algebraic EOS consistency of these written values; it does not establish closure of the chemistry at that same temperature or exact agreement with the intended hydrostatic average.

The earlier D5a memo separately quotes a discrepancy of `2.87e-8` in the particle count divided by mass, not 4.7%. Those two comparisons must not be conflated. The fresh ghost temperatures also differ from the saved ghosts by approximately `7.1e-7` relatively, which is relevant to the temperature leg but not an isolated cause of the mass-row change. This HeH9.7 check is not a repeat of D2b's `carrier_model_a_newton` experiment and does not establish that its reported temperature difference is absent.

**Recommendation:** replace the generic “4.7% inconsistent ghost” claim with a table comparing the reservoir level, face trace, and ghost averages at their own locations. Evaluate the EOS, hydrostatic continuation, molecular caloric relation, and chemistry-temperature consistency separately. Do not force a hydrostatically continued ghost to have the input level temperature solely to remove a percentage discrepancy.

## 6. R3 and R9: Termination and provenance must precede new campaigns

### 6.1 P4 contains a correctness defect, not just record cleanup

[run_case.sh](../LHS1140b/models/run_case.sh), lines 287–337, falls back from a missing outer stationary verdict to the last inner `JFNK`/`PTC` completion line. It then accepts `INFO=0` before checking whether a complete, finite state exists.

The production-function probes reproduced:

```text
inner success, no outer verdict, no state files -> solved
outer success, no state files                 -> solved
inner failure, no outer verdict               -> state_missing
```

The third case also should not be interpreted as a normal completed solve simply because an inner solver returned. These are controlled synthetic logs, not claims that every current catalog record is wrong.

The correction needs more than removing one fallback:

1. Determine the execution route and require that route's terminal event.
2. Keep process termination, solver completion, stored-state completeness, and certification as distinct facts.
3. Require a complete state before declaring a stored solve successful.
4. Bind certification and terminal events to the same run and state identity.
5. Report conflicting evidence explicitly, including final `info=0` with a refusing certificate.
6. Test truncated logs, timeout after a successful inner pass, missing output, nonfinite output, stale output, and normal success.

The publisher also needs attention: [publish_state.py](../LHS1140b/models/publish_state.py), `solver_from_log`, lines 1159–1176, independently returns the last inner `info` when the outer verdict is absent. Repairing only `run_case.sh` leaves inconsistent interpretations in generated manifests.

This does not demonstrate a bypass of every certificate-identity check in the publisher. It does demonstrate that the solve outcome and provenance can be mislabeled, affecting continuation and reporting decisions. Move this part of P4 before P1/P7 campaigns.

### 6.2 New solve generations can lose their parent too

The `g0004` manifest for the molecular HeH9.7 solve has null parent identity and a message saying it was not recorded. Its catalog memo names `g0003` as the seed. This is not one of the old imported evaluate-only records described by P4.

The normal solve publication at `run_case.sh`, lines 1181–1184, passes `${FIRST_GENERATION:-}` as the parent. That records an internal continuation when present, but does not necessarily record an explicitly selected external or existing-generation seed.

**Recommended correction:** carry the resolved seed generation identity from selection through solve publication, alongside the warm-restart policy. Record atomic-to-molecular conversion as a distinct transformation when relevant. Recover old facts into a separate provenance attachment or a new versioned record; do not mutate an immutable state manifest to pretend the information was present originally.

Matching an exact hydro/species pair is useful evidence of state identity, but does not alone identify a unique run: several attempts can write the same state. Preserve ambiguity if run ID, configuration, or time evidence cannot distinguish them.

## 7. R4 and R8: Split P5 into separate contracts

### 7.1 Certification capacity is a production safety issue

`add_entry` returns when `rep%n >= cert_max_entries`, while entry writers continue assigning `rep%e(rep%n)` ([certification.f90](../src/modules/time_step/certification.f90), lines 1897–1931).

The production-linked public-API probe obtained:

```text
entry_count=32
last_name=existing equation sentinel
last_measure=0
silent_overwrite=T
```

This can leave the old equation name attached to a new equation's measure. The artificial full-report input establishes the latent bug; this review did not establish that a normal present configuration reaches capacity.

Fail explicitly before any partial write, or return a checked entry index/status that callers must honor. Increasing 32 to another fixed number does not repair the contract. Test exactly-full and one-more-entry cases. This change is not production physics, but it is production certification logic and should not be described as a shell-harness-only change.

### 7.2 The current output-directory complaint is not established

The inspected `grid_and_gates/run.sh` contains 17 shell-script invocations. Thirteen explicitly forward `EXHALE_TEST_OUT`; four do not:

- `base_level_single_statement.sh`
- `restart_grid_guard.sh`
- `sed_coverage_stop.sh`
- `momentum_row_from_fluxes_only.sh`

All 17 child scripts read `EXHALE_TEST_OUT`. A value supplied in the driver's environment remains exported to its child commands even when the command adds only `EXHALE_EXE`. A direct shell check confirmed that inheritance.

Thus the historical observation may have described an earlier script version, but the plan's assertion that 14 current invocations ignore the requested directory is not supported by this tree. Reproduce an actual escaping write before declaring a remaining defect. Explicit forwarding can still improve clarity, but it should not be reported as a demonstrated physics or isolation fix.

### 7.3 Executable selection should be explicit and testable

The separate `EXHALE_RESID_EXE` override is real in both residual-determinism scripts. Standardize on `EXHALE_EXE`, or document a compatibility policy that refuses conflicting settings. Have the test print and check the actual executable identity.

A test of an override should use a sentinel executable or an unmistakable diagnostic from the selected build. A passing numerical test alone does not prove that it ran the requested executable.

### 7.4 The two Jacobian failures are not one kind of defect

The current `coupled_block_jacobian` probe returned:

| Assertion | Measured result | Interpretation |
| --- | ---: | --- |
| Carrier-direction action against central difference | Mean relative error `0.89540842`; threshold `0.01` | A substantial discrepancy in this test; cause unresolved |
| Energy-direction action against central difference | `0.0056023809`; threshold `0.01` | Pass for this direction and fixture |
| Species-row scale bitwise identity | Three unequal comparisons out of 13 | Printed differences are at the last digits; assess arithmetic error separately |
| Two-cell reconstruction entries | Zero missing out of 13 | Pass |
| Outer ghost rule equality | Zero unequal out of two | Pass |

The carrier sample's largest response in its support is only `4.7885355e-12`, and its moved-row cutoff is `4.7885355e-15`. A large relative error can therefore reflect unresolved finite differences, an inconsistent residual map, or an action implementation error. It is not sufficient evidence to choose among them. Repeat the perturbation ladder and measure the absolute noise floor before modifying the gate or the derivative.

The scale comparison should distinguish an intended algebraic identity from an exact floating-point operation sequence. If two equivalent evaluations necessarily round differently, define and test a defensible arithmetic bound. Do not use this argument to excuse the much larger carrier-action discrepancy.

“It fails on the control” establishes that the latest increment did not introduce the failure. It does not establish that the control is correct. Move the action investigation into the numerical diagnosis relevant to P2, while noting that this coupled fixture is not identical to P2's partitioned hydrodynamic block.

## 8. P2: Worth opening as a bounded investigation, with a stronger design

The proposal to investigate one frozen system rather than increase wall time is reasonable. However, a written hydro/species pair alone may not reproduce the exact stalled Newton operator. Capture or reconstruct and verify the radiation context, composition elimination state, branch decisions, scaling, pseudo-time shift, and reconstruction state used by that iteration.

Before estimating a condition number, establish:

1. The repeated residual is one map at fixed complete inputs.
2. Directional actions have an interval of perturbation lengths where they are resolved.
3. The direction does not cross an active bound or a contact branch unnoticed.
4. Row and unknown scalings are specified.

For a noisy residual, a forward-difference action contains both truncation error proportional to the step and evaluation noise divided by the step. An arbitrarily small perturbation can make the action worse, not more accurate. A condition number of an unresolved finite-difference matrix then describes neither a reliable Jacobian nor the nonlinear problem.

Measure separately:

- The unshifted, scaled Jacobian.
- The pseudo-time-shifted linear system actually solved.
- The preconditioned operator seen by Krylov iteration.
- The true linear residual, the nonlinear trial residual, and the individual refusing equation rows.

Record accepted step lengths, changes in the state and merit, branch changes, limiter/positivity restrictions, and every reason a counter resets. Comparing the x0.01 and x0.10 cases is useful, but their chemical/radiative states differ as well as the Mach number; it is not a one-variable proof of a low-Mach mechanism.

Possible later remedies should follow the measured cause: consistent inner closure, a boundary coupling correction in the preconditioner, a resolved directional-difference policy, or revised globalization. None is authorized by this review. In particular, reopening the held movement bound remains a solver decision.

The five P2 entries should not be described as a proven single failure mode. The evidence listed for `.L22/i3_alt` says it was published but not re-solved, whereas the other entries describe actual failed attempts. Separate unattempted states, boundary-row stagnation, interior-row stagnation, and slow carrier relaxation.

The published Kelley and Keyes analysis formulates pseudo-transient continuation through a shifted Jacobian system and distinguishes the linear solve from timestep control; this supports measuring those mechanisms separately, not assuming a larger iteration budget resolves them. See the local journal version, Section 1.1, Equation (1.3): [Kelley and Keyes (1998)](../../references/Kelley_1998SIAMJNA_35_508.pdf).

## 9. R6 and R7: Correct the remaining P6 scope

### 9.1 CFL: physical control volumes, including boundary-face information

`eval_dt.f90`, lines 30–46, computes the minimum over `1-Ng:N+Ng`. That fact is confirmed. It does not by itself prove an unstable or physically wrong timestep; extra ghost restrictions can be conservative.

The physical stability question concerns the evolved cells and the numerical fluxes at their faces. A boundary Riemann problem can contain a faster wave from the reservoir even when no ghost cell is evolved. Simply taking the same cell-centered formula over `1:N` can remove the only conservative representation of that boundary speed.

Derive the restriction for the actual flux, reconstruction, geometry, and integrator. Conceptually, each physical cell must be constrained by its volume, adjacent face areas, and bounding signal speeds, with the coefficient appropriate to the scheme. Include the external boundary trace in the face-speed estimate.

Also separate physical timestepping from pseudo-time initialization. `eval_dt` supplies the stationary route's initial timestep in `EXHALE_main.f90:5815` and the JFNK pseudo-time floor in `steady_newton.f90:16618`. A change can alter solver behavior, not only the final bits of a marching timestep. The no-step evaluate route does not call that solve-entry `eval_dt`; therefore this change cannot directly explain the already measured P1 evaluate refusal.

### 9.2 Line-search collapse does not always wait for 3000 iterations

The current JFNK implementation already has `n_no_descent_max = 12`, and at `steady_newton.f90:18070–18090` it stops after consecutive unsuccessful searches. An exhausted damped-model attempt can stop earlier.

A run can nevertheless reach the cap if tiny accepted steps repeatedly reset the counter, another control resets it, or a different path remains active. Diagnose which event occurs. Do not add a duplicate “stop after failed line searches” counter without first explaining why the existing one did not fire.

Stopping on stagnation, changing recovery behavior, and changing the returned iterate are distinct decisions. Such work belongs with P2's solver investigation rather than a miscellaneous optimization. As an interface example, PETSc explicitly distinguishes line-search failure, linear-solve failure, and iteration exhaustion in [SNESConvergedReason](https://petsc.org/release/manualpages/SNES/SNESConvergedReason/); this is a useful reporting distinction, not a recommendation to replace the solver.

### 9.3 Stationary spectra are not blocked solely by missing `_adv` products

`EXHALE_transit.py:78–102` already accepts `EXHALE_TRANSIT_STATE=solution`. Its invoked selector, `exhale_transit_lib.py:800–833`, returns `Hydro_ioniz.txt` and `Ion_species.txt` together. The selector was exercised directly in the audit.

Therefore the plan's assertion that the tool reads only `_adv` profiles is outdated. No new spectrum was computed here, so this is verification of the implemented selection path, not validation of a particular He 10830 result.

For transported ionization stages, using the solved composition is also the physically important distinction: an independent post-processing correction can solve a different transport approximation. Do not create `_adv` files merely to satisfy a filename expectation, and do not present a second discretization as the same certified state.

Replace the item with tests of the solution-state spectrum path, metastable population meaning, state-pair consistency, and provenance. Prefer a resolved generation directory or an immutable staged input pair over a moving `output/` directory. There is no need to generate unrelated post-processing products solely for this path.

### 9.4 Charge-exchange turnover guard

The missing guard in `cx_add_to_turnover` is confirmed at `charge_exchange.f90:627–641`. It reads thread-local `cx_kc` and `cx_metal_base` without the temperature check used by other entry points. Relevant callers include reaction-residual normalization and constrained chemical equilibrium.

This review did not reproduce a production call with stale rates. Distinguish that unproven failure from the confirmed missing API guard. Still, an incorrect turnover scale can change chemical acceptance and thus certification, so this is more consequential than a cosmetic diagnostic.

Supply the expected temperature/context or make rates explicit, check initialization, and test alternating cell temperatures on the same thread and multiple threads. Do not assume `threadprivate` alone establishes that a thread's rates belong to the current cell.

### 9.5 Elemental boundary flux diagnostics

`element_base_diffusive_flux` is assigned from `Jf(0)` in `binary_element_diffusion.f90:856`; it represents the flux from the last step's operator coefficients. Printing it after a different evaluation can therefore describe a different state unless its origin is retained.

Report the base advective and diffusive contributions, sign convention, units, and the associated state/evaluation. For a stationary evaluation, use the flux from that evaluation, not an uninitialized or historical step diagnostic. A report-only change should not alter the physical elemental budget. If a golden changes because a diagnostic record is added, identify that separately from a numerical flux change.

## 10. P3, P7, and the catalog counts

Keep the existing retention decision. Neither an old certificate becoming invalid nor a failed solve makes a historical state disposable: those states can be essential controls for P1/P2. No deletion is proposed here.

P7 should not claim that every molecular solve encounters P1. The direct evidence is one successfully solved molecular state refused on re-evaluation; other cases may stop earlier for different reasons. The listed `diffusion_check` cases include atomic cases, so their ability to be evaluated is not logically contingent on a molecular ghost fix.

A reasonable operational policy is to defer expensive molecular re-solves while the shared boundary changes. That is a resource and dependency decision, not proof that every case has the same defect. Cheap diagnostic evaluations and independent atomic checks can remain separate.

The catalog population also needs a reproducible definition. Direct enumeration during this review found:

- 137 `state_index.json` files recursively under `LHS1140b/models`, with 83 non-null `latest_certified` references. This includes older closure rungs and archival directories.
- Using the current `status.py` `GROUPS`/`Case.product_dir` selection: 93 indexed products, with 83 non-null certified references.

These counts do not reproduce the plan's unqualified 88/74. They do not prove that the historical counts were wrong: the plan may refer to a narrower case set or an earlier snapshot. Publish the selected case list, selection rule, timestamp, binary criterion, and treatment of archives/rungs with the counts. A non-null reference is also not a fresh revalidation of its certificate under every current operator setting.

The plan's 18 golden comparisons, 882 parsed inputs, historical EW differences, and longer solve timings were not remeasured. They remain reported results of their cited records, not outcomes of this review.

## 11. Recommended revised order and decisions

| Order | Work | Completion criterion |
| --- | --- | --- |
| 1 | Freeze case membership and source/binary identities; preserve existing products | Reproducible evidence scope |
| 2 | Repair termination classification and the parallel publisher interpretation; record the selected seed | Interrupted or incomplete solves cannot appear completed |
| 3 | Repair certification capacity and the molecular electron diagnostic; verify executable/output isolation | Trustworthy measurement and certification paths |
| 4 | Run P1's controlled attribution, including the temperature/rate closure | Identify which inputs and closure defects move the base row |
| 5 | Implement the selected complete boundary-closure accuracy contract | Same physical state and declared branch give consistent, adequately resolved boundary residuals |
| 6 | Run targeted molecular reload and boundary tests, plus atomic controls that share changed code | Physical and numerical acceptance, with no claim about unperformed tests |
| 7 | If authorized, perform bounded P2 operator/globalization diagnosis | Report cause-specific evidence without changing solver policy |
| 8 | Re-solve affected catalog cases and refresh affected numerical references once | Results tied to the corrected operator and validated state identity |
| 9 | Complete historical provenance recovery and revisit retention only on user direction | Evidence retained and ambiguous history labeled |

The solution-state transit verification can proceed independently. A CFL redesign should be justified and validated on its own affected paths, not inserted as a prerequisite for interpreting a no-step evaluation.

### Recommended decisions for the user

- **P1:** Authorize diagnosis and a consistent ghost thermochemical closure, with P1a reformulated as an error-controlled boundary solve. Do not approve P1b merely to recover the old certificate. Choose any physical closure change explicitly if the investigation reveals that the current constraints are incompatible.
- **P2:** Opening a bounded investigation is worthwhile. Require replay/noise checks before condition-number conclusions, and no solver-control changes without a subsequent decision.
- **P3:** Keep all generations under the existing rule.
- **P4:** Promote termination correctness and new-run parent identity ahead of campaign work; historical record recovery can remain later.
- **P5:** Separate harness interfaces, a latent production certification defect, and an unresolved Jacobian-action issue.
- **P6:** Fold the temperature leg into P1, correct the current-code descriptions, and require an explicit stability argument for CFL changes.
- **P7:** Triage by actual failure mechanism and model dependency instead of imposing a universal molecular-refusal explanation.

## 12. Source and literature limitations

Implementation evidence takes precedence over old measurements and comments. The full plan, relevant source paths, saved generation records, D5a/D5b-2/D2b/D9fix discussions, and the catalog continuation record were inspected to the extent cited above.

For the numerical-method context, the local final journal version of Kelley and Keyes (1998), *SIAM Journal on Numerical Analysis* 35, 508–523, was inspected, including the algorithm and shifted linear system in Section 1.1. An ADS lookup was attempted first for its bibliographic check but failed with `URLError`; no successful ADS verification is claimed. The local published PDF supplied the bibliographic and equation information. PETSc's official nonlinear-solver documentation was consulted for the separation of convergence and failure reasons, not as evidence of EXHALE's implementation.

The accepted D2b model is the starting point of this review. This document does not silently replace it with the earlier Codex preference for B. Any proposed correction must first satisfy the physical and numerical contract of the model the user actually selected.
