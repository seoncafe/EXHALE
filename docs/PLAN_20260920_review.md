# Review of PLAN_20260920

Review date: 2026-09-20 (KST). Scope: the P2, M1, and M2 investigations proposed in `PLAN_20260920.md`, checked against the current implementation, stored generations, diagnostic drivers, and relevant published methods.

## 1. Recommendation

**Proceed with a revised investigation, not with the plan's present decision rules.** Separating reproducibility, derivative accuracy, linear convergence, nonlinear acceptance, and continuation is useful. Keeping the weak-XUV atomic cases separate from molecular composition failures is also appropriate. However, several stated instruments either change the algorithm or do not run on the proposed system, and several conclusions are stronger than their measurements can support.

The most important corrections are:

1. Resolve the identity of the state before interpreting any failing row. The P2 table combines live case-log results with indexes that select different, older generations.
2. Separate observation-only diagnostics from changes to scaling, Krylov termination, projection, and the solver route. Some Stage C diagnostics require the trust-region route and an additional diagnostic gate; the proposed row-attribution tools inspect species rows, not P2's mass and energy rows.
3. Correct the M1 claim that nothing drives the carrier residual. The production outer iteration explicitly relaxes carrier and elemental balances even when they are not Newton unknowns.
4. Replace categorical diagnoses with discriminating tests. A small linear residual does not validate a Jacobian; a missing finite-difference plateau does not uniquely identify noise; failure of parameter continuation does not demonstrate a physical limit.
5. Treat H3+ domain counters as evaluation history, not a map of invalid physical cells. A focused compiled test also reproduced two endpoint-counting defects.
6. Rebase M2 on its existing September 20 evaluation and the actual conversion contract. A converted atomic root is not expected to remain a root of the molecular equations.

The review does **not** establish that any of the remaining states has, or lacks, a stationary physical solution. It identifies defects in the investigation design and in diagnostic accounting, and gives a narrower route to deciding which production changes are justified.

Section 16 adds further recommendations requested after the initial review. Its highest-priority additions are a two-dimensional closure/derivative accuracy experiment and an independent discrete conservation audit. These are proposals, not completed atmosphere tests.

## 2. Evidence and limits of this review

Evidence labels used below:

- **INSPECTED:** current executable statements and their callers.
- **READ:** numerical values already present in a stored certificate or log. Reading them again is not a new atmospheric calculation.
- **EXECUTED:** a check actually run for this review.
- **DERIVED:** a mathematical consequence or a counterexample, not an EXHALE measurement.

The plan was read in full. Principal implementation paths inspected were `src/EXHALE_main.f90`, `src/modules/time_step/steady_newton.f90`, `src/modules/time_step/certification.f90`, `src/modules/files_IO/input_read.f90`, `src/modules/lower_atmosphere/diffusive_photochemistry.f90`, `src/modules/lower_atmosphere/h3p_cooling.f90`, `src/modules/radiation/util_ion_eq.f90`, and `src/modules/init/molecular_seed_from_atomic_state.f90`. The existing `coupled_block_jacobian` diagnostic scripts and linked-object driver were also inspected.

Snapshot identities:

- Plan MD5: `1481dfbb045ceec458868d70c94ab12f`.
- Current checkout HEAD at inspection: `93eed86`.
- `steady_newton.f90` SHA256: `aaddd643ed80c6a1c314aae5bec492e62ac7a025cb6ac4b5b03a59f23bb5f1cc`.
- Production H3+ source MD5: `a0095be5c0a7d30c46453937a45e1ad1`.
- Existing root executable MD5: `f89569405648f9ef8ff8dc78b1ab6810`.
- The inventory used by the plan names catalog binary `75d55d9d4fd0e748cd01d6e35713e34c`.

`make -q` returned 1, and a dry build showed scheduled compilation. Consequently, this review does not assert that the existing root executable is a fresh build of every inspected source. No production rebuild or atmospheric solve was performed. The H3+ test below directly compiled the current source into its own audit directory, avoiding that ambiguity.

New test sources, execution logs, and reproduction instructions are retained in [audit_20260905/plan_20260920_review](audit_20260905/plan_20260920_review/README.md). No case input, stored state, production source, or production executable was changed.

## 3. High-priority correction: identify the state, not just the case

### 3.1 The P2 discrepancy is already partly explained by the records

The plan correctly recognizes that conflicting row locations must be resolved. It does not yet distinguish a generation selected by the index from a later iterate printed in the case's `run.log`.

The following are **READ** from `state_index.json`, the selected generation's `manifest.json`, and its own `certification.txt`; the audit script calculated the component hashes directly.

| Case suffix | Selected latest complete generation | Binary prefix | Stored result |
|---|---|---|---|
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | `g0004_20260919T094145Z_259fe9c3` | `2c3b0acc` | Ending is `hydrodynamic_refusal`, not a wall stop; maximum mass measure 0.6925 at cell 9; tolerance-normalized mass verdict 0.5880 at cell 11; energy 1.949 at cell 7 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | `g0002_20260919T004834Z_d59ec9ff` | `7670f310` | An `evaluate` generation; mass 0.9589 at cell 2; energy 1.112 at cell 1 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `g0002_20260919T004843Z_f7485b14` | `7670f310` | An `evaluate` generation; mass 0.9806 at cell 2; energy 1.071 at cell 1 |

These are not the three states described by the plan's opening table. In particular, the `HeH9.7` scalar case's `run.log:1784` contains energy 0.9978 at cell 248, and `run.log:1807` reports `info=2`, `||R||=0.9978`, and flux spread `4.424e-12`. Its selected generation still contains energy 1.071 at cell 1. The case log also identifies the newer `...ghost_fixed_point_seed_reservoir_row_v3` boundary, whereas the selected old scalar generation identifies `...contact_upwind_v2`.

Thus a case name does not identify one numerical experiment. At least the generation, boundary implementation, and execution phase differ here. The evidence does not support the stronger statement that the same state, measured by the same operator, has contradictory worst cells. Nor does the existence of a log entry establish that its exact iterate was saved and is reloadable.

**Required revision:** insert Stage A0 before the existing Stage A. For each proposed checkpoint, record the immutable state pair and hashes, generation or checkpoint origin, complete input set and spectrum identity, binary/build identity, boundary model, row registry, and the exact log block belonging to that state. If the interesting stalled iterate exists only in a log, report that limitation and decide whether a bounded reproduction is necessary. Do not silently substitute the index-selected parent.

### 3.2 Maximum row measure is not always the certification bottleneck

The mass tolerance depends on the cell. Therefore, `argmax |R_j|/S_j` need not equal `argmax (|R_j|/S_j)/tol_j`. The photochemical generation above explicitly gives different cells for those two questions.

The certified molecular `HeH2.13` generation provides another concrete example: its largest mass measure is `1.369e-9` at cell 1, but the decisive distance is `0.9384` at cell 191, with measure `2.261e-11` and that cell's own tolerance. The certificate states that the tolerance at cell 1 is `7.8e-9`.

Record both locations, their radii, signed residuals, normalizers, local tolerances, and distances. Do not compare a displayed maximum at one cell with a tolerance selected at another. Also, the reported `du` in the cited JFNK output is the flux-spread diagnostic; it is not the norm of the state update requested by Stage D.

### 3.3 What is justified about the newly certified molecular state

The index and certificate support removing `molecular_scalar_gj1132_kzz1e9/HeH2.13` from the unresolved numerical-certification list. Generation `g0008_20260920T053924Z_43df1674` is an evaluation of the named parent, and its report records zero accepted and zero attempted physical steps.

The justified conclusion is that **this reevaluated state needs no additional hydrodynamic solve to satisfy the current certification contract**. The sentence that it "never needed" a solve is broader than that evidence. Certification is also explicitly limited by the wind gate for transport equations and by the validity of the adopted physical closures; it is not an unconditional validation of the entire atmosphere.

## 4. The proposed diagnostic interface needs correction

### 4.1 Options that alter the numerical experiment

The plan's opening promise that A-E change no solver control is not true if all the listed options are enabled as proposed.

| Option | Current behavior, INSPECTED | Consequence for the plan |
|---|---|---|
| `EXHALE_JUDGED_ROWS` | Changes row scaling; reader at `steady_newton.f90:7707` | A controlled scaling variant, not an observation-only switch |
| `EXHALE_JUDGED_MERIT` | Changes the merit used by the solve; reader near line 7718 | Can change acceptance and the trajectory |
| `EXHALE_JV_COLUMN_SCALE`, `EXHALE_JV_PROBE_ARC` | Change the directional probe length; readers at 7653-7673, uses near 10850 and 10915 | Appropriate for a frozen perturbation sweep, but a changed production action if used during a solve |
| `EXHALE_GM_TRUE_RESIDUAL` | Selects the returned Krylov step and termination from sampled true residuals; `pgmres` near 12787 and 13158 | Not just a true-residual monitor |
| `EXHALE_SPECIES_BOUND_HOLD`, `EXHALE_SPECIES_BOUND_PROJECT` | Change active-bound treatment and projection; reader at 7520-7555 | These are algorithm comparisons, not Stage D logging |
| `EXHALE_KRYLOV_TOL_ABS` | Boolean switch tightening the forcing term to the judged tolerance scale, only with judged row scaling and no species rows; 7791-7797 and 17403-17409 | Does not accept a numerical residual-noise floor and is not the proposed general absolute stopping threshold |

For the last option, the code uses

```text
eta_forcing = min(eta_forcing, 1 / max(beta_judged, 1))
```

under its stated guards. It can ask for a **smaller** residual, whereas section 9 proposes it as a way to stop above an observed noise floor. That recommendation could make the oversolving problem worse.

Maintain two experiment groups: an unchanged-route observation and an explicitly named numerical variant. Even a diagnostic that discards its computed steps must restore the complete evaluation context and must not consume the same wall budget without accounting for the additional evaluations.

### 4.2 Several instruments do not inspect P2's default system

The spectral, band-difference, and Krylov-size diagnostics are called inside `trust_region_step` at `steady_newton.f90:15599` onward. Each also requires `elem_diag_here`. That gate comes from `EXHALE_ELEM_DIAG` and selected Newton iterations, including the first two (`16908` onward).

The route selection is explicit at `16661`:

```text
use_tr = (nspec_row > 0) or trust_region_on_the_hydrodynamic_rows
```

For the uncoupled, three-unknown hydrodynamic route, merely setting `EXHALE_PRECOND_SPECTRUM`, `EXHALE_BAND_DIFFERENCE`, or `EXHALE_KRYLOV_SIZE_SCAN` does not reach those calls. Enabling the trust-region route to reach them changes the experiment. A diagnostic run must assert that the requested output was actually produced, rather than interpreting silence as a good result.

There are two additional coverage mismatches:

- `EXHALE_FRONT_ROW` calls `the_binding_species_row`; it returns when `nspec_row <= 0` and inspects slots above the three hydrodynamic variables (`14369` onward). It does not attribute P2's mass or energy row.
- `EXHALE_RESID_ROW_CONDITION` loops over `ie = 1, nspec_row` (`9139` onward). Its displayed cancellation measures concern carrier and element rows, not the proposed hydrodynamic rows. Registering extra unknowns simply to obtain output changes the system under study.

For P2, use the existing hydrodynamic `write_residual_breakdown` output and the hydrodynamic linear-row diagnostic `EXHALE_LINEAR_ROWS` as starting points. The latter evaluates the actual returned direction in the three-unknown solve near `17426`. Extend a frozen operator measurement only where existing entry points cannot expose the needed quantity; do not replace the route merely to activate a monitor.

### 4.3 Existing drivers that the plan should reuse carefully

`EXHALE_RESID_DETERMINISM` already replays the full residual vector, including intervening hydrodynamic and species trials, a discarded evaluation, and a directional probe. Its entry is in `EXHALE_main.f90:1375`, not just in `steady_newton.f90`. This is a better starting point for Stage A than repeated printed maxima.

There is also a concrete executable driver, `src/tests/coupled_block_jacobian/coupled_block_linear_system.f90`, and a build script `run_linear_system.sh` that links the production objects except the main program. It samples a frozen state and directional actions. It is not a missing capability.

However, that driver explicitly calls `set_transported_species_rows(.true.)` at line 143 and carries an old fixed noise threshold of `1.6e-12` near line 95. It is useful for M1's coupled operator, but cannot automatically serve as the unchanged P2 three-unknown control. Select the correct row registry and measure a state-specific floor before reusing its pass/fail rules.

Finally, `EXHALE_RESID_BRANCH_REPORT` is a report on an eliminated-composition seed family (`branch_and_closure_of_the_seed_family`, line 8688 onward), not a general trace of every boundary, limiter, and reconstruction branch. `EXHALE_RESID_ROW_CELLS` selects cells for a report; it does not itself produce all the desired terms.

## 5. Stage A: test reproducibility at the right level

The physical stationary residual is a function of the physical state, fixed parameters, and the chosen closure. A scaled or shifted linear model additionally depends on frozen numerical controls. These should be named separately:

```text
physical residual:        F(Y; physical parameters, closure)
linear model:            Dr^-1 [J(Y) + shift] D
merit/certificate:       separate functions of F, scales, and tolerances
```

Do not let changes in row scaling or pseudo-time shift masquerade as changes to the physical equations. Conversely, if a closure intentionally selects a root using an additional branch label, that label is part of the model definition and must be recorded.

The proposed equality of row norms is too weak: two different residual vectors can have identical maxima. Compare every residual component, report the number of bit differences, the largest absolute and scaled differences, and the affected cells. Repeat after unrelated trial states as well as immediately after the first call. Ensure that a fresh process reads exactly the same serialized state as the in-process control.

The statement that an OpenMP difference *is* a reduction-order effect should be removed. Reduction order is one explanation; shared mutable state, races, initialization differences, and branch-sensitive chemical iteration are other explanations. Fix BLAS threading as well as OpenMP threading. Compare repeated executions at each thread count before comparing counts. A stable rounding difference and a nondeterministic race require different remedies.

Failure of exact replay under genuinely identical inputs is a reason to investigate reproducibility before trusting sensitive derivative estimates. It is not, by itself, proof that the physical model lacks a state function. Also, a discontinuous closure can still be a deterministic function of state: differentiability is the separate Stage B question.

## 6. Stage B and the unresolved Jacobian tests

### 6.1 Missing plateaus do not identify a unique cause

For a smooth residual with evaluation uncertainty, a forward difference has a competing truncation and cancellation error. But a sweep can also fail to show a plateau because it crosses a bound or branch, uses a poorly scaled direction, compares different closures, or samples too narrow an interval. Central differences are not a reference derivative if either side crosses a nonsmooth transition or is clipped differently.

The quoted response `4.79e-12` and relative discrepancy `0.895` in section 13 do not establish a noise floor by their magnitude alone. The relevant quantities include the unperturbed residual, its evaluation error, the actual perturbation, and the row units. These values were not remeasured in this review and remain **READ** historical evidence, not an exoneration of the implementation.

For each direction, record `h`, the actual state displacement after bound treatment, `F(Y+h v)-F(Y)`, an independently measured uncertainty estimate, and the forward/central differences in identical units. Add a Taylor remainder test over feasible directions:

```text
E(alpha) = || Dr^-1 [F(Y + alpha s) - F(Y) - alpha J_ref s] ||.
```

For a sufficiently smooth map and an accurate reference action, the quadratic regime should be visible before the uncertainty floor; an inaccurate reference adds a linear contribution. Measure hydrodynamic, carrier, and element blocks separately. Include the actual solver direction and localized directions at the blocking cells.

### 6.2 Nonadditivity is quantitative evidence, not a universal impossibility proof

At finite step length, the action `(F(Y+h v)-F(Y))/h` is generally not exactly linear in `v`, even for a smooth nonlinear `F`. The useful question is whether the action error is small enough relative to the forcing accuracy and the nonlinear reduction being requested. Significant nonadditivity can invalidate the current Arnoldi model, but any nonzero defect does not imply that no Krylov method can converge.

Measure both additivity and homogeneity, with fixed context and feasible directions, and compare their defects with the achieved linear residual. Do not demand a Krylov reduction below the accuracy with which its operator is represented. PETSc's primary [SNES documentation on inexact Newton methods](https://petsc.org/release/manual/snes/#inexact-newton-like-methods) likewise separates approximate linear solution accuracy from the nonlinear iteration; it is a methodological reference, not a proposed EXHALE dependency.

### 6.3 One existing test has a demonstrably invalid arithmetic assumption

The scale comparison at `steady_newton.f90:6836-6853` uses exact equality between `F0(i)/tscale_code` and a newly evaluated dimensional carrier residual. The production residual was formed by multiplication by `tscale_code` near line 4374. The shell test counts these Boolean comparisons at `src/tests/coupled_block_jacobian/run.sh:188-201`.

Multiplication followed by division need not recover the original floating-point number exactly. **EXECUTED arithmetic example:** `(0.1*3)/3 == 0.1` is false, with difference `1.3877787807814457e-17`. Thus the bitwise assertion can fail even if the dimensional operator was identical. The plan's general concern about the test is justified, but the inspected operation is a unit-conversion round trip, not demonstrated nondeterministic summation.

Proposed repair: compare the dimensional residual before conversion, or compare both sides after the identical conversion operation; if a conversion round trip must be tested, specify its rounding bound. Preserve independent tests of the operator, dimensions, and ghost rule. Do not loosen all derivative assertions because this one arithmetic assertion is too strict.

## 7. Stage C: separate the linear model from its correctness

### 7.1 Small linear residual does not validate the Jacobian

The inference in sections 6 and 9 that a small linear residual proves the linearization right is false. **DERIVED and executed as a toy check:** for `F(x)=x` at `x=1`, take an incorrect approximate Jacobian `J_tilde=-1`. The step `s=1` solves `J_tilde*s+F=0` exactly, while `F(x+s)=2`; every positive step along that direction increases the residual. The issue is not merely an excessively large trust radius.

A small linear residual validates the solve of the supplied model. A Taylor test or an independent action comparison is needed to validate that model against the nonlinear residual. In a pseudo-time system, the shifted and unshifted residuals also differ by the known shift contribution, which must be reported explicitly.

Use this ordering:

1. Verify the action over a resolved interval.
2. Measure the residual of the returned direction against the same frozen action, including the actual shift and scaling.
3. Compare model prediction with nonlinear changes over progressively smaller feasible steps.
4. Only then classify model inaccuracy, insufficient linear reduction, globalization failure, or a changed active branch.

### 7.2 Ritz values alone cannot diagnose conditioning or row scaling

The Jacobian is not established to be normal. For the simple matrix `[[1,K],[0,1]]`, both eigenvalues equal one for every `K`, while its singular-value conditioning becomes poor for large `K`. Conversely, widely spread eigenvalues do not uniquely identify bad row units. A finite Arnoldi projection also produces approximations whose accuracy needs to be checked.

Report Ritz residuals, scaling conventions, relevant singular-value or inverse-norm estimates where practical, and the measured effectiveness of the preconditioner on the actual right-hand side. The species and hydrodynamic compressions are not Schur complements: they omit cross-block feedback. A well-conditioned isolated block can participate in a poorly conditioned coupled system.

If the derivative is accurate and the band is demonstrably missing important coupling, a block or Schur-complement approximation, or a targeted treatment of the omitted radiative coupling, is a reasonable next experiment. It is not justified solely by a wide spectral plot.

### 7.3 A persistent local residual is not proof of incompatible physics

One row can remain large because of an omitted derivative, a bad state, a constrained direction, an inappropriate outer closure, source cancellation, or insufficient spatial resolution. Scaling experiments do not exclude those explanations. Replace the decision-table statement that no step control helps with: "localize the remaining imbalance and test its operator, admissibility, and discretization before deciding whether continuation is the next useful experiment."

## 8. Stage D and E: make the comparisons controlled

A useful accepted-step ledger should name the outer composition pass and inner Newton iteration separately, and contain:

- The state identity, active unknowns, unscaled row measures, and tolerance-normalized distances.
- The actual scaled state displacement, step length or trust radius, predicted reduction, and observed reduction with fixed comparison weights.
- The reduced Krylov residual and independently evaluated residual of its returned step.
- Bound hits, projection displacement, limiter/contact events, and changes to the composition or radiation closure.
- Termination status and counter-reset reasons, distinguishing an inner return, a completed outer route, and external interruption.

The current additivity hook is scheduled for Newton iterations 1, 20, and the configured final iteration; the jump scan is scheduled at the first iteration (`steady_newton.f90:16918`). These are not automatically the one selected outer pass in section 16. Name the checkpoint and invocation policy explicitly, and require a bounded diagnostic exit. An external timeout is a safety ceiling, not a substitute for that exit.

The x0.10 states are useful controls, provided that they are reevaluated with the same operator and configuration conventions as the x0.01 test. The plan correctly rejects attribution to Mach number alone. Also distinguish a certified state, on which the solver may exit before any Krylov work, from a representative nonzero-residual iterate. A preconditioner or step-acceptance comparison needs the latter. A controlled feasible perturbation of a certified state can supply it without pretending that the already-certified root exercises a full solve.

## 9. Stage F: continuation is worth trying, but failure is not a physical result

Adaptive continuation across the factor-of-ten XUV gap is a reasonable conditional experiment. Its initial purpose is to test reachability from a nearby certified branch, not to establish nonexistence. Retain the same base model, composition, grid policy, and definition of the scaled radiation band, and archive the actual spectrum of every point.

Start from a state certified under the measurement binary, change `log(XUV)` in moderate increments, and reduce the increment when correction fails. Use a predictor from previous successful states when available. Record branch markers, the minimum successful increment, residual structure, and the reason each correction ended. If a turning point is suspected, consider an augmented continuation equation rather than assuming a monotonic parameterization remains valid. Arc-length continuation adds the load parameter to the system and constrains the joint update; this is distinct from physical time evolution, as described in the primary [PETSc continuation documentation](https://petsc.org/release/manual/snes/#newton-with-arc-length-continuation).

The second outcome in section 10 should not call refusal a physical limit. It contradicts the caution in section 11. Neither a smoothly growing failing row nor repeated failure below a parameter value excludes an unresolved branch, a fold, a discretization artifact, or a numerical failure.

Recommended result language: **"lowest certified XUV on the tested branch, with the stated grid, boundary law, and numerical method."** A claim about loss of physical stationarity additionally needs a reliable branch/stability investigation and checks of domain and closure validity. An accurately integrated physical-time calculation can test time-dependent behavior; a failed or oscillating pseudo-time iteration cannot establish it. Kelley and Keyes' published PTC algorithm explicitly separates the shifted linear correction from a fully solved implicit time step; see section 13 below.

## 10. M1: distinguish a solved hydrodynamic subsystem from a solved atmosphere

### 10.1 Carrier and element equations are already being driven

Section 14.2 is correct that `Coupled carrier solve: False` leaves these quantities outside the Newton vector. Its next inference is incorrect.

The actual outer route in `EXHALE_main.f90` calls:

- `relax_element_composition` near line 7373 when element diffusion is active, using the face mass flux of the current state.
- `relax_photochemical_composition` near line 7438 when transported rows exist and `block_now` is false.
- Thermodynamic and chemical refreshes around these updates, followed by joint certification in the outer iteration.

The composition is therefore not an unsolved spectator. It is a separate participant in an inexact nonlinear block iteration. A small hydrodynamic norm and a large carrier residual are consistent, but they indicate incomplete joint convergence, not that no algorithm ever addressed the carrier balance.

Call the two `HeH0.083` cases "hydrodynamically within tolerance at the inspected state," not "solved." Their atmospheric stationarity claim still fails the active transport equations. When composition changes, the hydrodynamic rows can become unsatisfied again through the EOS, pressure, rates, and heating/cooling.

### 10.2 A flat sequence does not prove a non-joint fixed point

For an exact, mutually consistent block Gauss-Seidel map, a fixed point satisfies each block equation. A flat maximum residual in this code can instead indicate an inexact inner solve, a movement bound, rejected updates, a cycle, or changes in the cell attaining the maximum. Flat scalar diagnostics do not prove convergence of the full iterate.

Measure both the actual composition displacement and the unconstrained requested displacement, the inner residual and ending, the post-refresh joint residual, and whether successive state vectors converge or cycle. If admissible updates shrink while a nonzero balance remains, identify the restriction responsible before assigning the failure to splitting itself.

### 10.3 `On stall` has a narrower trigger than the plan states

The switch at `EXHALE_main.f90:7877` requires all of the following: the option is enabled, no previous handover, `outer_no_progress`, transported rows present, the route still uncoupled, the current carrier relaxation ending on `carrier_relax_movement_bound`, and enough consecutive such endings (`n_bound_endings >= outer_no_fall_max`).

It does not switch for every flat row history, element refusal, wall stop, or exhausted pass budget. Before assigning an `On stall` trial, inspect whether the selected checkpoint's route can meet those conditions. If not, its failure to switch is expected behavior, not evidence against the coupled method.

`True` is worth a bounded controlled experiment, but not yet a demonstrated answer. It changes the unknown set and normally changes the globalization route to the coupled trust-region method. The unresolved coupled-action tests in section 13 matter directly to this experiment and should precede an expensive coupled campaign on the relevant states.

Running `True` and `On stall` on both `HeH0.083` states means four solves, not one solve for each case; an unchanged `False` control adds two more if fresh controls are needed. A cheaper design is one promising case with a verified coupled action, a bounded `True` comparison, and an `On stall` comparison only if its trigger is relevant. Judge success by the refreshed full certificate, not the hydrodynamic norm alone.

### 10.4 Outer-cell failures warrant a boundary experiment, not a boundary verdict

The outer carrier continuation is explicit: `carrier_outflow_ghost` copies the last physical carrier value into the outer ghost cells (`diffusive_photochemistry.f90:4181-4188`); the associated face rule is defined just above it. This is a concrete boundary closure to test.

Nevertheless, a large residual at cells 499-500 does not prove that the closure is wrong or that the wind is irrelevant. Advection, diffusion, chemical loss, a reconstruction stencil, and the relative normalization can all contribute. A small abundance may also leave a large relative measure whose absolute balance must still be reported.

Use term-by-term carrier balances over the outer interval, including the boundary flux and its derivative. For a radius experiment, keep the inner mesh fixed and add an outer extension instead of rebuilding a different mixed grid everywhere. Compare physical fluxes, inner profiles, integrated carrier budgets, and the radiation column, not just the index of the worst cell. Check the assumed collisional/fluid regime and outward-flow condition before treating a farther boundary as automatically more physical. If an imposed incoming characteristic or an unmodeled incoming carrier flux is needed, changing the boundary model requires a physical decision.

A residual that tracks the outer boundary supports a localization hypothesis. It does not, on its own, establish which boundary law is physically correct or permit those rows to be removed from certification.

### 10.5 H3+ cooling: domain evidence is not yet causal evidence

The current rate used by `util_ion_eq.f90:3078-3113` comes from the actual H3+ module. It implements the LTE emission fit, an H2-collider departure factor, a low-density continuation, temperature clamps, and explicit joins. Miller et al.'s published calculation supplies the LTE fits and the non-LTE table; it does not validate every extension or the optional external-IR closure used here.

In particular, collider density below `1e6 cm^-3` is treated by an explicitly stated linear low-density limit, not by holding the low-density table edge. This is physically distinct from clamping an out-of-range temperature. A count of that limit's use is not automatically evidence of erroneous cooling.

The source counters are incremented on function calls. There is no production call to `h3p_reset_domain_records`; the calls found by the audit are in the focused physics test. Certification copies and sums those counters at `certification.f90:2070`. A cell may contribute to more than one category, be evaluated repeatedly, and be visited during a rejected trial. A net-cooling evaluation can also evaluate the emission at the radiation temperature. Accordingly, 334, 428, and 526 are neither unique affected-cell counts nor a clean measure of the final state's physical invalidity.

There is an immediate counterexample to ranking failure by this total: the newly certified molecular `HeH2.13` report records **1320** activations, including **219** above the fit temperature range. These are **READ** values, not evidence that those extrapolations are physically acceptable. They show that the count alone does not explain numerical refusal.

For the failing energy cells, record `T`, `n(H2)`, `n(H3+)`, the actual H3+ energy term, all competing energy terms, fit/interpolation branch, and directional derivatives. Determine whether the state is outside a validated physical domain and whether a join or clamp actually affects the attempted correction. Co-location alone does not identify the solver's failure mechanism. Conversely, a physically unsupported expression must be repaired or restricted even when its numerical contribution is small; effect size is not the acceptance criterion.

## 11. Confirmed diagnostic defects and proposed repairs

### 11.1 Valid lower endpoints are counted as outside the domain

**INSPECTED and EXECUTED:** `h3p_emission_lte` tests `.not. (T > 30)` and increments `h3p_n_below_fit_T`; `h3p_nonlte_factor` tests `.not. (T > 300)` and increments `h3p_n_outside_nonlte_T`. Equality therefore counts as outside at both valid lower endpoints.

The audit probe compiled the unchanged module with `gfortran -O0 -fopenmp -fcheck=all` and produced:

```text
T=30 K: below-fit count=1
T=300 K: outside-non-LTE count=1
One evaluation: below-collider, outside-T=1 1
Repeated evaluation: below-collider, outside-T=2 2
```

The first two results are an endpoint-classification defect. They do not demonstrate an incorrect rate at the endpoint. The last two establish the overlapping, cumulative semantics used by the plan's statistical argument.

Proposed repair: separate finiteness handling from physical interval membership, count only strict excursions below or above the tabulated domain, and classify nonfinite inputs separately. Add tests at the endpoint and the adjacent representable values. This review did not implement the repair.

### 11.2 Validity of the final state and history of attempted evaluations are conflated

Cumulative evaluation counts are legitimate performance/history data. Copying them into a final-state validity report without that distinction is misleading, especially when A-C deliberately add many discarded probes. Identical final states can acquire different counts simply because they were reached or inspected differently.

Proposed repair: retain a clearly named process-history record, and separately compute a final-state domain map in a controlled evaluation. Distinguish unique cells from evaluations, distinguish analytic-limit use from unsupported continuation, and keep rejected trials separate. Avoid resetting global counters opportunistically inside arbitrary residual calls; the evaluation/report ownership must be explicit.

### 11.3 The Jacobian test's unit-conversion equality is too strict

The arithmetic-contract problem in section 6.3 is a confirmed weakness of the test, not a confirmed explanation of all historical coupled-action failures. Repair the comparison itself while retaining the unresolved derivative discrepancy as an open investigation.

No new production error in the hydrodynamic flux, molecular reaction network, or carrier transport equations was established by a fresh atmosphere experiment in this review. The diagnostic defects above are not presented as the cause of the unresolved winds.

## 12. M2: audit the conversion contract and the present state separately

### 12.1 The September 18 result is not the latest available evaluation

The current index selects `g0005_20260920T053700Z_f59de41d`, evaluated with binary `75d55d9d`, with parent `g0003_20260917T222518Z_c55563e3`. Its stored certificate reports:

- Mass measure `1.041` at cell 2, local tolerance `3.4e-10`.
- Momentum measure `1.450e-3` at cell 1.
- Energy measure `1.021` at cell 1.
- Gated H2 measure `1.000` at cell 432.
- Gated He/H transport measure `1.494e-4` at cell 217.

These are **READ**, not a new solve, and they are not necessarily a freshly converted seed. They should nevertheless appear in the plan before the older `3146d11b` attempt is used to diagnose current behavior. Establish whether each measurement describes the conversion output, a later relaxed state, or an evaluated descendant.

The quoted September 18 mass value `0.2254` divided by `2.2e-11` is `1.0245e10`, about **10.01 decades**, not two decades. This arithmetic was executed by the audit script. That old mass row is at cell 2; the later assertion that its refusal is at cell 1 also needs correction. The underlying `not_solved.md` repeats the two-decade wording, so copying it is not independent verification.

### 12.2 Changing the chemical model changes the equations

The conversion routine explicitly supports two thermodynamic contracts:

```text
invariant=p: preserve rho, v, p; reconstruct T from the new particle count
invariant=T: preserve rho, v, T; reconstruct p from the new particle count
```

See `molecular_seed_from_atomic_state.f90:52-64` and the executable branch at 638. The species transfer checks nucleus and mass budgets, and the converted energy is closed with the production EOS. Neither option claims to preserve the source model's chemical or energy equilibrium.

Consequently, comparing an atomic residual before conversion with a molecular residual afterward changes both the state and the operator. Every added residual is not automatically a conversion error. Even a pressure-preserving map can change thermal rates, cooling, caloric response, boundary closure, and numerical flux response. A stationary atomic certificate is not a molecular certificate.

First test the promised conversion invariants, composition admissibility, charge/nucleus conservation, EOS round trip, level and reservoir metadata, and the base closure reconstructed by the target operator. Then decompose the target molecular residual into physical terms. Classify a failed invariant as a conversion defect, and a correctly converted but off-equilibrium state as an initialization/continuation problem.

The plan's categorical sentence that the seed is the origin of everything is unsupported. In fact, the cited September 18 memo already notes that the nucleus, helium, and mass budgets were exact and distinguishes a valid conversion from failure to reach a molecular wind.

### 12.3 Reuse the actual local-chemistry initialization capability

`EXHALE_MOLECULAR_SEED_X2=local` already uses the cell's thermochemical estimate together with `carrier_h2_chemical_root`; the latter is called near line 554. This is not a full transported stationary solution, but it is more than assigning the base molecular fraction to the whole column. The new plan should state which mode produced each old seed and compare existing modes before proposing an undefined additional relaxation.

"Relax chemistry at fixed hydrodynamics" also needs an explicit thermodynamic invariant. Fixed conserved energy, fixed pressure, and fixed temperature are different experiments when molecular abundances and heat capacities change. Use the production EOS throughout and report which quantity is allowed to move.

Refusing a metal-free scalar restart as though it were an identical metal-bearing photochemical state is appropriate. An explicit conservative remapping could create a new initial guess, but it would not be a compatible warm restart or inherit the source certificate. Any such conversion should be named, validated, and approved as a separate experiment. It is unnecessary to settle every M1 question before conducting M2's read-only state and conversion audit.

## 13. Physical and numerical references checked

The local final journal versions were read with `pdftotext -layout`; the relevant methodological text, rather than abstracts alone, was inspected.

1. **Miller, Stallard, Tennyson, and Melin (2013), "Cooling by H3+ Emission," J. Phys. Chem. A 117, 9770-9777.** Local file: `references/Miller_2013_JPCA_117_9770.pdf` at the workspace root. Tables 5 and 6 and the discussion of equations 8-10 distinguish LTE emission from the H2-dependent departure factor. The non-LTE calculation uses the proton-hopping rate and explicitly discusses its uncertainty and possible upper-limit character at low density. The low-density continuation, interpolation joins, and out-of-range temperature handling must therefore be identified as code choices, not assigned to the paper. [Published DOI](https://doi.org/10.1021/jp312468b).
2. **Kelley and Keyes (1998), "Convergence Analysis of Pseudo-Transient Continuation," SIAM J. Numer. Anal. 35, 508-523.** Local file: `references/Kelley_1998SIAMJNA_35_508.pdf`. Section 1, especially Algorithm 1.1 and equations 1.2-1.4, defines the shifted correction and its relation to an implicit step; section 1.2 also warns that small steps need not imply small residuals. This supports separating the physical residual, pseudo-time model, and state-update diagnostic. It does not certify any particular EXHALE stopping rule. [Published journal copy](https://www.jstor.org/stable/2587140).
3. **PETSc primary SNES documentation**, consulted for the distinction between an inexact Newton solve, Jacobian verification, and arc-length continuation. Its methods are examples of the relevant numerical distinctions, not a recommendation to replace EXHALE with a library in this task. [SNES manual](https://petsc.org/release/manual/snes/).

ADS was attempted first for the citation check, but the shell request failed at DNS resolution. The publisher DOI request for Miller through the web tool was also unavailable. Bibliographic and content checks for the papers therefore relied on the local published journal PDFs, not on an arXiv version or an unverified search summary. No new external literature result is claimed from the failed ADS request.

## 14. Revised execution order and decision gates

| Phase | Minimum useful work | Gate before proceeding |
|---|---|---|
| A0: identity | Select immutable checkpoints, log phases, input/SED hashes, binary and boundary identities; include M2's September 20 evaluation | No mixed state/measurement rows in the case table |
| A: repeatability | Full-vector immediate and intervening-state replay; fresh-process replay; controlled OpenMP/BLAS comparison | Reproducible operator or a quantified and explained evaluation uncertainty |
| B: action accuracy | Feasible perturbation ladders, Taylor remainder, additivity/homogeneity, actual active-set and branch data | Resolved directional information on the rows being investigated; unresolved tests remain open |
| M1 local audit | Fixed-wind carrier/element outcomes and unaccepted requested movement; outer boundary terms; state-specific cooling-domain map | Identify whether movement limits, inner accuracy, outer closure, or a physical-domain problem is actually present |
| C: frozen linear comparison | Correct route and row registry; shifted/unshifted residuals; baseline action, band, scaling, and actual RHS | Separate action error from linear-solve/preconditioner failure |
| D/E: bounded comparisons | Accepted-step ledger; matched nearby control with a meaningful nonzero residual | A reproducible mechanism, not just a correlation |
| M1 coupled trial | One justified `True` comparison first; `On stall` only with a relevant trigger; preserve unchanged control | Full refreshed certification or a localized, recorded failure |
| M2 seed work | Verify conversion contracts and available local initialization; define the thermodynamic invariant | An admissible target-model initial state with traceable residual terms |
| F: continuation | Adaptive XUV path, with separate authorization and immutable results | Report reachable certified branch, not unproved physical nonexistence |

Predeclare resource limits as residual-evaluation, Krylov-product, inner-iteration, and outer-pass limits in addition to a wall ceiling. Estimate time from one measured frozen evaluation on the actual checkpoint. The plan's unmeasured "minutes" estimates should not be promises, especially when eliminated chemistry and spectral diagnostics dominate the cost.

The requested decisions should be narrow:

- Approve observation-only A0/A and state-local audits first.
- Permit numerical variants only as labeled, frozen or bounded comparisons; do not alter the baseline implicitly through diagnostic flags.
- Defer a broad coupled solve campaign until its relevant action has passed the accuracy checks.
- Keep changes to the physical boundary law, cooling closure, certification gate, and continuation campaign separate from solver instrumentation.

Do not relax a tolerance or exclude a failing equation merely to obtain a certificate. If a physical approximation is unsupported, define and implement the physically justified model or restrict its stated applicability; if the difficulty is arithmetic or discretization, measure and address that difficulty explicitly.

## 15. Validation performed and work left open

**Executed:** the generation/manifest/certificate audit and component hashes; the M2 tolerance-ratio calculation; two small arithmetic counterexamples; direct compilation and execution of the unchanged H3+ module to test endpoint and cumulative counters. The corrected H3+ build and probe exited successfully. All new probe code and logs are retained in the audit directory linked above.

**Not executed:** fresh residual evaluations of the P2/M1/M2 atmospheres, the full coupled Jacobian suite, OpenMP reproducibility measurements, a new atmosphere solve, radius/grid experiments, or XUV continuation. This review used targeted checks rather than carrying out the proposed campaign; the state/binary and diagnostic-route ambiguities should be resolved before broader measurements are interpreted. No physical-model changes were requested or made. No full regression suite was run because production behavior was not changed.

**Bottom line:** the investigation is worth doing. Its current strongest claims are not yet earned by the evidence. Correcting state identity, diagnostic reachability, M1's interpretation, and the causal decision rules will make a short targeted investigation substantially more informative than another broad solver campaign.

## 16. Additional recommendations

Added on 2026-09-20 after the user's request for further suggestions. The paths below were inspected for this addition; no new atmosphere calculation or production change was made.

### 16.1 Measure closure bias separately from repeatability and derivative error

**Priority: immediate, before interpreting a successful Stage B test.** A deterministic residual can be an inaccurate approximation to the intended fully closed residual. Replaying it bit for bit does not measure that bias.

The distinction exists explicitly in the implementation. `eval_residual` repeats the eliminated chemistry/radiation update and allows either a convergence test or a prescribed number of sweeps (`steady_newton.f90:4175-4278`). The Newton loop separately forms `F_jac` with `n_eq_sweeps_fixed=n_eq_sweeps_model` near line 17037. Directional probes use the fixed-count interface as well. This is sensible for comparing the same finite-iteration map, but does not prove that its derivative is already that of the converged closure.

Add a two-dimensional experiment at one frozen physical state:

1. Hold the seed, row/column scales, reconstruction, and branch fixed. Increase the closure sweep count, for example through a short doubling sequence, and save `F_k(Y)`, the closure increments, the relevant chemical balance, and the temperature/heat/cooling changes. Stop when the available evidence resolves the intended accuracy, not merely when two printed norms coincide.
2. For each count, repeat the directional perturbation ladder. Compare `J_k v` both across perturbation lengths and across closure counts. Within each difference, use the same count and seed at both endpoints.
3. Compare the production convergence-controlled residual with the fixed-count map used by the production Jacobian. Record `||Dr^-1(F - F_jac)||` and its row distribution at the actual iterate.

Two distinct outcomes matter. If the derivative stabilizes with perturbation length but moves with closure count, the finite-difference rule is not the primary issue: the inner map is not yet resolved. If both stabilize while the linear solve still fails, stronger preconditioning becomes a better-founded candidate. A change in the selected chemical branch is a separate outcome, not closure error to average away.

For a forward difference in consistent scaled coordinates, a useful error model is `O(h) + O(epsilon_F/h)`, where `epsilon_F` is uncertainty in evaluating the relevant map, not automatically machine epsilon. Deterministic closure bias also requires its own convergence check; it is not estimated by immediate replay. PETSc's [matrix-free function-error interface](https://petsc.org/release/manualpages/Mat/MatMFFDSetFunctionError/) explicitly relates its perturbation control to estimated function-evaluation error. This supports measuring that error rather than assuming that `sqrt(machine epsilon)` is always the correct physical probe scale.

Only after this experiment should an adaptive inner-accuracy policy be considered. Such a policy should keep the error induced by the inner closure below the accuracy required by the outer step and final certificate. Raising every inner accuracy globally would be expensive and would not repair a branch-selection problem.

### 16.2 Add an independent discrete conservation audit

**Priority: immediate for the energy failures and the outer carrier rows.** The question is not just whether an individual normalized residual is small, but whether the discrete equations transfer conserved quantities consistently across the entire domain.

The necessary quantities already have concrete owners:

- `RK_integration` stores the numerical `face_flux` assembled by `RK_rhs`.
- `steady_residual.f90:267-277` constructs mass, momentum, and energy residuals from flux divergence, source terms, radiative heating/cooling, and optional viscous/conductive sources.
- `carrier_steady_residual` exposes dimensional residual and term arrays (`diffusive_photochemistry.f90:3630` onward). The carrier row dump also records the two advective faces, the two diffusive faces, and chemical production/loss separately near lines 5854-5868.

For each conserved quantity, form the finite-volume identity from the **actual assembly**, schematically

```text
sum_j V_j R_j
  = A_outer F_outer - A_base F_base - sum_j V_j S_j.
```

Use the cell volumes, face areas, units, and source signs of that operator. If an exported quantity already includes area or inverse volume, undo or retain that factor consistently rather than applying it twice. Do not replace face fluxes with a cell-centered product merely because both approach the same continuum expression.

Report the signed identity defect, the sum of absolute local defects, and the largest local defect. Signed cancellation alone can hide large opposing errors. Repeat the budget over subdomains ending near the failing cells; the flux crossing the certification gate remains part of the physical budget even when equations below it do not gate a certificate.

For elemental budgets, weight the species by their nucleus counts. Chemical source terms should cancel for the appropriate conserved element after all participating species are included; H2 alone is not conserved. A source-free individual-species flux is constant only in a deliberately inert test, not through a dissociation front.

For energy, first state the implemented energy convention. The hydrodynamic physical flux uses `v*(E+p)` (`Num_Fluxes.f90:538-543`), while the residual carries other energy sources separately. The caloric module describes bound H2 rotational/vibrational energy and places dissociation energy in chemical source accounting rather than in that heat capacity. Therefore, do not append a chemical binding-energy flux or gravitational-energy correction by analogy without deriving its relation to the implemented variables. Audit reaction heating, radiative losses, and any energy transported by diffusing material against one explicit convention. This is a proposed accounting check, not a finding that such a term is currently missing.

The budget should be reconstructed independently from exported face and source data, then compared with the summed residual. Summing the residual through the same reporting function twice is not an independent check of assembly or units.

### 16.3 Use a small diagnostic hierarchy to localize the responsible coupling

**Priority: before another full campaign.** A new method should first be tested on the smallest operator that still reproduces the defect. Use existing boundary, carrier, and caloric test entry points where they match that operator; do not create a second implementation of the physics for convenience.

| Controlled problem | What remains active | Main question |
|---|---|---|
| Compatible hydrostatic state | Actual gravity, EOS, reconstruction, and boundary law; a prescribed consistent source balance | Does the weak-flow residual and its derivative respect the intended equilibrium? |
| Inert carrier in a prescribed wind | Actual carrier advection/diffusion and face rules; chemical sources disabled only for this diagnostic problem | Do the boundary flux and discrete transport budget agree? |
| Local reacting parcel | Actual reaction network and EOS with a stated thermodynamic invariant; no spatial transport | Are the local chemical root, energy accounting, and directional derivatives mutually consistent? |
| Coupled physical checkpoint | The full original model on a selected state | Which discrepancy appears when the previously checked components interact? |

These are proposed validation problems, not substitutes for the original atmosphere. A hydrostatic test must use a compatible reservoir and pressure profile; otherwise it tests an imposed mismatch rather than well balancing. An inert test must state how its sources were changed and must never be published as a physical atmospheric solution. A coarse grid can be useful for direct linear algebra, but only if it retains the failing branch/front; otherwise it does not reproduce the original defect.

There are already relevant test directories, including `src/tests/carrier_outer_boundary`, `src/tests/carrier_boundary_jacobian`, `src/tests/carrier_returned_state_acceptance`, and `src/tests/h2_level_ladder`. Their existence is not a claim that they all implement this hierarchy or pass on the current build. Select and inspect the needed caller before running or extending one.

### 16.4 Check whether eliminating a chemical variable is physically justified

**Priority: high for M1/M2, but distinct from choosing coupled versus alternating transport solves.** Moving H2 into the Newton vector does not automatically validate the instantaneous local equilibrium used for the remaining species. A perfectly converged stationary solver can still solve an inappropriate closure.

The code already contains `chemical_decay_against_the_transport_times` (`steady_newton.f90:9177-9252`). It reads local molecular reaction rates and compares them with transport over a resolved abundance-gradient length. This is a useful starting point and should be included at the molecular front, the outer failing cells, and representative interior cells.

However, interpret its output carefully:

- An undefined advective or diffusive time is printed as zero. In a stagnant cell, zero advection does not mean an instantaneous transport time. New reports should mark the mathematical infinite-time limit or unavailable coefficient explicitly, rather than treating the printed zero as a measured timescale.
- `measure_molecular_decay_rates` forms a reaction Jacobian and then uses `abs(wr(i))` to summarize its eigenvalues (`ionization_equilibrium.f90:4896`). The reported magnitudes do not distinguish relaxation from growth, and exactly zero real parts are skipped. This is an inspected limitation of the diagnostic summary, **not evidence that an unstable or neutral mode occurs in these states**.
- A diagonal decay estimate for one species does not replace the slow modes of the constrained coupled network. The supplied calculation also holds a local background; radiation feedback and thermal coupling can change the full system's modes.

Before using fast chemistry as a reason to eliminate a variable, retain the signed real parts of the physical kinetic eigenvalues, check eigenvalue accuracy and near-neutral modes, and compare the relevant stable relaxation times with the physical transport and thermal times. Use the sign convention of `dot x = P-L`, not that of an oppositely signed stationary residual. If a slow or growing mode is found reliably, reconsider the physical closure and which variables are transported. Merely increasing Newton effort cannot correct the wrong set of equations.

### 16.5 A physics-based derivative and preconditioner are conditional development options

**Priority: only after the action/closure experiments establish a resolved map.** If finite differencing the nested chemistry remains costly or inaccurate, there is a more specific direction than repeatedly increasing the Krylov dimension.

Let `Y` be the retained hydrodynamic and transported variables, and `z` the variables determined by a closure `C(Y,z)=0`. For a differentiable, locally unique branch with invertible `C_z`, implicit differentiation gives

```text
z_Y = -C_z^-1 C_Y
dF/dY = F_Y - F_z C_z^-1 C_Y.
```

This is a mathematical proposal, not a statement that these blocks are already assembled by EXHALE. It shows why differentiating a source while holding its eliminated composition fixed omits a physical feedback. A useful initial implementation experiment would compare this sensitivity with resolved finite differences on a small local system, including its EOS and energy terms.

For the fully coupled hydrodynamic/transport system, write its Jacobian as `[[A,B],[C,D]]`. A composition-block preconditioner that approximates the feedback in `D-C A^-1 B` can be more informative than treating the isolated blocks as unrelated matrices. The [primary PETSc block-preconditioning documentation](https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/) gives this Schur-complement structure. Citing it does not imply that a PETSc migration is needed.

The applicability limits matter. Optical-depth dependence couples cells, so a local chemical derivative is only part of the full action. Near a chemical branch fold, `C_z` may be poorly conditioned and elimination may cease to be a useful formulation. A changing limiter or active bound can invalidate a smooth derivative model. These are reasons to measure the omitted feedback and possibly retain more physical unknowns, not to declare an approximate local derivative exact.

Adopt any such method only if it improves the independently checked linear residual and nonlinear prediction at matched checkpoints, then reaches the unchanged final certification requirements. A smaller Krylov iteration count by itself is insufficient, especially if a more expensive closure is hidden inside each product.

### 16.6 Separate equation accuracy, physical validity, and scientific output accuracy

**Priority: part of the final acceptance design.** Keep three different claims separate:

1. The computed state satisfies the implemented discrete equations.
2. The chosen equations and closures are applicable in the state being modeled.
3. The scientific quantities of interest are resolved under changes in grid/domain and defensible physical uncertainty.

A numerical certificate addresses the first under its stated gate; it cannot establish the other two. Conversely, apparently stable mass loss or a nearly unchanged line profile does not excuse a failed conservation equation or an invalid closure.

After a case is numerically certified, compare quantities such as mass flux, location of the molecular transition, relevant species columns, and integrated heating/cooling on a small, justified refinement study. Use conservative restriction/prolongation when comparing cell averages on different grids. For a radius study, retain the inner mesh as proposed in section 10.4. If a spectral observable is the eventual target, propagate only the selected certified comparison states into that calculation; do not regenerate all catalog products during solver diagnosis.

Predeclare the accuracy needed for the intended scientific conclusion. Report numerical changes and physical-model uncertainties separately. Neither should be converted into an unmotivated relaxation of the equation tolerances.

### 16.7 The next experiment should answer one decision

My preferred first addition is **one frozen checkpoint with a two-dimensional closure/probe sweep and a discrete budget report**, not simultaneous changes to several solver controls. Choose an energy-failing molecular checkpoint if the immediate goal is M1; choose an identity-resolved weak-XUV checkpoint if the immediate goal is P2. Use a matched certified control only for the same measurement, not as proof that the failing model is equivalent.

Before the run, state which result would lead to which action:

- Closure-count dependence: resolve or reformulate the inner closure before optimizing Krylov work.
- A face/source accounting mismatch: investigate the discrete operator or units before changing the step method.
- Accurate closure and action, but poor returned linear residual: test a targeted preconditioner change.
- Accurate local prediction at small feasible steps, but poor acceptance of larger steps: test globalization or the constraint treatment.
- A physically invalid equilibrium assumption: revise the physical model, even if its effect appears small.

Archive one compact evidence package containing the checkpoint identities, full configuration, actual enabled diagnostics, requested and realized perturbations, signed row terms, and explicit exit status. Keep diagnostic runtime separate from ordinary solver cost. Treat "no diagnostic output," "timeout before the checkpoint," and "measurement unresolved" as distinct results; none is a numerical pass.

These additions extend the review's proposed experiments. They do not add new completed validation claims to section 15, and they do not authorize changes to the boundary law, chemical network, or published case states.
