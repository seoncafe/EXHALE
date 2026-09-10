# Critical review of PLAN_20260909

Date: 2026-09-09, KST

## 1. Recommendation

**Accept the overall direction, but revise the plan before authorizing its execution.** The division between numerical repair, constrained chemistry, stationary solutions, restart/observable handling, and physical integration is sensible. The refusal to excuse a nonzero balance row at a carrier bound is particularly important and should remain unchanged.

The current document is not yet a sufficiently precise execution contract. Its main weaknesses are:

- N16 omits an existing production-level defect: the computed integration error is not connected to step rejection.
- The step-doubling path can compare states at different times after a retry and can advance the clock by a different interval from the retained state.
- Several proposed tests target the wrong layer or are algebraic demonstrations rather than tests of the production decision being repaired.
- The constraint map does not yet specify consistently which species are fixed, rescaled, or eliminated in each calculation.
- Stage B's written gate omits completion of N4b, while N8's final tolerance study depends on solutions that the stage ordering has not yet produced.
- Universal bitwise requirements and repeated planetary runs conflict with affected-path validation and physical correctness as the acceptance criterion.
- Separate build directories do not isolate simultaneous edits to the same source file or API.
- Several decisions unnecessarily block an entire stage instead of the specific interface change that requires approval.

I recommend authorizing a bounded first increment only after these points are incorporated: N0's revised tests, N3/N1/N2 under coordinated source ownership, and the first two iterations of N7. In parallel, the physical-step defects should be diagnosed and repaired before any new claim of time-dependent accuracy. Full planetary reruns, broad unknown-layout changes, and schema changes remain separate decisions.

This review is not permission to lift the code pause. Only this review file was created; production source, running calculations, reference outputs, and the original plan were not modified.

## 2. Evidence and inspection scope

The plan was read completely, including its item list, dependency diagram, test matrix, execution rules, and decisions 13-18. Its statements were compared with current implementations, not accepted solely because they appeared in the previous review.

Evidence labels used below:

- **SOURCE:** established from the current code and the relevant caller or consumer.
- **REPORTED:** a result stated in the plan, issue inventory, update log, or earlier review; not remeasured in this turn.
- **DERIVED:** an analytical consequence of the equations or inspected control flow.
- **PROPOSED:** a correction, decision, or test that has not been implemented here.

No numerical test, new build, planetary integration, or regression refresh was performed in this turn. In particular, the earlier audit harness was inspected but not rerun. Its previous measurements remain REPORTED here.

Principal inspected paths include:

| Subject | Implementation and consumer |
| --- | --- |
| Stationary unknown registry and trust region | `steady_newton.f90`: `set_transported_species_rows`, `solve_steady_jfnk`, `trust_region_step`, `dogleg_step`, `pgmres` |
| Elemental budgets | `species_table.f90`; `diffusive_photochemistry.f90`: `carrier_state`, `carrier_source`, `limit_to_element_budget`, carrier write-back |
| Physical error estimation and adoption | `EXHALE_main.f90`: attempted-step loop and `err_pass`; `attempted_step.f90`: error estimate, reduced-step policy, checkpoint handling |
| Restart | `load_IC.f90`, `write_output.f90`, and the main-program restart/steady entry |
| Post-process and transit | `T_equation.f90`, `post_process_adv.f90`, `EXHALE_transit.py`, `exhale_transit_lib.py`, and selected spectrum readers |
| Tests and build routing | Makefile; `steady_species_rows`, `steady_selfconsistent_residual`, and `attempted_step` test sources; the existing audit harness |

The local `../references/` collection was checked. The final ApJ identity of the local Huang paper was confirmed from its PDF header. No new literature claim or rate recommendation is introduced here; the plan's shell-physics note does not require reopening the rate audit. For numerical API distinctions, current official [KINSOL documentation](https://sundials.readthedocs.io/en/latest/kinsol/Mathematics_link.html) and the [PETSc reduced-space variational-inequality solver description](https://petsc.org/release/manualpages/SNES/SNESVINEWTONRSLS/) were consulted. Neither library's existence proves that EXHALE already implements the required method.

### 2.1 Identity of the reviewed files

HEAD: `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`. The worktree contains extensive existing changes, so HEAD is not an adequate source identity by itself.

| File | SHA-256 |
| --- | --- |
| `docs/PLAN_20260909.md` | `d88b225d02e060a1cc56791878c453978e5d3c83061e20652a4b4fb6f6996239` |
| `src/modules/time_step/steady_newton.f90` | `0678dca1af85b9de93f56662e12348c382554f961837501f52a54dd1e0984e28` |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | `04c2e7e74254897ce3d0deb2ba6a2d52e7b4d8416f6652dd9ab11ed94cbdc8bf` |
| `src/EXHALE_main.f90` | `71d20728d256854d95ef3cbe486193a743fde56b8b6355eda5513cbf5443809a` |
| `src/modules/time_step/attempted_step.f90` | `a9137543ffef21bd555ef7e23ff75f159553c695e71402d36a22113699d8f725` |
| `docs/ISSUES_20260909.md` | `e4d737465c49392474aa35c1a0e21255eed4c9376ab4258202b02c817908187d` |

The issue inventory now includes completed B5l results. Its reported carrier runs with GMRES sizes 40, 80, and 160 should not be treated as experiments that have never been made. A corresponding atomic-element experiment in N7 is a different question and may still be useful.

The named `build_lwv/` and `EXHALE_lwv.x` exist and have the reported September 9 01:04 timestamps. This confirms the paths and timestamps, not exact source equivalence. Before an execution gate, identify the actual binary by hash, compiler flags, linked libraries, and a source manifest.

The paths `scratchpad/plan_rev2_brief_common.md` and `scratchpad/atomic_elem/` were not found under the inspected workspace. They may refer to an external working location; this review cannot verify their contents. The plan must provide resolvable paths or durable fixtures before making them prerequisites.

## 3. Important additional findings: N16 is not just a species-norm extension

### F1. Integration error is calculated but does not reject a production step

**Severity: critical for physical-time accuracy. Evidence: SOURCE.**

In `EXHALE_main.f90`, approximately lines 2806-2824, the sampled step-doubling path:

1. calls `attempted_step_error_estimate`;
2. prints `e_err`;
3. restores the stored full-step state;
4. proceeds to the adoption boundary.

There is no comparison of `e_err` against one on that path. A source-wide search shows the distinction clearly:

- `as_reject_int_error` is declared and imported;
- `attempted_step_reduced_dt` contains a controller formula for that reason;
- unit tests call the controller formula directly;
- the production caller does not issue that rejection reason from the newly calculated error.

Thus, a valid routine and a passing policy test exist, but the error-controlled adoption path is not connected. The sampled estimate is currently a diagnostic even for the hydrodynamic variables it measures.

This finding strengthens the previous review's R9. That review identified species coverage and sampling frequency, but did not explicitly identify the missing production rejection connection. The new plan should correct that omission rather than simply repeat the earlier wording.

Required N16 amendment:

- Add an end-to-end gate between error estimation and physical adoption.
- Reject when the finite normalized estimate exceeds its acceptance limit; reject or stop on a nonfinite estimate.
- Restore the complete state at the beginning of the macrostep, reduce the proposed interval, and recompute both comparison trajectories.
- Keep physical clocks, accepted counters, and physical energy/source accumulations unchanged during rejected attempts and diagnostic branches.
- Test the main-program route, not only the function that computes a smaller `dt`.

A test should force an excessive temporal error while all local nonlinear solves and positivity checks succeed. The expected result is rejection, unchanged physical time, a shorter retry, and eventual acceptance only after the error requirement is met.

### F2. Retried comparison branches can cover different intervals and corrupt the clock

**Severity: critical when step doubling and a shortened retry occur together. Evidence: SOURCE and DERIVED control-flow consequence; not replayed in a full run.**

The main program saves the proposed `dt_step` before entering `err_pass`. An outer retry can reduce `dt`. However:

- the half-step branch is initialized from `0.5*dt_step`, not necessarily half the interval actually covered by the accepted full-step branch;
- the third pass continues with the current `dt` rather than independently certifying coverage of the remaining target interval;
- after restoring `as_full`, the code resets `dt = dt_step`;
- the physical clock then advances by that reset value.

Relevant source locations are approximately lines 1706-1734, 2797-2799, 2817-2820, and 2837 of `EXHALE_main.f90`.

A conditional example follows directly from those assignments:

```text
Proposed macrostep: H
Full branch: rejected at H, then accepted after H/2
Comparison branches: H/2 + H/2, assuming neither is rejected
Compared endpoints: t + H/2 versus t + H
Retained state: full branch at t + H/2
Recorded clock increment after restoration: H
```

The difference is no longer a temporal discretization estimate at a common time, and the retained state is labeled with the wrong elapsed interval. Shortening a half-step branch creates related coverage problems.

Recommended design: the error controller owns a proposed macro-interval. Each comparison branch must either reach that exact endpoint or return a named failure. If a branch cannot do so under the chosen algorithm, reject the complete macrostep and restart from its common beginning with a shorter interval. Do not silently compare partially completed intervals.

Required tests include rejection in the full branch, rejection in the first half, rejection in the second half, and exhaustion in each branch. Check the endpoint, restored state, accepted interval, clock, and physical accumulations together. A correct checkpoint of a substep is not automatically the correct checkpoint of the complete error-control transaction.

### F3. State selection and estimator normalization need an explicit contract

**Evidence: DERIVED; proposed numerical design requirement.**

The code retains the full-step trajectory. For a method with leading local error `C H^(p+1)`, the full-versus-two-half-step difference is, to leading order,

```text
Delta = C H^(p+1) (1 - 2^(-p)).
```

Therefore, the error estimate for the retained full step is `Delta/(1 - 2^(-p))`; for the retained two-half-step result it is `Delta/(2^p - 1)`. These estimates refer to different accepted states.

For a first-order complete split update, the full-step correction factor is two. A raw difference is not automatically the error estimate of the state the code retains. The order must be established for the complete update, not inherited from the RK hydrodynamic substep alone.

N16 should state which trajectory is accepted, its effective order, the estimator convention, and the scales for every transported quantity. It should also require inner nonlinear error to be below the temporal-error budget. These requirements do not justify increasing an unvalidated physical tolerance.

Recommendation: split N16 into:

- **N16a:** common endpoint, rollback, clock, and actual error-driven rejection;
- **N16b:** species-aware norm, accepted-state normalization, and temporal convergence;
- **N16c:** measured performance and any proposed intermittent sampling policy.

Move N16a into the immediate physical-safety work. N16b is a prerequisite for time-dependent molecular claims and N6's dynamical branch study. N16c may remain later. A strictly stationary/initialization-only project can proceed with Stage A while physical-trajectory claims remain explicitly withheld.

## 4. Numerical repair and test design

### F4. N0 must distinguish a reproducer from a correction test

**Evidence: SOURCE inspection of the audit harness.**

The earlier harness deliberately asserts that the old failure occurs. For example, its zero-operator GMRES test succeeds when the returned solution is nonfinite. That is appropriate for recording a defect, but the assertion cannot be copied unchanged into a permanent correctness suite.

There are also different levels of coverage:

- `pgmres`, `dogleg_step`, and the enthalpy diagnostic are extracted from production source;
- the Cauchy-length example calculates the old and correct formulas directly in the harness;
- the shell and shared-budget examples demonstrate geometry and algebra; they do not invoke the production census or future constraint map.

Consequently, N0 should preserve the historical reproducer and add production-linked correctness tests. It should not describe every existing row as an end-to-end production regression.

Specific changes to the test matrix:

| Planned row | Required target and expectation |
| --- | --- |
| Zero radius with nonzero dogleg legs | Test radius initialization/recovery in the caller. A geometric routine should not return a nonzero step inside a ball of radius zero. A named invalid-radius result is acceptable if that is its API contract. |
| Approximate-gradient Cauchy example | Invoke the actual production step-length calculation or trust-region path. A copied scalar formula will not detect a later regression in production. |
| Zero operator in GMRES | Expect finite output, correct nonconvergence status, and a nonzero actual residual for a nonzero right-hand side. |
| Ratio 0.5 is not stationary | Test the stationarity/adoption decision, not a demand that the diagnostic ratio change. The ratio calculation can remain correct. |
| Outer shell is counted | Invoke the same shell-selection/census function used by the spectrum code. |
| Shared budget is infeasible | Invoke the canonical constraint evaluator, not a second hard-coded sum in the test alone. |
| Error changes the next step | Exercise the production rejection, restoration, and clock boundary, not just `attempted_step_reduced_dt`. |

The raw zero-radius behavior is evidence of a caller-level absorbing state, not proof that the elementary dogleg geometry itself is wrong. This distinction should be explicit in the revised plan and preserves the intended correction without introducing an invalid nonzero step.

Keep all existing audit files under `docs/audit_20260905/`. Add the maintained production tests separately. The user's request to retain audit code should not be interpreted as permission to move away or overwrite its historical reproducer.

### F5. N1/N3 need a clearer linear-solve outcome contract

The proposed named outcomes are an improvement, but `truncated` is ambiguous. Reaching the requested subspace size with insufficient reduction is different from failing to sample a direction, a preconditioner failure, or a singular reduced problem.

Recommended outcomes distinguish at least:

- requested linear tolerance reached;
- subspace budget exhausted with a finite approximate step;
- Jacobian action unavailable;
- preconditioner failure;
- singular/inconsistent reduced problem;
- nonfinite arithmetic.

A small Arnoldi remainder does not by itself imply either success or failure. A finite useful step from a rank-deficient least-squares problem may be retained as an approximate direction, but its true residual and model reduction must be reported. A zero right-hand side must terminate successfully before unnecessary factorization or normalization.

N1 should specify a finite positive radius in scaled coordinates, how it recovers after a failed initialization, and how a lack of feasible descent is reported. An arbitrary positive radius cannot create a missing direction. N2 should retain the actual line minimum and the distinction between an approximate descending direction and an exact-gradient Cauchy guarantee.

These are numerical correctness decisions, not reasons to change physical certification. The official [KINSOL stopping description](https://sundials.readthedocs.io/en/latest/kinsol/Mathematics_link.html) likewise distinguishes residual convergence from small-step stagnation.

### F6. The plan misidentifies the Cauchy path of the atomic control

**Evidence: SOURCE.**

N2 says that the atomic route uses the Cauchy leg on the coupled default. In the inspected implementation:

```fortran
use_tr = (nspec_row .gt. 0)
```

The registry adds rows only for activated elemental or molecular transport. The checked `backup/regression/wasp_full_newton/input.inp` has no activated transport rows; the corresponding parameter defaults are false. It is the three-unknown control, not evidence that the species trust-region path was exercised.

Correct the explanation: atomic gas **with elemental transport** can use the trust region; the atomic three-unknown case does not. Merely setting a coupled-solve flag when the registry remains empty also does not create those rows.

N2 needs a species-row fixture to validate its changed behavior. N3 changes shared GMRES code, so both the three-unknown and species-row callers are relevant. An unchanged atomic control is useful, but it cannot certify the repaired Cauchy branch.

## 5. Constraint physics, reproducibility, and stage gates

### F7. N4a must define full inventory and conditional availability separately

The proposed shared hydrogen/oxygen/carbon constraints are necessary. However, adopting one expression named `available` everywhere without specifying the chemical closure could preserve an existing inconsistency under a new module name.

SOURCE findings:

- `species_table.f90` already contains hydrogen, helium, oxygen, carbon, and charge information for the relevant species, together with metal-stage metadata.
- `carrier_state` constructs `nO_free` and `nC_free` from totals that include atomic ion stages. They are not simply neutral reservoirs.
- `carrier_source` calculates atomic oxygen after subtracting `cbg_nOion` and the molecular carriers.
- the carrier write-back redistributes the remaining oxygen and carbon among their atomic stages.
- hydrogen availability changes when H+ is itself transported, and HeH+ removes inventory from both H and He accounting.

The plan must distinguish:

1. the full physical constraint, including every species in the current state;
2. the conditional constraint of a solve that actually holds some densities fixed;
3. the constraint after those species are allowed to re-equilibrate or are consistently redistributed.

For example, when oxygen ions are held fixed in a local source solve, the molecular oxygen budget must account for them. When they are eliminated consistently in a global solve, freezing the previous ion inventory would impose an artificial restriction. The correct response is not to subtract an old ion density indiscriminately at every call site.

Required deliverable for N4a: a table assigning every species to transported unknowns, conditionally eliminated unknowns, or fixed data in each context, with one stoichiometric source of truth. Build the constraint map from existing species metadata rather than a new independent reaction-species list.

The plan's requirement that every golden state be feasible should be changed to: **evaluate and classify every selected reference state under the physical map.** A formerly accepted state may fail the corrected constraint. That should expose a defect in the state or its earlier closure; it must not cause the new map to be weakened.

N4a can first be diagnostic on the marching path. In that mode it does not yet establish that trial construction, chemistry, and certification all enforce the same constraints. The stage gate must not claim completion of that later integration.

### F8. Stage B's gate is incomplete, and N8 needs two milestones

The dependency sketch includes N4b and N8, but the written Stage B gate explicitly names only N4a and N5. It is therefore possible to satisfy its stated tests without implementing the constrained step that the next stage relies on.

Require all of the following for the relevant species set:

- N4a's full and conditional inventories are defined and tested;
- N4b can move along the coupled constraint surface and handles infeasible initial data;
- N5 establishes replay determinism and conditional closure accuracy;
- N8a establishes usable row scales, boundary derivatives, and an inner-error budget;
- a nonzero physical balance at a bound is still refused by certification.

Split N8:

- **N8a, before the solution campaign:** manufactured states, scaling, boundary equations, derivative tests, and error amplification.
- **N8b, after candidate solutions:** grid refinement and the final empirical basis for carrier/element certification tolerances.

Requiring final physical-grid tolerance calibration before the first candidate solution would create a dependency loop. Conversely, the first candidate should not immediately establish the final tolerance by definition.

Split Stage C into C1 atomic element transport, C2 H2 transport, and C3 oxygen carriers. Give each its own entry conditions and acceptance evidence. The wording “one solution” followed by all three families currently leaves the stopping point unclear.

### F9. N5 must not demand global seed independence across genuine physical branches

Replaying the same declared state must be deterministic. But the proposed test varies every admissible fast-chemistry seed while N6 simultaneously allows multiple self-consistent radiative/chemical roots.

These requirements are compatible only after stating the branch contract. Separate:

- **replay determinism:** identical complete inputs and branch data return the same result regardless of unrelated prior probes;
- **conditional closure convergence:** seeds within a stated branch converge to its same fast solution;
- **branch discovery:** other admissible seeds may identify distinct roots and must be recorded, not silently averaged or forced onto a preferred root.

Where the eliminated system is genuinely multivalued, either its branch/history becomes explicit data or the slow mode must remain an unknown. A test that forces all seeds to one answer can hide real bistability.

The B5h case is already refused in its unsafe public configuration, and the issue inventory reports much smaller seed sensitivity with H2 transported. Therefore, do not label the proposed production N5 test “RED on current source” solely from the earlier forbidden configuration. State which internal diagnostic reproduces the old issue and test the current supported configuration separately.

N6's first local fixed-radiation study can start early as a diagnostic. Its time-dependent stability stage must wait for N16a/N16b. Decision 18 must not block an unrelated atomic element solution.

### F10. Stage D is only partly independent

N9's round-trip investigation, N12's caloric correction, and N13's shell selection are largely independent of Stage A. N10's validation on a trustworthy stationary state is not: it requires a reproducible residual and a certified fixture. N11's physical acceptance policy also depends on a valid stationary evaluator.

Distinguish implementation availability from the evidence needed to certify its physical use. Unit tests may proceed independently; claims about a molecular stationary correction may not.

For N12, passing the upstream specific internal energy is preferable to adding only another H2-fraction argument. That directly represents the quantity transported and covers a current atomic cell with a molecular upstream neighbor. Keep the thermal/formation-energy convention consistent; do not add a dissociation energy term twice.

For N11, replacing the current term ratio with the correct mass residual is necessary but not sufficient. A post-process can preserve an input density/velocity field while changing pressure or species fluxes. State whether the product is a one-way diagnostic on that fixed wind or a fully self-consistent solution, and do not let a successful scalar temperature solve imply the latter.

## 6. Execution rules that should change

### F11. Replace universal bitwise requirements with affected-path validation

Section 7 requires every item to measure `wasp_full_newton`, including unrelated restart metadata, transit masks, and documentation. This is unnecessarily broad and can consume most of the execution budget without testing the changed path.

Recommended scope:

| Change | Relevant validation |
| --- | --- |
| N1/N2 trust-region behavior | numerical step tests and a short species-row fixture |
| N3 shared GMRES behavior | numerical tests, three-unknown caller, species-row caller, targeted bounds/FPE checks |
| N4/N5 constraints and residual | selected carrier/element configurations and replay/derivative tests |
| N9 restart metadata or reader | round trip, compatibility/failure tests, no new hydrodynamic evolution |
| N12 caloric residual | thermodynamic unit tests and a small composition-gradient profile |
| N13 census or spectral metadata | synthetic ray geometry and reader tests using fixed existing arrays |
| N16 controller | forced-error/retry integration tests and temporal convergence |
| N17 documentation | source and link verification only |

Bitwise comparison remains useful where exactly the same serialized state or untouched arithmetic is intended. It is not a universal numerical acceptance criterion. In particular, corrected stable norms, rotations, or reorthogonalization in shared GMRES may change rounding while improving correctness.

The rule never to refresh a reference for a change below 0.1% is also unsuitable as a general policy. It does not define which quantity is measured and can conflict with a bitwise comparator. Preserve historical references, document changed physics or numerics, and choose reference updates from the validated change and comparison contract, not a universal percentage. No update should conceal a failed physical test.

### F12. Private build directories do not isolate source edits

N1/N3 and N2 all change `steady_newton.f90`; the `pgmres` result contract is consumed inside the trust-region routine. Two private object directories still compile the same live source file. A build can observe one worker's API change and another worker's unfinished caller.

Recommended first implementation: one owner for the coordinated N3 -> N1 -> N2 change, with separate tests/review work if explicitly authorized. Alternatively, use truly isolated source snapshots and an explicit integration sequence. Do not rely on `OBJDIR=build_<item>` as source isolation.

Private outputs also need isolation. Selected test scripts use fixed working paths or regression inputs; inspect those paths before launching simultaneous tests. The Makefile performs generated-source work while being parsed, so even a nominal inspection command should not be assumed to have no file effects.

No agents were delegated work in this review. The plan's proposed worker assignments are part of the document being reviewed, not an instruction to launch them now.

### F13. Bitwise physical restart is an overstrong default gate

An unchanged-state round trip is the correct goal, but its comparison must distinguish serialization from reconstruction. The writer uses dimensional quantities; the loader reconstructs quantities and applies unit conversions and inventory calculations. Even round-trippable decimal scalars do not guarantee bitwise identity after those arithmetic operations.

Use separate tests:

- stored scalar fields round-trip under a specified format;
- conserved state and inventories differ only within a measured floating-point reconstruction budget;
- no chemistry, remapping, boundary reinitialization, or physical-time advance occurs accidentally;
- the stationary residual remains within its independently established evaluation budget;
- a second round trip does not cause secular drift.

If bitwise restart of internal conservative arrays is required, design and approve that checkpoint format explicitly. Do not force a derived human-readable profile file to satisfy an unstated binary-checkpoint contract.

The existing `heh_dev > 1e-6` check should not simply be loosened. N9 correctly begins with the actual failing file and the origin of its inventory difference.

### F14. State and artifact ownership need explicit limits

N17 is described as both “now” and a first action after the pause lifts, while the plan says nothing has started. Clarify whether documentation edits are authorized separately from code. This review request authorizes the review file; it does not itself authorize editing the issue inventory or code-status document.

Resolve the missing scratch paths before execution. A reproducer should identify its input files, restart files, source/binary identity, controls, and expected diagnostic entry point. Preserve the historical audit files and keep future review/test code in the requested audit directory when it is not part of a maintained production suite.

## 7. Recommendations for the user decisions

These are recommendations, not recorded user approvals.

| Decision | Recommendation | Scope and reason |
| --- | --- | --- |
| **13: lift the code pause** | Approve a bounded first increment after plan revision | N0, coordinated N3/N1/N2, and the short N7 diagnostic. Keep current running jobs untouched. Add N16a as immediate work if physical-mode calculations are to be relied upon. Do not interpret this as approval for all stages or planetary regeneration. |
| **14: constrained step or new parameterization** | Start with route (i), a coupled linearized constraint set on the present unknowns | This directly tests the diagnosed feasible-set problem without simultaneously changing every unknown representation. Approval should require the closure/inventory contract, nonlinear feasibility checks, and physical residual certification. Keep route (ii) open if conditioning or zero-abundance handling makes it the technically better formulation. |
| **15: restart metadata and intent** | Approve a versioned restart contract in principle; approve exact keys only after a concise interface design | Store sufficient state/provenance and separate trajectory continuation, relaxation, and stationary evaluation. Prefer one coherent intent selector if needed, not three independent Boolean switches. Define consistency with existing `Run mode`, `Solver`, and load controls. |
| **16: temperature/composition status** | Approve explicit separate status fields with a versioned schema | Require unambiguous meanings for corrected, retained, failed, unsupported, and not evaluated. Preserve simultaneous reasons where relevant. Old files without the fields must read as unknown legacy validity, not as successful correction. |
| **17: transit metadata** | Approve structured comment metadata while keeping numerical spectral columns unchanged initially | The writer already emits a comment header and selected consumers use `np.loadtxt`. Additional comments are a smaller change than a numerical-column redesign, but check custom readers. Do not regenerate atmospheric or spectral data solely to add metadata. |
| **18: two stable branches** | Defer the final physical choice until the branch evidence exists; approve the diagnostic study | Initial/boundary history is meaningful input, not permission to select the preferred output. Record both branches and their domains if both survive. This decision blocks only branch-dependent interpretation, not all of Stage C. |
| **Existing regeneration gate** | Keep unchanged | No automatic reruns of planets, benchmarks, paper, or poster products. |

### 7.1 Conditions on decision 14

Route (i) should not become another coordinate-box projection. The step model must include shared elemental constraints, density-dependent totals, and derivatives through logarithmic carriers. Re-evaluate the nonlinear constraints at the candidate state and restore feasibility when needed.

Use dimensionless log coordinates such as `log(n/n_ref)` with a documented reference. Species at exactly zero elemental abundance need an explicit absent-species treatment; a positive logarithmic floor must not invent nuclei. Neither a new parameterization nor an active-set method removes the need to handle this case.

Also avoid replacing the equation `F(Y)=0` with a variational inequality merely because a library supports bounds. PETSc's [reduced-space VI method](https://petsc.org/release/manualpages/SNES/SNESVINEWTONRSLS/) solves a bounded variational-inequality problem; its convergence must not be equated with zero EXHALE balance rows at active bounds. The user's decision 12 c remains the physical acceptance rule.

Choose route (ii) if the evidence supports it, not merely if route (i) fails an arbitrary iteration budget. Before approval, show independent degrees of freedom, treatment of shared H/O/C/He species, zeros and trace abundances, derivative conditioning, and restart impact.

### 7.2 Conditions on decisions 16 and 17

Changing a status schema does not by itself certify a result. Store the input state's stationary certification, source identity, post-process model restrictions, and the meaning of each status alongside the output.

An optical-depth share of refused material is a contribution diagnostic, not an error bar on transit depth. Transmission depends nonlinearly on total optical depth; a saturated contribution does not produce a proportional change in depth. Label the quantity accordingly and use controlled profile perturbations or validated alternate states for uncertainty estimates.

The radius-mask bug in N13 can be corrected and tested without deciding a new metadata format. Similarly, the upstream-energy correction in N12 need not wait for a restart interface decision. Approval should block the dependent write/schema action, not independent analysis or numerical correction after the code pause is lifted.

### 7.3 Add a physical-time policy decision

The plan leaves an important policy implicit: is physical mode intended to produce a trajectory with controlled numerical error, or only to expose a diagnostic estimate?

Recommended new decision **19**:

> Physical-mode results intended for scientific time-dependent interpretation must use error-driven rejection at a common endpoint, include transported species, and validate the accepted-state estimator. Until that contract is implemented and tested, existing diagnostic estimates must not be described as error-controlled integration.

Every-step control is the recommended initial reference implementation. An intermittent policy should be considered later only with an explicit validated error argument. Enabling the current estimator every step alone does not fix F1 or F2.

The extra cost should be measured. Three attempted advances for step doubling do not imply a universal wall-time multiplier because their nonlinear costs can differ. Keep initialization mode separate; it need not pay for physical-trajectory error control when it makes no trajectory claim.

## 8. Revised dependency structure

The following structure makes the necessary dependencies explicit without forcing unrelated tasks to wait:

```text
Plan corrections and artifact/permission audit
    |
    +-- N0 correctness tests -> N3 -> N1 -> N2 -> short N7 diagnostic
    |
    +-- N4a inventory/closure definition -> decision 14 -> N4b
    |                                       |
    +-- N5 replay and branch contract ------+-- Stage B gate
    +-- N8a scaling and boundary derivatives-+        |
    |                                               v
    |                                    C1 atomic -> C2 H2 -> C3 oxygen
    |                                               |
    |                                      N8b tolerance/refinement study
    |
    +-- N9 reader/writer tests -> decision 15 -> N10 interface/routing
    |                                                |
    |                                  N5 + certified fixture -> N10 gate
    |
    +-- N12 caloric correction
    +-- N13 shell census -> decision 17 -> metadata integration
    +-- N11 validity design -> decision 16 -> schema and reader tests
    |
    +-- N16a endpoint/rejection/clock -> N16b species/time accuracy
                                                   |
                       N6 local/radiative branches -> N6 dynamic stability

N5 + source accuracy references -> N14 performance experiment
N15 reference checks run before and alongside N14, not only afterward
```

N7 may stop after a meaningful diagnostic if the complete nonlinear solve is not yet ready. Its first result should identify the failure mode, not be forced to produce a converged planet. Stage C's physical certification remains the later criterion.

## 9. Minimum revised acceptance matrix

| Area | Required evidence before completion |
| --- | --- |
| Radius initialization | Recovery from a zero first direction; finite bounds; zero-residual success; named failure when no feasible descent exists |
| Cauchy/model step | Actual production length tested with exact and approximate gradients; projected-step model agreement; insufficient Krylov reduction handled explicitly |
| GMRES | Identity, zero RHS, zero operator, compatible singular and incompatible singular systems, near breakdown, unusable Jacobian action, and preconditioner failure |
| Element constraints | Complete species census; shared budgets; fixed versus eliminated ion contexts; zero/trace inventories; current-density dependence; inward motion from a constraint surface |
| Residual | Replay independence from discarded probes; conditional seed convergence; explicit branch handling; inner-error and finite-difference interval measurements |
| Stationary solution | Complete active row registry; dimensional and normalized balances; conserved fluxes; independent re-evaluation and initialization; subsequent spatial/tolerance refinement |
| Restart | Correct scalar serialization and conservative reconstruction; no unrequested state evolution; invariant/energy checks; clock and intent semantics; no secular round-trip drift |
| Post-process | Upstream mixture energy; actual discrete stationarity test; explicit one-way/model restrictions; status schema round trip |
| Transit | All contributing shells and correct interpolation support; physical-cell selection; metadata readers; contribution diagnostic not mislabeled as uncertainty |
| Physical controller | Error above threshold actually rejects; common endpoint after every retry; whole-transaction restoration; correct clock; species-aware norm; temporal convergence |

For tests not executed yet, use “expected to expose the defect” rather than declaring RED as a measured result. Documentation corrections and new diagnostics do not need an artificial failing numerical test. Record a failure mechanism before repair where practical, and a directly relevant correctness test afterward.

Targeted bounds and floating-point checks belong at the numerical/controller gates. They need not become a full planetary regression campaign. A passing unit test for a policy function is insufficient when the production caller can bypass that policy, as F1 demonstrates.

## 10. Suggested approval statement and final assessment

A concise approval, if the user chooses to proceed, would be:

> Revise the plan to include the production error-rejection and common-endpoint defects, corrected test targets, explicit Stage B dependencies, and affected-path validation. Then implement N0 and N3/N1/N2 with coordinated source ownership and run only the short N7 diagnostic. Keep running calculations untouched. Prepare the N4a/N4b and restart/status interface designs for the stated decisions. Do not regenerate scientific products or treat the remaining stages as automatically authorized.

If physical trajectories are needed immediately, add N16a/N16b to that authorized scope and withhold quantitative trajectory claims until their end-to-end tests pass.

The plan already identifies the right physical destination: species transport, energy closure, reproducible residuals, and certification of the equations actually solved. Its revision should focus on the missing production connections and precise gates, rather than adding more broad studies. The most consequential new finding is that an integration-error routine and its controller formula can both exist and be unit-tested while the main program still accepts the step without applying them.
