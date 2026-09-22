# Review of `PLAN_20260920_rev1.md`

Review date: 2026-09-20 (KST).

## Overall assessment

Revision 1 is a substantial improvement over the original plan. It correctly separates stored generations from log-only iterates, removes the unsupported physical-limit conclusions, corrects the M1 interpretation, identifies the actual reachability of the diagnostic flags, and corrects the M2 order-of-magnitude error.

The plan is now suitable as a review-driven investigation design, but it is not yet ready to execute unchanged. Four points require correction before expensive runs:

1. `On stall` is indeed unreachable in the present pass ordering, but repairing it requires an explicit handover-state contract, not only a saved outcome token.
2. The proposed "current binary and current boundary" rerun is too categorical. A loaded old generation can still be a valid experiment of the current operator if the boundary is deliberately rebuilt and all differences are labeled.
3. The closure-sweep proposal must control the actual sweep count. `n_eq_sweeps_model` is initialized to one but is overwritten by measured sweep counts during the solve; a statement that the production Jacobian always uses one sweep is not generally correct.
4. A few diagnostic-output claims need an explicit file/exit gate. A requested diagnostic that produces no file must be recorded as "not executed," not interpreted as a clean result.

No production code was changed and no new atmospheric solve was run for this review. The conclusions below are based on direct source inspection and stored records. The earlier audit artifacts remain in `docs/audit_20260905/plan_20260920_review/`.

## 1. What rev1 corrected successfully

### 1.1 State identity and row locations

The revised two-table presentation is correct and important. The selected generations and the 2026-09-20 log-only iterates do not describe the same numerical state. They differ in binary, boundary identity, generation phase, and sometimes in the cell attaining the largest row measure. A report must not combine a row value from a log with a state file selected by `state_index.json`.

The rev1 table also correctly distinguishes a row maximum from a tolerance-normalized verdict. This is required because the mass tolerance is cell dependent. The displayed maximum residual measure and the cell that is closest to the certification threshold can be different.

One wording should be tightened. The table labels some values as "READ" while the radius values are measured by indexing the generation file. Use separate labels such as `READ: certificate value` and `MEASURED: radius looked up in the selected generation`. This prevents later readers from assuming the entire row was recomputed.

### 1.2 Diagnostic reachability

The conclusions about the Stage C diagnostics match the implementation:

- `EXHALE_PRECOND_SPECTRUM`, `EXHALE_BAND_DIFFERENCE`, `EXHALE_KRYLOV_SIZE_SCAN`, and the species binding-row report are inside `trust_region_step` and require the element diagnostic gate.
- `EXHALE_TRUST_REGION=1` changes both the step controller and the merit scaling.
- `EXHALE_KRYLOV_TOL_ABS` tightens a forcing term under restricted guards; it is not an absolute tolerance above an evaluation-noise floor.
- `EXHALE_LINEAR_ROWS` is the relevant existing diagnostic for the default three-unknown route.

The distinction between observation and numerical variant is now correct. Keep it in every run directory and manifest, not only in prose: record the complete environment-variable set, including variables explicitly set to `0`.

### 1.3 M1 interpretation

The revised description correctly identifies the default path as an alternating block iteration. `relax_element_composition` and `relax_photochemical_composition` do act on the carrier and elemental equations outside the Newton unknown vector. Therefore a small hydrodynamic norm together with a large carrier norm means that the joint state is not certified; it does not mean that the carrier equation was never advanced.

The proposed measurements are appropriate: requested versus adopted composition displacement, the relaxation ending token, post-refresh joint residual, and possible cycling must all be recorded. A scalar maximum alone is insufficient.

### 1.4 M2 interpretation

The revision correctly treats atomic-to-molecular conversion as a change of model and state, not as a root-preserving transformation. It correctly documents the `p` and `T` invariant modes and the existing `handoff`, `local`, and explicit-fraction seed choices. The correction from approximately two decades to 10.01 decades and the distinction between mass cell 2 and energy/momentum cell 1 are also correct.

## 2. Important remaining correction: `On stall` handover

### 2.1 The reachability diagnosis is correct

The current ordering is:

1. The progress logic can set `outer_ending = outer_no_progress`.
2. The pass then resets `carrier_outcome = carrier_relax_nothing_to_advance`.
3. The composition update is guarded by `outer_ending == outer_running`, so no relaxation runs on that pass.
4. The handover condition later requires both `outer_ending == outer_no_progress` and `carrier_outcome == carrier_relax_movement_bound`.
5. The no-progress branch returns from the stationary loop.

Thus the current-pass outcome cannot satisfy the handover condition. This is a real control-flow defect, not merely an untested option.

### 2.2 The proposed repair is incomplete

Rev1 proposes retaining the last relaxation outcome. That is necessary but not sufficient. A safe handover needs to preserve, from the same last relaxation pass:

- the composition state actually adopted;
- the fixed-wind state and face mass-flux array used by that relaxation;
- the carrier and elemental residual vectors, scales, and worst-cell identities;
- the bound value and the exact reason for the movement-bound ending;
- the outer progress reference and the number of consecutive no-fall passes;
- the boundary/reconstruction state that the coupled solve will consume.

Otherwise the coupled block can be entered with a residual measured on one composition and a bound/outcome token belonging to another. That would recreate the state/path inconsistency the plan is trying to diagnose.

The handover should be implemented as an explicit state transition, for example `alternation_stalled_at_bound -> coupled_block_entry`, with one immutable entry snapshot. The transition must occur before the no-progress return, or the no-progress decision must defer return when `carrier_newton_on_stall` is enabled and the saved last-relaxation contract is satisfied.

Required tests before using the option:

1. A synthetic or archived case that reaches the movement bound and no-progress condition.
2. A log assertion that handover occurs exactly once.
3. A check that the coupled block's first residual equals a fresh residual of the handover snapshot.
4. A control with `On stall` disabled that follows the old refusal path.
5. A case where no movement-bound ending occurred, proving that no handover is triggered by a merely flat scalar norm.

Do not report `On stall` as a solver comparison until these tests pass.

## 3. Current boundary versus old restart state

Rev1 says that the "honest choice" is to rerun all three P2 cases under the current binary and boundary before measuring anything. This is too strong as a general rule.

`load_IC` explicitly detects a boundary-model mismatch, rebuilds the boundary from the physical column and current reservoir, and states that the residual is not the residual under which the restart was written. That makes an old generation a legitimate, clearly labeled **current-operator-on-old-state** experiment. It can answer whether the current boundary law evaluates the old state consistently, whether the residual is reproducible, and which terms change when the boundary is rebuilt.

It cannot answer whether the current boundary law has a converged state on the current branch. That question requires a current-model solve or a current-model checkpoint. The plan should therefore define two valid modes:

- **Mode R:** reload the old state, rebuild the current boundary, and measure only. This is read-only and useful for A0/A/B diagnostics.
- **Mode C:** run or continue the current binary and current boundary until a publishable checkpoint exists. This is required for claims about current-branch convergence, continuation, or solver performance.

The report must never call Mode R a reproduction of the historical certificate. The boundary model, reservoir, input hash, and whether a physical step was taken must be printed beside every result.

This distinction can save a costly campaign while preserving scientific honesty. A rerun is required only when the question depends on a current-branch state, not for every reproducibility or operator-localization test.

## 4. Closure-sweep experiment: correct the control variable

The proposed closure-bias experiment is valuable, but its implementation description is currently too simple.

`read_composition_elimination_controls` initializes `n_eq_sweeps_model = 1`. During the solve, however, `n_eq_sweeps_model` is updated from `n_eq_sweeps_last` after the adopted evaluation and can be increased after accepted trials. The Newton-loop Jacobian call uses `n_eq_sweeps_fixed=n_eq_sweeps_model`, so the map is fixed for an individual Jacobian evaluation but is not necessarily a one-sweep map throughout the run.

The experiment must therefore record, for every base and perturbed evaluation:

- requested fixed sweep count;
- actual sweep count completed;
- convergence criterion and whether it was reached;
- seed composition identity;
- number of chemical roots without convergence;
- closure increment, temperature change, and heating/cooling change.

Do not implement the experiment by merely changing `EXHALE_RESID_SC_MAX`; that variable is a cap and does not necessarily force the exact fixed count required for a comparable map. Use the existing `n_eq_sweeps_fixed` interface in a dedicated diagnostic path or add an explicit, documented diagnostic override. Any added override must be disabled by default and must not alter the production route.

The correct comparison is two-dimensional:

```text
probe length h  x  fixed closure count k
```

For each `k`, use identical seeds and endpoint evaluation order. Compare the action with the same `k` at both endpoints, then compare the resulting actions across `k`. A stable finite-difference plateau at one `k` does not show that the converged closure derivative is correct. Conversely, a change with `k` does not prove a finite-difference defect.

The plan should also state which residual is being judged. `EXHALE_RESIDUAL` can use the loaded-state diagnostic route, while the Newton Jacobian uses the fixed-count `eval_residual` map. These are related but not automatically identical maps. The closure experiment must not compare their numbers without recording the route and sweep policy.

## 5. Carrier row-term diagnostic: useful, but add an execution gate

The claim that `EXHALE_CARRIER_ROW_TERMS=1` requires no new code is correct in the current source. The option enables the row dump, and certification writes `output/carrier_row_terms.txt`; the dump contains face terms, transport terms, production/loss, reaction channels, scales, and the residual.

However, the output is conditional. It is written by the certification/report path and represents the last assembled carrier rows. A run that terminates before certification, times out, or never assembles a carrier row may produce no file. The plan must require:

- file existence and nonzero size;
- header and row-count validation;
- matching state/generation identity;
- a statement of whether the rows came from the final state, a discarded trial, or an intermediate Newton evaluation;
- explicit exit class.

The dump is not by itself an independent conservation audit. It records the production operator's terms. An independent budget still needs a separate summation of face fluxes, cell volumes, and sources, with signed and absolute defects.

## 6. Stage B/C numerical interpretation

The revised decision table is much safer, but two refinements are needed.

### 6.1 Actual perturbation after projection

The plan correctly asks for requested and realized perturbations. Add the full active-set mask and the projection displacement in every direction. If a direction is clipped, its finite difference is a derivative of a projected map, not the unconstrained residual. A plateau in such data should not be used to certify the unconstrained Jacobian.

Also record whether the reconstruction mode and limiter weights were frozen. The production action and trial residual use different reconstruction settings in parts of `steady_newton.f90`; the comparison must state whether it measures the production action, a frozen model, or the nonlinear trial operator.

### 6.2 `EXHALE_LINEAR_ROWS` interpretation

The linear-row diagnostic computes the residual of the shifted model for the returned direction. This is useful, but it is not the nonlinear residual of the accepted or rejected state. Keep three quantities separate:

```text
r_model = -F - J_action(s) - shift(s)
r_trial = F(Y+s)
r_state = F(Y_adopted)
```

The first validates the solve of the supplied action; the second tests prediction; the third is the accepted-state quantity. A small `r_model` with large `r_trial` is not sufficient to classify the Jacobian as correct or incorrect.

### 6.3 Linear-model row scaling

The plan should state that row scales are formed from the latest residual evaluation. In the coupled route, the species row scales and constraint boxes are derived from the current evaluation and can be stale if a diagnostic evaluates another state without restoring all module products. Every Stage C run should record the state at which `Drow` and active bounds were formed, and verify that it is the same base state used for the action.

## 7. M1 additions

The M1 redesign is sound, but it should add two safeguards.

First, compare the full state vector between outer passes, not only carrier displacement and scalar joint distance. A moving worst cell can make a maximum appear flat while the state is changing; conversely, a flat maximum can hide a two-cycle. Save a compact norm of `Y_{k+1}-Y_k`, `Y_{k+1}-Y_{k-1}`, and the composition vector displacement.

Second, verify that `Coupled carrier solve: True` is compatible with the selected transport options. The input reader rejects the combination with `Ionization transport: True`, because the variables and bounds are different. This restriction must be checked before selecting a case; otherwise a failed launch is a configuration incompatibility, not a solver result.

The proposed one-case coupled trial should include a matched `False` run with the same current checkpoint and budget. A previously certified hydrodynamic state is not enough if the comparison begins after one route has already refreshed composition and the other has not.

## 8. M2 additions

The invariant audit should include the actual base ghost and face state after `Apply_BC`, not only interior density, pressure, and composition. The conversion routine's pressure or temperature invariant is a cell-state property; the lower boundary reconstructs ghost quantities from composition and reservoir data. A conversion can preserve cell pressure while changing the ghost pressure/particle count used by the base flux.

For the energy residual, use the code's explicit convention. `steady_residual.f90` assembles

```text
R_energy = dF_energy - S_energy - (heat - cool)
```

and `Num_Fluxes.f90` uses the energy flux `v*(E+p)`. The audit should report each term in this convention and should not add molecular binding energy a second time. The caloric EOS documents bound H2 internal energy separately from chemical source accounting; this is precisely why the invariant audit must be term based.

## 9. Additional physical/numerical recommendation

The plan should add a **residual-route equivalence test**. There are at least three related paths: loaded-state `EXHALE_RESIDUAL`, Newton `eval_residual`, and the marching/RK operator. They intentionally differ in whether a state is advanced, how many composition sweeps are used, and which reconstruction mode is selected. For one checkpoint, evaluate all three at the same physical state where possible and compare:

- mass, momentum, and energy row vectors;
- base and outer ghost states;
- face mass fluxes;
- heating/cooling arrays;
- composition and particle counts;
- carrier and element row vectors.

If the paths are intended to represent the same stationary equations, their differences should be explained term by term. If they intentionally represent different maps, the plan must name the map used by each certificate and derivative test. This check is more informative than comparing only the final norm and can expose a route mismatch before a long continuation campaign.

## 10. Revised execution order

I recommend the following order for rev1:

1. A0 identity and Mode R/Mode C selection.
2. Read-only residual-route equivalence and full-vector repeatability on one checkpoint.
3. Closure-count versus perturbation-length sweep on the same checkpoint.
4. Independent carrier and hydrodynamic conservation budgets.
5. M1 full-state alternation audit and one matched coupled trial only if the configuration is admissible.
6. M2 invariant and base-ghost audit.
7. P2 Stage C/D diagnostics and only then continuation.

This order is more efficient because it checks whether the maps and balances are coherent before spending time on spectral diagnostics or a new XUV ladder.

## 11. Final verdict

`PLAN_20260920_rev1.md` is substantially reasonable and is a good basis for controlled investigation. It should be accepted with the following mandatory edits:

- implement and test `On stall` as an explicit state transition, or remove it from the experimental options until repaired;
- replace the unconditional current-boundary rerun requirement with clearly labeled old-state/current-operator and current-branch modes;
- specify an exact closure-sweep control and record actual versus requested sweep counts;
- add output-existence and exit-status gates for carrier row diagnostics;
- add residual-route equivalence and full-state/cycle checks;
- keep the independent conservation audit separate from production row dumps.

After these edits, the plan's strongest scientifically defensible conclusion remains: it can identify whether a particular stored/current operator and numerical route are reproducible and where their balance fails. It cannot, by itself, prove nonexistence of a stationary atmosphere or establish that a remaining failure is physical.

## 12. Documentation and reproducibility corrections

The opening statement that every review claim was checked "by three independent verifications" is not itself evidence. Rev1 gives source locations and stored values, but it does not provide three independently rerunnable commands or three independent result files for each claim. Replace that sentence with a narrower statement such as "the decisive claims were checked against the cited source paths and stored records." Where a claim is computed, retain the command, input hash, output hash, and exit status in the audit directory.

The statement that the earlier H3+ probe "was not rebuilt by me" is also not useful in the plan unless the document identifies which results are inherited and which were independently inspected. The audit directory already contains a successful direct compilation and execution log from the preceding review. Rev1 should either cite that artifact as inherited evidence or omit the author-specific qualification.

Section 10.1 contains a duplicated phrase, "Proposed repair," in consecutive sentences. Consolidate it into one precise repair statement and include a regression test for the pass-order transition. This is editorial, but ambiguity in a control-flow repair can create a materially different handover behavior.

Finally, each `MEASURED` result in the plan should have a stable artifact reference. A number calculated interactively and copied into prose is not enough for a later code review, especially when the executable hash is changing. The minimum record is: command, working directory, binary hash, source-manifest hash, environment, input/generation identity, output path, and exit status. A timeout must never be converted into a normal measured ending.
