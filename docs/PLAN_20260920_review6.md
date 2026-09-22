# Review of PLAN_20260920_rev6.md

Review date: 2026-09-20 (KST).

## Decision

Accept rev6 as the investigation plan and proceed with its initial diagnostics. The substantive additions from review5 agree with the current implementation: all four momentum branches are represented, cancellation is restricted to shared conservative fluxes, and reconstruction continuation has an explicit diagnostic policy. I found no new blocking physical or numerical error in these revised equations.

Two statements still need correction: the old Mode R/C definitions were not actually removed, and the claim that the momentum row cannot be reconstructed from "the final arrays at all" is broader than the code supports. Neither requires another complete planning revision. The recommendations below specify how to finish the document and implement its audit without expanding the solver work.

## Scope and evidence

This review compares rev6 with rev5 and review5 and inspects the relevant executable statements in `RK_rhs.f90`, `Num_Fluxes.f90`, `steady_residual.f90`, and `EXHALE_main.f90`. It also distinguishes assembled rows from diagnostic term arrays and restored reconstruction flags. Source locations below are relative to `src/`.

The conclusions are based on source and document inspection, not new solver runs. Historical certificate values, timings, H3+ probes, and finite-difference measurements were not reproduced. No production code, inputs, or solution products were changed. This document is the only new artifact.

## 1. The physical and numerical corrections are satisfactory

### Momentum branch table

Section 13.4 correctly reproduces the assembly at `modules/time_step/RK_rhs.f90:231` through `:270`:

- Ordinary PLM uses the area-weighted numerical momentum flux difference and subtracts its explicit source.
- Ordinary WENO3 additionally includes `(pR-pL)/dr`.
- The well-balanced PLM branch replaces that expression with the area-weighted flux and equilibrium-pressure departures.
- The well-balanced WENO3 branch replaces it with the stored flux divergence plus the departure gradient.

The well-balanced block replaces the earlier expression; it does not add another ordinary WENO pressure gradient. Rev6 states this correctly. The viscous momentum source must subsequently be subtracted where active, as the plan requires.

### Internal-face cancellation

The restriction in section 13.2 is correct. `modules/flux/Num_Fluxes.f90:167` through `:168` explicitly gives `q_dn = q_up - dp_eq`. The two departures reference different equilibria, so demanding their equality would reject intended arithmetic. The complete momentum balance is the appropriate check.

### Energy

The retained energy identity and the separate gas-energy check remain appropriate. `RK_rhs.f90:273` through `:275` carries gravitational energy work inside `dF(3,:)`; it must be included even though the explicit energy source is zero. The potential-weighted mass residual must remain in the augmented-energy identity at a nonstationary checkpoint.

### Reconstruction continuation

Section 13.5 correctly identifies the missing branch provenance in a consumer that reads only the last pressure-departure arrays and the restored method flags. `modules/time_step/steady_residual.f90:163` through `:211` evaluates the branches, blends their outputs, and restores the caller's flags. Capturing both branch ledgers before blending is a valid complete solution. Explicitly limiting the first diagnostic to pure endpoints is also a reasonable initial scope.

These findings close the main technical requests from review5. They do not establish that the proposed export or its consumer has already been implemented or tested.

## 2. Narrow the statement about unavailable momentum reconstruction

The introduction says that the momentum row inside continuation "cannot be rebuilt from the final arrays at all." Verdict-table item 3 similarly refers broadly to "the final arrays." Section 13.5 is more precise: it identifies the final departure arrays as insufficient.

The broader statement is incorrect. The continuation routine returns blended `dF` and `S` and explicitly blends the following diagnostic arrays:

```text
momentum_ram_divergence
momentum_pressure_gradient
momentum_gravity
equilibrium_pressure_force
```

Evidence: `steady_residual.f90:185` through `:205`. These outputs exist after the blend; the assembled momentum row is available from `dF(2,:)-S(2,:)`, followed by the transport-source subtraction. The stored term arrays also support attribution, although their interpretation must respect the balance option.

This does not make summing the production row an independent audit, and it does not recover the two original branch face ledgers. The actual limitation is insufficient information for an independent reconstruction from only the final face-pressure/departure data and restored flags.

Recommended replacement in the introduction and verdict table:

> Inside reconstruction continuation, the final pressure-departure arrays and restored method flags do not contain the two branch ledgers needed for an independent face-level momentum reconstruction. The blended residual and blended momentum diagnostics remain available.

This qualification also applies to any overly broad reading of review5. The physical problem has not become unauditable; the proposed face-level audit needs additional records.

## 3. The mode cleanup was claimed but not completed

Section 2 says "The obsolete mode sentences are removed," but section 3 still contains:

- Mode R: "with the current binary and measure only," followed later by a table allowing one restored disposable update.
- Mode C: "until a publishable checkpoint exists," followed later by text admitting a short failed trajectory.
- "Anything read out of an existing `run.log` is a historical audit," while section 4.1 correctly says that putting a new diagnostic result into a log does not change the experiment's identity.

Replace the first two definitions with the following concise definitions:

> Mode R evaluates the current operator on an immutable checkpoint, optionally including one explicitly identified disposable update followed by restoration. Record the input state, actual operations, and output identity.

> Mode C advances a retained trajectory under the current solver, with a declared budget and recorded ending. A published or certified checkpoint is not required for the trajectory to provide diagnostic evidence.

For log evidence, distinguish the act of reading a stored artifact from the identity of the experiment it records. Retain that experiment's executable, state, and mode metadata. A filename does not establish whether it used the current operator.

Also change phase 1's label to "historical audit for stored identity records; R for new evaluations" and qualify section 16's assertion that phases 1 through 5 are all Mode R. Section 4.1 already has the correct distinction. These are consistency edits, not reasons to delay the experiments.

## 4. Concrete implementation choices

I recommend the assembly-point export in section 13.2 and the following initial scope:

1. Support the four pure momentum branches first, with the energy and mass identities already specified.
2. Record whether reconstruction continuation is enabled and its lambda. Refuse intermediate lambda in the first consumer unless both branch ledgers have been captured.
3. Derive the effective endpoint from the continuation control when it is enabled. Even at a pure endpoint, the routine restores the caller's flags after using the selected branch; those flags alone do not establish which branch ran.
4. Use the existing isolated-output and evaluation-identity protocol. Preserve the complete exported terms at the required precision.
5. Validate the export against the assembled rows, then separately apply the physical stationarity criteria. Do not treat agreement with the production arithmetic as independent validation of every source model.

The endpoint point follows directly from `steady_residual.f90:163` through `:171` and `:209` through `:211`: lambda selects PLM or WENO3 before evaluation, and the saved flags are restored afterward. This is an implementation detail within the already required effective-configuration record.

Update phase 5's reference from sections 13.2-13.4 to sections 13.2-13.5 so that its execution gate includes continuation handling. No change to the intended physical model is implied.

## 5. Existing defects remain separate work

The `On stall` reachability problem remains visible in the inspected code. `EXHALE_main.f90:7355` resets the outcome, `:7356` skips the relaxation when the ending is no progress, and `:7876` requires the reset outcome to indicate a movement-bound ending. Rev6 correctly keeps that option out of comparisons until its repair is tested.

The inherited H3+ and coupled scaling-test findings were not rerun here. Their inclusion in an accepted plan does not mean their code repairs are complete. The charge assertion remains a validation improvement, not evidence that the neutral transfer changes charge.

The practical next step is execution of the initial diagnostics, followed by the scoped conservation export. Apply the wording corrections locally; the technical plan does not need another broad redesign before that work begins.
