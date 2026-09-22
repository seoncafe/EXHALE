# Review of PLAN_20260920_rev4.md

Review date: 2026-09-20 (KST).

## Decision

Accept rev4 as the investigation framework, with the targeted corrections below. The previous review's principal recommendations are now present: correct carrier geometry semantics, adequate audit precision, a properly qualified map defect, neutral-transfer charge bookkeeping, a defined Stage D, and coupled-action validation before the coupled trial. Initial identity and frozen-state diagnostics can proceed without another comprehensive planning revision.

The most important remaining technical issue is the energy audit's treatment of gravity. The implementation places gravitational energy work inside the flux-difference assembly, while the explicit energy source is zero. An export that includes only the numerical energy faces and the array named `S` would omit that work. This is a specification gap in the proposed audit, not a newly demonstrated error in the production energy equation.

Other corrections concern inconsistent mode assignments, output provenance, and a historical numerical value presented next to a different current state.

## Verification scope

I read rev4, compared its changes with rev3 and review3, and inspected the relevant executable statements in the following files:

- `src/modules/time_step/RK_rhs.f90` and `hydrodynamic_rows_body.inc`;
- `src/modules/states/Source.f90`;
- `src/modules/time_step/steady_residual.f90`, `viscous_conduction.f90`, `steady_newton.f90`, and `certification.f90`;
- `src/modules/lower_atmosphere/diffusive_photochemistry.f90`;
- `src/modules/init/molecular_seed_from_atomic_state.f90`;
- `src/EXHALE_main.f90`.

INSPECTED below means current source inspection. DERIVED means an algebraic consequence of the inspected equations. The two logarithmic ratios in section 5 were calculated directly with `awk` from the values printed in the plan; they are not new measurements of the atmosphere. Historical certificates and solver timings were not rerun. No source code, case input, or solution product was changed, and no regression suite or solver campaign was executed.

## 1. Energy conservation must include the gravity term embedded in dF

Priority: required before implementing or accepting the hydrodynamic conservation audit.

Sections 11.3 and 13.2 correctly require the code's signs, units, energy zero, and source decomposition. However, specifying "energy faces plus gravitational source" is insufficient without identifying where gravity actually enters the energy row.

INSPECTED: `src/modules/states/Source.f90:44` returns all zeros under the well-balanced option. Outside that branch, `:49` and `:51` set `S(1)=0` and `S(3)=0`. This does not mean the energy equation lacks gravity.

INSPECTED: `src/modules/time_step/RK_rhs.f90:273` forms:

```text
dF3p = A_R*Fmass_R*(phi_R - phi_cell)
     - A_L*Fmass_L*(phi_L - phi_cell)
dF_E = (A_R*Fenergy_R - A_L*Fenergy_L + dF3p)/V
```

Here the potentials are `Gphi_i` at faces and `Gphi_c` at the cell center; V is `spherical_cell_volume(j)`. The generic arithmetic implementation contains the same term at `hydrodynamic_rows_body.inc:1413`. The stationary row then subtracts heating minus cooling and transport energy sources (`steady_residual.f90:269` through `:278`).

A consumer that reconstructs energy from `face_flux(3,:)`, `S(3,:)`, heating, cooling, and `Sene` alone will miss `dF3p/V`. Its mismatch would be an audit implementation error.

Recommended export contract:

1. Export the numerical mass and energy face fluxes, exact face geometry and cell volume, `Gphi_i`, and `Gphi_c`.
2. Export `dF3p/V` as a separately identified gravitational energy contribution, or reconstruct it from those quantities in the independent consumer.
3. State that it enters the residual with a positive sign inside `dF_E`; do not additionally subtract a second gravitational work term.
4. Export radiative/chemical heating and cooling and the transport contributions in the existing energy convention.

DERIVED: because the mass source is zero, a useful independent identity for a static prescribed potential is

```text
Rmass = (A_R*Fmass_R - A_L*Fmass_L)/V
Q = heat - cool + Sene

Renergy + phi_cell*Rmass
  = [A_R*(Fenergy_R + phi_R*Fmass_R)
     - A_L*(Fenergy_L + phi_L*Fmass_L)]/V - Q.
```

This follows directly by expanding `dF3p`. After volume weighting, internal faces telescope using the energy flux augmented by potential-energy transport. The `phi_cell*Rmass` term must be retained at the refusing checkpoints, where mass balance is not necessarily small. Dropping it would incorrectly use a stationary mass identity while auditing a nonstationary state.

Use this as an assembly-consistency check in addition to the gas-energy row check. It is not a claim that every physical source model has been validated. If a future model makes the prescribed potential time dependent, its energy identity must be reconsidered.

The momentum export needs the analogous implementation awareness already requested by rev4: WENO adds a face-pressure difference outside the stored momentum flux, and the well-balanced branch uses pressure departures (`RK_rhs.f90:235`, `:263`). Raw `face_flux(2,:)` and `S(2,:)` alone are not a complete momentum ledger in those branches.

## 2. Mode assignments remain internally inconsistent

Priority: correct before scheduling the M1 trajectory experiment.

Section 3 now explicitly assigns the M1 alternation audit to Mode C. This is correct when the experiment advances a sequence of outer passes. However:

- Section 15 still labels phase 6, "M1 alternation audit," as Mode R.
- Section 16 still says that everything in phases 1 through 6 and phase 8 is Mode R.
- Section 3 introduces "Three things are called Mode R" and then includes a current-solver trajectory that the same paragraph assigns to Mode C.
- The first definition of Mode C says to run until a publishable checkpoint exists, whereas the later definition correctly allows a short controlled failed solve with preserved evidence.

Use one consistent classification:

| Activity | Required label |
|---|---|
| Reading historical logs and stored products | Historical audit, with the historical executable identity |
| Current operator evaluation or a restored disposable update from an immutable checkpoint | Mode R, stating whether an update occurred |
| Advancing a sequence under the current solver | Mode C, including short diagnostic trajectories |

Change phase 6 to C for newly executed trajectories. If it only reads existing logs, label it a historical audit. Replace section 16's blanket claim accordingly. Define Mode C by the experiment performed, not by whether it successfully publishes a checkpoint.

The coupled Stage B prerequisite is now scheduled explicitly as phase 6b. That resolves the previous ordering problem and should be retained.

## 3. Output identity needs an operational implementation

Priority: required when running diagnostics in a directory that might contain prior output.

Section 4.2 already requires output existence, row count, state identity, and exit class. That is a sound requirement. The current carrier file does not itself supply enough identity to meet it automatically.

INSPECTED: `certification.f90:1492` writes `output/carrier_row_terms.txt` inside `carrier_rows_of_state`, before restoring its isolated module workspace. `carrier_row_terms_write` returns immediately if recording is disabled or its arrays are unallocated, and otherwise opens that fixed filename with `status='replace'` (`diffusive_photochemistry.f90:2788` through `:2790`). Its header describes "the last assembly" and has no generation or evaluation identifier.

Consequences:

- A failed new run may leave a valid-looking file from a previous run.
- Multiple certification evaluations may overwrite the same path.
- Restoring the numerical workspace does not restore or remove the file; the certification routine explicitly treats file output separately.

An existence check is therefore necessary but not sufficient. Implement the existing provenance requirement with an isolated output directory for each diagnostic invocation, preserving existing products. Archive the file with a sidecar identifying the invocation, input state, executable, evaluation stage, and exit status. For a run with multiple evaluations, either emit identifiers at the writer or capture the intended evaluation explicitly. Modification time alone is weak evidence of which state produced the rows.

If execution ends before the intended evaluation, record that outcome even if another evaluation produced a file. A timeout is an external termination, while the last completed inner solve may still be a separately identified observation. Neither should be promoted into a completed outer result.

## 4. Stage D is defined correctly, but instrumentation gaps should be explicit

The new section 7.2 appropriately separates outer passes from inner Newton iterations and requests observed and predicted reductions under fixed weights. Its hook schedule matches the inspected conditions at `steady_newton.f90:16911` through `:16923`: the element diagnostic uses the first two iterations and `elem_diag_third`; additivity uses iterations 1, 20, and `maxit_used`; the jump scan activates at iteration 1.

Those sparse hooks do not by themselves produce a complete record for every iteration. The plan should identify fields already available in logs and fields requiring additional diagnostic capture. Where an outer pass has no Krylov solve or a rejected state has no accepted step, record "not applicable" rather than a zero that might be read as convergence.

There are two simple cross-reference errors: the introduction and verdict-table item 8 say Stage D is defined in section 8. It is defined in section 7.2; section 8 remains closure bias. Fix those references without changing the protocol.

## 5. The M2 decade count mixes historical and current values

Priority: correct before reporting the severity of the current M2 state.

Section 11.1 prints a current mass measure of `1.041e0` against `3.4e-10`, then says the mass row is 10.01 decades above tolerance using `0.2254 / 2.2e-11`. These are different pairs of numbers.

Direct calculation from the displayed values gives:

| Values used | Ratio | Decimal logarithm of ratio |
|---|---:|---:|
| `1.041 / 3.4e-10` | approximately `3.0618e9` | approximately 9.49 |
| `0.2254 / 2.2e-11` | approximately `1.0245e10` | approximately 10.01 |

Calculation command, executed from the project directory:

```sh
awk 'BEGIN {a=1.041/3.4e-10; b=.2254/2.2e-11; printf "current %.12g %.12f\nhistorical %.12g %.12f\n",a,log(a)/log(10),b,log(b)/log(10)}'
```

The original review, `PLAN_20260920_review.md`, explicitly associates the second pair with the September 18 value. Rev4 should preserve that historical qualifier. State approximately 9.49 decades for the displayed current pair, subject to checking the unrounded certificate, and retain 10.01 only as the correction to the historical memo. This review recalculated the printed arithmetic; it did not remeasure the checkpoint or validate its current index selection.

This is a reporting error. Both pairs indicate a very large refusal, but their similarity in interpretation does not justify mixing identities.

## 6. Previous findings that should remain closed or explicitly deferred

The geometry, precision, and map-defect corrections are now satisfactory at the conceptual level. The carrier writer confirms that its signed face contributions already include geometry and that its output uses `ES14.7`; the new audit precision requirement is appropriate. The H2 channel count applies specifically to the H2 `chan` record, not to the new hydrodynamic export. Section 9.3 should replace the ambiguous phrase "That record is 15 fields" with "The H2 channel record contains 15 channel values."

The charge classification is also correct. The neutral transfer changes only H I and H2 (`molecular_seed_from_atomic_state.f90:735` through `:742`), while the checks at `:610` through `:633` omit the computed charge census. Preserve the distinction between an algebraically conserved quantity and an explicit assertion that has not been implemented. Do not reopen this as a demonstrated physical charge violation.

The `On stall` reachability defect remains present on direct inspection: the no-progress assignment at `EXHALE_main.f90:7323`, unconditional outcome reset at `:7355`, skipped relaxation at `:7356`, and required movement-bound outcome at `:7876` remain incompatible. Keep it out of solver comparisons until repaired and tested. Rev4's coherent entry-state requirement is appropriate.

This review did not rerun the inherited H3+ endpoint probe or the coupled floating-point equality test. Their repairs remain separately recorded work; adopting rev4 is not evidence that those repairs have been applied.

## Recommended next action

Proceed with the existing initial diagnostics. Before the affected phases, make four local amendments: include embedded gravitational energy work in the audit contract, correct the mode table, implement output provenance, and separate the historical M2 ratio from the current one. Select the assembly-point export proposed in section 13.2 as the default design, with a separate consumer performing the checks above. The remaining cross-reference and record-label corrections are editorial.

No additional broad planning cycle is necessary to begin the clearly specified work. The unresolved scientific question remains whether a state satisfies the current discrete equations and whether that model is physically applicable; solver refusal alone answers neither.
