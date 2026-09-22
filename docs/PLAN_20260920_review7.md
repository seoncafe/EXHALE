# Review of PLAN_20260920_rev7.md

Review date: 2026-09-21 (KST).

## Decision and scope

Rev7 resolves the main document issues from review6. The corrected mode definitions, narrower continuation statement, and effective reconstruction identity are appropriate. Initial identity, restoration, and production-operator diagnostics can proceed.

The new section 13.6 is useful, but its generic-precision experiment needs additional conditions before it can serve as an acceptance gate. Source inspection identifies a conditional mismatch between the generic assembly's continuation endpoint and the scheme used to reconstruct its momentum diagnostics. The quadruple path also rounds rows and faces separately to double, which limits what a consumer can establish from the returned arrays. Finally, sharing `assemble_residual` does not prove that two callers selected the same assembly branch.

This review compares rev7 with rev6 and review6, inspects the newly cited paths and their downstream consumers, and checks the relevant existing test definitions and build inclusion. Source references are relative to `src/`. No solver campaign, compiled test, or historical numerical result was rerun. The endpoint finding below is established from control flow and formulas; its occurrence or numerical impact in the listed atmosphere cases has not been measured. No source code, input, or solution product was modified.

## 1. Changes that are correct

- The old Mode R "measure only" and Mode C "until a publishable checkpoint exists" definitions have actually been replaced. The distinction between disposable probes, retained trajectories, and stored evidence is now usable.
- `reconstruction_continuation_rhs` returns blended `dF` and `S` and retains blended momentum diagnostics. Rev7 correctly limits the missing information to the independent reconstruction from branch face data.
- Recomputing a numerical flux from blended reconstructed states does not generally reproduce the blend of the two numerical fluxes. Section 13.5 correctly requires the actual operator fluxes. Recomputing both original branches from an immutable state and then blending is another valid implementation, provided their complete configuration is reproduced.
- The assembly selector really is read on each call to `generic_precision_rows_selected`, and it is consulted by `assemble_residual` only outside the marching state kind. The plan is right to add this identity.
- The generic implementation is included in the build (`Makefile:144`), has callable double and quadruple wrappers, and is invoked by the stationary residual dispatcher. It is an existing implementation, not just a proposed diagnostic.

One qualification to section 13.5: `face_q_up` and `face_q_dn` contain the last WENO departures after a blended production evaluation only when the well-balanced branch actually writes them. With that option off, those arrays are unused and may retain earlier values or zeros. Mark them inapplicable rather than calling them current WENO departures unconditionally.

## 2. Conditional code defect: generic endpoint assembly and momentum diagnostics can select different schemes

Priority: resolve or exclude the triggering configuration before using generic mode 1 or 2 to judge momentum scales or normalized residuals.

The problem is in the new experiment's permitted endpoint configurations, not in the already refused open continuation interval.

### Evidence chain

1. `modules/time_step/hydrodynamic_rows_body.inc:1464` through `:1475` chooses the generic flux assembly from lambda when continuation is enabled: lambda at or below zero selects PLM, and lambda at or above one selects WENO3. The calls pass an explicit logical scheme selector; they do not change the global scheme flags.
2. `modules/time_step/hydrodynamic_rows.f90:289` through `:309` refuses intermediate lambda, low-Mach damping, and internally inconsistent scheme flags. It does not reject an otherwise consistent caller flag set that disagrees with the selected endpoint.
3. After the generic call, `modules/time_step/steady_residual.f90:263` through `:266` invokes `equilibrium_pressure_force_of_state` and `momentum_row_terms_of_state`.
4. Those routines select their formulas from global `use_plm` and `use_weno3`, not the generic call's effective endpoint. See `RK_rhs.f90:326` onward for the equilibrium force, `:450` onward for the pressure-gradient weight, and `:395` onward for the pressure/departure decomposition.
5. `steady_residual.f90:375` onward uses those terms to form `momentum_largest_term`, which is downstream diagnostic and scaling data.

For example, continuation enabled with lambda zero and a consistent global WENO3 flag set passes the guard. The generic assembly evaluates PLM, but subsequent momentum diagnostic reconstruction uses WENO3. The reverse mismatch is possible at lambda one with global PLM flags. The production continuation temporarily sets the effective flags while running `RK_rhs` and calculating its momentum terms; the generic wrapper does not follow that pattern for the subsequent diagnostic calls.

This is a conditional implementation inconsistency: the raw generic row can be formed with the intended endpoint while its pressure-gradient and gravity diagnostics use another scheme. It can affect the largest-term momentum scale and therefore normalized comparisons. It is not evidence that every generic evaluation is wrong, nor that the current atmosphere checkpoints trigger the condition.

### Recommended handling

For the first planned comparison, require either continuation disabled or caller flags consistent with the effective endpoint. Record that condition explicitly. Report other endpoint configurations as unsupported pending repair.

The proper repair is to derive the effective reconstruction once and use it consistently for flux assembly, equilibrium-force evaluation, and momentum-term reconstruction. Passing the effective scheme explicitly is a clear design; a scoped save/set/restore implementation must preserve the same invariant. Do not change the physical endpoint merely to make a diagnostic agree.

Add targeted checks for both endpoints with both matching and mismatching caller flags. Compare raw rows, pressure-gradient diagnostics, equilibrium-force diagnostics, and momentum scales against the production endpoint. Verify restoration of the caller configuration afterward.

Existing tests are relevant but do not establish coverage of this condition: the inspected well-balanced generic diagnostic test in `tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90:3513` disables continuation and explicitly sets matching PLM/WENO3 flags around `:3559`. This review inspected those tests; it did not run them.

## 3. Quadruple internal arithmetic is returned through double arrays

Priority: required when interpreting an audit of `EXHALE_RESID_QUAD=1`.

Section 13.6 correctly names the quadruple assembly, but it should also name the precision of each returned quantity.

INSPECTED: the wrapper `hydrodynamic_rows_in_quadruple_precision` in `hydrodynamic_rows.f90:208` declares its returned states, differences, sources, and face data as `real*8`. The generic implementation converts them separately at `hydrodynamic_rows_body.inc:1504` through `:1511`:

```text
dF8 = real(dF, kind(1.0d0))
S8  = real(S,  kind(1.0d0))
ff8 = real(ff, kind(1.0d0))
fp8 = real(fpr, kind(1.0d0))
```

The assembled difference is calculated before this rounding. A consumer that subtracts the separately rounded double face fluxes may obtain a different small residual. Seventeen decimal digits preserve those double values, but cannot recover the discarded working-precision face values.

Consequently:

- A face reconstruction mismatch at the double cancellation scale is not automatically an error in the quadruple assembly.
- Record internal arithmetic precision, returned-array precision, and serialization precision separately.
- For the initial production/double audit, the existing double precision contract is appropriate.
- For a genuinely working-precision audit, capture faces, geometry, sources, and row differences before conversion, and perform the consumer arithmetic at adequate precision. Otherwise propagate the conversion uncertainty and label smaller differences unresolved.

The generic routine offers optional `dF_at_working_kind`, but that alone is not a complete high-precision face/source ledger. Also, the surrounding chemistry, returned residual arrays, and subsequent heating/cooling subtraction remain double in the inspected route. Describe mode 1 as higher-precision hydrodynamic assembly, not an entirely quadruple coupled residual.

## 4. A shared entry point does not prove an identical evaluated operator

Section 5.1 says that a difference between the loaded-state and Newton routes "cannot be a different row assembly" because both call `assemble_residual`. That inference needs qualification, particularly now that section 13.6 documents a dispatcher inside that routine.

INSPECTED: `steady_residual.f90:239` through `:247` chooses production, generic double, or generic quadruple according to both the state kind and selector. The same environment value can still produce production rows in the marching state kind and generic rows in a stationary state kind. Reconstruction, boundary data, mixture state, and transport activation also remain inputs or module state read by the routine.

Recommended replacement:

> Both routes share the residual dispatcher. Compare their effective assembly kind and prepared inputs first. When those agree, differences should be investigated in preparation, closure, mutable state, or output transformation rather than attributed to separate dispatcher implementations.

At each evaluation, record the actual `rows_kind`, `ieq_sweep_state_kind`, effective reconstruction, transport activation, and relevant input identities. Recording only `EXHALE_RESID_QUAD` does not prove which branch executed. In particular, a supposed mode-2 experiment that remains in the marching state kind has not exercised generic double assembly.

## 5. Use generic double as a conditional comparison, not a presumed pass

The proposed mode-2 comparison is useful. The bitwise statement in the selector's comments is an intended invariant, and tests exist for it. It is not a measured result for the new checkpoint or proof of all supported configurations.

Section 13.6 should define success as observed agreement after matching the complete input and effective configuration, including the endpoint condition in section 2 above. Compare the unnormalized hydrodynamic row arrays first, then face data and diagnostic scales. If the raw rows agree but the scales differ, investigate the diagnostic reconstruction rather than classifying a chemistry or flux failure.

The guard has three refusal classes, not just the two listed in rev7: inconsistent reconstruction flags/method are also refused (`hydrodynamic_rows.f90:302` through `:309`). List that precondition. Preserve unsupported results as configuration outcomes; do not silently change damping or reconstruction in a comparison intended to preserve the operator.

A successful comparison supports implementation agreement at the sampled state. It cannot replace the conservation consumer or physical model validation because two implementations can share a formula error. Likewise, "one extra evaluation" is only an incremental budget when a matching production control and restoration already exist; record actual counts rather than treating this as a general runtime estimate.

## 6. Recommended disposition

Accept the main plan and proceed with production-operator identity and repeatability checks. Before the newly proposed generic-precision experiment:

1. Add the conditional endpoint restriction or repair the effective-scheme mismatch and test it.
2. Record the actual selected assembly and state kind, not only the requested environment value.
3. State the double-return limitation of quadruple assembly and adjust the conservation comparison accordingly.
4. Qualify the shared-dispatcher inference and make bitwise agreement a measured gate.

The previous momentum branch table, gravity identity, output provenance, and mode definitions need no further redesign. This review introduces a narrowly identified code issue and precision limits relevant to rev7's new experiment; it does not establish a failure of the default production atmosphere calculations.
