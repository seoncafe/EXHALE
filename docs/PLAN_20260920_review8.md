# Review of PLAN_20260920_rev8.md

Review date: 2026-09-21 (KST).

## Decision

Accept rev8 for execution within its stated configuration restrictions. It incorporates the substantive requirements of review7: the generic endpoint mismatch is recorded, unsupported configurations are excluded, precision conversion is explicit, actual assembly selection must be recorded, and generic-double agreement must be measured.

No new blocking defect was found in the revised plan during this inspection. Three local clarifications remain: limit the statement about unaffected calculations to the default production assembly, distinguish rejection tests from tests of a complete repair, and avoid treating a skipped generic branch as an evaluation that did nothing. These do not require another complete planning cycle.

## Verification scope

I compared rev8 with rev7 and review7 and re-inspected the changed claims against the current dispatcher, generic wrappers, endpoint selection, momentum diagnostic construction, downstream scaling, build inclusion, and relevant existing test definitions.

Evidence below is source inspection. No atmospheric solution, benchmark, historical certificate, or compiled test was rerun. In particular, the numerical impact of defect 10.5 on the named checkpoints remains unmeasured. Only this review document was created; production source, input files, and solution products were not changed.

## 1. Findings from review7 are adequately addressed

| Review7 finding | Rev8 disposition | Current implementation evidence |
|---|---|---|
| Generic endpoint and momentum diagnostics can select different reconstructions | Correctly recorded as defect 10.5; section 13.6 excludes mismatching endpoints until addressed | `hydrodynamic_rows_body.inc:1464-1475` selects the endpoint explicitly, while `steady_residual.f90:263-266` invokes diagnostics that read global scheme flags |
| Quadruple assembly returns double arrays | Correctly distinguished from a fully quadruple coupled residual | `hydrodynamic_rows.f90:208-213`; separate conversions at `hydrodynamic_rows_body.inc:1504-1511` |
| A common residual routine does not guarantee a common assembly | Corrected to a dispatcher comparison | `steady_residual.f90:239-247` branches on the state kind and selected arithmetic |
| Generic guard has three refusal classes | All three are now listed | `hydrodynamic_rows.f90:289-309` |
| Face departures are not meaningful when the well-balanced option does not populate them | Correctly marked inapplicable | Assignments in `RK_rhs.f90:193-199` are conditional |
| Generic-double bitwise agreement is an intended property, not a new measurement | Correctly made an observed gate | Section 13.6 requests raw rows, faces, and scales separately |

Source paths in this table are under `src/modules/time_step/`, except where a full path is given elsewhere. The diagnostic implementation exists and is included by `Makefile:144`; accepting this plan does not mean its new experiments have already been executed.

The precision contract is now particularly useful: the internal working precision, precision of returned arrays, and text serialization precision are distinct. Preserve that distinction when reporting a small face-reconstruction discrepancy. Seventeen printed digits cannot restore information already lost when a working-precision face was converted to double.

## 2. Qualify the impact statement for defect 10.5

Section 10.5 concludes:

> None on the production atmosphere calculations, which never take this path; the defect is confined to the experiment of section 13.6.

The first clause is defensible only if "production atmosphere calculations" explicitly means evaluations selecting `ROWS_PRODUCTION`. The code does not confine the generic path to a dedicated test driver. `assemble_residual` invokes it for stationary evaluations when the selector requests it, and Newton calls that dispatcher.

The affected quantities also have downstream consumers beyond printed diagnostics. Source inspection shows:

- `steady_residual.f90:377` forms `momentum_largest_term` from the reconstructed terms.
- `steady_residual.f90:747` uses that value as the momentum scale.
- `steady_newton.f90:2775` reads `residual_row_scale` in forming model row scaling.
- `steady_newton.f90:5792` uses the row scale in a normalized measure.

Thus a stationary solve that enables generic arithmetic and reaches the mismatching endpoint configuration can consume inconsistent scales. The default production assembly avoids this particular defect, but the affected code is callable from a solver run, not exclusively from the newly planned audit.

Recommended replacement:

> The defect does not affect evaluations using `ROWS_PRODUCTION`. It can affect momentum diagnostics and their consumers in stationary evaluations that select generic arithmetic and satisfy the mismatching-endpoint condition. Its occurrence and numerical impact in the historical atmosphere runs have not been established.

Also replace "what moves is the scale" with "the scale can change." A different component formula does not guarantee a different maximum for every state: an unaffected term may dominate, or a special state may make the formulas coincide. This does not make the inconsistency acceptable; it simply avoids claiming an unmeasured effect size.

## 3. Separate the tests for restriction and for repair

Section 10.5 reasonably describes both a complete effective-scheme repair and a smaller guard extension. Its test paragraph should distinguish their expected outcomes.

| Implementation choice | Matching endpoint flags | Mismatching endpoint flags |
|---|---|---|
| Immediate guard extension | Evaluation succeeds and agrees with the corresponding production endpoint | Evaluation terminates with the specified unsupported-configuration result before the inconsistent diagnostics are formed |
| Effective-scheme repair | Evaluation succeeds and agrees with the corresponding production endpoint | Evaluation succeeds using the endpoint's scheme throughout; rows, diagnostic terms, and scales agree with the production endpoint |

For successful evaluations, check that the caller's configuration is restored. For an expected `error stop`, run the refusal check in a separate process and assert its exit status and diagnostic reason; there is no returned caller state to inspect.

The current test at `src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90:3513` disables continuation, and the loop near `:3559` sets matching PLM/WENO3 flags. It remains useful coverage but does not demonstrate either mismatching-endpoint behavior. This review inspected that test and did not rerun it.

The guard extension is a valid immediate restriction. It should not be described as full support for the previously inconsistent configuration. The complete repair is the consistent effective-scheme evaluation that the plan already recommends.

## 4. Clarify what happens when the generic selector is bypassed

Section 13.6 says that a mode-2 attempt remaining in the marching state kind "has exercised nothing at all." The implementation does perform a production residual evaluation; it has simply not exercised the requested generic branch.

At `steady_residual.f90:239-247`, `rows_kind` starts as `ROWS_PRODUCTION`, and the generic selector is consulted only outside the marching state kind. Replace that sentence with:

> A mode-2 request evaluated in the marching state kind exercises the production assembly, not the generic-double assembly. Record the generic comparison as not performed.

This distinction matters for interpreting a misleadingly identical pair: comparing production against production is not a test of generic-double equivalence. Rev8's requirement to record actual `rows_kind` and state kind already supplies the right prevention mechanism.

## 5. Recommended execution decision

Proceed with the initial production identity, restoration, repeatability, and closure diagnostics. For the first generic-double comparison, use continuation disabled or an endpoint with matching effective flags, require that the actual selected branch is generic double, and compare raw rows before faces and scales. Apply the existing isolation and provenance requirements.

Defect 10.5 remains a source-level finding, not an applied repair. The other previously listed defects likewise do not become fixed because the plan has been accepted. Their implementation and focused verification can proceed as separate tasks under the plan.

The remaining amendments concern accurate scope and test outcomes. Rev8's main physical and numerical investigation structure can now be used without another broad revision.
