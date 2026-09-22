# Review of PLAN_20260920_rev5.md

Review date: 2026-09-20 (KST).

## Decision and verification scope

Rev5 is suitable for starting the investigation. It incorporates the substantive corrections from review4: gravitational energy work, output provenance, the M1 trajectory classification, and the distinction between current and historical M2 ratios. No further comprehensive planning revision is needed before the initial diagnostics.

Before implementing the momentum conservation audit, make its reconstruction rules explicit and limit the internal-face cancellation test to quantities that are actually shared conservative fluxes. A few older sentences also remain inconsistent with the revised mode definitions. These are local corrections to the plan, not evidence of a newly demonstrated production-code failure.

I read all of rev5, compared it with rev4, and inspected the current residual assembly, numerical flux, reconstruction continuation, source, certification, output writer, and relevant solver paths. INSPECTED below denotes executable source inspection; DERIVED denotes algebraic consequences. Historical certificates, timings, and inherited probes were not rerun. No numerical campaign or regression suite was executed, and no source or solution data was modified. The only new artifact is this review.

## 1. Review4 findings that rev5 resolves

| Finding | Rev5 assessment |
|---|---|
| Gravitational work omitted from the proposed energy export | Resolved in sections 13.2 and 13.4. The export now includes potentials and the embedded energy term. |
| Current M1 trajectories labeled Mode R in the execution table | Resolved in phase 6: newly advanced sequences are Mode C; existing logs are historical evidence. |
| A fixed diagnostic filename cannot identify its generating evaluation | Resolved as a protocol in section 4.2: isolated output, an identifying sidecar, and explicit capture of the intended evaluation. Implementation remains future work. |
| Historical 10.01-decade ratio attached to current M2 rows | Resolved in section 11.1. The displayed current pair is associated with approximately 9.49 decades. These values remain document-derived here. |
| Sparse hooks presented beside a complete Stage D record | Resolved conceptually: section 7.2 requires identification of missing instrumentation and uses "not applicable" for undefined fields. |

INSPECTED: `certification.f90:1492` still writes the carrier file before restoring its workspace; `diffusive_photochemistry.f90:2788` through `:2790` still returns without writing when unavailable and otherwise replaces the fixed path. Thus rev5's output protocol addresses an actual implementation limitation. Its presence in the plan does not mean the diagnostic writer now supplies identifiers.

## 2. The new energy identity is correct

INSPECTED: `src/modules/time_step/RK_rhs.f90:273` through `:275` includes

```text
dF3p = A_R*Fmass_R*(phi_R-phi_cell)
     - A_L*Fmass_L*(phi_L-phi_cell)
dF_E = (A_R*Fenergy_R-A_L*Fenergy_L+dF3p)/V.
```

`src/modules/states/Source.f90:44` through `:51` makes the explicit mass and energy sources zero. `steady_residual.f90:269` through `:278` subtracts heating minus cooling and the transport energy source. The generic arithmetic implementation has the same gravitational energy construction at `hydrodynamic_rows_body.inc:1413`.

DERIVED: the identity in section 13.4 follows directly:

```text
Renergy + phi_cell*Rmass
  = divergence(Fenergy + phi_face*Fmass) - (heat-cool+Sene).
```

Keeping `phi_cell*Rmass` at a refusing state is essential. Rev5 correctly distinguishes this discrete assembly identity from validation of the physical source models. Retain the gas-energy check alongside the augmented-energy check; agreement in one combined equation should not hide compensating mistakes in the mass and energy rows.

One explanatory refinement: the quoted comment in `Source.f90` describes the cancellation of equilibrium pressure and gravity in the momentum equation. The executable `dF3p` statements establish how gravity enters energy. State these as separate reasons rather than using the momentum comment as the explanation for every zero source component.

## 3. Specify all momentum branches in the audit

Priority: required for the momentum audit consumer; does not block frozen-state identity or repeatability checks.

Section 13.4 discusses the well-balanced `face_q` expressions, but it should also explicitly include ordinary WENO3. INSPECTED: `RK_rhs.f90:235` adds `(pR-pL)/dr` when WENO3 is active, even when the well-balanced option is disabled. In that case the required pressure information is `face_p`, not `face_q_up` or `face_q_dn`.

The following table is a direct transcription of the relevant assembly, with A denoting face area, V cell volume, and F the stored numerical momentum flux. The final momentum residual subtracts `Smom` in every branch where transport is active.

| Reconstruction and balance option | dF_momentum | Explicit S_momentum |
|---|---|---|
| PLM, ordinary | `(A_R*F_R-A_L*F_L)/V` | Gravity plus the geometric pressure term from `Source.f90` |
| WENO3, ordinary | `(A_R*F_R-A_L*F_L)/V + (p_R-p_L)/dr` | Gravity from `Source.f90` |
| PLM, well-balanced | `[A_R*(F_R+q_up,R)-A_L*(F_L+q_dn,L)]/V` | Zero |
| WENO3, well-balanced | `(A_R*F_R-A_L*F_L)/V + (q_up,R-q_dn,L)/dr` | Zero |

Evidence: `RK_rhs.f90:231` through `:270`, `Source.f90:44` through `:54`, and `steady_residual.f90:277`. The ordinary WENO pressure term also appears in `hydrodynamic_rows_body.inc:1380`.

Export the actual branch, `dr_j`, the necessary pressure/departure arrays, geometry, explicit source, and viscous momentum term. Do not add an explicit gravity source to the well-balanced row a second time. Conversely, do not treat ordinary WENO pressure as already contained in the stored momentum flux.

### The internal-face cancellation test needs a narrower statement

Section 13.2 says to verify cancellation between the outward contribution of one cell and the inward contribution of its neighbor. This is appropriate for a shared conservative transport rate. It does not apply indiscriminately to all pressure-departure terms.

INSPECTED: `src/modules/flux/Num_Fluxes.f90:167` through `:168` explicitly constructs `q_dn = q_up - dp_eq`. The two quantities at one face refer to different cells' equilibrium pressures and need not be equal. Other numerical-flux branches also return these distinct departures. Their failure to cancel alone is not evidence of broken conservation.

Recommended wording:

> Verify internal-face cancellation for shared conservative fluxes after the specified volume weighting. Verify momentum pressure, geometry, and gravity terms through the complete branch-specific discrete identity; do not require the two equilibrium-referenced pressure departures at a face to be equal.

This also corrects an overly broad reading of the cancellation recommendation in earlier reviews. Spherical radial momentum with gravity is a balance law with source terms; a zero source-free global sum is not its acceptance condition.

## 4. Reconstruction continuation must be represented explicitly

Priority: conditional on auditing a state with `0 < recon_lambda < 1`.

INSPECTED: `reconstruction_continuation_rhs` in `steady_residual.f90:117` evaluates PLM and WENO3 separately and blends their `dF`, `S`, face fluxes, and stored momentum terms. At `:209` it restores the caller's reconstruction flags. Thus the flags seen after the call need not describe one pure reconstruction that can reproduce the blended operator.

The existing code explicitly blends `momentum_ram_divergence`, `momentum_pressure_gradient`, and `momentum_gravity`. It does not construct an analogous blended pair of `face_q_up` and `face_q_dn` in that routine. Exporting the final departure arrays and a restored method flag is therefore insufficient to reconstruct the momentum row in the open continuation interval.

For this case, record lambda and capture the two branch ledgers before blending, then verify their weighted sum. Alternatively, state that the first audit supports only the pure endpoints and reject an intermediate lambda explicitly. That restriction should be reported as an unsupported diagnostic configuration, not as failure of the physical equations.

The energy identity remains algebraically compatible with blending at a fixed state and potential because its flux expression is linear in the blended mass and energy fluxes. This observation does not supply the missing momentum decomposition.

## 5. Mode definitions retain a few obsolete sentences

The revised classification table is usable, but the paragraphs around it should be normalized:

- Section 3 still opens Mode C with "until a publishable checkpoint exists." Remove that stopping condition; later text correctly permits a short failed trajectory with preserved evidence.
- The first Mode R paragraph says "measure only," while the table allows a restored disposable update. Refer explicitly to that limited probe rather than leaving the two descriptions contradictory.
- Section 16 assigns phase 6 to Mode C "as soon as it advances a pass." Clarify that this means advancing a retained trajectory, since the table allows a single restored disposable update in Mode R.
- Stage A0 includes inspecting existing manifests and logs. Label that portion a historical audit; only its new operator evaluations are Mode R. Storing newly generated diagnostic results in a log does not itself make the underlying experiment a historical-operator test.

These are operational definitions, not new physics requirements. Use the classification table as the authority and amend the older prose accordingly.

## 6. Implementation status and next action

The previously identified `On stall` defect remains present: the outcome reset at `EXHALE_main.f90:7355` and skipped update at `:7356` conflict with the movement-bound outcome required by the handover condition at `:7876`. Rev5 correctly excludes that option until repair and verification. The review did not rerun the H3+ endpoint probe or coupled scaling test, so their recorded repairs remain unverified work.

I recommend accepting rev5 and beginning its initial diagnostics. Before building the conservation consumer, add the momentum branch table, qualify the cancellation test, and specify treatment of intermediate reconstruction continuation. Correct the remaining mode wording locally. These additions make the current contract executable without requiring a broader solver redesign or another complete planning cycle.
