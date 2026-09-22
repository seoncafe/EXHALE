# Review of PLAN_20260920_rev3.md

Review date: 2026-09-20 (KST).

## Assessment and scope

Rev3 is a sound basis for the investigation. Its corrections to operator identity, closure counts, diagnostic reachability, residual comparisons, and final-state certification should be retained. Initial identity and frozen-state diagnostics can proceed. Before implementing the conservation audit or interpreting the alternation experiment, correct the geometry/output contract and the meaning of the proposed map defect. Before the coupled trial, correct the experiment dependency order.

This review compares the plan and review2 with executable statements in the current source. Source inspection is labeled INSPECTED below; mathematical deductions are identified as deductions. Historical certificate values and timing claims in the plan were not remeasured. No solver campaign, compiled probe, or regression suite was run. This is a document review; the source, case inputs, and stored solutions were not modified. The remaining questions can be resolved from the implementation and elementary discrete identities; no new literature claim is introduced.

## 1. Corrections that rev3 gets right

- The H2 channel record always has 15 fields: 13 molecular channels followed by two oxygen chemistry channels. INSPECTED: `src/modules/lower_atmosphere/diffusive_photochemistry.f90:1300`, `:2841`, and `:2851`. Review2's wording "up to 15" was unnecessarily conditional. Rev3 correctly distinguishes a fixed output layout from whether the oxygen terms are active.
- Effective parsed settings belong in the experiment identity. INSPECTED: `src/modules/files_IO/input_read.f90:2193` derives `carrier_transport` from `thereis_oxychem` when the former was not stated.
- A residual-route comparison requires state restoration and explicit closure policies. The loaded diagnostic calls boundary preparation and residual assembly; `assemble_residual` also updates stored row terms. Comparing sequential evaluations without restoration can compare different states or workspaces.
- The hydrodynamic conservation export is now correctly described as work that must be implemented. INSPECTED: `src/EXHALE_main.f90:1299` exports three residuals and mass faces, while `src/modules/time_step/steady_residual.f90:355` retains mass faces and largest-term scales rather than a complete conservation record.
- The closure experiment now records requested/completed counts, seeds, admissibility, and refusal reasons. This is appropriate for separating evaluation uncertainty from changes in the eliminated chemistry map.
- A controlled evaluation of the final state is the correct source for H3+ domain diagnostics. Historical call counts cannot identify affected cells in that final state.

These items are sufficiently specified conceptually. They do not need another round of broad redesign.

## 2. Required correction: the carrier file does not export raw faces and geometry

Section 9.3 says that the carrier output contains four faces "with their areas and the cell volume." That wording is misleading if interpreted as separately available quantities.

INSPECTED: `carrier_row_terms_write`, `diffusive_photochemistry.f90:2784`, writes a main row, a `face` row, and an H2 `chan` row. The face header is:

```text
face cell carrier fdif_in fdif_out fadv_in fadv_out Frho_in Frho_out
```

There are no separate cell-volume or face-area fields. The writer explicitly describes the first four quantities as contributions that already include geometry. Their construction at `:5754` and `:5855` is:

```text
Kj = 1 / (cv(j) * R0cb)
sL = fa(j-1) * R0sq
sR = fa(j)   * R0sq
fdif_in  = -Kj * sL * Jf(j-1,ic)
fdif_out =  Kj * sR * Jf(j,ic)
```

Thus these are signed volumetric divergence contributions, not raw flux densities. The advective entries have the same divergence interpretation in the file header. Applying area/volume factors again would be dimensionally wrong. Summing the volumetric entries without volume weights would also fail to produce an integrated conservation balance on a nonuniform grid.

Recommended replacement for section 9.3:

> The file writes signed face contributions to each volumetric divergence, with area and volume factors already applied. It does not separately export the geometry or raw carrier face fluxes. These terms support local attribution; an integrated budget additionally requires exact cell volumes and a documented reconstruction of face transport rates.

For the new audit, export the exact geometry and a shared face array, preferably with one record for each face. State whether a quantity is a flux density, an area-integrated transport rate, or a contribution divided by cell volume. State the spherical normalization, including whether the common solid-angle factor is omitted. Verify that the outward contribution of cell j and the inward contribution of cell j+1 cancel after volume weighting.

This is a remaining documentation error, not evidence that the production carrier operator applies the geometry incorrectly.

## 3. Required correction: output precision must support the claimed budget

INSPECTED: the carrier writer uses `ES14.7` for terms, faces, and channels (`diffusive_photochemistry.f90:2829`, `:2835`, `:2841`). This provides eight significant decimal digits. The loaded hydrodynamic profile uses `ES16.8` (`EXHALE_main.f90:1305`), providing nine.

Deduction: subtracting large terms reconstructed from these rounded files can lose a much smaller physical residual. A discrepancy in that reconstructed balance can arise from serialization even when the in-memory arithmetic is consistent. A small printed residual does not restore the digits discarded from its constituent terms.

Use at least 17 significant decimal digits for binary64 audit values, such as `ES25.16E3`, or a lossless binary representation with a documented schema. Export geometry at matching precision. If inspecting the existing files, propagate their rounding uncertainty and label any cancellation below that uncertainty as unresolved.

Separate two checks:

1. Assembly consistency: independently reconstruct the residual from exported terms and compare it with the production residual.
2. Stationarity: determine whether that reconstructed residual satisfies the physical acceptance criterion.

A nonstationary checkpoint may pass the first and fail the second. Agreement between a writer and its consumer does not establish that the physical flux or source formula is correct. Rev3 already recognizes this limitation; the acceptance gates should say it explicitly.

## 4. Required correction: define the alternating map before using its defect

Section 9.2 asks for an "unconstrained map defect" and says that only this metric indicates whether a fixed point is near. That conclusion is too strong. Review2 introduced this recommendation without a sufficiently precise definition; it should be corrected rather than carried forward as an established requirement.

For a deterministic implemented update `Y_next = G(Y)`, the map defect `G(Y)-Y` is exactly the state displacement. Giving the two different names does not create an independent diagnostic. They differ only if the proposed G omits a constraint, damping, rejection, or part of the implemented update. In that case G is a different map and must be defined.

INSPECTED: the carrier relaxation rejects a trial, restores thermochemical state, reduces `grow`, and can end on `carrier_relax_movement_bound` (`diffusive_photochemistry.f90:8006` through `:8053`). A trial before that rejection is not automatically a complete alternating-map application. Removing its movement limit also does not imply that chemistry remains admissible or that the fixed-wind approximation remains appropriate.

Recommended experiment:

- Record the requested carrier trial, the accepted carrier change, and the reason for any reduction or rejection.
- Record the full accepted outer-state displacement and the difference spanning two outer passes, with fixed documented block scales.
- Evaluate the refreshed joint equation residual as the stationarity criterion.
- If an additional unconstrained map is evaluated, specify its inner stopping criteria, sweep policy, admissibility checks, and which constraint is removed. Run it from a restored snapshot and report an undefined or inadmissible result when it cannot be evaluated.

Deduction: even a genuine map defect does not bound distance to a solution without a suitable stability or contraction estimate. Locally, `G(Y)-Y` is approximately `(DG-I)(Y-Y*)`; a direction with an eigenvalue of DG close to one can have a small defect and a large state error. Therefore neither a small displacement nor a small map defect should replace the joint physical residual.

## 5. Charge: a missing assertion, not a demonstrated violation

Rev3 correctly finds that `nchg0` and `nchg1` are computed and not compared (`src/modules/init/molecular_seed_from_atomic_state.f90:450`, `:610`). However, its introduction calls this a further defect, and section 16 leaves open whether charge should be an invariant. The implementation allows a more precise conclusion.

INSPECTED: `transfer_h2` restores `f_atomic`, computes the neutral pool, and changes only H I and H2 (`molecular_seed_from_atomic_state.f90:704` through `:748`). The changed entries are:

```text
f_pool = f_HI + 2*f_H2
f_H2   = 0.5*dfc
f_HI   = f_pool - dfc
```

INSPECTED: both species have charge zero in `src/modules/init/species_table.f90:152`. The census sums species charge weights and metal ion stages (`src/modules/functions/element_census.f90:152` through `:214`); it does not subtract an independent electron density.

Deduction: for finite admissible inputs, this neutral redistribution preserves the census charge algebraically because all charged entries remain unchanged. The missing comparison is a validation gap; this review has not found a charge-changing operation in the transfer.

Recommendation: add the charge assertion when the conversion checks are next modified. Define it as preservation of the ion charge density before and immediately after the neutral transfer. If quasineutrality is also checked, compare the independently obtained electron density with that census and acknowledge when electron density is constructed from the same census, making the check an identity rather than an independent test. Use an absolute-plus-relative tolerance that remains defined for neutral cells.

Keep the subsequent boundary reconstruction separate. A reservoir may reset ghost composition; preserving each old ghost charge value across `Apply_BC` is not the invariant of the interior neutral transfer. The boundary must instead satisfy its own reservoir and charge closure.

## 6. Execution order and Mode R need small but necessary edits

Section 9.4 requires the coupled action to pass Stage B before the coupled trial. Section 15 schedules the coupled trial as phase 7 and Stage B/C/D as phase 9. The closure-count experiment in phase 4 does not automatically validate the coupled action, whose unknown set includes carrier variables.

Make coupled Stage B an explicit prerequisite immediately before phase 7. Perform it with the coupled configuration on an immutable checkpoint. Then run the matched coupled and uncoupled solves with the same physical starting state and refresh stage. Since this comparison changes both the unknown set and globalization, its result concerns the combined configuration; it cannot by itself attribute success to either change separately. An additional trust-region control is useful only if that attribution becomes necessary.

Section 3 defines Mode R as reload and measure only, whereas section 9.2 asks for sequences of outer passes and an application of the alternating map. Resolve this in the plan:

- Reading existing logs is a historical trajectory audit and retains the historical executable identity.
- A restored one-step probe can be a disposable diagnostic derived from an immutable checkpoint; record that it performed an update even if no state was published.
- Advancing a sequence under the current solver to study progress is a current-solver experiment and belongs with Mode C controls and budgets.

Mode C need not imply a long campaign or an already certified checkpoint. A short controlled failed solve still supplies useful evidence if its initial state, budget, and ending are preserved. Publication of a diagnostic checkpoint must remain distinct from certification of a solution.

Stage D and "the ledger" in phase 9 also need an explicit reference or a short definition. They are not defined elsewhere in rev3. The original plan has an accepted-step record in section 7, but also has stronger causal interpretations elsewhere that later reviews corrected. Import only the intended recording protocol, not the old interpretation table implicitly.

## 7. Conservation export: recommended design and limits

Choose section 13.2 option 1: a diagnostic export at assembly, followed by a separate consumer that performs the sums and comparisons. This gives a concrete implementation boundary and avoids reconstructing a subtly different operator from rounded checkpoint products.

The claim that "only the writer is missing" should be softened. INSPECTED: `assemble_residual` has `Smom` and combined `Sene`, while `viscous_conduction_sources` forms the latter as `vel*Smom + qv + Qc` (`src/modules/time_step/viscous_conduction.f90:609`). Separate viscous and conductive energy terms are local to that routine. Exporting them separately requires diagnostic capture or an extended diagnostic interface. Likewise, preserving the pressure/gravity decomposition requires the actual reconstruction convention, not simply labeling an aggregate source as gravity.

Specify the following before implementation:

- Exact numerical face fluxes, cell volumes, face areas, sources, residuals, physical units or conversion factors, and binary64 output precision.
- The state and evaluation identity, reconstruction/precision settings, and whether boundary preparation preceded the export.
- Integrated budgets for the physical cells, subdomains ending at the failing rows, and the face at each certification gate.
- A consistency tolerance derived from arithmetic and serialization uncertainty, separate from the physical stationarity tolerance.

The discrete energy convention should be stated explicitly. INSPECTED: `caloric_eos.f90:485` includes translational energy and H2 rovibrational excitation through `1.5*T + x2*u_rv/T0`. This expression is not an H2 binding-energy term measured relative to separated atoms. The plan's warning against adding dissociation energy twice is appropriate; spell out the energy zero when documenting the complete source budget. This inspection does not independently validate every chemical heating/cooling channel.

## 8. Status of the four previously reported defects

| Item | Current review | Appropriate action |
|---|---|---|
| `On stall` cannot fire | Confirmed from the reset at `EXHALE_main.f90:7355`, skipped update at `:7356`, and conjunction at `:7876` after `outer_no_progress` | Keep disabled in comparisons until the transition and snapshot tests pass. |
| H3+ endpoints counted outside the domain | Confirmed from the strict lower comparisons in `h3p_cooling.f90:180` and `:245` | Separate finite-input checks from strict interval excursions; preserve the valid endpoint evaluation. |
| Historical H3+ calls appear in the final certificate | `certification.f90:2069` reads accumulated counters; the values are informational and do not by themselves invalidate certification | Distinguish historical calls, unique affected cells, and final-state domain status. Do not imply that these counters currently reject a solution. |
| Exact equality after scaling and unscaling | Confirmed at `steady_newton.f90:6839` and `:6849` | Compare identical scaled quantities or use a stated rounding allowance; retain derivative checks. |

The `On stall` design should also identify one coherent entry state. The composition after relaxation and the wind used before relaxation are provenance records that may describe different thermodynamic stages. Preserve both if useful, but mark which complete state is the actual input to the coupled solver. Re-evaluating its residual must not silently combine fields from different snapshots.

## 9. Recommended disposition

Accept rev3's overall structure. Start the identity and frozen-state checks, and amend the relevant sections before their dependent experiments:

1. Correct the carrier face/geometry description and specify sufficient output precision.
2. Define any additional alternating map and retain the joint physical residual as the acceptance criterion.
3. Classify the unused charge census as a missing assertion and state the neutral-transfer invariant precisely.
4. Put coupled-action validation before the coupled solve; distinguish historical logs, disposable probes, and current solver trajectories.
5. Specify the conservation export, including the transport-source decomposition, and define Stage D explicitly.

One editorial correction is also warranted: section 13.2's "MEASURED against the current source" is source inspection and should be labeled INSPECTED. None of the findings in this review establishes that the atmosphere lacks a stationary solution, and none requires another broad planning cycle before the already specified initial diagnostics can begin.
