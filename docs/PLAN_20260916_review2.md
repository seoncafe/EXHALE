# Review of PLAN_20260916_rev2.md

Date: 2026-09-16 (KST)

## 1. Overall judgment

**Revision 2 is suitable for starting bounded diagnostic work.** The principal objections to the previous designs have been addressed. In particular, L26 no longer prescribes an invalid zero-flux repair, L22 distinguishes pseudo-time iteration from physical subcycling, and L25 preserves the complete physical residual while separating accuracy from nonlinear convergence.

A further wholesale rewrite is unnecessary. Before implementing the affected parts, clarify five points:

1. A carrier relaxation that makes admissible progress and reaches its excursion bound is an incomplete inner solve, not automatically a failed outer solve.
2. Failure of the proposed mask experiment cannot exclude all interval-control designs.
3. Conservation tolerances must refer to the transported inventory in consistent units, not simply the wind-window mass flux.
4. The HLLC/Roe/corrected-Roe comparison must allow both the flux-family change and the Mach correction to contribute.
5. L9 must make the thermodynamic definition and certification of the chemistry-refreshed work state explicit.

These are targeted implementation-contract corrections. They do not invalidate the revised research program or justify delaying unrelated metadata, reproducibility, and data-inventory work.

## 2. What was verified

The complete rev2 plan was read, with focused source checks of:

- `src/EXHALE_main.f90`: carrier progress, outer-loop outcomes, and stationary evaluation.
- `src/modules/lower_atmosphere/diffusive_photochemistry.f90`: local intervals, cumulative bounds, rollback, and the fixed-conserved-state EOS calculation.
- `src/modules/states/base_boundary.f90`: the wind-window dependence and the circular dependency of a directly reused base flux.
- `src/modules/files_IO/load_IC.f90` and `src/modules/time_step/certification.f90`: imported versus evaluated certification.
- `src/modules/flux/Num_Fluxes.f90`: Roe acoustic coefficients.
- `src/modules/post_process/post_process_adv.f90`: input state and derived-output path.

Direct file checks gave:

```text
MD5
c2e9c9990b9f14f1be8cd77abca68945  EXHALE.x
8e818aa4699a2c8708e2da0f7ed77e84  src/EXHALE_main.f90
b445869ac32112a8498f1117c080af8e  src/modules/states/base_boundary.f90
ca4db5847c20101622a819c84f4673dc  src/modules/lower_atmosphere/diffusive_photochemistry.f90

SHA-256
c23d4274b7b2275e3727e06671983dc0a90346e66023215f4e7a3dab2c600ae1  docs/PLAN_20260916_rev2.md
```

The selected source identities have not changed since the previous review. This is a review of a revised plan, not evidence that its proposed corrections have already been implemented. No production executable, simulation, regression suite, or new numerical experiment was run for this review. The reported campaign and regression numbers remain values quoted from prior records.

No further literature search was needed: rev2 does not introduce a new paper-dependent numerical method beyond those checked previously against the supplied published Rieper papers. The remaining findings follow from the current implementation and mathematical contracts. The earlier boundary-limit probe remains relevant but was not rerun or modified.

## 3. Previous findings that can now be closed at the planning level

| Topic | Assessment of rev2 |
|---|---|
| L26 zero-window repair | Correctly withdrawn; multidirectional probes and a later user decision are appropriate |
| L26 fully static state | Correctly distinguishes coincident weight limits from the general discontinuity |
| Carrier progress | Correctly states that the default carrier measure already uses a steady residual |
| Element progress | Correctly separates an undamped finite-map displacement from a residual at the returned state |
| Local intervals | Correctly recognizes that `dt_code(j)` already exists |
| Mask definition | Empty-set handling and explicit masks improve the experiment |
| L25 damping sign | The negative-semidefinite operator and fixed-M qualification are now stated correctly |
| Low-Mach flux | Three-way comparison and separate accuracy/regularity/convergence measurements are appropriate |
| R15 energy model | Data inventory, prompt/internal partition, and unverified all-level quenching are correctly distinguished |
| L12 stage transport | Corrected eddy-flux decomposition and all-stage face identity are appropriate |
| L9 output provenance | Naming the loaded, refreshed, and advection-derived states is a useful correction |
| Metadata and fixtures | Imported claims remain distinct; the original difficult fixture is preserved |

“Closed at the planning level” means that the design objection is addressed. It does not mean that the corresponding physics or implementation has passed validation.

## 4. L22: distinguish incomplete inner relaxation from failure

Rev2 says that reaching the cumulative bound returns control to the hydro solve, but also says that a row which stops updating with a nonzero residual fails the pass. These statements need separate outcome definitions.

The current implementation saves the latest accepted state before each trial. A rejected trial restores that saved state; if the interval multiplier reaches its minimum, the routine returns a reason such as `carrier_relax_movement_bound` (`diffusive_photochemistry.f90:5958–5998`). Earlier accepted substeps are retained. Therefore a bound-limited return can be a useful, admissible partial advance even though the carrier is not stationary yet.

Requiring every inner relaxation to satisfy the carrier residual before another hydro update would defeat the purpose of bounding a fixed-hydro composition update. The composition may need a refreshed wind before further admissible progress is possible.

Use distinct outcomes, without weakening final certification:

| Inner outcome | Interpretation | Outer action |
|---|---|---|
| Residual satisfied | Carrier block closed at the current background | Continue the coupled check |
| Accepted progress, then excursion bound | Admissible partial update; not stationary | Return to hydro and reassess the joint residual |
| No admissible advance | No useful update was obtained | Report the exact rejection condition; apply a bounded recovery or stop policy |
| Nonfinite/invalid returned state | The update is unusable | Restore an admissible state and report failure |
| Repeated coupled stagnation | The complete alternating iteration is not improving under the tested controls | Stop or change the algorithm under an explicit policy |

The statement about a frozen row should be rewritten as: a nonzero physical row cannot pass final certification; an incomplete inner return may still advance the outer iteration. The new element residual should likewise be evaluated on the composition actually returned after damping and required state refreshes.

### Finite is not the same as non-NaN

The R23 verification-table row says that `carrier_close` takes `crc_max` when finite. The actual check at `EXHALE_main.f90:6794` is only `crc_max == crc_max`. Under ordinary IEEE semantics, this rejects NaN but accepts positive and negative infinity. Rev2 already proposes proper unavailable/nonfinite handling in L22 step 1; correct the evidence wording and implement an actual finiteness check. Also require nonnegativity for a norm. This is a verified guard limitation, not evidence that a current campaign produced infinity.

## 5. L22: define what the causal comparison can conclude

The following conclusion remains too strong: if the far-cell residual does not decrease in either selected experiment, the obstruction is elsewhere and no interval design answers it.

A negative result only rejects the particular mask, starting state, update length, and observation horizon tested. It does not exclude a local controller, a different redistribution of pseudo-time, or improvement that becomes visible only after the coupled hydro refresh. Likewise, a positive result may occur because a formerly limiting region was allowed to leave the small-response regime, not because the new strategy is suitable for production.

For each comparison record:

- The exact mask and the cells omitted from interval selection.
- Accepted intervals and all rejection reasons.
- Whether the original particle-count bound is still checked on omitted cells.
- Changes in pressure, temperature, source terms, and physical residuals throughout the domain.
- The residual immediately after the carrier update and after a complete coupled outer update, if that comparison is authorized.

If omitted cells still veto the same trial through the old global bound, the experiment may be unchanged despite a different proposed multiplier. If the bound is temporarily waived there, identify it explicitly as a diagnostic experiment and enforce positivity, element/charge feasibility, finite thermodynamics, and a separate safety stop. Excluding cells from interval selection is not a sufficiently complete description of trial acceptance.

Recommended conclusion wording: “No improvement was observed for the tested masks and intervals; test the named competing explanations before selecting a production design.” Keep the bound/controller as a leading hypothesis rather than a uniquely proven cause until an intervention distinguishes it.

## 6. L22: use an inventory-based conservation tolerance

Rev2 improves the physical-subcycling test by requiring volumes, integrated face transfers, and sources. However, its parenthetical returns to a face-flux mismatch relative to the wind-window mass flux. That is not a general scale for a species or elemental inventory:

- A number flux and a mass flux need a species mass conversion before comparison.
- An inventory defect and an instantaneous flux have different time dimensions.
- The distant net wind can vanish while internal diffusive exchange remains nonzero.
- A trace species can have a large relative conservation error hidden by the bulk mass scale.

For each conserved elemental inventory a, define over a common synchronization interval:

```text
E_a = (N_a,end - N_a,start)
      + outward_boundary_transfer_a
      - integrated_source_a.
```

All three terms must use the same units and geometry. Internal reaction sources for a conserved element should cancel by stoichiometry; included external sources should be stated explicitly. Require

```text
abs(E_a) <= absolute_tolerance_a + relative_tolerance_a * inventory_scale_a.
```

Construct the scale from inventories and/or absolute integrated transfers, with a declared zero-inventory policy. Supplement this with the face transfer cancellation test, rather than replacing the inventory test by it. The wind-window mass flux can remain a diagnostic but should not be the sole acceptance denominator.

For stationary pseudo-time iteration, rev2 correctly keeps final acceptance on the unmodified conservative residual and elemental/charge identities, without imposing a synchronized physical-time inventory equation on intermediate iterates.

## 7. L25: correct the attribution rule, retain the experiment

The three-way comparison is correct, but its attribution sentence is not: rev2 says that an improvement is attributed to the Mach factor only if corrected and unmodified Roe differ **and HLLC and unmodified Roe do not**.

Both changes can contribute. If HLLC performs poorly, Roe improves it, and corrected Roe improves it further, the latter difference is still the incremental effect of the correction under the tested Roe configuration. Equality between the first two is unnecessary.

Report two contrasts separately:

- Unmodified Roe minus HLLC: the flux-family change, including its existing entropy and admissibility behavior.
- Corrected Roe minus unmodified Roe: the incremental Mach-factor change under the same remaining settings.

Use the separate accuracy, Jacobian-action, and nonlinear-work measurements specified in rev2. Different trajectories can enter different basins; certification and grid checks still matter when comparing converged states. The acoustic coefficient locations in `Num_Fluxes.f90:324–331` support this experiment, but do not establish its outcome.

The remaining L25 changes are acceptable at the diagnostic-design level. A base-only mask is properly distinguished from Rieper's published factor. The discrete damping energy check has the correct negative sign and fixed-M condition. The plan no longer substitutes a hydrostatic-only certification row. Do not require a new flux to preserve the old finite-resolution root when physical consistency calls for a change.

## 8. L26: proceed with probes; preserve a single-valued residual

The revised L26 is now appropriately scoped. It recognizes the directional limit problem, does not call RMS an automatic repair, and puts the physical boundary decision after probes. No additional formula should be chosen merely to make the plan appear complete.

One implementation constraint should be explicit. `base_boundary.f90:148–157` explains that the base face flux is produced by the Riemann solve which consumes the boundary state. Reading the previously stored face flux to choose the new boundary can introduce evaluation-history dependence.

A local face-flux solution should therefore be a deterministic local boundary/flux solve, or another construction shown to be a function of the current state and reservoir alone. Its closure accuracy must be adequate for Jacobian differencing. Include repeated evaluations and changed call order in the boundary tests. A stale flux must not become a hidden additional state variable of the stationary residual.

This is a guard on the proposed local approach, not a rejection of it. The separation between probes and a later user decision on the physical repair is appropriate.

## 9. L9: specify the thermodynamics and meaning of the work-state certificate

Three named states improve provenance, but a chemistry-refreshed work state must also be thermodynamically defined.

The current evaluation path derives pressure from the loaded conserved state, calls `ioniz_eq`, refreshes particle densities, and recomputes temperature from pressure (`EXHALE_main.f90:5274–5293`). The molecular relaxation has a different explicit contract: `pressure_and_temperature_at_fixed_conserved_state` recomputes pressure using the current caloric EOS at fixed `u` (`diffusive_photochemistry.f90:5580–5603`). A composition change can alter heat capacity as well as particle count.

For the new work state, decide which quantities are held fixed:

- If conserved density, momentum, and total energy are retained, derive pressure and temperature from the refreshed composition and EOS, then recompute consistent source terms.
- If pressure or temperature is intentionally retained instead, reconstruct and label the resulting change in conserved energy. That is a different state, not a raw evaluation of the loaded one.

One chemistry sweep does not automatically produce a closed chemical/thermal fixed point. If the work state is not closed, report the failed closure or residuals; naming it “refreshed” must not imply success. This review does not establish a new numerical failure of the current evaluator; it identifies a necessary contract for the planned extension.

Retain two distinct answers where appropriate: whether the file's original stationary claim reproduces, and whether the explicitly refreshed work state passes its own evaluation. A successful refreshed state must not hide failure of the original claim.

Finally, the advection-derived product must not inherit the work state's certificate merely because it shares the output writer. `post_process_adv.f90:1807` calls `write_output` for derived quantities. Trace the metadata consumed by that call and mark the input certificate as provenance, not as certification of the derived composition. Rev2 states this separation correctly; a negative test should enforce it.

## 10. Other implementation notes and decisions

**L23:** The explicit metadata policy is adequate. Validate reason length before copying into a short destination; the current parser uses fixed-length tokens. Distinguish an omitted legacy field from an explicitly conflicting value when comparing the two files. Test conflicting claims and duplicate keys, and preserve unknown tokens only as provenance. No numerical model change is required for this item.

**L15:** The revised localization is appropriate. “Goldens: none” should describe expected scope, not forbid a justified single-thread numerical change if a shared-state defect also affects that path. Same-setting repetition, thread-count comparison, and exclusion of shared writes are complementary evidence. Instrument only the relevant path and preserve the first divergent checkpoint.

**L7g and L12:** The inventory-first molecular model and the corrected stage-flux identity are acceptable starting points. Implementation remains conditional on actual data coverage and explicit boundary/matching conditions. No new literature-derived coefficient should be inferred from the acceptance wording alone.

**L24 and bookkeeping:** Preserve the original fixture and numerical products as planned. No broad regression or data regeneration is warranted for path-only changes.

### Recommended authorization boundaries

| Work | Review decision |
|---|---|
| L23 metadata implementation | Ready after the small parser-policy clarifications |
| L15 localization | Ready for bounded diagnostics |
| L26 probes | Ready; physical repair remains a later decision |
| L22 measurements | Ready after trial-acceptance and outcome definitions are made explicit |
| L22 production solver change | Conditional on the measurements and conservation contract |
| L25 three-way experiment | Ready after correcting the attribution sentence |
| L7g inventory and L12 derivation | Ready; do not infer implementation completion |
| L9 implementation | Define fixed quantities and certification/output semantics first |

The next useful step is targeted measurement, not another expansive plan revision. The remaining corrections can be made locally in the affected sections without reopening the conclusions already settled.

## 11. Deliverable and limitations

Only this review document was created. Existing reviews, production sources, the executable, campaign products, and expected regression outputs were left untouched. No diagnostic code was created in this turn, so there is no new program to archive under `docs/audit_20260905/`.

The review verifies source behavior where stated and critiques proposed contracts analytically. It does not claim that the planned algorithms converge, that the catalog failures are solved, or that any unperformed test passed.
