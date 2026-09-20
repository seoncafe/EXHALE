# Second review of PLAN_20260918: revision 1

Date: September 18, 2026 (KST)

Reviewed document: [PLAN_20260918_rev1.md](PLAN_20260918_rev1.md).
Previous review: [PLAN_20260918_review.md](PLAN_20260918_review.md).
The filename of this second review is `PLAN_20260918_review1.md`, as requested.

## 1. Executive judgment

Revision 1 substantially improves the plan. It withdraws the incorrect characteristic argument, preserves direct evaluation of returned states, distinguishes chemical turnover from transport accuracy, includes the helium inventory defect, and expands D5 into a coupled boundary-state problem. Those are meaningful corrections, not merely wording changes.

However, the new plan should not be approved as one undifferentiated implementation package. D7 exposes an additional defect in an existing reference implementation; D8 does not yet define a safe checkpoint transaction or a state-consistent ranking; and D9 turns one recorded convergence history into an unsupported universal pass-budget rule. The molecular seed restrictions added to D5 are also stronger than either conservation or the loader actually requires.

Most importantly, an isolated execution of the unchanged production charge-exchange module confirms that its generic helium contribution has the wrong orientation when inserted directly into the He I balance used by `System_HeH_TR_metals`. Reusing this routine without an explicit equation-basis conversion would reproduce that error in the transported sources. Agreement between two implementations is not sufficient when both use the same wrong sign.

| Item | Judgment on revision 1 | Remaining requirement |
|---|---|---|
| D1 | Direction approved with small corrections. | Use a safe numeric format; pin historical inputs before changing a global default; scope validation to affected grids. |
| D2 | Previous blocking premise corrected. | Distinguish interior, characteristic-face, and numerical-flux velocities; present concrete boundary-model choices before requesting a decision. |
| D3 | Diagnostic direction approved. | Include the new atomic branch and all source corrections; separate measurement from enabling a new certification gate. |
| D4 | Revised contract investigation approved. | Evaluate the unmodified returned state with isolated diagnostic state; do not equate a chemical residual with abundance error. |
| D5 | Core design direction approved; new seed restriction rejected as written. | Separate restart compatibility from construction of a new molecular seed; clarify the dependency on D2. |
| D6 | Conservation correction approved. | Include the excited-level admissibility condition and a clearly owned frozen molecular inventory. |
| D7 | Missing-source diagnosis supported; implementation prescription incomplete. | Add the interim refusal, then correct row orientations, dimensions, and thread-local rate initialization before enabling metal-bearing stage transport. |
| D8 | Important, but not implementation-ready. | Distinguish stale generations from mismatched file halves; define immutable generations, same-state certification, and crash-safe publication. |
| D9 | More time may be justified; the quantitative rule is not. | Replace the 3.3-cells/40-passes rule with measured, case-specific budgets and honest incomplete outcomes. |

Recommended immediate work is D7's unsupported-configuration guard and focused sign tests, D6's admissibility tests/correction, D5a/D4a diagnostics, and preservation of run generations. Do not refresh reference outputs merely to match an implementation whose source signs and state identity are still under review.

## 2. Evidence and validation scope

This review reread all 743 lines of revision 1, the relevant L34c/L36e records, and the actual source paths discussed below. The code has changed since the previous review: the carrier source now contains the L36e atomic branch and explicit H/He charge exchange. Earlier line numbers and earlier source hashes were not assumed to describe the new tree.

Three types of evidence are distinguished:

- **Executed production-module check:** a small program compiled the actual `parameters.f90`, `species_table.f90`, and `charge_exchange.f90`, and called their existing routines. It did not compile or run the full atmosphere solver.
- **Executed workflow check:** existing Python functions were called with controlled mocked inputs, or the actual nested mapper selector was extracted and called. The campaign pipeline check is a deliberately reduced shell example, not a campaign execution.
- **Recorded campaign result:** convergence numbers from L34c/L36e were read and checked for consistency, not remeasured in a new atmospheric run.

The retained checks are:

- [plan_20260918_rev1_source_probe.f90](audit_20260905/plan_20260918_rev1_source_probe.f90)
- [plan_20260918_rev1_workflow_checks.py](audit_20260905/plan_20260918_rev1_workflow_checks.py)
- [run_plan_20260918_rev1_review_checks.sh](audit_20260905/run_plan_20260918_rev1_review_checks.sh)
- [plan_20260918_rev1_review_checks.log](audit_20260905/plan_20260918_rev1_review_checks.log)

The script builds into a new `/tmp` directory and does not replace `EXHALE.x` or the project's normal objects. All checks completed with GNU Fortran 16.2.0. Production source, catalog inputs/products, and regression references were not modified.

## 3. D7: source completeness is necessary, but the equation orientation is also wrong in an existing caller

### 3.1 The missing metal contribution is confirmed

The current `carrier_source` starts near `diffusive_photochemistry.f90:4716`. It calls `mol_heh_rows` for molecular gas or `heh_tr_rows` for atomic gas, and now adds `he_h_cx_fvec` with a branch-dependent helium sign. It still does not call the generic metal charge-exchange assembly or an equivalent physical source evaluator.

The local metal systems do call `cx_add_to_fvec`, including `System_HeH_metals.f90:162`, `System_HeH_TR_metals.f90:125`, and `System_HeH_mol_metals.f90:259`. The transported and locally evaluated source equations therefore differ for metal-bearing configurations when the corresponding reactions are active.

The proposed interim refusal of `Ionization transport: True` with metals is justified. It is **not implemented yet** in the inspected `input_read.f90:2189–2263`: that block checks helium availability, molecular carrier transport, and incompatibility with the coupled route, but does not reject metals.

Make this refusal a distinct, small task rather than leaving it implicit in the eventual source addition. A later precise capability check could be narrower, but an unsupported combination must not silently execute meanwhile.

### 3.2 Correct the description of reaction groups

The D7 statement that groups C and D enter helium rows 2 and 3 is not consistent with `charge_exchange.f90`'s reaction descriptors:

- Group A exchanges charge between metals and H/H+.
- Group B is the separately evaluated H/He pair; it is excluded from the generic active set.
- Group C exchanges charge between metals and He/He+. Its listed helium transitions are neutral/singly ionized, not He II/He III transitions.
- Group D exchanges charge between metals. It does not directly add an H or He source, although its effect on metal populations can matter indirectly.
- The additional group E reaction, O III + H I to O II + H II, contributes to the proton source and is active through its own setting independently of `cx_full`.

D7 should include the actual active reaction set, especially the default group E contribution, rather than restricting its checklist to the plan's group description.

The published [Huang et al. (2023), ApJ 951, 123](https://doi.org/10.3847/1538-4357/accd5e), section 2.4.2 and Table 4, was checked in the local journal PDF `references/Huang_2023_ApJ_951_123.pdf` at the workspace level. It treats charge exchange alongside ionization/recombination and lists the Si + He+ reaction used below. The defect demonstrated here concerns stoichiometry and equation orientation, not a proposed change to the tabulated rate coefficient.

### 3.3 New confirmed finding: the generic helium contribution has the wrong sign in the triplet-system He I row

`cx_add_to_fvec` at `charge_exchange.f90:337–368` adds a positive reaction rate to the donor's lower-to-upper transition row and a negative rate to the acceptor's transition row. `cx_fvidx` maps the neutral/singly ionized helium transition directly to row 2. This is an ionization-positive convention.

However, `heh_tr_rows` at `ion_residual_core.f90:143–153` writes row 2 as the **He I gain** balance. A reaction destroying He II and creating He I must add a positive contribution to this row. `System_HeH_TR_metals` inserts `cx_add_to_fvec` directly after `heh_tr_rows`, without changing that generic contribution's sign. The dedicated H/He pair already has a `he_row_sign=-1` argument for this reason; the generic metal contribution has no equivalent conversion at this call site.

Consider only

\[
\mathrm{Si\ I}+\mathrm{He\ II}\longrightarrow
\mathrm{Si\ II}+\mathrm{He\ I}.
\]

At `T=10000 K`, with both reactant number densities set to `1 cm^-3` and the other reactants zero, the executed production-module check gives:

| Quantity | Measured value, cm^-3 s^-1 | Required physical sign |
|---|---:|---|
| Si II production from the generic routine | `+7.561e-11` | Positive |
| Generic increment to row 2 | `-7.561e-11` | Negative only for an ionization-positive helium row |
| Required increment to the TR He I balance | `+7.561e-11` | Positive |
| He II source inferred using the TR conversion `-f(2)-f(3)` | `+7.561e-11` | Must instead be negative |
| Total ionic-charge source under that inconsistent interpretation | `+1.5122e-10` | Must be zero for this exchange |

The kernel outputs were executed; the TR interpretation follows the inspected caller and the explicit source transformation. A complete atmospheric solution or the full nonlinear TR system was not executed. Thus this confirms the local source-assembly inconsistency, not its quantitative effect on an observed spectrum or an entire wind solution.

The affected existing route is not limited to `Ionization transport: True`. A local atomic TR+metals system with the relevant group C reactions enabled has this sign mismatch as well. Blocking transported metals alone does not correct that local-system defect.

### 3.4 Do not append the generic call to the present small array

There are two further implementation hazards:

1. `carrier_source` currently declares `fv(10)`. The generic routine also writes metal rows using `cx_metal_base + 2*(element-1) + stage`. A full metal layout exceeds ten rows. Calling it with the present small array would risk an out-of-bounds write, even for a reaction whose rate is zero, because the routine still addresses the indexed rows.
2. `cx_kc` and `cx_metal_base` are thread-local state. `cx_set_cell(T)` must initialize the rates for the actual current cell and temperature. `carrier_source` is called inside an OpenMP loop near `diffusive_photochemistry.f90:5280`; inheriting whichever rates a worker last used in another sweep would give a schedule-dependent source. Thread-local storage prevents a race but does not guarantee that stored values belong to the current cell.

Consequently, D7 needs explicit ownership of the reaction coefficients, metal-stage densities, row dimensions, and equation basis. “The same routine” is not a complete implementation specification.

### 3.5 Preferred repair and acceptance tests

Prefer a physically named reaction-source interface that computes stoichiometric number-density source increments for H and He stages and the affected metal stages. Each solver then maps those physical sources to its own residual equations. Alternatively, retain the existing transition-row evaluator but accumulate into a correctly sized separate array and perform an explicit, tested equation-basis transformation before adding it.

Do not negate a complete assembled row after the fact: that would also negate unrelated reaction terms. Do not negate the generic helium row globally, since its current orientation is appropriate for other equation bases.

Required tests:

- One active reaction at a time, in both available directions: species creation/destruction signs, each element's nucleus balance, and charge conservation.
- Standard atomic, atomic TR, and molecular equation bases, including zero-abundance limits.
- Group B counted once, group E included under its own switch, and groups C/D correctly selected by `cx_full`.
- Equality between the physical sources represented by each solver, not equality of raw residual arrays whose bases differ.
- Finite-difference derivatives and analytic derivatives where used, after the same basis transformation.
- Several cells at different temperatures evaluated in changing order and with one versus several OpenMP threads. Rate initialization must be independent of the previous cell on a worker.
- Bounds-checking builds for the enlarged/source-specific arrays.

The changes may correctly move key-off TR+metals results. Preserving an existing wrong sign is not an acceptable byte-identity requirement. This also means that the eventual reference refresh may involve more than D9's four key-on cases.

## 4. D8: distinguish run provenance, physical snapshots, and exact continuation

### 4.1 The runner's continuation policy is indeed too broad

`run_case.sh:398–428` launches the high-pseudo-time continuation whenever `INFO != 0` and both output files exist. The check for hydrodynamic rows being within tolerance occurs afterward, near lines 440–445, and only influences trying another seed.

Therefore, a composition-only refusal can already have triggered the continuation before the script identifies that condition. Moving the classification before the continuation is justified.

However, “hydrodynamic refusal” alone is still too broad. Missing/nonfinite state, a failed physical-admissibility condition, exhausted outer passes, and a flat hydrodynamic solve are not equivalent failures. Use an explicit ending reason and validated state availability, not just a nonzero aggregate `info`. Continue only for endings the chosen continuation is intended to address.

The L34c change from `4.13e-2` to `7.37e-2` is a **recorded result**. It supports preserving both endpoints. It does not establish that every composition-only continuation is harmful, or that a smaller single printed row defines the uniquely best state.

### 4.2 A seed and a product are not the two halves of one pair

The two halves of an output state are `Hydro_ioniz.txt` and `Ion_species.txt`. The input pair uses the corresponding two `_IC.txt` filenames. An old product beside a new input seed is a stale-generation problem, but it does not by itself prove that either individual hydro/species pair has mismatched halves.

L34c describes the old product/new seed situation. Revision 1 should not present this as a direct measurement of a partially written hydro/species pair.

A real partial-write risk nevertheless exists independently: `write_output.f90:239–327` writes and closes the hydro file before opening the species file. An interruption between these operations can leave mixed output generations. The revised plan should address both problems and test them separately.

### 4.3 There is an existing reader path that can deliberately mix suffixes

`map_state_to_grid.py:520–523` chooses the source filename independently for each component: use `name.txt` if it exists, otherwise use `name_IC.txt`.

The executed check of that actual nested selector, with only the hydro product present, selects:

```text
Hydro_ioniz.txt
Ion_species_IC.txt
```

This proves that filename selection can combine a product and a seed. Subsequent radius or reservoir checks can detect some inconsistencies, but are not a generation-identity check and need not reject two states on the same grid and reservoir.

Resolve a complete generation first, then read both components from it. If the selected product generation is incomplete, either reject it or explicitly select a complete checkpoint generation; do not choose each component independently.

### 4.4 “Best state” needs a same-state measure and a non-destructive policy

Define what is being ranked. A useful diagnostic ordering can prioritize physical admissibility and complete evaluability, then certification, then a named normalized distance such as the largest active residual/tolerance ratio. Preserve the full vector of row measures, scales, gates, and model provenance beside that score. A state closest under one diagnostic norm is not necessarily the best basin for a subsequent nonlinear solve.

Keep the latest complete checkpoint, the best measured candidate, and every certified generation distinct. Selecting a best product for reporting should not silently reset the solver trajectory to that product.

There is a specific current-code trap. In the outer loop:

1. Certification is computed near `EXHALE_main.f90:6858`.
2. Composition can subsequently change in `relax_photochemical_composition` near line 7214 and in the element update.
3. The summary near line 7418 prints the earlier `sp_meas`, `sp_cell`, and related certification measures.

A checkpoint written after the update cannot simply inherit the earlier certificate or use that score to declare itself best. Either checkpoint the exact state that was evaluated, or evaluate the completed updated state before attaching its score. Store the state identifier in both the certificate and the snapshot. This is directly related to D5 and prevents D8 from repeating the same different-state comparison in a new file layout.

### 4.5 Two filenames cannot be made one atomic transaction by separate renames

The requirement to write “the state pair atomically” needs an actual publication protocol. A safe design is:

1. Write a new generation directory with both data files and a manifest. Do not alter the currently published generation.
2. Include a unique run/attempt identifier, iteration identity, source/binary/configuration identifiers, file hashes, dimensions, completion status, and the certificate's state identifier.
3. Close and validate all files. Apply the filesystem synchronization required by the intended durability guarantee.
4. Publish one pointer or manifest to that completed generation using an atomic operation within the same filesystem.
5. Readers resolve that pointer once and use the selected generation for both components. Preserve the previous completed generation until the new one is usable.

Two separate renames leave an interruption window. Matching pass numbers alone are insufficient: a restarted run can repeat the same pass number, or a file can be truncated while retaining its header. A completion manifest and content validation are needed.

Because this workspace is on NFS, distinguish process-interruption safety from durability across a server failure. State which guarantee is implemented and test it on the relevant filesystem; do not infer both from a successful local rename.

This is a data-layout and reader-interface change. Its intended shape should be agreed before implementation, consistent with the plan's approval rule for D5b.

### 4.6 A physical-state checkpoint is not an exact solver checkpoint

`trust_pass` is initialized from `carrier_trust` near `EXHALE_main.f90:6756` and can shrink during the solve. Restarting only from the hydro/species pair resets this and other numerical progress history. Saving the pair after every pass does not by itself prevent the changed trajectory observed when the continuation starts again with a larger movement bound.

Decide explicitly between:

- **Warm restart:** preserve the physical state, rebuild derived quantities, and restart the solver controls under a documented policy.
- **Exact algorithmic continuation:** also preserve the relevant trust bound, relaxation parameters, progress history, reconstruction state, iteration counters, and other state that changes the iteration map. Verify against an uninterrupted control.

The first is a valid and simpler contract; it must not be described as exact continuation. Neither implies physical time evolution: these are stationary solver iterates.

Do not call the complete `write_output` path at every outer pass without checking its cost and state dependencies. It emits additional diagnostics and products, not merely a compact restart pair. Prefer a focused checkpoint writer whose inputs are a single synchronized snapshot and whose execution does not change subsequent solver behavior.

### 4.7 Reader and campaign completion rules are part of D8

Two additional checks found relevant current behavior:

- `pick_seed.py::certified` trusts a `REPRODUCE.md` certification line before reading a log, without binding that verdict to hashes or a generation identifier of the current data files. The controlled function call accepts such a record. This is not proof that a particular catalog seed is mislabeled, but it is an interface capable of trusting a stale certificate.
- `run_campaign.sh` captures `run_case.sh` through `tail -n 1`, then prints a line from a shell without `pipefail`. A reduced execution of the same command structure returns exit status zero after the task fails. The campaign's final message therefore cannot serve as a machine-readable success verdict. Its last-line capture also normally keeps the re-measurement sentence rather than the preceding detailed `DONE` line.

In addition, `run_case.sh` temporarily edits `input.inp` for evaluation and restores it only after the binary returns. D8's interruption tests should cover this phase. Prefer a separate evaluation directory/configuration so a killed runner cannot leave the main solve input transformed into an evaluation-only input.

Add `pick_seed.py`, `status.py`, `write_reproduce.py`, the mapper, and campaign outcome handling to D8's ownership/consumer list. The writer and runner alone do not establish a reliable run-directory contract.

### 4.8 Required D8 tests

- Composition-only budget exhaustion, supported hydrodynamic continuation, nonfinite refusal, and missing final verdict are classified separately.
- Interrupt before any output, during the hydro write, between the two component writes, during manifest publication, and during post-processing.
- After each interruption, a reader obtains either the previous completed generation or the new completed generation, never a mixed pair.
- An old certificate with new files, identical pass numbers from different attempts, and a truncated file are rejected.
- The certificate used for best-state ranking belongs to the exact stored state.
- One failed campaign case produces a non-success aggregate status while preserving all case outcomes.
- Enabling checkpoint output does not modify the next numerical update; if exact continuation is claimed, compare interrupted/resumed and uninterrupted iteration histories.

## 5. D9: the recorded success is useful; the proposed budget formula is not

### 5.1 The arithmetic and the record do not support the stated rule

L36e reports the following binding-cell sequence for its atomic fiducial:

| Outer pass | Binding cell | Worst gated species row |
|---:|---:|---:|
| 1 | 217 | `3.51e-2` |
| 2 | 416 | `2.05e-2` |
| 10 | 444 | `1.91e-2` |
| 20 | 475 | `1.80e-2` |
| 25 | 491 | `1.76e-2` |

The executed arithmetic checks give:

```text
(491 - 217)/(25 - 1) = 11.4166666667 cells each pass
(491 - 416)/(25 - 2) =  3.2608695652 cells each pass
500/3.3              = 151.5151515152 passes
```

The approximately 3.3 value describes the later movement of the binding cell, not a traversal of all 500 cells from the start. Revision 1's “at least 40 passes before any gated row can begin to fall” is contradicted by the table, where the maximum falls immediately, and by its own statement that the outer boundary is reached at pass 26.

The reported acceptance at pass 56 is evidence for that configuration, initial state, algorithm, and tolerances. It does not derive a general lower pass bound for every transported-ionization run.

### 5.2 A moving maximum is not itself a measured front speed

The binding cell is an `argmax` of a normalized residual. It can jump between different cells or species as scales and residuals change. Pass 1 to pass 2 already shows a large jump. A physical or numerical front should be identified from its profile, a consistent stage-fraction threshold, or another separately defined feature, not only from the location of the largest certification measure.

Moreover, outer iterations are not physical timesteps, and the production grid is nonuniform. A fixed number of cells in one part of that grid is not the same radial or chemical distance elsewhere. The observed progression can depend on seed, grid, transport coefficients, coupling, trust controls, and reaction stiffness.

The front interpretation is plausible but should remain an interpretation until those profiles are checked. The source-level progress logic and the record do not establish a universal finite-propagation speed for the solver.

### 5.3 Revised budget policy

For the measured fiducial, a bounded repeat with a budget such as the recorded 90-pass allowance is a defensible test of reproducibility. Name it as an empirical allowance, not a law inferred from `N/3.3`. For other configurations, use a stated pass/wall-time ceiling with progress evidence and explicitly authorized extensions.

Record the residual profile, binding species and radius, relevant scales, physical admissibility, and solver-control history. Budget exhaustion while required rows remain above tolerance is **uncertified and incomplete**. It need not mean that no solution exists or that the solver can never converge, but must not be relabeled as a successful solution because a front is suspected.

A modest matrix of different seeds, resolutions, and chemical/transport regimes would determine which configurations need larger budgets. No full catalog campaign is needed to establish that policy initially.

### 5.4 Seed ranking is reasonable only after compatibility checks

The proposed preference for a nearby certified transported state is sensible. The actual `Physics` class in `pick_seed.py` does not yet parse or compare the ionization-transport key. The retained check constructs two inputs differing only in that key; `Physics.matches` returns `True`.

Adding a preference in one tier is insufficient:

- The same-case tier obtains an existing certified state before the other compatibility ranking. It must also inspect the state generation's resolved options, not just the current directory's input file.
- Compatibility must precede distance ranking: reservoir inventory, available species, molecular chemistry, diffusion model, boundary/radiation configuration, and the permitted conversion must be known.
- A local-equilibrium state may still be a legitimate initialization when transported states are unavailable. Record it as an explicit model-option transition with the appropriate restart permission, not as identical physics.
- No certificate from the source problem transfers to the target problem merely because that source is a good seed.

This work depends on D8's generation identity and the D5 seed distinction below.

### 5.5 Reference refresh is acceptance work, not harmless bookkeeping

The four key-on differences reported by L36c/L36d/L36e should be preserved as observations. They are not, alone, authorization to replace all reference values. D6/D7 and boundary/closure work can still change the correct result, and the newly confirmed TR+metals sign issue can affect a key-off route.

Freeze a coherent corrected implementation, pass the independent physical/source tests, and then refresh affected reference outputs once with provenance and measured differences. Keep unchanged paths out of unnecessary regeneration. A larger iteration budget and a new reference baseline do not fix missing reaction terms, wrong signs, or stale state evaluation.

## 6. D5's added molecular-seed restriction is too strong

Revision 1 correctly preserves the refusal of an incompatible direct restart. It then concludes that a composition step inside a molecular family is impossible whenever HeH+ is present and that the only entry to a molecular well-mixed case is conversion of an atomic well-mixed wind. Neither conclusion follows.

### 6.1 The current implementation refuses a particular rescaling, not every possible molecular seed

`load_IC.f90:650–667` refuses a He/H mismatch requiring rescaling when HeH+ is present and helium diffusion is disabled. A compatible molecular well-mixed state already having the target composition does not enter that mismatch branch.

The mapper also already contains a conservation check near `map_state_to_grid.py:709–728`. It tests the resulting reservoir ratio and rejects the single-factor transformation when molecules prevent it from satisfying the requested inventory. D5 should distinguish this implemented check from the additional compatibility checks it still needs, including the difference between a diffused column and a well-mixed target.

### 6.2 Conservation does not make a new initialization impossible

For example, hold total hydrogen nuclei at one and HeH+ at `0.1`. A target He/H ratio of `0.5` leaves `0.4` helium nuclei available to the atomic helium pool. This is a nonnegative, inventory-consistent state with nonzero HeH+. The arithmetic check verifies the resulting total He/H of `0.5`.

More generally, constructing a new seed can solve a constrained inventory problem using the species stoichiometric matrix. It must choose what to preserve, such as pressure, temperature, or conserved thermal energy, and must recompute density, electrons, and the EOS consistently. If a target budget cannot support the retained molecules, changing those molecules requires an explicit repartition rule. This is not a passive restart or a one-factor rescale.

The new state is uncertified initialization and must be solved under the target equations. A direct restart refusal remains appropriate until such a separate conversion is designed and authorized. The plan should say “unsupported by the present simple conversion,” not “physically impossible.” A compatible molecular seed and a suitable cold initialization are also distinct possible routes; the current runner's preferred route is not a physical uniqueness result.

### 6.3 Clarify the D2/D5 dependency

D5a can diagnose current evaluation order without a new boundary model. D5b can also make a deliberately frozen current model internally consistent. But revision 1 explicitly requires D5b's inflow/outflow composition rule to follow the D2 model, which is not yet chosen.

Either put the necessary D2a specification before that portion of D5b, or state that D5b first preserves the current boundary model and that any later physical change is a separate D2 increment. Do not leave both implications in the dependency diagram.

The D5 physical particle-count comparison should also remain a sensitivity observation rather than a rigorous bound on every effect of composition refresh. A molecular caloric state and radiation response depend on species partition, not only total particle count. The newly proposed detailed D5a state measurements are the right way to isolate the cause.

## 7. Remaining refinements to D1–D4 and D6

### 7.1 D1: numeric format and actual migration order

The direction of resolved-input disclosure is good. However, the suggested `ES24.17E3` is not a safe general signed-real format. The executed compiler check prints a positive `2e-4`, but prints 24 asterisks for `-2e-4`. `ES26.17E3` prints the negative value successfully. Seventeen digits after the decimal in scientific notation also mean 18 significant digits, not 17.

The listed grid quantities are positive, so this is not evidence that this proposed format would fail for the current grid width. It is a concrete reason not to use it as a generic resolved-real writer. Use a sufficiently wide explicit format, or an explicitly defined round-trip serializer, and test both signs, zero, representative exponents, and parsing back to the original binary value.

A global default cannot distinguish old and new input files by their age. “New configurations only” becomes true only after historical configurations are pinned or otherwise identified before the default changes. Unrecorded omitted-key inputs remain affected. State this migration precondition explicitly.

D1a should not regenerate numerical products. D1b tests should cover affected Mixed-grid defaults, explicit old widths, explicit new widths, and restart compatibility; do not require a complete unrelated output regeneration merely because a reporting field changed.

### 7.2 D2: one factual correction remains, and the decision needs concrete alternatives

The revised physical reasoning is much better. The opening still says “face Mach zero” for the zero-momentum fixture. L35 section 4.3 lists a nonzero numerical-flux-based local face Mach, `6.048816e-6` for the fiducial rest perturbation. An interior velocity, the characteristic boundary-state velocity, and a velocity inferred from a numerical Riemann flux are different quantities. Name the intended one instead of calling them all the face Mach.

Also, the table's `3.79361%` density difference uses the remote density as denominator. It is not literally the same percentage as the remote density's deficit relative to the reservoir density. This is a reporting correction, not a reason to select either closure.

For D2a, present a small set of complete physical choices and their incoming-characteristic requirements. Asking the user to specify an abstract contact model without those alternatives shifts too much of the design burden to the decision. Retain the manufactured-equilibrium, acoustic, weak-flow, and contact tests already added.

### 7.3 D3: include the new atomic branch and correct source diagnostics

D3a mentions exposing gross channels through `mol_heh_rows`. Since L36e, atomic stage transport also calls `heh_tr_rows`; that branch needs the same diagnostic coverage if the investigation is meant to inform stage certification generally.

There is a newly relevant discrepancy in the molecular proton diagnostic:

- `mol_heh_rows` returns `pHp` and `lHp`.
- `carrier_source` then changes `fv(1)` through `he_h_cx_fvec`.
- `src(Hp)` receives the corrected `fv(1)`, but molecular `sprod(Hp)` and `sloss(Hp)` still receive the original `pHp` and `lHp` near lines 5005 and 5031.

Consequently, `sprod(Hp)-sloss(Hp)` need not equal the actual assembled proton source. The retained production-kernel check supplies a pure H/He correction of `-1` while the earlier base contributions can remain zero. This is a source-ledger defect, not evidence that the final proton numerator omits that correction.

Require channel sums to reproduce the final assembled source after all corrections, for every gas branch and relevant switch. The complete stoichiometric ledger proposed for D7 can serve D3 as well.

The statement about normalized convergence order should remain conditional: a face-magnitude scale behaves like `1/h` where a smooth nonzero flux dominates it. The review did not remeasure the dimensional order of the actual L36d case. Do not turn the possible extra power of `h` into a newly measured dimensional order without executing the dimensional test.

Finally, D3a step 4 enables a formerly report-only stage-sum gate. That is an acceptance-interface change, not just measurement. Split the measurement/bound derivation from enabling the gate, add `certification.f90` to ownership, and explicitly approve that interface change. The algebraic identity still deserves a gate; its implementation should be honestly classified.

### 7.4 D4: preserve diagnostic isolation and conditioning

The revised D4 is substantially correct. Retain both possible outcomes: a production stopping-criterion defect or an unjustified abundance-movement bound in the test. Neither has been eliminated merely by identifying a small failing number.

An extra sweep used only to measure a closure increment changes module state as well as the array passed to it. Save and restore the rates, radiation-dependent context, caloric state, boundary caches, and ledgers needed by the next real operation, or evaluate on a properly isolated snapshot. Otherwise the diagnostic can itself settle the state production will later use and conceal the defect.

For coupled ion/electron balances, use the relevant constrained Jacobian or a controlled perturbation study to connect residuals to abundance error. The scalar estimate in D3 is useful for intuition but is not a universal error bound for that coupled network.

### 7.5 D6: add the excited-level constraint to source admissibility

The HeH+ inventory correction remains valid and is not yet present in the inspected projection. Extend its tests to the fact that He 2^3S is a subset of neutral helium:

\[
0\leq n_{\mathrm{He}(2^3S)}\leq n_{\mathrm{HeI}}.
\]

`carrier_source` takes `n_heiTR` from the frozen background, recomputes total neutral helium from the trial ion fractions, and uses `max(n_hei-n_heiTR,0)` for the singlet. If the trial leaves less neutral helium than the frozen triplet population, this clip hides an inconsistent trial state even after reserving HeH+.

If the triplet population is frozen during a trial, the admissible set must also reserve that neutral subset. If it is recomputed or repartitioned, define the corresponding derivative and closure consistently. Do not count the triplet as an additional helium nucleus in the census. The later write-back rescales it with neutral helium, but that does not repair a source already evaluated at an impossible trial partition.

This is a source-level admissibility observation; no catalog frequency or spectral impact was measured here. It belongs in D6's targeted tests because the plan already promises consistent admissibility during source evaluation and differentiation, not just after write-back.

## 8. Baseline, ownership, and implementation order

### 8.1 The proposed baseline procedure can interfere with other work

Section 0.2 requires `make -q` to return zero and then a build from clean objects in the current tree. A dirty, actively edited tree may correctly be out of date. Requiring an up-to-date result before building is not a baseline definition.

More importantly, this Makefile executes build-stamp shell commands while being parsed, including creation and possible replacement of files under `OBJDIR`. `make -q` is therefore not a strictly read-only probe in this project. A shared clean build can also collide with the work the plan says is still running.

Use a captured source snapshot or an otherwise frozen build interval with a content manifest, and a separate object directory and executable. This Makefile accepts command-line `OBJDIR` and `EXE` overrides. Record compiler, effective flags, linked libraries, source hashes, and the resolved input/data identity. Verify the source manifest remained the same throughout the build. An executable hash names the artifact; it does not establish which unrecorded input data it consumed.

The current review used a narrow isolated build and did not run `make clean`, replace the main binary, or demand that the active tree be up to date.

### 8.2 Update the ownership table and dependency graph

The table omits D7, D8, and D9 despite scheduling D8/D9 first. It also omits `certification.f90` from the D3a gate change. At minimum, add:

- D7: `charge_exchange.f90`, relevant local equation systems, `carrier_source` and frozen-background construction, input validation, source/Jacobian tests, and parallel consistency tests.
- D8: runner, campaign, writer/checkpoint interface, readers including the mapper and loader, seed selector, status/reproduction reporting, and the outer-loop snapshot point.
- D9: campaign budget configuration, seed ranking across all tiers, outcome reporting, and affected reference-generation scripts.
- D3: both molecular and atomic reaction kernels and the certification interface when a gate is enabled.

Calling D8 and D9 mere bookkeeping understates their scope. A checkpoint format and interruption protocol change data interfaces; a new retry policy changes the algorithm followed; a reference refresh changes what regression accepts. None changes a physical equation by itself, but all require explicit contracts and tests.

### 8.3 Recommended staged order

1. Preserve current run generations and freeze targeted experimental snapshots without altering catalog products.
2. Add the D7 unsupported-combination refusal; reproduce and correct the existing TR+metals sign inconsistency with isolated reaction tests. Plan the general stage-source addition after its interface is defined.
3. Implement D6's complete inventory/level admissibility and its independent counterexamples.
4. Execute D5a and D4a measurements, with diagnostics that do not change later behavior.
5. Execute D3a source/conditioning measurements on valid supported configurations. Metal-bearing measurements depend on the relevant D7 correction; atomic and molecular branches are labeled separately.
6. Agree on D2a and the D5b state interface/dependency. Implement the selected consistency and physical changes in separate measured increments.
7. Implement the agreed D8 generation/checkpoint protocol and reader updates. A minimal preservation-only runner improvement can occur earlier, but it should not claim exact continuation or attach stale certificates.
8. Measure D9 budgets on a small declared matrix and add compatibility-first seed ranking. Refresh only affected references after the relevant physics and interface tests pass.
9. Perform D1 disclosure independently when its files are available; change the default only after historical input pinning is established.

## 9. Acceptance matrix for the revised plan

| Concern | Minimum evidence before declaring the item complete |
|---|---|
| D7 source completeness/orientation | Isolated reaction signs and conservation, correct residual-basis conversion, source-ledger identity, dimensions checked, thread-order independence. |
| D6 admissibility | H/He/mass census and level bounds on actual trial/returned species, including HeH+ and nearly exhausted neutral helium. |
| D5 same-state evaluation | Same physical state and inputs produce consistent boundary, caloric, cached-face, flux, and signed residual quantities across all relevant call paths. |
| D4 closure | Direct chemical/thermal measurements on the unmodified returned state, with diagnostic isolation and justified residual-to-error interpretation. |
| D3 norms and stage-sum gate | Dimensional residuals and conditioning measured; algebraic error bound tested; any gate activation explicitly identified as an interface change. |
| D8 snapshot correctness | Certificate and snapshot have one state identity; crash injection cannot publish mixed/incomplete generations; consumers reject stale certificates. |
| D8 continuation | Supported ending reasons only; warm restart versus exact continuation stated; best/latest/certified products kept distinct. |
| D9 budget | A recorded case-specific convergence history, honest incomplete outcomes at ceilings, and no certification inferred from movement of an argmax alone. |
| D1 migration | Resolved values round-trip; old configurations are pinned before the default changes; only intentional new-grid seeds are mapped. |

## 10. Limitations and reproducibility

No full regression matrix, new catalog solve, long transported-front experiment, or full TR+metals atmosphere was executed. The L34c/L36e pass histories and reference differences remain recorded results. The narrow checks establish the local source-orientation issue and selected workflow behavior without attributing a numerical size to their effect on a complete production solution.

Literature verification used the local final journal version of Huang et al., including its methods subsection and Table 4. ADS was attempted first but the shell request failed on DNS resolution; opening the DOI through the web tool was blocked. The citation and content comparison therefore rely on the local published PDF, not a successful online metadata response or a preprint. No claim is made that every rate in the table was independently revalidated in this review.

The repository HEAD was `3c73905ca8a7fe2af92a5c2b014c225b2feedff3`, with substantial existing changes. Selected SHA-256 identifiers from this review interval are:

```text
docs/PLAN_20260918_rev1.md
063d2facb6370de88de69f8d209e82645e24d3b32238c0ea90b9e5872574ad62

src/EXHALE_main.f90
effc018fed7950c494c0fdf27cd78ab1837800f770e6850a488aa86c0d582f6a

src/modules/radiation/charge_exchange.f90
c24e394f26dad00df18adb12c92145fe63bb9402895541579e8a14ef4db3cbbf

src/modules/nonlinear_system_solver/System_HeH_TR_metals.f90
d70185484e1da140a371be7dc42705ad1508d766b8cf4f0c61970efd960c02c5

src/modules/lower_atmosphere/diffusive_photochemistry.f90
a9d7c0a1eaae367b8c4550ce5ecbd71400daa8448ba3e46ec621dc35d22c5884

src/utils/map_state_to_grid.py
672fcb62f0c7b17d4e17d249e55dfa33e66af91e88cee5b39272a79758f6261b

LHS1140b/models/run_case.sh
8d9819fe214bfba7f7be8e33cb93cde049ff5ad429de2ca440773c6b2ab623e2

LHS1140b/models/run_campaign.sh
8d9c170c3561ea9f029a2eaf191a569d1fb508ab0ed904f4d3b510d257468cbe

LHS1140b/models/pick_seed.py
adfd70e30c8707f0aab0c6b9485b11994846abdbded98a29c368ce21fc3b7e8e
```

Run the retained targeted checks from the project root with:

```bash
bash docs/audit_20260905/run_plan_20260918_rev1_review_checks.sh
```

These checks intentionally report current defects. Their successful execution means the observations were reproduced, not that the affected production behavior is correct. The changed artifacts of this review are this report and the retained audit programs/log; no production correction has been applied.

## Conclusion

Revision 1 is a much stronger plan than its predecessor. D1–D6 largely move in the right direction, with the qualifications above. The principal remaining blockers are now concrete: D7 must not copy an existing wrong helium-row sign, D8 needs transactional generations and certificates tied to the exact saved state, and D9 needs measured budgets rather than the unsupported 40-pass rule. Correct those points and the molecular-seed overstatement before treating the entire plan as ready for implementation.
