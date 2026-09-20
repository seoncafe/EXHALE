# Review of PLAN_20260918, including D5

Date: September 18, 2026 (KST)

## 1. Overall judgment

The plan identifies worthwhile problems, especially the molecular restart inconsistency in D5. However, it should not be executed unchanged. D2 starts from an incorrect characteristic argument, D4 would replace a returned-state acceptance test with a test that repairs the state, and D3 has not established that the refusing row is impossible to satisfy. D5 has the right architectural direction but understates the coupled boundary problem and overstates what its proposed implementation would guarantee.

There is also a separate physical admissibility defect in the inspected source: the helium ion-stage projection does not reserve the helium already held in HeH+. This must be addressed before a different residual normalization is used to accept transported helium states.

| Item | Recommendation | Required correction |
|---|---|---|
| D1 | Approve precision disclosure; revise migration strategy. | Preserve old grids with explicit resolved inputs. A default change does not require mapping every archive. |
| D2 | Do not approve the proposed correction yet. | Specify the physical entropy/composition boundary at zero flow and reversal. An acoustic characteristic still leaves the domain at rest. |
| D3 | Approve a diagnostic investigation, not a new acceptance rule. | Measure the actual transport/reaction balance and its conditioning; do not choose a denominator simply because it makes the fixture pass. |
| D4 | Reject the test-only replacement. | Test the state production returns. Diagnose and, if necessary, repair the production thermochemical closure or explicitly revise an unjustified test contract. |
| D5 | Approve a staged boundary/restart consistency investigation and redesign. | Treat hydrodynamic ghosts, composition, caloric state, cached face states, and initialization order together. Do not promise bit identity or certification of every historical state in advance. |
| Additional finding | Prioritize a targeted conservation fix. | Include HeH+ in the helium admissible set and verify actual species inventories after write-back. |

The suggested implementation order is therefore not “D4 test change, D5, D2, D3.” First establish coherent state evaluation and valid element budgets, then decide which failures are solver errors and which require changes to physical boundary conditions or acceptance norms.

## 2. Scope and evidence

This review covers the updated 533-line plan containing D5. It checks the implementation and call sites, not only the descriptions in the plan. The working tree already contains substantial changes from ongoing work. No production source, catalog product, regression reference, input, or tolerance was changed for this review.

The inspected paths include:

- Grid defaults, input parsing, construction, restart radius checks, setup reports, and resolved configuration output.
- `characteristic_base_face_state`, `Apply_BC`, `Rec_BC`, reconstruction, stationary residual assembly, and the stationary restart evaluation route.
- Carrier state construction, physical residual normalization, helium source terms, stage projection, write-back, and certification.
- The fixed-conserved-state thermochemical closure and its returned-state tests.
- The ion-stage manufactured-column test and the L35, L36c, L36d, and L37 records relevant to the plan.

Evidence labels used below have distinct meanings:

- **Source finding:** directly established by the inspected expressions and their callers.
- **Arithmetic check:** independently executed on a small counterexample; not execution of the EXHALE implementation.
- **Recorded result:** a value read from an existing experiment record, not remeasured here.
- **Hypothesis:** a possible explanation requiring a controlled production experiment.

The retained arithmetic checks are [plan_20260918_review_checks.py](audit_20260905/plan_20260918_review_checks.py), with output in [plan_20260918_review_checks.log](audit_20260905/plan_20260918_review_checks.log). They demonstrate the precision issue, the discontinuity of a change confined to zero, the distinction between chemical turnover and transport accuracy, the helium inventory counterexample, and nonmonotone convergence of an iteration.

## 3. D1: real grid-configuration defect, unnecessarily expensive migration

### 3.1 What the implementation establishes

`src/modules/init/parameters.f90:143` and `src/modules/files_IO/input_read.f90:451` assign the default with `2.0e-4`, whereas the input key is read into a double-precision variable near `input_read.f90:1174`.

With the binary32/binary64 behavior assumed by the current build, the executed arithmetic check gives:

```text
default-real literal widened to binary64: 0.00019999999494757503
decimal input parsed as binary64:        0.00020000000000000001
relative difference in width:            2.5262124836444855e-08
```

The old value is recovered exactly by the explicit string `0.00019999999494757503`. The maximum relative radius displacement of `6.58e-9` is a **recorded L35 result**, not a newly executed grid measurement. It is distinct from the relative difference in the width itself.

The restart radius test at `load_IC.f90:2221–2244` uses a relative threshold of `1e-10`. Keeping that compatibility check is appropriate. Increasing it to accommodate a silently different grid would conflate loading a state with interpolating a state.

### 3.2 Correct the claimed scope

`define_grid.f90:31–35` assigns `drc = dr_base`, but uses this base-spacing construction in the `Mixed` branch. The claim that changing the literal moves “every grid” is too broad. The affected set is Mixed grids that depend on this default; explicitly specified widths and other grid constructions must be distinguished.

More importantly, an old archive does not have to be mapped merely because a new default is introduced. It can retain its original grid by specifying its resolved old width and base-cell count explicitly. This preserves the archived numerical problem rather than changing it and then trying to recover its solution.

Recommended migration:

1. Record the resolved grid type, width, base-cell count, total cell count, outer radius, and relevant geometry version.
2. Add explicit historical values to reproduction inputs for archives whose configuration is established.
3. Use the corrected double-precision default for newly defined calculations after the user approves that default change.
4. Map a state only when intentionally moving it to a different grid. Such a mapped state is a seed, not a stationary solution certified on the new operator.
5. For archives without enough provenance, report that the resolved configuration is unknown; do not silently assign the present default.

A full catalog re-solve may be useful for other reasons, but is not a mathematical prerequisite for correcting the default. Likewise, L30 sensitivity to a `4e-13` seed perturbation does not predict the response to a grid/operator change. The plan should remove that quantitative forecast.

### 3.3 Improve the proposed disclosure

There is already a parsed-configuration dump: `write_setup_report.f90:732` writes `dr_base`. Its `put_r` format near line 858 is `ES23.15E3`, which supplies 16 significant digits, not a universal 17-digit round-trip representation for binary64.

There is also `write_resolved_config` near line 883, producing a machine-readable resolved configuration. Extend that record with the grid controls and sufficient precision, instead of making a reproduction utility depend on human-oriented report text alone.

Reporting “default” versus “key” requires recording whether the key was supplied. Value equality cannot establish its origin. The D1 ownership list must therefore include the input reader and the declaration of that provenance flag. A single named default should replace the duplicated initialization literals.

Tests should compare parsed values and generated coordinates directly. The L35 radius discrepancy should remain a case-specific observation, not a universal hardcoded expected discrepancy for all grids.

### 3.4 Separate this from physical grid convergence

The current Mixed geometry also places the base face according to the ghost/interior center construction. In the uniform base region, `r(0)=1`, `r(1)=1+dr_base`, and the intervening face is at `1+dr_base/2`. Refining the width therefore changes the location of the numerical boundary relative to the reservoir reference level.

This is a distinct issue from literal precision. A resolution study must state whether it holds the physical boundary face fixed. Changing that geometry would require its own design and validation; neither printing 17 digits nor replacing the literal resolves it.

## 4. D2: the stated rest-boundary argument is physically incorrect

### 4.1 A stationary fluid still has acoustic characteristics

The normal Euler characteristic velocities are

\[
v-c_s,\qquad v,\qquad v+c_s.
\]

At rest these are `-c_s`, `0`, and `+c_s`. At the lower boundary of a domain extending toward increasing radius, one acoustic characteristic leaves the domain and one enters it. The contact/entropy characteristic is stationary. The plan's assertion that no characteristic leaves at rest is therefore false.

This is consistent with the characteristic construction in the final published version of [Thompson, 1987, Journal of Computational Physics 68, 1–24](https://www.sciencedirect.com/science/article/pii/0021999187900416), also available locally as `references/Thompson_1987JCP_68_1.pdf` relative to the workspace root. The paper does not establish that reservoir density must replace every zero-velocity face state.

Whether the zero-speed contact receives reservoir entropy, interior entropy, or a regularized transition is an additional boundary-model decision. It must follow the physical reservoir/contact problem. This observation does not prove that the current arithmetic average is the best boundary model; it invalidates the proposed proof that the average is necessarily wrong.

### 4.2 Zeroing a wind velocity is not a hydrostatic equilibrium test

`characteristic_base_face_state` in `src/modules/states/base_boundary.f90:618–816` continues the interior state to the face and uses the outgoing acoustic relation

\[
v_b=v_i+\frac{p_b-p_i}{\rho_i c_i}.
\]

Consequently, setting the interior velocity to zero does not imply `v_b=0`. A pressure mismatch can generate a legitimate acoustic response.

The rest fixture in `src/tests/grid_and_gates/base_closure_candidates_probe.f90:300–320` removes momentum and the corresponding kinetic energy from a wind state. It does not construct a new density/pressure profile satisfying the discrete hydrostatic equations and the same reservoir boundary condition.

L35 records a nonzero base mass flux in that fixture. Its section 4.4 gives `1.2751005420e-5` for the production closure and `1.3790806400e-5` for the local-reservoir candidate. These are **recorded results**. The candidate increases that flux by about 8.2%; it does not demonstrate restoration of a stationary equilibrium.

Well-balanced preservation must be tested on the discrete equilibrium of the actual discretization. The published [Käppeli and Mishra, 2016, Astronomy & Astrophysics 587, A94](https://www.aanda.org/articles/aa/pdf/2016/03/aa27815-15.pdf), particularly the discrete balance and preservation argument around equations (18)–(20), makes this distinction explicit. Its particular discretization should not be mistaken for a direct proof about EXHALE's molecular boundary.

### 4.3 A smooth change cannot affect only the zero point

The current branch weight is continuous and equals `0.5` at zero. If the new weight equals the old weight at every nonzero point approaching zero, continuity requires its value at zero to remain `0.5`.

Thus a proposed weight that equals zero at rest, is continuously differentiable, and leaves every arbitrarily weak nonzero flow unchanged cannot exist. A smooth multiplier can vanish at zero, but then changes a neighborhood of zero. Its width becomes a model/numerical parameter requiring justification.

The retained arithmetic check shows the divergent difference quotient if only the zero value is replaced: for arguments `1e-2`, `1e-4`, and `1e-6`, the quotient grows from about `49` to `5e3` to `5e5`.

Byte identity may hold for a selected set of saturated flowing controls. It cannot be promised for all flowing states. A singular discriminant whose argument is sent to infinity at zero also needs a complete limiting and finite-arithmetic analysis; absence of an `if` statement is not a differentiability criterion.

### 4.4 A better D2 decision and experiment

Specify first:

- Which thermodynamic quantities the lower atmosphere prescribes, and at which reference radius.
- Whether entropy and molecular composition are imposed only for inflowing material or also through a separate thermal/contact boundary model.
- How a stationary contact with different interior and reservoir entropies is treated.
- How the boundary transitions through weak inflow, weak outflow, and reversal without imposing outgoing information.

Then test:

1. A manufactured discrete hydrostatic state consistent with the molecular EOS and reservoir, checking face velocity, mass flux, momentum balance, and energy balance.
2. A stationary pressure-balanced contact with differing entropy/composition, with its expected behavior stated in advance.
3. Acoustic perturbations of a stationary background, retaining the outgoing acoustic information.
4. Weak flows of both signs and directional derivative tests around zero local and window flux, including a zero window mean with nonzero spatial fluctuations.
5. Reservoir-limited and interior-limited branches, comparing composition, pressure, caloric energy, and characteristic count together.

Keep D2 separate from D5 initially. D5 should make a declared boundary model reproducible; it should not silently select a different entropy boundary to obtain a passing restart test.

## 5. D5: high-priority consistency issue, but not just file precision

### 5.1 What the recorded experiment supports

L37 contains strong evidence that the lower ghost state affects molecular restart residuals. In its swap experiment, a settled state supplied with the archived ghost composition reproduces the archived first-evaluation mass row. The first-evaluation discrepancies are materially larger than the existing rounding allowance.

It is reasonable to reject widening `c_round` to 100 as the first response. A changed derived boundary is not the same perturbation as unavoidable rounding in one evaluation of a fixed operator.

However, the record also reports disagreement between two certifications inside one run. That makes initialization order, cached state, and closure consistency relevant even without serialization. The evidence does not isolate insufficient printed precision as the sole cause.

### 5.2 Correct the factual summary before using it as acceptance evidence

These are comparisons with the tables in L37, not new EXHALE measurements:

| Claim in the D5 summary or supporting summary | What the detailed L37 record actually shows |
|---|---|
| All molecular states settle by the third or fourth re-entry. | Section 4.2 says four of six settle; `heh213_pre` still changes and `heh055_new` has continuing changes in cells 7 and 8. |
| Atomic states are settled on the first reading. | `at213` has a period-two sequence with a stated amplitude of about 1.4%; only the other two atomic states repeat immediately throughout the tested cells. |
| Atomic response to the composition perturbation is zero at printed precision. | Section 4.4 changes `at213` from `1.1852e-9` to `9.5494e-10`, and `at97` from `1.8517e-10` to `1.8526e-10`. The response plateaus with the tested perturbation size; it is not zero. |
| Ghost density drift is of order `1e-8`. | Section 4.1 gives relative density changes `9.07e-16` and `1.36e-16`. The `4.835e-8` number refers to ghost H I, not density. |
| The conserved state is established to round-trip bit for bit. | The displayed comparison is principally of `rho`, `v`, `p`, and species. It does not directly compare stored/reconstructed total energy. |

The distinction in the last row matters for molecules. `write_output.f90:298–319` writes density, velocity, pressure, temperature, and source quantities, not the conserved total energy as an independent column. When caloric energy depends on composition, matching printed `rho`, `v`, and `p` does not prove identical `u(3)` after decoding and EOS conversion.

Also, L37's six molecular snapshots were not all produced by the evaluation binary. Its provenance table includes `3146d11b`, older `c2e9c999` states, and an L22 state. Same-operator round-trip consistency and validity of a historical certification under changed equations are different questions. Historical certification must not be treated as an unconditional expected pass under a new operator.

### 5.3 Source finding: the evaluation path mixes closure stages

The current call sequence in `stationary_state_of_the_loaded_restart` is important:

1. `EXHALE_main.f90:5483–5486` installs composition-derived quantities and calls `Apply_BC(u)`.
2. `Apply_BC_W` computes and stores `base_face_W`, `base_ghost_W`, and `base_face_lower_W` in `Apply_BC.f90:175`.
3. The restart route calls `ioniz_eq` near line 5502, updating composition and chemical state.
4. It refreshes species/EOS quantities and derives pressure at fixed conserved energy near lines 5517–5524.
5. It calls `assemble_residual` near line 5594 without another `Apply_BC` between that refresh and the assembly.
6. `assemble_residual` calls reconstruction and the hydrodynamic RHS. `reconstruction_continuation_rhs` does not rebuild the boundary. `Rec_BC` inserts the previously cached `base_face_W` at `Apply_BC.f90:343`.
7. Later, `stationary_face_mass_flux` does refresh composition-derived state and calls `Apply_BC` on a copy before reconstructing its reported fluxes (`stationary_operator.f90:155–173`).

Therefore, the residual and the later face-flux report can use boundary states installed at different closure stages. This is a directly identifiable consistency risk; its numerical contribution to each L37 discrepancy still needs measurement. It must be tested before concluding that changing only the loader resolves D5.

`Apply_BC(u)` currently has no composition argument. It updates hydrodynamic ghosts using the already installed caloric mixture. Its interior-preserving behavior is an improvement and should be retained, but it does not make the complete hydrodynamic/composition boundary idempotent.

The state-installation semantics also matter for diagnostics. `stationary_face_mass_flux` restores reconstruction selection, but explicitly leaves composition-derived module state installed. Do not assume that saving the caller's conserved array makes this diagnostic free of relevant side effects. Repeated and interleaved diagnostic calls must be included in the consistency tests.

### 5.4 Source finding: ghost chemistry currently retains an iteration history

In `ionization_equilibrium.f90:1999–2010`, the imposed lower-boundary molecular fraction is computed as

\[
x_{\mathrm{H_2,all}}=q_{\mathrm{H_2,base}}(1-x_{\mathrm{ion,entry}}),
\]

where the ionized hydrogen-nucleus fraction uses the composition entering the sweep, including H II, H2+, H3+, and HeH+. The code explicitly describes reaching the joint condition by alternation.

This explains why deleting a ghost-file read is not sufficient: one chemical sweep from two different admissible starting ghost compositions can impose two different H2 constraints. A complete boundary reconstruction must either solve this coupled condition to its declared accuracy or define and carry any genuinely independent reservoir state explicitly.

The ionization sweep also depends on the radiation environment and thermochemical quantities. The boundary is not generally determined by the first physical cell alone. D5 should describe it as a function of the full declared physical state, prescribed reservoir data, and radiation/boundary inputs, not of unspecified cached state.

There is a further consistency question to resolve explicitly. `input_read.f90:2116–2117` initializes the stored reservoir pressure, temperature, and particle count once. The sweep later updates `ntot_bc` and `dp_bc` near `ionization_equilibrium.f90:3168–3195`, but this is not another call to `set_base_reservoir`. D5 must identify which count belongs to the prescribed reservoir and which belongs to a solved ghost. Updating the molecular EOS while retaining a different reference particle count is not automatically a consistent entropy/density prescription. Conversely, refreshing a prescribed reservoir pressure from a trial state would change the boundary problem. The intended physical definition must decide the synchronization rule.

### 5.5 Recommended design

Prefer option (1) in principle: ghost values that are derived boundary quantities should not be independent restart degrees of freedom. But replace “make `Apply_BC` form composition on every call” with a defined boundary-state operation having an explicit input/output contract.

The operation should conceptually take the physical conserved state and composition, reservoir specification, radiation context, and model options; it should return a mutually consistent ghost composition, ghost conserved state, caloric state, and reconstructed boundary face data.

Requirements:

- Physical cells remain unchanged during a boundary reconstruction.
- Hydrogen and helium inventories, charge, species nonnegativity, and excited-level bookkeeping are respected. He 2^3S is a level within He I, not an extra helium nucleus.
- The molecular partition and ionization prescription are satisfied simultaneously to a justified boundary-closure tolerance.
- The same composition defines pressure, temperature, heat capacity, sound speed, and conversion between primitive and conserved ghost variables.
- The inflow/outflow composition rule follows the declared physical boundary model; a file-read repair must not impose reservoir material on an outgoing characteristic by accident.
- Every residual consumer obtains boundary data corresponding to the same installed state. A cached face state must have explicit validity or be recomputed with its inputs.
- If a reservoir quantity is physically independent and cannot be reconstructed, serialize that quantity with its physical meaning and versioned provenance. Do not hide it in an arbitrary ghost row.

A coupled boundary solve may contain internal iterations. This is fundamentally different from repeatedly writing and rereading the entire state until certification passes: the former solves a specified boundary equation with a reported residual and failure condition; the latter changes the evaluation contract without identifying the equation being closed.

Changing the boundary interface can be the correct redesign. Its intended state layout and ownership should be agreed before implementation because it reaches initialization, marching stages, steady residuals, checkpoints, restart evaluation, and output. The present D5 file list is too narrow.

### 5.6 Stronger tests than an “atomic band”

The ratio `E1/in-run` is unstable when the reference residual is close to zero, discards sign, and can mask changes in the normalization. It is not a reliable equality criterion for two evaluations of one state.

Use the following targeted tests instead:

1. **Boundary independence:** hold all physical cells and declared inputs fixed; replace old ghost data with several different finite, admissible states. The reconstructed boundary must agree within a justified arithmetic/closure bound.
2. **Boundary idempotence:** reconstruct twice, then compare composition, conserved ghosts, caloric quantities, cached face states, and boundary-closure residuals. Also evaluate state A, then B, then A to expose stale global state.
3. **Same-state operator agreement:** compare the signed mass residual, both adjacent face mass fluxes, normalization terms, and rounding floor between in-run certification and stationary evaluation.
4. **Serialization:** compare physical conserved variables, including energy, and species before writing and after loading. Separately measure boundary reconstruction and residual assembly. This separates representation error from closure error.
5. **Actual call paths:** exercise a marching-stage boundary, steady trial, restored checkpoint, stationary restart, and the face-flux diagnostic. The same named state must lead to the same boundary/operator.
6. **Model limits:** include atomic and molecular states, molecular handoff on/off, transported ionization on/off, weak inflow/outflow, and a pressure-balanced rest fixture. Radiation inputs must be held fixed or recomputed by the same declared rule.
7. **Repeated evaluations:** test once, twice, and in a fresh process. No physical cells may drift merely because a diagnostic was requested.

For residual comparison, use an absolute difference in dimensional residuals or normalize the signed difference by a common recorded scale and a justified arithmetic floor. Report any change in that scale separately. Bit identity is a useful diagnostic where the same arithmetic is intended, but must not be promised before accounting for EOS round trips and nonlinear closure.

### 5.7 Catalog handling

Do not begin with in-place catalog evaluation or `status.py --write`. Work on immutable copies with recorded binary/configuration identifiers. Distinguish:

- reproduction under the same operator;
- measurement under a corrected operator;
- a new stationary solve required because the corrected boundary changes the root.

A corrected operator is allowed to reject a formerly certified state. Requiring both catalog states to certify immediately could turn a physical correction into pressure to weaken the test. Update catalog products only after the result is understood and the relevant update is authorized. Keeping a backup is useful, but does not make premature replacement an appropriate validation method.

## 6. D3: establish the balance and its conditioning before changing the gate

### 6.1 The actual denominator is already identifiable

The relevant call chain is:

```text
certification_evaluate
  -> carrier_rows_of_state
     -> carrier_steady_residual
        -> carrier_residual
  -> carrier_row_entry
```

`carrier_rows_of_state` starts near `certification.f90:1438`. The physical carrier ratio is formed at `diffusive_photochemistry.f90:3546` with `row_terms_phys`. The inventory/signal scale computed nearby is used for the optional legacy diagnostic, not this certification denominator.

In `carrier_residual`, near lines 5282–5354, the physical scale contains:

\[
S_j=\frac{A_R|J_R|+A_L|J_L|}{V_j}
    +S_{\mathrm{adv},j}+|\mathcal S_j|+S_{\mathrm{absent},j}.
\]

Here `J` denotes the corresponding diffusive contribution, `S_adv` the absolute advective contributions, and `mathcal S` the net chemical source. The time term is excluded from `row_terms_phys`. The current geometry uses the shell volume from the shared geometry routine; the previously reviewed alternative volume formula should not be reported as a current defect.

Thus neither of the plan's simplified alternatives is a complete description: the denominator is not just the net chemical source, and it is not the whole-element signal-speed scale.

### 6.2 “No composition can certify” is not demonstrated

Let `D` be the total transport divergence, and let the chemical source be `P-L`. The steady equation is

\[
R=D-(P-L)=0.
\]

At a discrete steady root, this residual is zero regardless of which positive scale is used. A denominator containing `|P-L|` can be badly conditioned near local equilibrium, particularly if transport is zero and the denominator tracks the same net source as the numerator. That is a legitimate concern, not a proof that a transported steady root cannot pass.

A locally equilibrated control satisfies `P=L`, not necessarily `D=P-L`. It is not an independent exact reference for a transported problem when `D` is nonzero. Close agreement in abundance, including the plan's recorded small composition discrepancy, can coexist with a large residual in a stiff reaction.

Before restricting the investigation to normalization, measure the worst cells with the actual transported fractions held consistently and determine:

- `P`, `L`, `P-L`, advection, diffusion, the unscaled residual, and every denominator component;
- the reaction sensitivity to stage fractions and electron density;
- the abundance correction required to balance transport, approximately `D/(d(P-L)/dx)` near local equilibrium;
- the floating-point cancellation floor, checked with a higher-precision source evaluation where practical;
- whether the residual changes smoothly under perturbations that respect the element budgets;
- whether the same composition and radiation field are used by the equilibrium sweep and transported source evaluation.

These checks may identify a conditioning issue, a genuine solver defect, an inconsistent source evaluation, or a combination. The present evidence does not justify ruling out solver work in advance.

### 6.3 Existing helium diagnostics do not expose gross production and loss

`carrier_source` near lines 4985–5006 explicitly fills the helium diagnostics as

```text
sprod = max(src, 0)
sloss = max(-src, 0)
```

These are the positive and negative parts of the **net** source. They are not the sums of individual production and destruction channels. Their sum cannot measure fast cancellation between photoionization and recombination.

D3a therefore needs access to a true reaction-channel decomposition. Prefer optional diagnostic outputs from the same chemical-network expressions used by `mol_heh_rows`, including subsequent source corrections. Do not reproduce the chemistry independently in a test and risk validating a different equation.

The claim that D3a touches only test programs and a memo is consequently incomplete unless an existing lower-level API already exposes all the required channel rates. That API must be identified and its values checked against the assembled source.

### 6.4 Gross turnover is useful but is not sufficient for transport accuracy

The executed arithmetic counterexample uses

```text
P = 100000000, L = 99999999, correct D = 1.
```

The correct balance gives zero residual. For an incorrect divergence `D=2`, the residual is one. With absolute face contributions of three:

```text
|R| / (3 + |P-L|) = 0.25
|R| / (3 + P + L) = 4.99999995e-9
```

The turnover-normalized row passes `1e-5` despite an order-one error relative to the transport balance in this example. A turnover norm measures chemical backward error; it does not by itself certify a correct transport balance or an accurate trace abundance.

Use complementary criteria where stiffness requires them: chemical backward error, transport/global inventory balance, and an estimate of the abundance or observable error implied by the remaining residual. When comparing proposed norms, compare the allowed dimensional residual `tolerance × scale`, not just the numerical size of the normalized residual or threshold.

If a fast stage is intentionally eliminated by a local-equilibrium approximation, state its timescale validity and validate the resulting reduced model. Omitting a required balance from certification while retaining an unrestricted “certified” label is not an acceptable substitute.

### 6.5 The present manufactured test does not establish a universal stiff-stage tolerance

The manufactured-column test in `src/tests/ionization_stage_flux/ionization_stage_flux_tests.f90:1711–1906` is useful for the transport discretization. It prescribes smooth profiles, computes an analytic advective/diffusive divergence, and compares the numerical divergence with that manufactured source.

It does not exercise a stiff photoionization/recombination cancellation with independently large production and loss. Its principal convergence measure also excludes the two boundary cells. It therefore cannot alone anchor the full helium reaction/transport gate or the D5 boundary behavior.

There is an additional interpretation issue. Near lines 1823–1829, its scale contains absolute flux magnitudes divided by cell volume. For a smooth nonzero flux, that scale grows as `1/h`. If the dimensional truncation residual scales as `h^p`, the reported normalized residual can scale as `h^(p+1)`.

Accordingly, the recorded normalized order `2.886` is not automatically a third-order dimensional truncation result. It can be consistent with a dimensional order near `1.886` in the regime just described. The assertion requiring normalized order at least two also does not, by itself, prove second-order dimensional accuracy.

Report dimensional residual convergence, solution error where available, and convergence under a grid-independent reference scale alongside the current acceptance-normalized result. Extend the manufactured tests to stiff reaction/transport balance, several chemical-to-transport timescale ratios, the Mixed grid, and the actual boundary operator before extrapolating a tolerance to the production He II problem.

### 6.6 Structural negative tests are a separate obligation

Charging helium to the hydrogen flux and projecting helium stages independently are valuable deliberately incorrect controls. They test structural conservation and admissibility, not simply the magnitude of a species steady residual. A denominator should not be selected merely because those controls happen to fail under it.

The production stage-nucleus-sum entry is currently report-only: `certification.f90:1260–1280` sets `within_tol = .true.` with no gating tolerance. Therefore, a failing structural test does not imply that the corresponding production certification would refuse the state.

The sum of stage fluxes is an algebraic identity, not a spatial truncation approximation. Its gate should have a floating-point error bound consistent with the number and magnitudes of summed terms, with execution-based validation. It should not wait for a spatial manufactured-solution tolerance intended for a different equation.

## 7. D4: test the returned state, not a later repaired state

### 7.1 The proposed assertion is not the current assertion

The current test near `carrier_retry.f90:1263–1300` asks whether one further sweep moves the composition returned by production by more than the chosen bound. The proposed replacement permits up to five extra sweeps and asks whether a later composition has stopped moving.

In mathematical terms, the current assertion concerns `T(x_returned)-x_returned`; the new one concerns `T(x_k)-x_k` after the test has advanced the state to `x_k`. Those are different contracts. Keeping the old test name would conceal the change.

The recorded third-sweep movement of approximately `1e-10` is useful diagnostic evidence. It does not show that the original returned state already satisfies the original assertion. A test-only loop cannot improve the composition that the production caller received.

### 7.2 There is a real production stopping-criterion gap to investigate

`equilibrate_chemistry_at_fixed_conserved_state` in `diffusive_photochemistry.f90:7014–7124` repeatedly performs a chemical sweep and reconstructs pressure/temperature at fixed conserved energy. It rejects nonfinite and thermodynamically inadmissible states and reports failure when its cycle budget is exhausted. Those safeguards should be preserved.

However, its successful termination condition is the relative temperature increment `dT < chem_cycle_tol`. It does not separately require the eliminated chemical abundances or their reaction residuals to be converged at the final returned temperature.

The nearby comment asserts that a temperature change below the reaction-residual tolerance cannot move the composition beyond that tolerance. There is no general implication of this kind. A trace abundance can change with almost no caloric effect, and a temperature-sensitive reaction can amplify a small thermal change. The closure also retains rates and source quantities evaluated during a sweep at the preceding temperature, whose consistency with the returned temperature must be bounded.

This is an implementation-level closure weakness, not proof that the particular `1.034627e-6` test value is a scientifically significant error. The effect on the actual fixture was not remeasured here. Its physical importance and the appropriate norm must be established, rather than inferred from the trace abundance or the small amount by which the test failed.

### 7.3 The existing test norm also needs a clear justification

`ieq_res_tol` is a reaction-residual tolerance. It is not automatically a bound on the relative abundance movement of every species after an additional nonlinear sweep. A residual-to-abundance error estimate requires conditioning information.

Thus two possibilities remain open:

- Production promises a closed thermochemical state and exits prematurely; repair its stopping condition.
- The abundance-movement test assumes a stronger bound than the production residual contract justifies; revise and rename that assertion using a defensible norm, while adding a direct returned-state residual check.

Neither case is resolved by allowing the test to iterate until it passes. Similarly, `f_sp > 1e-20` is a numerical presence rule in mass-scaled number fractions, not a universal physical importance criterion. A low-mass species can dominate electrons or a reaction pathway. The recorded L36d electron contribution of H3+ is a concrete reason not to dismiss the molecular chain simply because HeH+ is small.

### 7.4 Recommended replacement task

1. Save the unmodified production return and all state needed to evaluate it consistently.
2. Measure its chemical residuals at the returned temperature and electron density, with transported fractions held according to their actual ownership.
3. Separately measure one complete fixed-conserved-energy closure-map increment, including the EOS/temperature update, rather than only repeated chemistry at frozen temperature.
4. Define successful closure using thermal consistency and justified chemical residual/increment criteria. Record the binding cell and species and fail explicitly if the iteration budget is exhausted.
5. Verify that failed trials restore composition, rates, caloric state, and boundary context used by the next trial.
6. Keep additional settling sweeps as a separately named diagnostic. Do not overwrite the direct returned-state assertion with their result.

Strict decrease of every successive maximum norm is not a general property of a convergent coupled iteration. The retained arithmetic example has spectral radius `sqrt(0.1) < 1`, yet successive infinity-norm changes begin `1, 10, 0.1, 1`. An unconditional monotonicity requirement could reject legitimate convergence. If contraction is required, state the norm and assumptions that establish it.

## 8. Additional source defect: HeH+ is missing from the helium stage bound

### 8.1 Evidence and counterexample

The following inspected expressions do not describe one common admissible set:

1. `carrier_state` obtains total helium nuclei from `element_nucleus_counts` near `diffusive_photochemistry.f90:2875`. That census includes HeH+ through the species stoichiometry.
2. `carrier_stage_nucleus_density` and `carrier_nucleus_reference`, near lines 4193–4220, normalize He II and He III by that total helium inventory.
3. `carrier_ionization_stage_projection`, near lines 6109–6154, imposes only `x_HeII + x_HeIII <= 1`.
4. `limit_to_element_budget` subsequently checks the carbon, oxygen, and hydrogen carrier constraints but does not reserve the helium in HeH+ from the helium ion-stage sum.
5. `carrier_source` near line 4760 subtracts HeH+ when reconstructing neutral helium, then clips a negative remainder to zero.
6. `carrier_write_back`, near lines 6245–6268, preserves HeH+, writes the ionized stages, and also clips the remaining neutral helium to zero.

For a normalized total helium inventory of one, take

```text
HeH+ helium inventory = 0.10
x(He II)             = 0.60
x(He III)            = 0.35
```

The ion-stage simplex accepts `0.95 <= 1`. The write-back expressions then give neutral helium `max(1-0.60-0.35-0.10, 0)=0`, leaving a total helium inventory of `1.05`.

This counterexample was executed in the retained arithmetic script. It proves that the inspected projection/write-back expressions do not enforce the claimed inventory invariant. It does **not** establish that a catalog run currently reaches this state, or that this defect caused the recorded L36 failure. No production reachability test was run in this review.

### 8.2 Physical correction

When stage fractions use total helium nuclei as their denominator, the constraint must include helium held outside the atomic stages:

\[
x_{\mathrm{HeII}}+x_{\mathrm{HeIII}}
\leq 1-\frac{n_{\mathrm{HeH^+}}}{n_{\mathrm{He,nuc}}}.
\]

Use the actual stoichiometric inventory, including any additional helium-bearing molecules introduced later. A joint projection or a consistently defined available atomic-helium budget can enforce this. If the normalization is changed to an atomic pool instead, its flux and source definitions must change consistently; changing only the limiter denominator is not sufficient.

Apply the same constraint in trial admissibility, finite-difference perturbations, source evaluation, projection, and write-back. Do not repair a negative neutral abundance by clipping while leaving excess positive species untouched.

Tests must measure total H nuclei, total He nuclei, total mass, and charge consistency from the actual returned species. Include nonzero HeH+, nearly exhausted neutral helium, physical cells, and ghosts according to their boundary ownership. An algebraically closed sum of an artificial “remaining stage” and two ion fluxes does not establish these actual species inventories.

## 9. Revised execution plan and decisions

### 9.1 Recommended sequence

1. **Record a stable baseline.** Identify the exact source and executable used for each targeted experiment. The active working tree and historical L37 binaries are not interchangeable.
2. **D5 diagnosis and state contract.** Add comparisons around boundary construction, chemical refresh, EOS installation, residual assembly, and the later face-flux report. Separate stale cached faces from serialization and incomplete ghost chemistry.
3. **Helium admissibility repair and tests.** Establish a common molecular/atomic element budget before accepting altered He II residual norms.
4. **D4 returned-state closure investigation.** Define and test production's actual closure contract; repair production or justify a revised test norm, not extra unreported settling in the assertion.
5. **D3 measurements.** Obtain true reaction-channel production/loss, actual transport balances, conditioning, and stiff manufactured cases. Diagnostic work can proceed earlier when file ownership allows it; the acceptance decision should follow coherent state evaluation.
6. **D2 physical decision.** Choose and document the zero-flow/reversal entropy and composition model. Implement only after its mathematical properties and valid tests are stated.
7. **D1 disclosure and migration.** Precision/provenance work can proceed independently when its shared files are available. Correcting the default and changing the physical base-face geometry are separate decisions.

After a boundary/closure change, run focused caller and downstream-consumer tests first. Since a shared boundary operation reaches several execution routes, those routes deserve broader coverage than a single catalog evaluation. Unrelated products need not be regenerated.

### 9.2 Decisions to present to the user

| Decision | Recommended choice | Evidence required before adoption |
|---|---|---|
| Historical grid reproduction | Pin resolved old grid inputs; use corrected defaults for new configurations after approval. | Coordinate equality and restart compatibility on representative affected grids. |
| D2 boundary at zero/reversed flow | Decide the physical contact/reservoir model before selecting a blending formula. | Characteristic count, discrete hydrostatic and contact tests, weak-flow derivative checks. |
| D3 certification norm | Retain the current gate while diagnosing; introduce a documented multi-criterion rule only if needed. | Dimensional balances, cancellation floor, conditioning, stiff manufactured tests, and structural conservation tests. |
| D4 closure | Require a directly evaluated returned-state contract; retain settling as a separate diagnostic. | Chemical/thermal residuals at the returned state and correct rollback. |
| D5 restart state | Derive boundary ghosts from explicit physical state and inputs; persist only genuinely independent reservoir variables. | Boundary idempotence, same-state operator agreement, and controlled serialization measurements. |
| Catalog replacement | Validate copied products first; re-solve if corrected equations require it. | Same-operator versus changed-operator provenance and a reviewed result, not merely an old certification flag. |

### 9.3 Correct the plan's internal bookkeeping

- The title, opening, “The four,” and final exclusion paragraph still describe four items despite D5.
- D5 allows a reference-product refresh, while the final paragraph says none of the four items refreshes a reference. State the actual policy and the numerical reason for each intentional change.
- Add D5 to the ownership table. Its likely scope includes composition/EOS installation, ionization, initialization, stationary evaluation, and state snapshots, not just `Apply_BC.f90`, `base_boundary.f90`, and `load_IC.f90`.
- D4 is not established to be a one-hour test-only task. D3a is not established to require only test files.
- Sequence shared closure/boundary edits with the ongoing work, not merely D2 against D5. A source file becoming available does not prove that its interface assumptions have remained unchanged.
- Replace promises that old states “now certify” or all flowing controls remain byte-identical with physical acceptance conditions and measured outcomes.

## 10. Validation performed, limitations, and provenance

### Performed

- Direct inspection of the implementation paths cited above, including callers and consumers of the affected values.
- Comparison of the updated D5 summary with the detailed L37 tables and provenance.
- Execution of the retained Python arithmetic checks. These are counterexamples and representation checks, not an EXHALE regression suite.
- Reading of the relevant published Thompson and Käppeli–Mishra material from the local reference PDFs. The publisher records/PDF links are given at the associated physical claims. ADS was attempted first but could not be reached from the shell environment; that limitation was not treated as a successful citation check through ADS.

### Not performed

- No new production executable was built and no full regression suite was run.
- No catalog state was evolved, remapped, rewritten, or re-certified.
- L35/L36/L37 numerical values were not remeasured. They are explicitly treated as recorded results, including where their summaries disagree with their tables.
- No production test of the helium counterexample's reachability or attribution of D5's discrepancy among its candidate causes was run.
- No claim is made that an unperformed test passed or that a current source snapshot reproduces an archived executable's numbers.

This scope is appropriate for a plan review in a changing working tree: it establishes incorrect premises and concrete source-level risks without silently changing the system being reviewed. The implementation tasks above include the targeted execution needed before accepting a fix.

### Snapshot identifiers

The repository HEAD at inspection was `3c73905ca8a7fe2af92a5c2b014c225b2feedff3`, with extensive existing tracked and untracked changes. HEAD alone does not identify the inspected implementation.

Selected SHA-256 values checked during the review:

```text
PLAN_20260918.md
be4184564f966bd5c38cad42c5617571d953768c7a1d2a86d1057c1444e16e69

src/modules/states/base_boundary.f90
a2c46b4f1fbfc16879aa3bda684a303d2ffa8d6e0ee16bfc562eebd86fdfcc6d

src/modules/lower_atmosphere/diffusive_photochemistry.f90
08eef7af73f259baede1774fddc45caffb10f5cd96c00fa453630e38c4358601

src/modules/time_step/certification.f90
0dcfbd3ce079891a6aefc41706d3eeafba110722c89278c1224b82debc5755f9

src/tests/ionization_stage_flux/ionization_stage_flux_tests.f90
51ff8eda62b386f19298fec616a79bd5e6c25241b74390fe947377252e0bee2b
```

`carrier_retry.f90` changed during the review. Its returned-state assertion was reread after the change and was still present; the subsequent hash was `4407b5205d0defc40ae3937c376866e47e616a93cebf378710150581de09e166`. Line numbers in this review are navigation aids for that working interval, not an assertion that other ongoing work has stopped.

## Conclusion

D5 is worth addressing promptly, but as a coupled boundary-state and evaluation-consistency problem, not a guaranteed two-file restart repair. D2 needs a corrected physical premise; D3 needs a demonstrated balance and conditioning analysis; D4 must continue to judge the state that production actually returns. The newly identified helium molecular-inventory constraint belongs ahead of any relaxation of helium acceptance. These changes would make the plan a defensible path toward physical and numerical correctness rather than a sequence aimed at recovering existing pass/fail outcomes.
