# Review of session_handoff_20260915_rev1.md

Review date: 2026-09-16 (KST)

## 1. Conclusions and scope

The failed cases deserve further investigation, but the handoff does not yet establish their causes. The most useful next step is a short, instrumented comparison at a saved failing state, not an unconditional extension of the iteration budget.

The current implementation has two confirmed consistency problems: the stationary and physical-time routes use different base-boundary operators, and the outer progress controller treats a smaller constrained composition update as evidence of contraction. There is also a boundary derivative omission in the approximate carrier Jacobian under specific flow/boundary conditions. These findings are distinct from the unproven claim that an outer-boundary condition causes the well-mixed failures.

This review inspected the current source, callers, selected existing run logs, and local published papers. It did not rebuild an executable, start a simulation, rerun a regression suite, or modify production code. Numerical results below are **read from existing logs**, not independently reproduced solver measurements. Equations and limiting arguments are analytical checks. Proposed experiments have not been performed.

The repository HEAD was `43bc28cef58772560bae019559f3592335922212`, with extensive uncommitted work. HEAD alone does not identify the inspected implementation. The campaign executable identified in the handoff is `db87b88d1ce53facf1d61084fa535ca5` (MD5); several current source corrections postdate it. Current-source findings must not automatically be attributed to that executable. Conversely, a correction present in the tree does not repair its archived results.

## 2. Findings ranked by consequence

| ID | Finding | Evidence status | Consequence |
|---|---|---|---|
| R1 | Base boundary depends on whether the stationary route is active | Confirmed in implementation and callers | Certification is not certification of the physical-time operator |
| R2 | A smaller limited composition update is used as a progress surrogate | Confirmed control logic; effect on a specific refusal not isolated | False progress or premature refusal is possible |
| R3 | Global particle-count bound can couple distant carrier regions | Confirmed algorithm; dominance in these cases unproven | Local stiffness can restrict useful advance elsewhere |
| R4 | Approximate carrier Jacobian omits copied-ghost dependence in certain boundary branches | Confirmed algebraic omission; triggering flow in failing cases unmeasured | Less accurate boundary linearization, potentially worse convergence |
| R5 | Outer-boundary diagnosis exceeds the available evidence | Logs locate a maximum, not its cause | Changing the boundary prematurely can change the physical problem |
| R6 | R23 photon subtraction is present, but internal excitation is still deposited immediately | Confirmed model approximation | A physically justified thermalization regime remains necessary |
| R7 | Eight-thread nonreproducibility has not been localized to a reduction | Existing handoff reports divergence; cause unverified | A race, statefulness, and rounding amplification remain distinct possibilities |

No evidence in this review shows that R2 or R4 bypasses final certification. A poor progress decision and a false successful certification are different defects.

## 3. R1: stationary and physical-time base boundaries differ

### Implementation evidence

`src/modules/states/base_boundary.f90:539–558` explicitly describes the distinction, and the executable branch at approximately lines 572–583 implements it. The distant wind-flux discriminant is used only when `base_branch_on_wind_flux`, `base_branch_stationary_route`, and `have_F` are all true. Otherwise the local discriminant controls the entropy blend.

This is not an unused option. `src/EXHALE_main.f90:1585–1589` and `3899–3903` enable and disable the stationary-route flag around stationary solves. Diagnostic/certification contexts also enable it, including lines 1169 and 5263.

The newer implementation blends branch **weights** smoothly. The previous review's discontinuous hard-switch finding should therefore not be repeated as an unchanged current defect. Smoothness, however, does not make the two boundary operators identical.

### Physical judgment

If the stated objective is the steady state of the same physical model used for time evolution, selecting a different boundary condition by solver route violates that contract. Let the two spatial residuals be `R_stat(U)` and `R_time(U)`. A successful check of `R_stat(U*) = 0` does not establish `R_time(U*) = 0`.

At a continuum steady state, conserved mass flux can justify using a distant flux as information about the base. It does not prove equivalence of the two discrete entropy selectors when the local cell velocity contains the very oscillatory mode that motivated the change. The source acknowledges this distinction.

The stationary boundary can instead be treated as an explicitly different model. That is an honest interim interpretation, but its certification must be labeled accordingly. Small changes in mass loss or agreement in an observable do not establish physical equivalence.

### Recommendation

First evaluate both boundary constructions on exactly the same saved state. Record ghost density, pressure, temperature, entropy source, base mass/energy flux, and the first few residual rows. No long time integration is required for this comparison.

For a durable correction, define a single physical reservoir/characteristic boundary and use it in both routes. Address low-Mach face-velocity errors in the discrete flux or boundary reconstruction itself. A distant flux may inform initialization or a preconditioner, but should not silently replace the physical boundary according to solver choice. A redesign of the boundary API requires a separate implementation decision.

The distant-flux dependence also introduces nonlocal Jacobian coupling. If retained, the preconditioner should approximate that coupling, for example through a small explicit boundary block. This is a numerical recommendation, not a justification of the boundary model.

## 4. R2–R3: the composition update is not a fixed-point defect

### What the code measures

In `src/EXHALE_main.f90:6636–6674`, `composition_closing` can become true when the hydro distance is smaller than the joint distance and `comp_drift_last < comp_drift_before`. This can prevent the no-progress counter from increasing even when the joint residual does not fall. The controller also changes `comp_omega` and the carrier movement bound. The drift history is updated around lines 6779–6780.

In `relax_photochemical_composition`, the returned carrier drift is the largest accepted change in a solved carrier fraction, divided by `x_ref` (`diffusive_photochemistry.f90:5869–5882`). It is a displacement between the entry and exit compositions of a bounded, finite-budget relaxation. It is not the residual of the unlimited fixed-point map.

For example, with `delta f = alpha d`, reducing `alpha` reduces `||delta f||` even when the direction `d` and physical residual do not improve. Thus the current comparison cannot by itself establish contraction. Conversely, an active movement bound can hold the update size nearly constant during useful redistribution, making displacement an unreliable stagnation test as well.

### Global bound and its diagnostic distinction

`carrier_particle_count_change` at lines 5324–5340 computes

`max_j |n_particles,new(j) - n_particles,entry(j)| / n_particles,entry(j)`.

The reference is the composition at entry to the relaxation, not the preceding accepted substep. At lines 5800–5855 a violating trial is rejected, the trial state is restored, and a single scalar `grow` is reduced for the entire grid. The local base intervals differ, but all receive this same multiplier. This is a global excursion budget, not merely a bound on each successive step.

That design can protect the fixed-hydrodynamic-state approximation. It can also make an outer chemical region advance only as far as a different region permits. The source establishes the possibility; the inspected logs do not isolate it as the dominant cause.

Importantly, `carrier drift worst at cell ...` identifies the largest carrier-fraction displacement (`pct_drift_j`). It is **not** necessarily the cell controlling the particle-count bound (`bound_last_j`), and neither is necessarily the maximum-residual cell. These three indices must be recorded separately.

The existing debug probes are fixed at radii 1.20, 1.36, and 1.60 (`diffusive_photochemistry.f90:5727,5890–5905`). They are insufficient to diagnose a row at the outer edge without additional dynamically selected probes.

### Recommended control changes

1. Measure residual change with consistent component scales before and after a complete coupled update. Log the unscaled balance terms as well.
2. Record movement limits, accepted intervals, and whether each update is constraint-limited. Do not interpret a reduction caused by a tighter limit as physical contraction.
3. Compare fixed-point defects only at a common relaxation definition and budget. Residual-based actual/predicted reduction is preferable when a usable local model exists.
4. If the global bound is demonstrated to bind at an unrelated location, test regional relaxation or a coupled carrier/thermal block. Shared face fluxes must remain conservative; independent cell updates without matching face fluxes are not an acceptable repair.
5. Retain strict final residual, positivity, elemental, charge, and mass-closure checks. Improve the path to a solution rather than weakening acceptance.

## 5. R4–R5: what the outer carrier boundary actually does

### Flux accounting

`carrier_residual` initializes `Jf = 0` and fills only faces `1:N-1` (`diffusive_photochemistry.f90:4289–4306`). `carrier_face_coefficients` likewise fills interior faces. Thus the outer **diffusive** flux is zero. The interior expression contains both a gradient contribution and a gravitational drift contribution (`carrier_face_flux`, lines 4154–4177).

The advective contribution is separate. `carrier_advective_divergence` (lines 3452–3488) constructs mass fractions, calls `species_face_fraction`, forms species face fluxes from `Frho`, and takes their divergence. `species_face_fraction` in `src/modules/flux/species_face_flux.f90:106–151` chooses reconstructed donor values according to the sign of the mass flux. Outer carrier ghosts copy the last physical carrier fraction, including line 4702.

Therefore the implementation does **not** impose zero total molecular flux at the outer boundary. Advective escape is allowed. Zero diffusive flux is also not interchangeable with zero mixing-ratio gradient in the presence of the drift term.

For a carrier with number density `n_c`, the relevant steady balance is schematically

`div(n_c v + J_c) = P_c - L_c`.

A large last-cell residual can arise from its advective flux difference, the transition to the imposed diffusive face flux, chemistry, radiation/shielding, or inadequate relaxation. Its index alone cannot distinguish these causes. Nor does a moving maximum prove the propagation of a physical front: the identity of the largest row can change among several spatial structures.

### Conditional Jacobian omission

The carrier matrix deliberately uses a donor-cell approximation to the full reconstructed advective operator (`diffusive_photochemistry.f90:4610–4653`). Omitting higher-order limiter derivatives is documented and is not itself an exact-Jacobian bug.

There is a more specific issue. For `j=N` with `Frho(N)<0`, the matrix inserts neither the usual diagonal donor term nor a neighboring term, because the latter requires `j<N`. Yet the outer ghost fraction is copied from `f_c(N)`. At fixed `msum`, its donor mass fraction is `m_c f_c(N)/msum(N+1)`, so its derivative with respect to `f_c(N)` is nonzero. This dependence is missing even from the intended first-order approximation.

Similarly, `carrier_face_mass_fraction` (lines 3524–3537) assigns `Y(ghost)=Y(1)` when no carrier base composition is imposed. In that branch, a positive inner face flux depends on cell 1, whereas the matrix's inner-boundary comment treats the ghost as independent data and supplies a diagonal only for negative inner flux.

These omissions should be repaired by deriving boundary derivatives from the same imposed-versus-copied ghost rule used by the residual. They need not change the root of a fully converged solve. This review did not measure the sign of `Frho(N)` or establish activation of the unimposed-base branch in the reported failures; neither omission is claimed as their demonstrated cause.

A focused validation should compare directional derivatives for positive and negative face fluxes, imposed and copied base composition, and constant/linear carrier profiles. Test the donor-cell contribution separately before interpreting differences caused by the intentionally deferred reconstruction terms.

## 6. Case-by-case reassessment

### 6.1 Well-mixed He/H = 0.083

The existing `LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.083/run.log` directly records:

| Outer pass | Maximum gated species row | Cell | Mass row | Energy row |
|---|---:|---:|---:|---:|
| 21 | 2.42e-2 | 500 | 1.25e-8 | 3.25e-8 |
| 22 | 2.43e-2 | 500 | 6.36e-9 | 1.66e-8 |
| 23 | 2.45e-2 | 500 | 2.49e-9 | 1.15e-8 |
| 24 | 2.46e-2 | 500 | 1.06e-8 | 3.11e-8 |

These entries occur at lines 4188, 4427, 4610, and 4827. Line 4833 records refusal after three passes without the required progress. Nearby displacement diagnostics identify cell 231 at radius 1.257. That is evidence that the displacement and residual maxima differ, not proof that cell 231 controls the particle-count bound.

The residual is slightly increasing in this interval. A longer identical run is not justified by an assumed monotone decrease. First instrument cell 500, its neighbor, the actual bound-controlling cell, and the displacement maximum. Dump both face fluxes, production/loss, the row denominator, accepted physical interval, temperature, radiation rates, and the boundary mass-flux sign.

Only after this decomposition should a domain extension or physically justified outer diffusive boundary be compared. Keep the shared interior mesh and physical reservoir fixed when possible. Do not remove cell 500 from certification or impose local chemical equilibrium there merely to make the residual smaller.

### 6.2 Well-mixed He/H = 0.55

The corresponding existing log records passes 31–35 at lines 9154, 9329, 9500, 9664, and 9828. The species maxima are `0.0227, 0.0227, 0.0229, 0.0230, 0.0232`, at cells `294, 295, 297, 298, 299`.

The maximum moves outward, but the residual is flat and then rising. These cells are not the last cell. The addendum's description of a feature passing through the outer boundary is therefore not established by these diagnostics. Full profiles are needed to distinguish transport of a feature, a switching maximum, and insufficient coupled relaxation.

Apply the same instrumented comparison as above. Extend the pass budget only if a fixed-definition residual or a resolved profile diagnostic shows meaningful improvement, not merely because the maximum's index increases.

### 6.3 Kzz = 1e9, He/H = 9.7

The handoff reports first-pass hydrodynamic residuals of order unity and an inadmissible element-relaxation attempt after mapping a molecular He/H = 2.13 seed. Those specific numbers are quoted from the handoff here, not remeasured.

Trying the newly certified atomic state at the same elemental abundance is a reasonable seed comparison: it separates a large change in elemental inventory from activation of molecular chemistry. It is not guaranteed to solve the molecular branch. The previously successful 0.3 dex ladder spacing is an empirical continuation choice, not a physical admissibility limit.

Before any further continuation, identify the exact element-relaxation rejection condition and verify density, conserved energy, element budgets, charge, and thermal admissibility after mapping. If a homotopy is needed, use consistent stoichiometry and energy accounting, and certify only the complete target model.

### 6.4 Kzz = 1e9, He/H = 0.083

The handoff identifies elemental transport, not H2 transport, as the later limiting row. This requires a different diagnostic emphasis. The initial 1.4 dex abundance change is a plausible contributor, not an established explanation of late-time convergence.

An own-abundance atomic seed and intermediate molecular abundance steps provide useful controlled comparisons. Check elemental flux balance and the element under-relaxation response before altering H2 boundaries.

The quoted residual sequence includes an increase from pass 5 to pass 10 and should not be called globally monotone. Even a constant 2% decrease from `8.35e-4` would require approximately `ln(1e-5/8.35e-4)/ln(0.98)`, or 219 **additional** passes. This is an extrapolation, not a forecast; it is particularly unreliable when the active row or movement bound changes.

## 7. Physical chemistry: corrected energy accounting is not complete thermalization physics

The R23 photon subtraction is present in `src/modules/lower_atmosphere/molecular_reaction_heat.f90:418–449`, and the heating term uses the corrected quantity. It would be incorrect to report the old full reaction-energy deposition as unchanged in the current source.

The published Boehringer and Arnold paper was inspected from `../references/Bohringer_1986JCP_84_1459.pdf`, including the experimental method and discussion on pages 1459–1461. It measures the reaction in a selected-ion drift tube; the discussion of the approximate 153 nm photon and preferential vibrational excitation summarizes a proposed radiative mechanism. It is not a measurement of thermalization efficiency in an escaping planetary atmosphere. The similarly named file ending in `2097.pdf` concerns a different reaction and is not the source for this claim.

The current source immediately deposits the post-photon remainder and explicitly acknowledges unresolved H2+ vibrational excitation. This still needs a physical validity criterion. Removing photon energy is necessary, but internal excitation becomes translational heat only through the relevant collisional/reactive channels. Rapid destruction can transfer that energy to other products rather than directly to thermal motion.

Use an energy ledger separating escaping radiation, product translation, internal excitation, and subsequent transfer. A reduced model may use an effective thermalization fraction if its collisional, radiative, reactive, and advective timescale conditions are stated and checked. The analogous R15 formation-state question remains open; an enthalpy identity alone does not determine the immediate heat fraction.

The R23-corrected case agreeing in an observable to four digits, as reported by the handoff, is an impact observation. It is not a physical justification for neglecting excitation or a validation of all molecular cases.

The chemistry closure routine also requires careful terminology. `equilibrate_chemistry_at_fixed_conserved_state` (approximately lines 5510–5610) rejects nonfinite states and tests thermal admissibility and temperature stabilization. It deliberately does not reject every off-simplex root-search event. That is reasonable when the returned state is repaired and its balance is separately checked. However, a finite, temperature-stable returned state is not by itself proof of a chemical root. Keep accepted-state reaction residuals and closure status distinct in diagnostics and final certification; do not treat either the search-excursion count or its absence as decisive.

## 8. Low-XUV pseudo-time behavior

The handoff's improvement after changing `EXHALE_PTC_DTAU0` from 1 to 1e8 is evidence about those recorded runs, not proof that the explicit stability interval is the cause of the earlier stalls.

Kelley and Keyes, *SIAM Journal on Numerical Analysis* 35, 508–523 (1998), was inspected in its published local PDF, particularly section 1.1 and Algorithm 1.1. Its pseudo-transient correction has the form

`(V/delta + F'(x)) s = -F(x)`.

Small and large pseudo-time intervals change the regularization relative to the Jacobian. They are not interchangeable physical timesteps, and a large interval can help an already useful seed by approaching a Newton correction. The paper does not justify a universal value of 1e8 for this implementation or diagnose its cell-500 energy row.

A focused comparison should keep seed, executable, thread count, residual definition, and boundary fixed while varying the initial interval. Record the diagonal shift relative to the local Jacobian/preconditioner scales, accepted line-search factors, Krylov work, and unscaled energy terms. Attribute the improvement only after separating these effects. The x0.01 failures remain failures of attempted routes, not a proof that the physical solution does not exist or that one particular flux redesign is uniquely necessary.

## 9. Parallel reproducibility and regression interpretation

The handoff's phrase identifying an unlocated order-dependent reduction is too specific. Inspection of the ionization sweep's OpenMP clauses (`src/modules/radiation/ionization_equilibrium.f90:1702–1706`) shows reductions involving counters, timings, and maxima. Their presence does not identify a floating-point sum that perturbs the physical state. This is not an exhaustive race audit or proof of thread safety.

Start with repeated residual evaluations of the same frozen state, first at one thread and then at eight, with state restoration and explicit cache policy. Compare returned physical arrays before comparing nonlinear iteration histories. If the residual already differs, localize its earliest changed component. If it agrees but Krylov steps differ, inspect linear algebra ordering and thread configuration. Also inspect saved scratch state, initialization, and mutable caches; do not presume harmless rounding until shared writes have been excluded.

Re-evaluation of a saved solution at one thread, repeatability of the solve at eight threads, and agreement of independent converged solutions are three different tests. None substitutes for the others. Exact iteration-history identity is not the physical acceptance criterion, but unexplained state differences require investigation.

Regression differences should be attributed on an identified executable before replacing expected outputs. A zero-versus-nonzero ghost species value needs its boundary and loader semantics checked, not merely a comparator-tolerance change. A relative difference of one against a zero reference is not automatically a 100% error in a physically important quantity; report absolute value, units, floor, and physical constraint. Conversely, a small difference does not excuse incorrect physics.

## 10. Recommended bounded experiment sequence

These are proposals, not experiments performed for this review.

| Order | Comparison | Main question | Decision criterion |
|---|---|---|---|
| 1 | Same saved state, stationary versus physical-time boundary | Are the certified and evolved operators consistent? | Compare base fluxes and residuals against the declared physical boundary |
| 2 | Frozen well-mixed failing state, full row decomposition | What actually dominates the carrier residual? | Separate advection, diffusion, source, scaling, and limiting cell |
| 3 | Repeated frozen-state residuals at one and eight threads | Does nonreproducibility start inside residual evaluation? | Localize the first changed physical array before new long solves |
| 4 | Boundary directional-derivative probes | Are copied-ghost derivatives represented? | Check first-order terms separately from deferred reconstruction |
| 5 | Fixed seed, controlled composition-bound comparison | Does global limiting suppress the relevant region? | Require residual reduction and conservation, not smaller updates alone |
| 6 | Own-abundance seed versus mapped seed for Kzz failures | Is the dominant issue initialization? | Compare rejection reasons and consistent full residuals |
| 7 | Domain/boundary sensitivity only if step 2 implicates the edge | Is truncation physically affecting the solution? | Compare shared-domain profiles, fluxes, and observables |

Keep the campaign executable and archived products untouched. Any future diagnostic source should be retained under `docs/audit_20260905/`, with exact inputs, executable identity, thread settings, and output locations. Private build directories are appropriate only after implementation/testing is authorized. This review created no diagnostic code requiring archival.

The first implementation priorities should be a consistent physical boundary contract, honest residual-based progress diagnostics, and the boundary linearization corrections. A stronger coupled solve is a subsequent option if targeted measurements show that partitioned carrier/thermal coupling is the actual bottleneck. Increasing pass counts, loosening tolerances, suppressing boundary rows, or changing physical rates to obtain convergence are not substitutes.

## 11. Evidence inventory and limitations

Primary reviewed paths:

- `docs/session_handoff_20260915_rev1.md`, especially sections 12 and 16.
- `src/EXHALE_main.f90`: stationary boundary context and partitioned progress control.
- `src/modules/states/base_boundary.f90`: entropy-branch discriminants and blend.
- `src/modules/lower_atmosphere/diffusive_photochemistry.f90`: carrier fluxes, matrix entries, copied ghosts, movement bound, and chemistry closure.
- `src/modules/flux/species_face_flux.f90`: reconstructed advective species flux.
- `src/modules/lower_atmosphere/molecular_reaction_heat.f90`: R23 heat recipient.
- `src/modules/radiation/ionization_equilibrium.f90`: relevant parallel reduction declarations.
- The two well-mixed run logs cited in section 6.
- `LHS1140b/models/BINARY_MANIFEST_db87b88d1ce5.txt`: retrospective source manifest and explicit later-edit warning.

SHA-256 identities read during this review:

```text
d5ebed732dd5dc30383f96b0aa5e6007721ece0075fc2c3f8bc20a3881d41ce1  src/EXHALE_main.f90
f4f66400b5e119d5d9a10689ab02c20b9a063f4ee0e1ea13df861e7299e2787b  src/modules/states/base_boundary.f90
f26cc47bb75ef5f30ee035c3aa8091819159a0314a1ed71dc78600b492f56755  src/modules/lower_atmosphere/diffusive_photochemistry.f90
b343039eab5a1f4bb19e03f9022fb3750045256ad9bfd02adb299ce21a1757f1  src/modules/lower_atmosphere/molecular_reaction_heat.f90
```

Line references apply to that inspected working tree and can move as concurrent development continues. Existing logs may also continue to grow. The manifest was written after the campaign build; its timestamp-based explanation is useful provenance but is not equivalent to a complete immutable prebuild source snapshot and build transcript.

Local published references read for the specific checks above:

- Boehringer and Arnold (1986), *Temperature and pressure dependence of the reaction of He+ ions with H2*, JCP 84, 1459–1462, DOI `10.1063/1.450490`; `../../references/Bohringer_1986JCP_84_1459.pdf` relative to this document.
- Kelley and Keyes (1998), *Convergence Analysis of Pseudo-Transient Continuation*, SIAM J. Numer. Anal. 35, 508–523; `../../references/Kelley_1998SIAMJNA_35_508.pdf` relative to this document.

This is a focused source and evidence review, not a new end-to-end validation or an exhaustive audit of the complete molecular network. No unperformed experiment is counted as passed, and no source-only finding establishes the precise contents or behavior of an older executable.
