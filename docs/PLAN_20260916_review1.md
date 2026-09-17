# Second review: PLAN_20260916_rev1.md

Date: 2026-09-16 (KST)

## 1. Judgment

Revision 1 is substantially better than the original plan. It correctly withdraws the modified momentum-certification proposal, separates steady pseudo-time iteration from physical subcycling, recognizes the limits of the parallel evidence, and puts the molecular energy model behind a data-inventory gate. The metadata and output-path contracts are also much clearer.

**Proceed with the diagnostic work, but revise L26 before implementing its proposed repair.** Defining the zero-window weight as one half and replacing the range by a normalized RMS do not, in general, make the boundary continuous. This is the principal remaining mathematical defect in the revision.

Other required clarifications concern the L22 diagnostic partition and progress metric, the sign of a replacement dissipation operator, the comparison design for switching from HLLC to Roe, and acceptance tests that distinguish shared-face bookkeeping from full inventory conservation.

### Evidence and scope

This review re-read the revision and the implementing paths in the current main program, base boundary, carrier relaxation, element relaxation, numerical flux, dissipation module, and output/post-processing routines. The previous literature check against the supplied published Rieper papers remains applicable. No additional external literature is needed to establish the mathematical issues identified here.

The executable MD5 is still `c2e9c9990b9f14f1be8cd77abca68945`. Directly checked source identities are:

```text
8e818aa4699a2c8708e2da0f7ed77e84  src/EXHALE_main.f90
b445869ac32112a8498f1117c080af8e  src/modules/states/base_boundary.f90
ca4db5847c20101622a819c84f4673dc  src/modules/lower_atmosphere/diffusive_photochemistry.f90
```

The reviewed plan SHA-256 is `1afde61b0293b2cff0861bb44615ea1600bac8f6a65491074c49ccc122994e7f`. These checks identify selected files; they are not a complete build reconstruction.

One small analytical probe was executed and retained at `docs/audit_20260905/plan_20260916_rev1_boundary_limit.py`. It evaluates the proposed normalized-RMS boundary formula, not EXHALE. No production executable was rebuilt or run, no campaign was restarted, and no numerical products or expected regression outputs were changed.

## 2. Correction to my previous review: carrier progress already uses a residual

The revision's qualification that the carrier half has already been corrected is supported by the current code. My previous review described the progress path too broadly as still using the bounded carrier displacement. That description is incorrect for the default carrier path in the reviewed build and should not be carried forward.

`src/EXHALE_main.f90:6762–6796` distinguishes the diagnostic `carrier_drift` from `carrier_close`. By default, when `bg_ready` is true, it saves the carrier module state, calls `carrier_steady_residual` on the returned composition, restores the module state, and assigns that residual to `carrier_close`. The old displacement behavior is selected by `EXHALE_CARRIER_DRIFT_IS_DISPLACEMENT=1`. The carrier residual scale excludes the pseudo-time term (`diffusive_photochemistry.f90:3198–3223`).

The element path is different. In `binary_element_diffusion.f90:1246–1266`, `drift = element_composition_distance(Ynow,Ypass,Yres)` is formed **before** applying `omega`. Reducing `omega` therefore does not directly scale this diagnostic down. Nevertheless, it is the displacement toward the result of a finite inner relaxation, not necessarily the residual of the state actually returned after damping. The routine can exit with a step-budget outcome, so its relaxed endpoint is not always a solved fixed point.

Consequently, the remaining L22 task is not simply to replace all bounded updates by an unbounded update:

- Distinguish carrier residual, element inner-map displacement, actual accepted displacement, and the full coupled residual in both names and logs.
- Decide whether the element progress measure should be a physical transport residual at the returned state, or a fixed-map defect with an explicitly fixed inner accuracy and budget.
- Normalize the quantities before taking their maximum; two dimensionless values with different meanings and reference scales are not automatically comparable progress measures.
- Keep the physical residual scales independent of iteration length. If scales depend on the state, also report absolute terms or a fixed-reference progress measure so scale movement is visible.
- A missing background or nonfinite diagnostic must be reported as unavailable/invalid, not silently interpreted through the old displacement measure. The current default assignment to `carrier_drift` and the `crc_max == crc_max` check deserve an explicit policy.

The header-level statement in rev1 that the wind response is independent of particle count should also retain its condition: `p=(gamma-1)e_th` is independent of particle count **at fixed gamma and fixed thermal energy**. The implemented molecular caloric EOS changes heat capacity with composition, so gamma is not generally fixed. Monitoring pressure, temperature, sources, and the coupled residual remains the correct recommendation.

## 3. L26: the proposed zero-flux repair is not sufficient

### 3.1 There is no unique limit to assign at the zero window

The source computes a relative spread from the distant window fluxes and blends a local branch weight `w_i` with a wind branch weight. In schematic notation,

```text
w = s(d) w_wind(mean(F)) + (1-s(d)) w_i,
d = spread(F) / abs(mean(F)).
```

The proposed RMS changes the numerator but leaves the scale-free dependence on the direction of the flux vector. Hold the local state fixed with `w_i != 1/2`, and approach the zero vector along two paths:

1. `F = epsilon (1,1,...,1)`: RMS spread is zero; the window weight is one. As epsilon tends to zero, `w` tends to `1/2`.
2. A nonuniform vector `F = epsilon a` whose relative spread is above the rejection threshold: the window weight is zero, independently of epsilon. Now `w` tends to `w_i`.

Both paths approach the same zero window. Unless the local weight happens to equal one half, their limits differ. Changing the formula at the single zero point cannot repair that discontinuity.

The retained probe uses two representative flux values, the proposed relative RMS, the source's cubic branch function, threshold `0.01`, Mach blend width `1e-8`, and fixed `w_i=0`. Its directly executed output is:

| epsilon | Uniform path `(epsilon,epsilon)` | Nonuniform path `(epsilon,2 epsilon)` |
|---:|---:|---:|
| 1e-10 | 0.492500250000 | 0.000000000000 |
| 1e-12 | 0.499925000000 | 0.000000000000 |
| 1e-14 | 0.499999250000 | 0.000000000000 |
| 1e-16 | 0.499999992500 | 0.000000000000 |

This is a direct check of an analytical counterexample, not a simulation or a measured density jump in a catalog case. A density jump additionally requires the local and reservoir entropy constructions to give different densities.

### 3.2 A completely static state is a special case

Rev1's phrase that the boundary is discontinuous on the class containing a hydrostatic initial condition is too broad. At a completely static local base, `M_i=0` and `w_i=1/2`, the two weight limits above coincide. The demonstrated discontinuity is at a zero **distant window** with a different local branch weight. That is still a relevant admissible transient or nonlinear iterate, but it should not be confused with every fully static initial condition.

A fully static state also needs derivative tests; continuity of the value at that particular state is not a proof of differentiability in its neighborhood.

### 3.3 RMS is not automatically smooth

A Euclidean RMS norm is not differentiable at its zero vector. Dividing it by the absolute mean retains a singular normalization. Squared variance avoids the square-root cusp, but dividing by a vanishing squared mean remains problematic. Smoothness of a scalar cubic blend alone does not settle smoothness of the composite state-dependent operator.

Replacing a range by RMS also changes the meaning of the calibrated threshold. An isolated discrepant cell contributes approximately its departure divided by the square root of the number of cells to an RMS, whereas the range detects its full departure. The previously measured separation of catalog values cannot simply be reused as an equivalent physical criterion.

### 3.4 Recommended revision

Make the L26 repair a design decision after the probes, not the predetermined assignment of a half weight at zero:

- Prefer a base closure based on the physical reservoir and local characteristic information, with an appropriate face reconstruction or local boundary solve. This is the long-term physical repair.
- If a nonlocal diagnostic is temporarily retained, require its influence to vanish smoothly and uniformly as the magnitude of the entire window flux tends to zero, not only along the uniform path. Specify any regularization scale and test sensitivity to it.
- A regularized squared coherence measure and an amplitude weight can be investigated mathematically, but are not by themselves a physical justification for distant control of a transient boundary.
- Test directional limits along several nonuniform shapes, local/window sign disagreement, zero mean with nonzero variance, a fully static state, and the window extrema transitions of the current operator.
- Require a transient test for the stated physical boundary contract. Rev1 mentions a transient probe but lists only tests (i)–(iii) in its test suite; include (iv) explicitly.

The statement that regression bases must remain byte-identical because they lie outside the old gate is conditional. A new norm or amplitude criterion can activate different states. First specify the correct boundary, then verify where the new gate is inactive. Output identity is not a reason to preserve an incorrect zero-flux branch.

## 4. L22: remaining problems in the proposed diagnostic and solver contract

### 4.1 The proposed excluded region can contain the whole grid

Step 2 defines the frozen front region by `x2 > 1e-2`, while explicitly noting that in one target state `x2` never falls below that value. The complement is then empty. The multiplier chosen by the tightest cell outside the front is undefined. Calling the partition diagnostic does not remove this problem.

Instead, identify the actual bound-controlling cells from the rejected trial, freeze a specified mask for the controlled experiment, and report its membership. Compare the full set, a narrow mask around the active region, and a separate far-wind mask. If no cell remains in the comparison set, stop with an explicit diagnostic outcome rather than an arbitrary unlimited interval. Keep excluded cells out of the interval-selection experiment only, not out of final physical certification.

The evidence paragraph should also separate the far-wind bound location for the Kzz case from the two well-mixed cases. Do not collapse three different profiles into a single front/wind picture.

### 4.2 Local base intervals already exist

`relax_photochemical_composition` already computes a cell array `dt_code(j)` from diffusion and advection scales (`diffusive_photochemistry.f90:5876–5890`) and passes `dt_code*grow` to the transport routine. The new work is independent control of the multiplier/acceptance, not introducing cell-local intervals from nothing.

The cumulative bound references `nsum0` at relaxation entry. Once a cell reaches that excursion boundary, even an arbitrarily short step may point outside it. Local pseudo-time adjustment alone does not enlarge the feasible set. The plan must define what happens then: return to the hydro solve, accept another globally admissible coupled step, or use a coupled correction. Repeatedly shrinking that cell's interval toward zero must not be mistaken for convergence.

In particular, the statement that positive local pseudo-time factors preserve roots assumes the factors do not become a substitute for testing the original residual. An effectively frozen row can stop updating while its physical residual remains nonzero. Require the original carrier and coupled residuals to pass, and distinguish solver exhaustion from convergence.

### 4.3 Conservation tests need more than a shared-face comparison

The proposed mismatch divided by the larger face flux is undefined when both fluxes vanish. Define a zero case or a physically scaled absolute tolerance. If one shared flux object is used, comparing its two identical uses is almost a bookkeeping tautology.

For physical subcycling, test the inventory equation with cell volumes, accumulated boundary transfers, and integrated sources over a common interval. For stationary pseudo-time iteration, test the unmodified conservative residual and elemental/charge identities at convergence; do not impose a physical-time inventory constraint on intermediate pseudo-time updates.

Failure of the first local pseudo-time experiment also does not automatically justify physical subcycling. These are different algorithms for different objectives. Compare a coupled block or preconditioner improvement against subcycling using measured stiffness/coupling before selecting a larger implementation.

## 5. L25: the revised direction is sound, with four additions

### 5.1 Use a three-way comparison when changing flux families

The original failures use HLLC, while the proposed Rieper correction is on the existing Roe branch. Compare:

1. Unmodified HLLC.
2. Unmodified Roe.
3. Roe with the velocity-jump correction.

Keep boundary, reconstruction, seed, source terms, and solver settings fixed. Otherwise an improvement attributed to the Mach factor may instead be caused by changing flux family, entropy treatment, or admissibility handling.

The two acoustic expansion coefficients at `Num_Fluxes.f90:324–331` are indeed the relevant place to inspect. As verified previously in Rieper (2011), this is an accuracy correction, not a proof that the stiff nonlinear stationary problem will converge.

### 5.2 State the gate without a dimensional ambiguity

The proposed criterion compares face Mach with truncation of a pressure departure. Mach is dimensionless; a pressure error is not until a reference is supplied. If the intended comparison is with `delta_p/p_ref`, justify the asymptotic power of Mach being compared and derive it for the gravitational equilibrium and caloric EOS.

Also separate the published `min(Ma_face,1)` factor, which operates throughout subsonic regions, from an additional base-layer mask. The latter is a new algorithmic choice that needs its own consistency and resolution study. The revised plan should not imply that the published formula itself is confined to the base.

### 5.3 A composition of negative-semidefinite operators has the wrong sign unless specified

Rev1 proposes a weighted composition of negative-semidefinite second differences if the existing stress fails the energy check. The product of two such second-difference operators can be positive-semidefinite: in the simple symmetric case, `L^2` has nonnegative eigenvalues when `L` has nonpositive eigenvalues.

Write the actual dissipative evolution operator. For example, with positive volume/mass matrix `M` and nonnegative weight matrix `K`,

`dv/dt = -M^{-1} B^T K B v`

gives `d(v^T M v/2)/dt = -(Bv)^T K (Bv) <= 0` when M is held fixed for that check and boundary terms are handled consistently. Derive the associated conservative momentum and energy fluxes for the actual grid; do not insert this illustrative operator without that derivation.

The previous review's description of an undivided third difference at a face was correct. Its divergence gives the fourth-difference effect in the cell update. Calling the method a fourth-difference stress does not invalidate the variable-coefficient warning. This terminology does not need further debate; the discrete energy test does.

### 5.4 Do not use the old state as the unchanged-solution oracle

A corrected flux can have a different finite-resolution root. The certified 0.03 control should be re-solved under a changed operator when necessary, and the physical residual, conservation, and grid convergence evaluated. An old state failing the new residual tolerance is not automatically a regression in physics.

Retain the unchanged-state check for genuinely unchanged operator branches. Also retain the full mass, momentum, and energy certification; rev1 correctly withdraws the hydrostatic-only acceptance shortcut.

## 6. Other items: largely acceptable, with bounded clarifications

### L15: good localization design; avoid proof by finite repetition

The revised trajectory-based experiment is appropriate. The wording that five repeated evaluations exclude residual dependence even at that state should be softened to: no difference was observed in the compared outputs during those trials. Finite repetition does not prove absence of a schedule-sensitive race.

For restart of the first divergent step, include all mutable state actually read, including thread-local caches and active-set/scaling history. A fresh process and a restored in-process evaluation should be compared if these states cannot all be serialized reliably. Keep test instrumentation outside the numerical reduction order where possible.

The new attribution of the old 2.07 speed-up to 16 threads differs from the earlier plan's eight-thread wording. Rev1 labels it as a memory record, not a fresh measurement. Treat it as unverified until a specific benchmark record or a new scoped timing establishes thread count, hardware, workload, and baseline. It should not be an acceptance threshold for the fix.

### L23: metadata contract is now appropriate

Preserving imported certification separately from a new evaluation is the right design. Add a policy for reasons longer than 32 characters and conflicting metadata between the two state files. Avoid silently attaching a valid-looking imported claim to a mismatched pair. Unknown tokens should remain provenance or be explicitly rejected, not treated as newly verified certification.

### L7g: the data-first gate is appropriate, but the current validity claim is unproved

The code comment at `molecular_reaction_heat.f90:463–476` requires density above the critical density of every relevant level, then supports the claim only with a v=1 critical density. That does not establish the all-level condition. Rev1 should call this an unverified validity assertion pending the inventory, not a verified high-density regime for the whole cascade.

The translational share is prompt kinetic energy of the products. Immediate deposition into a single thermal gas also assumes thermalization faster than transport and relevant reactions; verify that condition where the approximation is used. This qualification does not restore the old practice of multiplying all of D0 by a vibrational fraction.

The prompt/internal energy partition, R6 molecular-product audit, collider inventory, and double-counting tests are otherwise well framed. If a reduced uncertain model is needed, derive parameter bounds from data or explicit energy limits rather than adjusting them to recover a desired base temperature.

### L12: the corrected decomposition is a sound starting point

The additional stage flux `-n_el K grad(x)` avoids counting the element eddy flux twice, provided the same convention and face interpolation are used consistently. The discrete sum-of-stage-fluxes identity should include all neutral and ionized stages, not only the transported subset.

Make the flux boundary conditions explicit: inflowing stage fractions, outgoing characteristics, zero/reversed element flux, and matching to an inner equilibrium approximation. Include neutral-background limits, nonuniform element abundance, and helium fractions near simplex faces. No new physical acceptance shortcut is needed; the full coupled residual remains the criterion.

### L9: protect in-memory state, not just the input files

The explicit post-processing call is necessary and the new contract is much better. However, testing that input files are unchanged is insufficient: the existing loader already reads them without needing to overwrite them. Compare the loaded in-memory conserved state and composition before and after evaluation, including relevant module state that the chemistry and post-process use.

The output must say whether a product describes the original loaded state, a chemistry-refreshed work state, or an advection-derived composition. Decide which of these the post-process consumes. A certification of one must not be attached ambiguously to another. An empty product directory or stale-file sentinel should still test that the required `_adv` producer actually ran.

### L24 and records

Keeping the difficult fixture and adding a named large-pseudo-time variant is appropriate. Do not let a successful variant erase the original failure characterization. The revised naming and path-only verification policy is also appropriate.

## 7. Recommended decisions and revised gates

| Work | Decision | Gate before implementation or acceptance |
|---|---|---|
| L23 | Proceed | Metadata length, conflicting-pair, and re-evaluation semantics |
| L15 | Proceed with diagnostics | Compare actual perturbed states and all relevant mutable state |
| L26 probes | Proceed first | Include multidirectional zero-window limits and transients |
| L26 proposed repair | Do not adopt as written | A unique continuous limit and physically justified local boundary behavior |
| L22 experiments | Proceed after clarifying masks | Nonempty diagnostic sets; identify the metric already implemented |
| L22 solver change | Conditional | Define exit at an active cumulative bound and full-residual acceptance |
| L25 | Proceed as controlled experiments | HLLC/Roe/corrected-Roe comparison; correctly signed damping; defined gate |
| L7g inventory | Proceed | Do not presume the all-level high-density limit |
| L12 derivation | Proceed | Discrete flux identity and explicit boundary/matching conditions |
| L9 | Proceed after output-state choice | Verify memory-state preservation and actual downstream products |

The revision is now a useful basis for work. It does not need another broad rewrite. It needs a targeted correction of L26, clearer experimental definitions for L22/L25, and accurate descriptions of what the current code already implements.

## 8. Validation record

- Executed: selected executable/source hashes; the retained analytical boundary-limit script, including its assertions.
- Inspected: implementing call paths and current formulas cited above; revision and prior review.
- Not executed: EXHALE boundary tests, nonlinear solves, molecular/atomic campaign runs, parallel experiments, or regression suites. This turn reviews a plan and does not change the production model, so those runs would not validate an implemented correction.
- New files: this review and `docs/audit_20260905/plan_20260916_rev1_boundary_limit.py`.
- The earlier review is retained as history. Section 2 of this document explicitly corrects its carrier-progress claim; it should not be used as evidence that the default carrier path still measures only bounded displacement.
