# Review of the September 9 EXHALE issues and current implementation

Date: 2026-09-09, KST

Scope: `ISSUES_20260909.md`, the current source paths that implement its open issues, and downstream consumers of their results. This review proposes solutions and corrections; it does not modify production source or regenerate planetary products.

## 1. Executive assessment

The main development priority should be to establish a well-defined, constrained stationary residual and a reliable numerical method for reducing it. Increasing the marching budget, relaxing certification tolerances, or accepting a constrained least-squares minimum as a physical stationary state would not address the current problems.

The issue inventory is useful, especially its withdrawal of the earlier molecular convergence claims after discovering residual-history dependence. However, several explanations are too categorical, some entries lag behind the current source, and several additional defects need attention.

The most important findings are:

1. The stationary solver can initialize its trust-region radius to zero and then repeatedly produce a zero step even when a subsequent direction is available.
2. Its approximate Cauchy direction uses a step-length formula valid for an exact gradient. The current formula can increase the quadratic model even when the direction has a negative initial slope.
3. The production GMRES routine lacks a safe singular-Arnoldi termination. An isolated execution of the current routine returned a nonfinite solution while reporting a relative residual of zero for a zero operator.
4. Independent carrier ceilings are not the complete elemental feasible set. Shared hydrogen and oxygen budgets impose coupled constraints, and the current box also freezes bounds that depend on density and composition.
5. The `_adv` enthalpy-term test is a useful warning but is not a certification of stationarity. It accepts nonzero mass divergence whenever its ratio is below one.
6. The molecular branch of the post-process energy residual evaluates the upstream energy using the current cell's H2 fraction. This is not the energy difference of a compositionally varying mixture.
7. The transit warning counts only shells with `r <= Rib`, although its ray integrals also use shells with `r > Rib`. The diagnostic can underreport the atmospheric material carrying refused corrections.
8. The restart claim about exact comparison of printed He/H digits is not supported by the current implementation. The actual threshold is a relative difference of `1e-6`; the reported pair differs by approximately `8.83e-5`.
9. The physical integration error monitor compares hydrodynamic variables only. It does not directly control the integration error of transported composition, and its default sampling interval is 20 accepted steps.

Do not change the user's decision that a carrier bound does not excuse a nonzero physical balance row. A bound-constrained optimizer can stop at a least-squares minimum with a nonzero residual; that is numerical stagnation or an incompatible closure, not a stationary atmosphere.

## 2. Evidence, scope, and source identity

### 2.1 Evidence categories

- **SOURCE:** verified by reading the current implementation and its caller or consumer.
- **MEASURED:** executed during this review.
- **REPORTED:** a value in `ISSUES_20260909.md` or the update log, not independently reproduced here.
- **INFERENCE:** a proposed explanation or consequence that needs a targeted production experiment.

The issue document's statement that its numbers are measured does not make them new measurements in this review. In particular, the quoted planetary residuals, iteration counts, execution times, and spectral changes remain REPORTED.

No full regression suite, long planetary integration, or new stationary planetary solution was run. The affected numerical routines were exercised in an isolated harness, and the surrounding production paths were inspected. This scope avoids modifying the working executable or outputs while other source work is in progress.

### 2.2 Source and executable are different evidence

The repository HEAD is `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`, with extensive existing modifications. HEAD alone does not identify the reviewed source.

At inspection, `EXHALE.x` was dated September 8 at 19:03 KST and `build/steady_newton.o` at 15:40 KST. The current `steady_newton.f90` was dated September 9 at 00:01 KST. Therefore, executing the existing binary would not establish the behavior of every current source addition.

The Makefile includes `steady_newton.f90` and links LAPACK. Symbol inspection of the existing executable found `solve_steady_jfnk`, `pgmres`, `dogleg_step`, and `enthalpy_flux_term_ratio`, together with LAPACK references. This verifies linked entry points, not source equivalence. The isolated harness compiled selected routines directly from current source.

SHA-256 identities recorded during this review:

| File | SHA-256 |
| --- | --- |
| `docs/ISSUES_20260909.md` | `8d4768bb75ec4ca559c1004722241bff40597380d9d5f5c4bdb24ff22fb95683` |
| `src/modules/time_step/steady_newton.f90` | `0678dca1af85b9de93f56662e12348c382554f961837501f52a54dd1e0984e28` |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | `04c2e7e74254897ce3d0deb2ba6a2d52e7b4d8416f6652dd9ab11ed94cbdc8bf` |
| `src/EXHALE_main.f90` | `71d20728d256854d95ef3cbe486193a743fde56b8b6355eda5513cbf5443809a` |
| `src/modules/files_IO/load_IC.f90` | `de914886a055f5535d3957f3076207b8031a419d7270c865be238408a1fa9e87` |
| `src/modules/post_process/post_process_adv.f90` | `f328f158c4c6da0dfa8a5c3cfd9c18d7462142cb11368fdd6e7cef75404ccdc0` |
| `src/modules/nonlinear_system_solver/T_equation.f90` | `43b9c3c1cf43f7bf485bebb33bbef7ad9b577d0295ff2bb854801dd3f0dd7802` |
| `EXHALE.x` | `e1c191a0b0d95b1ac7c1a20513523373b1c93892a7024615f39f16c1ad0bab0b` |

Line numbers below refer to this inspected source. Routine names are the more durable locator.

### 2.3 Literature handling

The local `../references/` directory was inspected before literature lookup. The local Huang et al. paper is the final ApJ version; its methods, including Section 2.4.1, were inspected with `pdftotext -layout`. Its publication record was checked through NASA ADS: `2023ApJ...951..123H`, DOI `10.3847/1538-4357/accd5e`.

The file named `Schulik_Owen_2025_MNRAS_542_927.pdf` identifies itself internally as an arXiv preprint with placeholder journal metadata. It was not used as final-publication evidence. A journal-like local filename does not establish the version of its contents.

For numerical method recommendations, the official [SUNDIALS KINSOL mathematical documentation](https://sundials.readthedocs.io/en/latest/kinsol/Mathematics_link.html) was checked. It distinguishes variable scaling, residual scaling, inexact linear solutions, globalization, and successful residual convergence from small-step stagnation. These distinctions are relevant regardless of whether EXHALE continues using its own solver.

## 3. Corrections to the issue inventory

| Issue statement | Current assessment | Required correction |
| --- | --- | --- |
| B5l's carrier upper face is pending | SOURCE: `freeze_species_unknown_box` is present and called from the JFNK iteration | Describe implementation and validation separately. No new converged solution is established by routine existence. |
| GMRES size is a literal 40 | SOURCE: the caller still passes 40, but `solve_steady_jfnk` reads `EXHALE_GM_M` and reports the achieved linear residual | Keep 40 as the default description; do not state that changing it requires editing source. |
| Headroom availability has only a `huge` sentinel | SOURCE: `carrier_headroom_known()` also exists | Report the remaining magnitude-based checks, not absence of an availability query. |
| Molecular restart refuses exact equality of printed He/H digits | SOURCE: `heh_dev > 1e-6` is the branch condition | Identify the actual abundance discrepancy and its origin before changing the loader. |
| B5f supplies a stationary state accurate to `1e-8` for the reload-kick argument | REPORTED: the same record withdraws B5f as evidence of convergence because its residual was history-dependent | Call it a previously reported low-residual state; re-certify it with the corrected evaluator before using it as a stationary preservation test. |
| Seed retention proves an eigenvalue exactly equal to one | Not established by a finite seed-retention measurement alone | Distinguish slow physical chemistry, slow iteration, truncation, and multiple roots. |
| Coupled-loop pass cost is essentially fixed | Supported only for the tested fixed-point iteration and modifications | State that scope. A different nonlinear solve can change the iteration matrix and its contraction rate. |
| Every accepted `_adv` correction is exact | SOURCE: the equation has improved, but the validity screen is heuristic, the discretization is finite resolution, and the molecular closure is restricted | Say which equation and approximation were solved, to what residual, and under what validity conditions. |

The latest user decisions remain the starting policy. In particular, the amnesty cap is already selected, carrier logarithms are already selected, and no bound-based certification exemption is allowed. This review does not treat those decisions as unanswered questions.

## 4. Confirmed numerical defects and correction proposals

### R1. A zero initial trust-region radius is an absorbing failure state

**Priority: immediate. Evidence: SOURCE and MEASURED isolated behavior.**

Location: `solve_steady_jfnk`, approximately lines 6075-6081; `dogleg_step`, lines 4931-4969; `trust_region_step`, approximately lines 5289-5297.

The initial radius is constructed as:

```fortran
call pgmres(..., dZ, ...)
tr_delta = 1.0d-2*sqrt(sum(dZ*dZ))
tr_dmax  = 1.0d2*tr_delta
```

There is no requirement that this initialization call produced a nonzero, finite direction. `pgmres` legitimately returns a zero step when it cannot sample its first direction. The initialization condition is `tr_delta < 0`, so a zero result prevents a later retry of initialization.

With `delta = 0`, the current `dogleg_step` returns a zero step even when supplied with a nonzero Cauchy point and a nonzero Newton point. The outer method then reports a zero dogleg step without recovering the radius. The isolated test reproduced this behavior using the production routine.

This is a concrete mechanism consistent with the reported B5m radius-zero symptom. It does **not** yet establish why the first B5m initialization call produced zero; that needs its own trace. Counts of admitted residual evaluations do not answer whether the first scaled direction was usable.

Recommended correction:

1. Return a named linear-solve status, finite-step flag, achieved residual, and step norm from the initialization call.
2. Initialize from a finite, feasible descent direction. If the Krylov direction is unavailable, compute a valid projected gradient direction before choosing the radius.
3. Keep a positive, dimensionless lower initialization bound tied to meaningful variable scales. Do not derive both the initial radius and every future maximum radius from a failed solve.
4. If no direction exists, stop with the actual linear/model failure. Repeating zero-radius iterations is not recovery.
5. Let an adequate model grow the radius beyond an arbitrary ceiling inherited from its first step, subject to physically meaningful feasibility constraints and model quality.

Validation: zero initial Krylov direction followed by an available direction; nonfinite initialization; an initially tiny valid direction; a genuinely stationary zero-residual state; and one short atomic element-row restart after the isolated tests pass.

### R2. The approximate Cauchy point has the wrong line-minimizing length

**Priority: immediate. Evidence: SOURCE and MEASURED algebraic counterexample.**

Location: `trust_region_step`, approximately lines 5158 and 5220-5265.

Let the scaled residual be `r`, the full scaled Jacobian be `A`, and the approximate band matrix be `B`. The implementation forms:

```text
g = B^T r
Ag = A g
tau_current = (g^T g)/(Ag^T Ag)
sU = -tau_current g
```

For the model `m(s) = 0.5 ||r + A s||^2`, the minimizer along the actual direction `-g` is:

```text
tau_line = (r^T A g)/(Ag^T Ag), provided r^T A g > 0.
```

The numerators are equal only when `g` is the exact gradient, or in a special coincidental case. Checking `r^T A g > 0` rejects an initially ascending direction but does not make the old length valid. With the old length, even model decrease requires `2 r^T A g > g^T g`.

The isolated example `r = 1`, `A = 1`, `B = 10` gives:

| Quantity | Measured value |
| --- | ---: |
| Initial model merit | 0.5 |
| Merit at the current Cauchy point | 40.5 |
| Merit at the line-minimizing point | 0 |

A smaller trust radius can sometimes hide the overshoot, and the actual-reduction test can reject it. Neither repairs the claimed Cauchy decrease property.

The source comment defending the old expression refers to B5e's 66-iteration convergence. The issue record subsequently invalidates that convergence as history-dependent. It is not a sound reason to retain a mathematically incorrect line minimizer.

Recommended correction:

- Use the actual line minimizer when the band-derived direction is descending.
- Treat it as an approximate descent direction, not an exact-gradient Cauchy point.
- When it is not descending, obtain a direction with a verified decrease guarantee or report model failure. Dropping the direction does not retain the classical trust-region guarantee by itself.
- Verify the achieved residual of the Krylov leg before assuming that every truncated Newton leg reduces the model.
- Recompute predicted decrease for the actual projected step, retaining the correction already present in the source.

Validation: exact and approximate gradients, positive and negative slopes, a truncated Krylov leg, active constraints, and a direction whose unconstrained line minimum lies beyond the trust radius.

### R3. GMRES can report zero residual while returning a nonfinite solution

**Priority: immediate. Evidence: MEASURED execution of the current production routine.**

Location: `pgmres`, lines 4751-4913, especially Givens rotations and back substitution.

The routine replaces a zero Givens denominator with one, but then both rotation coefficients can be zero. In that case `gg(j+1)` becomes zero and is interpreted as convergence. Back substitution still divides by the zero diagonal entry of the Hessenberg matrix.

With an identity preconditioner, `A = 0`, and `b = 1`, the review harness measured:

```text
zero_operator_gmres_solution_finite=F
zero_operator_reported_relative_residual=0.0000000E+00
```

An identity-operator control returned the correct solution, one. The test uses production `pgmres`, with a deliberately simple Jacobian action and preconditioner. It proves the numerical failure mode, not that a planetary case has already reached this exact operator.

Recommended correction:

1. Distinguish successful Arnoldi termination from singular or inconsistent breakdown. A vanishing next basis vector is not sufficient evidence that the linear system is solved.
2. Check the reduced-system rank before triangular division; return a named failure or solve the rank-deficient least-squares problem safely.
3. Check LAPACK status and finiteness after every preconditioner application and after construction of the returned step.
4. Verify the true scaled linear residual when the reduced system is nearly singular, Arnoldi orthogonality is suspect, or convergence is claimed at breakdown.
5. Use stable norms and rotations, and add reorthogonalization when its measured loss warrants it.

The current preconditioner is fixed within a cycle, so ordinary right-preconditioned GMRES is conceptually appropriate. Flexible GMRES is needed only if that preconditioner is changed within a cycle. Replacing the solver library is not a prerequisite for correcting this defect.

## 5. Species constraints and the stationary formulation

### R4. A carrier box is a useful local model, not the full physical feasible set

**Priority: high. Evidence: SOURCE; production impact requires measurement.**

Relevant paths:

- `steady_newton.f90`: `freeze_species_unknown_box`, `fix_active_species_bounds`, `eval_residual`, and `write_species_rows_into_composition`.
- `diffusive_photochemistry.f90`: `carrier_steady_residual`, `carrier_element_headroom`, `hydrogen_available_to_carriers`, and `carrier_source`.

The current implementation has progressed beyond the issue document: upper faces are constructed, trial vectors are projected, and the model action is recomputed after projection. These are useful changes. However, the physical feasible set is not generally an independent box.

For example, with normalized available hydrogen equal to one and both H2 and H+ transported, the separate conditions are:

```text
n(H2) <= 0.5
n(H+) <= 1
```

The state at both individual ceilings contains two hydrogen nuclei in these normalized units. The actual shared constraint is:

```text
2 n(H2) + n(H+) + n(OH) + 2 n(H2O) <= available hydrogen.
```

Likewise, OH, H2O, and CO share an oxygen constraint; CO also consumes carbon. Helium-hydride and the eliminated molecular ions must be included consistently in the elemental inventory, according to which species remain fixed in the closure.

`carrier_source` already exposes the shared sums: it constructs atomic H and O from the remainder and clamps that remainder to zero. This means that individual box feasibility cannot by itself establish that the physical composition is valid. Other chemical or final certification checks may reject an excessive combined inventory; this review does not claim they all fail. The defect is the incomplete feasible-set model used to choose and sample directions.

There is a second issue. The current box freezes the density ceiling and headroom at the outer iterate. Density and elemental fractions are also Newton unknowns, so physically admissible simultaneous changes can require a moving ceiling. A frozen ceiling can exclude an inward or outward direction of the true coupled constraint surface.

A third concern is explicit in `eval_residual`: headroom admission compares the number of violating cells with the iterate's count. Equal counts do not bound violation magnitudes or guarantee that the same cells violate. `freeze_species_unknown_box` also raises a frozen headroom face to the iterate when necessary. These can be legitimate continuation devices, but neither is a physical feasibility certificate.

Recommended formulation:

1. Write one elemental constraint map from the species stoichiometry and the current density/element unknowns.
2. Evaluate physical feasibility from that map at every candidate state.
3. Construct the local linearized constraint set from its derivatives, including density and elemental-fraction derivatives. For logarithmic carrier coordinates, include `dn/d(log n) = n`.
4. Use a constrained step that can move along shared element-budget surfaces. Coordinate clipping alone does not represent these directions.
5. Keep an explicit restoration phase for an infeasible initial state, measuring the magnitude of constraint violation rather than only a count.
6. Require the original physical balance rows and elemental constraints at certification. Do not replace either with optimizer stationarity conditions.

An alternative parameterization can enforce element conservation by construction. It must preserve independent degrees of freedom, handle trace and zero-abundance species, and include its Jacobian consistently. A broad change to unknown layout should be presented for approval before implementation.

For the next B5l measurement, compare the old and new methods at the same saved state using constraint violations, actual merit reduction, achieved linear residual, and elapsed time. Do not compare only Newton iteration counts or whether a component touches a bound.

### 5.1 Shielded H2: retain transport, but diagnose the source of seed dependence correctly

The refusal of the currently unsafe eliminated-H2 coupled configuration is justified. A Newton residual must be reproducible from its declared unknowns, not from a previously visited chemical state.

However, a finite sweep retaining part of a seed perturbation does not prove that local equilibrium is mathematically nonunique or that a fixed-point eigenvalue is exactly one. Possible explanations are:

- a nearly conserved chemical mode with a long physical relaxation time;
- a poorly contracting iteration for a unique equilibrium;
- several self-consistent radiative/chemical roots;
- a stopping metric insensitive to an important residual consumer;
- incomplete restoration of an intermediate state or changed radiation columns.

The physically relevant comparison is the chemical relaxation time against transport times. Obtain chemical decay rates from the reaction Jacobian after removing elemental conservation modes. Compare them with an advective time based on the resolved gradient length and a diffusion time based on that same length. A cell width alone is not automatically the physical gradient scale.

Transported H2 is appropriate when the slow molecular mode is controlled by flow, diffusion, or boundary supply. The fast species can remain eliminated if their conditional root is sufficiently unique, stable, and accurate. There is no need to promote every fast intermediate into the global unknown vector by default.

Recommended residual contract:

```text
declared state Y
  -> current elemental totals and transported carriers
  -> conditional fast chemistry and radiation closure
  -> hydrodynamic, elemental, and carrier balance rows
```

Repeated evaluations with the same `Y` and different admissible fast-chemistry seeds must agree below a measured error budget of the outer residual. If not, either improve that closure or promote the unresolved slow mode into `Y`.

### 5.2 Two ionization-front roots require a branch study, not a residual-only choice

The reported different front locations with small local chemical residuals are important, but the evidence does not yet identify physical bistability.

At identical hydrodynamic unknowns, changed chemistry also changes optical depths and radiation. Thus, identical `Y` does not necessarily mean identical local rate coefficients. The study must distinguish a local reaction system at fixed radiation from the full self-consistent radiative/chemical system.

Recommended sequence:

1. Freeze temperature, elemental totals, and the complete radiation field. Solve the same local reaction system from multiple seeds and verify positivity, charge, and elemental conservation.
2. Release the radiation closure and repeat at fixed hydrodynamics. Record incident boundary data, columns, rates, and cell-averaging choices for both solutions.
3. Repeat on successively refined front grids, including the attenuation integral within the front cell. Compare front position and integrated columns, not just a cell index.
4. Continue irradiation or boundary abundance gradually from a simple solution. Reverse the continuation to test hysteresis.
5. Where physically admissible time integration is available, perturb each branch and measure its stability. Include the relevant chemical and transport time scales.

If two stable branches remain, the initial state and boundary history become physical model inputs. If one branch disappears with closure tolerance or grid refinement, the original observation was a numerical branch ambiguity. Do not choose the branch merely because it gives a preferred mass-loss rate or spectrum.

### 5.3 Atomic element-row stagnation: the next diagnostic should be short and ordered

The atomic case is valuable because it separates the element-row problem from H2 transport. The reported 11-variable layout is consistent with the registry: three hydrodynamic quantities plus eight element rows. Those element rows are helium and the seven trace elements C, O, N, Mg, Ca, Na, and Fe.

The current source already gives a reservoir equation to the first element cell and an identity-like absent-element row. Do not restart the investigation with an unsupported assertion that the base row is still empty.

Use the existing saved atomic state and inspect the first two Newton iterations:

- minimum/maximum column and row scales, with the associated physical quantity and cell;
- zero, nonfinite, or effectively unsampled Jacobian columns;
- band-factorization status, pivot magnitudes, and a condition estimate;
- initialization GMRES status, returned norm, and achieved residual;
- counts of bound-limited directions separately from residual-evaluation failures;
- projected gradient norm, `r^T A g`, predicted decrease, and trust radius before any trial.

First correct R1-R3. Then study whether the preconditioner needs more local chemistry coupling or a treatment of nonlocal radiation. A larger GMRES space cannot repair a zero radius or an invalid Cauchy point.

Use the existing `EXHALE_GM_M` control to compare 40, 80, and, if warranted, a larger space at identical states. Select based on achieved linear residual and total residual-evaluation cost, not a requirement that every cycle use fewer iterations.

### 5.4 Base amplification and certification scales

The reported base energy amplification is plausible: `eval_residual` updates `n_part_cell1` before applying the boundary condition precisely because the ghost temperature depends on the current composition. It is a coupled boundary sensitivity, not grounds for exempting the first physical cell from energy balance.

Recommended checks:

1. Determine whether each first-cell equation is a physical control-volume balance or an imposed reservoir condition. Use the corresponding equation consistently in the residual and Jacobian.
2. Differentiate through composition, particle count, caloric EOS, and the boundary reconstruction together.
3. Resolve cancellation between energy flux, pressure work, radiation, and boundary terms using dimensional values as well as normalized residuals.
4. Measure directional-derivative error against perturbation size and inner-closure tolerance at the base.
5. Use the measured map from closure error to outer residual to set the inner tolerance. A universal `eq_sweep_reltol * merit` estimate is not a rigorous bound when the boundary amplifies composition error.

Keep nonlinear row scaling separate from certification normalization. Improving a row scale to condition the linear algebra is allowed; changing the physical acceptance tolerance to hide a troublesome row is not.

The `1e-8` element and carrier tolerances are currently provisional. Anchor them through repeated residual evaluation, derivative convergence, a known steady manufactured solution, grid refinement, and error in conserved fluxes and observables. Different precision limits in wind and diffusion-dominated regimes should be reported explicitly rather than averaged into one permissive threshold.

## 6. Restart and stationary-state preservation

### R5. The reported He/H failure is not an exact-digit comparison

**Evidence: SOURCE and MEASURED arithmetic; original failing file not replayed.**

`load_IC.f90`, lines 273-315, reconstructs hydrogen and helium nuclei including molecular contributions, checks physical cells rather than ghosts, and forms:

```text
heh_dev = max_j |(nHe/nH)_j - HeH_input| / max(HeH_input, 1e-30).
```

The branch is entered when `heh_dev > 1e-6`. A nonzero HeH+ population then prevents a simple uniform rescaling when helium diffusion is disabled.

The quoted values `7.9307E-02` and `0.0793` differ relatively by `8.8272383e-5`, approximately 88 times the threshold. Moreover, `write_output` uses list-directed double-precision output for the species columns, not the short diagnostic format used in the error message.

Consequently, changing equality to an approximate comparison is not a fix: an approximate comparison already exists. Nor should the tolerance simply be raised to accommodate this pair.

Recommended investigation and correction:

1. Identify the exact input and restart pair, their source and constant versions, and whether the atmosphere was allowed to develop a varying elemental composition.
2. Perform a writer-to-loader round trip without evolution or chemistry updates, comparing density, elemental totals, charge, molecular abundances, and thermal energy.
3. If the writer and loader reconstruct different inventories, fix that inconsistency.
4. If the snapshot itself violates its stated fixed elemental abundance, correct the path that created the state rather than silently changing it on restart.
5. Store the input reservoir abundance, species schema, physical radius grid, constants, model options, and source identity in restart metadata at sufficient precision.
6. Separate an exact restart from an intentional composition-changing initialization. The latter needs an explicit conservative mapping and is not an unchanged restart.

For changed Jupiter radius conventions, retaining the grid guard is correct. A changed physical radius is not merely a renamed coordinate. A conservative remap may be possible, but it must preserve the intended physical boundary and integrated inventories, and it requires a separate authorized workflow. Do not regenerate all stored products as part of this review.

### 6.1 Do not march a stationary restart merely to reach Newton

SOURCE: the ordinary Newton finish is invoked from inside the time loop, after an attempted marching step. A direct steady entry also exists through the `EXHALE_PTC` route and can select JFNK. Therefore, the issue is not absence of a callable stationary solver; it is routing and initialization semantics.

A normal restart workflow should distinguish:

- continue a physical trajectory;
- continue relaxation toward stationarity;
- evaluate or improve a stationary state directly.

For the third case, load, rebuild required derived quantities without changing conserved state, evaluate the actual stationary residual, and either return the certified state or enter Newton immediately. Do not unconditionally take CFL steps or re-equilibrate a transported quantity first.

The claimed `1e-8 -> 1.226` jump must be remeasured from a state certified by the corrected, seed-independent evaluator. The old B5f state is not currently that reference.

Separately, splitting can disturb a steady balance because individual transport and source updates are nonzero even when their sum vanishes. Whether the full-step displacement starts at order `dt` or `dt^2` depends on the actual composition reset and splitting formula. Measure the preserved-state defect with `dt`, `dt/2`, and `dt/4`; do not infer its order solely from one two-step restart.

## 7. Advection post-processing and transit consumers

### R6. The enthalpy ratio is not a stationarity certificate

**Priority: high for physical interpretation. Evidence: SOURCE and MEASURED isolated function.**

For an internal-energy equation without additional transport terms, the steady conservative form can be expanded as:

```text
div(rho v e) + p div(v) = Q
rho v de/dr - (p/rho) v d(rho)/dr + h div(rho v) = Q
h = e + p/rho.
```

Restoring the enthalpy term is algebraically correct. However, without a mass source, a physically stationary atmosphere must also satisfy `div(rho v) = 0`. Solving an energy equation on a prescribed nonstationary mass profile does not make the complete atmosphere stationary.

The current `enthalpy_flux_term_ratio` compares the omitted-mass-balance contribution with two energy terms. A ratio below one can still represent a large inconsistency. The production function returns 0.5 for the harness example with a nonzero mass divergence, and the present threshold accepts it. This is a counterexample to the interpretation, not a claim about the magnitude in a particular planet.

Also, the post-process uses center-to-center `rho v r^2` differences. The production hydrodynamic balance uses numerical face fluxes. The continuum form is related, but these are not identical discrete stationarity tests. A finite-volume steady solution can have a small numerical mass residual without exactly constant center-sampled `rho v r^2`.

Recommended correction:

1. Use the actual stationary certification and the same numerical mass-flux operator for any claim that the input is stationary.
2. Retain the enthalpy ratio as a diagnostic of sensitivity and energy-term size, not as the sole physical validity test.
3. If approximate local use is allowed, state a tolerance tied to the required error and describe the output as a conditional correction, not an exact steady state.
4. Test momentum, composition, and energy compatibility where relevant; continuity alone is necessary, not sufficient.
5. Record separate temperature and composition validity flags. One integer representing the first refusal cannot fully describe every partially corrected row.

A profile assembled from corrected and retained cells also need not satisfy interface balances. This limitation must follow the profile into the observable calculation.

### R7. Upstream molecular caloric energy is evaluated at the wrong composition

**Priority: high for molecular `_adv` use. Evidence: SOURCE.**

Locations: `T_equation.f90`, approximately lines 255-260; `post_process_adv.f90`, approximately lines 647-650 and the energy-solve caller.

The molecular energy difference uses:

```text
E(x_H2,j, T_j) / mu_j - E(x_H2,j, T_{j-1}) / mu_{j-1}.
```

The second term should use the upstream H2 particle fraction, not the current cell's fraction. The issue is not repaired by using the upstream mean molecular weight: molecular rotational and vibrational heat capacity also depend on molecular abundance.

Across a varying H2 front, this omits part of the advected caloric-energy change. The atomic branch can likewise miss the upstream molecular energy when the current cell has zero H2 and the upstream cell does not.

Recommended correction: pass the upstream specific internal energy, or the complete upstream caloric state, into the residual. Prefer evaluating the conservative energy transport through the same thermodynamic definitions used by the main solver. Test an isothermal composition gradient and a stationary flow across a molecular dissociation layer.

The writer already warns that molecular and oxygen columns are not advection-corrected. That restriction should remain. Correcting this energy difference alone will not turn the post-process into a fully coupled molecular transport solve; elemental, particle-count, charge, chemical-energy, and species-flux consistency must all be established first. Until then, molecular `_adv` products should not be presented as fully physical molecular stationary solutions.

### R8. Transit refusal counts omit contributing outer shells

**Priority: high for output reliability. Evidence: SOURCE and MEASURED geometry example.**

Location: `EXHALE_transit.py`, lines 399-423; compare the ray construction around lines 578-605 and the resonance-line integral in `exhale_transit_lib.py` around lines 483-490.

The warning selects:

```python
_in_range = (r >= 1.0) & (r <= Rib)
```

Its comment argues that the impact-parameter interval `1 <= b <= Rib` excludes atmospheric shells outside the same radial interval. This is geometrically incorrect.

The actual ray integral includes every available shell satisfying `r >= b`, with line-of-sight coordinate `x = sqrt(r^2-b^2)`. For example, when `Rib = 10`, the shell at `r = 12` intersects the ray at `b = 2`, at `|x| = sqrt(140) = 11.8321596`. The harness measured this geometry. That shell can contribute absorption despite being omitted from the refusal count.

Recommended correction:

- Construct the diagnostic from the same shell/ray selection used by the optical-depth integration.
- Count all contributing physical shells; do not impose `r <= Rib` on a radial-shell census.
- Where possible, report the optical-depth contribution of retained or refused material, resolved by line and wavelength. Counts alone do not measure spectral uncertainty.
- Carry the validity information into saved spectral metadata, not only console output.

This finding concerns the refusal diagnostic. The inspected ray selection already includes outer shells, so this review does not claim that the corresponding optical-depth integral incorrectly truncates them at `Rib`.

## 8. Coupled-loop cost and physical time integration

### 8.1 COST5 identifies one iteration's contraction, not an algorithmic lower bound

For a fixed-point iteration, a relation between pass count, logarithmic tolerance reduction, and contraction rate is expected. It explains the measured pass count for that map. It does not prove that the equations intrinsically require those passes.

For example, a scalar fixed point `x = theta*x + c` can require many Picard iterations when `theta` is close to one, while Newton solves the corresponding linear residual in one exact step. A coupled chemistry/temperature Newton solve changes the iteration map; so can acceleration or a better block elimination.

The reported unsuccessful ordering changes and extrapolation experiments are useful negative results. Keep extrapolation off unless the actual error and execution-time measurements justify it. Do not repeat those experiments without a new reason.

Recommended next experiment, after residual correctness:

1. Isolate the molecular source block at fixed density and transported inventories.
2. Solve temperature and conditional chemical composition together, including derivatives of particle count, electron density, third-body terms, and caloric energy.
3. State how columns and radiation participate. A local block Newton does not automatically eliminate global shielding modes.
4. Use the accepted-state error and energy closure as criteria, and measure residual evaluations and elapsed time.

### 8.2 COST7's factor-three estimate is empirical, not a universal error bound

SOURCE: `coupled_pair_geometric_error_estimate` uses successive increment vectors, a scalar mode fit, model-agreement guards, and a safety factor multiplying `|theta/(1-theta)|` times the latest increment.

For a genuine single contracting mode, that is the correct geometric-tail relation. The code's guards are useful. Nevertheless, multiple modes, changing norms, changing active chemistry, or a small unresolved component with strong downstream sensitivity can invalidate a scalar estimate.

Keep periodic overconverged references, including failed or unusually slow cases. Compare errors in temperature, important carriers, particle count, heating/cooling, and assembled residuals. Near an inner iteration floor, report an unattained tolerance rather than equating a small change with an accurate fixed point.

### R9. Physical integration error control does not directly cover transported species

**Priority: high before quantitative time-dependent molecular claims. Evidence: SOURCE.**

`attempted_step_error_estimate` in `attempted_step.f90`, lines 1254-1280, receives only the three hydrodynamic conservative arrays and compares full and two-half-step results in those variables. `attempted_step_configure` sets the physical-mode sampling interval to 20 by default. Thus, many accepted physical steps are not individually checked with step doubling, and no transported-species error appears directly in that estimate.

The carrier substep acceptance checks a solved time-discrete equation. That is not the same as controlling the integration error of the chemical trajectory. A backward-Euler chemical step can solve its nonlinear equations accurately while having a large temporal error.

Recommended correction:

1. Include transported carrier and elemental variables in the full-versus-half-step error comparison, with abundance-aware absolute and relative scales.
2. Include charge and thermodynamic consistency in the comparison where changes in trace species strongly affect opacity or chemistry despite a small hydrodynamic response.
3. For accuracy-critical physical trajectories, use an error controller on every accepted step, or provide a validated bound for an intermittent policy. A diagnostic check every 20 steps is not such a bound.
4. Test temporal convergence of a chemically evolving state whose hydrodynamic variables change little, as well as a coupled thermal/chemical case.

The hydrodynamic row called `hydro_row_max` also deserves precise wording: it compares the RK update with a single initial Euler operator evaluation, and is deliberately diagnostic rather than decisive in the current certification routine. Simply enforcing its current threshold would not create a correct RK residual test. Use the actual staged discrete equation or the complete step error estimate.

The present separation of initialization and physical modes is useful. Non-root classes 4 and 6 are refused in physical mode by the inspected certification path, while the selected amnesty policy remains a continuation allowance. Preserve that distinction; a green relaxation matrix is not evidence of a validated physical trajectory.

## 9. Recommended implementation order and acceptance tests

### Stage A: recover a reliable numerical step

Correct R1-R3 first. Add named GMRES outcomes, finite checks, a positive recoverable radius initialization, and the correct line minimum for the direction actually used.

Acceptance:

- isolated identity and singular linear tests behave correctly;
- no accepted or reported-converged linear step is nonfinite;
- an available descent direction is not suppressed by a zero radius;
- predicted decrease agrees with an explicitly evaluated small linear model;
- failure reasons distinguish linear breakdown, infeasibility, and nonlinear rejection.

### Stage B: establish the constrained residual

Audit and then correct shared elemental constraints, moving bounds, and residual state ownership. Use one small atomic element case and one small H2 transport case before large atmosphere integrations.

Acceptance:

- repeated `F(Y)` evaluations are independent of prior probes and seeds within a measured budget;
- linearity and finite-difference convergence of `Jv` are measured over a useful perturbation interval;
- the same stoichiometric constraints govern trial construction, chemistry closure, and certification;
- a state at a bound with a nonzero balance row remains uncertified;
- a feasible direction along a shared elemental constraint is not discarded merely because a coordinate bound was frozen.

### Stage C: obtain the first trustworthy species-row stationary solution

Start with the atomic element case to separate linear algebra from molecular closure. Then solve H2 transport, followed by oxygen carriers and the more complete chemistry.

Acceptance should include the complete row registry, dimensional and normalized residuals, elemental/charge closure, discrete mass flux, thermal balance, and reproducibility from independently constructed nearby initial states. Then perform grid and inner-tolerance refinement.

A successful single low-residual iteration history is not sufficient if restarting or re-evaluating the same state changes its residual materially.

### Stage D: repair restart and observable interpretation

Complete unchanged-state restart tests, direct Newton routing, R6-R8, and explicit validity metadata. Keep old products as historical results with their source and physics versions. Do not silently relabel them as current predictions.

Only after the user authorizes the affected planetary work should the stale SED/physics products be regenerated. A file-layout correction or warning fix by itself does not require regenerating physical profiles.

### Stage E: validate physical trajectories and optimize cost

Complete species-aware time-error control and temporal refinement before using physical integration to decide front stability or hysteresis. Then benchmark the coupled source improvements against those validated solutions.

This order prevents a faster implementation from accelerating an ambiguous residual, and prevents a new spectrum from being interpreted before its underlying state is physically certified.

## 10. Tests performed and remaining limitations

The following review code is retained under the requested audit directory:

- [Fortran harness prefix](audit_20260905/issues_review_prefix_20260909.f90)
- [Fortran test program](audit_20260905/issues_review_suffix_20260909.f90)
- [Execution script](audit_20260905/run_issues_review_20260909.sh)

Run from `EXHALE_v1.00`:

```bash
bash docs/audit_20260905/run_issues_review_20260909.sh
```

The script extracts the current production `pgmres`, `dogleg_step`, and `enthalpy_flux_term_ratio` into a compiler input stream. It supplies an isolated Jacobian action and identity preconditioner, compiles with GNU Fortran 16.2.0, `-O0 -fcheck=all`, and does not change production source, the normal executable, or planetary output files. Temporary compiled products are placed in a uniquely named `/tmp/exhale-issues-review-20260909.*` directory.

Measured results:

| Check | Result | Meaning |
| --- | --- | --- |
| Nonzero dogleg directions with radius zero | Step norm zero | Confirms the radius failure mechanism |
| Approximate-gradient Cauchy example | Merit 0.5 to 40.5 | Current length does not guarantee model decrease |
| Correct line minimum in the same example | Merit zero | Establishes the mathematical correction |
| Identity operator through production GMRES | Solution one | Control test |
| Zero operator through production GMRES | Nonfinite solution; reported relative residual zero | Confirms singular-Arnoldi handling defect |
| Nonzero mass divergence in production enthalpy test | Ratio 0.5 | Below-one ratio does not imply stationarity |
| Quoted He/H pair | Relative difference `8.8272383e-5` | Not an exact-equality problem |
| H2 and H+ at individual normalized ceilings | Hydrogen inventory two, available one | Independent boxes do not express the shared budget |
| Shell `r=12`, ray `b=2`, projected limit ten | Chord coordinate `11.8321596` | Outer shell contributes despite warning mask |

All harness assertions passed, meaning the stated failure modes and controls were reproduced. This does **not** mean the production defects were corrected.

Not performed:

- reproduction of B5m's original complete run or proof of the first zero direction's cause;
- a new converged atomic element or molecular carrier atmosphere;
- replay of the exact molecular restart file that triggered the He/H error;
- measurement of molecular `_adv` energy or transit-warning errors in a planet;
- full bounds/FPE regression of the current complete source tree;
- regeneration of any user-gated planet, benchmark, paper, or poster result;
- a fresh audit of every reaction rate and radiative correction listed as resolved.

The resolved radiative work should still be judged against physical closure, not only regression identity. For example, the published [Huang et al. methods](https://doi.org/10.3847/1538-4357/accd5e) explicitly combine outer- and inner-shell photoionization cross sections. That supports checking a shell-complete opacity treatment, but does not by itself validate EXHALE's selected effective charge-stage and heating approximation. The current user-selected approximation should retain an explicit applicability statement and its own charge/energy audit.

The immediate recommendation is therefore narrow and actionable: correct the demonstrated numerical defects, establish the coupled elemental feasible set and residual reproducibility, and obtain one independently certified species-row solution before changing tolerances or launching extensive production reruns.
