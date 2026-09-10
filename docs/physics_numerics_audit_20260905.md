# EXHALE v1.00: physical and numerical audit, including molecular atmospheres

Date: 2026-09-05 (KST). Inspected revision: `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`.

This report analyzes the current implementation and proposes a development sequence. It does not implement the proposed changes. Production source, input data, existing simulation products, and regression reference products were not modified.

## 1. Main assessment

EXHALE contains substantial molecular physics already: H2/H2+/H3+/HeH+ chemistry, an OH/H2O oxygen cycle, CO storage and transport, H2 rovibrational internal energy, chemical heating, ultraviolet photolysis, infrared emission, and transport of selected composition variables. The central problem is consistency among these components, rather than the absence of a molecular module.

The most consequential issues found in this audit are:

1. The marching composition update changes thermal energy at fixed temperature, separately from the explicit physical energy sources. This is a conservation defect in the update map and affects atomic as well as molecular composition changes.
2. H2 double photoionization spends 51.4 eV in the heating calculation, but the asymptotic chemical energy of its represented products is only about 31.675 eV. The difference, associated with fragment kinetic energy in the vertical-ionization picture, has no identified recipient in the inspected energy calculation.
3. The LLF signal speed omits the negative acoustic characteristic. Two additional errors occur in the optional ROE wave-speed calculation.
4. The transported-proton Jacobian and residual floor use the oxygen abundance as their reference scale. The H+ carrier was added without extending these two species-dependent branches.
5. H3+ cooling has a nonzero low-collider-density floor and discontinuous temperature joins. The discontinuities are present even when the published polynomial coefficients are transcribed correctly.
6. Chemistry and the caloric EOS use different H2 internal-state models. The resulting partition function and internal energy do not satisfy a common thermodynamic identity.
7. The oxygen cycle has chemistry and photolysis heating but lacks a complete collisional reaction-energy budget. CO could be removed by an equilibrium ceiling without a kinetic destruction rate; that half is resolved (P6, items B3b-CO and CEILING-DEL, 2026-09-06).
8. The steady solver cannot solve the transported-proton system, while its existing composition consistency check is weaker than convergence of the full chemical state.

These findings change the recommended priority in `code_status_20260905.md`. Refreshing regression reference products should not be the first scientific milestone. First establish conservation, correct the demonstrated implementation errors, and define a common set of equations for evolution and steady convergence. A close match to a published profile cannot compensate for an incorrect equation.

The atomic branch has valuable reported validation history, but it should not be described as generally physically validated across all regimes. The molecular branch has enough implemented physics to support a systematic verification program, but current fixed-step snapshots cannot establish a converged molecular escape solution.

## 2. Evidence, scope, and limitations

### 2.1 What was checked

The supplied status document was read and compared with the code. Inspection followed these paths:

| Area | Implementation inspected | Role in the calculation |
|---|---|---|
| Execution and configuration | `EXHALE_main.f90`, `input_read.f90`, `Makefile` | Invocation order, option restrictions, compiled routines |
| Hydrodynamics | `Num_Fluxes.f90`, `speed_estimate_ROE.f90`, `RK_rhs.f90`, `eval_dt.f90`, boundary and grid routines | Fluxes, positivity repair, geometry, time stepping |
| Thermodynamics | `caloric_eos.f90`, `UW_conversions.f90`, `composition.f90`, `mol_rates.f90` | Pressure, internal energy, heat capacity, partition functions |
| Chemistry and radiation | `ionization_equilibrium.f90`, `util_ion_eq.f90`, `System_HeH_mol.f90`, molecular and oxygen rate modules | Species sources, reaction heat, radiation accounting |
| Composition transport | `diffusive_photochemistry.f90`, relevant binary-element diffusion APIs | Transport variables, coefficients, constraints, source Jacobian |
| Steady convergence | `steady_newton.f90`, `steady_residual.f90` | Unknowns, residuals, acceptance conditions |
| Lower boundary and products | Lower-profile reader/schema, `lower_column.f90`, `post_process_adv.f90`, collision diagnostic | Handoff meaning, downstream composition, validity limits |

Inspection concentrated on the expressions and their callers that determine these conclusions. It was not a line-by-line audit of every atomic coefficient, every scientific-data row, every post-processing line profile, or the complete Wind-AE implementation.

The Makefile includes the molecular and flux routines discussed below. `nm EXHALE.x` also showed linked symbols for `photochemical_transport_step`, `h3p_cooling_rate`, `keq_H_H_to_H2`, `molecular_chemical_heating`, `oxygen_carrier_rows`, and `speed_estimate_ROE`. Thus these are executable components, not merely unused source files. The existence of a linked routine is distinguished here from validation of its physical approximation.

The working tree was clean at the start. The supplied document instead describes HEAD `6d07d48` and 175 outstanding files. That is a historical description, not the state inspected here. Its planet-specific residuals, temperatures, density ratios, and mass-loss rates were not remeasured in this audit.

### 2.2 Evidence labels

- **Confirmed implementation defect:** a wrong expression or inconsistent update is present in a reachable code path. A reduced execution is included where practical.
- **Confirmed model limitation:** the implemented equations omit a process or impose a restricted approximation; this is not necessarily a typographical mistake.
- **Verification risk:** the present checks do not establish the required property. An additional test is specified rather than an unmeasured failure asserted.
- **Historical claim:** a result reported in the supplied status document, not independently reproduced here.

Severity describes the affected property and scope. It does not claim an unmeasured change in a planet's mass-loss rate.

### 2.3 Executed diagnostics

The attached [driver](audit_20260905/audit_probe.f90) compiles and calls the current production routines. The [script](audit_20260905/run_checks.sh) builds in a new temporary directory and also executes the existing H2-channel and oxygen-kinetics diagnostic drivers.

```bash
cd EXHALE_v1.00
bash docs/audit_20260905/run_checks.sh
```

The compiler used was GNU Fortran 16.2.0, with `-O0 -g -fcheck=all -fbacktrace -fopenmp`. Final outputs for this audit were written to `/tmp/exhale-physics-checks-xNoZEQ/`. The script is retained so the evidence does not depend on that temporary directory surviving.

These are diagnostic programs. A zero exit status means execution completed; it does not mean every printed physical property passed. The driver intentionally exposes failures. It does not replace the full atmospheric calculation.

| Diagnostic | Directly measured result | What it establishes |
|---|---|---|
| LLF negative flow | Implemented speed 1; required acoustic bound 3 | Missing negative characteristic |
| ROE rarefaction, density and pressure multiplied by 10 | Left star sound speed 1.21221780 becomes 0.80090436 | Unphysical dependence on density normalization |
| ROE shock, same scaling | Left star sound speed 4.90216311 becomes 15.50200089 | Missing upstream-density factor |
| H3+ at 1000 K, n(H3+)=1 cm^-3 | Cooling is 5.45570188e-15 erg cm^-3 s^-1 for n(H2)=0, 1, and 1e6 cm^-3 | Low-density clamp prevents the collisional cooling limit |
| H3+ temperature joins | Relative jumps +0.4304733 at 300 K, -0.0011301 at 800 K, -0.0235334 at 1800 K | Discontinuous cooling function |
| H2 internal energy at 4000 K | Chemistry-derived rovibrational energy is 5.1400% below the EOS value | Thermodynamic models differ |
| H2 double-channel energy | 51.4 - [2 I(H) + D0(H2)] = 19.72505589 eV | Energy absent from the represented chemical products |
| Infrared thermal equilibrium at 700 K | Net/emission: H2 +5.04e-5; H2O -1.14e-5; CO -9.96e-5 | Separate table interpolation does not exactly preserve zero net exchange |
| Existing H2-channel diagnostic | Cross-section closure error <=4.441e-16; zero element/charge violations in four channels | Stoichiometry and channel sum work; energy closure is a separate issue |
| Existing oxygen-kinetics diagnostic | Largest printed OH scaled residual 6.120609e-9 at 300 K; <=6.509e-15 over 864-3000 K | Reduced equilibrium algebra is well resolved over the latter range |

The oxygen diagnostic reconstructs reduced balance expressions. It does not run the full atmospheric chemistry/transport integration. The first compile attempt lacked the water-photolysis dependency; after including its actual dependency chain, compilation and execution completed. No substitute physics module was used.

No full regression matrix, new planet convergence calculation, transit regeneration, or full-project rebuild was performed. The selected compilation checks the directly probed production routines; it does not certify unrelated source files. This scope is appropriate for an analysis-only request, while leaving the proposed physical acceptance tests explicit.

## 3. Confirmed numerical and conservation defects

### N1. Fixed-temperature composition projection changes energy without a physical source

**Priority: P0. Confirmed implementation defect.**

Evidence: [EXHALE_main.f90](../src/EXHALE_main.f90), lines 1162-1228; [composition.f90](../src/modules/functions/composition.f90); [UW_conversions.f90](../src/modules/functions/UW_conversions.f90).

After the hydrodynamic update, the marching loop obtains temperature from the current pressure and composition. It then transports carriers and solves ionization at that temperature. It calculates the new particle count, reconstructs pressure at the unchanged temperature, and calls `W_to_U` before applying heating and cooling.

For the monatomic limit,

\[
u_{\rm th}=\frac32 n_{\rm part}k_BT,
\qquad
\Delta u_{\rm projection}=\frac32 k_BT
       (n_{\rm part,new}-n_{\rm part,old}).
\]

This contribution is not obtained by integrating a physical heating term. The molecular EOS adds a change in the H2 rovibrational reservoir as well. Nuclear conservation does not remove either energy change.

A simple example is a hydrogen gas at fixed hydrogen-nucleus density. Increasing its ionized fraction from zero to 0.5 changes the particle count from n_H to 1.5 n_H. Keeping T fixed increases thermal energy by 50%, before the specified photoheating is added. The attached driver uses the production conversion routine to demonstrate the corresponding code-unit change from 1.5 to 2.25. This is an isolated demonstration of the projection arithmetic, not a claim that a particular planet changes its ionized fraction by 0.5 in one step.

For real photoionization, the new electron shares the available kinetic energy with the gas. Creating another thermal particle does not provide its equilibrium kinetic energy for free. Similarly, dissociation changes particle count and binding energy together; it cannot be represented by a temperature-preserving composition reset plus an independently applied heat source without a consistent energy derivation.

**Required correction:** update species and energy consistently. Either evolve thermal energy with all reaction-energy sources and invert the new EOS after the composition change, or include formation energies in the evolved total internal energy. Do not both include formation energy in the conserved variable and add its reaction heat again.

**Acceptance:** a closed homogeneous reacting cell conserves total thermal plus chemical energy; a prescribed photoionization event deposits exactly the specified photon energy after accounting for stored chemical energy and escaping radiation. The result must converge under time-step refinement. Measure the existing projection contribution separately before replacing it.

### N2. LLF dissipation misses the negative acoustic characteristic

**Priority: P0 for robustness. Confirmed implementation defect, reproduced.**

Evidence: [Num_Fluxes.f90](../src/modules/flux/Num_Fluxes.f90), line 287:

```fortran
a1 = max(abs(vL+aL),abs(vR+aR))
```

For Euler flow the relevant bound is

\[
\alpha=\max(|v_L|+a_L,\ |v_R|+a_R).
\]

The current expression omits `v-a`. At v=-2 and a=1 on both sides, it returns 1 instead of 3. With states `(rho,v,p)=(1,-2,0.6)` and `(2,-2,1.2)`, the measured LLF mass flux is -3.5; the corrected speed would produce -4.5.

This is not confined to users selecting `Numerical flux: LLF`. `positivity_limited_fluxes` in `RK_rhs.f90` calls the same routine to repair an inadmissible high-order stage. Its claimed positivity protection therefore does not have the necessary wave-speed bound in inward-flow cells.

**Required correction:** use the full acoustic bound. Revisit the positivity argument with spherical geometry, sources, the caloric EOS, and actual time-step limits included.

**Acceptance:** test mirrored left/right flow, inward supersonic advection, near-vacuum rarefaction, and the stage-repair caller. Do not infer that this defect caused the documented base oscillation until that case is measured.

### N3. Two errors in the optional ROE star-state estimate

**Priority: P1; confined to the ROE path. Confirmed implementation defects, reproduced.**

Evidence: [speed_estimate_ROE.f90](../src/modules/flux/speed_estimate_ROE.f90), lines 76-77 and 101-104; its live caller is the ROE branch of `Num_flux`.

The two-rarefaction estimate defines z=(gamma-1)/(2 gamma). Its pressure should be

\[
p_* = \left[
\frac{a_L+a_R-\frac12(\gamma-1)(v_R-v_L)}
     {a_Lp_L^{-z}+a_Rp_R^{-z}}
\right]^{1/z}.
\]

The code takes the square root of the bracket. The bracket has pressure units raised to z, so this is also a dimensional inconsistency. For gamma=5/3, the required exponent is 5, not 1/2.

The shock estimate calculates a compression ratio but assigns it directly to `rhoL_star` and `rhoR_star`. Each ratio must be multiplied by its respective upstream density before constructing a sound speed.

The driver multiplies both upstream densities and pressures by 10 while preserving velocity and p/rho. Physical sound speeds should remain unchanged. Both affected branches fail this invariant, with the measured numbers listed in Section 2.3.

There is an additional **model limitation** in applying constant-gamma Roe eigenvectors and an averaged gamma to a composition-dependent caloric gas. Correcting these two algebraic errors does not establish a thermodynamically consistent general-EOS Roe solver. A conservative flux on each side alone is insufficient to validate its dissipation matrix.

**Acceptance:** normalization invariance, exact constant-gamma Riemann problems, entropy behavior, and molecular contact/shock tests. Using HLLC for molecular production calculations is a reasonable interim choice, subject to its own verified bounds.

### N4. H+ inherits oxygen-based finite-difference and residual scales

**Priority: P1. Confirmed implementation defect by source inspection.**

Evidence: [diffusive_photochemistry.f90](../src/modules/lower_atmosphere/diffusive_photochemistry.f90), lines 1680-1688 and 1830-1838.

The carrier set includes `ic_Hp=5`, but the Jacobian step selection distinguishes only H2 and CO. Every other carrier uses `nO_free`:

```fortran
if (kc .eq. ic_H2) then
   nref = nH_free(j)
else if (kc .eq. ic_CO) then
   nref = min(nO_free(j), nC_free(j))
else
   nref = nO_free(j)
endif
dn = max(1.0d-6*abs(nc(kc)), 1.0d-12*nref, 1.0d-300)
```

Thus the proton derivative uses oxygen as its element reference. The residual's absolute floor repeats the same classification. In an oxygen-free H/He atmosphere with initially zero H+, the derivative perturbation can collapse to 1e-300 instead of a resolvable fraction of the hydrogen density. Subtracting reaction rates evaluated at such nearby states can erase the derivative. This is particularly relevant to the very H/He molecular benchmark for which proton transport was introduced.

**Required correction:** define the element reference for every transported species in metadata and use that definition for both derivative steps and residual scales. H+ uses hydrogen. An analytic or automatically differentiated source Jacobian is preferable once the reaction representation is unified.

**Acceptance:** evaluate the actual source Jacobian in metal-free gas at zero, tiny, and moderate H+; compare directional derivatives with independently scaled perturbations. This audit did not execute that private source-Jacobian path in a full carrier step, so its planet-level convergence impact remains unmeasured.

### N5. The semi-implicit energy update returns cooling at the preceding iterate

**Priority: P1. Confirmed implementation defect and convergence risk.**

Evidence: [energy_semi_implicit.f90](../src/modules/time_step/energy_semi_implicit.f90), lines 211-295.

Two temperature iterations are always performed. Cooling is recalculated only when `iter == 1`; the second iteration changes temperature again, but the routine returns `cool_trial` from the first updated temperature. The returned `W` and `u` therefore belong to a different temperature from the returned cooling, except when the last correction vanishes.

The scalar nonlinear equation also has no final residual acceptance test. Replacing dC/dT with its absolute value changes the iteration matrix; after only two corrections this cannot be assumed to solve the implicit equation. In molecular gas, the quantity named heating contains temperature-dependent collisional reaction terms, including endothermic dissociation, but it is held fixed during this temperature update.

**Required correction:** solve the chosen energy equation to a stated tolerance, evaluate all returned rates at the returned state, and reject or subdivide a failed source step. Use the derivative of net heating at the appropriate coupled state. A safeguarded scalar solve is useful when composition is explicitly fixed; a coupled species-energy solve is needed when reactions are stiff.

**Acceptance:** stiff cooling and dissociation cell problems, falling cooling branches, multiple thermal roots, and decreasing time steps. Report positivity repair separately from energy convergence.

### N6. The carrier residual is reported before the final composition limiter

**Priority: P1. Confirmed state/diagnostic mismatch.**

Evidence: `diffusive_photochemistry.f90`, lines 1889-1893. `pct_newton_resid` is set before `limit_to_element_budget` modifies the solution. That limiter also applied the CO equilibrium ceiling until it was deleted on 2026-09-06 (item CEILING-DEL); the element rescale it still applies makes the same point.

A small reported residual is consequently not necessarily the residual of the returned carrier state. Separately, the line search can accept its final attempted damping factor without an improvement, and reaching the iteration cap does not produce a caller-visible success/failure result comparable to a time-step rejection.

**Required correction:** evaluate the residual and constraints on the state actually returned; return a convergence status. Prefer constrained steps that stay inside the admissible composition region instead of silently projecting after the solve. A projection used for round-off repair should be measured and bounded separately.

**Acceptance:** a case that deliberately activates each limiter must either satisfy the final physical residual or reject the step. Element closure alone is not an adequate acceptance criterion.

### N7. Steady acceptance needs composition convergence, not only a stable residual norm

**Priority: P1. Verification risk with an existing partial safeguard.**

The code does have a final composition-consistency loop: `steady_newton.f90`, lines 2636-2655. It repeatedly feeds the output composition into a new residual evaluation at fixed hydrodynamic unknowns. Therefore it would be incorrect to claim there is no attempt to evaluate the final state's own radiation and composition.

However, the loop stops when the scalar residual norm changes little. It does not explicitly require convergence of all species, temperature, optical depths, and diffuse-radiation rates. A stable norm can coexist with a moving weak species or an oscillating composition. During a single `eval_residual`, rates use the temperature derived from the seed composition; the post-sweep composition can change the EOS and particle count before hydrodynamic residual assembly. The final loop is intended to close that difference, so its stopping criterion matters physically.

**Required correction:** test the state change and physical algebraic residuals as well as the hydrodynamic norm. A coupled residual with its own thermodynamic state avoids treating these conditions as loosely related stopping rules.

**Acceptance:** repeated evaluation at a fixed final state must reproduce species, optical depth, temperature, heating, and cooling within their tolerances. Test deterministic evaluation and self-consistency separately: a deterministic map can still represent an incompletely converged elimination.

## 4. Molecular energy and rate consistency

### P1. H2 double ionization has an unassigned fragment-energy term

**Priority: P0 for the molecular energy budget. Confirmed accounting defect for the represented channel.**

Evidence: [util_ion_eq.f90](../src/modules/radiation/util_ion_eq.f90), lines 932-968; [h2_photo_channels.f90](../src/modules/functions/h2_photo_channels.f90); `parameters.f90`, `e_th_H2_dd`; [molecular_reaction_heat.f90](../src/modules/lower_atmosphere/molecular_reaction_heat.f90).

The represented final state is

\[
\mathrm{H_2}+h\nu\rightarrow 2\mathrm{H^+}+2e^-.
\]

Yan et al. describe 51.4 eV as a **vertical** threshold, in their published Section 4, Equation (20). A vertical threshold is not the energy stored in two infinitely separated ground-state protons relative to H2. The latter is

\[
\Delta\epsilon_{\rm chem}=D_0(\mathrm{H_2})+2I(\mathrm H)
                         =31.67494411\ \mathrm{eV}.
\]

The driver measures the difference from the code's threshold as 19.72505589 eV. In the vertical-ionization picture the repulsive nuclear configuration converts this excess into fragment kinetic energy. The code treats only h nu - 51.4 eV as available kinetic energy, applies an electron-degradation fraction to it, and does not identify an additional ion-fragment energy term. The molecular chemical-heating module deliberately excludes photon-driven reactions. Subsequent atomic recombination cannot restore the missing 19.725 eV through the represented chemical ledger.

The conclusion is an energy-accounting error for this final-state representation; it does not establish a precise electron/ion energy distribution from the threshold alone. In a collisional gas the fragment energy must thermalize or enter an explicit nonthermal reservoir. Where particles escape, their escaping kinetic energy must instead be recorded. No such recipient was found in the inspected consumers.

At 80 eV the current `chung80` split gives a double-event fraction of 0.02253775. That is a channel-weight measurement, not an atmospheric heating correction: the latter requires the local spectrum, optical depth, and ion abundance.

**Required correction:** distinguish the absorption threshold, chemical product energy, electron energy distribution, fragment kinetic energy, internal excitation, and radiated energy. Apply electron degradation only to electron energy. Include fragment heating with a collision/escape validity criterion. The approximation that reuses the single-dissociative electron-degradation row also needs replacement or a quantified range; its error is not guaranteed to be second order merely because an energy dependence is described as slow.

**Acceptance:** each photoevent must close its energy ledger, as well as H nuclei and charge. Check thresholds, intermediate photon energies, and the high-energy limit. Do not change the cross-section threshold to 31.675 eV simply to close the ledger: that would confuse absorption physics with energy partition.

Source: [Yan et al. (1998), published article](https://doi.org/10.1086/305420).

### P2. H3+ cooling violates the low-collider-density limit

**Priority: P1. Confirmed model defect outside the table's density range.**

Evidence: [h3p_cooling.f90](../src/modules/lower_atmosphere/h3p_cooling.f90), lines 134-135. The H2 density is clamped to the first table column, 1e6 cm^-3. The input density can be zero and the same non-LTE factor is still returned.

The published Table 6 density range is 1e12-1e20 m^-3. The conversion to 1e6-1e14 cm^-3 in the code is correct. The problem is extending the lowest tabulated factor as a constant below that range.

With only H2 collisional excitation represented, collisional thermal cooling must approach zero when n(H2) approaches zero. In the low-density limit its leading behavior should be proportional to n(H3+) n(H2) at fixed temperature. The measured constant value in Section 2.3 fails that limit.

Radiative or chemical pumping could produce H3+ emission without H2 collisions, but that emission is not automatically cooling of the thermal gas. It requires its own excitation source and energy ledger. It cannot justify a constant floor in a function used as collisional cooling.

**Required correction:** derive a controlled low-density continuation or solve statistical equilibrium with the relevant colliders and sources. Record out-of-domain evaluations. Replacing one published correction table with another does not by itself correct the missing asymptote.

### P3. Published H3+ polynomial joins create nonphysical discontinuities

**Priority: P1. Confirmed numerical/model defect, reproduced.**

The published Table 5 coefficients agree with the inspected code. Nevertheless, direct evaluation gives a 43.047% jump across 300 K and a 2.353% downward jump across 1800 K. These are not compiler failures or transcription mistakes. They result from using the rounded piecewise fits exactly as separate branches.

A physical emission coefficient does not have a finite jump at an arbitrary fitting boundary. Such a function can create apparent thermal barriers, distort derivative-based iterations, and make a finite-difference Jacobian depend on whether its perturbation crosses a join. Keeping the discontinuity because it occurs in a published fit does not make it physically acceptable.

**Required correction:** obtain a consistent line/level calculation or reconstruct a smooth positive cooling function constrained by the published emission data and their uncertainty. Document any discrepancy between the polynomial and the underlying table. Do not erase genuine features of a cooling curve by indiscriminate smoothing.

The Miller non-LTE calculation is itself approximate and discusses uncertainty in the excitation rate. Matching a temperature profile by selecting a correction table is weaker evidence than verifying the emitted power, collider dependence, and limiting behavior separately.

Source for P2-P3: [Miller et al. (2013), published article, Tables 5-6 and the non-LTE method](https://doi.org/10.1021/jp312468b).

### P4. H2 chemical equilibrium and the caloric EOS use different state spectra

**Priority: P1. Confirmed thermodynamic inconsistency, reproduced.**

`mol_rates::q_rovib_H2` constructs an approximate vibrational/rotational ladder using spectroscopic constants and a dissociation cutoff. The caloric EOS instead uses the level energies and degeneracies in `molecular_infrared_data`.

A single thermodynamic model must satisfy

\[
u_{\rm rv}/k_B=T^2\frac{d\ln Z_{\rm rv}}{dT}.
\]

The diagnostic differentiates the actual chemistry partition function and compares it with the actual EOS API:

| T [K] | EOS u_rv/k_B [K] | Chemistry T^2 d ln Z/dT [K] | Relative difference |
|---:|---:|---:|---:|
| 300 | 274.80745 | 273.07437 | -0.63065% |
| 1000 | 1005.61451 | 987.34034 | -1.81721% |
| 2000 | 2386.88008 | 2312.46102 | -3.11784% |
| 4000 | 6273.37129 | 5950.91964 | -5.14001% |
| 8000 | 15760.44654 | 14454.50527 | -8.28619% |

These percentages refer to the rovibrational contribution, not the total gas internal energy or an atmospheric temperature shift. The direct measurements nevertheless establish that one common thermodynamic potential does not generate both pieces.

The oxygen reverse rates introduce another H2 thermodynamic representation through Shomate data. Pairwise detailed balance within one module is useful, but it does not establish detailed balance around cycles that cross modules.

**Required correction:** choose a common internal-state/thermochemical representation. Derive partition functions, sensible energy, heat capacity, equilibrium constants, and reverse rates from it. Keep formation energies distinct from thermal excitation. Specify whether ortho/para H2 follows equilibrium populations or a frozen ratio; LTE spin populations require a physical conversion timescale.

**Acceptance:** the identity above, positive heat capacity, EOS round trips, common equilibrium constants, and closed-network cycle consistency over the intended temperature range. At densities too low for rovibrational LTE, replace the LTE assumption with an explicit justified approximation or level populations.

### P5. Oxygen chemistry lacks a closed energy budget

**Priority: P1; becomes essential for oxygen-rich or metal-rich gas. Confirmed model limitation.**

The oxygen-carrier rows and their thermochemical reverse rates are implemented and called. FUV photolysis heat is also added. However, the collisional chemical-heating routine only covers the H/He molecular network; it does not include the oxygen reactions.

`ionization_equilibrium.f90`, around lines 2094-2104, explicitly describes retained O(1D) excitation energy and its later chemical release as absent from the heat sum. The inspected energy-source paths provide no general oxygen reaction-energy contraction. The reduced oxygen network also eliminates O(1D), so its lifetime approximation and the energy released through its assumed sink must be handled together.

In addition, the associative branch involving He(2^3S) and HeH+ has an admitted omitted heat term in `molecular_reaction_heat.f90`. A small contribution in one old atmosphere is not a validity condition for every future molecular atmosphere.

**Required correction:** construct an energy ledger for every retained reaction and every eliminated intermediate. Derive it from the same species energies and reaction rates that produce composition. Record escaping photons, collisional thermalization, and excited product energy explicitly.

**Acceptance:** complete H/O/He reaction-cycle ledgers and a homogeneous energy-conservation test with photolysis disabled. With photolysis enabled, close incident/escaped/stored/thermal energy. Quantify the effect separately for solar abundance and enhanced oxygen abundance.

### P6. CO transport includes an equilibrium ceiling without a destruction rate

**Priority: P1 for C/O chemistry. Confirmed model limitation. RESOLVED
2026-09-06.** The required correction below was carried out: `carrier_source`
no longer sets `src(ic_CO) = 0` but destroys CO by He+ + CO -> C+ + O + He
(UMIST RATE22 entry 4068) and by shielded photodissociation on the
Lyman-Werner beam (Visser et al. 2009), with the domain in which a
destruction-only row is legitimate measured cell by cell and reported (item
B3b-CO); the equilibrium ceiling and its counters were then deleted from
`limit_to_element_budget`, which now applies conservation only (item
CEILING-DEL). MEASURED at the deletion: on `oxygen_chemistry` at 1000 steps
the ceiling never fired, so every output column was byte-identical across it.
The finding as written follows.

Evidence: `carrier_source` sets `src(ic_CO)=0`, while `limit_to_element_budget` reduces CO to `co_equilibrium_density` when the computed ceiling is exceeded. CO is therefore an inert transported species until an algebraic operation destroys it.

An equilibrium abundance does not determine how quickly a molecule approaches equilibrium. A parcel can retain CO above the local equilibrium abundance when transport is faster than destruction. The ceiling removes precisely this possible quenched state. Conversely, the absence of a production rate prevents relaxation toward the same equilibrium from below through that mechanism.

The code comments acknowledge this restriction. It should remain a restricted CO reservoir approximation, not be presented as a solved CO photochemical network. The change in CO binding energy when the ceiling acts also requires an explicit accounting choice.

**Required correction:** implement a physically selected CO formation/destruction mechanism, photodissociation, and shielding, or limit the supported science to conditions in which the frozen reservoir approximation is justified. A table of equilibrium CO is not a kinetic substitute.

**Acceptance:** equilibrium and quenched limits must emerge from reaction and transport timescales, with elemental and energy budgets closed. Do not use agreement with one prescribed CO profile as the only gate.

### P7. Infrared emission, thermodynamics, and radiation exchange need a common validity range

**Priority: P1/P2 depending on molecular abundance. Confirmed limitations and a measured interpolation defect.**

The implementation has real infrared calculations: H2 line emission, H2O/CO band emission and absorption, and an imposed radiation temperature with geometric dilution. It is inaccurate to call this an absent infrared treatment.

However:

- H2O and CO use LTE band means with a limited cross-section temperature range. Evaluating the surrounding emission table to higher temperature does not extend the validity of the opacity data.
- H2O, CO, OH, H2+, H3+, and HeH+ do not acquire their own rovibrational heat capacities in the current caloric EOS. Treating them as translational particles is a trace-species approximation, not a general molecular EOS.
- Cooling by spontaneous emission is equal to thermal cooling only when the excitation energy comes from thermal collisions in the assumed balance. Radiative pumping, formation pumping, and optically thick trapping change that balance.
- A band-mean cross section is not generally sufficient to obtain transmitted flux through a line-rich, optically thick column: averaging opacity and exponentiating do not commute.
- At Tgas=Trad and dilution W=1, the independently interpolated emission and absorption tables leave the small nonzero rates measured in Section 2.3. Their size is small, but the exact equilibrium identity should be built into the representation.

**Required correction:** compute a shared net-exchange expression that preserves detailed balance; implement collider-dependent statistical equilibrium where LTE fails; add optical-depth treatment appropriate to the bands and column. Expand molecular heat capacities when those species carry a material fraction of the particles. Define runtime validity diagnostics for T, density, abundance, and opacity.

**Acceptance:** zero net thermal radiation exchange in an equilibrium cavity, the optically thin limit, the diffusion limit where applicable, LTE/non-LTE limits, and integrated radiative energy conservation. Check these before comparing a planet's final temperature.

## 5. Transport, radiation, and domain limitations

### 5.1 The transported-proton steady system is genuinely unavailable

`input_read.f90`, lines 1845-1881, explicitly rejects `Ionization transport: True` without molecular chemistry/carrier transport, with `Coupled carrier solve: True`, or with `Solver: Newton`. The existing steady unknown packing in `steady_newton.f90` supports three hydrodynamic quantities and optionally one H2 density; it does not pack a proton continuity unknown.

This confirms the specific limitation in the supplied status document. It does **not** mean there is no Newton solver, no chemistry solver, or no molecular steady-solver entry point. All exist, but they do not solve the requested transported-H+ system.

The appropriate next design is a steady system with the actual active species-continuity equations, not a new special proton correction applied after convergence. H+ should work in an atomic atmosphere as well. Molecules should not be required to activate a physically atomic transport equation.

A long chemical relaxation time does not prove that a physical steady state does not exist. Nor does an unsuccessful marching calculation prove that every possible steady solver must fail. Distinguish a demonstrated failure of the present algorithm, nonexistence of a stationary solution under the chosen boundaries, and an unstable stationary solution. A bifurcation/stability study is needed to discriminate the latter two.

### 5.2 Select equilibrium species using chemical relaxation, not abundance alone

He+, molecular ions, and metal stages remain in local equilibrium in the relevant hydrodynamic paths. This can be justified for some fast intermediates. It is not justified merely because a species is rare: a rare ion can dominate electrons, catalytic H2 destruction, or radiative cooling.

For species coupled through a reaction network, the useful timescales come from the constrained chemical Jacobian. Let

\[
J_{ss'}=\frac{\partial S_s}{\partial n_{s'}},
\qquad
\tau_k\simeq\frac{1}{|\operatorname{Re}\lambda_k|},
\]

where elemental and charge-conservation null modes are excluded, and stable decaying modes are considered. Compare relevant relaxation times with transport, temperature change, and radiation-change times. A negative or complex mode structure must be interpreted rather than reduced indiscriminately to a positive lifetime.

For an advection scale L and velocity v,

\[
\mathrm{Da}_k=\frac{L/|v|}{\tau_k},
\qquad
\tau_{\rm diff}\sim\frac{L^2}{D+K_{zz}}.
\]

Large Da supports local chemical relaxation, provided temperature and irradiation vary slowly too. Small Da requires transport/history. The status document's P r/v criterion is a useful simple illustration for a photoionized species, but it is not a general relaxation criterion for a coupled molecular network. At a stagnation point use the actual forcing and diffusion times rather than interpreting an infinite r/v literally.

**Development consequence:** transport H+ and He+ where required, then retain or eliminate H2+, H3+, HeH+, and excited states according to demonstrated fast-mode limits. Algebraic elimination must preserve nuclei, charge, energy, and derivatives. The same decision should apply to atomic and molecular atmospheres.

### 5.3 The diffusion model is a restricted mixture approximation

The existing binary-element diffusion code already has neutral, polarization, Coulomb, and electric-field related routines. It is not appropriate to propose all of that as missing work. However, the neutral molecular-carrier operator specifically calls `hard_sphere_pair_diffusion` in its Blanc mixture estimate. For H+ it sets molecular diffusion to zero and retains bulk advection and eddy mixing. These are different implementations and should not be conflated.

The neutral approximation can be useful in an H2-rich, weakly ionized gas. It needs reassessment across an ionization front, for polar H2O, and whenever a transported molecule ceases to be a trace constituent. H2 itself can be most of the gas; a trace-solute interpretation of its diffusion is then inadequate without a mixture derivation.

For a physically closed mixture in the barycentric frame, diffusive fluxes must satisfy

\[
\sum_s m_s J_s=0,
\]

and charged transport must satisfy the intended current closure, with the electric field obtained consistently. Charge neutrality of cell populations is not a substitute for consistency of charge flux. The gas energy equation also needs the corresponding enthalpy transport, schematically \(\sum_s h_s J_s\), and any work terms required by the selected closure.

The inspected carrier update adjusts local atomic reservoirs to preserve element totals. That is a useful invariant check, but it does not prove that the compensating species and energy fluxes follow a multicomponent momentum balance. Derive that reduction explicitly or use a Maxwell-Stefan/Burgers-type mixture formulation with a clearly defined trace limit.

### 5.4 Species advection should use the hydrodynamic mass flux consistently

The hydrodynamic scheme uses a face Riemann mass flux and the spherical cell volume `(r_right^3-r_left^3)/3`. The carrier residual uses a cell-centered derivative of a mass-normalized fraction, with a limited correction; its diffusion divergence uses `r_center^2*dr` as the volume factor. At the base, the inflow direction is inferred from a wind-region average of `rho*v*r^2`, while the magnitude uses a local cell quantity.

The continuum equation for a mass-normalized fraction is legitimate when derived using mass continuity. The unresolved issue is whether these separate discrete operators preserve the same integrated species budget and the intended boundary species flux during nonsteady motion and grid stretching. This audit did not measure a whole-column conservation failure, so this is a verification risk rather than a reported numerical loss.

An apparent density-normalization error was specifically checked and ruled out: although the formal argument of `carrier_geometry` is named `ntot`, both the marching and steady callers pass `nrho`. Its advection coefficient is therefore based on the mass-normalization density. Reading only the callee would have produced an incorrect diagnosis.

**Recommended design:** form advective species fluxes from the same face mass flux as hydrodynamics, reconstruct the species mass fractions consistently, and use a common spherical volume. The total species mass flux should agree with the hydrodynamic mass flux. Any relative diffusion then enters as an additional explicitly balanced flux. This also makes base boundary conditions interpretable without substituting a distant average for an actual face flux.

### 5.5 Spectral conservation requires cell and band ledgers

EXHALE reads SEDs, integrates photoionization channels, computes attenuation, and has separate LW and oxygen photolysis treatments. Those paths exist. The next task is to establish that they share a consistent photon budget at overlapping wavelengths and at optical depths that change rapidly within a cell.

For a beam crossing a cell, use the actual attenuation difference to check the total absorbed photon count. When local coefficients are integrated at cell centers, refinement must recover that same count. Absorber competition must distribute the absorbed photons without creating a second copy of the beam for another process.

The LW code contains both a level-resolved shielding factor for dissociation and a band-equivalent-width treatment for transmission. They encode different observables. They need a joint photon-conservation test under the actual spectrum, temperature, and columns; a shielding-table lookup by itself does not establish a conserved radiation field. Rapid velocity gradients can also alter line overlap and shielding relative to a static column.

Important tests are monochromatic absorption, a single optical-depth jump, competing H/H2/He absorbers, the LW/oxygen overlap bands, and refinement around ionization thresholds. A spectrum that resolves the absorbed **energy** does not automatically resolve the absorbed **photon count**, and neither alone establishes the photoelectron-energy integral.

Case-A/Case-B recombination choices and helium recombination-photon coupling should be checked together with photon escape. An imposed recombination coefficient and an explicit local photon-reabsorption correction must not count the same recycling twice. The present coupling is implemented; a full closure audit of all its lines and continua was not completed here.

### 5.6 Lower-atmosphere input is more than a file-format problem

The repository has an analytic lower column and profile adapters/readers for external atmosphere models. The Python schema tracks several molecules, including CO2, CH4, NH3, HCN, and N2, as diagnostics and elemental contributors. The main EXHALE species table does not evolve all of those molecules.

Therefore reading a profile containing a molecular column does not mean that EXHALE computes that molecule's upper-atmosphere abundance. Preserving its elemental nuclei while replacing the molecular carrier by an atomic reservoir also does not preserve molecular opacity, reaction energy, or heat capacity.

At a lower/upper interface require explicit consistency of:

- Radius and reference pressure, gravity convention, and geometric area.
- Temperature, pressure, density, and the definition of particle versus nucleus mixing ratio.
- Elemental abundances in the gas, including the effect of condensed reservoirs when relevant.
- Species fluxes or a mathematically justified reduced boundary condition, not abundance matching alone.
- Radiative/enthalpy/chemical energy crossing the interface.
- Stellar irradiation, eddy diffusion, and any photochemical boundary flux imposed by the lower model.

Do not independently impose every species density and flux at the same boundary. The number of constraints must match the transport equations and incoming characteristics. A matching-pressure study should show an overlap region in which the coupled solution is insensitive to moving the interface.

For warm H/He escape, prioritize a verified H2/He network. For carbon- or oxygen-rich atmospheres, add a selected C/H/O network with CO/CO2 and relevant intermediates. For cooler reducing atmospheres, evaluate CH4/NH3/N2/HCN pathways, H- chemistry, and radiative association according to reaction flux and timescale. A long species list without verified rates, energy, and transport is not an adequate extension.

### 5.7 Continuum, common-temperature, and geometric assumptions need runtime checks

`src/utils/collisional_validity.py` already provides collision/Knudsen diagnostics; collision checking is not wholly absent. Its documented omissions and duplicated Python/Fortran coefficients should be considered when using it. The inspected main solver does not replace its fluid closure with a kinetic escape boundary when that diagnostic fails.

Track at least Knudsen number, ion-neutral/electron-heavy-particle thermal equilibration, chemical relaxation, infrared critical densities, and relevant optical depths. A common bulk velocity and temperature are valid only when momentum and thermal coupling support them. Molecular ions may require different relaxation checks from neutral H2.

Likewise, a spherically symmetric wind with a Roche potential does not become a multidimensional Roche-overflow solution. Geometry, irradiation averaging, area expansion, and the outer boundary must be stated consistently. An outflow boundary beyond a collisional sonic point is a different problem from a subsonic or collisionless upper boundary. Magnetic confinement and stellar-wind interaction are additional physical models when the target demands them, not corrections that should be hidden in a heating fraction.

### 5.8 Molecular post-processing is not self-consistent throughout the molecular layer

The advection post-processing module explicitly excludes molecular species from its solved H/He background and associated particle counts. It can carry atomic/metal advection corrections on that background, while molecular output columns describe a different approximation.

For an atomic, ionized line-forming region this may be acceptable after a quantified locality check. For a molecular base, the hydro composition, `_adv` populations, electron density, and temperature should not be presented as one fully solved state without a consistency test.

**Required correction:** obtain line-forming species from the same non-equilibrium atmosphere, or clearly restrict a post-processing calculation to a domain where neglected molecules and feedback are demonstrably irrelevant. Validate spectra only after the atmospheric state is converged under the chosen physics.

## 6. Recommended governing system

### 6.1 One species description and one reaction representation

Extend the useful existing species metadata into the owner of species mass, elements, charge, formation energy, internal-state model, and transport properties. Preserve distinct treatment of an excited level already included in its parent atomic population; do not double count He(2^3S).

Represent reactions using a stoichiometric matrix \(\nu_{s\alpha}\) and reaction progress rates \(\mathcal R_\alpha\):

\[
S_s=\sum_\alpha\nu_{s\alpha}\mathcal R_\alpha.
\]

For every reaction verify elemental and charge invariants:

\[
A\nu=0,\qquad z^T\nu=0.
\]

The same reaction object should determine species sources, reaction-energy terms, derivatives, and rate-validity diagnostics. Keep special scientific rate routines where needed; the change is that several manually maintained source/heat lists should no longer define separate versions of the chemistry.

Rate metadata should state reaction, units, collider, temperature range, source, uncertainty, reversible/irreversible interpretation, and whether a coefficient assumes LTE internal populations. Below- or above-range evaluation should be reported. An endpoint clamp is a numerical continuation, not a validated high- or low-temperature physical limit.

### 6.2 Species continuity and constraints

For each independent transported species use

\[
\frac{\partial n_s}{\partial t}
+\frac{1}{r^2}\frac{\partial}{\partial r}
 \left[r^2(n_s v+J_s)\right]=S_s.
\]

Choose independent species unknowns after accounting for elemental and charge constraints. Either eliminate dependent species analytically in an admissible domain or use a constrained nonlinear method. Avoid solving a singular set consisting of every species equation plus duplicate conservation equations.

Enforce positivity and elemental budgets during the nonlinear step. Logarithmic variables can help strictly positive concentrations, but exact zeros, nearly exhausted neutral reservoirs, and multiple elements require explicit treatment. A generic log transform alone does not solve the constraint problem.

Electron density may be eliminated through charge neutrality in a sufficiently quasineutral single-fluid model. If charged diffusion is active, determine electric-field/current closure consistently. Specify when electron temperature can be identified with the heavy-particle temperature.

### 6.3 A single energy convention

The preferred long-term formulation includes chemical formation energies in total gas energy:

\[
E=\frac12\rho v^2
 +\sum_s n_s\left[\epsilon_s^0+u_s(T,\mathbf p_s)\right],
\]

where \(\mathbf p_s\) denotes internal populations when they are not assumed to be in LTE. The pressure is obtained from the same species state. Gravity can be handled with a consistent potential-energy flux as the current hydrodynamic discretization already attempts.

In this convention, nonradiative chemistry transfers energy among components without an additional formation-energy heat source. Radiation, nonthermal particles, conduction, viscosity, and species enthalpy diffusion enter their corresponding budgets explicitly.

A thermal-energy-only formulation is also physically valid and could be the smaller initial repair. Then chemical energy must appear through exactly matched source terms, and changes of particle count/internal populations must be handled in the energy inversion. The two conventions must give the same closed-cell and radiative-event results.

This is a substantial representation change. The design and restart/input consequences should be reviewed before implementation; this report does not choose a new public API or silently redefine stored energy.

### 6.4 A shared evolution and steady residual

Define a spatial operator \(\mathcal F(U,\mathbf n,\mathcal J)\) from the common hydro/species fluxes, reaction network, energy ledger, and radiation state. Evolution integrates this operator; steady calculations solve its zero. Sources, constraints, boundary states, and radiation closure should not change when a Newton finish begins.

Use an implicit or IMEX integration for stiff chemical/thermal terms, with a well-defined source-step convergence criterion. RK3 for the homogeneous hydro stages does not make a Lie-split composition/energy scheme third order in time. Establish the actual temporal order with a coupled manufactured or reference solution.

For steady convergence, add all required species rows to the nonlinear system. A block preconditioner can exploit local chemical coupling and nearest-neighbor transport. Attenuation couples distant radii; its response must enter the residual/Jacobian action or a radiation iteration driven to a stated tolerance. A banded preconditioner need not reproduce every nonlocal derivative to be useful, but the true residual must contain the physics being solved.

Residual scaling must be attached to the species/element represented, with both absolute and relative tolerances. Report maximum local imbalance and integrated imbalance; a volume-weighted norm alone can hide a narrow chemically important front.

### 6.5 Sound speeds and the dissociation front

The present caloric EOS uses a frozen-composition heat capacity to form a sound speed. That is the relevant fast-wave limit if chemistry cannot relax during the acoustic disturbance. An equilibrium reacting gas has a different effective compressibility because ionization and dissociation absorb energy and change particle count.

An explicit hyperbolic step followed by implicit reactions can be consistent with frozen characteristic speeds, but its stiff relaxation limit must recover the equilibrium response. Check the frozen and equilibrium limits and their compatibility. Do not substitute an effective gamma into constant-gamma shock formulas and assume the resulting Riemann solver is exact for a reacting mixture.

## 7. Development sequence and concrete acceptance criteria

### Stage A: correct demonstrated errors and establish conservation

Deliverables:

1. Correct the LLF speed and both ROE algebraic errors; qualify the molecular ROE path.
2. Fix H+ derivative/reference scales through species metadata.
3. Return energy and carrier diagnostics at the actual accepted state; expose failed solves.
4. Remove the unaccounted fixed-temperature composition-energy projection through a derived conservative source update.
5. Add the H2 photoevent energy ledger, including double-ionization fragments.
6. Add missing oxygen/associative energy terms or explicitly restrict the affected options until their treatment is complete.

Acceptance: source-level invariants and reduced execution tests pass, integrated cell energy closes, and no accepted state relies on unreported floors/projections. Run only the regression cases touching each correction, then expand coverage when a shared energy or flux change reaches more branches. Preserve old outputs as historical evidence; numerical identity with them is not the physical acceptance condition.

### Stage B: unify molecular thermodynamics and source integration

Deliverables: a common H2 state model; compatible oxygen thermochemistry; explicit rate validity; continuous H3+ emission and its low-density limit; a converged coupled species-energy source solve.

Acceptance: common thermodynamic identities, energy/charge/element invariants, correct reaction equilibrium and frozen limits, and systematic convergence with source tolerances and time step. Verify short-lived eliminated intermediates against an explicitly integrated reference cell network.

### Stage C: converge the non-equilibrium escape equations

Deliverables: a shared active-species residual; H+ and necessary He+ transport in atomic and molecular gas; Newton/PTC integration with the same species and energy state; nonlocal radiation closure and composition convergence gates.

Acceptance: a residual-converged molecular wind from multiple initial conditions where the physical branch is stable, including independent tolerances for mass, momentum, energy, charge, and species. A result must survive a fresh final-state residual evaluation. If no steady solution is found, report whether the obstacle is nonlinear convergence, boundary incompatibility, or physical instability; investigate instead of declaring a fixed-step snapshot converged.

### Stage D: resolve boundaries and discrete transport consistency

Deliverables: common face mass/species fluxes and spherical volumes; boundary conditions matching incoming characteristics and specified composition flux; a verified hydrostatic balance; consistent diffusion/enthalpy fluxes.

Acceptance: zero-irradiation hydrostatic atmospheres remain stationary to truncation error, tracer/element column budgets close, and base behavior converges as the grid and time step are refined. A smaller local pseudo-time step can accelerate or stabilize a steady iteration, but it does not by itself prove that the boundary equation is correct. Use global physical time steps for transient claims.

### Stage E: extend molecular radiation and composition for the target class

Deliverables: a selected kinetic CO/C/O network; additional species only as required by reaction flux and the lower-atmosphere composition; molecular heat capacities, non-LTE cooling, and line/band transfer within verified ranges.

Acceptance: equilibrium, photochemical, advective-quench, optically thin, and optically thick tests. A lower/upper coupling benchmark must conserve elements and energy across the interface and be insensitive to a reasonable change of matching pressure within an overlap region.

### Stage F: external atmospheric validation and observable predictions

Only after the corresponding conservation and convergence gates should published planet profiles serve as the acceptance target. Compare temperature, velocity, neutral and ion densities, energy terms, species source balances, and mass flux together. Spectra and transit depths are downstream tests, not substitutes for these balances.

Track discretization error, iterative error, rate/data uncertainty, and model discrepancy separately. One set of tuned heating parameters should not absorb all four.

## 8. Verification matrix

The following are proposed tests, not tests reported as executed in this audit. Numerical tolerances should be justified against the scales below; the suggested values are starting targets rather than universal physical constants.

| Test | Setup | Required property | Main implementation reached |
|---|---|---|---|
| Reaction invariants | Every retained reaction | Exact integer element/charge closure; consistent energy recipient | Species metadata, reaction source assembly |
| Thermodynamic identity | Temperature grid covering the declared range | u=T^2 k_B d ln Z/dT, positive heat capacity, compatible K_eq | EOS and reverse rates |
| Isolated reacting cell | No transport/radiation; several initial molecular fractions | Constant total thermal+chemical energy; correct equilibrium | Coupled source update |
| Photon-event test | One absorber/channel, prescribed photon energy | Photon energy equals heat+stored+escaped energy | Photo channels and degradation |
| Finite-volume absorption | Thin/thick cells and competing absorbers | Absorption equals incoming minus outgoing photons | SED, opacity, optical-depth integration |
| Infrared equilibrium cavity | Tgas=Trad, W=1 | Net exchange zero to interpolation/round-off tolerance by construction | Molecular IR and H3+ field treatment |
| H3+ low-density limit | Fixed T,n(H3+), decreasing n(H2) | Correct collider scaling; no artificial floor | H3+ non-LTE treatment |
| Frozen passive scalar | Smooth flow and a transported composition feature | Correct speed/order and integrated tracer conservation | Hydro/species flux coupling |
| Reaction-advection test | Prescribed wind, analytic or independently integrated kinetics | Local equilibrium at large Da; frozen/advected state at small Da | Transport plus chemistry |
| Charged diffusion | Neutral-dominated and ionized mixtures | Correct limiting diffusion, zero chosen current, no net diffusive mass | Mixture diffusion/electric field |
| Riemann and mirrored-flow tests | Atomic shock/rarefaction/contact; reversed velocities | Wave-speed bounds, conservation, entropy and positivity | HLLC/LLF/ROE |
| Molecular contact/shock | Jump in H2 fraction and heat capacity | No spurious pressure creation; consistent shock energy | EOS, reconstruction, fluxes |
| Hydrostatic column | No irradiation or outflow drive | Gravity-pressure balance and no growing base mode | Grid, sources, boundaries |
| Source-step refinement | dt,dt/2,dt/4 at fixed spatial grid | Measured temporal convergence; accepted nonlinear residual | Split/IMEX integration |
| Spatial refinement | At least three grids, independently refined fronts/base | Quantified profile and integral convergence | Spatial scheme and boundaries |
| Final-state residual | Reload converged solution; rebuild rates/columns | Same physics, consistent composition and small residual | Solver acceptance and restart |
| Lower/upper interface | Move match pressure through overlap region | Flux/energy continuity and stable upper solution | Profile handoff and boundary solve |
| Molecular benchmark | Fully stated Model A-like inputs | Converged profiles and budgets; documented remaining model differences | Complete molecular wind |

For discrete element/charge checks, exact stoichiometry should normally give round-off-level closure. Integrated conservation targets should include boundary flux and source integration error rather than comparing density sums alone. The existing 1e-5 steady residual target can be retained provisionally, but it must specify normalization, the included equations, and whether the maximum local error also passes.

A convergence study should report errors in quantities that control the science: mass loss, front location, integrated absorption, thermal peak, key ion columns, and line-forming populations. Global L1 error may look small while an unresolved narrow layer controls the answer. Refine the base scale height, chemical front, and radiation absorption length separately where possible.

Acceptance tests should fail with a nonzero exit status when their stated criterion is violated. Several existing diagnostic drivers print balances without enforcing all of them. Preserve their useful diagnostics, but distinguish a numerical report from an automated gate. Make the required tests discoverable from the build/test entry points instead of depending on obsolete compilation instructions or private scratch directories.

## 9. How to interpret the Koskinen comparison

The published Koskinen et al. (2022) method was inspected directly, including the upper-atmosphere description and Appendix B. It transports the specified H/He molecular/ionic species and describes photoelectron heating using h nu - I, as well as H3+ cooling with a non-LTE treatment. Its Appendix B exists in the final journal version.

The status document reports a factor of 2.5-3 discrepancy in inferred heating at similar ionization rate and attributes it to the external implementation. This audit has not reconstructed that comparison numerically. The currently inspected EXHALE implementation also has its own unresolved energy issues. Therefore the evidence does not justify treating the entire remaining discrepancy as a confirmed error in the published code.

Reconstruct the comparison in this order:

1. Fix geometry, reference radius, area averaging, and the pressure/temperature/composition at the lower boundary.
2. Use the same dimensional incident spectrum and integration limits; compare absorbed photon counts and energy before solving a wind.
3. At identical prescribed n_s(r), T(r), and radiation, compare photoionization, electron energy, fragment energy, chemical heating, and each cooling term.
4. Compare the reaction network and rate coefficients at those states, including the handling of internal excitation and reverse reactions.
5. Solve a prescribed-wind transport problem before comparing two fully coupled winds.
6. Converge the final atmosphere and demonstrate spatial/source-step convergence.

The paper states that ionization-potential energy is lost through subsequent recombination or chemistry in its adopted treatment. EXHALE's collisional molecular heating is a different physical choice when nonradiative molecular reactions return part of that energy to the gas. A physically improved model may intentionally differ from the benchmark. Describe that difference explicitly rather than treating agreement as the sole standard.

Using a fraction of full photon energy can be a diagnostic bracket, but it is not a general correction to an unknown energy ledger. Likewise, forcing the He+ profile or H3+ cooling to fall within digitized points is not a substitute for verifying the transport and level-population equations. Observed/model profile ratios combine abundance, temperature, radiation, and normalization; they do not identify one coefficient uniquely.

Source: [Koskinen et al. (2022), final journal article](https://doi.org/10.3847/1538-4357/ac4f45). The status document's quoted profile ratios and mass-loss comparisons remain historical claims in this report.

## 10. What is already useful, and what this audit does not establish

The existing code has several sound foundations that should be retained while removing inconsistent interfaces:

- Species metadata and element-census checks distinguish excited levels from independent nuclei.
- Four explicit H2 absorption channels close nuclei and charge in direct execution.
- H2 association/dissociation coefficients have an implemented detailed-balance relation, and oxygen reverse rates are constructed from thermochemistry rather than independently guessed expressions.
- The caloric EOS computes H2 internal energy and heat capacity and provides energy/pressure conversions; it is more than a single altered gamma.
- H3+ and other molecular radiation routines are linked and callable, with identifiable data origins.
- The carrier operator, constrained equilibrium routines, JFNK/PTC solver, residual diagnostics, and lower-profile interfaces provide useful components for a common coupled system.
- The final Newton path attempts composition consistency, and the collision diagnostic already exists. These should be strengthened rather than inaccurately described as missing.

These positive checks do not validate every rate coefficient, every physical approximation, or every planetary solution. This audit also does not establish the unique cause of the reported base limit cycle, the final effect of any defect on mass loss, the existence or stability of every requested molecular steady state, or the correctness of the comparison authors' unpublished implementation.

The immediate scientific milestone should be a molecular atmosphere whose species, radiation, and energy budgets are internally closed and whose final state satisfies the same equations used during evolution. Extending that verified core to more molecules and observables is then a controlled physical development, with a stated range of validity.

## 11. Literature and reproducibility notes

NASA ADS was queried before general web access for the relevant published records. The ADS results confirmed the journal bibcodes and DOIs below. Publisher access through the web tool was unsuccessful; the final journal PDFs already present in the shared `references/` directory were read using `pdftotext -layout`. No arXiv version was used as the content source.

| Reference | ADS bibcode | Published source inspected | Use in this report |
|---|---|---|---|
| Koskinen et al. (2022) | `2022ApJ...929...52K` | `references/Koskinen_2022_ApJ_929_52.pdf` | Upper-atmosphere method, Appendix B, comparison interpretation |
| Miller et al. (2013) | `2013JPCA..117.9770M` | `references/Miller_2013_JPCA_117_9770.pdf` | H3+ cooling polynomials, density units, non-LTE assumptions |
| Yan et al. (1998) | `1998ApJ...496.1044Y` | `references/Yan_1998_ApJ_496_1044.pdf` | Vertical double-ionization threshold and channel distinction |

The `references/` paths in this table are relative to the shared ExoAtmosphere directory, one level above EXHALE_v1.00. ADS also returned the Hollenbach and McKee (1979) record, but this report does not treat its detailed heating derivation as independently checked. Other provenance statements embedded in rate-module comments were not automatically promoted to verified literature claims.

The retained diagnostic script compiles current source directly, rather than trusting the timestamp or revision label of an old executable. It writes only new temporary build/log products. Its output is intended to be read alongside the severity and scope qualifications in this report. All proposed physics changes and extended validation remain future work.
