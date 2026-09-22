# Hydrodynamic Methods in Frelikh, Yabe, and Kuramoto: Assessment for EXHALE

Date: September 21, 2026

Scope: hydrodynamic discretization, spherical geometry, gravity, boundary conditions, conservation, numerical stability, computational cost, and an adoption strategy. Atomic and molecular physics models, reaction rates, radiation physics, and their scientific adequacy are intentionally excluded. An equation of state and prescribed energy sources appear only as interfaces required to close and test the fluid equations.

## 1. Recommendation

**A controlled experiment is worthwhile; an immediate replacement of EXHALE's hydrodynamic solver is not justified by these papers.**

The most useful ideas are conservative semi-Lagrangian mass transport, compact interpolation with additional moments, staggered pressure-velocity coupling, and a pressure-implicit treatment of acoustic stiffness. These are separate ideas, not one interchangeable "CIP solver."

My recommended order is:

1. Establish a hydrodynamics-only comparison against the current EXHALE operator, with identical geometry, gravity, thermodynamic closure, forcing, and physical boundary conditions.
2. Implement and test spherical CIP-CSL2 transport in isolation.
3. Build a staggered, explicit hydrodynamic comparison following the Kuramoto strategy. Label its conservation limitations explicitly if momentum and energy use original CIP.
4. If explicit acoustic steps are the measured bottleneck, test a conservative pressure-implicit extension inspired by CCUP. Also consider applying a pressure correction to the existing finite-volume operator: adopting CIP advection is not a prerequisite for addressing acoustic stiffness.
5. Promote a new method only after conservation, hydrostatic balance, positivity, boundary sensitivity, and accuracy at matched computational cost have passed explicit acceptance gates.

For a production replacement, I recommend **conservative mass, compatible momentum and total-energy updates, and balanced pressure-gravity coupling**. That is a proposed extension of the cited methods, not a claim that Frelikh or Kuramoto already implemented that exact combination.

For speed, explicit CIP offers a possible reduction in transport cost and numerical diffusion, but it does not remove acoustic stiffness. A pressure-implicit method has a stronger potential benefit in a very subsonic atmosphere, at the cost of a more substantial derivation and validation effort. Neither method establishes that a difficult steady solution exists or is dynamically stable.

## 2. Evidence and limits of this review

### 2.1 Sources actually examined

The supplied PDFs were found in the workspace-level directory named references, rather than in EXHALE_v1.00/docs:

| Label | Published paper | Material used |
|---|---|---|
| [F] | Frelikh and Murray-Clay, 2026, ApJ 996, 96 | Fluid equations; Sections 2.5-2.8; Appendix A, including rendered pages 28-29 |
| [Y91] | Yabe and Aoki, 1991, Computer Physics Communications 66, 219-232 | CIP interpolation and gradient evolution; Sections 7-8 on hydrodynamics and geometry |
| [Y01M] | Yabe, Tanaka, Nakamura, and Xiao, 2001, Monthly Weather Review 129, 332-344 | Original CIP, conservative integral constraints, CIP-CSL4, CIP-CSL2, and characteristic tracing |
| [Y01J] | Yabe, Xiao, and Utsumi, 2001, Journal of Computational Physics 169, 556-593 | Primitive-variable splitting, shock treatment, Section 4.2 on CCUP, Section 6 on conservative transport |
| [K13] | Kuramoto, Umemoto, and Ishiwatari, 2013, Earth and Planetary Science Letters 375, 312-318 | Section 2 numerical method; Sections 3.1-3.2 validation and energy balance |

The PDFs were read through layout-preserving text extraction. The suspect formulas in Frelikh's Appendix A were also checked against rendered journal pages, not inferred from text extraction alone. NASA ADS returned publication metadata and DOI records for all five papers. Attempts to open three DOI landing pages through the web tool failed; the supplied final journal PDFs remain the evidence for the equations and methods discussed here.

### 2.2 Implementation examined

The review inspected the current working files, not just plans or earlier review documents. The repository HEAD was 93eed86, but the working tree contained extensive existing changes, including new conservation-reporting code. Consequently, HEAD alone does not identify everything reviewed.

Important implementation evidence includes:

- [Main integration and solver selection](../src/EXHALE_main.f90): three hydrodynamic RK stages around lines 2334-2430; steady JFNK/PTC calls around lines 7122-7125.
- [Conservative and primitive conversions](../src/modules/functions/UW_conversions.f90).
- [Reconstruction](../src/modules/states/Reconstruction.f90) and [numerical fluxes](../src/modules/flux/Num_Fluxes.f90).
- [Spherical flux divergence and gravitational energy work](../src/modules/time_step/RK_rhs.f90), especially lines 190-275.
- [Geometric and gravitational sources](../src/modules/states/Source.f90).
- [Explicit time-step calculation](../src/modules/time_step/eval_dt.f90).
- [Base face closure](../src/modules/states/base_boundary.f90), especially characteristic_base_face_state.
- [Ghost and reconstructed boundary states](../src/modules/states/Apply_BC.f90).
- [Stationary operator selection](../src/modules/states/stationary_operator.f90), [residual assembly](../src/modules/time_step/steady_residual.f90), and [certification](../src/modules/time_step/certification.f90).
- [Viscosity and conduction](../src/modules/time_step/viscous_conduction.f90).
- [Artificial low-Mach stress](../src/modules/flux/low_mach_dissipation.f90).
- [Conservation-budget writer](../src/modules/time_step/conservation_budget.f90).
- [Build source list](../Makefile), which includes the active hydrodynamic, steady-solver, and transport modules.

No production source code was modified. No EXHALE atmospheric simulation, full regression suite, runtime profile, or new CIP solver was executed. The only executed numerical check was a small independent algebra test of published formulas, described in Section 11. Performance comparisons below are therefore predictions and proposed measurements, not benchmark results.

## 3. What the papers actually implement

### 3.1 Original CIP is not CIP-CSL2, and neither implies CCUP

| Method | Additional information carried | Main advantage | What it does not establish |
|---|---|---|---|
| Original CIP | Point values and spatial derivatives | Compact, accurate characteristic interpolation with relatively low transport diffusion | Exact conservation of cell-integrated mass, momentum, or total energy |
| CIP-CSL2 | Point values and cell integrals; a quadratic transported profile obtained from a cubic integral profile | Conservative transport of the selected integral quantity | Conservation of other equations that are still advanced with original CIP |
| CIP-CSL4 | Point values, derivatives, and cell integrals | Higher-degree constrained interpolation with conservative integral transport | Automatic positivity, entropy stability, or low cost in every application |
| CCUP | A pressure solve coupled to density, velocity, and thermal variables | Implicit acoustic coupling and a route toward low-Mach calculations | Automatic conservation, correct shock energetics, or unlimited accurate time steps |

The number in CSL2 describes the quadratic profile construction; it should not be used as a blanket statement that a complete split, nonuniform-grid fluid solver has a particular convergence order.

For a scalar equation

\[
\partial_t f+v\partial_r f=G,
\]

original CIP uses both f and its gradient g. The gradient satisfies

\[
\partial_t g+v\partial_r g
=-(\partial_r v)g+\partial_rG.
\]

Interpolation alone is not the entire method: the source-gradient and stretching terms, departure-point calculation, and boundary data must all be consistent. Replacing the evolved gradient by an unrelated finite difference after every step produces a different method.

### 3.2 Frelikh: staggered, explicit, primitive-variable CIP

Frelikh evolves density, radial velocity, and specific internal energy, together with gradients. Density and energy occupy the regular grid, while velocity occupies staggered locations. The hydrodynamic update is split between advection and explicit nonadvection terms. Artificial pressure is applied in compression. These statements follow Sections 2.5-2.7 and Appendix A of [F].

The fiducial grid has 500 physical regular cells, with four regular-grid ghosts; the cited base spacing and geometric stretching ratio are case settings, not general recommendations for EXHALE. The paper's outer-boundary argument assumes a supersonic outflow with the boundary sufficiently beyond the sonic point. Its lower boundary fixes density information and pressure, while extrapolating velocity.

**Do not describe this implementation as a demonstrated exactly conservative CIP-CSL2 solver or as an implicit acoustic solver.** Those conclusions do not follow from the published method.

The paper reports that artificial pressure vanishes in its final monotonically accelerating, expanding winds. That is a property of those solutions and the compression sensor, not a theorem about all stagnant, reversing, shocked, or weak-wind atmospheres.

### 3.3 Kuramoto: a mixed conservative/nonconservative strategy

Section 2 of [K13] makes an important distinction:

- CIP-CSL2 is used for the continuity equation.
- Original CIP is used for advection in the momentum and energy equations to reduce computational expense.
- Time integration is explicit.
- Thermal diffusion uses centered spatial differences.

The paper writes conservative spherical Euler equations, including total energy. However, writing a conservative differential equation does not by itself prove that every discrete update conserves its cell integral. In particular, the authors identify continuity as the equation receiving CIP-CSL2.

The published isothermal tests use escape parameters 5, 15, and 25. Section 3.1 reports mass-flux errors below 2% with 1000 grid points. Section 3.2 reports integrated energy-budget discrepancies below 3% over the cases studied. **These are quoted published results, not measurements reproduced in this review.** They support the usefulness of that implementation but are not an exact energy-conservation guarantee or an accuracy bound for EXHALE.

### 3.4 Yabe: the broader method family

[Y91] provides the primitive-variable hydrodynamic split, artificial viscosity, and extensions to cylindrical and spherical geometry. [Y01M] explains how integral constraints yield exact conservation of a transported scalar. [Y01J] treats pressure-based acoustic coupling as a separate extension.

Their published demonstrations justify investigating these tools. They do not establish a universal superiority over a contemporary, balanced, contact-resolving finite-volume solver on EXHALE's particular problem.

## 4. What EXHALE currently does, and why that changes the comparison

### 4.1 The current method already carries conservative information

The actual conversion routines store

\[
U=(\rho,\rho v,E),\qquad E=\rho e+\tfrac12\rho v^2,
\]

and obtain primitive states as density, velocity, and pressure. RK_rhs computes shared face fluxes, then forms their spherical divergence with

\[
A_f=r_f^2,\qquad
V_j=\frac{r_{j+1/2}^3-r_{j-1/2}^3}{3}.
\]

The common solid-angle factor is omitted consistently. For a full spherical mass-loss rate, it must be restored as 4 pi; comparisons must not mix a flux integrated over solid angle with a flux normalized by solid angle.

The main loop invokes three SSP-RK stages. Reconstruction includes PLM and WENO3, and numerical flux choices include LLF, HLLC, and Roe, with an HLLE branch for inadmissible Roe states. The existing positivity repair replaces shared face fluxes and can reject a step. This is substantially different from the older Lax-Friedrichs comparison motivating [K13].

The momentum discretization requires care: pressure is placed differently in the PLM and WENO3 branches, and the well-balanced branch combines equilibrium pressure and gravity algebraically. "Replace Num_flux" is therefore not a sufficient description of replacing this hydrodynamic operator.

### 4.2 Existing low-Mach and hydrostatic treatments must be included

Code inspection confirms:

- A well-balanced option uses hydrostatic pressure departures in reconstruction, numerical fluxes, and the momentum row.
- A Roe-specific low-Mach velocity-jump option changes acoustic dissipation. It does not change the explicit acoustic time-step restriction.
- An optional artificial stress acts through momentum and energy fluxes. It is distinct from the compressive artificial pressure used by Frelikh.

These features must be documented in the benchmark configuration. Their existence does not show that they are enabled in every run or that every low-Mach problem is solved.

Likewise, a staggered arrangement prevents the elementary collocated central-pressure checkerboard described in Appendix A.3 of [F], but this does not prove that EXHALE's observed alternating modes have that same cause. EXHALE uses reconstructed Riemann fluxes, not merely the central pressure derivative used in that illustrative argument.

### 4.3 Gravity and energy already have a coupled discrete treatment

RK_rhs adds

\[
D_{\Phi,j}
=A_+F_{\rho,+}(\Phi_+-\Phi_j)
-A_-F_{\rho,-}(\Phi_--\Phi_j)
\]

inside the energy divergence. The source routine returns zero for its explicit energy source. Thus, gravitational work is carried by the mass flux rather than added a second time as a separate energy source.

For a time-independent potential, adding Phi_j times the mass row to the energy row gives a conservative flux for E + rho Phi. This algebraic property is useful and should survive a redesign. It applies to the hydrodynamic rows; it is not a claim that every split stage, boundary operation, or local pseudo-time update automatically satisfies a physical-time global budget.

### 4.4 Physical-time marching and stationary solving are different targets

The main program calls steady JFNK or PTC solvers, and stationary_operator explicitly selects WENO3 for stationary evaluation. A CIP state is therefore not automatically a solution of the current stationary operator.

There are two legitimate uses:

1. CIP produces an initial state, which is then solved and certified with the existing EXHALE operator.
2. CIP defines a new discretization, which has its own discrete fluxes, residual, and certification.

A small old-operator residual is an additional cross-check for the second use, not its definition of success. Conversely, a small change between successive CIP steps does not certify the existing WENO3 equations.

The current time-step routine also supports local pseudo-time intervals. With different intervals in adjacent cells, their shared instantaneous flux need not cancel in the sum of finite state changes over one iteration. That can be a useful stationary acceleration, but it is not conservative evolution over one common physical interval. Conservative semi-Lagrangian remapping must use a common interval unless a conservative asynchronous flux exchange is separately derived.

## 5. Advantages and disadvantages for EXHALE

### 5.1 Advantages worth testing

**Lower contamination of very small advective fluxes.** In an illustrative first-order diffusive flux, a density jump produces a contribution of order c Delta-r |d rho/dr|. Its ratio to rho |v| scales roughly as Delta-r/(M H), where H is a density scale and M = |v|/c. This explains why diffusion can overwhelm a weak wind. It is not the error formula for the current well-balanced Roe/WENO3 implementation.

A CIP-CSL2 comparison can determine whether EXHALE's weak-flow transport error is reduced at a fixed resolution. A valid comparison must distinguish the numerical face mass flux from a cell-centered product rho v r squared. Conservation of the former does not guarantee accurate reconstruction of the latter, and disagreement between them is not automatically a mass leak.

**Direct pressure-velocity coupling on a staggered grid.** Adjacent pressures determine face acceleration without the simplest even-odd decoupling of collocated central differences. This is a plausible benefit for nearly stagnant layers, subject to testing of the full spherical operator and boundaries.

**Independent discretization as a diagnostic.** Agreement between two genuinely different methods after refinement is more informative than agreement between two parameter settings of one method. It can separate a spatial-discretization problem from a nonlinear stationary-solver problem.

**Potentially efficient transport.** Compact interpolation avoids a full Riemann calculation during the advection substep and may reduce the number of expensive state conversions. This is a performance hypothesis, not a measured advantage.

**A route to implicit acoustics.** CCUP directly addresses pressure coupling. In a low-Mach atmosphere, this can matter more for runtime than changing the interpolation polynomial.

### 5.2 Costs and risks

**Conservation can regress.** Original CIP transports point information accurately without enforcing the required cell integrals. Mass-only CSL2 does not protect momentum and total energy. Long relaxation times, large gravitational energy changes, and small net energy imbalances make this significant even if a short test looks smooth.

**Artificial viscosity remains part of the method.** Original CIP is not universally stable without shock treatment. When artificial pressure q is added, momentum and energy must use compatible pressure work. In an internal-energy formulation, the compressional contribution includes -q div(v); using q only to damp momentum omits the corresponding irreversible heating.

**Low diffusion can reveal or retain unwanted modes.** A nearly stationary contact is almost unmoved by advection. Accurate interpolation alone cannot damp a pressure-entropy-velocity mode at v approximately zero. Gradient overshoots, boundary reflections, or inadequate source coupling can remain.

**Staggering is a state-layout change.** Cell-averaged momentum cannot generally be identified with rho at a center times an arithmetic average of neighboring face velocities. A definition of momentum control volumes and kinetic energy is required. Repeated conversions between layouts can create energy errors or erase the information CIP was intended to preserve.

**Nonuniform-grid details are substantial.** Hermite coefficients, source derivatives, gradients, characteristic tracing, cell integrals, and boundary extrapolation must use actual coordinates. Reusing a uniform-grid formula with an ambiguous Delta-r is not sufficient.

**Higher interpolation degree does not imply high global order.** A first-order source split, inaccurate departure points, or first-order boundaries can dominate a cubic reconstruction. Temporal refinement is required separately from spatial refinement.

**It does not cure an inappropriate boundary model.** A solver change cannot determine the correct reservoir condition, select a physically absent transonic branch, or replace a required exterior pressure.

## 6. Computational speed: what is plausible and what must be measured

### 6.1 Separate the cost of a step from the number of steps

For a physical duration T, a useful cost model is

\[
t_{\rm wall}\simeq
N_{\rm accepted}\left(C_{\rm transport}+C_{\rm pressure}
+C_{\rm sources}+C_{\rm boundaries}\right)
+C_{\rm rejected}+C_{\rm diagnostics}.
\]

For a steady solution, the number of accepted steps is also affected by relaxation time, boundary reflections, physical instabilities, splitting error, and the stopping criterion. A cheap update that requires many more steps is not a faster solver.

The measured comparison should be wall time to a specified error or certified stationary residual, not wall time for an equal number of iterations.

### 6.2 Explicit CIP: possible cheaper updates, the same acoustic bottleneck

The current hydrodynamic step makes three reconstruction/Riemann/RHS evaluations. A split CIP implementation may do less work in its transport stage. However, it also evolves gradients, evaluates source-gradient corrections, interpolates between staggered locations, and performs boundary and positivity checks. A nominal comparison of "three RK stages versus one interpolation" misses much of this work.

More importantly, both Frelikh's and Kuramoto's stated hydrodynamic integrations are explicit in their pressure-related evolution. Their complete fluid solver remains subject to an acoustic stability scale of order

\[
\Delta t_{\rm acoustic}\sim
\min_j \frac{\Delta r_j}{|v_j|+c_j},
\]

with geometry and scheme-dependent constants. EXHALE currently uses spherical volume-to-face-area weights, not simply a Cartesian cell width.

The large-Courant capability discussed in [Y01M] concerns semi-Lagrangian advection with appropriate multi-cell departure tracing. It does **not** authorize an explicit pressure update to take arbitrarily large acoustic steps.

Consequently, explicit CIP should initially be expected to offer, at most, an unmeasured improvement in update cost or resolution efficiency. An orders-of-magnitude runtime claim would require evidence of a much smaller required grid, much faster relaxation, or a genuinely different treatment of the stiff terms.

### 6.3 CCUP: stronger speed potential, additional work and restrictions

For the same cell and length scale in a low-Mach region, the ratio between an advective interval and an acoustic interval can approach c/|v|. That ratio is a possible source of speed improvement for an implicit acoustic method, not a forecast for the entire atmosphere.

Actual improvement can be smaller because:

- A different cell or the sonic region may determine the global interval.
- Pressure accuracy can require a smaller interval than stability alone permits.
- Gravity, thermal transport, prescribed source variation, and boundary motion impose additional scales.
- Nonlinear pressure corrections may require multiple iterations.
- Rejected updates and conservation corrections have a cost.
- A stationary solve may already bypass physical acoustic marching.

A one-dimensional pressure equation with nearest-neighbor coupling normally gives a tridiagonal linear system. A linear solve costs O(N), which is favorable. However, nonlinear iterations, coefficient reconstruction, boundary coupling, and global constraints can change that simple cost picture. The pressure solve must be profiled as part of the complete update.

CCUP should therefore be compared both with explicit EXHALE and with its existing stationary JFNK/PTC route when the scientific target is a steady atmosphere.

### 6.4 Cost outside hydrodynamics limits any total speedup

Let f be the measured fraction of runtime that a change accelerates, and s its speedup. With everything else unchanged, the total speedup is

\[
S=\frac{1}{(1-f)+f/s}.
\]

For illustration only, accelerating a component that occupies 10% of runtime by a factor of ten gives S approximately 1.10, not ten. These are hypothetical inputs, not timings measured in EXHALE.

This review excludes the physical content of the other operators, but their execution time still matters. Profile their aggregate cost without changing their equations during the method comparison.

### 6.5 Memory and parallel execution

Original CIP adds derivative fields to the advected point variables. CSL2 instead carries integral and point information. The increase in total program memory is not a universal factor of two because other state arrays and solver workspaces already exist. Measure actual peak memory.

The current RK_rhs already uses separate OpenMP face and cell loops with a synchronization point between them. A new local remap also parallelizes naturally, provided every shared face transfer is computed once and reused.

Potential limitations include:

- Extra gradient and moment arrays increase memory traffic.
- Several short parallel regions can cost more than their work on a small radial grid.
- A standard tridiagonal pressure solve is sequential even though its assembly is parallel.
- A large advection Courant number expands the departure stencil and, with domain decomposition, communication requirements.
- Boundary state, gradient updates, and global acceptance checks need deterministic ordering.

For radial grids with hundreds to a few thousand cells, test one, two, four, and additional useful thread counts rather than assuming that more threads help. Parallel execution of independent atmospheres may be more efficient than adding threads to a small pressure solve; that is a scheduling hypothesis to measure, not a result established here.

### 6.6 Required timing report

Report at least:

- Wall time to a common physical end time and, separately, to a common stationary acceptance target.
- Accepted and rejected updates, characteristic traces, pressure solves, and nonlinear correction iterations.
- Time in reconstruction/remapping, pressure coupling, other source operators, boundaries, and diagnostics.
- Density, velocity, temperature, mass-loss, and energy-budget errors at that time.
- Grid size, time-step policy, thread count, compiler options, and peak memory.

Use repeated timed runs after correctness is established, report variability, and keep diagnostic output settings identical. No such timing experiment has been completed in this review.

## 7. Stability: distinct questions that need separate answers

| Issue | Explicit original CIP / Kuramoto-style hybrid | Pressure-implicit extension | Acceptance requirement |
|---|---|---|---|
| Acoustic stiffness | Remains | Can be substantially reduced | Time-step sweep with phase and damping errors |
| Advection over several cells | Possible only with proper departure tracing and integral remap | Still requires accurate tracing | Constant- and variable-velocity tests |
| Negative density or pressure | Polynomial overshoot is possible | An implicit pressure solve does not guarantee positivity | Conservative limiting or rejection; no silent floors |
| Shocks | Requires suitable dissipation and energy work | Still requires entropy-producing shock treatment | Shock speed, jumps, and total-energy budget |
| Hydrostatic equilibrium | Not automatically preserved | Not automatically preserved | Pressure-gravity cancellation on the actual grid |
| Nearly stationary alternating modes | Staggering may help; advection alone is insufficient | Pressure coupling may help; boundary modes can remain | Mode-growth tests about representative steady states |
| Conduction | Centered spatial differences alone do not remove a diffusive time-step limit | Needs its own consistent treatment | Thermal boundary and positivity tests |
| Reversal and boundary transitions | Need correct characteristic data and derivative states | Pressure boundary rows must represent the same physics | Reversing-flow tests and domain-extension tests |

### 7.1 Conservation, positivity, and stability are not synonyms

Conservative CSL2 can still overshoot. A positive result can still have the wrong energy budget. A linearly stable implicit update can still strongly damp a physically important mode or converge to a splitting-dependent steady state.

For a limited CSL reconstruction, preserve its cell integral while constraining the profile. Clipping point values or cell averages independently after remapping can destroy either the integral constraints or the relation between moments. If a repair changes a conserved quantity, record that change as a budget defect rather than concealing it.

### 7.2 Artificial pressure must use the spherical compression rate

Use

\[
\theta=\frac{1}{r^2}\partial_r(r^2v),
\]

not just dv/dr, to define volumetric compression in a radial flow. A possible artificial pressure has the schematic form

\[
q=\rho\left(C_1c\,\ell\,[-\theta]_+
+C_2\ell^2[-\theta]_+^2\right).
\]

Here ell is an explicitly defined local length. This is a schematic dimensional form, not a transcription of all coefficients in [F]. The force and pressure work must be derived together.

At a smooth expanding steady solution q can vanish. During compression, reversal, acoustic oscillation, or shock formation it can remain active. Measure its heating and force contributions and their refinement dependence; "the reference paper used 0.75" is not an acceptance test.

### 7.3 Hydrostatic balance is a first-order design requirement

The continuous equilibrium relation is dp/dr = -rho dPhi/dr. Discretely, a small mismatch between two large terms can drive a velocity larger than the wind being studied.

Staggering changes where that mismatch occurs but does not eliminate it automatically. Construct a discrete equilibrium on the proposed primary and momentum grids, and evolve departures from the same equilibrium used in the source and pressure operators. Require a rest state to remain at rest to the designed discrete tolerance.

### 7.4 Implicit does not mean nonoscillatory

Backward Euler is a useful initial pressure-relaxation prototype because it damps unresolved stiff linear modes. That damping is also an accuracy cost. A second-order method should be introduced only with a demonstrated temporal-convergence test.

Crank-Nicolson is not L-stable: very stiff modes need not decay rapidly and can alternate in sign. The existing EXHALE transport routine uses a Crank-Nicolson-type tridiagonal update with admissibility checks. Reusing its linear solver is reasonable; assuming its time integrator is automatically the right pressure integrator is not.

### 7.5 Splitting can create a false appearance of steady convergence

If the split step map is S_dt, a fixed point S_dt(U) = U need not be a zero of the intended unsplit residual at finite dt. This matters especially when pressure work, gravity, and prescribed heating nearly cancel.

Check both the step defect and a correctly defined balance residual. For a split map, repeat at smaller dt and establish that the stationary profiles and fluxes approach a common limit. Local pseudo-time intervals make this distinction still more important.

### 7.6 Stability of the atmosphere itself is a separate question

An explicit physical-time method need not converge to an unstable steady root that a Newton solver can locate. Conversely, extra numerical damping can make an unstable state look attractive.

If one solver oscillates and another converges, compare the physical-time response to controlled perturbations and track how growth rates change with grid size, time step, and artificial viscosity. Do not identify convergence alone with physical correctness.

### 7.7 A stability qualification in the existing EXHALE comparison

The active code in contact_mode_dissipation_flux multiplies a third velocity difference by a spatially varying positive coefficient, then supplies the paired energy flux v_face times the momentum correction. That pairing supports total-energy accounting, but does not by itself make the correction dissipative in kinetic energy.

At frozen density, the relevant interior quadratic form has the structure

\[
\dot K_{\rm art}=d^TWLd,
\]

where d contains neighboring velocity differences, L is the negative second-difference operator, and W contains positive spherical face weights and damping coefficients. For constant W this is nonpositive; for variable W, the symmetric part of WL need not be negative semidefinite.

For example, the abstract two-component choice L = [[-2,1],[1,-2]], W = diag(1,100), and d = (10,1) gives d^TWLd = 610, a positive value. This is an algebraic counterexample to unconditional dissipativity, not a measured EXHALE atmosphere. The source file itself discusses this limitation; the implemented loop still has the variable-coefficient third-difference form.

Therefore, benchmark this optional stress separately from the undamped operator. If it remains necessary, derive a conservative stress whose volume-weighted energy form is nonpositive for variable coefficients and stretched geometry, with matching total-energy work. Do not infer that either CIP's removal of a numerical oscillation or the current stress's smoothing of a profile proves a correct energy transfer. The existing stress tests were not rerun in this review.

## 8. Published formulas that should not be copied literally

The following are defects or ambiguities in the displayed equations of [F], confirmed against journal pages 28-29. **They are not proof of defects in the authors' executable code**, which was not available for this review.

| Location | Problem in the printed expression | Required treatment |
|---|---|---|
| A9 | The density source uses a velocity difference divided by r_i squared, without the face-radius-squared factors needed to discretize A2 | Differentiate r squared times v, with consistent face coordinates |
| A10 | The displayed planetary gravity is outward; the stellar term is written with a/r_i cubed instead of r/a cubed as in A3 | Derive acceleration from the selected potential; for the paper's stated approximation it is -GM_p/r squared + 3GM_star r/a cubed |
| A4 and A11 | The specific-energy conduction term lacks the density and spherical-divergence normalization required for a conductivity; A11 also adds rate-like terms to an energy without the corresponding displayed time factors | Derive the energy update from units and the governing equation, not these displays |
| A12 | The left derivative is printed with respect to r where the nonadvection evolution requires time | Use the gradient evolution equation consistently |
| A16-A18 | Spatial source/stretch derivatives contain a denominator 2 Delta-t where a spatial separation is required | Re-derive from A6 and use actual nonuniform coordinates |
| A27 | The coefficient contains Y_(i+1) minus 2Y_i, whereas the Hermite constraints require Y_(i+1) plus 2Y_i in that bracket | Use coefficients satisfying both endpoint values and gradients |
| Section 2.5 versus Appendix A.2 | The order of advection and nonadvection steps is described differently | Specify the implemented composition and verify its temporal order |

For signed interval h from one interpolation node to its neighbor, consistent cubic Hermite coefficients are

\[
F(\xi)=a\xi^3+b\xi^2+g_i\xi+f_i,
\]

\[
a=\frac{g_i+g_k}{h^2}
-\frac{2(f_k-f_i)}{h^3},\qquad
b=\frac{3(f_k-f_i)}{h^2}
-\frac{2g_i+g_k}{h}.
\]

The departure point determines the appropriate interval; it is not always the interval on the increasing-radius side.

The stretching term also deserves an implementation note. Appendix A discusses it in the nonadvection gradient equation and again when differentiating the advective equation. A consistent algorithm can allocate it to one part of the split or account for it through a differentiated characteristic map. Applying the full contribution in both places would double count it. The published description alone is insufficient to establish which treatment the authors' code uses.

These findings do not invalidate the CIP concept or the paper's numerical results. They make a literal equation-by-equation transcription unsafe.

## 9. A physically consistent adoption design

### 9.1 Fix the mathematical problem before changing the algorithm

For an initial hydrodynamics-only comparison, use a specified ideal-gas closure with constant gamma, a fixed potential Phi(r), and prescribed heating and conductivity. This removes changes in other physical models as a confounding factor.

The radial equations to be reproduced are

\[
\partial_t\rho+\frac{1}{r^2}\partial_r(r^2\rho v)=0,
\]

\[
\partial_t(\rho v)
+\frac{1}{r^2}\partial_r[r^2(\rho v^2+p)]
=\frac{2p}{r}-\rho\Phi',
\]

\[
\partial_tE+\frac{1}{r^2}\partial_r[r^2v(E+p)]
=-\rho v\Phi'+H-\frac{1}{r^2}\partial_r(r^2q_{\rm th}),
\qquad q_{\rm th}=-\kappa\partial_rT.
\]

Viscosity, when enabled, requires its matched momentum stress and energy work. Do not compare a viscous calculation with an inviscid one and attribute the difference solely to transport interpolation.

### 9.2 Use the correct conserved measure for CSL2

Two mathematically consistent choices for continuity are:

1. In radial coordinates, transport w = r squared times rho through
   \[
   \partial_t w+\partial_r(vw)=0.
   \]
   The conserved moment is the integral of w over dr.
2. In the volume coordinate x = r cubed/3, transport rho through
   \[
   \partial_t\rho+\partial_x(a\rho)=0,\qquad a=r^2v.
   \]
   The conserved moment is the integral of rho over dx.

The second choice makes the moment equal to the existing spherical cell mass without a geometric quadrature factor. It does not make the radial momentum and gravity equations Cartesian; those must be transformed separately.

My preference for the first experiment is a radial-coordinate formulation with exact spherical mass moments, because EXHALE's geometry and gravity are already expressed in r. A volume-coordinate prototype is also valid. Choose one and derive all point and moment updates consistently rather than mixing their transport velocities.

Do not conserve the unweighted integral of rho over dr and call it spherical mass conservation.

### 9.3 Compute shared, time-integrated transfers

For each face, compute one swept mass transfer over the accepted interval and use it with opposite signs in its two adjacent cells. For a large Courant number, include every fully crossed cell and the partial-cell contribution, with an appropriate velocity trajectory.

An instantaneous flux and a time-integrated transfer are different quantities. Existing consumers of face_flux must not receive a transfer without a documented conversion and time-level convention.

A passive tracer with constant initial fraction should remain constant when it uses that mass transfer. This is a transport-consistency test, not an analysis of atomic or molecular physics.

### 9.4 Decide the momentum and energy contract explicitly

For a reference reproduction, the Kuramoto-style mass-only CSL2 approach is a reasonable comparison, provided momentum and energy defects are measured and reported.

For a production method, I recommend a conservative remap of momentum and total energy, followed by compatible pressure work and gravity. The final energy budget should close for E + rho Phi when the potential is static, including boundary fluxes and external energy sources.

On a staggered grid, define:

- The control volume of each velocity or momentum degree of freedom.
- The density used in face acceleration.
- The discrete kinetic energy.
- The pressure-work pairing between mass cells and momentum volumes.
- The conversion back to the public cell-averaged state.

The pressure gradient and velocity divergence should satisfy an appropriate discrete integration-by-parts relation on their volume weights. That is a stronger criterion than using apparently symmetric averages.

### 9.5 Add pressure implicitness as a separate experiment

For an illustrative frozen-coefficient acoustic substep without gravity or other sources,

\[
v^{n+1}=v^*-\Delta t\,\frac{1}{\rho^*}\partial_r p^{n+1},
\]

\[
p^{n+1}-p^*=-K^*\Delta t\,D v^{n+1},
\qquad
K^*=\rho^*c_{\rm ad}^{*2},\quad
Dv=\frac{1}{r^2}\partial_r(r^2v).
\]

Elimination gives

\[
\frac{p^{n+1}-p^*}{K^*\Delta t^2}
-D\left(\frac{1}{\rho^*}\partial_rp^{n+1}\right)
=-\frac{Dv^*}{\Delta t}.
\]

This is a schematic acoustic linearization related to [Y01J], Section 4.2, not a complete proposed EXHALE equation. Gravity, nonlinear coefficients, boundary relations, energy consistency, and spherical weights must still be supplied. In a stratified atmosphere, solve for a pressure departure about a matching hydrostatic state so that implicit pressure and gravity do not generate a new imbalance.

**Avoid double-counting compression.** A conservative CSL2 continuity update already includes compression through the divergence of the mass flux. Appending the primitive CCUP density change to that completed update can count compression twice. Either derive a conservative advection/acoustic split with complementary fluxes or compute a pressure correction to the mass transfer. Apply the same accounting to pressure work in energy.

A useful alternative is to build the pressure correction around EXHALE's existing conservative state and use CIP-CSL2 only if it independently improves transport accuracy. This separates the two scientific questions and can avoid a premature full layout rewrite.

### 9.6 Preserve physical boundary semantics, not just ghost-array sizes

The current base routine uses a reservoir and an outgoing acoustic relation. It selects reservoir density at rest/inflow and an interior isentropic density on reversal, with a wind-window criterion or matched face velocity determining the branch. It also contains configurable acoustic matching and a velocity guard.

That is not identical to Frelikh's fixed-density/fixed-pressure boundary with extrapolated velocity. A method comparison should implement the same physical boundary model in both layouts before exploring a different one.

At an outer boundary:

- Supersonic outward flow has no incoming Euler characteristic.
- Subsonic outward flow has one incoming acoustic characteristic and needs an exterior closure.
- Inward flow changes the characteristic count again.
- Thermal conduction needs a thermal boundary condition even when Euler waves all leave.

The current free_outflow_ghost routine actually imposes an isothermal hydrostatic continuation of density and pressure and copies velocity. It has no explicit switch implementing a separately specified exterior acoustic condition. Its introductory comment that no characteristic enters an outflowing face is valid only for supersonic outflow. For subsonic use, this continuation is an implicit numerical boundary closure whose reflection and physical interpretation require testing. Neither keeping it silently nor replacing it with a paper's supersonic Neumann condition resolves that issue.

Store and reconstruct gradient/moment boundary information from the chosen physical closure. Independently extrapolating those fields can contradict the primitive face state.

### 9.7 Keep the experiment isolated and restartable

Use a distinct solver selection and state type. Do not overload a reconstruction label to mean a different arrangement of variables, time integrator, and conservation law.

Recommended responsibilities, named for their calculations, include spherical_mass_transport, staggered_pressure_acceleration, pressure_work, and acoustic_pressure_correction.

Exact restart of an original-CIP method requires the independent gradients, and CSL methods require the independent moments. A file containing only physical cell averages can instead support a documented warm restart that reinitializes the additional information. These are different guarantees.

Rejected steps must restore all gradients, moments, boundary data, diagnostic accumulators, and pressure-solver state that influence the next evaluation. A rollback of only rho, momentum, and energy is insufficient.

Changing public state layout or restart interfaces should follow a separate design decision before implementation. This report proposes that design; it does not implement or authorize it.

## 10. Staged implementation and acceptance plan

### Stage 0: establish a measured baseline

Freeze the test configurations and identify the bottleneck. Measure the existing explicit and stationary routes with their actual flux, reconstruction, well-balanced, damping, boundary, and time-step settings.

Use the current conservation-budget export to inspect the spatial rows, but do not mistake an instantaneous stationary-row export for a time-integrated budget of an accepted physical step. Add the latter to the experiment if it is needed.

Deliverable: baseline errors, budgets, time-step restrictions, and runtime breakdown. No new hydro method is needed for this stage.

### Stage 1: verify interpolation and spherical transport

Implement a small separate CIP/CSL2 test program, not a replacement in the main atmospheric loop.

Test:

- Constant and linear profiles, including both velocity signs.
- Smooth transport on uniform and stretched grids.
- A nonuniform velocity field with a known characteristic map.
- Transport across multiple cells in one interval.
- Discontinuous positive profiles and their limiting.
- Spherical mass conservation and a constant passive fraction.

Acceptance: polynomial constraints hold, integral updates telescope to boundary transfers, no unreported negative state or mass repair occurs, and refinement establishes the claimed accuracy. Store test code and results under docs/audit_20260905.

### Stage 2: construct an explicit hydrodynamic comparison

Implement the pressure, gravity, energy, and boundary updates with a common global interval. Start with fixed ideal-gas parameters and prescribed forcing.

Compare original CIP and mass-conservative CSL2 variants where feasible. This identifies whether the benefit comes from moment conservation, staggering, or interpolation. If only the hybrid is implemented, do not infer what the untested variants would have done.

Acceptance: the tests below pass before any production atmosphere is used as evidence.

### Stage 3: test pressure implicitness only if warranted

Add a pressure-implicit version if Stage 0 or Stage 2 shows that acoustic steps dominate the cost. Begin with a robust first-order acoustic update, then establish the accuracy of a higher-order composition separately.

Acceptance: larger steps remain physically admissible; slow-flow accuracy and budgets are preserved; pressure iterations converge; and wall time decreases at a matched error. Merely surviving a large time step is not a pass.

### Stage 4: compare with EXHALE at fixed physical inputs

Use identical domain, potential, thermodynamic service, forcing, and physical boundaries. Keep the excluded physical models unchanged. Compare profiles, sonic structure, mass loss, boundary fluxes, and integrated budgets under spatial and temporal refinement.

A CIP state can also be transferred as an initial guess to the current steady solver. Measure whether it reduces total time to the existing certificate, including the time spent generating that initial state.

Acceptance: a new discretization has its own named certificate; an initializer succeeds only when the existing target operator is subsequently solved.

### Stage 5: production decision

Promote the method only if it offers a reproducible benefit at a specified accuracy and no unresolved conservation or boundary defect. Retain a comparison route until independent benchmarks cover its intended regime.

Stop development of a particular variant if its energy or boundary error cannot be controlled by refinement, if positivity requires persistent nonconservative repairs, or if its measured advantage disappears at matched accuracy. A different method may still be useful as a diagnostic even when it is not a suitable production replacement.

### Benchmark matrix

| Test | What it isolates | Measurements |
|---|---|---|
| Uniform stationary state | Consistency and unintended source terms | Drift, positivity, conservation |
| Isothermal and adiabatic hydrostatic columns | Pressure-gravity and boundary balance | Spurious velocity, force residual, boundary drift |
| Weak smooth flow over a Mach-number sequence | Low-Mach transport and stiffness | Flux error, profile error, accepted interval, runtime |
| Stationary and moving contacts | Diffusion, overshoot, pressure consistency | Width, extrema, pressure error |
| Linear acoustic pulse | Stability versus accuracy | Phase error, damping, reflection |
| Shock tube and strong expansion | Shock energetics and positivity | Jump conditions, shock speed, total-energy error |
| Isothermal Parker wind | Transonic selection and small basal velocity | Mass-loss error, sonic radius, resolution dependence |
| Spherical wind with prescribed smooth heating | Full hydrodynamic energy balance | Boundary enthalpy/kinetic/potential fluxes and source integral |
| Conduction with prescribed thermal boundaries | Diffusive coupling | Analytic or manufactured solution, thermal budget |
| Stagnation and flow reversal | Branch changes and characteristic closure | Admissibility, reflection, flux continuity in time |
| Extended outer domain | Dependence on the exterior closure | Interior profiles, incoming waves, mass loss |
| Interrupted and rejected updates | State completeness | Restart and rollback reproducibility |

For a first comparison, 250, 500, and 1000 cells provide a useful proposed refinement sequence, not a guarantee that the finest grid is sufficient. Keep the physical domain fixed and reduce local spacing consistently; increasing the number of outer cells while leaving the base resolution unchanged is not a base-layer convergence test.

Use multiple time steps on each useful grid. For isothermal Parker tests, energy is deliberately constrained by the isothermal model, so they cannot validate the unconstrained total-energy equation. Include a separate energy test.

Conservation tolerances should distinguish roundoff-level algebraic cancellation, pressure-solve tolerance, quadrature error, spatial truncation, and time splitting. Near zero mass flux, use absolute as well as relative errors. A paper's 2% or 3% result is a reference point, not a universal EXHALE acceptance threshold.

## 11. Checks completed in this review

The independent script [check_published_formulas.py](audit_20260905/cip_hydrodynamics_20260921/check_published_formulas.py) was executed successfully with Python 3. It tests normalized algebraic examples, not an atmospheric simulation.

| Check | Literal printed expression | Consistent expression |
|---|---|---|
| A27, linear profile with h = 1, endpoint values 2 and 3, both slopes 1 | Endpoint value error 4; endpoint gradient error 8 | Both errors zero |
| A9, constant radial velocity 0.5 at radius 2 with density 1 | Density source zero | Density source -0.5 |
| A10, planetary gravity alone with GM = 1 at radius 2 | Acceleration +0.25 | Acceleration -0.25 |

These directly measured results confirm that the selected printed formulas cannot be implemented literally. They do not measure the accuracy, speed, or stability of the authors' code or EXHALE.

Other completed checks:

- Read the active main-loop calls, conservative conversions, reconstruction/flux/source paths, boundaries, and stationary selection.
- Checked the listed implementation modules in the Makefile.
- Read the paper descriptions of explicit CIP, mass-only CSL2, and CCUP rather than treating them as one method.
- Checked the dimensions and spherical geometry of the relevant displayed equations.
- Verified publication metadata through ADS.

Not performed: compilation of a new solver, new EXHALE runtime experiments, code-to-code benchmarks, full regressions, authors' source-code verification, or analysis of atomic/molecular physics. A full regression was not run because no production code was changed and it would not establish the accuracy of a solver that has not yet been implemented.

## 12. Final judgment

The papers provide a credible reason to test a different hydrodynamic discretization, particularly for weak winds on steeply stratified grids. They do not show that explicit original CIP is a universally faster, more stable, or more conservative replacement for the present EXHALE solver.

The best near-term experiment is **spherical CIP-CSL2 transport plus a carefully matched staggered hydrodynamic test**, with conservation and hydrostatic balance measured from the beginning.

If runtime is dominated by acoustic marching, **pressure implicitness is the more directly relevant improvement**, and it should be evaluated independently of the choice of CIP interpolation. For production, preserve conservative mass and total energy, derive compatible pressure work and gravity, and state the boundary physics explicitly.

The main practical warning is simple: **do not copy Frelikh's displayed Appendix A equations directly, and do not identify Kuramoto's mass-conservative hybrid with an all-equation conservative implicit solver.**

## References

[F]: ../../references/Frelikh_2026_ApJ_996_96.pdf
[Y91]: ../../references/Yabe_1991CPC_66_219.pdf
[Y01M]: ../../references/Yabe_2001MWRv_129_332.pdf
[Y01J]: ../../references/Yabe_2001JCP_169_556.pdf
[K13]: ../../references/Kuramoto_2013EPSL_375_312.pdf

- [Frelikh and Murray-Clay, published PDF][F]. 2026, ApJ, 996, 96. DOI: [10.3847/1538-4357/ae18d5](https://doi.org/10.3847/1538-4357/ae18d5).
- [Yabe and Aoki, published PDF][Y91]. 1991, Computer Physics Communications, 66, 219-232. DOI: [10.1016/0010-4655(91)90071-R](<https://doi.org/10.1016/0010-4655(91)90071-R>).
- [Yabe, Tanaka, Nakamura, and Xiao, published PDF][Y01M]. 2001, Monthly Weather Review, 129, 332-344. DOI: 10.1175/1520-0493(2001)129&lt;0332:AECSLS&gt;2.0.CO;2. [ADS publication record](https://ui.adsabs.harvard.edu/abs/2001MWRv..129..332Y/abstract).
- [Yabe, Xiao, and Utsumi, published PDF][Y01J]. 2001, Journal of Computational Physics, 169, 556-593. DOI: [10.1006/jcph.2000.6625](https://doi.org/10.1006/jcph.2000.6625).
- [Kuramoto, Umemoto, and Ishiwatari, published PDF][K13]. 2013, Earth and Planetary Science Letters, 375, 312-318. DOI: [10.1016/j.epsl.2013.05.050](https://doi.org/10.1016/j.epsl.2013.05.050).
