# Review of the solver-approach analysis dated September 10, 2026

Review date: September 10, 2026  
Reviewed source revision: `a0f4292b39e68cb83c1578a63cf4d50bd081242f`  
Scope: [solver_approach_analysis_20260910.md](solver_approach_analysis_20260910.md), the current implementation, relevant experiment records, and published numerical methods.

## 1. Overall assessment

The document makes a reasonable practical recommendation to stop extending the current species-row JFNK experiment without a new hypothesis. It also correctly identifies time integration and partitioned iteration as serious alternatives. Its account of the historical CETIMB time-integration method is substantially supported by the published papers.

However, its central mathematical explanation is too strong and partly incorrect. The experiments establish failure of particular implementations, initial states, preconditioners, and iteration budgets. They do not establish that a coupled stationary root cannot be obtained, that its Newton basin does not exist, or that spatial discretization error prevents a smaller algebraic residual. Several proposed changes already have substantial implementations, with limitations that the analysis omits.

My recommendation is to retain the experimental record and the pause, but revise the diagnosis and decision criteria before using this document to select a new architecture or relax acceptance requirements. In particular, do not replace the existing wind species gate with a report-only outcome on the strength of these arguments.

### Principal findings

| Priority | Finding | Consequence |
| --- | --- | --- |
| Critical | Spatial truncation error is treated as a lower bound on the residual of the discrete equations. | The proposed tolerance justification and the main method-selection test are invalid. |
| High | Ratios of Ritz-value magnitudes, including masked compressions, are interpreted as conditioning and near-singularity of the coupled physical system. | The measurements support a difficult current linear solve, not the claimed impossibility result. |
| High | Partitioning is described as removing coupling and requiring a new outer loop. | Coupling remains in the outer iteration; a callable outer loop and relaxation routines already exist. |
| High | Removing species from Newton followed by a hydrodynamic finish is presented as the CETIMB method. | That finish is an EXHALE proposal, not established by the cited papers, and can invalidate the species balance. |
| High | The well-balanced option is said to remove the rounding amplifier, although N37 reports that the measured floor did not fall. | The summary contradicts its own experiment record and overlooks an active normalization defect. |
| Medium | Similarity and matched column scaling are generalized into claims that scaling cannot affect convergence. | Only narrower algebraic statements are valid. |
| Medium | Success of isolated components and a small source contribution are treated as evidence that the complete molecular transport problem is sound. | Chemical stiffness, admissibility, closure, and coupled boundary conditions remain separate questions. |

## 2. Evidence and limits of this review

Three kinds of evidence are distinguished below:

- **Code-verified:** the current routine body and its relevant callers were inspected. This is not a claim that a fresh production run passed.
- **Recorded:** a number was read from an existing experiment record. The expensive atmosphere experiments were not repeated in this review.
- **Executed here:** small, independent algebraic examples were run to check implications used in the analysis. These do not measure EXHALE's convergence.

The principal code paths inspected were:

| Implementation | Relevant location and purpose |
| --- | --- |
| [EXHALE_main.f90](../src/EXHALE_main.f90) | Calls at lines 1487, 3765, and 5214; `steady_wind_with_element_diffusion`, lines 5989–6235. |
| [steady_newton.f90](../src/modules/time_step/steady_newton.f90) | Species registry, lines 1866–1946; returned-state certification, lines 5116–5144 and 5546–5612; masked Arnoldi diagnostics, lines 11991–12239; trust-region selection, line 14855; final acceptance, lines 16048–16066. |
| [binary_element_diffusion.f90](../src/modules/functions/binary_element_diffusion.f90) | `relax_element_composition`, lines 815–1021; fixed-background transport and composition relaxation. |
| [diffusive_photochemistry.f90](../src/modules/lower_atmosphere/diffusive_photochemistry.f90) | `relax_photochemical_composition`, lines 4718–4863; bounded carrier advancement, not an unconditional stationary solve. |
| [ionization_equilibrium.f90](../src/modules/radiation/ionization_equilibrium.f90) | Transported molecular and proton fractions imposed on the chemical solve, around lines 1885–1917. |
| [input_read.f90](../src/modules/files_IO/input_read.f90) | Configuration restrictions on ionization transport and uncoupled Newton, lines 2094–2175. |
| [steady_residual.f90](../src/modules/time_step/steady_residual.f90) | `store_row_terms`, lines 280–302; momentum normalization. |
| [Source.f90](../src/modules/states/Source.f90) and [RK_rhs.f90](../src/modules/time_step/RK_rhs.f90) | The well-balanced source cancellation and momentum assembly. |
| [certification.f90](../src/modules/time_step/certification.f90) | Wind tolerances, species gates, and independent transport-row evaluation. |
| [certification_contexts.f90](../src/tests/certification/certification_contexts.f90) | Synthetic-column relaxation and restriction to coarser grids, especially lines 860–965. |

The build source list includes the relevant modules, and the main-program calls establish that the outer iteration is not merely an unused design sketch. I did not rebuild or inspect the symbol tables of a particular executable; its build provenance should be checked before any new production measurement.

Supporting records include [the tolerance anchoring memo](certification_tolerance_anchoring_20260910.md), [decisions 22 and 23](To_be_determined_by_user_20260906.md), and N35–N38 in [Update_EXHALE.md](Update_EXHALE.md). Their reported run results are not relabeled as new measurements here.

No production source, configuration, restart, or reference output was changed. No atmosphere was advanced. The existing instruction to stop after N37/N38 was respected.

## 3. The most important correction: two different errors

### 3.1 A discrete root can have a residual far below the spatial error

Let the continuous stationary equations be

\[
\mathcal F(U)=0,
\]

and let the implemented spatial operator on a fixed mesh be

\[
F_h(U_h)=0.
\]

Three quantities must not be conflated:

1. The algebraic residual of the current iterate, \(F_h(U_h^{(k)})\).
2. The truncation defect obtained by inserting a sampled continuous solution, \(F_h(I_hU)\).
3. The error in the converged numerical solution, \(U_h-I_hU\).

The first measures how accurately the discrete equations have been solved. The second and third characterize the discretization. A spatial error of order \(h^p\) does not prevent the first quantity from approaching floating-point accuracy. Indeed, algebraic error normally must be sufficiently smaller than discretization error to make a grid-convergence study interpretable.

It is reasonable to stop an iterative solve when its remaining error is demonstrably negligible for the intended physical accuracy. That is an efficiency decision supported by error estimates. It is not a mathematical lower bound on the residual, and two dimensionless numbers are not automatically comparable when their normalizations differ.

This distinction invalidates the statements in target sections 3, 5.1, and 7 that a roughly `5e-4` spatial error at `N = 500` prevents `1e-5` certification of the discrete equations. It also invalidates using the same argument to establish that the coupled system has no suitable Newton basin.

### 3.2 The existing anchoring experiment does not measure such a floor

The synthetic-column test first relaxes a fine-grid column, then selects its values on coarser nested grids and evaluates the coarse operators. The code does not independently solve the stationary problem on each coarse grid before computing those residuals.

The anchoring memo records a fine-grid residual of approximately `7.90e-13`, followed by coarse-grid residuals of approximately `2.47e-6` and `7.43e-6`, with an inferred order near `1.59`. Those are useful consistency measurements on a restricted profile. They are not the minimum residuals attainable by solving the coarse discrete equations. The fine-grid result itself demonstrates that a nonuniform species column can satisfy its discrete balance much more accurately than the quoted spatial-error estimate.

The candidate remapping experiment has a related limitation. Evaluating an unconverged candidate after remapping from `N = 500` to `N = 1000`, and applying a Richardson-style formula with an assumed order, does not establish that a stated fraction of the original residual is unavoidable discretization error. A nearly exact round trip through interpolation does not bound the fine-grid derivative error: interpolation and restriction can cancel in a round trip while still changing gradients on the intermediate grid.

The estimated order from a smooth synthetic column also cannot be transferred without testing to the molecular front, a different limiter regime, or a different row normalization. In the anchoring memo, several carrier entries are explicitly unmeasured. An elemental estimate therefore does not justify a numerical floor for every carrier equation.

**Required revision:** describe the remapping results as sensitivity or consistency diagnostics. Reserve a spatial-accuracy claim for independently converged meshes, or for a manufactured-solution experiment with known forcing and a controlled norm.

### 3.3 Executed counterexample

The archived script solves the centered discrete Poisson equation analytically for a sine forcing and evaluates the residual in double precision. Its directly executed output was:

| Interior cells | Residual normalized by the forcing amplitude | Maximum solution error against the continuous sine |
| --- | --- | --- |
| 31 | `4.573807e-14` | `8.035777e-4` |
| 63 | `2.031103e-13` | `2.008218e-4` |
| 127 | `7.540370e-13` | `5.020092e-5` |

This is a counterexample to the proposed general implication, not a prediction of the residual EXHALE can attain.

### 3.4 Rounding, derivative error, and acceptance also need separation

Repeated evaluation, input quantization, cancellation in flux differences, approximate chemical elimination, finite-difference derivative error, and spatial truncation error are different quantities. They need different tests.

For example, a finite-difference action has the schematic error

\[
\frac{\widetilde F(Y+\epsilon v)-\widetilde F(Y)}{\epsilon}
=Jv+O(\epsilon)+O(\delta_F/\epsilon).
\]

This explains why both overly small and overly large probes can be inaccurate. It does not, without further analysis, identify a lower bound on the final residual. Likewise, a relative mismatch between a banded approximate Jacobian and another derivative estimate is not itself an absolute residual floor. A consistently imperfect derivative can permit repeated linear convergence; in this code the band also serves as a preconditioner, not as the exact matrix defining every Krylov product.

N31–N34 provide substantial evidence that the current matrix-free actions become unreliable on particular directions. That is important evidence against the current derivative/preconditioner combination. It does not prove that every derivative formulation or nonlinear method must have the same limitation.

Decision 22 has already selected `1e-5` wind gates and report-only layer species rows. This review does not change that decision or recommend silently restoring another value. It recommends correcting its numerical justification and keeping the existing decision explicit until a separate, properly supported acceptance decision is made.

## 4. What the Krylov diagnostics establish, and what they do not

### 4.1 The reported ratios are not measured condition numbers

The diagnostic performs Arnoldi iteration and calls `dgeev` on the resulting Hessenberg matrix. It reports the largest and smallest Ritz-value magnitudes and their ratio. For the species and hydrodynamic measurements, a mask is applied to the initial vector and again after every operator action.

Thus the species diagnostic samples a compression of the preconditioned operator onto the species subspace. It is not an explicit singular-value decomposition of the full coupled Jacobian. It is also not the Schur complement obtained by eliminating the hydrodynamic variables.

For a generally nonnormal matrix,

\[
\frac{\max |\lambda|}{\min |\lambda|}
\ne \kappa_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.
\]

Moreover, finite-dimensional Ritz approximations need convergence checks before being identified with extreme eigenvalues of the full operator. Direction-dependent finite-difference errors further weaken a literal spectral interpretation.

The N35 records are credible evidence that the current preconditioned action has difficult species-associated directions. They are not sufficient evidence for the statement that the original coupled physical problem is near-singular in the claimed sense.

### 4.2 Compression is not elimination

Write a coupled linearization as

\[
J=\begin{pmatrix} A&B\\ C&D\end{pmatrix}.
\]

If the hydrodynamic variables are eliminated, the relevant species operator is

\[
S=D-CA^{-1}B,
\]

not simply the species principal block or a masked compression of a right-preconditioned operator.

For a small exact counterexample,

\[
J=\begin{pmatrix}1&1\\1&0\end{pmatrix}
\]

has a zero species principal block but determinant `-1` and a 2-norm condition number of approximately `2.618`. Its species Schur complement is `-1`. The archived script evaluates these quantities. This does not claim that EXHALE's operator is well conditioned; it shows why the diagnostic needs a more limited interpretation.

### 4.3 Small eigenvalues and one failed GMRES cycle do not prove impossibility

The comment above `ritz_values_of_the_preconditioned_operator` also overstates the relationship between the eigenvalue spectrum and GMRES. For a diagonalizable nonnormal matrix, bounds generally involve the eigenvector conditioning as well as the polynomial on the spectrum. More generally, nonnormality and the right-hand side matter. A small eigenvalue alone does not force an arbitrarily long Krylov iteration: a system with only two distinct nonzero eigenvalues can be solved in at most two exact-arithmetic GMRES steps despite a very large eigenvalue-magnitude ratio.

A best step found in a particular finite Krylov basis is the best step in that sampled space, under its measured action and objective. It is not the best step that the operator admits among all directions. The target's section 2 should make that restriction explicit.

Similarly, N35's small sampled differences between full and band actions are encouraging for local stencil coverage. They do not establish that the band misses nothing. For right preconditioning,

\[
JM^{-1}=I+(J-M)M^{-1}.
\]

A small difference in sampled columns need not be small after multiplication by a poorly conditioned inverse, nor on untested directions. Constraint projections, chemical elimination, and radiation responses should be included in identifying the operator actually tested.

### 4.4 The scaling statements require qualifications

There is a valid narrow statement about matched column scaling. In exact arithmetic, if the same invertible scaling `E` is applied consistently to the operator and preconditioner,

\[
(JE)(ME)^{-1}=JM^{-1}.
\]

That explains why a particular matched column-scaling experiment may leave the right-preconditioned cycle unchanged. It does not cover changes to variables that also change constraints, finite-difference probes, admissibility limits, globalization, or the construction accuracy of the preconditioner.

For matched row scaling `S`, the resulting preconditioned operator is similar:

\[
(SJ)(SM)^{-1}=SJM^{-1}S^{-1}.
\]

Similarity preserves eigenvalues, not the Euclidean residual norm, singular values, or finite-step GMRES behavior. The executed example `[[1,a],[0,2]]` with the corresponding right-hand side `(0,1)` gives a one-step relative residual of `0.999800060` for `a = 100`, but `0.447213595` for the diagonally similar matrix with `a = 1`.

The correct conclusion is that the particular N36/N38 choices did not help under their tested conditions. This is not a recommendation to resume unguided scaling experiments. It is a correction to the claim that their failure follows generally from similarity.

### 4.5 Last-bit sensitivity does not locate or exclude a Newton basin

Sensitivity can arise from inaccurate actions, ill-conditioned solves, limiter transitions, constraints, acceptance branches, or accumulated differences outside a local convergence region. It does not identify which mechanism dominates, prove physical chaos, or prove nonexistence of a smooth root nearby.

Likewise, an improved initial state after 2000 marching steps shows that transient evolution or continuation was useful in that experiment. It does not distinguish a poor preconditioner from nonlinear front motion, an inaccurate derivative, a distant root, or an incompatible boundary specification.

## 5. The proposed segregated route already has substantial code

### 5.1 What exists and is invoked

`steady_wind_with_element_diffusion` already alternates a stationary wind solve with composition updates. It reads an outer-pass limit, element under-relaxation, and a carrier trust setting. It calls:

1. `solve_steady_jfnk` or `solve_steady_ptc`;
2. the state and chemical refresh;
3. `relax_element_composition`, when element diffusion is enabled;
4. `relax_photochemical_composition`, when carriers are transported outside the coupled Newton system;
5. a carrier residual diagnostic and outer stopping conditions.

The element routine holds the hydrodynamic background and the supplied face mass flux fixed, advances composition on a transport time scale, and updates its coefficients. This is substantial reuse for section 5.1, not a missing frozen-background concept.

The carrier routine is different: its implementation exits when a composition-movement threshold is exceeded, when a small incremental change is observed, or when its step budget is exhausted. It is a bounded advancement procedure. Its existence does not establish an available solver that returns a certified stationary molecular state on every frozen background.

The historical behavior discussed in that routine's comments includes saturation at an elemental ceiling and difficulty with under-relaxation. Those numerical values were not remeasured here, but the current bounded-advance logic is code-verified. The target should address this existing experience before recommending a fresh full stationary carrier relaxation.

**Required revision:** describe section 5.1 as redesigning or validating the existing alternation, with a stricter equation and acceptance contract. Do not estimate its work as writing a previously absent outer loop.

### 5.2 Partitioning moves coupling; it does not remove it

Let the two nonlinear blocks be `H(u,c) = 0` and `G(u,c) = 0`. For exact successive block solves, linearization near a common root gives

\[
e_c^{k+1}\simeq D^{-1}CA^{-1}B\,e_c^k.
\]

Under-relaxation changes this iteration matrix to

\[
(1-\omega)I+\omega D^{-1}CA^{-1}B.
\]

The coupling remains in this matrix. Convergence is not guaranteed by successful individual block solves. Nor does positive under-relaxation stabilize every unstable mode. The archived example `J = [[1,2],[2,1]]` has trivial scalar subsolves but an undamped composition error multiplier of `4`; damping gives `1 + 3 omega`, still greater than one for every positive `omega`.

Therefore the comparison-table entry should be changed from “removes the failing coupling” to “moves the coupled iteration outside the Krylov solve.” Local linear convergence is conditional, not an established property of the proposed route.

### 5.3 Current stopping logic is not a joint certification contract

The outer loop uses an element drift threshold and, for carriers, a volume-weighted residual condition. It modifies composition after the hydrodynamic solve. At a successful outer exit, it does not itself perform another complete hydrodynamic residual evaluation after all composition changes.

Conversely, `stationary_rows_of_the_returned_state` intentionally checks the equations represented in the current Newton registry. With the registry empty, species entries are evaluated and reported by certification but do not decide that solver's success flag. A hydrodynamic `info = 0` therefore cannot be used as evidence that the outer composition problem has converged.

This is not evidence that the final writer necessarily hides the failure; the certification report can still expose it. It is a verified mismatch between the existing subsolve-success contract and the joint-success contract section 5.1 needs.

A revised outer iteration needs:

- A final evaluation of all active hydrodynamic, transport, chemical-closure, and conservation conditions on exactly the returned state.
- Residual magnitudes below the selected gates, not merely residuals or abundances that have stopped changing.
- Explicit failure when the pass budget ends without satisfying that joint test.
- Consistent face mass flux, density, temperature, composition, and radiation background throughout each block solve and after its update.
- Element and charge conservation and nonnegative abundances through every proposed composition blend or projection.

### 5.4 Removing the registry is not a complete molecular implementation

The chemistry routines already contain explicit fixed-fraction controls for transported `H2` and, when enabled, `H+`. This is important existing functionality; the review does not assume that every chemical sweep resets all transported species to local equilibrium.

Nevertheless, the input reader explicitly refuses `Ionization transport: True` with `Solver: Newton` unless `Coupled carrier solve: True` is also selected. The proposed universal switch to hydro-only Newton therefore does not currently support every transported configuration. Bypassing that refusal would require a validated contract for the fixed transported fractions, not deletion of the check alone.

The carrier chemical source also couples species through reaction rates, free-element budgets, electron density, temperature, and shielding. A scalar tridiagonal transport operator at fixed coefficients is not proof that the full molecular subproblem is a set of independent scalar tridiagonals. A block formulation or an internally converged chemical/transport iteration may be necessary.

Finally, the recorded small chemical contribution to one H2 row at one front cell does not establish weak chemical stiffness. The source value and its derivative are different quantities. Large production and loss derivatives can nearly cancel in the net source, and different cells or species can have very different time scales.

### 5.5 A concrete mismatch in the retained hydrodynamic solver description

The target says to retain the three-unknown JFNK “with its trust region.” The current implementation sets

```fortran
use_tr = (nspec_row .gt. 0)
```

and its environment override can disable that selection. Consequently, the default three-unknown path does not use this species trust-region branch. The retained component should be described as the existing hydrodynamic JFNK path with its actual globalization, not as a verified three-unknown trust-region solver.

## 6. What the published atmospheric papers support

Citation metadata were checked through NASA ADS before general web lookup. The local Koskinen and Huang PDFs have final journal publication headers. The numerical-method and boundary-condition passages were read with `pdftotext -layout`, rather than relying on search summaries.

### 6.1 Koskinen et al. (2013)

Section 2.1.3 describes splitting nonadvective and advective updates, updating the former before advection, van Leer advection, Crank–Nicolson treatment of viscosity and conduction, a 1 s time step, and periodic Shapiro filtering. It judges steady behavior using flux constancy and approximate energy-flux conservation. This supports the target's historical description of the time-integration approach. It does not establish an EXHALE species-residual tolerance or show that a final hydro-only Newton correction preserves transported abundances in steady balance. [Published paper, DOI 10.1016/j.icarus.2012.09.027](https://doi.org/10.1016/j.icarus.2012.09.027).

The lower-boundary discussion is **section 2.1.1**. Section 2.1.2 concerns the **upper** boundary. The references in target sections 5.3 and 8 should be corrected accordingly. The paper's molecular-validity qualifications also mean that its lower-boundary choice is not a universal argument for excluding EXHALE's molecular layer. [Local published PDF](../../references/Koskinen_2013Icarus_226_1678.pdf).

### 6.2 Koskinen et al. (2022)

Section 2.5 describes the time-dependent upper-atmosphere model, its Kinetic PreProcessor chemistry, semi-implicit diffusion, and flux-limited advection. Appendix B gives species continuity with chemical sources and diffusion velocities, including the zero net diffusive mass-flux condition in equation B10. It also describes the radiation update cadence and boundary conditions. Importantly, the discussion associated with Figures 15 and 16 explicitly examines transport–chemistry balance for H and H2. The paper does not treat stationary species balance as physically optional simply because its numerical presentation does not use EXHALE's normalized residual gate. [Published paper, DOI 10.3847/1538-4357/ac4f45](https://doi.org/10.3847/1538-4357/ac4f45); [local published PDF](../../references/Koskinen_2022_ApJ_929_52.pdf).

### 6.3 Huang et al. (2023)

Section 2.3 identifies the model as time-dependent CETIMB and refers to the 2022 Appendix B for general details. It describes a molecular/ionized atmosphere with transport and additional thermal and chemical processes, and an outflow upper boundary outside the sonic point for its application. This supports using time integration as a credible approach, not the specific EXHALE workflow “march species, then finish hydrodynamics once, then waive species acceptance.” [Published paper, DOI 10.3847/1538-4357/accd5e](https://doi.org/10.3847/1538-4357/accd5e); [local published PDF](../../references/Huang_2023_ApJ_951_123.pdf).

### 6.4 Consequences for section 5.2

The proposed method should be named an **EXHALE time-integration route informed by CETIMB**, with any hydrodynamic stationary correction identified separately.

Total mass-flux constancy is necessary for a stationary source-free mass equation but does not establish every species balance. Reactions redistribute species while conserving elemental nuclei. A physically valid molecular solution must satisfy

\[
\frac{1}{r^2}\frac{d}{dr}
\left[r^2(n_s v+\Phi_s^{\mathrm{diff}})\right]=P_s-L_s
\]

for each independent transported species, with appropriate boundary conditions. Agreement in total mass flux can coexist with incorrect H2 dissociation, ionization, or elemental partitioning.

A final hydrodynamic correction changes the velocity, density, temperature, and possibly the radiation field that determine the species residual. The species must be checked again afterward. If they no longer satisfy their gates, another coupled update is needed; the hydro correction is then part of an outer iteration, not a terminal certification step by itself.

### 6.5 A split fixed point need not solve the unsplit equations

Time integration with splitting can approach a fixed point of the discrete step map while retaining an unsplit stationary defect at finite time step. This is not a claim that all split schemes necessarily have such a defect; steady-state-preserving constructions can avoid it.

The executed scalar example splits `dx/dt = 1 - 2x` into two exactly integrated stages, `dx/dt = 1 - x` followed by `dx/dt = -x`. The continuous steady solution is `x = 0.5`; the split map instead has the fixed point `1/(exp(dt) + 1)`.

| Time step | Fixed point of the split map | Residual of the unsplit equation |
| --- | --- | --- |
| `0.1` | `0.475020813` | `0.049958375` |
| `0.05` | `0.487502604` | `0.024994793` |
| `0.025` | `0.493750326` | `0.012499349` |

Thus a longer march alone cannot distinguish incomplete relaxation from a persistent splitting defect. A time-step study, an unsplit residual measurement, and checks on radiation/chemistry refresh cadence are needed. Filters should likewise have an explicit physical and numerical justification; a historical filter is not permission to suppress an unresolved physical instability.

## 7. Well-balanced treatment: useful, but the current summary is inaccurate

### 7.1 N37 contradicts the claim that the measured floor was removed

N37 records substantial improvement in preserving the scheme's own discrete hydrostatic equilibrium. It separately records almost unchanged non-smoothness floors: for the energy row of cell 1, approximately `1.009e-11` with the option and `1.046e-11` without it. These are existing recorded measurements, not new executions.

Therefore section 5.4 must not state without qualification that the implemented option removes the rounding amplifier from the hydrodynamic rows. It removes an algebraically canceling equilibrium contribution in the constructed balance. Rounded pressure inputs and other row contributions still limit the measured action. Exact preservation of a specified discrete equilibrium and a lower floating-point floor for arbitrary perturbed states are different properties.

Likewise, N34/N35 do not measure the effectiveness of the subsequent N37 implementation. N37's own results must be used when assessing that option.

### 7.2 The momentum reference-scale defect is visible in the current code

With `well_balanced` enabled, `Source.f90` returns `S = 0`, because the equilibrium source has been absorbed analytically into the flux construction. This is not physically missing gravity.

However, `store_row_terms` still computes

```fortran
momentum_largest_term(j) = max(abs(dF(2,j)), abs(S(2,j)), abs(Smom(j)))
```

When `S(2) = 0` and other momentum sources are absent, the remaining momentum residual is also its own reference scale. Above the absolute floor, the normalized magnitude is then one, regardless of how small the dimensional imbalance becomes. The independent script reproduces this normalization identity; N37 records the corresponding production behavior.

This is a genuine code-interface defect for the enabled option. It is already acknowledged in N37, rather than newly discovered here. It prevents interpreting its current stationary-solver comparison as a clean test of the well-balanced discretization alone.

**Suggested correction, not implemented:** carry physically meaningful reference terms alongside the cancellation-free residual. For example, retain the reconstructed equilibrium pressure force or gravitational weight as a scale even though those contributions cancel analytically in the numerator. Include dynamic and other source terms and a controlled zero-gravity/near-zero-force limit. The numerator must remain the actual imbalance, not a restored cancellation-prone expression.

Validation should test the residual and its normalization separately on hydrostatic equilibria, perturbations, moving flows, and zero-gravity cases. Enabling the option by default should wait for that acceptance-interface correction and its targeted validation. Agreement with an old output is not the physical acceptance criterion.

### 7.3 Literature verification limit

The local Kappeli–Mishra files are explicitly institutional report versions, not verified final publisher PDFs. N37 itself records that the publisher versions remain wanted. I checked the journal metadata through ADS and attempted to retrieve the final 2016 A&A PDF, but the publisher returned HTTP 403. Consequently, I do not independently certify the target's detailed equation-number or section-number attribution to the final 2014/2016 papers.

The relevant publication records are [Kappeli and Mishra (2014), DOI 10.1016/j.jcp.2013.11.028](https://doi.org/10.1016/j.jcp.2013.11.028) and [Kappeli and Mishra (2016), DOI 10.1051/0004-6361/201527815](https://doi.org/10.1051/0004-6361/201527815). Those links identify the publications; the review's claims about EXHALE's implemented cancellation and normalization above follow from its own code and algebra, not from unverified final-paper passages.

The implementation should be judged against the discrete equilibrium it actually constructs. Preserving that equilibrium does not imply exact preservation of every irradiated, composition-stratified atmosphere, or eliminate boundary-driven drainage.

## 8. Moving the lower boundary or changing coordinates

### 8.1 These are changes to the modeled problem, not just solver remedies

A lower boundary near one microbar can be appropriate when an external lower-atmosphere calculation supplies consistent radius, pressure, temperature, and composition. It is not automatically appropriate for all EXHALE configurations or for an objective that includes molecular formation, dissociation, and transport within the removed region.

Before moving the boundary, specify what remains physically represented:

- Reservoir abundances and elemental supply, including H2 and any transported ions.
- Molecular shielding and the radiation column below or near the interface.
- Thermal conduction and the energy carried across the interface.
- Diffusive separation, eddy mixing, and the flux conditions for each independent element or carrier.
- The validity of fluid transport and the upper boundary relative to the sonic point and exobase.

The existence of `base.inp` and lower-atmosphere profile input means that some changes of boundary location may be achievable with existing interfaces. It is premature to say that every such change requires rewriting the grid and every reference result. A mass-coordinate reformulation is a much larger change than selecting an already supported boundary location.

### 8.2 A new coordinate does not generally linearize the whole atmosphere

For hydrostatic spherical gas,

\[
\frac{dp}{dr}=-\rho g_{\mathrm{eff}}(r).
\]

For a column mass coordinate increasing inward, `dm/dr = -rho`, one obtains `dp/dm = g_eff(r)`. This is a constant-coefficient linear relation only under additional assumptions, such as constant gravity. Radius-dependent gravity, Roche terms, composition-dependent density, and the coordinate mapping remain coupled.

For logarithmic pressure,

\[
\frac{dr}{d\ln p}=-\frac{p}{\rho g_{\mathrm{eff}}},
\]

which retains the equation of state and the temperature/composition dependence. Changing variables redistributes metric factors; it does not eliminate floating-point cancellation by itself. The transformed flux and source discretization still needs a consistent balance.

Finally, the target's association of layer nonstationarity solely with thin radial cells is too strong. N37 records that its mechanical-column drainage persisted and was controlled by the boundary rather than cured by equilibrium reconstruction. Physical forcing, boundary compatibility, and stability must be distinguished from grid spacing before removing the layer.

## 9. An omitted middle ground: partitioned components inside a coupled method

The choice is not limited to the existing band-preconditioned JFNK or a completely separate outer fixed point. Existing transport solves can also be used in a block preconditioner or nonlinear acceleration while retaining the same complete stationary residual as the target.

The final published Knoll–Keyes survey discusses pseudo-transient continuation in section 2.4.2, physics-based and split preconditioning in section 3.4, and nonlinear preconditioning in section 3.6. Its methods explicitly allow segregated components to serve a coupled Newton–Krylov method. This supports keeping such designs conceptually available; it does not predict that one will succeed on the present EXHALE fixtures. [Published survey, DOI 10.1016/j.jcp.2003.08.010](https://doi.org/10.1016/j.jcp.2003.08.010); [final journal PDF read for this review](https://www.inf.ufes.br/~luciac/comcie/JFNK-Knoll.pdf).

Possible hypotheses for a separately authorized follow-up are:

- Use the transport and reaction blocks to approximate the physically relevant coupled response, including a Schur-complement approximation where justified.
- Improve the derivative of transport rows without changing their physical residual, while retaining consistent treatment of limiter and constraint behavior.
- Accelerate a demonstrated contractive or nearly contractive outer iteration, always accepting against the original joint residual.
- Examine a coupled pseudo-transient formulation on the constrained state manifold, with differential variables and algebraic constraints treated appropriately.

EXHALE already has PTC and continuation-related code. The last item is not a claim that PTC is absent. It asks whether the implemented time shift, variable scaling, species registry, and constraint treatment correspond to the proposed coupled continuation problem.

None of these possibilities warrants an open-ended new solver campaign now. Their relevance is that N31–N38 have not exhausted the entire class of coupled methods. Stopping the current experiment is reasonable without claiming a mathematical impossibility that the evidence does not establish.

## 10. Recommended decision framework

### 10.1 Decisions that should remain separate

| Decision | Recommendation | Reason |
| --- | --- | --- |
| Resume the existing N0–N38 campaign? | Keep the pause. | A new experiment should answer a new question, not merely increase budgets again. |
| Accept the target's impossibility diagnosis? | No; replace it with a bounded implementation-specific conclusion. | Ritz compressions, rounding diagnostics, and trial failures do not prove it. |
| Change the selected `1e-5` wind gate? | Do not change it through this review. Correct the justification and revisit only with explicit approval and evidence. | An unconverged candidate's residual cannot define its own acceptance threshold. |
| Declare layer species behavior certified? | No; preserve the explicit report-only scope. | Wind acceptance does not make the layer stationary or physically irrelevant. |
| Try a partitioned stationary approach? | Consider a bounded evaluation of the existing route after defining joint acceptance. | Much code exists, but convergence and molecular closure are not established. |
| Use time integration as the physical reference? | Yes, conditionally, with time-step, history, and unsplit-balance checks. | It can establish dynamic behavior that a stationary root alone cannot. |
| Add a hydro-only finish to a march? | Only if the returned state's species and closure are rechecked. | Hydrodynamic changes can undo the species balance. |
| Enable well-balanced treatment by default? | Not before correcting the momentum reference scale and validating the intended equilibrium. | The current enabled-option acceptance measure is degenerate. |
| Remove the deep layer? | Decide by the physical scope and boundary model, not solver inconvenience. | It changes molecular supply, shielding, and thermal coupling. |

### 10.2 A useful first measurement is a baseline, not a binary decision

The proposed residual evaluation of an existing marched state is worthwhile. It does not by itself decide between time integration and partitioned stationarity, and a new long march is not necessarily cheap when slow diffusive modes remain.

If follow-up work is authorized, the first stage should read available saved states and evaluate, without advancing them:

1. The unsplit hydrodynamic and species residuals on each state's own composition and background.
2. Elemental and total mass fluxes on the numerical faces used by the operators.
3. Mass and charge closure, positivity, any active abundance constraints, and the validity of the chemical closure.
4. The same quantities after any proposed terminal hydrodynamic correction.
5. Separate maxima and locations in the wind, transition region, and deep layer, using the existing regime definitions rather than a newly convenient boundary.

Recorded states must have compatible configuration and build provenance. A different molecular network, boundary reservoir, or transport ownership is not a comparison of solvers for the same equations.

### 10.3 Then isolate the uncertainty

**If the unsplit residual decreases with longer physical evolution:** estimate the remaining transport time scales and compare time steps on a matched physical interval. Step count alone is not elapsed relaxation time. Check that skipped or incomplete carrier advancement has not invalidated the history.

**If the split map stops changing while the unsplit residual remains appreciable:** test time-step dependence and update cadence before attributing the plateau to spatial discretization.

**If frozen-background species relaxation is considered:** first establish that its returned state satisfies the same stationary species operator and constraints. Compare its entry and exit residuals, not just its drift or internal step count. Then measure the response of the wind–composition outer iteration and require joint certification.

**If the coupled linearization remains the focus:** use a small controlled configuration where explicit derivative comparisons or singular-value diagnostics are affordable. Preserve the same boundary and constraint treatment. Distinguish the full operator, its preconditioner, its masked compressions, and its Schur response. Do not use a single Ritz-value ratio as the deciding condition estimate.

**If a spatial-accuracy claim is needed:** solve compatible problems independently on multiple meshes, with algebraic errors controlled on each mesh. Compare physical profiles, front positions, integrated fluxes, and observables in common coordinates. Merely remapping one stalled state is insufficient.

### 10.4 Acceptance must describe both physics and numerical status

A converged algebraic root is not automatically the physical atmosphere: it can be dynamically unstable, violate a modeling approximation, or use incompatible boundaries. Conversely, failure of this Newton implementation does not prove physical unsteadiness.

Report at least these statuses separately:

- Which discrete equations meet their stationary gates.
- Which regions are outside that certification scope.
- Whether transport history and conservation constraints are valid.
- Whether the solution is stable or evolving under the intended physical time equations.
- What spatial and temporal accuracy has actually been demonstrated.

This prevents both errors: certifying an unphysical numerical root, and waiving a missing physical balance because a solver plateau appears small.

## 11. Suggested replacement for the central conclusion

The following wording preserves the useful result without overstating it:

> The tested species-row JFNK configurations have not delivered a state that satisfies the selected stationary gates. The measurements identify inaccurate matrix-free actions in some directions and poor progress of the present right-preconditioned Krylov solve, with difficult directions associated with transported species. They do not establish nonexistence of a discrete stationary root or a residual floor set by spatial truncation error. Further work should begin with a common full-system residual and compare the existing partitioned and time-integration routes under that contract. Any change in certification scope, boundary physics, or discretization requires its own justification and decision.

## 12. Validation performed and retained artifacts

The independent checks are retained in [solver_approach_checks_20260910.py](audit_20260905/solver_approach_checks_20260910.py), as requested for review test code.

Run from `EXHALE_v1.00`:

```bash
python3 docs/audit_20260905/solver_approach_checks_20260910.py
```

All illustrative assertions passed. The checks demonstrate:

1. A discrete residual can be far smaller than the spatial solution error.
2. Diagonal similarity can change finite-step GMRES behavior despite identical eigenvalues.
3. A singular block compression does not imply a singular full coupled operator.
4. Individually solvable blocks need not yield a convergent damped alternation.
5. A fixed point of a split step can have a nonzero unsplit residual.
6. Using a residual as its own scale makes its normalized magnitude insensitive to improvement.

These are mathematical counterexamples, not EXHALE regression results. No full regression suite, production restart, long integration, or new grid sequence was run: the deliverable is an implementation- and literature-based review, and production behavior was not changed. Historical numerical results remain labeled as recorded. Final-paper verification was completed for the atmospheric method passages and the Knoll–Keyes survey; final publisher-text verification of the two Kappeli–Mishra papers remains incomplete because of access restrictions.

The practical outcome is not a recommendation to abandon stationary methods or to accept a larger residual. It is a recommendation to correct the interpretation of the evidence, preserve the physical equations and explicit acceptance scope, and make any next experiment discriminate among mechanisms that remain genuinely unresolved.
