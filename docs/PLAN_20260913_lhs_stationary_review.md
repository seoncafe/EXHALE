# Review of the LHS 1140 b Stationary-Solver Plan

Date: 2026-09-13

## Assessment

The well-balanced stationary route and an explicit atomic-to-molecular initialization product are reasonable directions. The plan is not ready to execute unchanged. Its molecular construction contains concrete inconsistencies with the implementation, its identity test is not a valid physical invariant, and the subsonic outer boundary deserves a separate acceptance item.

Recommended scope: finish a controlled L4d experiment, investigate L6 before accepting metals-on results, and develop a small molecular initialization pilot together with the ghost-cell closure and caloric-domain diagnostics. Do not require the complete atomic campaign before investigating the molecular blockers. One suitable certified atomic seed is sufficient for that investigation.

Evidence: direct inspection of the current working tree at HEAD `43bc28cef58772560bae019559f3592335922212`, including the uncommitted changes to `steady_newton.f90`. Existing modifications and products were preserved. All atmospheric numbers quoted by the plan remain **READ**, not independently remeasured in this review. No solver run, compilation, or campaign was performed. Formulas below are **DERIVED** from the stated definitions and source; code constants and behavior are **READ** from implementation. This report is a review, not an implementation authorization or a success certificate.

## 1. L4d is a useful experiment, not yet a sufficient solution

The working-tree diff already implements both proposed changes in `src/modules/time_step/steady_newton.f90`:

- Around lines 15870–15895, a residual-ratio forcing term replaces the fixed Krylov tolerance on the relevant non-trust-region path.
- Around lines 16180–16200, a lower merit resets the stagnation counter even if the best judged distance does not improve.

These are no longer merely absent features proposed in the plan. Their effectiveness still needs the failing-rung reproduction and controls. Existing ongoing edits should remain under their current owner's control.

The original stagnation rule can indeed stop an iteration that is making useful nonlinear progress. However, the new comparison `f2 < f2_best` has two limitations:

1. `cell_state_scales` and, when relevant, `cell_row_scales` are rebuilt around lines 15540–15562. The merit is subsequently `norm(F/D)` or `norm(F/Drow)`. Historical merits with changing denominators do not necessarily measure improvement in the same norm. A declining merit can partly reflect rescaling.
2. Any strict decrease, including insignificant numerical changes, resets the counter. That avoids premature termination but can consume the iteration budget without meaningful progress.

Use a fixed diagnostic scaling for the stagnation history, or reevaluate retained residuals under the current scaling. Keep the physical certification gates separate. Define meaningful progress relative to the measured residual-evaluation uncertainty and a rolling window; do not select a threshold solely to recover an existing pass count.

The new forcing range is READ as `1e-4` to `1e-1`. Tighter global Krylov tolerance may help, but a Euclidean residual dominated by energy does not guarantee adequate continuity accuracy. Also, the new safeguard `0.9*eta_prev**1.618... > 0.1` cannot activate while `eta_prev <= 0.1`: its maximum is about 0.0217 (DERIVED). This does not invalidate the method, but the safeguard is ineffective under the selected cap. The paper attribution in the new comment was not independently verified here.

Measure the actual unpreconditioned linear residual separately for mass, momentum, and energy, as well as the nonlinear trial residual. Distinguish failure to reach the requested Krylov tolerance from failure of the nonlinear model. If continuity remains the obstruction, investigate row/block scaling or a continuity-constrained step rather than indefinitely tightening a single scalar tolerance.

The HD 209458 b control is valuable, but success on one planet does not prove that the route's only LHS obstruction is its seed. Conditioning, coupling, boundary conditions, and closure validity change with the physical regime. Replace that absolute conclusion in section 1 with a narrower statement.

## 2. L7(A) contains three concrete construction errors

### 2.1 The stated partition is not what `set_IC` computes

`src/modules/init/set_IC.f90:211–225` obtains `x2_ic` and then uses:

```fortran
dfHI_ic = min(x2_ic/mass_per_H, f_sp(j,isp_HI))
f_sp(j,isp_H2) = 0.5d0*dfHI_ic
f_sp(j,isp_HI) = f_sp(j,isp_HI) - dfHI_ic
```

Under that initializer's abundance normalization, this allocates a fraction of the total hydrogen budget, limited by available neutral H. The plan instead proposes `n(H2)=x2*n(HI)/2`, which allocates a fraction of the neutral budget. They agree for fully neutral hydrogen, but not in partially ionized layers.

Choose the intended seed prescription explicitly. For a total-hydrogen target, a density-space construction is

```text
delta_nHI = min(x2 * nH_nuclei, nHI_available)
nH2_new = nH2_old + delta_nHI/2
nHI_new = nHI_old - delta_nHI
```

A neutral-weighted seed is also possible as an initialization choice, but must not be described as an exact reuse of `set_IC`. With element diffusion, use each cell's actual element census rather than assuming the base He/H applies at every radius. Distinguish an imposed boundary partition from an arbitrary extension of that partition through the entire wind.

### 2.2 The plan applies the mixing-ratio ceiling to the wrong variable

`composition.f90:298–353` distinguishes the molecular mixing ratio `q` from the hydrogen-nuclei fraction `x2`. With `y=He/H`,

```text
q_max = 1/(1+2*y)
x2 = 2*q*(1+y)/(1+q)
0 <= q <= q_max  corresponds to  0 <= x2 <= 1.
```

Section 5.2 instead requires `x2` to lie below `q_max`. This wrongly rejects valid molecular fractions. For example, `y=2.13`, `q=0.1` gives `q_max=0.190114...` and `x2=0.569091...` (DERIVED): the input is admissible, but the proposed check would reject it.

Check `q` against its ceiling, `x2` against unity, and the actual density transfer against the available hydrogen budget.

### 2.3 Copying the hydro file unchanged contradicts the loader contract

`load_IC.f90:1162–1171` compares the metadata fields of both restart files and refuses disagreement. Updating the species file's options/species metadata while copying the hydro file unchanged therefore fails before the intended molecular solve on schema-bearing inputs.

Generate mutually consistent target metadata in **both** files, including the target species declaration and applicable reservoir fields. Derive options from the target configuration; do not hard-code `molbase=T carrier=T` for every molecular group. Keep the atomic certification as provenance only, and write the new state as an uncertified initialization seed.

## 3. Molecular conversion must choose a thermodynamic invariant

The plan says that the hydrodynamic state does not know about the partition. The implementation contradicts this: `load_IC` reconstructs mass density from species, but reads pressure and temperature from the hydro file; it returns primitive `W=(rho,v,p)` around lines 834–837. The caloric EOS uses the molecular composition, and the conserved/primitive conversion routines call it.

Converting two H atoms to one H2 molecule preserves nuclei and mass but decreases the number of particles. With initially unchanged electrons and other species,

```text
n_particles,new = n_particles,old - delta_nHI/2
p = n_particles * k_B * T.
```

Consequently, rho, v, p, and T cannot all remain unchanged. Molecular internal energy also changes the caloric relation. A seed need not conserve the source state's energy, but it must state that choice and be internally consistent.

For this nearly hydrostatic wind, a useful first pilot is to preserve `rho`, `v`, and `p`, recompute T from the new particle census, and recompute energy with the production EOS. This preserves the initial pressure-force profile, but changes reaction rates and does not preserve chemical/energy equilibrium. An alternative preserves `rho`, `v`, and T and recomputes pressure and energy, followed by a pressure-balance adjustment. Compare these alternatives on one seed; do not promise that only chemistry will move.

Prefer an explicit conversion executable that calls the production census, EOS, boundary, and metadata routines. A Python utility can orchestrate it, but independently reproducing these physical rules in Python introduces a second implementation to maintain. This remains an explicit initialization product, not permission to weaken ordinary restart validation.

Option B does not become numerically better simply because option A's seed is distant: identical partitions and thermodynamics produce the same distance regardless of whether construction occurs offline or in `load_IC`.

## 4. Replace the L7 acceptance tests

**T-L7-1 is not a valid molecular identity test as written.** Zero initial H2 does not disable the molecular network. The molecular source equations contain formation terms that can create molecules from atomic material; see `System_HeH_mol.f90`, including the density-dependent terms near lines 252–297. Molecular equilibrium is not required to reproduce the atomic solution with an H2 column below the proposed threshold.

Separate the tests:

1. Conversion identity: zero requested transfer preserves source species/nuclei and the explicitly chosen primitive invariants before chemistry is enabled.
2. EOS closure: the written composition, pressure, temperature, and energy describe the same state.
3. Loader acceptance: both headers and all target species columns are valid, and the file is treated as initialization rather than continuation.
4. Chemical validity: reaction residuals, charge closure, and element budgets after equilibration satisfy their physical criteria. Do not require H2 to remain zero.
5. End-to-end convergence: one LHS seed reaches joint hydro/species certification under the target physics.

The hot-Uranus test must define how a molecular state is first projected into an atomic seed. Removing columns alone loses nuclei in molecular ions and, for HeH+, helium. Specify ion charge/electron handling and the thermodynamic invariant before defining a round-trip result.

A failed first chemical sweep also does not prove that no root exists. Preserve the failing cell and distinguish an inconsistent reservoir constraint, a poor initial guess, an inaccurate Jacobian/scaling, and genuinely incompatible equations. Examine multiple admissible starting states and reaction residuals before attributing failure to low-temperature rate validity.

## 5. Resolve ghost chemistry and the caloric domain as independent issues

The molecular-cell test in `states/caloric_eos.f90:215–230` activates the molecular EOS whenever H2 density is positive. `time_step/energy_semi_implicit.f90:790–798` then applies the molecular ceiling, READ as 50000 K. Thus even a trace positive H2 abundance can select that ceiling.

The plan's statement that a hot cell has no H2 worth considering is not sufficient to remove the ceiling. Transported or nonequilibrium H2 need not vanish instantly. An arbitrary abundance cutoff would introduce a discontinuity and could violate energy consistency.

Inspect the energy bracket with composition held exactly as the affected operator holds it. Then decide whether to extend the molecular caloric model consistently, change the coupled chemistry/energy solve, or explicitly refuse states outside a documented validity domain. Any high-temperature extension must have consistent internal energy and heat capacity and avoid counting dissociation energy twice.

The ghost-cell failure should be reproduced with the actual boundary-imposed composition and thermodynamics, not with an unconstrained interior equilibrium solve. A reservoir prescribing a molecular partition is not automatically in local chemical equilibrium with the irradiated ghost state. Determine which quantities the boundary is allowed to prescribe and which chemical rows it must solve. Do not simply exempt all ghosts from checks or tighten the iteration tolerance without this analysis.

## 6. Add an outer-boundary acceptance item

`states/Apply_BC.f90:192–313` implements an isothermal hydrostatic continuation into outer ghosts. Its opening claim that an outflow has no incoming characteristic is only true for supersonic outflow. For subsonic outward Euler flow, the characteristic speeds are `v-c`, `v`, and `v+c`, and `v-c` enters the domain.

The code later acknowledges that its continuation matters at subsonic edges, but the plan does not establish the external condition it represents for this LHS solution. A converged state can therefore be a valid discrete solution of the selected boundary prescription without being a boundary-independent physical mass-loss prediction.

State the intended external-pressure or incoming-characteristic condition. Test sensitivity to outer radius and that physical boundary condition, and distinguish a breeze selected by the environment from a transonic escape solution. Extending the domain is appropriate only while the physical model remains valid. Do not infer physical uniqueness from agreement among several seeds at one fixed boundary.

Well-balanced reconstruction is still appropriate. The local published copy `references/Kappeli_2016A&A_587_A94.pdf` describes preserving hydrostatic balance through pressure reconstruction combined with the gravitational source and suitable fluxes. That supports the numerical direction, not the physical uniqueness of the LHS boundary-value problem. A detailed new validation of that paper's method was not performed here.

## 7. L6, post-processing, and the order of work

L6 should remain a gate for accepting metals-on results. A photon grid reaching below the C I edge establishes wavelength support, not the actual photoionization rate. Audit the local attenuated rate, recombination, charge exchange, density convention, and the cooling contributions at the same state. Compare the implemented excitation/de-excitation and radiative treatment within their density and temperature domains. A large neutral fraction is a reason to investigate, not by itself proof of an error.

`EXHALE_main.f90:3928` exits for `do_only_pp` late in the time loop, after state updates and output handling. The tiny-CFL recipe suppresses a change rather than making post-processing read-only. The source also activates secondary ionization for this mode. Therefore explicitly verify the physics configuration used by solve, certification, and output; do not assume the output operation preserves the certified state.

Create a direct diagnostic/output path that evaluates the held state without a time update. Whether it may re-equilibrate chemistry should be a separate explicit choice. Until that exists, treat the tiny-CFL recipe as a temporary workaround and recheck the resulting state; it is not a mathematical identity operation.

Recommended sequence:

1. Finish the already-started L4d failing-rung experiment and its named controls. Track fixed-scale progress and actual linear residuals.
2. Diagnose L6 before accepting metals-on campaign results; run boundary and grid sensitivity on a small representative atomic subset.
3. Correct the L7 specification and implement one explicit, thermodynamically closed initialization pilot. Inspect its ghost state before a long run.
4. Address the ghost closure and caloric-domain problems using isolated reproductions. Use adaptive continuation in an explicitly stated physical boundary/composition parameter if the direct seed remains too distant, and certify only the final target problem.
5. Expand to the campaign only after the representative configurations satisfy joint physical/numerical acceptance and output preserves that meaning.

Suggested user decision: approve the explicit initialization product plus investigation of the two molecular blockers as one bounded pilot. Do not yet authorize a general restart-layout exception or a broad solver rewrite. A later continuity-constrained stationary formulation may be justified if L4d fails, but it should preserve the actual numerical face mass flux, not merely impose constant cell-centered `rho*v*r^2`.

Only this review document was added. Existing source edits, model outputs, and ongoing work were not changed.
