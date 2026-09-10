# Review of the EXHALE development plan dated 2026-09-05

Review date: 2026-09-05 (KST).

Reviewed document: [development_plan_20260905.md](development_plan_20260905.md).
Implementation checked: revision `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`.

## 1. Assessment

The plan is a useful development outline, but it is not yet a sufficiently consistent implementation specification. Its main direction is reasonable: fix demonstrated conservation and numerical defects before claiming a validated molecular wind, unify thermodynamics, couple species transport to the steady solver, and compare published atmospheres only after checking their constituent processes.

Several details require correction before implementation. The most important are the energy-update design and its dependencies, an incorrect description of the existing steady residual, and an H3+ extrapolation prescription that does not satisfy the physical limit it claims to enforce. Other corrections concern incomplete derivatives, missing test coverage, and work dependencies that would cause repeated redesign.

The recommended disposition is to retain the overall scope and most of the task inventory, revise the findings below, and then implement bounded tasks against explicit acceptance conditions. This review does not recommend rejecting the entire plan or restarting the project design from nothing.

Production source, inputs, simulation products, and regression reference products were not changed during this review. The reviewed plan itself was not edited.

## 2. What was independently checked

The complete 430-line plan was read. Targeted source inspection covered:

- The marching composition and energy update in `EXHALE_main.f90`.
- `steady_newton::eval_residual`, its final consistency check, and the actual Jacobian-vector product.
- The semi-implicit energy residual, derivative, iteration sequence, and returned cooling.
- H3+ non-LTE table values and their density coordinates.
- Molecular reaction-energy exclusions and existing cooling consumers.
- Relevant input defaults and the configurations of the regression cases cited by the plan.
- The Makefile `check` target and the regression script's case-selection interface.

Two limited calculations were executed:

1. Numerical arrays in `Hydro_ioniz.txt` and `Ion_species.txt` were read for the twelve named cases and compared with their existing reference arrays. All twelve cases differed in both files. This verifies a difference in numerical contents, beyond timestamps or comment headers. It does not establish the cause, magnitude, or provenance of each difference.
2. The exponent implied by the proposed two-column H3+ power-law continuation was calculated from the actual 1000 K table entries. The result is 0.459740201590, rather than the required low-density exponent 1.

No atmosphere was rerun. No new measurement of mass loss, base oscillations, full-network Jacobians, or physical convergence was made. Statements about those effects below are either source-based conclusions or proposed tests, not newly demonstrated planetary results.

## 3. Proposals that are reasonable

The following parts can be retained with their intended scope:

| Proposal | Assessment |
|---|---|
| A1: correct the LLF acoustic bound | Correct. Check both direct LLF use and HLLC stage repair. |
| A2: correct the two ROE formulas | Correct for the constant-gamma formulas. Restricting an unverified molecular ROE path is reasonable. |
| A3: one species-dependent reference scale and a final-state residual | Correct direction. Its error status and CO-ceiling behavior need the dependency changes below. |
| A5-A6: explicit photon and reaction energy recipients | Necessary. Energy closure must be distinguished from a uniquely known energy partition. |
| A7: residual-based source convergence and cooling at the returned temperature | Correct. The residual and derivative must use the same heating function. |
| B1: common molecular thermodynamics | Necessary. The standard-state free energy requires more than an internal partition function. |
| C1: variable species unknowns in the steady system | Better than another special post-convergence proton correction. |
| C2: equilibrium elimination based on relaxation times | Physically appropriate, provided all relevant forcing times are considered. |
| C3: species/radiation convergence as well as a scalar norm | Necessary and well motivated by the current implementation. |
| D1-D3: common transport fluxes and verified boundary constraints | Important. Their mathematical design should precede the steady-system implementation. |
| Tests with explicit failure status | A substantial improvement over diagnostics that only print numbers. |
| Term-by-term reconstruction of the Koskinen comparison | Appropriate. The discrepancy should not be attributed to either code before those checks. |

The plan also correctly treats the thermal-energy and formation-energy-inclusive conventions as alternative valid formulations. Retaining a thermal-energy variable can be a reasonable initial choice. The issue is whether the resulting discretization conserves energy, not whether one particular storage convention is universally preferable.

## 4. Findings requiring revision

### R1. A4 does not yet define an energy-conserving reacting update

**Priority: high. References: plan lines 113-144; B3 at lines 240-249.**

A4 proposes to solve composition at the old temperature, retain thermal energy, and invert the EOS to obtain a new temperature. It then adds reaction-energy sources. Retaining thermal energy during an otherwise source-free EOS conversion removes the old particle-count energy injection, but this operation alone is not conservation of thermal plus chemical energy during a reaction.

For a homogeneous cell, the necessary accounting identity is

\[
E_{\rm th}^{n+1}-E_{\rm th}^{n}
=-\sum_s\epsilon_s^0(n_s^{n+1}-n_s^n)
 + Q_{\rm nonchemical},
\]

where the last term is the integrated energy supplied or removed through explicitly defined external/radiative channels. If the species change is generated by finite-rate chemistry, its reaction progress and the integrated chemical energy term must agree. If composition is imposed through instantaneous equilibrium instead, energy conservation and equilibrium must be solved together. A finite composition jump cannot generally be reconciled by subsequently multiplying an unrelated set of final-state rates by a time step.

There are three problems in the current description:

- The composition returned by an equilibrium solve at the old temperature is generally not equilibrium at the newly inverted temperature. A further coupled solve or a controlled finite-rate update is required.
- The phrase saying that binding energies are partly represented by the old projection is physically misleading. That projection changes sensible energy through particle counts and molecular thermal excitation; it is not a representation of ionization potentials or molecular bond energies.
- Adding the binding energy of every composition change on top of existing collisional-ionization cooling and molecular reaction heat can count the same contribution twice. The current `molecular_reaction_heat` exclusions and `eval_cool` consumers must be reconciled through an explicit reaction ledger.

**Revision:** move the mathematical design of B3 into the prerequisite for A4. A4 and B3 may be implemented incrementally, but Phase 1 must not claim a closed reacting-cell energy gate using an update that still lacks the required coupling.

Separate three kinds of operation in the derivation: a change of thermodynamic coordinates at fixed physical state, a chemical reaction, and spatial transport of species and their energy. The same formula cannot be assigned to all three without specifying the associated fluxes and sources.

The closed-cell test also needs a precise definition. Turning off incident radiation does not disable radiative recombination or spontaneous emission. Either turn off those radiative processes for an isolated material-energy test, retain the radiation energy inside the closed system, or include escaping radiation in the measured total budget. Conservation error should be compared with nonlinear tolerances and round-off; time-step convergence of the trajectory is a separate test.

### R2. The plan's description of `eval_residual` is incorrect

**Priority: high. References: plan lines 39-44 and 131-144.**

The plan states that the steady residual follows the same temperature-from-pressure, sweep, pressure-from-temperature sequence as the marching loop. It does not.

The actual paths are:

| Marching update | `eval_residual` |
|---|---|
| Derive temperature from current pressure/composition | Derive temperature from the seed composition and prescribed U |
| Transport composition and call `ioniz_eq` | Call `ioniz_eq` |
| Obtain new particle counts | Obtain new particle counts |
| Call `comp_p_from_T`, then `W_to_U` | Reapply ghost boundaries and call `assemble_residual` with the interior U retained |

Evidence: [EXHALE_main.f90](../src/EXHALE_main.f90), lines 1162-1228, and [steady_newton.f90](../src/modules/time_step/steady_newton.f90), lines 562-606. `Apply_BC` rewrites ghost states; it does not perform the claimed interior energy reconstruction.

This matters for A4: there is no corresponding interior `comp_p_from_T` operation in `eval_residual` to replace mechanically. The required steady-side work is self-consistent evaluation of temperature, composition, radiation, and the physical energy residual at the specified unknowns.

At an exactly self-consistent steady state, the equilibrium sweep can return its seed composition. That observation is reasonable, but it is not sufficient to assert that the implemented finite-step marching map has no projection contribution or that every atomic solution will remain unchanged after A4-A7. Operator splitting and incomplete closure can affect that conclusion.

**Revision:** correct the source description and state separately what is being held fixed in evolution, residual evaluation, and final composition closure. Treat preservation of a previously converged atomic solution as a hypothesis to test for the same physical equations, not an established consequence of the current implementation.

### R3. B2's suggested power-law continuation violates its stated low-density limit

**Priority: high. Reference: plan lines 231-236.**

For the collisional H3+ cooling component, the plan requires cooling proportional to n(H3+) n(H2) at low H2 density. It then permits a power-law continuation inferred from the first two table columns. Those are not equivalent prescriptions.

At 1000 K the actual table gives s=0.0013 at n(H2)=1e6 cm^-3 and s=0.0108 at 1e8 cm^-3. A power law fitted through those points has

\[
\alpha=\frac{\ln(0.0108/0.0013)}{\ln(10^8/10^6)}
      =0.459740201590.
\]

It approaches zero, but as n(H2)^0.45974 rather than n(H2). Reducing density by a factor of 1000 gives a cooling factor of 0.0417618 with this continuation, compared with 0.001 for a linear asymptote with the same boundary normalization. Thus the proposed test that cooling merely tends to zero would accept the wrong density scaling.

**Revision:** enforce the physical exponent 1 in the asymptotic collisional limit. Determine its normalization from a valid low-density rate model or a stated matching approximation, and use a smooth transition to the tabulated regime. Do not assume the first two tabulated columns already lie in that asymptote.

For temperature joins, require continuity by construction and inspect derivative continuity where needed by the nonlinear solver. A relative jump below 1e-3 still allows a finite jump; it should not be the definition of a continuous physical coefficient. Fit accuracy against published emission data and mathematical continuity are separate acceptance conditions.

### R4. A7 must change the heating residual together with its derivative

**Priority: high. Reference: plan lines 180-187.**

The plan offers adding the derivative of molecular chemical heating or documenting that heating is frozen. These are different numerical models, not interchangeable ways to complete the same implicit solve.

The existing energy residual uses a fixed `heat(j)`. If the new residual is still defined that way, adding dH/dT changes only the iteration matrix and does not solve an equation containing H(T). If the intended residual contains H(T), then H itself must be recalculated at the trial state:

\[
F(T)=u(T)-u(T_{\rm old})-\Delta t[H(T)-C(T)].
\]

The derivative then follows this same function, including composition derivatives when chemistry is coupled. Holding chemical heating fixed can remain a documented split approximation only with an accuracy/step-size criterion. A header describing a lag is not sufficient to demonstrate stability or physical accuracy in a stiff dissociation layer.

Retaining absolute-value damping as a trial iteration matrix can be acceptable with a reliable residual test and globalization. It does not guarantee convergence to the intended physical root, particularly with multiple thermal roots.

**Revision:** explicitly specify the residual being solved, its derivative or approximate iteration matrix, the physical root-selection policy, and failure/retry behavior. Coordinate this work with A4 and B3.

### R5. The dependency graph postpones equations that the steady solver needs

**Priority: high. References: A3, C1, D1, and plan lines 400-407.**

The plan first implements a larger steady system using the current carrier operator, then replaces that operator's spatial fluxes and geometry in D1. This is possible, but the stated dependency D1-after-C1 is not necessary and is likely to duplicate Jacobian, residual, and boundary work.

There is also a direct dependency between A3 and the CO limitation. A3 requires that a state modified by the CO ceiling satisfy the original carrier residual or fail the step. However, CO has no kinetic source in that residual, while its ceiling imposes an extra destruction operation. Repeatedly subdividing a step cannot be relied on to reconcile a kinetic equation with an incompatible algebraic constraint. It can instead produce repeated rejection or impractically small steps.

**Revision:** settle the species-flux, cell-volume, boundary, and CO-constraint equations before finalizing C1's residual and preconditioner. The complete CO chemistry can remain later work, but Phase 1 must define a valid restricted configuration or an explicit constraint model for any case it intends to accept.

A6 and A7 also depend on the energy ledger selected by A4/B3; they are not fully independent tasks. Independent algebraic fixes and read-only diagnostics can proceed separately, but shared energy/source interfaces need one agreed derivation first.

### R6. The named regression cases do not exercise transported H+

**Priority: high for verification. References: A3 and C1-C2.**

The plan predicts movement in `mol_sec_ion` from the carrier reference-scale correction. Its current input enables molecular chemistry and immediate secondary ionization, but does not enable molecular carrier transport or ionization transport. Those switches default to false in this configuration. Secondary ionization and advection of H+ are distinct processes.

The `mol_carrier` case enables H2 carrier transport, but does not enable `Ionization transport`. A search of the inspected `backup/regression/*/input.inp` files found no explicit transported-H+ option. These cases are useful controls, but they cannot replace a test of the proton row being corrected.

**Revision:** add dedicated H+ configurations with zero oxygen and zero/tiny initial proton abundance, plus a mixed H2/H+ front. The proposed local directional-derivative test is good and should call the actual source assembly; supplement it with a transported-proton integration test.

C1 should compare a converged steady solution with a time-converged reference using the same physical and spatial equations. It should not require equality to an old finite-step snapshot. Replace the general P r/v criterion with the network's relevant relaxation times, and include diffusion, temperature variation, and irradiation changes in the forcing times. C2's comparison with digitized He+ points is an external validation check, not a substitute for a continuity-equation test.

### R7. C4 incorrectly calls the JFNK action exact

**Priority: medium. Reference: plan lines 303-306.**

The actual code uses

```fortran
Jv = (Fp - F0)/eps
```

after evaluating `eval_residual(Y + eps*v, ...)`; see `steady_newton.f90`, lines 1376-1385. This is a finite-difference approximation, with perturbation-size, rounding, and inner-solve errors. Calling the full residual does not make this action exact.

Furthermore, differentiating a residual that uses a frozen seed or incomplete radiation/composition closure is different from differentiating the fully eliminated, self-consistent physical system. The existing photoionization calls precede the local chemistry sweep, and the steady path has a separate composition-consistency iteration. Their roles need to be specified.

**Revision:** state which variables are independent, frozen, and eliminated. Verify directional derivatives over a perturbation-size range and verify inner closure separately. The banded preconditioner may intentionally omit long-range coupling; the physical residual must still contain the required response.

### R8. B1 needs full thermodynamic normalization and a justified spin policy

**Priority: medium. Reference: plan lines 209-222.**

Using one level model is reasonable, but the instruction to replace the Shomate H2 entry with Z_rv is incomplete. A dimensionless internal partition function is not itself a standard molar Gibbs energy. The construction must include translation, standard pressure, formation-energy reference, statistical-weight convention, and the same molecular mass convention used elsewhere.

The selected line-derived level list also needs a completeness/range check before it is made the source for chemical equilibrium up to 10000 K. A useful infrared data set is not automatically a complete thermodynamic state sum over that full range.

Equilibrium ortho/para populations should be supported by a conversion-timescale argument in the intended gas. The plan currently makes equilibrium the default unless a contrary argument is written. For physical acceptance the burden is to state and check the conditions under which that equilibrium assumption holds.

**Revision:** derive standard-state thermodynamic functions from the selected spectrum and references, then test cycles spanning modules. Describe the expected change in atmospheric outputs as unmeasured: the 0.6-8% difference in rovibrational energy from the audit does not predict a similar percentage change in a molecular snapshot.

### R9. A5 closes an energy deficit but does not determine a unique energy distribution

**Priority: medium. Reference: plan lines 150-160.**

The proposed separation of chemical, electron, and fragment energy is correct in purpose. However, 19.725 eV is the difference between the adopted vertical threshold and the chemical energy of the represented asymptotic products. It is not, by itself, a measurement establishing the same fragment kinetic energy for every event and photon energy.

Likewise, evaluating an electron-degradation fraction at the mean electron energy is an approximation when that fraction is nonlinear in energy. The average of a nonlinear response need not equal the response at the average energy.

**Revision:** first impose the exact event-energy sum. Then state the adopted distribution or simplified partition and its validity range. Compare alternatives where the distribution is uncertain. Thermalization of energetic fragments should use their energy-dependent stopping/collision scale, rather than assuming a thermal-gas mean free path alone establishes local deposition.

This qualifies the implementation prescription; it does not remove the need to correct the missing energy recipient identified in the audit.

### R10. Baseline management is useful but should not block scientific work

**Priority: medium. References: plan lines 34-38, 196-198, and 355-357.**

The twelve current/reference numerical differences were confirmed during this review. Establishing a reproducible baseline before a sequence of edits is reasonable. It is not necessary to overwrite the old scientific record before every independent unit test, source inspection, or algebraic correction can proceed.

Also, an existing output directory is not automatically a verified reference for the current revision. Some inspected provenance headers identify an older dirty revision. Associate each reference with input contents, source identity, toolchain, enabled physics, stopping reason, and test purpose. Preserve the old baseline separately if its numerical history is useful.

The plan's command `make check` does not select affected cases: the current Makefile invokes `run_check.sh check` without arguments. To run a bounded case set, use the script's actual case arguments or add an explicit filtering interface when test tooling is changed.

Finally, numerical agreement within 0.1% is not identity. Distinguish byte equality, numerical equality, and agreement within a declared tolerance. A physical defect cannot be accepted on the basis of a small output difference. The claim that A7 necessarily moves every reference output is also unverified and should be replaced by an affected-path list and measurements.

### R11. Some necessary physical closures disappear from the implementation phases

**Priority: medium; essential to the stated scope of a physically correct molecular code.**

The plan includes many audit items, but several are reduced to brief mentions without an owner or a closure test:

- Species enthalpy diffusion and, when charged transport is included, the associated electric-field/current and energy closure. Zero net diffusive mass alone does not establish energy conservation.
- Momentum and thermal coupling validity, including a supported response when a single-temperature fluid approximation becomes inappropriate.
- Lower/upper atmosphere matching of species/element fluxes and energy, including a matching-pressure sensitivity test.
- Consistency of molecular composition and electron/temperature fields in `_adv` post-processing and observables.
- Optical-depth-dependent photon conservation and overlap between continuum absorption and the LW/oxygen bands.

**Revision:** assign these to explicit stages or declare a restricted supported regime that excludes the unresolved physics. They need not all become large new implementations immediately, but their validity checks cannot be omitted from a claim of physical completeness.

The atomic transport decision should distinguish an implementation schedule from physical necessity: where atomic ionization cannot equilibrate on the flow timescale, turning transport off is a stated approximation, not an equivalent physical option. Similarly, labeling an oxygen energy budget incomplete is useful disclosure, but does not satisfy Phase 1's physical closure criterion for oxygen-rich production calculations.

### R12. Parallel independence and session estimates are overstated

**Priority: planning. References: plan lines 65-67 and 410-415.**

The reaction ledger, source integration, residual evaluation, and species-flux design share equations and interfaces. Their implementation should not proceed as independent mutations before those interfaces are agreed. In contrast, LLF/ROE formula checks, coverage inspection, rate-range inventory, and baseline comparisons are bounded independent tasks.

The session estimates have no measured implementation basis. A6 includes a complete oxygen energy budget and eliminated excited states; B2 may require a new physical fit or level calculation. Neither is reliably a one-session task from the evidence given.

**Revision:** estimate work after each derivation/prototype has a concrete scope. Use milestones with acceptance tests and identify uncertainty in data availability and nonlinear convergence. Avoid making a calendar estimate the basis for retaining incomplete physics.

## 5. Recommended revised order

The dependency graph should follow the equations that downstream work consumes:

1. **Establish evidence and focused gates.** Preserve baseline provenance; convert the existing diagnostics into tests; add the missing H+ configurations. This does not require replacing every old reference first.
2. **Implement bounded corrections.** Correct LLF/ROE, carrier reference scales, and returned-state diagnostics. Define safe behavior for a failed carrier solve and an active CO ceiling.
3. **Derive the common physical system.** Specify species variables, energy convention, reaction/photon ledgers, common face transport, cell volume, and boundary constraints. A thermal-energy formulation remains a valid option.
4. **Implement thermodynamics and the coupled source update.** Coordinate A4, A6, A7, B1, and B3. Require closed-cell conservation and temporal accuracy before accepting the energy correction as complete. Correct H3+ asymptotes and continuity with separate physical and numerical gates.
5. **Complete the spatial operator and its validity conditions.** Establish consistent hydro/species fluxes, diffusion and enthalpy transport, and the intended radiation closure. A staged implementation can start in a restricted H/He regime.
6. **Extend the steady solver to that operator.** Add H+ and required additional species, composition/radiation convergence, derivative checks, and final-state residual acceptance. Expand atomic transport support according to the physical regime.
7. **Validate coupled atmospheres and extend composition.** Test boundaries, matching pressure, spatial convergence, molecular benchmarks, and observables. Add broader C/O/N chemistry within explicitly verified thermodynamic and transport ranges.

This order moves the design of D1/D3 ahead of C1, and treats B3 as part of completing A4 rather than an optional later improvement. It does not require all molecular extensions to be finished before a restricted H/He system can be validated.

## 6. Minimum revisions before using the plan as an implementation specification

The following edits are the minimum necessary:

1. Correct the description of `eval_residual` and remove the unsupported inference that atomic converged results must be preserved.
2. Replace A4's informal energy recipe with a species-energy derivation, making B3's necessary coupling a prerequisite for closure.
3. Make A7's residual, rate evaluations, derivative, and acceptance criterion describe one mathematical problem.
4. Replace B2's unconstrained two-column power law with a verified linear collisional asymptote and stronger continuity tests.
5. Add actual H+ transport tests and correct the stated coverage of `mol_sec_ion` and `mol_carrier`.
6. Move transport/boundary equation design before steady-system implementation and resolve the CO-ceiling/rejection conflict.
7. Replace the claim of an exact JFNK action with the actual finite-difference and inner-closure requirements.
8. Assign explicit validity or implementation gates to mixture energy transport, lower-atmosphere matching, and molecular post-processing.

The overall direction is sound after these revisions. At present, the plan is best treated as a task inventory with promising priorities; its energy and solver prescriptions still need mathematical consolidation before they can guide a physically reliable implementation.
