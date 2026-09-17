# Review of PLAN_20260916.md

Date: 2026-09-16 (KST)

## 1. Overall judgment

The plan identifies important remaining work and generally separates diagnostics from implementation. However, it is **not ready to implement unchanged**. The main revisions needed are:

1. Correct the conservation design of L22 before implementing regional carrier intervals.
2. Treat L25 as a discrete-operator and boundary investigation, not as permission to replace the momentum certification equation.
3. Investigate the newly revised wind-window boundary, including its zero-flux discontinuity and nonlocal transient behavior.
4. Keep perturbed residual evaluations inside the L15 investigation; the frozen-state test does not exclude them.
5. Specify the energy recipients and missing molecular data before promising a level-resolved R15 cascade.
6. Correct the inherited L12 eddy-flux decomposition and validate the new transport equations independently of the approximate post-process.
7. Expand L9 beyond a runner input change: the current evaluation route stops before the advection post-process and modifies chemistry working arrays.

I recommend authorizing bounded diagnostics and small metadata repairs after their contracts are clarified. The regional transport design, physical boundary/flux changes, and formation-state model should each have a revised design gate before implementation.

### Verification scope

This review read the plan, relevant current source and callers, the L15 and L7e diagnostic records, the L12 design, and selected literature. It did not rebuild, run a simulation, replace expected regression outputs, modify production code, or reproduce the campaign. Numerical campaign results quoted below remain reported results, not new measurements.

Direct file hashing found `EXHALE.x` MD5 `c2e9c9990b9f14f1be8cd77abca68945`, matching the plan. The inspected main program, base boundary, and carrier module also match their entries in `LHS1140b/models/BINARY_MANIFEST_c2e9c9990b9f.txt`. This checks those files, not every linked object or the full build environment. HEAD is `43bc28cef58772560bae019559f3592335922212` with substantial uncommitted work.

The reported 17/17 regression result was not rerun. Agreement with refreshed expected files is useful reproducibility evidence, not independent validation of the underlying physics.

## 2. Decisions by item

| Item | Recommendation | Required change to the plan |
|---|---|---|
| L23 | Proceed after metadata semantics are specified | Preserve the imported certification claim separately from newly evaluated certification |
| L15 | Proceed with revised localization | Include perturbed residuals, histories, and the failing trajectory, not only Krylov arithmetic |
| L22 | Proceed with diagnostics; redesign implementation | Distinguish local pseudo-time preconditioning from conservative physical subcycling; retain nonlinear acceptance |
| L24 | Proceed with controlled experiments | Preserve the original fixture; add a successful-start variant instead of replacing its historical purpose |
| L25 | Rewrite before adopting a treatment | Audit the new boundary; preserve the full physical residual; assess low-Mach consistency rather than unchanged mass loss |
| L7g | Proceed first as a data/model study | Separate prompt translation from internal excitation; do not reuse grain formation or v=1 data as general state-resolved data |
| L12 | Raise its scientific priority | Correct flux partitioning, test full coupling, and define the inner-domain approximation |
| L9 | Proceed as a distinct output-path change | Explicitly invoke the post-process and protect the loaded state; test actual downstream files |
| L16/naming | Keep bounded and separate | Do not infer a historical cause or regenerate data for a name-only change |

## 3. L23: a real missing field, but not a one-line certification fix

### What the implementation shows

`src/modules/files_IO/load_IC.f90:1770–1829`, `parse_coupling_header`, reads `certified` into `ic_certified` but has no `cert_reason` case. This supports the plan's finding.

However, `src/modules/functions/utilities.f90:31–38`, `set_state_certified`, sets **both** the Boolean and the reason. Calling `set_state_certified(.false., reason)` indiscriminately while parsing a reason is not a neutral metadata restoration. Moreover, `src/modules/time_step/certification.f90:1358–1369` can assign a qualification token to a **true** certification, not just a reason to a false one.

### Recommended contract

Read the file's Boolean and reason into imported-state metadata first. Keep these distinct from the certification of the state under the current executable and configuration. A raw roundtrip may preserve the imported pair; an actual evaluation must write its newly determined pair. Do not let token order determine the result, and do not preserve a stale success qualification after the state or model changes.

Tests should cover true and false claims, absent/unknown reasons, qualification tokens, reordered header fields, and re-evaluation that changes the verdict. Check reason length against the current 32-character destination. The existing test is `src/tests/grid_and_gates/restart_intent_and_metadata.sh`, not a separate directory with that name.

A metadata-only change warrants these focused tests and the roundtrip script. A full physics matrix is unnecessary unless the implementation alters state handling beyond metadata. Replacing a stored roundtrip header and saying that no expected output changes should not be presented as the same statement.

## 4. L15: the stage-1 conclusion is too broad

The stage-1 memo evaluates one certified atomic state repeatedly with a particular private executable and compares formatted files and logs. Its limitations section is more careful than the plan's assertion that the residual evaluation and ionization sweep are excluded.

### Why this does not localize the failure to Krylov arithmetic

`src/modules/time_step/steady_newton.f90:10186` defines `jv_product`, which samples admissible perturbations of the state. `pgmres` around lines 11820–11925 applies the banded preconditioner, calls `jacobian_action_of_direction`, and then performs Arnoldi operations.

Even if repeated `F(U*)` evaluations agree at a certified state, `F(U + epsilon v)` can exercise different chemistry branches, active bounds, failed-root handling, or mutable state. The failing `atomic_elem_newton` trajectory is also different from that certified fiducial. Finite-difference cancellation can amplify a small residual discrepancy by a factor proportional to `1/epsilon`.

The shown Arnoldi inner products are ordinary Fortran `sum` expressions inside a serial loop. The inspected `steady_newton.f90` contains no explicit OpenMP directives. An intrinsic `sum` is not, by itself, evidence of an OpenMP reduction. Called routines and linked libraries still need inspection.

Finally, identical formatted state files do not establish bitwise equality of all internal arrays or the complete residual vector. Values can differ below the output precision or in quantities not written.

### Revised experiment

1. Reproduce two eight-thread trajectories from the actual failing seed, with the same executable and BLAS settings. Locate the first divergent iteration with compact array checksums or exact binary dumps.
2. At the preceding checkpoint, record `U`, composition, scaling, active constraints, cached radiation/chemistry state, preconditioner data, and the pseudo-time shift.
3. Repeat both `F(U)` and the actual perturbed evaluations used by the Jacobian action. Restore state between evaluations; also test changed evaluation order to detect hidden history dependence.
4. Compare preconditioner inputs/outputs, each Jacobian action, then Arnoldi coefficients and the merit calculation. Dump full basis vectors only around the first divergence.
5. Include at least one molecular configuration after the atomic cause is understood.

A barrier alone cannot repair a race in which threads overwrite shared scratch storage while needing distinct values. Ownership or local storage must be corrected. Likewise, adding thread partial sums in thread order does not generally give the same result at one and eight threads because the partitions differ. Thread-count-independent arithmetic requires a fixed partition/reduction tree independent of thread count, or another explicitly reproducible summation method.

The proposed 200-iteration repetition is useful for the fixture that first diverged around iteration 21. Starting ten cases at their already certified states may exercise almost no Krylov work and therefore is not a sufficient stress test. Include difficult saved iterates or controlled admissible perturbations. Printed-digit agreement is an output criterion, not proof that a race is absent.

Annotating affected reproducibility records is sensible, but the wording should be precise: the solve path was not demonstrated repeatable in the stated configuration. Do not attribute a race or reduction mechanism before localization.

## 5. L22: the proposed regional scheme is underspecified

### 5.1 Stronger evidence, but still a causal experiment is needed

The updated L7e section 23 provides row decomposition and a trajectory ending repeatedly on the movement bound. This supports prioritizing the controller. It does not prove that removing that restriction will produce a jointly admissible solution. The planned bound-exclusion experiment is useful precisely because that causal question remains open.

The current carrier Jacobian now includes copied-ghost derivatives and removes the extra carrier-mass division in the advective coefficients (`diffusive_photochemistry.f90`, assembly around lines 4700–4770). These older defects should not be listed as unfixed causes. The updated memo also states that the failing cases remain after the corrections.

The claim that four cases are held by the same bound should distinguish directly failing trajectories from `wellmixed/HeH2.13`, which is listed as lacking a seed. Absence of a seed is not a measured bound failure.

### 5.2 A maximum bound and bounds in every cell are equivalent

For nonnegative changes `d_j`,

`max_j d_j <= trust` if and only if `d_j <= trust` for every cell.

Thus the defect is not the mathematical use of a maximum. It is the use of one update multiplier and one cumulative excursion budget for regions with different responses. Merely moving the maximum into separate regions does not enlarge the feasible set if every cell must still remain inside the same entry-referenced bound. It can only improve the path within that set.

The description of wind response as local also needs qualification. The equation of state is local, but hydrodynamic pressure coupling and attenuation/shielding are not generally local. At fixed conserved thermal energy, particle count is not even a complete pressure-response measure: for an ideal monatomic mixture, `p = (gamma-1)e_th` is independent of particle count while temperature changes. Molecular heat capacities and reaction heating introduce further composition dependence.

The current `pressure_and_temperature_at_fixed_conserved_state` at lines 5580–5603 uses the caloric EOS at fixed `u`; this is the relevant implementation to differentiate. Monitor changes in pressure, temperature, sources, and the coupled residual in addition to particle count.

### 5.3 Local pseudo-time is not automatically nonconservative

There are two distinct algorithms:

**Steady local pseudo-time preconditioning.** Solve `M_j delta U_j / delta_tau_j + R_j(U) = 0`, where `R` retains one common physical face flux. Different positive pseudo-time intervals do not change the root `R(U)=0`. Such iterations need not conserve a physical-time inventory at every intermediate update, because they are not a synchronized physical evolution. They must preserve admissibility, and the final unmodified conservative residual must vanish.

**Conservative physical subcycling.** Neighboring cells must exchange the same time-integrated interface flux with opposite signs over a common synchronization interval. This needs accumulated flux accounting and appropriate source integration.

The plan mixes these contracts. Taking the shorter interval at an interface is insufficient if adjacent cells still multiply the same instantaneous flux by different cell intervals. Their inventory changes are then `-dt_L A F` and `+dt_R A F`, which do not cancel. Conversely, replacing the physical face flux by a face-dependent interval-weighted flux inside the steady balance can change the equation being solved unless the preconditioning and normalization are derived consistently.

Choose and document one contract before introducing a region table. For a stationary solver, cell-local pseudo-time preconditioning of the unchanged residual is worth testing before adding permanent front/wind categories. For actual subcycling, define the common interval and accumulated interface transfers explicitly.

The `x2 > 1e-2` division is empirical and can have no crossing, multiple crossings, or move between iterations. The diagnostic record itself reports a case whose H2 fraction never falls below that value. Two regions are therefore not a general physics-based decomposition. Freeze any diagnostic partition during its controlled comparison and test sensitivity to the partition separately.

### 5.4 Coupling does not eliminate the need for globalization

Candidate b2 says that adding a linearized wind response removes the need for a bound. This is too strong. If hydro variables remain fixed, a local pressure derivative is not a full wind response. A coupled Newton block or Schur complement must include the relevant hydro/carrier derivatives. Even a fully coupled formulation still needs positivity, feasible element/charge fractions, and nonlinear step acceptance.

The existing outer progress logic still compares `comp_drift_last` with `comp_drift_before` (`EXHALE_main.f90:6629–6633,6809–6810`). A reduction caused by a smaller movement limit can still masquerade as contraction. Fix that diagnostic/control issue within L22 rather than only changing the spatial intervals.

Acceptance should include residual decrease under fixed scales, elemental and charge balances, admissibility, and independence from iteration controls once converged. Agreement between two solutions is useful only under the same physical operator and identified branch. A residual tolerance is not automatically an error bound on state variables. Specify a normalized conservation metric for the proposed `1e-14` threshold; a dimensionless tolerance cannot be applied to an unspecified dimensional inventory.

## 6. L25: revise the physical and numerical premises

### 6.1 Low Mach number does not mean zero wind

A base with Mach number `1e-8` may still supply a finite steady mass flux because its density is large. Hydrostatic balance is a leading-order momentum balance, not a license to drop the finite-flow equations or their certification.

The expectation that a base-local change must move mass loss by less than `1e-4` because the sonic point determines it is unjustified. Reservoir entropy, heating, and the subsonic structure help select the transonic solution. If the sonic point is outside the computational domain, its regularity condition is not automatically enforced by the truncated-domain calculation.

### 6.2 Candidate (a) is not yet an all-speed flux

`Num_Fluxes.f90:140–191` computes HLLC star states from `S_star` and selects the star flux by its sign. Replacing a wave speed, averaging states, averaging star fluxes, and using a Rusanov flux are different operations. The proposal does not define which conservative flux would result.

The existing Rusanov routine around lines 368–405 uses an acoustic signal speed `max(|v|+c)`. Smoothing a switch does not by itself remove low-Mach acoustic-scale dissipation; the `abs`, `min`, and `max` operations also need attention when differentiability is claimed. A smooth central blend is not automatically stable, positivity-preserving, well-balanced, or asymptotically correct.

Require an explicit flux formula, its low-Mach scaling, preservation of hydrostatic balance, boundary compatibility, and tests of stationary contacts, a resolved low-speed flow, and a compressible control. Define the spatial extent by a resolution-aware criterion rather than cells 1–20. Report solution movement after judging physical correctness; unchanged old mass loss is not the acceptance rule.

### 6.3 Candidate (b) should inspect the existing implementation first

`src/modules/flux/low_mach_dissipation.f90` already defines `contact_mode_dissipation_flux` and a Mach-dependent artificial stress. It pairs momentum and energy flux contributions and is invoked by `RK_rhs.f90:165–166`; the production stationary residual also calls `RK_rhs` (`steady_residual.f90:151–172`). The separate generic/quadruple diagnostic route in `hydrodynamic_rows.f90:297–301` explicitly refuses active low-Mach damping, so that route must not be assumed to include it. The plan should compare this existing option before inventing a separate one-cell viscosity.

There is an additional caution in that module. The flux uses a third difference of velocity multiplied by a spatially varying density, wave speed, and gate (lines 222–232). The header's nonnegative integrated heating argument is not established for arbitrary varying coefficients and the nonuniform grid by merely dropping boundary terms. Variable-coefficient integration introduces extra terms. Conservative total-energy pairing alone does not prove a nonnegative kinetic-energy dissipation operator.

Before relying on this option, check its discrete energy quadratic form with the actual volume weights, grid, and gate. A construction such as a properly weighted negative-semidefinite second-difference composition provides a clearer route to a damping proof. This review identifies an unsupported guarantee, not a measured unstable run.

### 6.4 Reject candidate (c) as an acceptance shortcut

An algebraically equivalent, well-balanced evaluation of the **full** momentum residual is appropriate. Certifying only the hydrostatic part and omitting finite-flow terms is not. Recalibrating a tolerance does not repair a change in the governing equation. Also, `|v| c` has units of velocity squared and cannot be compared directly with a pressure-gradient truncation error without a specified normalization.

If a distinct hydrostatic reduced model is wanted, state its domain, matching conditions, approximation error, and equations explicitly. It must not be introduced by changing only the certification row of a wind calculation. The reported failures include the mass row, so a momentum-only certification change does not even target every stated symptom.

### 6.5 Update after reading the supplied Rieper papers

The user supplied three published Rieper PDFs during this review. The 2011 paper's methods, sections 2 and 3, relevant numerical-test descriptions in section 4, and conclusion in section 5 were inspected. The 2010 paper's characteristic analysis, sections 1–4, was also inspected. The 2009 file was identified, but its detailed derivation was not reviewed and is not used as evidence here.

**The actual 2011 correction is different from candidate (a).** Section 3.1, equations (3:15)–(3:16), scales the normal velocity jump entering the Roe dissipation by `min(Ma_face,1)`. In the paper's two-dimensional definition, `Ma_face = (abs(U_Roe) + abs(V_Roe))/a_Roe`. Eigenvalues and eigenvectors are not replaced by an averaged HLLC star state. Section 3.2 assumes the local and reference Mach numbers have comparable asymptotic order. This is a controlled change to numerical dissipation, not simply smoothing the contact-selection branch.

There is a direct code location for a bounded Roe-based diagnostic: `Num_Fluxes.f90:320–330` forms `dvel = vR-vL` and the acoustic coefficients

```text
a1 = (dp - rho_avg*a_avg*dvel)/(2*a_avg**2)
a3 = (dp + rho_avg*a_avg*dvel)/(2*a_avg**2)
```

The velocity-jump contribution is the candidate to examine, while preserving the physical central flux, the pressure-departure treatment `dp_wb`, and the existing admissibility logic. This is a proposed experiment, not a completed derivation for EXHALE. A one-dimensional analogue of the Mach factor is straightforward to write, but its behavior when neighboring velocities cancel in their Roe average deserves an explicit probe.

**Do not claim a convergence cure from an accuracy correction.** Section 5 expressly distinguishes the improved low-Mach accuracy from stiffness: explicit timesteps remain acoustically limited, and implicit algorithms can still converge slowly. LMRoe also contains absolute values and a minimum, so the paper is not a proof of globally C1 residuals for JFNK. L25 should separately measure spatial accuracy, Jacobian-action regularity, and nonlinear convergence.

**The papers do not justify replacing HLLC with Rusanov.** Rieper (2010), section 3.2, derives excessive HLL dissipation of contact/shear modes relative to slow transport; section 4 gives characteristic-wave tests to diagnose a flux. Rieper (2011), section 5, explicitly distinguishes schemes such as Rusanov and HLL, whose additional entropy/shear diffusion is not cured by this velocity-jump correction. The latter paper identifies HLLC as a possible extension, but its demonstrated method and analysis use Roe. An HLLC adaptation therefore needs its own derivation and validation.

**The gravitational problem needs an additional consistency argument.** Rieper's 2011 asymptotics use the source-free Euler equations, with spatially uniform leading pressure, and the reported low-Mach tests use an ideal gas with fixed gamma. EXHALE instead has a stratified hydrostatic pressure gradient, gravity, spherical geometry, radiation/chemistry, and a composition-dependent caloric EOS. Apply the low-Mach analysis to departures from the appropriate gravitational equilibrium, not by flattening the physical background pressure. The multidimensional pressure-checkerboard results also do not establish that EXHALE's radial alternating-velocity mode has the same cause.

**Revised practical recommendation:** first characterize the new boundary and the existing damping option; then compare a precisely defined Roe velocity-jump correction on controlled fixed-background tests, followed by gravity-balanced and chemically coupled tests. Use the Roe experiment to distinguish dissipation error from branch/conditioning problems. Adopt an HLLC variant only after deriving its wave contributions. Neither paper supports changing the physical certification equation or prescribing that corrected mass loss must remain unchanged.

## 7. Additional prerequisite: audit the new base-window boundary

The previous review's route-dependent boundary split has been removed. The current `base_boundary.f90` uses the same wind-window blend without `base_branch_stationary_route`. This is a real implementation change and should not be described as the old split still being present.

However, the replacement deserves a new item in the plan:

- `wind_window_mass_flux` (lines 484–526) measures `(max F - min F)/abs(mean F)` over a distant window and rejects exactly zero mean flux.
- The default blend (lines 660–677) trusts a sufficiently uniform window regardless of the absolute flux magnitude.
- The header acknowledges that a transient base can have a different flux from that distant window. A uniform distant flux does not establish the direction of characteristics at the base during such a transient.

### A concrete continuity counterexample

Hold the local base state fixed with local branch weight `w_i != 1/2`. Set all window fluxes to the same nonzero value `epsilon`. Their relative spread is zero, so the window receives full weight. As `epsilon` tends to zero, the wind Mach tends to zero and the smoothstep branch weight tends to `1/2`. At exactly zero, `have_F` becomes false and the branch weight returns to `w_i`.

Thus the default boundary is discontinuous at that class of states. The code comment acknowledges the zero-window exception, but calling it a measure-zero case does not make it harmless to a Newton perturbation or a hydrostatic initial condition. The resulting density jump is nonzero when the two entropy constructions differ.

In addition, a cubic smoothstep of the relative spread is not generally C1 as a function of the full state: `max` and `min` have kinks when the cells attaining the extrema change inside the active blending interval. The weight is smooth in its scalar argument, not necessarily in `U`.

Add frozen-state boundary continuity and directional-derivative probes to L25 **before** attributing every kink to HLLC. Test constant-window flux approaching zero from both signs, local/window sign disagreement, extrema switching, and transients. The durable objective remains a boundary determined by the physical reservoir and local characteristics; agreement between solver routes is necessary but not sufficient for that physical contract.

## 8. L7g: the formation-state proposal needs different inputs

The present R15 term multiplies the full reaction energy by `h2_vibrational_heat_fraction` (`molecular_reaction_heat.f90:469–496`). That function uses one effective radiative decay rate and H/H2 v=1 quenching fits (`h2_vibrational_relaxation.f90:189–230`). It is not a table of rates for each vibrational/rotational transition.

The published Hollenbach and McKee (1979) text was checked locally, including its chemical-heating treatment near equations 6.42–6.43. That treatment distinguishes grain formation from the H-minus route and assigns energy among excitation and translation. It is not a state distribution for gas-phase `H + H + M` association. A grain-formation distribution cannot be transferred to R15 without a physical derivation.

For three-body formation, write an energy partition such as

`D0 = E_translation + E_internal(H2) + E_internal(M)`

with any change in reactant/product kinetic energies handled consistently. A reduced thermal deposition can then be written as prompt thermalized translation plus the collisionally thermalized share of internal excitation. Multiplying the **whole** binding energy by a vibrational branching fraction incorrectly suppresses any prompt translation in the low-density limit. Formation probabilities, third-body identity, destruction, and radiative escape matter.

The IR module contains level energies and radiative line information (`molecular_infrared_data.f90:44–51`). Existence of those arrays does not establish a complete connected transition network, state-resolved collisional rates, or a cascade solver. Before committing to the proposed three-to-four-day implementation, inventory coverage and acquire justified formation and collision data. Helium collisions matter especially for helium-rich configurations; the current effective routine explicitly omits them.

There is also a concrete error in the plan's proposed exemption for dissociative recombination. The source reaction table at `molecular_reaction_heat.f90:261–267` contains **R6: H3+ + e -> H2 + H**. Its products are not all atoms, and H2 can carry internal excitation. Audit energy recipients channel by channel rather than categorizing R5/R6/R7/R16 together as fragments without excited molecules. This does not prescribe a particular branching fraction without data.

Recommended acceptance tests include zero- and high-collision limits, normalized product distributions, event-level energy conservation, appropriate equilibrium/detailed-balance checks where the adopted network supports them, and a check that IR cooling and the caloric EOS do not count the same energy twice. Three density samples alone cannot validate a level model. If data are inadequate, retain an explicitly bounded approximation and report uncertainty rather than inventing a detailed cascade from v=1 coefficients.

## 9. L12: correct the inherited flux equation and raise priority

Transporting ion fractions on their element fluxes is a sensible direction. Current atomic ionization systems already call `impose_transported_ionization_fractions`, for example `System_HeH_TR.f90:79`; the remaining work is not the invention of all imposition infrastructure. It is the invoked atomic transport path, conservative flux partition, coupled state update, and certification.

### Eddy-flux double counting in the supporting design

The plan relies on `lhs1140b_stationary_L12a_design_20260913.md`, section 1.3. That design uses the full element flux, including diffusion, for transport of the stage fraction, then adds `-n_tot Kzz grad(x n_el/n_tot)` as the stage eddy term. The implemented element flux already contains eddy mixing: `binary_element_diffusion.f90:2009` includes `Df + Kf`.

Let `y = n_el/n_tot` and `x` be a stage fraction. In the number-fraction convention used by that design,

`J_el,eddy = -n_tot K grad(y)`.

The full stage eddy flux is

`J_stage,eddy = -n_tot K grad(x y) = x J_el,eddy - n_el K grad(x)`.

If `x J_el` already includes `x J_el,eddy`, the additional relative flux must be `-n_el K grad(x)`, not the complete product gradient again. Otherwise summing stage fluxes double-counts the elemental eddy contribution. The actual implementation must make the same derivation using its own mass/number convention; this is an error in the proposed decomposition, not evidence that the unimplemented stages already run with it.

A mandatory test is `sum_stages F_stage = F_element` on the actual discrete faces, with a nonuniform element abundance, both constant and varying stage fractions, and nonzero Kzz. Also test charge consistency and the helium simplex. Neglect of drift between charge stages remains a physical approximation that needs its collisional/electric-field regime stated.

### Acceptance and sequencing

The `_adv` post-process solves an approximate conditional problem, not the same equations as the proposed coupled transport. Agreement with it is a comparison, not a correctness oracle. The plan's stage-B wording mixes an ion-profile comparison with a `1e-2` mass-row tolerance; define independent composition, transport-residual, and hydro criteria.

First validate fixed-background transport with known limits or manufactured balances. Then recouple opacity, electron density, heating/cooling, pressure, and hydrodynamics. Accept the fully coupled residual, not a profile held to the old approximate answer. Ion transport can change the hydro state, so later residual and observable movement may be physically required.

Do not hide unsatisfied inner transport rows solely because a radius gate excludes them. If an inner equilibrium approximation is retained, validate its timescale regime and flux matching across the transition. A certification limited to the outer wind must be described as such.

Given the plan's own reported discrepancy between equilibrium and `_adv` line predictions, L12 is a prerequisite for a self-consistent scientific interpretation of the atomic observables, not merely the seventh cleanup item. Its operator design and data contracts can proceed alongside single-threaded L22/L25 diagnostics. L15 should block claims of eight-thread repeatability, not all single-threaded scientific progress.

## 10. L9: evaluation and post-processing are not the same entry point

The runner currently takes a tiny step in `LHS1140b/models/run_case.sh:485–514`, then extracts the steady mass-loss text and uses `_adv` output. `EXHALE_transit.py:79–80` reads `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` explicitly.

In contrast, the `stationary_evaluate_only` branch (`EXHALE_main.f90:5331–5355`) writes equilibrium output and breakdowns and then stops. The ordinary `post_process_adv` call is around line 4238, and the mass-loss report is later in that ordinary output path. Merely changing `Restart intent` will not establish that these downstream products are generated.

Also, the evaluation path calls `ioniz_eq` on `f_sp` and reconstructs diagnostics after that sweep. It takes no time step but is not automatically a byte-preserving evaluation of every loaded composition. The explicit `stationary_equilibrate_loaded` branch controls additional equilibration there; `EXHALE_RELOAD_EQ` is not a substitute for tracing that specific path.

Define two operations: evaluation on copied work state, with the original retained, and optional derived advection post-processing. Expose or call the existing post-process explicitly without reentering time integration. Do not rely on stale `_adv` files left from earlier runs.

Tests must verify the actual required files, mass-loss value, certification metadata, and unchanged input state. A stale-output sentinel or an empty test output directory can reveal a skipped producer. Comparing outputs only to the certification tolerance is insufficient: a residual tolerance does not directly bound an observable or a file's state error.

## 11. L24, L16, and acceptance rules

For L24, a single-threaded pseudo-time-start comparison is appropriate. Failure at an explicit-scale interval after 200 iterations does not prove a particular startup defect; preconditioner quality, admissibility, state-dependent residual branches, and coupling remain possible. Preserve the original difficult fixture and add a named large-initial-interval variant if useful. Replacing the original fixture would lose coverage of the behavior it was designed to expose.

For L16, annotate unreconstructable historical numbers honestly. Later localization of a threading problem does not retrospectively prove that it caused every historical discrepancy.

For directory renaming, update paths and verify resolution without regenerating unchanged numerical products. Names should express the physical configuration or comparison quantity, not only a permanent control/measured role.

Replace the blanket default-off rule with two categories: experimental algorithms may remain optional while evaluated; confirmed corrections to physical equations should become the supported behavior, with old behavior retained only for an explicitly identified diagnostic need. Conversely, enabling a method because three catalog cases converge is not sufficient if its conservation or low-Mach consistency is unproved.

Byte comparisons are valuable for untouched numerical paths and deterministic fixtures. They should not force preservation of a physically incorrect answer. Tolerances, conservation checks, solution uncertainty, and regression snapshots have different purposes and should be stated separately.

## 12. Revised order and evidence to retain

1. Record the exact inspected build and make the small metadata contract/test changes separately.
2. Localize L15 along the failing trajectory while continuing other controls at one thread.
3. Add zero-flux and extrema-switch probes for the new boundary; correct progress diagnostics before expanding relaxation complexity.
4. Run the bounded L22 causal comparisons and L24 initial-interval comparison. Choose the pseudo-time versus physical-subcycling contract from those results.
5. Develop L25 with the full conservative residual unchanged for certification; assess existing damping and a defined low-Mach flux separately.
6. Begin the L12 flux derivation and R15 data inventory without waiting for every catalog case to solve. Their implementation estimates depend on what these steps establish.
7. Complete a true no-step evaluation/post-processing path with downstream product tests.
8. After each physically coherent change series, run its affected tests and report movement before updating expected outputs. Do not routinely rerun unrelated broad tests.

Retain future diagnostic programs under `docs/audit_20260905/` with input state, executable identity, thread configuration, residual definitions, and raw outputs. This review created no diagnostic program and changed no numerical products.

## 13. Sources and limits

The primary evidence is the implementation and paths cited above. The plan's SHA-256 at review was `d1520b473e5796f82e3877ef6e393990d19b2c7d2c7daf9a43850f09ea2f78f2`. Line references can move with subsequent development.

Local published literature inspected for the formation-energy discussion: Hollenbach and McKee (1979), ApJS 41, 555, `../../references/Hollenbach_1979_ApJS_41_555.pdf`, particularly the formation/chemical-heating treatment around equations 6.42–6.43. Its reaction mechanisms do not supply a three-body product distribution for R15.

An ADS API lookup was attempted first, but the execution environment could not resolve the ADS host. A subsequent publisher search verified Rieper (2011), *A low-Mach number fix for Roe's approximate Riemann solver*, JCP 230, 5263–5287, [publisher record](https://www.sciencedirect.com/science/article/pii/S0021999111001689). The user then supplied the published PDFs, which were read locally with `pdftotext -layout`; section 6.5 supersedes the initial literature-access limitation.

- Rieper (2011): `../../references/Rieper_2011JCP_230_5263.pdf`, DOI `10.1016/j.jcp.2011.03.025`. Methods and scope checked as described in section 6.5.
- Rieper (2010), *On the dissipation mechanism of upwind-schemes in the low Mach number regime: A comparison between Roe and HLL*, JCP 229, 221–232: `../../references/Rieper_2010JCP_229_221.pdf`. Characteristic analysis in sections 1–4 checked.
- Rieper and Bader (2009), JCP 228, 2918–2933: `../../references/Rieper_2009JCP_228_2918.pdf`. Publication identified only; not used to claim a verified result about EXHALE's grid.

No new state-resolved formation or collision dataset was verified. No new race, failing trajectory, or numerical convergence rate was experimentally established here. Confirmed implementation facts, analytical counterexamples, reported prior measurements, and proposed tests have been kept distinct throughout.
