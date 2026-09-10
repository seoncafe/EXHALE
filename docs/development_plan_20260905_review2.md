# Second review of the EXHALE development plan

Review date: 2026-09-05 (KST).

Reviewed document: [development_plan_20260905_rev1.md](development_plan_20260905_rev1.md).
Previous review: [development_plan_20260905_review.md](development_plan_20260905_review.md).
Implementation inspected: `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`.

## 1. Assessment

Revision 1 is substantially more reasonable than the original plan. Its revised order, governing-equation derivation, species-energy coupling, thermodynamic normalization, and independent convergence gates should be retained. The document now distinguishes several previously conflated operations and withdraws unsupported predictions about atmospheric results.

It still requires corrections before it can serve as an implementation specification. The most serious new design problem is that its explicit energy residual combines a chemical-energy increment with heating from which that chemical energy has already been subtracted. An isolated photoionization event exposes the resulting double subtraction. The CO constraint is also insufficiently defined, and the charged diffusion and thermal-energy flux conventions remain incomplete.

Reading the published papers in the workspace's `references/` directory revealed a separate, confirmed source error that both the original audit and the first review missed: the H2 channel module interprets Chung et al.'s measured detector-signal ratio as a proton-number ratio. Their instrument records one pulse for an unresolved double-ionization event. The distinction changes the reconstructed channel probabilities while leaving element, charge, and total-cross-section checks satisfied. This correction must be added to A5 and the E1 tests.

A new direct execution also found a nonvacuum rarefaction for which the current ROE routine returns NaN sound speeds. The two algebraic changes listed in A2 do not reach that branch. Thus even Phase 1 needs a more complete acceptance specification.

Recommended disposition: retain the phase structure, revise the issues below, and use the corrected D0 derivation to settle the shared equations. A1 and the reference-scale portion of A3 can proceed as bounded work. A2 needs the additional branch tests, and the CO model should leave the category of corrections requiring no derivation.

This review adds this document and preserves the additional diagnostic sources and execution instructions in `docs/audit_20260905/`, as requested. Production source, input files, the reviewed plan, and simulation products were not modified. No atmosphere or full regression matrix was run.

## 2. Evidence and verification scope

The complete revision 1 plan was read and compared with the first review. Source inspection followed the marching and steady call paths, the EOS conversions, the carrier residual and limiter, H2 photoabsorption and heating, H3+ cooling, relevant initialization and input restrictions, and regression tooling.

The following published papers were inspected from the existing local collection. Reading covered the methods and the specific equations, tables, or discussion supporting the findings; this was not a claim to have audited every result in every paper.

| Published source in the local collection | Material checked and its role |
|---|---|
| [Chung et al. 1993](../../references/Chung_1993JCP_99_885.pdf), JCP 99, 885 | Experiment; Results and Discussion on pp. 886-888; detector response, Eqs. (1)-(5), Table II, double-ionization fraction, excited neutral fragments. The relevant page was also rendered and visually inspected. |
| [Yan et al. 1998](../../references/Yan_1998_ApJ_496_1044.pdf), ApJ 496, 1044 | Section 4; the vertical double-ionization threshold and the double-to-single ionization ratio. |
| [Miller et al. 2013](../../references/Miller_2013_JPCA_117_9770.pdf), JPCA 117, 9770 | Cooling-curve calculation and non-LTE method; Eqs. (6)-(10), Tables 4-6, the scaled and unscaled emission choices and collision-rate uncertainty. |
| [Dalgarno et al. 1999](../../references/Dalgarno_1999_ApJS_125_237.pdf), ApJS 125, 237 | Gas-state assumptions, electron energy degradation, and Section 5.3's distinction between vibrational radiation and collisional conversion into heat. |
| [Glassgold et al. 2012](../../references/Glassgold_2012_ApJ_756_157.pdf), ApJ 756, 157 | Construction of the energy inventory from electron degradation and subsequent chemical processes. |
| [Koskinen et al. 2022](../../references/Koskinen_2022_ApJ_929_52.pdf), ApJ 929, 52 | Upper-atmosphere method and Appendix B; species transport, common momentum and energy equations, diffusion constraints, and limitations of molecular photodissociation. |
| [Schulik and Booth 2023](../../references/Schulik_2023MNRAS_523_286.pdf), MNRAS 523, 286 | Hydrodynamic discretization, Section 3.4's photoionization update, and Section 4.2's hydrostatic tests. |

NASA ADS metadata confirmed the seven published identities: `1993JChPh..99..885C`, `1998ApJ...496.1044Y`, `2013JPCA..117.9770M`, `1999ApJS..125..237D`, `2012ApJ...756..157G`, `2022ApJ...929...52K`, and `2023MNRAS.523..286S`. Content checks used the local journal PDFs. No preprint was substituted for these published versions.

The existing [audit driver](audit_20260905/run_checks.sh) was executed again in a temporary build directory with GNU Fortran 16.2.0 and `-O0 -g -fcheck=all -fbacktrace -fopenmp`. It compiled and ran the actual production routines, the E1 diagnostic, and the oxygen kinetics diagnostic. Additional small programs called the ROE, H3+, and H2 channel routines directly. Simple arithmetic checks evaluated the proposed energy equation; these are counterexamples to a proposed formula, not executions of a future implementation.

Temporary outputs are in `/tmp/exhale-physics-checks-O3Wzma/`. That location is evidence from this session, not a permanent project dependency. The important inputs and results are recorded below so that the conclusions do not depend on retaining it.

The successful exit of a diagnostic program does not certify every printed quantity. In particular, the existing E1 checks establish certain algebraic invariants but do not check the experimental detector interpretation.

## 3. Disposition of the first review

| First-review item | Assessment of revision 1 |
|---|---|
| R1: conserving species-energy update | The coupling and double-counting concerns are acknowledged, but the new explicit photon-energy residual is inconsistent. See Q2. |
| R2: steady residual call path | The major factual correction is adopted. Its description of the temperature and seed still needs the qualification in Q9. |
| R3: H3+ density continuation | The invalid exponent prescription is withdrawn. Full cooling normalization, the high-density join, and independent fit data remain to be specified. See Q7. |
| R4: heating residual and derivative | Substantially addressed: the residual, derivative, returned state, and failure test now describe the same intended solve. Its energy terms must first be corrected under Q2. |
| R5: transport and CO dependencies | Moving the spatial operator before the steady solver is correct. Merely replacing a CO residual by distance to a ceiling does not finish the constrained physics. See Q4. |
| R6: transported-H+ coverage | The new test inventory is appropriate. Initialization must verify that the desired seed survives to the actual carrier call. See Q10. |
| R7: finite-difference JFNK action | Corrected. Perturbation-size and inner-closure checks are appropriate. |
| R8: thermodynamic normalization and spin | Substantially addressed. The internal-state assumptions must also agree with the radiative and electron-degradation models. See Q3 and Q12. |
| R9: energy closure versus energy distribution | Improved, but the channel data themselves were misinterpreted, and neutral excited fragments remain important. See Q1 and Q3. |
| R10: regression provenance | Partly addressed. The refresh order and the use of “round-off” still need correction. See Q11. |
| R11: missing physical closures | Most now have owners. Electrical-current closure is still absent from the actual D1 specification, and some required local tests are delayed too far. See Q5 and Q12. |
| R12: dependencies and effort | Removing calendar estimates and defining D0 is reasonable. Calling all of A3 independent contradicts its new constraint-model work. |

These are assessments of the proposed work. Acceptance of an explanation in a plan does not establish that a source defect has been repaired.

## 4. Findings and required revisions

### Q1. The H2 branching construction misreads Chung's detector observable

Priority: critical for the molecular photoionization model. Status: confirmed source error, supported by the published experimental method and direct execution. Relevant plan tasks: A5, E1, D0 reaction and photon ledgers.

In [h2_photo_channels.f90](../src/modules/functions/h2_photo_channels.f90), the module header and `h2_channel_cross_sections` interpret Chung's ratio as

\[
r=\frac{N_S+2N_D}{N_M},
\]

where M is molecular single ionization, S is dissociative single ionization, and D is double ionization. The source then says the paper's cross-section normalization overcounts the double events.

The published method states that the two ions from a double event arrive almost simultaneously and cannot be resolved. The text on p. 886 includes the exact phrase “detector responds with a single pulse”. The signal ratio therefore corresponds to

\[
r_{\rm detector}=\frac{N_S+N_D}{N_M}.
\]

The normalization in the paper must be interpreted with that detector response. Its equality between the sum of the two measured ion-channel cross sections and the ionization-event cross section is not contradicted by the existence of two physical protons in a D event. [Chung et al., published version, pp. 885-888](../../references/Chung_1993JCP_99_885.pdf).

With \(q_D=N_D/(N_S+N_D)\) and an ionizing absorption cross section \(\sigma_i\), the corresponding reconstruction is

\[
\sigma_M=\frac{\sigma_i}{1+r},\qquad
\sigma_S=(1-q_D)\frac{r\sigma_i}{1+r},\qquad
\sigma_D=q_D\frac{r\sigma_i}{1+r}.
\]

The physical proton and electron yields are subsequently computed from stoichiometry: S contributes one, and D contributes two. Detector response belongs in the data inversion; particle multiplicity belongs in the reaction source.

The paper's statement that about 20% of its H+ cross section at 80 eV comes from double ionization must be interpreted consistently with that event observable. It cannot be assigned directly to the fraction `2*N_D/(N_S+2*N_D)` used by `frac_H2_double_of_protons`.

A direct test called `h2_channel_cross_sections(80.d0,1.d0,sig)` with `h2_double_ionization_model='chung80'`. The comparison below retains the code's existing interpolation/rounding of the measured ratio; it is not a new fit to Table II.

| Quantity | Actual routine | Reconstruction using the detector interpretation |
|---|---:|---:|
| Molecular channel M / total absorption | 0.797160208148 | 0.779590000000 |
| Single dissociative channel S / total absorption | 0.180302037202 | 0.176328000000 |
| Double channel D / total absorption | 0.022537754650 | 0.044082000000 |
| `(S+D)/M` | 0.254452981696 | 0.282725535217 |
| `D/(S+D)` | 0.111111111111 | 0.200000000000 |

The input ratio recovered from `frac_H2_dissociative_ionization` is 0.282725535217. The current routine reproduces it as `(S+2*D)/M`, which is the wrong observable. Its four channel fractions still sum to one. Thus the measured E1 channel-sum error of at most `4.441e-16` and zero stoichiometric violations do not detect this error.

Required revision: reconstruct the M/S/D probabilities from the experimental detector model, redefine the double fraction after deciding what it represents, and rederive the `yan_rho` conversion consistently. A published event ratio \(\rho_D=N_D/(N_M+N_S)\) implies \(\sigma_D/\sigma_i=\rho_D/(1+\rho_D)\); combine this with the detector ratio and check nonnegative S. Recheck the underlying ratio interpolation against Table II. Add tests that reconstruct both the measured detector ratio and the physical particle yields.

The first review did not identify this problem. Its E1 invariant results remain valid as algebraic measurements, but they must not be used as evidence that the four probabilities match the experiment.

### Q2. The proposed energy residual subtracts the photoionization energy twice

Priority: critical. Status: inconsistency in the proposed equation. Relevant plan sections: D0 item 3 and Section 4.3.

Revision 1 proposes

\[
F_E=\Delta u_{\rm th}+\sum_s\epsilon_s^0\Delta n_s
-\Delta t\,[H_{\rm rad}-C]=0,
\]

but defines `H_rad` as photoelectron, fragment, and photolysis heating. Such heating is already the part of the photon energy left after chemical storage and other destinations have been removed.

Consider a single H photoionization, with no secondary ionization, escaping radiation, or transport. Take a photon energy of 20 eV and a chemical ionization energy of 13.598434599 eV. The physical increments are

\[
\Delta E_{\rm chem}=13.598434599\ {\rm eV},\qquad
\Delta u_{\rm th}=6.401565401\ {\rm eV}.
\]

Putting the physical heat into the proposed right-hand side instead gives

\[
\Delta u_{\rm th}=6.401565401-13.598434599
=-7.196869198\ {\rm eV}.
\]

The chemical threshold has been charged twice. At the photoionization threshold the same formula would cool the gas by the ionization energy even though the incoming photon supplies that energy completely.

For the represented H2 double-ionization products, the newly executed audit gives a chemical increment of 31.67494411407 eV. At 80 eV, assuming all nonchemical energy thermalizes locally, the correct thermal increment is 48.32505588593 eV. Substituting that heating into revision 1's equation gives only 16.65011177186 eV. These values are evaluations of the written equation, not atmospheric predictions.

There are two consistent formulations:

1. Use the chemical increment of every reaction and put the full net energy transferred from radiation and external agents on the right-hand side. For a photoevent this includes the energy going into chemical storage. For radiative recombination, escaping photon energy includes the release of stored chemical energy as well as the captured electron's energy.
2. Keep already partitioned thermal photoheating and thermal radiative cooling, and include chemical terms only for processes whose chemical energy is not already accounted for by those thermal source definitions. This needs a reaction-by-reaction derivation, especially for finite composition jumps and eliminated species.

The first formulation gives the clearer D0 accounting identity:

\[
F_E=\Delta u_{\rm th}+\Delta E_{\rm chem}
-\int_{t_n}^{t_{n+1}}Q_{\rm external,net}\,dt=0.
\]

Here `Q_external,net` is energy entering the represented material system minus energy leaving it. If fast electrons or excited states have independent stored energies, include their increments in the left-hand side or describe their consistent elimination. The variable should not be called thermal heating if it includes chemical storage.

This also changes the recombination migration. In [util_ion_eq.f90](../src/modules/radiation/util_ion_eq.f90), `eval_cool` calls the existing recombination-cooling coefficients and separately adds ionization-potential cooling for collisional ionization. Removing the latter does not automatically convert the former into full escaping photon energy. Derive each mapping; retaining a thermal recombination fit while adding the entire chemical release can create artificial heat.

The general need to carry composition derivatives in material internal energy is also explicit in the [Cantera reactor energy derivation](https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html#energy-equation). This review's photon counterexamples follow directly from the event ledger, rather than from adopting another code's implementation.

Required gates: isolated photoionization at threshold and above threshold; radiative recombination with escaping energy; collisional ionization and its reverse; a photon-driven chemical cycle with unchanged final composition; and the closed reacting cell. A zero net species increment does not mean zero reaction throughput or zero radiation exchange.

### Q3. The photon ledger must include excited fragments and the energy-degradation assumptions

Priority: high. Status: source omission and remaining model specification. Relevant tasks: A5, A6, B1, B4.

The revised distinction between the H2 double-ionization vertical threshold and the chemical energy of separated products is appropriate. However, 51.4 eV is a vertical threshold in Yan et al.'s Section 4, while Chung et al. discuss an observed double-ionization appearance near 47-48 eV. These represent different physical quantities. The plan should label a hard zero below 51.4 eV as part of the adopted vertical approximation and investigate the threshold region before treating it as a measured onset. A repair to the energy sum alone does not validate the cross-section threshold model. [Yan et al. 1998, Section 4](../../references/Yan_1998_ApJ_496_1044.pdf); [Chung et al. 1993, p. 886](../../references/Chung_1993JCP_99_885.pdf).

The same paper reports excited neutral dissociation products and discusses a contribution producing H(2p)+H(2s) on p. 887. By contrast, the current neutral channel stores only two ground-state H atoms, and `util_ion_eq.f90` assigns the entire photon energy minus D0(H2) to heat. That cannot represent an excited-fragment branch without accounting for its excitation and subsequent radiation or collisional quenching. The total neutral cross section alone does not specify all final-state branching fractions. The dissociative single-ionization channel likewise permits an excited H fragment in the experimental description.

Required revision: extend the event metadata to distinguish translational energy, internal excitation, and escaping or reabsorbed photons. Where detailed state branching is unavailable, record a bounded approximation and its supporting evidence. Do not infer complete local thermalization from the existence of a conserved ground-state stoichiometric row.

Dalgarno et al.'s degradation model assumes a particular H2 internal population and local deposition; Section 5.3 distinguishes vibrational radiation from heat produced by collisional quenching. Glassgold et al. additionally account for heat from later chemical reactions. Consequently, electron degradation, subsequent reaction heat, and the chosen thermal EOS must be reconciled without assigning one excitation energy twice. A scalar fraction of electron energy is not automatically the complete eventual heat of the irradiated gas. [Dalgarno et al. 1999](../../references/Dalgarno_1999_ApJS_125_237.pdf); [Glassgold et al. 2012](../../references/Glassgold_2012_ApJ_756_157.pdf).

The plan's stopping-length criterion for energetic fragments is useful. Extend it to the primary and secondary electrons where local deposition is assumed. A calculation at the mean electron energy remains an approximation to an integral over the electron energy distribution; the proposed label should identify that distribution and its domain of validity.

### Q4. A distance to the CO ceiling does not define the constrained reaction-transport problem

Priority: high. Status: incomplete proposed closure. Relevant tasks: A3, D0, A4+B3, D1.

The actual `carrier_source` sets `src(ic_CO)=0`. After solving the carrier system, `limit_to_element_budget` changes CO to its equilibrium ceiling and may also rescale OH, H2O, H2, and H+. Revision 1 proposes to report distance to the CO ceiling wherever the ceiling binds and the transport residual elsewhere.

This can give zero for a clipped value without establishing a physically admissible reaction source. For a one-sided destruction model, one possible mathematical statement is

\[
\partial_t n_{\rm CO}+\nabla\cdot F_{\rm CO}=-\lambda,
\quad \lambda\geq0,\quad n_{\rm CO}\leq n_{\rm eq},
\quad \lambda(n_{\rm eq}-n_{\rm CO})=0.
\]

All terms have to satisfy this system. At the active bound, evaluate the implied sink and verify its sign. A boundary value at the ceiling requiring a net creation source cannot be accepted as instantaneous destruction. The active set and its derivatives must follow the same equations as the returned residual.

CO destruction must also release the appropriate carbon and oxygen inventory and consume the corresponding dissociation energy unless an external channel supplies it. The ceiling depends on temperature and available elements, so the constraint, other species, and energy are coupled. A numerical report describing CO as a reservoir does not supply these equations or a physical relaxation-time argument.

The first review's statement about step subdivision should be understood conditionally: an already clipped state generally fails the original zero-source equation. A new solution at a smaller step might avoid activating the bound, so step halving is not universally incapable of finding an admissible step. It nevertheless does not justify an unexplained sink whenever the ceiling is active.

Required revision: either derive and validate the restricted instantaneous-destruction model, including its energy and applicability conditions, or exclude configurations that require this ceiling until a justified closure exists. Move that derivation out of the independent Phase 1 list. The reference-scale correction and returned-state diagnostic can remain bounded tasks, but their acceptance must not certify this incomplete CO physics.

### Q5. The thermal-energy flux convention and electrical-current closure remain unspecified

Priority: high. Status: incomplete governing system. Relevant tasks: D0 items 3 and 5, D1, Section 6's closure table.

The plan requires `sum(m_s*J_s)=0` and an enthalpy flux `sum(h_s*J_s)` without stating whether `h_s` includes formation energy. This matters because the plan evolves energy excluding formation energy.

For number fluxes J_s in the barycentric frame, define

\[
E_{\rm chem}=\sum_s\epsilon_s^0 n_s,
\qquad F_{\rm chem}=E_{\rm chem}v+\sum_s\epsilon_s^0J_s.
\]

Species continuity gives

\[
\partial_t E_{\rm chem}+\nabla\cdot F_{\rm chem}
=\sum_s\epsilon_s^0S_s.
\]

Subtracting this equation from the formation-energy-inclusive material equation yields the thermal convention. Its diffusive enthalpy contribution uses \(h_s^{\rm th}=h_s^{\rm full}-\epsilon_s^0\), and its chemical source uses \(-\sum_s\epsilon_s^0S_s\). Adding full enthalpy diffusion to a thermal-energy equation without this subtraction counts the transport of chemical energy inconsistently. The finite-volume identity must also use the same chemical face flux as the species rows.

In the existing source, `U(3)` stores kinetic plus thermal energy, as shown by [UW_conversions.f90](../src/modules/functions/UW_conversions.f90). D0 should distinguish this variable from thermal internal energy alone. Keeping the restart column does not eliminate the need to state the energy convention in restart metadata and verify old-state loading.

Mass balance alone also does not close charged diffusion. A quasineutral mixture needs electron transport and a specified current constraint or electric-field equation. A common isolated-column approximation is zero net current, \(\sum_s z_s e J_s=0\), with the corresponding electric field; it is not implied by \(\sum_s m_sJ_s=0\). Thermal and momentum exchange must be consistent with that choice.

The current `carrier_diffusivities` explicitly sets proton molecular diffusion to zero and continues bulk advection; its comments identify the missing ambipolar treatment. This is not an implementation of charged diffusion. Koskinen et al.'s Appendix B explicitly includes the electric term in its diffusion equation and the mass constraint separately in Eqs. (B9)-(B10). [Koskinen et al. 2022, Appendix B](../../references/Koskinen_2022_ApJ_929_52.pdf).

Required revision: specify either a justified common-velocity approximation with negligible charged drift, or an explicit charged transport closure. Add separate tests for net diffusive mass, electrical current, and energy flux. A neutral diffusion budget does not certify the charged model.

### Q6. A2 misses a ROE branch that returns NaN for a valid rarefaction

Priority: high. Status: newly reproduced source defect. Relevant task: A2.

The proposed exponent and upstream-density corrections are both necessary. However, [speed_estimate_ROE.f90](../src/modules/flux/speed_estimate_ROE.f90) only enters its rarefaction/shock alternatives when `Q = p_max/p_min` exceeds 2. An equal-pressure expansion stays in the initial PVRS estimate regardless of whether that estimate is admissible.

Directly executed input:

```text
gamma = 5/3
WL = (rho, v, p) = (1, -1, 1)
WR = (rho, v, p) = (1,  1, 1)
```

This is a nonvacuum symmetric rarefaction. The exact constant-gamma result is

\[
c=\sqrt{5/3},\quad
p_* = \left[1-\frac{(\gamma-1)(v_R-v_L)}{4c}\right]^{2\gamma/(\gamma-1)},
\quad c_*=c-\frac{\gamma-1}{2}.
\]

The independently evaluated pressure is `0.224614296367211`; the star sound speed is `0.957661115402472`. The production routine returns `u_star=0` and `cL_star=cR_star=NaN`. Here its initial pressure estimate is negative while the estimated star density remains positive. A2's two edited formulas are inside a branch that this input never enters.

The old audit's much stronger separating state returns finite sound speeds, but finiteness alone is insufficient: negative pressure and negative density can produce a positive ratio inside a square root. That diagnostic must not become a passing vacuum test merely because both returned values are finite.

Required revision: include branch selection, positive star-state admissibility, and an explicit vacuum limit in A2. Compare symmetric expansions, asymmetric shocks and rarefactions, stationary contacts, and density/pressure rescalings with suitable reference solutions. The routine uses approximate star estimates, so comparison tolerances should reflect that approximation; exact agreement is appropriate only in cases where the approximation reduces to the exact solution. Retain the proposed restriction to the verified constant-gamma regime.

### Q7. H3+ needs a complete emission choice and tests of both interpolation variables

Priority: medium to high. Status: remaining design gaps and directly verified source behavior. Relevant task: B2.

The revised low-density exponent and continuity requirements are good. The executed production probe reproduces the previously reported nonzero cooling at zero H2: at T=1000 K and n(H3+)=1 cm^-3, the routine returns `5.455701880206e-15 erg s^-1 cm^-3` for n(H2)=0, 1, and 1e6 cm^-3. This is an error for the isolated collisional component in the stated H2-collider model.

Smoothing the LTE emission alone leaves another discontinuity. At T=5000 K, the actual non-LTE routine returns

```text
s(1e14 * (1 - 1e-8) cm^-3) = 0.998499999743115
s(1e14 cm^-3)               = 1.000000000000000
```

The table's last entry is 0.9985, but the code switches immediately to one at the last density. B2 must cover the full function of temperature and collider density, including high-density continuation and derivative behavior at table boundaries.

The published tables have distinct roles. Table 4 contains emission values, Table 5 contains polynomial coefficients, and Table 6 contains non-LTE emission and factors. Table 4 distinguishes the unscaled and scaled choices; at 5000 K its quoted values are `3.4119e-18` and `4.7688e-18 W molecule^-1 sr^-1`. Table 6's high-density emission approaches the scaled value. The current LTE routine returns `3.42707344482506e-18` at 5000 K, following the unscaled fit. A proposal to fit “the published data” must name the selected emission definition and use compatible non-LTE normalization. [Miller et al. 2013, Eqs. (6)-(10), Tables 4-6](../../references/Miller_2013_JPCA_117_9770.pdf).

Comparing a new smooth function only with Table 5's discontinuous polynomial outputs is not independent physical validation at their joins. Use the emission data, a justified reconstruction, or the underlying state/transition calculation where available, and report the uncertainty of interpolation between sparse samples. The published non-LTE calculation itself has an uncertain collision-rate normalization.

The linear low-density gate applies to collisionally powered cooling under the selected population assumptions. Radiation pumping, formation pumping, and other colliders require their own accounting; they must not be forced to disappear merely because H2 density vanishes.

### Q8. The time discretization and steady operator must be specified together

Priority: high for numerical claims. Status: missing specification. Relevant tasks: A4+B3, C1, source-step and coupled convergence gates.

The written species and energy residuals are backward-Euler equations. Their time accuracy is first order unless the algorithm is changed. A third-order hydrodynamic Runge-Kutta step does not make a split backward-Euler chemistry update third order. State the intended physical time integrator, the ordering of transport, chemistry, radiation, and energy updates, and the expected global order.

The phrase “built from the Phase 3 source step” is also ambiguous for C1. A physical steady equation contains instantaneous source rates together with spatial flux divergence. A finite-step map, for example `(Phi_dt(x)-x)/dt`, can depend on its splitting and step size when combined with transport. A fixed point of separated updates need not exactly satisfy the unsplit stationary equations at finite dt.

Required revision: define one continuous source function and one spatial operator. Derive the time discretization and the stationary residual from them. If a finite-step map is used to accelerate convergence, test the physical unsplit residual independently. Separate a common physical time step from local pseudo-time stepping used only to seek a steady state.

Step rejection must restore the entire attempted state: species, energy, cached radiation/columns, and any mutable background used in rate evaluation. Otherwise the retried step solves a different problem. Source substeps need their own accumulated radiation exchange and reaction progress for the energy budget.

An equilibrium elimination that is solved together with the energy constraint can conserve material energy; it need not wait for finite-rate integration of every species. Therefore revision 1's statement that no increment can claim any closed-cell gate before increment (iii) is too categorical. An earlier increment may pass the gate for its explicitly defined equilibrium model. Time-resolved chemical trajectory convergence remains a different claim.

### Q9. The revised steady-path description is mostly correct but still compresses different states

Priority: medium. Status: source-description precision and acceptance-gate consistency.

The main correction in Section 3.1 is verified: `eval_residual` retains interior U after `ioniz_eq` and refreshes the ghosts before assembly. However, the phrase describing the flux temperature as p divided by the new particle count needs to specify which p is meant.

`assemble_residual` calls reconstruction, which derives primitive states from U using the current caloric composition. Its viscous/conduction branch also calls `U_to_W` again before computing `Tc=W(3,:)/n_part`. With a molecular caloric EOS, this recomputed pressure can differ from the pressure derived using the seed. For a constant-gamma gas, pressure at fixed U is composition independent. These cases must not be conflated. See [steady_newton.f90](../src/modules/time_step/steady_newton.f90), [steady_residual.f90](../src/modules/time_step/steady_residual.f90), and [UW_conversions.f90](../src/modules/functions/UW_conversions.f90).

Likewise, saying all rates use the seed composition is too broad. The radiation calculation is lagged in specific ways, but the local nonlinear chemistry evaluates trial populations, and the final heating/cooling assembly uses the returned composition at the sweep temperature. Name the seed-dependent optical depths and coefficients separately from the nonlinear local state.

C1's “atomic cases with the keys off byte-identical” gate also needs a reference point. It can be a useful isolation check against the immediately preceding, physically corrected implementation for an untouched option-off path. It cannot mean identity to the pre-Phase-3 results after shared energy and spatial equations have changed. Sections 3.1 and 9 correctly admit those results may move.

### Q10. The new tests need construction and acceptance details

Priority: medium to high. Status: incomplete validation design.

The three proposed H+ configurations are relevant, but a filename and intended abundance do not establish that the actual solver receives that abundance. `input_read` requires molecular chemistry and carrier transport for transported H+, and rejects the present Newton route. In `ionization_equilibrium.f90`, the cold-start sweep can initialize carriers from equilibrium before `bg_ready` is true; the transported H+ constraint is applied under the background/restart condition. A naive zero-seed input can therefore become a nonzero equilibrium state before the code path under test is reached.

Construct the test through a controlled driver or a verified restart. Assert H+ and the other relevant populations immediately before the actual source/Jacobian assembly. Check the implied electron density, element budgets, and lower-boundary state as well. Keep the tiny-abundance derivative test separate from a full transport-front test.

The proposed zero-irradiation hydrostatic column is not automatically stationary if radiative cooling, chemical relaxation, or heat conduction remains active. For the mechanical test, use a fixed composition and a thermodynamically stationary state with compatible boundaries, or a stated balancing source. Then test radiative-chemical equilibrium separately. Zero incident radiation alone is not a thermal equilibrium condition.

D3 also needs a count for the full differential system. Incoming hyperbolic characteristics already include composition modes in a multispecies Euler system; they should not be counted twice by adding an independent condition for every advected species. Diffusion, viscosity, and conduction introduce additional boundary requirements. Count those after naming independent species and constraints, and cover flow reversal and the outer boundary.

Finally, changing local time steps cannot repair a spatial pressure-gravity imbalance at a stationary state. D2 should first evaluate the discrete hydrostatic residual. If a residual remains, correct the spatial balance before diagnosing a time-integration instability. The published AIOLOS hydrostatic tests provide a useful example of this distinction, without requiring EXHALE to copy that solver. [Schulik and Booth 2023, Sections 3.1 and 4.2](../../references/Schulik_2023MNRAS_523_286.pdf).

### Q11. Baseline handling and invocation diagnostics still contain factual problems

Priority: medium for provenance; low for the counter duplication. Status: source-verified tooling behavior and terminology.

Section 3.2 describes refreshing goldens and moving old goldens afterward. The actual `run_check.sh golden` branch copies files directly into `golden/<case>/`, overwriting the previous reference. Archive and verify the old reference first. If provenance is unavailable, mark it unknown rather than inventing a toolchain, stopping reason, or a clean source identity. Copying an output does not recover its provenance.

Do not promote current outputs solely because all twelve differ from older references. Validate their identity and purpose first. The current review did not repeat the earlier twelve-case numeric comparison because none of those products was changed and the present issues are local.

The expression “round-off (below 0.1%, from flags or constant unification)” is incorrect. A 0.1% threshold is a declared relative tolerance. Round-off depends on precision, conditioning, accumulation, and the numerical path. Changing a physical constant is a model/input change, even if its effect happens to be below that tolerance. The script currently defaults `REGRESSION_REL_TOL` to `1e-3`; a successful comparison at that setting must not be described as round-off equality without further evidence.

A1 asks to add a positivity-repair counter, but `RK_rhs.f90` already has `n_faces_flux_positivity_limited`, and `EXHALE_main.f90` reports it when positive. Decide whether the missing information is attempted calls, accepted repaired interfaces, a zero-valued summary, or direct LLF use. Reuse the existing quantity where appropriate. A counter can demonstrate that an affected path was reached; it cannot predict the magnitude of the resulting atmospheric change.

### Q12. Phase ownership must include the physical validity of the molecular state

Priority: medium to high. Status: remaining planning dependencies.

The common H2 partition function is a sensible thermodynamic foundation. It must remain consistent with the internal populations assumed by chemistry, electron degradation, and infrared transfer. A gas can retain substantial H2 above the density at which some internal levels cease to be thermalized. In that regime, assigning an LTE internal energy and an unrelated non-LTE emission correction requires a derivation or a clearly limited approximation. An ortho/para policy by itself does not settle vibrational equilibrium.

The plan also postpones heat capacities of other molecules and ions until Phase 6 while proposing a general H/O/He energy ledger in Phase 3. This is acceptable only if the Phase 3 accepted configurations explicitly satisfy the trace-species assumption used by the EOS. A generalized reaction ledger does not make an incomplete caloric model valid for water-rich gas. Add a runtime/gate condition based on the actual internal-energy and heat-capacity contribution or restrict the accepted composition domain.

Photon conservation across cells is assigned only to Phase 6, but the Phase 3 energy source already consumes those absorption rates. The basic monochromatic cell balance should accompany the source/radiation interface: absorbed photons equal entering minus leaving photons, and each absorbed event enters the material ledger once. Complex band overlaps and the coupled atmosphere can remain later work. Similarly, local energy conservation should be tested at the earliest increment that claims it, rather than being deferred solely by phase numbering.

Validity flags are useful evidence. They do not authorize production interpretations outside the stated domain. When a Knudsen or thermal-coupling diagnostic fails, the report must identify which conclusions remain supported. A finite-range assumption can be deliberate, but its justification must be physical. If the equations require ionization transport, a user preference can select a restricted experiment; it cannot establish that equilibrium elimination remains accurate there.

## 5. Recommended changes to the development sequence

Retain Phases 0-6 with the following adjustments:

| Stage | Revision needed before accepting its result |
|---|---|
| Evidence and test construction | Preserve references before refresh; add detector-observable H2 tests and the equal-pressure ROE rarefaction; verify H+ seeds at the call being tested. |
| Bounded source corrections | A1 and the A3 reference-scale change remain bounded. A2 includes admissible branch selection. CO constraint physics is removed from this category. |
| D0 derivation | Fix the full photon/material energy convention, excitation recipients, CO treatment, thermal versus full enthalpy, charged-current closure, independent species, boundary count, and physical time integration. |
| Thermodynamics and sources | Implement from the corrected ledger. Test isolated photoionization, recombination, excited products, and a closed reacting cell. Verify the elementary finite-volume photon balance here. |
| Spatial operator | Use the shared species face fluxes and the matching energy convention. Verify hydrostatic balance, mass and element budgets, diffusive energy, and the stated electrical constraint. |
| Steady solver | Solve the unsplit stationary equations; verify deterministic inner closure and distinguish the result from a finite-step split fixed point. |
| Coupled validation and extension | Retain term-by-term Koskinen reconstruction, grid/time convergence, lower/upper matching, CO kinetics, molecular infrared transfer, and observables. Compare models only after their different physics is explicitly identified. |

The local Koskinen paper is already sufficient to begin the equation and process comparison. Its Appendix B states the molecular photodissociation limitations and diffusion formulation. A missing original spectrum file can limit an exact numerical reproduction, but does not prevent those comparisons. Any request to the authors remains a separate user-authorized action.

## 6. Minimum acceptance tests added by this review

| Test | Required assertion |
|---|---|
| Chung detector inversion | Reconstructed `(S+D)/M` matches the adopted measured signal ratio; the double-event fraction uses the same definition; physical proton/electron multiplicities remain correct. |
| H photoionization | A 20 eV event stores 13.598434599 eV chemically and leaves 6.401565401 eV thermally under the stated simple local-deposition assumptions. |
| Recombination and photon-driven cycle | Stored chemical energy and escaping/recycled photons close the budget even when net species change is zero. |
| Excited molecular fragments | Excitation energy is assigned to a represented level, a controlled eliminated level, radiation, or collisional quenching exactly once. |
| CO active bound | The implied source has the permitted sign, product species close element budgets, and energy includes the destruction process. |
| Equal-pressure ROE expansion | The nonvacuum state in Q6 produces finite admissible estimates and passes the specified physical accuracy criterion. |
| H3+ cooling surface | Low-density collisional scaling, both temperature and density joins, high-density approach, and the selected published emission normalization are checked independently. |
| Diffusion-only column | Species fluxes close mass, elements, the specified current constraint, and thermal plus chemical energy with the chosen boundary fluxes. |
| Physical time refinement | The complete coupled update converges at its declared order under a common physical time step. Local pseudo-time iterations are not used as physical trajectories. |
| Hydrostatic column | The actual discrete pressure-gravity residual and spurious velocity decrease with the stated spatial scheme and compatible energy/boundary conditions. |
| H+ initialization | The asserted zero/tiny seed is present immediately before the production call, rather than only in an input or an earlier initialization state. |

These tests do not require a complete planetary atmosphere for every correction. They target the affected equations before a coupled model can hide the cause of a failure.

## 7. Reproduction notes and limitations

The diagnostics can be reproduced from the project root with:

```bash
bash docs/audit_20260905/run_checks.sh
bash docs/audit_20260905/run_review2_checks.sh
```

The new ROE case calls the actual routine with `WL=[1.d0,-1.d0,1.d0]`, `WR=[1.d0,1.d0,1.d0]`, and `gamma_ad`. It links `parameters.f90` and `speed_estimate_ROE.f90`. The H3+ check additionally links `h3p_cooling.f90` and evaluates the two density arguments printed in Q7. Its preserved source is [roe_equal_pressure.f90](audit_20260905/roe_equal_pressure.f90).

The H2 test links `parameters.f90`, `cross_sec.f90`, and `h2_photo_channels.f90`; sets the model to `chung80`; and calls the cross-section routine at 80 eV with unit total cross section. Its preserved source is [h2_detector_ratio.f90](audit_20260905/h2_detector_ratio.f90). The detector inversion equations and all comparison values are included in Q1. The arithmetic checks are preserved in [photoevent_energy_check.py](audit_20260905/photoevent_energy_check.py). The [diagnostic README](audit_20260905/README.md) describes the complete file inventory and execution scope.

The executed tests establish the local defects and algebraic balances described above. They do not measure their effect on mass loss, the base limit cycle, full-network Jacobians, or the existence/stability of a transported molecular steady wind. Those claims remain assigned to the revised gates. No existing atmospheric output or golden reference was regenerated during this review.
