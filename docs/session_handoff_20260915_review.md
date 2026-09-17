# Review of the September 15 Session Handoff

Date: 2026-09-15  
Document reviewed: [session_handoff_20260915.md](session_handoff_20260915.md)

## 1. Verdict and verification scope

The series contains useful corrections, particularly the distinction between conserved density and composition normalization, improved diagnostics, and recognition that a collocated velocity is not the numerical face flux. However, the conclusion that the molecular heat budget is physically settled is too strong. The radiative charge-transfer channel R23 is currently assigned its entire formation-energy difference as gas heating. The cited experimental papers explicitly identify a radiative mechanism. Energy bookkeeping alone cannot justify this thermal assignment.

The revised base boundary also has two unresolved structural problems: its branch depends on a distant diagnostic window, and switching between that window and cell 1 introduces a discontinuity despite the smooth branch-weight function. These matter especially for the low-Mach stationary difficulty under investigation.

| Finding | Assessment | Priority |
| --- | --- | --- |
| R23 radiative energy treated entirely as heat | Confirmed mismatch between implemented thermal term and stated radiative mechanism | High |
| Base entropy branch selected by distant wind data | Conditional stationary approximation applied without a physical-time restriction | High |
| Hard switch between base-branch discriminants | Source-established discontinuity for admissible opposing discriminants | High |
| Density authority inferred from discrepancy size | Incomplete restart-validation contract; substantial inconsistencies can change density instead of being refused | Medium |
| Unbounded chemical mass projection | Correct idea for roundoff, insufficient safeguards for larger closure failures | Medium |
| Pseudo-time/CFL interpretation | Stronger claims than the algorithm or cited theory establishes | Medium |
| Thread reproducibility and low-XUV conclusions | Need uncertainty and diagnostic qualification before scientific acceptance | Medium |

I inspected the working-tree implementation, not just the handoff. HEAD remains `43bc28cef58772560bae019559f3592335922212`, but numerous relevant files have uncommitted changes. A direct `md5sum EXHALE.x` returned **`db87b88d1ce53facf1d61084fa535ca5`**, matching the handoff's campaign binary. This does not prove that every currently edited source file is represented by that binary. The review concerns the source inspected during this session; a frozen source manifest is needed to associate each finding with an exact executable.

No build, atmospheric integration, regression refresh, campaign control, or source modification was performed. Atmospheric values from the handoff and item memos are **READ**, not independently remeasured here. Algebraic examples below are **DERIVED**, not atmospheric test results. No runtime race detector or complete heat-budget experiment was run. Only this review file was added.

## 2. R23: chemical energy conservation is not a thermalization model

### Implementation and actual consumer

`src/modules/lower_atmosphere/mol_rates.f90:950–966` identifies

```text
H2 + He+ -> H2+ + He
```

as radiative charge transfer, with the constant coefficient `7.2e-15 cm^3 s^-1`.

`src/modules/lower_atmosphere/molecular_reaction_heat.f90:454` nevertheless includes

```fortran
+ k23*nh2*nheii(j) *q(ir_R23)
```

in `gamma_chem`, where `q` is the difference of formation energies. The caller in `src/modules/radiation/util_ion_eq.f90:2121–2123` places that result in heating channel 15 and adds it to total gas heating. The inspected contraction contains no R23 photon-energy subtraction or radiative thermalization factor. The similarly named `R23` in `Cool_coeff.f90` is an unrelated transition-rate variable in a level-population calculation, not compensation for this molecular channel.

### Published evidence

[Schauer et al. (1989), JCP 91, 4593](../../references/Schauer_1989JCP_91_4593.pdf), introduction and experimental/analysis sections, distinguish dissociative and radiative charge transfer and explain the product-ion inference. [Boehringer and Arnold (1986), JCP 84, 1459](../../references/Bohringer_1986JCP_84_1459.pdf), experimental section and discussion on p. 1461, describe the H2+ channel as involving radiation; their discussion estimates an emitted wavelength near 153 nm and a vibrationally excited product. These are published journal PDFs, not preprints. The wavelength estimate is a model description, not a measured universal photon energy to insert unquestioned into the code.

Thus a formation-energy difference of about 9.16 eV does not imply 9.16 eV of immediate local gas heat. The energy must be divided among radiation, product internal excitation, and kinetic energy. Radiation can subsequently heat gas only through a specified absorption process. Neither high gas density nor telescoping reaction energies proves local absorption and complete thermalization of the emitted photon.

This is not necessarily double counting of a formation energy. It is an incorrect or unstated assignment of the energy recipient. The categorical L7g claim that its construction makes an energy-budget defect impossible should be withdrawn. Correct stoichiometric bookkeeping is necessary but not sufficient.

### Proposed correction and tests

Represent each relevant channel with separate chemical-energy release, prompt kinetic heating, internal excitation, and radiative emission. For R23, first establish a defensible product/radiation partition. Follow emitted photons with an explicit transfer or escape/reabsorption approximation whose domain is stated. Do not simply set the complete reaction heat to zero without accounting for the nonradiative remainder.

Add reaction-level tests that distinguish gas energy from gas-plus-radiation energy. Then measure R23's contribution at the same held state used for the L7g budget before re-solving an atmosphere. Only after that can the effect on the reported molecular base temperature be quantified. This review does not establish that R23 explains the full hot-base discrepancy, or even that it dominates it.

## 3. L21: useful stationary diagnosis, incomplete general boundary correction

### Distant cell averages are still not the numerical base flux

`wind_window_mass_flux` in `src/modules/states/base_boundary.f90:392–420` averages `W(1,j)*W(2,j)*r(j)^2` over `j_flux:N`. It does not evaluate the base Riemann flux, or even average numerical face fluxes. The LHS stationary profiles can make this a good approximation to the common wind flux, but the equality is conditional.

`base_boundary_states` calls this routine before constructing the characteristic base state. `characteristic_base_face_state` then uses the distant average to choose the entropy source. No physical-time versus stationary gate appears in that branch.

Consequences:

- Changing the convergence window can change a physical boundary condition.
- During a transient, material can leave the base while the distant flow reverses, or vice versa. A large distant flux is not proof of the local characteristic direction.
- The base residual acquires dependence on distant cells. A banded local preconditioner does not represent that dependence merely because a matrix-free residual evaluation sees it.

The successful L21 stationary comparison is valuable evidence of the old sign problem, but it does not establish that the replacement is a general characteristic boundary condition. Unchanged short control runs are also not proof of correctness for arbitrary transients.

### The selection guard introduces a discontinuity

The code around lines 507–512 is:

```fortran
M_branch = M_i
if (base_branch_on_wind_flux .and. have_F) then
   M_wind = F_wind/(rho_i*rb*rb)/c_i
   if (abs(M_wind) .gt. base_face_mach_blend) M_branch = M_wind
endif
w_rev = characteristic_branch_weight(M_branch/base_face_mach_blend)
```

Let `b` be the positive blend width and consider `M_i=-2*b`. For `M_wind` just below `+b`, the selected reversal weight is 1. For `M_wind` just above `+b`, it is 0. An arbitrarily small change in the distant flux can therefore switch between two different entropy states. This is a direct algebraic consequence of the implemented branch; the cubic weight being C1 does not make the composite selection C1.

The jump requires different candidate entropy states to produce a state discontinuity; that is precisely the kind of situation described by the L21 reservoir/interior temperature discrepancy. A reduced blend width does not remove the structural problem.

### Better direction

Determine inflow/outflow from a self-consistent local boundary-face solution, with the outgoing characteristic supplied by the interior and the appropriate incoming data supplied by the reservoir. Entropy should follow the locally resolved contact/mass transport. If this requires a small implicit boundary solve, include its derivative or a consistent numerical action in the stationary operator.

As a temporary stationary-only device, a continuation flux can be an explicitly named algorithmic parameter, but it should not silently become a nonlocal physical-time boundary law. Merely smoothing the distant-window switch may aid iteration but does not resolve the physical dependency.

Required focused tests are local flow reversal with a distant flow of the opposite sign, rest equilibrium, variation through both selection thresholds, and independence from a diagnostic-window change. Diagnose the L17 low-XUV residual again after separating the boundary-selection discontinuity from the HLLC contact-selection kink.

## 4. L18/L19: preserve conserved density, but do not infer intent from error magnitude

Preserving the loaded conserved density for a same-problem restart is the right direction. However, `load_IC.f90:884–912` decides which density is authoritative by whether the maximum discrepancy exceeds `restart_density_agreement_tol`, READ as `1e-8`. Below it, species are rescaled to the hydro density; above it, species determine the density.

This conflates two distinct situations: an explicitly authorized composition transformation and an inconsistent or damaged restart pair. A large discrepancy does not prove that a transformation was requested, and a small intentional transformation is not impossible. The handoff's statement that conserved density is now the authority is therefore only conditionally true.

Use explicit transformation provenance and restart intent:

- Same-equation restart: preserve conserved quantities and reject inconsistency beyond a stated roundoff allowance.
- Explicit initialization conversion: construct the target composition and thermodynamics under the conversion's chosen invariants, recording that it is a new seed.
- Small roundoff correction: report its magnitude and location, and verify the affected invariants.

The sweep projection at `ionization_equilibrium.f90:2870–2910` likewise forms `s_close=n_in_dim/m_of_comp` and rescales all species without a maximum correction check. The chemistry ledger is finalized before that projection. For roundoff-scale changes this is a reasonable normalization operation. For a larger failure it can hide a defect and leave acceptance diagnostics referring to the pre-projection state.

Add finite/positive checks, a maximum correction diagnostic, and a refusal or explicit recovery status beyond the rounding regime. Validate conserved element budgets, not only element ratios: uniform rescaling preserves ratios but changes all absolute nucleus densities. Reevaluate relevant returned-state reaction residuals when a nontrivial correction is made. Do not present the observed small corrections in the current campaign as a proof that larger corrections cannot occur.

## 5. Molecular rates and vibrational heating: what the papers do and do not establish

### R17 and removal of the old R20

The published Schauer experiment supports rejecting the earlier use of its low-temperature result as a large, independently measured HeH+ formation coefficient. Its product detection constrains combinations of pathways, and the paper explains that interpretation. This is a sound reason to correct the attribution and network assignment. It is not a proof that the removed pathway is identically zero at every possible collision energy.

The current R17 implementation subtracts the separately carried radiative channel from a measured low-temperature total and adds an Arrhenius term. The low-temperature total agrees in form with Boehringer and Arnold's fit. Their experiment does not uniquely determine the entire temperature-dependent branching prescription. The Johnsen drift-tube study and its effective collision-energy treatment should not be conflated with direct thermal measurements over the entire hot-wind range. The source marks important extrapolations; retain them as model uncertainty rather than calling the full interpolation experimentally established.

For the omitted three-body process, an applicability test should use the actual collider density and temperature, for example the ratio of three-body to two-body loss. A historical representative density is not a domain guard. If the ratio is outside the intended approximation range, include the channel with justified products or flag the model as out of scope.

### Black's H2+ + He coefficient is conditional, not a universal upper bound

[Black (1978), ApJ 222, 125](../../references/Black_1978ApJ_222_125.pdf), p. 126 and its rate table, obtains the approximate coefficient from collision cross sections averaged over a thermal H2+ vibrational population. The surrounding methods discussion explicitly assumes that population is thermalized and discusses possible abundance overestimation if competing losses deplete excited levels.

This supports a conditional thermal-population approximation. It does not establish a universal upper bound on the reaction coefficient for a gas with arbitrary pumping and product-state populations. The current `mol_rates.f90` commentary upgrades that conditional caution to a stronger bound. Replace it with the actual assumption and test formation, relaxation, reaction, and destruction times for H2+.

### R15's heat fraction needs a formation-state model

`molecular_reaction_heat.f90:430–448` applies `h2_vibrational_heat_fraction` to the entire R15 association energy. The relaxation routine uses an effective collisional-versus-radiative competition for excited H2, with H and H2 colliders. This is not by itself a distribution of three-body formation energy among translation, rotation, and vibration.

The local published copies of Hollenbach and McKee (1979) and Burton et al. (1990) motivate collisional quenching of excitation, but using a quenching fraction for a specified excitation does not establish that the whole chemical binding-energy release entered that excitation. Treat formation heating, fluorescent pumping, and thermal rovibrational cooling as distinct processes that share level populations and rates where appropriate.

A useful model separates prompt translational energy from product internal excitation and then follows relaxation or radiation of the latter. Include helium colliders if reliable data are available; otherwise state the missing channel and test its possible significance. Check the low-temperature validity again with the corrected 226 K reservoir rather than relying on comments describing an older 900–1350 K layer.

The small change produced by toggling the present IR modules does not establish that the opacity or heat partition is correct. It only measures sensitivity to those implemented modules under their assumptions.

## 6. Pseudo-time control and the low-XUV wall

The local published Kelley and Keyes (1998) paper, section 1.2, describes residual-based switched evolution relaxation and cautions that small steps do not by themselves imply small residuals. Mulder and van Leer (1985), Coffey et al. (2003), and Gropp et al.'s NASA report provide context for continuation and robust stationary solution methods. They do not make the exact current line-search-factor cut a general convergence theorem.

In `steady_newton.f90` around lines 17005–17015, the comment claims that a pseudo-time no larger than a multiple of the explicit CFL interval makes the shift dominate every Jacobian entry and makes the step explicit Euler. That conclusion needs additional conditions. For

```text
(M/dtau + J) dY = -F,
```

the explicit-like limit requires `dtau*M^{-1}*J` to be small in a relevant operator sense. A hydrodynamic CFL estimate alone does not establish this for a coupled source/radiation/chemistry operator or its conditioning. A finite multiple of a stability limit is also not automatically an asymptotically small parameter.

Keep the counter as an empirical control if useful, but name that status accurately. Measure the implicit step's difference from the explicit-like direction, actual linear residual, and row-block contributions at the suspected floor. Do not use proximity to a CFL-derived interval as proof that further nonlinear progress is impossible.

Likewise, L17 establishes failure of the attempted routes, not a proof that the state cannot be solved without one particular flux replacement. Candidate investigations include a consistent local boundary treatment, directional residual scans through contact switches, a semismooth or regularized nonlinear treatment, and a verified all-speed well-balanced flux. Do not remove physical pressure terms or change conservation merely to smooth a residual.

## 7. Reproducibility and record corrections

L15 should not be labeled unconditionally nonblocking. Thread-dependent differences can be harmless final-bit changes, or they can change acceptance, solver branch, and observables. Establish their range on frozen inputs and a frozen executable before interpreting a campaign as one reproducible solution set. One-thread controls alone do not establish the behavior of the eight-thread campaign.

For a representative subset, repeat eight-thread runs, compare with one thread, and compare physical arrays, final residuals, certification, mass loss, and synthetic observables. Separate scheduling/log differences from physics differences. Reevaluate stored results with a deterministic read-only certification path. L9 is consequently a valuable prerequisite for trustworthy handoff validation, not merely cosmetic output work.

The handoff itself needs several corrections:

- Section 2's case counts sum to **101**, not 95 (DERIVED directly from the seven listed rows). Later scheduling counts use a different breakdown. Regenerate one catalog-based status table and distinguish cases, closure iterations, exclusions, and attempts.
- The opening says log Mdot moves by -1.8 percent. The L21 table reports Mdot itself decreasing by about 1.77 percent; the listed logarithms change from 7.8756 to 7.8679, approximately -0.0077 dex. These are different statements.
- Sections 7–8 contain pre-rerun operational instructions, while section 11 says the campaign is running and forbids rebuilding. Mark the older steps as historical/superseded rather than leaving executable-looking instructions in conflict.
- Do not refresh golden files solely because changes have been attributed. First establish the corrected physical contract and independent invariants; then a deliberate reference update can be considered by the user. This review authorizes no such update.
- A matching executable hash is useful but should accompany a frozen source/flags/dependency manifest, especially while the working source continues to change.

## 8. Recommended order of work

1. Reopen the thermal-recipient part of L7g for R23. Audit a held state before expensive reruns and quantify the affected energy terms.
2. Reopen the generality and smoothness of L21. Derive a local face-consistent entropy condition and test the switch independently of HLLC's interior contact selection.
3. Harden L18/L19 with explicit transformation intent and bounded, observable projection.
4. Complete a bounded L15 reproducibility study and a read-only certification/output path before accepting campaign-derived inference.
5. Then assess state-resolved molecular relaxation, rate uncertainties, and all-speed solver improvements on representative cases rather than tuning the entire campaign at once.

The useful outcome of this series is not merely a larger number of certified runs. It is a state whose discrete equations, boundary conditions, energy recipients, and stored representation describe the same physical problem. Numerical certification should remain necessary, but it cannot replace checking those model assumptions.

## 9. Literature actually consulted

The following existing local PDFs were read with `pdftotext -layout`; the review used the relevant experimental, methods, or algorithm sections rather than relying on the handoff's citation claims. Some review was sectional, not a page-by-page audit of every paper or all numerical values in its figures.

- [Schauer et al. 1989, JCP 91, 4593](../../references/Schauer_1989JCP_91_4593.pdf).
- [Johnsen et al. 1980, JCP 72, 3085](../../references/Johnsen_1980JCP_72_3085.pdf).
- [Boehringer and Arnold 1986, JCP 84, 1459](../../references/Bohringer_1986JCP_84_1459.pdf).
- [Black 1978, ApJ 222, 125](../../references/Black_1978ApJ_222_125.pdf).
- [Hollenbach and McKee 1979, ApJS 41, 555](../../references/Hollenbach_1979_ApJS_41_555.pdf).
- [Kelley and Keyes 1998, SIAM J. Numer. Anal. 35, 508](../../references/Kelley_1998SIAMJNA_35_508.pdf).
- [Mulder and van Leer 1985, JCP 59, 232](../../references/Mulder_1985JCP_59_232.pdf).
- [Coffey et al. 2003, SIAM J. Sci. Comput. 25, 553](../../references/Coffey_2003SIAMJSC_25_553.pdf).
- [Gropp et al., NASA CR-1998-208435](../../references/Gropp_1998_NASA_CR_208435.pdf), a technical report rather than a journal article.

Burton et al. (1990) was identified locally and its use inspected in the source; its full coefficient comparison was not independently completed. No new literature search or ADS citation lookup was necessary to locate the existing PDFs. In particular, this review does not certify every quoted coefficient or published-fit comparison in the handoff.
