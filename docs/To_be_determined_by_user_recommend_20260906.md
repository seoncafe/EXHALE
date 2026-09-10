# Recommendations on the September 6 decisions

Date: September 6, 2026, KST.

## 1. Recommendation in brief

Do not approve the five recommendations as one package. Approve bounded failure
recovery, reject non-root chemistry as an accepted physical state, and distinguish
an audit of the existing implementation from approval of the equations that will
replace it. The proposed WASP-121 b composite is a reasonable provisional model,
but neither its normalization nor its helium photoionization accuracy is
established by the decision document.

| Decision | Recommendation | Condition for approval |
|---|---|---|
| 1. Carrier Newton failure | Choose (b), preferably an outer attempted-step controller for the intended coupled model. | Complete rollback, bounded retries, final-state residual checks, and explicit failure after recovery is exhausted. |
| 2. Ionization relaxation amnesty | Reject (a) as a physical acceptance policy. Choose the rejection branch of (b), with a separate initialization/continuation procedure. | An old composition is a restored state or an initial guess, not automatically a solution at the new density, temperature, and radiation field. |
| 3. D0 acceptance | Accept a corrected D0 as an implementation audit; do not yet accept it as the completed Phase 2 governing-system specification. | Resolve the conservation equations, independent variables, transport closure, energy accounting, and acceptance tests required by rev3 Section 4.2. |
| 4. WASP-121 b SED | Prefer candidate C over A or B as a provisional composite, subject to spectral validation. | Label the proxy and unresolved absolute normalization; do not call it an exact Huang reproduction or an unbiased helium input. |
| 5. Loaded spectra | Approve conversion of individually validated configurations; keep scientific reruns separate. | Match each benchmark to its actual reference, preserve existing outputs, and do not certify unresolved WASP-121 b inputs as reference reproductions. |

The most urgent additional implementation findings are:

1. The carrier failure check can be bypassed when the cell with the largest
   returned residual was constrained, without establishing that every other
   unconstrained row converged.
2. The steady solver counts acceptance class 5 as uncertified even though the
   ionization solver now defines it as a verified constrained-continuation root.
3. The final steady-state gate does not itself require certified chemistry, and
   the marching caller disables its carrier residual condition.
4. The energy update returns after two iterations without a final residual
   check and can impose a temperature floor without accounting for its energy.
5. The current SED source, documentation, and validation assumptions disagree
   about active thresholds and wavelength coverage.

These are correction or validation requirements, not choices between equally
valid physical models.

## 2. Scope and evidence

This review addresses
[the decision document](To_be_determined_by_user_20260906.md),
[development plan rev3](development_plan_20260905_rev3.md), and
[D0](d0_governing_system_20260906.md). It follows the relevant production callers
and consumers instead of treating those documents as proof of implementation.

The working tree contains ongoing changes by another developer. The inspected
HEAD was `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`; HEAD alone does not identify
the reviewed source. A snapshot of selected file hashes is recorded in Section 10.
Line numbers refer to the inspected working tree and may move during development.

Evidence labels used below:

- **Source-confirmed:** executable statements and their relevant callers were
  inspected. This does not establish the frequency or atmospheric impact of a
  branch unless it was also executed.
- **Measured:** computed directly during this review from the stored SED files.
- **Published:** checked in a local final journal PDF, not just a project summary.
- **Reported:** a result stated by another document, not reproduced here.
- **Proposed:** a design or test requirement, not an implemented capability.

The published methods were read first in the local Huang (2023), Salz (2019),
and Bourrier (2020) PDFs under `references/`. ADS was used to confirm the Huang
and Bourrier publication records. The official MAST MUSCLES description was also
checked to distinguish observations, reconstructions, and synthetic components.

No production source or scientific input was changed, no executable was rebuilt,
and no atmospheric simulation or full regression suite was run. The decision
document's 141 Newton iterations, residual of `5.2e4`, electron-density jump,
sixteen/seventeen regression cases, binary identifier, and bounds-check results
are **reported**, not new measurements from this review.

The only new numerical diagnostic is preserved in
[decision_sed_check_20260906.py](audit_20260905/decision_sed_check_20260906.py).
It reads existing data and does not regenerate spectra or simulation products.

## 3. Decision 1: failure recovery for carrier chemistry and transport

### 3.1 Choose rejection and retry, but do not promise that it always succeeds

**Physical judgment:** a failed nonlinear solve is not an admissible substitute
for the species balance equations. Increasing an iteration budget can be useful,
but it cannot replace a residual-based acceptance test.

Choose option (b). The statement that a rootless step necessarily becomes
solvable at sufficiently small `dt`, so that stopping disappears, is too strong.
For a smooth finite-rate system and an admissible initial state, a sufficiently
small backward-Euler step often has a nearby solution. That statement requires
regularity and admissibility assumptions; it does not prove that the particular
Newton algorithm will find it. Inconsistent boundary conditions, incompatible
constraints, discontinuous closures, invalid rates, and non-finite evaluations
can survive arbitrary step reduction.

An iteration-cap exit also does not establish that a root is absent. Keep
separate failure reasons for exhausted iterations, failed line search,
inadmissible state, non-finite rates, and an invalid physical configuration.

### 3.2 Prefer a controller around the complete attempted update

**Source-confirmed:** `EXHALE_main.f90:1055-1171` contains the existing hydro
retry loop. Element diffusion, carrier transport, excited hydrogen, ionization,
and the energy update occur after that loop. Reusing its name or reducing `dt`
inside a chemistry routine does not extend its rollback scope.

For the intended coupled species-energy development, prefer this structure:

1. Save the accepted state at the beginning of a physical step.
2. Execute hydro, transport, radiation/chemistry, and energy as trial operations.
3. Evaluate positivity, element and charge constraints, equation residuals, and
   the energy budget on the returned trial state.
4. Adopt the entire state only if every required condition passes.
5. Otherwise restore the saved state, reduce the physical interval, and retry.
6. Stop with diagnostics after the retry limit or minimum time interval is reached.

The checkpoint must cover more than `u`: species fractions, element abundances,
carrier backgrounds, excited populations, thermodynamic caches, photon columns,
reaction-rate caches, and any accumulated quantities that influence later
physics. Attempt counters may accumulate, but accepted-step conservation and
constraint statistics must not include rejected trials.

A carrier-only substep controller remains a valid intermediate option if it
holds an explicitly defined background fixed and the successful substeps cover
the original interval. It does not repair a bad ionization sweep or an
inconsistent energy update outside its scope. Do not present that option as
complete recovery for the coupled system.

### 3.3 Additional defect: a constrained cell can exempt a failed solve

**Source-confirmed:** in
`src/modules/lower_atmosphere/diffusive_photochemistry.f90:742-764`, a failed
`solve_status` stops only when
`carrier_cell_is_constrained(pct_worst_j)` is false. The cell identifier is the
location of the largest residual after the limiter. The solver status, however,
is established before the limiter (`solve_carriers`, approximately lines
2188-2223).

Consequently, a failed solve whose largest returned residual is in a constrained
cell can reach `carrier_write_back`. This condition does not establish that all
other unconstrained cells satisfied their equations. Nor does a constraint
acting somewhere in a cell establish the validity of every carrier row there.
The review did not execute a constructed atmosphere that triggers this bypass;
the insufficient acceptance condition is visible directly in the source.

**Required correction:** retain solver convergence separately from constraint
activity, and check every unconstrained row. A physically derived constrained
model needs its own residual and admissibility conditions. A CO ceiling that has
not been physically justified is a diagnostic modification, not an exception
that certifies a failed chemistry solve.

**Focused validation:** inject a failure after mutating a trial state; verify
complete restoration, no failed write-back, correct elapsed physical time, and
agreement with an explicitly subdivided integration. Include a case with a
large constrained residual and a smaller but still unacceptable unconstrained
residual. This last case targets the existing masking condition directly.

## 4. Decision 2: remove non-root acceptance from physical evolution

### 4.1 A cap near one is not a chemical-equilibrium tolerance

**Source-confirmed:** `ionization_equilibrium.f90:266` sets
`ieq_res_tol = 1.0d-6`. Class 4 is explicitly a non-root, and the molecular and
atomic acceptance branches can install that state into the returned species
profiles. `nonroot_streak_update` reports such events and stops only after a
persistent streak; the stated limit is 1000 sweeps.

The proposed residual cap of order one is about six orders of magnitude above
the current normalized root tolerance. A histogram of previously successful
relaxations does not make such an iterate satisfy the chemistry equations.
Reporting an event is valuable but does not validate its state.

Reject option (a) as a production acceptance policy. Choose the step-rejection
part of option (b). The alternative wording within (b), "always keeps the
previous state," is insufficient unless it means a full rejected-step restore.
After hydro has changed density and temperature, the old composition generally
is not equilibrium under the new state or radiation field.

### 4.2 Separate three different uses of an iterate

| Use | Permitted behavior | What must not be claimed |
|---|---|---|
| Physical time integration | Accept a finite-rate step that meets its discretized species-energy equations, or a properly constrained equilibrium update within its stated validity range. | A failed equilibrium iterate is not a finite-rate chemistry solution. |
| Initialization or pseudo-time continuation | Use intermediate non-root guesses while solving an explicitly identified numerical problem. | An intermediate guess is not a completed physical step or a converged atmosphere. |
| Jacobian or line-search evaluation | Evaluate trial states without modifying the adopted physical state; control inner-solve errors. | Finite residual entries alone do not prove that the trial lies on the intended equilibrium closure. |

For equilibrium chemistry, `dt` is not generally an argument of the local
algebraic equations. Reducing the outer step may improve the input state and
initial guess, but does not directly cure an intrinsically invalid equilibrium
closure. Use robust constrained solves or continuation for initialization. Where
reaction times are not short compared with transport and thermal times, use
finite-rate species evolution instead of enforcing equilibrium.

### 4.3 Additional defect: inconsistent interpretation of acceptance class 5

**Source-confirmed:** `ionization_equilibrium.f90:1624-1647` defines class 5 as
a root of the constrained element-conserving continuation solve. The reporting
routine and streak logic also treat it as a root. In contrast,
`steady_newton.f90:630` sets:

```fortran
n_uncertified_last = sweep%acc_n(4) + sweep%acc_n(5)
```

The steady solver then compares this count with the current iterate's count to
accept or reject trial states. Thus a valid class-5 solution can be classified
as uncertified and can contribute to rejecting a trial. Its comments also
describe class 5 as chemistry that the code cannot certify, contrary to the
current producer of that status.

This is an interface-contract inconsistency, not a user preference. Unify the
meaning of the acceptance classification, and check all consumers rather than
changing only one sum. If the constrained solver certifies a different equation
set in some configuration, represent that limitation explicitly and recheck the
requested equations; do not infer invalidity from the algorithm that found a root.

**Focused validation:** provide states certified through both the fraction solve
and constrained continuation; require identical root-validity decisions. Then
provide a genuine class-4 state and verify that it cannot pass final physical
acceptance. No atmospheric impact of this mismatch was measured here.

## 5. Decision 3: D0 is an audit, not yet the completed target specification

### 5.1 Conditional acceptance is appropriate

D0 states explicitly that it records what the code currently solves, proposes
no fix, and does not derive the charged transport closure, CO model, or complete
independent stationary species space. That is a useful audit deliverable.

However, rev3 Section 4.2 requires those derivations, thermal and chemical
finite-volume equations, exact invariants, consistent thermodynamics, a rollback
contract, and a verification matrix before D0 is complete. Accepting the present
audit as though it fulfilled that gate would change the meaning of the plan.

Keep a clear separation between:

- the corrected description of the existing implementation;
- the approved physical and numerical system to be implemented.

They may be separate documents or explicitly separated sections. Each known
inconsistency needs a target equation, an owner, and an acceptance test. Physical
derivation work can continue while these choices are settled; do not implement
an unresolved public state layout or closure by assumption.

### 5.2 Correct misleading statements before approval

**The stationary system does have an optional H2 row.** D0 Section 7.1 and the
decision document summarize the steady solver as having no species rows. The
hydrodynamic assembler alone has three rows, but
`steady_newton.f90:613-625` calls `carrier_steady_residual` and installs an H2
equation in the fourth component when `nvar_jac >= 4`. D0 Sections 7.2 and 7.3
actually describe that path, so its own summary contradicts its detailed account.

The correct statement is that the global solver has either three hydrodynamic
unknowns or those three plus H2, not a complete independent species-transport
system. H+ transport is rejected with the direct Newton configuration by
`input_read.f90:1919-1940`. This proves an unsupported direct-solver combination;
it does not prove that time marching can never approach a stationary solution.
Certification of that stationary limit is a separate issue.

**The mass-row alternatives are not physically equivalent.** D0 question 3
offers either removing chemistry's density overwrite or adding its production
term to the steady residual. In this nonrelativistic chemical model, reactions
redistribute conserved nuclei and mass; they do not justify a new bulk mass
source. Do not make a nonconservative update acceptable by copying it into the
stationary equations. Keep hydro density fixed in a local chemistry step and
verify the species mass sum against it. A density reconstruction that agrees to
round-off can be a diagnostic, not an independent physical update. The actual
mass discrepancy in a running atmosphere was not measured here.

**Small Lyman-alpha escape probability does not by itself invalidate full
collisional-excitation cooling.** Rev3's local energy-ownership rule is important:
collisional excitation removes thermal energy from the local gas even when the
eventual photon does not immediately leave the atmosphere. Radiative scattering,
collisional de-excitation, and true absorption must be accounted for consistently.
Do not simply multiply the existing excitation cooling by an escape probability
while retaining separately computed de-excitation heating. That can count the
same trapping effect twice. C7 identifies a coupling question, but its wording
must not imply that an escape-probability multiplier is the established fix.

**Already implemented changes are not open design questions.** The inspected
`sed_read` builds threshold-split bins; D0 C5 marks this closed, while Section 9
still asks whether to implement it. The decision document also describes all
planet inputs as power laws, whereas `WASP-52b/input.inp` already selects a loaded
spectrum. Correct the inventory at the level of individual configurations.

### 5.3 Recommended answers to the remaining governing-system choices

| Topic | Recommended decision |
|---|---|
| Energy variable | Retain thermal plus kinetic energy as the evolved hydrodynamic variable for the planned source update, but require an explicit formation/excitation-energy ledger and a conservative coupled temperature-composition solve. A formation-inclusive evolved variable is also valid if all fluxes, EOS inversions, and radiation terms are consistently reformulated. |
| Composition projection | Remove the temperature-preserving reset as an independent update. Solve temperature and composition under the chosen energy constraint. Keeping thermal energy as a variable does not justify keeping the reset. |
| Mass density | Chemistry must preserve the mass supplied to its local source step. Use reconstruction for consistency checking; do not add a fictitious chemical mass source to the steady equations. |
| Charged transport | Prefer a derived ambipolar closure with charge neutrality and an explicitly stated current condition. For a model without imposed current, zero net current is a defensible closure to derive, not a universal identity for every magnetized atmosphere. A common-velocity approximation is acceptable only within a quantified coupling regime. |
| CO | Until destruction rates, products, energy, and validity are justified, exclude any ceiling-active atmosphere from validated physical results. Keep its whole-atmosphere diagnostic label and cumulative activation record. **Overtaken 2026-09-06:** the rates, products, energy and domain argument were supplied (items B3b-CO, B3b-CO2) and the ceiling was then deleted (item CEILING-DEL), so there is no ceiling-active atmosphere to exclude; a run reports the destruction model's domain instead. |
| H3+ emission | Choose the model from the published fit definitions and their internal-population convention, not by which scaling preserves old cooling. Document collider dependence, temperature limits, and treatment of extrapolation; validate the low-collider limit. The disputed fit choice was not independently rederived in this decision review. |
| Oxygen energy accounting | Hold oxygen-enabled configurations outside the validated set until their reaction and excitation-energy ledger is complete. Do not block independently valid atomic configurations merely because oxygen work is unfinished. |
| Atomic H+ transport | Make atomic hydrogen ionization transport a physically independent capability rather than requiring molecules. Do not enable it for scientific use until charged closure, species equations, and its stationary acceptance test exist. |
| `_adv` outputs | Use the same physical source assembly as the main solution. Until a consistent molecular correction exists, identify uncorrected molecular columns explicitly or refuse the unsupported corrected product. Never silently mix incompatible states. |
| Duplicate cooling assemblies | Use one definition of the physical source terms with interfaces appropriate to the time-dependent and stationary equations. A test between duplicated formulas is useful during transition, but is not the preferred permanent ownership model. |
| Trace-metal diffusion | In the target conservative model, include the background response so the complete diffusive mass flux sums to zero. Derive it with the charged closure. Retain a trace approximation only with explicit omitted-order estimates and a stated domain. |
| Loaded-SED quadrature | Keep threshold splitting, but also define what the input samples mean and verify band energy, photon number, and reaction rates through the production consumers. Threshold placement alone is not a complete quadrature validation. |

For a fixed-volume local source step without transport or mechanical work, the
energy requirement can be written schematically as

```text
Delta u_thermal + Delta u_formation/excitation = integral Q_external dt.
```

Here `Q_external` is the net energy exchanged with radiation and other external
reservoirs under the declared local boundary. Energy already held in the EOS
must not also appear in the separate formation/excitation reservoir. In the full
flow, chemical and thermal transport contributions must sum to the material
energy flux using the same species face fluxes. Neither energy-variable choice
permits a photon to supply its full energy twice or an untracked reaction to
create thermal energy.

## 6. Decision 4: WASP-121 b requires a provisional model label

### 6.1 Candidate A matches relative band shape, not absolute band flux

The final Huang et al. (2023) paper, Section 2.2, states an adopted integrated
XUV flux of `1.6e6 erg cm^-2 s^-1`, a solar-based spectrum below 1700 A, an
LLmodels spectrum above it, and a Lyman-alpha flux of `1.0e5 erg cm^-2 s^-1`.
It also gives a Lyman-continuum flux of `2.69e5` and the three narrow-band values
listed below. These statements were checked in the
[local published PDF](../../references/Huang_2023_ApJ_951_123.pdf), with the
[journal record](https://doi.org/10.3847/1538-4357/accd5e) confirmed through ADS.

**Measured from the currently stored candidate A:**

| Band, A | Published Huang flux | Stored candidate A flux | A / published |
|---|---:|---:|---:|
| 700-800 | 1.4000e4 | 6.87981e4 | 4.914 |
| 800-912 | 3.6000e4 | 1.75716e5 | 4.881 |
| 912-1170 | 6.3000e4 | 3.02335e5 | 4.799 |

Flux units are `erg cm^-2 s^-1`. These measurements integrate a piecewise-linear
`F_lambda` with interpolated band endpoints; they are not measurements of the
production `F_E` quadrature or attenuated atmospheric rates.

Dividing the three ratios by the middle-band ratio gives `1.00679`, `1`, and
`0.98319`. Thus the statement about agreement within approximately 3 percent is
reasonable for the *relative band shape*. It cannot be used to imply agreement
in absolute irradiation. Candidate C inherits the same discrepancy below
1700 A if it uses A unchanged there.

The stored A header and its builder already describe this tension, but the
short decision document omits its importance. A match at the 1700 A join does
not resolve it. The adjacent-node flux ratio across A's join is **1.17944**, so
even the description "continuous" is approximate, not exact.

**Recommendation:** retain `1.6e6` as an explicitly chosen high-irradiation
scenario if that is the intended physical experiment. Do not substitute
`2.69e5` automatically, and do not declare the paper internally inconsistent for
every possible solar-based reconstruction. Clarify wavelength definitions,
spectral construction, geometric factors, and the model variant being compared.
The paper also discusses redistribution using a zenith angle and a flux factor;
these must be distinguished from the incident orbital SED. An author-supplied
spectrum or a verified reconstruction is needed for a strong reproduction claim.

### 6.2 Candidate C improves one assumption, not every uncertainty

The measured integrated fluxes over the decision document's 1700-2583 A band are:

- A: `1.72517e8 erg cm^-2 s^-1`.
- B: `9.56730e7 erg cm^-2 s^-1`.
- A/B: `1.80319`.

This supports the reported difference between the blackbody and the scaled
WASP-17 spectrum over that band. It does **not** establish an equal factor in
the He 2^3S photoionization rate or helium transit depth. The relevant rate is

```text
Gamma_23S = integral F_lambda(lambda) sigma_23S(lambda)
                    lambda/(h c) d(lambda),
```

with wavelength units handled consistently, the actual threshold, and the
appropriate attenuation in an atmosphere. Photons below 1700 A also contribute.
No such rate or atmospheric helium response was computed in this review.

The MUSCLES broadband products combine observations with reconstructed and
synthetic components; they are not direct measurements at every wavelength.
The [official MAST product description](https://archive.stsci.edu/hlsp/muscles)
identifies EUV scaling/DEM models, reconstructed Lyman-alpha, HST UV data, and
PHOENIX photospheric components. In addition, WASP-17 is a proxy for WASP-121,
not a measurement of the target star.

The current B builder normalizes its long-wavelength component using the ratio
of the target bolometric flux to the proxy's `BOLOFLUX`. That fixes an integrated
scale but does not transform the proxy's effective temperature, gravity,
metallicity, or UV line blanketing into those of WASP-121. The decision document's
claim that C has "no He 2^3S bias" should be removed.

### 6.3 Recommended spectral hierarchy

1. Prefer the original Huang input for reproducing that calculation, or a
   target-star atmosphere model with documented stellar parameters and observed
   NUV constraints for a new physical model.
2. If those are unavailable, construct C as a clearly labeled provisional
   solar-XUV plus scaled-WASP-17 long-wavelength composite.
3. Keep A as a blackbody comparison. Keep B for diagnosing the consequences of
   its spectral construction, not as an equally credible reference spectrum.
4. Verify component provenance and normalization before joining them. Do not
   force continuity by arbitrarily renormalizing a component that has an
   independent physical normalization.
5. Compare cross-section-weighted ionization and heating integrals, not just
   broad-band energy totals. Resolve the absolute XUV discrepancy separately
   from the photospheric choice.

The measured adjacent-node ratio across B's 1700 A join is `0.016393`, a drop
by about 61. This confirms the stated discontinuity. Candidate C has not been
built or tested here, so its own join and rates remain unverified.

## 7. Decision 5: convert validated inputs without overwriting scientific history

### 7.1 Separate the reference experiment from the best available stellar model

For **HD 189733 b**, the two files answer different questions:

- Use the Salz-normalized eps Eri proxy for a deliberately defined comparison
  with Salz et al. (2016), while labeling that its spectral shape is a proxy,
  not the original Salz spectrum. Equal integrated luminosities do not establish
  equal photoionization or heating profiles.
- Prefer the Bourrier reconstruction for a model intended to represent the
  observed HD 189733 stellar environment described by that study. Select the
  relevant visit for an epoch-specific comparison. Use the B-E mean only when
  the averaging assumption is part of the experiment, with variability treated
  as an input uncertainty.

Bourrier et al. (2020), Section 4.4, reconstructs the XUV spectrum using
observational constraints and a differential emission measure model. It should
not be described as a directly measured EUV spectrum. The local final PDF and
Table 3 also use an EUV interval of 62-912 A, different from the 100-912 A
division used elsewhere in the project. Compare like bands. See the
[published paper](https://doi.org/10.1093/mnras/staa256) and
[local PDF](../../references/Bourrier_2020MNRAS_493_559.pdf).

The two stored files have identical sampled long-wavelength flux in the
1700-2583 A diagnostic band, while their measured 100-912 A fluxes are
`1.51308e4` and `2.07046e4 erg cm^-2 s^-1`, respectively. Their 912-1110 A
fluxes are `2.63238e3` and `1.84130e3`. These are checks of the stored tables,
not a determination that either is the correct spectrum for every observation.

For **WASP-121 b**, permit a provisional C experiment only with the qualifications
in Section 6. Do not hold all other planet configurations while resolving that
normalization, and do not label the unresolved case a completed literature match.

### 7.2 Conditions for configuration changes

Before switching each configuration, record its purpose, reference model,
SED identifier and checksum, wavelength coverage, flux units at the orbit,
distance/normalization assumptions, active absorbers, and geometric dilution.
Check both low- and high-energy coverage; a file opening successfully is not a
check of the radiation actually used by its consumers.

Keep the previous output products intact and mark them as generated with the
previous input and SED. Changing a power law to a loaded spectrum changes the
physical experiment, even if only a configuration file changes. Do not present
old outputs as results of the new configuration. Atmospheric reruns remain
separate from the configuration update, as the user previously required.

Preserve `VULCAN/atm/stellar_flux/sflux-HD189_B2020.txt`. Its removal is unnecessary
for this work and could affect another code. If its provenance is unresolved,
do not select it as EXHALE's reference input. The decision document's claimed
21 percent discrepancy for that file was not remeasured in this review.

## 8. Additional important issues to resolve

### 8.1 Final steady-state certification is weaker than its description

**Source-confirmed:** `steady_residual.f90:746-797` tests a hydrodynamic residual,
a mass-flux spread, and an optional carrier residual. It receives no chemical
root classification. The marching caller at `EXHALE_main.f90:1419` explicitly
passes `.false., 0.0d0` for the carrier condition.

The steady trial policy in `steady_newton.f90:647-650` allows the number of
uncertified cells to remain unchanged. That can be a continuation heuristic,
but it is not final chemical certification: an unchanged count can also hide
larger residuals or different failed cells. The class-5 mismatch in Section 4.3
further undermines the meaning of that count.

**Recommendation:** immediately before declaring a stationary physical result,
reevaluate the active equations on the exact state to be written. Require valid
chemical closures, residuals for every active independent transport equation,
energy balance, and the relevant invariants. A transported species needs a
stationary residual even when it is not a global Newton unknown. Do not use the
last backward-Euler residual as a substitute for its steady equation.

This review establishes the missing conditions in the inspected gate and
callers. It does not claim a newly executed false-convergence case.

### 8.2 Energy floors and two-iteration solves require explicit failure handling

**Source-confirmed:** `energy_semi_implicit.f90:210-295` takes two updates,
refreshes cooling after the first, and returns the second temperature with the
previous cooling. Heating and the cooling derivative are frozen. The derivative
uses `abs(dC_dT)` as a stabilizing iteration coefficient rather than the exact
Newton derivative on a falling cooling branch.

Such an iteration can be a valid solver strategy only when its final equation
residual is checked. The fact that a modified iteration has the same formal
fixed point does not prove that two iterations reach it. The same routine clamps
`T_trial` to `0.01` in code units and then constructs the returned energy.
Counting that event does not account for the associated energy modification.

**Recommendation:** first implement a residual-controlled temperature solve
with a physical bracket or safeguarded iteration, evaluate sources at the
returned temperature, and reject an unacceptable update. Next implement the
coupled temperature-composition source step required by the plan. If a temperature
floor represents external heating, specify and budget that heating; otherwise
the floor cannot certify the physical step.

### 8.3 Input coverage must follow the current physics, not old constants

**Source-confirmed:** the current `parameters.f90:768,791` uses
`e_th_HI = 13.598434599 eV` and `e_th_HeTR = 4.767775 eV`. The current
`J_inc.f90` defines the hydrogen n=2 threshold as `e_th_HI/4`, and the photon
grid can extend to that threshold when Balmer coupling is active.

The SED README and `inputdata/sed/build/sed_coverage_and_order_check.py` still
describe 4.80 eV / 2583 A for the triplet and a metal floor near 2856 A as
sufficient for every listed configuration. The actual triplet edge is near
2600 A; active Balmer photoionization requires coverage near 3647 A. The reader
and Balmer routine now contain coverage checks for the latter, so the missing
piece is not simply an absent Balmer guard: it includes stale inventory and test
assumptions. A table ending near 2995 A cannot supply the full active Balmer band.

**Recommendation:** derive coverage tests from the actual active threshold
inventory, and check each configuration. Also include the upper photon-energy
limit when comparing integrated XUV values: flux in 1-10 A is not included in a
run whose upper energy cutoff corresponds to 10 A.

The diagnostic in this review deliberately retains 2583 A for comparison with
the decision document's stated NUV band. It is not a production triplet-coverage
test and does not endorse that value as the current threshold.

### 8.4 Threshold splitting does not by itself establish one radiation field

**Source-confirmed:** `sed_read.f90:367-466` converts stored `F_lambda` samples
to `F_E`, assigns histogram values to bins with geometric energy edges, and
splits those bins at thresholds. `J_inc::stellar_flux_eV` instead interpolates
the table nodes in log-log space where the flux is positive. The Balmer
photoionization and heating routines integrate that interpolated field on their
own 400-point grid.

These constructions can approximate the same smooth spectrum as resolution
improves, but they are not identical finite-resolution representations. Splitting
a histogram bin preserves its existing energy contribution; it does not prove
that the histogram reproduces the original wavelength-bin energy, photon number,
or a different consumer's interpolated integral. Table nodes and averages over
finite wavelength bins also need distinct meanings.

**Recommendation:** specify whether the input values are point samples or bin
averages, preserve supplied edges when necessary, and use a consistent spectral
reconstruction with controlled quadrature errors. Test a narrow line, a sharp
threshold, a coarse bin, and a steep continuum through the actual ionization,
heating, and Balmer consumers. Require convergence of rates as resolution
increases. No production quadrature error was measured in this review.

### 8.5 Radiation, chemistry, and EOS must share energy recipients

Do not let solver-policy work obscure the remaining molecular physics. D0 still
identifies mismatched H2 internal populations between EOS and chemistry,
unassigned energy in H2 photochannels, H3+ low-density behavior, and an incomplete
oxygen reaction ledger. These were not all rederived or retested here and remain
open audit items, not newly confirmed numerical measurements.

For Phase 3, each reaction/photoevent must identify the initial and final
chemical states, excitation energy, prompt kinetic energy, electron degradation,
and emitted radiation. A value assigned to one recipient cannot also be charged
as another source. Apply the same rule to Lyman-alpha excitation and de-excitation
rather than inferring the thermal loss from escape probability alone.

### 8.6 Do not encode geometry twice in the SED and transfer calculation

The published Huang methods distinguish the incident spectrum at the planet from
the atmospheric irradiation geometry. A radial model's scalar dilution and
column prescription are also not automatically the same as a calculation with
an explicit oblique ray. State which factors belong to the orbital SED, which
multiply the incident beam, and which change optical depth. Check them once
through the actual consumer before using a normalization discrepancy to select
another SED. No new geometry equivalence or atmospheric correction factor was
derived here.

## 9. Proposed implementation order and acceptance checks

The following work order is recommended, not authorization to make all of these
changes in the current review task.

1. **Repair acceptance contracts:** reconcile class 5, remove the constrained-cell
   failure exemption, and separate a numerical intermediate from an accepted
   physical state. Establish independent final chemical and species gates.
2. **Approve the target governing system:** settle the energy ledger, invariant
   species space, current/mass-flux closure, CO exclusion, boundary conditions,
   and the complete rollback state. Correct D0's contradictory descriptions.
3. **Implement bounded recovery and the residual-controlled thermal solve:**
   verify state restoration and final-state energy residuals before expanding
   the source coupling.
4. **Implement the coupled species-energy update and consistent stationary
   equations:** use the same physical terms and compatible face fluxes.
5. **Validate and select spectra independently:** prepare C as provisional if
   desired, resolve the Huang normalization separately, and update only the
   configurations whose experiment and coverage have been checked.
6. **Run affected physical cases after their prerequisites pass:** keep old
   products, record the input/model differences, and request the separately
   scoped scientific reruns when needed.

| Focused check | Required result | Status in this review |
|---|---|---|
| Failed trial with mutated carriers and caches | No rejected state survives; accepted physical time is correct. | Proposed, not run. |
| Constrained worst cell plus unconstrained failed row | Cannot certify the failed carrier solve. | Deficient branch confirmed in source; injected test not run. |
| Valid constrained-continuation root | Same validity classification as another root of the requested equations. | Class-5 contract mismatch confirmed in source; runtime test not run. |
| Small hydro residual with failed chemistry or drifting species | Cannot be called a stationary physical solution. | Gate limitations confirmed in source; constructed run not performed. |
| Closed reacting cell | Conserves thermal plus chemical/excitation energy under the declared model. | Proposed; no source integration run here. |
| Returned scalar temperature | Sources and residual belong to the returned state; an unbudgeted floor cannot pass. | Current two-iteration behavior confirmed in source. |
| Stored SED diagnostic bands | Reproducible integrals and documented differences from published values. | Executed successfully; numerical results in Sections 6 and 7. |
| Production photon integrals | Consistent band coverage, energy, photon number, rates, and resolution convergence. | Not run; independent table integration is not a substitute. |
| Candidate C atmosphere | Validated input and a converged physical solution under that input. | C was not constructed; no atmospheric result claimed. |

## 10. Reproducibility and source records

Run the read-only measurement from the project root:

```bash
python3 docs/audit_20260905/decision_sed_check_20260906.py
```

The script uses the Python standard library, checks wavelength ordering and
requested-band coverage, interpolates the band endpoints, and prints SHA-256
hashes of every spectrum it reads. It computes diagnostic integrals without
modifying the original data. It completed with exit code 0 after its join-index
selection was checked against the actual 1699/1700 A rows of candidate B.

The measured SED identifiers were:

```text
wasp121b_solar_huang2023.txt
e504432bc16e51f8944a19c76be07db7de85b2f8dc2d484614c6460ab88e830a
wasp121b_wasp17_huang2023.txt
1e20037210a908cd6ef7aba27331bb3800bd3fdf106b7ccc5bff62e7913bdc02
hd189733b_epseri_salz2016.txt
dd0764ffadd986d969f2005ca061e4247812e7e4c48c804f663b0a43cc7b0cd4
hd189733b_bourrier2020.txt
da76f1fd69920f6aba33fb4f5503ebec04459eca19a34dcb981d5bee2bf3a1af
```

Selected document/source hashes, sampled at September 6, 2026, 08:18:14 KST:

```text
docs/To_be_determined_by_user_20260906.md
c61438c4b0a60befed0a6646b4b94052a8aa822987e5777ef6187939d3f8d18a
docs/development_plan_20260905_rev3.md
1febaa0017e9f3374719b57586d7bf07cdd26d1d2e4a822616a6e3a635410c8c
docs/d0_governing_system_20260906.md
3339ca2eeca362e38a10d44b08cde6ae902b4de2d9c733729a73739597d16f33
src/modules/time_step/steady_newton.f90
4efff857e10f0c022b9db599c98625b0007d20bd9cc2175dead00ad5a624587d
src/modules/lower_atmosphere/diffusive_photochemistry.f90
43be19aa481a875c3d30e2d4b966974f064dc279c92a250a2b7c5a69118ecbca
src/modules/radiation/ionization_equilibrium.f90
8e5dec6fa2a5af22ee2f7a9d9658b3a927bcd2524c8d528400b02699b0fc4914
```

These hashes identify the inspected working files, not a new build or a
regression result. Recheck affected findings if the ongoing implementation
changes them.

## 11. Suggested decision statement

Approve bounded rejection and retry with a complete rollback contract. Do not
accept non-root chemical states as completed physical steps or final stationary
solutions. Reconcile the root-status interface and require independent chemical,
species, and energy acceptance conditions. Accept corrected D0 as an audit;
require its missing physical derivations before declaring the Phase 2 gate met.
Retain thermal plus kinetic energy only with a conservative chemical-energy
ledger and coupled source update. Exclude physically unvalidated CO-ceiling (superseded 2026-09-06: the ceiling
is deleted and CO carries published rates) and incomplete oxygen-energy
configurations from validated results. Permit a
provisional WASP-121 b composite C with explicit proxy and normalization
uncertainties, while withholding an exact Huang-reproduction claim. Select
HD 189733 b spectra according to the intended reference or observational epoch.
Switch validated inputs without regenerating or overwriting existing scientific
products, and leave unrelated VULCAN data untouched.
