# Implementation Review of the LHS 1140b Lower-Atmosphere Plan

- **Review date:** 2026-08-22
- **Document reviewed:** `docs/lhs1140b_lower_atmosphere_plan.md`
- **Review basis:** Current EXHALE, Photochem, and p-winds source trees in this
  workspace

> [2026-08-27] Kept unchanged as the record of what the trees held on
> 2026-08-22.  Two of the modules it inspects have since moved:
> `src/modules/functions/species_diffusion.f90` was deleted on 2026-08-25 and
> replaced by `src/modules/functions/binary_element_diffusion.f90` (Phase D),
> and the lower-atmosphere handoff gained the profile route of Phase E
> (`src/modules/files_IO/lower_atmosphere_profile.f90`).  The plan of record
> is `lhs1140b_lower_atmosphere_plan_new.md`, whose Phase D/E/F entries carry
> the status of each item.

## 1. Scope and method

This review evaluates whether the proposed LHS 1140b workflow is supported by
the current implementation and whether the proposed work packages are
physically sufficient. It does not assess the plan only from its description.
The relevant definitions, callers, parsers, output consumers, and existing run
records were inspected.

The review followed these distinctions:

- A routine can exist without being connected to a user-facing execution path.
- A parser can recognize a quantity without enforcing its physical
  consistency.
- A lower-atmosphere model can produce a state without determining the escape
  flux needed to close the coupled problem.
- A transit routine can calculate a line profile without calculating the
  observational metrics named in the plan.
- An existing diffusion term can be valid in a trace limit without being valid
  at arbitrary composition.

The paper values quoted in the plan were not independently remeasured from the
published PDFs during this implementation review. They are treated here as
quoted target values. No new EXHALE, Photochem, or p-winds science simulation
was run.

## 2. Executive assessment

The scientific objective is well motivated, but the plan is not yet ready to
serve as an implementation specification. The largest problem is that the
proposed chain is not closed:

```text
lower atmosphere -> H production -> diffusive separation -> wind
```

The current handoff supplies a prescribed He/H ratio to EXHALE. It does not
carry an elemental hydrogen flux, and it does not match the lower model's upper
boundary flux to the wind's lower boundary flux. Therefore, He/H cannot become
a self-consistent output merely by adding a Photochem adapter.

Three other issues must be corrected before science runs:

1. The current normalization is not globally based on hydrogen nuclei as the
   plan states. The hydrodynamic density scale already uses H+He nuclei.
2. Replacing the hydrogen background in the existing trace-species diffusion
   kernel is not a valid arbitrary-composition diffusion model. A conservative
   binary or multicomponent formulation is required.
3. `EXHALE_transit.py` does not consume the base radius resolved through
   `base.inp`, and it does not calculate the astrophysical FWHM or the blue/red
   triplet ratio requested by WP4.

The external SED path is already exercised by a converged benchmark, contrary
to WP3. The LHS 1140b SED still needs new validation because metastable helium
requires spectral coverage down to 4.8 eV, beyond the usual XUV-only band.

## 3. Finding 1: the normalization premise is inaccurate

### 3.1 What the plan states

Section 3, Gap 1 states that EXHALE "normalizes everything to hydrogen nuclei"
and identifies

```text
mass_per_H = 1 + 4*HeH
```

as the relevant scaling. This description is incomplete for the current code.

### 3.2 What the implementation does

`src/modules/files_IO/input_read.f90` reads `He/H number ratio` into `HeH`, but
the derived normalization uses the composition module:

```fortran
mass_per_H = comp_mass_per_H()
ntot_bc    = comp_ntot_bc()
rho_bc     = comp_rho_bc()
```

The comments at `input_read.f90:961-975` explicitly define `n0` as an H+He
nuclei scale and state

```text
n_H = n0/(1 + HeH).
```

`src/modules/functions/composition.f90:188-191` defines the mass-density
factor as

```fortran
comp_rho_bc = comp_mass_per_H()/(1.0d0 + HeH)
```

The numerator also includes configured metal mass when the metal equation of
state is active. Therefore, `1 + 4*HeH` is not the whole composition policy.

For `HeH = 1000`, direct evaluation gives

```text
n_H/n0              = 9.99000999e-4
mass_per_H           = 4001.0       (H/He only)
rho/(m_H n0)         = 3.997003
```

The large value of `mass_per_H` does not make the physical density scale
diverge because `rho_bc` also divides by `1 + HeH`.

### 3.3 What still requires a trace-hydrogen audit

The correction above does not make WP0 unnecessary. Several species fractions,
initial guesses, and diagnostic quantities are still expressed relative to the
hydrogen-nuclei density. Examples of unguarded warm-start divisions include:

- `src/modules/radiation/ionization_equilibrium.f90:596-606`
- `src/modules/radiation/ionization_equilibrium.f90:623-632`
- `src/modules/post_process/post_process_adv.f90:438-448`
- `src/modules/post_process/post_process_adv.f90:489-491`

At `HeH = 1000`, these divisions are not divisions by zero solely because of
the abundance ratio. They can nevertheless become poorly conditioned if a
cell's computed H density approaches a numerical floor. The audit must follow
all paths that construct a Newton initial guess, apply the advection correction,
or report an abundance relative to H.

### 3.4 Cooling and electron accounting already present

The plan describes the cooling as H-centric and lists helium recombination,
free-free emission, and triplet lines as channels that may need to close the
budget. The current cooling implementation is broader than that description.
`src/modules/radiation/util_ion_eq.f90:604-705` includes:

- H II, He II, and He III recombination cooling;
- H I, He I, He II, and He 2^3S collisional-ionization cooling;
- free-free cooling from H and He ions, including the He III charge factor;
- H I, He I, and He II collisional-excitation cooling;
- He 10830 A and singlet-conversion terms when the triplet state is active.

The required question is therefore not simply whether helium cooling routines
exist. It is whether their sum closes the energy equation in the helium-rich,
low-XUV solution and whether assumptions such as photon escape remain valid.

### 3.5 Required revision to WP0

WP0 should be reframed as an audit of the remaining H-referenced algebra rather
than an audit based on a globally H-based hydrodynamic normalization.

Recommended acceptance criteria are:

1. Run the atomic H/He model at a sequence such as `HeH = 1, 10, 100, 1000`
   under floating-point and bounds checking.
2. Confirm finite Newton initial guesses and converged species fractions in
   every cell.
3. Confirm charge neutrality from the actual H, He, and electron densities.
4. Confirm that the reported heating and cooling channel sums reconstruct the
   energy source term.
5. Report the dominant cooling and heating channels rather than assuming them.
6. Compare a run with and without the advection correction to isolate any
   failure specific to `post_process_adv.f90`.
7. Record every applied density floor and demonstrate that the converged result
   is insensitive to a reasonable change in that floor.

## 4. Finding 2: the proposed diffusion generalization is underdefined

### 4.1 Current transport equation

`src/modules/functions/species_diffusion.f90:4-16` defines the transported
element ratio as

```text
fHe = n_He/n_H
```

and advances He relative to a hydrogen background. The conservative face flux
implemented at `species_diffusion.f90:136-156` is

```text
J_He = n_He v - D n_H [d(n_He/n_H)/dr + (n_He/n_H) G].
```

The H density is lagged from the current species state. After the solve, the
code rescales H and He fractions to restore the local mass normalization. A
reservoir cap enforces `n_He/n_H <= HeH`.

The same kernel is reused for each metal element at
`species_diffusion.f90:197-255`. Each metal is treated independently relative
to H and does not feed back on the H background.

This is an explicit trace-species formulation. It is not an
arbitrary-composition binary diffusion solver.

### 4.2 Why the two proposed alternatives are not equivalent

WP1 proposes making the background either the dominant species or, "cleaner,"
the total gas. These are not interchangeable implementation choices.

Switching between H and He as the named background preserves a trace-species
form only in the corresponding dilute limit. Substituting total gas into the
existing ratio and face coefficients does not enforce the equal and opposite
diffusive fluxes required by a binary mixture. It can change the bulk mass flux
and need not reproduce either dilute limit correctly.

A physically valid arbitrary-composition formulation should evolve a conserved
element variable, such as one independent mole or mass fraction, with a binary
Maxwell-Stefan flux or an equivalent formulation. It must enforce:

```text
sum_i rho_i w_i = 0
```

for the diffusive velocities relative to the chosen bulk velocity. Gravity,
the ambipolar electric field, thermal diffusion, and eddy diffusion must be
defined consistently with that bulk velocity.

### 4.3 Required diffusion tests

The proposed isothermal profiles in both orientations are necessary but not
sufficient. WP1 should require:

1. zero-flux diffusive equilibrium in an isothermal hydrostatic column;
2. recovery of the current He-in-H result as `n_He/n_H -> 0`;
3. recovery of the H-in-He trace limit as `n_H/n_He -> 0`;
4. conservation of total H nuclei and total He nuclei in a closed column;
5. zero net diffusive mass flux at every face;
6. preservation of a spatially uniform mixture when gravity and imposed
   gradients are absent;
7. convergence with grid spacing and diffusion timestep;
8. a homopause test with both molecular and eddy diffusion;
9. a moving-wind test that verifies the total elemental face flux rather than
   only the local abundance profile;
10. tests for neutral, partially ionized, and strongly ionized limits if the
    ambipolar correction remains active.

### 4.4 Crossover mass and heavy-element drag do not follow automatically

The plan states that crossover mass and O/C/N drag would become code results
after WP1. That conclusion is too strong.

The current metal loop treats each element as a trace constituent relative to
H. It does not solve a coupled multicomponent momentum system, and the response
of the dominant H/He mixture to metal drag is absent. Reversing the H/He
orientation alone cannot establish a physically correct crossover mass for O,
C, or N in a helium-dominated wind.

The plan should either:

- limit WP1 to validated H/He binary separation and retain crossover mass as a
  separate analytic diagnostic; or
- define a later multicomponent transport package that includes H, He, and the
  relevant heavy elements with a tested momentum or diffusion closure.

### 4.5 Molecular chemistry is currently incompatible with He diffusion

`src/modules/files_IO/input_read.f90:895-905` rejects a run when both molecular
chemistry and `He_diffusion` are enabled. This directly conflicts with the
planned sequence through a molecular lower boundary and diffusive separation.

WP1 or WP2 must define how this incompatibility is removed. Possible designs
include a single transport formulation valid throughout the matching region or
a carefully defined transition above the molecular layer. Merely completing
the atomic H/He kernel will not make the full chain executable.

## 5. Finding 3: the lower-atmosphere handoff does not close the problem

### 5.1 Current analytic lower column

`src/modules/lower_atmosphere/lower_column.f90` implements an isothermal
hypsometric column with the H2/H/He equilibrium fit described in its header.
The source states that the fit was validated from 1000 to 2500 K. Its current
purpose is to estimate the base radius and H2/H/He partition for hot or warm
H/He atmospheres.

The plan is correct that this routine is not a physically adequate climate and
cold-trap model for a 226 K atmosphere. The precise limitation is its
temperature range, isothermal structure, equilibrium H2 chemistry, and lack of
water condensation and photochemical transport, rather than only its use for
hot Jupiters.

### 5.2 Current automatic lower-atmosphere execution path

`src/modules/files_IO/input_read.f90:1228-1280` recognizes only:

```text
Lower atmosphere: analytic
Lower atmosphere: vulcan
Lower atmosphere: none
```

There is no Photochem execution adapter in the EXHALE source tree. The existing
Tier-3-style automated path runs VULCAN and converts its result with
`src/utils/vulcan_to_base.py`.

The existence of the adjacent `photochem` source tree is not a callable EXHALE
entry point. A new driver, configuration generator, result reader, failure
policy, and reproducibility record are required.

### 5.3 Current `base.inp` contract

`read_base_inp` recognizes only:

```text
T_base
r_base
HeH_base
Kzz_base
q_H2_base
p_base
```

Unknown keys are ignored. Existing VULCAN output files may include H2O, CO,
CH4, CO2, NH3, and HCN in comments, but EXHALE does not consume those values.

The parser also states that `q_H2_base` and `HeH_base` must originate from the
same lower-atmosphere solution, but it only prints this requirement. It does
not reject an inconsistent pair (`input_read.f90:1206-1217`).

The current contract cannot carry:

- the water abundance at the cold trap;
- the atomic-H production rate;
- an upward elemental-H flux;
- a vertical composition profile across the matching region;
- a temperature or Kzz profile;
- condensation or rainout information;
- the lower model's top boundary condition;
- uncertainty information needed to propagate the paper's input constraints.

### 5.4 Photochem does not automatically supply the claimed climate solution

The adjacent Photochem tree includes climate capability, but the chemistry
atmosphere initialization APIs require supplied temperature and eddy-diffusion
profiles. For example,
`../photochem/photochem/cython/EvoAtmosphere.pyx:153-205` accepts `temperature` and
`edd` arrays as inputs.

The current settings parser in
`../photochem/src/photochem_types_create.f90:102-107` rejects the obsolete
`evolve-climate` setting. The same file states that old water-fixing and
water-condensation settings are unsupported and that H2O condensation must be
represented through ordinary boundary conditions and an H2O particle in the
reaction mechanism (`photochem_types_create.f90:144-158`). Gas rainout also
requires a specified tropopause altitude (`photochem_types_create.f90:160-195`).

Thus, WP2 must explicitly define:

1. which climate API or external radiative-convective model produces T(P);
2. how that profile is passed to Photochem chemistry;
3. the reaction mechanism and H2O particle treatment;
4. the lower water inventory and surface boundary conditions;
5. the tropopause and rainout configuration;
6. the Kzz profile;
7. the photochemical upper boundary condition;
8. how the chemistry solution is iterated with the escaping wind.

### 5.5 Why He/H remains an input without flux matching

The plan expects the lower model to determine how much hydrogen reaches the
wind and expects the resulting EXHALE calculation to return He/H. A steady
lower-atmosphere chemistry model cannot determine this value independently of
its upper boundary flux. Conversely, the wind's lower composition depends on
the material supplied through that boundary.

The coupled problem therefore needs at least one flux-continuity condition. A
possible iterative structure is:

```text
1. Supply a trial elemental-H flux to the Photochem upper boundary.
2. Solve climate and chemistry to the matching pressure.
3. Pass the matching state and elemental flux to EXHALE.
4. Solve the wind and measure its lower-boundary H and He fluxes.
5. Update the Photochem upper boundary until elemental fluxes agree.
```

The exact algorithm requires a separate design decision, but some equivalent
closure is necessary. Without it, `HeH_base` remains a prescribed boundary
composition and the WP2 acceptance criterion is circular.

### 5.6 Recommended handoff quantities

Before changing a public input format, the project should define the matching
surface and the ownership of each quantity. A complete contract is likely to
need:

- matching pressure and radius;
- temperature and total number density;
- elemental H and He abundances;
- molecular and atomic H partition;
- upward elemental H and He fluxes;
- electron or charge information if the matching layer is ionized;
- Kzz at the match and, if needed, a short profile around it;
- selected molecular abundances needed by the EXHALE network;
- a source-model version and configuration fingerprint.

Whether this should remain a scalar `base.inp` file or become a profile file is
an API and data-layout decision that should be confirmed before implementation.

## 6. Finding 4: the molecular-network description needs correction

`src/modules/lower_atmosphere/mol_rates.f90:132-173` defines reactions R16
through R23. However,
`src/modules/nonlinear_system_solver/System_HeH_mol.f90:21-25` states that R21
and R22 are excluded from the molecular residual to preserve the atomic limit.
The system instead imports the shared H-He charge-exchange implementation from
`charge_exchange`.

The molecular coefficient setup stores R16-R20 and R23, not R21 or R22
(`System_HeH_mol.f90:55-90`). Therefore, the proposed audit should cover:

- R16-R20 for HeH+ and He+/H2 chemistry;
- R23 for H2 + He+ charge exchange;
- the shared H-He atomic charge-exchange path;
- He 2^3S + H2 Penning ionization;
- the electron closure when He supplies most electrons;
- the H-nucleus closure when H-bearing species are trace constituents.

Referring only to "R16-R21" would miss R23 and would incorrectly imply that
R21 is directly evaluated in the molecular residual.

## 7. Finding 5: the external SED path is already exercised

### 7.1 Existing implementation and execution evidence

`src/modules/radiation/sed_read.f90` reads a two-column numerical spectrum,
checks the selected energy coverage, requires at least two points, and reports
missing low-energy coverage. `src/modules/init/set_energy_vectors.f90:127-145`
calls this reader when `Spectrum type` requests a file.

The WASP-52b benchmark already selects

```text
Spectrum type: Load from file..
Spectrum file: ../../WASP-52b/eps_eri_sed_fxuv1p0.txt
```

in `benchmarks/wasp52/input.inp:11-12`. Its existing
`run_20260815_kb.log:17-18` records successful numerical-spectrum reading, and
the setup report records an external SED and a converged solution.

WP3 should therefore not state that all current regression cases use a power
law or that LHS 1140b will be the first production exercise of this path.

### 7.2 Spectral range required by metastable helium

When He 2^3S is enabled, `sed_read.f90:37-38` lowers the minimum selected photon
energy to the triplet threshold, 4.8 eV. This corresponds to a wavelength near
2590 A. Ground-state He and H additionally require coverage across their own
ionization edges.

An input described only as an XUV SED may reproduce the high-energy heating
while omitting photons that photoionize metastable helium. That omission can
directly change the 10830 A population even if the integrated XUV flux is
correct.

### 7.3 Recommended WP3 checks

WP3 should require:

1. wavelength coverage from at least the active low-energy threshold through
   the configured X-ray upper bound;
2. explicit treatment of gaps and overlap among the GJ 1132, GJ 699, and XMM
   components;
3. a documented absolute normalization at the LHS 1140b orbit;
4. consistency between integrated X-ray/EUV luminosities and the setup report;
5. monotonic wavelength ordering and sufficient sampling near the H I, He I,
   He II, and He 2^3S thresholds;
6. a sensitivity test for the 4.8-13.6 eV band;
7. preservation of the exact resolved SED file in the run record.

## 8. Finding 6: the transit validation cannot currently evaluate all targets

### 8.1 The transit calculation does not read the resolved base radius

`EXHALE_transit.py:80-93` obtains planet parameters through
`exhale_transit_lib.read_input_params`, which reads only `input.inp`.
`exhale_transit_lib.py:307-331` returns the original planet radius and
equilibrium temperature. It does not read `base.inp` or
`EXHALE_setup.out`.

The wind solver reads `base.inp` before constructing derived normalization
constants. Therefore, a lower-atmosphere handoff can change the radius used by
the wind without changing the radius used by the transit calculation.

A direct check of the existing WASP-52b benchmark measured:

```text
Quantity             input.inp       base.inp
Planet radius        1.27000 R_J     1.40492 R_J
Temperature          1304.00 K       1181.20 K
He/H                  0.020408        0.095919
```

The stale `T0` value is currently printed but is not otherwise used by the
transit script. The radius mismatch is physical: `Rp` enters radial gradients,
line-of-sight chord lengths, absorbing areas, orbital scaling, and rotational
velocity calculations throughout `EXHALE_transit.py`.

Before using a lower-atmosphere handoff, the transit script must read the same
resolved planet radius used by the wind. The preferred source is a structured
resolved-configuration output written by EXHALE, not a second independent
interpretation of both input files.

### 8.2 Requested observational metrics are not implemented

The code calculates the minimum transmission of the full He triplet and saves
the theoretical, instrument-convolved, and rotation-plus-instrument curves.
It does not calculate:

- the astrophysical FWHM of the modeled red feature;
- the blue-component depth;
- the blended-red depth under a defined window;
- the blue/red depth ratio;
- uncertainties or grid sensitivity for those metrics.

`FWHM_HeTR` at `EXHALE_transit.py:398` is the Gaussian line-spread width
derived from instrument resolving power. It is not the FWHM of the modeled
absorption feature.

WP4.3 requires a new metric extractor with documented continuum treatment,
wavelength windows, interpolation, and instrument resolution. The same
definition used for the published observational values must be applied to the
model curve.

## 9. Finding 7: the p-winds oracle is not yet a reproducible case

The adjacent p-winds tree contains the public calculation routines, but no LHS
1140b driver, frozen input configuration, or reproduction artifact was found.
The three scalar values named in WP4.1--mass-loss rate, isothermal temperature,
and H/He--are insufficient to reproduce the published profile.

A useful oracle case also needs:

- the exact p-winds revision or the changes used for the paper;
- the lower and outer radii and the radius grid;
- lower-boundary ion fractions and density normalization;
- the exact stellar spectrum used for ionization and triplet destruction;
- stellar and planetary parameters;
- the line-profile and instrument-convolution settings;
- any wavelength or velocity shifts applied before comparison;
- a saved reference profile and numerical tolerances.

The first oracle acceptance criterion should be reproduction of the published
p-winds line profile and reported metrics with a frozen local driver. Only then
should EXHALE profiles be compared against it.

## 10. Finding 8: the fluid-validity conclusion is too strong

No Knudsen-number or exobase implementation was found in the EXHALE source
path. WP4.5 therefore describes a new diagnostic, not the use of an existing
one.

For a partially ionized helium-rich wind, the mean free path is not defined by
one unspecified cross section. The diagnostic should consider, as appropriate:

- neutral He-He collisions;
- neutral H-He collisions;
- ion-neutral momentum transfer;
- Coulomb coupling among ions and electrons;
- species coupling times compared with the flow time.

The physically correct conclusion when the sonic point lies above the exobase
is that the continuum hydrodynamic solution is not validated through its
critical point. It does not prove that the hydrodynamic mass-loss rate is an
overestimate. The sign and magnitude of the error require a kinetic or
transitional-flow comparison.

Recommended acceptance criteria are:

1. define the collision model and cross-section sources;
2. calculate species-relevant Knudsen numbers across the heating and
   acceleration region;
3. compare collision, ionization, recombination, and flow times;
4. locate the exobase and critical point with stated definitions;
5. classify the hydrodynamic mass-loss rate as validated or unvalidated rather
   than assigning an unsupported error direction;
6. use a kinetic escape estimate or model if the critical region is not
   collisional.

## 11. Recommended replacement work structure

The original order allows several later tasks to begin before their inputs and
metrics are defined. The following order reduces that ambiguity.

### Phase A: freeze the target and reference configuration

- Verify the paper values against the published article and supplement.
- Recover or reconstruct the exact p-winds configuration.
- Define the observational line metrics and instrument settings.
- Construct the full stellar spectrum, including the 4.8-13.6 eV band.
- Record all source files, revisions, scaling rules, and uncertainties.

**Acceptance:** A frozen p-winds driver reproduces the reference line profile
and the chosen depth, FWHM, and blue/red metrics within stated tolerances.

### Phase B: make the wind and transit configuration consistent

- Write a structured resolved-configuration product from EXHALE.
- Make `EXHALE_transit.py` consume the resolved base radius and other relevant
  parameters.
- Implement observational FWHM and blue/red metric extraction.
- Add a focused test in which `base.inp` changes the planet radius.

**Acceptance:** The wind and transit reports show the same resolved radius, and
the metric extractor recovers synthetic profiles with known answers.

### Phase C: audit the helium-rich atomic wind

- Execute the revised WP0 abundance sequence under runtime checking.
- Audit the Newton initialization and advection-correction divisions.
- Verify charge and energy closure from channel-resolved outputs.
- Determine whether the triplet and line-cooling escape assumptions remain
  valid.

**Acceptance:** A converged `HeH = 1000` atomic run has finite state variables,
closed charge and energy budgets, documented active floors, and stable results
under resolution and floor changes.

### Phase D: implement conservative H/He separation

- Derive the arbitrary-composition binary transport equation.
- Define the bulk velocity and diffusive mass-flux constraint.
- Implement and verify the neutral and ionized limits.
- Resolve the current molecular-chemistry incompatibility.

**Acceptance:** The implementation passes conservation, zero-flux equilibrium,
both dilute limits, homopause, and moving-wind flux tests.

### Phase E: design and implement the Photochem coupling

- Select the climate calculation and define T(P), condensation, rainout, and
  Kzz.
- Create a reproducible LHS 1140b Photochem configuration.
- Define the matching surface and handoff data layout.
- Couple elemental fluxes iteratively across the matching surface.
- Propagate the lower-atmosphere state into the molecular and atomic wind
  regions.

**Acceptance:** The lower and upper models agree on elemental H and He fluxes at
the match, the cold-trap water abundance is documented, and He/H is insensitive
to the initial trial flux after coupled convergence.

### Phase F: science validation and physical applicability

- Run the prescribed-composition EXHALE/p-winds comparison.
- Run the coupled-composition case.
- Compare density, velocity, temperature, ionization, and metastable-He
  profiles.
- Evaluate the 2024 line metrics and the 2025 XUV sequence.
- Calculate the collisional-validity diagnostics.

**Acceptance:** Every reported observable traces to one resolved configuration,
and the applicability statement follows from the calculated collision and flow
scales.

## 12. Proposed acceptance matrix

| Area | Required result | Current support | Missing work |
|---|---|---|---|
| He-rich normalization | Finite and converged trace-H solution | H+He density normalization exists | Stress tests and H-referenced path audit |
| Energy balance | Closed channel-resolved budget | H and He channels exist | He-rich validation and radiative-assumption checks |
| H/He diffusion | Conservative arbitrary-composition transport | He-in-H trace kernel exists | Binary formulation, conservation tests, H-in-He limit |
| Molecular transition | Chemistry and separation in one chain | Molecular network exists | Molecular chemistry and He diffusion are incompatible |
| Climate and cold trap | Reproducible T(P) and H2O supply | Adjacent climate capability exists | Explicit climate-chemistry driver and condensation setup |
| Lower/upper coupling | Matched elemental fluxes | Scalar base-state parser exists | Flux contract and iterative closure |
| Stellar spectrum | Full active photoionization band | External SED reader exists and has run | LHS spectrum, 4.8 eV coverage, stitching tests |
| Transit geometry | Same resolved radius as wind | He triplet forward model exists | Read resolved base parameters |
| Line comparison | Depth, FWHM, blue/red ratio | Peak triplet absorption exists | Astrophysical metric extractor |
| p-winds reference | Reproducible paper configuration | Public routines exist | LHS driver and frozen reference outputs |
| Fluid validity | Defined collisional regime | Sonic-point tools exist | Knudsen, exobase, and coupling-time diagnostics |
| Heavy-element drag | Physically coupled O/C/N transport | Independent trace-metal diffusion exists | Multicomponent transport or limited analytic claim |

## 13. Specific edits recommended for the original plan

The following statements should be corrected before the plan is used for
implementation:

1. Replace "normalizes everything to hydrogen nuclei" with a description of
   the H+He hydrodynamic scale and the remaining H-referenced species algebra.
2. Replace "dominant species (or total gas)" with a requirement for a
   conservative binary formulation and a defined bulk velocity.
3. Remove the claim that the H/He diffusion rewrite alone calculates crossover
   mass and O/C/N drag.
4. Replace the statement that Photochem automatically performs the required
   climate and cold-trap calculation with an explicit climate-chemistry
   coupling task.
5. State that the current automatic lower-atmosphere route is VULCAN or the
   analytic column and that a Photochem adapter does not yet exist.
6. Replace the scalar-composition acceptance criterion with a flux-continuity
   criterion at a defined matching surface.
7. Correct the molecular reaction audit from R16-R21 to the reactions and
   shared charge-exchange path actually invoked by `System_HeH_mol`.
8. Remove the claim that all current regression cases use power-law spectra.
9. Add the 4.8-13.6 eV spectral band to the stellar-input requirements.
10. Add a prerequisite that the transit calculation consume the wind's
    resolved base radius.
11. State that FWHM and blue/red metrics require new code.
12. Replace the assertion that an exobase below the sonic point proves an
    overestimated mass-loss rate with a statement that the hydrodynamic result
    is then unvalidated without a kinetic comparison.
13. Change "No code has been written for this yet" to "No LHS 1140b-specific
    driver, configuration, or validation case has been written." Substantial
    supporting infrastructure already exists.

## 14. Validation performed for this review

The following read-only checks were performed:

- inspected the full plan document;
- followed `HeH`, `mass_per_H`, `ntot_bc`, and `rho_bc` through input parsing
  and composition normalization;
- inspected the H/He and metal diffusion face fluxes and rescaling;
- inspected the molecular-chemistry compatibility gate;
- inspected the analytic lower-column validity range and outputs;
- inspected the automatic lower-atmosphere dispatch and `base.inp` parser;
- inspected the current Photochem atmosphere initialization and water settings;
- inspected the molecular reaction definitions and their actual invocation;
- followed the external SED reader through its caller;
- verified an existing external-SED benchmark configuration and execution log;
- followed transit parameter parsing into geometry calculations;
- searched the transit output path for the requested line metrics;
- searched the EXHALE source for Knudsen-number and exobase diagnostics;
- directly compared the WASP-52b values read by the transit parser with the
  values in its `base.inp`;
- directly evaluated the H/He-only normalization factors at `HeH = 1000`.

No repository file was modified as part of the analysis that preceded this
review document. This document itself is the only new file created for the
requested write-up.

## 15. Validation not performed

The following work was intentionally not performed because it would exceed the
scope of a code-and-plan review:

- no EXHALE build or new hydrodynamic simulation;
- no `HeH = 1000` runtime-checking experiment;
- no Photochem build or LHS 1140b chemistry calculation;
- no p-winds reproduction run;
- no new SED construction;
- no transit-spectrum generation;
- no independent remeasurement of the paper values from the published PDFs;
- no full regression suite.

These items are not reported as passing. They belong to the revised work phases
and acceptance criteria above.
