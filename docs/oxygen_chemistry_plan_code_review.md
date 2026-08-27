# EXHALE oxygen chemistry plan: code review and revised proposal

> **2026-08-19, later the same day:** the plan of record is now
> `docs/oxygen_chemistry_new_plan.md`, which folds this review into the merged
> plan after re-verifying its code claims against the source (all checked
> claims held). Section 9 below (the recommended order) is superseded by that
> document's roadmap. This file remains the review record.

- Review date: 2026-08-19
- Document reviewed: `docs/oxygen_chemistry_options.md`
- Evidence used: the current EXHALE Fortran and Python implementation, plus the
  Photochem source tree in this workspace
- Scope: this is an implementation review and revised plan, not a record of
  completed chemistry changes

## 1. Main conclusions

The original document correctly separates the problem into two areas:

1. O/OH/H2O chemistry that controls the H2/H partition near 1 microbar
2. Observability and completeness of O I/O II/O III, which EXHALE already treats
   as atomic ions

Reading the code changes the priority and expected cost of several options.

1. **`q_H2_base` does not pin the H2 abundance at the lower boundary.**
   The value in `base.inp` passes through
   `composition::h2_bound_fraction` and affects only `ntot_bc`. It changes the
   particle count used by the base pressure and EOS normalization, but it is not
   a species boundary condition. `set_IC.f90` still initializes H2 with the
   chemical equilibrium fit, after which the local steady-state molecular solve
   in `ionization_equilibrium.f90` determines H2 again. Option A0 should therefore
   describe `q_H2_base` as an external photochemical EOS anchor, not an H2 pin.

2. **Photochem already provides a gas giant and hot Jupiter chemistry path.**
   `photochem/examples/GasGiants.ipynb` is a gas giant tutorial and explicitly
   models SO2 photochemistry on the hot Jupiter WASP-39b. It uses
   `zahnle_rx_and_thermo_files` to construct an H/He/N/O/C/S network and solves it
   with `EvoAtmosphereGasGiant`. Saying that Photochem's Zahnle mechanism is only
   for rocky planets is therefore too narrow. A more accurate description is
   that the source mechanism is derived from `zahnle_earth.yaml`, while Photochem
   supplies an official gas giant subset workflow and a hot Jupiter example.

3. **A2 is not a simple addition of two or three species.**
   Atomic oxygen is currently owned by the generic metal block, where total
   oxygen is defined as O I + O II + O III. Adding OH and H2O requires coupled H
   and O conservation, EOS changes, opacity, restart and output support,
   diffusion decisions, and stronger nonlinear-root checks. There must be one
   total oxygen budget shared by O/O+/O++/OH/H2O. Atomic O must not be duplicated
   in a second molecular reservoir.

4. **B1 requires more than adding the O I 1302/1304/1306 line entries.**
   The transit tool can read total O I density, but the triplet components start
   from different fine-structure levels of the O I ground term. The cooling code
   solves a three-level statistical equilibrium problem internally, but it does
   not export the level populations. Applying total O I to all three lines would
   count the same atoms repeatedly. B1 first needs a lower-level population
   policy, followed by a combined triplet spectrum and an observational treatment
   of the integration band, ISM absorption, and geocoronal contamination.

5. **Most of the Huang Table 4 audit proposed in B2 has already been done.**
   `charge_exchange.f90` states that Huang Table 4 contains only neutral and
   singly ionized charge transfer for C/N/O, with no O2+ reaction. The older
   Kingdon-Ferland O2+ + H term was removed intentionally. B2 should therefore ask
   whether an O2+ reaction outside Table 4 is scientifically necessary and which
   modern rate source should be used, rather than asking whether Table 4 contains
   such a reaction.

The safest order is to integrate and compare the external photochemistry first.
A2 should begin only if those results show that finite-rate oxygen chemistry
inside the EXHALE wind is required.

## 2. Chemistry structure in the current EXHALE code

### 2.1 Species layout and ownership of oxygen

`n_species` is fixed at 37 in `global_parameters`.

| Block | `f_sp` columns | Implementation |
|---|---:|---|
| H/He/He 2^3S | 1-6 | Fixed layout |
| C/O/N/Mg/Si/Ca/Na/K/S/Fe ion stages | 7-33 | Generic metal metadata in `species_table.f90` |
| H2, H2+, H3+, HeH+ | 34-37 | Used only when `thereis_mol` is true |

Oxygen is already represented as one metal element:

- `im_OI=4`, `im_OII=5`, and `im_OIII=6`
- The elemental abundance is `X_O = n_O/n_H`
- The current cell total is
  `nm_tot(:,O) = n(OI) + n(OII) + n(OIII)`
- The generic metal residual fixes total O, solves two unknown fractions for O+
  and O++, and assigns the remainder to O I

Adding an independent molecular oxygen reservoir while retaining the present O I
reservoir would duplicate total oxygen. The closure must instead be

```text
n_O,tot = n(O) + n(O+) + n(O++) + n(OH) + n(H2O) + ...
```

O I can remain the neutral stage in the atomic metal block, but its remainder
must exclude oxygen bound in OH and H2O.

### 2.2 Meaning of the current molecular solve

`System_HeH_mol.f90` has the following seven unknowns, with one more when He 2^3S
is enabled:

```text
H+, He+, He++, 2H2/H, 2H2+/H, 3H3+/H, HeH+/H, [He 2^3S]
```

When metals are enabled, `System_HeH_mol_metals.f90` appends two X+/X++ unknowns
for each metal. The current coupling between the molecular and metal blocks is
limited to:

- a shared electron density
- H/He/metal charge exchange
- H2 and He ion chemistry

There are no reactions that connect an oxygen-bearing molecule to atomic oxygen.

The solve is also a **local steady-state calculation in each cell**. Molecular
species do not have continuity or advection equations, and `post_process_adv.f90`
explicitly excludes H2/H2+/H3+/HeH+. Adding OH and H2O to this residual would not
turn it into a general transport photochemistry model. A2 needs a timescale gate:

```text
tau_chem << tau_adv  -> a local steady-state A2 may be justified
tau_chem ~ tau_adv   -> species transport or external photochemistry is needed
```

### 2.3 Radiation model

Current molecular photoprocesses use two radiation treatments:

- H2 photoionization is evaluated on the existing XUV energy grid, above its
  15.4 eV threshold.
- H2 Lyman-Werner photodissociation uses a separate band-integrated stellar LW
  flux with H2 self-shielding.

Important H2O photolysis wavelengths lie mostly in the FUV below the H I
ionization edge. The standard EXHALE input grid starts at 13.6 eV, and the
numerical SED is extended to lower energies only for He 2^3S or low-ionization
potential metals. A2 therefore needs more than new rate coefficients and cross
sections. It requires one of the following:

1. Extend the numerical SED into the H2O photolysis band and lower `e_low`
   automatically when water chemistry is active.
2. Add a dedicated FUV band or wavelength-resolved photolysis module, analogous
   in interface to the current LW treatment.

The design must also account for OH and H2O columns and overlapping absorption by
other FUV species. Water photolysis with self-shielding cannot be represented by
one extra term in the current radiation functions.

### 2.4 What `base.inp` actually passes to EXHALE

`read_base_inp` recognizes exactly six keys:

```text
T_base, r_base, HeH_base, Kzz_base, q_H2_base, p_base
```

The consumer path for `q_H2_base` is

```text
q_H2_base
  -> h2_mixing_ratio_base()
  -> h2_bound_fraction()
  -> comp_ntot_bc()
  -> base pressure/EOS normalization
```

The cold-start H2 profile follows a different path:

```text
set_IC.f90
  -> q_h2_equilibrium(local p, T0)
  -> f_sp(:,H2)
```

`q_H2_base` does affect the wind solution, but indirectly through the base
particle count. The documentation and any publication based on this interface
should state that distinction explicitly.

### 2.5 Output, restart, transit, and diffusion

- `write_output.f90` writes a fixed layout containing 27 metal ions and four
  molecular species.
- `load_IC.f90` maps schema labels back to `f_sp` indices and reconstructs rho
  using both molecular and metal masses.
- `EXHALE_transit.py` explicitly reads only Mg II, Ca II, and Na I from
  `Ion_species_adv.txt`. O I is already present in NumPy column 10 of the same
  file, but the transit calculation does not use it.
- `He_metal_diffusion` can diffuse atomic oxygen as it does other metal elements.
  Input validation, however, forbids all combinations of molecular chemistry and
  `He_diffusion`. B4 can therefore be tested directly only in atomic runs. Using
  diffusion with A2 would require an oxygen flux shared by atomic and molecular
  carriers.

## 3. Reassessment of options A0 to A4

### A0: retain the current treatment

This remains a scientifically possible option, but it should be renamed and
described more precisely.

- Current function: use a photochemical H2 mixing ratio to determine the base
  EOS and particle count
- Not currently done: pin H2 at the boundary, transport H2, or initialize the H2
  species profile from `q_H2_base`
- Small possible improvement: use `q_H2_base` in the base-cell initial seed while
  stating clearly whether it is only an initial guess or a boundary constraint

Recommended name: **A0: external photochemical EOS anchor**.

### A1: expand the handoff vector

The original plan suggests that writing O, OH, H2O, and CO values to `base.inp`
would make the lower boundary more self-consistent. The current parser ignores
such keys, and EXHALE has no physics that consumes those molecular abundances. A1
should be divided into three layers.

1. **A1a, provenance:** store `chemistry_code`, `chemistry_version`,
   `network_id/hash`, `profile_id`, and `p_base` as machine-readable values.
2. **A1b, elemental handoff:** connect `O_H_base`, `C_H_base`, `N_H_base`, and
   `S_H_base` to `melem_ab`. This is physically meaningful because EXHALE already
   transports and redistributes elemental reservoirs.
3. **A1c, molecular-species handoff:** treat OH/H2O/CO values as diagnostic
   metadata until the EOS, chemistry, opacity, initial condition, or a boundary
   condition actually consumes them.

Adding species keys without these consumers would widen the interface without
changing the model.

### A2: minimal O/OH/H2O catalytic network

The direction is reasonable, but the unit of implementation must be a
**single-owner C/H/O element-budget extension**, not simply two new species.

At minimum, the design must include:

- one total O budget containing O/O+/O++/OH/H2O
- subtraction of one H nucleus for OH and two for H2O from the H budget
- reactions that couple the water family back to the atomic O reservoir
- FUV flux and column attenuation for OH and H2O photolysis
- photochemical heating and reaction-energy bookkeeping
- molecular mass, particle-count, and electron contributions to the EOS
- root physicality checks and an element-simplex check
- IC, restart, and output schema support
- a Damkohler-number test of the local-equilibrium assumption

The two reactions in the original document may not close the connection between
the OH/H2O family and atomic oxygen. Before implementation, the reference network
should be reduced at the reaction level and should at least examine these
families:

```text
H2O + photon -> OH + H
OH + H2       -> H2O + H
O + H2        <-> OH + H
OH + photon   -> O + H        [include and branch only if supported by the audit]
```

This is a conservation skeleton, not a proposed final network. Rates, reverse
reactions, and photolysis branches should be selected only after a direct
comparison of the relevant VULCAN and Photochem hot Jupiter mechanisms.

### A3: reduced C/H/O photochemistry inside EXHALE

The main risk is not the amount of code but the change in model meaning. EXHALE
uses local algebraic steady state, whereas Photochem and VULCAN are kinetics
models with vertical transport. Expanding the residual to 15 to 25 species
without transport would produce a large local-equilibrium residual, not a compact
version of Photochem.

Since Photochem already supplies an official gas giant workflow, A3 should remain
low priority unless eliminating the external dependency is itself a clear science
requirement.

### A4: two-way offline iteration

This is possible in principle, but the feedback variables must be defined first.

- Photochem/VULCAN to EXHALE: temperature, radius, elemental abundances, and
  species composition
- EXHALE to the lower-atmosphere model: upper-boundary species fluxes, total
  escape flux, irradiation or attenuation, and possibly temperature or pressure

The current handoff implements only an EOS anchor in the first direction. It has
neither a reverse species-flux schema nor a defined overlap region in which both
models are valid. The interface physics must therefore be specified before an A4
driver is written.

## 4. Missing point: Photochem chemistry for hot Jupiters

The omission noted in the original plan is real. The Photochem source in this
workspace provides three relevant layers.

1. `photochem/utils/_format.py::zahnle_rx_and_thermo_files`
   - Uses H/He/N/O/C/S as the default element set.
   - Selects species and reactions from `zahnle_earth.yaml` to generate gas giant
     reaction and thermodynamic files.
2. `photochem/extensions/gasgiants.py::EvoAtmosphereGasGiant`
   - Provides equilibrium initialization and a quench-based initial guess.
   - Supports pressure-based temperature and Kzz profiles.
   - Applies gas giant boundary policies and manages the top-of-atmosphere
     pressure.
3. `photochem/examples/GasGiants.ipynb`
   - Explicitly identifies WASP-39b as a hot Jupiter.
   - Generates an H/He/N/O/C/S mechanism.
   - Calculates SO2 photochemistry for WASP-39b.

The existing EXHALE comparison used

```python
zahnle_rx_and_thermo_files(
    atoms_names=['H', 'He', 'N', 'O', 'C'],
    remove_reaction_particles=True,
)
```

This was not an accidental use of an unrelated rocky-planet network. It follows
the same generation path as Photochem's official gas giant tutorial, with sulfur
omitted. Three qualifications still matter:

- The historical scope of the source mechanism is not the same as the scope of
  the hot Jupiter benchmarks.
- An official example does not by itself validate H/H2 at 1 microbar for
  HD 189733b.
- The previously reported factor of 7.2 mixes network differences with upper
  boundary and domain differences. A reaction-level audit and a matched-domain
  rerun are still needed.

The Photochem repository exists at `./photochem`, but this workspace currently
lacks the compiled extension and `photochem_clima_data`. As a result,
`from photochem.extensions import gasgiants` fails. The source supports the
workflow, but the environment is not yet ready to run it as an EXHALE production
pre-step.

## 5. Reassessment of B1 to B4

### B1: O I 1302/1304/1306 transit spectrum

B1 has high scientific value, but it is not only a line-list change. The required
sequence is:

1. Read total O I from `Ion_species_adv.txt`.
2. Determine the populations of the O I ground-term `3P2/3P1/3P0` levels.
   - A fast first implementation could use a documented LTE/Boltzmann partition.
   - The preferred implementation would reuse the three-level statistical
     equilibrium physics used by cooling, as a function of T, ne, nHI, and the
     radiation field, and export the resulting fractions.
3. Assign the correct lower-level density to each 1302/1304/1306 component.
4. Calculate the full triplet spectrum and the same integration band used by the
   observation.
5. Decide which parts of ISM absorption, geocoronal contamination or masking,
   and the intrinsic stellar profile belong in the forward model or the
   observation-comparison layer.

The generic `resonance_spectrum` routine can be reused for planetary absorption,
but it is not sufficient by itself to reproduce the published Vidal-Madjar band
depth.

### B2: oxygen charge exchange

The O-related Huang Table 4 pairs currently present are:

- A13/A14: O + H+ / O+ + H, active by default
- C5/C6: O + He+ / O+ + He, active with `cx_full`
- D17/D18: Fe + O+ / Fe+ + O, active with `cx_full`

The `cases 14, 33, 48` labels in the original table do not match the current
`select case` numbering. The present case numbers are 13/14 for A13/A14, 30/31
for C5/C6, and 48/49 for D17/D18.

The direction of the endothermic barrier in the O/H pair has already been audited
and corrected in `HUANG2023_TABLE4_OXYGEN_ERRATUM.md` and in the code comments.
The revised B2 should answer two questions:

1. Does a reaction outside Table 4, such as O2+ + H -> O+ + H+, materially change
   the O III profile at EXHALE temperatures and ion fractions?
2. If so, which modern source and reverse-rate or detailed-balance policy should
   be adopted?

### B3: cooling and radiative pumping

The original document is correct that O I and O II cooling are implemented.
However, the exact ground-term statistical equilibrium is an internal cooling
calculation, not an exported level population. A useful shared deliverable for B1
and B3 is therefore an **O I ground-term level-fraction helper and output field**.
If Ly-beta pumping is added, cooling and transit should use the same populations.

### B4: oxygen diffusion

An atomic-oxygen experiment is possible without code changes under these
conditions:

- `He_diffusion: True`
- `He_metal_diffusion: True`
- `Molecular chemistry: False`, because the parser currently forbids the
  combination
- Element diffusion is computed from the sum of the O ion stages

After A2, oxygen in OH and H2O would have to participate in the same element flux,
so the present B4 implementation could not be used unchanged.

## 6. Revised roadmap

### P0: correct the documentation and interface semantics

1. Describe `q_H2_base` as a base EOS anchor, not a composition pin.
2. Add code and network provenance keys to `base.inp` or its successor.
3. Document Photochem's official gas giant and hot Jupiter workflow.
4. Divide A1 into provenance, elemental abundances, and molecular species.

This phase removes ambiguity in result interpretation without adding chemistry.

### P1: evaluate Photochem as a peer external pre-step to VULCAN

1. Install the compiled Photochem package and data package in a reproducible
   environment.
2. Implement a separate `Lower atmosphere: photochem <R_1bar>` adapter.
3. Compare the following runs with matched T(p), Kzz, spectrum, gravity,
   abundances, and vertical domain:
   - VULCAN with the reference hot Jupiter network
   - Photochem with the same converted network
   - Photochem with its official gas giant H/He/N/O/C(/S) mechanism
4. Do not rely only on a case such as HD 189733b where EXHALE itself does not
   reach steady state. Measure the EXHALE-level effect on at least one converged
   planet.
5. Export the production and loss budgets of the H2O/OH catalytic reactions near
   1 microbar.

These results should define the A2 reaction list and its validation targets.

### P2: define an explicit handoff contract

Every handoff field should belong to one of the following categories:

| Semantics | Examples | Required EXHALE consumer |
|---|---|---|
| Provenance | code/network/version/hash | Setup report and output header |
| Hydro/EOS boundary | T, r, p, q_H2 | `composition` and base boundary condition |
| Elemental reservoir | He/H, O/H, C/H, N/H, S/H | `melem_ab`, IC, and restart |
| Species initial guess | H2, OH, H2O | `set_IC` |
| Species boundary constraint | H2, OH, H2O | Chemistry boundary or residual |

The current `q_H2_base` belongs only to the hydro/EOS row. An initial guess and a
boundary constraint need separate, explicit semantics.

### P3: complete the atomic-oxygen observable first

1. Redefine B2 as a sensitivity audit of omitted O2+ rates.
2. Add an O I ground-term population helper and output.
3. Implement the O I triplet planetary spectrum.
4. Compare HD 209458b with matching observational bands and explicit ISM and
   geocoronal treatment.
5. Run atomic-oxygen diffusion A/B tests.

B1 still offers a high return, but the rate audit and level-population prerequisite
should be closed before quoting agreement with the observations.

### P4: A2 go or no-go decision

Start in-code A2 only if all of the following conditions are met:

1. A reaction budget from an external hot Jupiter network shows that the H2O/OH
   cycle dominates H2 destruction near 1 microbar.
2. The chemical time in that layer is sufficiently shorter than the advective
   time.
3. Feedback between the EXHALE wind and the lower model creates a measurable
   science requirement that an EOS-only handoff cannot meet.
4. The project accepts the validation cost of refactoring O/H conservation and
   the FUV radiation treatment.

Otherwise, an explicit handoff from a Photochem or VULCAN external pre-step is
more accurate and easier to validate.

## 7. Code areas affected by A2

| Area | Main files | Required changes |
|---|---|---|
| Species metadata | `parameters.f90`, `species_table.f90` | OH/H2O indices and H/O nuclei, mass, and charge metadata |
| Nonlinear residual | `System_HeH_mol_metals.f90`, `ion_residual_core.f90` | One total-O closure and coupled H/O rows |
| Cell driver | `ionization_equilibrium.f90` | Include molecules in total O, plus seeding, extraction, and root validation |
| Rates | `mol_rates.f90` or a new O chemistry module | Thermal reactions with documented sources and validity ranges |
| Photolysis | `cross_sec.f90`, `sed_read.f90`, `util_ion_eq.f90`, or a new FUV module | OH/H2O cross sections, branching, columns, and attenuation |
| EOS | `composition.f90`, `utilities.f90` | OH/H2O mass and particle counts without double counting |
| Energy | `util_ion_eq.f90` and cooling/heating diagnostics | Photolysis heating and reaction energetics |
| IC/restart/output | `set_IC.f90`, `load_IC.f90`, `write_output.f90` | Schema labels and conservation-preserving reload |
| Diffusion | `binary_element_diffusion.f90` (this row named `species_diffusion.f90`, deleted 2026-08-25) | One element flux for atomic and molecular O if required |
| Transit | `EXHALE_transit.py` | Separate total and atomic O, with level populations |
| Tests and docs | Regression matrix, input schema, setup report | Zero-O limit, old-default identity, and element budgets |

A2 is therefore substantially larger than adding two or three unknowns. The
existing merged molecular and metal residual is still the right starting point.
Its element bookkeeping should be generalized rather than replaced by a new
solver.

## 8. Minimum validation gates

### External Photochem and VULCAN comparison

- Repeat the code comparison with the same network and vertical domain.
- Record the species count, reaction count, and mechanism hash for the official
  Photochem gas giant mechanism.
- Store H2O/OH production and loss budgets at 1 microbar, not only H2/H.
- Preserve all handoff provenance in the setup report and final output.

### A1 handoff

- Classify every key as provenance, EOS, element, initial guess, or boundary
  constraint.
- Mark any unconsumed species key explicitly as diagnostic.
- Reproduce the previous result with the legacy six-key `base.inp`.
- Check closure of the H, He, C, N, O, and S budgets after elemental handoff.

### A2 in-code chemistry

- Recover the existing molecular run at zero oxygen abundance.
- Recover the existing O I/O II/O III result in the atomic limit where OH and H2O
  vanish.
- Close the H and O element budgets within solver tolerance in every cell.
- Output the chemical-time to advection-time profile.
- Match absorbed energy and heating when the FUV opacity is expanded.
- Preserve H/O mass and species in an IC write/read round trip.
- If molecular post-processing remains molecule-free, either prohibit `_adv`
  output for these runs or emit an explicit warning.

### O I transit

- Verify that the three lower-level populations sum to total O I.
- Test both an optically thin column and a saturated line.
- Compare band-integrated depths over the observational band, not only line-center
  absorption.
- Record the ISM and geocoronal treatment in result metadata.

## 9. Recommendation

Replace the original short-term order `B1 -> B2 -> A1 -> A2` with:

```text
P0  Correct the documentation, q_H2 semantics, and provenance
 -> P1  Run matched comparisons including Photochem's official hot Jupiter path
 -> P2  Define the handoff contract and pass elemental abundances
 -> B2' Test omitted O2+ rates and add an O I level-population helper
 -> B1  Build the O I triplet observational forward model
 -> B4  Run atomic-O diffusion A/B tests
 -> A2  Implement only after the reaction-budget and Damkohler gates pass
```

The first question is not simply whether to use Photochem or VULCAN. The model
must record which hot Jupiter mechanism and vertical domain produced the lower
boundary state, and which part of that state EXHALE actually consumes. The current
plan mixes these issues, which makes both the role of `q_H2_base` and the need for
A2 appear stronger than the code supports.
