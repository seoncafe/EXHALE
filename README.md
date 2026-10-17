# EXHALE v1.10

**EXoplanetary Hydrodynamic Atmospheric Loss and Escape** (EXHALE) is a 1-D
radiation-hydrodynamics code that simulates photoionization-driven atmospheric
mass loss from irradiated exoplanets. It solves the spherically symmetric Euler
equations together with the ionization and thermal balance of a hydrogen/helium
atmosphere carrying trace metals and, optionally, molecules, and returns the
mass-loss rate, the wind structure, and the transit transmission spectra that
follow from it. EXHALE started from the ATES code (Caldiroli et al. 2021;
Biassoni et al. 2024), of which it keeps the outline of the time-marching
hydrodynamic loop; almost everything else was written for
EXHALE, and the code is essentially a new one: the metal chemistry and
cooling, the molecular and lower-atmosphere layers, the diffusive separation
of the elements, the Ly-alpha radiative transfer, the well-balanced
reconstruction and the base boundary, the certification of a steady state,
the Newton-Krylov steady-state solver with its partitioned composition
iteration, the heat conduction, and the transit synthesis.
Every extension that adds a species, a transport process or a field is
opt-in; corrections to the physics of the ATES problem itself (the He 2^3S
triplet, secondary ionization, case A/B recombination, the He/H charge
exchange) are on by default, so a bare `input.inp` solves the atomic H/He
problem with that physics rather than reproducing the legacy ATES model bit
for bit.

---

This README starts from the v1.00 README and adds the v1.10 implementation
changes below. The detailed physics description and bibliography are retained
from v1.00. Numerical benchmarks quoted in that inherited description are
historical v1.00 measurements, not renewed validation of v1.10. The v1.10
implementation and limitations in this section take precedence over historical
recipes. References to documents available only in the neighboring v1.00 tree
are labeled accordingly; that tree is needed to open those local links.

## Changes in v1.10

### Regional stabilization and stationary convergence

Two distinct outcomes are used:

| Outcome | Required evidence |
| --- | --- |
| **Regional stabilization** | Upper-atmosphere profiles remain stable and lower-atmosphere changes stay within the configured bounds throughout a finite confirmation interval |
| **Stationary convergence** | The enabled fluid, chemistry, radiation, and boundary equations satisfy their stated stationary residual and consistency tolerances |

Regional stabilization is implemented in
[`regional_convergence.f90`](src/modules/time_step/regional_convergence.f90)
and invoked by the accepted-step path in
[`EXHALE_main.f90`](src/EXHALE_main.f90). It does not require JFNK and does not
start a stationary solve when its criteria are met. A small change over time
alone does not establish stationary convergence.

The monitor is disabled unless `Regional convergence: True` is present.
Its supported configuration is atomic, global physical-time marching:
`Run mode: phys`, with local time stepping, molecular chemistry, Newton solving,
post-processing-only operation, and stationary restart intent disabled.
Unsupported combinations stop at startup. Use a complete physical input file;
the two files below are **settings fragments**, not standalone planet inputs:

- [Huang regional settings](examples/regional_convergence/Huang.inp): lower
  boundary between regions at 1.3 R_p, upper at 1.85 R_p; window
  87965.04804471003 s.
- [Salz regional settings](examples/regional_convergence/Salz.inp): both region
  boundaries at 1.2 R_p, with no transition region; window 57600 s.

Merge the appropriate fragment into `input.inp`, replacing duplicate keys.
The proposed default thresholds in the current source are:

| Quantity | Upper region | Lower and transition regions |
| --- | ---: | ---: |
| Window-mean change and span | 0.01 | 0.10 |
| Full range | 0.05 | 0.20 |
| Regional mass and internal-energy change and range | 0.01 | 0.05 |
| Mean actual face mass-flow change and span | 0.01 | 0.05 |

All configured criteria must pass for eight consecutive windows
(`Regional validation windows: 4` times `Regional confirmation factor: 2`),
with at least 40 accepted endpoints in each window. The monitor samples every
accepted step; rejected attempts do not advance its physical confirmation time.
It uses the numerical face mass flux combined with the Runge-Kutta weights,
and regional internal energy excludes kinetic energy. Velocity and species
normalization uses sound speed and hydrogen-nucleus density where appropriate;
these limits are not relative tolerances on every rare species or on a nearly
zero velocity.

All limits, region boundaries, durations, and sampling requirements are input
settings. The source defaults place both region boundaries at 1.2 R_p.
`Regional window seconds: 0` selects twice the initial outward acoustic travel
time. `Regional maximum seconds: 0` selects the settling duration plus two
complete confirmation intervals; it does not mean unlimited execution.
`Regional settling seconds` defaults to zero. Explicit regional inventory and
mass-flow limits override the common `Regional inventory limit` and
`Regional mass flow limit` values independently of input-line order. Unknown
`Regional ...` keys are rejected.

These are proposed finite-duration criteria, not universal physical accuracy
bounds or a proof that oscillations remain bounded indefinitely. Neither Huang
nor Salz is accepted by this README update. See
[the definitions](md/CONVERGENCE_LEVELS.md),
[the default-setting rationale](md/PROPOSAL_20261016_upper_focused_regional_defaults.md),
and [the completed-results assessment](md/REVIEW_20261016_completed_results_convergence_conditions.md).

### Progress output and plots

The regional monitor writes:

| File | Contents |
| --- | --- |
| `output/regional_convergence.csv` | Regional criterion history |
| `output/regional_convergence_status.txt` | Current status, effective settings, and achieved convergence level |
| `output/regional_profile_history.txt` | Accepted profile records with region identifiers |

Profiles/status are recorded at `Regional profile interval` (40 accepted steps
by default), with initial/final profile records. The stopping statistics still
use every accepted step. On successful completion, the status includes:

```text
target_convergence_level=regional_stabilization
convergence_level=regional_stabilization
stationary_root_claim=False
status=regional_limits_satisfied
```

An incomplete regional run has `convergence_level=none`; the main program's
normal finalization exits with status 4 if regional monitoring was enabled but
its criteria were not met. Other startup/runtime failures can terminate earlier.
Read the status and numerical records as well as the exit code.

From the calculation directory, visualize existing output with:

```bash
python3 /path/to/EXHALE_v1.10/EXHALE_plots.py --regional --live 60
python3 /path/to/EXHALE_v1.10/EXHALE_plots.py --regional --save output/regional
```

The first command refreshes figures every 60 seconds; the second writes PDF
figures. The regional plots separate lower, transition, and upper profiles and
show the evolution of basic physical quantities and convergence criteria.
Plotting reads existing outputs and does not advance the solver. `--profiles`
selects the ordinary profile view when automatic regional detection is active.

### Stationary boundary chemistry

[`stationary_boundary_chemistry.f90`](src/modules/states/stationary_boundary_chemistry.f90)
and [`reservoir_composition.f90`](src/modules/states/reservoir_composition.f90)
provide the lower-boundary thermochemical calculation used by the stationary
residual and stored-state certification. `steady_newton.f90` invokes
`stationary_ghost_state`; `certification.f90` queries its boundary result.
An unclosed boundary cannot establish a certified stationary solution.
This implementation does not establish stationary solutions for all historical
Yan, Huang, or Salz states.

The boundary design and validation record are in the
[design](md/DESIGN_20261012_stationary_boundary_closure.md),
[integration report](md/IMPLEMENTATION_20261012_stationary_boundary_chemistry.md),
and [boundary analysis](md/stationary_boundary_thermochemistry_efficient_analysis.md).
Optional fixed-state stellar-rate reuse is implemented in `util_ion_eq.f90`;
it is a computational option with restricted reuse conditions, not new physics.

### Build, controls, and validation differences

- Build the current source explicitly. Existing executables and workflow
  snapshots may represent earlier source revisions; their presence does not
  show that the latest source was compiled.
- `make check` requires `backup/regression/run_check.sh`. That harness is absent
  from this tree, and the target returns failure when it is unavailable.
  `make test` runs the available scientific suites and can execute model cases;
  it is not a documentation-only check.
- `Max steps` and `EXHALE_MAXSTEPS` currently set separate limits. Likewise,
  `EXHALE_MARCH_FLUX_TH` controls marching stabilization, whereas `Flux spread
  tol` controls stationary acceptance. Do not treat them as interchangeable.
- `EXHALE_MONITOR_STEPS` must be positive. The source rejects zero and negative
  values before the monitor-window calculation.
- The material created for error testing was removed from the source on
  2026-10-16: the fault-injection environment controls and test hooks, the
  expected-failure and invalid-input test rows, and their drivers and
  fixtures. Real failure detection, retry, rollback and the rejection of
  invalid inputs remain. The proposed input migrations are not implemented.
  The [options and cleanup review](md/environment_variable_options_survey_20261016.md)
  lists the removed names and the decisions; check the actual reader for the
  interface supported by a particular executable. Workflow snapshots and
  executables built before that date still contain the removed controls.

The development tree originated from the October 12, 2026 v1.00 copy, recorded
in `provenance/source_copy.json`. Historical source, input, and executable
records are retained under `provenance/`; calculations use `runs/` and workflow
records use `workflows/`. These local development artifacts and the linked
`md/` reports may not be distributed with a source checkout. The
[handoff](md/HANDOFF_20261015_marching_goal_active.md) records the user-directed
calculation hold; the word "active" in its historical filename is not an
instruction to resume calculations.


---

## Features

**Chemistry and radiation**

- Trace metals (C, N, O, Mg, Si, Ca, Na, K, S, Fe) solved self-consistently
  inside the coupled ionization system, with Badnell RR+DR recombination
  (Huang et al. 2023 fits for Fe I/II and Shull & Van Steenberg for Ca I, whose
  isoelectronic sequences the Badnell project does not reach), Voronov
  collisional ionization, and the charge-exchange reactions of Huang et al.
  (2023) Table 4, extended by the `O2+ + H0 -> O+ + H+` and
  `N2+ + H0 -> N+ + H+` electron captures of Barragan et al. (2006) that
  Table 4 omits (on by default; `metals.inp: cx_O2p_H 0` or `cx_N2p_H 0`
  leaves either out); each rate is taken from its original source and each
  pair held to detailed balance with internal partition functions (the
  source's direction kept, the other derived, and for the three forwards whose
  sources name their product levels the reverse follows the product
  populations of a five-level statistical equilibrium at the local electron
  density), the energy defects heat the gas (`heat_charge_exchange`), and the
  recombination and collisional-ionization energies of the metals enter the
  cooling
- Metal-line cooling from CHIANTI v11 (collision strengths descaled with the
  natural spline CHIANTI prescribes, observed level energies): closed-form
  fits for Mg II and Na I, (T, n_e) statistical-equilibrium tables for Ca II,
  Mg I, Fe II and the levels of C/N/O above their ground terms, and the Fe I
  table in its LTE limit, with
  the split ground terms of C I, C II, N II and O I solved in exact statistical
  equilibrium at the local `(n_e, n_HI)` instead of the coronal limit, which
  gives the density-dependent saturation of the [C II] 158 um and [O I] 63 um
  floors
- Line trapping in the eight ground-term fine-structure lines of C I, C II,
  N II and O I (the line-center escape probability from the column above each
  cell enters the statistical-equilibrium solution as `A_ul -> beta*A_ul`),
  and a smooth cutoff of the coronal fits below their
  10^3 K validity floor, leaving only the explicit fine-structure terms
- H and He recombination with case A where the gas is thin and case B where
  it is thick (`H_rec_escape`, `He_rec_coupling`; ground captures from the
  Milne relation on the transfer's cross sections), the He II and He III
  recombination radiation absorbed on the spot, the He I 584 A resonance line
  scattered by He I and partly converted into the 2^1S two-photon continuum
  (Wood, Mathis & Ercolano 2004), and every recombination cooling term
  `k T alpha (3/2 + dln alpha/dln T)` on the coefficient the balance uses
- He I 2^3S metastable triplet in the coupled solver, with the
  electron-impact 2^3S -> 1^1S de-excitation (the detailed-balance reverse of
  the excitation), the collisional feed through the higher triplets, the
  collision strengths of Bray et al. (2000), temperature-dependent He(2^3S)+H
  and +H2 destruction rates
  (Garcia Munoz 2025, continuous closed forms with a 0.9/0.1
  Penning/associative branching)
- Updated photoionization data: He I ground state from Verner et al. (1996),
  a He I 2^3S cross section extended past 60 eV against TOPbase, and the
  inner shells of the metal ions (Verner & Yakovlev 1995) with their Auger
  and fluorescence relaxation (Kaastra & Mewe 1993): the Auger electrons heat
  the gas and a neutral that loses two or more electrons jumps two stages
- He+ + H charge exchange in its radiative channel (Stancil, Lepp & Dalgarno
  1998, row 19) and its non-radiative one (the H+ + He cross section of
  Loreau, Ryabchenko & Vaeck 2014 by reciprocity, taken as
  `sigma_r(15) sqrt((x - x_th)/(15 - x_th))` between the 13.655 eV/u
  threshold and the 15 eV/u where their calculation starts), its
  detailed-balance reverse and the radiative He2+ + H capture in every
  helium-bearing system (`He_H_charge_exchange`, default on)
- One stellar spectrum type builds every band of the photon grid, the XUV and
  the part below 13.6 eV where the He 2^3S metastable and the low-IP metals
  absorb: `Spectrum type: Power-law` is the power law everywhere,
  `Spectrum type: Planck` is the photospheric blackbody
  `pi B_nu(T_eff) (R_star/a)^2` built from `Stellar Teff` and `Stellar radius`,
  and `Spectrum type: Load` is the table everywhere. A loaded table that stops
  above the lowest threshold of an active absorber (4.768 eV = 2600 A with
  `Include He23S? True`, or a neutral metal's threshold) stops the run rather
  than leaving that absorber without a field; the startup report states the
  type, the source of the band below 13.6 eV, and the integrated grid flux
- Secondary ionization by fast photoelectrons: the ionization budget of
  Shull & van Steenberg (1985), with the heating fraction, its dependence on
  the photoelectron energy and the H2 channels from Dalgarno, Yan & Liu (1999),
  above 30 eV (`Secondary_ionization: True`, the default, applies the
  coupling from the first step; `Staged` keeps the earlier default, which
  switched it on only once the marching had settled without it; `False`
  turns it off). Below 30 eV the photoelectrons are deposited whole as heat by
  default; the opt-in key `Low-energy electron partition: Furlanetto2010`
  partitions those of 10.2-30 eV into heat, H I and He I ionizations and
  escaping H I line radiation by the Monte Carlo tables of Furlanetto &
  Stoever (2010), in cells without H2 (a cell with any n(H2) > 0 keeps the
  heat-only treatment: the tables are those of an atomic H/He gas;
  `src/modules/radiation/low_energy_electron_degradation.f90`). The
  photoelectrons of the recombination radiation absorbed on the spot pass,
  at their representative energy, through the same bands and partition as
  the stellar ones
- Non-LTE H(n=2) populations and Ly-alpha radiative transfer, from either a
  fast Neufeld escape-probability closure or a field imported from the LaRT
  Monte Carlo code, including the Ly-alpha emitted in situ within the wind
- Molecular chemistry: H2, H2+, H3+ and HeH+ in the coupled ionization
  equilibrium, with H2 photoionization opacity/heating, H3+ infrared cooling
  (the LTE cooling function of the MiZATeP line list, Mizus et al. 2017,
  times the non-LTE departure factor of Miller et al. 2013), and H2
  photodissociation in the Lyman-Werner bands
  from a line-by-line calculation over every Lyman and Werner transition,
  tabulated on our own (T, n_H, N_H2) grid -- the dissociation cross section
  itself, with the line self-shielding and the line overlap inside it and the
  trapping of the fluorescent decay photons taken from level-resolved CLOUDY
  runs, rather than a published closed-form fit
- Oxygen chemistry: OH, H2O and CO in the same coupled system, with the FUV
  photolysis of H2O and OH in four bands, so the base H2/H partition is
  computed rather than imported. Rates from Baulch et al. (2005) and the IUPAC
  evaluations, reverse rates by detailed balance against a NIST-JANAF Shomate
  table, CO carried as an oxygen reservoir. The first band is the 912-1201 A
  Lyman-Werner interval, where H2, H2O and OH share one beam: the H2 lines and
  the H2O/OH continuum each attenuate what the other sees, so the interval has
  one incident flux (`Stellar LW flux`). With oxygen chemistry, molecular
  carrier transport is enabled by default (`Molecular carrier transport`);
  without oxygen chemistry it is disabled by default. The transported
  carriers are H2, H2+, H3+, HeH+, and, with oxygen chemistry, OH, H2O and
  CO. They use an
  implicit diffusion-advection solve coupled to their chemistry. A stationary
  solve can take those balances as its own unknowns (`Coupled carrier solve`,
  three-valued: `False` by default, `True`, or `On stall`, which runs the
  alternation and hands the state to the coupled block at the pass where the
  alternation gives up; a word the key has no meaning for stops the run); on a
  molecular configuration that needs
  `Molecular carrier transport: True` and is refused at startup without it,
  because an eliminated H2 content is one the local equilibrium cannot
  determine and the stationary residual is then not a function of its
  unknowns. The
  ionization stages can ride the same operator (`Ionization transport:
  True`, default off), in an atomic gas as well as a molecular one: x(H II)
  per hydrogen nucleus and x(He II), x(He III) per helium nucleus are carried
  on their own element's nucleus face flux with the local ionization balance
  as their source, the neutral stage closing each simplex, and every
  equilibrium sweep is handed the transported fractions instead of putting
  them back on their local root. It needs helium in the mixture, in a
  molecular gas also `Molecular carrier transport: True`. An explicit
  `Coupled carrier solve: True` is accepted in a molecular gas, but refused
  in an atomic gas; `On stall` is refused with ionization transport. The local root is wrong where
  a parcel leaves its shell faster than it ionizes (`P r/|v| < 1`), which is
  what the Koskinen et al. (2022) comparison of
  `../EXHALE_v1.00/md/k22_electron_density_excess.md` needed
- He 2^3S transport (`He 2^3S transport: True`, default off; needs
  `Include He23S: True` and `Ionization transport: True`): the metastable
  population n(2^3S)/n(He nuclei) is carried as a fourth helium state of the
  same operator, on the helium nucleus flux, with the level balance of the
  coupled system as its source (the steady advection balance of Oklopcic &
  Hirata 2018 and Lampon et al. 2020, with the eddy term of the stage flux),
  instead of the local root of that balance in every cell; its certification
  entry is `carrier balance He2^3S`
- Thermal infrared field of the atmosphere below the base, so the molecular
  and fine-structure coolants return the net rate rather than the vacuum limit
- Molecular infrared bands (`Molecular IR bands`): the H2 quadrupole and
  magnetic dipole line spectrum (Roueff et al. 2019) and the H2O and CO
  vibration-rotation bands (HITEMP), in LTE and exchanging with that same
  field, so the layer below the H2 -> H front settles on a radiative
  equilibrium temperature instead of radiating itself away
- Composition-dependent caloric equation of state: below the H2 -> H front the
  gas stores its energy in the H2 rotational and vibrational ladder as well as
  in translation, so the internal energy is built from the observed 302-level
  H2 ladder (Roueff et al. 2019) rather than from a constant gamma = 5/3. An
  atomic run is byte-identical
- Chemical heat of the molecular network (`Molecular reaction heat`, default
  on; the key exists to turn it off): in a molecular gas the ionization energy
  a photon spends comes back to the GAS through the dissociative recombination
  of H3+ and H2+, not out of it as the Lyman photon of a radiative
  recombination, and the code now deposits it -- 81 percent of the total
  heating rate at the base of a converged hot Uranus. Each channel deposits
  what its measured product state leaves to translation and no more: the
  dissociative recombinations of H2+ and HeH+ put one H atom in n = 2 and
  return 0.749 and 1.554 eV to the gas, the H3+ recombination leaves 2.4795 eV
  as internal energy of the H2 fragment, and the three-body association's
  4.478 eV branches through the vibrational quench fraction, which is formed
  from the all-level radiative rate of the 302-level ladder against published
  collisional rates for each of the three colliders H, H2 and He. The
  association itself uses the collider-resolved rate coefficient, not the
  M = H2 coefficient on the total density. `EXHALE_REACTION_HEAT_RECIPIENTS=0`
  restores the older prompt-deposit behavior for comparison
- Diffusive separation of helium and metals: hydrogen and helium are
  transported as a two-component mixture, with bulk advection, binary
  diffusive settling in the computed ambipolar field, and an optional eddy
  term; each trace metal can diffuse independently. The friction is resolved
  by ionization stage -- hard sphere, polarization and Coulomb -- so an ion
  is held to the protons instead of settling at a neutral rate. The energy
  equation carries the enthalpy the moving elements take with them, the
  interdiffusion enthalpy flux `q_d = sum_s h_s J_s` (Cook 2009, eqs. 11-13),
  in the marching update and the stationary energy row alike
  and the enthalpy the transported molecular carriers (H2, OH, H2O, CO) carry
  relative to the rest of their elements, and, with `Ionization transport:
  True`, the enthalpy the carried ionization stages move relative to their
  own element
  (`Interdiffusion enthalpy flux:`, on by default whenever `He_diffusion` is,
  a molecular carrier is transported or the stages are carried;
  `False` only to reproduce published models that omit it; 13-26 percent of
  the local energy budget at 1.1-5 R_p on the LHS 1140 b diffusion states)

**Hydrodynamics and solvers**

- Finite-volume marching with PLM or WENO3 reconstruction and a stated
  interface flux (`Numerical flux:` is mandatory; `HLLC`, `ROE` or `LLF`),
  usable as a two-stage PLM -> WENO3 sequence. `Low Mach velocity jump:`
  (default off) scales the Roe branch's velocity jump by `min(Ma, 1)` after
  Rieper (2011); it is kept for study and is NOT recommended on a
  near-hydrostatic column, where it was measured to remove the damping of the
  odd-even velocity mode at exactly the faces that mode lives on and to prevent
  a solve unmodified Roe completes
- Jacobian-free Newton-Krylov steady-state solver with PTC warm-up, SER ramp,
  a non-monotone (Grippo) line search and a scaled trust region. The
  stationary system carries a row and an unknown for every transported
  balance the configuration activates, not only the hydrodynamic triple
- With regional monitoring disabled, the ordinary marching criterion stops
  when the NUMERICAL face mass flux `r^2 F_rho` is flat
  over `r >= r_flux` to a spread below 4e-4 (`march_flux_spread_th`) AND the
  state has stopped changing: over one monitor window of 500 steps no cell's
  rho or total energy moves by 1e-4 relative, nor its momentum by 1e-4 of
  `sqrt(rho E)` (`march_change_th`). A `[stop]` line prints both numbers
  at every monitor step. The cell-centered spread `du` of `rho v r^2` is
  printed as a diagnostic only (it has a resolution floor next to the outer
  ghost) and still drives the PLM -> WENO3 switch of `du_th [PLM,WENO3]:`.
  `Solver: Newton [switch]` hands the state to the JFNK solve once that face
  spread falls below `switch` (default 0.1; it was a `du` threshold before
  2026-10-10). `EXHALE_MARCH_FLUX_TH`, `EXHALE_MARCH_CHANGE_TH` and
  `EXHALE_MONITOR_STEPS` override the two thresholds and the window for
  measurement
- The Mach number of the outermost cell is printed with every Mdot, with a
  note below 1.1, "NOT A QUANTITATIVE MDOT" (also in `EXHALE_setup.out`): at
  a subsonic outer face one acoustic characteristic enters through the outer
  ghost condition, which can then set Mdot, and no other subsonic-face
  treatment exists (decision 2026-10-12). For a
  quantitative Mdot put the outer face past the sonic point with
  `Outer radius [R_p]:` (on `wasp_full_newton`, measured before update log
  stage 3 section 138 changed the columns, the default Roche L1 face,
  1.574 R_p at Mach 0.95, moves rho and v by up to 15-17% near 1.31 R_p and
  log10 Mdot by 0.01 against an outer face at 2 R_p;
  `../EXHALE_v1.00/md/outer_boundary_subsonic_review_20261009.md` section 5)
- `Composition update holds: pressure` (default `energy`): where the
  stationary solve alternates a hydrodynamic solve with the fixed-wind
  element and carrier relaxations, the composition update keeps the density,
  velocity and pressure of every cell and rebuilds the total energy from the
  caloric equation of state of the new composition, instead of keeping the
  conserved energy and letting the pressure move. At a base many decades
  subsonic a pressure step of 1e-4 against the reservoir drives a face mass
  flux decades above the wind's into the next solve; holding the pressure
  removes that step. A route through the iteration, not an equation set: at
  a fixed point both values leave the state unchanged, and the energy each
  update adds is printed at every pass. On the LHS 1140 b molecular states it
  removes the base step but does not by itself make the next solve converge
  (`../EXHALE_v1.00/md/lhs1140b_composition_update_pressure_20260924.md`)
- `Well balanced: True` (default False) carries the departure from each
  cell's own local hydrostatic equilibrium through the reconstruction, the
  Riemann jumps and the pressure force, so a discrete hydrostatic
  equilibrium is preserved to rounding (Kappeli and Mishra 2014, 2016). A
  measurement option: turning it on moves every result, and the
  species-row stationary solves do not improve under it
- Stationary certification reports the active balance checks: one evaluator states
  one condition per active balance and the run prints, and the state file
  records, which conditions were met and which row and cell refused. The
  tolerances are anchored by measurement, not inherited from the solver's own
  target: mass 3e-12, momentum 1e-8, energy 1e-6, the closure and the level
  balances at the tolerance the cell sweep already judges a root by, and the
  element and carrier rows at 1e-5 for `r >= 1.20 R_p`, the rows below that
  reported and not gating, because the same balance is a cancellation of
  advective terms in the wind and of eddy and settling terms in the layer and
  is not resolvable to one number in both. The certification is a
  diagnostic of the specified stationary criteria. Regional stabilization is
  a separate finite-duration outcome, as defined above; neither classification
  alone supplies an observational error budget. A state that is not certified
  carries its failing rows. The legacy non-regional exit statuses include 0
  for a certified state or a run that made no stationary claim (the marching
  stop, whose state is a relaxation snapshot, a step cap, a stall, a NaN); 2
  for a state the run declared stationary (the residual gate met while
  marching, a JFNK finish with info = 0) that is not certified (a flag, the state is written in full; a gfortran runtime
  error also exits 2, and is told apart by the certification's closing line
  `(certification) the state was written in full; the run exits with status 2
  because it is not certified.`); 3 when the stationary solve of a steady
  route failed (the solver's last iterate is written with
  `cert_reason=stationary_solve_failed`)
- Molecular transport: the Navier-Stokes viscous force, its dissipation, and
  heat conduction, integrated Crank-Nicolson and entering the steady residual
  with the same operator. `Conduction: True` (default off) uses the
  conductivity of the mixture, `kappa = (n_e 1.2e-6 T^2.5 + n_HI 379 T^0.69 +
  n_HeI 299 T^0.69 + n_H2 k_H2(T))/n` erg cm^-1 s^-1 K^-1 (Salz et al. 2015;
  Sutton et al. 2015; an Eucken-form fit to the tabulated H2 conductivity,
  valid 200-2000 K), with the base face held at the bath temperature of the
  lower atmosphere (the profile's own temperature when a lower-atmosphere
  profile is handed over, the reservoir's isentrope otherwise). The heat the
  wind conducts into the lower atmosphere through that face is measured on
  the certified state and written to `EXHALE_resolved.out`
  (`base_conductive_flux_cgs`, positive into the lower atmosphere, with
  `base_conduction_T_bath_K` and `base_conduction_T_cell1_K`).
  `EXHALE_CONDUCTION_SCALE=s` (0 <= s <= 1, unset = 1) multiplies the
  conductivity, a continuation factor for reaching a state with conduction
  from one solved without it, carried in the state's `cond` token
- Low-Mach damping, a gated fourth-difference stress for a shell that has
  stopped flowing, and an optional Shapiro filter for a breathing base
- Roche-potential geometry, in either a spherical or a Roche-lobe domain;
  `Outer radius [R_p]:` carries a Roche domain past L1 with the full tidal
  potential (as Huang et al. 2023 and Yan et al. 2022 do). Past the ridge the
  cold IC holds the density at its ridge value, the ghosts of a supersonic
  outermost cell are the steady isothermal wind through it, and the absorbing
  columns start at the outer face `r_edg(N)` in every domain (the outer
  ghosts carry no column)
- Dual-energy method for cold, highly supersonic expanding atomic gas
  (`rho e < 1e-2 E`, M >~ 13): the internal energy from the advected adiabat
  `p/rho^gamma` instead of `E - rho v^2/2` (Ryu et al. 1993; Bryan et al.
  1995); the run log prints the reset cell-steps
- Choice of initial condition: cold hydrostatic, warm Parker seed, an
  automatic selector, or a warm start from an included Wind-AE solver
  (Murray-Clay et al. 2009 / Broome et al. 2025)
- Lower-atmosphere connection: an analytic Koskinen et al. (2022) column, an
  EOS-only molecular-base correction, and a `base.inp` handoff generated
  either analytically or from a VULCAN photochemistry run that EXHALE launches
  itself. A handoff that states its own level (`p_base`) fixes the base level
  of the run, `n0 = p_base/(k_B T0 ntot_bc)`, so `Log10 lower boundary number
  density` is then unnecessary, and a pair that disagrees by more than 1% is
  refused at startup rather than one of the two silently winning
- OpenMP parallelization of the cell sweeps. The historical v1.00
  thread-comparison result was bitwise identical: the suite `src/tests/thread_count_identity` runs 200 steps of the
  metals-on regression case at 1 and 8 threads and compares every data line
  (it found, and the tree fixed on 2026-09-26, a shared scratch variable of
  the metal sweep that made metals-on runs between 2026-09-25 and that fix
  depend on the thread count)

**Tooling**

- Runtime configuration by file presence, not compile flags: `metals.inp`
  present means metals on, likewise `opacity.inp` and `base.inp`
- Runtime grid size: `Grid cells: <N>` in `input.inp` (default 500), so
  base-refinement studies run without a rebuild
- Shells beyond the constructed grid, off by default: `Outer shells
  [r_face,cells]: <r_face> <cells>` appends `<cells>` cells outside the
  outer face of the grid, the last face landing on `r_face` [R_p], with the
  cells, faces and volumes below the old outer face unchanged bit for bit, so
  a domain-size test moves only the outer boundary (`../EXHALE_v1.00/md/input_schema.md`,
  K33d)
- `EXHALE_transit.py` transmission post-processor: He I 10830 A, Ly-alpha,
  H-alpha, H-beta, the metal resonance doublets Mg II h&k, Ca II H&K and
  Na I D, and the O I 1302/1304/1306 triplet out of its three resolved
  ground-term fine-structure levels, with impact-parameter Voigt integration,
  instrument and rotation convolution, and an optional triaxial Roche geometry.
  Every `tpm_*.txt` carries a `# transit_schema 1` metadata block above its
  column line: which `_adv` profile it was built from, that profile's validity
  schema, certification verdict and provenance, the census of rows the line
  took its depth from that were refused a correction, and the
  `EXHALE_TRANSIT_*` overrides in effect, so a curve says what it stands on
- A model line compared with an observed transit profile is fitted to the
  observed equivalent width with its width matched by an added Gaussian
  velocity kernel, for every planet: the excess absorption of the
  spectrum convolved with the planet's rotation and the instrument (column 4,
  `T_rot+instr`) is convolved with a Gaussian of FWHM `f` in
  velocity, `f` chosen for each model so that the extractor used on the data
  returns the observed FWHM, and the EW is the window EW of that curve on the
  observation's samples (trapezoid weights). The kernel stands for a velocity
  field the 1-D model does not contain, not a computed physics term, and every
  result that uses it reports `f`. For LHS 1140 b He 10830 the code is
  `../EXHALE_v1.00/LHS1140b/he10830_equivalent_width.py` (`broaden_excess`,
  `width_matching_fwhm_kms`, `width_matched_product`, CLI `--width-matched`)
  and `../EXHALE_v1.00/LHS1140b/models/.P1/fiducial_conduction/ew_of_state.sh`, whose
  `EW_NONTHERMAL_KERNEL=width_match` adds `EW_kernel=` and `fwhm_kms=` to
  `EW.txt` (`EW=` stays the unbroadened value); the script evaluates the state
  with `Restart intent: stationary evaluate`, stages its
  `conserved_state.txt` and records the written state's certificate as
  `written_state`
- Lower atmosphere from Photochem with a two-way coupling:
  `src/utils/photochem_to_lower_profile.py` writes the lower-atmosphere profile
  and, with `--conducted-heat-flux F` (default 0, profile unchanged), solves
  the energy balance of the column below the matching level under the heat
  the wind conducts into it (`src/utils/conducted_heat_column_balance.py`:
  mixture conductivity, band emission of CH4, NH3, H2O, C2H2, CO and CO2 from
  the climate model's k-distributions, published V-T rates for the
  two-level bands; `--column-energy-treatment nlte|nlte_partial|lte|rotational`,
  `--nh3-h2-rate-analog cd4_bending|water_bending_overtone`,
  `--column-band-transfer thin|cooling_to_space|cooling_to_space_exchange`
  with its validity bound `--column-exchange-shift-max`), re-solving the
  chemistry on the new temperature; the profile header records
  `conducted_heat_flux_cgs`, `column_energy_treatment`, `column_T_match_K`,
  `nh3_h2_rate_analog` (treatment `nlte`), `column_band_transfer` (off the
  default `thin`) or, under `thin`, the `thin_layer_*` record of the
  optically thin validity check (the column is kept and flagged past its
  10 percent limit). `src/utils/element_flux_closure.py` iterates wind and
  chemistry to one elemental flux at the matching level (for hydrogen the
  hydrogen-nucleus flux `F_Hnuc` of `EXHALE_resolved.out`, not the
  component-1 mass flux `F_H`, which carries the metals too; the column
  models take the stellar energy the wind's input states, `--face-geometry
  wind`, checked against the wind's `dayside_dilution`) and, with
  `energy_coupling` on, passes the wind's `base_conductive_flux_cgs` to the
  next column (damped by `energy_coupling_omega`), converging only when the
  elemental fluxes and the conducted flux agree within `tol` and `T_match`
  moves less than `T_match_tol` (default 6 K)
- An independent conservation audit of the discrete rows, off by default.
  `EXHALE_CONSERVATION_BUDGET=<n>` exports the next `n` stationary residual
  assemblies to `output/conservation_budget_<nnnn>.txt`: per cell the face
  areas, the cell volume, both faces of all three conserved variables, the
  face potentials and the cell potential, the gravitational work the energy
  flux difference carries inside it, the face pressures or the equilibrium
  pressure departures the momentum branch in force uses, the explicit and
  transport sources, the radiative heating and cooling, and the assembled
  rows, every real at seventeen significant decimal digits. Schema 3 appends,
  for the upper face of each row, the left and right face states (rho, v, p)
  the Riemann problem was solved on, the well-balanced pressure data
  (`wb_dp_eq`, `wb_dev_L`, `wb_dev_R` and the jump `dp_WB` the Roe flux uses),
  the central mass transport 0.5 (F_L + F_R), the part of the numerical mass
  flux the pressure jump carries (the production `Num_flux` with the
  production jump minus the same call with the jump set to zero, an
  attribution and not a flux; NaN where the flux is not Roe or the option is
  off) and the recomputed `Num_flux` mass flux; per cell |R_mass|, the scale
  the stationary solve divides it by, their ratio and the signed running sum
  of V R_mass over the physical cells; and a header block for the base: the
  reservoir, the interior state the base condition reads at the face (the
  first cell's own constant-density equilibrium under `Well balanced: True`,
  its hydrostatic isentrope otherwise), the characteristic face velocity, both states of face 0, the numerical base
  mass flow and the wind-window flux. The export is a measurement: the run's
  outputs are bitwise the same with it on or off. Absent or `0`
  nothing is written and no file is opened.
  `EXHALE_conservation_budget.py` reads that file, rebuilds each row from
  the exported terms alone, and reports two things separately: assembly
  consistency, whose tolerance is arithmetic, and stationarity, whose
  tolerance is the certification's
- Python loaders ([`examples/exhale_io.py`](examples/exhale_io.py)) read the
  output column schemas. The Makefile runs the available executable
  `src/tests/*/run.sh` scripts plus its explicit scientific test targets.
  The v1.00 regression matrix and historical suite counts are not a statement
  of available or passing validation in v1.10. See the build limitations above.
- A restart is a contract, not a file copy. Both state files carry a
  `restart_schema 1` metadata block (reservoir, species columns, grid,
  constants, twenty-six option switches, physical time, source), and a load whose
  grid, reservoir, constants or options disagree with the input is refused by
  name rather than silently accepted. `Restart intent:` says what the loaded
  state is (`trajectory`, `relaxation`, `stationary`, which enters the steady
  solver with no time step, or `stationary evaluate`, which certifies the stored
  pair as loaded, takes no step at all and writes the state, the
  advection-corrected profiles and the mass-loss line from that pair, its exit
  status being the stored pair's verdict; one composition sweep from the
  loaded state is then taken and its movement and certification reported
  beside the verdict, never written. The same held-composition evaluation,
  the temperature, sources and lower ghost of the stored composition from one
  sweep that returns the interior composition bit for bit, is what every
  stationary outer pass and every written file is certified on, so a file's
  certificate is the one its reload measures), and `Restart option change:` names the
  option tokens a deliberate ladder is allowed to differ in. The two
  continuation-factor tokens, `cond` (`EXHALE_CONDUCTION_SCALE`) and
  `interdiff_enth` (`EXHALE_INTERDIFF_ENTH_SCALE`), carry a factor other than
  1 as `ES23.16E3`, so two factors are two states; a five-decimal token
  written before 2026-10-01 stands for the interval [v - 5e-6, v + 5e-6], and
  a requested factor inside it loads the state as a seed with its claim
  dropped (`legacy_factor_token`), one outside it is a physics change to be
  named. A file with no
  block is legacy, loaded as before and marked `provenance unknown`, and that
  mark is inherited by everything the run writes. The two files state the mass
  density twice -- the `rho` column, and the species densities that weigh it --
  and the conserved variable is the authority: the density is read from its own
  column and the loaded species are projected onto it, so the state a restart
  evaluates is the state the file names to the last bit, whatever mass closure
  the written composition carries. The departure is reported on every restart.
  The grid is part of the contract: since 2026-09-27 the first cell face of
  every grid type is the base level `r = 1` (the reservoir the lower
  atmosphere states is then the face state itself, with no half-cell slab of
  reservoir gas under the domain), so every cell center moved and a state
  written before then is refused; it enters a run only as a seed mapped with
  `src/utils/map_state_to_grid.py <src> <target_grid_file> <out> --ic`
  (`--reservoir El/H <value>` rescales an element, `--uniform` sets it to one
  ratio in every row; a target cell whose center lies between the base face
  and the source's first cell center takes that cell's values, never an
  interpolation with the reservoir ghost)
- Beside every equilibrium state pair (`Hydro_ioniz.txt`, `Ion_species.txt`,
  including the pass-state generations and the `EXHALE_DUMP_IC=1` dump) the
  run writes `output/conserved_state.txt`, the exact code-unit state: the
  conserved variables (rho, rho v, E) and every species fraction of every
  cell, ghosts included, as `ES25.17E3` text, which reproduces each binary64
  value exactly. It is copied with the pair to `conserved_state_IC.txt`, and
  a restart that finds it starts from the bits the run ended on instead of
  converting the dimensional pair (the conversion from physical units, the
  species fractions recovered from densities and the energy rebuilt from `p`
  move the last digits: 146 of 500 conserved rows and all 500 composition rows
  of a 20-step `wasp_he23off` state). The file records the format version, the
  shape (`N`, `Ng`, number of species), the rolling checksums of the two
  dimensional files it was written with, the reservoir and options fields of
  the restart metadata, the base boundary model, the species order, the
  normalization constants, checksums of the physical input files (spectrum,
  lower-atmosphere profile or `base.inp`, `metals.inp`, the opacity
  parameters and tables) and of the gravity and irradiation fields, the cell
  and face radii, and a checksum of its own data. It is used only when all of
  them agree with the run; otherwise, or when the load applied a `Restart
  option change:`, migrated a legacy factor token or read a pair without
  metadata, it is ignored with a printed notice and the restart proceeds from
  the dimensional pair, which stays the state of record for analysis and the
  restart contract. A molecular-seed conversion and a seed projected with
  `EXHALE_SEED_HYDROSTATIC_CELLS` never read it. Physical cells and the upper
  ghosts are restored; the lower ghosts are rebuilt by the boundary from the
  run's reservoir. A state that cannot be written exactly (a nonfinite value,
  a nonpositive density or internal energy, a negative species fraction) is
  not written, any older file at that path is removed, a warning is printed
  and the run's exit status is unaffected; a write of a pair without the
  conserved state (the molecular seed) also removes the older file. The
  graphical interface copies the file with the pair and deletes a stale
  `conserved_state_IC.txt` when the state has none;
  `src/utils/map_state_to_grid.py` refuses a destination that already holds
  one, since a mapped state is a new state
- A stationary run that alternates a wind solve with a composition relaxation
  publishes the state of every completed outer pass as one generation,
  `output/pass_state/<state_id>/` (`Hydro_ioniz.txt`, `Ion_species.txt`,
  `conserved_state.txt` when it can be written, `manifest.txt`), written by
  the final-state writers, with
  `<state_id> = r<run identity>_p<pass>`; the one-line file
  `output/pass_state/current` names the generation to take, and it is
  replaced by a single rename only after that generation is complete, so a
  run stopped from outside at any moment leaves `current` naming a whole
  state. The generation `current` named before is kept and older ones are
  removed. Copy the pair `current` names (and its `conserved_state.txt`) to
  the `_IC` names to continue from that pass. The pair is never a certified state: both halves say
  `certified=F cert_reason=pass_snapshot_p<n> ... state_id=<state_id>`
  (`load_IC` refuses halves of two generations), and a `# pass_snapshot`
  line records the carrier movement bound and the element under-relaxation
  factor of the next pass, which a restart does not inherit
  (`EXHALE_CARRIER_TRUST`, `EXHALE_DIFF_OMEGA` pass them back). Until
  2026-09-29 the pair was `output/{Hydro_ioniz,Ion_species}_last_pass.txt`
  (`../EXHALE_v1.00/md/Update_EXHALE_stage3.md` section 33)
- The advection-corrected `_adv` profiles the analysis and transit tools read
  say row by row what they are: two validity fields, `adv_T_status` and
  `adv_comp_status`, for the temperature and the composition separately, and
  `adv_mass_row`, the measure both were decided by. The post-process solves
  the steady ionization and energy equations along the recorded flow, second
  order in the cell width (variable-step BDF2 since 2026-09-27; the
  first-order recursion before it put 1-4 percent into the He 10830
  equivalent width of the LHS 1140 b states), each element's stages stepped
  on the velocity of that element's nuclei (`v_H`, `v_He`, which differ from
  the bulk velocity when `He_diffusion` moves the elements), with the heat
  conduction and the element, carrier and stage enthalpy flux divergences of
  the run kept in the energy equation it re-solves (the conduction couples
  neighboring cells, so the corrected column is solved together by a Newton
  iteration and accepted on its residual), and the
  flow it integrates along changes its face mass flux by a fraction of itself
  across each cell; the correction is first order in that fraction, so a row
  at or below the `# adv_conditional_tol` of the file is corrected and is a
  CONDITIONAL correction accurate to that fraction of itself, while a row
  above it, or one whose local radiative balance rather than the flow sets
  its temperature, or one the gas flows into, or one whose ionization is
  already equilibrated, or one whose energy balance is carried by the
  molecular channels the post-process omits, keeps the run's own state and
  says so. Whether the
  whole input state passed the stationary certification is a separate
  statement in the same header (`# adv_input_certified`), and the product is
  declared there as a one-way correction on a fixed density and velocity
  field, so a spectrum built on a breathing base is labeled as such. The
  `# adv_derived_state:` line of the same header records how the solves that
  produced the rows ended as a whole: the column energy solve of the
  transport terms (rejected, with its reason and triggering cell, or
  `converged_on_residual`: every unknown's energy row, formed again at the
  final profile, within 10 eps times the sum of the parts it is assembled
  from + 1e-6 of its largest term, with
  the worst cell and residual in the record), the chemistry, the row counts,
  the cells whose own marching energy step found no root (unknowns of the
  column by default) and `outer=unverified`; a rejected column writes every
  one of its cells `failed` and keeps the marching profile in the numeric
  columns. `unverified` is the best a record says: the column is solved with
  its transport coefficients formed on the input composition and the
  fine-structure escape probabilities of the start of the pass, and neither
  the full energy residual of the derived state, nor the convergence of the
  ten outer passes, nor the momentum balance of the corrected pressure is
  tested. A row where H+, He+ or He(2^3S) loses faster to H2 than to the
  channels the molecule-free advection systems carry keeps the run's state
  (status 3), and every refused field carries the run's own numbers exactly
  (the temperature of a failed row of a rejected column excepted). The
  `stationary evaluate` route exits 7 when that record is rejected (2 stays
  the certification refusal of the input), and `EXHALE_transit.py` refuses an
  `adv` spectrum from a rejected record, or from a file with no record (written
  before the record existed: an UNKNOWN derived state), with exit 7 unless
  `EXHALE_TRANSIT_DIAGNOSTIC=1`

---

## Requirements and build

| Component | Version |
|-----------|---------|
| Fortran compiler | `gfortran` >= 9.3 or `ifort`/`ifx` >= 2021 |
| Python 3 | >= 3.8; packages `numpy`, `scipy`, `matplotlib`, `tkinter`, `astropy` |
| MINPACK | included in `src/modules/nonlinear_system_solver/` |

Build from the v1.10 source directory after installing the required compiler
and Python dependencies.

```bash
cd /path/to/EXHALE_v1.10
make                                    # gfortran; make FC=ifort or FC=ifx
# make check requires the v1.00 regression harness, absent from this tree
make test                               # the assertion suites (physics, grid and gates, ...)
```
The Makefile takes the compiler from `PATH` and, when that compiler's prefix
carries an OpenBLAS (the conda-forge gfortran 16.2 of this machine does), links
that prefix's LAPACK and records it in the rpath -- so the binary and its
LAPACK share one libgfortran runtime; otherwise it links the system `-llapack`
(`make LAPACK_LIBS='...'` overrides). OpenBLAS runs its own thread pool,
independent of OpenMP and by default as wide as the machine: at startup the
binary sets it to one thread unless `OPENBLAS_NUM_THREADS` is stated, and the
setup report says which library it found and what it set (the band
factorizations of the stationary solve are too small to gain from a BLAS team:
MEASURED 16.33 s against 16.32 s on the `wasp_full_newton` reload at 1 and 8
BLAS threads, 2026-09-13). The compiler's path and version are part
of the rebuild stamp: objects of one gfortran are never linked by another.
The Intel compilers (`make FC=ifx`, `FC=ifort`) link MKL (`-qmkl=sequential`)
by the same rule; `make FC=ifx OBJDIR=build_ifx EXE=EXHALE_ifx.x` builds a
second binary beside the gfortran one.
```
cd /path/to/prepared_run
/path/to/EXHALE_v1.10/EXHALE.x  # reads ./input.inp, writes ./output/
```

The binary always reads `./input.inp` and writes `./output/` relative to the
current working directory, so a run is simply a directory. The v1.10
`examples/regional_convergence/` files are partial settings. Complete historical
configurations are available in [the v1.00 examples](../EXHALE_v1.00/examples/README.md);
check their physics, paths, and supported options before adapting them. `./run_EXHALE.sh` instead opens a Tk interface
that fills in the planetary parameters and builds and runs for you.

For a stationary-solve attempt, the inherited two-stage reconstruction and
Newton recipe is shown below. It is separate from regional stabilization, which
requires Newton solving to be disabled. Choose the outer radius for the actual
sonic point; the illustrative value 2.0 is not sufficient for every planet:

```
Reconstruction scheme:  PLM+WENO3
du_th [PLM,WENO3]:      0.5 1.0e-3
Solver:                 Newton
Outer radius [R_p]:     2.0
```

The first number of `du_th [PLM,WENO3]:` is the `du` at which the PLM stage
hands over to WENO3; the second only has to be smaller than the first (that
selects the two-stage mode) and is no longer a stop: the marching stop and
the JFNK hand-off read the face mass-flux spread (see Features). Check the line `Mach number of the
outermost cell` in the results: below 1 the outer ghost condition can set
Mdot, and the outer radius must be moved out. The inherited run recipes with
inputs, commands, run times and outputs are in
`../EXHALE_v1.00/examples/recipe_*/README.md`.

This recipe does not guarantee certification. A run that carries a species row,
an element row from `He_diffusion` or a carrier row from the molecular
carriers, is solved by alternating the wind solve with the composition
relaxations over outer passes; the v1.00 record reports certified states (the LHS 1140 b
molecular and atomic diffusion states among them, `cert_reason=certified_in_wind`,
the species rows gated at `r >= 1.20 R_p`), and some do not yet, what limits
them being measured and recorded. The `# coupling:` line of every state file
says which: quote a mass-loss rate only from a state that says `certified=T`,
or with that qualification.

The inherited compiler variants, physics keys, output-file schemas,
convergence recipes, and post-processing instructions are described in
[`README_HOWTO.md` (v1.00)](../EXHALE_v1.00/README_HOWTO.md).

---

## Documentation

The v1.10 user manual and physics document below contain current additions
before or alongside inherited v1.00 material. Links explicitly labeled v1.00
refer to historical notes and examples in the neighboring source tree.


**Start here**

- [`README_HOWTO.md` (v1.00)](../EXHALE_v1.00/README_HOWTO.md), task-oriented recipes: one entry per
  task, with the exact input lines, the expected output, and where the full
  documentation lives. All the details this file used to carry are there
- [`docs/EXHALE_user_manual.pdf`](docs/EXHALE_user_manual.pdf), the reference
  manual: every input key, every output column, the physics and the solver
- [`docs/EXHALE_physics_and_algorithms.pdf`](docs/EXHALE_physics_and_algorithms.pdf),
  the physical equations and the numerical algorithms as the code applies them
  (hydrodynamics, radiation, ionization, heating and cooling, molecules and
  metals, the marching scheme, boundary and initial conditions, the stationary
  solve), written from the source at the depth of the ATES paper; the section
  files are `docs/physics_overview/*.tex`
- [`../EXHALE_v1.00/docs/Update_EXHALE_stage3.pdf` (v1.00)](../EXHALE_v1.00/docs/Update_EXHALE_stage3.pdf): the current
  update log (stage 3, from 2026-09-23; source `../EXHALE_v1.00/md/Update_EXHALE_stage3.md`).
  It carries the code-size appendix against the original ATES and the list of
  source inherited unchanged from it (`../EXHALE_v1.00/docs/Update_EXHALE_appendix.tex`,
  measured 2026-10-05, after section 102, with `src/utils/codesize.py`: 236
  core Fortran files, 109,539 code lines, 31 times the 3,514 of ATES);
  [`../EXHALE_v1.00/docs/Update_EXHALE_stage2.pdf` (v1.00)](../EXHALE_v1.00/docs/Update_EXHALE_stage2.pdf): sections
  1-17 (2026-09-05 to 09-22);
  [`../EXHALE_v1.00/docs/Update_EXHALE_stage1.pdf` (v1.00)](../EXHALE_v1.00/docs/Update_EXHALE_stage1.pdf): sections 1-171, the dated changelog against
  the original ATES
- [`examples/README.md` (v1.00)](../EXHALE_v1.00/examples/README.md): what each of the numbered
  example configurations demonstrates, and the exact lines it adds
- [`README_photochem.md` (v1.00)](../EXHALE_v1.00/README_photochem.md): the Photochem build the
  lower-atmosphere profile handoff runs on: where the source comes from, the
  four corrections applied to it, and how the source and its environment are
  rebuilt. Neither is in this repository

**Physics and numerics notes**

- [`../EXHALE_v1.00/docs/cooling_formulas.pdf` (v1.00)](../EXHALE_v1.00/docs/cooling_formulas.pdf): the analytic
  CHIANTI cooling-coefficient fits and their accuracy
- [`../EXHALE_v1.00/docs/photoion_cross_sections.pdf` (v1.00)](../EXHALE_v1.00/docs/photoion_cross_sections.pdf) and
  [`../EXHALE_v1.00/docs/recombination_coefficients.pdf` (v1.00)](../EXHALE_v1.00/docs/recombination_coefficients.pdf):
  the H/He/He 2^3S atomic data and its benchmarks
- [`../EXHALE_v1.00/docs/lower_atmosphere_coupling.pdf` (v1.00)](../EXHALE_v1.00/docs/lower_atmosphere_coupling.pdf),
  the lower-atmosphere connection: analytic column, molecular chemistry,
  Lyman-Werner photodissociation, base infrared field, the H2/H2O/CO infrared
  bands, VULCAN handoff
- [`../EXHALE_v1.00/docs/molecular_hydrogen_treatment.pdf` (v1.00)](../EXHALE_v1.00/docs/molecular_hydrogen_treatment.pdf),
  how H2 is treated: it is one species with no (v,J) resolution, so every
  level distribution is an assumption; the assumptions, their sources and
  their validity ranges collected in one place
- [`../EXHALE_v1.00/docs/version_compare.pdf` (v1.00)](../EXHALE_v1.00/docs/version_compare.pdf), diffusive
  separation of He and metals, and its measured effect on He 10830
- [`../EXHALE_v1.00/docs/transmission_spectrum.pdf` (v1.00)](../EXHALE_v1.00/docs/transmission_spectrum.pdf): how the
  transit spectra are computed from the wind profiles
- [`../EXHALE_v1.00/docs/EXHALE_BC_and_IC.pdf` (v1.00)](../EXHALE_v1.00/docs/EXHALE_BC_and_IC.pdf) and
  [`../EXHALE_v1.00/docs/steady_solver_memo.pdf` (v1.00)](../EXHALE_v1.00/docs/steady_solver_memo.pdf): boundary and
  initial conditions, the convergence criteria, and the Newton-Krylov design
- [`../EXHALE_v1.00/docs/wind_ae_solver.pdf` (v1.00)](../EXHALE_v1.00/docs/wind_ae_solver.pdf): the included Wind-AE
  solver behind `IC mode: windae`
- [`../EXHALE_v1.00/md/viscosity_conduction.md` (v1.00)](../EXHALE_v1.00/md/viscosity_conduction.md), molecular
  viscosity and heat conduction: derivation and where they matter
- [`../EXHALE_v1.00/docs/code_comparison.pdf` (v1.00)](../EXHALE_v1.00/docs/code_comparison.pdf) and
  [`../EXHALE_v1.00/docs/methodology_comparison.pdf` (v1.00)](../EXHALE_v1.00/docs/methodology_comparison.pdf):
  comparison with ATES, Salz, Kubyshkina, Murray-Clay, AIOLOS, Taylor, Xing

`../EXHALE_v1.00/docs/` holds the inherited typeset notes on individual
investigations. The local `md/` directory holds v1.10 development reports.

---

## References

1. [Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., Rauscher,
   E. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...655A..30C). *Irradiation-driven escape of
   primordial planetary atmospheres I. The ATES photoionization hydrodynamics
   code.* A&A, 655, A30; and [(2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...663A.122C),
   *... II. Evaporation efficiency of sub-Neptunes through hot Jupiters.*
   A&A, 663, A122.  (The code EXHALE forks.)

2. [Biassoni, F., Caldiroli, A., Gallo, E., Haardt, F., Spinelli, R., Borsa,
   F. (2024)](https://ui.adsabs.harvard.edu/abs/2024A%26A...682A.115B). *Self-consistent modeling of
   metastable helium exoplanet transits.* A&A, 682, A115.

3. [Huang, C., Koskinen, T., Lavvas, P., Fossati, L. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...951..123H).
   *A Hydrodynamic Study of the Escape of Metal Species and Excited Hydrogen
   from the Atmosphere of WASP-121b.* ApJ, 951, 123.

4. [Koskinen, T. T., Lavvas, P., Huang, C., et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...929...52K).
   *Mass Loss by Atmospheric Escape from Extremely Close-in Planets.* ApJ,
   929, 52.  (Lower-atmosphere column and molecular rate coefficients.)

5. [Taylor, A. R., Koskinen, T., et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...989...68T).
   *A Multispecies Atmospheric Escape Model with Excited Hydrogen and Helium:
   Application to HD209458b.* ApJ, 989, 68.  (Diffusive-separation reference
   model; its He(2^3S)+H Penning fit was replaced by the Garcia Munoz 2025
   rate.)

6. [Xing, L., Yan, D., Guo, J. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...953..166X). *The Mass
   Fractionation of Helium in the Escaping Atmosphere of HD 209458b.* ApJ,
   953, 166.  (Multi-fluid He/H fractionation reference.)

7. [Verner, D. A., Ferland, G. J., Korista, K. T., Yakovlev, D. G.
   (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V). *Atomic Data for Astrophysics. II. New
   Analytic Fits for Photoionization Cross Sections of Atoms and Ions.* ApJ,
   465, 487.

8. [Murray-Clay, R. A., Chiang, E. I., Murray, N. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...693...23M).
   *Atmospheric Escape From Hot Jupiters.* ApJ, 693, 23; and
   [Broome, M. I., Murray-Clay, R., Vissapragada, S., et al.
   (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...995..198B). *Wind-AE: A Fast, Open-source 1D
   Photoevaporation Code with Metal and Multifrequency X-Ray Capabilities.*
   ApJ, 995, 198.  (The steady Parker-wind solver behind `IC mode: windae`.)

9. [Tsai, S.-M., Lyons, J. R., Grosheintz, L., et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJS..228...20T).
   *VULCAN: An Open-source, Validated Chemical Kinetics Python Code for
   Exoplanetary Atmospheres.* ApJS, 228, 20.  (Cite together with
   [Tsai et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...923..264T), ApJ, 923, 264, when the
   VULCAN pre-step is used.)

---

## Author

Kwang-Il Seon (KASI / UST)

Last updated: 2026-10-17 21:03 KST
