# EXHALE

**EXoplanetary Hydrodynamic Atmospheric Loss and Escape** (EXHALE) is a 1-D
radiation-hydrodynamics code that simulates photoionization-driven atmospheric
mass loss from irradiated exoplanets. It solves the spherically symmetric Euler
equations together with the ionization and thermal balance of a hydrogen/helium
atmosphere carrying trace metals and, optionally, molecules, and returns the
mass-loss rate, the wind structure, and the transit transmission spectra that
follow from it. EXHALE is a heavily extended fork of the ATES code (Caldiroli
et al. 2021; Biassoni et al. 2024): the metal chemistry and cooling, the
molecular and lower-atmosphere layers, the diffusive separation of helium, the
Ly-alpha radiative transfer, and the Newton-Krylov steady-state solver are new.
Every extension is opt-in, so a bare `input.inp` still reproduces the legacy
ATES model.

---

## Features

**Chemistry and radiation**

- Trace metals (C, N, O, Mg, Si, Ca, Na, K, S, Fe) solved self-consistently
  inside the coupled ionization system, with Badnell RR+DR recombination
  (Huang et al. 2023 fits for Fe I/II and Shull & Van Steenberg for Ca I, whose
  isoelectronic sequences the Badnell project does not reach), Voronov
  collisional ionization, and Huang et al. (2023) charge exchange, extended by
  the `O2+ + H0 -> O+ + H+` electron capture of Barragan et al. (2006) that
  Table 4 omits (on by default; `metals.inp: cx_O2p_H 0` restores the
  Table-4-only reaction set)
- Metal-line cooling as closed-form analytic fits to CHIANTI v11 (C I/II,
  N I/II, O I/II, Mg I/II, Ca II, Na I, Fe II; 0.1-3% accuracy), with
  the split ground terms of C I, C II, N II and O I solved in exact statistical
  equilibrium at the local `(n_e, n_HI)` instead of the coronal limit, which
  gives the density-dependent saturation of the [C II] 158 um and [O I] 63 um
  floors
- Line trapping in the eight ground-term fine-structure lines of C I, C II,
  N II and O I (the line-center escape probability from the column above each
  cell enters the statistical-equilibrium solution as `A_ul -> beta*A_ul`),
  and a smooth cutoff of the coronal fits below their
  10^3 K validity floor, leaving only the explicit fine-structure terms
- He I 2³S metastable triplet in the coupled solver, with
  temperature-dependent He(2³S)+H and +H₂ destruction rates
  (García Muñoz 2025, continuous closed forms with a 0.9/0.1
  Penning/associative branching)
- Updated photoionization data: He I ground state from Verner et al. (1996),
  and a He I 2³S cross section extended past 60 eV against TOPbase
- One stellar spectrum type builds every band of the photon grid, the XUV and
  the part below 13.6 eV where the He 2³S metastable and the low-IP metals
  absorb: `Spectrum type: Power-law` is the power law everywhere,
  `Spectrum type: Planck` is the photospheric blackbody
  `pi B_nu(T_eff) (R_star/a)^2` built from `Stellar Teff` and `Stellar radius`,
  and `Spectrum type: Load` is the table everywhere. A loaded table that stops
  above the lowest threshold of an active absorber (4.80 eV = 2583 Å with
  `Include He23S? True`, or a neutral metal's threshold) stops the run rather
  than leaving that absorber without a field; the startup report states the
  type, the source of the band below 13.6 eV, and the integrated grid flux
- Secondary ionization by fast photoelectrons (Shull & van Steenberg 1985)
- Non-LTE H(n=2) populations and Ly-alpha radiative transfer, from either a
  fast Neufeld escape-probability closure or a field imported from the LaRT
  Monte Carlo code, including the Ly-alpha emitted in situ within the wind
- Molecular chemistry: H2, H2+, H3+ and HeH+ in the coupled ionization
  equilibrium, with H2 photoionization opacity/heating, Miller et al. (2013)
  H3+ infrared cooling, and H2 photodissociation in the Lyman-Werner bands
  from a level-resolved CLOUDY calculation tabulated on our own (T, n_H,
  N_H2) grid -- the dissociation cross section itself, with the line
  self-shielding and the trapping of the fluorescent decay photons inside it,
  rather than a published closed-form fit
- Oxygen chemistry: OH, H2O and CO in the same coupled system, with the FUV
  photolysis of H2O and OH in four bands, so the base H2/H partition is
  computed rather than imported. Rates from Baulch et al. (2005) and the IUPAC
  evaluations, reverse rates by detailed balance against a NIST-JANAF Shomate
  table, CO carried as an oxygen reservoir. The first band is the 912-1201 A
  Lyman-Werner interval, where H2, H2O and OH share one beam: the H2 lines and
  the H2O/OH continuum each attenuate what the other sees, so the interval has
  one incident flux (`Stellar LW flux`). The molecular carriers H2, OH, H2O and
  CO are transported by default (`Molecular carrier transport`), an
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
  hydrogen ionization state can ride the same operator (`Ionization
  transport: True`, default off): H+ becomes a fifth carrier and the sweep is
  handed the transported fraction where a parcel leaves its shell faster
  than it ionizes (`P r/|v| < 1`), which is what the Koskinen et al. (2022)
  comparison of `docs/k22_electron_density_excess.md` needed
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
  is held to the protons instead of settling at a neutral rate

**Hydrodynamics and solvers**

- Second-order marching with PLM or WENO3 reconstruction and a stated
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
- `Well balanced: True` (default False) carries the departure from each
  cell's own local hydrostatic equilibrium through the reconstruction, the
  Riemann jumps and the pressure force, so a discrete hydrostatic
  equilibrium is preserved to rounding (Kappeli and Mishra 2014, 2016). A
  measurement option: turning it on moves every result, and the
  species-row stationary solves do not improve under it
  (`docs/well_balanced_flux_difference_design_20260910.md`)
- A steady state is not "converged", it is CERTIFIED: one evaluator states
  one condition per active balance and the run prints, and the state file
  records, which conditions were met and which row and cell refused. The
  tolerances are anchored by measurement, not inherited from the solver's own
  target: mass 3e-12, momentum 1e-8, energy 1e-6, the closure and the level
  balances at the tolerance the cell sweep already judges a root by, and the
  element and carrier rows at 1e-5 for `r >= 1.20 R_p`, the rows below that
  reported and not gating, because the same balance is a cancellation of
  advective terms in the wind and of eddy and settling terms in the layer and
  is not resolvable to one number in both
- Molecular transport: the Navier-Stokes viscous force, its dissipation, and
  heat conduction, integrated Crank-Nicolson and entering the steady residual
  with the same operator
- Low-Mach damping, a gated fourth-difference stress for a shell that has
  stopped flowing, and an optional Shapiro filter for a breathing base
- Roche-potential geometry, in either a spherical or a Roche-lobe domain
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
- OpenMP parallelization of the cell ionization sweep, bitwise identical to
  the serial result

**Tooling**

- Runtime configuration by file presence, not compile flags: `metals.inp`
  present means metals on, likewise `opacity.inp` and `base.inp`
- Runtime grid size: `Grid cells: <N>` in `input.inp` (default 500), so
  base-refinement studies run without a rebuild
- `EXHALE_transit.py` transmission post-processor: He I 10830 Å, Ly-alpha,
  H-alpha, H-beta, the metal resonance doublets Mg II h&k, Ca II H&K and
  Na I D, and the O I 1302/1304/1306 triplet out of its three resolved
  ground-term fine-structure levels, with impact-parameter Voigt integration,
  instrument and rotation convolution, and an optional triaxial Roche geometry.
  Every `tpm_*.txt` carries a `# transit_schema 1` metadata block above its
  column line: which `_adv` profile it was built from, that profile's validity
  schema, certification verdict and provenance, the census of rows the line
  took its depth from that were refused a correction, and the
  `EXHALE_TRANSIT_*` overrides in effect, so a curve says what it stands on
- An independent conservation audit of the discrete rows, off by default.
  `EXHALE_CONSERVATION_BUDGET=<n>` exports the next `n` stationary residual
  assemblies to `output/conservation_budget_<nnnn>.txt`: per cell the face
  areas, the cell volume, both faces of all three conserved variables, the
  face potentials and the cell potential, the gravitational work the energy
  flux difference carries inside it, the face pressures or the equilibrium
  pressure departures the momentum branch in force uses, the explicit and
  transport sources, the radiative heating and cooling, and the assembled
  rows, every real at seventeen significant decimal digits. Absent or `0`
  nothing is written and no file is opened.
  `EXHALE_conservation_budget.py` reads that file, rebuilds each row from
  the exported terms alone, and reports two things separately: assembly
  consistency, whose tolerance is arithmetic, and stationarity, whose
  tolerance is the certification's
- Python loaders (`examples/exhale_io.py`) driven by the `# columns` schema
  header every output file carries; a bitwise regression harness over a
  sixteen-case physics matrix (`make check`); and 33 assertion suites under
  `src/tests/` plus three standalone test programs (`make test`), which print
  one `PASS|FAIL <name> measured= reference= tol=` line per assertion
- A restart is a contract, not a file copy. Both state files carry a
  `restart_schema 1` metadata block (reservoir, species columns, grid,
  constants, twenty option switches, physical time, source), and a load whose
  grid, reservoir, constants or options disagree with the input is refused by
  name rather than silently accepted. `Restart intent:` says what the loaded
  state is (`trajectory`, `relaxation`, `stationary`, which enters the steady
  solver with no time step, or `stationary evaluate`, which measures the state
  and takes no step at all and writes the advection-corrected profiles and the
  mass-loss line from that measurement), and `Restart option change:` names the
  option tokens a deliberate ladder is allowed to differ in. A file with no
  block is legacy, loaded as before and marked `provenance unknown`, and that
  mark is inherited by everything the run writes. The two files state the mass
  density twice -- the `rho` column, and the species densities that weigh it --
  and the conserved variable is the authority: the density is read from its own
  column and the loaded species are projected onto it, so the state a restart
  evaluates is the state the file names to the last bit, whatever mass closure
  the written composition carries. The departure is reported on every restart
- The advection-corrected `_adv` profiles the analysis and transit tools read
  say row by row what they are: two validity fields, `adv_T_status` and
  `adv_comp_status`, for the temperature and the composition separately, and
  `adv_mass_row`, the measure both were decided by. The post-process solves
  the steady ionization and energy equations along the recorded flow, and the
  flow it integrates along changes its face mass flux by a fraction of itself
  across each cell; the correction is first order in that fraction, so a row
  at or below the `# adv_conditional_tol` of the file is corrected and is a
  CONDITIONAL correction accurate to that fraction of itself, while a row
  above it, or one whose local radiative balance rather than the flow sets
  its temperature, or one the gas flows into, or one whose ionization is
  already equilibrated, keeps the run's own state and says so. Whether the
  whole input state passed the stationary certification is a separate
  statement in the same header (`# adv_input_certified`), and the product is
  declared there as a one-way correction on a fixed density and velocity
  field, so a spectrum built on a breathing base is labeled as such

---

## Requirements and build

| Component | Version |
|-----------|---------|
| Fortran compiler | `gfortran` >= 9.3 or `ifort`/`ifx` >= 2021 |
| Python 3 | >= 3.8; packages `numpy`, `scipy`, `matplotlib`, `tkinter`, `astropy` |
| MINPACK | included in `src/modules/nonlinear_system_solver/` |

There is no installation step beyond cloning.

```bash
git clone https://github.com/seoncafe/EXHALE
cd EXHALE
make                                    # gfortran; make FC=ifort or FC=ifx
make check                              # bitwise regression over the case matrix
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
cd examples/tutorial && ../../EXHALE.x  # reads ./input.inp, writes ./output/
```

The binary always reads `./input.inp` and writes `./output/` relative to the
current working directory, so a run is simply a directory. `examples/` holds
nineteen ready-made configurations, one for each physics or solver option;
copy one as a starting point. `./run_EXHALE.sh` instead opens a Tk interface
that fills in the planetary parameters and builds and runs for you.

For a quantitative mass-loss rate, converge with the two-stage scheme and the
Newton finish:

```
Reconstruction scheme:  PLM+WENO3
du_th [PLM,WENO3]:      0.5 1.0e-3
Solver:                 Newton
```

An atomic run finished this way certifies. A run that carries a species row,
an element row from `He_diffusion` or a carrier row from the molecular
carriers, does not: no such configuration has yet reached its certification
tolerance on any route, and what limits them is measured and recorded in
[`docs/code_status_20260910.md`](docs/code_status_20260910.md) section 3.3.
Quote a mass-loss rate from such a run only with that qualification.

Everything else (compiler variants, all opt-in physics keys, output-file
schemas, convergence recipes, post-processing) is in
[`README_HOWTO.md`](README_HOWTO.md).

---

## Documentation

**Start here**

- [`docs/code_status_20260910.md`](docs/code_status_20260910.md): the state of
  the code physics by physics (implemented / verified / limited / missing),
  what is certified and what limits the rest, the ordered problems and the
  task list. Read this one first in a new session
- [`docs/ISSUES_20260909.md`](docs/ISSUES_20260909.md): every problem the
  2026-09 work exposed, resolved (a table) and open (with its evidence and
  its next item)
- [`README_HOWTO.md`](README_HOWTO.md), task-oriented recipes: one entry per
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
- [`docs/Update_EXHALE_stage2.md`](docs/Update_EXHALE_stage2.md): the current update log (stage 2, from 2026-09-05), whose PDF carries the
  code-size appendix against the original ATES and the list of source inherited unchanged from it
  (`docs/Update_EXHALE_appendix.tex`, re-measured with `src/utils/codesize.py`);
  [`docs/Update_EXHALE_stage1.pdf`](docs/Update_EXHALE_stage1.pdf): sections 1-171, the dated changelog against
  the original ATES
- [`examples/README.md`](examples/README.md): what each of the numbered
  example configurations demonstrates, and the exact lines it adds
- [`README_photochem.md`](README_photochem.md): the Photochem build the
  lower-atmosphere profile handoff runs on: where the source comes from, the
  four corrections applied to it, and how the source and its environment are
  rebuilt. Neither is in this repository

**Physics and numerics notes**

- [`docs/cooling_formulas.pdf`](docs/cooling_formulas.pdf): the analytic
  CHIANTI cooling-coefficient fits and their accuracy
- [`docs/photoion_cross_sections.pdf`](docs/photoion_cross_sections.pdf) and
  [`docs/recombination_coefficients.pdf`](docs/recombination_coefficients.pdf):
  the H/He/He 2³S atomic data and its benchmarks
- [`docs/lower_atmosphere_coupling.pdf`](docs/lower_atmosphere_coupling.pdf),
  the lower-atmosphere connection: analytic column, molecular chemistry,
  Lyman-Werner photodissociation, base infrared field, the H2/H2O/CO infrared
  bands, VULCAN handoff
- [`docs/molecular_hydrogen_treatment.pdf`](docs/molecular_hydrogen_treatment.pdf),
  how H2 is treated: it is one species with no (v,J) resolution, so every
  level distribution is an assumption; the assumptions, their sources and
  their validity ranges collected in one place
- [`docs/binary_diffusion_design.md`](docs/binary_diffusion_design.md),
  [`docs/design_hehe_diffusion.md`](docs/design_hehe_diffusion.md) and
  [`docs/version_compare.pdf`](docs/version_compare.pdf), diffusive
  separation of He and metals, and its measured effect on He 10830
- [`docs/transmission_spectrum.pdf`](docs/transmission_spectrum.pdf): how the
  transit spectra are computed from the wind profiles
- [`docs/EXHALE_BC_and_IC.pdf`](docs/EXHALE_BC_and_IC.pdf) and
  [`docs/steady_solver_memo.pdf`](docs/steady_solver_memo.pdf): boundary and
  initial conditions, the convergence criteria, and the Newton-Krylov design
- [`docs/wind_ae_solver.pdf`](docs/wind_ae_solver.pdf): the included Wind-AE
  solver behind `IC mode: windae`
- [`docs/viscosity_conduction.md`](docs/viscosity_conduction.md), molecular
  viscosity and heat conduction: derivation and where they matter
- [`docs/code_comparison.pdf`](docs/code_comparison.pdf) and
  [`docs/methodology_comparison.pdf`](docs/methodology_comparison.pdf):
  comparison with ATES, Salz, Kubyshkina, Murray-Clay, AIOLOS, Taylor, Xing

- [`docs/PLAN_20260909_rev1.md`](docs/PLAN_20260909_rev1.md) and
  [`docs/worker_rules.md`](docs/worker_rules.md): the plan the open items
  belong to, and the rules every worker on this tree follows
- [`docs/certification_tolerance_anchoring_20260910.md`](docs/certification_tolerance_anchoring_20260910.md):
  the five measured anchors behind the element and carrier tolerances, and
  why the wind and the layer cannot share one number
- [`docs/steady_solver_design.md`](docs/steady_solver_design.md): the
  stationary solver section by section, ending with what the design became

`docs/` holds roughly forty further memos on individual investigations;
[`docs/TO_BE_DONE.md`](docs/TO_BE_DONE.md) is the open-items list.

---

## References

1. [Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., Rauscher,
   E. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...655A..30C). *Irradiation-driven escape of
   primordial planetary atmospheres I. The ATES photoionization hydrodynamics
   code.* A&A, 655, A30; and [(2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...663A.122C),
   *… II. Evaporation efficiency of sub-Neptunes through hot Jupiters.*
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
   Application to HD209458b.* ApJ, 989, 68.  (Temperature-dependent Penning
   rate; diffusive-separation reference model.)

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

Last updated: 2026-09-22 12:07
