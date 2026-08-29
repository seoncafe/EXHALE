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
  collisional ionization, and Huang et al. (2023) charge exchange
- Metal-line cooling as closed-form analytic fits to CHIANTI v11 (C I/II,
  N I/II, O I/II, Mg I/II, Ca II, Na I, Fe II; 0.1–3% accuracy), with
  the split ground terms of C I, C II, N II and O I solved in exact statistical
  equilibrium at the local `(n_e, n_HI)` instead of the coronal limit, which
  gives the density-dependent saturation of the [C II] 158 um and [O I] 63 um
  floors
- Line trapping in the eight ground-term fine-structure lines of C I, C II,
  N II and O I (the line-center escape probability from the column above each
  cell enters the statistical-equilibrium solution as `A_ul -> beta*A_ul`),
  and a smooth cutoff of the coronal fits below their
  10^3 K validity floor, leaving only the explicit fine-structure terms
- He I 2³S metastable triplet in the coupled solver, with a
  temperature-dependent He(2³S)+H Penning-ionization rate (Taylor et al. 2025)
- Updated photoionization data: He I ground state from Verner et al. (1996),
  and a He I 2³S cross section extended past 60 eV against TOPbase
- Secondary ionization by fast photoelectrons (Shull & van Steenberg 1985)
- Non-LTE H(n=2) populations and Ly-alpha radiative transfer, from either a
  fast Neufeld escape-probability closure or a field imported from the LaRT
  Monte Carlo code, including the Ly-alpha emitted in situ within the wind
- Molecular chemistry: H2, H2+, H3+ and HeH+ in the coupled ionization
  equilibrium, with H2 photoionization opacity/heating, Miller et al. (2013)
  H3+ infrared cooling, and H2 photodissociation in the Lyman-Werner bands
  with Draine & Bertoldi (1996) self-shielding
- Thermal infrared field of the atmosphere below the base, so the molecular
  and fine-structure coolants return the net rate rather than the vacuum limit
- Diffusive separation of helium and metals: hydrogen and helium are
  transported as a two-component mixture, with bulk advection, binary
  diffusive settling in the computed ambipolar field, and an optional eddy
  term; each trace metal can diffuse independently. The friction is resolved
  by ionization stage -- hard sphere, polarization and Coulomb -- so an ion
  is held to the protons instead of settling at a neutral rate

**Hydrodynamics and solvers**

- Second-order marching with PLM or WENO3 reconstruction and HLLC fluxes,
  usable as a two-stage PLM -> WENO3 sequence
- Jacobian-free Newton-Krylov steady-state solver with PTC warm-up, SER ramp,
  and non-monotone (Grippo) line search
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
  itself
- OpenMP parallelization of the cell ionization sweep, bitwise identical to
  the serial result

**Tooling**

- Runtime configuration by file presence, not compile flags: `metals.inp`
  present means metals on, likewise `opacity.inp` and `base.inp`
- Runtime grid size: `Grid cells: <N>` in `input.inp` (default 500), so
  base-refinement studies run without a rebuild
- `EXHALE_transit.py` transmission post-processor: He I 10830 Å, Ly-alpha,
  H-alpha, H-beta, and the metal resonance doublets Mg II h&k, Ca II H&K and
  Na I D, with impact-parameter Voigt integration, instrument and rotation
  convolution, and an optional triaxial Roche geometry
- Python loaders (`examples/exhale_io.py`) driven by the `# columns` schema
  header every output file carries, and a bitwise regression harness over a
  five-case physics matrix

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
cd examples/tutorial && ../../EXHALE.x  # reads ./input.inp, writes ./output/
```

The binary always reads `./input.inp` and writes `./output/` relative to the
current working directory, so a run is simply a directory. `examples/` holds
sixteen ready-made configurations, one for each physics or solver option;
copy one as a starting point. `./run_EXHALE.sh` instead opens a Tk interface
that fills in the planetary parameters and builds and runs for you.

For a quantitative mass-loss rate, converge with the two-stage scheme and the
Newton finish:

```
Reconstruction scheme:  PLM+WENO3
du_th [PLM,WENO3]:      0.5 1.0e-3
Solver:                 Newton
```

Everything else — compiler variants, all opt-in physics keys, output-file
schemas, convergence recipes, post-processing — is in
[`README_HOWTO.md`](README_HOWTO.md).

---

## Documentation

**Start here**

- [`README_HOWTO.md`](README_HOWTO.md) — task-oriented recipes: one entry per
  task, with the exact input lines, the expected output, and where the full
  documentation lives. All the details this file used to carry are there
- [`docs/EXHALE_user_manual.pdf`](docs/EXHALE_user_manual.pdf) — the reference
  manual: every input key, every output column, the physics and the solver
- [`docs/Update_EXHALE.pdf`](docs/Update_EXHALE.pdf) — dated changelog against
  the original ATES, with a code-size appendix
- [`examples/README.md`](examples/README.md) — what each of the sixteen
  example configurations demonstrates, and the exact lines it adds
- [`README_photochem.md`](README_photochem.md) — the Photochem build the
  lower-atmosphere profile handoff runs on: where the source comes from, the
  four corrections applied to it, and how the source and its environment are
  rebuilt. Neither is in this repository

**Physics and numerics notes**

- [`docs/cooling_formulas.pdf`](docs/cooling_formulas.pdf) — the analytic
  CHIANTI cooling-coefficient fits and their accuracy
- [`docs/photoion_cross_sections.pdf`](docs/photoion_cross_sections.pdf) and
  [`docs/recombination_coefficients.pdf`](docs/recombination_coefficients.pdf)
  — the H/He/He 2³S atomic data and its benchmarks
- [`docs/lower_atmosphere_coupling.pdf`](docs/lower_atmosphere_coupling.pdf) —
  the lower-atmosphere connection: analytic column, molecular chemistry,
  Lyman-Werner photodissociation, base infrared field, VULCAN handoff
- [`docs/binary_diffusion_design.md`](docs/binary_diffusion_design.md),
  [`docs/design_hehe_diffusion.md`](docs/design_hehe_diffusion.md) and
  [`docs/version_compare.pdf`](docs/version_compare.pdf) — diffusive
  separation of He and metals, and its measured effect on He 10830
- [`docs/transmission_spectrum.pdf`](docs/transmission_spectrum.pdf) — how the
  transit spectra are computed from the wind profiles
- [`docs/EXHALE_BC_and_IC.pdf`](docs/EXHALE_BC_and_IC.pdf) and
  [`docs/steady_solver_memo.pdf`](docs/steady_solver_memo.pdf) — boundary and
  initial conditions, the convergence criteria, and the Newton-Krylov design
- [`docs/wind_ae_solver.pdf`](docs/wind_ae_solver.pdf) — the included Wind-AE
  solver behind `IC mode: windae`
- [`docs/viscosity_conduction.md`](docs/viscosity_conduction.md) — molecular
  viscosity and heat conduction: derivation and where they matter
- [`docs/code_comparison.pdf`](docs/code_comparison.pdf) and
  [`docs/methodology_aiolos_taylor_xing.pdf`](docs/methodology_aiolos_taylor_xing.pdf)
  — comparison with ATES, Salz, Kubyshkina, Murray-Clay, AIOLOS, Taylor, Xing

`docs/` holds roughly forty further memos on individual investigations;
[`TO_BE_DONE.md`](TO_BE_DONE.md) is the open-items list.

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

Last updated: 2026-08-30 00:18
