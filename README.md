# EXHALE

**EXoplanetary Hydrodynamic Atmospheric Loss and Escape** (EXHALE) is a 1-D
radiation-hydrodynamics code that simulates photoionization-driven atmospheric
mass loss from irradiated exoplanets. It solves the spherically symmetric Euler
equations together with the ionization and thermal balance of a hydrogen/helium
atmosphere carrying trace metals and, optionally, molecules, and returns the
mass-loss rate, the wind structure, and the transit transmission spectra that
follow from it.

EXHALE started from the ATES code (Caldiroli et al. 2021; Biassoni et al.
2024), of which it keeps the outline of the time-marching hydrodynamic loop;
almost everything else was written for EXHALE: the metal chemistry and
cooling, the molecular and lower-atmosphere layers, the diffusive separation
of the elements, the Ly-alpha radiative transfer, the well-balanced
reconstruction and the base boundary, the certification of a steady state,
the Newton-Krylov steady-state solver, the heat conduction, and the transit
synthesis. Every extension that adds a species, a transport process or a
field is opt-in. Corrections to the physics of the ATES problem itself (the
He 2^3S triplet, secondary ionization, case A/B recombination, the He/H
charge exchange) are on by default, so a bare `input.inp` solves the atomic
H/He problem with that physics rather than reproducing the legacy ATES model
bit for bit.

---

## Features

One line each; the input keys, defaults, validity ranges and measurements are
in the [user manual](docs/EXHALE_user_manual.pdf), the recipes in
[`README_HOWTO.md`](README_HOWTO.md).

**Chemistry and radiation**

- Ten trace metals (C, N, O, Mg, Si, Ca, Na, K, S, Fe) solved inside the coupled
  H/He ionization system, with radiative and dielectronic recombination,
  collisional ionization, and charge exchange held to detailed balance
- Metal-line cooling from CHIANTI 11, with the ground-term fine-structure
  lines of C I, C II, N II and O I in statistical equilibrium and line trapping
- Case A/B hydrogen and helium recombination by local escape, with the
  recombination radiation and the He I 584 A line absorbed on the spot
- He I 2^3S metastable triplet in the coupled solver (the He 10830 line)
- Photoionization including the metal inner shells and their Auger electrons
- He/H charge exchange in every helium-bearing system
- Stellar spectrum as a power law, a photospheric blackbody, or a loaded table
  covering every band of the photon grid
- Secondary ionization by fast photoelectrons, and an optional partition of
  the 10.2-30 eV photoelectrons
- Non-LTE H(n=2) populations and Ly-alpha transfer (escape probability or an
  imported field)
- Molecular chemistry (H2, H2+, H3+, HeH+) with H3+ infrared cooling and
  Lyman-Werner photodissociation of H2
- Oxygen chemistry (OH, H2O, CO) with FUV photolysis, so the base H2/H
  partition is computed rather than imported
- Transport with the flow of the molecular carriers, the ionization stages and
  the He 2^3S population (each optional)
- Thermal infrared field of the lower atmosphere and the H2, H2O and CO
  infrared bands
- Composition-dependent caloric equation of state built from the H2 level
  ladder, and the chemical heat of the molecular network
- Diffusive separation of helium and the metals, with the interdiffusion
  enthalpy flux in the energy equation

**Hydrodynamics and solvers**

- Second-order marching (PLM or WENO3; HLLC, Roe or local Lax-Friedrichs
  interface flux), with physical-time (`Run mode: phys`) and initialization
  modes
- Optional Jacobian-free Newton-Krylov steady-state solver and a
  certification of every claimed steady state, with documented exit statuses
- Viscosity and heat conduction, low-Mach damping, an optional Shapiro filter,
  an optional well-balanced reconstruction, and a dual-energy update for cold
  supersonic gas
- Spherical or Roche-potential domain; the outer face can be placed past the
  sonic point, and the outer Mach number is reported with every mass-loss rate
- Initial conditions: cold hydrostatic, warm Parker seed, an automatic
  selector, or a warm start from the included Wind-AE solver
- Lower-atmosphere connection: an analytic column, a `base.inp` handoff, a
  lower-atmosphere profile, or a VULCAN or Photochem photochemistry model
- OpenMP parallelization of the cell sweeps, bitwise identical to the serial
  result

**Configuration, restart and analysis**

- Runtime configuration: plain-text input keys, and `metals.inp`,
  `opacity.inp`, `base.inp` switched on by their presence; runtime grid size
- Checked restarts with an exact code-unit state file, and a restartable
  state after every stationary outer pass
- Advection-corrected `_adv` profiles that state the validity of every row
- `EXHALE_transit.py`: He I 10830, Ly-alpha, H-alpha, H-beta, Mg II, Ca II,
  Na I D and O I 1302 transit spectra, with instrument and rotation
  convolution and optional Roche geometry
- Python loaders for every output file (`examples/exhale_io.py`), a live
  profile plotter, and an independent conservation audit of the discrete rows

---

## Requirements and build

| Component | Version |
|-----------|---------|
| Fortran compiler | `gfortran` >= 9.3 or `ifort`/`ifx` >= 2021 |
| Python 3 | >= 3.8; packages `numpy`, `scipy`, `matplotlib`, `tkinter`, `astropy` |
| MINPACK | included in `src/modules/nonlinear_system_solver/` |

```bash
git clone https://github.com/seoncafe/EXHALE
cd EXHALE
make                     # ifx (default); make FC=gfortran or FC=ifort to override
```

Compiler variants, BLAS/LAPACK selection and the other `make` targets:
[`README_HOWTO.md`, Build](README_HOWTO.md#build).

---

## Quick start

```bash
cd examples/tutorial && ../../EXHALE.x  # reads ./input.inp, writes ./output/
```

`EXHALE.x [input_file]` reads `input.inp` in the working directory when no
file is named. Files named in the input (`metals.inp`, `base.inp`,
`opacity.inp`, the lower-atmosphere profile, the spectrum) are found beside
it. Output goes to `output/` for `input.inp` (another input file writes to a
directory named after its stem), or to the directory set by
`Output directory: <path>`; `EXHALE_setup.out` there echoes the configuration
the run used. `examples/` holds nineteen ready-made configurations, one for
each physics or solver option; `./run_EXHALE.sh` instead opens a Tk interface
that fills in the planetary parameters, builds and runs. Details:
[Run](README_HOWTO.md#run).

**Before quoting a result.** Marching profile convergence is the default
route to a result. For a physical-time result set `Run mode: phys`, keep
profile snapshots, and check the temperature, density and velocity bands over
the final three flow-crossing times
([recipe](README_HOWTO.md#run-a-physical-time-marching-model)); the default
`Run mode: init` is initialization or continuation and makes no such claim.
`Solver: Newton` adds an optional stationary finish; failing to find a
stationary root does not invalidate a marching result that passes its
criterion. Treat the mass-loss rate as quantitative only when the outer face
is supersonic: below Mach 1.1 the run labels it "NOT A QUANTITATIVE MDOT".

---

## Inputs and outputs

- Input keys: [user manual](docs/EXHALE_user_manual.pdf), section "Input
  files" (`input.inp` core block and optional keys, `metals.inp`,
  `opacity.inp`, the lower-atmosphere files)
- Output files and their columns: user manual, section "Output files";
  summary in [`README_HOWTO.md`, Output files](README_HOWTO.md#output-files)
- Restarts: [Restart a run](README_HOWTO.md#restart-a-run-load-ic-restart-intent-restart-option-change)
- Python access: [Reading output in Python](README_HOWTO.md#reading-output-in-python)
- Transit spectra: [Transmission spectra](README_HOWTO.md#transmission-spectra-exhale_transitpy)

---

## Directory layout

```
Makefile, run_EXHALE.sh   build; optional GUI launcher
src/                      Fortran sources (EXHALE_main.f90, modules/) and utils/
examples/                 ready-made configurations, tutorial, exhale_io.py
inputdata/, cooling_data/ opacity and SED samples; cooling-fit scripts
docs/, md/                manuals, physics notes, update logs; Markdown memos
EXHALE_transit.py, EXHALE_plots.py, roche_recon.py   analysis tools
```

Full layout: [Directory layout](README_HOWTO.md#directory-layout).

---

## Documentation

| Question | Document |
|----------|----------|
| How do I do X (exact lines and commands)? | [`README_HOWTO.md`](README_HOWTO.md) |
| What does each input key, output column and tool do? | [`docs/EXHALE_user_manual.pdf`](docs/EXHALE_user_manual.pdf) |
| Which equations and algorithms does the code apply? | [`docs/EXHALE_physics_and_algorithms.pdf`](docs/EXHALE_physics_and_algorithms.pdf) |
| What changed, when and why? | [`docs/Update_EXHALE_stage3.pdf`](docs/Update_EXHALE_stage3.pdf) (current), [stage 2](docs/Update_EXHALE_stage2.pdf), [stage 1](docs/Update_EXHALE_stage1.pdf) |
| What does each example demonstrate? | [`examples/README.md`](examples/README.md) |
| How is the Photochem lower atmosphere built? | [`README_photochem.md`](README_photochem.md) |

Topic notes (cooling fits, atomic data, lower atmosphere, H2 treatment,
diffusion, transit spectra, boundary conditions, steady solver, Wind-AE, code
comparisons) are listed in
[Where results and documents live](README_HOWTO.md#where-results-and-documents-live).

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

Last updated: 2026-11-13 20:43 KST
