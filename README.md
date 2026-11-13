# EXHALE

**EXoplanetary Hydrodynamic Atmospheric Loss and Escape** (EXHALE) is a 1-D
radiation-hydrodynamics code that simulates photoionization-driven atmospheric
mass loss from irradiated exoplanets. It solves the spherically symmetric Euler
equations together with the ionization and thermal balance of a hydrogen/helium
atmosphere carrying trace metals and, optionally, molecules, and returns the
mass-loss rate, the wind structure, and the transit transmission spectra that
follow from it. EXHALE started from the ATES code (Caldiroli et al. 2021;
Biassoni et al. 2024), of which it keeps the outline of the time-marching
hydrodynamic loop; almost everything else was written for EXHALE: the metal
chemistry and cooling, the molecular and lower-atmosphere layers, the diffusive
separation of the elements, the Ly-alpha radiative transfer, the well-balanced
reconstruction and the base boundary, the certification of a steady state, the
Newton-Krylov steady-state solver, the heat conduction, and the transit
synthesis.

Every extension that adds a species, a transport process or a field is
opt-in; corrections to the physics of the ATES problem itself (the He 2^3S
triplet, secondary ionization, case A/B recombination, the He/H charge
exchange) are on by default, so a bare `input.inp` solves the atomic H/He
problem with that physics rather than reproducing the legacy ATES model bit
for bit.

## This branch: v1.10

This is EXHALE v1.10 (git branch `v1.10`). This README,
[`README_HOWTO.md`](README_HOWTO.md) and the user manual are the v1.00 files
with what differs in v1.10 added where it applies, in a sentence that begins
"v1.10:" or names the version. Where a v1.10 sentence and the inherited text
disagree, the v1.10 sentence holds. Numbers quoted in the inherited text are
v1.00 measurements and were not measured again on v1.10 unless a v1.10
sentence says so. The changes against v1.00 (items P1-P9: the thermochemical
closure of the stationary lower boundary, certification of the state as
written, transport of the molecular ions, the selectable certification region,
regional stabilization, and others) are listed in section 1.5 of the user
manual, "Source changes against v1.00".

Material that exists only in the v1.00 branch (most examples, the update logs,
the physics notes) is linked below on GitHub; HOWTO and manual name it by its
path in a side-by-side checkout, `../EXHALE_v1.00/...`.

## Requirements and build

A Fortran compiler (`ifx` by default; `ifort` >= 2021 or `gfortran` >= 9.3),
`make`, and, for the GUI and the analysis tools, Python 3 (>= 3.8) with
`numpy`, `scipy`, `matplotlib`, `tkinter` and `astropy`. MINPACK is included.
There is no installation step beyond cloning:

```bash
git clone https://github.com/seoncafe/EXHALE
cd EXHALE
git checkout v1.10
make                                    # ifx; make FC=gfortran or FC=ifort to override
```

Compiler variants, the BLAS/LAPACK linking and a separate GNU Fortran build:
[`README_HOWTO.md`, Build](README_HOWTO.md#build).

## Quick start

```bash
cd examples/regional_convergence && ../../EXHALE.x Salz.inp   # reads ./Salz.inp, writes ./Salz/
```

This runs HD 209458 b after Salz et al. (2016) as a finite-duration
physical-time integration (its run time was not measured for this README).
In general, run `EXHALE.x [input_file]` from the run directory; with no
argument it reads `input.inp` there. Output goes to `Output directory:` when
the input names one, else to `output/` for `input.inp` and to a directory
named after the input's stem otherwise. `./run_EXHALE.sh` instead opens a Tk
interface that fills in the planetary parameters, builds and runs.
Ready-made configurations: [`README_HOWTO.md`, Example
configurations](README_HOWTO.md#example-configurations) (the tutorial and
most numbered examples are in the v1.00 branch).

**What counts as a result.** Marching profile convergence is the default
scientific route, in v1.00 and v1.10 alike. Set `Run mode: phys` (the
executable's default, `init`, is initialization or continuation and
establishes no physical elapsed time), keep profile snapshots, and accept the
state when every physical cell stays within 0.004 dex in T and rho and within
1 percent in velocity over the final K = 3 flow-crossing times; ordinary
solver stops do not apply this test. A stationary JFNK solution
(`Solver: Newton`) and the v1.10 regional stabilization test
(`Regional convergence: True`) are optional and answer different questions;
failing to obtain a stationary root does not invalidate a marching result
that passes its criterion. Before quoting a mass-loss rate, check that the
outer face is supersonic. Recipe: [`README_HOWTO.md`, Run a physical-time
marching model](README_HOWTO.md#run-a-physical-time-marching-model); criteria:
manual sections 1.1 and 3.4.

## Inputs and outputs

A run is a directory. `input.inp` holds the keys; `metals.inp`, `opacity.inp`
and `base.inp` switch their physics on by being present. Every input key, with
its default and meaning, is in the user manual (section "Input files"); every
output file and column is in the manual (section "Output files") and, in
short, in [`README_HOWTO.md`, Output files](README_HOWTO.md#output-files).
`examples/exhale_io.py` loads every output file by the `# columns` header it
carries.

## Features

**Chemistry and radiation**

- Trace metals (C, N, O, Mg, Si, Ca, Na, K, S, Fe) inside the coupled ionization system, with recombination, collisional ionization and detailed-balanced charge exchange
- Metal-line cooling from CHIANTI v11, with the C/N/O fine-structure ground terms in statistical equilibrium and line trapping
- H and He recombination switching between case A and case B with optical depth; He I 584 A scattering; on-the-spot absorption of He recombination radiation
- He I 2^3S metastable triplet in the coupled solver, for the He 10830 line
- Updated photoionization data, metal inner shells with Auger and fluorescence relaxation
- He+ + H charge exchange, radiative and non-radiative, with its reverse
- One stellar spectrum type (power law, Planck or a table) for every band of the photon grid
- Secondary ionization by fast photoelectrons; optional partition of 10.2-30 eV photoelectrons
- Non-LTE H(n=2) and Ly-alpha transfer (escape probability or a LaRT field)
- Molecular chemistry: H2, H2+, H3+, HeH+, H3+ cooling and Lyman-Werner photodissociation
- Oxygen chemistry (OH, H2O, CO) with FUV photolysis, so the base H2/H partition is computed
- Transport of the molecular carriers, of the ionization stages and of the He 2^3S population
- Thermal infrared field of the lower atmosphere and molecular infrared bands
- Composition-dependent caloric equation of state and the chemical heat of the molecular network
- Diffusive separation of helium and metals, with the interdiffusion enthalpy flux

**Hydrodynamics and solvers**

- Second-order marching with PLM or WENO3 reconstruction and an HLLC, ROE or LLF flux
- Optional Jacobian-free Newton-Krylov stationary solver for all active balances
- Certification of a stationary state, one condition per active balance, with exit statuses
- v1.10: regional stabilization, a finite-duration test for physical-time marching
- Viscosity and heat conduction, low-Mach damping, an optional Shapiro filter
- Well-balanced hydrostatic option; dual-energy update for cold supersonic gas
- Spherical or Roche-potential geometry, with the outer face past L1 if needed
- Initial conditions: cold hydrostatic, warm Parker seed, automatic, or a Wind-AE warm start
- Lower-atmosphere connection: analytic column, `base.inp` handoff, or a VULCAN/Photochem profile
- OpenMP parallelization, bitwise identical to the serial result

**Tooling**

- Runtime configuration by input keys and file presence; runtime grid size
- Restarts with a metadata contract, an exact code-unit state file and outer-pass snapshots
- Advection-corrected `_adv` profiles with row-by-row validity fields
- `EXHALE_transit.py`: He 10830, Ly-alpha, H-alpha, H-beta, Mg II, Ca II, Na I D and O I 1302 transit spectra
- Fitting a model line to an observed one by equivalent width with a width-matching velocity kernel
- `src/utils/element_flux_closure.py`: wind and lower-atmosphere chemistry iterated to one elemental flux
- Conservation audit of the discrete rows (`EXHALE_conservation_budget.py`)
- Python loaders (`examples/exhale_io.py`) and a live plotter (`EXHALE_plots.py`)

Each item is described, with its keys, defaults and validity range, in the
user manual and, as a recipe, in [`README_HOWTO.md`](README_HOWTO.md).

## Directory layout

```
EXHALE/
|-- Makefile, run_EXHALE.sh   # build; optional GUI launcher
|-- src/                      # EXHALE_main.f90, modules/, utils/
|-- examples/                 # input configurations, exhale_io.py
|-- inputdata/, cooling_data/ # spectra and opacity tables; cooling-fit scripts
|-- docs/                     # user manual, physics document, notes
|-- EXHALE_transit.py, EXHALE_plots.py, EXHALE_conservation_budget.py
```

The full tree, with the third-party codes fetched on demand:
[`README_HOWTO.md`, Directory layout](README_HOWTO.md#directory-layout).

## Documentation

| Question | Document |
|---|---|
| How do I do X (exact lines, commands, expected output)? | [`README_HOWTO.md`](README_HOWTO.md) |
| What does each input key, output column and option mean? | [`docs/EXHALE_user_manual.pdf`](docs/EXHALE_user_manual.pdf) (section 1: what changed in v1.10) |
| Which equations and algorithms does the code apply? | [`docs/EXHALE_physics_and_algorithms.pdf`](docs/EXHALE_physics_and_algorithms.pdf), sources in `docs/physics_overview/` |
| What do the numbered example configurations add? | [`examples/README.md` (v1.00)](https://github.com/seoncafe/EXHALE/blob/v1.00/examples/README.md) |
| How is Photochem obtained and built for the profile handoff? | [`README_photochem.md` (v1.00)](https://github.com/seoncafe/EXHALE/blob/v1.00/README_photochem.md) |
| What changed, and when? | update logs (v1.00): [stage 3](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/Update_EXHALE_stage3.pdf) (from 2026-09-23), [stage 2](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/Update_EXHALE_stage2.pdf), [stage 1](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/Update_EXHALE_stage1.pdf); v1.10 keeps no update log of its own |

Physics and numerics notes (v1.00 branch):
[cooling formulas](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/cooling_formulas.pdf),
[photoionization cross sections](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/photoion_cross_sections.pdf),
[recombination coefficients](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/recombination_coefficients.pdf),
[lower-atmosphere coupling](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/lower_atmosphere_coupling.pdf),
[molecular hydrogen treatment](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/molecular_hydrogen_treatment.pdf),
[diffusive separation and He 10830](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/version_compare.pdf),
[transmission spectrum](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/transmission_spectrum.pdf),
[boundary and initial conditions](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/EXHALE_BC_and_IC.pdf),
[steady solver](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/steady_solver_memo.pdf),
[Wind-AE solver](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/wind_ae_solver.pdf),
[viscosity and conduction](https://github.com/seoncafe/EXHALE/blob/v1.00/md/viscosity_conduction.md),
comparison with other escape codes ([code](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/code_comparison.pdf),
[methodology](https://github.com/seoncafe/EXHALE/blob/v1.00/docs/methodology_comparison.pdf)).
The list with one line on each is in
[`README_HOWTO.md`, Where results and documents live](README_HOWTO.md#where-results-and-documents-live).

Notes in this branch: [`docs/jfnk_march_jfnk_lhs1140b.pdf`](docs/jfnk_march_jfnk_lhs1140b.pdf)
and [`docs/stationary_vs_marching_lhs1140b.pdf`](docs/stationary_vs_marching_lhs1140b.pdf)
(LHS 1140 b, stationary and marching solutions; sources `.tex` with class
`docs/my_memo.cls`), `docs/compare_vulcan_photochem.py`, and
`docs/defect_fixes_20261027/`.

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

Last updated: 2026-11-13 20:48 KST
