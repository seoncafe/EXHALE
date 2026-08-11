# EXHALE

**EXoplanetary Hydrodynamic Atmospheric Loss and Escape** (EXHALE) is a
1-D radiation-hydrodynamics code that simulates photoionization-driven
atmospheric mass loss from irradiated exoplanets.  It is a heavily extended
fork of the ATES code (Caldiroli et al. 2021; Biassoni et al. 2024), adding:

- Trace metals (C, N, O, Mg, Si, Ca, Na, K, S, Fe) solved self-consistently
  inside the coupled ionization system, with Badnell RR+DR recombination,
  Voronov collisional ionization, and Huang et al. (2023) charge exchange
- Metal-line cooling as **closed-form analytic formulas fitted to CHIANTI
  v11** (C I/II, N I/II, O I/II, Mg I/II, Ca II, Na I, Fe II; 0.1–3%
  accuracy), with density-dependent saturation of the [C II] 158 um /
  [O I] 63 um fine-structure floors and a 2-D statistical-equilibrium
  Fe II coefficient — see `docs/cooling_formulas.pdf`
- **Line trapping in the two fine-structure coolants**: the [O I] 63 um and
  [C II] 158 um line-center escape probability is computed from the column
  above each cell (Hollenbach & McKee 1979 / de Jong, Boland & Dalgarno 1980
  form, renormalized so beta(0) = 1) and applied as `A_ul -> beta*A_ul`
  *inside* the two-level solution, so the subcritical limit stays independent
  of beta. beta -> 1 in the thin wind; it matters only at a cold dense base.
  Resonance metal lines (Mg II h&k) are still treated as thin
- **The CHIANTI coronal fits are switched off below their 10^3 K validity
  floor** by a smooth cutoff (value and dT-slope continuous, bit-identical at
  and above 10^3 K), leaving only the explicit two-level fine-structure terms,
  which are valid there. Below the floor the coronal fits were dominated by
  exponentials whose excitation temperatures match no ground-term splitting of
  the ion. `cno_cool 0` (legacy AIOLOS branch) is deliberately not cut off
- He I 2³S metastable triplet state (coupled solver), with a
  **temperature-dependent He(2³S)+H Penning-ionization rate** (Taylor et
  al. 2025; replaces the classic 5e-10 constant)
- Updated photoionization data: **He I ground state = Verner et al. (1996)**
  by default (legacy ATES 2-term fit via `ATES_photoionization_rate: True`),
  and a **He I 2³S cross section extended past 60 eV** (two Verner-form
  wings + power-law bridge fitted to TOPbase/Opacity-Project data) — see
  `docs/photoion_cross_sections.pdf` and `docs/recombination_coefficients.pdf`
- **Diffusive separation of helium and metals** (opt-in `He_diffusion: True`):
  the He/H element ratio is transported with bulk advection + molecular-
  diffusion settling (Banks & Kockarts binary D, eddy `He_Kzz`, ambipolar-
  corrected settling mass, optional thermal diffusion `He_alphaT`), and each
  trace metal can diffuse independently (`He_metal_diffusion: True`) — He/H
  declines with altitude as in Taylor et al. (2025) / Xing et al. (2023);
  before/after impact on He 10830: `docs/version_compare.pdf`
- Non-LTE H(n=2) and Ly-alpha radiative transfer, computed either with a fast
  Neufeld escape-probability closure or by coupling to the **LaRT** Monte Carlo
  Ly-alpha radiative-transfer code, which sets the n=2 population from the full
  scattered field — including the Ly-alpha emitted in situ within the wind by
  recombination and collisional excitation
- A Jacobian-free Newton-Krylov (JFNK) steady-state solver with PTC warm-up,
  SER ramp, and non-monotone (Grippo) line search
- **Molecular transport** (opt-in `Viscosity: True` / `Conduction: True`): the
  Navier-Stokes viscous force, its dissipation, and heat conduction
  (`kappa(T) = 4.45e4 (T/1000 K)^0.7`, Watson et al. 1981; `mu(T)` tied to it
  by the monatomic Chapman-Enskog relation), integrated Crank-Nicolson and
  entering the steady residual with the same operator — see
  `docs/viscosity_conduction.md`, which also records that on the cases tested
  they change nothing measurable except at a cold molecular base
- Roche-potential geometry (spherical or Roche-lobe domain modes)
- **Lower-atmosphere connection** (opt-in): an analytic
  Koskinen+2022 lower column (`Lower column:` key reports the derived 1-ubar
  base radius and base H2/H/He), an EOS-only molecular-base correction
  (`Molecular base:`), **full molecular chemistry** (`Molecular chemistry:
  True` — H2/H2+/H3+/HeH+ in the coupled ionization equilibrium, Yan+1998 H2
  photoionization opacity/heating, Miller+2013 H3+ IR cooling; hot Jupiters
  develop a sharp H2->H front above a thin molecular base), and a `base.inp`
  handoff with two generators: `src/utils/run_lower.py` (analytic column,
  isothermal or Guillot T(p)) and `src/utils/vulcan_to_base.py` (converts a
  **VULCAN** photochemistry output — the photochemical H2/H state, which on
  HD 189733 b dissociates ~11x more H than equilibrium at 1 ubar) — see
  `docs/lower_atmosphere_coupling.pdf`
- **`EXHALE_transit.py`** transmission post-processor: transit spectra
  for He I 10830 Å, Ly-alpha, H-alpha, H-beta, and the metal resonance
  doublets Mg II h&k, Ca II H&K, and Na I D

For a complete description of the physics, solver, and all input parameters
see **`docs/EXHALE_user_manual.pdf`**; the full dated changelog (with a
code-size appendix vs. the original ATES) is **`docs/Update_EXHALE.pdf`**.
A task-oriented quick reference ("how do I run X?") is
**[`README_HOWTO.md`](README_HOWTO.md)**.

---

## Requirements

| Component | Version |
|-----------|---------|
| Fortran compiler | `gfortran` >= 9.3 or `ifort`/`ifx` >= 2021 |
| Python 3 | >= 3.8; packages: `numpy`, `scipy`, `matplotlib`, `tkinter`, `astropy` (`astropy` is used by `EXHALE_transit.py` for `astropy.convolution`) |
| MINPACK | included in `src/modules/nonlinear_system_solver/` |

---

## Installation

```bash
git clone https://github.com/seoncafe/EXHALE
cd EXHALE
```

No additional installation step is required.

---

## Directory layout

```
EXHALE/
├── Makefile               # incremental build (FC=gfortran default)
├── run_EXHALE.sh            # optional GUI launcher (writes input.inp, calls make)
├── src/
│   ├── EXHALE_main.f90    # program entry point
│   ├── modules/           # Fortran source modules (flux, init, radiation, …)
│   └── utils/             # Python GUI (EXHALE_interface_main.py), fortdep.py
├── VULCAN/                # third-party VULCAN (+FastChem), NOT in this repo;
│                          #   fetched by src/utils/setup_vulcan.sh — see below
├── inputdata/             # opacity / SED table samples (*.opa, Jlya.txt, …)
├── cooling_data/          # CHIANTI cooling-formula fit scripts + notebooks
├── examples/
│   ├── 01_legacy_marching/ … 15_molecular/  # ready-made input configs (15 folders)
│   │                          #   (solver stages, metals, He 2³S, Balmer/Lya,
│   │                          #    Wind-AE IC, lower atmosphere, He/metal
│   │                          #    diffusion, full molecular chemistry)
│   ├── README.md          # one-line description of each config folder
│   ├── exhale_io.py         # Python loaders for all output files
│   ├── EXHALE_analysis.ipynb
│   └── tutorial/          # minimal worked example (generic hot Jupiter)
├── docs/
│   ├── EXHALE_user_manual.pdf   # full reference manual
│   ├── Update_EXHALE.pdf      # dated changelog + code-size appendix vs ATES
│   ├── cooling_formulas.pdf   # analytic cooling-coefficient reference
│   ├── photoion_cross_sections.pdf   # H/He/He2³S cross-section benchmarks + TOPbase ext.
│   ├── recombination_coefficients.pdf # H/He recombination-data review
│   ├── design_hehe_diffusion.md      # He/H + metal diffusive-separation design/validation
│   ├── version_compare.pdf    # v1.0 vs current: diffusion impact on He 10830 (2 planets)
│   ├── methodology_aiolos_taylor_xing.pdf # AIOLOS/Taylor/Xing methodology comparison
│   ├── EXHALE_BC_and_IC.pdf     # boundary- and initial-condition reference
│   ├── code_comparison.pdf    # BC/IC/solver vs ATES, Salz, Kubyshkina, Murray-Clay
│   ├── steady_solver_memo.pdf # Newton-Krylov design notes
│   ├── newton_scaling_and_base_wall.md # JFNK diagonal scaling, line-search merit, watchdog
│   ├── viscosity_conduction.md # molecular viscosity + heat conduction: derivation, gates
│   ├── wind_ae_solver.pdf     # included Wind-AE solver (IC mode: windae)
│   ├── lower_atmosphere_coupling.pdf # lower-atmosphere connection: analytic column, molecular chemistry, VULCAN
│   ├── code_review_20260702.md # full-code review report (fixes + recommendations)
│   └── …
├── observational_data/    # digitized observational comparison data
├── EXHALE_transit.py                 # transmission spectrum post-processor
├── EXHALE_plots.py          # live / static profile plotter
├── roche_recon.py         # Roche-lobe geometry helper
└── eta_approx.py          # analytic heating-efficiency approximation
```

Build artifacts land in `build/`; `make clean` removes object/module files
while keeping `EXHALE.x`.

---

## Building

```bash
make                    # gfortran (default)
make FC=ifort           # ifort
make FC=ifx             # ifx
make clean              # remove build/ objects and modules (keeps EXHALE.x)
make distclean          # remove build/ and EXHALE.x
make wind_ae_ic         # optional: standalone Wind-AE IC generator (./wind_ae_ic.x)
```

---

## Running

### Option A — GUI (writes input.inp automatically)

```bash
chmod +x run_EXHALE.sh
./run_EXHALE.sh           # gfortran
./run_EXHALE.sh --ifort   # ifort
./run_EXHALE.sh --ifx     # ifx
```

The Tk interface opens, you fill in the planetary parameters, press **Done**,
and the code builds and runs.  The system preset list is stored in
`src/utils/params_table.txt` and can be extended with the **Add planet**
button.

The GUI writes the original ATES-format `input.inp`, which runs with the
legacy convergence behavior (single-stage marching at `du < 1e-3`, no
Newton finish).  All EXHALE extensions are opt-in keys appended to
`input.inp` (`du_th [PLM,WENO3]`, `Solver: Newton`, `Domain mode`, ...)
or separate runtime files (`metals.inp`, `opacity.inp`), so to use them
add the lines by hand or start from `examples/` (Option B).

### Option B — direct (recommended for scripted or repeated runs)

1. Pick or create a run directory containing an `input.inp` (copy one from
   `examples/` as a starting point; `examples/tutorial/` is a ready-to-run demo).
2. If metals are required, place a `metals.inp` in the same directory.
3. Build once at the repo root, then run the binary from inside the run directory:

```bash
make                                    # build ./EXHALE.x at the repo root
cd examples/tutorial && ../../EXHALE.x  # reads ./input.inp, writes ./output/
```

### Recommended convergence workflow

For robust convergence on typical hot-Jupiter / sub-Neptune problems, run the
two-stage scheme and **leave the Shapiro filter off**:

```
# in input.inp
Reconstruction scheme:    PLM+WENO3    # two-stage; PLM (or WENO3) alone is single-stage
du_th [PLM,WENO3]:        0.5 1.0e-3   # PLM until du<0.5, then WENO3 until du<1e-3
# Solver:  Newton                      # OPTIONAL extra residual-tightening finish
# Shapiro filter:  -1                  # OFF by default; opt-in only for breathing cases
```

**Convergence criterion (flux-based).** Convergence is judged on the *flux*
criterion of the reference codes: the fractional spread of the mass flux,
`dMdot/Mdot < du_th` (in the code `du` *is* the radial spread of `rho*v*r^2`).
This is the ATES test (Caldiroli 2021, `< 1e-3`) and is equivalent to the CETIMB
requirement (Koskinen 2013a) that `F_c = rho*v*r^2` be constant with altitude.
So the two lines above are the whole recipe: `Reconstruction scheme: PLM+WENO3`
enables the two-stage run, PLM switches to WENO3 at `du < 0.5`, and the run
**stops when `du < 1e-3`** (flux-converged). (With `Reconstruction scheme: PLM`
or `WENO3` the run is single-stage and only the first `du_th` value is used.)
The `Solver: Newton` key is **optional** and does *not* change this criterion.

**About `Solver: Newton` (optional).** `du` (a mass-flux flatness) and the steady
residual `||R|| = ||du/dt||` (mass+momentum+energy) are *different* quantities:
`du` can reach `1e-3` while `||R||` is still larger (an operator-split / energy
imbalance). When `Solver: Newton` is set, the WENO3 stage does **not** stop at
`du < 1e-3`; instead it warms up until the flux metric `du < newton_du_switch`
(default `1e-2`, the du-keyed hand-off) and then the JFNK solver drives the
**residual** `||R||` (not `du`) toward `1e-3` -- a stricter, separate check, not
a repeat of the `du` test. The hand-off point is adjustable via an optional
third token on the key: `Solver: Newton <du_switch>` (e.g. `Solver: Newton 5e-3`
hands off later; bare `Solver: Newton` keeps the `1e-2` default). `||R||` is
otherwise computed and reported
**for reference only** (volume-weighted by default, since the L-inf max is
dominated by the small near-base cells) and gates the stop only if you set
`Resid tol:`. A run that stops on `du` without reaching `||R|| < 1e-3` is still
flux-converged by the reference standard; it is not a failure, but the Mdot then
carries the `du`-stop path dependence of a few percent, so quote a
Newton-converged value where one is available. Full-physics cases
(He 2^3S + metals, secondary ionization, a `base.inp` handoff) that used to
abort at `info = 2` around `||R|| ~ 2.8e-3` now reach `info = 0`; the cause was
the solver's diagonal scaling and its stagnation watchdog, not the lower
boundary condition -- see `docs/newton_scaling_and_base_wall.md`,
`docs/EXHALE_BC_and_IC.pdf` (convergence-criterion and test-matrix sections)
and `docs/steady_solver_memo.pdf`.

### Base grid resolution (`Base grid [dr,cells]`)

With `Grid type: Mixed` the grid is `cells` uniform cells of size `dr` (in R_p)
stacked on the lower boundary, followed by `N - cells` stretched cells out to
`r_max` (`N = 500`, a compile-time constant). The default reproduces the
historical hardcoded grid:

```
# in input.inp -- optional; this line IS the default
Base grid [dr,cells]:  2.0e-4 50     # 50 cells of 2e-4 R_p = 0.01 R_p uniform
Base grid [dr,cells]:  5.0e-5 200    # the 4x refinement, same 0.01 R_p extent
```

The two numbers belong on one line because their product is the extent of the
uniform region: refine at fixed extent by dividing `dr` and multiplying `cells`
by the same factor. The key is ignored by `Grid type: Uniform` and `Stretched`.

**Why it matters.** `dr` has to resolve the base density scale height
`H = kT/(mu g)`. Where `H/dr` is only a few cells the discretization carries a
*stationary* 2-cell entropy (contact) mode that nothing damps: HLLC resolves a
zero-speed contact exactly, the gravity source is cell-local, and the WENO3
pressure gradient sees only interface pressures. Measured amplitudes are ~1e-2
in `ln rho` for `H/dr < 5` and ~1e-4 for `H/dr > 100`. High-gravity planets have
the least margin, and a cool shielded base shrinks `H` further.
`EXHALE_setup.out` echoes the value in effect as
`Base scale-height resolution: H(T_eq)/dr` and warns below 10 cells. Costs: the
CFL step scales with the smallest cell (a 4x finer base means ~4x more steps for
the same physical time), and since `N` is fixed the stretched region gives up
those cells and coarsens. Full investigation:
`docs/hd189_base_checkerboard.md`.

### Enabling metal chemistry

Place a `metals.inp` in the run directory listing the total elemental
abundances n_X/n_H by number (the equilibrium solver distributes each
element over its ionization stages):

```
# metals.inp — solar abundances (Asplund+2009)
C     2.69e-4
N     6.76e-5
O     4.90e-4
Mg    3.98e-5
Ca    2.19e-6
Na    1.74e-6
Fe    3.16e-5
```

Optional control keys: `pp_metals 0|1|2` (metal treatment in the advection
post-process), `cx_full 0|1` (full Huang+2023 charge-exchange network),
`cno_cool 0|1` (C/N/O cooling source: `1` = CHIANTI fits including
N I/N II, the default; `0` = legacy AIOLOS fits).  No recompile is
needed; remove `metals.inp` to run without metals.  A template with all
ten elements is in `inputdata/metals.inp.example`.

### Diffusive separation of He (and metals)

By default the He/H ratio is frozen at the input value at all radii.  To let
helium diffusively separate from hydrogen (settling vs. wind drag, so He/H
declines with altitude and the He 10830 line weakens on gentle escapers), add
to `input.inp`:

```
He_diffusion:        True     # transport He/H (advection + settling); default False
He_metal_diffusion:  True     # optional: each metal diffuses with its own mass/D
# He_Kzz:            1.0e9    # eddy-diffusion coefficient [cm^2/s]
# He_ambipolar:      False    # turn OFF the ambipolar settling correction (default on)
# He_alphaT:         0.15     # thermal-diffusion factor (default 0 = off)
```

All flags default off, so standard runs are unchanged.  Physics, numerics,
and validation: `docs/design_hehe_diffusion.md`; quantitative before/after
comparison on HD 209458 b and WASP-121 b: `docs/version_compare.pdf`.

### Molecular chemistry (H2, H2+, H3+, HeH+)

For warm Neptunes / sub-Neptunes (or to *verify* the atomic base of a hot
Jupiter), enable the molecular network:

```
Molecular chemistry: True   # coupled H2/H2+/H3+/HeH+ equilibrium + H3+ cooling
Molecular base:      True   # recommended companion (consistent base pressure)
```

`Ion_species*.txt` gains four columns (`H2 H2p H3p HeHp`). Requires He/H > 0;
v1 is exclusive with `metals.inp` and `He_diffusion` (the parser refuses the
combinations). Local-equilibrium caveats in
`docs/lower_atmosphere_coupling.pdf` §4.

### Lower-atmosphere pre-step: VULCAN as a subroutine

EXHALE can generate its own lower-boundary conditions before the wind solve.
One line in `input.inp` is enough:

```
Lower atmosphere: vulcan 1.36     # run VULCAN photochemistry (auto-fetched if absent)
#Lower atmosphere: analytic 1.36  # or: fast chemical-equilibrium column
#  (no key at all = classic base; the VULCAN step is fully optional)
```

The number is the 1-bar (transit) radius in R_J. On startup EXHALE invokes
`src/utils/vulcan_driver.py`, which (if `VULCAN/` is not yet present) fetches
it via `src/utils/setup_vulcan.sh`, copies the tree into `<run_dir>/vulcan_work/`, builds the planet's T(p)/Kzz atmosphere (Guillot
2010) and picks a stellar UV spectrum by the host Teff, compiles FastChem if
needed, runs VULCAN to steady state (**hours** on the first run; later runs
reuse the cached `.vul`), converts the result to `base.inp` (photochemical
H2/H/He state, base temperature and radius), and then proceeds with the wind
solve on that base — i.e. VULCAN acts as a subroutine of EXHALE. An existing
`base.inp` always wins (delete it to regenerate); `EXHALE_ROOT` overrides the
code-root path for relocated installs. Manual invocation and finer control:

```bash
python3 src/utils/vulcan_driver.py <run_dir> --r1bar 1.36 [--force] [--sflux F] [--atm F]
python3 src/utils/run_lower.py     <run_dir> --r1bar 1.36 [--guillot]
```

VULCAN provides H/C/N/O(/S) composition only — metal abundances stay in
`metals.inp`. The quick one-line consistency check without any handoff is
`Lower column: <transit radius>` in `input.inp`.

## Obtaining VULCAN and FastChem (third-party; not in this repo)

The lower-atmosphere pre-step uses the **VULCAN** photochemical-kinetics code,
which ships the **FastChem** equilibrium-chemistry code inside it. These are
third-party open-source codes and are **not** committed to the EXHALE
repository (`VULCAN/` is git-ignored). Fetch and prepare them with one command:

```bash
src/utils/setup_vulcan.sh          # clones VULCAN into EXHALE/VULCAN/ + patches + builds FastChem
```

(EXHALE also runs this automatically the first time you use `Lower atmosphere:
vulcan …` and `VULCAN/` is missing.) The script clones VULCAN, applies the two
modifications EXHALE needs, builds FastChem, and removes stale run products.

**Download sources / required citations** (please cite in any publication that
uses the pre-step):

- VULCAN — https://github.com/exoclime/VULCAN (mirror
  https://github.com/shami-EEG/VULCAN). Cite **Tsai et al. 2017, ApJS 228, 20**
  and **Tsai et al. 2021, ApJ 923, 264**.
- FastChem — https://github.com/NewStrangeWorlds/FastChem (shipped inside
  VULCAN). Cite **Stock et al. 2018, MNRAS 479, 865** and, for FastChem 2,
  **Stock et al. 2022, MNRAS 517, 4070**.

**What `setup_vulcan.sh` changes in the cloned VULCAN tree** (i.e. the only
modifications you need to make if you set it up by hand):

| File (in the VULCAN clone) | Modification |
|---|---|
| `make_chem_funs.py` | add `encoding=None` to the `np.genfromtxt(vulcan_cfg.com_file, …)` call so the element-conservation check does not crash on Python 3 (bytes-vs-str) |
| `vulcan_cfg.py` | set `use_photo = True` and `use_live_plot = False` (this file is the driver's template) |
| `output/*.vul` | delete (run products; regenerated per run in `<run_dir>/vulcan_work/`) |
| `fastchem_vulcan/` | `make` to build the `fastchem` binary |

**Added / modified on the EXHALE side** (these *are* in this repository):

| File | Change |
|---|---|
| `src/utils/setup_vulcan.sh` | **added**: fetch VULCAN + apply the above patches + build FastChem |
| `src/utils/vulcan_driver.py` | **added**: subroutine-style driver (planet cfg, cached VULCAN run, `base.inp`) |
| `src/utils/vulcan_to_base.py` | **added**: `.vul` → `base.inp` converter |
| `src/modules/files_IO/input_read.f90`, `src/modules/init/parameters.f90` | **modified**: the `Lower atmosphere: vulcan\|analytic <R_1bar>` key and the `run_lower_atm_prestep` invocation |

### Legacy atomic-data switch

The He I ground-state photoionization cross section defaults to the Verner
et al. (1996) fit; add `ATES_photoionization_rate: True` to revert to the
original ATES 2-term fit.  (The He 2³S cross section — TOPbase-extended
past 60 eV — and the temperature-dependent Penning rate are always on;
see `docs/photoion_cross_sections.pdf`.)

### Wind-AE warm-start initial condition (`IC mode: windae`)

For a planet that is hard to launch from the default cold/auto initial
conditions, EXHALE can build the initial condition from an included 1-D
steady-state Parker-wind solver (a Fortran port of **Wind-AE**; Murray-Clay
et al. 2009 / Broome et al. 2025, under `src/modules/wind_ae/`).  There are
two ways to use it.

**In-process (recommended)** — add one line to `input.inp`:

```
IC mode:                  windae
Load IC?                  False
Solver:                   Newton    # always Newton-finish for a quantitative Mdot
```

On startup EXHALE maps the planet's parameters to Wind-AE, ramps a shipped
seed solution to the planet, solves the steady wind, writes the IC onto the
EXHALE grid (`output/*_IC.txt`), and loads it automatically — the whole
input → Wind-AE solve → IC → run sequence is a single `./EXHALE.x` (or
`./run_EXHALE.sh`) invocation.  It is a *warm start*: EXHALE re-solves the
ionization (and any metals) from the first step, so the Wind-AE IC need not
be exactly self-consistent.  The shipped seed and spectrum live in
`inputdata/windae_seed.csv` and `inputdata/windae_spectrum.inp`.

To instead bootstrap from the **nearest** converged solution in the full
Broome et al. (2025a) grid (~1000 solutions; `pick_nearest_seed`), download the
grid database separately and place it under `inputdata/windae_grid/` — see
[`inputdata/README_windae_grid.md`](inputdata/README_windae_grid.md) for the
download links and the expected layout.

This works for hot Jupiters close to the shipped seed
(`examples/12_windae_ic_hd209/`, HD 209458 b) and, via the self-consistent-BC
continuation (re-converging the base boundary conditions, and turning the
molecular layer off when the base sinks into the wind), for strongly-bound,
far-from-seed planets too — including HD 189733 b (`examples/11_windae_ic/`),
whose Wind-AE ramp now converges and writes a valid IC.  (If a ramp ever
fails, the run prints a message advising `IC mode: auto`.)  Note that whether
EXHALE then *time-integrates* a given planet cleanly is a separate question
from the Wind-AE IC: HD 189733 b, for instance, hits an EXHALE-side
base-breathing instability near 1.07 R_p regardless of the IC source.

**Standalone generator** — the same solver also builds an IC out of process:

```bash
make wind_ae_ic     # builds ./wind_ae_ic.x (separate from EXHALE.x)
./wind_ae_ic.x <input.inp> <seed.csv> <spectrum.inp> \
               <IC_dump_grid> <outdir> [rhoscale]
```

It reads an EXHALE `input.inp`, ramps from the seed, and writes
`<outdir>/{Hydro_ioniz,Ion_species}_IC.txt` on the grid given by an
`IC_dump.txt` from a prior EXHALE run with the same domain settings.  Then
run EXHALE on the result with `Load IC?  True` and `Solver: Newton`.

The full algorithm, continuation ramp, and boundary/initial conditions are
documented in `docs/wind_ae_solver.pdf`.

---

## Output files

All output is written to `output/` in the run directory.

| File | Contents |
|------|----------|
| `Hydro_ioniz.txt` | Radius, number density, velocity, pressure, temperature, heating rate, cooling rate (columns vs. radius) |
| `Ion_species.txt` | Number densities of H I, H II, He I, He II, He III, He 2³S, and the metal ionization stages (33 species; zero columns when a species is off) |
| `Hydro_ioniz_adv.txt` | Post-processed version of `Hydro_ioniz.txt` (advection-corrected) |
| `Ion_species_adv.txt` | Post-processed version of `Ion_species.txt` |
| `Cooling_breakdown.txt` | Radiative cooling by channel vs. radius |
| `Excited_H.txt` | Non-LTE H(n=2) populations (when the Balmer/Ly-alpha physics is on) |

Every file starts with a `# columns ...` schema header, so analysis tools
adapt to the column layout automatically.

When **Load IC** is enabled the previous outputs are copied to `*_IC.txt`
and read back as initial conditions for a restart run.

Full column definitions are in `docs/EXHALE_user_manual.pdf` §4.

---

## Reading output in Python

```python
import sys
sys.path.insert(0, 'examples/')
import exhale_io

run = exhale_io.Run('output/')   # load all output files
print(run.r)        # radius [Rp]
print(run.T)        # temperature [K]
print(run.Mdot)     # mass-loss rate [g/s]
```

See `examples/exhale_io.py` for the full API and `examples/EXHALE_analysis.ipynb`
for a worked example.

---

## Live plot during a run

```bash
python3 EXHALE_plots.py          # plot current output (static)
python3 EXHALE_plots.py --live 4 # refresh every 4 s
```

The plotter reads the `# columns` headers, overlays the post-processed
`*_adv` profiles (dashed) when present, and adds a metal-ion-density
figure automatically for a metals-on run.

---

## Transmission spectra (`EXHALE_transit.py`)

After a converged run, compute the transit transmission spectra with:

```bash
MPLBACKEND=Agg python3 EXHALE_transit.py
```

`EXHALE_transit.py` reads `input.inp` and the
`*_adv.txt` profiles in `output/`, and
produces the model transmission curves (theoretical, instrument-convolved,
and instrument+rotation-convolved) — and, on request, figures as PNG +
vector PDF — for **He I 10830 Å,
Ly-alpha, H-alpha, H-beta** and the metal resonance doublets **Mg II h&k,
Ca II H&K, Na I D** (skipped automatically for a metals-off run).  The whole
system comes from the run directory's `input.inp`: the stellar radius and
effective temperature from `Stellar radius [R_sun]` / `Stellar Teff [K]`, and
the planet spin period as the tidally-locked orbital period built from
`Orbital distance` and `Parent star mass` + `Planet mass`.  Each is resolved as
`EXHALE_TRANSIT_*` environment override > `./input.inp` > built-in default, and
the resolved values and their sources are printed at startup.  Instrument
resolving powers and wavelength windows stay in the script's header block.
Each line is written into the run directory under its own key as
`tpm_<line>.txt` (`He10830`, `Lya`, `Halpha`, `Hbeta`, `MgII`, `CaII`, `NaI`);
the figures carry the same keys and are saved only when
`EXHALE_TRANSIT_FIG_PREFIX` is set.  A 3-D
Roche-equipotential geometry is available via `geometry = 'triaxial'`
(`roche_recon.py`).  Full description in `docs/transmission_spectrum.pdf`
and the manual's transmission-spectra section.

---

## Example configurations

Ready-made `input.inp` templates covering all physics/solver combinations
are in `examples/` (folders `01`–`11` are HD 189733 b so each option can be
isolated; `12`, `14`, `15` are HD 209458 b, where the Wind-AE warm start,
diffusion and molecular options are validated; `13` spans four planets).
See `examples/README.md` for the exact lines each one adds:

| Folder | Description |
|--------|-------------|
| `01_legacy_marching/` | PLM-only, single du threshold (simplest) |
| `02_two_stage/` | PLM → WENO3 two-stage |
| `03_newton/` | Two-stage + Newton finish (recommended default) |
| `04_newton_from_state/` | Resume from saved state with Newton |
| `05_metals/` | Metals on (solar C/N/O/Mg/Ca/Na/Fe) |
| `06_he23s/` | He I 2³S triplet included |
| `07_balmer_lya/` | Balmer + Ly-alpha RT |
| `08_full/` | Full physics (metals + He 2³S + Balmer/Lya) |
| `09_spherical/` | Extended spherical domain |
| `10_warm_seed_ic/` | Warm-seed initial condition |
| `11_windae_ic/` | In-process Wind-AE IC for HD 189733 b — the self-consistent-BC ramp converges and writes the IC; EXHALE's own base-breathing instability (separate from the IC) then limits the warm start |
| `12_windae_ic_hd209/` | In-process Wind-AE warm-start IC that works (HD 209458 b, seed-adjacent) |
| `13_lower_atmosphere/` | Lower-atmosphere connection for four planets: analytic 1-µbar base column + `base.inp` handoff (isothermal / Guillot T(p)) |
| `14_diffusion/` | Diffusive separation of He and metals (HD 209458 b): He/H declines with altitude, each metal settles independently, reshaping He 10830 |
| `15_molecular/` | Full molecular chemistry (HD 209458 b): H2/H2+/H3+/HeH+ coupled equilibrium; sharp H2→H front above a thin molecular base (metals/diffusion off) |

---

## References

1. Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., Rauscher, E.
   (2021). *Irradiation-driven escape of primordial planetary atmospheres I.
   The ATES photoionization hydrodynamics code.* A&A, 655, A30.

2. Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., Rauscher, E.
   (2022). *Irradiation-driven escape of primordial planetary atmospheres II.
   Evaporation efficiency of sub-Neptunes through hot Jupiters.* A&A, 663, A122.

3. Biassoni, F., Caldiroli, A., Gallo, E., Haardt, F., Spinelli, R., Borsa, F.
   (2024). *Self-Consistent Modeling of Metastable Helium Exoplanet Transits.*
   A&A, 682, A115.

4. Huang, C., Koskinen, T., Lavvas, P., Fossati, L. (2023). *A Hydrodynamic
   Study of the Escape of Metal Species and Excited Hydrogen from the
   Atmosphere of WASP-121b.* ApJ, 951, 123.

5. Taylor, A. R., Koskinen, T., et al. (2025). *A Multispecies Atmospheric
   Escape Model with Excited Hydrogen and Helium: Application to HD209458b.*
   ApJ, 989, 68.  (Temperature-dependent Penning rate; diffusive-separation
   reference model.)

6. Xing, L., Yan, D., Guo, J. (2023). *The Mass Fractionation of Helium in
   the Escaping Atmosphere of HD 209458b.* ApJ, 953, 166.  (Multi-fluid
   He/H fractionation reference.)

7. Verner, D. A., Ferland, G. J., Korista, K. T., Yakovlev, D. G. (1996).
   *Atomic Data for Astrophysics. II.* ApJ, 465, 487.  (Photoionization
   cross sections.)

---

## Author

Kwang-Il Seon (KASI / UST)

Last updated: 2026-08-11 19:00 KST
