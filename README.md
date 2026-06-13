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
- He I 2³S metastable triplet state (coupled solver)
- Non-LTE H(n=2) and Ly-alpha radiative transfer via the Neufeld
  escape-probability method
- A Jacobian-free Newton-Krylov (JFNK) steady-state solver with PTC warm-up,
  SER ramp, and non-monotone (Grippo) line search
- Roche-potential geometry (spherical or Roche-lobe domain modes)
- **TPM** (Transmission Probability Module) post-processor: transit spectra
  for He I 10830 Å, Ly-alpha, H-alpha, H-beta, and the metal resonance
  doublets Mg II h&k, Ca II H&K, and Na I D

For a complete description of the physics, solver, and all input parameters
see **`docs/EXHALE_user_manual.pdf`**.

---

## Requirements

| Component | Version |
|-----------|---------|
| Fortran compiler | `gfortran` >= 9.3 or `ifort`/`ifx` >= 2021 |
| Python 3 | >= 3.8; packages: `numpy`, `scipy`, `matplotlib`, `tkinter` |
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
├── EXHALE_main.f90          # program entry point
├── Makefile               # incremental build (FC=gfortran default)
├── run_EXHALE.sh            # optional GUI launcher (writes input.inp, calls make)
├── src/
│   ├── modules/           # Fortran source modules (flux, init, radiation, …)
│   └── utils/             # Python GUI (EXHALE_interface_main.py), fortdep.py
├── inputdata/             # opacity / SED table samples (*.atesopa, Jlya.txt, …)
├── cooling_data/          # CHIANTI cooling-formula fit scripts + notebooks
├── examples/
│   ├── 01_legacy_marching/ … 12_windae_ic_hd209/  # ready-made input configs
│   ├── README.md          # one-line description of each config folder
│   ├── exhale_io.py         # Python loaders for all output files
│   ├── EXHALE_analysis.ipynb
│   └── tutorial/          # minimal worked example (generic hot Jupiter)
├── docs/
│   ├── EXHALE_user_manual.pdf   # full reference manual
│   ├── cooling_formulas.pdf   # analytic cooling-coefficient reference
│   ├── EXHALE_BC_and_IC.pdf     # boundary- and initial-condition reference
│   ├── code_comparison.pdf    # BC/IC/solver vs ATES, Salz, Kubyshkina, Murray-Clay
│   ├── steady_solver_memo.pdf # Newton-Krylov design notes
│   ├── wind_ae_solver.pdf     # bundled Wind-AE solver (IC mode: windae)
│   └── …
├── observational_data/    # digitized observational comparison data
├── TPM.py                 # transmission spectrum post-processor
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

1. Edit `input.inp` (copy from `examples/` as a starting point).
2. If metals are required, place a `metals.inp` in the same directory.
3. Build and run:

```bash
make
./EXHALE.x
```

### Recommended convergence workflow

For robust convergence on typical hot-Jupiter / sub-Neptune problems:

```
# in input.inp
Reconstruction:           PLM
du_th [PLM,WENO3]:        0.5 1.0e-3   # two-stage: coarse PLM then fine WENO3
Solver:                   Newton        # JFNK finish once ||R||_inf < 0.05
```

The two-stage key `du_th [PLM,WENO3]` replaces the old single-value `du_th`.
EXHALE switches from PLM to WENO3 automatically when du falls below the first
threshold, then hands off to the Newton-Krylov solver near the fixed point.
See `docs/steady_solver_memo.pdf` for details.

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

### Wind-AE warm-start initial condition (`IC mode: windae`)

For a planet that is hard to launch from the default cold/auto initial
conditions, EXHALE can build the initial condition from a bundled 1-D
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

This works best for hot Jupiters close to the shipped seed (see the working
example `examples/12_windae_ic_hd209/`, HD 209458 b).  A strongly-bound,
far-from-seed planet can stall the static-BC ramp — the run then prints a
message advising `IC mode: auto` (see the deliberately non-converging
`examples/11_windae_ic/`, HD 189733 b).

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
| `Cooling_breakdown.txt` | Per-channel radiative cooling vs. radius |
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

## Transmission spectra (TPM)

After a converged run, compute the transit transmission spectra with:

```bash
MPLBACKEND=Agg python3 TPM.py
```

TPM reads `input.inp` and the `*_adv.txt` profiles in `output/`, and
produces spectrum figures (PNG + vector PDF; theoretical, instrument-
convolved, and instrument+rotation-convolved curves) for **He I 10830 Å,
Ly-alpha, H-alpha, H-beta** and the metal resonance doublets **Mg II h&k,
Ca II H&K, Na I D** (skipped automatically for a metals-off run).  Stellar
parameters (`R_star`, `rot_period`, `T_star` for the Balmer lines) and the
output figure names are set in the script's header block.  A 3-D
Roche-equipotential geometry is available via `geometry = 'triaxial'`
(`roche_recon.py`).  Full description in `docs/transmission_spectrum.pdf`
and the manual's TPM section.

---

## Example configurations

Ready-made `input.inp` templates covering all physics/solver combinations
are in `examples/` (folders `01`–`11` are HD 189733 b; `12` is HD 209458 b).
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
| `11_windae_ic/` | In-process Wind-AE warm-start IC — deliberate non-converging case (HD 189733 b stalls the static-BC ramp; use `IC mode: auto`) |
| `12_windae_ic_hd209/` | In-process Wind-AE warm-start IC that works (HD 209458 b, seed-adjacent) |

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
