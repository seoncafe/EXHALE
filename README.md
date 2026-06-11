# EXHALE

**EXoplanetary Hydrodynamic Atmospheric Loss and Escape** (EXHALE) is a
1-D radiation-hydrodynamics code that simulates photoionization-driven
atmospheric mass loss from irradiated exoplanets.  It is a heavily extended
fork of the ATES code (Caldiroli et al. 2021; Biassoni et al. 2024), adding:

- Trace metals (C, N, O, Fe, Mg, Ca, Na) solved self-consistently inside
  the MINPACK ionization system, with Badnell RR+DR recombination, Voronov
  collisional ionization, and Kingdon & Ferland charge-transfer with H
- Metal-line cooling from CHIANTI (Fe II, Mg II, Ca II, Na I, ...) and an
  optional AIOLOS/Black two-level fine-structure channel
- He I 2³S metastable triplet state (coupled solver)
- Ly-alpha radiative transfer via the Neufeld core-skipping escape-probability
  method
- A Jacobian-free Newton-Krylov (JFNK) steady-state solver with PTC warm-up,
  SER ramp, and non-monotone (Grippo) line search
- Roche-potential geometry (spherical or Roche-lobe domain modes)
- **TPM** (Transmission Probability Module) post-processor: transit spectra
  for He I 10830 Å, Ly-alpha 1215.67 Å, H-alpha 6562.8 Å, and H-beta 4861 Å

For a complete description of the physics, solver, and all input parameters
see **`docs/ATES_user_manual.pdf`**.

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
├── ATES_main.f90          # program entry point
├── Makefile               # incremental build (FC=gfortran default)
├── run_ATES.sh            # optional GUI launcher (writes input.inp, calls make)
├── src/
│   ├── modules/           # Fortran source modules (flux, init, radiation, …)
│   └── utils/             # Python GUI (ATES_interface_main.py), fortdep.py
├── inputdata/             # opacity / SED table samples (*.atesopa, Jlya.txt, …)
├── cooling_data/          # metal cooling tables (metal_cooling_chianti.txt, …)
├── examples/
│   ├── inputs/            # 10 ready-made HD 189733 b input configurations
│   ├── ates_io.py         # Python loaders for all output files
│   ├── ATES_analysis.ipynb
│   └── tutorial/          # minimal worked example (no metals)
├── docs/
│   ├── ATES_user_manual.pdf   # full reference manual
│   ├── steady_solver_memo.pdf # Newton-Krylov design notes
│   └── …
├── observational_data/    # digitized observational comparison data
├── TPM.py                 # transmission spectrum post-processor
├── ATES_plots.py          # live / static profile plotter
├── roche_recon.py         # Roche-lobe geometry helper
└── eta_approx.py          # analytic heating-efficiency approximation
```

Build artifacts land in `build/`; `make clean` removes object/module files
while keeping `ATES.x`.

---

## Building

```bash
make                    # gfortran (default)
make FC=ifort           # ifort
make FC=ifx             # ifx
make clean              # remove build/ objects and modules (keeps ATES.x)
make distclean          # remove build/ and ATES.x
```

---

## Running

### Option A — GUI (writes input.inp automatically)

```bash
chmod +x run_ATES.sh
./run_ATES.sh           # gfortran
./run_ATES.sh --ifort   # ifort
./run_ATES.sh --ifx     # ifx
```

The Tk interface opens, you fill in the planetary parameters, press **Done**,
and the code builds and runs.  The system preset list is stored in
`src/utils/params_table.txt` and can be extended with the **Add planet**
button.

### Option B — direct (recommended for scripted or repeated runs)

1. Edit `input.inp` (copy from `examples/inputs/` as a starting point).
2. If metals are required, place a `metals.inp` in the same directory.
3. Build and run:

```bash
make
./ATES.x
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

Place a `metals.inp` in the run directory listing the trace-metal abundances
(relative to solar):

```
# metals.inp — example
CI   1.0
NI   1.0
OI   1.0
```

No recompile is needed.  Remove `metals.inp` to run without metals.

---

## Output files

All output is written to `output/` in the run directory.

| File | Contents |
|------|----------|
| `Hydro_ioniz.txt` | Radius, density, velocity, pressure, temperature, heating rate, cooling rate, heating efficiency (columns vs. radius) |
| `Ion_species.txt` | Number densities of H I, H II, He I, He II, He III and (if metals active) metal ionization states |
| `Hydro_ioniz_adv.txt` | Post-processed version of `Hydro_ioniz.txt` (advection-corrected) |
| `Ion_species_adv.txt` | Post-processed version of `Ion_species.txt` |

When **Load IC** is enabled the previous outputs are copied to `*_IC.txt`
and read back as initial conditions for a restart run.

Full column definitions are in `docs/ATES_user_manual.pdf` §4.

---

## Reading output in Python

```python
import sys
sys.path.insert(0, 'examples/')
import ates_io

run = ates_io.Run('output/')   # load all output files
print(run.r)        # radius [Rp]
print(run.T)        # temperature [K]
print(run.Mdot)     # mass-loss rate [g/s]
```

See `examples/ates_io.py` for the full API and `examples/ATES_analysis.ipynb`
for a worked example.

---

## Live plot during a run

```bash
python3 ATES_plots.py          # plot current output (static)
python3 ATES_plots.py --live 4 # refresh every 4 s
```

---

## Transmission spectra (TPM)

After a converged run, compute the transit transmission spectrum with:

```bash
python3 TPM.py
```

TPM reads `input.inp` and the `*_adv.txt` profiles in `output/`, and produces
PNG figures for He I 10830 Å, Ly-alpha, H-alpha, and H-beta.  Stellar
parameters (`T_star`, `R_star`) must be set in `input.inp` for the Balmer
lines.  Full description in `docs/Halpha_transmission.pdf`.

---

## Example configurations

Ready-made `input.inp` templates for HD 189733 b covering all physics/solver
combinations are in `examples/inputs/`:

| Folder | Description |
|--------|-------------|
| `01_legacy_marching/` | PLM-only, single du threshold (simplest) |
| `02_two_stage/` | PLM → WENO3 two-stage |
| `03_newton/` | Two-stage + Newton finish (recommended default) |
| `04_newton_from_state/` | Resume from saved state with Newton |
| `05_metals/` | Metals on (C/N/O) |
| `06_he23s/` | He I 2³S triplet included |
| `07_balmer_lya/` | Balmer + Ly-alpha RT |
| `08_full/` | Full physics (metals + He 2³S + Balmer/Lya) |
| `09_spherical/` | Extended spherical domain |
| `10_warm_seed_ic/` | Warm-seed initial condition |

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
