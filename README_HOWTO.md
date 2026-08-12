# EXHALE HOWTO — task-oriented quick reference

One entry per task: the exact lines/commands, the expected output, and where
the full documentation lives. Everything here is opt-in; a bare `input.inp`
runs the legacy ATES-compatible model. (Reference manual:
`docs/EXHALE_user_manual.pdf`; changelog: `docs/Update_EXHALE.pdf`.)

## Run a standard converged model

```
# input.inp
Reconstruction scheme:    PLM+WENO3
du_th [PLM,WENO3]:        0.5 1.0e-3
Solver: Newton
```
```bash
make                                   # build ./EXHALE.x once at the repo root
cd examples/tutorial && ../../EXHALE.x # run in a self-contained run dir (input.inp + output/)
```
`Reconstruction scheme: PLM+WENO3` enables the two-stage PLM->WENO3 marching
(PLM to the first `du_th` value, WENO3 to the second); `Reconstruction scheme:
PLM` (or `WENO3`) alone is single-stage and uses only the first `du_th` value.
Then the JFNK Newton finish runs. If the flux metric
plateaus just above the hand-off threshold (seen with He diffusion), the
hand-off now fires on the plateau automatically. Quantitative Mdot always
needs the Newton finish. -> manual §2.5–2.6.

Check `EXHALE_setup.out` for `Base scale-height resolution: H(T_eq)/dr`. Below
about 10 cells (the code warns) the base carries a stationary cell-to-cell
entropy mode that the inviscid HLLC fluxes do not damp, and the converged
profile shows an alternating density/temperature pattern at constant pressure.
Refine at fixed extent — halve `dr`, double the cell count:

```
# input.inp   (Mixed grid only; default is 2.0e-4 50)
Base grid [dr,cells]:  1.0e-4 100
```

Of the four production planets only HD 189733 b needs this (27 cells on the
default grid, against 102 for WASP-121 b). The cells come out of the stretched
region, so a finer base is a coarser wind and a smaller CFL step; 2x is the
step that has been carried to convergence. -> `docs/hd189_base_checkerboard.md`.

## Damp the base with viscosity and heat conduction

```
# input.inp
Viscosity:  True     # radial viscous force + its dissipation q_mu
Conduction: True     # heat conduction, kappa(T) = 4.45e4 (T/1000 K)^0.7
```
Adds the Navier-Stokes molecular transport that CETIMB carries and the inviscid
HLLC scheme lacks, integrated Crank-Nicolson so the stiff base cells do not
limit the step. Both keys default off and a run without them is unchanged. On
the hot-Jupiter cases tested the terms change nothing measurable; on a cold
molecular base (HD 209458 b with a `base.inp` handoff) conduction moves base
`T` by up to 4% and `rho` by up to 13%.
`mu(T)` follows `kappa(T)` through the monatomic Chapman-Enskog relation
(Prandtl 2/3); `Viscosity: <mu0> [<s>]` instead sets a diagnostic power law
`mu = mu0*T^s` in code units. Coefficients are the neutral atomic-hydrogen
values — the ionized-wind (Spitzer) conductivity is not included. ->
`docs/viscosity_conduction.md`.

## Add trace metals

Put a `metals.inp` in the run directory (abundances n_X/n_H by number; template
`inputdata/metals.inp.example`). Remove the file to turn metals off — no
rebuild. Optional keys inside: `pp_metals`, `cx_full`, `cno_cool`, `eos_metals`.
-> manual §3.2, README "Enabling metal chemistry".

## He I 10830 / metastable triplet

```
Include He23S? True
```
`EXHALE_transit.py` computes the 10830 line from the `*_adv` output:
```bash
MPLBACKEND=Agg python3 EXHALE_transit.py          # He 10830, Lya, Halpha, Hbeta, metal lines
# wide window when the line is broad:
# widen the He window (auto by default; override with EXHALE_TRANSIT_HE_LMIN/LMAX/N):
EXHALE_TRANSIT_HE_LMIN=10827.5 EXHALE_TRANSIT_HE_LMAX=10832.5 python3 EXHALE_transit.py
```
-> manual §5.2.

### Where the star and planet parameters come from

`EXHALE_transit.py` takes the whole system from the `input.inp` of the run
directory, so the spectra always describe the same system as the simulation.
Each parameter is resolved as

```
EXHALE_TRANSIT_* environment override  >  ./input.inp  >  built-in default
```

| Quantity | `input.inp` label | Override |
|----------|-------------------|----------|
| `R_p`, `M_p`, `T_eq`, `a`, `M_star` | `Planet radius`, `Planet mass`, `Equilibrium temperature`, `Orbital distance`, `Parent star mass` | - |
| `L_EUV`, day-night `xi` | `Log10 of EUV luminosity`, `2D approximate method` | `xi_override` in the script |
| `R_star` | `Stellar radius [R_sun]` | `EXHALE_TRANSIT_RSTAR_RSUN` |
| `T_star` | `Stellar Teff [K]` | `EXHALE_TRANSIT_TSTAR` |
| planet spin period | none - computed as `P_orb` from `a` and `M_star + M_p` (tidal locking) | `EXHALE_TRANSIT_ROTP` |

The script prints every resolved parameter and its source at startup, and warns
when `Stellar radius`/`Stellar Teff` are missing from `input.inp`. Read that
block first when a spectrum looks wrong: `R_star` sets the transit
normalization, so a wrong stellar radius rescales every absorption depth.
Instrument resolving powers, wavelength windows, and geometry stay
script/environment settings - they describe the observation, not the system.

### Where the spectra go

Every line carries one key - `He10830`, `Lya`, `Halpha`, `Hbeta`, `MgII`,
`CaII`, `NaI` - and both products of a line are named from it, in the run
directory:

| Product | Name | Default |
|---------|------|---------|
| model curve | `<path>/tpm_<line>.txt` | always written |
| figure | `<path>/<prefix><line>.png` (and `.pdf`) | off; set `EXHALE_TRANSIT_FIG_PREFIX` |

```bash
MPLBACKEND=Agg EXHALE_TRANSIT_PATH=WASP-121b \
  EXHALE_TRANSIT_FIG_PREFIX=tpm_ python3 EXHALE_transit.py
```

The curves are written for whatever the run contains, so a metals-off run
simply has no `tpm_MgII/CaII/NaI.txt`. Each file has columns
`lambda[A]  T_theo  T_instr  T_rot+instr`; excess absorption in percent is
`(1 - T)*100`. `EXHALE_TRANSIT_SAVE_PREFIX` decorates the curve name
(`<path>/<prefix>tpm_<line>.txt`); both prefixes are resolved against the run
directory, so give an absolute prefix to write somewhere else.

### Instrument resolution and rotation (env overrides)

The instrument resolving power for each line `R = lambda/Delta-lambda` and the azimuthal
sampling of the rotation integral are run-time overridable:

```bash
# resolving power per line (defaults in parentheses):
#   RES_HETR (8e4, He 10830)   RES_HI  (5e4, Lya)    RES_HA (1.15e5, Halpha)
#   RES_HB   (1.15e5, Hbeta)   RES_MGII(3e4, Mg II)  RES_CAII/RES_NAI (=RES_HA)
EXHALE_TRANSIT_RES_HETR=5e4 EXHALE_TRANSIT_RES_HI=1.14e5 python3 EXHALE_transit.py

# planet rotation period [days] and azimuthal samples of the exact disk integral.
# ROTP overrides the tidally-locked default (the orbital period built from
# `Orbital distance` and `Parent star mass` + `Planet mass`); use it only for a
# planet that is not tidally locked, or to test the sensitivity to the spin.
EXHALE_TRANSIT_ROTP=2.2185 EXHALE_TRANSIT_ROT_NPHI=64 python3 EXHALE_transit.py
```

Rotation is computed as the exact projected-disk integral (each chord Doppler-shifted by
its local solid-body velocity `Omega*b*cos(phi)` and averaged over azimuth), not a
single-velocity Gaussian convolution; it conserves each line's equivalent width.
-> manual §5.2.

## He/H (and metal) diffusive separation

```
He_diffusion: True          # He element transported (advection + settling)
He_metal_diffusion: True    # optional: each metal with its own mass/D
# He_Kzz: 1.0e9   He_ambipolar: True   He_alphaT: 0.0
```
Default off. With `Solver: Newton` the code co-converges the diffused He/H
field with the steady wind (outer JFNK<->diffusion iteration). Physics and
two-planet impact: `docs/design_hehe_diffusion.md`,
`docs/version_compare.pdf`. -> manual §3.6.

## Check the base radius (analytic lower column)

```
Lower column: 1.36     # the 1-bar (transit) radius [R_J]
```
Startup report: derived r(1 ubar) as an [equilibrium, fully-atomic] bracket,
base H2/H/He fractions and mu, next to the input "Planet radius". Rule of
thumb from the four worked planets (`examples/13_lower_atmosphere/`):
HD 189733 b input radius falls inside its bracket; WASP-52 b's transit-radius
shortcut is 0.17–0.35 R_J too deep (biases Mdot ×1.5, He 10830 unchanged).
-> `docs/lower_atmosphere_coupling.pdf` §4.1, §5.

## Molecular chemistry (warm Neptunes / sub-Neptunes)

```
Molecular chemistry: True
Molecular base:      True
```
H2/H2+/H3+/HeH+ join the coupled ionization equilibrium; H2 photoionization
opacity/heating (Yan+1998) and H3+ IR cooling (Miller+2013) are included;
`Ion_species*.txt` gains `H2 H2p H3p HeHp` columns. Requires He; `He_diffusion`
is refused. A `metals.inp` may be present -- the metal stages are then solved in
the same system as the molecular network, which they share the free electron
density with (`examples/16_molecular_metals`). Hot Jupiters: thin molecular base, sharp
H2->H front, atomic wind above (the atomic assumption becomes a result).
Caveat: local equilibrium (no molecular advection).
-> `docs/lower_atmosphere_coupling.pdf` §4.3.

## Hand off a lower-atmosphere model (`base.inp`)

```bash
# analytic column (isothermal Teq; --guillot for the semi-grey T(p)):
python3 src/utils/run_lower.py <run_dir> --r1bar 1.36
```
writes `<run_dir>/base.inp` with `T_base r_base HeH_base Kzz_base`; EXHALE
reads it at startup and echoes every override (absent file = strict no-op).
The VULCAN converter (`src/utils/vulcan_to_base.py`, below) adds `q_H2_base`
and `p_base`: with `Molecular base: True` the photochemical H2 mixing ratio
then replaces EXHALE's chemical-equilibrium fit in the base particle count.
-> `docs/lower_atmosphere_coupling.pdf` §4.4.

## Use VULCAN photochemistry for the base state (subroutine-style)

```
# input.inp — one line; EXHALE runs VULCAN itself on startup (auto-fetched):
Lower atmosphere: vulcan 1.138      # arg = 1-bar (transit) radius [R_J]
```
First run takes hours (VULCAN to steady state in `<run_dir>/vulcan_work/`);
later runs reuse the cached `.vul`. `Lower atmosphere: analytic 1.138` uses
the fast equilibrium column instead; omit the key to skip the pre-step
entirely. Manual control: `python3 src/utils/vulcan_driver.py <run_dir>
--r1bar 1.138 [--force]`. The photochemical base differs from equilibrium
(HD 189733 b: q_H = 0.23 at 1 ubar vs 0.020 — 11x more dissociation).
VULCAN is H/C/N/O(/S) only: metal abundances stay in `metals.inp`.
VULCAN+FastChem are third-party codes fetched by
`src/utils/setup_vulcan.sh` into `VULCAN/` (not committed; EXHALE also
fetches them automatically on first use) — see README "Obtaining VULCAN and
FastChem" for the download URLs and required citations.
-> `docs/lower_atmosphere_coupling.pdf` §4.4 and Fig. 2.

## Warm-start a hard planet (Wind-AE IC)

```
IC mode:  windae
Load IC?  False
Solver:   Newton
```
-> manual §2.7, `docs/wind_ae_solver.pdf`.

## Post-process into transmission spectra

`MPLBACKEND=Agg python3 EXHALE_transit.py` in the run directory (reads `input.inp` +
`output/*_adv.txt`). Metal doublets (Mg II, Ca II, Na I D) appear
automatically for metals-on runs; `geometry='triaxial'` enables the
Roche-equipotential geometry (`roche_recon.py`). -> manual §5.2,
`docs/transmission_spectrum.pdf`.

## Regression / hygiene

```bash
make check                              # golden regression (= backup/regression/run_check.sh check)
./backup/regression/run_fcheck.sh       # runtime-checked build (bounds/mem), periodic
./backup/regression/test_roundtrip.sh   # restart round-trip
```
The harness lives in `backup/regression/`, which is a working-copy directory
and is not in the git remote. Default matrix (four cases, each bitwise against
its golden `Hydro_ioniz.txt` / `Ion_species.txt`):

| case | what it guards |
|---|---|
| `wasp_full` | WASP-121 b with He 2³S **and** metals |
| `wasp_he23off` | the same with He 2³S off (the HeITR-off branch) |
| `mol_base_handoff` | hot-Uranus Tier-2 gate: molecular chemistry + a `base.inp` handoff whose `q_H2_base` drives the photochemical base particle count; 12000-step snapshot |
| `mol_metals` | the same gate + solar C/N/O/Mg/Ca/Na/Fe: the molecular and metal networks in one system; 12000-step snapshot |

Any other case directory can be named on the command line. `run_check.sh golden`
does **not** re-run — it snapshots whatever `output/` sits in each case
directory, so the order after a code change is `check`, `golden`, `check`.
`run_fcheck.sh` rebuilds with `-fcheck`, runs a bounded HD 209458 b case,
fails on any runtime trap, then restores the production build (it caught a
real out-of-bounds read on WASP-121b in the 2026-07-02 review).

## Where results/documents live

- `docs/EXHALE_user_manual.pdf` — full reference (inputs, outputs, physics)
- `docs/Update_EXHALE.pdf` — dated changelog + code-size appendix vs ATES
- `docs/lower_atmosphere_coupling.pdf` — lower-atmosphere connection: survey,
  implementation, 4-planet examples, figures
- `docs/newton_scaling_and_base_wall.md` — JFNK diagonal scaling, line-search
  merit and stagnation watchdog; why the base momentum row is not the blocker
- `examples/` — ready-made configs 01–16; the planet directories `HD209458b/`,
  `HD189733b/`, `WASP-121b/`, `WASP-52b/` sit at the repo root and are
  self-contained (own `input.inp`, output, notebooks)
