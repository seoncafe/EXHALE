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
make && ./EXHALE.x
```
`Reconstruction scheme: PLM+WENO3` enables the two-stage PLM->WENO3 marching
(PLM to the first `du_th` value, WENO3 to the second); `Reconstruction scheme:
PLM` (or `WENO3`) alone is single-stage and uses only the first `du_th` value.
Then the JFNK Newton finish runs. If the flux metric
plateaus just above the hand-off threshold (seen with He diffusion), the
hand-off now fires on the plateau automatically. Quantitative Mdot always
needs the Newton finish. -> manual §2.5–2.6.

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

### Instrument resolution and rotation (env overrides)

The per-line instrument resolving power `R = lambda/Delta-lambda` and the azimuthal
sampling of the rotation integral are run-time overridable:

```bash
# resolving power per line (defaults in parentheses):
#   RES_HETR (8e4, He 10830)   RES_HI  (5e4, Lya)    RES_HA (1.15e5, Halpha)
#   RES_HB   (1.15e5, Hbeta)   RES_MGII(3e4, Mg II)  RES_CAII/RES_NAI (=RES_HA)
EXHALE_TRANSIT_RES_HETR=5e4 EXHALE_TRANSIT_RES_HI=1.14e5 python3 EXHALE_transit.py

# planet rotation period [days] and azimuthal samples of the exact disk integral:
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

## Check the base radius (Tier-1 lower column)

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
`Ion_species*.txt` gains `H2 H2p H3p HeHp` columns. Requires He; v1 refuses
`metals.inp` and `He_diffusion`. Hot Jupiters: thin molecular base, sharp
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
./regression/run_check.sh        # golden regression
./regression/run_fcheck.sh       # runtime-checked build (bounds/mem), periodic
./regression/test_roundtrip.sh   # restart round-trip
```
`run_fcheck.sh` rebuilds with `-fcheck`, runs a bounded HD 209458 b case,
fails on any runtime trap, then restores the production build (it caught a
real out-of-bounds read on WASP-121b in the 2026-07-02 review).

## Where results/documents live

- `docs/EXHALE_user_manual.pdf` — full reference (inputs, outputs, physics)
- `docs/Update_EXHALE.pdf` — dated changelog + code-size appendix vs ATES
- `docs/lower_atmosphere_coupling.pdf` — lower-atmosphere connection: survey,
  implementation, 4-planet examples, figures
- `examples/` — ready-made configs 01–13 (+ per-planet folders `HD209458b/`,
  `HD189733b/`, `WASP-121b/`, `WASP-52b/` with analysis notebooks)
