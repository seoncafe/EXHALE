# p-winds oracle for LHS 1140b — build and validation record (Phase A3)

Built 2026-08-22; rebuilt the same day on the **v23** SEDs (see
`../sed/README.md` — the paper-era Mega-MUSCLES release recovered from
local holdings).  Driver: `run_pwinds.py` (frozen; edit only with a
recorded reason).  p-winds: the workspace clone `../../../p-winds` at
commit `bd01d2f` (2024-08-16).  Inputs: `../system_parameters.md`,
`../sed/`.  Metrics: `../../he_line_metrics.py` (repo root; the paper's
three-Gaussian definition; selftest PASS).

## Cases and formal metrics (instrument-convolved, R = 68,000)

| Case | SED | T [K] | Mdot [g/s] | red depth | blue depth | red/blue | FWHM [A] |
|---|---|---|---|---|---|---|---|
| published_gj1132 | gj1132 v23 | 5160 | 2.03e8 | 2.27% | 1.00% | 2.3 | 0.46 |
| **matched_gj1132** | gj1132 v23 | 6100 | 2.03e8 | **1.27%** | **0.44%** | **2.9** | 0.47 |
| published_gj699 | gj699 v23 | 5160 | 2.03e8 | 3.45% | 1.72% | 2.0 | 0.49 |
| observed 2024 | — | — | — | 1.24 +0.22/-0.23% | 0.25 +0.14/-0.12% | 6.7 +12.7/-3.1 | 0.86 +0.15/-0.27 |

At the paper's exact best fit (T = 5160 K) the forward model overshoots the
red depth 1.8x (GJ 1132 SED).  Moving T to 6100 K — inside the paper's own
explored range (they span 5160-6400 K) — reproduces the observation: red
0.1 sigma, blue ~1.3 sigma, ratio ~1.2 sigma.  The model FWHM is ~1.5 sigma
narrower than observed in every case; the measured 23.9 km/s exceeds the
thermal+LSF width and our model carries no bulk line-of-sight wind, which
matches the paper's own discussion of the width.

## Residual difference from the published fit, and what was ruled out

With the v23 SEDs the input provenance is closed; the remaining 1.8x at
T = 5160 K belongs to unpublished configuration details of the paper's
p-winds setup (their exact revision, grids, boundary fractions,
convolution) — the underdetermination the plan's review predicted.
Sensitivity levers measured along the way (v25-era scans, retained in this
directory; the conclusions carry over):

- **T is the strong axis**: base density is exponentially sensitive;
  red depth 2.12% -> 1.20% across 5160 -> 6100 K at fixed Mdot (v23 SED).
- **Mdot is weak** (`scan_mdot_fxuv33.txt`, v25-era): the lines are
  saturated near the published point.
- **FUV destruction band is a null lever**: tripling the 911-2600 A flux
  moves n_He3 by 0.2% — the triplet sink is radiative decay
  (A31 = 1.27e-4 /s), not photoionization.
- **The SED version moves the EUV production band 4.2x** (v23 vs v25,
  `../sed/README.md`), worth ~0.5% in red depth at 5160 K.

## What the oracle is for

Downstream acceptance (Phases C/F) is **EXHALE vs p-winds with identical
inputs** — same SED file, same system parameters, same (Mdot, T, H:He).
`matched_gj1132` is the reference case for comparisons against the observed
2024 line; `published_gj1132` is the reference at the paper's parameter
point; `published_gj699` is the SED-sensitivity alternate.

## Files

- `run_pwinds.py` — the frozen driver (three cases above)
- `tspec_<case>.txt` — vacuum wavelength [A], transmission, excess [%],
  convolved excess [%]
- `profile_<case>.txt` — r [R_p], v [km/s], rho [g/cm3], f_HII, f_He3,
  n_He3 [cm^-3]
- `scan_mdot.py`, `scan_mdot_fxuv33.txt` — Mdot scan (v25-era record; the
  fxuv33 SED it references was removed with the v23 rebuild)
- `run_pwinds.log` — console record of the production run

Model configuration: isothermal Parker wind with self-consistent mu from
`hydrogen.ion_fraction` (relax_solution, mu ~ 4.0); He populations with
all-singlet initial state; transit geometry b = 0.658 (i = 89.96 deg),
grid 121 px supersampled 10x; NIST triplet properties from
`p_winds.lines`; air->vacuum by the line-0 ratio; Gaussian LSF R = 68,000.

## Broadening-treatment test (2026-08-24)

Question: the model line at the paper's retrieval point is much narrower than
the measured one (FWHM 0.47 vs 0.84 A) -- is that an option of the radiative
transfer? Measured with `broadening_test.py` (matched_gj1132 physics, only
the `radiative_transfer_2d` call changed):

| treatment | red | blue | red/blue | FWHM |
|---|---|---|---|---|
| `average` (the oracle's) | 1.265 % | 0.437 % | 2.90 | 0.465 A |
| `formal` | 1.265 % | 0.436 % | 2.90 | 0.466 A |
| `formal` + turbulence | 1.265 % | 0.436 % | 2.90 | 0.466 A |
| observed (released data) | 1.254 % | 0.198 % | 6.35 | 0.841 A |

No: the three treatments are indistinguishable here. The width and the
saturated red/blue ratio are intrinsic to the Parker model at these
parameters.

Note on the paper's figures: the curves drawn over the spectrum in Fig. 3
are the **three-Gaussian line fits** (their caption: "best-fitting
absorption line models" with fitted depth and Doppler shift; the purple
family is MCMC draws of that fit), not p-winds. The p-winds physical model
is overlaid on the spectrum once, in **Fig. 4** -- the violet "2024 best
fit (T_wind = 5160 K)" curve, which the caption miscalls "dark blue". That
overlay shows the same character as this oracle's reproduction: its blue
component reaches about 0.6 % where the data show 0.1-0.2 % (the saturated
red/blue ratio), and its red feature is slightly narrower than the measured
one. Fig. 3's apparently perfect model-data agreement is the Gaussian fit,
not the physical model.

## 2026-08-24: the EW test that identified the paper's SED normalization

At the Fig. S5 posterior medians (log Mdot = 8.31, log T = 3.7127,
log H:He = -2.99, v_wind = 2260 m/s -- identical to the main-text values),
measured in one aperture (air frame, red-trough 10829.6-10831.2 A) against
the paper's Fig. 4 best-fit curve digitized from the published figure
(depth 1.323 %, FWHM 0.708 A, EW 0.912 %.A):

- catalog-A SED: depth 2.150, FWHM 0.571, EW 1.128 (+24 %) -- no kernel and
  no parameter vector inside the 1-sigma box reaches the plotted EW;
- integral-matched SED: depth 1.817, FWHM 0.511, EW 0.906 (+0.7 %) -- a
  ~13 km/s FWHM Gaussian kernel then maps depth and width onto the plotted
  curve simultaneously.

Conclusion: the retrieval consumed the integral-matched SED. The ~13 km/s
kernel between the computed and plotted model remains undocumented in the
paper (candidates: model evaluated on a coarse wavelength grid, or an LSF
applied to the model beyond the R = 68,000 already included); it conserves
EW, so it does not bias the retrieval of a saturated line.

## 2026-08-24 (authors' code obtained): the "kernel" identified

The paper's own p-winds script (LHS1140b_zenodo/lhs1140b_2024B_pwinds.py,
Zenodo 20723095) settles the reconciliation. Reproduced verbatim
(`reproduce_authors_fig4.py`):

- SED normalization is const-res v23 x Lx-ratio **0.59** (catalog ratio) --
  NOT the integral-matched 0.35 this repo switched to. The EW test had
  pointed at integral matching only because our oracle was missing the two
  broadening switches below, which the EW test could not see. Corrected: the
  paper uses 0.59.
- the extra width we modeled as a ~13 km/s Gaussian kernel is
  **turbulence_broadening=True** (Lampon+2020) plus **5-phase transit
  averaging**, both explicit in their `transmission_model`. No undocumented
  kernel exists.

The reproduction is now run at the **MCMC medians** (Mdot 2.03e8 g/s,
T 5160 K, H:He 1.01e-3), which is the vector the paper's Fig. 4 legend states
for the plotted curve. With the authors' exact settings this gives red
1.376 %, blue 0.534 %, FWHM 0.575 A, against the paper's Fig. 4 curve
(1.323 %, ~0.71 A) -- a 3 % match in depth and a line ~20 % too narrow.

The width residual was previously attributed to the plotted curve being drawn
at the medians rather than at the grid argmax. That was tested and is wrong:
the grid argmax (1.83e8, 5131.6, 8.338e-4) gives red 1.366 %, FWHM 0.581 A,
so the two vectors differ by less than 0.01 A in width. What produces the
extra width of the published curve has not been identified.

Bulk-velocity convention: `tspec_authors_*.txt` are written **without** any
line-of-sight bulk velocity applied inside the model; v_wind = 2.26 km/s is
applied at plot time, as it is for every other model file here
(`tspec_published_*`, the EXHALE `tpm_*` outputs). This is exact for
`wind_broadening_method='average'`, where `bulk_los_velocity` is a rigid
shift: with and without it, red/blue depths and FWHM agree to the third
decimal and only the fitted shift moves (+0.080 A vs -0.002 A). Earlier
versions of these files carried the shift internally *and* had it added again
at plot time, so Figures 1(a) and 4 of the memo were drawn ~0.08 A too red.

Consequence for our own oracle: `run_pwinds.py` used turbulence OFF and a
single phase, so its lines were too narrow -- the model comparisons in the
notebook should use the authors' two switches. The SED should revert to the
0.59 (catalog) normalization.
