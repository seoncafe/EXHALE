# Observational comparison data for the EXHALE model planets

Working catalog of the measured atmosphere data available under
`~/Exoplanetary_Atmosphere/` and inside this workspace, assembled for the paper
(Task B) comparison of EXHALE transmission spectra against observations.

EXHALE transit spectra come from `EXHALE_transit.py` as `tpm_*.txt`
(columns: wavelength [Angstrom] vs transit absorption, some as `R_eff/R_star`),
covering He I 10830, Ly-alpha 1215.67, H-alpha 6562.8, H-beta 4861, and the
metal doublets Mg II h&k, Ca II H&K, Na I D. **Every measured file below uses a
slightly different y-axis convention (fractional excess absorption, percent, or
normalized in/out flux) — none is `R_eff/R_star`. Normalize the EXHALE output to
each file's own convention before overplotting.**

Primary model planets: **HD209458b, HD189733b, WASP-52b, WASP-121b**.
Secondary/validation: GJ1214b, GJ436, HAT-P-11b, WASP-107b (+ off-list WASP-69b).

Root paths: `A` = the local collection of digitized measurements,
`W` = the workspace holding EXHALE and the reference codes,
`E = W/EXHALE` = this repository.

---

## 0. In-repo digitized observations (already formatted for comparison)

`E/observational_data/` and the planet folders hold hand-digitized
observations already cast in the EXHALE comparison convention: column 2 is
`dF/F` (relative flux change), so transmission `T = 1 + dF/F`. These are
lower-fidelity eye-traces, but they are ready to overplot and cover lines the
raw data files in Sections 1-2 lack (Ly-alpha for HD209458b, H-alpha for
HD189733b). Provenance is in each file header and `E/observational_data/README.md`.

| Planet | Line | File | x / y | Reference |
|---|---|---|---|---|
| **HD209458b** | Ly-alpha | `E/observational_data/vidalmadjar2003_HD209458b_Lya.txt` | velocity [km/s] / dF/F; core +/-40 km/s omitted (geocorona/ISM), 22 rows | Vidal-Madjar et al. 2003 (Nature 422,143) / Ehrenreich et al. 2008 |
| **HD209458b** | H-alpha | `E/HD209458b/jensen2012_HD209458b_Ha.txt` (+ `_coarse`) | wavelength [A] / dF/F, 96 rows; noisy trace of Christie+2013 Fig.14 | Jensen et al. 2012 (ApJ 751,86) |
| **HD189733b** | H-alpha | `E/observational_data/cauley2015_HD189733b_Ha.txt` | wavelength [A] / dF/F, 23 rows | Cauley et al. 2015 (ApJ 810,13) |
| **HD189733b** | H-alpha | `E/HD189733b/jensen2012_HD189733b_Ha.txt` | wavelength [A] / dF/F, 21 rows; clean -2.2% dip, trace of Christie+2013 Fig.6 | Jensen et al. 2012 (ApJ 751,86) |
| **HD189733b** | He I 10830 | `E/observational_data/salz2018_HD189733b_He.txt` | wavelength [A] / dF/F, 18 rows | Salz et al. 2018 (A&A 620,A97) |
| WASP-69b | He I 10830 | `E/observational_data/nortmann2018_WASP69b_He.txt` | wavelength [A] / dF/F, 17 rows | Nortmann et al. 2018 (Science 362,1388) |

`E/observational_data/create_observational_data.py` regenerates these; the
`Halpha_compare_*.ipynb` notebooks in the planet folders already overlay the
Jensen+2012 traces on the EXHALE model. **Not observations:** the
`paper_tpm_*.txt` files in `E/benchmarks/*/` and the planet folders are saved
EXHALE transit model curves (columns `T_theo / T_instr / T_rot+instr`), used as
regression/reference snapshots — not measured data.

---

## 1. Ready-to-use observations, by model planet

| Planet | Line | File | y-axis convention | Coverage | Instrument / reference |
|---|---|---|---|---|---|
| **WASP-52b** | H-alpha | `A/WASP-52b/WASP52b_Halpha_transpec_unbinned.txt` | `TS` = F_in/F_out - 1 (fraction) + `e_TS`; `TS_model` col is all zeros | 6559.76-6565.81 A, ~0.01 A, 606 rows | VLT/ESPRESSO, Chen et al. 2020 (A&A 635, A171) |
| **WASP-52b** | H-alpha (binned) | `A/WASP-52b/WASP52b_Halpha_transpec_0p10_AA_binned.txt` | same convention, 0.10 A bins | 6559.85-6565.75 A, 60 rows | as above |
| **WASP-52b** | He I 10830 | `A/WASP-52b/wasp52_He10830_trans_spec.dat` | excess absorption in **percent** (+/- ~0.5%) + 1-sigma | 10827.94-10837.76 A, planet rest frame, 95 rows | Keck/NIRSPEC, Kirk et al. 2022 (AJ 164, 24) |
| **WASP-121b** | He I 10833 | `A/WASP-121b/WASP-121b_3D-Models_RT_code/Transmission_spec_observation.dat` | normalized flux (in/out ratio ~1.0; absorption = 1 - flux) + err + weight | 10830.08-10837.95 A, 74 rows | VLT/CRIRES+, Czesla et al. 2024 (A&A 692, A230) |
| **HD189733b** | He I 10830 | `A/SPIRou/sp/HD189733_b.dat` | normalized flux + 1-sigma (**not** a ratio yet) | order #71, 10703.9-10953.4 A, 3057 rows | CFHT/SPIRou, Masson et al. 2024 (A&A 688, A179); confirmed He detection |
| **HD209458b** | He I 10830 | `A/SPIRou/sp/HD209458_b.dat` | normalized flux + 1-sigma (**not** a ratio yet) | order #71, 10704.8-10963.1 A, 3184 rows | as above; tentative He detection |

Notes:
- The WASP-52b files carry the full provenance (Yan et al. 2022, ApJ 936, 177,
  shared by Dongdong Yan; original data from Chen+2020 and Kirk+2022). The
  `WASP-52b_YanDongdong/` folder is byte-identical to `WASP-52b/` plus the
  provenance email — treat as one dataset.
- WASP-121b: the same folder is mostly the Czesla et al. 2024 **3D model**
  (Athena++ + Monte Carlo RT) with 84+84 synthetic He 10833 spectra in
  `Spectra_1.0/` and `Spectra_0.5/`; only `Transmission_spec_observation.dat`
  is the measurement. Use the observation for validation, the model spectra as a
  model-to-model cross-check.
- SPIRou spectra are **reduced normalized-flux** spectra over the He order, not
  pre-formed excess absorption — build the in/out-of-transit ratio (or read the
  absorption straight from Masson+2024) before comparing.

## 2. Secondary / validation observations (He I 10830)

| Planet | File | y-axis | Reference |
|---|---|---|---|
| GJ1214b | `A/CARMENES/J1214b_TS.dat` | excess absorption, percent | CARMENES GTO (ref to resolve) |
| GJ1214b | `A/SPIRou/sp/Gj1214_b.dat` | normalized flux | Masson+2024 (nondetection) |
| HAT-P-11b | `A/DACE/HAT-P-11_b_He.dat` | excess absorption, percent | p-winds tutorial dataset |
| HAT-P-11b | `A/SPIRou/sp/HAT-P-11_b.dat` | normalized flux | Masson+2024 (confirmed) |
| WASP-107b | `A/DACE/WASP-107_b_He.dat` | excess absorption, percent | p-winds tutorial dataset |
| WASP-69b | `A/p-winds_data/WASP-69_b.dat` | normalized flux + err | off-list bonus |
| GJ436 | `A/CARMENES/GJ436.zip` (raw), `A/SPIRou/sp/Gj436_b.dat` | raw FITS / normalized flux | needs reduction |

Also `A/p-winds_data/jisan/wasp-52b*.txt` hold a separate (earlier/rougher)
WASP-52b He 10830 reduction (29 and 209 rows) — cross-check only; prefer the
provenance-clean `A/WASP-52b/` files above.

## 3. Stellar SED / EUV inputs (model drivers, not comparison data)

| Star / planet | File | Coverage |
|---|---|---|
| HD209458 | `W/p-winds/data/HD209458b_spectrum_lambda.dat` (= `p-winds_org/`) | 1-1200 A XUV/FUV; host-star irradiation at the planet from the X-exoplanets database (INTA-CSIC), photons->erg and scaled to the HD209458b orbit (p-winds `quickstart.ipynb`). Stellar SED driver, not a transmission observation. |
| GJ1214 | `A/CARMENES/J1214_EUV.dat`; `A/MUSCLES/gj1214/` | 1-1200 A; panchromatic |
| GJ436 | `A/MUSCLES/gj436/` | panchromatic MUSCLES |
| HAT-P-11 | `A/DACE/HAT-P-11_spec.dat` | 4.8-2600 A at planet |
| WASP-69 / HD97658 | `A/p-winds_data/real_data_wasp-69.txt`, `hd97658_...txt` | X-ray-FUV |

`A/MUSCLES/` (STScI HLSP MUSCLES/Mega-MUSCLES) has **none of the four primary
hosts** — only low-mass K/M dwarfs plus a few hot-Jupiter hosts (HAT-P-12/26,
WASP-17/43/77a/127). HD189733, WASP-52, and WASP-121 SEDs are not on disk here;
those runs use the literature-pinned inputs already in their planet folders.

## 4. Rival models for code-to-code cross-checks (not observations)

- `A/WASP-121b/.../Spectra_{0.5,1.0}/` — Czesla et al. 2024 3D He 10833 spectra (WASP-121b).
- `A/Nail2025/WASP-52b/` — Athena++ + Cloudy RT reproduction package (WASP-52b);
  `A/Nail2025/HAT-P-67b/` — 3D He 10830 model spectra (HAT-P-67b, off-list).
- `A/Linssen2025_sunset/` — sunbather Parker-wind transit-spectrum grid (~200
  planets incl. GJ1214b, GJ436b, HAT-P-11b; no primary planet).
- `A/Linssen2024_sunbath/`, `A/MacLeod2022/`, `A/MacLeod2025/` — escape-code
  reproduction packages / 3D wind simulations.
- `W/p-winds_org/HeI_test/` — LaRT Monte Carlo **synthetic** He 10830 image
  cubes for HD209458b (filenames say `_obs` but are LaRT observer-frame output,
  not a measurement). Small `HeI_limb*.fits.gz` hold the ready 1-D emergent
  spectrum (`wavelength`, `Jout`, `Jin`); absorption = 1 - Jout/Jin.

## 5. Coverage summary

| Planet | He 10830 obs | H-alpha obs | Ly-alpha obs | Stellar SED on disk |
|---|---|---|---|---|
| HD209458b | SPIRou (tentative, raw) | Jensen+2012 (digitized) | Vidal-Madjar+2003 (digitized) | yes (p-winds) |
| HD189733b | Salz+2018 (digitized) + SPIRou (confirmed, raw) | Cauley+2015 & Jensen+2012 (digitized) | - | no |
| WASP-52b | Kirk+2022 (raw, ready) | Chen+2020 (raw, ready) | - | no |
| WASP-121b | Czesla+2024 (raw, ready) | - | - | no |

Two data tiers: **raw reductions** (Sections 1-2, higher fidelity, y-convention
varies per file) and **in-repo digitized traces** (Section 0, `dF/F`, ready but
eye-traced). No H-beta, Mg II, Ca II, or Na I D transmission spectra are on disk
for any planet (only a Na/K mention in the WASP-52b ESPRESSO paper, Chen+2020,
not stored as data). Richest tests: **WASP-52b** (raw H-alpha + He 10830),
**HD189733b** (digitized H-alpha + He, plus raw SPIRou He), and **HD209458b**
(digitized H-alpha + Ly-alpha, plus raw SPIRou He). **WASP-121b** has one ready
raw He 10830 spectrum.
