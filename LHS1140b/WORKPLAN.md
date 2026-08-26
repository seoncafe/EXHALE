# LHS 1140b — execution board

Plan of record: `../docs/lhs1140b_lower_atmosphere_plan_new.md` (Phases A-F).
This file tracks execution; every product lives under this directory.
Target paper: Cherubim et al. (2026, Science), PDFs in `../../references/`.

## Layout

```
LHS1140b/
  WORKPLAN.md          this board
  LHS1140b_analysis.ipynb  results notebook, Phases A-C incl. the He:H
                       profile comparison (executed; rerun with
                       jupyter nbconvert --execute from this folder)
  system_parameters.md Phase A1: verified system + retrieval parameters
  sed/                 Phase A2: stellar-spectrum construction
  pwinds_oracle/       Phase A3: frozen p-winds reference case, plus the
                       H:He scan of the He 10830 line (scan_hhe.py)
  (wind run dirs)      Phases C/F: EXHALE runs, one directory per case

Line metrics: `../he_line_metrics.py` (repo root, shared with
EXHALE_transit.py since Phase B).
```

## Status

| Step | What | State |
|---|---|---|
| A1 | Scaffold + parameter record with provenance (`system_parameters.md`; Cadieux 2024 PDF added to `references/`) | **done 2026-08-22** |
| A2 | LHS 1140 spectrum from Mega-MUSCLES proxies at the b orbit (`sed/`). Rebuilt on **v23** (the paper-era release, recovered from `~/Exoplanetary_Atmosphere/MUSCLES/` after MAST withdrew it); C confirmed to be the catalog X-ray flux ratio; F_XUV = 39.0 vs the paper's 33; EUV band differs 4.2x between v23/v25 (measured) | **done 2026-08-22** |
| A3 | p-winds oracle on the v23 SEDs (`pwinds_oracle/`): published point (T = 5160 K) overshoots red 1.8x — remaining gap is the paper's unpublished p-winds configuration; `matched_gj1132` (T = 6100 K, inside their range) reproduces the 2024 line (red 0.1 sigma); levers measured (T strong, Mdot weak/saturated, FUV null — A31-dominated sink) | **done 2026-08-22** |
| A4 | Metric extractor per the paper's three-Gaussian definition, selftest PASS (`../he_line_metrics.py`) | **done 2026-08-22** |
| B | `EXHALE_resolved.out` written by the wind (post-base.inp values), consumed by `EXHALE_transit.py` (fallback: input.inp); metric fit + `tpm_He10830_metrics.txt` in the transit tool; override + fallback verified on a scratch run; `make check` **5/5 byte-identical PASS** (`docs/Update_EXHALE.md` sec 64) | **done 2026-08-22** |
| C | He-rich audit, **rescoped 2026-08-24** onto the converged LHS 1140 b solutions themselves (`exhale/`, He/H = 0.55 ... 10000 on both SEDs) instead of the tutorial-planet series in `heh_series/`, whose He/H = 100 marching never reached a usable state in any of three variants. `heh_series/audit_checks.py` now takes case directories as arguments; record: `exhale/audit_summary.md`. All ten cases finite, charge neutrality to 4.1e-16, heating/cooling channel closure to 3.4e-15/4.2e-16; a JFNK finish of He/H = 1000 under `-fcheck=bounds,do,mem` (`exhale/heh1000_fcheck/`) exits info=0 with no runtime check firing; the H/He path carries no density floors, only division guards sitting >=1e29x below the solution (n_HI reaches an exact zero where H is fully ionized, and nothing divides by it). Not covered: resolution sensitivity, and the triplet/line-cooling photon-escape assumptions of the plan's fourth bullet. The one thing the trapping run did expose -- `0/0` in the primitive conversion of not-yet-filled ghost cells (`Apply_BC`, JFNK line search), harmless to results -- is fixed (`docs/Update_EXHALE.md` sec. 66; `make check` 5/5 byte-identical) and the JFNK path now runs clean under `-ffpe-trap` at both ends of the composition axis | **done 2026-08-24 (rescoped)** |
| D | Conservative binary H/He diffusion (`../docs/binary_diffusion_design.md`, through decision D3). On LHS 1140 b the operator with the default `He_Kzz = 0` removes the helium the observation is about, so the eddy coefficient had to be argued: literature survey (`kzz_literature.md`) + seven-point scan (`kzz_decision.md` sec. 3, cases `exhale/heh0p55_diff_kzz*`). **`He_Kzz = 1.0e9` adopted 2026-08-25 (user decision, `kzz_decision.md` sec. 0)** | **done 2026-08-25** |
| D2 | Composition rescanned under diffusion at the adopted `He_Kzz = 1e9` (`exhale/heh{1,2,2p13,4}_diff_kzz1e9`, all `info = 0`): the He 10830 equivalent width crosses the measured 1.108 +/- 0.030 %A at **He/H = 2.09** (bracketed by runs at 2.00 and 2.13), against 0.55 in the well-mixed limit. Record: `kzz_decision.md` sec. 6. `../docs/lhs1140b_exhale_vs_pwinds.tex` still states 0.55 and takes this at its next revision | **done 2026-08-25** |
| E | Photochem adapter + flux-continuity coupling | pending |
| F | Science runs (prescribed H:He first, then coupled) + collisional validity | pending |

## Decisions taken here (report-level, not user decisions)

- The fiducial configuration follows the paper's GJ 1132 SED (their fiducial
  F_XUV = 0.033 W/m^2 is quoted for that SED); the GJ 699 SED is the
  sensitivity alternate, as in the paper.
- Directory convention copies the other planet folders: self-contained,
  notebooks/scripts read `./` relative paths.
