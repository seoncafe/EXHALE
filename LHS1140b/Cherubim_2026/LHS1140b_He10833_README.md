# LHS 1140b He I 10833 A spectrum

Target: Cherubim et al. (2026), Figure 3B, 2024 in-transit spectrum.

The attached paper identifies Figure 3B black points as the data after Gaussian-process correlated-noise correction. Instead of estimating marker centers from the raster, the CSV uses the numerical spectrum released by the authors in the Zenodo record cited as reference 49 of the paper (DOI: 10.5281/zenodo.20723095). The reconstruction follows the authors' Figure 3 notebook: `fig3b_excess_absorption_percent = raw_excess_absorption_percent + gp_correction_applied_percent`.

## CSV columns

- `wavelength_vacuum_planet_rest_A`: Figure 3B x coordinate.
- `wavelength_air_planet_rest_A`: Figure 4 x-coordinate convention.
- `raw_excess_absorption_percent`: pre-GP-correction spectrum (Figure 2 quantity).
- `gp_correction_applied_percent`: signed GP correction from the released file.
- `fig3b_excess_absorption_percent`: Figure 3B black-point ordinate. Absorption is negative.
- `absorption_depth_positive_percent`: the same signal with positive absorption convention.
- `uncertainty_1sigma_percent`: 1-sigma uncertainty plotted in Figure 3B.
- `fig4_normalized_flux`: `1 + fig3b_excess_absorption_percent / 100`.
- `fig4_normalized_flux_uncertainty_1sigma`: normalized-flux 1-sigma uncertainty.

There are 129 samples spanning 10830.0220606-10834.6442331 A (vacuum, planetary rest frame).
The QA overlay uses red rings; their alignment with Figure 3B black markers verifies the reconstruction. The paper reports a fitted line depth of 1.24 (+0.22/-0.23) percent and a Doppler shift of 2.0 (+2.0/-2.2) km/s.

---

## Checks run on this reconstruction (2026-08-23)

**Our own line-metric extractor reproduces the authors' fit.** Running
`../../he_line_metrics.py` -- the three-Gaussian fit with shared width and
shift that we apply to every EXHALE and p-winds spectrum -- on
`absorption_depth_positive_percent`:

| | this extractor | paper |
|---|---|---|
| blended red depth | 1.254 % | 1.24 (+0.22/-0.23) % |
| blue depth | 0.198 % | 0.25 (+0.14/-0.12) % |
| red/blue amplitude | 6.35 | 6.7 (+12.7/-3.1) |
| FWHM | 0.841 A | 0.86 (+0.15/-0.27) A |
| Doppler shift | +0.103 A | +0.072 A (2.0 +2.0/-2.2 km/s) |

Every value falls inside the published uncertainty. The models can therefore
be compared to the measurement through one extractor rather than against
numbers fitted by a different method.

**The earlier figure digitization agrees.** `fig2_digitization/` recovered
Fig. 2 from the PDF raster before this file was available; its deepest point
is -1.488 % at 10833.426 A against -1.497 % at 10833.416 A here. It is kept
as an independent check and is not the reference.

## Which column to compare a model to

`fig3b_excess_absorption_percent` (equivalently
`absorption_depth_positive_percent`) is the spectrum the paper fits and
quotes 1.24 % from. `raw_excess_absorption_percent` is Fig. 2, before the
correlated-noise correction; fitting it gives 1.339 %, so the choice of
column moves the target by about 7 %.
