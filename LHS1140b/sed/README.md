# LHS 1140 stellar spectrum at the b orbit: build record (Phase A2)

Built 2026-08-22 by `make_sed.py`; **rebuilt the same day on Mega-MUSCLES
v23** after the paper-era files were found in the local holdings
`~/Exoplanetary_Atmosphere/MUSCLES/` (MAST now serves only v25).  Recipe:
Cherubim et al. (2026) Science Supplement: proxy SEDs scaled by the X-ray
normalization C, then moved to the b orbit by (d/a)^2 with d = 14.96 pc,
a = 0.0946 AU.  Numbers behind the recipe: `../system_parameters.md`.

## Products (v23, the release the paper used)

The numbers are MEASURED 2026-09-21 on the files as they stand, by
trapezoidal integration of each file's own two columns.

| File | A | F_X(0.25-2 keV) at b | F_XUV(10-1300 A) | F_EUV(100-911 A) |
|---|---|---|---|---|
| `lhs1140_sed_gj1132_at_b.txt` (PRIMARY) | 0.5900 | 4.407 | 33.449 | 10.051 |
| `lhs1140_sed_gj699_at_b.txt` (alternate) | 0.0560 | 3.492 | 169.931 | 70.181 |

(erg s^-1 cm^-2.  Paper comparison: F_X 2.94 from their APEC fit; their
GJ 1132 fiducial F_XUV = 0.033 W/m^2 = 33, and the shipped product's 33.449
is that fiducial.  Both products are the `const-res` v23 release times the
catalog X-ray flux ratio, which is the combination the 2026-08-24 sections
below settle; the older table here quoted the `adapt-var-res` generation and
its 39.0 and has been re-measured.)

Format: two columns, bin-center wavelength [A] increasing, flux at the
planet [erg cm^-2 s^-1 A^-1], read directly by EXHALE `sed_read` and
p-winds `make_spectrum_from_file`.  Coverage ~8 A to 3 um: the He 2^3S
threshold band (4.8-13.6 eV) and the H I / He I / He II edges are inside
the tabulated range.

## What the normalization constant C actually is

The Supplement says "we normalize the panchromatic SEDs ... so their
integrated X-ray flux equals that measured for LHS 1140", and quotes
C_699 = 0.056, C_1132 = 0.59.  Tested here on both v23 and v25: the SEDs'
own X-ray integrals do **not** equal the proxies' cataloged fluxes
(GJ 699 1.15x, GJ 1132 1.78x, identical in both versions, so the X-ray
segments did not change between releases), and only
**C = F_X(LHS 1140)/F_X(proxy) from the cataloged measurements** reproduces
the quoted constants (0.5926 = 3.2/5.4, 0.0559 = 2.7/48.3).  `make_sed.py`
uses that ratio AS THE PAPER ROUNDS IT, `A = 0.59` and `A = 0.056`, which is
what the authors' released script states; the file headers carry the value
used.

## The version difference that matters (measured)

v25 replaced the GJ 1132 EUV estimate (euv-scaling in v23) with a DEM
model: F_EUV(100-911 A) at b is **10.1 (v23) vs 42.3 (v25)**, 4.2x in the
helium-ionizing band, while the X-ray is unchanged.  Built on v25, the
forward-modeled He 10830 line at the paper's retrieval point is
correspondingly deeper and the paper's fiducial F_XUV cannot be recovered
(79.6 vs 33).  This was originally documented here as an irreducible
provenance gap with an `fxuv33` rescaled product as a control; **both are
obsolete** now that the v23 inputs exist, and the `fxuv33` product was
removed.

## Inputs kept

`muscles/hlsp_muscles_multi_multi_{gj1132,gj699}_broadband_v23_{const-res,adapt-var-res}-sed.fits`,
the paper-era release, copied 2026-08-22 from
`~/Exoplanetary_Atmosphere/MUSCLES/{gj1132,gj699}/` (withdrawn from MAST).
`make_sed.py` reads the `const-res` pair; the `adapt-var-res` pair is what the
superseded products were built from.
`muscles/...v25...` is the current MAST release (downloaded 2026-08-22),
kept for the version-difference record above.  "adapt" = negative-flux bins
adaptively rebinned; residual negatives (36 / 5 rows in v23) are clipped at
build time.

## 2026-08-24: the paper's exact input identified

The paper's SED is pinned by requiring two of its published numbers at once:

- Supplement: normalization constants A_1132 = 0.59, A_699 = 0.056 (the
  catalog X-ray flux ratios; the prose describes integral-matching, but no
  v23 product integrates to its catalog flux, so the two definitions differ
  by 1.78 for GJ 1132);
- main text: fiducial XUV flux at the orbit 0.033 W/m2 (10-1300 A).

Sweeping every v23 GJ 1132 product x normalization: only the plain
`const-res` product x A = 0.59 lands on the fiducial (33.6 measured).
The `adapt-*` products x A give 36-39; integral-matched gives 20-22.
`make_sed.py` now builds from `const-res` x catalog ratio; the superseded
products are kept as `*_catalogC.txt` (adapt-var-res x A) and
`*_integralC.txt` (adapt-var-res, integral-matched).

Validation of the result (paper-matched SED, at the orbit): F_XUV = 33.60
(paper 33), F_X(0.25-2 keV) = 4.43 (paper's own APEC measurement is 2.94 --
the proxy SED and the APEC fit disagree in the X-ray band, and any model
built on the proxy SED inherits that), He-2^3S-ionizing band 4.8-13.6 eV = 26.5
(halved from the adapt-var-res product's 58.8, the largest physical change).

## 2026-08-24 (later): which normalization the retrieval actually consumed

The two conventions cannot both describe one input, and the paper carries
both: the quoted constants (A = 0.59/0.056) and the fiducial XUV flux
(0.033 W/m2) belong to the catalog-ratio convention, while the Supplement's
stated procedure is integral matching. The equivalent-width test against
the paper's own plotted best-fit curve (Fig. 4 purple, digitized: depth
1.323 %, FWHM 0.708 A, EW 0.912 %.A) decides it:

| SED fed to p-winds at the Fig. S5 medians | EW [%.A] | vs 0.912 |
|---|---|---|
| const-res x catalog A = 0.59 | 1.128 | +24 % |
| const-res x integral-matched A_eff = 0.35 | 0.906 | +0.7 % |

That test read the retrieval as consuming the integral-matched SED, and the
shipped products were switched to it for a time.

**THE SWITCH WAS REVERTED, and the shipped products carry the catalog ratio
A = 0.59 / 0.056.**  VERIFIED 2026-09-21 by reading the files: the header of
`lhs1140_sed_gj1132_at_b.txt` states `A = 0.5900 (catalog X-ray flux ratio;
authors' released script)` and `make_sed.py` sets `C_eff = c['A']` with
`A = 0.59`.  The reason is in `make_sed.py`'s own docstring: the
equivalent-width test above was confounded by two broadening switches missing
from our oracle (turbulence and phase averaging), and the authors' released
script (`LHS1140b_zenodo/lhs1140b_2024B_pwinds.py`) normalizes the `const-res`
v23 product by the catalog X-ray flux RATIO, not by an integral match.  The
paragraph that said the shipped products are integral-matched described a
state the directory left on 2026-08-24 and has been corrected here.

What the other files are: `*_fiducialA.txt` is the same `const-res` product at
the unrounded ratio A = 0.5926 (MEASURED: it differs from the shipped file by
0.44 per cent, which is 0.5926/0.59); `*_catalogC.txt` and `*_integralC.txt`
are the superseded `adapt-var-res` generations.

A residual ~13 km/s FWHM broadening kernel between our forward model and
the plotted curve remains unexplained (it reshapes at fixed EW; see
`../pwinds_oracle/README.md`).

## The XUV-scaled grid, and the script that writes it

`lhs1140_sed_gj1132_at_b_xuv<factor>.txt` is the whole flux column of
`lhs1140_sed_gj1132_at_b.txt` times the factor, on the same wavelength grid:
0.01, 0.015, 0.02, 0.03, 0.05, 0.06, 0.065, 0.07, 0.10, 0.15, 0.16, 0.18,
0.20, 0.25, 0.26, 0.28, 0.30, 0.33.  They are the XUV axis of the model
catalog and of the continuations that reach its low end.

The scaling is ACHROMATIC and redefines no band, so every band integral is the
reference integral times the factor and the series is a one-parameter family
in the irradiation strength alone.  The flux a run sees is the file's: with
`Spectrum type: Load from file..` the code overwrites `LX` and `LEUV` with the
file's own band integrals (`sed_read.f90`, "Calculate LEUV and LX
luminosities"), so the luminosity lines of an `input.inp` do not enter.

`scale_sed.py <factor> <tag> "<why>"` writes one, and refuses to overwrite an
existing spectrum.  `scale_sed.py --check` MEASURED 2026-09-21 that the rule
reproduces every scaled file on the tree byte for byte, all 30003 lines each.
