# LHS 1140 b He 10830: is the measured line width what we think it is?

Audit of the *measurement* side of the He I 10830 width discrepancy
(observation 0.841 A against models 0.233-0.271 A). It asks two questions:
does anything on the observing or fitting side inflate the published width,
and are the model and the measurement quantities defined the same way?

Sources read directly: `../../references/Cherubim_2026Science.pdf` and
`../../references/Cherubim_2026Science_Supplement.pdf` (via `pdftotext
-layout`), the authors' released spectrum and MCMC posterior in
`../LHS1140b/Cherubim_2026/`, and the extractor and transit code in this
repository. Every number below was recomputed here; none is quoted from a
previous memo.

Note on quotations: the Science PDF text layer renders math italics as
substitute glyphs (`R` appears as a script character, minus signs as arrows,
sigma as a script `N`). Quoted passages restore those symbols inside square
brackets and are otherwise character-for-character.

---

## Verdict

**The measurement side explains essentially none of the missing width.** Of
the ~22 km/s of extra broadening the models need, the instrument line spread
function, the triplet fine structure, stellar rotation, planetary rotation
and orbital smearing together account for **at most about 0.4 km/s**, and
most of that is already removed when the comparison is done consistently.
The like-for-like requirement, with both sides deconvolved from their own
instrument kernels, is **21.96 km/s FWHM (sigma = 9.33 km/s)**, against the
22.06 km/s obtained by the naive quadrature of the two convolved widths.
The audit moves the requirement by 0.1 km/s, not by a factor.

Two secondary results run the other way, i.e. they make the tension slightly
worse or the comparison slightly cleaner rather than closing the gap:

1. Our EXHALE transit spectra for this target were convolved at **R = 80,000**
   (the CARMENES default in `EXHALE_transit.py:126`) instead of WINERED's
   **R = 68,000**. Corrected -- and it has since been corrected, see item 3
   of section 3 -- the model blend width at `flux_closure/heh10p7/k01` is
   0.2815 A rather than 0.2684 A, so the naive gap shrinks by 0.9 %, and the
   deconvolved requirement is unchanged.
2. The published FWHM uncertainties appear **transposed**. The released
   posterior gives FWHM = 0.865 +0.261/-0.155 A (23.93 +7.24/-4.30 km/s);
   the paper prints 0.86 +0.15/-0.27 A (23.9 +4.2/-7.5 km/s). The 1 sigma
   lower bound on the measured width is therefore 0.709 A, not 0.59 A, so
   narrow models are excluded more firmly than the printed bar suggests.

Accounting, in the units the requirement is stated in:

| Term | Size | Effect on the ~22 km/s requirement |
|---|---|---|
| Instrument LSF, R = 68,000 | 4.41 km/s FWHM (0.1593 A) | removed by deconvolution; 22.4 -> 22.0 km/s |
| Triplet splitting, 0.0897 A | 0.7 % of the observed blend width, 8.2 % of the model's | present on both sides, nearly cancels; < 0.1 km/s |
| Our R = 8e4 instead of 6.8e4 | model 0.2684 -> 0.2815 A | 0.2 km/s on the naive comparison, 0 on the deconvolved one |
| Stellar rotation (P_rot = 131 d) | v_eq = 0.083 km/s, broadening ~0.14 km/s | < 0.05 km/s in quadrature |
| Planet rotation (tidally locked) | 0.032 km/s at R_p, 0.097 km/s at 3 R_p | negligible |
| Orbital RV smearing across transit | 0.93 km/s if uncorrected; corrected by the authors | residual 0.022 km/s (ephemeris) |
| Smearing within a 300 s exposure | 0.037 km/s | negligible |
| **Residual not explained by the measurement** | | **21.9 - 23.2 km/s FWHM** |

The range in the last row is the full spread over the convention choices
examined below (our extractor vs. the paper's posterior, deconvolved or
not), not an error bar.

---

## 1. What the published width actually is

**Finding.** The quoted 0.86 A is `2.3548 x sigma` of a **single Gaussian
component** whose width is **shared by all three triplet lines**, fitted to
the excess-absorption spectrum in the planetary rest frame. It is not a
red-component-specific width, and it is not the width of the blended red
pair. Whether it is deconvolved from the instrument profile could not be
settled from the paper alone; the released posterior weakly indicates it is
not, and the distinction is worth 1.7 % in width either way.

**Evidence.**

The main text states the model and its free parameters:

> "helium is expected to produce a triplet of closely spaced absorption
> lines, which we modeled with three Gaussian profiles at 10,832.057,
> 10,833.217 and 10,833.306 Å [rest wavelengths in vacuum (19)]; the latter
> two lines are blended at the resolution of the WINERED spectra."

> "The MCMC process followed a Bayesian retrieval framework with five free
> parameters: the three peak amplitudes, a shared peak width, and a shared
> Doppler shift (18)."

and reports

> "The full-width at half-maximum (FWHM) of the blended helium absorption
> lines is 0.86+0.15−0.27 Å, corresponding to 23.9+4.2−7.5 km s−1."

The authors' released MCMC chain,
`../LHS1140b/Cherubim_2026/LHS1140b_He10833_reproduction_files/lhs1140b_flat_samples_GP_IT_vac_3lines.pickle`,
is a (46400, 5) array whose columns are the five parameters above. Its
medians and 16/84 percentiles are

```
col0 (blue amplitude)  4.733  +3.844 -3.024
col1 (red amplitude 1) 15.346 +11.463 -10.481
col2 (red amplitude 2) 16.193 +10.886 -10.894
col3 (shared width)    0.36714 +0.11105 -0.06600   [A]
col4 (shared shift)    0.07175 +0.07769 -0.07273   [A]
```

col4 reproduces the published shift `0.072 +0.080/-0.073 A` exactly, which
identifies the chain. Then

* `2.3548 x 0.36714 = 0.8646 A`, and `0.8646/10833 x c = 23.93 km/s`. Both
  published numbers, to the digits printed.
* The width of the blended red pair, measured numerically per sample, is
  `0.879 +0.259/-0.162 A` instead. That is **not** what the paper printed.

So despite the wording "of the blended helium absorption lines", the
published number is the component FWHM. The two conventions differ by 1.7 %
at the observed width, so nothing downstream depends on it, but the
description is loose.

**A second width convention inside the same paper.** The Supplement's
stellar-variability test reports

> "We used the same MCMC process as before to determine the posterior
> probability distributions of the peak amplitudes, Gaussian standard
> deviations (corresponding to a FWHM), and Doppler shifts. The FWHM of the
> stellar helium line (FWHM_stellar) is 0.22176+0.0027−0.00094 Å the FWHM of
> the in-transit feature (FWHM_in−transit) is 0.323+0.012−0.011 Å and the
> FWHM of the pre-ingress feature (FWHM_pre−ingress) is 0.392+0.018−0.017 Å"

Those numbers are **sigma, not FWHM**. Fitting the same three-Gaussian model
with our extractor to the same in-transit spectrum without a GP (the
`raw_excess_absorption_percent` column of the released file, which is the
pre-GP spectrum the test describes) returns `sigma = 0.3241 A`, matching
their 0.323 to 0.3 %; the corresponding FWHM is 0.763 A, matching neither.
Converted, the Supplement's three widths are FWHM 0.522, 0.761 and 0.923 A.
The `~9 sigma` separation the paragraph concludes with is a difference
between two widths and is unaffected by the mislabeling; the absolute values
are.

**Rest frame and bulk velocity.** The line-of-sight motion of the planet is
removed before averaging:

> "The planetary transmission spectrum was determined by shifting the time
> series spectra into the planetary rest frame, assuming a circular orbit,
> using the previously measured orbital inclination of 89.96 ± 0.04 degrees
> and stellar mass of 0.1844 ± 0.0045 Solar masses (L_Sun) (13), and taking
> the mean of all in-transit spectra (fig. 2)."

Any residual bulk atmospheric velocity is absorbed by the shared shift
parameter (0.072 A = 2.0 km/s), not by the width. The width parameter
therefore measures the velocity *dispersion* along the line of sight, which
is the quantity the models are being asked to produce.

---

## 2. Is the observed width just the instrument profile?

**Finding. No.** The line is resolved by a factor 5.4. This candidate is
dead, and it was the one that would have made the rest unnecessary.

The Supplement fixes the kernel:

> "To incorporate instrumental broadening, we assumed a Gaussian line spread
> function with FWHM of 10,833 Å/[R], where [R] = 68,000 is the resolving
> power of the WINERED spectrograph in HIRES-Y mode. We validated this
> kernel by comparing it to nearby telluric features, finding that it
> reproduced the observed widths. The observed helium absorption feature has
> a larger FWHM than the instrument broadening profile, therefore thermal
> broadening dominates over instrumental and natural (Lorentzian)
> broadening."

and the same kernel is used in their p-winds retrieval:

> "the instrumental line spread function which we modeled as a Gaussian with
> FWHM = 4.4 km s−1 (equivalent to [R] = 68,000)."

Numbers:

```
LSF FWHM  = 10833/68000 = 0.15931 A = 4.409 km/s   (sigma = 0.06765 A)
observed  = 0.8646 A = 23.93 km/s                  -> 5.43 x LSF
observed, 16th percentile of the released posterior
          = 0.7092 A = 19.63 km/s                  -> 4.45 x LSF
```

Sampling is 0.03611 A per point (129 points over 4.622 A), i.e. 4.4 samples
per LSF FWHM, consistent with the Supplement's "~5 pixels per resolution
element"; a rectangular pixel of that size adds 0.3 % in quadrature to the
LSF and nothing measurable to a 0.86 A line.

Deconvolving the LSF from our own fit to the released spectrum:

```
sigma_obs (as measured)  = 0.35473 A
sigma_obs (LSF removed)  = sqrt(0.35473^2 - 0.06765^2) = 0.34822 A = 9.637 km/s
```

The instrument is 3.6 % of the observed width in quadrature.

---

## 3. Triplet structure and extractor conventions

**Finding.** Our extractor and the paper measure widths on the same model
and, to within 1.7 %, the same quantity; this is not an apples-to-oranges
comparison. The 0.0897 A splitting of the red pair is present in both the
observed and the model widths and very nearly cancels: it contributes 0.7 %
of the observed blend width but 8.2 % of the model's, and correcting for it
changes the broadening requirement by less than 0.1 km/s.

**Atomic data.** The transit code carries the NIST triplet with oscillator
strengths `f = 2.9958e-1, 1.7974e-1, 5.9902e-2` for 10830.33977,
10830.25010 and 10829.09114 A in air (`exhale_transit_lib.py:56-76`), i.e.
exactly 5 : 3 : 1, the `2J+1` weights of the `1s2p 3P_{2,1,0}` upper levels.
The red pair is the `J = 2` (stronger, redder) and `J = 1` (weaker) lines,
separated by 0.08967 A in air and 0.089 A in vacuum (10833.306 - 10833.217),
and the blue line carries 1/9 of the total. The paper's optically thin
expectation of 8 for the red-to-blue ratio is this 8 : 1.

**What each side measures.** `he_line_metrics.py:69-75` builds the two red
Gaussians on a fine grid, takes the maximum of their sum, and measures the
FWHM of that sum numerically, so our `fwhm_A` is the **blend** width. The
paper prints `2.3548 sigma` of one component (item 1). Both are computed
from the same five-parameter fit; the difference between them is

```
observed  (sigma = 0.3547):  blend 0.8409 A vs 2.3548 sigma 0.8353 A   (+0.7 %)
paper     (sigma = 0.3671):  blend 0.8708 A vs 2.3548 sigma 0.8645 A   (+0.7 %)
model     (sigma = 0.1048):  blend 0.2684 A vs 2.3548 sigma 0.2467 A   (+8.2 %)
```

The splitting matters eight times more for the narrow model than for the
broad observation, because at `sigma = 0.35 A` the pair is a single hump
whose width the 0.09 A offset barely perturbs. That asymmetry works in the
*wrong* direction for closing the gap: our extractor already credits the
model with 0.022 A of width it does not get from gas motion.

**How much of the observed width is instrument plus splitting alone.** Two
delta components in the physical 3 : 5 ratio, convolved with the R = 68,000
kernel and nothing else, give a blend FWHM of 0.1963 A = 5.43 km/s, which is
23 % of the observed 0.841 A linearly and 3 % of it in quadrature. The
remaining 0.818 A = 22.6 km/s has to come from gas.

**The like-for-like requirement.** Removing each side's own instrument
kernel and comparing component widths:

```
sigma_obs, intrinsic  = sqrt(0.35473^2 - 0.06765^2) = 0.34822 A = 9.637 km/s
sigma_mod, intrinsic  = sqrt(0.10478^2 - 0.05749^2) = 0.08760 A = 2.425 km/s
sigma_extra required  = sqrt(0.34822^2 - 0.08760^2) = 0.33702 A = 9.327 km/s
FWHM_extra required   = 0.7936 A = 21.96 km/s
```

against 22.06 km/s from the naive `sqrt(0.841^2 - 0.268^2)`. Using the
paper's own posterior sigma instead of our fit to its released points, the
requirement is 23.23 km/s (if the published sigma is already deconvolved) or
22.81 km/s (if it is not). Every route lands between 22 and 23.2 km/s.

**One caveat on the fitted amplitudes.** In an unresolved pair the split of
the blend between its two components is not constrained: our fit to the
released spectrum returns `a(10833.217) : a(10833.306) = 1.50`, the paper's
posterior returns 0.95, and the atomic ratio is 0.60. The blend peak and the
shared width are the meaningful outputs; the individual red amplitudes are
not, and neither the depth nor the width above depends on which partition is
used (forcing the physical 3 : 5 split changes the observed blend FWHM from
0.8409 to 0.8414 A).

**Model side, checked in the code.** `EXHALE_transit.py:405-406, 737-738,
772` convolve the disk-averaged and rotation-broadened profiles with a
normalized Gaussian of FWHM `1.0830e-6/Instr_res_HeTR`, and
`he_line_metrics.py` is run on the last column of `tpm_He10830.txt`, which is
the rotation-plus-instrument curve. The metrics files under
`../LHS1140b/exhale/` therefore do describe an instrument-convolved model
measured exactly as the observation is. The one defect is the value of
`Instr_res_HeTR`: every LHS 1140 b `transit.log` in the tree reported
`He 8e+04`, the CARMENES default at `EXHALE_transit.py:126`, where WINERED
HIRES-Y is 68,000.

**Corrected after this audit.** The transit post-processing of every LHS
1140 b run this memo quotes was re-run at `EXHALE_TRANSIT_RES_HETR=68000`
(winds untouched), and the setting is now a convention rather than a
per-command flag: `LHS1140b/winered_hires_y.sh` exports it and every LHS
1140 b run script sources it. The built-in default is deliberately left at
8e4, because it is shared with the other planets in this repository. The
measured effect, at `flux_closure/heh10p3/k01`, is a blend width of 0.2805
against 0.2674 A and a red depth of 3.526 against 3.696 per cent; the
red-pair equivalent width changes by less than 1e-4 %A anywhere, so no
crossing moves. The required broadening falls by 0.12 km/s on every arm.
The pre-correction curves are kept beside the new ones as
`tpm_*_R80k.txt`.

---

## 4. Stellar contamination

**Finding.** Stellar rotation cannot broaden anything: `v_eq = 0.083 km/s`
for this 131 d rotator, so the Rossiter-McLaughlin and center-to-limb
residuals that widen lines on fast rotators are three orders of magnitude
too small here. A stellar residual of a different kind, imperfect division
of a variable chromospheric line, is not excluded by the paper on width
grounds and is not excluded here either, but it cannot supply the observed
width on its own: the star's own He line is narrower than the transit
feature.

**Rotation.** With `R_star = 0.2159 R_sun` and `P_rot = 131 d` (both from
`../LHS1140b/system_parameters.md`, sourced to Cadieux 2024 and Cherubim
2026):

```
v_eq = 2 pi R_star / P_rot = 2 pi (1.502e5 km) / (1.132e7 s) = 0.0834 km/s
rotational broadening ~ 1.7 v sin i = 0.14 km/s = 0.0051 A
```

That is 3 % of the instrument FWHM and 0.6 % of the observed width. The
Rossiter-McLaughlin distortion scales with the same `v sin i` and is
correspondingly irrelevant. Neither "Rossiter", "McLaughlin", "limb" nor
"CLV" appears anywhere in the paper or its Supplement (checked by text
search over both PDFs), so the authors do not discuss these; on this system
there is nothing to discuss.

**What the paper does argue.** Their exclusion of stellar origin is about
depth and time behavior, not width:

> "For the TLS effect to produce our observed excess absorption depth of
> 1.24%, the transit chord would need to be approximately 3.3× brighter at
> 10,833 Å relative to the average flux at this wavelength across the
> stellar disk, given the white-light transit depth of 0.55% (36). This
> would require intense, localized helium emission along the transit chord
> or a perfectly helium-absorbing, unocculted area covering 70% of the
> stellar surface."

and, on flares,

> "LHS 1140 flares infrequently at a rate of one major flare every 5,770
> days (16). Hence it is unlikely that LHS 1140 produced a flare during the
> six and a half hours of observations."

**Width scale of a stellar residual.** From item 1, the star's own
metastable He line has `sigma = 0.22176 A`, i.e. FWHM 0.522 A = 14.5 km/s
observed, 0.497 A = 13.8 km/s with the LSF removed. So an imperfectly
divided stellar line injects a residual of ~14 km/s FWHM, wide compared with
our 0.27 A model line but still well short of the 0.86 A observed feature. A
stellar residual therefore cannot be the whole width, and a mixture of a
narrow planetary line with a displaced stellar residual is not something
this audit can rule out or quantify from the released data.

**A related observation, not a stellar-contamination argument.** The
pre-ingress feature, when the planet occults nothing, has depth 1.01 % and
`sigma = 0.392 A` (FWHM 0.923 A = 25.6 km/s), i.e. it is at least as broad
as the in-transit feature. Whatever produces the broad width is not tied to
the transit chord. The paper reads this as a leading tail.

---

## 5. Planetary rotation and orbital motion

**Finding.** All terms are below 0.5 km/s, and the largest is removed by the
authors' rest-frame correction. Nothing here is at the 20 km/s scale.

With `R_p = 1.730 R_earth`, `P = 24.73723 d`, `a = 0.0946 AU`,
`R_star = 0.2159 R_sun`, `b = 0.23`:

```
tidally locked surface speed  v = 2 pi R_p / P = 0.0324 km/s
  at 1.5 R_p  0.049 km/s ; at 3 R_p  0.097 km/s
orbital speed                 v_orb = 2 pi a / P = 41.60 km/s
transit duration              T14 = (P/pi)(R_star/a) sqrt((1+k)^2 - b^2)
                                  = 7570 s = 2.10 h
LOS RV change, half transit   41.60 x 2 pi x 3785/P = 0.463 km/s
  full transit                0.926 km/s  (equivalent Gaussian FWHM 0.63 km/s)
LOS RV change, one 300 s exposure           0.037 km/s
LOS RV change, 178 s ephemeris uncertainty  0.022 km/s
```

Even entirely uncorrected, the full-transit smear would raise 23.93 km/s to
23.93 km/s (the change is 0.008 km/s, because 0.63 adds in quadrature to
23.9). It is in any case corrected: the spectra are shifted into the
planetary rest frame before averaging (quoted in item 1), leaving only the
ephemeris residual of 0.022 km/s. Planetary rotation is already in the model
curve we measure (`EXHALE_transit.py:772`, with `ROTP` defaulting to the
tidally locked `P_orb`), and at 0.03-0.10 km/s it changes nothing.

---

## Residual

```
required extra broadening, like for like : 21.96 km/s FWHM (sigma 9.33 km/s)
explained by measurement-side effects    : <= 0.4 km/s, most of it already
                                            removed by deconvolving the LSF
residual                                 : 21.9 - 23.2 km/s FWHM
```

The width discrepancy is not a measurement artifact and is not a definition
mismatch between our extractor and the paper. It survives every convention
change examined: component width or blend width, deconvolved or not, our fit
to the released spectrum or the authors' own posterior.

---

## Two findings recorded in passing

**The paper's own fiducial model does not reproduce its own quoted width,
and is not optically thin.** Our reproduction of the authors' p-winds
setting at their retrieved medians (`Mdot = 2.03e8 g/s`, `T = 5160 K`,
`H:He = 1.01e-3`), in `../LHS1140b/pwinds_oracle/tspec_authors_noturb.txt`
and `tspec_authors_turb.txt`, both convolved at R = 68,000, measured with
the same extractor:

| | red depth | red/blue | FWHM |
|---|---|---|---|
| p-winds, turbulence off | 1.447 % | 2.36 | 0.458 A |
| p-winds, turbulence on | 1.376 % | 2.57 | 0.575 A |
| EXHALE `flux_closure/heh10p3/k01` | 3.684 % | 6.98 | 0.268 A |
| measurement | 1.254 % | 6.35 | 0.841 A |

Two things follow. First, the retrieval that produced the published
parameters lands at 0.46-0.58 A, so the 0.86 A width is not reproduced by
the paper's model either; the gap is a shared problem, not an EXHALE-specific
one. Second, the p-winds model reaches its width partly through saturation:
its red-to-blue ratio is 2.36 against the optically thin limit of 8, whereas
our model sits at 6.98 and the measurement at 6.35. The premise that the line
is optically thin, so that column density cannot be traded for width, holds
for our solutions and for the data, but not for the comparison model. Purely
thermal broadening at their retrieved 5160 K is 7.71 km/s FWHM = 0.279 A,
which is close to our model width and far from 0.86 A, so their extra width
comes from the velocity field and from saturation, not from temperature.

**Transposed uncertainties on the published FWHM.** From the released chain,
`FWHM = 2.3548 sigma` has percentiles `0.8645 +0.2615/-0.1554 A`, i.e.
`23.93 +7.24/-4.30 km/s`; the paper prints `0.86+0.15−0.27 Å` and
`23.9+4.2−7.5 km s−1`. The magnitudes match the posterior after swapping the
two, to 2-4 %. This is stated as an apparent transposition rather than a
certainty, since we cannot rule out that the printed bars came from a
slightly different chain than the deposited one. It matters for our purposes
only in one direction: the lower 1 sigma bound on the measured width is
0.709 A, not 0.59 A.

---

## What was not determined

* **Whether the published sigma is deconvolved from the instrument profile.**
  The Supplement says the model incorporates the LSF, which implies the
  fitted sigma is intrinsic, but two checks point weakly the other way: our
  own LSF-free fit to the same points without a GP returns `sigma = 0.3241 A`
  against their 0.323 (a deconvolving fit should have returned 0.3170), and
  refitting the released posterior curve to the released data with a free
  amplitude scale prefers no additional LSF convolution by `delta chi^2 =
  6.6` over 129 points. Neither test is clean, and the answer is worth 1.7 %
  in width, so it was not pursued further.
* **Whether a displaced stellar residual blended with a narrow planetary
  line could mimic the observed width.** The width scale of a stellar
  residual is pinned above (14 km/s FWHM), but the released data do not
  include the stellar template, so the mixture was not tested.
* **The unit of the amplitude parameters in the released chain.** A free
  scale fit gives 1 chain unit = 0.0400 % excess absorption; the origin of
  that factor was not traced. It does not enter any width.
