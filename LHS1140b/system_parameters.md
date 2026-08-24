# LHS 1140 system parameters — verified record (Phase A1)

Every value below was read from the named published PDF (all in
`../../references/`), not from memory or third-party summaries.
C26 = Cherubim et al. (2026, Science, aea9708; + Supplement),
C24 = Cadieux et al. (2024, ApJL 960, L3).

## Star (LHS 1140 = GJ 3053)

| Quantity | Value | Source |
|---|---|---|
| Spectral type | mid-M dwarf (M4.5, C24 abstract) | C24 |
| M_star | 0.1844 +/- 0.0045 M_sun | C24 Table (also quoted in C26 Supplement) |
| R_star | 0.2159 +/- 0.0030 R_sun | C24 |
| L_bol | 0.0038 +/- 0.0003 L_sun | C26 (citing C24) |
| Distance | 14.96 +/- 0.01 pc (14.97 in the Supplement X-ray section) | C26 |
| Rotation period | 131 d; age >~ 3.1 Gyr | C26 |
| L_X (0.25-2.0 keV) | 7.8 +1.1/-1.3 x 10^25 erg/s | C26 Supplement |
| L_X/L_bol | 6.49 +1.34/-1.25 x 10^-6 | C26 Supplement |

## Planet b

| Quantity | Value | Source |
|---|---|---|
| M_p | 5.60 +/- 0.19 M_Earth | C24 (abstract), C26 |
| R_p | 1.730 +/- 0.025 R_Earth | C24 (abstract), C26 |
| P | 24.73723 +/- 0.00002 d | C24 |
| a | 0.0946 +/- 0.0017 AU | C24 |
| Inclination | 89.96 +/- 0.04 deg | C26 Supplement (citing C24) |
| Irradiation | 42% of Earth's | C26 |
| T_eq (A=0) | 226 +/- 4 K | C26 |
| F_X at orbit (0.25-2.0 keV) | 2.94 +0.61/-0.57 erg/s/cm^2 | C26 Supplement |

Consistency checks run for this record (not published values):
`a` from Kepler's third law with M_star and P gives 0.0945 AU (matches);
L_bol/a^2 gives 42.5% of Earth's irradiation (matches);
(1.52 R_p / R_star)^2 = 1.25% reproduces the paper's excess-depth <->
equivalent-opaque-radius conversion, confirming R_star = 0.2159 R_sun is the
radius behind the C26 numbers.

Note the convention that check implies. The paper's 1.52 R_p follows from
delta = (R_eff/R_star)^2. The plotted spectrum, however, has its continuum
at zero, so the white-light transit is already divided out and the measured
delta is the *excess* area: delta = (R_eff^2 - R_p^2)/R_star^2, which for
delta = 1.24% gives R_eff = 1.82 R_p. The two differ by 20% in radius. Model
comparisons here are made on depth, which is unambiguous; where an effective
radius is quoted, the convention is stated with it.

## 2024 helium detection (C26) — the Phase A/F comparison targets

The numerical spectrum behind these is in `Cherubim_2026/`, reconstructed
from the authors' Zenodo deposit (10.5281/zenodo.20723095), not read off a
figure. Our own extractor (`../he_line_metrics.py`) fitted to it returns
red 1.254%, blue 0.198%, ratio 6.35, FWHM 0.841 A -- every value inside the
published uncertainty -- so models and measurement can be compared through
one extractor.

| Quantity | Value |
|---|---|
| Excess depth, blended red pair | 1.24 +0.22/-0.23 % |
| Excess depth, blue line | 0.25 +0.14/-0.12 % |
| Blended-red / blue amplitude ratio | 6.7 +12.7/-3.1 (optically thin expectation 8) |
| FWHM of the blended feature | 0.86 +0.15/-0.27 A = 23.9 +4.2/-7.5 km/s |
| Doppler shift | +0.072 +0.080/-0.073 A = 2.0 +2.2/-2.0 km/s |
| Equivalent opaque radius | 1.52 R_p |
| Triplet rest wavelengths (vacuum) | 10832.057, 10833.217, 10833.306 A |
| 2025 upper limit | ~0.6% (no detection) |

## p-winds retrieval point (C26) — the Phase A3 oracle target

| Quantity | Value |
|---|---|
| Mdot | 2.03 +0.58/-0.67 x 10^8 g/s |
| T_wind | ~5160 K (best fit; grid explored to 6400 K) |
| H:He number ratio | <~ 1e-3 (fit requires low H) |
| Fiducial F_XUV (10-1300 A, at the orbit) | 0.033 W/m^2 = 33 erg/s/cm^2 |
| Instrument LSF | Gaussian FWHM 4.4 km/s (R = 68,000) |
| Energy-limited check | (6.2-29) x 10^7 g/s |
| Survival bound on the mean rate | < 5 x 10^9 g/s |

## Stellar-spectrum recipe (C26 Supplement) — the Phase A2 specification

Proxy stars of matched type and rotation, panchromatic Mega-MUSCLES SEDs
(X-ray to IR, EUV reconstruction included), normalized so the integrated
X-ray flux equals the XMM measurement for LHS 1140:

| Proxy | Type / P_rot | X-ray flux at Earth (band) | LHS 1140 in that band | Scale C |
|---|---|---|---|---|
| GJ 699 | M4 / 145 d | (4.83 +/- 0.52) x 10^-14 (0.3-10 keV) | (2.7 +/- 0.5) x 10^-15 | **0.056** |
| GJ 1132 | M3.5 / 125 d | (5.4 +/- 1.4) x 10^-15 (0.2-2.4 keV) | (3.2 +/- 0.6) x 10^-15 | **0.59** |

Predicted LHS 1140 flux at Earth: C x F_proxy(Earth); at the b orbit,
multiply by (d/a)^2 with d = 14.96 pc, a = 0.0946 AU.  The paper's fiducial
is the **GJ 1132** SED (the 0.033 W/m^2 and the fiducial Mdot are quoted for
it); GJ 699 is the sensitivity alternate.  EXHALE use additionally requires
coverage down to 4.8 eV (He 2^3S threshold in `sed_read`), which the
panchromatic SEDs provide.
