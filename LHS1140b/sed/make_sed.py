#!/usr/bin/env python3
"""Build the LHS 1140 stellar spectrum at the LHS 1140b orbit (Phase A2).

Recipe (Cherubim et al. 2026, Science Supplement, "X-ray observations and
analysis"): take the Mega-MUSCLES panchromatic SEDs of the proxies GJ 699
and GJ 1132 (flux at Earth), multiply by the X-ray normalization C, and
scale from Earth distance to the planet orbit by (d/a)^2.

C is the ratio of the MEASURED X-ray fluxes,
    C = F_X(LHS 1140) / F_X(proxy)   in the proxy's measured band:
    C_699  = 2.7e-15 / 4.83e-14 = 0.0559   (0.3-10 keV)
    C_1132 = 3.2e-15 / 5.4e-15  = 0.5926   (0.2-2.4 keV)
(the paper quotes the rounded 0.056 / 0.59).  Verified here: the SEDs' own
X-ray integrals do NOT equal the catalog fluxes (1.15x / 1.78x, identical
in v23 and v25), so the paper's "normalize the integrated X-ray flux"
sentence operationally means the catalog-ratio C above — only that choice
reproduces their quoted constants.

SED VERSION: the paper used the Mega-MUSCLES release of its time, v23 —
withdrawn from MAST (v25 only there), but held locally at
~/Exoplanetary_Atmosphere/MUSCLES/ and copied into muscles/.  v23 is the
PRIMARY input; with it and C above, F_XUV(10-1300 A) at the b orbit comes
out 38.8 erg/s/cm2, consistent with the paper's fiducial "0.033 W/m^2"
(= 33) to ~18% (their rounding/procedure details unrecoverable).  The
decisive version difference is the EUV reconstruction:
F_EUV(100-911 A) at b is 10.0 (v23 euv-scaling) vs 42.3 (v25 DEM) for
GJ 1132 — a 4.2x change in the helium-ionizing band.

Outputs (two columns: bin-center wavelength [A], flux at the b orbit
[erg cm^-2 s^-1 A^-1], increasing wavelength — read by both EXHALE
`sed_read` and p-winds `make_spectrum_from_file`):

  lhs1140_sed_gj1132_at_b.txt   PRIMARY (paper fiducial proxy, v23)
  lhs1140_sed_gj699_at_b.txt    sensitivity alternate (v23)
"""
import numpy as np
from astropy.io import fits

PC_AU = 206264.806
D_PC  = 14.96               # distance [pc]      (system_parameters.md)
A_AU  = 0.0946              # b semi-major axis [AU]
SCALE_ORBIT = (D_PC*PC_AU/A_AU)**2

HC_KEV_A = 12.398425

# proxy -> (v23 SED, C = catalog X-ray flux ratio)
CASES = {
    # The combination the paper's own numbers pin down: the plain const-res
    # v23 product scaled by the catalog X-ray flux ratio A (Supplement:
    # A_1132 = 0.59, A_699 = 0.056). Cross-check: this and only this
    # combination reproduces the paper's fiducial XUV flux at the orbit,
    # 0.033 W/m2 = 33 erg/s/cm2 (10-1300 A, GJ 1132 proxy) -- measured 33.6
    # here, against 39.0 for adapt-var-res x A and 22.0 for the
    # integral-matched normalization the Supplement prose describes.
    # Superseded products are kept as *_catalogC.txt (adapt-var-res x A) and
    # *_integralC.txt (adapt-var-res, integral-matched).
    # The authors' released script (LHS1140b_zenodo/lhs1140b_2024B_pwinds.py)
    # normalizes const-res v23 by the catalog X-ray flux RATIO, not by an
    # integral match: Lx_ratio_1132 = 0.59, Lx_ratio_699 = 0.056. This repo
    # briefly switched to integral matching (A_eff = 0.35) after an EW test,
    # but that test was confounded by two broadening switches missing from
    # our oracle (turbulence + phase averaging); with the authors' code in
    # hand, 0.59 is confirmed.
    'gj1132': dict(fits='muscles/hlsp_muscles_multi_multi_gj1132_broadband_v23_const-res-sed.fits',
                   A=0.59),
    'gj699':  dict(fits='muscles/hlsp_muscles_multi_multi_gj699_broadband_v23_const-res-sed.fits',
                   A=0.056),
}

W_MAX = 30000.0             # keep native start (~8 A) .. 3 um

def band_flux(w, f, wlo, whi):
    m = (w >= wlo) & (w <= whi)
    return np.trapz(f[m], w[m])

for tag, c in CASES.items():
    with fits.open(c['fits']) as h:
        d = h[1].data
        w = np.asarray(d['WAVELENGTH'], dtype=float)
        f = np.asarray(d['FLUX'], dtype=float)
    keep = (w <= W_MAX)
    w, f = w[keep], f[keep]
    # The adapt-var-res release stitches instrument segments and repeats a
    # wavelength where two of them meet. A repeated abscissa leaves the
    # energy-bin width undefined, and EXHALE's SED reader refuses the file,
    # so collapse each repeat to one bin carrying the mean flux.
    wu, inv = np.unique(w, return_inverse=True)
    ndup = len(w) - len(wu)
    if ndup:
        f = np.bincount(inv, weights=f)/np.bincount(inv)
        w = wu
    neg = int((f < 0).sum())
    f = np.clip(f, 0.0, None)

    # Normalization: const-res v23 x catalog X-ray flux ratio A, exactly as
    # the authors' released script does (see README). An earlier integral-
    # matched normalization was reverted once the authors' code confirmed A.
    C_eff = c['A']              # catalog X-ray flux ratio (authors' script)
    fb = f*C_eff*SCALE_ORBIT

    out = f'lhs1140_sed_{tag}_at_b.txt'
    with open(out, 'w') as fh:
        fh.write(f"# LHS 1140 SED at the LHS 1140b orbit, proxy {tag.upper()}\n"
                 f"# Built by make_sed.py from Mega-MUSCLES v23 const-res-sed\n"
                 f"# (the release the paper used; local copy, withdrawn from MAST)\n"
                 f"# A = {C_eff:.4f} (catalog X-ray flux ratio; authors' released script),\n"
                 f"# d = {D_PC} pc, a = {A_AU} AU, (d/a)^2 = {SCALE_ORBIT:.6e},\n"
                 f"# negative input fluxes clipped: {neg}\n"
                 f"# repeated wavelengths merged: {ndup}\n"
                 f"# col 1: bin center wavelength [A] (increasing)\n"
                 f"# col 2: flux at the planet [erg cm^-2 s^-1 A^-1]\n")
        np.savetxt(fh, np.column_stack([w, fb]), fmt='%.6e')

    print(f"{tag} (v23, C_eff = {C_eff:.4f}): rows={len(w)}, clipped_neg={neg}")
    print(f"  F_X(0.25-2 keV)  at b = {band_flux(w,fb,HC_KEV_A/2.0,HC_KEV_A/0.25):8.3f}"
          f" erg/s/cm2  (paper APEC: 2.94)")
    print(f"  F_XUV(10-1300 A) at b = {band_flux(w,fb,10.0,1300.0):8.3f}"
          f" erg/s/cm2  (paper fiducial, GJ1132: 33)")
    print(f"  F_EUV(100-911 A) at b = {band_flux(w,fb,100.0,911.0):8.3f} erg/s/cm2")
    print(f"  F(4.8-13.6 eV)   at b = {band_flux(w,fb,911.0,12398.425/4.8):8.3f}"
          f" erg/s/cm2  (He 2^3S ionizing band)")
