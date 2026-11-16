#!/usr/bin/env python3
"""He 10830 line metrics, defined as in Cherubim et al. (2026) (Phase A4).

The paper fits the excess-absorption spectrum with THREE Gaussians at the
fixed vacuum rest wavelengths of the He triplet, with five free parameters:
the three peak amplitudes, one shared width, and one shared Doppler shift
(Science main text, "Interpretation" section; Supplement "Model fitting").
The two long-wavelength lines are blended at WINERED resolution.

Reported metrics (this module reproduces each):

  red_depth   excess absorption [%] at the position of the blended
              long-wavelength pair = maximum of (G1 + G2)
  blue_depth  amplitude [%] of the short-wavelength Gaussian (A0)
  red_blue    ratio of the blended red PEAK amplitude to the blue peak
              amplitude = red_depth / A0   (paper value 6.7; optically thin
              expectation 8)
  fwhm_A      full width at half maximum [A] of the blended red feature
              (G1 + G2), measured numerically on the fitted model
  shift_A     shared Doppler shift [A] relative to the rest wavelengths

Usage as a module:

  from he_line_metrics import fit_metrics
  m = fit_metrics(wavelength_vacuum_A, excess_percent)

`python3 he_line_metrics.py <file>` runs on a two-or-more-column text file
(col 1 vacuum wavelength [A]; excess [%] taken from the LAST column, which
in the oracle `tspec_*.txt` files is the instrument-convolved excess).
`python3 he_line_metrics.py --selftest` validates the extractor on synthetic
profiles with known answers.
"""
import sys
import numpy as np
from scipy.optimize import curve_fit

# He I 10830 triplet rest wavelengths [A].  Vacuum: Cherubim 2026, citing
# Drake 1996.  Air: the values EXHALE_transit.py / exhale_transit_lib.py use
# (l_He3_3/2/1 there, increasing wavelength here).
LAM_VAC = np.array([10832.057, 10833.217, 10833.306])
LAM_AIR = np.array([10829.09114, 10830.25010, 10830.33977])


def _three_gauss(w, a0, a1, a2, sigma, shift, lam):
    g = np.zeros_like(w)
    for a, l in zip((a0, a1, a2), lam):
        g += a*np.exp(-0.5*((w - l - shift)/sigma)**2)
    return g


def fit_metrics(wl_A, excess_pct, p0_sigma=0.4, frame='vacuum'):
    """Fit the paper's three-Gaussian model and return the metric dict.

    frame: 'vacuum' (default; the paper's convention) or 'air' (the
    EXHALE_transit.py wavelength convention).  Reported widths and the
    shift are in the input frame; the metric definitions are identical.
    """
    LAM_VAC_LOCAL = {'vacuum': LAM_VAC, 'air': LAM_AIR}[frame]
    w, y = np.asarray(wl_A, float), np.asarray(excess_pct, float)
    a_guess = max(y.max(), 1e-3)
    p0 = [0.2*a_guess, 0.4*a_guess, 0.6*a_guess, p0_sigma, 0.0]
    bounds = ([0, 0, 0, 0.05, -2.0], [np.inf, np.inf, np.inf, 5.0, 2.0])
    model = lambda ww, a0, a1, a2, sigma, shift: \
        _three_gauss(ww, a0, a1, a2, sigma, shift, LAM_VAC_LOCAL)
    popt, pcov = curve_fit(model, w, y, p0=p0, bounds=bounds)
    a0, a1, a2, sigma, shift = popt

    # blended red feature on a fine grid
    wf = np.linspace(LAM_VAC_LOCAL[1] - 5, LAM_VAC_LOCAL[2] + 5, 20001)
    red = (a1*np.exp(-0.5*((wf - LAM_VAC_LOCAL[1] - shift)/sigma)**2) +
           a2*np.exp(-0.5*((wf - LAM_VAC_LOCAL[2] - shift)/sigma)**2))
    red_depth = red.max()
    half = red_depth/2.0
    above = wf[red >= half]
    fwhm = above[-1] - above[0] if len(above) > 1 else np.nan

    return dict(red_depth=red_depth, blue_depth=a0,
                red_blue=red_depth/a0 if a0 > 0 else np.inf,
                fwhm_A=fwhm, shift_A=shift, sigma_A=sigma,
                amplitudes=(a0, a1, a2))


def _selftest():
    """Synthetic profiles with known answers."""
    w = np.linspace(10826.0, 10838.0, 1200)
    ok = True

    # 1) pure three-Gaussian input -> exact recovery
    truth = dict(a0=0.25, a1=0.45, a2=0.85, sigma=0.35, shift=0.07)
    y = _three_gauss(w, *truth.values(), LAM_VAC)
    m = fit_metrics(w, y)
    # expected blended-red peak measured from the same model
    wf = np.linspace(10828.0, 10838.0, 20001)
    red = (truth['a1']*np.exp(-0.5*((wf-LAM_VAC[1]-truth['shift'])/truth['sigma'])**2)
           + truth['a2']*np.exp(-0.5*((wf-LAM_VAC[2]-truth['shift'])/truth['sigma'])**2))
    exp_red = red.max()
    for name, got, exp, tol in [
            ('red_depth', m['red_depth'], exp_red, 1e-6),
            ('blue_depth', m['blue_depth'], truth['a0'], 1e-6),
            ('shift_A', m['shift_A'], truth['shift'], 1e-6),
            ('sigma_A', m['sigma_A'], truth['sigma'], 1e-6)]:
        good = abs(got - exp) < tol
        ok &= good
        print(f"  [{'OK' if good else 'FAIL'}] {name}: {got:.6f} vs {exp:.6f}")

    # FWHM sanity: single-dominant-line limit -> 2.3548 sigma
    y2 = _three_gauss(w, 0.0, 0.0, 1.0, 0.30, 0.0, LAM_VAC)
    m2 = fit_metrics(w, y2)
    exp_fwhm = 2.0*np.sqrt(2.0*np.log(2.0))*0.30
    good = abs(m2['fwhm_A'] - exp_fwhm) < 1e-3
    ok &= good
    print(f"  [{'OK' if good else 'FAIL'}] fwhm single-line: "
          f"{m2['fwhm_A']:.4f} vs {exp_fwhm:.4f}")

    # 2) noisy realization -> recovery within tolerance
    rng = np.random.default_rng(42)
    yn = y + rng.normal(0, 0.02, w.size)
    mn = fit_metrics(w, yn)
    good = abs(mn['red_depth'] - exp_red) < 0.05
    ok &= good
    print(f"  [{'OK' if good else 'FAIL'}] noisy red_depth: "
          f"{mn['red_depth']:.3f} vs {exp_red:.3f} (tol 0.05)")

    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
        sys.exit(_selftest())
    dat = np.loadtxt(sys.argv[1])
    m = fit_metrics(dat[:, 0], dat[:, -1])
    print(f"file: {sys.argv[1]}")
    print(f"  red_depth  = {m['red_depth']:.3f} %   (obs 1.24 +0.22/-0.23)")
    print(f"  blue_depth = {m['blue_depth']:.3f} %   (obs 0.25 +0.14/-0.12)")
    print(f"  red/blue   = {m['red_blue']:.2f}      (obs 6.7 +12.7/-3.1)")
    print(f"  FWHM       = {m['fwhm_A']:.3f} A   (obs 0.86 +0.15/-0.27)")
    print(f"  shift      = {m['shift_A']:+.3f} A   (obs +0.072)")
