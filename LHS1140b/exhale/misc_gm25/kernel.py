#!/usr/bin/env python3
"""The added Gaussian broadening one He I 10830 curve needs to reach the
measured red-pair width, and the line metrics before and after it.

Same construction as ../../make_memo_figures.py (`matched_kernel`): the
three-Gaussian extractor of ../../../he_line_metrics.py on the excess
absorption, the red-pair equivalent width over the vacuum window
10832.60-10834.20 A, and the measured width FWHM = 0.841 A.

usage: ./kernel.py <label>=<tpm_He10830.txt> ...
"""
import sys, os
import numpy as np
from scipy.optimize import brentq
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00')
from he_line_metrics import fit_metrics

AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20
C_KMS = 2.99792458e5
FWHM_OBS = 0.841


def curve(p):
    s = np.loadtxt(p)
    return s[:, 0], (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100


def broaden(lam, exc, f, lam0=10830.0):
    if f <= 0:
        return exc
    dl = np.median(np.diff(lam))
    return gaussian_filter1d(exc, f/C_KMS*lam0/2.35482/dl)


def metrics(lam, exc, f):
    e = broaden(lam, exc, f)
    d = fit_metrics(lam, e, frame='air')
    lv = lam*AIR
    m = (lv >= EW_LO) & (lv <= EW_HI)
    d['ew'] = np.trapz(e[m], lv[m])
    return d


def matched_kernel(lam, exc):
    return brentq(lambda f: metrics(lam, exc, f)['fwhm_A'] - FWHM_OBS,
                  1.0, 60.0, xtol=1e-3)


if __name__ == '__main__':
    for a in sys.argv[1:]:
        lab, path = a.split('=', 1)
        if not os.path.isfile(path):
            print('%-28s MISSING %s' % (lab, path))
            continue
        lam, exc = curve(path)
        d0 = metrics(lam, exc, 0.0)
        f = matched_kernel(lam, exc)
        d = metrics(lam, exc, f)
        print('%-28s kernel=%6.2f km/s (sigma %5.2f) | as solved red=%.4f '
              'blue=%.4f FWHM=%.4f EW=%.4f | at kernel red=%.4f blue=%.4f '
              'ratio=%.2f EW=%.4f'
              % (lab, f, f/2.35482, d0['red_depth'], d0['blue_depth'],
                 d0['fwhm_A'], d0['ew'], d['red_depth'], d['blue_depth'],
                 d['red_blue'], d['ew']))
