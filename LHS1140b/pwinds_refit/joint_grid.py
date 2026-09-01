#!/usr/bin/env python3
"""Control experiment: joint likelihood grids with the broadening fixed and
free, in the style of the authors' own grid search.

Cherubim et al. (2026) retrieved on a 20 x 20 x 20 grid in (log Mdot, T,
log h_fraction) with the broadening set by the model.  Here the same kind of
grid is evaluated twice: once at their broadening (v_nt = 0) and once with
v_nt as a fourth axis.  Profile likelihoods along each axis then give the
maximum-likelihood parameters and Delta chi^2 = 1 intervals, which is what the
comparison with their published values needs.

Outputs: joint_grid.npz and a summary on stdout.

Usage: OMP_NUM_THREADS=1 python3 joint_grid.py [nproc]
"""
import sys
import time
import numpy as np
from multiprocessing import Pool

import refit_lib as R

NPROC = int(sys.argv[1]) if len(sys.argv) > 1 else 64

X, Y, YERR = R.load_data()
KERNEL = R.instrument_kernel(X)
R._ray_tracing()

# the authors' grid support; finer in Mdot and T than their 20 nodes over the
# full range, because the interesting region is narrow
LOG_MDOT = np.linspace(7.0, 10.0, 25)
T_VALS = np.array([4500., 4700., 4900., 5100., 5300., 5500., 5800., 6100.,
                   6500., 7000., 7500.])
LOG_H = np.arange(-4.5, -0.49, 0.25)   # carried past the authors' -1.5 edge
V_NT = np.array([0., 2., 4., 6., 8., 9., 10., 11., 12., 14., 16., 20.]) * 1e3


def chi2_of(args):
    lm, T, lh, v = args
    m = R.cascading_model((lm, np.log10(T), 10 ** lh, v), X, KERNEL)
    return np.sum((Y - m) ** 2 / YERR ** 2)


def run(v_grid, tag):
    jobs = [(lm, T, lh, v) for lm in LOG_MDOT for T in T_VALS
            for lh in LOG_H for v in v_grid]
    t0 = time.time()
    with Pool(NPROC) as pool:
        out = pool.map(chi2_of, jobs, chunksize=4)
    c = np.array(out).reshape(len(LOG_MDOT), len(T_VALS), len(LOG_H),
                              len(v_grid))
    print(f'{tag}: {c.size} models in {(time.time()-t0)/60:.1f} min')
    return c


def profile(chi2, axis_values, axis, name, unit=''):
    """Profile chi^2 along one axis; print the minimum and Delta chi^2 = 1."""
    other = tuple(i for i in range(chi2.ndim) if i != axis)
    p = np.min(chi2, axis=other)
    k = int(np.argmin(p))
    lo, hi = None, None
    for i in range(k, -1, -1):
        if p[i] - p[k] > 1.0:
            lo = np.interp(p[k] + 1.0, [p[i + 1], p[i]],
                           [axis_values[i + 1], axis_values[i]])
            break
    for i in range(k, len(p)):
        if p[i] - p[k] > 1.0:
            hi = np.interp(p[k] + 1.0, [p[i - 1], p[i]],
                           [axis_values[i - 1], axis_values[i]])
            break
    lo_s = f'{lo:.4g}' if lo is not None else '< grid'
    hi_s = f'{hi:.4g}' if hi is not None else '> grid'
    print(f'{name:22s} best {axis_values[k]:11.4g} {unit:8s} '
          f'delta chi2 = 1 interval [{lo_s}, {hi_s}]')
    return axis_values[k], lo, hi, p


def report(chi2, v_grid, tag):
    print()
    print(f'=== {tag} ===')
    i = np.unravel_index(np.argmin(chi2), chi2.shape)
    print(f'best chi2 {chi2[i]:.3f} for {len(X)} points at '
          f'log Mdot {LOG_MDOT[i[0]]:.3f} (Mdot {10**LOG_MDOT[i[0]]:.4g}), '
          f'T {T_VALS[i[1]]:.0f} K, log h {LOG_H[i[2]]:.2f} '
          f'(He/H {(1-10**LOG_H[i[2]])/10**LOG_H[i[2]]:.0f}), '
          f'v_nt {v_grid[i[3]]/1e3:.1f} km/s')
    profile(chi2, LOG_MDOT, 0, 'log10 Mdot [g/s]')
    profile(chi2, T_VALS, 1, 'T [K]', 'K')
    profile(chi2, LOG_H, 2, 'log10 h_fraction')
    if len(v_grid) > 1:
        profile(chi2, v_grid / 1e3, 3, 'v_nt [km/s]', 'km/s')


def main():
    c_fixed = run(np.array([0.0]), 'fixed broadening')
    report(c_fixed, np.array([0.0]), 'fixed broadening (the C26 setup)')
    c_free = run(V_NT, 'free broadening')
    report(c_free, V_NT, 'free broadening')
    print()
    print(f'delta chi2 (fixed -> free) = '
          f'{np.min(c_fixed) - np.min(c_free):.2f} for one extra parameter')
    np.savez('joint_grid.npz', log_mdot=LOG_MDOT, T=T_VALS, log_h=LOG_H,
             v_nt=V_NT, chi2_fixed=c_fixed, chi2_free=c_free)


if __name__ == '__main__':
    main()
