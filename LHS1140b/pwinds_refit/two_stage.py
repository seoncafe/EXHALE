#!/usr/bin/env python3
"""Two-stage retrieval: match the red-component equivalent width first, then
fit the line shape with the non-thermal broadening free.

Stage 1 -- the equivalent width selects the absorbing column.  The equivalent
width of the blended red pair is one number, so with the broadening held at the
Cherubim et al. (2026) value it defines a surface in (Mdot, T, h_fraction)
rather than a point.  The surface is mapped by fixing (T, h_fraction) on a grid
and solving for the log Mdot that reproduces the measured equivalent width.

Stage 2 -- the shape.  On that surface the non-thermal dispersion v_nt is
freed.  Because the line is saturated the equivalent width itself grows with
v_nt, so Mdot is re-solved at every (T, h_fraction, v_nt) node: the equivalent
width stays pinned to the measurement and only the shape is being fitted.  The
residual coupling between the two stages is exactly the amount by which Mdot
has to move to keep the equivalent width fixed as v_nt grows, and it is
reported.

Two grids:

  authors_prior   T and h_fraction inside the support of the authors' own
                  grid search (T in [4500, 7500] K, log h in [-4.5, -1.5],
                  i.e. He/H >= 30.6);
  extended        the same construction with the composition axis carried to
                  He/H ~ 0.8 and the temperature axis down to 4000 K, to see
                  where the fit goes when their prior no longer bounds it.

Outputs: two_stage_<grid>_stage{1,2}.txt and a summary on stdout.

Usage: OMP_NUM_THREADS=1 python3 two_stage.py [grid] [nproc]
"""
import sys
import time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import brentq, minimize_scalar

import refit_lib as R

GRID = sys.argv[1] if len(sys.argv) > 1 else 'authors_prior'
NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 34

AIR2VAC = 10832.057 / 10829.09114
# Red aperture: the blended long-wavelength pair (10833.217, 10833.306 A in
# vacuum), in the aperture used by ../pwinds_oracle (air 10829.6-10831.2 A).
RED_LO_VAC, RED_HI_VAC = 10829.6 * AIR2VAC, 10831.2 * AIR2VAC

X, Y, YERR = R.load_data()
KERNEL = R.instrument_kernel(X)
WL_VAC = X * 1e10 * AIR2VAC
DLAM = np.median(np.diff(WL_VAC))
RED = (WL_VAC >= RED_LO_VAC) & (WL_VAC <= RED_HI_VAC)

# C26 published medians
MDOT_PUB, T_PUB, H_PUB = 2.03e8, 5160.0, 1.01e-3


def equivalent_width(flux):
    """Red-component equivalent width [A]: integral of (1 - F) over the red
    aperture, on the data grid."""
    return np.sum((1.0 - flux)[RED]) * DLAM


EW_OBS = equivalent_width(Y)
# The error bars are those of individual pixels; the released spectrum also
# carries a GP component
# whose correlations are not propagated here, so this is a lower bound on the
# uncertainty of the equivalent width.
EW_ERR = np.sqrt(np.sum(YERR[RED] ** 2)) * DLAM


def model_flux(log_mdot, T, h_fraction, v_nt):
    return R.cascading_model((log_mdot, np.log10(T), h_fraction, v_nt),
                             X, KERNEL)


def chi2_of(flux):
    return np.sum((Y - flux) ** 2 / YERR ** 2)


def solve_mdot(T, h_fraction, v_nt, lo=7.0, hi=10.0, n_coarse=13):
    """Mass-loss rates that reproduce the observed red equivalent width.

    The equivalent width is not monotonic in Mdot -- it rises, peaks, and falls
    again as the denser wind loses its metastable fraction -- so the measured
    value is typically reached twice, on either side of the peak, and the two
    solutions can fit the line profile very differently.  All crossings are
    returned and the caller picks the one with the better full-spectrum chi^2.

    Returns (list of log Mdot, matched) where ``matched`` says whether the
    observed equivalent width is attainable at all at this node; if it is not,
    the single returned value is the Mdot of closest approach.
    """
    grid = np.linspace(lo, hi, n_coarse)
    ew = np.array([equivalent_width(model_flux(g, T, h_fraction, v_nt))
                   for g in grid])
    d = ew - EW_OBS
    crossings = np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]
    if len(crossings):
        roots = []
        for i in crossings:
            try:
                roots.append(brentq(lambda g: equivalent_width(
                    model_flux(g, T, h_fraction, v_nt)) - EW_OBS,
                    grid[i], grid[i + 1], xtol=1e-4, rtol=1e-8, maxiter=60))
            except Exception:
                roots.append(grid[i])
        return roots, True
    j = int(np.argmin(np.abs(d)))
    a = grid[max(j - 1, 0)]
    b = grid[min(j + 1, n_coarse - 1)]
    res = minimize_scalar(
        lambda g: abs(equivalent_width(model_flux(g, T, h_fraction, v_nt))
                      - EW_OBS), bounds=(a, b), method='bounded',
        options={'xatol': 1e-3})
    return [res.x], False


def node(args):
    T, h, v_nt = args
    try:
        roots, matched = solve_mdot(T, h, v_nt)
    except Exception:
        return T, h, v_nt, np.nan, np.nan, np.nan, 0
    best = None
    for lm in roots:
        flux = model_flux(lm, T, h, v_nt)
        c2 = chi2_of(flux)
        if best is None or c2 < best[1]:
            best = (lm, c2, equivalent_width(flux))
    return T, h, v_nt, best[0], best[1], best[2], int(matched)


HEADER = ('%s\n'
          'T [K]  log10 h_fraction  v_nt [km/s]  log10 Mdot (EW-matched)  '
          'chi2 (full spectrum)  EW [A]  reached (1 = the observed EW is '
          'attainable at this node)')


def pack(res):
    return np.array([[r[0], np.log10(r[1]), r[2] / 1e3, r[3], r[4], r[5], r[6]]
                     for r in res])


def main():
    print(f'red aperture (vacuum): {RED_LO_VAC:.3f}-{RED_HI_VAC:.3f} A, '
          f'{RED.sum()} pixels of {len(X)}')
    print(f'observed red equivalent width: {EW_OBS:.6f} +/- {EW_ERR:.6f} A '
          f'({EW_OBS/EW_ERR:.1f} sigma)')
    f_pub = model_flux(np.log10(MDOT_PUB), T_PUB, H_PUB, 0.0)
    print(f'C26 published vector: EW {equivalent_width(f_pub):.6f} A '
          f'({(equivalent_width(f_pub)-EW_OBS)/EW_ERR:+.2f} sigma), '
          f'chi2 {chi2_of(f_pub):.2f}')

    if GRID == 'authors_prior':
        T_grid = np.array([4500., 4800., 5160., 5500., 6000., 6500., 7000.])
        logh_grid = np.arange(-4.5, -1.49, 0.25)
        v_grid = np.arange(0.0, 18.001, 1.5) * 1e3
    elif GRID == 'extended':
        T_grid = np.array([4000., 4300., 4600., 4800., 5160., 5500., 6000.,
                           6500.])
        logh_grid = np.arange(-4.5, -0.24, 0.5)
        v_grid = np.array([0., 1.5, 3., 4.5, 6., 7.5, 9., 10.5, 12., 15.,
                           18.]) * 1e3
    else:
        raise SystemExit(f'unknown grid {GRID}')
    h_grid = 10 ** logh_grid

    t0 = time.time()
    # ---- stage 1: the broadening at the C26 value ----------------------
    jobs = [(T, h, 0.0) for T in T_grid for h in h_grid]
    with Pool(NPROC) as pool:
        res1 = pool.map(node, jobs)
    arr1 = pack(res1)
    np.savetxt(f'two_stage_{GRID}_stage1.txt', arr1, header=HEADER % 'stage 1: '
               'equivalent width matched at the C26 broadening (v_nt = 0)')
    ok = np.isfinite(arr1[:, 4])
    b1 = arr1[ok][np.argmin(arr1[ok][:, 4])]
    print()
    print('=== stage 1: equivalent-width-matched family, v_nt = 0 ===')
    print(f'{int(arr1[:, 6].sum())} of {len(arr1)} nodes reach the observed '
          f'equivalent width inside log Mdot = [7, 10]')
    print(f'best full-spectrum chi2 on the family: T {b1[0]:.0f} K, '
          f'log h {b1[1]:.2f} (He/H {(1-10**b1[1])/10**b1[1]:.0f}), '
          f'log Mdot {b1[3]:.3f}, EW {b1[5]:.6f} A '
          f'({(b1[5]-EW_OBS)/EW_ERR:+.1f} sigma), chi2 {b1[4]:.2f}')
    print('the family at T = 5160 K:')
    sel = arr1[np.isclose(arr1[:, 0], 5160.)]
    print(f'{"log h":>8s} {"He/H":>9s} {"log Mdot":>9s} {"EW [A]":>10s} '
          f'{"EW-EW_obs":>10s} {"matched":>8s} {"chi2":>9s}')
    for row in sel:
        h = 10 ** row[1]
        print(f'{row[1]:8.2f} {(1-h)/h:9.1f} {row[3]:9.3f} {row[5]:10.6f} '
              f'{(row[5]-EW_OBS)/EW_ERR:+10.2f} {int(row[6]):8d} '
              f'{row[4]:9.2f}')

    # ---- stage 2: free the broadening on the EW-matched family ---------
    jobs = [(T, h, v) for T in T_grid for h in h_grid for v in v_grid]
    with Pool(NPROC) as pool:
        res2 = pool.map(node, jobs)
    arr2 = pack(res2)
    np.savetxt(f'two_stage_{GRID}_stage2.txt', arr2, header=HEADER % 'stage 2: '
               'equivalent width kept matched while v_nt is free')
    ok = np.isfinite(arr2[:, 4])
    b2 = arr2[ok][np.argmin(arr2[ok][:, 4])]
    print()
    print('=== stage 2: shape fitted with v_nt free, EW held ===')
    print(f'best: T {b2[0]:.0f} K, log h {b2[1]:.2f} '
          f'(He/H {(1-10**b2[1])/10**b2[1]:.0f}), log Mdot {b2[3]:.3f} '
          f'(Mdot {10**b2[3]:.3g} g/s), v_nt {b2[2]:.2f} km/s, '
          f'EW {b2[5]:.6f} A ({(b2[5]-EW_OBS)/EW_ERR:+.1f} sigma), '
          f'chi2 {b2[4]:.2f}')
    print(f'stage 1 -> stage 2: delta chi2 = {b1[4] - b2[4]:.2f}')

    print()
    print('profile of the best chi2 against v_nt (minimized over T and h):')
    print(f'{"v_nt [km/s]":>12s} {"chi2":>9s} {"T [K]":>7s} {"log h":>7s} '
          f'{"He/H":>8s} {"log Mdot":>9s} {"EW-EW_obs":>10s}')
    for v in np.unique(arr2[:, 2]):
        s = arr2[(arr2[:, 2] == v) & np.isfinite(arr2[:, 4])]
        if len(s) == 0:
            continue
        r = s[np.argmin(s[:, 4])]
        h = 10 ** r[1]
        print(f'{v:12.2f} {r[4]:9.2f} {r[0]:7.0f} {r[1]:7.2f} '
              f'{(1-h)/h:8.1f} {r[3]:9.3f} {(r[5]-EW_OBS)/EW_ERR:+10.2f}')

    print()
    print('residual coupling: Mdot needed to hold the equivalent width, '
          'as v_nt grows (T = 5160 K, log h = -3.00)')
    s = arr2[np.isclose(arr2[:, 0], 5160.) & np.isclose(arr2[:, 1], -3.0)]
    if len(s) == 0:                     # -3.00 is not on the 0.25 grid
        j = np.argmin(np.abs(logh_grid + 3.0))
        s = arr2[np.isclose(arr2[:, 0], 5160.)
                 & np.isclose(arr2[:, 1], logh_grid[j])]
    print(f'{"v_nt [km/s]":>12s} {"log Mdot":>9s} {"Mdot/Mdot(0)":>13s} '
          f'{"chi2":>9s}')
    ref = None
    for row in s[np.argsort(s[:, 2])]:
        if not np.isfinite(row[3]):
            continue
        if ref is None:
            ref = row[3]
        print(f'{row[2]:12.2f} {row[3]:9.3f} {10**(row[3]-ref):13.3f} '
              f'{row[4]:9.2f}')
    print(f'\nelapsed {(time.time()-t0)/60:.1f} min')


if __name__ == '__main__':
    main()
