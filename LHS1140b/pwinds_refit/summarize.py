#!/usr/bin/env python3
"""Tables and figures for the free-broadening re-retrieval.

Reads the two-stage grids (two_stage_<grid>_stage{1,2}.txt) and, if present,
the joint control grid (joint_grid.npz), and writes summary.log content to
stdout plus fig1_spectrum_fits.pdf and fig2_broadening.pdf.

Usage: OMP_NUM_THREADS=1 python3 summarize.py
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00')
import refit_lib as R
from he_line_metrics import fit_metrics

CLIGHT = 2.99792458e5
LAM0_VAC = 10832.057
AIR2VAC = 10832.057 / 10829.09114

# Cherubim et al. (2026), GJ 1132 SED: value, +1 sigma, -1 sigma
C26 = {'mdot': (2.03e8, 0.67e8, 0.58e8), 'T': (5160., 46., 50.),
       'log_h': (np.log10(1.01e-3), 0.261, 0.291)}

X, Y, YERR = R.load_data()
KERNEL = R.instrument_kernel(X)
WL_VAC = X * 1e10 * AIR2VAC
DLAM = np.median(np.diff(WL_VAC))
RED = (WL_VAC >= 10829.6 * AIR2VAC) & (WL_VAC <= 10831.2 * AIR2VAC)
EW_OBS = np.sum((1 - Y)[RED]) * DLAM
EW_ERR = np.sqrt(np.sum(YERR[RED] ** 2)) * DLAM

COLS = 'T log_h v_nt log_mdot chi2 ew reached'.split()


def load(grid, stage):
    f = f'two_stage_{grid}_stage{stage}.txt'
    return np.loadtxt(f) if os.path.exists(f) else None


def he_h(log_h):
    h = 10 ** log_h
    return (1 - h) / h


def model_of(T, log_h, v_nt, log_mdot):
    return R.cascading_model((log_mdot, np.log10(T), 10 ** log_h, v_nt),
                             X, KERNEL)


def describe(row, tag):
    T, lh, v, lm, c2, ew = row[:6]
    m = model_of(T, lh, v * 1e3, lm)
    d = fit_metrics(WL_VAC, (1 - m) * 100, frame='vacuum')
    print(f'{tag:26s} T {T:6.0f}  He/H {he_h(lh):8.1f}  '
          f'Mdot {10**lm:9.3g}  v_nt {v:5.2f}  chi2 {c2:8.2f}  '
          f'chi2/dof {c2/(len(X)-4):6.3f}  EW {ew:.6f} '
          f'({(ew-EW_OBS)/EW_ERR:+5.1f} s)  red {d["red_depth"]:5.3f}%  '
          f'FWHM {d["fwhm_A"]:5.3f} A = {d["fwhm_A"]/LAM0_VAC*CLIGHT:5.2f} km/s')
    return m, d


def sigma_shift(val, ref, ref_p, ref_m):
    return (val - ref) / (ref_p if val > ref else ref_m)


def main():
    print(f'data: {len(X)} pixels; red aperture {RED.sum()} pixels')
    print(f'observed red equivalent width {EW_OBS:.6f} +/- {EW_ERR:.6f} A '
          '(errors of individual pixels only)')
    print(f'chi2 of a flat continuum: {np.sum((Y-1)**2/YERR**2):.1f}')
    print()

    th_pub = (np.log10(C26['mdot'][0]), np.log10(C26['T'][0]), 1.01e-3)
    m_pub = R.cascading_model(th_pub, X, KERNEL)
    d_pub = fit_metrics(WL_VAC, (1 - m_pub) * 100, frame='vacuum')
    chi2_pub = np.sum((Y - m_pub) ** 2 / YERR ** 2)
    ew_pub = np.sum((1 - m_pub)[RED]) * DLAM
    print('=== reference points ===')
    print(f'{"C26 published vector":26s} T {C26["T"][0]:6.0f}  '
          f'He/H {he_h(C26["log_h"][0]):8.1f}  Mdot {C26["mdot"][0]:9.3g}  '
          f'v_nt  0.00  chi2 {chi2_pub:8.2f}  chi2/dof '
          f'{chi2_pub/(len(X)-3):6.3f}  EW {ew_pub:.6f} '
          f'({(ew_pub-EW_OBS)/EW_ERR:+5.1f} s)  red {d_pub["red_depth"]:5.3f}%'
          f'  FWHM {d_pub["fwhm_A"]:5.3f} A = '
          f'{d_pub["fwhm_A"]/LAM0_VAC*CLIGHT:5.2f} km/s')
    d_obs = fit_metrics(WL_VAC, (1 - Y) * 100, frame='vacuum')
    print(f'{"observed":26s} {"":41s}                      '
          f'EW {EW_OBS:.6f}            red {d_obs["red_depth"]:5.3f}%  '
          f'FWHM {d_obs["fwhm_A"]:5.3f} A = '
          f'{d_obs["fwhm_A"]/LAM0_VAC*CLIGHT:5.2f} km/s')

    best = {}
    spec = {}
    print()
    print('=== two-stage results ===')
    for grid in ('authors_prior', 'extended'):
        a1, a2 = load(grid, 1), load(grid, 2)
        if a1 is None or a2 is None:
            continue
        r1 = a1[np.argmin(a1[:, 4])]
        r2 = a2[np.argmin(a2[:, 4])]
        best[grid] = {'stage1': r1, 'stage2': r2}
        spec[grid] = {}
        spec[grid]['stage1'] = describe(r1, f'{grid}: stage 1')[0]
        spec[grid]['stage2'] = describe(r2, f'{grid}: stage 2')[0]
        print(f'{"":26s} stage 1 -> stage 2: delta chi2 = '
              f'{r1[4]-r2[4]:.1f}')

    print()
    print('=== displacement from the published values (two-stage best) ===')
    for grid in best:
        r = best[grid]['stage2']
        print(f'{grid:14s} Mdot {10**r[3]:9.3g} '
              f'({sigma_shift(10**r[3], *C26["mdot"]):+5.2f} sigma_pub)  '
              f'T {r[0]:6.0f} ({sigma_shift(r[0], *C26["T"]):+6.1f} '
              f'sigma_pub)  log h {r[1]:6.2f} '
              f'({sigma_shift(r[1], *C26["log_h"]):+5.2f} sigma_pub)')

    # composition profile on the equivalent-width-held surface
    for grid in best:
        a2 = load(grid, 2)
        print()
        print(f'=== {grid}: composition profile, equivalent width held ===')
        print(f'{"log h":>7s} {"He/H":>10s} {"chi2 (v_nt=0)":>14s} '
              f'{"chi2 (v_nt free)":>17s} {"v_nt":>6s} {"T":>7s} '
              f'{"log Mdot":>9s}')
        for lh in np.unique(a2[:, 1]):
            s = a2[a2[:, 1] == lh]
            r = s[np.argmin(s[:, 4])]
            s0 = s[s[:, 2] == 0]
            c0 = s0[:, 4].min() if len(s0) else np.nan
            print(f'{lh:7.2f} {he_h(lh):10.1f} {c0:14.2f} {r[4]:17.2f} '
                  f'{r[2]:6.2f} {r[0]:7.0f} {r[3]:9.3f}')

        print()
        print(f'=== {grid}: v_nt profile, equivalent width held ===')
        p = prof_v(a2)
        print(f'{"v_nt":>6s} {"chi2":>9s} {"T":>7s} {"He/H":>10s} '
              f'{"log Mdot":>9s}')
        for r in p:
            print(f'{r[2]:6.2f} {r[4]:9.2f} {r[0]:7.0f} {he_h(r[1]):10.1f} '
                  f'{r[3]:9.3f}')

        print()
        print(f'=== {grid}: residual coupling, Mdot that holds the '
              f'equivalent width (T = 5160 K) ===')
        print(f'{"v_nt":>6s}' + ''.join(
            f'{"He/H=" + f"{he_h(lh):.0f}":>14s}' for lh in (-4.0, -3.0, -2.0)))
        vs = np.unique(a2[:, 2])
        for v in vs:
            line = f'{v:6.2f}'
            for lh in (-4.0, -3.0, -2.0):
                s = a2[np.isclose(a2[:, 1], lh) & np.isclose(a2[:, 0], 5160.)
                       & (a2[:, 2] == v)]
                line += f'{10**s[0, 3]:14.3g}' if len(s) else f'{"-":>14s}'
            print(line)

    # joint control grid
    if os.path.exists('joint_grid.npz'):
        g = np.load('joint_grid.npz')
        print()
        print('=== joint control grid ===')
        for tag, c in (('fixed broadening', g['chi2_fixed']),
                       ('free broadening', g['chi2_free'])):
            i = np.unravel_index(np.argmin(c), c.shape)
            v = g['v_nt'][i[3]] / 1e3 if c.shape[3] > 1 else 0.0
            print(f'{tag:18s} chi2 {c[i]:8.2f}  Mdot '
                  f'{10**g["log_mdot"][i[0]]:9.3g}  T {g["T"][i[1]]:6.0f}  '
                  f'He/H {he_h(g["log_h"][i[2]]):8.1f}  v_nt {v:5.2f}')
            # restricted to the authors' composition prior
            mask = g['log_h'] <= -1.5
            cc = c[:, :, mask, :]
            j = np.unravel_index(np.argmin(cc), cc.shape)
            vv = g['v_nt'][j[3]] / 1e3 if cc.shape[3] > 1 else 0.0
            print(f'{"  within log h <= -1.5":18s} chi2 {cc[j]:8.2f}  Mdot '
                  f'{10**g["log_mdot"][j[0]]:9.3g}  T {g["T"][j[1]]:6.0f}  '
                  f'He/H {he_h(g["log_h"][mask][j[2]]):8.1f}  v_nt {vv:5.2f}')
        print()
        print('joint grid, profile chi2 along each axis (free broadening; '
              'fixed broadening in brackets)')
        cf, cv = g['chi2_fixed'], g['chi2_free']
        for lab, vals, ax in (('log h / He-H', g['log_h'], 2),
                              ('T [K]', g['T'], 1),
                              ('log Mdot', g['log_mdot'], 0),
                              ('v_nt [km/s]', g['v_nt'] / 1e3, 3)):
            print(f'  {lab}')
            oth = tuple(i for i in range(4) if i != ax)
            pv = np.min(cv, axis=oth)
            pf = (np.min(cf, axis=oth) if ax != 3
                  else np.full_like(pv, np.min(cf)))
            for k, x in enumerate(vals):
                extra = (f'  (He/H {he_h(x):9.1f})' if ax == 2 else '')
                print(f'    {x:9.3f}  chi2 free {pv[k]:9.2f}   '
                      f'fixed {pf[k]:9.2f}{extra}')

    make_figures(best, spec, m_pub)


def prof_v(a2):
    """chi2 minimized over T and h at each v_nt, with the best node."""
    vs = np.unique(a2[:, 2])
    out = []
    for v in vs:
        s = a2[a2[:, 2] == v]
        out.append(s[np.argmin(s[:, 4])])
    return np.array(out)


def make_figures(best, spec, m_pub):
    plt.rcParams.update({'font.size': 9})
    c_pub, c_s1, c_s2 = '0.45', '#1f77b4', '#d62728'

    # ---- Figure 1 -------------------------------------------------------
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.6, 4.8), sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.07})
    ax.errorbar(WL_VAC, (Y - 1) * 100, yerr=YERR * 100, fmt='o', ms=2.4,
                color='0.3', elinewidth=0.7, capsize=0,
                label=r'WINERED 2024B, in transit (GP corrected)')
    ax.plot(WL_VAC, (m_pub - 1) * 100, '-', lw=1.5, color=c_pub,
            label=r'C26 published vector (fixed broadening)')
    g = 'extended' if 'extended' in spec else 'authors_prior'
    ax.plot(WL_VAC, (spec[g]['stage1'] - 1) * 100, '--', lw=1.5, color=c_s1,
            label=r'stage 1: equivalent width matched, $v_{\rm nt}=0$')
    ax.plot(WL_VAC, (spec[g]['stage2'] - 1) * 100, '-', lw=1.8, color=c_s2,
            label=r'stage 2: $v_{\rm nt}$ free, equivalent width held')
    ax.axhline(0, color='0.75', lw=0.6)
    ax.set_ylabel(r'excess absorption [\%]')
    ax.legend(frameon=False, fontsize=8, loc='lower left')
    for m, c in ((m_pub, c_pub), (spec[g]['stage1'], c_s1),
                 (spec[g]['stage2'], c_s2)):
        axr.plot(WL_VAC, (Y - m) / YERR, '-', lw=1.0, color=c)
    axr.axhline(0, color='0.75', lw=0.6)
    axr.set_ylabel(r'residual [$\sigma$]')
    axr.set_xlabel(r'vacuum wavelength [\AA]')
    axr.ticklabel_format(axis='x', useOffset=False, style='plain')
    fig.savefig('fig1_spectrum_fits.pdf', bbox_inches='tight')
    fig.savefig('fig1_spectrum_fits.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # ---- Figure 2 -------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2))
    style = {'authors_prior': ('#1f77b4', '-', r"inside the C26 prior"),
             'extended': ('#d62728', '-', r'composition axis extended')}

    ax = axes[0]
    for grid in best:
        a2 = load(grid, 2)
        p = prof_v(a2)
        col, ls, lab = style[grid]
        ax.plot(p[:, 2], p[:, 4], ls, color=col, marker='o', ms=3, label=lab)
    ax.axvline(9.33, color='k', ls='--', lw=1.0)
    ax.annotate(r'width memo: 9.33', xy=(9.33, 40), xytext=(3, 0),
                textcoords='offset points', fontsize=7.5, rotation=90,
                va='bottom', ha='left')
    ax.set_xlabel(r'$v_{\rm nt}$ [km s$^{-1}$]')
    ax.set_ylabel(r'$\chi^2$ (equivalent width held)')
    ax.set_yscale('log')
    ax.legend(frameon=False, fontsize=7.5)
    ax.set_title(r'(a) where the broadening goes', fontsize=9)

    # (b) the composition axis, with and without the free broadening
    ax = axes[1]
    a2 = load('extended', 2)
    if a2 is None:
        a2 = load('authors_prior', 2)
    lhs = np.unique(a2[:, 1])
    free = np.array([np.min(a2[a2[:, 1] == lh][:, 4]) for lh in lhs])
    fixed = np.array([np.min(a2[(a2[:, 1] == lh) & (a2[:, 2] == 0)][:, 4])
                      for lh in lhs])
    ax.plot(he_h(lhs), fixed, '-o', ms=3, color='#1f77b4',
            label=r'$v_{\rm nt}=0$ (C26 broadening)')
    ax.plot(he_h(lhs), free, '-o', ms=3, color='#d62728',
            label=r'$v_{\rm nt}$ free')
    ax.axvline(he_h(C26['log_h'][0]), color='k', ls=':', lw=1.0)
    ax.axvspan(he_h(C26['log_h'][0] + C26['log_h'][1]),
               he_h(C26['log_h'][0] - C26['log_h'][2]), color='k', alpha=0.10)
    ax.annotate(r'C26 He/H', xy=(he_h(C26['log_h'][0]), 700), xytext=(3, 0),
                textcoords='offset points', fontsize=7.5, rotation=90,
                va='top', ha='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'He/H number ratio')
    ax.set_ylabel(r'$\chi^2$ (equivalent width held)')
    ax.legend(frameon=False, fontsize=7.5, loc='upper center')
    ax.set_title(r'(b) what is left of the composition', fontsize=9)

    # (c) residual coupling between the two stages
    ax = axes[2]
    for lh, col, lab in ((-3.0, '#d62728', r'He/H $=999$ (C26)'),
                         (-2.0, '#1f77b4', r'He/H $=99$'),
                         (-4.0, '#2ca02c', r'He/H $=10^4$')):
        s = a2[np.isclose(a2[:, 1], lh) & np.isclose(a2[:, 0], 5160.)]
        if len(s) == 0:
            continue
        s = s[np.argsort(s[:, 2])]
        ax.plot(s[:, 2], 10 ** s[:, 3], '-o', ms=3, color=col, label=lab)
    ax.axhline(C26['mdot'][0], color='k', ls=':', lw=1.0)
    ax.set_yscale('log')
    ax.set_xlabel(r'$v_{\rm nt}$ [km s$^{-1}$]')
    ax.set_ylabel(r'$\dot{M}$ that holds the equivalent width [g s$^{-1}$]')
    ax.set_title(r'(c) coupling left between the stages, $T=5160$ K',
                 fontsize=9)
    ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig('fig2_broadening.pdf', bbox_inches='tight')
    fig.savefig('fig2_broadening.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('\nfigures written: fig1_spectrum_fits.pdf, fig2_broadening.pdf')


if __name__ == '__main__':
    main()
