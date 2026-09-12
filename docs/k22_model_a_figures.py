#!/usr/bin/env python3
"""Figures of docs/koskinen2022_model_a_comparison.tex: an EXHALE hot-Uranus
run against Koskinen et al. (2022, ApJ 929, 52) Model A, a = 0.05 au.

Usage: k22_model_a_figures.py <run_output_dir> [<second_run_output_dir>] [<third_run_output_dir>]
       (default: benchmarks/koskinen2022_model_a/matched_hnu_minus_I/output, the
       matched run; benchmarks/koskinen2022_model_a/output is the superseded
       Roche-domain state of the morning of 2026-09-05, see the README there)

The Model A curves are the digitized readings of docs/p23_published_profiles.md
section 5.1 (Figure 8 densities, Figure 7 T and v; log10 n in cm^-3), on
r/r_base = (r/R_p)/1.34, r_base being their 1 microbar lower boundary; ours is
r/R_p, since our R_p IS the lower boundary (same absolute radius to 2e-4).
Writes PDFs into docs/figures/k22_model_a/ (or $K22_FIG_DIR when set).
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'examples'))
import exhale_io as aio

# --- Koskinen et al. (2022) Model A, digitized (p23_published_profiles.md 5.1) ---
# r/R_p, r/r_base, log10 n [cm^-3]: H2, H, He, H+, e-, He+, H3+, H2+, HeH+, T [K], v [km/s]
NAN = np.nan
K22 = np.array([
 [1.34, 1.000, 12.69, 11.21, 11.89,  4.59, 4.82, 0.17, 4.41, -2.83, -3.39, 1067,  NAN],
 [1.40, 1.045, 11.84, 10.42, 11.08,  4.84, 6.71, 1.87, 3.43, -0.84, -1.67, 1096,  NAN],
 [1.50, 1.119, 10.86,  9.59, 10.11,   NAN, 6.72, 2.32, 3.86,  0.70, -1.47, 1611,  NAN],
 [1.60, 1.194, 10.12,  9.22,  9.36,   NAN, 6.96, 3.15, 3.94,  1.78, -0.51, 1884,  NAN],
 [1.75, 1.306,  9.31,  8.88,  8.60,   NAN, 7.01, 3.80,  NAN,  2.73,  0.14, 2286, 0.09],
 [2.00, 1.493,  8.38,  8.38,  7.77,   NAN, 6.77, 4.07, 4.03,  3.28,  0.26, 2855, 0.30],
 [2.30, 1.716,  7.82,  7.94,  7.24,   NAN, 6.48, 4.05, 3.97,  3.55,  0.22, 3517, 0.93],
 [2.70, 2.015,  7.31,  7.56,  6.81,   NAN, 6.24, 4.04, 3.86,  3.61,  0.14, 4141, 1.91],
 [3.20, 2.388,  6.92,  7.26,  6.45,  6.03, 6.07, 4.04, 3.61,  3.55,  0.09, 4488, 2.92],
 [4.00, 2.985,  6.51,  6.92,  6.09,  5.87, 5.89, 4.05, 3.29,  3.48,  0.03, 4740, 4.36],
 [5.00, 3.731,  6.13,  6.61,  5.79,  5.69, 5.73, 4.04, 3.00,  3.36,   NAN, 4924, 5.69],
 [6.00, 4.478,  5.86,  6.39,  5.53,  5.56, 5.61, 4.01, 2.77,  3.27, -0.11, 5087, 6.84],
 [7.50, 5.597,  5.55,  6.12,  5.24,  5.40, 5.43, 3.95, 2.52,  3.16, -0.23, 5292, 8.21],
 [9.00, 6.716,   NAN,  5.90,  5.03,  5.26, 5.30, 3.88, 2.34,  3.08, -0.33, 5477, 9.31],
 [9.60, 7.164,  5.19,  5.82,  4.96,  5.20, 5.25, 3.87, 2.27,  3.05, -0.37, 5527, 9.70]])
KC = {'r': 1, 'H2': 2, 'H': 3, 'He': 4, 'Hp': 5, 'e': 6, 'Hep': 7, 'H3p': 8, 'H2p': 9, 'HeHp': 10, 'T': 11, 'v': 12}

def kcol(name):
    c = K22[:, KC[name]]; m = np.isfinite(c)
    return K22[m, KC['r']], c[m]

def load(outdir):
    h = aio.load_hydro(os.path.join(outdir, 'Hydro_ioniz.txt'))
    r, ion = aio.load_ions(os.path.join(outdir, 'Ion_species.txt'))
    nH = ion['HI'] + ion['HII'] + 2*(ion['H2'] + ion['H2p']) + 3*ion['H3p'] + ion['HeHp']
    ne = ion['HII'] + ion['HeII'] + 2*ion['HeIII'] + ion['H2p'] + ion['H3p'] + ion['HeHp']
    return dict(r=r, T=h['T'], v=h['v']/1e5, H2=ion['H2'], H=ion['HI'], Hp=ion['HII'],
                He=ion['HeI'], Hep=ion['HeII'], H3p=ion['H3p'], H2p=ion['H2p'], HeHp=ion['HeHp'],
                e=ne, fH2=2*ion['H2']/nH, xHp=ion['HII']/nH, nH=nH)

def main():
    args = sys.argv[1:]
    out = args[0] if args else os.path.join(HERE, '..', 'benchmarks', 'koskinen2022_model_a',
                                            'matched_hnu_minus_I', 'output')
    A = load(out); B = load(args[1]) if len(args) > 1 else None; C = load(args[2]) if len(args) > 2 else None
    fdir = os.environ.get('K22_FIG_DIR') or os.path.join(HERE, 'figures', 'k22_model_a'); os.makedirs(fdir, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.labelsize': 11, 'legend.fontsize': 8.5})
    lab_ours = os.environ.get('K22_LABEL1', 'EXHALE'); lab_prev = os.environ.get('K22_LABEL2', 'EXHALE, second run'); lab_k = 'Koskinen et al. (2022) Model A'; lab_third = os.environ.get('K22_LABEL3', 'EXHALE, third run'); col3 = 'C1'
    rmax = 7.3

    # 1. temperature and velocity
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.4))
    for a, key, yl in ((ax[0], 'T', r'$T$ [K]'), (ax[1], 'v', r'$v$ [km s$^{-1}$]')):
        a.plot(A['r'], A[key], 'k-', lw=1.6, label=lab_ours)
        if B is not None: a.plot(B['r'], B[key], '-', color='0.6', lw=1.0, label=lab_prev)
        if C is not None: a.plot(C['r'], C[key], '--', color=col3, lw=1.4, label=lab_third)
        rk, yk = kcol(key); a.plot(rk, yk, 'o', ms=4.5, mfc='none', color='C3', label=lab_k)
        a.set_xlabel(r'$r/r_{\rm base}$'); a.set_ylabel(yl); a.set_xlim(1, rmax); a.grid(alpha=0.3)
    ax[0].set_ylim(0, 6000); ax[1].set_ylim(-0.5, 11); ax[0].legend(loc='lower right')
    fig.tight_layout(); fig.savefig(os.path.join(fdir, 'temperature_velocity.pdf')); plt.close(fig)

    # 2. number densities
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.6))
    for a, keys in ((ax[0], (('H2', r'H$_2$', 'C0'), ('H', 'H', 'C2'), ('He', 'He', 'C1'))),
                    (ax[1], (('e', r'$e^-$', 'k'), ('Hp', r'H$^+$', 'C0'), ('Hep', r'He$^+$', 'C1'),
                             ('H3p', r'H$_3^+$', 'C3'), ('H2p', r'H$_2^+$', 'C4')))):
        for key, lab, col in keys:
            a.semilogy(A['r'], A[key], '-', color=col, lw=1.5, label=lab)
            if B is not None: a.semilogy(B['r'], B[key], ':', color=col, lw=1.0)
            if C is not None: a.semilogy(C['r'], C[key], '--', color=col, lw=1.1)
            rk, yk = kcol(key); a.semilogy(rk, 10.0**yk, 'o', ms=4, mfc='none', color=col)
        a.set_xlabel(r'$r/r_{\rm base}$'); a.set_ylabel(r'$n$ [cm$^{-3}$]'); a.set_xlim(1, rmax); a.grid(alpha=0.3)
        a.legend(loc='upper right', ncol=1)
    ax[0].set_ylim(1e4, 1e13); ax[1].set_ylim(1e0, 1e8)
    ax[0].set_title('neutrals (solid: ' + lab_ours.split(',')[0] + '; circles: Model A' + ('; dotted/dashed: second/third run' if B is not None else '') + ')', fontsize=8)
    ax[1].set_title('ions and electrons', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(fdir, 'densities.pdf')); plt.close(fig)

    # 3. fractions: H nuclei in H2, and ionized H
    rk, lH2 = kcol('H2'); _, lH = kcol('H'); rH2 = rk
    # Model A f(H2) needs n(H nuclei) ~ 2 n(H2) + n(H) (+ ions, negligible where readable)
    rr = np.intersect1d(kcol('H2')[0], kcol('H')[0])
    lh2 = np.interp(rr, *kcol('H2')); lh = np.interp(rr, *kcol('H'))
    fk = 2*10**lh2/(2*10**lh2 + 10**lh)
    rk2, lHp = kcol('Hp'); lh2b = np.interp(rk2, *kcol('H2')); lhb = np.interp(rk2, *kcol('H'))
    xk = 10**lHp/(2*10**lh2b + 10**lhb + 10**lHp)
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.4))
    ax[0].plot(A['r'], A['fH2'], 'k-', lw=1.6, label=lab_ours)
    if B is not None: ax[0].plot(B['r'], B['fH2'], '-', color='0.6', lw=1.0, label=lab_prev)
    if C is not None: ax[0].plot(C['r'], C['fH2'], '--', color=col3, lw=1.4, label=lab_third)
    ax[0].plot(rr, fk, 'o', ms=4.5, mfc='none', color='C3', label=lab_k)
    ax[0].set_ylabel(r'$f({\rm H}_2) = 2n({\rm H}_2)/n_{\rm H}$'); ax[0].set_ylim(0, 1.02); ax[0].legend(loc='lower left')
    ax[1].semilogy(A['r'], A['xHp'], 'k-', lw=1.6)
    if B is not None: ax[1].semilogy(B['r'], B['xHp'], '-', color='0.6', lw=1.0)
    if C is not None: ax[1].semilogy(C['r'], C['xHp'], '--', color=col3, lw=1.4)
    ax[1].semilogy(rk2, xk, 'o', ms=4.5, mfc='none', color='C3')
    ax[1].set_ylabel(r'$x({\rm H}^+) = n({\rm H}^+)/n_{\rm H}$'); ax[1].set_ylim(1e-8, 1)
    for a in ax: a.set_xlabel(r'$r/r_{\rm base}$'); a.set_xlim(1, rmax); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(fdir, 'fractions.pdf')); plt.close(fig)
    # 3b. heating and cooling rates against the digitized Figure 9
    f9 = os.path.join(HERE, '..', 'benchmarks', 'koskinen2022_model_a', 'model_a_fig9_digitized.txt')
    if os.path.exists(f9):
        d9 = np.loadtxt(f9)
        hb = np.loadtxt(os.path.join(out, 'Heating_breakdown.txt')); cb = np.loadtxt(os.path.join(out, 'Cooling_breakdown.txt'))
        fig, a = plt.subplots(figsize=(5.2, 3.6))
        a.semilogy(hb[:, 0], hb[:, 3]*0.1, 'r-', lw=1.6, label=lab_ours + ' stellar heating')
        a.semilogy(cb[:, 0], cb[:, 3]*0.1, 'b--', lw=1.4, label=lab_ours + ' radiative cooling')
        if B is not None and len(args) > 1 and os.path.exists(os.path.join(args[1], 'Heating_breakdown.txt')):
            hb2 = np.loadtxt(os.path.join(args[1], 'Heating_breakdown.txt')); cb2 = np.loadtxt(os.path.join(args[1], 'Cooling_breakdown.txt'))
            a.semilogy(hb2[:, 0], hb2[:, 3]*0.1, 'r:', lw=1.2, label=lab_prev + ' heating')
            a.semilogy(cb2[:, 0], cb2[:, 3]*0.1, 'b:', lw=1.0, label=lab_prev + ' cooling')
        if C is not None and len(args) > 2 and os.path.exists(os.path.join(args[2], 'Heating_breakdown.txt')):
            hb3 = np.loadtxt(os.path.join(args[2], 'Heating_breakdown.txt')); cb3 = np.loadtxt(os.path.join(args[2], 'Cooling_breakdown.txt'))
            a.semilogy(hb3[:, 0], hb3[:, 3]*0.1, '--', color=col3, lw=1.3, label=lab_third + ' heating')
            a.semilogy(cb3[:, 0], cb3[:, 3]*0.1, '--', color='C9', lw=1.0, label=lab_third + ' cooling')
        a.semilogy(d9[:, 1], d9[:, 2], 'o', ms=4.5, mfc='none', color='r', label='Model A heating (Fig. 9)')
        a.semilogy(d9[:, 1], d9[:, 4], 's', ms=4, mfc='none', color='b', label='Model A radiative cooling (Fig. 9)')
        a.semilogy(d9[:, 1], d9[:, 3], '^', ms=4, mfc='none', color='0.4', label='Model A adiabatic (Fig. 9)')
        a.set_xlim(1, rmax); a.set_ylim(1e-12, 1e-7); a.set_xlabel(r'$r/r_{\rm base}$'); a.set_ylabel(r'rate [W m$^{-3}$]'); a.grid(alpha=0.3); a.legend(loc='upper right', fontsize=6.5)
        fig.tight_layout(); fig.savefig(os.path.join(fdir, 'heating_cooling.pdf')); plt.close(fig)

    # 4. the table of the memo: nine radii, ours (and the earlier run) against Model A
    RB = [1.05, 1.10, 1.15, 1.20, 1.30, 1.50, 2.00, 3.00]
    KT = {'f': [0.982, 0.976, 0.963, 0.936, 0.843, 0.688, 0.508, 0.403],
          'h3': [3.1e3, 5.0e3, 7.0e3, 9.1e3, 1.2e4, 1.1e4, 7.2e3, 1.9e3],
          'T': [1180, 1570, 1720, 1900, 2270, 2870, 4120, 4740],
          'ne': [5.0e6, 4.4e6, 7.2e6, 9.6e6, 1.0e7, 5.7e6, 1.8e6, 7.8e5]}
    def at(D, x):
        j = int(np.argmin(abs(D['r'] - x))); return j
    rows = []
    for k, x in enumerate(RB):
        j = at(A, x); cells = [f'{x:.2f}']
        for key, fmt in (('fH2', '{:.3f}'), ('H3p', '{:.1e}'), ('T', '{:.0f}'), ('e', '{:.1e}')):
            v = A[key][j]; cells.append(fmt.format(v))
            if B is not None: cells.append(fmt.format(B[key][at(B, x)]))
            if C is not None: cells.append(fmt.format(C[key][at(C, x)]))
            kk = {'fH2': 'f', 'H3p': 'h3', 'T': 'T', 'e': 'ne'}[key]
            cells.append(fmt.format(KT[kk][k]))
        rows.append(' & '.join(cells).replace('e+0', r'\times10^{').replace('e-0', r'\times10^{-') + r' \\')
    # exponents: 'a.bx10^{c' -> close the brace; only the 1-digit exponents occur here
    fixed = []
    for row in rows:
        out = ''; i = 0
        while i < len(row):
            if row.startswith(r'\times10^{', i):
                j = i + len(r'\times10^{'); sign = ''
                if row[j] == '-': sign = '-'; j += 1
                out += r'\times10^{' + sign + row[j] + '}'; i = j + 1
            else: out += row[i]; i += 1
        fixed.append('$' + out.replace(' & ', '$ & $').replace(r' \\', r'$ \\'))
    ncol = 2 + (B is not None) + (C is not None)
    head = (r'\begin{tabular}{l' + ('c'*ncol)*4 + '}' + '\n' + r'\toprule' + '\n'
            + r'$r/r_{\rm base}$ & \multicolumn{%d}{c}{$f({\rm H}_2)$} & \multicolumn{%d}{c}{$n({\rm H}_3^+)$ [cm$^{-3}$]} & \multicolumn{%d}{c}{$T$ [K]} & \multicolumn{%d}{c}{$n_e$ [cm$^{-3}$]} \\' % ((ncol,)*4) + '\n'
            + ' & ' + ' & '.join([(r'$h\nu - I$ & $h\nu$ & $0.55\,h\nu$ & A' if C is not None else (r'$h\nu - I$ & $h\nu$ & A' if B is not None else r'ours & A'))]*4) + r' \\' + '\n' + r'\midrule' + '\n')
    with open(os.path.join(fdir, 'comparison_table.tex'), 'w') as fh:
        fh.write(head + '\n'.join(fixed) + '\n' + r'\bottomrule' + '\n' + r'\end{tabular}' + '\n')
    print('figures and comparison_table.tex written to', fdir)

if __name__ == '__main__':
    main()
