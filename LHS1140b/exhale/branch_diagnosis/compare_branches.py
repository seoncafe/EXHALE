#!/usr/bin/env python3
"""T1: where in radius do the cool and hot converged states differ?

Overlays the two K_zz = 1e7, He/H = 3.4 solutions of kzz_profile_scan and
reports the relative difference |dq|/q of T, v, rho, x(HII), x(HeII) and
n(He 2^3S) as a function of radius, together with the share of the total
integrated difference that falls inside the He 10830 line-forming region
against the share that falls beyond the exobase.

Usage:  python3 compare_branches.py <cool_output_dir> <hot_output_dir> <label>
"""
import sys
import os
import numpy as np

HYD_COLS = ['r', 'rho', 'v', 'p', 'T']
# Ion_species column order after r
ION_NAMES = ['HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR',
             'CI', 'CII', 'CIII', 'OI', 'OII', 'OIII', 'NI', 'NII', 'NIII',
             'MgI', 'MgII', 'MgIII', 'SiI', 'SiII', 'SiIII',
             'CaI', 'CaII', 'CaIII', 'NaI', 'NaII', 'KI', 'KII',
             'SI', 'SII', 'FeI', 'FeII', 'FeIII']


def load(outdir, adv=True):
    suf = '_adv' if adv else ''
    h = np.loadtxt(os.path.join(outdir, 'Hydro_ioniz%s.txt' % suf))
    i = np.loadtxt(os.path.join(outdir, 'Ion_species%s.txt' % suf))
    d = {'r': h[:, 0], 'rho': h[:, 1], 'v': h[:, 2], 'p': h[:, 3], 'T': h[:, 4],
         'heat': h[:, 5], 'cool': h[:, 6]}
    for k, nm in enumerate(ION_NAMES):
        d['n_' + nm] = i[:, 1 + k]
    nH = d['n_HI'] + d['n_HII']
    nHe = d['n_HeI'] + d['n_HeII'] + d['n_HeIII']
    d['x_HII'] = d['n_HII'] / np.maximum(nH, 1e-300)
    d['x_HeII'] = d['n_HeII'] / np.maximum(nHe, 1e-300)
    d['x_HeIII'] = d['n_HeIII'] / np.maximum(nHe, 1e-300)
    d['n_He23S'] = d['n_HeITR']
    d['nH'] = nH
    d['nHe'] = nHe
    return d


def main():
    cool_dir, hot_dir, label = sys.argv[1], sys.argv[2], sys.argv[3]
    c = load(cool_dir)
    h = load(hot_dir)

    # the two runs share the grid by construction; verify
    assert len(c['r']) == len(h['r']), 'grid length differs'
    dr = np.max(np.abs(c['r'] - h['r']) / np.maximum(c['r'], 1e-30))
    r = c['r']

    quantities = ['T', 'v', 'rho', 'x_HII', 'x_HeII', 'n_He23S', 'p',
                  'n_HI', 'n_HeI', 'heat', 'cool']

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'T1_profile_diff_%s.txt' % label)
    with open(out, 'w') as f:
        f.write('# T1: cool vs hot converged state, relative difference vs radius\n')
        f.write('# cool = %s\n# hot  = %s\n' % (cool_dir, hot_dir))
        f.write('# max relative grid mismatch = %.3e\n' % dr)
        f.write('# columns: r[Rp] ' +
                ' '.join('%s_cool %s_hot reldiff_%s' % (q, q, q)
                         for q in quantities) + '\n')
        for j in range(len(r)):
            row = ['%18.10e' % r[j]]
            for q in quantities:
                a, b = c[q][j], h[q][j]
                den = max(abs(a), abs(b), 1e-300)
                row.append('%14.6e %14.6e %12.4e' % (a, b, abs(b - a) / den))
            f.write(' '.join(row) + '\n')

    # --- where does the difference live? ---
    print('=== %s ===' % label)
    print('grid mismatch (max relative) = %.3e' % dr)
    print('r range = %.4f .. %.4f R_p, N = %d' % (r[0], r[-1], len(r)))
    print()

    bands = [('base 1-1.2',   1.0, 1.2),
             ('line 1-3',     1.0, 3.0),
             ('line 1-5',     1.0, 5.0),
             ('mid 3-10',     3.0, 10.0),
             ('outer 10-20',  10.0, 20.0),
             ('exo >20',      20.0, 1e9)]

    print('%-14s %8s' % ('band', 'Ncells') +
          ''.join('%12s' % q for q in quantities))
    print('%-14s %8s' % ('', '') +
          ''.join('%12s' % 'max|dq|/q' for q in quantities))
    for nm, lo, hi in bands:
        m = (r >= lo) & (r < hi)
        if m.sum() == 0:
            continue
        line = '%-14s %8d' % (nm, m.sum())
        for q in quantities:
            a, b = c[q][m], h[q][m]
            den = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)
            line += '%12.3e' % np.max(np.abs(b - a) / den)
        print(line)
    print()

    # share of the L1-norm of the difference, weighted by cell radial extent,
    # that falls in each band (normalized per quantity)
    w = np.gradient(r)
    print('share of the radius-weighted integrated |dq|/q in each band [%]')
    print('%-14s' % 'band' + ''.join('%12s' % q for q in quantities))
    tot = {}
    for q in quantities:
        a, b = c[q], h[q]
        den = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)
        tot[q] = np.sum(np.abs(b - a) / den * w)
    for nm, lo, hi in bands:
        m = (r >= lo) & (r < hi)
        if m.sum() == 0:
            continue
        line = '%-14s' % nm
        for q in quantities:
            a, b = c[q][m], h[q][m]
            den = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)
            s = np.sum(np.abs(b - a) / den * w[m])
            line += '%12.2f' % (100.0 * s / max(tot[q], 1e-300))
        print(line)
    print()

    # point values at a few radii
    print('point values')
    for rq in [1.02, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 25.0]:
        if rq > r[-1]:
            continue
        j = int(np.argmin(np.abs(r - rq)))
        print('  r = %6.3f : T %8.1f / %8.1f K   v %10.3e / %10.3e   '
              'rho %10.3e / %10.3e   x_HII %7.4f / %7.4f   '
              'x_HeII %8.5f / %8.5f   n23S %10.3e / %10.3e'
              % (r[j], c['T'][j], h['T'][j], c['v'][j], h['v'][j],
                 c['rho'][j], h['rho'][j], c['x_HII'][j], h['x_HII'][j],
                 c['x_HeII'][j], h['x_HeII'][j],
                 c['n_He23S'][j], h['n_He23S'][j]))
    print()
    print('T_max: cool %.1f K at r = %.3f ; hot %.1f K at r = %.3f'
          % (c['T'].max(), r[np.argmax(c['T'])],
             h['T'].max(), r[np.argmax(h['T'])]))
    print('wrote %s' % out)


if __name__ == '__main__':
    main()
