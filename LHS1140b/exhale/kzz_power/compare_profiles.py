#!/usr/bin/env python3
"""Level-by-level comparison of two lower-atmosphere profiles."""
import sys
import numpy as np


def read(path):
    names, rows = None, []
    hdr = {}
    for line in open(path):
        if line.startswith('#'):
            if line.startswith('# columns:'):
                names = line.split(':', 1)[1].split()
            elif line.startswith('# notes'):
                hdr['notes'] = line[7:].strip()
            continue
        rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    return {n: a[:, i] for i, n in enumerate(names)}, hdr


def at(c, p_bar, key):
    lp = np.log(c['p'])
    o = np.argsort(lp)
    return float(np.interp(np.log(p_bar), lp[o], np.asarray(c[key])[o]))


if __name__ == '__main__':
    ref, hr = read(sys.argv[1])
    new, hn = read(sys.argv[2])
    levels = [('deep     16.73 bar', 16.7309812445485271),
              ('troppause 2.636 bar', 2.636),
              ('1 bar             ', 1.0),
              ('1 mbar            ', 1.0e-3),
              ('match     1 ubar  ', 1.0e-6)]
    keys = ['T', 'Kzz', 'X_He', 'X_C', 'X_N', 'X_O', 'q_H2', 'q_H', 'q_H2O',
            'q_CO', 'q_CO2', 'q_CH4', 'q_NH3', 'q_HCN', 'q_N2', 'q_OH']
    print('%-20s %-8s %12s %12s %10s' % ('level', 'field', 'const 1e9',
                                         'K_zz(p)', 'new/ref'))
    for lab, p in levels:
        for k in keys:
            if k not in ref or k not in new:
                continue
            a, b = at(ref, p, k), at(new, p, k)
            r = b/a if a != 0 else float('nan')
            print('%-20s %-8s %12.5E %12.5E %10.4G' % (lab, k, a, b, r))
        print()
