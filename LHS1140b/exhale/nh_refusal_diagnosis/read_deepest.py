#!/usr/bin/env python3
"""Deepest-level El/H departure from the input abundances, read out of any
lower_atmosphere_profile.dat.  This is exactly the quantity the adapter's
conservation check (lower_profile_schema.write_handoff) forms: the file's
deepest row IS ratios[el][0] when no --p-deep clip is asked for."""
import sys
import numpy as np

SOLAR = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5,
             O=6.06178e-4, S=1.31826e-5)


def read(path):
    cols, hdr = None, {}
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                if line.startswith('# columns:'):
                    cols = line.split(':', 1)[1].split()
                else:
                    w = line[1:].split(None, 1)
                    if w:
                        hdr[w[0]] = w[1].strip() if len(w) > 1 else ''
                continue
            rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    return cols, hdr, a


def deviations(path, heh):
    cols, hdr, a = read(path)
    want = dict(SOLAR)
    want['He'] = heh
    out = {}
    for el in ('He', 'C', 'N', 'O'):
        key = 'X_' + el
        if key not in cols:
            continue
        got = a[0, cols.index(key)]
        out[el] = (got, abs(got/(want[el]/want['H']) - 1.0))
    return cols, hdr, a, out


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        path, _, heh = spec.partition('@')
        cols, hdr, a, out = deviations(path, float(heh))
        print('%-58s p_deep=%.4e T=%.1f' % (path, a[0, 0], a[0, 2]))
        for el, (got, dev) in out.items():
            print('    %-3s got %.9e  dev %.3e' % (el, got, dev))
