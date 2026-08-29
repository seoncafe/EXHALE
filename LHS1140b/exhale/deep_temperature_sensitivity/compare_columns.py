#!/usr/bin/env python3
"""Two handoff columns, compared at fixed level index and at fixed pressure.

The distinction is the whole point.  Photochem's grid is an ALTITUDE grid
whose extent is itself part of the solution (`var.top_atmos` is solved for so
that the average top-of-atmosphere pressure meets `TOA_pressure_avg`), so the
pressure of level `i` is an output, not an axis.  Comparing two columns row by
row therefore mixes a difference of state with a difference of grid.  EXHALE
reads the file as a table over PRESSURE, so the second comparison -- both
columns interpolated onto the same pressures -- is the one that says what the
wind actually inherits.

usage: compare_columns.py <dirA> <dirB>
"""
import os
import sys

import numpy as np


def read(d):
    path = d if d.endswith('.dat') else os.path.join(
        d, 'lower_atmosphere_profile.dat')
    cols, data, hdr = None, [], {}
    for line in open(path):
        if line.startswith('#'):
            if line.startswith('# columns:'):
                cols = line.split(':', 1)[1].split()
            else:
                w = line[1:].split()
                if len(w) >= 2:
                    hdr[w[0]] = w[1]
            continue
        data.append([float(x) for x in line.split()])
    return cols, np.array(data), hdr


ca, a, ha = read(sys.argv[1])
cb, b, hb = read(sys.argv[2])
print('A = %s' % sys.argv[1])
print('B = %s' % sys.argv[2])
print('levels %d / %d   top %.6e / %.6e bar   deep %.6e / %.6e bar'
      % (len(a), len(b), a[:, 0].min(), b[:, 0].min(),
         a[:, 0].max(), b[:, 0].max()))

if len(a) == len(b):
    dp = np.abs(b[:, 0]/a[:, 0] - 1.0)
    i = int(np.argmax(dp))
    print('\nat fixed LEVEL INDEX')
    print('  pressure: worst %.3e relative at level %d '
          '(%.6e vs %.6e bar)' % (dp[i], i, a[i, 0], b[i, 0]))
    worst = []
    for j, c in enumerate(ca):
        if c == 'p':
            continue
        x, y = a[:, j], b[:, cb.index(c)]
        m = x != 0.0
        if not m.any():
            continue
        r = np.abs(y[m]/x[m] - 1.0)
        worst.append((r.max(), c, int(np.argmax(r))))
    for r, c, k in sorted(worst, reverse=True)[:8]:
        print('  %-8s worst %.3e relative' % (c, r))

print('\nat fixed PRESSURE (log-linear interpolation onto the common range)')
lo = max(a[:, 0].min(), b[:, 0].min())
hi = min(a[:, 0].max(), b[:, 0].max())
pq = np.logspace(np.log10(lo), np.log10(hi), 200)


def interp(cols, arr, name, p):
    j = cols.index(name)
    x = np.log10(arr[:, 0])[::-1]
    y = arr[:, j][::-1]
    pos = y > 0
    if pos.all():
        return 10.0**np.interp(np.log10(p), x, np.log10(y))
    return np.interp(np.log10(p), x, y)


rows = []
for c in ca:
    if c == 'p':
        continue
    ya = interp(ca, a, c, pq)
    yb = interp(cb, b, c, pq)
    m = ya != 0.0
    if not m.any():
        continue
    r = np.abs(yb[m]/ya[m] - 1.0)
    rows.append((r.max(), c, float(pq[m][int(np.argmax(r))])))
for r, c, p in sorted(rows, reverse=True):
    print('  %-8s worst %.3e relative at p = %.3e bar' % (c, r, p))

print('\nat the MATCHING LEVEL (the node at 1e-6 bar)')
ia = int(np.argmin(np.abs(a[:, 0] - 1.0e-6)))
ib = int(np.argmin(np.abs(b[:, 0] - 1.0e-6)))
print('  p_A = %.6e   p_B = %.6e bar' % (a[ia, 0], b[ib, 0]))
for j, c in enumerate(ca):
    x, y = a[ia, j], b[ib, cb.index(c)]
    rel = abs(y/x - 1.0) if x != 0.0 else float('nan')
    print('  %-8s %-22.12e %-22.12e  %.3e' % (c, x, y, rel))
