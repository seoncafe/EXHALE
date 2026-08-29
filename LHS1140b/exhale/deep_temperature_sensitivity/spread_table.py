#!/usr/bin/env python3
"""The spread of the handoff column over the runs of `runs/`.

Every run is the same physical problem: same planet, same reservoir, same
stellar flux, same mechanism, same K_zz, same imposed escape fluxes.  Only
`--climate-t-deep-guess` differs, and the guess is not physics -- it is where
MINPACK starts looking for the deep boundary temperature.  The spread down each
column is therefore the reproducibility of the stored profile and nothing else.

Three readings, because they answer different questions:
  index    -- row i against row i.  Photochem's grid is an altitude grid whose
              extent is part of the solution, so this mixes state with grid.
  pressure -- both columns interpolated onto common pressures.  This is what
              EXHALE reads.
  match    -- the single level EXHALE takes its base state from.

usage: spread_table.py [run-name prefix ...]   (default: the `g4` guess ladder)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'runs')


def read(d):
    cols, data = None, []
    for line in open(os.path.join(d, 'lower_atmosphere_profile.dat')):
        if line.startswith('#'):
            if line.startswith('# columns:'):
                cols = line.split(':', 1)[1].split()
            continue
        data.append([float(x) for x in line.split()])
    return cols, np.array(data)


def steady(d):
    t, top = float('nan'), float('nan')
    for line in open(os.path.join(d, 'adapter.log')):
        if 'STEADY' in line:
            t = float(line.split('t = ')[1].split(' s')[0])
        if line.startswith('wrote ./lower'):
            top = float(line.rsplit('..', 1)[1].split()[0])
    return t, top


KEEP = sys.argv[1:] or ['g4']
names = sorted(os.listdir(RUNS))
names = [n for n in names
         if any(n.startswith(k) for k in KEEP)
         and os.path.exists(os.path.join(RUNS, n,
                                         'lower_atmosphere_profile.dat'))]
tab = {n: read(os.path.join(RUNS, n)) for n in names}

print('# runs')
print('# %-12s %-12s %-14s %-12s' % ('run', 't_steady[s]', 'p_top[bar]',
                                     'exit'))
for n in names:
    t, top = steady(os.path.join(RUNS, n))
    print('  %-12s %-12.3e %-14.6e %-12s'
          % (n, t, top, 'equilibrium_time' if t > 1.0e17 else 'conv_longdy'))

ref = names[0]
cr, ar = tab[ref]

print('\n# spread over the %d runs, at fixed LEVEL INDEX '
      '(max over levels of max-min / mean)' % len(names))
for c in cr:
    stack = np.array([tab[n][1][:, tab[n][0].index(c)] for n in names])
    m = np.abs(stack).mean(axis=0)
    good = m > 0
    if not good.any():
        continue
    rel = (stack.max(axis=0) - stack.min(axis=0))[good]/m[good]
    print('  %-8s %.3e' % (c, rel.max()))

lo = max(tab[n][1][:, 0].min() for n in names)
hi = min(tab[n][1][:, 0].max() for n in names)
pq = np.logspace(np.log10(lo), np.log10(hi), 400)


def interp(cols, arr, name, p):
    x = np.log10(arr[:, 0])[::-1]
    y = arr[:, cols.index(name)][::-1]
    if (y > 0).all():
        return 10.0**np.interp(np.log10(p), x, np.log10(y))
    return np.interp(np.log10(p), x, y)


print('\n# spread over the %d runs, at fixed PRESSURE' % len(names))
for c in cr:
    if c == 'p':
        continue
    stack = np.array([interp(tab[n][0], tab[n][1], c, pq) for n in names])
    m = np.abs(stack).mean(axis=0)
    good = m > 0
    if not good.any():
        continue
    rel = (stack.max(axis=0) - stack.min(axis=0))[good]/m[good]
    print('  %-8s %.3e  at p = %.3e bar'
          % (c, rel.max(), pq[good][int(np.argmax(rel))]))

print('\n# spread over the %d runs, at the MATCHING LEVEL (1e-6 bar)'
      % len(names))
for c in cr:
    vals = []
    for n in names:
        cc, aa = tab[n]
        i = int(np.argmin(np.abs(aa[:, 0] - 1.0e-6)))
        vals.append(aa[i, cc.index(c)])
    v = np.array(vals)
    if v.mean() == 0.0:
        continue
    print('  %-8s %-22.12e %.3e' % (c, v.mean(),
                                    (v.max() - v.min())/abs(v.mean())))
