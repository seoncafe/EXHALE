#!/usr/bin/env python3
"""Count genuine discontinuities in n(2^3S), independent of grid stretching.

A fixed cell-ratio window such as [0.9, 1.12] is not scale free: on the
outer stretched grid (dr ~ 0.4 R_p) a smooth exponential decline trips it.
Here the cell ratio is compared with the LOCAL TREND, the geometric mean of
the ratios of the two flanking cell pairs; a cell is flagged when its own
ratio departs from that trend by more than 10%.
"""
import numpy as np, os, sys
B='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/exhale/'
def flags(r, n, rmin=1.05, tol=0.10):
    q = np.full_like(n, np.nan)
    q[1:] = n[1:]/np.maximum(n[:-1], 1e-300)
    out = []
    for i in range(2, len(n)-2):
        if r[i] < rmin or n[i] <= 0 or n[i-1] <= 0: continue
        trend = np.sqrt(q[i-1]*q[i+1]) if (q[i-1] > 0 and q[i+1] > 0) else np.nan
        if not np.isfinite(trend) or trend <= 0: continue
        if abs(q[i]/trend - 1.0) > tol:
            out.append((i, r[i], q[i], q[i]/trend))
    return out
def load(run, which):
    D = os.path.join(B, run, 'output')
    h = open(os.path.join(D, 'Ion_species.txt')).readlines()[1].split()[2:]
    a = np.loadtxt(os.path.join(D, 'Ion_species%s.txt' % ('' if which=='eq' else '_adv')))
    hy = np.loadtxt(os.path.join(D, 'Hydro_ioniz%s.txt' % ('' if which=='eq' else '_adv')))
    hn = open(os.path.join(D, 'Hydro_ioniz.txt')).readlines()[1].split()[2:]
    return a[:,0], a[:,h.index('HeITR')], hy[:,hn.index('T[K]')]
runs = ['heh0p55','heh2p13_diff_kzz1e9','flux_closure/heh10p3/k01',
        'flux_closure/heh11p1/k01','flux_closure/hi/k06',
        'heh0p55_diff_kzz1e10','heh0p55_diff_kzz1e9','heh0p55_diff_kzz1e8',
        'heh0p55_diff_ctrl']
for run in runs:
    for w in ('eq','adv'):
        r,n,T = load(run,w)
        f = flags(r,n)
        tag = lambda i: ('Penning T=4000' if (T[i-1]-4000)*(T[i]-4000) <= 0 else
                         ('outer boundary' if r[i] > 29 else 'OTHER'))
        print('%-26s %-4s  %d flagged: %s' % (run.split('/')[-1], w, len(f),
              '; '.join('r=%.4f q=%.3f [%s]' % (x[1],x[2],tag(x[0])) for x in f[:8])))
