#!/usr/bin/env python3
"""The one-key-at-a-time control ladder of docs/lhs1140b_exhale_vs_pwinds.tex
Sect. sec:basemetals, rung by rung, on the current binary.

For each rung: the escape rate from pp.log, the He 10830 red-pair depth from
he_line_metrics.py at both kernels (R = 68,000, the WINERED HIRES-Y kernel the
memo measures through, and R = 80,000, the kernel the stored ladder of
Update_EXHALE Sect. 78 was measured at), and the red-pair equivalent width on
the memo's window with ../crossings_gm25/measure.py's own definition.

usage: ./ladder_table.py <label>=<dir> ...
"""
import os, sys
import numpy as np

AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20


def red_ew(path):
    if not os.path.isfile(path):
        return float('nan')
    s = np.loadtxt(path)
    lam = s[:, 0]*AIR
    dep = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(dep[m], lam[m])


def metric(path, key='red_depth'):
    if not os.path.isfile(path):
        return float('nan')
    for line in open(path):
        if line.startswith(key):
            return float(line.split()[1])
    return float('nan')


def mdot(d):
    p = os.path.join(d, 'pp.log')
    v = float('nan')
    if os.path.isfile(p):
        for line in open(p):
            if 'steady-state Mdot' in line:
                v = float(line.split('=')[-1].split()[0])
    return v


rows = []
print('%-14s %10s %9s %9s %9s %9s' % ('rung', 'log10Mdot', 'red68k',
                                      'red80k', 'EW68k', 'EW80k'))
for a in sys.argv[1:]:
    lab, d = a.split('=', 1)
    r = dict(lab=lab,
             md=mdot(d),
             r68=metric(os.path.join(d, 'tpm_He10830_metrics.txt')),
             r80=metric(os.path.join(d, 'tpm_He10830_metrics_R80k.txt')),
             e68=red_ew(os.path.join(d, 'tpm_He10830.txt')),
             e80=red_ew(os.path.join(d, 'tpm_He10830_R80k.txt')))
    rows.append(r)
    print('%-14s %10.4f %9.4f %9.4f %9.5f %9.5f'
          % (lab, r['md'], r['r68'], r['r80'], r['e68'], r['e80']))

print()
print('step to step (each rung against the one above it):')
for a, b in zip(rows[:-1], rows[1:]):
    print('  %-12s -> %-12s  dlog10Mdot %+6.3f   red68k %+7.2f%%   '
          'red80k %+7.2f%%   EW68k %+7.2f%%'
          % (a['lab'], b['lab'], b['md']-a['md'],
             (b['r68']/a['r68']-1)*100, (b['r80']/a['r80']-1)*100,
             (b['e68']/a['e68']-1)*100))
