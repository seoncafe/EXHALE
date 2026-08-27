#!/usr/bin/env python3
"""One row per reservoir He/H of the flux-closed ladder (Phase F item 2).

Reads each arm's `closure_history.txt` for the converged iterate, the
converged iterate's `EXHALE_resolved.out` and profile header for the
composition and the cold trap, `tpm_He10830_metrics.txt` for the line, and
`tpm_He10830.txt` for the red-pair equivalent width on the same vacuum
window `make_memo_figures.py` uses.
"""
import os, re, sys
import numpy as np

AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20
EW_OBS, EW_ERR = 1.108, 0.030


def red_ew(path):
    if not os.path.isfile(path):
        return float('nan')
    s = np.loadtxt(path)
    lam = s[:, 0]*AIR
    dep = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(dep[m], lam[m])


def metrics(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        k, v = line.split()
        out[k] = float(v)
    return out


def last_row(case):
    p = os.path.join(case, 'closure_history.txt')
    rows = [l.split() for l in open(p) if l.strip() and not l.startswith('#')]
    return rows[-1]


def cold_trap_OH(kdir):
    p = os.path.join(kdir, 'lower_atmosphere_profile.dat')
    for line in open(p):
        if line.startswith('# notes'):
            m = re.search(r'O/H gas ([0-9.eE+-]+)', line)
            return float(m.group(1)) if m else float('nan')
    return float('nan')


HDR = ('reservoir  k  F_H[g/s]      F_He[g/s]     He/H@match   log10Mdot'
       '  red[%]  blue[%] FWHM[A]  EW[%A]  EW/obs  O/H_coldtrap')
print(HDR)
for case, resv in sys.argv[1:] and [(a.split('=')[0], a.split('=')[1])
                                    for a in sys.argv[1:]] or []:
    r = last_row(case)
    k = int(r[0])
    kdir = os.path.join(case, 'k%02d' % k)
    mt = metrics(os.path.join(kdir, 'tpm_He10830_metrics.txt'))
    ew = red_ew(os.path.join(kdir, 'tpm_He10830.txt'))
    print('%-9s %2d  %.5E  %.5E  %.6f  %8.3f  %6.3f  %6.3f  %6.4f  %6.4f'
          '  %6.3f  %.4e'
          % (resv, k, float(r[3]), float(r[4]), float(r[11]), float(r[12]),
             mt.get('red_depth', float('nan')),
             mt.get('blue_depth', float('nan')),
             mt.get('fwhm_A', float('nan')), ew, ew/EW_OBS,
             cold_trap_OH(kdir)))
