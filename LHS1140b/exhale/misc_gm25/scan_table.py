#!/usr/bin/env python3
"""One row per arm of the GJ 699 composition scan, in the columns the memo's
table carries: log10 Mdot, T_max and where it sits, v(20 R_p), the peak
metastable density and where it sits, and the red-pair depth, FWHM and
equivalent width of the He I 10830 line.

T, v and the metastable are read the way the memo reads them: the
hydrodynamic quantities and the metastable from the equilibrium
`Hydro_ioniz.txt` and `Ion_species.txt`.  The line comes from the arm's own
`tpm_He10830.txt` and `tpm_He10830_metrics.txt`.

usage: ./scan_table.py <label>=<dir> ...
"""
import os, sys
import numpy as np

AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20


def cols(p):
    return open(p).readlines()[1].split()[2:]


def red_ew(p):
    if not os.path.isfile(p):
        return float('nan')
    s = np.loadtxt(p)
    lam = s[:, 0]*AIR
    dep = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(dep[m], lam[m])


def metrics(p):
    out = {}
    if os.path.isfile(p):
        for line in open(p):
            if line.startswith('#') or not line.strip():
                continue
            k, v = line.split()
            out[k] = float(v)
    return out


def mdot(d):
    p = os.path.join(d, 'pp.log')
    if os.path.isfile(p):
        for line in open(p):
            if 'steady-state Mdot' in line:
                return float(line.split('=')[-1].split()[0])
    return float('nan')


if __name__ == '__main__':
    print('%-16s %8s %9s %8s %8s %11s %8s %8s %8s %8s'
          % ('label', 'He/H', 'log10Mdot', 'T_max', '(r/Rp)', 'v(20Rp)',
             'n23S_max', '(r/Rp)', 'red[%]', 'FWHM[A]'), end='')
    print(' %8s' % 'EW[%A]')
    for arg in sys.argv[1:]:
        lab, d = arg.split('=', 1)
        o = os.path.join(d, 'output')
        hn = cols(o+'/Hydro_ioniz.txt')
        h = np.loadtxt(o+'/Hydro_ioniz.txt')
        r, T, v = h[:, 0], h[:, hn.index('T[K]')], h[:, hn.index('v[cm/s]')]
        sn = cols(o+'/Ion_species.txt')
        a = np.loadtxt(o+'/Ion_species.txt')
        n23 = a[:, sn.index('HeITR')]
        heh = float('nan')
        for line in open(os.path.join(d, 'input.inp')):
            if line.startswith('He/H number ratio:'):
                heh = float(line.split(':')[1])
        i, j = int(np.argmax(T)), int(np.argmax(n23))
        mt = metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
        print('%-16s %8.4g %9.4f %8.0f %8.2f %11.2f %8.4g %8.2f %8.3f %8.4f %8.4f'
              % (lab, heh, mdot(d), T[i], r[i],
                 float(np.interp(20.0, r, v))/1e5, n23[j], a[j, 0],
                 mt.get('red_depth', float('nan')),
                 mt.get('fwhm_A', float('nan')),
                 red_ew(os.path.join(d, 'tpm_He10830.txt'))))
