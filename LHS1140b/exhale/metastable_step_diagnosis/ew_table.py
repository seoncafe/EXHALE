#!/usr/bin/env python3
"""He 2^3S column and He 10830 red-pair EW of the as-run and counterfactual runs.

EW definition copied from ../kzz_scan_table.py: the instrument-convolved
curve (column index 2) of tpm_He10830.txt normalized to its own continuum,
in percent, integrated over the vacuum window 10832.60-10834.20 A.
Column: line-of-sight integral of n(2^3S) along the sub-planetary ray,
    N = int n_23S dr  over the whole domain (cm^-2), from Ion_species_adv.txt.
"""
import numpy as np, os
SC = ('/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/'
      'a3b4f6bc-eb2d-4715-9dd2-73142c6adf66/scratchpad')
RJ = 7.1492e9
AIR_TO_VAC = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20
def ew(p):
    a = np.loadtxt(p); lam = a[:,0]*AIR_TO_VAC
    dep = (a[:,2].max() - a[:,2])/a[:,2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(dep[m], lam[m])
def depth(p):
    for line in open(p):
        if line.startswith('red_depth'): return float(line.split()[1])
    return float('nan')
def column(d):
    D = os.path.join(d, 'output')
    h = open(os.path.join(D,'Ion_species_adv.txt')).readlines()[1].split()[2:]
    a = np.loadtxt(os.path.join(D,'Ion_species_adv.txt'))
    Rp = 0.0
    for line in open(os.path.join(d,'input.inp')):
        if line.startswith('Planet radius'): Rp = float(line.split(':')[1])*RJ
    return np.trapz(a[:,h.index('HeITR')], a[:,0]*Rp)
runs = ['heh2p13_diff_kzz1e9','heh0p55','flux_closure/heh10p3/k01',
        'flux_closure/heh11p1/k01','heh0p55_diff_kzz1e10','heh0p55_diff_ctrl']
vars_ = ['base','lo','blendA','blendB']
print('%-24s %-7s %11s %8s %11s %8s %11s'
      %('run','variant','N23S[cm^-2]','dN/N','EW [%A]','dEW/EW','red depth[%]'))
for r in runs:
    N0 = e0 = None
    for v in vars_:
        d = os.path.join(SC,'cf_'+v,'exhale',r)
        N = column(d); E = ew(os.path.join(d,'tpm_He10830.txt'))
        Dp = depth(os.path.join(d,'tpm_He10830_metrics.txt'))
        if v == 'base': N0, e0 = N, E
        print('%-24s %-7s %11.4e %+7.2f%% %11.5f %+7.2f%% %11.5f'
              %(r.split('/')[-1] if v=='base' else '', v, N,
                100*(N/N0-1), E, 100*(E/e0-1), Dp))
    print()
