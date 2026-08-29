#!/usr/bin/env python3
"""The measurement record of the K_zz(p) test (see closure.json comments)."""
import math, os, sys
import numpy as np
sys.path.insert(0, '../flux_closure')
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    import closure_heh_table as T   # printing its own header at import
sys.path.insert(0, '.')
from compare_profiles import read, at

ARMS = [('11.1', '../flux_closure/heh11p1/k01', 'heh11p1_ref/k00',
         'heh11p1_chem/k00', 'heh11p1_full/k00'),
        ('2.09', '../flux_closure/hi/k06', 'hi_ref/k00',
         'hi_chem/k00', 'hi_full/k00')]
LAB = ('recorded (other seed)', 'reference K_zz = 1e9',
       'A chemistry only', 'B chemistry + wind')


def wind(d):
    a = np.loadtxt(d + '/output/element_flux_profile.txt', usecols=(1, 2))
    r, X = a[:, 0], a[:, 1]
    heh = (X/4.0)/(1.0 - X)
    h = np.loadtxt(d + '/output/Hydro_ioniz.txt')
    lg = float('nan')
    for l in open(os.path.join(d, 'pp.log'), errors='replace'):
        if 'steady-state Mdot' in l:
            lg = float(l.split('=')[1].split()[0].replace('D', 'E'))
    return (float(np.interp(2.0, r, heh)), float(np.interp(10.0, r, heh)),
            h[:, 4].max(), lg)


print('K_zz(p) = 1.0e9 (p / 1.1e-8 bar)^-1/2, against the constant 1.0e9')
print('anchored at the profile top so the wind inherits the same value above')
print('the file.  Implied: 1.05e8 at 1 ubar (the match), 1.05e5 at 1 bar,')
print('2.56e4 at 16.7 bar.\n')
for arm, rec, ref, cha, ful in ARMS:
    print('=== reservoir He/H = %s ===' % arm)
    print('%-22s %8s %8s %9s %8s %9s %9s %8s %8s'
          % ('run', 'red[%]', 'FWHM[A]', 'EW[%A]', 'EW/obs', 'He/H@2Rp',
             'He/H@10Rp', 'Tmax[K]', 'lgMdot'))
    ews = {}
    for lab, d in zip(LAB, (rec, ref, cha, ful)):
        mt = T.metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
        ew = T.red_ew(os.path.join(d, 'tpm_He10830.txt'))
        h2, h10, tmax, lg = wind(d)
        ews[lab] = ew
        print('%-22s %8.4f %8.4f %9.4f %8.4f %9.4f %9.4f %8.1f %8.3f'
              % (lab, mt.get('red_depth', np.nan), mt.get('fwhm_A', np.nan),
                 ew, ew/T.EW_OBS, h2, h10, tmax, lg))
    e0 = ews[LAB[1]]
    print('  same profile, other seed: EW %+.2f %% -- the JFNK path spread'
          % (100*(ews[LAB[0]]/e0 - 1.0)))
    for lab in LAB[2:]:
        d = math.log(ews[lab]/e0)
        print('  %-20s EW %+.2f %%, crossing 11.73 -> %.2f (%+.2f %%)'
              % (lab, 100*(ews[lab]/e0 - 1.0), 11.73*math.exp(-d/0.3819),
                 100*(math.exp(-d/0.3819) - 1.0)))
    prof_ref, prof_new = read(rec + '/lower_atmosphere_profile.dat')[0], \
        read(ful + '/lower_atmosphere_profile.dat')[0]
    print('  at the match (1 ubar), K_zz(p)/const:', end=' ')
    for k in ('X_He', 'X_C', 'X_N', 'X_O'):
        print('%s %.4f' % (k, at(prof_new, 1e-6, k)/at(prof_ref, 1e-6, k)),
              end='  ')
    print()
    print('  carriers at the match, K_zz(p)/const:', end=' ')
    for k in ('q_CH4', 'q_CO', 'q_CO2', 'q_NH3', 'q_N2', 'q_HCN', 'q_H2O',
              'q_H'):
        print('%s %.3G' % (k, at(prof_new, 1e-6, k)/at(prof_ref, 1e-6, k)),
              end='  ')
    print('\n')
