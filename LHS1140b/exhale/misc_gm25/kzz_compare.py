#!/usr/bin/env python3
"""The numbers Fig. lhs1140b_kzz_profiles.pdf and its caption quote, read off
four solutions that differ only in the eddy coefficient at a fixed reservoir.

Read exactly as ../../make_memo_figures.py reads them: the hydrodynamic
panels (T, v, rho) from the equilibrium `Hydro_ioniz.txt`, the composition and
metastable panels from `Ion_species_adv.txt`, and the elemental helium as
HeI + HeII + HeIII divided by the base reservoir He/H = 0.55.  HeITR is a
level of He I and is already inside the HeI column (species_table.f90).

usage: ./kzz_compare.py <label>=<dir>[:K_zz] ...
"""
import os, re, sys
import numpy as np

HEH0 = 0.55


def cols(p):
    return open(p).readlines()[1].split()[2:]


def read(d):
    o = os.path.join(d, 'output')
    hn = cols(o+'/Hydro_ioniz.txt')
    h = np.loadtxt(o+'/Hydro_ioniz.txt')
    sn = cols(o+'/Ion_species_adv.txt')
    a = np.loadtxt(o+'/Ion_species_adv.txt')
    nh = a[:, sn.index('HI')] + a[:, sn.index('HII')]
    nhe = sum(a[:, sn.index(s)] for s in ('HeI', 'HeII', 'HeIII'))
    return dict(r=h[:, 0], T=h[:, hn.index('T[K]')], v=h[:, hn.index('v[cm/s]')],
                rho=h[:, hn.index('rho[mH/cm3]')], r_adv=a[:, 0],
                ratio=nhe/np.maximum(nh, 1e-99)/HEH0,
                n23s=a[:, sn.index('HeITR')])


def homopause(d, kzz):
    if kzz <= 0.0:
        return float('nan')
    for rel in ('output/element_flux_profile.txt', 'diffusion_faceflux.txt'):
        p = os.path.join(d, rel)
        if os.path.isfile(p):
            break
    else:
        return float('nan')
    icol = 8
    for line in open(p):
        if not line.startswith('#'):
            break
        if 'columns:' in line:
            head = line.split('columns:', 1)[1].strip()
            toks = (re.split(r'\s{2,}', head) if re.search(r'\s{2,}', head)
                    else head.split())
            if 'D_eff[cm2/s]' in toks:
                icol = toks.index('D_eff[cm2/s]')
            break
    a = np.loadtxt(p, usecols=(1, icol))
    s = a[:, 1] - kzz
    for j in range(len(s)-1):
        if s[j]*s[j+1] < 0.0:
            f = s[j]/(s[j]-s[j+1])
            return a[j, 0] + f*(a[j+1, 0]-a[j, 0])
    return float('nan')


def mdot(d):
    p = os.path.join(d, 'pp.log')
    if os.path.isfile(p):
        for line in open(p):
            if 'steady-state Mdot' in line:
                return float(line.split('=')[-1].split()[0])
    return float('nan')


if __name__ == '__main__':
    print('%-10s %9s %10s %11s %9s %11s %10s %8s %10s'
          % ('label', 'log10Mdot', 'v(10Rp)', 'rho(2Rp)', 'T(1.5Rp)',
             'He/H(1.05)', 'peak n23S', 'r_peak', 'homopause'))
    rows = []
    for arg in sys.argv[1:]:
        lab, d = arg.split('=', 1)
        kzz = 0.0
        if ':' in d:
            d, k = d.rsplit(':', 1)
            kzz = float(k)
        p = read(d)
        j = int(np.argmax(p['n23s']))
        row = dict(lab=lab, md=mdot(d),
                   v10=float(np.interp(10.0, p['r'], p['v']))/1e5,
                   rho2=float(np.interp(2.0, p['r'], p['rho'])),
                   T15=float(np.interp(1.5, p['r'], p['T'])),
                   rat=float(np.interp(1.05, p['r_adv'], p['ratio'])),
                   n23s=p['n23s'][j], rpk=p['r_adv'][j], hp=homopause(d, kzz))
        rows.append(row)
        print('%-10s %9.4f %10.4g %11.4g %9.1f %11.3f %10.4g %8.3f %10.4f'
              % (lab, row['md'], row['v10'], row['rho2'], row['T15'],
                 row['rat'], row['n23s'], row['rpk'], row['hp']))
    if len(rows) > 1:
        f = lambda k: max(x[k] for x in rows)/min(x[k] for x in rows)
        print('\nspread: v(10Rp) %.1f%%, rho(2Rp) %.1f%%, log10Mdot %.3f dex; '
              'T(1.5Rp) %.0f -> %.0f K; peak n(2^3S) %.3g -> %.3g'
              % (100*(f('v10')-1), 100*(f('rho2')-1),
                 max(x['md'] for x in rows)-min(x['md'] for x in rows),
                 rows[0]['T15'], rows[-1]['T15'],
                 rows[0]['n23s'], rows[-1]['n23s']))
