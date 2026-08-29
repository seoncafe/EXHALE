#!/usr/bin/env python3
"""Every number Table tab:thermostat and its two paragraphs quote, measured
on one run directory.

All structure quantities are read from the advection-corrected profiles
(`output/Hydro_ioniz_adv.txt`, `output/Ion_species_adv.txt`), as the memo
states; the heating and cooling channel breakdowns are the equilibrium files
(`Heating_breakdown.txt`, `Cooling_breakdown.txt`), which is where the code
writes them.

Definitions follow the memo and the code:
  mean mass per particle = (n_H + 4 n_He)/(n_H + n_He + n_e), nuclei counted
      as HI+HII and HeI+HeII+HeIII (`calc_ne`: metals in n_e when eos_metals
      is on; here the H/He electron count is used, as the memo's formula is
      written with n_e of the H/He mixture);
  n_tot follows `calc_ntot` -- one particle per species, electrons excluded,
      He 2^3S already inside He I;
  the 2^3S budget is `tr_triplet_row` of ion_residual_core.f90 with the rate
      coefficients of Cool_coeff.f90 (rec_HeII_23S, coex_HeI_1S_23S,
      coex_HeI_23S_21S, coex_HeI_23S_21P, ci_HeI23S, ioniz_HeI23S_H, and
      A31 = 1.272e-4).  The metastable photoionization rate g_heiTR is not
      written to any output file, so it is obtained as the residual of the
      balance; the reported closure error is the check that the rest is right.

usage: ./thermostat.py <label>=<dir> ...
"""
import os, sys
import numpy as np

AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20
KB_EV = 8.617333262e-5
A31 = 1.272e-4
E_TH_HETR = 4.8                       # eV
ERG2EV = 6.2415091e11                 # eV per erg (EXHALE's erg2eV)


def cols(path):
    return open(path).readlines()[1].split()[2:]


def load(d):
    o = os.path.join(d, 'output')
    hn = cols(os.path.join(o, 'Hydro_ioniz.txt'))
    sn = cols(os.path.join(o, 'Ion_species.txt'))
    h = np.loadtxt(os.path.join(o, 'Hydro_ioniz_adv.txt'))
    a = np.loadtxt(os.path.join(o, 'Ion_species_adv.txt'))
    s = {c: a[:, i] for i, c in enumerate(sn)}
    s['r'] = a[:, 0]
    s['T'] = h[:, hn.index('T[K]')]
    s['rho'] = h[:, hn.index('rho[mH/cm3]')]
    s['v'] = h[:, hn.index('v[cm/s]')]
    s['heat'] = h[:, hn.index('heat[erg/cm3/s]')]
    s['cool'] = h[:, hn.index('cool[erg/cm3/s]')]
    s['metal_names'] = [c for c in sn if c not in
                        ('r[Rp]', 'HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR')]
    return s


def rates(T):
    """The coefficients of Cool_coeff.f90, vectorized."""
    ups13 = 0.06703*np.exp(-2.673e-6*T) - 0.03584*np.exp(-3.880e-4*T)
    q13 = 2.10e-8*np.sqrt(13.60/(KB_EV*T))*np.exp(-19.81/(KB_EV*T))*ups13
    upsa = 2.847*np.exp(-1.252e-5*T) - 1.953*np.exp(-3.558e-4*T)
    q31a = 2.10e-8*np.sqrt(13.60/(KB_EV*T))*np.exp(-0.80/(KB_EV*T))*upsa/3.0
    upsb = 1.185*np.exp(-4.749e-6*T) - 0.9131*np.exp(-1.669e-4*T)
    q31b = 2.10e-8*np.sqrt(13.60/(KB_EV*T))*np.exp(-1.40/(KB_EV*T))*upsb/3.0
    a23s = 2.10e-13*(T/1.0e4)**(-0.778)
    lnT = np.log(T)
    Q31 = 1.0e-9*np.exp(-8.64804e1/T - 2.86766e-1*lnT
                        + 8.68445e-2*lnT**2 - 5.73001e-3*lnT**3)
    b23s = 6.41e-21*np.sqrt(T)*np.exp(-55338.0/T)/(E_TH_HETR/ERG2EV)
    return dict(q13=q13, q31a=q31a, q31b=q31b, a23s=a23s, Q31=Q31, b23s=b23s)


def at(r, y, r0):
    return float(np.interp(r0, r, y))


def column(r, n, lo, hi):
    m = (r >= lo) & (r <= hi)
    rr = np.concatenate(([lo], r[m], [hi]))
    nn = np.concatenate(([np.interp(lo, r, n)], n[m], [np.interp(hi, r, n)]))
    return float(np.trapz(nn, rr))


def red_ew(path):
    if not os.path.isfile(path):
        return float('nan')
    s = np.loadtxt(path)
    lam = s[:, 0]*AIR
    dep = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return float(np.trapz(dep[m], lam[m]))


def report(label, d):
    s = load(d)
    r = s['r']
    nH = s['HI'] + s['HII']
    nHe = s['HeI'] + s['HeII'] + s['HeIII']
    ne_hhe = s['HII'] + s['HeII'] + 2.0*s['HeIII']
    nmet = sum(s[c] for c in s['metal_names']) if s['metal_names'] else np.zeros_like(r)
    ne_met = np.zeros_like(r)
    for c in s['metal_names']:
        if c.endswith('III'):
            ne_met += 2.0*s[c]
        elif c.endswith('II'):
            ne_met += s[c]
    ne = ne_hhe + ne_met
    ntot = nH + nHe + nmet
    mu = (nH + 4.0*nHe)/(nH + nHe + ne_hhe)
    n23 = s['HeITR']

    print('=' * 68)
    print('%s   (%s)' % (label, d))
    print('=' * 68)
    print('  T_max                       %.1f K   at r = %.4f Rp'
          % (s['T'].max(), r[s['T'].argmax()]))
    print('  T(2 Rp)                     %.1f K' % at(r, s['T'], 2.0))
    print('  mean mass per particle(2Rp) %.4f mH' % at(r, mu, 2.0))
    print('  n(HI) at 2 Rp               %.4e cm^-3' % at(r, s['HI'], 2.0))
    print('  He/H by nuclei at 2 Rp      %.4f' % at(r, nHe/nH, 2.0))
    print('  He/H by nuclei at 5 Rp      %.4f' % at(r, nHe/nH, 5.0))
    print('  He/H by nuclei at 10 Rp     %.4f' % at(r, nHe/nH, 10.0))
    imax = n23.argmax()
    print('  max n(2^3S)                 %.4f cm^-3  at r = %.4f Rp'
          % (n23[imax], r[imax]))
    print('  int n(2^3S) dr  1-10 Rp     %.4f' % column(r, n23, 1.0, 10.0))
    print('     of which     1- 2 Rp     %.4f' % column(r, n23, 1.0, 2.0))
    print('     of which     2-10 Rp     %.4f' % column(r, n23, 2.0, 10.0))
    print('  EW red pair                 %.5f %%A' % red_ew(
        os.path.join(d, 'tpm_He10830.txt')))

    print('  -- composition / ionization at 2 Rp (calc_ntot convention) --')
    print('     nH/ntot %.4f   nHe/ntot %.4f   nmet/ntot %.3e'
          % (at(r, nH/ntot, 2.0), at(r, nHe/ntot, 2.0), at(r, nmet/ntot, 2.0)))
    print('     x(HII)  %.4f   x(HeII)  %.4f   x(HeIII) %.3e'
          % (at(r, s['HII']/nH, 2.0), at(r, s['HeII']/nHe, 2.0),
             at(r, s['HeIII']/nHe, 2.0)))
    m30 = (r >= 1.0) & (r <= 30.0)
    print('     max nmet/ntot over 1-30 Rp  %.3e   (max x(HeIII) inside 3 Rp %.3e)'
          % ((nmet/ntot)[m30].max(),
             (s['HeIII']/nHe)[(r >= 1.0) & (r <= 3.0)].max()))

    print('  -- 2^3S budget at the metastable peak (r = %.4f Rp, T = %.1f K) --'
          % (r[imax], s['T'][imax]))
    k = rates(s['T'])
    nSI = s['HeI'] - n23
    src_rec = ne*s['HeII']*k['a23s']
    src_cex = ne*nSI*k['q13']
    src = src_rec + src_cex
    lo_pen = s['HI']*k['Q31']
    lo_dee = ne*(k['q31a'] + k['q31b'])
    lo_ci = ne*k['b23s']
    lo_A = A31
    lo_known = lo_pen + lo_dee + lo_ci + lo_A
    lo_tot = src/np.maximum(n23, 1e-300)      # steady state: src = n23 * sink
    lo_pi = lo_tot - lo_known                 # residual = photoionization
    i = imax
    print('     source  rec %.4e  coll-exc %.4e  (coll-exc share %.3e)'
          % (src_rec[i], src_cex[i], src_cex[i]/src[i]))
    print('     sink [s^-1] total %.5e  Penning %.5e  e-deexc %.5e'
          % (lo_tot[i], lo_pen[i], lo_dee[i]))
    print('                 A31 %.5e  coll-ion %.5e  photoion(resid) %.5e'
          % (lo_A, lo_ci[i], lo_pi[i]))
    print('     shares: Penning %.4f  e-deexc %.4f  A31+photoion+collion %.4f'
          % (lo_pen[i]/lo_tot[i], lo_dee[i]/lo_tot[i],
             (lo_A + lo_ci[i] + lo_pi[i])/lo_tot[i]))
    print('     lifetime 1/sink            %.4f s' % (1.0/lo_tot[i]))
    print('     source rate at peak        %.5e cm^-3 s^-1' % src[i])
    print('     n(HI) at peak              %.5e cm^-3' % s['HI'][i])
    print('     n(2^3S) at peak            %.5f    src/sink check %.6f'
          % (n23[i], src[i]/(n23[i]*lo_tot[i])))

    # heating / cooling
    o = os.path.join(d, 'output')
    hb = np.loadtxt(os.path.join(o, 'Heating_breakdown.txt'))
    hbn = open(os.path.join(o, 'Heating_breakdown.txt')).readlines()[2]
    cb = np.loadtxt(os.path.join(o, 'Cooling_breakdown.txt'))
    rh = hb[:, 0]
    heat_tot, heat_HI, heat_met = hb[:, 3], hb[:, 4], hb[:, 9]
    cool_tot = cb[:, 3]
    cmet = cb[:, 11:].sum(axis=1)
    print('  -- heating / cooling (equilibrium breakdown) --')
    for rr in (1.2, 1.5, 2.0):
        ht = at(rh, heat_tot, rr)
        print('     r=%.1f Rp: cool/heat %.4f  metal share of cool %.4f  '
              'metal cool %.4e' % (rr, at(rh, cool_tot, rr)/ht,
                                   at(rh, cmet, rr)/max(at(rh, cool_tot, rr), 1e-300),
                                   at(rh, cmet, rr)))
    m12 = (rh >= 1.2) & (rh <= 2.0)
    print('     cool/heat over 1.2-2 Rp: %.4f - %.4f'
          % ((cool_tot/heat_tot)[m12].min(), (cool_tot/heat_tot)[m12].max()))
    print('     heating per particle at 2 Rp  %.5e erg/s'
          % (at(rh, heat_tot, 2.0)/at(r, ntot, 2.0)))
    print('     heating_HI per neutral H at 2 Rp %.5e erg/s'
          % (at(rh, heat_HI, 2.0)/at(r, s['HI'], 2.0)))
    print('     heat_Hpe[excitedH] max over domain %.3e' % np.abs(hb[:, 10]).max())
    if s['metal_names']:
        names = cols(os.path.join(o, 'Ion_species.txt'))
        cn = [c for c in names if c.startswith(('CI', 'CII', 'CIII'))]
        cidx = [11 + [c for c in names if c not in ('r[Rp]', 'HI', 'HII', 'HeI',
                'HeII', 'HeIII', 'HeITR')].index(c) for c in cn]
        big = cb[:, 11:]
        mnames = [c for c in names if c not in ('r[Rp]', 'HI', 'HII', 'HeI',
                  'HeII', 'HeIII', 'HeITR')]
        j = np.argmin(np.abs(rh - 1.2))
        top = np.argsort(big[j])[::-1][:3]
        print('     largest metal coolants at 1.2 Rp: '
              + ', '.join('%s %.3e' % (mnames[t], big[j, t]) for t in top))
    print()
    return dict(label=label, T=s['T'], r=r, n23=n23, mu=mu, nH=nH, nHe=nHe,
                HI=s['HI'], ne=ne, src=src, sink=lo_tot, pen=lo_pen,
                dee=lo_dee, imax=imax,
                ew=red_ew(os.path.join(d, 'tpm_He10830.txt')))


if __name__ == '__main__':
    out = []
    for a in sys.argv[1:]:
        lab, d = a.split('=', 1)
        out.append(report(lab, d))
