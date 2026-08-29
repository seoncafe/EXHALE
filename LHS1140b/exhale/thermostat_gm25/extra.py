#!/usr/bin/env python3
"""The remaining quantities the thermostat paragraphs quote: the two
definitions of heating per particle, the ionization-budget numbers at 2 Rp,
and the He III profile the compprof caption bounds."""
import os, sys
import numpy as np

def cols(p): return open(p).readlines()[1].split()[2:]

def at(r, y, r0): return float(np.interp(r0, r, y))

def go(label, d):
    o = os.path.join(d, 'output')
    hn = cols(os.path.join(o, 'Hydro_ioniz.txt'))
    sn = cols(os.path.join(o, 'Ion_species.txt'))
    h = np.loadtxt(os.path.join(o, 'Hydro_ioniz_adv.txt'))
    a = np.loadtxt(os.path.join(o, 'Ion_species_adv.txt'))
    s = {c: a[:, i] for i, c in enumerate(sn)}
    r = a[:, 0]
    T = h[:, hn.index('T[K]')]
    heat_adv = h[:, hn.index('heat[erg/cm3/s]')]
    nH = s['HI'] + s['HII']; nHe = s['HeI'] + s['HeII'] + s['HeIII']
    met = [c for c in sn if c not in ('r[Rp]','HI','HII','HeI','HeII','HeIII','HeITR')]
    nmet = sum(s[c] for c in met) if met else np.zeros_like(r)
    ntot = nH + nHe + nmet
    ne = s['HII'] + s['HeII'] + 2*s['HeIII']
    for c in met:
        if c.endswith('III'): ne = ne + 2*s[c]
        elif c.endswith('II'): ne = ne + s[c]
    hb = np.loadtxt(os.path.join(o, 'Heating_breakdown.txt'))
    rh, heat_eq, heat_HI = hb[:, 0], hb[:, 3], hb[:, 4]
    heat_HeI, heat_Herec = hb[:, 5], hb[:, 12]
    Rp_cm = None
    # planet radius from input.inp (R_J)
    for line in open(os.path.join(d, 'input.inp')):
        if line.startswith('Planet radius'):
            Rp_cm = float(line.split(':')[1])*7.1492e9
    # tau(13.6 eV) of the H I column above 2 Rp
    sig = 6.30e-18
    m = r >= 2.0
    NHI = np.trapz(s['HI'][m], r[m]*Rp_cm)
    # case-B recombination, Penning and He+ + H charge exchange at 2 Rp
    T2 = at(r, T, 2.0)
    # rec_HII_B = alphaB_HII_new: Badnell RR minus the Mao & Kaastra n=1 term
    def rr_badnell(T, A, B, T0, T1, Cc, T2_):
        bp = B + Cc*np.exp(-T2_/T); tt = np.sqrt(T/T0)
        return A/(tt*(1.0+tt)**(1.0-bp)*(1.0+np.sqrt(T/T1))**(1.0+bp))
    def rr_mao(T, a0, b0, c0, a1, b1, a2, b2):
        te = T*8.617333262e-5
        return a0*1e-10*te**(-b0-c0*np.log(te))*(1.0+a2*te**(-b2))/(1.0+a1*te**(-b1))
    aB = (rr_badnell(T2, 8.318e-11, 0.7472, 2.965, 7.001e5, 0.0, 1.0)
          - rr_mao(T2, 2.3390e-2, 1.2310, 1.0380e-2, 1.6430e1, 7.5080e-1,
                   8.9970e-2, 3.5560e-1))
    kcx = 1.25e-15*(T2/300.0)**0.25            # He+ + H0 -> He0 + H+ (row B2)
    lnT = np.log(T2)
    Q31 = 1.0e-9*np.exp(-8.64804e1/T2 - 2.86766e-1*lnT
                        + 8.68445e-2*lnT**2 - 5.73001e-3*lnT**3)
    f_pen = 0.9
    nHI2, nHII2, ne2 = at(r,s['HI'],2.0), at(r,s['HII'],2.0), at(r,ne,2.0)
    n232 = at(r, s['HeITR'], 2.0)
    rec = aB*nHII2*ne2
    pen = f_pen*n232*nHI2*Q31
    print('%-10s  T(2Rp)=%.1f' % (label, T2))
    print('   heat/particle @2Rp  eq-breakdown %.4e   adv-column %.4e erg/s'
          % (at(rh, heat_eq, 2.0)/at(r, ntot, 2.0),
             at(r, heat_adv, 2.0)/at(r, ntot, 2.0)))
    print('   heat_HI per neutral H @2Rp %.4e ; (heat_HI+heat_Herec)/nHI %.4e'
          % (at(rh, heat_HI, 2.0)/nHI2,
             (at(rh, heat_HI, 2.0)+at(rh, heat_Herec, 2.0))/nHI2))
    print('   heat_HeI/nHI @2Rp %.4e   heat_He_recomb/nHI @2Rp %.4e'
          % (at(rh, heat_HeI, 2.0)/nHI2, at(rh, heat_Herec, 2.0)/nHI2))
    print('   N(HI) above 2Rp %.4e cm^-2   tau(13.6 eV) %.4f' % (NHI, sig*NHI))
    print('   alpha_B(T) n_e @2Rp %.4e s^-1   case-B recomb rate %.4e cm^-3 s^-1'
          % (aB*ne2, rec))
    print('   Penning H+ source @2Rp %.4e  = %.4f of case-B removal'
          % (pen, pen/rec))
    print('   Penning per neutral H @2Rp %.4e s^-1' % (f_pen*n232*Q31))
    nHeII2 = at(r, s['HeII'], 2.0)
    cx = kcx*nHeII2*nHI2
    print('   He+ + H charge exchange @2Rp %.4e  = %.4f of case-B removal'
          % (cx, cx/rec))
    print('   remainder (photoionization) share %.4f'
          % (1.0 - pen/rec - cx/rec))
    x3 = s['HeIII']/nHe
    m3 = (r >= 1.0) & (r <= 3.0)
    print('   x(HeIII): @2Rp %.3e  @3Rp %.3e  max inside 3Rp %.3e at r=%.3f'
          % (at(r, x3, 2.0), at(r, x3, 3.0), x3[m3].max(), r[m3][x3[m3].argmax()]))
    print()

for a in sys.argv[1:]:
    lab, d = a.split('=', 1); go(lab, d)
