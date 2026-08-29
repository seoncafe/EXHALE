#!/usr/bin/env python3
"""Rewrite the He 2^3S column of Ion_species_adv.txt with the He(2^3S)+H
   Penning rate coefficient made CONTINUOUS in T, so that the T = 4000 K
   discontinuity of the Taylor et al. (2025) two-branch fit is removed.

   The as-run rate (Cool_coeff.f90 :: penning_HeI_23S) is
       Q31 = 1.9e-9*(300/T)**0.07   (T <= 4000 K)
       Q31 = 9.1e-9*(300/T)**0.50   (T >  4000 K)
   whose branches differ by 1.572 at the join.

   Counterfactuals:
     lo         low branch at every T (upper bound on n_23S above 4000 K)
     blendA     log-log blend of the two branches over 3000-5000 K
     blendB     log-log blend of the two branches over 2000-8000 K

   n_23S = S/L with S the (continuous) source and
     L = Q31*n_HI + A31 + (q31a+q31b)*n_e + P_HeITR + a_ion_HeITR*n_e,
   so n_23S(cf) = n_23S(as run) * L(as run)/L(cf).  P_HeITR and
   a_ion_HeITR are not written to any output file and are dropped from L;
   since they only add a continuous positive term, dropping them makes
   |L0/Lcf - 1| larger, i.e. the rescaling is an UPPER bound on |dn_23S|.

   usage: make_counterfactual.py <run_dir> <lo|blendA|blendB>
"""
import numpy as np, os, sys
kb_eV = 8.617333262e-5
A31   = 1.272e-4
def Q31_lo(T): return 1.9e-9*(3.0e2/T)**0.07
def Q31_hi(T): return 9.1e-9*(3.0e2/T)**0.50
def Q31(T):    return np.where(T <= 4.0e3, Q31_lo(T), Q31_hi(T))
def Q31_blend(T, T1, T2):
    w = np.clip((np.log(T) - np.log(T1))/(np.log(T2) - np.log(T1)), 0.0, 1.0)
    return np.exp((1.0 - w)*np.log(Q31_lo(T)) + w*np.log(Q31_hi(T)))
def q31a(T):
    u = 2.847*np.exp(-1.252e-5*T) - 1.953*np.exp(-3.558e-4*T)
    return 2.10e-8*np.sqrt(13.60/(kb_eV*T))*np.exp(-0.80/(kb_eV*T))*u/3.0
def q31b(T):
    u = 1.185*np.exp(-4.749e-6*T) - 0.9131*np.exp(-1.669e-4*T)
    return 2.10e-8*np.sqrt(13.60/(kb_eV*T))*np.exp(-1.40/(kb_eV*T))*u/3.0

def rescale(d, mode):
    D   = os.path.join(d, 'output')
    p   = os.path.join(D, 'Ion_species_adv.txt')
    hdr = open(p).readlines()[:2]
    names = hdr[1].split()[2:]
    itr, ihi = names.index('HeITR'), names.index('HI')
    a  = np.loadtxt(p)
    hn = open(os.path.join(D, 'Hydro_ioniz_adv.txt')).readlines()[1].split()[2:]
    T  = np.loadtxt(os.path.join(D, 'Hydro_ioniz_adv.txt'))[:, hn.index('T[K]')]
    ne = np.loadtxt(os.path.join(D, 'Heating_breakdown.txt'))[:, 2]
    nHI = a[:, ihi]
    Lq  = A31 + (q31a(T) + q31b(T))*ne
    L0  = Q31(T)*nHI + Lq
    Qcf = {'lo':     Q31_lo(T),
           'blendA': Q31_blend(T, 3000.0, 5000.0),
           'blendB': Q31_blend(T, 2000.0, 8000.0)}[mode]
    Lcf = Qcf*nHI + Lq
    n3  = a[:, itr].copy()
    a[:, itr] = n3*L0/Lcf
    with open(p, 'w') as f:
        f.writelines(hdr)
        for row in a:
            f.write(''.join('%26.17E' % v for v in row) + '\n')
    ch = a[:, itr]/np.maximum(n3, 1e-300)
    return ch.min(), ch.max()

if __name__ == '__main__':
    lo, hi = rescale(sys.argv[1], sys.argv[2])
    print('%-46s %-7s n3 rescale min/max = %.4f / %.4f'
          % (os.path.basename(sys.argv[1]), sys.argv[2], lo, hi))
