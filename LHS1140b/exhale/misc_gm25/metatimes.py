#!/usr/bin/env python3
"""He 2^3S source and sink timescales on a solved profile (memo Table 4).

  t_rec     = 1/(n_e alpha_2^3S),  alpha = 2.10e-13 (T/1e4)^-0.778
  t_destr   = 1/(n_HI Q31),        Q31 the He(2^3S)+H total ionization rate
  t_life    = 1/(1/t_destr + A31), A31 = 1.272e-4 s^-1
  L         = 9.34 km/s * t_life  [R_p]

Q31 is either the Garcia Munoz (2025) closed form now in Cool_coeff.f90
(`gm25`) or the retired Taylor et al. (2025) Table 2 two-branch fit
(`taylor`), so that the same profile can be read on either coefficient.

usage: ./metatimes.py <run_dir> [gm25|taylor] [output_subdir]
"""
import os, re, sys
import numpy as np

A31 = 1.272e-4
RP_CM = 0.157692*7.1492e9      # LHS 1140 b radius, R_J -> cm
V_KMS = 9.34


def q31_gm25(T):
    c, d1, d2, d3 = -8.64804e1, -2.86766e-1, 8.68445e-2, -5.73001e-3
    lt = np.log(T)
    return 1e-9*np.exp(c/T + d1*lt + d2*lt**2 + d3*lt**3)


def q31_taylor(T):
    return np.where(T <= 4000.0, 1.9e-9*(300.0/T)**0.07,
                    9.1e-9*(300.0/T)**0.50)


def cols(p):
    return open(p).readlines()[1].split()[2:]


def table(d, rate='gm25', sub='output'):
    o = os.path.join(d, sub)
    hn, sn = cols(o+'/Hydro_ioniz_adv.txt'), cols(o+'/Ion_species_adv.txt')
    h, a = np.loadtxt(o+'/Hydro_ioniz_adv.txt'), np.loadtxt(o+'/Ion_species_adv.txt')
    r, T = h[:, 0], h[:, hn.index('T[K]')]
    nh = h[:, hn.index('rho[mH/cm3]')]        # not used; kept for reference
    nhi = a[:, sn.index('HI')]
    # electrons: every ionization stage the file carries, with its charge
    ne = np.zeros(a.shape[0])
    for i, c in enumerate(sn):
        m = re.match(r'^([A-Z][a-z]?)(I+)$', c)
        if m and m.group(2) != 'I':
            ne += (len(m.group(2)) - 1)*a[:, i]
    q = q31_gm25(T) if rate == 'gm25' else q31_taylor(T)
    trec = 1.0/(ne*2.10e-13*(T/1e4)**-0.778)
    tdes = 1.0/(nhi*q)
    tlife = 1.0/(1.0/tdes + A31)
    out = []
    for rr in (1.2, 2.0, 3.0, 5.0, 8.0, 20.0):
        f = lambda y: float(np.interp(rr, r, y))
        tl = f(tlife)
        out.append((rr, f(T), f(trec), f(tdes), tl, V_KMS*1e5*tl/RP_CM))
    return out


if __name__ == '__main__':
    d = sys.argv[1]
    rate = sys.argv[2] if len(sys.argv) > 2 else 'gm25'
    sub = sys.argv[3] if len(sys.argv) > 3 else 'output'
    print('%s  [%s, %s]' % (d, rate, sub))
    print('%6s %8s %12s %12s %12s %10s'
          % ('r[Rp]', 'T[K]', 't_rec[s]', 't_destr[s]', 't_life[s]', 'L[Rp]'))
    for row in table(d, rate, sub):
        print('%6.1f %8.0f %12.3g %12.3g %12.3g %10.4g' % row)
