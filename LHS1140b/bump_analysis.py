#!/usr/bin/env python3
"""Look for the Schulik & Owen (2025) adiabatic-cooling He 2^3S population
bump in the converged EXHALE solutions for LHS 1140 b, and measure the
velocity at which the bump sits.

Run from LHS1140b/.  Writes a table to stdout.  make_memo_figures.py
imports this module directly for its bump figure, so a change to RUNDIR or
to load() below moves that figure with it.

The rate coefficients below are transcribed from the EXHALE source so that
the reference curves use the code's own temperature dependence:

  alpha(2^3S)  rec_HeII_23S        src/modules/radiation/Cool_coeff.f90:985
  q13          coex_HeI_1S_23S     src/modules/radiation/Cool_coeff.f90:898
  q31a         coex_HeI_23S_21S    src/modules/radiation/Cool_coeff.f90:920
  q31b         coex_HeI_23S_21P    src/modules/radiation/Cool_coeff.f90:943
  b_23S        ci_HeI23S           src/modules/radiation/Cool_coeff.f90:1091
  Q31          ioniz_HeI23S_H      src/modules/radiation/Cool_coeff.f90:1040
  Q31(H2)      ioniz_HeI23S_H2     src/modules/radiation/Cool_coeff.f90:1076
  A31 = 1.272e-4 s^-1              src/modules/radiation/util_ion_eq.f90:1235

and the steady-state 2^3S balance they enter is tr_triplet_row,
src/modules/nonlinear_system_solver/ion_residual_core.f90:81.
"""
import os
import sys

import numpy as np

# The EXHALE profile files carry two GHOST rows at each end -- a fixed base
# state below, a zero-gradient / WENO3 extrapolation above -- and they are not
# solution cells; exhale_io.loadtxt_cells drops them.
sys.path.insert(0, os.path.join('..', 'examples'))
from exhale_io import loadtxt_cells                          # noqa: E402

KB_EV = 8.617333262e-05          # parameters.f90:380
ERG2EV = 6.241509075e11          # parameters.f90:387
E_TH_HETR = 4.80                 # eV, parameters.f90:414
A31 = 1.272e-4                   # s^-1, util_ion_eq.f90:1235


def alpha_23S(T):
    return 2.10e-13*(T/1.0e4)**(-0.778)


def _ups(T, a, b, c, d):
    return a*np.exp(b*T) + c*np.exp(d*T)


def q13(T):
    ups = _ups(T, 0.06703, -2.673e-6, -0.03584, -3.880e-4)
    return 2.10e-8*np.sqrt(13.60/(KB_EV*T))*np.exp(-19.81/(KB_EV*T))*ups


def q31a(T):
    ups = _ups(T, 2.847, -1.252e-5, -1.953, -3.558e-4)
    return 2.10e-8*np.sqrt(13.60/(KB_EV*T))*np.exp(-0.80/(KB_EV*T))*ups/3.0


def q31b(T):
    ups = _ups(T, 1.185, -4.749e-6, -0.9131, -1.669e-4)
    return 2.10e-8*np.sqrt(13.60/(KB_EV*T))*np.exp(-1.40/(KB_EV*T))*ups/3.0


def ci_23S(T):
    return 6.41e-21*np.sqrt(T)*np.exp(-55338.0/T)/(E_TH_HETR/ERG2EV)


# Branching of the He(2^3S) + neutral total ionization between the Penning
# and the associative channel, Garcia Munoz (2025) "an average 0.9:0.1"
# (Cool_coeff.f90:552).  The metastable balance takes the TOTAL -- both
# channels quench it -- so the sink below is unscaled; the factor is kept
# here because it is what the heating and proton-source terms carry.
F_PENNING_HEI23S = 0.9


def ioniz_HeI23S_H(T):
    """TOTAL He(2^3S) + H ionization rate coefficient [cm^3 s^-1], Penning
    plus associative, as the Garcia Munoz (2025, A&A 698, A199) closed form
    of the Movre & Meyer (1997) cross sections.  This is Q31 in
    tr_triplet_row."""
    lnT = np.log(T)
    return 1.0e-9*np.exp(-8.64804e1/T - 2.86766e-1*lnT
                         + 8.68445e-2*lnT**2 - 5.73001e-3*lnT**3)


def ioniz_HeI23S_H2(T):
    """TOTAL He(2^3S) + H2 ionization rate coefficient [cm^3 s^-1], Cohen &
    Lane (1977) cross sections through the Garcia Munoz (2025) network file.
    Not exercised by the LHS 1140 b runs analyzed here, which carry no
    molecular chemistry, but transcribed with the atomic channel so the two
    stay together."""
    return 5.408222e-12 * T**6.75388e-1 * np.exp(-6.96275e2/T)


def load(tag, adv):
    """Return a dict of the converged profile.  adv=True reads the
    advection-corrected files that EXHALE_transit.py uses."""
    sfx = '_adv' if adv else ''
    base = os.path.join(RUNDIR, tag, 'output')
    hy = loadtxt_cells(os.path.join(base, 'Hydro_ioniz%s.txt' % sfx))
    fi = os.path.join(base, 'Ion_species%s.txt' % sfx)
    cols = [l for l in open(fi) if l.startswith('# columns')][0].split()[2:]
    k = {n: j for j, n in enumerate(cols)}
    io = loadtxt_cells(fi)
    d = dict(r=hy[:, 0], v=hy[:, 2], T=hy[:, 4])
    for s in ('HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR'):
        d[s] = io[:, k[s]]
    # No metals.inp / no molecules in these runs, so the electron budget is
    # exactly the H + He one the solver uses.
    d['ne'] = d['HII'] + d['HeII'] + 2.0*d['HeIII']
    d['nHe'] = d['HeI'] + d['HeII'] + d['HeIII']
    # nhei written out is total neutral He including the metastable
    # (ionization_equilibrium.f90:891), so the ground singlet is the remainder.
    d['HeI_SI'] = np.maximum(d['HeI'] - d['HeITR'], 0.0)
    d['f3'] = d['HeITR']/np.maximum(d['nHe'], 1e-300)
    return d


def gamma_tr(d):
    """Back out the 2^3S photoionization rate [s^-1] from the equilibrium
    solution: it is the only term of tr_triplet_row not reconstructible from
    the output files, and the equilibrium file satisfies that row."""
    T, ne = d['T'], d['ne']
    src = ne*(d['HeII']*alpha_23S(T) + d['HeI_SI']*q13(T))
    sink_known = (A31 + d['HI']*ioniz_HeI23S_H(T)
                  + ne*(q31a(T) + q31b(T) + ci_23S(T)))
    with np.errstate(divide='ignore', invalid='ignore'):
        g = src/np.maximum(d['HeITR'], 1e-300) - sink_known
    return g


def n23_reference(d, g_tr, T_q):
    """Steady-state 2^3S density with the collisional and recombination rates
    evaluated at T_q instead of the local temperature (S&O's 'const q'
    reference when T_q = T_max); the ionization terms are kept local."""
    ne = d['ne']
    src = ne*(d['HeII']*alpha_23S(T_q) + d['HeI_SI']*q13(T_q))
    sink = (A31 + d['HI']*ioniz_HeI23S_H(T_q) + g_tr
            + ne*(q31a(T_q) + q31b(T_q) + ci_23S(T_q)))
    return src/np.maximum(sink, 1e-300)


def bump_metrics(d, rmin=1.05, rmax=25.0):
    """Locate the 2^3S maximum and characterize the rise in f3."""
    r, T, v = d['r'], d['T'], d['v']/1e5
    m = (r >= rmin) & (r <= rmax)
    i = np.where(m)[0][np.argmax(d['HeITR'][m])]
    j = np.where(m)[0][np.argmax(d['f3'][m])]
    iT = np.where(m)[0][np.argmax(T[m])]
    # Does f3 rise outward over the cooling leg (from T_max out to the 2^3S
    # peak)?  That rise is the bump in S&O's sense.
    leg = slice(iT, i + 1)
    f3rise = (d['f3'][i]/max(d['f3'][iT], 1e-300)) if i > iT else np.nan
    return dict(i=i, r_peak=r[i], n_peak=d['HeITR'][i], T_peak=T[i],
                v_peak=v[i], j=j, r_f3=r[j], f3_max=d['f3'][j],
                T_f3=T[j], v_f3=v[j],
                T_max=T[iT], r_Tmax=r[iT], iT=iT, f3_rise=f3rise,
                f3_rise_j=d['f3'][j]/max(d['f3'][iT], 1e-300),
                f3_monotonic=bool(np.all(np.diff(d['f3'][leg]) > -1e-30)))


def velocity_distribution(d, rmin=1.0, rmax=25.0):
    """Absorber-weighted distribution of |v_r| over the 2^3S column,
    w = n_23S r^2 dr (a radial-projection upper bound: the true line-of-sight
    projection of a radial field is smaller)."""
    r, v = d['r'], np.abs(d['v'])/1e5
    m = (r >= rmin) & (r <= rmax)
    r, v, n3 = r[m], v[m], d['HeITR'][m]
    w = n3*r**2
    dr = np.gradient(r)
    w = w*dr
    o = np.argsort(v)
    v, w = v[o], w[o]
    cum = np.cumsum(w)/np.sum(w)
    frac = {t: float(np.sum(w[v > t])/np.sum(w)) for t in (2.0, 5.0, 10.0)}
    q = {p: float(np.interp(p, cum, v)) for p in (0.5, 0.9, 0.99)}
    vmean = float(np.sum(w*v)/np.sum(w))
    vrms = float(np.sqrt(np.sum(w*v**2)/np.sum(w)))
    return dict(v=v, cum=cum, frac=frac, vmean=vmean, vrms=vrms, q=q,
                vmax=float(v.max()))


# The well-mixed cases of the model tree of record (MODELS.md), read as
# `<group>/HeH<value>` under `models/`.
RUNDIR = 'models'

_WM = 'atomic_scalar_gj1132_wellmixed'
_G699 = 'atomic_scalar_gj699_wellmixed'

CASES = [(_WM + '/HeH0.55', 'GJ 1132 SED, He/H = 0.55'),
         (_WM + '/HeH0.083', 'GJ 1132 SED, solar He/H'),
         (_WM + '/HeH1000', 'GJ 1132 SED, He/H = 1000'),
         (_G699 + '/HeH0.050', 'GJ 699 SED, He/H = 0.050'),
         (_G699 + '/HeH0.083', 'GJ 699 SED, solar He/H'),
         (_G699 + '/HeH1000', 'GJ 699 SED, He/H = 1000')]


def main():
    out = {}
    hdr = ('%-40s %-4s %7s %8s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s'
           % ('case', 'prof', 'r_n23', 'n_peak', 'T_n23', 'v_n23',
              'r_f3', 'T_f3', 'v_f3', 'T_max', '<|v|>', '>2', '>5', '>10'))
    print(hdr)
    print('-'*len(hdr))
    for tag, _ in CASES:
        for adv in (True, False):
            d = load(tag, adv)
            b = bump_metrics(d)
            vd = velocity_distribution(d)
            g = gamma_tr(d) if not adv else None
            print('%-40s %-4s %7.2f %8.3g %7.0f %7.2f %7.2f %7.0f %7.2f '
                  '%7.0f %7.2f %6.1f%% %6.1f%% %6.1f%%'
                  % (tag, 'adv' if adv else 'eq', b['r_peak'], b['n_peak'],
                     b['T_peak'], b['v_peak'], b['r_f3'], b['T_f3'],
                     b['v_f3'], b['T_max'], vd['vmean'],
                     100*vd['frac'][2.0], 100*vd['frac'][5.0],
                     100*vd['frac'][10.0]))
            out['%s_%s' % (tag, 'adv' if adv else 'eq')] = (d, b, vd, g)
    print()
    print('f3 rise over the cooling leg (T_max -> 2^3S peak), '
          'and mdot_23S drop past the peak:')
    for tag, _ in CASES:
        for adv in (True, False):
            d, b, vd, g = out['%s_%s' % (tag, 'adv' if adv else 'eq')]
            r, v, n3 = d['r'], d['v'], d['HeITR']
            mdot = 4*np.pi*r**2*v*n3
            i = b['i']
            outer = (r > b['r_peak']) & (r <= 25.0)
            drop = (mdot[outer].min()/mdot[i]) if mdot[i] > 0 else np.nan
            print('  %-40s %-3s  f3(peak)/f3(T_max) = %6.2f  '
                  'monotonic rise: %-5s  mdot_23S(min past peak)/mdot(peak) '
                  '= %6.3f' % (tag, 'adv' if adv else 'eq', b['f3_rise_j'],
                               b['f3_monotonic'], drop))
    print()
    print('velocity percentiles of the 2^3S column [km/s] '
          '(w = n_23S r^2 dr, |v_r|):')
    for tag, _ in CASES:
        for adv in (True, False):
            d, b, vd, g = out['%s_%s' % (tag, 'adv' if adv else 'eq')]
            print('  %-40s %-3s  median %5.2f  90%% %5.2f  99%% %5.2f  '
                  'max %5.2f  rms %5.2f'
                  % (tag, 'adv' if adv else 'eq', vd['q'][0.5], vd['q'][0.9],
                     vd['q'][0.99], vd['vmax'], vd['vrms']))
    print()
    print('const-q reference test (equilibrium profiles; rates frozen at '
          'T_max vs local T):')
    for tag, _ in CASES:
        d, b, vd, g = out['%s_eq' % tag]
        n_var = d['HeITR']
        n_const = n23_reference(d, g, np.full_like(d['T'], b['T_max']))
        n_loc = n23_reference(d, g, d['T'])
        m = (d['r'] > b['r_Tmax']) & (d['r'] <= 25.0)
        print('  %-40s local-T reference reproduces solved n23 to %.2e '
              '(max rel. dev.);  n_var/n_const at 2^3S peak = %6.2f, '
              'max over r > r(T_max) = %6.2f; gamma_TR(peak) = %.3e s^-1'
              % (tag,
                 np.max(np.abs(n_loc[m]/np.maximum(n_var[m], 1e-300) - 1)),
                 n_var[b['i']]/max(n_const[b['i']], 1e-300),
                 np.max(n_var[m]/np.maximum(n_const[m], 1e-300)),
                 g[b['i']]))
    return out


if __name__ == '__main__':
    main()
