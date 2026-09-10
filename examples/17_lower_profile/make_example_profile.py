#!/usr/bin/env python3
"""Write the synthetic lower-atmosphere profile of this example.

The file this writes is NOT a photochemistry result: it is a hand-built,
hydrostatically integrated column with a smooth H2 -> H front, whose only
purpose is to exercise and document the schema of
docs/phase_e_flux_closure_design.md section 2 (the reader, the K_zz profile
path, the elemental reservoirs and the refusal rules).  The production
producer is the Photochem adapter of milestone E2.

Everything below is stated explicitly so the file can be regenerated and
diffed:

  planet      HD 209458 b, the planet of examples/14_diffusion
  match       p_match = 1e-6 bar, where the example places the EXHALE base
              and where r = 1.401 R_J, the planet radius of that example
  H2 front    the H nuclei fraction bound into H2 goes through 1/2 at
              p = 10^-5.7 bar
  T(p)        a monotone radiative profile through T = 1450 K at the match
  K_zz(p)     rising outward as p^-0.35, 1e8 cm^2/s at 1e-3 bar
  elements    solar C/N/O and He/H = 0.0833, constant with height, because
              photochemistry moves nuclei between molecules without creating
              or destroying them
  F_H, F_He   a constant trial elemental flux (the first arm of the closure
              holds helium at zero)
"""
import hashlib
import os
import sys

import numpy as np

# ---- constants (CODATA / IAU, cgs) ------------------------------------- #
# The Jupiter radius and mass are EXHALE's own values, imported from the one
# Python definition (examples/exhale_io.py, which carries the IAU 2015
# nominal values of parameters.f90), so that a radius written here means the
# same length the code reads back.  The directory is APPENDED to sys.path so
# nothing in it can shadow a standard-library module.
_EXAMPLES_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
if _EXAMPLES_DIR not in sys.path:
    sys.path.append(_EXAMPLES_DIR)
from exhale_io import RJ_CM, MJ as MJ_G                     # noqa: E402

KB = 1.380649e-16
MAMU = 1.66053906660e-24
GNEWT = 6.67430e-8
RJ = RJ_CM
MJ = MJ_G

# ---- the configuration this file is a solution of ---------------------- #
CFG = dict(
    planet='HD209458b',
    Mp_MJ=0.720,
    r_match_RJ=1.401,
    p_match_bar=1.0e-6,
    p_deep_bar=1.0e-3,
    p_top_bar=1.0e-8,
    nlev=51,
    HeH=0.083333333,
    X_C=2.70e-4,
    X_N=6.80e-5,
    X_O=4.90e-4,
    T_match_K=1450.0,
    trial_flux_H=5.0e10,
    trial_flux_He=0.0,
    iteration=0,
)

M_HE = 4.002602
M_H = 1.008


def h2_fraction(x):
    """Fraction of the H nuclei bound into H2, x = log10(p/bar)."""
    return 0.5*(1.0 + np.tanh((x + 5.7)/0.5))


def temperature(x):
    return CFG['T_match_K'] + 120.0*np.tanh((x + 6.0)/1.2)


def eddy(x):
    return 1.0e8*10.0**(-0.35*(x + 3.0))


def build():
    x = np.linspace(np.log10(CFG['p_deep_bar']),
                    np.log10(CFG['p_top_bar']), CFG['nlev'])
    p_bar = 10.0**x
    f_h2 = h2_fraction(x)
    T = temperature(x)
    kzz = eddy(x)

    heh = CFG['HeH']
    # per H nucleus: n_tot/N_H and mass/N_H
    per_H = 1.0 - 0.5*f_h2 + heh
    q_H2 = (0.5*f_h2)/per_H
    q_H = (1.0 - f_h2)/per_H
    mu = (M_H + heh*M_HE)/per_H

    n_tot = p_bar*1.0e6/(KB*T)
    rho = n_tot*mu*MAMU

    # hydrostatic radius, integrated outward and inward from the match
    imatch = int(np.argmin(np.abs(p_bar - CFG['p_match_bar'])))
    r = np.zeros_like(p_bar)
    r[imatch] = CFG['r_match_RJ']*RJ
    lnp = np.log(p_bar)
    mp = CFG['Mp_MJ']*MJ
    for i in range(imatch + 1, len(x)):          # outward (p falling)
        rr = r[i-1]
        for _ in range(4):
            Hs = KB*0.5*(T[i] + T[i-1])*rr*rr/(0.5*(mu[i] + mu[i-1])
                                               * MAMU*GNEWT*mp)
            rr = r[i-1] - Hs*(lnp[i] - lnp[i-1])
        r[i] = rr
    for i in range(imatch - 1, -1, -1):          # inward (p rising)
        rr = r[i+1]
        for _ in range(4):
            Hs = KB*0.5*(T[i] + T[i+1])*rr*rr/(0.5*(mu[i] + mu[i+1])
                                               * MAMU*GNEWT*mp)
            rr = r[i+1] - Hs*(lnp[i] - lnp[i+1])
        r[i] = rr
    r[imatch] = CFG['r_match_RJ']*RJ

    # The match radius is written EXACTLY as configured: (r_match*RJ)/RJ can
    # differ from r_match by one ulp, and the base radius the reader hands to
    # the code must be the number the example states.
    r_RJ = r/RJ
    r_RJ[imatch] = CFG['r_match_RJ']
    # State the match at a level of the table, exactly: a p_match_bar that
    # misses a node by one ulp is interpolated between its neighbours, which
    # is correct but no longer bit-for-bit the node's own value.
    p_bar[imatch] = CFG['p_match_bar']
    cols = {
        'p': p_bar,
        'r': r_RJ,
        'T': T,
        'n_tot': n_tot,
        'rho': rho,
        'Kzz': kzz,
        'q_H2': q_H2,
        'q_H': q_H,
        'X_He': np.full_like(p_bar, heh),
        'X_C': np.full_like(p_bar, CFG['X_C']),
        'X_N': np.full_like(p_bar, CFG['X_N']),
        'X_O': np.full_like(p_bar, CFG['X_O']),
        'F_H': np.full_like(p_bar, CFG['trial_flux_H']),
        'F_He': np.full_like(p_bar, CFG['trial_flux_He']),
        'q_H2O': 3.0e-4*f_h2,
    }
    return cols


def solution_id(cols):
    h = hashlib.sha256()
    h.update(repr(sorted(CFG.items())).encode())
    for name in sorted(cols):
        h.update(name.encode())
        h.update(cols[name].tobytes())
    return h.hexdigest()


def main():
    cols = build()
    names = list(cols)
    sid = solution_id(cols)
    imatch = int(np.argmin(np.abs(cols['p'] - CFG['p_match_bar'])))
    with open('lower_atmosphere_profile.dat', 'w') as f:
        f.write('# EXHALE lower-atmosphere profile'
                ' (docs/phase_e_flux_closure_design.md section 2)\n')
        f.write('# solution_id %s\n' % sid)
        f.write('# source_code analytic\n')
        f.write('# source_version examples/17_lower_profile/'
                'make_example_profile.py\n')
        f.write('# mechanism none (synthetic H2/H/He column,'
                ' no reaction network)\n')
        f.write('# stellar_flux none (no photochemistry was solved)\n')
        f.write('# p_match_bar %.17E\n' % CFG['p_match_bar'])
        f.write('# p_top_bar %.17E\n' % CFG['p_top_bar'])
        f.write('# p_deep_bar %.17E\n' % CFG['p_deep_bar'])
        f.write('# trial_flux_H %.17E\n' % CFG['trial_flux_H'])
        f.write('# trial_flux_He %.17E\n' % CFG['trial_flux_He'])
        f.write('# iteration %d\n' % CFG['iteration'])
        f.write('# reached_steady_state T\n')
        f.write('# notes synthetic example column; T(p) prescribed,'
                ' no climate solution, no cold trap\n')
        f.write('# units p[bar] r[R_J] T[K] n_tot[cm^-3] rho[g/cm^3]'
                ' Kzz[cm^2/s] q_*[-] X_*[El/H nuclei] F_*[g/s]\n')
        f.write('# columns: %s\n' % ' '.join(names))
        for i in range(len(cols['p'])):
            f.write(' '.join('%25.17E' % cols[n][i] for n in names) + '\n')
    print('wrote lower_atmosphere_profile.dat, %d levels, solution_id %s'
          % (len(cols['p']), sid))
    print('match level: p = %.3E bar, r = %.6f R_J, T = %.2f K, q_H2 = %.5f'
          % (cols['p'][imatch], cols['r'][imatch], cols['T'][imatch],
             cols['q_H2'][imatch]))


if __name__ == '__main__':
    main()
