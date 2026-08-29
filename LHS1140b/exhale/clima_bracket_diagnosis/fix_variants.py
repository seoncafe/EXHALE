"""Candidate repairs, measured in a replication of the same outer solve.

The replication reproduces the installed build: it drives the patched
residual -- log10(T_surf - T_trop), log10(T_trop), both residuals normalized
-- through the same MINPACK `hybrd`, and it reproduces the installed
build's converged deep temperature at the compositions that solve.  Nothing
in `photochem/` or in either installed environment is touched.

variants
  patched   what is installed: forward differences at the MINPACK default
            step, initial trust-region factor 100
  epsfcn    the same, with `epsfcn = 1e-4`, i.e. a relative difference
            step of 1e-2 in the log variable (MINPACK's `fdjac1` steps
            `sqrt(max(epsfcn, epsmch))*|x|`)
  factor    the same, with the initial trust-region bound at factor 1
  bounded   the trial surface temperature confined to
            (T_trop, T_max) by a logistic parameterization, T_max the
            upper end of the thermodynamic data
"""
import os
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

T_MAX = 6000.0          # measured ceiling of `make_profile`, T_ceiling.py

VARIANT = sys.argv[1]
HEHS = [float(a) for a in sys.argv[2:]]


def solve(c, P_i, bg, P_deep, variant):
    flux_scale = max(1.0e-6, 1.0e3*abs(c.rad.bolometric_flux())
                     + abs(c.surface_heat_flow))
    state = {'n': 0, 'err': None, 'T': np.nan, 'T_trop': np.nan}

    def unpack(x):
        T_trop = 10.0**x[1]
        if variant == 'bounded':
            u = 1.0/(1.0 + np.exp(-x[0]))
            T = T_trop + (T_MAX - T_trop)*u
        else:
            T = T_trop + 10.0**x[0]
        return T, T_trop

    def fcn(x):
        T, T_trop = unpack(x)
        c.T_trop = T_trop
        state['n'] += 1
        try:
            ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
        except Exception as e:
            state['err'] = str(e)
            raise StopIteration
        f1 = (ISR - OLR + c.surface_heat_flow)/flux_scale
        A = c.rad.wrk_sol.fup_n[-1]/c.rad.wrk_sol.fdn_n[-1]
        f2 = (c.rad.skin_temperature(A) - T_trop)/max(T_trop, 1.0)
        state['T'], state['T_trop'] = T, T_trop
        return np.array([f1, f2])

    T_g, T_t = bl.T_DEEP_GUESS, bl.T_TROP_GUESS
    if variant == 'bounded':
        u = (T_g - T_t)/(T_MAX - T_t)
        x0 = np.array([np.log(u/(1.0 - u)), np.log10(T_t)])
    else:
        x0 = np.array([np.log10(T_g - T_t), np.log10(T_t)])
    opts = {}
    if variant == 'epsfcn':
        opts['eps'] = 1.0e-4
    if variant == 'factor':
        opts['factor'] = 1.0
    try:
        s = root(fcn, x0, method='hybr', tol=np.sqrt(np.finfo(float).eps),
                 options=opts)
    except StopIteration:
        return None, state
    T, T_trop = unpack(s.x)
    ok = s.success and max(abs(s.fun)) < 1.0e-7
    return (T, T_trop, ok, state['n']), state


print('# variant = %s' % VARIANT)
for heh in HEHS:
    # one work directory per variant: the climate species/settings files are
    # rewritten on every build, so concurrent variants must not share a path
    wd = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'work_' + VARIANT, 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    c.T_trop = bl.T_TROP_GUESS
    c.solve_for_T_trop = True
    res, state = solve(c, P_i, bg, P_deep, VARIANT)
    if res is None:
        print('%-8s FAIL  %s' % (heh, state['err'][:70]))
    else:
        T, T_trop, ok, n = res
        print('%-8s %-4s T_deep=%10.5f  T_trop=%9.5f  nfev=%d'
              % (heh, 'OK' if ok else 'STALL', T, T_trop, n))
    sys.stdout.flush()
