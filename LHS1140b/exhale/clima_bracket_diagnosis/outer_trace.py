"""Which trial (T_surf, T_trop) does the outer solve hand to the bracketing?

`surface_temperature_bg_gas` is a two-variable MINPACK solve in
log10(T_surf - T_trop) and log10(T_trop); the bracket failure is raised from
inside one of its function evaluations, so the trial pair matters.  This
replays the patched outer residual through scipy's `hybr` (the same MINPACK
`hybrd`) with every evaluation logged, and reports the first trial whose
`make_profile_bg_gas` fails.
"""
import os
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1])
T_guess = float(sys.argv[2]) if len(sys.argv) > 2 else bl.T_DEEP_GUESS
T_trop_guess = float(sys.argv[3]) if len(sys.argv) > 3 else bl.T_TROP_GUESS

wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work',
                  'heh%s' % str(heh).replace('.', 'p'))
c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
c.T_trop = T_trop_guess
c.solve_for_T_trop = True

# flux_scale of the patch, evaluated before the solve as the Fortran does
flux_scale = max(1.0e-6, 1.0e3*abs(c.rad.bolometric_flux())
                 + abs(c.surface_heat_flow))
print('# He/H = %s  bg = %s  flux_scale = %.6e' % (heh, bg, flux_scale))
print('# %4s %12s %12s %14s %14s  %s'
      % ('eval', 'T_surf', 'T_trop', 'fvec1', 'fvec2', 'note'))

trials = []


def fcn(x):
    T_trop = 10.0**x[1]
    T = T_trop + 10.0**x[0]
    c.T_trop = T_trop
    try:
        ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
    except Exception as e:
        trials.append((T, T_trop, None, None, str(e)))
        print('  %4d %12.5f %12.5f %14s %14s  %s'
              % (len(trials), T, T_trop, '-', '-', str(e)[:80]))
        sys.stdout.flush()
        raise SystemExit(0)
    f1 = (ISR - OLR + c.surface_heat_flow)/flux_scale
    A = c.rad.wrk_sol.fup_n[-1]/c.rad.wrk_sol.fdn_n[-1]
    f2 = (c.rad.skin_temperature(A) - T_trop)/max(T_trop, 1.0)
    trials.append((T, T_trop, f1, f2, ''))
    print('  %4d %12.5f %12.5f %14.6e %14.6e' % (len(trials), T, T_trop, f1, f2))
    sys.stdout.flush()
    return np.array([f1, f2])


x0 = np.array([np.log10(T_guess - T_trop_guess), np.log10(T_trop_guess)])
sol = root(fcn, x0, method='hybr', tol=np.sqrt(np.finfo(float).eps))
print('# converged=%s  T_trop=%.5f  T_surf=%.5f'
      % (sol.success, 10.0**sol.x[1], 10.0**sol.x[1] + 10.0**sol.x[0]))
