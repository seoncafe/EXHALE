"""What the 0.8.4 build does at the same compositions.

Two things are measured: the solution `surface_temperature_bg_gas` returns
(so its background pressure can be checked against the interval the new
bracketing scans), and the trial temperatures its outer MINPACK solve
visits, using its own parameterization -- log10(T_surf), log10(T_trop),
unscaled flux and temperature residuals -- replayed through scipy's `hybr`.
"""
import os
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

MODE = sys.argv[1]                      # 'solve' or 'trace'
for heh in [float(a) for a in sys.argv[2:]]:
    wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work_old',
                      'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    if MODE == 'solve':
        try:
            T = c.surface_temperature_bg_gas(P_i, P_deep, bg,
                                             T_guess=bl.T_DEEP_GUESS)
            Pbg = c.f_i_surf[ind]*c.P_surf
            x = np.log10(Pbg)
            print('%-8s OK   T_deep=%10.5f  T_trop=%9.5f  P_bg=%.8e  '
                  'log10(P_bg)=%.6f  interval=[%.6f, %.6f]  inside=%s'
                  % (heh, T, c.T_trop, Pbg, x, np.log10(P_deep) - 12.0,
                     np.log10(P_deep), P_deep*1e-12 <= Pbg <= P_deep))
        except Exception as e:
            print('%-8s FAIL %s' % (heh, e))
    else:
        c.T_trop = bl.T_TROP_GUESS
        c.solve_for_T_trop = True
        n = [0]

        def fcn(x):
            T = 10.0**x[0]
            T_trop = 10.0**x[1]
            c.T_trop = T_trop
            n[0] += 1
            try:
                ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
            except Exception as e:
                print('  %4d %14.5f %12.5f  %s' % (n[0], T, T_trop, str(e)[:70]))
                raise SystemExit(0)
            f1 = ISR - OLR + c.surface_heat_flow
            A = c.rad.wrk_sol.fup_n[-1]/c.rad.wrk_sol.fdn_n[-1]
            f2 = c.rad.skin_temperature(A) - T_trop
            print('  %4d %14.5f %12.5f %14.6e %14.6e' % (n[0], T, T_trop, f1, f2))
            sys.stdout.flush()
            return np.array([f1, f2])

        print('# He/H = %s  old parameterization, unscaled residuals' % heh)
        x0 = np.array([np.log10(bl.T_DEEP_GUESS), np.log10(bl.T_TROP_GUESS)])
        s = root(fcn, x0, method='hybr', tol=np.sqrt(np.finfo(float).eps))
        print('# converged=%s  T_surf=%.5f  T_trop=%.5f'
              % (s.success, 10.0**s.x[0], 10.0**s.x[1]))
    sys.stdout.flush()
