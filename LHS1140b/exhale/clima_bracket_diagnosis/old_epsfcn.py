"""Does the same larger finite-difference step cure the 0.8.4 failures?

Runs 0.8.4's own outer solve -- log10(T_surf), log10(T_trop), unscaled
residuals -- over the composition grid, once at MINPACK's default step and
once at 1e-4, so the repair can be told apart from the parameterization the
patch changed.
"""
import os
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

rel = float(sys.argv[1])
for heh in [float(a) for a in sys.argv[2:]]:
    wd = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'work_old_eps%g' % rel, 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    c.solve_for_T_trop = True
    err = {}

    def fcn(x):
        T, T_trop = 10.0**x[0], 10.0**x[1]
        c.T_trop = T_trop
        try:
            ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
        except Exception as e:
            err['m'] = str(e)
            raise StopIteration
        A = c.rad.wrk_sol.fup_n[-1]/c.rad.wrk_sol.fdn_n[-1]
        return np.array([ISR - OLR, c.rad.skin_temperature(A) - T_trop])

    x0 = np.array([np.log10(bl.T_DEEP_GUESS), np.log10(bl.T_TROP_GUESS)])
    opts = {} if rel <= 0 else {'eps': rel}
    try:
        s = root(fcn, x0, method='hybr', tol=np.sqrt(np.finfo(float).eps),
                 options=opts)
    except StopIteration:
        print('%-8s FAIL  %s' % (heh, err['m'][:70]))
        sys.stdout.flush()
        continue
    ok = s.success and abs(s.fun[0]) < 1.0 and abs(s.fun[1]) < 1e-4
    print('%-8s %-5s T_deep=%10.5f  T_trop=%9.5f'
          % (heh, 'OK' if ok else 'STALL', 10.0**s.x[0], 10.0**s.x[1]))
    sys.stdout.flush()
