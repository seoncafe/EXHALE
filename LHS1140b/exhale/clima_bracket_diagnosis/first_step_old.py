"""The same first-step measurement on the 0.8.4 build, in its own variables.

0.8.4 solves in log10(T_surf) and log10(T_trop) with unscaled residuals.  If
its Jacobian is equally noisy at MINPACK's default step, then the defect
predates the patch and the patch only moved which compositions are unlucky.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

T0, T_trop0 = bl.T_DEEP_GUESS, bl.T_TROP_GUESS
STEPS = [np.sqrt(np.finfo(float).eps), 1.0e-4]

for heh in [float(a) for a in sys.argv[1:]]:
    wd = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'work_step_old', 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    c.solve_for_T_trop = True

    def F(x):
        T = 10.0**x[0]
        T_trop = 10.0**x[1]
        c.T_trop = T_trop
        ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
        A = c.rad.wrk_sol.fup_n[-1]/c.rad.wrk_sol.fdn_n[-1]
        return np.array([ISR - OLR, c.rad.skin_temperature(A) - T_trop])

    x0 = np.array([np.log10(T0), np.log10(T_trop0)])
    f0 = F(x0)
    print('# He/H = %-6s f = [%.6e, %.6e]' % (heh, f0[0], f0[1]))
    for rel in STEPS:
        J = np.empty((2, 2))
        for j in range(2):
            h = rel*abs(x0[j])
            xp = x0.copy(); xp[j] += h
            J[:, j] = (F(xp) - f0)/h
        dx = np.linalg.solve(J, -f0)
        xn = x0 + dx
        print('  rel step %-10.3e J = [[%11.4e, %11.4e], [%11.4e, %11.4e]]'
              '  -> T_surf = %12.4g K, T_trop = %9.4f K %s'
              % (rel, J[0, 0], J[0, 1], J[1, 0], J[1, 1],
                 10.0**xn[0], 10.0**xn[1],
                 'T_surf < T_trop' if xn[0] < xn[1] else ''))
    sys.stdout.flush()
