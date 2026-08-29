"""The first Newton step of the outer solve, at two finite-difference steps.

The step is what decides whether the trial surface temperature stays in the
range of the thermodynamic data.  It is built from a 2x2 forward-difference
Jacobian in (log10(T_surf - T_trop), log10(T_trop)); the cross term
df1/dx2 is tens of times the diagonal term df1/dx1, so the step's x1
component is the cross term divided by a small, noisy number.
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
                      'work_step', 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    c.solve_for_T_trop = True
    flux_scale = max(1.0e-6, 1.0e3*abs(c.rad.bolometric_flux()))

    def F(x):
        T_trop = 10.0**x[1]
        T = T_trop + 10.0**x[0]
        c.T_trop = T_trop
        ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
        f1 = (ISR - OLR)/flux_scale
        A = c.rad.wrk_sol.fup_n[-1]/c.rad.wrk_sol.fdn_n[-1]
        f2 = (c.rad.skin_temperature(A) - T_trop)/max(T_trop, 1.0)
        return np.array([f1, f2])

    x0 = np.array([np.log10(T0 - T_trop0), np.log10(T_trop0)])
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
        T_trop_new = 10.0**xn[1]
        T_new = T_trop_new + 10.0**xn[0]
        print('  rel step %-10.3e J = [[%11.4e, %11.4e], [%11.4e, %11.4e]]'
              % (rel, J[0, 0], J[0, 1], J[1, 0], J[1, 1]))
        print('  %-21s dx = [%9.4f, %9.4f]  ->  T_surf = %12.4g K, '
              'T_trop = %9.4f K   %s'
              % ('', dx[0], dx[1], T_new, T_trop_new,
                 'OUTSIDE the 6000 K data range' if T_new > 6000.0 else ''))
    sys.stdout.flush()
