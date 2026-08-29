"""How large a temperature step does the flux residual need to be a derivative?

`surface_temperature_bg_gas` builds its Jacobian by forward differences at
MINPACK's default step, sqrt(eps)*|x| in log10 T, i.e. a fraction of a
millikelvin.  This measures the residual's response to a range of steps
around the starting point, so the step at which the difference stops being
round-off can be read off.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1]) if len(sys.argv) > 1 else 9.4
T0, T_trop = 400.0, 120.0
wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work_noise',
                  'heh%s' % str(heh).replace('.', 'p'))
c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
c.solve_for_T_trop = True
c.T_trop = T_trop
flux_scale = max(1.0e-6, 1.0e3*abs(c.rad.bolometric_flux()))


def f1(T):
    c.T_trop = T_trop
    ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
    return (ISR - OLR)/flux_scale


x0 = np.log10(T0 - T_trop)
base = f1(T0)
print('# He/H = %s  T_surf = %.1f  T_trop = %.1f  f1 = %.12e'
      % (heh, T0, T_trop, base))
h_minpack = np.sqrt(np.finfo(float).eps)*abs(x0)
print('# MINPACK default step in x1 = %.4e  ->  dT_surf = %.4e K'
      % (h_minpack, (10.0**(x0 + h_minpack) - 10.0**x0)))
print('# %12s %14s %16s %16s' % ('h(x1)', 'dT_surf [K]', 'df1', 'df1/dx1'))
for h in [h_minpack, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
    dT = 10.0**(x0 + h) - 10.0**x0
    df = f1(T0 + dT) - base
    print('  %12.4e %14.4e %16.6e %16.6e' % (h, dT, df, df/h))
    sys.stdout.flush()
