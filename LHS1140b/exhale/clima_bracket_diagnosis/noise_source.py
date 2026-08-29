"""Where does the outer residual's noise come from?

`make_profile_bg_gas` terminates on a normalized surface-pressure tolerance
of 1e-10, so the background pressure it returns is a quantized function of
the trial surface temperature.  This walks T_surf in millikelvin steps and
records both the background pressure the inner solve settles on and the
flux residual that follows, so a jump in one can be matched to a jump in
the other.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1]) if len(sys.argv) > 1 else 9.4
wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work_noise2',
                  'heh%s' % str(heh).replace('.', 'p'))
c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
c.solve_for_T_trop = True
flux_scale = max(1.0e-6, 1.0e3*abs(c.rad.bolometric_flux()))
T_trop = 120.0
T0 = 400.0

print('# He/H = %s  T_trop = %.1f' % (heh, T_trop))
print('# %14s %22s %22s %14s'
      % ('T_surf - 400 [K]', 'P_bg [dyn/cm2]', 'f1', 'd f1 / d T'))
prev = None
for dT in [0.0, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 1e-2, 1e-1]:
    T = T0 + dT
    c.T_trop = T_trop
    c.make_profile_bg_gas(T, P_i, P_deep, bg)
    Pbg = c.f_i_surf[ind]*c.P_surf
    ISR, OLR = c.TOA_fluxes_bg_gas(T, P_i, P_deep, bg)
    f1 = (ISR - OLR)/flux_scale
    if prev is None:
        slope = float('nan')
    else:
        slope = (f1 - prev)/dT
    print('  %14.6e %22.14e %22.14e %14.6e' % (dT, Pbg, f1, slope))
    if prev is None:
        prev = f1
    sys.stdout.flush()
