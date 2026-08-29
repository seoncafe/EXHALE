"""Where in (T_surf, T_trop) does the bracketing scan lose its sign change?

For each trial temperature the outer `surface_temperature_bg_gas` solve can
visit, replay the patch's 49-point scan and report: the residual at the top
of the interval, how many scan points `make_profile` refused, and whether a
bracket survives.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1])
T_trop = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
Ts = [float(a) for a in sys.argv[3:]] or list(np.arange(150.0, 901.0, 25.0))

wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work',
                  'heh%s' % str(heh).replace('.', 'p'))
c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
c.T_trop = T_trop
print('# He/H = %s  bg = %s  T_trop = %.4f  P_surf = %.6e' % (heh, bg, T_trop, P_deep))
print('# %9s %8s %16s %16s %s'
      % ('T_surf', 'n_invalid', 'f_last_valid', 'f_top', 'status'))
for T in Ts:
    rows = bl.scan(c, P_i, ind, P_deep, T)
    ninv = sum(1 for r in rows if not r[2])
    valid = [r for r in rows if r[2]]
    flast = valid[-1][1] if valid else float('nan')
    ftop = rows[-1][1] if rows[-1][2] else float('nan')
    st, det = bl.bracket_from_scan(rows)
    print('  %9.4f %8d %16.6e %16.6e %s' % (T, ninv, flast, ftop, st))
    sys.stdout.flush()
