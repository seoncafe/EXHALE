"""The residual of `make_profile_bg_gas` across the patch's scan interval.

Prints the 49-point table the patched Fortran walks, at a stated surface
temperature and tropopause temperature, so a missing sign change can be told
apart from a sign change the grid steps over and from one the invalid-state
skip removes.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1])
T_surf = float(sys.argv[2]) if len(sys.argv) > 2 else 400.0
T_trop = float(sys.argv[3]) if len(sys.argv) > 3 else 120.0
nscan = int(sys.argv[4]) if len(sys.argv) > 4 else bl.NSCAN

wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work',
                  'heh%s' % str(heh).replace('.', 'p'))
c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
c.T_trop = T_trop
print('# He/H = %s  bg = %s (index %d)  P_surf target = %.6e dyn/cm2'
      % (heh, bg, ind, P_deep))
print('# T_surf = %.4f  T_trop = %.4f  nscan = %d' % (T_surf, T_trop, nscan))
print('# mixing: ' + ', '.join('%s=%.4e' % kv for kv in mixing.items()))
print('# %3s %14s %16s %6s  %s' % ('i', 'log10(P_bg)', 'f', 'valid', 'note'))
rows = bl.scan(c, P_i, ind, P_deep, T_surf, nscan=nscan)
for i, (x, f, ok, msg) in enumerate(rows):
    print('  %3d %14.6f %16.6e %6s  %s'
          % (i, x, f, 'T' if ok else 'F', msg[:90]))
st, det = bl.bracket_from_scan(rows)
print('# bracket status: %s  %s' % (st, det))
