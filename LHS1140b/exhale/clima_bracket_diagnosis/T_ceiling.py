"""The surface temperature above which `make_profile` refuses.

The scan of `make_profile_bg_gas` calls `make_profile` at the trial surface
temperature the outer solve hands it.  Above the range of the thermodynamic
polynomials every scan point is invalid, so the interval holds no valid
point and the bracket cannot exist whatever the pressure grid is.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1]) if len(sys.argv) > 1 else 9.4
wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work',
                  'heh%s' % str(heh).replace('.', 'p'))
c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
c.T_trop = 157.31584

lo, hi = 100.0, 1.0e9
for T in [200., 500., 1000., 2000., 3000., 4000., 5000., 6000., 1e4, 1e5]:
    f, ok, msg = bl.residual(c, P_i, ind, P_deep, T, np.log10(P_deep))
    print('T_surf = %12.4g   make_profile %s  %s'
          % (T, 'ok' if ok else 'REFUSED', msg[:60]))
# bisect the ceiling at the top of the pressure interval
lo, hi = 200.0, 1.0e6
for _ in range(60):
    mid = 0.5*(lo + hi)
    _, ok, _ = bl.residual(c, P_i, ind, P_deep, mid, np.log10(P_deep))
    if ok:
        lo = mid
    else:
        hi = mid
print('ceiling at P_bg = P_surf : T_surf = %.6f K' % lo)
lo2, hi2 = 200.0, 1.0e6
xlow = np.log10(P_deep) - 12.0
for _ in range(60):
    mid = 0.5*(lo2 + hi2)
    _, ok, _ = bl.residual(c, P_i, ind, P_deep, mid, xlow)
    if ok:
        lo2 = mid
    else:
        hi2 = mid
print('ceiling at P_bg = P_surf*1e-12 : T_surf = %.6f K' % lo2)
