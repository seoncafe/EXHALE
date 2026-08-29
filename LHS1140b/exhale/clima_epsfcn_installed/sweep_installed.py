"""Does the installed environment solve the He/H grid of Update_EXHALE sec. 91?

Calls `AdiabatClimate.surface_temperature_bg_gas` on the build that is
installed in `EXHALE_v1.00/env/photochem`, over the two composition grids
sec. 91 scans: He/H 8.00-11.00 in steps of 0.05 (61 points) and
9.30-9.60 in steps of 0.01 (31 points).  The climate object is the one the
LHS 1140 b closure ladder builds; the setup is imported unchanged from
`../clima_bracket_diagnosis/bracket_lib.py`, so the only difference from the
sec. 91 sweeps is the installed library.

usage: python sweep_installed.py <first> <last> <step>
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'clima_bracket_diagnosis'))
import bracket_lib as bl  # noqa: E402

first, last, step = (float(a) for a in sys.argv[1:4])
n = int(round((last - first)/step)) + 1
tag = sys.argv[4] if len(sys.argv) > 4 else 'grid'

for i in range(n):
    heh = round(first + i*step, 6)
    wd = os.path.join(HERE, 'work_' + tag, 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    t = time.time()
    try:
        T = c.surface_temperature_bg_gas(P_i, P_deep, bg, T_guess=bl.T_DEEP_GUESS)
        print('%-8s OK   T_deep=%9.5f  T_trop=%8.4f  P_bg=%.6e  bg=%s  %.1fs'
              % (heh, T, c.T_trop, c.f_i_surf[ind]*c.P_surf, bg, time.time() - t))
    except Exception as e:
        print('%-8s FAIL %s  (%.1fs)' % (heh, e, time.time() - t))
    sys.stdout.flush()
