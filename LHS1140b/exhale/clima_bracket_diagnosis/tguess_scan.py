"""Does the deep-temperature guess move the failure?

`--climate-t-deep-guess` is the only user-facing knob that changes where the
outer MINPACK solve starts, so it decides which trial states it visits.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

heh = float(sys.argv[1])
guesses = [float(a) for a in sys.argv[2:]] or [300., 350., 380., 400., 420.,
                                               450., 500., 600.]
wd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'work',
                  'heh%s' % str(heh).replace('.', 'p'))
for g in guesses:
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    c.T_trop = bl.T_TROP_GUESS
    try:
        T = c.surface_temperature_bg_gas(P_i, P_deep, bg, T_guess=g)
        print('  T_guess=%7.1f  OK   T_deep=%10.5f  T_trop=%9.5f'
              % (g, T, c.T_trop))
    except Exception as e:
        print('  T_guess=%7.1f  FAIL %s' % (g, str(e)[:80]))
    sys.stdout.flush()
