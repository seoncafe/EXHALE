"""Does the installed build refuse He/H = 9.4 and 9.5 and accept 9.3/9.45?"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bracket_lib as bl

for heh in [float(a) for a in sys.argv[1:]]:
    wd = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'work', 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    t = time.time()
    try:
        T = c.surface_temperature_bg_gas(P_i, P_deep, bg, T_guess=bl.T_DEEP_GUESS)
        print('%-8s OK   T_deep=%9.4f  T_trop=%8.4f  P_bg=%.6e  bg=%s  %.1fs'
              % (heh, T, c.T_trop, c.f_i_surf[ind]*c.P_surf, bg, time.time()-t))
    except Exception as e:
        print('%-8s FAIL %s  (%.1fs)' % (heh, e, time.time()-t))
    sys.stdout.flush()
