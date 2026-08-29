"""The deep boundary the climate solve returns, on the two installed builds.

Same climate object as the sweeps.  Run once under
`env/photochem/bin/python` (photochem 0.9.0 with the epsfcn repair) and once
under `/home/kiseon/.conda/envs/photochem_090_fix/bin/python` (the same
photochem 0.9.0 source without it), to measure how far the repair moves a
composition that both builds already solve.

usage: <python> deep_boundary_two_builds.py <He/H> ...
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'clima_bracket_diagnosis'))
import bracket_lib as bl  # noqa: E402

tag = os.path.basename(os.path.dirname(os.path.dirname(sys.executable)))
for heh in [float(a) for a in sys.argv[1:]]:
    wd = os.path.join(HERE, 'work_deep_' + tag, 'heh%s' % str(heh).replace('.', 'p'))
    c, P_i, bg, ind, P_deep, mixing = bl.build(heh, wd)
    T = c.surface_temperature_bg_gas(P_i, P_deep, bg, T_guess=bl.T_DEEP_GUESS)
    print('%-18s %-6s T_deep=%.10f  T_trop=%.10f' % (tag, heh, T, c.T_trop))
    sys.stdout.flush()
