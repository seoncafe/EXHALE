#!/usr/bin/env python3
"""Move the He/H crossings by the counterfactual equivalent-width change.

d ln(He/H) = -d ln(EW)/s.  s = d ln EW / d ln(He/H), measured on the stored
scan runs by slopes.py, except for the flux-closure crossing where the value
of kzz_decision.md section 8.1 (0.382 over the 11.11-12.01 interval) is used.
dEW/EW from ew_table.py.
"""
import numpy as np
cases = [
 ('heh0p55                  well mixed, diffusion off', 0.55,  0.709,
  {'lo': +0.00863, 'blendA': -0.00281, 'blendB': -0.01407}),
 ('heh2p13_diff_kzz1e9      K_zz = 1e9 (adopted)',      2.09,  0.921,
  {'lo': +0.01008, 'blendA': -0.00297, 'blendB': -0.01202}),
 ('flux_closure/heh10p3/k01 flux closure, 10.3 arm',    11.73, 0.382,
  {'lo': +0.01471, 'blendA': -0.00227, 'blendB': -0.01070}),
 ('flux_closure/heh11p1/k01 flux closure, 11.11 arm',   11.73, 0.382,
  {'lo': +0.01387, 'blendA': -0.00219, 'blendB': -0.01038}),
]
print('%-52s %9s %6s  %-8s %9s %10s'
      % ('case', 'crossing', 'slope', 'variant', 'dEW/EW', 'He/H'))
for name, x0, s, d in cases:
    xs = {}
    for v in ('lo', 'blendA', 'blendB'):
        xs[v] = x0*np.exp(-d[v]/s)
        print('%-52s %9s %6s  %-8s %+8.2f%% %10.3f'
              % (name if v == 'lo' else '', ('%.2f' % x0) if v == 'lo' else '',
                 ('%.3f' % s) if v == 'lo' else '', v, 100*d[v], xs[v]))
    lo, hi = min(xs.values()), max(xs.values())
    print('%-52s %9s %6s  %-8s %9s %10s  span %.3f - %.3f  (%+.3f / %+.3f)\n'
          % ('', '', '', '', '', '', lo, hi, lo-x0, hi-x0))
