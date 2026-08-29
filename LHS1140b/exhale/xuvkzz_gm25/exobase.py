#!/usr/bin/env python3
"""Exobase and sonic point of the re-solved ladder arms.

The argument that discards the jump-seeded 12.01 solution rests on where the
two solutions differ lying outside the exobase, so that the fluid equations
do not hold there.  This re-evaluates that on the current coefficient with
src/utils/collisional_validity.py, which reads a run directory and changes
nothing in it.
"""
import os, sys
import numpy as np
EX = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
sys.path.insert(0, os.path.join(EX, 'src', 'utils'))
sys.path.insert(0, EX)
import collisional_validity as CV

os.chdir(os.path.dirname(os.path.abspath(__file__)))
CASES = [('2.09', 'L2p09'), ('3.0', 'L3p0'), ('5.0', 'L5p0'),
         ('8.0', 'L8p0'), ('9.7', 'L9p7'), ('10.3', 'L10p3'),
         ('10.7', 'L10p7'), ('11.1', 'L11p1'),
         ('12.0 one step', 'L12p0ctl'), ('12.0 one jump', 'L12p0jump')]
print('%-16s %12s %12s' % ('case', 'r_exobase', 'r_sonic'))
lo, hi = [], []
for lab, d in CASES:
    try:
        res = CV.collisional_diagnosis(d, adv=True)
    except Exception as exc:
        print('%-16s  FAILED  %s' % (lab, exc))
        continue
    re_, rs = res['r_exobase'], res['r_sonic']
    print('%-16s %12s %12s'
          % (lab, 'above domain' if re_ is None else '%.2f' % re_,
             'none' if rs is None else '%.2f' % rs))
    if re_ is not None:
        lo.append(re_)
if lo:
    print('exobase range over the ladder: %.1f - %.1f R_p' % (min(lo), max(lo)))
