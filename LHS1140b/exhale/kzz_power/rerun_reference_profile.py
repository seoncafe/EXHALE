#!/usr/bin/env python3
"""Re-solve the wind on a profile that is already in the tree.

Used to put the constant-K_zz reference arm and the K_zz(p) arms through the
same binary, so that the difference between them is the eddy coefficient and
not a code change made in between.  The profile and its paired base.inp are
copied verbatim; nothing about the lower atmosphere is recomputed.
"""
import os, shutil, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../../../src/utils'))
import element_flux_closure as efc                          # noqa: E402

srck, dstk, seed, cfgpath = sys.argv[1:5]
cfg = efc.read_configuration(cfgpath)
os.makedirs(os.path.join(dstk, 'output'), exist_ok=True)
for f in (cfg['profile_name'], 'base.inp'):
    shutil.copyfile(os.path.join(srck, f), os.path.join(dstk, f))
_, info, mdot = efc.solve_escape_wind(cfg, dstk, seed, print)
print('info = %s, log10 Mdot = %.3f' % (info, mdot))
sys.exit(0 if info == 0 else 1)
