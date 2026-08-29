#!/usr/bin/env python3
"""Does the cool branch really end at He/H = 3.8 at K_zz = 1e7?

`kzz_profile_scan/results.txt` result 2 says the cool state "ceases to exist,
it does not merely fail to converge" above He/H = 3.8: the He/H = 4.2
continuation leaves it.  The census shows what the 4.2 arm actually did --
it acquired the band of exactly zero C/N/O.

Here the 4.2 step is taken again from the clean 3.8 state, the deleted C, N
and O are refilled at the handoff reservoir ratio, and the wind is re-solved.
If it lands on the cool branch, the branch does not end at 3.8.

Usage:  python3 cool_branch_continuation.py
"""
import os
import sys
import shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '../../..'))
sys.path.insert(0, os.path.join(EX, 'src/utils'))
sys.path.insert(0, HERE)
import element_flux_closure as efc                          # noqa: E402
from metal_restore_causal_test import restore               # noqa: E402

SCAN = os.path.join(EX, 'LHS1140b/exhale/kzz_profile_scan/kzz1e7')
ARM = os.path.join(SCAN, 'heh4p2/k00')          # the handoff of the 4.2 arm
FIRST = os.path.join(HERE, 'metal_dropout_runs/heh4p2_from_heh3p8')


def zero_band(outdir):
    a = np.loadtxt(os.path.join(outdir, 'Ion_species.txt'))
    r, C = a[:, 0], a[:, 7] + a[:, 8] + a[:, 9]
    z = np.where(C <= 0.0)[0]
    return len(z), (None if not len(z) else (r[z[0]], r[z[-1]]))


def main():
    ratios = {}
    for ln in open(os.path.join(ARM, 'EXHALE_resolved.out')):
        for el in ('C', 'O', 'N'):
            if ln.startswith('abundance_%s ' % el):
                ratios[el] = float(ln.split()[1])
    print('handoff reservoir ratios: ' +
          ', '.join('%s/H %.4e' % (k, v) for k, v in sorted(ratios.items())))

    src = os.path.join(FIRST, 'output')
    n0, b0 = zero_band(src)
    print('the re-taken 3.8 -> 4.2 step left %d zero-carbon cells %s'
          % (n0, '' if b0 is None else '(%.4f - %.4f R_p)' % b0))

    d = os.path.join(HERE, 'cool_continuation')
    seed = os.path.join(d, 'seed')
    os.makedirs(seed, exist_ok=True)
    os.makedirs(os.path.join(d, 'output'), exist_ok=True)
    shutil.copyfile(os.path.join(src, 'Hydro_ioniz.txt'),
                    os.path.join(seed, 'Hydro_ioniz.txt'))
    filled = restore(os.path.join(src, 'Ion_species.txt'),
                     os.path.join(seed, 'Ion_species.txt'), ratios)
    print('cells refilled: ' +
          ', '.join('%s %d' % (k, v) for k, v in sorted(filled.items())))

    cfg = efc.read_configuration(os.path.join(os.path.dirname(ARM),
                                              'closure.json'))
    for f in ('lower_atmosphere_profile.dat', 'base.inp'):
        shutil.copyfile(os.path.join(ARM, f), os.path.join(d, f))
    outd, info, mdot = efc.solve_escape_wind(cfg, d, seed, print)
    n1, b1 = zero_band(outd)
    h = np.loadtxt(os.path.join(outd, 'Hydro_ioniz.txt'))
    print('info = %s   log10 Mdot = %.3f   zero-carbon cells %d %s'
          % (info, mdot, n1,
             '' if b1 is None else '(%.4f - %.4f R_p)' % b1))
    print('T max = %.1f K at r = %.3f R_p' % (h[:, 4].max(),
                                              h[np.argmax(h[:, 4]), 0]))
    print()
    print('for comparison, the cool branch at this K_zz:')
    print('  He/H 3.0  log10 Mdot 7.560   3.4  7.540   3.8  7.520')
    print('the arm as the scan recorded it: He/H 4.2, log10 Mdot 7.860,'
          ' window spread 0.093, EW unresolved')


if __name__ == '__main__':
    main()
