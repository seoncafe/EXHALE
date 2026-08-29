#!/usr/bin/env python3
"""Reproduce the creation of the C/N/O dropout band.

The census (`T4_dropout_census.txt`) shows the band is created inside a wind
solve, not inherited: `kzz1e7/heh4p2/k00` and `kzz1e7/heh4p500/k00` carry it
while the arm they were seeded from, `kzz1e7/heh3p8`, does not.  Here the
same pair of steps is re-run with `EXHALE_DIFFUSION_CHECK=1`, which makes
`binary_element_diffusion` report the range of the helium mass fraction each
step, so that the state of the composition operator when the metals disappear
is on record.

Usage:  python3 metal_dropout_runs_reproduce.py
"""
import os
import sys
import shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '../../..'))
sys.path.insert(0, os.path.join(EX, 'src/utils'))
import element_flux_closure as efc                          # noqa: E402

SCAN = os.path.join(EX, 'LHS1140b/exhale/kzz_profile_scan/kzz1e7')


def carbon_zero_band(outdir):
    a = np.loadtxt(os.path.join(outdir, 'Ion_species.txt'))
    r, C = a[:, 0], a[:, 7] + a[:, 8] + a[:, 9]
    z = np.where(C <= 0.0)[0]
    if len(z) == 0:
        return 0, None
    return len(z), (r[z[0]], r[z[-1]])


def main():
    root = os.path.join(HERE, 'metal_dropout_runs')
    os.makedirs(root, exist_ok=True)
    cases = [
        # tag, profile source (the arm whose handoff is imposed), seed arm
        ('heh4p2_from_heh3p8', os.path.join(SCAN, 'heh4p2/k00'),
         os.path.join(SCAN, 'heh3p8/k02')),
        ('heh4p500_from_heh3p8', os.path.join(SCAN, 'heh4p500/k00'),
         os.path.join(SCAN, 'heh3p8/k02')),
    ]
    for tag, pdir, sdir in cases:
        n0, b0 = carbon_zero_band(os.path.join(sdir, 'output'))
        print('=== %s : seed %s has %d zero cells' % (tag, sdir, n0))
        d = os.path.join(root, tag)
        os.makedirs(d, exist_ok=True)
        cfg = efc.read_configuration(
            os.path.join(os.path.dirname(pdir), 'closure.json'))
        cfg['exhale_env'] = dict(cfg['exhale_env'])
        cfg['exhale_env']['EXHALE_DIFFUSION_CHECK'] = '1'
        for f in ('lower_atmosphere_profile.dat', 'base.inp'):
            shutil.copyfile(os.path.join(pdir, f), os.path.join(d, f))
        try:
            outd, info, mdot = efc.solve_escape_wind(
                cfg, d, os.path.join(sdir, 'output'), print)
        except Exception as exc:                            # noqa: BLE001
            print('   solve stopped: %s' % exc)
            outd = os.path.join(d, 'output')
            info, mdot = None, float('nan')
        n1, b1 = carbon_zero_band(outd)
        print('   info=%s  log10 Mdot=%s' % (info, mdot))
        print('   carbon zero cells: seed %d -> solved %d  %s'
              % (n0, n1, '' if b1 is None else '(%.4f - %.4f R_p)' % b1))
        # what the diffusion operator reported
        lines = [ln.rstrip() for ln in
                 open(os.path.join(d, 'run.log'), errors='replace')
                 if 'diffusion) X min' in ln]
        if lines:
            def field(ln, key):
                f = ln.split()
                return float(f[f.index(key) + 1])
            xmax = [field(ln, 'max') for ln in lines]
            xmin = [field(ln, 'min') for ln in lines]
            print('   diffusion steps reported: %d,  X in [%.6f, %.9f],'
                  ' steps with X = 1 exactly: %d'
                  % (len(lines), min(xmin), max(xmax),
                     sum(1 for x in xmax if x >= 1.0)))
            for ln in lines[:3] + ['   ...'] + lines[-3:]:
                print('     %s' % ln)
        else:
            print('   the diffusion operator wrote no step report')
        print()


if __name__ == '__main__':
    main()
