#!/usr/bin/env python3
"""How close are the accepted arms to the clamp that deletes the metals?

`project_elements` empties a cell of its trace metals when the helium mass
fraction of that cell reaches the clamp `Xhe = 1` of
`binary_element_diffusion.f90:439`.  This re-solves a few arms that carry no
dropout, from their own converged state and with their own handoff, with
`EXHALE_DIFFUSION_CHECK=1`, and reports how near the composition operator
comes to that clamp.  The re-solve writes into this directory; the arms
themselves are not touched.

Usage:  python3 clamp_margin.py
"""
import os
import sys
import shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '../../..'))
sys.path.insert(0, os.path.join(EX, 'src/utils'))
import element_flux_closure as efc                          # noqa: E402

ARMS = [
    # tag, arm iteration directory (handoff and seed both taken from it)
    ('kzz1e9_heh11p1', 'flux_closure/heh11p1/k01'),
    ('kzz1e8_heh14p643', 'kzz_profile_scan/kzz1e8/heh14p643/k01'),
    ('kzz1e7_heh3p8', 'kzz_profile_scan/kzz1e7/heh3p8/k02'),
]


def x_range(runlog):
    mn, mx = [], []
    for ln in open(runlog, errors='replace'):
        if 'diffusion) X min' not in ln:
            continue
        f = ln.split()
        mn.append(float(f[f.index('min') + 1]))
        mx.append(float(f[f.index('max') + 1]))
    return np.array(mn), np.array(mx)


def main():
    root = os.path.join(HERE, 'clamp_margin')
    os.makedirs(root, exist_ok=True)
    print('%-20s %8s %14s %12s %10s' % ('arm', 'steps', 'X max', '1 - X max',
                                        'X = 1 hits'))
    for tag, rel in ARMS:
        arm = os.path.join(EX, 'LHS1140b/exhale', rel)
        d = os.path.join(root, tag)
        os.makedirs(os.path.join(d, 'output'), exist_ok=True)
        cfg = efc.read_configuration(
            os.path.join(os.path.dirname(arm), 'closure.json'))
        cfg['exhale_env'] = dict(cfg['exhale_env'])
        cfg['exhale_env']['EXHALE_DIFFUSION_CHECK'] = '1'
        for f in ('lower_atmosphere_profile.dat', 'base.inp'):
            shutil.copyfile(os.path.join(arm, f), os.path.join(d, f))
        try:
            outd, info, mdot = efc.solve_escape_wind(
                cfg, d, os.path.join(arm, 'output'), lambda *a: None)
        except Exception as exc:                            # noqa: BLE001
            print('%-20s solve stopped: %s' % (tag, exc))
            outd = os.path.join(d, 'output')
        mn, mx = x_range(os.path.join(d, 'run.log'))
        a = np.loadtxt(os.path.join(outd, 'Ion_species.txt'))
        nz = int(((a[:, 7] + a[:, 8] + a[:, 9]) <= 0.0).sum())
        if len(mx) == 0:
            print('%-20s %8s %14s %12s %10s' % (tag, 0, '-', '-', '-'))
            continue
        print('%-20s %8d %14.12f %12.3e %10d   zero-C cells after: %d'
              % (tag, len(mx), mx.max(), 1.0 - mx.max(),
                 int((mx >= 1.0).sum()), nz))


if __name__ == '__main__':
    main()
