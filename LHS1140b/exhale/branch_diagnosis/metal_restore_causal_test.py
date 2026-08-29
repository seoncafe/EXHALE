#!/usr/bin/env python3
"""Is the metal dropout what makes the hot state hot?

Takes the hot converged solution, refills the band in which its C, N and O
densities are exactly zero with the reservoir ratio of the handoff (all of it
in the neutral stage, which is how `load_IC` rebuilds an element the restart
file does not carry), and re-solves the wind from that seed with the SAME
handoff the hot arm used.  Nothing else is changed.

If the state falls to the cool branch, the dropout is the cause of the second
state and not a symptom of it.  If it stays hot, it is not.

The restored seed is a legitimate restart file: `load_IC` reconstructs the
mass density from the loaded species with the run's own mass policy, so
adding the metal nuclei back adds their mass in the same way a cold start
does.

Usage:  python3 metal_restore_causal_test.py
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
HOT = os.path.join(SCAN, 'heh3p399/k03')
COOL = os.path.join(SCAN, 'heh3p4/k02')

# Ion_species columns: 0 r, 1 HI, 2 HII, 3 HeI, 4 HeII, 5 HeIII, 6 HeITR,
# then CI CII CIII OI OII OIII NI NII NIII ...
COL = {'C': 7, 'O': 10, 'N': 13}


def restore(src_ion, dst_ion, ratios):
    """Refill every cell whose element total is zero with ratio * n_H, all of
    it in the neutral stage."""
    head = [ln for ln in open(src_ion) if ln.startswith('#')]
    a = np.loadtxt(src_ion)
    nH = a[:, 1] + a[:, 2]
    filled = {}
    for el, c in COL.items():
        tot = a[:, c] + a[:, c + 1] + a[:, c + 2]
        m = tot <= 0.0
        filled[el] = int(m.sum())
        a[m, c] = ratios[el] * nH[m]
        a[m, c + 1] = 0.0
        a[m, c + 2] = 0.0
    with open(dst_ion, 'w') as fh:
        fh.writelines(head)
        for row in a:
            fh.write(' '.join('%23.16E' % v for v in row) + '\n')
    return filled


def main():
    d = os.path.join(HERE, 'metal_restore')
    out = os.path.join(d, 'output')
    os.makedirs(out, exist_ok=True)

    # the reservoir ratios the hot arm's own handoff states
    ratios = {}
    for ln in open(os.path.join(HOT, 'EXHALE_resolved.out')):
        for el in ('C', 'O', 'N'):
            if ln.startswith('abundance_%s ' % el):
                ratios[el] = float(ln.split()[1])
    print('handoff reservoir ratios: ' +
          ', '.join('%s/H %.4e' % (k, v) for k, v in sorted(ratios.items())))

    seed = os.path.join(d, 'seed')
    os.makedirs(seed, exist_ok=True)
    shutil.copyfile(os.path.join(HOT, 'output/Hydro_ioniz.txt'),
                    os.path.join(seed, 'Hydro_ioniz.txt'))
    filled = restore(os.path.join(HOT, 'output/Ion_species.txt'),
                     os.path.join(seed, 'Ion_species.txt'), ratios)
    print('cells refilled: ' +
          ', '.join('%s %d' % (k, v) for k, v in sorted(filled.items())))

    cfg = efc.read_configuration(os.path.join(os.path.dirname(HOT),
                                              'closure.json'))
    for f in ('lower_atmosphere_profile.dat', 'base.inp'):
        shutil.copyfile(os.path.join(HOT, f), os.path.join(d, f))
    outd, info, mdot = efc.solve_escape_wind(cfg, d, seed, print)
    print('info = %s   log10 Mdot = %.3f' % (info, mdot))

    a = np.loadtxt(os.path.join(outd, 'Ion_species.txt'))
    r, C = a[:, 0], a[:, 7] + a[:, 8] + a[:, 9]
    z = np.where(C <= 0.0)[0]
    print('carbon zero cells after the re-solve: %d %s'
          % (len(z), '' if not len(z) else '(%.4f - %.4f R_p)'
             % (r[z[0]], r[z[-1]])))

    h = np.loadtxt(os.path.join(outd, 'Hydro_ioniz.txt'))
    print('T max = %.1f K at r = %.3f R_p' % (h[:, 4].max(),
                                              h[np.argmax(h[:, 4]), 0]))
    for tag, p in (('cool reference', COOL), ('hot reference', HOT)):
        g = np.loadtxt(os.path.join(p, 'output/Hydro_ioniz.txt'))
        print('  %-16s T max = %.1f K at r = %.3f R_p'
              % (tag, g[:, 4].max(), g[np.argmax(g[:, 4]), 0]))
    print('reference log10 Mdot: cool 7.540, hot 7.870')


if __name__ == '__main__':
    main()
