#!/usr/bin/env python3
"""Elemental closure of the equilibrium initializer at every saved deep state.

For each stored state under `scan/`, this reads the deepest row of the
profile the adapter wrote (its pressure and temperature), rebuilds the
elemental vector that run asked for, and calls the same routine the adapter
calls, `gasgiants.composition_at_metallicity`, at that single level.  The
reported residual is the one that decides acceptance:

    max_El | molfracs_atoms(El)/requested(El) - 1 |

taken over the TOTAL elemental vector, so a condensed phase does not appear
as a failure.  `mass_tol` is set to the value the gas-giant constructor uses
by default, so the number measured here is the number a production run gets.

The same script runs against an unmodified Photochem: there the equilibrium
solver scales every elemental residual by the largest elemental abundance and
`composition_at_metallicity` does not verify closure, so the residual it
reports is the one the adapter's handoff check has to absorb.

usage: closure_saved_states.py [--mass-tol X] [scan_dir]
"""
import argparse
import glob
import os
import sys

import numpy as np

import photochem
from photochem import equilibrate
from photochem.extensions import gasgiants

# Lodders solar vector, the one photochem_to_lower_profile.py starts from.
SOLAR = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5,
             O=6.06178e-4, S=1.31826e-5)
BAR = 1.0e6      # dyn/cm^2


def deepest(path):
    cols = None
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                if line.startswith('# columns:'):
                    cols = line.split(':', 1)[1].split()
                continue
            row = [float(x) for x in line.split()]
            return cols, row
    raise RuntimeError('no data rows in %s' % path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('scan_dir', nargs='?',
                    default=os.path.join(os.path.dirname(
                        os.path.abspath(__file__)), 'scan'))
    ap.add_argument('--mass-tol', type=float, default=None,
                    help='override the equilibrium mass tolerance; default '
                         'is whatever the gas-giant constructor uses')
    args = ap.parse_args()

    default_tol = None
    try:
        import inspect
        sig = inspect.signature(gasgiants.GasGiantData.__init__)
        p = sig.parameters.get('equilibrium_mass_tol')
        if p is not None:
            default_tol = p.default
    except Exception:
        pass

    print('# photochem %s, equilibrate %s'
          % (photochem.__version__, equilibrate.__version__))
    print('# gas-giant default equilibrium mass_tol: %s'
          % ('absent (unmodified build)' if default_tol is None
             else '%.1e' % default_tol))
    print('# residual = max over requested elements of '
          '|molfracs_atoms/requested - 1|')
    print('%-9s %-11s %-9s %-6s %-12s %-4s %-12s %-12s'
          % ('He/H', 'P[dyn/cm2]', 'T[K]', 'conv', 'worst_rel', 'el',
             'N_rel', 'He_rel'))

    worst_all, worst_at = 0.0, None
    rows = []
    for d in sorted(glob.glob(os.path.join(args.scan_dir, 'heh*'))):
        prof = os.path.join(d, 'lower_atmosphere_profile.dat')
        thermo = os.path.join(d, 'photochem_work',
                              'zahnle_HHeNOC_thermo.yaml')
        if not (os.path.exists(prof) and os.path.exists(thermo)):
            continue
        heh = float(os.path.basename(d)[3:].replace('p', '.'))
        cols, row = deepest(prof)
        P = row[cols.index('p')]*BAR
        T = row[cols.index('T')]

        gas = equilibrate.ChemEquiAnalysis(thermo)
        tol = args.mass_tol if args.mass_tol is not None else default_tol
        if tol is not None:
            gas.mass_tol = tol
        atoms = list(gas.atoms_names)
        ab = dict(SOLAR)
        ab['He'] = heh
        tot = sum(ab[a] for a in atoms)
        want = np.array([ab[a]/tot for a in atoms])
        gas.molfracs_atoms_sun = want.copy()

        try:
            gasgiants.composition_at_metallicity(
                gas, np.array([T]), np.array([P]), 1.0, 1.0)
            conv = 'yes'
        except RuntimeError as exc:
            conv = 'REFUSED'
            print('%-9.4f %-11.4e %-9.3f %-6s  %s'
                  % (heh, P, T, conv, str(exc).splitlines()[0]))
            rows.append((heh, float('nan'), 'refused'))
            continue

        got = np.asarray(gas.molfracs_atoms, dtype=float)
        pos = want > 0.0
        rel = np.abs(got[pos]/want[pos] - 1.0)
        names = [a for a, k in zip(atoms, pos) if k]
        j = int(np.argmax(rel))
        rN = rel[names.index('N')] if 'N' in names else float('nan')
        rHe = rel[names.index('He')] if 'He' in names else float('nan')
        print('%-9.4f %-11.4e %-9.3f %-6s %-12.4e %-4s %-12.4e %-12.4e'
              % (heh, P, T, conv, rel[j], names[j], rN, rHe))
        rows.append((heh, float(rel[j]), names[j]))
        if rel[j] > worst_all:
            worst_all, worst_at = float(rel[j]), (heh, names[j])

    good = [r[1] for r in rows if r[1] == r[1]]
    print('#')
    print('# states measured: %d' % len(rows))
    if good:
        print('# largest relative residual: %.9e  (He/H = %.4g, %s)'
              % (worst_all, worst_at[0], worst_at[1]))
        print('# median  relative residual: %.9e' % float(np.median(good)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
