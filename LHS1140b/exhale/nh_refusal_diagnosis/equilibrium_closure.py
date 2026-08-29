#!/usr/bin/env python3
"""Elemental closure of Photochem's equilibrium solver, alone.

`gasgiants.composition_at_metallicity` builds the initial composition of the
column by calling `equilibrate.ChemEquiAnalysis.solve` level by level and --
its own comment -- does "not enforce convergence".  The probe
`probe_initial_vs_steady.py` showed that the deepest-level El/H departure the
adapter refuses on is already present in that composition and does not move
by one bit through the photochemical solve, so this is where it is made.

This script calls the same solver at the deepest level of the LHS 1140 b
climate column and reports, for a sweep of reservoir He/H:

  - whether `solve` reported convergence;
  - El/H of the returned gas composition against the El/H asked for;
  - the same for the solver's own `molfracs_atoms_gas`.

usage: equilibrium_closure.py <thermo.yaml> <P dyn/cm2> <T K> <heh> [heh ...]
"""
import sys

import numpy as np
from photochem import equilibrate

SOLAR = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5,
             O=6.06178e-4, S=1.31826e-5)


def main():
    thermo, P, T = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    gas = equilibrate.ChemEquiAnalysis(thermo)
    atoms = list(gas.atoms_names)
    names = list(gas.gas_names)
    # element counts of each gas species, from the solver's own matrix if it
    # exposes one, otherwise by formula (these are all plain formulas).
    import re
    rex = re.compile(r'([A-Z][a-z]?)(\d*)')

    def counts(sp):
        out, pos = {}, 0
        while pos < len(sp):
            m = rex.match(sp, pos)
            if m is None or m.start() != pos:
                return {}
            el = m.group(1)
            if el not in atoms:
                return {}
            out[el] = out.get(el, 0) + int(m.group(2) or 1)
            pos = m.end()
        return out

    cnt = {sp: counts(sp) for sp in names}
    bad = [sp for sp in names if not cnt[sp]]
    if bad:
        print('# unparsed species (skipped): %s' % ', '.join(bad))

    print('# P = %.6e dyn/cm2, T = %.4f K, thermo %s' % (P, T, thermo))
    print('%-8s %-5s %-11s %-11s %-11s %-11s %-11s'
          % ('He/H', 'conv', 'dev_He', 'dev_C', 'dev_N', 'dev_O', 'dev_N_atomsgas'))
    for spec in sys.argv[4:]:
        heh = float(spec)
        ab = dict(SOLAR)
        ab['He'] = heh
        mf = np.array([ab[a] for a in atoms])
        mf = mf/mf.sum()
        want = {a: ab[a]/ab['H'] for a in atoms if a != 'H'}
        for eps in [0.0, 1.0e-12, -1.0e-12, 1.0e-8, -1.0e-8]:
            conv = gas.solve(P, T + T*eps, molfracs_atoms=mf)
            if conv:
                break
        y = np.asarray(gas.molfracs_species_gas)
        nuc = {}
        for j, sp in enumerate(names):
            for el, k in cnt[sp].items():
                nuc[el] = nuc.get(el, 0.0) + k*y[j]
        got = {el: nuc[el]/nuc['H'] for el in nuc if el != 'H'}
        ag = np.asarray(gas.molfracs_atoms_gas)
        iH = atoms.index('H')
        gotA = {a: ag[i]/ag[iH] for i, a in enumerate(atoms) if a != 'H'}
        print('%-8.4f %-5s %+11.3e %+11.3e %+11.3e %+11.3e %+11.3e'
              % (heh, conv,
                 got['He']/want['He'] - 1.0, got['C']/want['C'] - 1.0,
                 got['N']/want['N'] - 1.0, got['O']/want['O'] - 1.0,
                 gotA['N']/want['N'] - 1.0))


if __name__ == '__main__':
    main()
