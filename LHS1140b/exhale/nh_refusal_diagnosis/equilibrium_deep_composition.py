#!/usr/bin/env python3
"""Deep-level equilibrium composition against the reservoir He/H.

Companion to equilibrium_closure.py: the same solver, the same level, but
reporting the mixing ratios themselves.  It answers whether the reservoir
band the adapter refuses is a different chemical branch (a jump in the
composition) or the same branch with a larger solver residual (a jump in the
residual alone).
"""
import sys
import numpy as np
from photochem import equilibrate

SOLAR = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5,
             O=6.06178e-4, S=1.31826e-5)
SHOW = ('H2', 'H2O', 'CH4', 'CO', 'NH3', 'N2', 'HCN', 'He')

thermo, P, T = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
gas = equilibrate.ChemEquiAnalysis(thermo)
atoms, names = list(gas.atoms_names), list(gas.gas_names)
print('# P = %.4e dyn/cm2, T = %.2f K' % (P, T))
print('%-8s %s' % ('He/H', ' '.join('%-12s' % s for s in SHOW)))
for spec in sys.argv[4:]:
    heh = float(spec)
    ab = dict(SOLAR); ab['He'] = heh
    mf = np.array([ab[a] for a in atoms]); mf = mf/mf.sum()
    for eps in [0.0, 1.0e-12, -1.0e-12, 1.0e-8, -1.0e-8]:
        if gas.solve(P, T + T*eps, molfracs_atoms=mf):
            break
    y = np.asarray(gas.molfracs_species_gas)
    print('%-8.3f %s' % (heh, ' '.join(
        '%-12.5e' % (y[names.index(s)] if s in names else np.nan)
        for s in SHOW)))
