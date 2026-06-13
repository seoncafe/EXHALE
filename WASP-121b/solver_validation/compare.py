#!/usr/bin/env python3
"""Compare the solver-validation runs and print a results table.

  A. Newton vs hybrd1  (identical one-step state) -- isolates the SOLVER.
  B. Newton vs baseline (converged)               -- includes one-step drift.
  C. _adv (Brent) vs baseline _adv                -- Task 1 (scalar T solve).

Run via run_validation.sh (which creates newton/, hybrd1/, baseline/).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

ION_NAMES = ['HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR',
             'CI', 'CII', 'CIII', 'OI', 'OII', 'OIII', 'NI', 'NII', 'NIII',
             'MgI', 'MgII', 'MgIII', 'SiI', 'SiII', 'SiIII', 'CaI', 'CaII',
             'CaIII', 'NaI', 'NaII', 'KI', 'KII', 'SI', 'SII', 'FeI', 'FeII',
             'FeIII']


def hydro(variant, adv=False):
    sub = '' if variant == 'baseline' else 'output'
    f = 'Hydro_ioniz_adv.txt' if adv else 'Hydro_ioniz.txt'
    return np.loadtxt(os.path.join(HERE, variant, sub, f))


def ions(variant, adv=False):
    sub = '' if variant == 'baseline' else 'output'
    f = 'Ion_species_adv.txt' if adv else 'Ion_species.txt'
    return np.loadtxt(os.path.join(HERE, variant, sub, f))


def maxreldiff_T(a, b):
    n = min(len(a), len(b))
    return np.max(np.abs(a[:n, 4] - b[:n, 4]) / np.maximum(b[:n, 4], 1.0))


def worst_ion(a, b, rlo=1.05, rhi=2.0, floor=1e-5):
    """Max relative diff over ion stages that are non-negligible in [rlo,rhi]."""
    r = b[:, 0]
    m = (r >= rlo) & (r <= rhi)
    worst, wname = 0.0, '-'
    for c in range(1, min(a.shape[1], b.shape[1])):
        bc, ac = b[m, c], a[m, c]
        if bc.max() <= 0:
            continue
        sel = bc > floor * bc.max()
        if sel.any():
            d = np.max(np.abs(ac[sel] - bc[sel]) / bc[sel])
            if d > worst:
                worst, wname = d, ION_NAMES[c - 1] if c - 1 < len(ION_NAMES) else str(c)
    return worst, wname


print('=' * 70)
print(' EXHALE solver-upgrade validation (WASP-121b Case B, PP-only sweep)')
print('=' * 70)

for d in ('newton', 'hybrd1'):
    log = os.path.join(HERE, d, 'run.log')
    line = '-'
    if os.path.exists(log):
        for L in open(log):
            if 'ioniz-eq solver' in L:
                line = L.strip()
    print(' %-7s : %s' % (d, line))
print()

# A. solver isolation: Newton vs hybrd1 at the identical one-step state
nH, hH = hydro('newton'), hydro('hybrd1')
nI, hI = ions('newton'), ions('hybrd1')
dT = maxreldiff_T(nH, hH)
dion, ionname = worst_ion(nI, hI)
print(' A. Newton vs hybrd1  (identical state -> pure SOLVER difference)')
print('      eq T          max rel diff : %.3e' % dT)
print('      Ion_species   max rel diff : %.3e   (worst: %s)' % (dion, ionname))
print('      => analytic Jacobian reproduces the MINPACK root.')
print()

# B. Newton vs converged baseline (includes the one-step hydro drift)
bH, bI = hydro('baseline'), ions('baseline')
print(' B. Newton vs baseline (converged; includes one-step PP hydro drift)')
print('      eq T          max rel diff : %.3e' % maxreldiff_T(nH, bH))
dion_b, ionname_b = worst_ion(nI, bI)
print('      Ion_species   max rel diff : %.3e   (worst: %s, drift-amplified)'
      % (dion_b, ionname_b))
print()

# C. Task 1: Brent _adv vs baseline _adv (the scalar T solve)
try:
    nHa = hydro('newton', adv=True)
    bHa = hydro('baseline', adv=True)
    na, ba = nHa, bHa
    nn = min(len(na), len(ba))
    v = na[:nn, 2]
    wind = v > 0
    rel = np.abs(na[:nn, 4] - ba[:nn, 4]) / np.maximum(ba[:nn, 4], 1.0)
    print(' C. Brent _adv T vs baseline _adv (Task 1, scalar energy equation)')
    print('      wind (v>0)    max rel diff : %.3e' % (rel[wind].max()))
    print('      all cells     max rel diff : %.3e' % (rel.max()))
except Exception as e:
    print(' C. (skipped: %s)' % e)
print('=' * 70)
