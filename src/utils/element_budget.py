#!/usr/bin/env python3
"""Element budget of a finished EXHALE run.

Every element the code carries is a conserved reservoir: the ionization
solver and the molecular network move nuclei between stages and molecules,
they never create or destroy them. So in every cell

    n_El(cell) / n_H(cell) = (El/H)_resolved,

where n_El sums the ion stages of that element, n_H counts hydrogen NUCLEI
(atomic, molecular and HeH+ bound), and (El/H)_resolved is the abundance the
run actually used -- from metals.inp, or from the "<El>_H_base" key of a
base.inp handoff, or from the X_<El> column of a lower-atmosphere profile.
Hydrogen itself closes against the mass density:
rho [m_H/cm^3] = mass_per_H * n_H.

Helium is the exception, and only when the element-diffusion operator is on.
`He_diffusion` solves a helium continuity equation of its own, so He/H is a
SOLVED PROFILE and the resolved `HeH_number_ratio` is the reservoir the base
is held at, not a column invariant: on LHS 1140 b it falls from 2.09 at the
base to 0.167 at 30 R_p, which is the separation the operator exists to
produce. Testing it cell by cell would call the physics a failure. With
`he_diffusion T` in EXHALE_resolved.out the He row is therefore checked at
the BASE cell, where it is a boundary condition, and the column-wide spread
is printed beside it as the measured separation. Hydrogen follows: its mass
closure runs through `mass_per_H`, which is built from the base He/H, so it
is a base-cell statement under diffusion for the same reason. The trace
elements are unaffected -- nothing diffuses them -- and stay column-wide.

This is the gate of P2 in docs/oxygen_chemistry_new_plan.md: it says that a
handoff which changes an elemental reservoir changes the whole column and
nothing else.

Usage:
    python3 src/utils/element_budget.py [run_dir] [--tol 1e-8] [--adv]

Reads <run_dir>/EXHALE_resolved.out (abundances, mass_per_H) and
<run_dir>/output/Ion_species.txt (+ Hydro_ioniz.txt for the H closure);
--adv checks the advection-corrected profiles instead. Exit status is 1 if
any element misses its reservoir by more than --tol (relative).
"""
import argparse
import os
import sys

import numpy as np

# stages that carry one nucleus of the element, by element symbol; built from
# the Ion_species.txt schema header at run time (the header is the authority)
MOLECULES = {           # species: (H nuclei, He nuclei)
    'H2':   (2, 0),
    'H2p':  (2, 0),
    'H3p':  (3, 0),
    'HeHp': (1, 1),
}
ATOMIC_H  = ('HI', 'HII')
# HeI is the TOTAL He I density, the He 2^3S metastable level included, so
# HeITR is deliberately not a member: it is an excited level of He I, not an
# independent species (bsp_is_excited_level in species_table.f90).
ATOMIC_HE = ('HeI', 'HeII', 'HeIII')
ROMAN = ('I', 'II', 'III', 'IV')


def read_resolved(path):
    """Abundances and EOS factors of the run (EXHALE_resolved.out)."""
    val = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            w = line.split()
            if len(w) >= 2:
                val[w[0]] = w[1]
    ab = {k[len('abundance_'):]: float(v)
          for k, v in val.items() if k.startswith('abundance_')}
    return val, ab


def read_columns(path):
    """Column labels from the '# columns' header, and the data block."""
    labels = None
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if 'columns' in line:
                    labels = line.split('columns', 1)[1].split()
                continue
            break
    if labels is None:
        raise SystemExit('%s: no "# columns" header' % path)
    data = np.loadtxt(path, comments='#')
    if data.shape[1] != len(labels):
        raise SystemExit('%s: %d columns, %d labels'
                         % (path, data.shape[1], len(labels)))
    return {lab: data[:, i] for i, lab in enumerate(labels)}, data[:, 0]


def element_stages(labels, sym):
    """Ion-stage column names present for element `sym` (CI, CII, CIII...)."""
    return [sym + r for r in ROMAN if sym + r in labels]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir', nargs='?', default='.')
    ap.add_argument('--tol', type=float, default=1.0e-8,
                    help='max allowed relative departure from the reservoir')
    ap.add_argument('--adv', action='store_true',
                    help='check the *_adv.txt advection-corrected profiles')
    a = ap.parse_args()

    tag = '_adv' if a.adv else ''
    res_file = os.path.join(a.run_dir, 'EXHALE_resolved.out')
    ion_file = os.path.join(a.run_dir, 'output', 'Ion_species%s.txt' % tag)
    hyd_file = os.path.join(a.run_dir, 'output', 'Hydro_ioniz%s.txt' % tag)
    for f in (res_file, ion_file, hyd_file):
        if not os.path.exists(f):
            raise SystemExit('missing %s' % f)

    val, ab = read_resolved(res_file)
    heh = float(val['HeH_number_ratio'])
    diffusing = val.get('he_diffusion', 'F').upper().startswith('T')
    mass_per_H = float(val['mass_per_H_amu'])
    col, r = read_columns(ion_file)
    hyd, _ = read_columns(hyd_file)
    rho = hyd['rho[mH/cm3]']

    # hydrogen nuclei: atomic + molecular carriers (the convention of
    # write_output.f90's Lyman-Werner diagnostic)
    nH = sum(col[s] for s in ATOMIC_H if s in col)
    nHe = sum(col[s] for s in ATOMIC_HE if s in col)
    for m, (h, he) in MOLECULES.items():
        if m in col:
            nH = nH + h * col[m]
            nHe = nHe + he * col[m]

    print('# element budget: %s' % os.path.abspath(a.run_dir))
    print('# %d cells, profiles %s, tolerance %.1e (relative)'
          % (len(r), os.path.basename(ion_file), a.tol))
    if diffusing:
        print('# He_diffusion is on: the H and He rows are checked at the '
              'base cell,')
        print('#   where He/H is the boundary reservoir; the solved '
              'separation is reported below.')
    print('%-6s %14s %14s %12s %10s' %
          ('elem', 'reservoir', 'max |n_El/n_H', 'at r[Rp]', 'verdict'))
    print('%-6s %14s %14s' % ('', '(El/H)', '- res|/res'))

    fail = 0

    def worst(dev):
        """The cell the check is stated at: the base cell when the helium
        is diffusing, otherwise the worst cell of the column."""
        return 0 if diffusing else int(np.argmax(dev))

    # H closes against the mass density instead of against a ratio
    dev = np.abs(rho / mass_per_H - nH) / nH
    j = worst(dev)
    ok = dev[j] <= a.tol
    fail += 0 if ok else 1
    print('%-6s %14s %14.3e %12.5f %10s'
          % ('H', 'rho/mass_per_H', dev[j], r[j], 'ok' if ok else 'FAIL'))

    if heh > 0.0:
        ratio = nHe / nH
        dev = np.abs(ratio - heh) / heh
        j = worst(dev)
        ok = dev[j] <= a.tol
        fail += 0 if ok else 1
        print('%-6s %14.6e %14.3e %12.5f %10s'
              % ('He', heh, dev[j], r[j], 'ok' if ok else 'FAIL'))
        if diffusing:
            print('%-6s %14s %14s %12s %10s'
                  % ('', 'separation', 'He/H top', '%.5f' % r[-1],
                     '%.6e' % ratio[-1]))

    for sym, x in sorted(ab.items()):
        if x <= 0.0:
            continue
        stages = element_stages(col, sym)
        if not stages:
            print('%-6s %14.6e %14s %12s %10s'
                  % (sym, x, '-', '-', 'no columns'))
            continue
        nEl = sum(col[s] for s in stages)
        dev = np.abs(nEl / nH - x) / x
        j = int(np.argmax(dev))
        ok = dev[j] <= a.tol
        fail += 0 if ok else 1
        print('%-6s %14.6e %14.3e %12.5f %10s'
              % (sym, x, dev[j], r[j], 'ok' if ok else 'FAIL'))

    if fail:
        print('\n%d element(s) miss their reservoir by more than %.1e'
              % (fail, a.tol))
    else:
        print('\nall element budgets close within %.1e' % a.tol)
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
