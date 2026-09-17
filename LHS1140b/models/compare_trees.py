#!/usr/bin/env python3
"""Put two or more generations of the LHS 1140 b catalog side by side.

The catalog is re-solved whenever a fix reaches the equations, and each time
the previous answers are moved aside by `models/archive_results.sh` into a
tree of the same shape.  This reads any number of those trees and prints one
row per case, so what a fix moved is read off the row and not looked up in a
memo.

    python3 models/compare_trees.py [--markdown] [--only <substring>]
                                    [<label>=<tree> ...]

With no tree given it compares the three generations of the 2026-09-15/16
work: `models_20260914_preL21/` (the results on the base boundary of item L21
as it stood), `models_20260915_db87/` (the re-run on the corrected boundary,
binary db87b88d1ce5) and `models/` (the current tree).  Name them explicitly
to compare others:

    python3 models/compare_trees.py preL21=../models_20260914_preL21 \\
        db87=../models_20260915_db87 final=.

Every number is read the way `models/status.py` reads it, so a value here and
a value in `MODELS.md` sections 7 and 8 are the same number, with two added:

  log10 Mdot   formed from the state itself, 4 pi r^2 rho v at one physical
               cell, because `pp.log` prints it to two decimals and that is
               2 per cent in Mdot -- coarser than the movement being measured
  T1, rho1     the first physical cell, which is the row the base boundary
               decides and the row the He I 10830 line is most sensitive to

The molecular groups are left out: they are solved on their own schedule and
no preserved tree holds them.
"""

import argparse
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LHS = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import status                                               # noqa: E402
from make_models import GROUPS, Group                       # noqa: E402

DEFAULT_TREES = [('preL21', os.path.join(LHS, 'models_20260914_preL21')),
                 ('db87', os.path.join(LHS, 'models_20260915_db87')),
                 ('final', HERE)]

# CODATA 2018 hydrogen-atom mass, the unit the `rho` column is written in.
M_H = 1.67353e-24

# The physical cell the mass flux is read at.  The wind is subsonic to the
# outer boundary, so the flux is the same at every physical cell to the
# accuracy of the solve; the same index is read in every tree.
FLUX_ROW = 501


def state_dir(case_dir):
    """Where the products of this case are: the case itself, or a rung's last
    iterate."""
    if os.path.isdir(os.path.join(case_dir, 'k00')):
        return status.last_iterate(case_dir)
    return case_dir


def base_state(case_dir):
    """(T of the inner ghost, T of cell 1, rho of cell 1)."""
    d = state_dir(case_dir)
    if d is None:
        return (None, None, None)
    path = os.path.join(d, 'output', 'Hydro_ioniz.txt')
    if not os.path.isfile(path):
        return (None, None, None)
    rows = np.loadtxt(path, comments='#')
    if rows.ndim != 2 or rows.shape[0] < 4:
        return (None, None, None)
    return (rows[1, 4], rows[2, 4], rows[2, 1])


def mass_loss_rate(case_dir, row=FLUX_ROW):
    """log10 of 4 pi r^2 rho v [g/s] at one physical cell of the state."""
    d = state_dir(case_dir)
    if d is None:
        return None
    path = os.path.join(d, 'output', 'Hydro_ioniz.txt')
    if not os.path.isfile(path):
        return None
    R0 = None
    with open(path, errors='replace') as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            if line.startswith('# grid'):
                f = line.split()
                R0 = float(f[f.index('R0[cm]') + 1])
    rows = np.loadtxt(path, comments='#')
    if R0 is None or rows.ndim != 2 or rows.shape[0] <= row:
        return None
    r = rows[row, 0]*R0
    mdot = 4.0*np.pi*r*r*rows[row, 1]*M_H*rows[row, 2]
    return np.log10(mdot) if mdot > 0.0 else None


def wall_clock(case_dir):
    path = os.path.join(case_dir, 'REPRODUCE.md')
    if not os.path.isfile(path):
        return None
    m = re.search(r'wall clock \*\*(\d+m\d+s)\*\*',
                  open(path, errors='replace').read())
    return m.group(1) if m else None


def outer_passes(case_dir):
    d = state_dir(case_dir)
    if d is None:
        return None
    log = os.path.join(d, 'run.log')
    if not os.path.isfile(log):
        return None
    n = re.findall(r'outer pass +(\d+)', open(log, errors='replace').read())
    return int(n[-1]) if n else None


def read_tree(root, only=None):
    """Every case of the group table as `status.Case` reads it under `root`."""
    status.HERE = root
    out = {}
    for name, ladder in GROUPS:
        if name.startswith('molecular_'):
            continue
        if only and only not in name:
            continue
        g = Group(name)
        for heh in ladder:
            key = '%s/HeH%s' % (name, heh)
            cdir = os.path.join(root, name, 'HeH%s' % heh)
            c = status.Case(g, heh)
            there = os.path.isdir(cdir)
            tg, t1, r1 = base_state(cdir) if there else (None, None, None)
            out[key] = dict(case=c, group=name, heh=heh, T_ghost=tg,
                            T_cell1=t1, rho_cell1=r1,
                            mdot=mass_loss_rate(cdir) if there else None,
                            wall=wall_clock(cdir) if there else None,
                            passes=outer_passes(cdir) if there else None)
    return out


def crossings(tree, root):
    status.HERE = root
    out = {}
    for name, ladder in GROUPS:
        if name.startswith('molecular_'):
            continue
        heh_num, ew_num = [], []
        for heh in ladder:
            e = tree.get('%s/HeH%s' % (name, heh))
            if e is None:
                continue
            c = e['case']
            heh_num.append(c.closure[1] if c.closure is not None
                           else float(heh))
            ew_num.append(c.ew)
        if heh_num:
            out[name] = status.equivalent_width_crossing(heh_num, ew_num)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('trees', nargs='*', metavar='label=tree')
    ap.add_argument('--markdown', action='store_true')
    ap.add_argument('--only', default=None)
    a = ap.parse_args()

    trees = []
    for spec in a.trees:
        label, _, path = spec.partition('=')
        if not path:
            ap.error('give each tree as <label>=<path>')
        trees.append((label, os.path.abspath(
            path if os.path.isabs(path) else os.path.join(HERE, path))))
    if not trees:
        trees = DEFAULT_TREES
    read = [(label, root, read_tree(root, a.only)) for label, root in trees]

    head = ['case', 'state (%s)' % read[-1][0]]
    for label, _, _ in read:
        head += ['log10 Mdot %s' % label, 'EW %s' % label,
                 'depth %s' % label, 'T1 %s' % label]
    head += ['d Mdot [%]', 'd EW [%]', 'passes', 'wall']
    rows = [head]
    for key in sorted(read[0][2]):
        entries = [t[key] for _, _, t in read]
        if all(e['case'].state == 'none' for e in entries):
            continue
        row = [key, entries[-1]['case'].state]
        for e in entries:
            c = e['case']
            row += [status.fmt(e['mdot'], '.4f'), status.fmt(c.ew, '.4f'),
                    status.fmt(c.depth, '.3f'),
                    status.fmt(e['T_cell1'], '.1f')]
        first, last = entries[0], entries[-1]
        row += [pct(last['mdot'], first['mdot'], log=True),
                pct(last['case'].ew, first['case'].ew),
                str(last['passes']), str(last['wall'])]
        rows.append(row)

    w = [max(len(r[i]) for r in rows) for i in range(len(head))]
    for i, r in enumerate(rows):
        cells = [' %-*s ' % (w[j], r[j]) for j in range(len(head))]
        print(('|' + '|'.join(cells) + '|') if a.markdown else ' '.join(cells))
        if a.markdown and i == 0:
            print('|' + '|'.join('-'*(w[j] + 2) for j in range(len(head)))
                  + '|')

    summarize(read)
    agreement(read)
    print('')
    print('crossings He/H (measured EW = %.3f %%A)' % status.EW_OBS)
    cx = [(label, crossings(t, root)) for label, root, t in read]
    for name, _ in GROUPS:
        if name.startswith('molecular_') or (a.only and a.only not in name):
            continue
        if not any(name in c for _, c in cx):
            continue
        print('  %-46s %s' % (name, '  ->  '.join(
            '%s %s' % (label, '%.4f' % c[name][0] if c.get(name) else
                       'no bracket') for label, c in cx)))


def pct(new, old, log=False):
    if new is None or old is None:
        return '--'
    if log:
        return '%+.2f' % (100.0*(10.0**(new - old) - 1.0))
    if old == 0.0:
        return '--'
    return '%+.2f' % (100.0*(new/old - 1.0))


def agreement(read):
    """How closely each tree reproduces the one before it, case by case.

    Where two generations solve THE SAME equations -- a binary that carries
    only changes the wind does not see -- the two states are the same fixed
    point and the departure has to sit inside the certification tolerances
    (1e-6 on a hydrodynamic row, 1e-5 on a gated species row). This block is
    the check: the median and the worst relative departure of each measured
    quantity, over the cases both trees certify, and how many are identical
    to the last bit.
    """
    for (la, _, a), (lb, _, b) in zip(read, read[1:]):
        rows = []
        for k in sorted(b):
            if a[k]['case'].state != 'info=0 certified' \
               or b[k]['case'].state != 'info=0 certified':
                continue
            rows.append((k, a[k], b[k]))
        if not rows:
            continue
        print('')
        print('%s against %s, over the %d cases both certify'
              % (lb, la, len(rows)))

        def report(what, value):
            v, worst, exact = [], (None, -1.0), 0
            for k, x, y in rows:
                d = value(x, y)
                if d is None:
                    continue
                v.append(d)
                if d > worst[1]:
                    worst = (k, d)
                if d == 0.0:
                    exact += 1
            if not v:
                return
            v.sort()
            digits = ('exactly' if v[-1] == 0.0
                      else 'to %.1f decimal digits' % -math.log10(v[-1]))
            print('  %-14s median %8.2e  worst %8.2e (%s)'
                  % (what, v[len(v)//2], v[-1], worst[0]))
            print('  %-14s identical to the last bit in %d of %d; agrees %s'
                  % ('', exact, len(v), digits))

        def rel(f):
            def g(x, y):
                p, q = f(x), f(y)
                if p is None or q is None or p == 0.0:
                    return None
                return abs(q/p - 1.0)
            return g
        report('log10 Mdot', lambda x, y: None
               if (x['mdot'] is None or y['mdot'] is None)
               else abs(10.0**(y['mdot'] - x['mdot']) - 1.0))
        report('red EW', rel(lambda e: e['case'].ew))
        report('red depth', rel(lambda e: e['case'].depth))
        report('T of cell 1', rel(lambda e: e['T_cell1']))
        report('rho of cell 1', rel(lambda e: e['rho_cell1']))


def summarize(read):
    """The movement of every later tree against the first."""
    first = read[0][2]
    for label, _, t in read[1:]:
        dm, de, dd, dt, dr, pas, walls = [], [], [], [], [], [], []
        for key, e in t.items():
            o = first[key]
            if e['mdot'] is not None and o['mdot'] is not None:
                dm.append(100.0*(10.0**(e['mdot'] - o['mdot']) - 1.0))
            if e['case'].ew and o['case'].ew:
                de.append(100.0*(e['case'].ew/o['case'].ew - 1.0))
            if e['case'].depth and o['case'].depth:
                dd.append(100.0*(e['case'].depth/o['case'].depth - 1.0))
            if e['T_cell1'] and o['T_cell1']:
                dt.append(100.0*(e['T_cell1']/o['T_cell1'] - 1.0))
            if e['rho_cell1'] and o['rho_cell1']:
                dr.append(e['rho_cell1']/o['rho_cell1'])
            if e['passes'] is not None:
                pas.append(e['passes'])
            if e['wall']:
                m_, s_ = e['wall'].rstrip('s').split('m')
                walls.append(int(m_) + int(s_)/60.0)
        print('')
        print('movement of %s against %s' % (label, read[0][0]))
        for vals, what, unit in ((dm, 'Mdot', ' %'), (de, 'red EW', ' %'),
                                 (dd, 'red depth', ' %'),
                                 (dt, 'T of cell 1', ' %')):
            if vals:
                print('  %-24s %3d cases, %+.2f to %+.2f%s (median %+.2f%s)'
                      % (what, len(vals), min(vals), max(vals), unit,
                         sorted(vals)[len(vals)//2], unit))
        if dr:
            print('  %-24s %3d cases, %.3f to %.3f (median %.3f)'
                  % ('rho of cell 1, ratio', len(dr), min(dr), max(dr),
                     sorted(dr)[len(dr)//2]))
        if pas:
            print('  %-24s %3d cases, %d to %d (median %d)'
                  % ('outer passes', len(pas), min(pas), max(pas),
                     sorted(pas)[len(pas)//2]))
        if walls:
            print('  %-24s %3d cases, %.1f to %.1f min, %.1f h in all'
                  % ('wall clock', len(walls), min(walls), max(walls),
                     sum(walls)/60.0))


if __name__ == '__main__':
    sys.exit(main())
