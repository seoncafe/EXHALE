#!/usr/bin/env python3
"""State the base grid explicitly in every input that relied on the old default.

Until 2026-09-19 the default width of the uniform base cells of the `Mixed`
grid was the single-precision neighbor of 2e-4, the binary64 value
1.9999999494757503e-4 R_p, with 50 uniform cells. From 2026-09-19 the default
is the exact double 2.0e-4, which builds another grid (cell centers displaced
by 3.8e-9 to 6.6e-9 relative on the grids measured, above the 1e-10 at which `load_IC` refuses a
state). An input written before that date with no `Base grid [dr,cells]:` line
would therefore build a grid on which its own stored results do not load.

This script writes the old resolved grid into those inputs:

    # Base grid of the stored results: the Mixed-grid default before 2026-09-19
    Base grid [dr,cells]: 1.9999999494757503e-4 50

immediately after the `Grid type:` line. The decimal
1.9999999494757503e-4 is read by gfortran's list-directed read into a real*8
as the bits 3F2A36E2E0000000, the bits of the default-real literal 2.0e-4
widened to double (MEASURED). 50 is the only value `N_low_cells` resolves to
without the key: `input_read` sets it to 50 and only the `Base grid` key
changes it.

WHICH INPUTS. A file is a candidate when it is named `input.inp` (or
`input_template.inp`, from which a flux-closure rung writes the `input.inp`
of each iterate), states `Grid type: Mixed` and has no non-comment
`Base grid` line. Where it lies decides what happens to it:

  LIVE      pinned: directories whose stored results are reproduced or
            continued (the LHS 1140 b catalog, the regression cases, the
            planet folders, the benchmarks, vulcan_work, and the examples
            that hold a stored state or are seeded from one).
  FROZEN    not edited: preserved records. Each tree root carries a
            GRID_DEFAULT_NOTE.md that states the line to add instead.
  TEMPLATE  not pinned: an example with no stored state is a starting point
            for new configurations and takes the present default.
  TRANSIENT not pinned: build/ and scratchpad/ hold products of tests and
            earlier work directories, rewritten by whatever makes them.

Nothing under a `states/` or `runs/` directory is ever opened for writing
(published generations and run records are immutable).

Usage (from anywhere; paths are relative to the EXHALE_v1.00 root):
    python3 src/utils/pin_base_grid.py            # dry run, prints the plan
    python3 src/utils/pin_base_grid.py --execute  # writes the LIVE pins
    python3 src/utils/pin_base_grid.py --list     # also lists every file
A second --execute writes nothing: a pinned file carries the key.
"""
import argparse
import collections
import os
import re
import sys

PIN_WIDTH = '1.9999999494757503e-4'
PIN_CELLS = 50
PIN_LINE = 'Base grid [dr,cells]: %s %d' % (PIN_WIDTH, PIN_CELLS)
PIN_COMMENT = ('# Base grid of the stored results: the Mixed-grid default'
               ' before 2026-09-19')

NAMES = ('input.inp', 'input_template.inp')

# Preserved records: never edited, one GRID_DEFAULT_NOTE.md at each root.
FROZEN_ROOTS = (
    'LHS1140b/archive_20260830',
    'LHS1140b/models_20260914_preL21',
    'LHS1140b/models_20260915_db87',
    'LHS1140b/examples',
    'backup/HD209458b_test',
    'backup/phase_d_baseline',
    'backup/lart_runs',
    'backup/example_HD209458b',
    'backup/regression/_quarantined',
    'docs/audit_20260905',
    'docs/lower_atmosphere_figs',
    'docs/version_compare',
)
LIVE_ROOTS = (
    'LHS1140b/models',
    'backup/regression',
    'HD209458b', 'HD189733b', 'WASP-52b', 'WASP-121b',
    'benchmarks',
    'vulcan_work',
)
TRANSIENT_ROOTS = ('build', 'scratchpad')
# Examples without a stored state of their own that the examples/README.md
# recipe seeds from the stored state of a pinned example: 04 restarts from
# 03_newton/output/, 16 from 15_molecular/output/. A seed and the run it
# seeds must build one grid, so these follow their seed.
EXAMPLE_SEEDED = {'examples/04_newton_from_state': 'examples/03_newton',
                  'examples/16_molecular_metals': 'examples/15_molecular'}


def label(line, key):
    """lbl_match of input_read.f90: key at the start of the left-trimmed line,
    followed by the end of the line or a separator (':', '?', ' ', tab, '=')."""
    t = line.lstrip(' \t')
    if not t.startswith(key):
        return False
    rest = t[len(key):].rstrip('\r\n')
    return rest == '' or rest[0] in ':? \t='


def is_comment(line):
    return line.lstrip(' \t').startswith('#')


def grid_type(lines):
    for ln in lines:
        if is_comment(ln) or not label(ln, 'Grid type'):
            continue
        words = ln.split()
        return words[2] if len(words) > 2 else ''
    return None


def has_base_grid(lines):
    return any(label(ln, 'Base grid') for ln in lines if not is_comment(ln))


def under(rel, root):
    return rel == root or rel.startswith(root + '/')


def has_stored_state(d):
    for name in os.listdir(d):
        p = os.path.join(d, name)
        if name.startswith('output') and os.path.isdir(p):
            if any(f.startswith('Hydro_ioniz') and f.endswith('.txt')
                   for f in os.listdir(p)):
                return True
    return False


def classify(root, rel):
    """(class, area) of the file at root/rel."""
    parts = rel.split('/')
    if 'states' in parts[:-1] or 'runs' in parts[:-1]:
        return 'GENERATION', parts[0]
    for r in FROZEN_ROOTS:
        if under(rel, r):
            return 'FROZEN', r
    for r in TRANSIENT_ROOTS:
        if under(rel, r):
            return 'TRANSIENT', r
    for r in LIVE_ROOTS:
        if under(rel, r):
            return 'LIVE', r
    if under(rel, 'examples'):
        d = os.path.dirname(rel)
        if has_stored_state(os.path.join(root, d)):
            return 'LIVE', 'examples'
        if d in EXAMPLE_SEEDED:
            return 'LIVE', 'examples'
        return 'TEMPLATE', 'examples'
    return 'UNCLASSIFIED', parts[0]


def candidates(root):
    for dirpath, dirnames, filenames in os.walk(root):
        rel_dir = os.path.relpath(dirpath, root)
        # Published generations and run records are never walked into.
        dirnames[:] = [x for x in dirnames if x not in ('states', 'runs')]
        for name in filenames:
            if name in NAMES:
                rel = os.path.normpath(os.path.join(rel_dir, name))
                yield rel


def pinned_text(lines):
    """The file with the comment and the key after its Grid type line."""
    nl = '\r\n' if lines and lines[0].endswith('\r\n') else '\n'
    out = []
    done = False
    for ln in lines:
        out.append(ln)
        if not done and not is_comment(ln) and label(ln, 'Grid type'):
            if not ln.endswith('\n'):
                out[-1] = ln + nl
            out.append(PIN_COMMENT + nl)
            out.append(PIN_LINE + nl)
            done = True
    return ''.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--execute', action='store_true',
                    help='write the LIVE pins (default: dry run)')
    ap.add_argument('--list', action='store_true',
                    help='list every candidate with its class')
    ap.add_argument('--root', default=os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', '..')),
        help='the EXHALE_v1.00 directory (default: this script\'s tree)')
    a = ap.parse_args()
    root = a.root

    count = collections.Counter()      # (class, area) -> files needing a pin
    keyed = collections.Counter()      # (class, area) -> Mixed with the key
    other = collections.Counter()      # (class, area) -> not a Mixed grid
    templates = []
    unclassified = []
    written = 0
    for rel in sorted(candidates(root)):
        cls, area = classify(root, rel)
        path = os.path.join(root, rel)
        with open(path, errors='replace', newline='') as fh:
            lines = fh.readlines()
        gt = grid_type(lines)
        if gt != 'Mixed':
            other[(cls, area)] += 1
            continue
        if has_base_grid(lines):
            keyed[(cls, area)] += 1
            continue
        count[(cls, area)] += 1
        if a.list:
            print('%-12s %s' % (cls, rel))
        if cls == 'TEMPLATE':
            templates.append(rel)
        if cls == 'UNCLASSIFIED':
            unclassified.append(rel)
        if cls == 'LIVE' and a.execute:
            text = pinned_text(lines)
            with open(path, 'w', newline='') as fh:
                fh.write(text)
            written += 1

    mode = 'EXECUTE' if a.execute else 'DRY RUN'
    print('pin_base_grid (%s): line "%s" after "Grid type: Mixed"'
          % (mode, PIN_LINE))
    print('%-11s %-34s %9s %9s %9s' % ('class', 'area', 'no key',
                                       'has key', 'not Mixed'))
    keys = sorted(set(count) | set(keyed) | set(other))
    tot = collections.Counter()
    for k in keys:
        print('%-11s %-34s %9d %9d %9d'
              % (k[0], k[1], count[k], keyed[k], other[k]))
        tot[k[0]] += count[k]
    print('files with Mixed grid and no key, by class: ' +
          ', '.join('%s %d' % (c, tot[c]) for c in sorted(tot)))
    if templates:
        print('TEMPLATE inputs left on the present default:')
        for t in templates:
            print('   ' + t)
    if unclassified:
        print('UNCLASSIFIED inputs (not pinned; classify them first):')
        for t in unclassified:
            print('   ' + t)
    if a.execute:
        print('written: %d' % written)
    else:
        print('would write: %d (the LIVE row total); run with --execute'
              % tot['LIVE'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
