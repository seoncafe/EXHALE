#!/usr/bin/env python3
"""Numeric comparison of a regression output against its golden.

WHY THIS EXISTS.  The regression's primary verdict is bitwise: a refactor with
no intended physics change must not move a single digit, and that is what
run_check.sh asks first.  But some changes move results by an amount that is
below any threshold the physics can notice -- a compiler option, unifying a
constant that two files had defined 0.03 per cent apart -- and spending a
golden refresh on those destroys the reference the whole matrix is measured
against.  The standing rule (2026-09-05) is that a difference under 0.1 per
cent counts as identical and the golden stays.  This script is how that rule
is measured rather than assumed: it reports the worst RELATIVE difference
between the two files, and which column and row carry it, so a verdict of
"identical within tolerance" is always backed by a number.

Usage:  compare_within_tolerance.py <current> <golden> [rel_tol]
        rel_tol defaults to 1e-3 (0.1 per cent).

Exit 0  numerically equal within rel_tol (or bitwise identical)
Exit 1  a value exceeds rel_tol
Exit 2  the files cannot be compared at all (shape, unreadable, non-numeric)

The relative measure is |a-b| / max(|a|,|b|), so a column whose golden value is
zero and whose new value is not counts as a full difference rather than being
divided by zero.  Rows beginning with '#' are skipped, as in run_check.sh.
"""
import sys


def read(path):
    rows = []
    try:
        fh = open(path)
    except OSError as e:
        print('CANNOT READ %s (%s)' % (path, e.strerror))
        raise SystemExit(2)
    with fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                raise SystemExit(2)
    return rows


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        raise SystemExit(2)
    cur, gold = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e-3
    a, b = read(cur), read(gold)
    if len(a) != len(b):
        print('SHAPE rows %d vs %d' % (len(a), len(b)))
        raise SystemExit(2)
    worst, wrow, wcol, wa, wb = 0.0, 0, 0, 0.0, 0.0
    for i, (ra, rb) in enumerate(zip(a, b), 1):
        if len(ra) != len(rb):
            print('SHAPE row %d: %d vs %d columns' % (i, len(ra), len(rb)))
            raise SystemExit(2)
        for j, (x, y) in enumerate(zip(ra, rb), 1):
            if x == y:
                continue
            d = abs(x - y) / max(abs(x), abs(y), 1.0e-300)
            if d > worst:
                worst, wrow, wcol, wa, wb = d, i, j, x, y
    if worst == 0.0:
        print('IDENTICAL')
        raise SystemExit(0)
    where = 'row %d col %d (%.17g vs %.17g)' % (wrow, wcol, wa, wb)
    if worst <= tol:
        print('WITHIN %.1e  max rel %.3e at %s' % (tol, worst, where))
        raise SystemExit(0)
    print('EXCEEDS %.1e  max rel %.3e at %s' % (tol, worst, where))
    raise SystemExit(1)


if __name__ == '__main__':
    main()
