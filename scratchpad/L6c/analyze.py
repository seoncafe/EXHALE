#!/usr/bin/env python3
"""Column-resolved comparison of two EXHALE output files.

Prints, for each column, the worst relative difference |a-b|/max(|a|,|b|) and
the row that carries it, with the two values and the radius of that row.
Rows beginning with '#' are skipped, so the row index matches the one
compare_within_tolerance.py reports.
"""
import sys

def read(path):
    cols, rows = None, []
    for line in open(path):
        s = line.strip()
        if s.startswith('# columns'):
            cols = s.split()[2:]
            continue
        if not s or s.startswith('#'):
            continue
        rows.append([float(x) for x in s.split()])
    return cols, rows

def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    thr = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    ca, A = read(a_path)
    cb, B = read(b_path)
    n = min(len(A), len(B))
    ncol = min(len(A[0]), len(B[0]))
    names = ca if ca and len(ca) >= ncol else ['c%d' % (j+1) for j in range(ncol)]
    print('%-14s %-10s %-6s %-9s %-14s %-14s' % ('column', 'max rel', 'row', 'r[Rp]', 'a', 'b'))
    for j in range(ncol):
        worst, wi = 0.0, -1
        for i in range(n):
            x, y = A[i][j], B[i][j]
            if x == y:
                continue
            d = abs(x - y) / max(abs(x), abs(y), 1e-300)
            if d > worst:
                worst, wi = d, i
        if worst > thr:
            print('%-14s %-10.3e %-6d %-9.4f %-14.6g %-14.6g'
                  % (names[j], worst, wi + 1, A[wi][0], A[wi][j], B[wi][j]))
    print('(columns with max rel <= %.1e omitted)' % thr)

main()
