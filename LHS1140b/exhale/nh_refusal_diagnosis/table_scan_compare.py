#!/usr/bin/env python3
"""The reservoir scan on both Photochem builds, side by side.

Left half is `scan/`, written by the conda photochem 0.8.4 package whose
equilibrium solver scales every elemental residual by the largest elemental
abundance.  Right half is `scan_pc090/`, written by the corrected 0.9.0 build
whose Equilibrate convergence test is element-relative.  Both were run with
the handoff check disarmed (--abundance-tol 1.0), so a row exists wherever the
climate and photochemical solve completed and the departure is the one the
check would have seen.

The last column is the deep temperature difference, which is the only thing
the Clima patch can move on a composition that already had a solution.
"""
import glob
import os
import re

import numpy as np

SOLAR = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5, O=6.06178e-4)
HERE = os.path.dirname(os.path.abspath(__file__))


def read(d):
    prof = os.path.join(d, 'lower_atmosphere_profile.dat')
    log = os.path.join(d, 'adapter.log')
    if not os.path.exists(prof):
        why = 'no profile'
        if os.path.exists(log):
            m = re.search(r'ClimaException: (.*)', open(log).read())
            if m:
                why = 'clima: ' + m.group(1).strip()[:34]
            else:
                m = re.search(r'RuntimeError: (.*)', open(log).read())
                if m:
                    why = 'refused: ' + m.group(1).strip()[:34]
        return None, why
    cols, data = None, []
    for line in open(prof):
        if line.startswith('#'):
            if line.startswith('# columns:'):
                cols = line.split(':', 1)[1].split()
            continue
        data.append([float(x) for x in line.split()])
    return (cols, np.array(data)), ''


def worst_dev(cols, a, heh):
    want = dict(SOLAR)
    want['He'] = heh
    out = {}
    for el in ('He', 'C', 'N', 'O'):
        k = 'X_' + el
        if k in cols:
            out[el] = a[0, cols.index(k)]/want[el] - 1.0
    el = max(out, key=lambda e: abs(out[e]))
    return abs(out[el]), el, a[0, 2]


def main():
    tags = set()
    for base in ('scan', 'scan_pc090'):
        for d in glob.glob(os.path.join(HERE, base, 'heh*')):
            tags.add(os.path.basename(d))
    rows = []
    for t in tags:
        heh = float(t[3:].replace('p', '.'))
        rows.append((heh, t))
    rows.sort()

    print('%-8s | %-11s %-4s %-8s | %-11s %-4s %-8s | %-8s'
          % ('He/H', '0.8.4 worst', 'el', 'T_deep', '0.9.0 worst', 'el',
             'T_deep', 'dT[K]'))
    print('-'*8 + '-+-' + '-'*25 + '-+-' + '-'*25 + '-+-' + '-'*8)
    new_worst, new_all = 0.0, []
    for heh, t in rows:
        cells = []
        Ts = []
        for base in ('scan', 'scan_pc090'):
            r, why = read(os.path.join(HERE, base, t))
            if r is None:
                cells.append('%-25s' % why[:25])
                Ts.append(None)
                continue
            w, el, T = worst_dev(r[0], r[1], heh)
            cells.append('%-11.3e %-4s %-8.3f' % (w, el, T))
            Ts.append(T)
            if base == 'scan_pc090':
                new_all.append(w)
                new_worst = max(new_worst, w)
        dT = ('%+8.3f' % (Ts[1]-Ts[0])) if (Ts[0] and Ts[1]) else '%8s' % '-'
        print('%-8.4f | %s | %s | %s' % (heh, cells[0], cells[1], dT))
    print()
    print('0.9.0: %d states, largest deepest-level departure %.4e, '
          'median %.4e' % (len(new_all), new_worst, float(np.median(new_all))))


if __name__ == '__main__':
    main()
