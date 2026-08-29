#!/usr/bin/env python3
"""One row per scanned reservoir He/H: the deepest-level El/H departure and
the deep temperature the climate step returned."""
import glob
import os
import re
import sys

import numpy as np

SOLAR = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5,
             O=6.06178e-4)
HERE = os.path.dirname(os.path.abspath(__file__))

rows = []
for d in sorted(glob.glob(os.path.join(HERE, 'scan', 'heh*'))):
    heh = float(os.path.basename(d)[3:].replace('p', '.'))
    prof = os.path.join(d, 'lower_atmosphere_profile.dat')
    log = os.path.join(d, 'adapter.log')
    if not os.path.exists(prof):
        why = 'climate solve failed'
        txt = open(log).read() if os.path.exists(log) else ''
        m = re.search(r'ClimaException: (.*)', txt)
        if m:
            why = 'clima: ' + m.group(1).strip()[:44]
        rows.append((heh, None, None, why))
        continue
    cols, data = None, []
    for line in open(prof):
        if line.startswith('#'):
            if line.startswith('# columns:'):
                cols = line.split(':', 1)[1].split()
            continue
        data.append([float(x) for x in line.split()])
    a = np.array(data)
    want = dict(SOLAR)
    want['He'] = heh
    dev = {}
    for el in ('He', 'C', 'N', 'O'):
        k = 'X_' + el
        if k in cols:
            dev[el] = a[0, cols.index(k)]/want[el] - 1.0
    rows.append((heh, a[0, 2], dev, ''))

rows.sort()
print('%-8s %-8s %-11s %-11s %-11s %-11s  %-6s %s'
      % ('He/H', 'T_deep', 'dev_He', 'dev_C', 'dev_N', 'dev_O',
         'ratio', 'verdict at tol 1e-4'))
for heh, T, dev, why in rows:
    if dev is None:
        print('%-8.4f %-8s %s' % (heh, '-', why))
        continue
    r = abs(dev['N']/dev['He']) if dev['He'] else float('nan')
    worst = max(abs(v) for v in dev.values())
    print('%-8.4f %-8.1f %+11.3e %+11.3e %+11.3e %+11.3e  %-6.0f %s'
          % (heh, T, dev['He'], dev['C'], dev['N'], dev['O'], r,
             'REFUSED' if worst > 1.0e-4 else 'accepted'))
