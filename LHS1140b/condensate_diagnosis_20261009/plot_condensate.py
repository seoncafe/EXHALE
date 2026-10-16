#!/usr/bin/env python3
"""Figure of the condensate diagnosis: condensed and gas-phase H2O, NH3, CH4
and T against pressure for the variant columns.  Usage:
plot_condensate.py <out.png> <case_dir>:<label> ..."""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                   # noqa: E402

sys.path.insert(0, '/home/kiseon/.claude/jobs/f24ea156/tmp/condensate')
from analyze import column_report                    # noqa: E402

out = sys.argv[1]
cases = [a.split(':', 1) for a in sys.argv[2:]]
fig, axs = plt.subplots(2, 4, figsize=(15, 10), sharey=True)
ax = axs[0]
bx = axs[1]
ref = None
styles = ['-', '--', '-.', ':', (0, (5, 1, 1, 1))]
for (case, label), ls in zip(cases, styles):
    o = column_report(case)
    p = o['p']
    if ref is None:
        ref = o
    else:
        for k, (key, tag) in enumerate([('qv', 'H2O'), ('qNH3', 'NH3'),
                                        ('qCH4', 'CH4'), ('T', 'T')]):
            d = np.abs(o[key]/ref[key] - 1.0)
            bx[k].plot(np.maximum(d, 1e-12), p, ls=ls, label=label)
    qa = np.where(o['qa'] > 0, o['qa'], np.nan)
    ax[0].plot(qa, p, ls=ls, label=label)
    ax[1].plot(o['qv'], p, ls=ls, label=label)
    ax[2].plot(o['qNH3'], p, ls=ls, color='C0' if ls == '-' else None)
    ax[2].plot(o['qCH4'], p, ls=ls, color='C1' if ls == '-' else None)
    ax[3].plot(o['T'], p, ls=ls, label=label)
ax[0].set_xscale('log')
ax[0].set_xlim(1e-8, 1e4)
ax[0].set_xlabel(r'condensed H$_2$O, $n_{\rm c}/n_{\rm gas}$')
ax[0].axvline(2.34e-4, color='k', lw=0.6)
ax[0].text(3e-4, 15, r'deep $q_{\rm H_2O}$', fontsize=8)
ax[1].set_xscale('log')
ax[1].set_xlabel(r'gas-phase $q_{\rm H_2O}$')
ax[1].axvline(1.354e-7, color='k', lw=0.6)
ax[1].text(1.5e-7, 1e-4, r'climate cold trap', fontsize=8, rotation=90)
ax[2].set_xscale('log')
ax[2].set_xlabel(r'$q_{\rm NH_3}$ (lower), $q_{\rm CH_4}$ (upper)')
ax[3].set_xlabel(r'$T$ [K]')
for k, tag in enumerate([r'$q_{\rm H_2O}$', r'$q_{\rm NH_3}$',
                         r'$q_{\rm CH_4}$', r'$T$']):
    bx[k].set_xscale('log')
    bx[k].set_xlim(1e-12, 1e0)
    bx[k].set_xlabel(r'$|$' + tag + r'/' + tag + r'$_{\rm (a)} - 1|$')
bx[0].set_ylabel(r'$p$ [bar]')
bx[0].legend(fontsize=8, loc='upper right')
for a in list(ax) + list(bx):
    a.set_yscale('log')
    a.axhline(1.004, color='0.5', lw=0.6)
ax[0].set_ylim(20, 1e-6)
ax[0].set_ylabel(r'$p$ [bar]')
ax[0].legend(fontsize=8, loc='upper left')
ax[3].text(250, 0.8, r'energy floor 1.004 bar', fontsize=8)
fig.tight_layout()
fig.savefig(out, dpi=150)
print('wrote', out)
