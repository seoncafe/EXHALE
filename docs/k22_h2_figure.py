#!/usr/bin/env python3
"""The H2 extent of the hot Uranus: Koskinen et al. (2022) Model A against every
EXHALE state of that configuration (docs/koskinen2022_model_a_comparison.tex,
section "The H2 extent, 2026-09-12").

Model A: benchmarks/koskinen2022_model_a/model_a_fig8_digitized.txt (the
2026-09-12 pixel digitization of their Figure 8) and the 2026-09-05 readings of
docs/p23_published_profiles.md 5.1 (circles). EXHALE: the matched marching run
(benchmarks/koskinen2022_model_a/matched_hnu_minus_I/output; NOT the
superseded Roche-domain state in benchmarks/koskinen2022_model_a/output, see
the README there), the carrier reload fixture's pinned
state (backup/regression/carrier_elem_newton/IC) and the two stationary-route
states pinned beside the benchmark (carrier_reload_states/pass12, fixed_wind).
Writes docs/figures/k22_model_a/h2_extent.pdf.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
sys.path.insert(0, os.path.join(ROOT, 'examples'))
import exhale_io as aio

dig = np.loadtxt(os.path.join(ROOT, 'benchmarks/koskinen2022_model_a/model_a_fig8_digitized.txt'))
rb, lH2, lH, fH2 = dig[:, 1], dig[:, 2], dig[:, 3], dig[:, 5]
K05 = np.array([[1.50, 10.86, 9.59], [1.60, 10.12, 9.22], [1.75, 9.31, 8.88], [2.00, 8.38, 8.38],
                [2.30, 7.82, 7.94], [2.70, 7.31, 7.56], [3.20, 6.92, 7.26], [4.00, 6.51, 6.92],
                [5.00, 6.13, 6.61], [6.00, 5.86, 6.39], [7.50, 5.55, 6.12], [9.60, 5.19, 5.82]])

def load(d):
    hf = os.path.join(d, 'Hydro_ioniz.txt'); sf = os.path.join(d, 'Ion_species.txt')
    if not os.path.exists(hf):
        hf = os.path.join(d, 'Hydro_ioniz_IC.txt'); sf = os.path.join(d, 'Ion_species_IC.txt')
    r, ion = aio.load_ions(sf)
    nH = ion['HI'] + ion['HII'] + 2*(ion['H2'] + ion['H2p']) + 3*ion['H3p'] + ion['HeHp']
    return r, ion['H2'], ion['HI'], 2*ion['H2']/nH

RUNS = [('matched run (2026-09-05, spherical 7.24, Rate/4 + Mdot, hnu-I)', 'benchmarks/koskinen2022_model_a/matched_hnu_minus_I/output', 'k-', 1.8),
        ('carrier reload, pinned state (Roche 4.7, Rate/2 + Mdot/2)', 'backup/regression/carrier_elem_newton/IC', '-', 1.2),
        ('carrier reload, 12 bounded passes', 'benchmarks/koskinen2022_model_a/carrier_reload_states/pass12', '--', 1.2),
        ('carrier reload, fixed-wind fixed point', 'benchmarks/koskinen2022_model_a/carrier_reload_states/fixed_wind', ':', 1.6)]
COLS = ['k', 'C3', 'C1', 'C0']

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
m = rb > 1.08
ax[0].plot(rb[m], lH2[m], color='C2', lw=2.2, label=r'Model A H$_2$ (Fig. 8, 2026-09-12 digitization)')
ax[0].plot(rb[m], lH[m], color='C2', lw=1.0, ls='--', label='Model A H')
ax[0].plot(K05[:, 0]/1.34, K05[:, 1], 'o', mfc='none', color='C2', ms=5, label='2026-09-05 readings')
ax[0].plot(K05[:, 0]/1.34, K05[:, 2], 's', mfc='none', color='C2', ms=4)
ax[1].plot(rb[m], fH2[m], color='C2', lw=2.2, label='Model A')
for (lab, d, ls, lw), c in zip(RUNS, COLS):
    r, n2, nh, f = load(os.path.join(ROOT, d))
    ax[0].plot(r, np.log10(np.maximum(n2, 1e-30)), ls=ls.strip('k'), color=c, lw=lw, label=lab)
    ax[1].plot(r, f, ls=ls.strip('k'), color=c, lw=lw, label=lab)
ax[0].set_xlim(1.0, 7.2); ax[0].set_ylim(2, 13.5)
ax[0].set_xlabel(r'$r / r_{\rm base}$'); ax[0].set_ylabel(r'$\log_{10} n\ [{\rm cm}^{-3}]$')
ax[0].legend(fontsize=7, loc='upper right')
ax[1].set_xlim(1.0, 4.8); ax[1].set_ylim(0, 1.02)
ax[1].set_xlabel(r'$r / r_{\rm base}$'); ax[1].set_ylabel(r'$f({\rm H}_2) = 2 n({\rm H}_2) / n_{\rm H}$')
ax[1].legend(fontsize=7, loc='upper right')
for a in ax: a.grid(alpha=0.3)
fig.tight_layout()
fdir = os.environ.get('K22_FIG_DIR', os.path.join(HERE, 'figures', 'k22_model_a'))
os.makedirs(fdir, exist_ok=True)
fig.savefig(os.path.join(fdir, 'h2_extent.pdf'))

# the table of the memo
rows = [1.2, 1.5, 2.0, 3.0, 4.0]
def at(x, r, y): return np.interp(x, r, y)
print('r/r_base  ModelA_fH2  ' + '  '.join(lab.split(',')[0][:20] for lab, *_ in RUNS))
for x in rows:
    vals = [at(x, rb, fH2)] + [at(x, *load(os.path.join(ROOT, d))[::3]) for _, d, *_ in RUNS]
    print(f'{x:5.1f}  ' + '  '.join(f'{v:6.3f}' for v in vals))
