"""Project the density-dependent Fe II cooling onto the WASP-121b profile.

Reads the converged run's Cooling_breakdown.txt (r, T, ne, per-channel cooling).
For each cell it forms the multilevel statistical-equilibrium Fe II cooling
coefficient Lambda_SE(T, ne) (cooling_data/fe2_cooling.py) and the coronal
coefficient Lambda_cor(T), then rescales the tabulated Fe II cooling channel by
Lambda_SE/Lambda_cor (both share the same beta_esc * ne * n_FeII prefactor), so
no ATES re-run is needed to see the impact.  Shows how much the coronal table
overestimates Fe II cooling vs radius, and the resulting change in total cooling.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'cooling_data')))

from chianti_cooling import cooling_effective                       # noqa: E402
from fe2_cooling import _load_fe2, fe2_levels_cooling               # noqa: E402

COOLF = os.path.join(HERE, 'output', 'Cooling_breakdown.txt')
ORDER = ('CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII '
         'CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII').split()


def main():
    d = np.loadtxt(COOLF)
    r, T, ne, tot = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    col = {k: d[:, 8 + i] for i, k in enumerate(ORDER)}
    feii = col['FeII']

    lev, sc, A = _load_fe2()
    lam_cor = cooling_effective('fe', 2, T, pop='coronal')
    ratio = np.ones_like(T)
    # only need cells where Fe II actually contributes
    sig = feii > 1e-6 * np.nanmax(tot)
    for j in np.where(sig)[0]:
        lam_se = fe2_levels_cooling(T[j], max(ne[j], 1.0), lev, sc, A)
        ratio[j] = (lam_se / max(ne[j], 1.0)) / max(lam_cor[j], 1e-99)
    ratio = np.clip(ratio, 0.0, 1.0)

    feii_corr = feii * ratio
    tot_corr = tot - feii + feii_corr

    sel = (r >= 1.0) & (r <= 1.5)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    a1.plot(r[sel], feii[sel], 'r--', lw=2, label=r'Fe\,II coronal (current table)')
    a1.plot(r[sel], feii_corr[sel], 'b-', lw=2.2,
            label=r'Fe\,II density-correct (SE)')
    a1.plot(r[sel], tot[sel], 'k:', lw=1.5, label='total cooling (current)')
    a1.plot(r[sel], tot_corr[sel], 'g-', lw=1.3, label='total (Fe\\,II corrected)')
    a1.set_yscale('log')
    a1.set_xlim(1.0, 1.5)
    ymax = np.nanmax(tot[sel])
    a1.set_ylim(ymax * 1e-6, ymax * 3)
    a1.set_xlabel(r'$r/R_p$')
    a1.set_ylabel(r'cooling rate $\mathrm{[erg\,cm^{-3}\,s^{-1}]}$')
    a1.set_title(r'Fe\,II cooling vs radius: coronal vs density-correct')
    a1.legend(loc='upper right', frameon=False, fontsize=9)
    a1.grid(True, which='both', alpha=0.25)

    a2.plot(r[sel], 1.0 / np.clip(ratio[sel], 1e-12, None), 'b-', lw=2)
    a2.set_yscale('log')
    a2.set_xlim(1.0, 1.5)
    a2.set_xlabel(r'$r/R_p$')
    a2.set_ylabel(r'coronal / density-correct (overestimate factor)')
    a2.set_title(r'Fe\,II coronal overestimate along the profile')
    a2.grid(True, which='both', alpha=0.25)
    ax2b = a2.twinx()
    ax2b.plot(r[sel], T[sel], color='0.6', lw=1.2, ls='--')
    ax2b.set_ylabel(r'$T\ \mathrm{[K]}$ (grey dashed)', color='0.5')

    fig.tight_layout()
    out = os.path.join(HERE, 'FeII_density_correction.png')
    fig.savefig(out, dpi=130)
    print('wrote', out)

    # numeric summary
    print('\n# r/Rp    T[K]      ne         FeII_cor    FeII_SE     overest')
    for rt in (1.00, 1.02, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40):
        j = int(np.argmin(abs(r - rt)))
        print('%6.3f  %8.1f  %.2e  %.3e  %.3e  %8.1f' %
              (r[j], T[j], ne[j], feii[j], feii_corr[j],
               1.0 / max(ratio[j], 1e-12)))
    base = int(np.argmin(abs(r - 1.0)))
    print('\nTotal cooling at base r=1.0: current %.3e -> corrected %.3e (%.1fx lower)'
          % (tot[base], tot_corr[base], tot[base] / max(tot_corr[base], 1e-99)))


if __name__ == '__main__':
    main()
