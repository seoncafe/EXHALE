#!/usr/bin/env python3
"""Phase 3b (a): H(n=2) heating budget for the mode-2 escape-probability run with
collisional de-excitation heating on. Compare to Huang et al. 2023 Figs. 10/26.

Excited_H.txt: col9 Hpe (Balmer photoelectric), col10 Hdx (collisional de-exc).
Hydro_ioniz.txt: col6 = total heating [erg cm^-3 s^-1] (already physical).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

e = np.loadtxt('output_m2_deexc/Excited_H.txt')
h = np.loadtxt('output_m2_deexc/Hydro_ioniz.txt')
r = e[:, 0]
Hpe, Hdx = e[:, 8], e[:, 9]
Hbal = Hpe + Hdx
Htot = h[:, 5]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

ax[0].plot(r, Htot, 'k-', label=r'total heating')
ax[0].plot(r, Hbal, 'b-', label=r'H($n{=}2$) total (Hpe$+$Hdx)')
ax[0].plot(r, Hpe, 'g--', label=r'Hpe (photoelectric)')
ax[0].plot(r, Hdx, 'r:', label=r'Hdx (de-excitation)')
ax[0].set_yscale('log')
ax[0].set_ylim(1e-10, 3e-4)
ax[0].set_ylabel(r'heating [erg cm$^{-3}$ s$^{-1}$]')

frac = np.where(Htot > 0, 100.0 * Hbal / Htot, 0.0)
ax[1].plot(r, frac, 'b-')
ax[1].set_ylabel(r'H($n{=}2$) heating / total [\%]')

for a in ax:
    a.set_xlabel(r'$r/R_p$')
    a.set_xlim(1.0, 2.2)
    a.axvspan(1.1, 1.8, color='0.92', zorder=0)
ax[0].legend(fontsize=8)

fig.suptitle(r'WASP-121b: H($n{=}2$) heating budget, mode-2 escape-prob.\ $+$ de-excitation '
             r'(cf.\ Huang Fig.\ 10)')
fig.tight_layout()
fig.savefig('lya_heating_budget.png', dpi=130)
print('wrote lya_heating_budget.png')
