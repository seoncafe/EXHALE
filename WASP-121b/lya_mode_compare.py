#!/usr/bin/env python3
"""Phase 3b: compare the H(n=2) pumping between the Phase-3a parameterized
Ly-alpha field (jlya_mode=0, output_pre3b_baseline/) and the in-line
escape-probability RT (jlya_mode=2, output_m2/) on WASP-121b.

The Phase-3a caveat was that the parameterized field peaks n(2p) near 2 R_p and
under-pumps the 1.1-1.8 R_p band. This figure shows J_lya(r), n(2p)(r), the n=2
photoionization proton source, and its fraction of the ground-state source
(compare Huang et al. 2023 Figs. 11/27).

Excited_H.txt columns (1-indexed): 1 r/Rp, 2 T, 3 nHI, 4 ne, 5 Jlya, 6 n2s,
7 n2p, 8 Sproton(n=2), 9 Hpe, 10 Hdx, 11 tau_Lya, 12 Sgrnd_photoion,
13 Scoll_ion, 14 Srecomb.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

m0 = np.loadtxt('output_pre3b_baseline/Excited_H.txt')   # jlya_mode 0 (param.)
m2 = np.loadtxt('output_m2c/Excited_H.txt')              # jlya_mode 2, trapped stellar
r0, r2 = m0[:, 0], m2[:, 0]

fig, ax = plt.subplots(2, 2, figsize=(11, 8))

ax[0, 0].plot(r0, m0[:, 4], 'k--', label=r'mode 0 (parameterized)')
ax[0, 0].plot(r2, m2[:, 4], 'b-', label=r'mode 2 (escape-prob.)')
ax[0, 0].set_ylabel(r'$\bar J_{\rm Ly\alpha}$ [cgs]')
ax[0, 0].set_yscale('log')

ax[0, 1].plot(r0, m0[:, 6], 'k--', label=r'mode 0')
ax[0, 1].plot(r2, m2[:, 6], 'b-', label=r'mode 2')
ax[0, 1].set_ylabel(r'$n_{2p}$ [cm$^{-3}$]')
ax[0, 1].set_yscale('log')

ax[1, 0].plot(r0, m0[:, 7], 'k--', label=r'mode 0: $n{=}2$ source')
ax[1, 0].plot(r2, m2[:, 7], 'b-', label=r'mode 2: $n{=}2$ source')
ax[1, 0].plot(r2, m2[:, 11], 'r:', label=r'ground-state photoion.')
ax[1, 0].set_ylabel(r'proton source [cm$^{-3}$ s$^{-1}$]')
ax[1, 0].set_yscale('log')

ax[1, 1].plot(r0, m0[:, 7] / np.maximum(m0[:, 11], 1e-99), 'k--', label=r'mode 0')
ax[1, 1].plot(r2, m2[:, 7] / np.maximum(m2[:, 11], 1e-99), 'b-', label=r'mode 2')
ax[1, 1].set_ylabel(r'$S(n{=}2)\,/\,S_{\rm ground}$')

for a in ax.flat:
    a.set_xlabel(r'$r/R_p$')
    a.set_xlim(1.0, 2.2)
    a.legend(fontsize=8)
    a.axvspan(1.1, 1.8, color='0.9', zorder=0)

fig.suptitle(r'WASP-121b: H($n{=}2$) pumping --- parameterized vs escape-probability '
             r'Ly$\alpha$ (shaded $1.1$--$1.8\,R_p$)')
fig.tight_layout()
fig.savefig('lya_mode_compare.png', dpi=130)
print('wrote lya_mode_compare.png')
