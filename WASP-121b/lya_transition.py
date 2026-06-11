#!/usr/bin/env python3
"""Phase 3b (b): smoothing the stellar-penetration photosphere transition.
Flat stellar profile T_* = max(0, 1-x1/Xs) (kink) vs Gaussian
T_* = erfc(x1/(sqrt2 Xs)) (smooth). Both runs are mode-2 + de-excitation heating.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

m0 = np.loadtxt('output_pre3b_baseline/Excited_H.txt')   # mode 0 (reference)
f = np.loadtxt('output_m2_deexc/Excited_H.txt')          # mode 2, flat T*
g = np.loadtxt('output_m2_gauss/Excited_H.txt')          # mode 2, Gaussian T*
r = f[:, 0]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

ax[0].plot(r, m0[:, 4], 'k:', label=r'mode 0')
ax[0].plot(r, f[:, 4], 'r--', label=r'mode 2, flat $T_*$ (kink)')
ax[0].plot(r, g[:, 4], 'b-', label=r'mode 2, Gaussian $T_*$ (smooth)')
ax[0].set_ylabel(r'$\bar J_{\rm Ly\alpha}$ [cgs]')
ax[0].set_yscale('log')
ax[0].set_ylim(1e-12, 1e-7)

ax[1].plot(r, f[:, 7] / np.maximum(f[:, 11], 1e-99), 'r--', label=r'flat $T_*$')
ax[1].plot(r, g[:, 7] / np.maximum(g[:, 11], 1e-99), 'b-', label=r'Gaussian $T_*$')
ax[1].set_ylabel(r'$S(n{=}2)\,/\,S_{\rm ground}$')

for a in ax:
    a.set_xlabel(r'$r/R_p$')
    a.set_xlim(1.0, 1.6)
    a.legend(fontsize=8)
    a.axvline(1.2, color='0.8', lw=1, zorder=0)

fig.suptitle(r'WASP-121b: stellar-penetration transition --- flat vs Gaussian $T_*$ '
             r'(mode 2 $+$ de-excitation)')
fig.tight_layout()
fig.savefig('lya_transition.png', dpi=130)
print('wrote lya_transition.png')
