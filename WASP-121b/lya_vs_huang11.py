#!/usr/bin/env python3
"""Phase 3b calibration: my n=2 photoionization fraction S(n=2)/S_ground vs
Huang et al. 2023 Fig. 11 (H(n=2) dotted / Photoionization-H solid), case A.

Huang Fig. 11 ratio digitized by eye from the log plot (approximate):
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Huang Fig. 11: H(n=2)/total-photoionization ratio, careful by-eye digitization
# (dotted/solid vertical gap on the log plot; ~factor-2 uncertain from the figure).
hr = np.array([1.3, 1.4, 1.5, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0])
hf = np.array([0.10, 0.30, 0.45, 0.50, 0.40, 0.40, 0.35, 0.32, 0.30, 0.28, 0.30, 0.30, 0.40])

m = np.loadtxt('output_caseA/Excited_H.txt')   # mode 2, spherical Case A to 13 Rp
r = m[:, 0]
RH = 1.25*r; mine = m[:, 7] / np.maximum(m[:, 11], 1e-99)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for a in ax:
    a.plot(r, mine, 'b-', label=r'mode 2 (boost=12)')
    a.plot(hr, hf, 'ks--', mfc='none', label=r'Huang 2023 Fig.\ 11 (digitized)')
    a.set_xlabel(r'$r/R_p$')
    a.set_ylabel(r'$S(n{=}2)\,/\,S_{\rm ground}$')
    a.legend(fontsize=8)
ax[0].set_xlim(1.0, 2.2)
ax[0].set_title(r'wind (linear)')
ax[1].set_xlim(1.3, 13)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_title(r'full range (log--log)')

fig.suptitle(r'WASP-121b: H($n{=}2$) photoionization fraction vs Huang Fig.\ 11')
fig.tight_layout()
fig.savefig('lya_vs_huang11.png', dpi=130)
print('wrote lya_vs_huang11.png')
# quick numeric gap
for rad in [1.3, 1.5, 2.0, 3.0]:
    i = int(np.argmin(np.abs(r - rad)))
    h = np.interp(rad, hr, hf)
    print(f'  r={rad:.1f}  mine={mine[i]:.3f}  Huang~{h:.2f}  (mine/Huang={mine[i]/h:.2f})')
