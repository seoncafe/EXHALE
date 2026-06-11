#!/usr/bin/env python3
"""Phase-3a: H II proton budget vs radius (cf. Huang et al. 2023 Figs. 11/27).

Compares the n=2 Balmer photoionization proton source against the OTHER
proton sources/sinks in the H ionization balance, all as VOLUMETRIC rates
[cm^-3 s^-1], read directly from output/Excited_H.txt (computed in-code from
the converged state):

    S_n2     = gamma2_bal * n2          (n=2 Balmer photoionization)   col8
    S_ground = P_HI_ground * nHI        (ground-state photoionization) col12
    S_coll   = a_ion_HI * ne * nHI      (collisional ionization)       col13
    S_recomb = alpha_B * ne * nHII      (recombination sink)           col14

Surfaces NUMBERS + FIGURE only; validation is the user's call.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DIR = 'output'

EH = np.loadtxt(f'{DIR}/Excited_H.txt')
r        = EH[:, 0]
T        = EH[:, 1]
nHI      = EH[:, 2]
ne       = EH[:, 3]
Sn2      = EH[:, 7]    # n=2 Balmer photoionization
Sground  = EH[:, 11]   # ground-state photoionization
Scoll    = EH[:, 12]   # collisional ionization
Srecomb  = EH[:, 13]   # recombination sink

m = (r >= 1.0)

# total proton SOURCE (the three ionization channels)
Stot = Sground + Scoll + Sn2

print("# r/Rp    S_n2       S_ground   S_coll     S_recomb   "
      "S_n2/S_grnd  S_n2/S_tot")
for rr in [1.0, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0, 2.5, 3.0]:
    i = np.argmin(np.abs(r - rr))
    fg = Sn2[i]/Sground[i] if Sground[i] > 0 else np.nan
    ft = Sn2[i]/Stot[i] if Stot[i] > 0 else np.nan
    print(f"{r[i]:6.3f}  {Sn2[i]:.3e}  {Sground[i]:.3e}  {Scoll[i]:.3e}  "
          f"{Srecomb[i]:.3e}  {fg:10.3e}  {ft:10.3e}")

frac_g = (Sn2/np.where(Sground > 0, Sground, np.nan))[m]
frac_t = (Sn2/np.where(Stot > 0, Stot, np.nan))[m]
print(f"\nmax S_n2/S_ground (physical) = {np.nanmax(frac_g):.3e} "
      f"at r/Rp = {r[m][np.nanargmax(frac_g)]:.3f}")
print(f"max S_n2/S_total  (physical) = {np.nanmax(frac_t):.3e} "
      f"at r/Rp = {r[m][np.nanargmax(frac_t)]:.3f}")

# ----- figure -----
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))

ax[0].plot(r[m], Sground[m], 'k-',   label=r'ground-state photoion. $P_{\rm HI}n_{\rm HI}$')
ax[0].plot(r[m], Sn2[m],     'C0-',  label=r'$n=2$ Balmer photoion. $S_{\rm n2}$')
ax[0].plot(r[m], Scoll[m],   'C2-.', label=r'collisional ioniz. $C n_e n_{\rm HI}$')
ax[0].plot(r[m], Srecomb[m], 'C3--', label=r'recombination $\alpha_B n_e n_{\rm HII}$')
ax[0].set_yscale('log'); ax[0].set_xlabel(r'$r/R_p$')
ax[0].set_ylabel(r'rate [cm$^{-3}$ s$^{-1}$]')
ax[0].set_title('H proton budget (cf. Huang Figs. 11/27)')
ax[0].legend(fontsize=8); ax[0].set_xlim(1, min(4, r.max()))

ax[1].plot(r[m], (Sn2/np.where(Sground > 0, Sground, np.nan))[m], 'C0-',
           label=r'$S_{\rm n2}/S_{\rm ground}$')
ax[1].plot(r[m], (Sn2/np.where(Stot > 0, Stot, np.nan))[m], 'C4--',
           label=r'$S_{\rm n2}/S_{\rm total\ source}$')
ax[1].set_yscale('log'); ax[1].set_xlabel(r'$r/R_p$')
ax[1].set_ylabel(r'fraction')
ax[1].set_title(r'$n=2$ photoionization share of the proton source')
ax[1].legend(fontsize=9); ax[1].set_xlim(1, min(4, r.max()))
ax[1].axhline(1.0, color='gray', lw=0.6, ls=':')

plt.tight_layout()
plt.savefig('proton_budget.png', dpi=130)
print("\nwrote proton_budget.png")
