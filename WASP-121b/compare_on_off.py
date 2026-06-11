#!/usr/bin/env python3
"""Phase-3a ON vs OFF comparison: temperature and H ionization profiles.

ON  = excited-H coupling enabled  (output_on_clean/)
OFF = excited-H coupling disabled (output_off_clean/)
Both produced from the SAME build; the only difference is the presence of the
stellar Teff/radius tail lines in input.inp (use_excited_H toggle).

Surfaces NUMBERS + FIGURE only; validation is the user's call.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(d):
    HY = np.loadtxt(f'{d}/Hydro_ioniz.txt')
    ION = np.loadtxt(f'{d}/Ion_species.txt')
    r = HY[:, 0]
    T = HY[:, 4]                       # col5 = T [K]
    # Ion_species.txt: col1=r/Rp, col2=nHI, col3=nHII, col4=nHeI ...
    nHI = ION[:, 1]
    nHII = ION[:, 2]
    xHII = nHII / np.maximum(nHI + nHII, 1e-300)
    return r, T, nHI, nHII, xHII


rON, TON, nHI_ON, nHII_ON, xON = load('output_on_clean')
rOFF, TOFF, nHI_OFF, nHII_OFF, xOFF = load('output_off_clean')

# physical domain mask (drop ghosts)
m = (rON >= 1.0)

# interpolate OFF onto ON grid for differencing (grids should match, but be safe)
TOFF_i = np.interp(rON, rOFF, TOFF)
xOFF_i = np.interp(rON, rOFF, xOFF)
nHI_OFF_i = np.interp(rON, rOFF, nHI_OFF)

print("# r/Rp    T_ON[K]   T_OFF[K]  dT[K]     xHII_ON   xHII_OFF  "
      "nHI_ON      nHI_OFF")
for rr in [1.0, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0, 2.5, 3.0]:
    i = np.argmin(np.abs(rON - rr))
    print(f"{rON[i]:6.3f}  {TON[i]:8.1f}  {TOFF_i[i]:8.1f}  "
          f"{TON[i]-TOFF_i[i]:+8.1f}  {xON[i]:.3e}  {xOFF_i[i]:.3e}  "
          f"{nHI_ON[i]:.3e}  {nHI_OFF_i[i]:.3e}")

dT = TON - TOFF_i
print(f"\nmax |dT| (physical domain) = {np.max(np.abs(dT[m])):.1f} K "
      f"at r/Rp = {rON[m][np.argmax(np.abs(dT[m]))]:.3f}")
print(f"mean |dT| (physical domain) = {np.mean(np.abs(dT[m])):.1f} K")
dx = (xON - xOFF_i)
print(f"max |dx_HII| (physical domain) = {np.max(np.abs(dx[m])):.3e} "
      f"at r/Rp = {rON[m][np.argmax(np.abs(dx[m]))]:.3f}")

# ----- figure -----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# (1) temperature
ax[0].plot(rON[m],  TON[m],  'C0-',  label='ON (excited-H)')
ax[0].plot(rOFF,    TOFF,    'C3--', label='OFF')
ax[0].set_xlabel(r'$r/R_p$'); ax[0].set_ylabel(r'$T$ [K]')
ax[0].set_title('Temperature'); ax[0].legend(fontsize=9)
ax[0].set_xlim(1, min(6, rON.max()))

# (2) H ionization fraction
ax[1].plot(rON[m],  xON[m],  'C0-',  label='ON')
ax[1].plot(rOFF,    xOFF,    'C3--', label='OFF')
ax[1].set_xlabel(r'$r/R_p$'); ax[1].set_ylabel(r'$x_{\rm HII}=n_{\rm HII}/n_{\rm H}$')
ax[1].set_title('H ionization fraction'); ax[1].legend(fontsize=9)
ax[1].set_xlim(1, min(6, rON.max()))

# (3) neutral H density (log) + dT on twin axis
ax[2].plot(rON[m],  nHI_ON[m],  'C0-',  label=r'$n_{\rm HI}$ ON')
ax[2].plot(rOFF,    nHI_OFF,    'C3--', label=r'$n_{\rm HI}$ OFF')
ax[2].set_yscale('log'); ax[2].set_xlabel(r'$r/R_p$')
ax[2].set_ylabel(r'$n_{\rm HI}$ [cm$^{-3}$]')
ax[2].set_title('Neutral H density'); ax[2].legend(fontsize=9, loc='upper right')
ax[2].set_xlim(1, min(6, rON.max()))

plt.tight_layout()
plt.savefig('compare_on_off.png', dpi=130)
print("\nwrote compare_on_off.png")
