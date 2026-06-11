#!/usr/bin/env python3
"""Phase-3a validation dump: H(n=2) Balmer proton source + photoelectric
heating vs radius, for comparison against Huang et al. (2023) Figs. 11/27
(proton source) and Figs. 10/26 (heating budget).

Surfaces NUMBERS + FIGURES only; validation is the user's call.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ----- normalization q0 (erg cm^-3 s^-1) from input.inp parameters -----
RJ   = 6.9911e9
kb   = 1.38e-16
mu   = 1.673e-24
n0   = 10.0**12.49
R0   = 1.766*RJ
T0   = 2358.0
v0   = np.sqrt(kb*T0/mu)
q0   = n0*mu*v0**3/R0
print(f"q0 = {q0:.4e} erg cm^-3 s^-1")

# ----- read Excited_H.txt -----
EH = np.loadtxt('output_on_clean/Excited_H.txt')
r      = EH[:, 0]
T      = EH[:, 1]
nHI    = EH[:, 2]
ne     = EH[:, 3]
Jlya   = EH[:, 4]
n2s    = EH[:, 5]
n2p    = EH[:, 6]
Sprot  = EH[:, 7]   # Balmer n=2 photoionization proton source [cm^-3 s^-1]
Hpe    = EH[:, 8]   # Balmer photoelectric heating [erg cm^-3 s^-1]
Hdx    = EH[:, 9]
tau    = EH[:, 10]

# ----- read Hydro_ioniz.txt (total heating) -----
# NOTE: write_output already multiplies col6/col7 by q0, so these columns are
# ALREADY in physical erg cm^-3 s^-1. Do NOT multiply by q0 again.
HY = np.loadtxt('output_on_clean/Hydro_ioniz.txt')
rH     = HY[:, 0]
heat_t = HY[:, 5]   # total radiative heating [erg cm^-3 s^-1], includes Hpe
cool_t = HY[:, 6]

# ----- read Ion_species.txt (nHII for recombination) -----
ION = np.loadtxt('output_on_clean/Ion_species.txt')
# Ion_species.txt: col1=r/Rp, col2=nHI, col3=nHII, col4=nHeI ...
nHII = ION[:, 2]    # col3 = HII number density [cm^-3]

# Case-B recombination rate alpha_B(T) [cm^3 s^-1] (Draine 2011 fit, same as module)
t4   = np.maximum(T, 1.0)/1.0e4
aB   = 2.54e-13*t4**(-0.8163 - 0.0208*np.log(t4))
recomb = aB*ne*nHII        # recombination sink [cm^-3 s^-1]

# physical-domain mask (drop ghosts)
m = (r >= 1.0) & (r <= r.max())

# ----- numeric summary -----
def at(rr):
    i = np.argmin(np.abs(r - rr))
    return i

print("\n# r/Rp   T[K]      nHI        n2s+n2p    Jlya       tau_Lya    "
      "Sproton    recomb     Sprot/recomb  Hpe        heat_tot   Hpe/heat")
for rr in [1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0]:
    i = at(rr)
    n2 = n2s[i] + n2p[i]
    sr = Sprot[i]/recomb[i] if recomb[i] > 0 else np.nan
    hf = Hpe[i]/heat_t[i] if heat_t[i] > 0 else np.nan
    print(f"{r[i]:6.3f}  {T[i]:8.1f}  {nHI[i]:.3e}  {n2:.3e}  {Jlya[i]:.3e}  "
          f"{tau[i]:.3e}  {Sprot[i]:.3e}  {recomb[i]:.3e}  {sr:10.3e}  "
          f"{Hpe[i]:.3e}  {heat_t[i]:.3e}  {hf:8.3e}")

# radius where tau_Lya ~ 1 (Ly-alpha photosphere)
itau = np.argmin(np.abs(tau[m] - 1.0))
print(f"\ntau_Lya ~ 1 (Ly-alpha photosphere) near r/Rp = {r[m][itau]:.3f}")
print(f"max Sprot/recomb (physical domain) = "
      f"{np.nanmax((Sprot/np.where(recomb>0,recomb,np.nan))[m]):.3e}")
print(f"max Hpe/heat_tot (physical domain) = "
      f"{np.nanmax((Hpe/np.where(heat_t>0,heat_t,np.nan))[m]):.3e}")

# ----- figures -----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# (1) proton source vs recombination
ax[0].plot(r[m], Sprot[m],  'C0-', label=r'$n=2$ Balmer photoion. $S_{\rm p}$')
ax[0].plot(r[m], recomb[m], 'C3--', label=r'recombination $\alpha_B n_e n_{\rm HII}$')
ax[0].set_yscale('log'); ax[0].set_xlabel(r'$r/R_p$')
ax[0].set_ylabel(r'rate [cm$^{-3}$ s$^{-1}$]')
ax[0].set_title('Proton source (cf. Huang Figs. 11/27)')
ax[0].legend(fontsize=8); ax[0].set_xlim(1, min(4, r.max()))

# (2) heating budget
ax[1].plot(r[m], heat_t[m], 'k-',  label=r'total heating')
ax[1].plot(r[m], Hpe[m],    'C2-', label=r'Balmer photoelec. $H_{\rm pe}$')
ax[1].plot(r[m], cool_t[m], 'C1:', label=r'total cooling')
ax[1].set_yscale('log'); ax[1].set_xlabel(r'$r/R_p$')
ax[1].set_ylabel(r'rate [erg cm$^{-3}$ s$^{-1}$]')
ax[1].set_title('Heating budget (cf. Huang Figs. 10/26)')
ax[1].legend(fontsize=8); ax[1].set_xlim(1, min(4, r.max()))

# (3) n=2 population + Jlya + tau
ax2 = ax[2].twinx()
ax[2].plot(r[m], (n2s+n2p)[m], 'C4-', label=r'$n_2 = n_{2s}+n_{2p}$')
ax[2].set_yscale('log'); ax[2].set_ylabel(r'$n_2$ [cm$^{-3}$]', color='C4')
ax2.plot(r[m], tau[m], 'C5--', label=r'$\tau_{{\rm Ly}\alpha}$')
ax2.set_yscale('log'); ax2.set_ylabel(r'$\tau_{{\rm Ly}\alpha}$', color='C5')
ax2.axhline(1.0, color='gray', lw=0.6, ls=':')
ax[2].set_xlabel(r'$r/R_p$'); ax[2].set_title(r'$n=2$ population and Ly$\alpha$ depth')
ax[2].set_xlim(1, min(4, r.max()))

plt.tight_layout()
plt.savefig('excited_H_validation.png', dpi=130)
print("\nwrote excited_H_validation.png")
