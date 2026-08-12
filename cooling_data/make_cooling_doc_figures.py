"""Generate the vector-PDF figures for docs/cooling_formulas.tex.

Self-contained: the fitted coefficients are inlined (they are the formulas
being documented; regenerate with fit_cooling_formulas.py /
fit_cno_formulas.py / fit_fs_saturation.py). The CHIANTI reference curves
are recomputed live via chianti_cooling.py.

Writes to ../docs/figures/:
  cool_resonance_fits.pdf   resonance lines + Fe II: table vs formula + resid
  cool_cno_compare.pdf      C/N/O: CHIANTI vs new fits vs AIOLOS + deviation
  cool_fs_saturation.pdf    C I/C II/N II/O I ground-term FS saturation
  cool_feii_decomp.pdf      Fe II coronal 4-term decomposition
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from chianti_cooling import cooling_effective

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "docs", "figures")
os.makedirs(OUT, exist_ok=True)
EV_ERG = 1.602176634e-12
K_B = 1.380649e-16
KB_CODE = 1.38e-16          # code-wide kb_erg
PREF = 8.629e-6

T = np.logspace(3.0, 5.0, 201)

# ------------------------------------------------------------------ tables
tab = np.loadtxt("metal_cooling_chianti.txt", skiprows=10)
Ttab = tab[:, 0]
assert np.allclose(T, Ttab)
table = dict(zip(["MgI_2853", "MgII_hk", "CaII_HK", "NaI_D",
                  "FeII_coronal", "FeII_boltz"], tab[:, 1:7].T))

# ------------------------------------------------- fitted formula library
def lam_line(T, g_l, dE_eV, a, b, T0):
    dE = dE_eV*EV_ERG
    ups = a + b*np.log(1.0 + T/T0)
    return PREF/(g_l*np.sqrt(T))*ups*dE*np.exp(-dE/(K_B*T))

def multiexp(A, Ti, T):
    A, Ti = np.asarray(A), np.asarray(Ti)
    return np.sum(A[:, None]*np.exp(-Ti[:, None]/T[None, :]),
                  axis=0)/np.sqrt(T)

FIT = {
  "MgI_2853": lam_line(T, 1, 4.3381, 0.355792, 15.1223, 94126.5),
  "MgII_hk":  lam_line(T, 2, 4.2766, 16.2305, 20.3448, 121896.0),
  "CaII_HK":  lam_line(T, 2, 3.1438, 14.7415, 19.492, 36218.3),
  "NaI_D":    lam_line(T, 2, 2.1037, 36.0744, 0.0, 1e30),
  "FeII_coronal": multiexp([1.9889e-18, 1.6793e-17, 6.5362e-16, 7.4638e-16],
                           [1350.1, 9630.3, 58939.3, 139593.3], T),
  "FeII_boltz": multiexp([1.5019e-18, 9.7419e-18, 5.2485e-16, 4.7387e-16],
                         [1175.7, 8678.5, 57594.7, 140620.2], T),
}
LBL = {"MgI_2853": r"Mg\,I $\lambda$2853", "MgII_hk": r"Mg\,II h\&k",
       "CaII_HK": r"Ca\,II H\&K", "NaI_D": r"Na\,I D",
       "FeII_coronal": r"Fe\,II coronal", "FeII_boltz": r"Fe\,II Boltz."}

# ---------------------------------------------- Fig 1: resonance + Fe II
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7.6), sharex=True,
                               gridspec_kw=dict(height_ratios=[2, 1]))
for nm in table:
    ln, = ax1.plot(T, table[nm], lw=2.3, alpha=0.42, label=LBL[nm])
    ax1.plot(T, FIT[nm], "--", lw=1.0, color=ln.get_color())
    ax2.plot(T, (FIT[nm]/table[nm] - 1.0)*100, lw=1.1, color=ln.get_color())
ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_ylim(1e-30, 1e-16)
ax1.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
ax1.legend(fontsize=8, ncol=2, loc="lower right")
ax1.set_title("CHIANTI tables (solid) vs closed-form fits (dashed)")
ax2.axhline(0, color="k", lw=0.6); ax2.axhspan(-1, 1, color="0.86")
ax2.set_ylim(-3.5, 3.5); ax2.set_xlabel(r"$T$ [K]")
ax2.set_ylabel(r"residual [\%]")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cool_resonance_fits.pdf"),
            bbox_inches="tight")
plt.close(fig)
print("wrote cool_resonance_fits.pdf")

# ---------------------------------------------------- C/N/O comparison
CHI = {nm: cooling_effective(el, st, T, pop="ground_term")
       for nm, (el, st) in
       dict(CI=("c", 1), CII=("c", 2), NI=("n", 1), NII=("n", 2),
            OI=("o", 1), OII=("o", 2)).items()}

def aiolos(name, T):
    if name == "CI":
        return 1.0e-24 + 3.1e-20*np.exp(-15162.0/T)*(1 + (T/2.0e4)**1.5)
    if name == "CII":
        return 1.5e-23 + 3.1e-20*np.exp(-45162.0/T)*(1 + (T/0.75e4)**1.5)
    if name == "OI":
        return 5.5e-24 + 1.1e-20*np.exp(-30162.0/T)*(1 + (T/0.75e4)**0.5)
    if name == "OII":
        return 5.1e-20*np.exp(-35162.0/T)*(1 + (T/0.75e4)**0.5)
    return np.zeros_like(T)

CNO = {
 "CI":  multiexp([1.10625037e-20, 1.39158295e-18, 6.96187591e-18,
                  6.59348477e-17, 3.01715156e-16],
                 [2351.38, 15172.3, 25411.0, 72609.1, 177264.0], T),
 "CII": multiexp([3.00421499e-20, 2.41049257e-20, 3.90434784e-17,
                  3.33191635e-16, 1.36009298e-15],
                 [61.3684, 7183.05, 64103.6, 124809.0, 246136.0], T),
 "NI":  multiexp([1.51401524e-18, 4.87732307e-18, 1.15970782e-17,
                  1.57486268e-16, 3.17477537e-16],
                 [29686.2, 34693.8, 50740.7, 139916.0, 254192.0], T),
 "NII": multiexp([2.40595483e-20, 1.24400479e-20, 8.04015794e-18,
                  1.54493674e-17, 2.88420949e-16, 6.94381304e-16],
                 [49.914, 6521.42, 23770.8, 64487.2, 157103.0, 295894.0], T),
 "OI":  multiexp([3.23621357e-22, 6.57165017e-20, 1.80264414e-18,
                  6.76183709e-18, 3.37196426e-17],
                 [847.937, 19058.9, 31568.4, 70528.3, 178827.0], T),
 "OII": multiexp([1.76816973e-17, 8.07102200e-18, 2.34447515e-16,
                  8.78914372e-16],
                 [44128.6, 59732.3, 179091.0, 325669.0], T),
}
LBL2 = dict(CI="C I", CII="C II", NI="N I", NII="N II", OI="O I", OII="O II")

fig, axes = plt.subplots(4, 3, figsize=(11.5, 10), sharex=True,
                         gridspec_kw=dict(height_ratios=[2, 1, 2, 1]))
pos = [(0, 0), (0, 1), (0, 2), (2, 0), (2, 1), (2, 2)]
for k, nm in enumerate(["CI", "CII", "NI", "NII", "OI", "OII"]):
    r, c = pos[k]
    axc, axr = axes[r, c], axes[r + 1, c]
    dat, fit, av = CHI[nm], CNO[nm], aiolos(nm, T)
    axc.plot(T, dat, lw=2.3, alpha=0.45, label="CHIANTI v11")
    axc.plot(T, fit, "k--", lw=1.0, label="new fit")
    if av.max() > 0:
        axc.plot(T, av, "r:", lw=1.5, label="AIOLOS")
    axc.set_xscale("log"); axc.set_yscale("log")
    axc.set_ylim(dat.max()*1e-9, dat.max()*5)
    axc.set_title(LBL2[nm], fontsize=11)
    axc.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
    axc.legend(fontsize=6.5, loc="lower right")
    axr.axhspan(-5, 5, color="0.88"); axr.axhline(0, color="k", lw=0.6)
    axr.plot(T, (fit/dat - 1)*100, "k-", lw=1.0)
    if av.max() > 0:
        axr.plot(T, (av/dat - 1)*100, "r:", lw=1.5)
    axr.set_ylim(-60, 60); axr.set_ylabel(r"dev. [\%]")
    if r == 2:
        axr.set_xlabel(r"$T$ [K]")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cool_cno_compare.pdf"), bbox_inches="tight")
plt.close(fig)
print("wrote cool_cno_compare.pdf")

# ------------------------------------------------- FS saturation figure
# The ground-term statistical-equilibrium coefficients, evaluated straight
# from the generator so the figure and the Fortran share one set of numbers.
from fit_fs_saturation import (FS_IONS, KH_TAB, REM_GUESS, ground_term_data,
                               fit_upsilon, fit_kH, above_term_coronal,
                               fit_multiexp, w_fs)

FSDATA = {}
for _nm, (_el, _ion, _nlev, _amu) in FS_IONS.items():
    _g, _Ek, _A, _ups = ground_term_data(_el, _ion, _nlev)
    _uf = {k: fit_upsilon(tr)[0] for k, tr in _ups.items()}
    _Tt, _kt = KH_TAB[_nm]
    _kf = {k: fit_kH(_Tt, v)[0:3:2] for k, v in _kt.items()}
    _Ar, _Tr = fit_multiexp(above_term_coronal(_el, _ion, _nlev, T),
                            REM_GUESS[_nm])
    FSDATA[_nm] = (_g, _Ek, _A, _uf, _kf, _Ar, _Tr, _nlev)


def fs_sat(name, Tv, ne, nHI):
    """Lambda_eff(T, ne, nHI) per (n_e n_ion), beta = 1."""
    g, Ek, A, uf, kf, Ar, Tr, nlev = FSDATA[name]
    x = np.log10(np.clip(Tv, 1e3, 1e5)/1e4)
    C = {}
    for (l, u), c in uf.items():
        C[(u, l)] = ne*PREF*10.0**np.polyval(c, x)/(g[u]*np.sqrt(Tv))
    for (u, l), (c, lim) in kf.items():
        uu = np.log10(np.clip(Tv, lim[0], lim[1])/1e3)
        C[(u-1, l-1)] = C[(u-1, l-1)] + nHI*10.0**np.polyval(c, uu)
    W = w_fs(g, Ek, A, C, np.ones((nlev, nlev)), Tv)
    return W/max(ne, 1e-30) + multiexp(Ar, Tr, Tv)


fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0))
CASES = [(1e2, 0.99, r"wind: $n_{\rm H}{=}10^2$, $x_e{=}0.99$"),
         (1e6, 0.5, r"$n_{\rm H}{=}10^6$, $x_e{=}0.5$"),
         (1e10, 1e-3, r"base: $n_{\rm H}{=}10^{10}$, $x_e{=}10^{-3}$"),
         (1e13, 1e-4, r"deep base: $n_{\rm H}{=}10^{13}$, $x_e{=}10^{-4}$")]
for ax, nm, lab in zip(axes.ravel(), ["CI", "CII", "NII", "OI"],
                       ["C I", "C II", "N II", "O I"]):
    ax.plot(T, CHI[nm], lw=2.6, alpha=0.4, label="coronal (CHIANTI, e-only)")
    for nH, xe, cl in CASES:
        ax.plot(T, fs_sat(nm, T, nH*xe, nH*(1 - xe)), lw=1.2, label=cl)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$T$ [K]")
    ax.set_ylabel(r"$\Lambda_{\rm eff}$ [erg cm$^3$ s$^{-1}$]")
    ax.set_title(lab); ax.legend(fontsize=7, loc="upper left")
    ax.set_ylim(1e-29, 1e-18)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cool_fs_saturation.pdf"), bbox_inches="tight")
plt.close(fig)
print("wrote cool_fs_saturation.pdf")

# ------------------------------------------- Fe II decomposition figure
A = np.array([1.9889e-18, 1.6793e-17, 6.5362e-16, 7.4638e-16])
Ti = np.array([1350.1, 9630.3, 58939.3, 139593.3])
glab = [r"$T_1{=}1350$ K (a$^6$D IR, 0.12 eV)",
        r"$T_2{=}9630$ K (metastable, 0.83 eV)",
        r"$T_3{=}58939$ K (UV, 5.1 eV)",
        r"$T_4{=}139593$ K (tail, 12 eV)"]
fig, ax = plt.subplots(figsize=(6.6, 4.6))
ax.plot(T, table["FeII_coronal"], lw=2.5, alpha=0.4, label="CHIANTI table")
ax.plot(T, multiexp(A, Ti, T), "k--", lw=1.1, label="4-term formula")
for i in range(4):
    ax.plot(T, A[i]*np.exp(-Ti[i]/T)/np.sqrt(T), ":", lw=1.3, label=glab[i])
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-24, 1e-17)
ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
ax.set_title("Fe II coronal: effective two-level decomposition")
ax.legend(fontsize=7.5, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "cool_feii_decomp.pdf"), bbox_inches="tight")
plt.close(fig)
print("wrote cool_feii_decomp.pdf")
