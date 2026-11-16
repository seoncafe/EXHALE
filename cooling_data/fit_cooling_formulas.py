"""Fit closed-form analytic formulas to the CHIANTI cooling tables.

Goal: express the hard-coded log10(Lambda) tables in Cool_coeff.f90 as
physically-motivated formulas.

Physics: every single-channel resonance line obeys exactly

    Lambda(T) = (8.629e-6 / (g_l sqrt(T))) * Upsilon(T) * dE * exp(-dE/kT)

so the only non-trivial factor is the effective collision strength
Upsilon(T), which is slowly varying (logarithmic in T for allowed lines,
Burgess & Tully 1992 type 1).  We fit

    Upsilon(T) = a + b * ln(1 + T/T0)

with dE left free as well (the tables were generated with the CHIANTI
.scups header dE, not the observed-wavelength energy; letting dE
float recovers it and keeps the formula to 3+1 parameters).  For Mg I and
Ca II that energy is the observed one; for Mg II it is the theoretical
4.27 eV, and the Mg II fit printed here is superseded by the observed-energy,
all-channel fit of magnesium_ii_line_cooling.py.

Fe II (coronal / Boltzmann-metastable) is a sum over hundreds of
transitions with dE from ~0.05 eV (a6D fine structure) to ~5 eV (UV), so a
single-exponential form cannot work; we fit a small sum of effective
two-level terms

    Lambda(T) = T^{-1/2} * sum_i  A_i * exp(-T_i / T)

(i = 1..4; each term is an effective transition group: IR fine-structure,
optical metastable, UV).  Nonlinear least squares on log Lambda.

Accuracy is reported over the full table (1e3-1e5 K) and over the range
where each coolant is non-negligible (Lambda > 1e-3 * peak).
Run:  python3 fit_cooling_formulas.py   (writes cooling_formula_fits.png)
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares

EV_ERG = 1.602176634e-12
K_B = 1.380649e-16
COLL_PREF = 8.629e-6

# ---------------------------------------------------------------- load table
tab = np.loadtxt("metal_cooling_chianti.txt", skiprows=10)
T = tab[:, 0]
cols = dict(zip(
    ["MgI_2853", "MgII_hk", "CaII_HK", "NaI_D", "FeII_coronal", "FeII_boltz"],
    tab[:, 1:7].T))

LINES = {  # g_l and initial dE guess [eV]
    "MgI_2853": dict(g_l=1.0, dE0=4.346),
    "MgII_hk":  dict(g_l=2.0, dE0=4.43),
    "CaII_HK":  dict(g_l=2.0, dE0=3.14),
    "NaI_D":    dict(g_l=2.0, dE0=2.104),
}

results = {}
errors = {}


def report(name, lam_fit, lam_tab):
    rel = np.abs(lam_fit / lam_tab - 1.0)
    mask = lam_tab > 1e-3 * lam_tab.max()
    errors[name] = (lam_fit, rel)
    print(f"    max rel. err (full 1e3-1e5 K)            : {rel.max()*100:6.2f} %")
    print(f"    max rel. err (relevant, T>{T[mask][0]:6.0f} K)    : "
          f"{rel[mask].max()*100:6.2f} %")
    print(f"    rms rel. err (relevant range)            : "
          f"{np.sqrt(np.mean(rel[mask]**2))*100:6.2f} %")


# ------------------------------------------ 1. resonance lines
for name, p in LINES.items():
    g_l = p["g_l"]
    lam = cols[name]

    if name == "NaI_D":
        # Na I was built from Van Regemorter (Upsilon = const by construction):
        # fit only (Upsilon, dE).
        def model(Tv, a, dE_eV):
            dE = dE_eV * EV_ERG
            return np.log10(COLL_PREF / (g_l * np.sqrt(Tv)) * a * dE
                            * np.exp(-dE / (K_B * Tv)))
        popt, _ = curve_fit(model, T, np.log10(lam),
                            p0=[36.0, p["dE0"]], maxfev=50000)
        a, dE = popt; b, T0 = 0.0, np.inf
    else:
        def model(Tv, a, b, T0, dE_eV):
            dE = dE_eV * EV_ERG
            ups = a + b * np.log(1.0 + Tv / T0)
            return np.log10(COLL_PREF / (g_l * np.sqrt(Tv)) * ups * dE
                            * np.exp(-dE / (K_B * Tv)))
        popt, _ = curve_fit(model, T, np.log10(lam),
                            p0=[1.0, 2.0, 2e4, p["dE0"]], maxfev=50000)
        a, b, T0, dE = popt
    lam_fit = 10.0 ** model(T, *popt)
    print(f"\n--- {name} ---")
    print(f"    Lambda = 8.629e-6/({g_l:.0f}*sqrt(T)) * Ups(T) * dE * exp(-dE/kT)")
    if b == 0.0:
        print(f"    dE = {dE:.4f} eV ;  Ups = {a:.6g}  (constant)")
    else:
        print(f"    dE = {dE:.4f} eV ;  Ups(T) = {a:.6g} + {b:.6g}*ln(1 + T/{T0:.6g})")
    report(name, lam_fit, lam)
    results[name] = dict(kind="bt", a=a, b=b, T0=T0, dE=dE, g_l=g_l)


# ------------------------------------------ 2. Fe II: sum of two-level terms
def fit_multiexp(name, nterm, Tguess):
    lam = cols[name]
    y = np.log(lam * np.sqrt(T))

    def resid(theta):
        A = np.exp(theta[:nterm])
        Ti = np.exp(theta[nterm:])
        return np.log(np.sum(A[:, None] * np.exp(-Ti[:, None] / T[None, :]),
                             axis=0)) - y

    A0 = np.full(nterm, np.log(lam[-1] * np.sqrt(T[-1]) / nterm))
    th0 = np.concatenate([A0, np.log(np.asarray(Tguess, dtype=float))])
    sol = least_squares(resid, th0, method="lm", max_nfev=40000)
    A = np.exp(sol.x[:nterm]); Ti = np.exp(sol.x[nterm:])
    order = np.argsort(Ti); A, Ti = A[order], Ti[order]
    lam_fit = np.sum(A[:, None] * np.exp(-Ti[:, None] / T[None, :]),
                     axis=0) / np.sqrt(T)
    print(f"\n--- {name} ---")
    print(f"    Lambda = T^-1/2 * sum_i A_i*exp(-T_i/T)")
    for ai, ti in zip(A, Ti):
        print(f"      A = {ai:.4e}   T_i = {ti:9.1f} K"
              f"   (dE_eff = {ti*K_B/EV_ERG:.3f} eV)")
    report(name, lam_fit, lam)
    return dict(kind="multiexp", A=A, Ti=Ti)


results["FeII_coronal"] = fit_multiexp(
    "FeII_coronal", 4, [550.0, 4000.0, 25000.0, 60000.0])
results["FeII_boltz"] = fit_multiexp(
    "FeII_boltz", 4, [550.0, 4000.0, 25000.0, 60000.0])

# ---------------------------------------------------------------- summary
print("\n================== summary of fitted formulas ==================")
print("Resonance lines: Lambda(T) = 8.629e-6/(g_l*sqrt(T))*Ups(T)*dE*exp(-dE/kT)")
for name, r in results.items():
    if r["kind"] == "bt":
        if r["b"] == 0.0:
            print(f"  {name:13s} g_l={r['g_l']:.0f}  dE={r['dE']:.4f} eV  "
                  f"Ups = {r['a']:.6g}  (constant)")
        else:
            print(f"  {name:13s} g_l={r['g_l']:.0f}  dE={r['dE']:.4f} eV  "
                  f"Ups = {r['a']:.6g} + {r['b']:.6g}*ln(1+T/{r['T0']:.6g})")
        # the digits Cool_coeff.f90 carries (cool_MgI_2853_coronal,
        # cool_CaII_HK_coronal): a, b, T0, dE [erg], dE/k [K]
        print(f"  {name:13s} Fortran: a b T0 dE_erg dE_over_k = "
              f"{r['a']:.8g} {r['b']:.8g} {r['T0']:.8g} "
              f"{r['dE']*EV_ERG:.8e} {r['dE']*EV_ERG/K_B:.3f}")
for name, r in results.items():
    if r["kind"] == "multiexp":
        terms = " + ".join(f"{a:.4e}*exp(-{t:.1f}/T)"
                           for a, t in zip(r["A"], r["Ti"]))
        print(f"  {name:13s} Lambda = T^-1/2 * [ {terms} ]")
print("  freefree_Z1   already analytic: 1.9095e-25 * Z^2 * (T/1e4)^0.55")

# ---------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAMES = list(cols.keys())
LABELS = {
    "MgI_2853": r"Mg\,I $\lambda$2853", "MgII_hk": r"Mg\,II h\&k",
    "CaII_HK": r"Ca\,II H\&K", "NaI_D": r"Na\,I D",
    "FeII_coronal": r"Fe\,II coronal", "FeII_boltz": r"Fe\,II Boltzmann",
}
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                               gridspec_kw=dict(height_ratios=[2, 1]))
for name in NAMES:
    ln, = ax1.plot(T, cols[name], lw=2.2, alpha=0.45, label=LABELS[name])
    ax1.plot(T, errors[name][0], "--", lw=1.0, color=ln.get_color())
    ax2.plot(T, errors[name][1] * 100, lw=1.2, color=ln.get_color())
ax1.set_xscale("log"); ax1.set_yscale("log")
ax1.set_ylim(1e-30, 1e-16)
ax1.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
ax1.legend(fontsize=8, ncol=2, loc="lower right")
ax1.set_title("CHIANTI tables (solid) vs analytic formulas (dashed)")
ax2.set_xlabel(r"$T$ [K]"); ax2.set_ylabel(r"$|$rel. error$|$ [\%]")
ax2.set_ylim(0, 5); ax2.axhline(1.0, color="k", lw=0.6, ls=":")
fig.tight_layout()
fig.savefig("cooling_formula_fits.png", dpi=150)
print("\nwrote cooling_formula_fits.png")
