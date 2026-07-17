#!/usr/bin/env python3
"""
Compare the recombination rate coefficients used in EXHALE for
H I (H+ -> H), He II (He++ -> He+), He I total (He+ -> He), and the
He I singlet/triplet channels (1^1S, metastable 2^3S) against other
published fits, over the temperature range relevant to escaping
thermospheres (3e3 - 2e4 K).

EXHALE sources (src/modules/radiation/Cool_coeff.f90):
  - Since 2026-07-17 the EXHALE default H/He case-B coefficients are
    alpha_A(Badnell 2023) - alpha_1(Mao & Kaastra 2016), with Badnell DR
    added for He II (functions alphaB_*_new / rr_badnell / rr_mao). The H+/
    He++ curves plotted here are the legacy Hui & Gnedin (1997) fits, now
    the `legacy_hhe_rates` option; the new default differs by ~1% for H+ and
    He++ and is ~7% higher for the He+ total (Badnell dielectronic term), so
    the plotted curves still represent the H+/He++ default. Port alphaB_*_new
    to overlay the exact new default.
  - H+, He++ -> case-B fits of Hui & Gnedin (1997) [legacy_hhe_rates]
  - He+ -> He total -> case-B of Hui & Gnedin (1997) [legacy_hhe_rates]
  - He+ -> He(1^1S) and He+ -> He(2^3S) -> Benjamin, Skillman & Smits
    (1999) effective (cascade-summed) coefficients, as parametrized by
    Oklopcic & Hirata (2018); used when the He 2^3S network is active.

Reference fits:
  - Draine (2011, "Physics of the ISM") case-A/B hydrogen fits
  - Verner & Ferland (1996) / Storey & Hummer (1995) total (case-A)
    coefficients in the form adopted by the recent multi-species escape
    model of Taylor et al. (2025, arXiv:2506.08232)

All coefficients in cm^3 s^-1.  Outputs (vector PDF) -> docs/figures/:
  rec_H.pdf    -- H+ -> H : EXHALE vs Draine vs Taylor2025 (case A vs B)
  rec_He.pdf   -- He recombination: He++->He+ and He+->He (+1S/3S split)

Author: Kwang-Il Seon
"""
import os
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": False, "font.family": "serif",
    "font.serif": ["Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 12,
    "axes.titlesize": 12.5, "axes.labelsize": 12,
    "legend.fontsize": 9.5, "lines.linewidth": 2.0,
})
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUTDIR, exist_ok=True)


# ============================ EXHALE (Hui & Gnedin 1997 case B) ===========
def aB_HII_HG97(T):
    xl = 2.0 * 157807.0 / T
    return 2.753e-14 * xl ** 1.5 / (1.0 + (xl / 2.740) ** 0.407) ** 2.242


def aB_HeIII_HG97(T):           # He++ -> He+ (hydrogenic, Z=2)
    xl = 2.0 * 631515.0 / T
    return 2.0 * 2.753e-14 * xl ** 1.5 / (1.0 + (xl / 2.740) ** 0.407) ** 2.242


def aB_HeII_HG97(T):            # He+ -> He total, case B
    xl = 2.0 * 285335.0 / T
    return 1.26e-14 * xl ** 0.75


# EXHALE He+ -> He(1^1S) / He(2^3S):  Benjamin+1999 / Oklopcic & Hirata 2018
def a_HeI_1S(T):
    return 1.54e-13 * (T / 1.0e4) ** (-0.486)


def a_HeI_3S(T):
    return 2.10e-13 * (T / 1.0e4) ** (-0.778)


# ============================ Reference fits ==============================
def aB_H_Draine(T):             # Draine 2011, case B
    t4 = T / 1.0e4
    return 2.54e-13 * t4 ** (-0.8163 - 0.0208 * np.log(t4))


def aA_H_Draine(T):             # Draine 2011, case A (total)
    t4 = T / 1.0e4
    return 4.13e-13 * t4 ** (-0.7131 - 0.0115 * np.log(t4))


def aA_H_Taylor25(T):           # Taylor+2025 (Storey & Hummer 1995), total
    return 4.0e-12 * (300.0 / T) ** 0.64


def aA_HeIII_Taylor25(T):       # Taylor+2025 (Verner & Ferland 1996), total
    return 1.92e-11 * (300.0 / T) ** 0.61


def aB_HeIII_Draine(T):         # hydrogenic scaling of Draine case B (Z=2)
    return 2.0 * aB_H_Draine(T / 4.0)


# ============================ Figure 1: H+ -> H ===========================
def fig_H():
    T = np.logspace(np.log10(3e3), np.log10(2e4), 600)
    ex = aB_HII_HG97(T)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 6.2), sharex=True,
                                   gridspec_kw={"height_ratios": [2.4, 1.0]})
    ax1.loglog(T, ex, color="#d62728",
               label="EXHALE legacy: Hui & Gnedin 1997 (case B)")
    ax1.loglog(T, aB_H_Draine(T), color="k", ls="--", label="Draine 2011 (case B)")
    ax1.loglog(T, aA_H_Draine(T), color="0.45", ls=":", label="Draine 2011 (case A)")
    ax1.loglog(T, aA_H_Taylor25(T), color="#1f77b4", ls="-.",
               label="Taylor+2025 / S&H95 (total, case A)")
    ax1.set_ylabel(r"$\alpha_{\rm H\,II}$ [cm$^3$ s$^{-1}$]")
    ax1.set_title(r"H$^+\!\to$ H recombination")
    ax1.grid(True, which="both", alpha=0.25); ax1.legend(loc="upper right")
    ax2.semilogx(T, aB_H_Draine(T) / ex, color="k", ls="--", label="Draine B / EXHALE")
    ax2.semilogx(T, aA_H_Taylor25(T) / ex, color="#1f77b4", ls="-.",
                 label="Taylor25 A / EXHALE")
    ax2.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax2.set_ylabel("ratio to EXHALE"); ax2.set_xlabel(r"$T$ [K]")
    ax2.set_xlim(3e3, 2e4); ax2.set_ylim(0.8, 1.9)
    ax2.grid(True, which="both", alpha=0.25); ax2.legend(loc="upper right", ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "rec_H.pdf")); plt.close(fig)


# ============================ Figure 2: He recombination ==================
def fig_He():
    T = np.logspace(np.log10(3e3), np.log10(2e4), 600)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.8))

    # He++ -> He+
    axL.loglog(T, aB_HeIII_HG97(T), color="#d62728",
               label="EXHALE legacy: Hui & Gnedin 1997 (case B)")
    axL.loglog(T, aB_HeIII_Draine(T), color="k", ls="--",
               label="Draine 2011 hydrogenic (case B)")
    axL.loglog(T, aA_HeIII_Taylor25(T), color="#1f77b4", ls="-.",
               label="Taylor+2025 / V&F96 (total, case A)")
    axL.set_ylabel(r"$\alpha_{\rm He\,III}$ [cm$^3$ s$^{-1}$]")
    axL.set_xlabel(r"$T$ [K]")
    axL.set_title(r"He$^{++}\!\to$ He$^+$ recombination")
    axL.set_xlim(3e3, 2e4); axL.grid(True, which="both", alpha=0.25)
    axL.legend(loc="upper right")

    # He+ -> He : total case B, plus 1S / 3S split
    axR.loglog(T, aB_HeII_HG97(T), color="0.45", ls=":",
               label="Hui & Gnedin 1997 total (case B)")
    axR.loglog(T, a_HeI_1S(T) + a_HeI_3S(T), color="#d62728",
               label=r"Benjamin+1999: $\alpha_1+\alpha_3$ (EXHALE)")
    axR.loglog(T, a_HeI_1S(T), color="#2ca02c", ls="--",
               label=r"$\alpha_1$ (singlet $1^1S$)")
    axR.loglog(T, a_HeI_3S(T), color="#1f77b4", ls="-.",
               label=r"$\alpha_3$ (triplet $2^3S$)")
    axR.set_ylabel(r"$\alpha_{\rm He\,II}$ [cm$^3$ s$^{-1}$]")
    axR.set_xlabel(r"$T$ [K]")
    axR.set_title(r"He$^+\!\to$ He recombination (+ $1^1S/2^3S$ split)")
    axR.set_xlim(3e3, 2e4); axR.grid(True, which="both", alpha=0.25)
    axR.legend(loc="upper right")

    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "rec_He.pdf")); plt.close(fig)


def print_tables():
    Ts = [5e3, 1e4, 1.5e4, 2e4]
    print("=== H+ -> H  (1e-13 cm^3/s) ===")
    print(f"{'T[K]':>7} {'EXHALE(B)':>10} {'DraineB':>9} {'DraineA':>9} {'Tayl25A':>9}")
    for T in Ts:
        print(f"{T:>7.0f} {aB_HII_HG97(T)*1e13:>10.3f} {aB_H_Draine(T)*1e13:>9.3f} "
              f"{aA_H_Draine(T)*1e13:>9.3f} {aA_H_Taylor25(T)*1e13:>9.3f}")
    print("\n=== He++ -> He+  (1e-13 cm^3/s) ===")
    print(f"{'T[K]':>7} {'EXHALE(B)':>10} {'DraineB':>9} {'Tayl25A':>9}")
    for T in Ts:
        print(f"{T:>7.0f} {aB_HeIII_HG97(T)*1e13:>10.3f} {aB_HeIII_Draine(T)*1e13:>9.3f} "
              f"{aA_HeIII_Taylor25(T)*1e13:>9.3f}")
    print("\n=== He+ -> He  (1e-13 cm^3/s):  total(HG97,B), a1(1S), a3(3S), a1+a3 ===")
    print(f"{'T[K]':>7} {'HG97 tot':>9} {'a1':>8} {'a3':>8} {'a1+a3':>8} {'a3/(a1+a3)':>11}")
    for T in Ts:
        a1, a3 = a_HeI_1S(T), a_HeI_3S(T)
        print(f"{T:>7.0f} {aB_HeII_HG97(T)*1e13:>9.3f} {a1*1e13:>8.3f} {a3*1e13:>8.3f} "
              f"{(a1+a3)*1e13:>8.3f} {a3/(a1+a3):>11.3f}")


if __name__ == "__main__":
    fig_H()
    fig_He()
    print_tables()
    print(f"\nFigures written to {OUTDIR}")
