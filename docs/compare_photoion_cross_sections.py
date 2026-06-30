#!/usr/bin/env python3
"""
Compare the photoionization cross sections used in EXHALE for
H I, He I, He II, and the metastable He I 2^3 S triplet, and benchmark
two of them against external references:

  (a) EXHALE sigma_HeI (2-term ATES fit) vs the Verner et al. (1996)
      single-shell analytic fit (VFKY96), over 24.6-500 eV.
  (b) EXHALE sigma_HeI23S (broken power-law fit to Norcross 1971) vs the
      p-winds tabulation of the same Norcross (1971) data, over 4.77-60 eV.

All cross sections are returned in units of 1e-18 cm^2 = 1 Mb, matching the
convention of src/modules/functions/cross_sec.f90.  Photon energy E is in eV.

Outputs (vector PDF) -> docs/figures/:
  xsec_overview.pdf            -- all four species on one log-log plot
  xsec_HeI_vs_Verner.pdf       -- He I: EXHALE vs Verner (+ ratio)
  xsec_HeI23S_vs_pwinds.pdf    -- He 2^3 S: EXHALE vs p-winds/Norcross (+ ratio)

Author: Kwang-Il Seon
"""
import os
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "lines.linewidth": 2.0,
})

# ---- EXHALE physical constants (src/modules/init/parameters.f90) ----
HP_EV   = 4.1357e-15        # Planck constant [eV s]
C_LIGHT = 2.99792458e10     # speed of light  [cm/s]
PI      = 3.1415926536
HC_EVA  = HP_EV * C_LIGHT * 1e8   # hc [eV*Angstrom] = 12398.47

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUTDIR, exist_ok=True)


# =====================================================================
#  EXHALE cross sections (ported verbatim from cross_sec.f90), in Mb
# =====================================================================
def sigma_hydrogenic(E, Z):
    """Hydrogenic atoms (H I: Z=1, He II: Z=2).  cross_sec.f90:sigma."""
    E = np.asarray(E, dtype=float)
    E0 = 13.6 * Z * Z
    out = np.where(
        E > E0,
        # E > E_th branch
        np.where(
            E > E0,
            6.3 / (Z * Z)
            * (E0 / np.where(E > 0, E, np.nan)) ** 4.0
            * np.exp(4.0 - 4.0 * np.arctan(np.sqrt(np.abs(E / E0 - 1.0)))
                     / np.sqrt(np.abs(E / E0 - 1.0)))
            / (1.0 - np.exp(-2.0 * PI / np.sqrt(np.abs(E / E0 - 1.0)))),
            0.0),
        6.3 / (Z * Z))           # E <= E_th: constant value (then cut below)
    cut = 0.99999 * E0
    out = np.where(E < cut, 0.0, out)
    return out


def sigma_HeI(E):
    """He I ground state, 2-term ATES fit.  cross_sec.f90:sigma_HeI."""
    E = np.asarray(E, dtype=float)
    eth = 24.6 * 0.999
    s = 0.6935 / ((E * 1.0e-2) ** 1.82 + (E * 1.0e-2) ** 3.23)
    return np.where(E >= eth, s, 0.0)


def sigma_HeI23S(E):
    """He I 2^3 S metastable, broken power-law fit to Norcross (1971).
    cross_sec.f90:sigma_HeI23S."""
    E = np.asarray(E, dtype=float)
    x1 = np.log10(HC_EVA / 2593.01) * 0.9999
    x2 = np.log10(HC_EVA / 1655.63)
    x3 = np.log10(HC_EVA / 357.340)
    x4 = np.log10(HC_EVA / 271.940)
    x5 = np.log10(HC_EVA / 209.490) * 1.0001
    a1, a2, a3 = -0.8134, -1.772, -3.039
    c1 = 1.240
    c2 = c1 + x2 * (a1 - a2)
    c3 = 5.470
    y3 = a2 * x3 + c2
    y4 = a3 * x4 + c3
    m = (y4 - y3) / (x4 - x3)

    logE = np.log10(np.where(E > 0, E, np.nan))
    logs = np.full_like(logE, np.nan)
    logs = np.where(logE <= x2, a1 * logE + c1, logs)
    logs = np.where((logE > x2) & (logE <= x3), a2 * logE + c2, logs)
    logs = np.where((logE > x3) & (logE < x4), m * (logE - x3) + y3, logs)
    logs = np.where(logE >= x4, a3 * logE + c3, logs)
    s = 10.0 ** logs
    s = np.where((logE < x1) | (logE > x5), 0.0, s)
    return s


# =====================================================================
#  Verner et al. (1996) single-shell analytic fit (VFKY96), in Mb
#  sigma_0 quoted in Mb directly.
# =====================================================================
def sigma_VFKY96(E, Eth, E0, s0, ya, P, yw, y0, y1):
    E = np.asarray(E, dtype=float)
    x = E / E0 - y0
    z = np.sqrt(x * x + y1 * y1)
    Q = 5.5 - 0.5 * P
    Fy = ((x - 1.0) ** 2 + yw ** 2) * z ** (-Q) * (1.0 + np.sqrt(z / ya)) ** (-P)
    return np.where(E >= Eth, s0 * Fy, 0.0)


# Verner+1996 Table 1 ground-shell parameters
#   [Eth, E0, sigma0(Mb), ya, P, yw, y0, y1]
VERNER = {
    "HI":   (13.60, 0.4298, 5.475e4, 32.88, 2.963, 0.0,   0.0,    0.0),
    "HeI":  (24.59, 13.61,  949.2,   1.469, 3.188, 2.039, 0.4434, 2.136),
    "HeII": (54.42, 1.720,  1.369e4, 32.88, 2.963, 0.0,   0.0,    0.0),
}


def sigma_verner(species, E):
    return sigma_VFKY96(E, *VERNER[species])


# =====================================================================
#  p-winds tabulation of Norcross (1971) for He 2^3 S
#  data_array = [wavelength (A), differential oscillator strength]
#  sigma[cm^2] = 8.0670e-18 * (df/deps)  ->  sigma[Mb] = 8.0670 * (df/deps)
# =====================================================================
PWINDS_NORCROSS = np.array([
    [2593.01, 0.605],  [2528.27, 0.589],  [2275.74, 0.537],
    [2023.15, 0.501],  [1655.63, 0.435],  [1214.41, 0.247],
    [958.87,  0.1572], [792.18,  0.1138], [674.86,  0.0780],
    [587.81,  0.0620], [520.65,  0.0557], [467.27,  0.0461],
    [423.81,  0.0358], [387.75,  0.0310], [357.34,  0.0325],
    [331.36,  0.0520], [271.94,  0.343],  [271.21,  0.338],
    [256.70,  0.274],  [243.01,  0.231],  [230.71,  0.200],
    [219.59,  0.1750], [209.49,  0.1537],
])
PW_LAMBDA = PWINDS_NORCROSS[:, 0]
PW_SIGMA_MB = 8.0670 * PWINDS_NORCROSS[:, 1]      # Mb
PW_ENERGY = HC_EVA / PW_LAMBDA                    # eV (4.78 ... 59.2 eV)


def sigma_pwinds(E):
    """Log-log interpolation of the p-winds/Norcross table (Mb)."""
    order = np.argsort(PW_ENERGY)
    xe = PW_ENERGY[order]
    ys = PW_SIGMA_MB[order]
    E = np.asarray(E, dtype=float)
    logs = np.interp(np.log10(E), np.log10(xe), np.log10(ys),
                     left=np.nan, right=np.nan)
    return 10.0 ** logs


# =====================================================================
#  TOPbase / Opacity Project tabulation of He I 1s2s 3S (NZ=2, NE=2,
#  iSLP=300, iLV=1).  Native units: photon energy in Rydberg, sigma in Mb.
#  Retrieved from the CDS TOPbase server (com=dt, ent=p); see docstring of
#  load_topbase().  Includes the full autoionizing-resonance structure.
# =====================================================================
RYD_EV = 13.605693  # 1 Rydberg in eV

TOPBASE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "data", "topbase_HeI_2_3S_photo.dat")


def load_topbase():
    """Return (E_eV, sigma_Mb) for the TOPbase He I 2^3 S cross section.
    Column 1 of the file is the photon energy in Rydberg, converted to eV."""
    arr = np.loadtxt(TOPBASE_FILE)
    E_eV = arr[:, 0] * RYD_EV
    sigma = arr[:, 1]
    return E_eV, sigma


def sigma_topbase_smooth(E, E_tb, s_tb, e_max_smooth=23.0):
    """Log-log interpolation of TOPbase BELOW the first autoionizing
    resonances (E <= e_max_smooth eV), where the cross section is smooth
    and directly comparable to the Norcross/EXHALE background."""
    m = E_tb <= e_max_smooth
    xe, ys = E_tb[m], s_tb[m]
    E = np.asarray(E, dtype=float)
    logs = np.interp(np.log10(E), np.log10(xe), np.log10(ys),
                     left=np.nan, right=np.nan)
    return 10.0 ** logs


# =====================================================================
#  Figure 1 -- overview of all four species (EXHALE-adopted curves)
# =====================================================================
def fig_overview():
    E = np.logspace(np.log10(1.0), np.log10(1000.0), 2000)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.loglog(E, np.clip(sigma_HeI23S(E), 1e-6, None),
              color="#d62728", label=r"He I $2^3$S (broken PL, Norcross 71)")
    ax.loglog(E, np.clip(sigma_hydrogenic(E, 1.0), 1e-6, None),
              color="#1f77b4", label=r"H I (hydrogenic)")
    ax.loglog(E, np.clip(sigma_HeI(E), 1e-6, None),
              color="#2ca02c", label=r"He I (2-term fit)")
    ax.loglog(E, np.clip(sigma_hydrogenic(E, 2.0), 1e-6, None),
              color="#9467bd", label=r"He II (hydrogenic)")
    for Eth, txt, col in [(4.781, r"He I $2^3$S", "#d62728"),
                          (13.6,  "H I", "#1f77b4"),
                          (24.6,  "He I", "#2ca02c"),
                          (54.4,  "He II", "#9467bd")]:
        ax.axvline(Eth, color=col, ls=":", lw=1.0, alpha=0.6)
    ax.set_xlabel(r"photon energy $E$ [eV]")
    ax.set_ylabel(r"$\sigma$ [Mb $=10^{-18}\,$cm$^2$]")
    ax.set_title("EXHALE photoionization cross sections")
    ax.set_xlim(1, 1000)
    ax.set_ylim(1e-3, 2e1)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="lower left", frameon=True)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "xsec_overview.pdf"))
    plt.close(fig)


# =====================================================================
#  Figure 2 -- He I: EXHALE 2-term fit vs Verner (24.6-500 eV)
# =====================================================================
def fig_HeI_vs_verner():
    E = np.logspace(np.log10(24.6), np.log10(500.0), 1500)
    s_ex = sigma_HeI(E)
    s_ve = sigma_verner("HeI", E)
    ratio = s_ex / s_ve

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7.0, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0]})
    ax1.loglog(E, s_ex, color="#2ca02c", label=r"EXHALE sigma_HeI (2-term)")
    ax1.loglog(E, s_ve, color="k", ls="--", label="Verner+1996 (VFKY96)")
    ax1.set_ylabel(r"$\sigma_{\rm He\,I}$ [Mb]")
    ax1.set_title(r"He I ground-state photoionization: EXHALE vs Verner+1996")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="upper right")

    ax2.semilogx(E, ratio, color="#2ca02c")
    ax2.axhline(1.0, color="k", lw=0.8, ls=":")
    ax2.set_ylabel("EXHALE / Verner")
    ax2.set_xlabel(r"photon energy $E$ [eV]")
    ax2.set_xlim(24.6, 500)
    ax2.set_ylim(0.6, 1.4)
    ax2.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "xsec_HeI_vs_Verner.pdf"))
    plt.close(fig)
    return E, s_ex, s_ve, ratio


# =====================================================================
#  Figure 3 -- He 2^3 S: EXHALE broken-PL vs p-winds/Norcross (4.77-60 eV)
# =====================================================================
def fig_HeI23S_vs_pwinds():
    E = np.logspace(np.log10(4.781), np.log10(59.0), 1500)
    s_ex = sigma_HeI23S(E)
    s_pw = sigma_pwinds(E)
    ratio = s_ex / s_pw

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7.0, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0]})
    ax1.loglog(E, s_ex, color="#d62728",
               label=r"EXHALE sigma_HeI23S (broken PL)")
    ax1.loglog(E, s_pw, color="#1f77b4", ls="--",
               label="p-winds log-log interp.")
    ax1.loglog(PW_ENERGY, PW_SIGMA_MB, "o", ms=4.5, color="k",
               label="Norcross 1971 (p-winds table)")
    ax1.set_ylabel(r"$\sigma_{{\rm He}\,2^3S}$ [Mb]")
    ax1.set_title(r"He I $2^3$S photoionization: EXHALE vs p-winds (Norcross 1971)")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="lower left")

    ax2.semilogx(E, ratio, color="#d62728")
    ax2.axhline(1.0, color="k", lw=0.8, ls=":")
    ax2.set_ylabel("EXHALE / p-winds")
    ax2.set_xlabel(r"photon energy $E$ [eV]")
    ax2.set_xlim(4.781, 59.0)
    ax2.set_ylim(0.5, 1.7)
    ax2.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "xsec_HeI23S_vs_pwinds.pdf"))
    plt.close(fig)
    return E, s_ex, s_pw, ratio


# =====================================================================
#  Figure 4 -- He 2^3 S: EXHALE + p-winds/Norcross + TOPbase (Opacity Proj.)
#  Two panels: (a) full range showing TOPbase autoionizing resonances;
#              (b) smooth near-threshold region (4.6-25 eV).
# =====================================================================
def fig_HeI23S_topbase():
    E_tb, s_tb = load_topbase()
    Efit = np.logspace(np.log10(4.781), np.log10(60.0), 1500)
    s_ex = sigma_HeI23S(Efit)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.8))

    # (a) full range with resonances
    axA.loglog(E_tb, np.clip(s_tb, 1e-3, None), color="0.45", lw=0.8,
               label="TOPbase / OP (full, with resonances)")
    axA.loglog(Efit, s_ex, color="#d62728", lw=2.0,
               label="EXHALE sigma_HeI23S (broken PL)")
    axA.loglog(PW_ENERGY, PW_SIGMA_MB, "o", ms=4.0, color="#1f77b4",
               label="Norcross 1971 (p-winds)")
    axA.set_xlabel(r"photon energy $E$ [eV]")
    axA.set_ylabel(r"$\sigma_{{\rm He}\,2^3S}$ [Mb]")
    axA.set_title("(a) Full range: TOPbase resonance structure")
    axA.set_xlim(4.0, 330)
    axA.set_ylim(1e-2, 5e3)
    axA.grid(True, which="both", alpha=0.25)
    axA.legend(loc="upper right", fontsize=8.5)

    # (b) smooth near-threshold region
    m = E_tb <= 25.0
    axB.loglog(E_tb[m], s_tb[m], color="0.45", lw=1.4, label="TOPbase / OP")
    axB.loglog(Efit[Efit <= 25], s_ex[Efit <= 25], color="#d62728", lw=2.0,
               label="EXHALE broken PL")
    pm = PW_ENERGY <= 25.0
    axB.loglog(PW_ENERGY[pm], PW_SIGMA_MB[pm], "o", ms=5.0, color="#1f77b4",
               label="Norcross 1971 (p-winds)")
    axB.set_xlabel(r"photon energy $E$ [eV]")
    axB.set_ylabel(r"$\sigma_{{\rm He}\,2^3S}$ [Mb]")
    axB.set_title("(b) Smooth near-threshold region")
    axB.set_xlim(4.4, 25)
    axB.set_ylim(0.4, 8)
    axB.grid(True, which="both", alpha=0.25)
    axB.legend(loc="lower left", fontsize=9)

    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "xsec_HeI23S_topbase.pdf"))
    plt.close(fig)
    return E_tb, s_tb


def print_topbase_table(E_tb, s_tb):
    print("\n=== He 2^3 S smooth region: EXHALE vs Norcross/p-winds vs TOPbase (Mb) ===")
    print(f"{'E[eV]':>8} {'EXHALE':>10} {'p-winds':>10} {'TOPbase':>10} "
          f"{'EX/TB':>8} {'PW/TB':>8}")
    for Eq in [5.0, 6, 8, 10, 12, 16, 20]:
        a = float(sigma_HeI23S(Eq))
        b = float(sigma_pwinds(Eq))
        c = float(sigma_topbase_smooth(Eq, E_tb, s_tb))
        print(f"{Eq:>8.1f} {a:>10.4f} {b:>10.4f} {c:>10.4f} "
              f"{a/c:>8.3f} {b/c:>8.3f}")
    # smooth-region ratio statistics
    Eg = np.logspace(np.log10(5.0), np.log10(20.0), 400)
    ex = sigma_HeI23S(Eg); tb = sigma_topbase_smooth(Eg, E_tb, s_tb)
    r = ex / tb
    print(f"EXHALE/TOPbase over 5-20 eV: min={r.min():.3f}, "
          f"max={r.max():.3f}, median={np.median(r):.3f}")
    print(f"TOPbase threshold (first point): E={E_tb.min():.3f} eV, "
          f"sigma={s_tb[np.argmin(E_tb)]:.3f} Mb; "
          f"resonance peak sigma_max={s_tb.max():.1f} Mb")


# =====================================================================
#  Impact of the >59 eV truncation on the He 2^3 S photoionization RATE.
#  Rate ~ INT sigma(E) N(E) dE, with photon-number weight N(E)=J_inc/E and
#  J_inc ~ E^PLind in the EXHALE power-law SED (EUV band [e_low,e_mid]
#  normalized to L_EUV, X-ray band [e_mid,e_top] to L_X).  We report the
#  relative increase of the rate when sigma is extended above 59.18 eV
#  with the TOPbase tail.
# =====================================================================
E_LOW_SED, E_MID_SED, E_TOP_SED = 4.80, 124.0, 1240.0  # parameters.f90


def rate_increase_above_cutoff(PLind, Lrapp, E_tb, s_tb):
    """Percent increase of the He 2^3 S photoionization rate if sigma is
    extended above 59.18 eV with the TOPbase tail, for an EXHALE power-law
    SED of index PLind and X-ray/EUV luminosity ratio Lrapp."""
    el, em, et = E_LOW_SED, E_MID_SED, E_TOP_SED
    P1 = PLind + 1.0
    if PLind != -1.0:
        jeuv = P1 / (em ** P1 - el ** P1)
        jx = P1 / (et ** P1 - em ** P1)
    else:
        jeuv = 1.0 / np.log(em / el)
        jx = 1.0 / np.log(et / em)
    jeuv *= 1.0 / (1.0 + Lrapp)
    jx *= Lrapp / (1.0 + Lrapp)

    def Nw(E):  # photon-number weight ~ J_inc/E
        return np.where(E < em, jeuv * E ** PLind, jx * E ** PLind) / E

    def sig_tb(E):
        return 10.0 ** np.interp(np.log10(E), np.log10(E_tb),
                                 np.log10(np.clip(s_tb, 1e-6, None)))

    Ea = np.linspace(4.80, 59.18, 6000)        # current EXHALE coverage
    Rc = np.trapz(sigma_HeI23S(Ea) * Nw(Ea), Ea)
    Eb = np.linspace(59.18, 323.0, 8000)       # TOPbase tail
    dR = np.trapz(sig_tb(Eb) * Nw(Eb), Eb)
    return dR / Rc * 100.0


def print_rate_impact(E_tb, s_tb):
    print("\n=== Rate impact of extending sigma_HeI23S above 59 eV (TOPbase tail) ===")
    print("relative increase [%] of the He 2^3 S photoionization rate")
    print(f"{'PLind':>7} {'Lx/Leuv=0.1':>13} {'0.3':>9} {'1.0':>9}")
    for P in [-1.0, -1.3, -1.6, -2.0]:
        row = [rate_increase_above_cutoff(P, L, E_tb, s_tb)
               for L in (0.1, 0.3, 1.0)]
        print(f"{P:>7.1f} {row[0]:>12.2f}% {row[1]:>8.2f}% {row[2]:>8.2f}%")


# =====================================================================
#  Impact of the autoionizing resonances (40-55 eV) on the He 2^3 S RATE.
#  WARNING: the raw TOPbase table must NOT be integrated directly --- the
#  tallest resonances (sigma up to 3223 Mb) are sampled by only 1-2 grid
#  points, so trapezoidal quadrature assigns them an essentially arbitrary
#  area set by the mesh spacing, not by the physical (Lorentzian) width.
#  We therefore (i) integrate EXHALE and TOPbase on the SAME native grid
#  and (ii) test convergence by clipping the unresolved peaks; a robust
#  resonance treatment needs a properly resonance-averaged cross section.
# =====================================================================
def resonance_impact(E_tb, s_tb, PLind=-1.3):
    """Return the relative change of the He 2^3 S rate from using the full
    TOPbase (resonances) instead of the EXHALE broken-PL bump, as a
    function of an upper clip on sigma (to expose the spike sensitivity).
    Photon weight N(E) ~ E^(PLind-1); the 4.8-59 eV range is entirely in
    the EUV band so the EUV/X-ray ratio cancels."""
    order = np.argsort(E_tb)
    Es, ss = E_tb[order], s_tb[order]
    sel = (Es >= 4.80) & (Es <= 59.18)
    E = Es[sel]
    s_TB = ss[sel]
    s_BPL = sigma_HeI23S(E)

    def integ(sig):
        return np.trapz(sig * E ** (PLind - 1.0), E)

    bpl = integ(s_BPL)
    out = []
    for cap in [np.inf, 500.0, 200.0, 100.0, 50.0]:
        out.append((cap, integ(np.minimum(s_TB, cap)) / bpl,
                    int((s_TB > cap).sum())))
    return out


def print_resonance_impact(E_tb, s_tb):
    print("\n=== Autoionizing-resonance impact on the He 2^3 S rate (PLind=-1.3) ===")
    print("TOTAL rate ratio TOPbase/EXHALE vs upper clip on sigma (spike test)")
    print(f"{'cap[Mb]':>10} {'TB/BPL':>9} {'pts>cap':>9}")
    for cap, ratio, npts in resonance_impact(E_tb, s_tb, PLind=-1.3):
        cs = "inf" if not np.isfinite(cap) else f"{cap:.0f}"
        print(f"{cs:>10} {ratio:>9.3f} {npts:>9d}")
    print("=> the apparent ~13% (no clip) collapses to ~1-2% once the 2-4 "
          "under-resolved peaks are removed: the net resonance effect is small.")


def print_tables(heI, he23S):
    E, s_ex, s_ve, ratio = heI
    print("\n=== He I: EXHALE 2-term vs Verner+1996 (Mb) ===")
    print(f"{'E[eV]':>8} {'EXHALE':>10} {'Verner':>10} {'ratio':>8}")
    for Eq in [24.6, 30, 50, 100, 200, 300, 500]:
        a = float(sigma_HeI(Eq)); b = float(sigma_verner("HeI", Eq))
        print(f"{Eq:>8.1f} {a:>10.4f} {b:>10.4f} {a/b:>8.3f}")
    m = (E >= 24.6) & (E <= 500)
    print(f"ratio over 24.6-500 eV: min={ratio[m].min():.3f}, "
          f"max={ratio[m].max():.3f}, median={np.median(ratio[m]):.3f}")

    E2, e_ex, e_pw, r2 = he23S
    print("\n=== He 2^3 S: EXHALE broken-PL vs p-winds/Norcross (Mb) ===")
    print(f"{'E[eV]':>8} {'EXHALE':>10} {'p-winds':>10} {'ratio':>8}")
    for Eq in [4.781, 6, 8, 12, 20, 30, 45, 55, 59]:
        a = float(sigma_HeI23S(Eq)); b = float(sigma_pwinds(Eq))
        print(f"{Eq:>8.2f} {a:>10.4f} {b:>10.4f} {a/b:>8.3f}")
    mm = np.isfinite(r2)
    print(f"ratio over 4.78-59 eV: min={r2[mm].min():.3f}, "
          f"max={r2[mm].max():.3f}, median={np.median(r2[mm]):.3f}")


if __name__ == "__main__":
    fig_overview()
    heI = fig_HeI_vs_verner()
    he23S = fig_HeI23S_vs_pwinds()
    E_tb, s_tb = fig_HeI23S_topbase()
    print_tables(heI, he23S)
    print_topbase_table(E_tb, s_tb)
    print_rate_impact(E_tb, s_tb)
    print_resonance_impact(E_tb, s_tb)
    print(f"\nFigures written to {OUTDIR}")
