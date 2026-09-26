"""Figures comparing the coded metal line-cooling coefficients with the rebuilt ones.

Writes to ../docs/figures/:
  atomic_audit_metal_cooling_mg2.pdf      Mg II: coded fit (.scups energy, h&k)
                                          vs the observed-energy CHIANTI sum
                                          (all channels) and its new fit
  atomic_audit_metal_cooling_ca2.pdf      Ca II: coded resonance fit vs the
                                          (T, n_e) statistical-equilibrium table
  atomic_audit_metal_cooling_mg1.pdf      Mg I: the same, with the sensitivity
                                          to the 3P_J fine-structure collisions
                                          that CHIANTI omits
  atomic_audit_metal_cooling_fe1.pdf      Fe I: coded Boltzmann-manifold table
                                          vs the ground-level and ground-term
                                          limits
  atomic_audit_metal_cooling_wasp121b.pdf the Mg I, Mg II and Ca II cooling
                                          rates, coded and rebuilt, on the
                                          T, n_e and ion columns of an existing
                                          WASP-121 b output (no new run)

The Ca II and Mg I table values are read from
metal_cooling_density_resolved_tables.txt and interpolated exactly as
remainder_table_value in Cool_coeff.f90 does (bilinear in log10 T, log10 n_e,
edges clamped).  Inputs: run magnesium_ii_line_cooling.py,
metal_cooling_density_resolved.py and fe1_level_population_limits.py first.

RUN:  XUVTOP=<dbase> python3 atomic_audit_metal_cooling_figures.py
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "docs", "figures")
LOGT = np.linspace(3.0, 5.0, 41)
LOGNE = np.linspace(0.0, 14.0, 29)
NE_SHOWN = [1.0e2, 1.0e6, 1.0e8, 1.0e9, 1.0e10, 1.0e12]
OLD_COLOR = "0.45"


def ramp(k, n):
    """One-hue sequential ramp (light to dark) for the n_e curves."""
    return plt.cm.Blues(0.35 + 0.65 * k / max(n - 1, 1))


# coded closed forms, exactly as in Cool_coeff.f90 before the patch
def cool_mgii_old(T):
    ups = 16.230522 + 20.344762 * np.log(1.0 + T / 121896.19)
    return 8.629e-6 / (2.0 * np.sqrt(T)) * ups * 6.85184290e-12 * np.exp(-49627.696 / T)


def cool_caii_old(T):
    ups = 14.741545 + 19.491995 * np.log(1.0 + T / 36218.324)
    return 8.629e-6 / (2.0 * np.sqrt(T)) * ups * 5.03698606e-12 * np.exp(-36482.742 / T)


def cool_mgi_old(T):
    ups = 0.35579231 + 15.122288 * np.log(1.0 + T / 94126.484)
    return 8.629e-6 / np.sqrt(T) * ups * 6.95044482e-12 * np.exp(-50341.867 / T)


def read_table(name):
    s = open(os.path.join(HERE, "metal_cooling_density_resolved_tables.txt")).read()
    m = re.search(r"cool_logLrem_" + name + r"\(NCOOLT,NCOOLNE\)\s*=\s*reshape\(\s*\[(.*?)\]", s, re.S)
    v = np.array([float(x) for x in m.group(1).replace("&", " ").replace("d0", "")
                  .replace(",", " ").split()])
    return v.reshape(LOGT.size, LOGNE.size, order="F")


def remainder_table_value(logl, T, ne):
    """Python copy of remainder_table_value (Cool_coeff.f90)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ne = np.broadcast_to(np.asarray(ne, dtype=float), T.shape)
    post = (np.log10(np.maximum(T, 1.0)) - LOGT[0]) / 0.05 + 1.0
    posn = (np.log10(np.maximum(ne, 1.0)) - LOGNE[0]) / 0.5 + 1.0
    kt = np.clip(np.floor(post).astype(int), 1, LOGT.size - 1)
    ft = np.where(post <= 1.0, 0.0, np.where(post >= LOGT.size, 1.0, post - kt))
    ke = np.clip(np.floor(posn).astype(int), 1, LOGNE.size - 1)
    fn = np.where(posn <= 1.0, 0.0, np.where(posn >= LOGNE.size, 1.0, posn - ke))
    kt -= 1
    ke -= 1
    l0 = logl[kt, ke] + ft * (logl[kt + 1, ke] - logl[kt, ke])
    l1 = logl[kt, ke + 1] + ft * (logl[kt + 1, ke + 1] - logl[kt, ke + 1])
    return 10.0**(l0 + fn * (l1 - l0))


def figure_mg2():
    d = np.loadtxt(os.path.join(HERE, "mg2_line_cooling_observed_energies.txt"))
    T, tot, hk, old, new = d.T
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.2, 6.4), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1.2]))
    a1.plot(T, tot, color="k", lw=2.0, label=r"CHIANTI, observed energies, all channels")
    a1.plot(T, hk, color=plt.cm.Blues(0.75), lw=1.4, ls="-.", label=r"CHIANTI, observed energies, h\&k only")
    a1.plot(T, new, color=plt.cm.Oranges(0.7), lw=1.2, ls=":", label=r"new fit (patch)")
    a1.plot(T, old, color=OLD_COLOR, lw=1.4, ls="--", label=r"coded fit (.scups energy, h\&k)")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_ylim(1e-30, 3e-18)
    a1.set_ylabel(r"$\Lambda_{\rm Mg\,II}$ [erg cm$^3$ s$^{-1}$]")
    a1.legend(fontsize=8, loc="lower right")
    a1.set_title(r"Mg\,II line cooling per $n_e n_{\rm Mg\,II}$ (coronal limit)")
    a2.plot(T, old / tot, color=OLD_COLOR, lw=1.4, ls="--", label=r"coded / new table")
    a2.plot(T, new / tot, color=plt.cm.Oranges(0.7), lw=1.2, ls=":", label=r"new fit / new table")
    a2.axhline(1.0, color="k", lw=0.6)
    a2.set_yscale("log"); a2.set_ylim(0.5, 8.0)
    a2.set_xlabel(r"$T$ [K]"); a2.set_ylabel(r"ratio")
    a2.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "atomic_audit_metal_cooling_mg2.pdf"))
    plt.close(fig)


def figure_table_ion(tag, label, old_fn, fname, extra=None):
    logl = read_table(tag)
    T = np.logspace(np.log10(2.0e3), np.log10(3.0e4), 300)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.2, 6.6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1.3]))
    old = old_fn(T)
    for k, ne in enumerate(NE_SHOWN):
        lam = remainder_table_value(logl, T, ne)
        c = ramp(k, len(NE_SHOWN))
        a1.plot(T, lam, color=c, lw=1.5,
                label=r"table, $n_e = 10^{" + f"{np.log10(ne):.0f}" + r"}$ cm$^{-3}$")
        a2.plot(T, lam / old, color=c, lw=1.5)
    a1.plot(T, old, color=OLD_COLOR, lw=1.5, ls="--", label=r"coded fit (resonance line only)")
    if extra is not None:
        extra(a1, a2, T, old)
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
    a1.set_title(label + r" line cooling per $n_e n_{\rm ion}$")
    a1.legend(fontsize=7.5, loc="lower right")
    a2.axhline(1.0, color="k", lw=0.6)
    a2.set_yscale("log")
    a2.set_xlabel(r"$T$ [K]"); a2.set_ylabel(r"table / coded fit")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname))
    plt.close(fig)


def mg1_fine_structure_band(a1, a2, T, old):
    """Mg I at n_e = 1e8 with a constant Upsilon = 1 and 10 added to each
    3P_J - 3P_J' pair, which CHIANTI v11.0.2 does not provide."""
    from multilevel_statistical_equilibrium import load_model_ion, line_cooling
    m = load_model_ion("mg", 1)
    tt = np.logspace(np.log10(2.0e3), np.log10(3.0e4), 40)
    base = line_cooling(m, tt, [1.0e8], 0)[:, 0]
    for ups, ls in ((1.0, (0, (3, 1))), (10.0, (0, (1, 1)))):
        extra = [(1, 2, ups), (1, 3, ups), (2, 3, ups)]
        lam = line_cooling(m, tt, [1.0e8], 0, extra_upsilon=extra)[:, 0]
        a2.plot(tt, lam / base * np.interp(tt, T, remainder_table_value(read_table("MgI"), T, 1.0e8)) / np.interp(tt, T, old),
                color=ramp(2, len(NE_SHOWN)), lw=1.0, ls=ls,
                label=r"$n_e = 10^{8}$, $\Upsilon(^3P_J\,^3P_{J'}) = " + f"{ups:g}" + r"$")
    a2.legend(fontsize=7, loc="upper right")


def figure_fe1():
    d = np.loadtxt(os.path.join(HERE, "fe1_level_population_limits.txt"))
    T, boltz, term, ground = d.T
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.2, 6.0), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1.2]))
    a1.plot(T, boltz, color=OLD_COLOR, lw=1.5, ls="--", label=r"coded: Boltzmann over 17 metastable levels ($n_e\to\infty$)")
    a1.plot(T, term, color=plt.cm.Blues(0.6), lw=1.3, ls="-.", label=r"a$^5$D term only")
    a1.plot(T, ground, color=plt.cm.Blues(0.95), lw=1.6, label=r"a$^5$D$_4$ ground level only ($n_e\to 0$)")
    a1.set_xscale("log"); a1.set_yscale("log"); a1.set_ylim(1e-26, 3e-18)
    a1.set_ylabel(r"$\Lambda_{\rm Fe\,I}$ [erg cm$^3$ s$^{-1}$]")
    a1.set_title(r"Fe\,I line cooling per $n_e n_{\rm Fe\,I}$: limits of the lower-level populations")
    a1.legend(fontsize=7.5, loc="lower right")
    a2.plot(T, ground / boltz, color=plt.cm.Blues(0.95), lw=1.6, label=r"ground level / coded")
    a2.plot(T, term / boltz, color=plt.cm.Blues(0.6), lw=1.3, ls="-.", label=r"a$^5$D term / coded")
    a2.axhline(1.0, color="k", lw=0.6)
    a2.set_ylim(0.4, 1.3)
    a2.set_xlabel(r"$T$ [K]"); a2.set_ylabel(r"ratio")
    a2.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "atomic_audit_metal_cooling_fe1.pdf"))
    plt.close(fig)


def new_mgii_fit(T):
    d = open(os.path.join(HERE, "mg2_line_cooling_observed_energies.txt")).read()
    a, b, t0, c, tx = [float(x) for x in
                       re.search(r"fit a, b, T0, c, T_x = (.*)", d).group(1).split()]
    e_erg = 35730.366 * 1.98644586e-16
    t_hk = 35730.366 * 1.438776877
    ups = a + b * np.log(1.0 + T / t0)
    return (8.629e-6 / (2.0 * np.sqrt(T)) * ups * e_erg * np.exp(-t_hk / T)
            + c * np.exp(-tx / T) / np.sqrt(T))


def wasp121b_profile():
    """Coded and rebuilt Mg I, Mg II and Ca II cooling rates on an existing
    WASP-121 b output (T, n_e, ion densities as written by that run)."""
    base = os.path.join(HERE, "..", "WASP-121b", "output")
    cb = np.loadtxt(os.path.join(base, "Cooling_breakdown.txt"))
    ions = np.loadtxt(os.path.join(base, "Ion_species.txt"))
    hdr = [l for l in open(os.path.join(base, "Ion_species.txt")) if l.startswith("# columns")][0].split()[2:]
    col = {nm: k for k, nm in enumerate(hdr)}
    r, T, ne, cool = cb[:, 0], cb[:, 1], cb[:, 2], cb[:, 3]
    metal_names = "CI CII CIII OI OII OIII NI NII NIII MgI MgII MgIII SiI SiII SiIII CaI CaII CaIII NaI NaII KI KII SI SII FeI FeII FeIII".split()
    mcol = {nm: 11 + k for k, nm in enumerate(metal_names)}
    assert np.allclose(ions[:, 0], r)
    n = {nm: ions[:, col[nm]] for nm in ("MgI", "MgII", "CaII")}
    old = {"MgII": ne * n["MgII"] * cool_mgii_old(T),
           "CaII": ne * n["CaII"] * cool_caii_old(T),
           "MgI": ne * n["MgI"] * cool_mgi_old(T)}
    new = {"MgII": ne * n["MgII"] * new_mgii_fit(T),
           "CaII": ne * n["CaII"] * remainder_table_value(read_table("CaII"), T, ne),
           "MgI": ne * n["MgI"] * remainder_table_value(read_table("MgI"), T, ne)}
    written = {nm: cb[:, mcol[nm]] for nm in ("MgI", "MgII", "CaII")}
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 6.6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1.2]))
    colors = {"MgII": plt.cm.Blues(0.9), "CaII": plt.cm.Oranges(0.7), "MgI": plt.cm.Greens(0.7)}
    labels = {"MgII": r"Mg\,II", "CaII": r"Ca\,II", "MgI": r"Mg\,I"}
    a1.plot(r - 1.0, cool, color="k", lw=1.2, label=r"total cooling (run)")
    for nm in ("MgII", "CaII", "MgI"):
        a1.plot(r - 1.0, new[nm], color=colors[nm], lw=1.6, label=labels[nm] + r" rebuilt")
        a1.plot(r - 1.0, old[nm], color=colors[nm], lw=1.0, ls="--", label=labels[nm] + r" coded")
        a2.plot(r - 1.0, new[nm] / old[nm], color=colors[nm], lw=1.4, label=labels[nm])
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_ylim(1e-16, 1e-4)
    a1.set_ylabel(r"cooling [erg cm$^{-3}$ s$^{-1}$]")
    a1.set_title(r"WASP-121\,b output profile (coefficient swap, no new run)")
    a1.legend(fontsize=7, ncol=2, loc="lower left")
    a2.axhline(1.0, color="k", lw=0.6)
    a2.set_yscale("log")
    a2.set_xlabel(r"$r/R_p - 1$"); a2.set_ylabel(r"rebuilt / coded")
    a2.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "atomic_audit_metal_cooling_wasp121b.pdf"))
    plt.close(fig)
    # numbers for the memo
    chk = {nm: np.nanmax(np.abs(old[nm] / np.where(written[nm] > 0, written[nm], np.nan) - 1.0))
           for nm in old}
    print("coded coefficients reproduce the run's written columns: max rel. diff",
          {k: f"{v:.1e}" for k, v in chk.items()})
    sel = r - 1.0 > 0.0
    for nm in ("MgII", "CaII", "MgI"):
        frac_old = old[nm] / cool
        frac_new = new[nm] / cool
        k = int(np.argmax(frac_old))
        print(f"{nm}: max share of total cooling, coded {frac_old.max():.3e} at r-1={r[k]-1:.3e} "
              f"(T={T[k]:.0f} K, n_e={ne[k]:.2e}); rebuilt/coded there {new[nm][k]/old[nm][k]:.3f}; "
              f"rebuilt share max {frac_new.max():.3e}; ratio range {np.min(new[nm][sel]/old[nm][sel]):.3f}"
              f"-{np.max(new[nm][sel]/old[nm][sel]):.3f}")
    print(f"profile: T {T.min():.0f}-{T.max():.0f} K, n_e {ne.min():.2e}-{ne.max():.2e} cm^-3")
    return r, T, ne


def main():
    plt.rcParams.update({"font.size": 10})
    figure_mg2()
    figure_table_ion("CaII", r"Ca\,II", cool_caii_old, "atomic_audit_metal_cooling_ca2.pdf")
    figure_table_ion("MgI", r"Mg\,I", cool_mgi_old, "atomic_audit_metal_cooling_mg1.pdf",
                     extra=mg1_fine_structure_band)
    figure_fe1()
    wasp121b_profile()
    print("wrote docs/figures/atomic_audit_metal_cooling_{mg2,ca2,mg1,fe1,wasp121b}.pdf")


if __name__ == "__main__":
    main()
