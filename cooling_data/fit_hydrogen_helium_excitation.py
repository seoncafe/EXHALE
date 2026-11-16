"""Closed-form H I, He I and He II collisional-excitation data.

Every coefficient this script prints is transcribed into
src/modules/radiation/Cool_coeff.f90; the Fortran comments name this file as
their source.

What is fitted, and why in this split.

  He I, 1^1S -> every upper term n <= 5 (full sum).  The coronal cooling
      of ground-state helium when the 2^3S metastable is not a level of the
      network: every excitation is followed by a radiative decay, so the
      excitation energy is lost from the electron gas.  Collision strengths:
      bray2000_helium.adopted_upsilon, the Bray, Burgess, Fursa & Tully
      (2000, A&AS 146, 481) Table 2 values inside 10^3.75-10^5.75 K and the
      CHIANTI v11 he_1 temperature dependence anchored to them outside it,
      the strengths the metastable network uses; energies the g-weighted
      observed term energies.

  He I, 1^1S -> every upper term EXCEPT 2^3S.  When the network tracks the
      metastable, the 1^1S -> 2^3S excitation is the network's own rate q13,
      and the energy it moves is charged with that rate together with its
      detailed-balance reverse.  The sum used beside it must then leave that
      one transition out, or the 19.82 eV of every excitation would be
      charged twice.

  He II, 1s -> every upper level (full sum).

  H I, the effective collision strengths of 1s-2s and 1s-2p (the two 2p
      fine-structure levels summed).  These are the excitation rates of the
      H(n=2) model (hydrogen_n2_rates.f90) and the n = 2 part of the H I
      cooling, so the population model and the cooling use one rate set.

  H I, 1s -> n >= 3 (full sum over the upper levels 3s ... 5g).

Energies.  Both the Boltzmann factor and the radiated energy use the OBSERVED
level energies of the .elvlc file, not the .scups header energy.  For these
three ions the two agree to 0.1 per cent or better (printed below), so the
choice moves the result by at most 1.5 per cent at 3000 K, but the observed
energy is the physical one.

Collision strengths: CHIANTI v11.0.2 ($XUVTOP), Burgess & Tully descaling
done by chianti_cooling.upsilon.  Sources named in the CHIANTI files:
  h_1   Anderson, Ballance, Badnell & Summers (2000, J. Phys. B 33, 1255;
        revised 2002)
  he_1  Sawey & Berrington (1993, ADNDT 55, 81); Bray, Burgess, Fursa &
        Tully (2000, A&AS 146, 481); used by CHIANTI for 3.0 < log T < 5.75
        (for He I the CHIANTI sum is printed as a comparison only)
  he_2  Ballance (2003, R-matrix with pseudostates, private communication
        to CHIANTI)

Fit forms.
  Sums:        Lambda(T) = C T^-1/2 (a + ln(1 + T/T0)) exp(-Tex/T)
               [erg cm^3 s^-1 per (n_e n_lower)]; the fitted Tex is the
               effective threshold of the sum.
  Upsilon:     log10 Upsilon = c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4,
               x = log10(T/1e4 K), fitted over 1e3-1e5 K; the Fortran holds
               T inside that range (constant Upsilon outside it).

Run:  XUVTOP=/path/to/CHIANTI/dbase python3 fit_hydrogen_helium_excitation.py
Writes hydrogen_helium_excitation_fits.txt.
"""

import os
import sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chianti_cooling as cc
import bray2000_helium as bh

K_B = 1.380649e-16          # erg/K
HC_K = 1.438776877          # h c / k [cm K]
# h^2 / ((2 pi m_e)^1.5 k^0.5) [cm^3 s^-1 K^1/2], CODATA 2018, the constant
# coll_rate_prefactor of Cool_coeff.f90
COLL_PREF = 6.62607015e-27**2 / ((2.0 * np.pi * 9.1093837015e-28)**1.5 * np.sqrt(K_B))


def load(elem, ion):
    d = cc.ion_dir(elem, ion)
    lev = cc.read_elvlc(os.path.join(d, f"{elem}_{ion}.elvlc"))
    trs = cc.read_scups(os.path.join(d, f"{elem}_{ion}.scups"))
    return lev, trs


def ground_sum(elem, ion, T, exclude=(), only=None):
    """Sum over the transitions out of level 1 of q_1u * dE_1u, observed dE."""
    lev, trs = load(elem, ion)
    g1 = lev[1][0]
    tot = np.zeros_like(T)
    for t in trs:
        if t["ll"] != 1 or t["ul"] in exclude:
            continue
        if only is not None and t["ul"] not in only:
            continue
        de = lev[t["ul"]][1] * HC_K * K_B
        tot += (COLL_PREF / (g1 * np.sqrt(T)) * cc.upsilon(t, T)
                * np.exp(-de / (K_B * T)) * de)
    return tot


def helium_ground_sum(T, exclude_23S=False):
    """He I coronal cooling out of 1^1S with the adopted Bray et al. strengths
    [erg cm^3 s^-1 per (n_e n_1^1S)], every Table 2 term (g(1^1S) = 1)."""
    tot = np.zeros_like(T)
    for term in bh.upper_terms(1):
        if exclude_23S and term == (2, 3, "S"):
            continue
        de = bh.term_energy_K(term) * K_B
        tot += (COLL_PREF / np.sqrt(T) * bh.adopted_upsilon(1, term, T)
                * np.exp(-de / (K_B * T)) * de)
    return tot


def upsilon_sum(elem, ion, T, uppers):
    lev, trs = load(elem, ion)
    out = np.zeros_like(T)
    for u in uppers:
        t = [x for x in trs if x["ll"] == 1 and x["ul"] == u][0]
        out += cc.upsilon(t, T)
    return out


def log_sum_form(T, lnC, a, T0, Tex):
    return (lnC - 0.5 * np.log(T)
            + np.log(np.clip(a + np.log1p(T / np.abs(T0)), 1e-300, None))
            - Tex / T)


def fit_sum(label, T, lam, p0, out):
    p, _ = curve_fit(log_sum_form, T, np.log(lam), p0=p0, maxfev=400000)
    r = np.exp(log_sum_form(T, *p)) / lam
    out.append(f"{label}")
    out.append(f"  C   = {np.exp(p[0]):.8e}   [erg cm^3 s^-1 K^1/2]")
    out.append(f"  a   = {p[1]:.8e}")
    out.append(f"  T0  = {abs(p[2]):.8e}   K")
    out.append(f"  Tex = {p[3]:.6f}   K  ({p[3] * K_B / 1.602176634e-12:.5f} eV)")
    out.append(f"  fit / sum over {T[0]:.0e}-{T[-1]:.0e} K: "
               f"{r.min():.4f} - {r.max():.4f}")
    return p


def fit_upsilon(label, uppers, out):
    T = np.logspace(3.0, 5.0, 201)
    x = np.log10(T / 1.0e4)
    u = upsilon_sum("h", 1, T, uppers)
    c = np.polyfit(x, np.log10(u), 4)[::-1]
    r = 10.0 ** np.polyval(c[::-1], x) / u
    out.append(f"{label}  (CHIANTI levels {uppers})")
    out.append("  c0..c4 = " + ", ".join(f"{v:.8e}" for v in c))
    out.append(f"  fit / CHIANTI over 1e3-1e5 K: {r.min():.5f} - {r.max():.5f}")
    return c


out = ["# Closed-form H I, He I, He II collisional-excitation data",
       "# (CHIANTI v11, observed level energies).",
       "# Generated by fit_hydrogen_helium_excitation.py", ""]

# Energies: .scups header against observed .elvlc, transitions out of level 1
for elem, ion in (("h", 1), ("he", 1), ("he", 2)):
    lev, trs = load(elem, ion)
    ratios = [t["de"] / (lev[t["ul"]][1] / 109737.31568)
              for t in trs if t["ll"] == 1]
    out.append(f"{elem}_{ion}: .scups dE / observed dE over the ground "
               f"transitions: {min(ratios):.5f} - {max(ratios):.5f}")
out.append("")

T_he1 = np.logspace(np.log10(3.0e3), 5.0, 300)
for label, excl, p0 in (("He I ground-level sum, all upper terms", False,
                         [np.log(1.1e-16), 0.15, 1.0e5, 2.3e5]),
                        ("He I ground-level sum without 1^1S -> 2^3S", True,
                         [np.log(1.1e-16), 0.1, 1.0e5, 2.4e5])):
    lam = helium_ground_sum(T_he1, excl)
    p = fit_sum(label + " (Bray et al. 2000 strengths)", T_he1, lam, p0, out)
    Tq = np.array([3e3, 5e3, 1e4, 2e4, 5e4])
    ch = ground_sum("he", 1, Tq, exclude=(2,) if excl else ())
    out.append("  CHIANTI he_1 sum / this sum at 3e3, 5e3, 1e4, 2e4, 5e4 K: "
               + ", ".join(f"{v:.4f}" for v in ch / helium_ground_sum(Tq, excl)))
    out.append("  probe reference at 3e3, 5e3, 1e4, 2e4, 5e4 K: "
               + ", ".join(f"{v:.6e}" for v in helium_ground_sum(Tq, excl)))
    if not excl:
        ates = 1.1e-19 * Tq**0.082 * np.exp(-2.3e5 / Tq)
        out.append("  ATES coefficient 1.1e-19 T^0.082 exp(-2.3e5/T) / this sum: "
                   + ", ".join(f"{v:.3f}" for v in ates / helium_ground_sum(Tq, excl)))
    else:
        full = helium_ground_sum(Tq, False)
        out.append("  share of the 2^3S channel in the full sum: "
                   + ", ".join(f"{v:.3f}" for v in 1.0 - helium_ground_sum(Tq, True) / full))
    out.append("")
T_he2 = np.logspace(np.log10(5.0e3), np.log10(2.0e5), 300)
fit_sum("He II ground-level sum, all upper levels", T_he2,
        ground_sum("he", 2, T_he2), [np.log(2.3e-16), 0.5, 3.0e5, 4.7e5],
        out)
out.append("")
fit_upsilon("H I Upsilon(1s-2s)", [2], out)
out.append("")
fit_upsilon("H I Upsilon(1s-2p)", [3, 4], out)
out.append("")
T_h = np.logspace(np.log10(3.0e3), 5.0, 300)
fit_sum("H I ground-level sum over the n >= 3 upper levels", T_h,
        ground_sum("h", 1, T_h, exclude=(2, 3, 4)),
        [np.log(3.0e-17), 0.5, 2.5e4, 1.4e5], out)
out.append("")

open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  "hydrogen_helium_excitation_fits.txt"), "w").write(
    "\n".join(out) + "\n")
print("\n".join(out))
