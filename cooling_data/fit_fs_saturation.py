"""Density-dependent fine-structure treatment for [C II] 158um / [O I] 63um.

The cno_chianti coronal fits include the FS floor terms unsaturated; their
critical densities (n_e ~ 30-1e5 cm^-3, n_H ~ 1e4-1e6 cm^-3) are far below
the atmosphere base density, so the floors overestimate base cooling.
This script builds the saturating hybrid

    Lambda_eff(T,ne,nHI) = W_FS(T,ne,nHI)/ne + Lambda_rem(T)

where W_FS is the EXACT two-level FS cooling per ion [erg/s]

    W_FS = f_l(T) * hv * A * x * Cdex / (A + Cdex*(1+x)),
    x    = (g_u/g_l) exp(-E/kT),
    Cdex = ne*k_e(T) + nHI*k_H(T)            (de-excitation rate, s^-1)
    f_l  = g_l / Z_term(T)                   (ground-term Boltzmann fraction)

and Lambda_rem is a multi-exp refit of the CHIANTI ground-term curve with
the (1->2) FS channel removed.  Limits: ne,nHI -> 0 reproduces the CHIANTI
coronal floor channel exactly (same Upsilon, same f_l); ne or nHI -> inf
saturates at the LTE value f_l*hv*A*x/(1+x), which for the f_l-weighted
two-level form equals the exact multi-level LTE population result.

Atomic data (CHIANTI v11.0.2):
  C II 2P:  E21 = 63.397 cm^-1 (91.21 K),  A21 = 2.290e-6 s^-1, g=2,4
            (two-level is exact: the ground term has only two levels)
  O I  3P:  E21 = 158.265 cm^-1 (227.71 K), A21 = 8.542e-5 s^-1, g=5,3
            (63um channel; the weaker [O I] 146um (3P1-3P0, A=1.64e-5)
             stays coronal in the remainder -- a few-% piece of the floor)
  k_e from the CHIANTI .scups Upsilon (Burgess-Tully descaled), fitted as
            Ups(T) = a + b*ln(1+T/T0);
  k_H de-excitation (approximate, as in the legacy use_2lev_cool branch):
            O I 63um + H : 4.2e-11*(T/100)^0.67   (Draine 2011/Lique+2017)
            C II 158um + H: 4.0e-11               (Goldsmith+2012)
            !To Be Checked/AIOLOS tuning? -- flagged approximate upstream.

Run:  python3 fit_fs_saturation.py   (prints Fortran-ready coefficients)
"""

import os
import numpy as np
from scipy.optimize import curve_fit, least_squares
from chianti_cooling import (read_elvlc, read_scups, upsilon, ion_dir,
                             cooling_effective, _level_energy_cm1)

K_B = 1.380649e-16
HC_OVER_K = 1.438776877      # cm K
COLL_PREF = 8.629e-6
EV_ERG = 1.602176634e-12

T = np.logspace(3.0, 5.0, 201)

FS = {
    # name: (elem, ion, g_l, g_u, E_cm1(1->2), A21, kH_deex(T))
    "CII": ("c", 2, 2.0, 4.0, 63.397, 2.290e-6,
            lambda T: 4.0e-11 * np.ones_like(T)),
    "OI":  ("o", 1, 5.0, 3.0, 158.265, 8.542e-5,
            lambda T: 4.2e-11 * (T/100.0)**0.67),
}


def get_ups12(elem, ion):
    d = ion_dir(elem, ion)
    trans = read_scups(os.path.join(d, f"{elem}_{ion}.scups"))
    tr = [t for t in trans if t["ll"] == 1 and t["ul"] == 2][0]
    return upsilon(tr, T)


def term_partition(elem, ion):
    """Ground-term Boltzmann partition & level-1 fraction (0.4 eV cut,
    same convention as cooling_effective(pop='ground_term'))."""
    lev = _level_energy_cm1(read_elvlc(
        os.path.join(ion_dir(elem, ion), f"{elem}_{ion}.elvlc")))
    thr = 0.4*EV_ERG/(HC_OVER_K*K_B)
    lows = [(g, e) for i, (g, e) in sorted(lev.items()) if e <= thr]
    Z = np.zeros_like(T)
    for g, e in lows:
        Z += g*np.exp(-e*HC_OVER_K/T)
    return lows[0][0]/Z          # f_1(T)


def w_fs(name, Tv, ne, nHI, ups, f1):
    """Two-level FS cooling per ion [erg/s]."""
    el, ion, g_l, g_u, Ecm, A21, kH = FS[name]
    Estar = Ecm*HC_OVER_K                       # E/k [K]
    hv = K_B*Estar
    k_e = COLL_PREF*ups/(g_u*np.sqrt(Tv))       # de-excitation by e
    Cdex = ne*k_e + nHI*kH(Tv)
    x = (g_u/g_l)*np.exp(-Estar/Tv)
    return f1*hv*A21*x*Cdex/(A21 + Cdex*(1.0 + x))


def fs_coronal_channel(name, ups, f1):
    """The (1->2) e-collision coronal channel as cooling_effective counts
    it: f1 * q_lu * dE."""
    el, ion, g_l, g_u, Ecm, A21, kH = FS[name]
    Estar = Ecm*HC_OVER_K
    hv = K_B*Estar
    q_lu = COLL_PREF*ups/(g_l*np.sqrt(T))*np.exp(-Estar/T)
    return f1*q_lu*hv


def fit_multiexp(lam, Tguess):
    y = np.log(np.maximum(lam, 1e-300)*np.sqrt(T))
    n = len(Tguess)

    def resid(th):
        A = np.exp(th[:n]); Ti = np.exp(th[n:])
        return np.log(np.sum(A[:, None]*np.exp(-Ti[:, None]/T[None, :]),
                             axis=0)) - y

    th0 = np.concatenate([np.full(n, np.log(lam[-1]*np.sqrt(T[-1])/n)),
                          np.log(np.asarray(Tguess, float))])
    s = least_squares(resid, th0, method="lm", max_nfev=60000)
    A = np.exp(s.x[:n]); Ti = np.exp(s.x[n:])
    o = np.argsort(Ti)
    return A[o], Ti[o]


def multiexp(A, Ti, Tv):
    return np.sum(np.asarray(A)[:, None]
                  * np.exp(-np.asarray(Ti)[:, None]/Tv[None, :]),
                  axis=0)/np.sqrt(Tv)


REM_GUESS = {"CII": [7000.0, 64000.0, 125000.0, 246000.0],
             "OI":  [400.0, 19000.0, 32000.0, 70000.0, 179000.0]}

out = {}
for name, (el, ion, g_l, g_u, Ecm, A21, kH) in FS.items():
    ups = get_ups12(el, ion)
    f1 = term_partition(el, ion)

    # Upsilon fit: log-quadratic, Ups = 10^(c0 + c1*x + c2*x^2),
    # x = log10(T/1e4)  (the BT log form degenerates for these curves)
    x = np.log10(T/1e4)
    cu = np.polyfit(x, np.log10(ups), 3)        # [c3, c2, c1, c0]
    ups_fit = 10.0**np.polyval(cu, x)
    e_ups = np.abs(ups_fit/ups - 1).max()

    # remainder = full ground-term coronal curve minus the (1->2) channel
    total = cooling_effective(el, ion, T, pop="ground_term")
    rem = total - fs_coronal_channel(name, ups, f1)
    A, Ti = fit_multiexp(rem, REM_GUESS[name])
    rem_fit = multiexp(A, Ti, T)
    mask = rem > 1e-3*rem.max()
    e_rem = np.abs(rem_fit/rem - 1)[mask].max()

    # low-density consistency: hybrid(ne->0) vs CHIANTI total.
    # (n_crit,e of [C II] 158um is only ~20 cm^-3, so the limit must be
    #  taken at ne << 1 to avoid genuine saturation entering the check)
    ne0 = 1.0e-6
    lam0 = w_fs(name, T, ne0, 0.0, ups_fit, f1)/ne0 + rem_fit
    e_lo = np.abs(lam0/total - 1).max()

    print(f"--- {name} ---")
    print(f"  log10 Ups = {cu[3]:.8g} + {cu[2]:.8g}*x + {cu[1]:.8g}*x^2"
          f" + {cu[0]:.8g}*x^3,  x=log10(T/1e4)   (max err {e_ups*100:.2f}%)")
    print(f"  f1 example: f1(1e3)={f1[0]:.4f} f1(1e4)={f1[100]:.4f}")
    print(f"  remainder fit ({len(A)} terms), max err {e_rem*100:.2f}%:")
    for a, t in zip(A, Ti):
        print(f"    A={a:.8e}  Ti={t:.6g}")
    print(f"  low-density hybrid vs CHIANTI total: max dev {e_lo*100:.2f}%")
    out[name] = dict(cu=cu, A=A, Ti=Ti, ups_tab=ups, f1=f1)

    # f1(T) analytic check: exact ground-term partition is closed-form
    el2, ion2 = FS[name][0], FS[name][1]
    if name == "CII":
        f1_exact = 2.0/(2.0 + 4.0*np.exp(-91.21/T))
    else:
        f1_exact = 5.0/(5.0 + 3.0*np.exp(-227.71/T) + np.exp(-326.58/T))
    print(f"  f1 closed-form vs table: max dev "
          f"{np.abs(f1_exact/f1 - 1).max()*100:.3f}%")

# --------------------------- saturation behaviour table ------------------
print("\n=== saturation check: Lambda_eff = W_FS/ne + rem  [erg cm^3/s] ===")
print("    (xe = ne/nH; nHI = nH*(1-xe))")
for name in FS:
    r = out[name]
    print(f"\n{name}:  T=4000 K")
    i = np.argmin(abs(T - 4000.0))
    for nH, xe in [(1e4, 0.5), (1e8, 1e-2), (1e10, 1e-3), (1e12, 1e-3)]:
        ne, nHI = nH*xe, nH*(1 - xe)
        Tv = np.array([4000.0])
        upsv = f_ups(Tv, *r["ups"]) if False else \
            np.interp(Tv, T, r["ups_tab"])
        f1v = np.interp(Tv, T, r["f1"])
        w = w_fs(name, Tv, ne, nHI, upsv, f1v)[0]
        rem_v = multiexp(r["A"], r["Ti"], Tv)[0]
        cor = (w_fs(name, Tv, 1e-6, 0.0, upsv, f1v)/1e-6)[0] + rem_v
        print(f"  nH={nH:8.0e} xe={xe:7.0e}:  Lam_eff={w/ne+rem_v:11.4e}"
              f"   coronal={cor:11.4e}   suppression={(w/ne+rem_v)/cor:8.2e}")
