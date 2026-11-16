"""Ground-term fine-structure statistical equilibrium for the C/N/O coolants.

The cno_chianti coronal fits carry the ground-term fine-structure (FS) lines in
the optically thin, low-density limit.  Their critical densities are of order
1e0-1e5 cm^-3, far below the base density of an irradiated planetary
atmosphere, so the coronal form overestimates the FS cooling there by many
decades.  This script builds, for every C/N/O coolant whose ground term is
split, the saturating hybrid

    Lambda_eff(T,ne,nHI) = W_FS(T,ne,nHI,beta)/ne + Lambda_rem(T)

where W_FS [erg/s per ion] is the EXACT statistical-equilibrium emission of the
ground term,

    solve  sum_j f_j R_ji = f_i sum_j R_ij,   sum_i f_i = 1,
    R_ul = C_ul + beta_ul A_ul,   R_lu = C_ul (g_u/g_l) exp(-E_ul/kT),
    C_ul = ne k_e,ul(T) + nHI k_H,ul(T),
    W_FS = sum_{u>l} f_u beta_ul A_ul k_B E_ul,

and Lambda_rem is a multi-exponential refit of the CHIANTI ground-term cooling
curve with EVERY within-ground-term channel removed.

Ions treated (all the C/N/O coolants with a split ground term):

    C I   2p2 3P_0,1,2   [C I]  609.1 / 370.4 um
    C II  2p  2P_1/2,3/2 [C II] 157.7 um
    N II  2p2 3P_0,1,2   [N II] 205.3 / 121.8 um
    O I   2p4 3P_2,1,0   [O I]   63.2 / 145.5 / 44.1 um   (inverted term)

N I and O II have a single-level 4S ground term; Mg I/II, Ca II and Na I have a
single ground level; Fe I is built from permitted lines only and Fe II is
already density-dependent through its 2-D statistical-equilibrium table.  So
these four are the complete set.

Limits.  ne, nHI -> 0 reproduces the CHIANTI coronal curve (the SE populations
collapse onto the ground level, and for the LS-coupled transitions that make up
the remainder the Boltzmann and ground-only lower-level weightings coincide
because Upsilon_lu is proportional to g_l).  ne or nHI -> infinity saturates at
the exact multilevel LTE emission sum_u f_u^Boltz A_ul h nu_ul.

Atomic data
  levels, A values, electron Upsilon : CHIANTI v11.0.2 (elvlc / wgfa / scups)
  k_H (H-atom de-excitation)
    C I, N II : Yan & Babb (2023), MNRAS 518, 6004, Tables 1 and 2
                (relaxation rates, tabulated 10-1e4 K; every value checked
                against the published article)
    O I       : Abrahamsson, Krems & Dalgarno (2007), ApJ 654, 1171,
                as tabulated in LAMDA oatom.dat (20-1000 K)
    C II      : Barinovs, van Hemert, Krems & Dalgarno (2005), ApJ 620, 537,
                as tabulated in LAMDA c+.dat (20-2000 K)

Run:  python3 fit_fs_saturation.py   (prints Fortran-ready coefficients)
"""

import os
import numpy as np
from scipy.optimize import least_squares
from chianti_cooling import (read_elvlc, read_scups, read_wgfa, upsilon,
                             ion_dir, cooling_effective, _level_energy_cm1)

K_B = 1.380649e-16
HC_OVER_K = 1.438776877      # cm K
COLL_PREF = 8.629e-6
EV_ERG = 1.602176634e-12

T = np.logspace(3.0, 5.0, 201)          # the range the CHIANTI fits cover

# --------------------------------------------------------------------------
# H-atom de-excitation rate coefficients k_H(u->l) [cm^3 s^-1], tabulated.
# Keys are (u,l) with levels indexed 1..n in ENERGY order.
# --------------------------------------------------------------------------
KH_TAB = {
    # Yan & Babb (2023) MNRAS 518, 6004, Table 1 (C + H relaxation)
    "CI": (np.array([10., 20., 50., 70., 100., 200., 500., 700., 1000.,
                     2000., 5000., 7000., 10000.]),
           {(2, 1): np.array([.199, .182, .157, .153, .152, .159, .189, .205,
                              .227, .284, .400, .458, .527])*1e-9,
            (3, 1): np.array([.086, .080, .077, .081, .088, .110, .149, .167,
                              .189, .245, .346, .393, .451])*1e-9,
            (3, 2): np.array([.252, .251, .256, .265, .283, .339, .450, .504,
                              .572, .743, 1.061, 1.213, 1.397])*1e-9}),
    # Yan & Babb (2023) Table 2 (N+ + H relaxation)
    "NII": (np.array([10., 20., 50., 70., 100., 200., 500., 700., 1000.,
                      2000., 5000., 7000., 10000.]),
            {(2, 1): np.array([.353, .333, .330, .324, .314, .290, .285, .292,
                               .306, .345, .415, .455, .511])*1e-9,
             (3, 1): np.array([.169, .156, .171, .184, .201, .244, .313, .342,
                               .376, .460, .624, .700, .790])*1e-9,
             (3, 2): np.array([.458, .461, .513, .544, .583, .678, .848, .927,
                               1.022, 1.256, 1.698, 1.907, 2.158])*1e-9}),
    # Abrahamsson, Krems & Dalgarno (2007), LAMDA oatom.dat
    "OI": (np.array([20., 35., 50., 70., 100., 150., 200., 250., 300., 350.,
                     400., 450., 500., 600., 700., 800., 900., 1000.]),
           {(2, 1): np.array([2.6, 2.8, 3.0, 3.2, 3.6, 4.2, 4.7, 5.2, 5.6,
                              6.0, 6.4, 6.7, 7.0, 7.6, 8.1, 8.5, 8.8,
                              9.2])*1e-10,
            (3, 1): np.array([0.85, 2.6, 2.7, 2.9, 3.2, 3.6, 4.1, 4.5, 4.8,
                              5.1, 5.4, 5.7, 5.9, 6.4, 6.8, 7.1, 7.4,
                              7.6])*1e-10,
            (3, 2): np.array([1.8, 2.1, 2.4, 3.3, 4.4, 5.8, 6.8, 7.6, 8.2,
                              8.7, 9.1, 9.4, 9.6, 10.0, 10.0, 10.0, 11.0,
                              11.0])*1e-10}),
    # Barinovs et al. (2005), LAMDA c+.dat
    "CII": (np.array([20., 40., 60., 80., 100., 140., 200., 300., 400., 600.,
                      800., 1000., 1500., 2000.]),
            {(2, 1): np.array([5.96, 6.79, 7.19, 7.42, 7.58, 7.84, 8.17, 8.63,
                               9.02, 9.67, 10.2, 10.66, 11.58,
                               12.31])*1e-10}),
}

# ion -> (CHIANTI element, stage, number of ground-term levels, amu)
FS_IONS = {
    "CI":  ("c", 1, 3, 12.011),
    "CII": ("c", 2, 2, 12.011),
    "NII": ("n", 2, 3, 14.007),
    "OI":  ("o", 1, 3, 15.999),
}

# Multi-exponential starting guesses for the remainder fits [K]
REM_GUESS = {"CI":  [15000., 25000., 73000., 177000.],
             "CII": [7000., 64000., 125000., 246000.],
             "NII": [23000., 50000., 128000., 250000.],
             "OI":  [19000., 32000., 70000., 179000.]}


def ground_term_data(elem, ion, nlev):
    """(g, Ek[K], A[u,l], {(l,u): scups transition}) for the ground term."""
    d = ion_dir(elem, ion)
    lev = _level_energy_cm1(read_elvlc(os.path.join(d, f"{elem}_{ion}.elvlc")))
    idx = sorted(lev)[:nlev]
    g = np.array([lev[i][0] for i in idx])
    Ek = np.array([lev[i][1]*HC_OVER_K for i in idx])
    A = np.zeros((nlev, nlev))
    for (ll, ul, wl, gf, a) in read_wgfa(os.path.join(d, f"{elem}_{ion}.wgfa")):
        if ll in idx and ul in idx and ul != ll:
            A[idx.index(ul), idx.index(ll)] = a
    ups = {}
    for tr in read_scups(os.path.join(d, f"{elem}_{ion}.scups")):
        if tr["ll"] in idx and tr["ul"] in idx:
            ups[(idx.index(tr["ll"]), idx.index(tr["ul"]))] = tr
    return g, Ek, A, ups


UPS_DEG = 4


def fit_upsilon(tr):
    """Quartic in x = log10(T/1e4) for log10(Upsilon) (max error 4.8%; a cubic
    reaches 8.8% on the C I 3P0-3P2 collision strength)."""
    u = upsilon(tr, T)
    x = np.log10(T/1e4)
    c = np.polyfit(x, np.log10(u), UPS_DEG)    # [c4, c3, c2, c1, c0]
    err = np.abs(10.0**np.polyval(c, x)/u - 1).max()
    return c, err, u


def fit_kH(Ttab, ktab):
    """Quadratic in u = log10(T/1e3) for log10(k_H), fitted over the
    tabulated temperatures >= 50 K and evaluated clamped to that range."""
    m = Ttab >= 50.0
    u = np.log10(Ttab[m]/1e3)
    c = np.polyfit(u, np.log10(ktab[m]), 2)    # [c2, c1, c0]
    err = np.abs(10.0**np.polyval(c, u)/ktab[m] - 1).max()
    return c, err, (Ttab[m].min(), Ttab[m].max())


def kH_eval(c, lim, Tv):
    u = np.log10(np.clip(Tv, lim[0], lim[1])/1e3)
    return 10.0**np.polyval(c, u)


def se_populations(g, Ek, A, C, beta, Tv):
    """Exact ground-term SE level fractions; C[u,l] de-excitation rate [s^-1]."""
    n = len(g)
    Tv = np.atleast_1d(Tv)
    f = np.zeros((n, Tv.size))
    for it, Ts in enumerate(Tv):
        R = np.zeros((n, n))
        for u in range(n):
            for l in range(u):
                Cul = C[u, l][it] if np.ndim(C[u, l]) else C[u, l]
                R[u, l] += Cul + beta[u, l]*A[u, l]
                R[l, u] += Cul*(g[u]/g[l])*np.exp(-(Ek[u]-Ek[l])/Ts)
        M = R.T - np.diag(R.sum(axis=1))
        M[-1, :] = 1.0
        b = np.zeros(n); b[-1] = 1.0
        f[:, it] = np.linalg.solve(M, b)
    return f


def w_fs(g, Ek, A, C, beta, Tv):
    f = se_populations(g, Ek, A, C, beta, Tv)
    W = np.zeros(np.atleast_1d(Tv).size)
    n = len(g)
    for u in range(n):
        for l in range(u):
            if A[u, l] > 0:
                W += f[u]*beta[u, l]*A[u, l]*K_B*(Ek[u]-Ek[l])
    return W


def above_term_coronal(elem, ion, nlev, Tv):
    """Coronal cooling of every channel that LEAVES the ground term, i.e. the
    remainder the SE solution does not carry: sum over .scups transitions whose
    lower level is in the ground term and whose upper level is not, weighted by
    the ground-term Boltzmann fraction of the lower level (the convention of
    cooling_effective(pop='ground_term')).  Built additively rather than by
    subtraction so it cannot go negative: the elvlc level energies used by the
    SE solution and the .scups transition energies differ by up to 3% for the
    FS transitions (N II), which a subtraction would turn into a large relative
    error where the remainder is small."""
    d = ion_dir(elem, ion)
    lev = _level_energy_cm1(read_elvlc(os.path.join(d, f"{elem}_{ion}.elvlc")))
    idx = sorted(lev)[:nlev]
    g = np.array([lev[i][0] for i in idx])
    Ek = np.array([lev[i][1]*HC_OVER_K for i in idx])
    Z = np.sum(g[:, None]*np.exp(-Ek[:, None]/Tv[None, :]), axis=0)
    out = np.zeros_like(Tv)
    from chianti_cooling import RY_ERG
    for tr in read_scups(os.path.join(d, f"{elem}_{ion}.scups")):
        if tr["ll"] not in idx or tr["ul"] in idx:
            continue
        il = idx.index(tr["ll"])
        de = tr["de"]*RY_ERG
        q_lu = COLL_PREF*upsilon(tr, Tv)/(g[il]*np.sqrt(Tv))*np.exp(-de/(K_B*Tv))
        out += (g[il]*np.exp(-Ek[il]/Tv)/Z)*q_lu*de
    return out


def fit_multiexp(lam, Tguess):
    y = np.log(np.maximum(lam, 1e-300)*np.sqrt(T))
    n = len(Tguess)

    def resid(th):
        A = np.exp(th[:n]); Ti = np.exp(th[n:])
        return np.log(np.sum(A[:, None]*np.exp(-Ti[:, None]/T[None, :]),
                             axis=0)) - y

    th0 = np.concatenate([np.full(n, np.log(lam[-1]*np.sqrt(T[-1])/n)),
                          np.log(np.asarray(Tguess, float))])
    s = least_squares(resid, th0, method="lm", max_nfev=200000)
    A = np.exp(s.x[:n]); Ti = np.exp(s.x[n:])
    o = np.argsort(Ti)
    return A[o], Ti[o]


def multiexp(A, Ti, Tv):
    return np.sum(np.asarray(A)[:, None]
                  * np.exp(-np.asarray(Ti)[:, None]/np.atleast_1d(Tv)[None, :]),
                  axis=0)/np.sqrt(np.atleast_1d(Tv))


if __name__ == "__main__":
    store = {}
    for name, (el, ion, nlev, amu) in FS_IONS.items():
        g, Ek, A, ups = ground_term_data(el, ion, nlev)
        print(f"\n=================== {name} "
              f"({nlev}-level ground term, amu = {amu}) ===================")
        for i in range(nlev):
            print(f"  level {i+1}: g = {g[i]:.0f}  E/k = {Ek[i]:.4f} K")
        for u in range(nlev):
            for l in range(u):
                if A[u, l] > 0:
                    lam_um = 1e4*HC_OVER_K/(Ek[u]-Ek[l])*1e4
                    print(f"  A({u+1}->{l+1}) = {A[u,l]:.4e} s^-1"
                          f"   ({1e4*HC_OVER_K/(Ek[u]-Ek[l]):.2f} um)")

        # electron collision strengths
        ufit = {}
        for (l, u), tr in ups.items():
            c, err, uv = fit_upsilon(tr)
            ufit[(l, u)] = (c, uv)
            terms = " + ".join(f"{cc:.8g}*x**{UPS_DEG-k}" if UPS_DEG-k else
                               f"{cc:.8g}" for k, cc in enumerate(c))
            print(f"  log10 Ups({l+1},{u+1}) = {terms}   (max err "
                  f"{err*100:.2f}%)")

        # H collisions
        kfit = {}
        if name in KH_TAB:
            Ttab, ktab = KH_TAB[name]
            for (u, l), kv in ktab.items():
                c, err, lim = fit_kH(Ttab, kv)
                kfit[(u, l)] = (c, lim)
                print(f"  log10 kH({u}->{l}) = {c[2]:.8g} + {c[1]:.8g}*u"
                      f" + {c[0]:.8g}*u^2,  u = log10(T/1e3), T clamped to"
                      f" [{lim[0]:.0f},{lim[1]:.0f}] K   (max err "
                      f"{err*100:.2f}%)")

        # remainder: every channel out of the ground term into a higher term
        rem = above_term_coronal(el, ion, nlev, T)
        Ar, Tr = fit_multiexp(rem, REM_GUESS[name])
        rem_fit = multiexp(Ar, Tr, T)
        mask = rem > 1e-3*rem.max()
        print(f"  remainder fit ({len(Ar)} terms), max err "
              f"{np.abs(rem_fit/rem - 1)[mask].max()*100:.2f}%:")
        for a, t in zip(Ar, Tr):
            print(f"    A = {a:.8e}   Ti = {t:.6g}")
        # how much the remainder depends on the lower-level weighting: the
        # ground-only (thin) versus Boltzmann (dense) choice.  They agree
        # closely because Upsilon_lu is roughly proportional to g_l for the
        # LS-coupled transitions that leave the ground term.
        cor_only = cooling_effective(el, ion, T, pop="coronal")
        d = ion_dir(el, ion)
        from chianti_cooling import RY_ERG
        rem_g = np.zeros_like(T)
        for tr in read_scups(os.path.join(d, f"{el}_{ion}.scups")):
            if tr["ll"] != 1 or tr["ul"] <= nlev:
                continue
            de = tr["de"]*RY_ERG
            rem_g += COLL_PREF*upsilon(tr, T)/(g[0]*np.sqrt(T)) \
                * np.exp(-de/(K_B*T))*de
        m2 = rem > 1e-3*rem.max()
        print(f"  remainder, ground-only vs Boltzmann weighting: max dev "
              f"{np.abs(rem_g/rem - 1)[m2].max()*100:.1f}%")

        # low-density check: the SE ground-term emission must collapse onto the
        # coronal (ground-populated) within-term channel
        ne0 = 1e-8
        C0 = {(u, l): ne0*COLL_PREF*10.0**np.polyval(ufit[(l, u)][0],
                                                     np.log10(T/1e4))
              / (g[u]*np.sqrt(T))
              for u in range(nlev) for l in range(u) if (l, u) in ufit}
        beta1 = np.ones((nlev, nlev))
        cor_fs = np.zeros_like(T)
        for (l, u), tr in ups.items():
            dE = Ek[u]-Ek[l]
            if l != 0:
                continue
            cor_fs += COLL_PREF*upsilon(tr, T)/(g[0]*np.sqrt(T)) \
                * np.exp(-dE/T)*K_B*dE
        print(f"  low-density SE vs coronal within-term channel: max dev "
              f"{np.abs(w_fs(g, Ek, A, C0, beta1, T)/ne0/cor_fs - 1).max()*100:.2f}%")

        # LTE check: SE at huge density vs the exact Boltzmann emission
        Cbig = {k: 1e30*np.ones_like(T) for k in C0}
        Zt = np.sum(g[:, None]*np.exp(-Ek[:, None]/T[None, :]), axis=0)
        lte = np.zeros_like(T)
        for u in range(nlev):
            for l in range(u):
                if A[u, l] > 0:
                    lte += (g[u]*np.exp(-Ek[u]/T)/Zt)*A[u, l]*K_B*(Ek[u]-Ek[l])
        print(f"  high-density SE vs exact LTE emission: max dev "
              f"{np.abs(w_fs(g, Ek, A, Cbig, beta1, T)/lte - 1).max()*100:.2e}%")
        store[name] = dict(g=g, Ek=Ek, A=A, ufit=ufit, kfit=kfit,
                           rem=(Ar, Tr), amu=amu)
