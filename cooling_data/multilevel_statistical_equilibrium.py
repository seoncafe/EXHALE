"""Line cooling from the multilevel statistical equilibrium of a CHIANTI ion.

WHAT THIS COMPUTES.  For one ion, the fractional level populations f_i(T, n_e)
under electron-impact excitation and de-excitation and spontaneous radiative
decay (optically thin), and from them the line cooling

    Lambda(T, n_e) = ( sum_{u>l} f_u A_ul E_ul ) / n_e     [erg cm^3 s^-1],

the radiated power of one ion divided by n_e, i.e. the coefficient that the
cooling assembly multiplies by n_e n_ion.  In the n_e -> 0 limit this tends to
the coronal sum over excitations out of the ground level,

    Lambda_coronal(T) = sum_u q_1u(T) E_1u ,

which coronal_line_cooling() evaluates directly.

LEVEL ENERGIES.  Every energy that enters a rate or a radiated power is the
OBSERVED level energy of the CHIANTI .elvlc file (the theoretical energy only
where no observed one exists): the Boltzmann factor of each excitation rate,
the detailed-balance ratio of excitation to de-excitation, and the energy
E_ul carried by each photon.  The energy in the header of each .scups entry
is the Burgess & Tully (1992) scaling energy that CHIANTI used to scale the
temperature axis of that transition's collision strength; it is used here
for that descaling and nowhere else.  For Mg II the .scups scaling energy of
3s-3p is 4.27 eV against 4.43 eV observed, which is a factor 1.8 in
exp(-dE/kT) at 3000 K, and for the O II 4S-2D excitation it is 14 percent
above the observed gap (CHIANTI v11.0.2).

RATES.  With Upsilon_lu(T) the Maxwellian-averaged collision strength,

    q_ul = C_e Upsilon_lu / (g_u sqrt(T)),
    q_lu = (g_u/g_l) q_ul exp(-(E_u - E_l)/kT),
    C_e  = h^2 / ((2 pi m_e)^{3/2} k^{1/2}) = 8.6291e-6 cm^3 s^-1 K^{1/2},

where "u" is the level of higher OBSERVED energy, so detailed balance holds
exactly for every pair even where the observed and theoretical orderings of
two nearly degenerate levels differ.  Radiative decay rates A_ul are the
.wgfa values, summed over duplicate entries of one pair.  Electron collisions
only.  Of the ions this module is used for, C II and N II carry proton
(.psplups) collision strengths in CHIANTI v11.0.2; they are left out because
the proton density is not an argument of these tables (ChiantiPy takes it
from CHIANTI's coronal ionization equilibrium, which is not the gas the
tables are used for).  Neutral-hydrogen impact is not in CHIANTI.

VALIDITY.  Optically thin lines, no radiative excitation by an external
field, no recombination or ionization feeding of excited levels (their rates
are orders of magnitude below the collisional and radiative rates of bound
levels at the densities considered), steady state.  No autoionization: for
an ion whose model has levels above the ionization energy (Mg II: 129 of
161), restrict_to_bound_levels removes them before use.  The model atom is the
CHIANTI one; transitions it omits are absent here too (for Mg I there are no
electron collision strengths among the 3s3p 3P_J fine-structure levels).

WHY NOT ChiantiPy.  ChiantiPy 0.15.2 `ion()` sets up the level data only for
ions listed in $XUVTOP/masterlist/masterlist.ions; mg_1 is not listed in
CHIANTI v11.0.2, so `populate()` raises AttributeError ('Nlvls').  Calling
`ion.setup()` explicitly works around it.  This module reads the same three files with no ChiantiPy code
and descales every collision strength with chianti_cooling.upsilon (the
natural spline of CHIANTI's DESCALE_SCUPS.PRO), and it is the solver of
metal_cooling_density_resolved.py and magnesium_ii_line_cooling.py.
validate() compares it with ChiantiPy populate, with the descaling made
equal for that comparison (ChiantiPy 0.15.2 descales with a not-a-knot
spline).

RUN:  XUVTOP=<dbase> python3 multilevel_statistical_equilibrium.py
      prints the comparison with ChiantiPy populate for ca_2, mg_1, n_1, o_2, mg_2.
"""

import os
import numpy as np

from chianti_cooling import read_elvlc, read_wgfa, read_scups, upsilon, ion_dir

# CODATA 2018 exact / recommended constants (cgs)
H_PLANCK = 6.62607015e-27        # erg s
K_B_ERG = 1.380649e-16           # erg / K
M_E = 9.1093837015e-28           # g
C_LIGHT = 2.99792458e10          # cm / s
HC_ERG_CM = H_PLANCK * C_LIGHT   # erg cm: E[erg] = E[cm^-1] * HC_ERG_CM
C_E = H_PLANCK**2 / ((2.0 * np.pi * M_E)**1.5 * np.sqrt(K_B_ERG))  # 8.6291e-6


def observed_level_energies(lev):
    """{idx: (g, E_obs, E_th)} -> arrays over sorted indices: idx, g, E_cm1.

    E is the observed energy where the .elvlc file gives one (E_obs > 0, or
    the ground level), else the theoretical energy.
    """
    idx = np.array(sorted(lev.keys()))
    g = np.array([lev[i][0] for i in idx], dtype=float)
    e = np.array([lev[i][1] if (lev[i][1] > 0.0 or i == idx[0]) else lev[i][2]
                  for i in idx], dtype=float)
    return idx, g, e


def load_model_ion(elem, ion, dbase_dir=None):
    """Read the CHIANTI .elvlc / .wgfa / .scups files of one ion.

    Returns a dict with the level indices, statistical weights, observed
    energies [cm^-1], the radiative matrix A[u, l] (s^-1, u the upper, summed
    over duplicate .wgfa rows), and the .scups transitions (level pairs as
    array positions).  Only levels that appear in the .elvlc file are kept.
    """
    d = dbase_dir or ion_dir(elem, ion)
    name = f"{elem}_{ion}"
    lev = read_elvlc(os.path.join(d, name + ".elvlc"))
    idx, g, e = observed_level_energies(lev)
    pos = {int(i): k for k, i in enumerate(idx)}
    n = idx.size
    arad = np.zeros((n, n))
    n_dropped = 0
    for ll, ul, wl, gf, a in read_wgfa(os.path.join(d, name + ".wgfa")):
        if a <= 0.0 or ll not in pos or ul not in pos:
            continue
        i, j = pos[ll], pos[ul]
        up, lo = (j, i) if e[j] > e[i] else (i, j)
        if e[up] == e[lo]:
            n_dropped += 1
            continue
        arad[up, lo] += a
    trans = []
    for tr in read_scups(os.path.join(d, name + ".scups")):
        if tr["ll"] in pos and tr["ul"] in pos:
            trans.append((pos[tr["ll"]], pos[tr["ul"]], tr))
    return dict(name=name, idx=idx, g=g, e_cm1=e, arad=arad, trans=trans,
                n_rad_dropped=n_dropped)


def ionization_energy_cm1(z, stage):
    """Ionization energy [cm^-1] of the ion (atomic number z, spectroscopic
    stage, 1 = neutral) from $XUVTOP/ip/chianti.ip."""
    path = os.path.join(os.environ.get("XUVTOP", ""), "ip", "chianti.ip")
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) >= 3 and t[0] == str(z) and t[1] == str(stage):
                return float(t[2])
    raise KeyError(f"no ionization energy for Z={z}, stage {stage} in {path}")


def restrict_to_bound_levels(model, ip_cm1):
    """The same model without the levels above the ionization energy.

    Those levels (inner-shell excitations, e.g. Mg II 2p5 3s2 at 50 eV)
    autoionize far faster than they radiate, so an excitation into them is
    not radiated; this module carries no autoionization rates, so they are
    removed rather than left to decay radiatively."""
    keep = np.where(model["e_cm1"] < ip_cm1)[0]
    new_pos = {int(k): n for n, k in enumerate(keep)}
    trans = [(new_pos[a], new_pos[b], tr) for a, b, tr in model["trans"]
             if a in new_pos and b in new_pos]
    out = dict(model)
    out.update(idx=model["idx"][keep], g=model["g"][keep], e_cm1=model["e_cm1"][keep],
               arad=model["arad"][np.ix_(keep, keep)], trans=trans)
    return out


def collision_matrix(model, T, extra_upsilon=None):
    """Electron-impact rate coefficients per unit n_e, Q[j, i] = rate i -> j.

    T is a 1-D array; returns Q of shape (nT, n, n) [cm^3 s^-1].
    extra_upsilon, if given, is a list of (pos_a, pos_b, Upsilon) with a
    temperature-independent collision strength added for that pair (used only
    to measure the sensitivity to a transition the model atom omits).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    g, e = model["g"], model["e_cm1"]
    n = g.size
    q = np.zeros((T.size, n, n))
    pairs = [(a, b, upsilon(tr, T)) for a, b, tr in model["trans"]]
    if extra_upsilon:
        pairs += [(a, b, np.full(T.size, float(u))) for a, b, u in extra_upsilon]
    sqrt_t = np.sqrt(T)
    for a, b, ups in pairs:
        up, lo = (b, a) if e[b] >= e[a] else (a, b)
        de_erg = (e[up] - e[lo]) * HC_ERG_CM
        q_ul = C_E * ups / (g[up] * sqrt_t)
        q_lu = q_ul * (g[up] / g[lo]) * np.exp(-de_erg / (K_B_ERG * T))
        q[:, lo, up] += q_ul
        q[:, up, lo] += q_lu
    return q


def level_populations(model, T, ne, extra_upsilon=None):
    """Fractional populations f[iT, iN, level] on the T x n_e grid (sum = 1)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ne = np.atleast_1d(np.asarray(ne, dtype=float))
    q = collision_matrix(model, T, extra_upsilon)
    arad = model["arad"]                     # arad[u, l]: rate u -> l
    rad_in = arad.T                          # rad_in[l, u] = rate u -> l
    n = arad.shape[0]
    f = np.zeros((T.size, ne.size, n))
    diag = np.diag_indices(n)
    rhs = np.zeros(n)
    rhs[0] = 1.0
    for it in range(T.size):
        for ie in range(ne.size):
            m = ne[ie] * q[it] + rad_in      # m[j, i] = rate i -> j
            a = m.copy()
            a[diag] -= m.sum(axis=0)         # loss of level i: sum_j m[j, i]
            a[0, :] = 1.0                    # normalization replaces row 0
            sol = np.linalg.solve(a, rhs)
            f[it, ie] = np.clip(sol, 0.0, None)
    return f


def line_cooling(model, T, ne, n_ground_term=0, extra_upsilon=None,
                 return_populations=False):
    """Lambda(T, n_e) = sum_{u, l} f_u A_ul E_ul / n_e  [erg cm^3 s^-1].

    Only decays whose upper level lies outside the first n_ground_term levels
    (energy order) are summed, the convention of the tables in Cool_coeff.f90
    (n_ground_term = 0 or 1 for an ion with a single ground level).
    Returns an array of shape (nT, nN).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ne = np.atleast_1d(np.asarray(ne, dtype=float))
    f = level_populations(model, T, ne, extra_upsilon)
    e = model["e_cm1"]
    arad = model["arad"].copy()
    order = np.argsort(e)
    ground_term = order[:n_ground_term]
    arad[ground_term, :] = 0.0
    power = arad * ((e[:, None] - e[None, :]) * HC_ERG_CM)   # [u, l] erg/s
    lam = np.einsum("tnu,u->tn", f, power.sum(axis=1)) / ne[None, :]
    return (lam, f) if return_populations else lam


def coronal_line_cooling(model, T, upper_levels=None):
    """n_e -> 0 limit: sum over excitations out of the ground level of q_1u E_1u.

    upper_levels optionally restricts the sum to CHIANTI level indices.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    g, e, idx = model["g"], model["e_cm1"], model["idx"]
    lo = int(np.argmin(e))
    total = np.zeros_like(T)
    for a, b, tr in model["trans"]:
        if lo not in (a, b):
            continue
        up = b if a == lo else a
        if upper_levels is not None and int(idx[up]) not in upper_levels:
            continue
        de_erg = (e[up] - e[lo]) * HC_ERG_CM
        ups = upsilon(tr, T)
        total += C_E * ups / (g[lo] * np.sqrt(T)) * np.exp(-de_erg / (K_B_ERG * T)) * de_erg
    return total


def chiantipy_line_cooling(name, T, ne, n_ground_term=0):
    """The same sum from ChiantiPy populate (the route of
    metal_cooling_density_resolved.py), for the validation below.  T and ne are
    equal-length arrays of (T, n_e) pairs; ion.setup() is called explicitly so
    ions missing from the CHIANTI masterlist (mg_1) are set up as well."""
    import ChiantiPy.core as ch
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ne = np.atleast_1d(np.asarray(ne, dtype=float))
    ion = ch.ion(name, temperature=T, eDensity=ne)
    if not hasattr(ion, "Nlvls"):
        ion.setup()
    ion.populate()
    f = np.asarray(ion.Population["population"])
    upper = np.asarray(ion.Wgfa["lvl2"])
    avalue = np.asarray(ion.Wgfa["avalue"])
    wvl = np.asarray(ion.Wgfa["wvl"])
    keep = (wvl != 0.0) & (avalue > 0.0) & (upper > n_ground_term)
    e_ul = HC_ERG_CM / (np.abs(wvl[keep]) * 1.0e-8)
    return (f[:, upper[keep] - 1] * (avalue[keep] * e_ul)[None, :]).sum(axis=1) / ne


def validate():
    """Compare this solver with ChiantiPy populate on a (T, n_e) grid.

    ChiantiPy 0.15.2 descales the .scups data with FITPACK splrep (s = 0),
    the not-a-knot cubic spline, so for this comparison only the shared
    descaling spline is set to the same not-a-knot spline; what is compared
    is then the equilibrium solution, not the interpolation.  Ions with a
    .psplups file (C II, N II) are not in the list, since ChiantiPy adds
    proton impact there.
    """
    import chianti_cooling
    from scipy.interpolate import CubicSpline
    natural = chianti_cooling.scaled_upsilon_spline
    chianti_cooling.scaled_upsilon_spline = lambda xs, ys: CubicSpline(xs, ys)
    try:
        compare_with_chiantipy_populate()
    finally:
        chianti_cooling.scaled_upsilon_spline = natural


def compare_with_chiantipy_populate():
    logt = np.linspace(3.0, 5.0, 11)
    logne = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
    tt = np.repeat(10.0**logt, logne.size)
    nn = np.tile(10.0**logne, logt.size)
    for elem, ion, ngt in [("ca", 2, 1), ("mg", 1, 1), ("n", 1, 1), ("o", 2, 1),
                           ("mg", 2, 1)]:
        model = load_model_ion(elem, ion)
        mine = line_cooling(model, 10.0**logt, 10.0**logne, ngt)
        ref = chiantipy_line_cooling(f"{elem}_{ion}", tt, nn, ngt).reshape(
            logt.size, logne.size)
        good = ref > 1.0e-40
        r = mine[good] / ref[good]
        cor = coronal_line_cooling(model, 10.0**logt)
        print(f"{elem}_{ion}: {model['idx'].size} levels, "
              f"{len(model['trans'])} collision pairs; this solver / ChiantiPy "
              f"over logT 3-5, log ne 0-14 (Lambda > 1e-40): "
              f"min {r.min():.5f} median {np.median(r):.5f} max {r.max():.5f}; "
              f"SE(n_e=1)/coronal at 1e4 K "
              f"{mine[5, 0]/cor[5]:.5f}")


if __name__ == "__main__":
    validate()
