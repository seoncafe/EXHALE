"""Fe I line cooling in the two density limits of its lower-level populations.

The Fe I coefficient in Cool_coeff.f90 (cool_logL_FeI, built by
fe1_cooling.py) takes the lower levels of its 828 permitted lines in a
Boltzmann distribution over the even-parity metastable manifold (17 levels
below the first odd level z7D, 19351 cm^-1) at every density.  That is the
limit n_e -> infinity of the metastable populations (collisions faster than
their forbidden radiative decays).  The opposite limit, n_e -> 0, puts the
whole population in the a5D4 ground level, because the metastables are
radiatively drained (slowly, by forbidden lines) faster than they are
collisionally refilled.  The true coefficient lies between the two and depends
on n_e through the critical densities of the metastable levels.

This script evaluates, with the same lines, f-values and Van Regemorter
compact rate (Huang et al. 2023, Eq. 10) as fe1_cooling.py:

    boltzmann_manifold  Boltzmann over the 17 even levels (the coded table)
    ground_term         Boltzmann over the a5D term only (5 levels, <= 978 cm^-1),
                        the limit in which the fine structure is thermalized
                        but the higher metastables (a5F at 0.86 eV and above)
                        are empty
    ground_level        a5D4 only (the n_e -> 0 limit)

and prints the ratios, the check that boltzmann_manifold reproduces
cool_logL_FeI, and the proposed comment for cool_logL_FeI.

WHY NO n_e-DEPENDENT TABLE.  An n-level solve needs, besides the permitted
rates, (i) the forbidden radiative rates that drain the metastables and (ii)
the collision strengths that couple them.  fe1_nist_lines.tsv covers 1500 to
13000 A only, so the a5D fine-structure lines (24.0, 34.7 um) and the
a5F-a5D lines near 1.4 um are not in it (103 M1/E2 rows are, from higher
levels); and no electron (or H-atom) collision strength for a forbidden Fe I
transition exists in CHIANTI v11.0.2 (no Fe I model atom) or anywhere in this
tree.  The Van Regemorter approximation applies to permitted lines only.
Published partial data exist (Hollenbach & McKee 1989 for the fine-structure
and 1.36/1.44 um lines, which Huang et al. 2023 used) but not for the higher
metastable terms, so the populations between the two limits cannot be
computed from data available here.

RUN:  python3 fe1_level_population_limits.py
"""

import os
import re
import numpy as np

from chianti_cooling import van_regemorter_clu_compact, EV_ERG, HC_OVER_K
from fe1_cooling import load_levels, load_lines, cooling_FeI, ECUT_CM1_DEFAULT

OUTDIR = os.path.dirname(os.path.abspath(__file__))
A5D_TOP_CM1 = 978.074        # a5D0, the highest level of the ground term


def cooling_with_population(T, lower_set_cm1):
    """Fe I cooling per (n_e n_FeI) with the lower levels whose energies are in
    lower_set_cm1 in a Boltzmann distribution among themselves (every other
    level empty).  Same lines and rates as fe1_cooling.cooling_FeI."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    elev, glev = load_levels()
    sel = np.isin(np.round(elev, 3), np.round(lower_set_cm1, 3))
    elev, glev = elev[sel], glev[sel]
    u = np.zeros_like(T)
    for e, g in zip(elev, glev):
        u += g * np.exp(-e * HC_OVER_K / T)
    ei, de_ev, gi, f_lu = load_lines()
    keep = np.isin(np.round(ei, 3), np.round(elev, 3))
    ei, de_ev, gi, f_lu = ei[keep], de_ev[keep], gi[keep], f_lu[keep]
    total = np.zeros_like(T)
    for k in range(ei.size):
        f_l = gi[k] * np.exp(-ei[k] * HC_OVER_K / T) / u
        total += f_l * van_regemorter_clu_compact(T, de_ev[k], f_lu[k]) * (de_ev[k] * EV_ERG)
    return total, int(ei.size)


def coded_table():
    """cool_logL_FeI as written in Cool_coeff.f90 (41 values, logT 3..5)."""
    path = os.path.join(OUTDIR, "..", "src", "modules", "radiation", "Cool_coeff.f90")
    s = open(path).read()
    m = re.search(r"cool_logL_FeI\(NCOOLT\)\s*=\s*\[(.*?)\]", s, re.S)
    return np.array([float(x) for x in m.group(1).replace("&", " ").replace("d0", "")
                     .replace(",", " ").split()])


def main():
    logt = np.arange(3.0, 5.0 + 1e-9, 0.05)
    T = 10.0**logt
    elev, glev = load_levels()
    manifold = elev[elev < ECUT_CM1_DEFAULT]
    lam_b, nline_b, nlev_b = cooling_FeI(T, ECUT_CM1_DEFAULT)
    lam_m, nline_m = cooling_with_population(T, manifold)
    lam_t, nline_t = cooling_with_population(T, elev[elev <= A5D_TOP_CM1 + 1e-3])
    lam_g, nline_g = cooling_with_population(T, np.array([0.0]))
    coded = coded_table()
    dlog = np.abs(np.log10(lam_b) - coded)
    print(f"boltzmann_manifold: {nlev_b} even levels, {nline_b} lines; "
          f"max |log10(fe1_cooling) - cool_logL_FeI| = {dlog.max():.2e}")
    print(f"same through cooling_with_population: max rel. diff "
          f"{np.abs(lam_m/lam_b - 1).max():.1e}")
    print(f"ground_term: {nline_t} lines; ground_level: {nline_g} lines")
    print("   T [K]   ground_level/Boltzmann   ground_term/Boltzmann")
    rows = []
    for tc in (1e3, 2e3, 3e3, 5e3, 7e3, 1e4, 2e4, 3e4, 1e5):
        i = int(np.argmin(np.abs(T - tc)))
        rows.append((T[i], lam_g[i] / lam_b[i], lam_t[i] / lam_b[i]))
        print(f"  {T[i]:7.0f}   {lam_g[i]/lam_b[i]:10.3f}              {lam_t[i]/lam_b[i]:8.3f}")
    band = (T >= 2.0e3) & (T <= 2.0e4)
    r = lam_g[band] / lam_b[band]
    print(f"ground_level/Boltzmann over 2e3-2e4 K: min {r.min():.3f} max {r.max():.3f}")
    path = os.path.join(OUTDIR, "fe1_level_population_limits.txt")
    with open(path, "w") as fh:
        fh.write("# Fe I line cooling per (n_e n_FeI) [erg cm^3 s^-1] in the limits of the\n")
        fh.write("# lower-level populations; NIST ASD f-values + Van Regemorter compact form\n")
        fh.write("# (fe1_cooling.py lines), built by fe1_level_population_limits.py\n")
        fh.write("# columns: T_K  boltzmann_manifold(coded)  ground_term_a5D  ground_level_a5D4\n")
        for i in range(T.size):
            fh.write(f"{T[i]:12.5e} {lam_b[i]:14.6e} {lam_t[i]:14.6e} {lam_g[i]:14.6e}\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
