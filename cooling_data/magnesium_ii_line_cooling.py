"""Mg II collisional line cooling with the observed level energies, and its fit.

WHAT THIS BUILDS.  The coronal (n_e -> 0) line cooling of Mg II per
(n_e n_MgII), summed over every electron-impact excitation out of the 3s 2S1/2
ground level into the 32 bound levels of the CHIANTI v11.0.2 model ion (its
other 129 levels lie above the 121267.61 cm^-1 ionization energy and
autoionize, so an excitation into them is not radiated; 0.07 percent of the
sum at 1e5 K, under 1e-5 below 5e4 K),

    Lambda_MgII(T) = sum_u q_1u(T) E_1u            [erg cm^3 s^-1],

with the OBSERVED level energies of mg_2.elvlc in both the Boltzmann factor of
q_1u and the radiated energy E_1u (multilevel_statistical_equilibrium.py,
coronal_line_cooling).  The fit it replaces (cool_MgII_func in Cool_coeff.f90,
from fit_cooling_formulas.py) was made to a table built by
chianti_cooling.cooling_lambda, which puts the .scups header energy in both
places; for 3s-3p that energy is the Burgess-Tully scaling energy 0.3137 /
0.3143 Ry = 4.268 / 4.276 eV, the theoretical level energy, against the
observed 35669.301 / 35760.898 cm^-1 = 4.4224 / 4.4338 eV.

DENSITY.  Every excited level of Mg II decays by a permitted line, so the full
statistical equilibrium equals the coronal sum to 0.2 percent for
n_e <= 1e12 cm^-3 (2 percent at 1e13, 17 percent at 1e14; printed by main()).

FIT.  The h&k doublet keeps the two-level skeleton of the old fit,

    Lambda_hk = 8.629e-6/(2 sqrt T) Ups_hk(T) E_hk exp(-E_hk/kT),
    Ups_hk(T) = a + b ln(1 + T/T0),

with E_hk FIXED at the g-weighted observed 3p 2P energy (2 E_1/2 + 4 E_3/2)/6
= 35730.366 cm^-1 = 4.4300 eV instead of being a free parameter; a single
exponential in place of the two lines of the doublet costs 0.2 percent at
1e3 K (the doublet splitting is 0.011 eV).  The channels above 3p (4s at
8.65 eV, 3d at 8.86 eV, 4p and higher; 0.3 percent of the total at 1e4 K,
4.5 percent at 2e4 K, 39 percent at 1e5 K) enter as one effective term

    Lambda_x = c exp(-T_x/T) / sqrt(T).

Least squares on log10 Lambda over 1e3-1e5 K (201 points).  The fit error is
printed and written to the table header.

OUTPUT.  mg2_line_cooling_observed_energies.txt: the fit coefficients and
the h&k energy in its header (read from there into the Fortran function),
then T, the total, the h&k part, the old fit and the new fit.

RUN:  XUVTOP=<dbase> python3 magnesium_ii_line_cooling.py
"""

import os
import numpy as np
from scipy.optimize import least_squares

from multilevel_statistical_equilibrium import (load_model_ion, line_cooling,
                                                coronal_line_cooling,
                                                ionization_energy_cm1,
                                                restrict_to_bound_levels,
                                                HC_ERG_CM, K_B_ERG)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
COLL_PREF_FORTRAN = 8.629e-6      # the prefactor written in the Fortran form
HC_OVER_K = HC_ERG_CM / K_B_ERG   # cm K


def cool_mgii_fit_old(T):
    """cool_MgII_func as coded in Cool_coeff.f90 before this refit."""
    ups = 16.230522 + 20.344762 * np.log(1.0 + T / 121896.19)
    return 8.629e-6 / (2.0 * np.sqrt(T)) * ups * 6.85184290e-12 * np.exp(-49627.696 / T)


def hk_energy_cm1(model):
    """g-weighted observed energy of the 3p 2P term (levels 2 and 3)."""
    k = [int(np.where(model["idx"] == i)[0][0]) for i in (2, 3)]
    g, e = model["g"][k], model["e_cm1"][k]
    return float((g * e).sum() / g.sum())


def fit_form(T, p, e_hk_cm1):
    a, b, t0, c, tx = p
    e_erg = e_hk_cm1 * HC_ERG_CM
    t_hk = e_hk_cm1 * HC_OVER_K
    ups = a + b * np.log(1.0 + T / t0)
    return (COLL_PREF_FORTRAN / (2.0 * np.sqrt(T)) * ups * e_erg * np.exp(-t_hk / T)
            + c * np.exp(-tx / T) / np.sqrt(T))


def fit_mgii(T, lam, e_hk_cm1):
    y = np.log10(lam)

    def resid(th):
        p = [th[0], th[1], np.exp(th[2]), np.exp(th[3]), np.exp(th[4])]
        return np.log10(fit_form(T, p, e_hk_cm1)) - y

    th0 = [16.0, 20.0, np.log(1.2e5), np.log(1.0e-15), np.log(1.1e5)]
    sol = least_squares(resid, th0, method="lm", max_nfev=200000,
                        xtol=1e-15, ftol=1e-15)
    th = sol.x
    return [th[0], th[1], np.exp(th[2]), np.exp(th[3]), np.exp(th[4])]


def main():
    full = load_model_ion("mg", 2)
    model = restrict_to_bound_levels(full, ionization_energy_cm1(12, 2))
    T = np.logspace(3.0, 5.0, 201)
    total = coronal_line_cooling(model, T)
    hk = coronal_line_cooling(model, T, upper_levels={2, 3})
    e_hk = hk_energy_cm1(model)
    p = fit_mgii(T, total, e_hk)
    new = fit_form(T, p, e_hk)
    old = cool_mgii_fit_old(T)
    rel = np.abs(new / total - 1.0)
    wind = (T >= 2.0e3) & (T <= 3.0e4)
    err_full, err_wind = rel.max(), rel[wind].max()

    path = os.path.join(OUTDIR, "mg2_line_cooling_observed_energies.txt")
    with open(path, "w") as fh:
        fh.write("# Mg II coronal line cooling per (n_e n_MgII) [erg cm^3 s^-1], CHIANTI v11.0.2\n")
        fh.write("# observed .elvlc level energies in the Boltzmann factor and the radiated energy\n")
        fh.write("# built by magnesium_ii_line_cooling.py (multilevel_statistical_equilibrium.py)\n")
        fh.write(f"# h&k effective energy (g-weighted observed 3p 2P) = {e_hk:.3f} cm^-1\n")
        fh.write(f"# fit a, b, T0, c, T_x = {p[0]:.8g} {p[1]:.8g} {p[2]:.8g} {p[3]:.8e} {p[4]:.9g}\n")
        fh.write(f"# h&k E_erg, T_hk = {e_hk*HC_ERG_CM:.8e} {e_hk*HC_OVER_K:.3f}\n")
        fh.write(f"# fit max rel. error: {100*err_full:.2f} % (1e3-1e5 K), {100*err_wind:.2f} % (2e3-3e4 K)\n")
        fh.write("# columns: T_K  Lambda_all  Lambda_hk  old_fit_scups_energy  new_fit\n")
        for i in range(T.size):
            fh.write(f"{T[i]:12.5e} {total[i]:14.6e} {hk[i]:14.6e} {old[i]:14.6e} {new[i]:14.6e}\n")
    print(f"wrote {path}")
    print(f"bound levels {model['idx'].size} of {full['idx'].size}; coronal sum over "
          f"bound upper levels / over all at 5e4, 1e5 K: "
          + " ".join(f"{x:.6f}" for x in coronal_line_cooling(model, [5e4, 1e5])
                     / coronal_line_cooling(full, [5e4, 1e5])))
    print(f"h&k g-weighted observed energy {e_hk:.3f} cm^-1 = "
          f"{e_hk*HC_ERG_CM/1.602176634e-12:.4f} eV")
    print(f"fit: a={p[0]:.8g} b={p[1]:.8g} T0={p[2]:.8g} c={p[3]:.8e} T_x={p[4]:.6g}")
    print(f"fit max rel. error {100*err_full:.2f} % over 1e3-1e5 K, "
          f"{100*err_wind:.2f} % over 2e3-3e4 K")
    print("  T [K]   old/new-table  old/hk-table  total/hk   newfit/table")
    for tc in (1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 3e4, 5e4, 1e5):
        i = int(np.argmin(np.abs(T - tc)))
        print(f"  {T[i]:7.0f}   {old[i]/total[i]:8.4f}      {old[i]/hk[i]:8.4f}    "
              f"{total[i]/hk[i]:7.4f}    {new[i]/total[i]:8.5f}")
    # density check: full statistical equilibrium / coronal
    tt = np.array([2e3, 5e3, 1e4, 2e4, 3e4])
    nn = np.array([1e8, 1e10, 1e12, 1e13, 1e14])
    se = line_cooling(model, tt, nn, 0)
    cor = coronal_line_cooling(model, tt)
    print("full statistical equilibrium / coronal, rows T, columns n_e =", nn)
    for i in range(tt.size):
        print(f"  {tt[i]:7.0f} " + " ".join(f"{x:.4f}" for x in se[i] / cor[i]))


if __name__ == "__main__":
    main()
