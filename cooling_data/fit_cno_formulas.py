"""CHIANTI-based closed-form cooling formulas for C I/II, N I/II, O I/II.

The current C/O coolants in Cool_coeff.f90 are AIOLOS analytic fits
(photochem.cpp port, no in-code citation, likely Cloudy-calibrated);
N I/N II have no line cooling at all (return 0).  This script

  1. computes the CHIANTI v11.0.2 optically-thin line-cooling curves
     Lambda(T) per (n_e * n_ion) for the six ions,
  2. fits each with the same effective two-level multi-exponential form
     used for Fe II,   Lambda(T) = T^{-1/2} sum_i A_i exp(-T_i/T),
  3. compares the new fits and the AIOLOS formulas against CHIANTI.

Population model: lower levels Boltzmann-distributed over the ground-term
fine structure (within 0.4 eV: C I 3P, C II 2P, N II 3P, O I 3P); the
4S3/2 ground of N I and O II is a single level, so 'ground_term' reduces
to coronal there.  Every collisional excitation is assumed to radiate
(optically thin, low-density limit).

CAVEAT (document, like Fe II coronal): the fine-structure floor (e.g.
[C II] 158um, [O I] 63um) has very low critical density (n_e ~ 10-1e5
cm^-3) and saturates at the dense atmosphere base; both the AIOLOS
constant floors and these coronal fits overestimate it there.  The
optional use_2lev_cool treatment handles the saturation separately.

Writes:  cno_formula_coefficients.txt (Fortran-ready coefficients)
Run:     python3 fit_cno_formulas.py
"""

import numpy as np
from scipy.optimize import least_squares
from chianti_cooling import cooling_effective

T = np.logspace(3.0, 5.0, 201)

IONS = {
    "CI":  ("c", 1),
    "CII": ("c", 2),
    "NI":  ("n", 1),
    "NII": ("n", 2),
    "OI":  ("o", 1),
    "OII": ("o", 2),
}


def aiolos(name, T):
    """Current Cool_coeff.f90 analytic coolants (AIOLOS port)."""
    if name == "CI":
        return 1.0e-24 + 3.1e-20*np.exp(-15162.0/T)*(1.0 + (T/2.0e4)**1.5)
    if name == "CII":
        return 1.5e-23 + 3.1e-20*np.exp(-45162.0/T)*(1.0 + (T/0.75e4)**1.5)
    if name == "OI":
        return 5.5e-24 + 1.1e-20*np.exp(-30162.0/T)*(1.0 + (T/0.75e4)**0.5)
    if name == "OII":
        return 5.1e-20*np.exp(-35162.0/T)*(1.0 + (T/0.75e4)**0.5)
    return np.zeros_like(T)  # N I, N II


def fit_multiexp(lam, nterm, Tguess):
    """Lambda = T^-1/2 sum A_i exp(-T_i/T); LM on log Lambda."""
    y = np.log(lam*np.sqrt(T))

    def resid(th):
        A = np.exp(th[:nterm]); Ti = np.exp(th[nterm:])
        return np.log(np.sum(A[:, None]*np.exp(-Ti[:, None]/T[None, :]),
                             axis=0)) - y

    th0 = np.concatenate([np.full(nterm, np.log(lam[-1]*np.sqrt(T[-1])/nterm)),
                          np.log(np.asarray(Tguess, float))])
    sol = least_squares(resid, th0, method="lm", max_nfev=60000)
    A = np.exp(sol.x[:nterm]); Ti = np.exp(sol.x[nterm:])
    o = np.argsort(Ti)
    return A[o], Ti[o]


def evaluate(A, Ti, T):
    return np.sum(A[:, None]*np.exp(-Ti[:, None]/T[None, :]), axis=0)/np.sqrt(T)


GUESSES = {  # rough dE ladders [K]; FS / forbidden optical / UV
    "CI":  [300.0, 15000.0, 90000.0],
    "CII": [200.0, 50000.0, 120000.0],
    "NI":  [28000.0, 42000.0, 120000.0],
    "NII": [300.0, 22000.0, 60000.0, 150000.0],
    "OI":  [400.0, 23000.0, 110000.0],
    "OII": [39000.0, 58000.0, 150000.0],
}

results, chianti, fits = {}, {}, {}
print(f"{'ion':5s} {'nterm':>5s} {'maxerr_fit':>11s} {'maxerr_relevant':>16s}"
      f" {'AIOLOS/CHIANTI @1e4K':>22s} {'@2e4K':>8s}")
for name, (el, st) in IONS.items():
    lam = cooling_effective(el, st, T, pop="ground_term")
    chianti[name] = lam
    best = None
    for nterm in (len(GUESSES[name]), len(GUESSES[name]) + 1,
                  len(GUESSES[name]) + 2):
        g = GUESSES[name] + [2.0e5]*(nterm - len(GUESSES[name]))
        try:
            A, Ti = fit_multiexp(lam, nterm, g)
        except Exception:
            continue
        fit = evaluate(A, Ti, T)
        err = np.abs(fit/lam - 1.0).max()
        if best is None or err < best[0]:
            best = (err, A, Ti, fit)
    err, A, Ti, fit = best
    fits[name] = (A, Ti)
    mask = lam > 1e-3*lam.max()
    rel_err = np.abs(fit/lam - 1.0)[mask].max()
    i4 = np.argmin(abs(T - 1e4)); i2 = np.argmin(abs(T - 2e4))
    av = aiolos(name, T)
    r4 = av[i4]/lam[i4] if lam[i4] > 0 else np.nan
    r2 = av[i2]/lam[i2] if lam[i2] > 0 else np.nan
    print(f"{name:5s} {len(A):5d} {err*100:10.2f}% {rel_err*100:15.2f}%"
          f" {r4:22.2f} {r2:8.2f}")
    results[name] = dict(A=A, Ti=Ti)

print("\n=== fitted coefficients:  Lambda = T^-1/2 * sum_i A_i exp(-T_i/T) ===")
with open("cno_formula_coefficients.txt", "w") as fh:
    fh.write("# CHIANTI v11.0.2 ground-term-Boltzmann coronal cooling fits\n")
    fh.write("# Lambda(T) = T^-1/2 * sum_i A_i*exp(-T_i/T)  [erg cm^3 s^-1]\n")
    fh.write("# per (n_e * n_ion); fit range 1e3-1e5 K\n")
    for name, r in results.items():
        line = f"{name:5s} " + "  ".join(
            f"A={a:.8e} T={t:.6g}" for a, t in zip(r["A"], r["Ti"]))
        print(" ", line)
        fh.write(line + "\n")
print("\nwrote cno_formula_coefficients.txt")
