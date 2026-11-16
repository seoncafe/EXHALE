"""Error made by using the CHIANTI .scups header energy as the transition energy.

chianti_cooling.cooling_lambda and cooling_effective put the energy written in
the header of each .scups entry (the Burgess & Tully scaling energy of that
collision strength) into the excitation Boltzmann factor and into the radiated
energy.  Physically both are the observed level gap of the .elvlc file (the
theoretical gap only where no observed energy exists).  This script prints,
for every ion whose cooling the scripts in this directory fit, the ratio

    Lambda(.scups energy) / Lambda(observed energy)

of the quantity each coefficient was fitted to (ground level for Mg I, Mg II,
Ca II, Fe II, H I, He I, He II; Boltzmann ground term for C, N, O, as in
fit_cno_formulas.py), and the range of the energy ratio over the transitions.  The
collision strength is descaled with the .scups energy in both cases.

RUN:  XUVTOP=<dbase> python3 scaling_energy_error.py
"""

import os
import numpy as np

from chianti_cooling import (read_elvlc, read_scups, upsilon, _level_energy_cm1,
                             ion_dir, RY_ERG, K_B_ERG, HC_OVER_K)
from multilevel_statistical_equilibrium import C_E, HC_ERG_CM

T = np.array([1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 3e4, 5e4, 1e5])


def coronal_sum(elem, ion, lower_pop, upper=None, energy="observed"):
    """sum over .scups transitions l->u of f_l q_lu dE_lu, with f_l given by
    lower_pop(T) -> {level: fraction array}."""
    d = ion_dir(elem, ion)
    lev = _level_energy_cm1(read_elvlc(os.path.join(d, f"{elem}_{ion}.elvlc")))
    pop = lower_pop(lev)
    total = np.zeros_like(T)
    ratios = []
    for tr in read_scups(os.path.join(d, f"{elem}_{ion}.scups")):
        if tr["ll"] not in pop or (upper is not None and tr["ul"] not in upper):
            continue
        de_obs = (lev[tr["ul"]][1] - lev[tr["ll"]][1]) * HC_ERG_CM
        if de_obs <= 0.0:
            continue
        ratios.append(tr["de"] * RY_ERG / de_obs)
        de = tr["de"] * RY_ERG if energy == "scups" else de_obs
        q = C_E / (lev[tr["ll"]][0] * np.sqrt(T)) * upsilon(tr, T) * np.exp(-de / (K_B_ERG * T))
        total += pop[tr["ll"]] * q * de
    return total, np.array(ratios)


def ground_level(lev):
    return {1: np.ones_like(T)}


def ground_term(n):
    def pop(lev):
        z = sum(lev[i][0] * np.exp(-lev[i][1] * HC_OVER_K / T) for i in range(1, n + 1))
        return {i: lev[i][0] * np.exp(-lev[i][1] * HC_OVER_K / T) / z for i in range(1, n + 1)}
    return pop


CASES = [("Mg I 2853", "mg", 1, ground_level, {5}),
         ("Mg II all channels", "mg", 2, ground_level, None),
         ("Ca II H&K", "ca", 2, ground_level, {4, 5}),
         ("Fe II coronal", "fe", 2, ground_level, None),
         ("H I ground", "h", 1, ground_level, None),
         ("He I ground", "he", 1, ground_level, None),
         ("He II ground", "he", 2, ground_level, None),
         ("C I ground term", "c", 1, ground_term(3), None),
         ("C II ground term", "c", 2, ground_term(2), None),
         ("N I ground term", "n", 1, ground_term(1), None),
         ("N II ground term", "n", 2, ground_term(3), None),
         ("O I ground term", "o", 1, ground_term(3), None),
         ("O II ground term", "o", 2, ground_term(1), None)]


def main():
    print("Lambda(.scups energy)/Lambda(observed energy) at T [K] =",
          " ".join(f"{t:.0e}" for t in T))
    for label, el, io, pop, up in CASES:
        ls, r = coronal_sum(el, io, pop, up, "scups")
        lo, _ = coronal_sum(el, io, pop, up, "observed")
        print(f"{label:20s} dE ratio {r.min():.4f}-{r.max():.4f} ({r.size} tr.)  "
              + " ".join(f"{x:.4f}" for x in ls / lo))


if __name__ == "__main__":
    main()
