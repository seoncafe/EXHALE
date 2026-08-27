"""Export CHIANTI-derived collisional line-cooling tables for EXHALE.

Writes, for each metal coolant available in CHIANTI v11, the effective cooling
coefficient Lambda(T) per (n_e * n_ion) in erg cm^3 s^-1 on a shared log-T grid.
These tables are the auditable bridge into Cool_coeff.f90 (Phase 2 port).

Coolants and channels (Huang 2023 Section 2.5):
  Mg I   lambda2853  (3s^2 1S0 -> 3s3p 1P1)              ground -> {5}
  Mg II  h&k lambda2796,2804 (3s 2S -> 3p 2P)            ground -> {2,3}
  Ca II  H&K lambda3934,3968 (4s 2S -> 4p 2P)            ground -> {4,5}
  Fe II  coronal (ground only) and Boltzmann metastable (E < 4.8 eV)

Na I and Fe I are NOT in CHIANTI v11.0.2 (only ionization data); they are built
separately from NIST oscillator strengths + Van Regemorter (see the Fe I / Na I
sub-phase) and appended to these tables there.
"""

import os
import numpy as np
from chianti_cooling import (cooling_lambda, cooling_effective,
                             free_free_cooling, cooling_NaI)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
FE2_ECUT_CM1 = 38459.0  # lowest odd-parity 3d6.4p z6D (Huang's 4.8 eV cutoff)


def build_table(T):
    cols = {}
    cols["MgI_2853"], _ = cooling_lambda("mg", 1, T, lower_levels={1},
                                         upper_levels={5})
    cols["MgII_hk"], _ = cooling_lambda("mg", 2, T, lower_levels={1},
                                        upper_levels={2, 3})
    cols["CaII_HK"], _ = cooling_lambda("ca", 2, T, lower_levels={1},
                                        upper_levels={4, 5})
    cols["NaI_D"] = cooling_NaI(T)
    cols["FeII_coronal"] = cooling_effective("fe", 2, T, pop="coronal")
    cols["FeII_boltz"] = cooling_effective("fe", 2, T, pop="boltzmann",
                                           e_cut_cm1=FE2_ECUT_CM1)
    cols["freefree_Z1"] = free_free_cooling(T, Z=1.0)
    return cols


def main():
    T = np.logspace(3.0, 5.0, 201)  # 1e3 - 1e5 K
    cols = build_table(T)
    names = list(cols.keys())
    path = os.path.join(OUTDIR, "metal_cooling_chianti.txt")
    with open(path, "w") as fh:
        fh.write("# CHIANTI v11.0.2 collisional line-cooling coefficients for EXHALE\n")
        fh.write("# Lambda(T) per (n_e * n_ion)  [erg cm^3 s^-1], optically thin\n")
        fh.write("# Source: chianti_cooling.py (Burgess-Tully descaled .scups)\n")
        fh.write("# FeII_coronal = ground-level only; FeII_boltz = Boltzmann over\n")
        fh.write("#   metastable levels E < 4.8 eV (38459 cm^-1), upper depleted\n")
        fh.write("# freefree_Z1 = Huang Eq.9, 1.9095e-25 * Z^2 * (T/1e4)^0.55, Z=1\n")
        fh.write("# NaI_D = 3s-3p D doublet, rigorous Van Regemorter (gbar=0.2,\n")
        fh.write("#   calibrated to CHIANTI Mg II/Ca II); Na I not in CHIANTI v11\n")
        fh.write("# columns:\n")
        fh.write("{:>12}".format("T_K"))
        for nm in names:
            fh.write(" {:>14}".format(nm))
        fh.write("\n")
        for i in range(T.size):
            fh.write("{:12.5e}".format(T[i]))
            for nm in names:
                fh.write(" {:14.6e}".format(cols[nm][i]))
            fh.write("\n")
    print(f"wrote {path}  ({T.size} rows, {len(names)} coolant columns)")
    # quick sanity echo
    for Tc in (5000., 1e4, 2e4):
        i = int(np.argmin(abs(T - Tc)))
        print(f"  T={T[i]:9.1f}  " +
              "  ".join(f"{nm}={cols[nm][i]:.3e}" for nm in names))


if __name__ == "__main__":
    main()
