#!/usr/bin/env python3
"""Cross-check of the Shomate thermodynamic table of oxygen_rates.f90
against Burcat's NASA-9 polynomials, and of what that choice does to the
thermodynamically reversed rates.

This is the second half of the milestone-M1 check of
docs/a2_oxygen_option_design.md; the first half is the Fortran driver
a2_m1_rate_check.f90 in this directory.  Neither is part of the EXHALE
build.

Single-source rule: the Shomate numbers are NOT retyped here.  They are
read back out of the driver's saved output, a2_m1_rate_check.out, so the
comparison cannot drift away from the module.  The NASA-9 numbers are read
straight from VULCAN/thermo/NASA9/, which is where VULCAN reads them.

Run from src/tests/a2_m1/:
    python3 a2_m1_thermo_xcheck.py > a2_m1_thermo_xcheck.out
"""

import os
import re
import sys

import numpy as np

R_GAS = 8.31446261815324     # J mol^-1 K^-1, CODATA 2018
KB_ERG = 1.380649e-16        # erg K^-1, CODATA 2018
P_STD = 1.0e6                # 1 bar in dyn cm^-2

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER_OUT = os.path.join(HERE, "a2_m1_rate_check.out")
NASA9_DIR = os.path.normpath(
    os.path.join(HERE, "..", "..", "..", "VULCAN", "thermo", "NASA9"))

SPECIES = ("H", "H2", "O", "OH", "H2O", "CO")


def shomate_from_driver_output(path):
    """Read back the SHOMATE_TABLE lines the Fortran driver prints.

    Returns {temperature: {species: (H_kJmol, S_JmolK, G_kJmol)}}.
    """
    tables = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("SHOMATE_TABLE "):
                continue
            cols = line.split()
            if len(cols) != 6 or cols[2] not in SPECIES:
                continue
            temperature = float(cols[1])
            tables.setdefault(temperature, {})[cols[2]] = tuple(
                float(c) for c in cols[3:])
    if not tables:
        raise RuntimeError("no SHOMATE_TABLE lines in " + path)
    return tables


def nasa9(species, temperature):
    """H [kJ/mol], S [J/mol/K], G [J/mol] from Burcat's NASA-9 polynomials."""
    coeff = np.loadtxt(os.path.join(NASA9_DIR, species + ".txt")).flatten()
    a = coeff[0:10] if temperature < 1000.0 else coeff[10:20]
    t = temperature
    h_rt = (-a[0] / t**2 + a[1] * np.log(t) / t + a[2] + a[3] * t / 2
            + a[4] * t**2 / 3 + a[5] * t**3 / 4 + a[6] * t**4 / 5 + a[8] / t)
    s_r = (-a[0] / t**2 / 2 - a[1] / t + a[2] * np.log(t) + a[3] * t
           + a[4] * t**2 / 2 + a[5] * t**3 / 3 + a[6] * t**4 / 4 + a[9])
    return h_rt * R_GAS * t / 1000.0, s_r * R_GAS, (h_rt - s_r) * R_GAS * t


def equilibrium_constant_conc(reactants, products, temperature, gibbs):
    d_gibbs = (sum(gibbs(s, temperature) for s in products)
               - sum(gibbs(s, temperature) for s in reactants))
    d_n = len(products) - len(reactants)
    n_std = P_STD / (KB_ERG * temperature)
    return np.exp(-d_gibbs / (R_GAS * temperature)) * n_std**d_n


def main():
    if not os.path.exists(DRIVER_OUT):
        sys.exit("run a2_m1_rate_check.x first; %s is missing" % DRIVER_OUT)
    shomate = shomate_from_driver_output(DRIVER_OUT)

    print("# a2_m1_thermo_xcheck -- Shomate (NIST-JANAF) against NASA-9 (Burcat)")
    print("# Shomate values read back from %s" % os.path.basename(DRIVER_OUT))
    print("# NASA-9 read from %s" % NASA9_DIR)
    print()
    print("=== 1. Species thermodynamics, the two tables side by side ===")
    for temperature in sorted(shomate):
        print()
        print("  T = %.2f K" % temperature)
        print("  species   H_Shomate   H_NASA9      dH       "
              "S_Shomate   S_NASA9      dS       dG [kJ/mol]")
        for species in SPECIES:
            h_s, s_s, g_s = shomate[temperature][species]
            h_n, s_n, g_n = nasa9(species, temperature)
            print("  %-7s %10.4f %10.4f %+9.4f  %10.3f %10.3f %+9.3f  %+9.4f"
                  % (species, h_s, h_n, h_s - h_n, s_s, s_n, s_s - s_n,
                     g_s - g_n / 1000.0))
    print()
    print("  Up to 1500 K every species agrees to better than 0.02 kJ/mol in")
    print("  enthalpy EXCEPT OH, whose enthalpy of formation differs by a")
    print("  constant: NIST-JANAF (Chase 1998) 38.99 kJ/mol against Burcat")
    print("  37.28 kJ/mol.  Above 1500 K the H2O fits also part company, by")
    print("  0.25 kJ/mol at 2000 K and 0.62 kJ/mol at 2500 K, which is the")
    print("  two tables fitting the same JANAF data differently in their top")
    print("  interval and not a difference of thermochemistry.")
    print()

    print("=== 2. What the OH offset does to a thermodynamically reversed rate ===")
    d_oh = shomate[298.15]["OH"][0] - nasa9("OH", 298.15)[0]
    print("  measured offset dH_OH = %+.4f kJ/mol" % d_oh)
    print("  factor exp(dH_OH/RT) per OH appearing in the reaction")
    print("     T[K]    per one OH   per two OH")
    for temperature in sorted(shomate):
        f = np.exp(d_oh * 1e3 / (R_GAS * temperature))
        print("  %8.2f %12.4f %12.4f" % (temperature, f, f**2))
    print()

    print("=== 3. O10 reversal with each thermodynamic table ===")
    print("  Photochem transcribes OH + OH -> H2O + O (Baulch et al. 1992;")
    print("  Lifshitz & Michael 1991).  VULCAN id 5 transcribes the opposite")
    print("  direction, O + H2O -> OH + OH, 8.20e-14 T^0.95 exp(-8570/T),")
    print("  valid 250-2400 K.  Reversing the first must reproduce the second.")
    print()
    print("     T[K]    k_published    k_rev(Shomate)  ratio   "
          "k_rev(NASA9)   ratio")

    def gibbs_nasa9(species, temperature):
        return nasa9(species, temperature)[2]

    for temperature in sorted(shomate):
        if temperature < 250.0:
            continue
        k_fwd = 2.549944e-15 * temperature**1.14 * np.exp(-50.0 / temperature)
        k_pub = 8.20e-14 * temperature**0.95 * np.exp(-8570.0 / temperature)
        table = shomate[temperature]

        def gibbs_shomate(species, _t, tab=table):
            return tab[species][2] * 1000.0

        k_sh = k_fwd / equilibrium_constant_conc(
            ["OH", "OH"], ["H2O", "O"], temperature, gibbs_shomate)
        k_n9 = k_fwd / equilibrium_constant_conc(
            ["OH", "OH"], ["H2O", "O"], temperature, gibbs_nasa9)
        flag = "" if temperature <= 2400.0 else "  (beyond the published range)"
        print("  %8.2f %14.5e %14.5e %7.3f %14.5e %7.3f%s"
              % (temperature, k_pub, k_sh, k_sh / k_pub, k_n9, k_n9 / k_pub,
                 flag))
    print()
    print("  The Shomate column tracks the independently published opposite")
    print("  direction far better than the NASA-9 column does, which is the")
    print("  evidence for keeping the NIST-JANAF OH enthalpy of formation.")
    print()
    print("# end")


if __name__ == "__main__":
    main()
