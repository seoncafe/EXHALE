"""Arithmetic counterexamples to the proposed energy equation, not a solver test."""

import math

chi_h = 13.598434599
photon_h = 20.0
heat_h = photon_h - chi_h
print("H photo: expected heat, rev1 heat, ledger loss [eV]:",
      heat_h, heat_h - chi_h, chi_h)

# Chemical product energy measured by audit_probe.f90 at the reviewed source.
chemical_h2_double = 31.67494411407
photon_h2 = 80.0
print("H2 double: expected heat, rev1 heat [eV]:",
      photon_h2 - chemical_h2_double,
      photon_h2 - 2.0 * chemical_h2_double)

# Miller Table 6 entries at 1000 K; also present in h3p_cooling.f90.
print("Table6 exponent:", math.log(0.0108 / 0.0013) / math.log(100.0))
