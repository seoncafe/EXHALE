"""Arithmetic checks of revision 2 specifications; no atmosphere is evolved."""

chi_h = 13.598434599
photon_h = 20.0
print("Rev2 H photo thermal increment [eV]:", photon_h - chi_h)
print("Rev2 H photo material sum [eV]:", (photon_h - chi_h) + chi_h)

# Published Miller Tables 4 and 6, at 5000 K. The final Table 6 density is finite.
unscaled_lte = 3.4119e-18
scaled_lte = 4.7688e-18
table6_emission = 4.7587e-18
print("Table6 emission / Table4 unscaled LTE:", table6_emission / unscaled_lte)
print("Table6 emission / Table4 scaled LTE:", table6_emission / scaled_lte)

# A photon emitted by cell A is fully absorbed by cell B, with none escaping.
photon_exchange = 10.0
print("Two-cell radiation exchange, correct material sum [eV]:",
      -photon_exchange + photon_exchange)
print("If only domain-escaping emission is subtracted [eV]:", photon_exchange)
