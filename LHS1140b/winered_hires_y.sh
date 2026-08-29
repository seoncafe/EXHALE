# The spectrograph of the LHS 1140 b He I 10830 observation.
#
# Cherubim et al. (2026) observed this transit with WINERED in HIRES-Y mode.
# Their Supplement fixes the kernel: "we assumed a Gaussian line spread
# function with FWHM of 10,833 A/R, where R = 68,000 is the resolving power
# of the WINERED spectrograph in HIRES-Y mode", and their p-winds retrieval
# uses "FWHM = 4.4 km/s (equivalent to R = 68,000)".
#
# EXHALE_transit.py carries R = 8e4 (CARMENES) as its built-in He I 10830
# default, which is right for the other planets in this repository and wrong
# for this one. Every LHS 1140 b transit synthesis therefore sources this
# file before calling EXHALE_transit.py, so that the modelled line is
# convolved with the kernel the measurement was made through.
#
# Do not change the built-in default in EXHALE_transit.py: it is shared with
# HD 189733 b, HD 209458 b, WASP-52 b and WASP-121 b, whose stored results
# were produced at R = 8e4.
export EXHALE_TRANSIT_RES_HETR=68000
