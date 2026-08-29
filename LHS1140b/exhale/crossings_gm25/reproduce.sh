#!/bin/bash
# Every run of this directory, in the order it was made.  The binary is the
# one of 2026-08-29 04:56 (Update_EXHALE sections 84, 86, 87 and 88 all in),
# and every transit synthesis goes through ../winered_hires_y.sh, so the
# He I 10830 line is convolved at WINERED HIRES-Y, R = 68,000.
#
# Nothing outside this directory is written: ../heh0p55,
# ../heh2_diff_kzz1e9, ../heh2p13_diff_kzz1e9, ../heh1_diff_kzz1e9 and
# ../flux_closure/heh9p7 are read as seeds only.
set -u
cd "$(dirname "$0")" || exit 1

# ---- 1. well mixed: no diffusion, scalar base, no metals, GJ 1132 SED -----
# lineage W1, seeded from ../heh0p55 (the stored He/H = 0.55 solution)
for h in 0.40 0.41 0.42 0.425 0.43 0.435 0.44 0.45 0.46; do
  ./run_wellmixed.sh wm_heh$(echo $h | tr . p) $h
done
# lineage W2, seeded from wm_heh0p45 -- a solution of the current binary,
# which is what removes the one non-monotonic point of W1
for h in 0.40 0.41 0.42 0.425 0.43 0.435 0.44 0.46; do
  ./run_wellmixed.sh rs_heh$(echo $h | tr . p) $h wm_heh0p45
done

# ---- 2. K_zz = 1e9, scalar base ------------------------------------------
for h in 1.55 1.60 1.65 1.70 1.75; do
  ./run_kzz1e9.sh kz_heh$(echo $h | tr . p) $h            # seed ../heh2_diff_kzz1e9
done
# the same two bracket arms from the two other available seeds, to measure
# how much of the crossing is the seed
for h in 1.60 1.65; do
  ./run_kzz1e9.sh kzB_heh$(echo $h | tr . p) $h ../heh2p13_diff_kzz1e9
  ./run_kzz1e9.sh kzC_heh$(echo $h | tr . p) $h ../heh1_diff_kzz1e9
done

# ---- 3. elemental-flux closure -------------------------------------------
# One parent for the whole ladder: ../flux_closure/heh9p7/k03, whose
# reservoir is 9.71, so no arm is more than 9.71/8.0 = 1.21x away from its
# seed.  run_closure.sh defaults to that parent and to its converged fluxes.
for h in 8.0 8.3 8.45 8.75 8.9 9.0 11.0 11.5; do
  ./run_closure.sh fcA_heh$(echo $h | tr . p) $h
done
for h in 8.8 9.1 9.3 9.5 9.7 10.0; do
  ./run_closure.sh fc_heh$(echo $h | tr . p) $h
done
HERE=$PWD
./run_closure.sh fc_heh10p3 10.3 $HERE/fc_heh10p0/k00/output 4.616192E+06 2.011657E+07
./run_closure.sh fc_heh10p6 10.6 $HERE/fc_heh10p0/k00/output 4.616192E+06 2.011657E+07
# reservoir He/H in [8.5, 8.7] cannot be run: the adapter refuses the
# photochemical column because the deepest-level N/H departs from the input
# by 1.0e-4 to 1.9e-4 against its 1.0e-4 tolerance.  The refusal is the same
# for every trial flux tried, so it is a property of the reservoir value.

# the wind of every closure iterate, re-solved on its own output with the
# profile held fixed, so that arms that stopped at different k are compared
# at the same place -- their common wind fixed point
for s in "fcA_heh8p0/k02 rw_8p0"   "fcA_heh8p3/k02 rw_8p3" \
         "fcA_heh8p45/k02 rw_8p45" "fcA_heh8p75/k01 rw_8p75" \
         "fc_heh8p8/k01 rw_8p8"    "fcA_heh8p9/k01 rw_8p9" \
         "fcA_heh9p0/k01 rw_9p0"   "fc_heh9p1/k01 rw_9p1" \
         "fc_heh9p3/k00 rw_9p3"    "fc_heh9p5/k00 rw_9p5" \
         "fc_heh9p7/k00 rw_9p7"    "fc_heh10p0/k00 rw_10p0" \
         "fc_heh10p3/k00 rw_10p3"  "fc_heh10p6/k01 rw_10p6" \
         "fcA_heh11p0/k02 rw_11p0" "fcA_heh11p5/k02 rw_11p5"; do
  set -- $s; ./resolve_wind.sh $1 $2 3
done
# three arms re-closed from a different seed, to measure the seed systematic
./run_closure.sh fcB_heh8p75 8.75 $HERE/fcA_heh8p45/k02/output 5.3790386146E+06 1.9060985299E+07
./run_closure.sh fcB_heh9p0  9.0  $HERE/fc_heh10p0/k00/output  4.6161923967E+06 2.0116567097E+07
./run_closure.sh fcB_heh9p1  9.1  $HERE/fcA_heh8p9/k01/output  5.1204540716E+06 1.9639477396E+07
for s in "fcB_heh8p75/k00 rwB_8p75" "fcB_heh9p0/k01 rwB_9p0" "fcB_heh9p1/k00 rwB_9p1"; do
  set -- $s; ./resolve_wind.sh $1 $2 3
done

# ---- the table and the crossings -----------------------------------------
./measure.py wm0p42=rs_heh0p42 wm0p425=rs_heh0p425      # ... see results.txt
./crossing.py "well mixed" 0.42=1.10083 0.425=1.10880   # ... see results.txt
