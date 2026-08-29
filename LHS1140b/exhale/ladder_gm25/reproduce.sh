#!/bin/bash
# Every run of this directory, in the order it was made.  The binary is the
# one of 2026-08-29 04:56 (Update_EXHALE sections 84, 86, 87 and 88 all in),
# and every transit synthesis goes through ../winered_hires_y.sh, so the
# He I 10830 line is convolved at WINERED HIRES-Y, R = 68,000.
#
# Nothing outside this directory is written.  ../heh2p13_diff_kzz1e9 is the
# K_zz control input, ../heh0p55_gj699 the GJ 699 one, and the seeds --
# ../heh2p7_diff_kzz1e8, ../heh1p4_diff_kzz1e10, ../heh3p4_diff_kzz1e7,
# ../heh1_diff_kzz1e11, ../heh4p3_diff_kzz1e6, ../heh4p5_diff_kzz1e5,
# ../heh4p5_diff_kzz0, ../heh0p04_gj699, ../flux_closure/heh9p7/k03 and
# ../crossings_gm25/fcA_heh8p45/k02 -- are all read only.
set -u
cd "$(dirname "$0")" || exit 1
t(){ echo "$1" | tr . p; }

# ---- 1. the reservoir band the adapter used to refuse, 8.5 - 8.7 ---------
# It runs now because src/utils/lower_profile_schema.py --abundance-tol was
# raised from 1.0e-4 to 1.0e-3.  Same single parent as the ladder of record,
# ../flux_closure/heh9p7/k03 (reservoir 9.71, so 9.71/8.5 = 1.14x), and the
# same starting elemental fluxes, which run_closure.sh defaults to.
for h in 8.5 8.55 8.6 8.65 8.7; do ./run_closure.sh fcA_heh$(t $h) $h; done
for s in "fcA_heh8p5/k02 rw_8p5"   "fcA_heh8p55/k02 rw_8p55" \
         "fcA_heh8p6/k02 rw_8p6"   "fcA_heh8p65/k02 rw_8p65" \
         "fcA_heh8p7/k01 rw_8p7"; do
  set -- $s; ./resolve_wind.sh $1 $2 3
done
# two arms of record re-run here, to check that this directory reproduces them
for h in 8.45 8.75; do ./run_closure.sh rep_heh$(t $h) $h; done
./resolve_wind.sh rep_heh8p45/k02 rwrep_8p45 3
./resolve_wind.sh rep_heh8p75/k01 rwrep_8p75 3
# and the same three arms re-closed from the neighbour below, the seed test
HERE=$PWD
./run_closure.sh fcB_heh8p5  8.5  ../crossings_gm25/fcA_heh8p45/k02/output \
                 5.3790386146E+06 1.9060985299E+07
./run_closure.sh fcB_heh8p65 8.65 $HERE/fcA_heh8p5/k02/output \
                 5.3604612723E+06 1.9355450861E+07
./run_closure.sh fcB_heh8p7  8.7  $HERE/fcA_heh8p6/k02/output
for s in "fcB_heh8p5/k00 rwB_8p5" "fcB_heh8p65/k00 rwB_8p65" \
         "fcB_heh8p7/k01 rwB_8p7"; do set -- $s; ./resolve_wind.sh $1 $2 3; done

# ---- 2. the K_zz ladder, scalar base -------------------------------------
# Five arms per decade, spanning +/-8% of the crossing predicted by scaling
# the stored ladder of ../../kzz_decision.md section 6.1 by the ratio the
# K_zz = 1e9 crossing moved (1.6261/2.09 = 0.778).  Every seed is the nearest
# stored arm of the same K_zz, all of them within 1.35x in reservoir.
./run_kzz.sh kz1e8_heh2p05  2.05 1.0e8  ../heh2p7_diff_kzz1e8
./run_kzz.sh kz1e8_heh2p15  2.15 1.0e8  ../heh2p7_diff_kzz1e8
./run_kzz.sh kz1e8_heh2p23  2.23 1.0e8  ../heh2p7_diff_kzz1e8
./run_kzz.sh kz1e8_heh2p32  2.32 1.0e8  ../heh2p7_diff_kzz1e8
./run_kzz.sh kz1e8_heh2p41  2.41 1.0e8  ../heh2p7_diff_kzz1e8
for h in 1.06 1.10 1.15 1.20 1.25 1.29; do
  ./run_kzz.sh kz1e10_heh$(t $h) $h 1.0e10 ../heh1p4_diff_kzz1e10; done
for h in 2.70 2.82 2.94 3.06 3.18; do
  ./run_kzz.sh kz1e7_heh$(t $h)  $h 1.0e7  ../heh3p4_diff_kzz1e7;  done
for h in 0.795 0.83 0.865 0.90 0.93; do
  ./run_kzz.sh kz1e11_heh$(t $h) $h 1.0e11 ../heh1_diff_kzz1e11;   done
for h in 3.19 3.32 3.46 3.60 3.74; do
  ./run_kzz.sh kz1e6_heh$(t $h)  $h 1.0e6  ../heh4p3_diff_kzz1e6;  done
for h in 3.35 3.50 3.64 3.79 3.93; do
  ./run_kzz.sh kz1e5_heh$(t $h)  $h 1.0e5  ../heh4p5_diff_kzz1e5;  done
for h in 3.35 3.50 3.64 3.79 3.93; do
  ./run_kzz.sh kz0_heh$(t $h)    $h 0.0    ../heh4p5_diff_kzz0;    done

# where the eddy term stops acting: one composition, He/H = 3.79 (the
# crossing of the plateau), solved at each decade from one common seed so
# that the whole probe is a single lineage
for k in 0.0 1.0e1 1.0e2 1.0e3 1.0e4 1.0e5 1.0e6; do
  s=$(echo $k | sed -e 's/1\.0e/1e/' -e 's/^0\.0$/0/')
  ./run_kzz.sh probe${s}_heh3p79 3.79 $k kz0_heh3p79
done

# ---- 3. GJ 699, well mixed ----------------------------------------------
# lineage G1, seeded from the stored ../heh0p04_gj699
for h in 0.042 0.044 0.046 0.048 0.050; do
  ./run_wellmixed.sh g1_heh$(t $h) $h ../heh0p55_gj699 ../heh0p04_gj699; done
# lineage G2, seeded from g1_heh0p050 -- a solution of the current binary,
# which is what removes the non-monotonic point of G1 (as W1 -> W2 did in
# ../crossings_gm25)
for h in 0.042 0.044 0.046 0.048; do
  ./run_wellmixed.sh g2_heh$(t $h) $h ../heh0p55_gj699 g1_heh0p050; done

# ---- the tables and the crossings ---------------------------------------
./measure.py kz1e8_heh2p05=kz1e8_heh2p05 ...        # see results.txt
./crossing.py "K_zz=1e8" 2.05=1.01072 ...           # see results.txt
