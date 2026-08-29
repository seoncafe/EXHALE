#!/bin/bash
# Every run of this directory, in the order it was made.  The binary is the
# one of 2026-08-29 04:56 and every transit synthesis goes through
# ../winered_hires_y.sh, so the He I 10830 line is convolved at WINERED
# HIRES-Y, R = 68,000.
#
# Nothing outside this directory is written: ../heh0p55, ../heh0p55_gj699,
# ../heh2p13_diff_kzz1e9, ../heh2_diff_kzz1e9, ../heh4p5_diff_kzz0,
# ../flux_closure/heh9p7, ../crossings_gm25 and ../ladder_gm25 are read as
# references and seeds only.
#
# Each crossing is solved twice: from the nearest stored arm, and from the
# same common parent the arms of record were seeded from (the _P arms).
set -u
cd "$(dirname "$0")" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# ---- 1. the four scalar-base crossings -----------------------------------
./run_wellmixed.sh wm1132_0p4245   0.4245 ../heh0p55       ../crossings_gm25/rs_heh0p425
./run_wellmixed.sh wm1132_0p4245_P 0.4245 ../heh0p55       ../crossings_gm25/wm_heh0p45
./run_wellmixed.sh wm699_0p0479    0.0479 ../heh0p55_gj699 ../ladder_gm25/g2_heh0p048
./run_wellmixed.sh wm699_0p0479_P  0.0479 ../heh0p55_gj699 ../ladder_gm25/g1_heh0p050
./run_kzz.sh       kz1e9_1p6261    1.6261 1.0e9 ../crossings_gm25/kz_heh1p65
./run_kzz.sh       kz1e9_1p6261_P  1.6261 1.0e9 ../heh2_diff_kzz1e9
./run_kzz.sh       kz0_3p8103      3.8103 0.0   ../ladder_gm25/kz0_heh3p79
./run_kzz.sh       kz0_3p8103_P    3.8103 0.0   ../heh4p5_diff_kzz0

# ---- 2. the elemental-flux closure at its crossing reservoir -------------
# The single parent of the whole ladder of record, ../flux_closure/heh9p7/k03
# (reservoir 9.71), and its converged fluxes: run_closure.sh defaults to both.
# The interpreter is the repository's env/photochem (photochem 0.9.0); see
# results.txt section 5 for why the stored 0.8.4 path cannot be used.
./run_closure.sh   fc_heh9p048     9.048
OMP_NUM_THREADS=4 ./resolve_wind.sh fc_heh9p048/k01 rw_9p048 3

# ---- 3. the table --------------------------------------------------------
./measure.py G699=wm699_0p0479 G699p=wm699_0p0479_P \
             W1132=wm1132_0p4245 W1132p=wm1132_0p4245_P \
             K1e9=kz1e9_1p6261 K1e9p=kz1e9_1p6261_P \
             K0=kz0_3p8103 K0p=kz0_3p8103_P \
             C9p048=rw_9p048:9.048
