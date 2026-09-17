#!/bin/bash
export OMP_NUM_THREADS=8
export EXHALE_PTC_DTAU0=1.0e8
export EXHALE_OUTER_PASSES=40
export NOSEED=1
export FORCE=1
export SEED_ATTEMPTS=1
export RUN_NOTE="L25 step 3 item 3, grid convergence: the seed was mapped onto this grid beforehand and NOSEED=1 was set, because run_case.sh maps onto the 500-cell grid of models/current_grid_Hydro_ioniz.txt only."
exec /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/run_case.sh .L25/g1000_hllc
