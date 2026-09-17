#!/bin/bash
export OMP_NUM_THREADS=8
export EXHALE_PTC_DTAU0=1.0
export EXHALE_OUTER_PASSES=40
export SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.stopped/atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7_dtau0-1_20260916055251/output
export FORCE=1
export SEED_ATTEMPTS=1
export RUN_NOTE="L25 step 3: the same case at the default pseudo-time start, seeded from the stopped catalog attempt's written state, which is the configuration step 2 measured the Roe convergence in."
exec /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/run_case.sh .L25/roe_s010_HeH9.7_d1
