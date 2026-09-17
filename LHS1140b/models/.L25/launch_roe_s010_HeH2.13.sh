#!/bin/bash
export OMP_NUM_THREADS=8
export EXHALE_PTC_DTAU0=1.0
export EXHALE_OUTER_PASSES=40
export SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/output
export FORCE=1
export SEED_ATTEMPTS=1
export RUN_NOTE="L25 step 3: the catalog recipe with \`Numerical flux: ROE\` as the only change. Run directory under models/.L25/, not the catalog."
exec /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/run_case.sh .L25/roe_s010_HeH2.13
