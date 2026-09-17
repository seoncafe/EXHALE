#!/bin/bash
export OMP_NUM_THREADS=8
export EXHALE_PTC_DTAU0=1.0
export EXHALE_OUTER_PASSES=40
export SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9.7/k04/output
export FORCE=1
export SEED_ATTEMPTS=1
export RUN_NOTE="L25 step 3: the catalog recipe with \`Numerical flux: ROE\` as the only change. Run directory under models/.L25/, not the catalog."
exec /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/run_case.sh .L25/roe_pfid_HeH9.7
