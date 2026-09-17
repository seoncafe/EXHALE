#!/bin/bash
# One XUV continuation chain for the L14 ladder: each rung is solved directly
# at the raised pseudo-time start, from the state the previous rung wrote.
# The dtau0 = 1 stage of run_case.sh is skipped: on the He/H = 2.13 column of
# LHS 1140 b it falls into the ramp collapse of plan item L4h (lam 1e-6,
# dtau 6e-4) at every XUV step tried, while the mapped seed at dtau0 = 1e8 is
# in the Newton basin (MEASURED 2026-09-14, L14).
#
#   ladder.sh <HeH> <first-seed-dir> <rung> [<rung> ...]
set -u
M=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
EB=$EX/EXHALE_L14.x
HEH=$1; SEED=$2; shift 2
for RUNG in "$@"; do
   if [ "$RUNG" = x001 ]; then
      # The catalog case: solved by the runner, so it gets its
      # advection-corrected profiles, its transit spectrum and its
      # REPRODUCE.md, with the raised pseudo-time start this ladder needs.
      CASE=atomic_scalar_gj1132x0.01_kzz1e9/HeH$HEH
      echo "[$RUNG/$HEH] seed $SEED, through run_case.sh"
      ( cd "$M" && OMP_NUM_THREADS=8 EXHALE_BIN=$EB SEED=$SEED \
          EXHALE_PTC_DTAU0=1.0e8 ./run_case.sh "$CASE" )
      echo "[$RUNG/$HEH] $(grep -hE 'outer pass [0-9]+: ACCEPTED|returned info' $M/$CASE/run.log | tail -n 1)"
      continue
   fi
   D=$M/.L14/${RUNG}_HeH$HEH
   mkdir -p "$D/output"
   echo "[$RUNG/$HEH] seed $SEED"
   python3 $EX/src/utils/map_state_to_grid.py "$SEED" \
       $M/current_grid_Hydro_ioniz.txt "$D/output" --ic > "$D/seed.log" 2>&1 \
       || { echo "[$RUNG/$HEH] seed mapping FAILED"; exit 1; }
   ( cd "$D" && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
       EXHALE_OUTER_PASSES=40 "$EB" > run.log 2>&1 )
   if grep -q 'certified=T' "$D/output/Hydro_ioniz.txt" 2>/dev/null; then
      echo "[$RUNG/$HEH] CERTIFIED: $(grep -hE 'outer pass [0-9]+: ACCEPTED' $D/run.log | tail -n 1)"
      SEED=$D/output
   else
      echo "[$RUNG/$HEH] NOT CERTIFIED: $(grep -hE 'returned info|spent its budget|REFUSED --' $D/run.log | tail -n 1)"
      exit 2
   fi
done
echo "[$HEH] ladder finished"
