#!/bin/bash
# The XUV continuation chain of the lower-atmosphere-profile case
# atomic_photochem_gj1132x*_kzzprofile/HeH9.7, from the certified 0.10 case
# down to 0.01.  Each rung goes through run_case.sh, which builds the target
# grid header from the profile's own matching level (L13) and writes the
# advection-corrected profiles, the transit spectrum and REPRODUCE.md.
#
# The pseudo-time start is raised to 1e8 for every rung: MEASURED 2026-09-14
# (L14), the dtau0 = 1 solve of these cold, weakly irradiated columns falls
# into the ramp collapse of plan item L4h, while the mapped seed at 1e8 is in
# the Newton basin.  The 20-fold XUV step this case was first given
# (0.20 -> 0.01) leaves no seed in the basin at all, which is what the ladder
# replaces.
set -u
M=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models
EB=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L14.x
SEED=${SEED:-$M/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/output}
for RUNG in ${RUNGS:-p005 p003 p002 p0015 p001}; do
   if [ "$RUNG" = p001 ]; then
      CASE=atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7
   else
      CASE=.L14/${RUNG}_HeH9.7
   fi
   echo "[$RUNG] seed $SEED"
   ( cd "$M" && OMP_NUM_THREADS=8 EXHALE_BIN=$EB SEED=$SEED \
       EXHALE_PTC_DTAU0=1.0e8 FORCE=1 ./run_case.sh "$CASE" )
   # THE CERTIFICATION IS ON THE _IC PAIR, NOT ON Hydro_ioniz.txt. run_case.sh
   # copies the solved state to output/*_IC.txt and then takes the
   # advection-corrected pass, which rewrites output/Hydro_ioniz.txt with a
   # header that makes no stationary claim ("certified=F
   # cert_reason=no_stationary_claim"); the state itself moves by 4e-14,
   # since that pass runs at CFL 1e-12, but its header is not the verdict.
   if grep -q 'certified=T' "$M/$CASE/output/Hydro_ioniz_IC.txt" 2>/dev/null; then
      echo "[$RUNG] CERTIFIED: $(grep -hE 'outer pass [0-9]+: ACCEPTED' $M/$CASE/run.log | tail -n 1)"
      SEED=$M/$CASE/output
   else
      echo "[$RUNG] NOT CERTIFIED: $(grep -hE 'returned info|spent its budget|REFUSED --' $M/$CASE/run.log | tail -n 1)"
      exit 2
   fi
done
echo "photochem ladder finished"
