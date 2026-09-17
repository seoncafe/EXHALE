#!/bin/bash
# Run one regression case on a scratch copy with a named binary, then compare
# the four output files against the reference snapshot.
#   run_case.sh <case> <binary> <tag>
set -e
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
R=$EX/backup/regression
S=$EX/scratchpad/L6c
c="$1"; exe="$2"; tag="$3"
W=$S/runs/$tag/$c
rm -rf $W; mkdir -p $W/output
for f in "$R/$c"/*; do
  b=$(basename "$f")
  case "$b" in output|check.log|run.log|EXHALE.x) ;; IC) \cp -f "$f"/*_IC.txt $W/output/ ;; *) \cp -f "$f" $W/ ;; esac
done
cap=""
[ -f "$R/$c/maxsteps" ] && cap="EXHALE_MAXSTEPS=$(cat $R/$c/maxsteps)"
\cp -f "$exe" $W/EXHALE.x
( cd $W && env $cap OMP_NUM_THREADS=1 ./EXHALE.x > run.log 2>&1 ) || true
{
  echo "[$c/$tag] $(grep -E 'final:' $W/run.log | tail -n1)"
  for f in Hydro_ioniz.txt Ion_species.txt Hydro_ioniz_adv.txt Ion_species_adv.txt; do
    g=$R/golden/$c/$f
    if [ ! -f "$g" ]; then echo "     MISS $f"; continue; fi
    if cmp -s <(grep -v '^ *#' $W/output/$f) <(grep -v '^ *#' $g); then
      echo "     PASS $f (data identical)"
    else
      msg=$("$R/compare_within_tolerance.py" "$W/output/$f" "$g" 1e-3 2>&1) && \
        echo "     PASS $f ($msg)" || echo "     FAIL $f ($msg)"
    fi
  done
} > $S/runs/$tag/$c.report 2>&1
