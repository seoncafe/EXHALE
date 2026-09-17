#!/bin/bash
# One equilibrium solve + post-process at the case's GOLDEN state, with a
# named binary. The state is held fixed, so the cooling and heating written
# are those of the two builds evaluated on the same gas.
set -e
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
R=$EX/backup/regression; S=$EX/scratchpad/L6c
c="$1"; exe="$2"; tag="$3"
W=$S/pp/$tag/$c; rm -rf $W; mkdir -p $W/output
for f in "$R/$c"/*; do
  b=$(basename "$f")
  case "$b" in output|check.log|run.log|EXHALE.x|IC|EXHALE*.out|*.log|maxsteps) ;; *) \cp -f "$f" $W/ ;; esac
done
\cp -f $S/runs/C/$c/output/Hydro_ioniz.txt $W/output/Hydro_ioniz_IC.txt
\cp -f $S/runs/C/$c/output/Ion_species.txt $W/output/Ion_species_IC.txt
sed -i -e 's/^Do only PP: .*/Do only PP: True/' -e 's/^Load IC?: .*/Load IC?: True/' \
       -e 's/^Load IC? .*/Load IC? True/' -e '/^Restart intent:/d' -e '/^Solver:/d' $W/input.inp
\cp -f "$exe" $W/EXHALE.x
( cd $W && OMP_NUM_THREADS=1 ./EXHALE.x > run.log 2>&1 ) || echo "RUN FAILED $c/$tag (see $W/run.log)"
