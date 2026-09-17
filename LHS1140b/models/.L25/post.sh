#!/bin/bash
# The advection-corrected profiles and the transit spectrum of one solved
# state, by the recipe the catalog products were made with (the tree binary
# predates the "Restart intent: stationary evaluate" route that writes the
# _adv files, so that route measures the state and writes no profiles here).
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
D=$1
cd "$D" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
for f in Hydro_ioniz Ion_species; do
   [ -f output/$f.txt ] || { echo "$D: no output/$f.txt"; exit 1; }
   \cp -f output/$f.txt output/${f}_IC.txt
done
\cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
       -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
$EX/EXHALE.x > pp.log 2>&1
pprc=$?
\mv -f input.inp.solved input.inp
[ $pprc -ne 0 ] && [ $pprc -ne 2 ] && { echo "$D: pp exited $pprc"; exit 1; }
. $EX/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=$EX python3 $EX/EXHALE_transit.py > transit.log 2>&1
trc=$?
[ -f tpm_He10830_metrics.txt ] || { echo "$D: transit exited $trc"; exit 1; }
echo "$D: products written"
