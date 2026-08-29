#!/bin/bash
# Restart one eddy-diffusion case from its own converged output and take one
# more wind + composition outer pass.  The composition relaxation of a case
# whose helium profile moves far from the Kzz = 0 seed does not settle inside
# the outer-loop ceiling of a single invocation, so the loop is closed by
# repeating the invocation until the reported composition drift stops moving.
# usage: ./continue_kzz_case.sh <tag> [n_repeats]
set -u
tag=$1; nrep=${2:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
cd "$(dirname "$0")/$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}

for i in $(seq 1 "$nrep"); do
   \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
   \cp -f output/Ion_species.txt output/Ion_species_IC.txt
   sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
   EXHALE_DIFFUSION_CHECK=1 EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 \
      "$BIN" > "run_cont$i.log" 2> diffcheck.log
   echo "[$tag] cont $i: $(grep -o 'done info=[0-9]*' run_cont$i.log | tail -n 1)" \
        "$(grep -o '||R||= *[0-9.E+-]*' run_cont$i.log | tail -n 1)" \
        "$(grep -o 'composition drift = *[0-9.E+-]*' run_cont$i.log | tail -n 1)"
done

\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$?  $(grep 'steady-state Mdot' pp.log | tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"
