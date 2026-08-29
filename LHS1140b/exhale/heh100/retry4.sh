#!/bin/bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
cd "$(dirname "$0")"
export OMP_NUM_THREADS=8
\cp -f results_constres_fiducialA/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f results_constres_fiducialA/Ion_species.txt output/Ion_species_IC.txt
for i in 1 2 3 4; do
  EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$EX/EXHALE.x" > run.log 2>&1
  R=$(grep -oE 'done info=[0-9]+ +\|\|R\|\|= *[0-9.E+-]+' run.log | tail -n 1)
  echo "JFNK pass $i: $R"
  grep -q "done info=0" run.log && { echo CONVERGED; break; }
  \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
  \cp -f output/Ion_species.txt output/Ion_species_IC.txt
done
grep -q "done info=0" run.log || { echo CHAIN_FAILED; exit 1; }
\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$EX/EXHALE.x" > pp.log 2>&1
sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
echo "PP $(grep 'steady-state Mdot' pp.log | tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
grep -i 'He 10830 metrics' transit.log | tail -n 1
echo CHAIN_OK
