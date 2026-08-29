#!/bin/bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
cd "$(dirname "$0")"
export OMP_NUM_THREADS=6
for c in solar heh1 heh10 heh100 heh1000 heh10000; do
  ( cd $c
    \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
    \cp -f output/Ion_species.txt output/Ion_species_IC.txt
    sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
    "$EX/EXHALE.x" > pp.log 2>&1
    sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
    MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1 )
  echo "$c PP+transit done: $(grep 'steady-state Mdot' $c/pp.log|tail -1)  $(grep 'He 10830 metrics' $c/transit.log|tail -1|sed 's/.*metrics: //')"
done
echo ALL_DONE
