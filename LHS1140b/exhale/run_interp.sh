#!/bin/bash
# Fill in the He/H interval between solar and 1, where the observed red
# equivalent width falls. Each case is seeded from heh1's converged solution;
# load_IC rescales the loaded H/He onto this case's composition.
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd "$(dirname "$0")"
export OMP_NUM_THREADS=6
for c in heh0p55 heh0p6; do
  ( cd $c
    sed -i 's/^Load IC? False$/Load IC? True/' input.inp
    \cp -f ../heh1/output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
    \cp -f ../heh1/output/Ion_species.txt output/Ion_species_IC.txt
    EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$EX/EXHALE.x" > run.log 2>&1
    if grep -q "done info=0" run.log; then
      \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
      \cp -f output/Ion_species.txt output/Ion_species_IC.txt
      sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
      "$EX/EXHALE.x" > pp.log 2>&1
      sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
      MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
      mkdir -p tpm_turb
      EXHALE_TRANSIT_TURB=1 EXHALE_TRANSIT_SAVE_PREFIX=tpm_turb/ MPLBACKEND=Agg \
        PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > tpm_turb/transit.log 2>&1
    fi )
  echo "$c: $(grep -oE 'done info=[0-9]+' $c/run.log|tail -1)  $(grep 'steady-state Mdot' $c/pp.log 2>/dev/null|tail -1|sed 's/.*= //')  $(grep 'He 10830 metrics' $c/transit.log 2>/dev/null|tail -1|sed 's/.*metrics: //;s/ ->.*//')"
done
echo ALL_DONE
