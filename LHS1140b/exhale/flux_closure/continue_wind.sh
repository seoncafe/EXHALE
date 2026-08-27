#!/bin/bash
# Continue one profile-driven wind solve from its own output until the JFNK
# finish reports info=0, and leave the converged state where a closure run
# can take it as its --seed.  The stall the closure driver stops on is the
# base contact mode: restarting the solve on its own iterate lets the outer
# diffusion passes carry the composition the rest of the way.
#
# usage: ./continue_wind.sh <dir_with_input.inp_and_output> [n_passes]
set -u
d=$1; n=${2:-8}
BIN=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x
cd "$d" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
for i in $(seq 1 "$n"); do
  \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
  \cp -f output/Ion_species.txt output/Ion_species_IC.txt
  EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$BIN" > "cont$i.log" 2>&1
  info=$(grep -o 'done info=[0-9]*' "cont$i.log" | tail -n 1)
  echo "[$d] pass $i: $info $(grep -o '||R||= *[0-9.E+-]*' cont$i.log | tail -n 1)" \
       "$(grep -o 'composition drift = *[0-9.E+-]*' cont$i.log | tail -n 1)"
  [ "$info" = "done info=0" ] && break
done
