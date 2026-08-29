#!/bin/bash
# Re-solve one closure iterate's wind repeatedly on its own output, with the
# lower-atmosphere profile held fixed, and re-synthesize the line.  This
# separates the two things a closure iteration does: it tests whether the
# outer region an arm inherited from its seed (Update_EXHALE section 83) is
# shed by re-solving the wind alone, with the composition unchanged.
#
# usage: ./resolve_wind.sh <src_k_dir> <tag> [n_passes]
set -u
src=$1; tag=$2; n=${3:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"
BIN=$EX/EXHALE.x
cd "$(dirname "$0")" || exit 1
rm -rf "$tag"; mkdir -p "$tag/output"
cp "$src/input.inp" "$src/lower_atmosphere_profile.dat" "$tag/"
cp "$src/output/Hydro_ioniz.txt" "$src/output/Ion_species.txt" "$tag/output/"
cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
for i in $(seq 1 "$n"); do
  \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
  \cp -f output/Ion_species.txt output/Ion_species_IC.txt
  EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$BIN" > "cont$i.log" 2>&1
  echo "[$tag] pass $i: $(grep -o 'done info=[0-9]*' cont$i.log | tail -n 1)" \
       "$(grep -o '||R||= *[0-9.E+-]*' cont$i.log | tail -n 1)"
done
\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] done $(grep 'steady-state Mdot' pp.log | tail -n 1)"
