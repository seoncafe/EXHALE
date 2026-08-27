#!/bin/bash
# One point of the XUV grid on the flux-closed crossing solution: the
# converged lower-atmosphere profile is held fixed and only the stellar
# spectrum the wind sees is scaled, so the experiment isolates the wind's
# response.  (The chemistry that made the profile saw the fiducial
# spectrum; re-closing the loop at each scaling is a separate experiment.)
#
# usage: ./run_xuv_profile_case.sh <source_iteration_dir> <tag> <sed_file> [n]
set -u
src=$1; tag=$2; sed_file=$3; nrep=${4:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
here=$(dirname "$0"); cd "$here" || exit 1

mkdir -p "$tag/output"
sed -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    -e "s|^Spectrum file: .*|Spectrum file: $sed_file|" \
    -e "s/^Do only PP: True$/Do only PP: False/" \
    "$src/input.inp" > "$tag/input.inp"
\cp -f "$src/lower_atmosphere_profile.dat" "$tag/"
\cp -f "$src/base.inp" "$tag/" 2>/dev/null
\cp -f "$src/output/Hydro_ioniz.txt" "$tag/output/Hydro_ioniz_IC.txt"
\cp -f "$src/output/Ion_species.txt" "$tag/output/Ion_species_IC.txt"

cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
echo "[$tag] START from $src, SED=$sed_file  $(date +%H:%M:%S)"
for i in $(seq 1 "$nrep"); do
   [ $i -gt 1 ] && { \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
                     \cp -f output/Ion_species.txt output/Ion_species_IC.txt; }
   sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
   EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$BIN" > "run$i.log" 2>&1
   echo "[$tag] pass-set $i: $(grep -o 'done info=[0-9]*' run$i.log|tail -n 1)" \
        "$(grep -o '||R||= *[0-9.E+-]*' run$i.log|tail -n 1)"
done
\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$?  $(grep 'steady-state Mdot' pp.log|tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"
