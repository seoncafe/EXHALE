#!/bin/bash
# Re-solve one stored LHS 1140 b arm on its own output, composition unchanged,
# so that the solution carries the current He 2^3S + H ionization coefficient
# (Update_EXHALE section 87) and the current advection post-process
# (section 88).  Copied from crossings_gm25/resolve_wind.sh and generalized:
# the source arms sit one directory level higher, so the relative SED path
# gains one "..", and the lower-atmosphere profile is copied only if the
# source has one (these six arms use the scalar base).
#
# usage: ./resolve_wind.sh <src_dir> <tag> [n_passes]
set -u
src=$1; tag=$2; n=${3:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"
BIN=$EX/EXHALE.x
cd "$(dirname "$0")" || exit 1
rm -rf "$tag"; mkdir -p "$tag/output"
cp "$src/input.inp" "$tag/"
[ -f "$src/lower_atmosphere_profile.dat" ] && cp "$src/lower_atmosphere_profile.dat" "$tag/"
[ -f "$src/metals.inp" ] && cp "$src/metals.inp" "$tag/"
cp "$src/output/Hydro_ioniz.txt" "$src/output/Ion_species.txt" "$tag/output/"
cd "$tag" || exit 1
sed -i -e 's|^Spectrum file: \.\./\.\./sed/|Spectrum file: ../../../sed/|' \
       -e 's/^Do only PP: True$/Do only PP: False/' input.inp
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
echo "[$tag] START src=$src threads=$OMP_NUM_THREADS $(date +%H:%M:%S)"
for i in $(seq 1 "$n"); do
  \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
  \cp -f output/Ion_species.txt output/Ion_species_IC.txt
  EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$BIN" > "cont$i.log" 2>&1
  echo "[$tag] pass $i: $(grep -o 'done info=[0-9]*' cont$i.log | tail -n 1)" \
       "$(grep -o '||R||= *[0-9.E+-]*' cont$i.log | tail -n 1)" "$(date +%H:%M:%S)"
done
\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] done $(grep 'steady-state Mdot' pp.log | tail -n 1) $(date +%H:%M:%S)"
