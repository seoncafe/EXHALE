#!/bin/bash
# Re-solve one stored arm's wind on the current binary, with its lower
# atmosphere held fixed, and re-synthesize the He 10830 line.
#
# Same recipe as ../crossings_gm25/resolve_wind.sh: seed the JFNK pass from
# the arm's own output, repeat until the residual stops moving, then one
# post-process pass and one transit synthesis.  The composition, the XUV
# scaling and K_zz are whatever the source arm carried -- nothing here
# changes them.  The source directory is read only.
#
# usage: ./resolve_wind.sh <src_dir> <tag> [n_passes]
set -u
src=$1; tag=$2; n=${3:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"
BIN=$EX/EXHALE.x
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
[ -d "$src" ] || src="$EX/LHS1140b/exhale/$src"
rm -rf "$tag"; mkdir -p "$tag/output"
cp "$src/input.inp" "$tag/"
for f in lower_atmosphere_profile.dat base.inp metals.inp opacity.inp; do
  [ -f "$src/$f" ] && cp "$src/$f" "$tag/"
done
cp "$src/output/Hydro_ioniz.txt" "$src/output/Ion_species.txt" "$tag/output/"
cd "$tag" || exit 1
# the SED path is relative in some source trees; make it absolute so the run
# directory can sit at any depth
sed -i -E 's#^(Spectrum file: *)\.\./\.\./sed/#\1'"$EX"'/LHS1140b/sed/#' input.inp
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
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
