#!/bin/bash
# Re-solve one stored run's wind on the current binary, seeded from its own
# stored output, with everything else -- composition, spectrum, base, and any
# lower-atmosphere profile -- held exactly as the source directory carries it.
# The source directory is read only; nothing in it is written.
#
# usage: ./resolve_wind.sh <src_dir> <tag> [n_passes]
set -u
src=$1; tag=$2; n=${3:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"          # He I 10830 at R = 68,000
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
here=$(cd "$(dirname "$0")" && pwd)
srcabs=$(cd "$src" && pwd)
cd "$here" || exit 1
rm -rf "$tag"; mkdir -p "$tag/output"
cp "$srcabs/input.inp" "$tag/"
for f in lower_atmosphere_profile.dat base.inp metals.inp opacity.inp; do
  [ -f "$srcabs/$f" ] && cp "$srcabs/$f" "$tag/"
done
cp "$srcabs/output/Hydro_ioniz.txt" "$srcabs/output/Ion_species.txt" "$tag/output/"
cd "$tag" || exit 1
# the spectrum path is relative to the source directory; make it absolute
spec=$(grep '^Spectrum file:' input.inp | sed 's|^Spectrum file: *||')
case "$spec" in
  /*) ;;
  *) sed -i "s|^Spectrum file: .*|Spectrum file: $(cd "$srcabs" && cd "$(dirname "$spec")" && pwd)/$(basename "$spec")|" input.inp ;;
esac
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
mkdir -p tpm_turb
EXHALE_TRANSIT_TURB=1 EXHALE_TRANSIT_SAVE_PREFIX=tpm_turb/ MPLBACKEND=Agg \
  PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > tpm_turb/transit.log 2>&1
echo "[$tag] done $(grep 'steady-state Mdot' pp.log | tail -n 1)"
