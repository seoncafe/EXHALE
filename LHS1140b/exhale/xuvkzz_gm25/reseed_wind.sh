#!/bin/bash
# The same wind solve as ./resolve_wind.sh, seeded from a different solution
# instead of from the arm's own stored output.  The configuration -- the
# spectrum, the composition, K_zz, the lower atmosphere -- is the target
# arm's; only the initial condition differs.  Two lineages measured this way
# bound the seed dependence of a grid point.
#
# usage: ./reseed_wind.sh <config_src> <seed_src> <tag> [n_passes]
set -u
cfg=$1; seed=$2; tag=$3; n=${4:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"
BIN=$EX/EXHALE.x
here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
rm -rf "$tag"; mkdir -p "$tag/output"
cp "$cfg/input.inp" "$tag/"
for f in lower_atmosphere_profile.dat base.inp metals.inp opacity.inp; do
  [ -f "$cfg/$f" ] && cp "$cfg/$f" "$tag/"
done
cp "$seed/output/Hydro_ioniz.txt" "$seed/output/Ion_species.txt" "$tag/output/"
cd "$tag" || exit 1
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
