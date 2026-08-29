#!/bin/bash
# One arm of the well-mixed (no diffusion, scalar base, no metals, GJ 1132 SED)
# He/H bracket, solved with the current binary.  The configuration is
# ../heh0p55 in every respect but the He/H reservoir ratio; that directory is
# read only, as the seed of the restart.
#
# usage: ./run_wellmixed.sh <tag> <He/H> [seed_dir]
#   seed_dir defaults to ../heh0p55 (its converged output/)
set -u
tag=$1; heh=$2; seed=${3:-../heh0p55}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"          # He I 10830 at R = 68,000
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
cd "$(dirname "$0")" || exit 1
ref=../heh0p55

mkdir -p "$tag/output"
sed -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    -e "s|^He/H number ratio: .*|He/H number ratio: $heh|" \
    -e "s/^Do only PP: True$/Do only PP: False/" \
    -e "s|^Spectrum file: .*|Spectrum file: ../$(grep '^Spectrum file:' $ref/input.inp | sed 's|^Spectrum file: ||')|" \
    "$ref/input.inp" > "$tag/input.inp"
\cp -f "$seed/output/Hydro_ioniz.txt" "$tag/output/Hydro_ioniz_IC.txt"
\cp -f "$seed/output/Ion_species.txt" "$tag/output/Ion_species_IC.txt"

cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
echo "[$tag] JFNK START He/H=$heh seed=$seed  $(date +%H:%M:%S)"
EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 "$BIN" > run.log 2>&1
info=$(grep -o 'done info=[0-9]*' run.log | tail -n 1)
echo "[$tag] JFNK DONE $info $(grep -o '||R||= *[0-9.E+-]*' run.log | tail -n 1)"
grep -q 'done info=0' run.log || { echo "[$tag] JFNK FAILED"; exit 1; }

\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$? $(grep 'steady-state Mdot' pp.log | tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"
