#!/bin/bash
# One arm of the K_zz = 1e9 scalar-base He/H bracket, solved with the current
# binary.  The configuration is ../heh2p13_diff_kzz1e9 in every respect but
# the He/H reservoir ratio; that directory and the seed are read only.
#
# Like ../run_heh_diff_case.sh, the JFNK solve is repeated until the
# composition outer loop's drift stops moving, then one "Do only PP" pass and
# the transit synthesis.
#
# usage: ./run_kzz1e9.sh <tag> <He/H> [seed_dir] [n_pass_sets]
set -u
tag=$1; heh=$2; seed=${3:-../heh2_diff_kzz1e9}; nrep=${4:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"          # He I 10830 at R = 68,000
BIN=${EXHALE_BIN:-$EX/EXHALE.x}
cd "$(dirname "$0")" || exit 1
ref=../heh2p13_diff_kzz1e9

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
echo "[$tag] START He/H=$heh seed=$seed  $(date +%H:%M:%S)"
for i in $(seq 1 "$nrep"); do
   [ $i -gt 1 ] && { \cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
                     \cp -f output/Ion_species.txt output/Ion_species_IC.txt; }
   sed -i 's/^Do only PP: True$/Do only PP: False/' input.inp
   EXHALE_DIFFUSION_CHECK=1 EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 \
      "$BIN" > "run$i.log" 2> diffcheck.log
   echo "[$tag] pass-set $i: $(grep -o 'done info=[0-9]*' run$i.log | tail -n 1)" \
        "$(grep -o '||R||= *[0-9.E+-]*' run$i.log | tail -n 1)" \
        "$(grep -o 'composition drift = *[0-9.E+-]*' run$i.log | tail -n 1)"
done
grep -q 'done info=0' "run$nrep.log" || { echo "[$tag] JFNK FAILED"; exit 1; }

\cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
\cp -f output/Ion_species.txt output/Ion_species_IC.txt
sed -i 's/^Do only PP: False$/Do only PP: True/' input.inp
"$BIN" > pp.log 2>&1
echo "[$tag] PP DONE rc=$? $(grep 'steady-state Mdot' pp.log | tail -n 1)"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] TRANSIT DONE rc=$?  $(date +%H:%M:%S)"
