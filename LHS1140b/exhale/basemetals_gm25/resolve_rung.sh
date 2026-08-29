#!/bin/bash
# Solve one rung of the one-key-at-a-time control ladder of
# docs/lhs1140b_exhale_vs_pwinds.tex Sect. sec:basemetals on the current
# binary, seeded from the previous rung's converged output, and synthesize
# the He 10830 line twice: at R = 68,000 (the WINERED HIRES-Y kernel the
# memo measures through) and at R = 80,000 (the kernel the stored ladder of
# Update_EXHALE Sect. 78 was measured at, kept only for lineage).
#
# The rung directory must already carry input.inp and, where the rung has
# one, base.inp / lower_atmosphere_profile.dat.  No directory outside this
# one is written.
#
# usage: ./resolve_rung.sh <tag> <seed_dir> [n_passes]
set -u
tag=$1; seed=$2; n=${3:-4}
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
BIN=$EX/EXHALE.x
cd "$(dirname "$0")" || exit 1
seed=$(cd "$seed" && pwd)
mkdir -p "$tag/output"
# seeding from the rung's own output (a re-solve at a tighter tolerance) is a
# no-op copy, so it is skipped rather than refused
if [ "$seed" != "$(cd "$tag" && pwd)" ]; then
  cp "$seed/output/Hydro_ioniz.txt" "$seed/output/Ion_species.txt" "$tag/output/"
fi
cd "$tag" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
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
# R = 80,000: EXHALE_transit.py's built-in He I 10830 default, unset here
env -u EXHALE_TRANSIT_RES_HETR MPLBACKEND=Agg PYTHONPATH="$EX" \
    python3 "$EX/EXHALE_transit.py" > transit_R80k.log 2>&1
for f in tpm_He10830 tpm_He10830_metrics tpm_Halpha tpm_Hbeta tpm_Lya; do
  [ -f "$f.txt" ] && \cp -f "$f.txt" "${f}_R80k.txt"
done
# R = 68,000: the WINERED HIRES-Y kernel of the observation, the memo default
. "$EX/LHS1140b/winered_hires_y.sh"
MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > transit.log 2>&1
echo "[$tag] done $(grep 'steady-state Mdot' pp.log | tail -n 1)"
