#!/bin/bash
# One arm of the elemental-flux-closure reservoir bracket, run on the
# corrected Photochem build (photochem_090_fix) with the
# current binary.  The closure iteration is executed in full -- chemistry,
# wind, elemental-flux measurement, under-relaxed update -- and is not a
# re-convergence on a frozen profile.
#
# The fixed configuration is ../flux_closure/heh9p7 in every respect but the
# reservoir He/H; the seed and the starting elemental fluxes are that arm's
# converged iterate, so no reservoir step exceeds 9.71/8.8 = 1.10x (the
# 1.50x jump of Update_EXHALE section 83 is what this avoids).
#
# usage: ./run_closure.sh <tag> <reservoir He/H> [seed_output] [phi0_H] [phi0_He]
#   TOL=<value>  closure tolerance on the undamped elemental-flux residual
#                (default 0.05, the value the stored ladder was run at)
#   KMAX=<n>     iteration ceiling (default 8)
set -u
tag=$1; resv=$2
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
FC=$EX/LHS1140b/exhale/flux_closure
seed=${3:-$FC/heh9p7/k03/output}
phiH=${4:-4.7682069882E+06}
phiHe=${5:-2.0433648550E+07}
DRV=$EX/src/utils/element_flux_closure.py
. "$EX/LHS1140b/winered_hires_y.sh"          # He I 10830 at R = 68,000
cd "$(dirname "$0")" || exit 1
HERE=$PWD
d=$HERE/$tag
mkdir -p "$d"

sed -e "s/He=9.7/He=$resv/" \
    -e "s|$FC/heh9p7/input_template.inp|$d/input_template.inp|" \
    -e "s/He\/H = 9.7 as the STARTING/He\/H = $resv as the STARTING/" \
    -e "s|/photochem_cmp/bin/python|/photochem_090_fix/bin/python|" \
    "$FC/heh9p7/closure.json" > "$d/closure.json"
sed -e "s|^He/H number ratio: .*|He/H number ratio: $resv|" \
    -e "s/^Planet name: .*/Planet name: LHS1140b_$tag/" \
    "$FC/heh9p7/input_template.inp" > "$d/input_template.inp"

echo "[$tag] CLOSURE START reservoir=$resv seed=$seed  $(date +%H:%M:%S)"
resume=""; [ -f "$d/closure_history.txt" ] && resume="--resume"
python3 "$DRV" "$d" --phi0-H "$phiH" --phi0-He "$phiHe" \
    --config "$d/closure.json" --seed "$seed" $resume \
    --tol ${TOL:-0.05} --kmax ${KMAX:-8} > "$HERE/$tag.closure.log" 2>&1
rc=$?
last=$(\ls -d --color=never $d/k[0-9][0-9] 2>/dev/null | sort | tail -n 1)
echo "[$tag] CLOSURE rc=$rc last=$last"
[ $rc -eq 0 ] && [ -n "$last" ] || { echo "[$tag] CLOSURE FAILED"; exit 1; }
( cd "$last" && MPLBACKEND=Agg PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" \
    > transit.log 2>&1 )
echo "[$tag] TRANSIT rc=$? in $last  $(date +%H:%M:%S)"
