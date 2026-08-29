#!/bin/bash
# Re-synthesize the He I 10830 line of an arm with the turbulence term on,
# the p-winds setting of Lampon et al. (2020): v_turb^2 = (5/6) kT/m added
# in quadrature. The wind is untouched; only the line synthesis changes, so
# the result lands in tpm_turb/ beside the nominal one, as ../regen_turb.sh
# does for the earlier scan.
#
# usage: ./run_turb.sh <tag> [tag ...]
set -u
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
. "$EX/LHS1140b/winered_hires_y.sh"   # He I 10830 at WINERED HIRES-Y, R = 68,000
cd "$(dirname "$0")" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
for tag in "$@"; do
  ( cd "$tag" || exit 1
    mkdir -p tpm_turb
    EXHALE_TRANSIT_TURB=1 EXHALE_TRANSIT_SAVE_PREFIX=tpm_turb/ MPLBACKEND=Agg \
      PYTHONPATH="$EX" python3 "$EX/EXHALE_transit.py" > tpm_turb/transit.log 2>&1 )
  echo "$tag: $(grep 'He 10830 metrics' "$tag/tpm_turb/transit.log" | tail -1 | sed 's/.*metrics: //;s/ ->.*//')"
done
echo ALL_DONE
