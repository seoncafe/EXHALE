#!/bin/bash
# Re-synthesize the transit spectra with turbulence broadening on, matching
# the p-winds setting the authors use. Winds are untouched; only the line
# synthesis changes, so results land in tpm_turb/ beside the nominal ones.
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd "$(dirname "$0")"
export OMP_NUM_THREADS=6
for c in solar heh1 heh10 heh100 heh1000 heh10000; do
  ( cd $c
    mkdir -p tpm_turb
    EXHALE_TRANSIT_TURB=1 MPLBACKEND=Agg PYTHONPATH="$EX" \
      EXHALE_TRANSIT_SAVE_PREFIX=tpm_turb/ \
      python3 "$EX/EXHALE_transit.py" > tpm_turb/transit.log 2>&1 )
  echo "$c: $(grep 'He 10830 metrics' $c/tpm_turb/transit.log | tail -1 | sed 's/.*metrics: //')"
done
echo ALL_DONE
