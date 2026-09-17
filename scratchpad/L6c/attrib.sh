#!/bin/bash
# Attribution table for one case:
#   (a) control vs golden   -- does reverting the two cooling files restore it
#   (b) measured vs golden  -- the movement the consolidated log reports
#   (c) measured vs control -- the movement the two cooling files alone produce
#   (d) measured vs the consolidated run's own output (EXHALE_C.x)
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
S=$EX/scratchpad/L6c; R=$EX/backup/regression
c="$1"
for f in Hydro_ioniz.txt Ion_species.txt Hydro_ioniz_adv.txt Ion_species_adv.txt; do
  cg=$($R/compare_within_tolerance.py $S/runs/C/$c/output/$f $R/golden/$c/$f 1e-3 2>&1)
  mg=$($R/compare_within_tolerance.py $S/runs/M/$c/output/$f $R/golden/$c/$f 1e-3 2>&1)
  mc=$($R/compare_within_tolerance.py $S/runs/M/$c/output/$f $S/runs/C/$c/output/$f 1e-3 2>&1)
  mx=$($R/compare_within_tolerance.py $S/runs/M/$c/output/$f $S/checkC_outputs/$c/$f 1e-3 2>&1)
  printf '%-22s %-24s\n' "$c/$f" ""
  printf '    control vs golden : %s\n' "$cg"
  printf '    measured vs golden: %s\n' "$mg"
  printf '    measured vs contrl: %s\n' "$mc"
  printf '    measured vs EXHALE_C.x run: %s\n' "$mx"
done
