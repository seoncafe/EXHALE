#!/bin/bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
S=$EX/scratchpad/L6c
M=$S/tree_M/EXHALE_L6c_M.x
C=$S/tree_C/EXHALE_L6c_C.x
# longest first; at most 6 at a time
JOBS="wasp_full:C wasp_full:M wasp_he23off:C wasp_he23off:M wasp_full_newton:C wasp_full_newton:M \
      lower_profile:C lower_profile:M mol_metals:C mol_metals:M mol_ir_bands:C mol_ir_bands:M \
      oxygen_chemistry:C oxygen_chemistry:M"
for j in $JOBS; do
  c=${j%%:*}; t=${j##*:}
  while [ $(jobs -rp | wc -l) -ge 6 ]; do sleep 10; done
  if [ $t = C ]; then e=$C; else e=$M; fi
  $S/run_case.sh $c $e $t &
done
wait
echo ALLDONE > $S/runs/.done
