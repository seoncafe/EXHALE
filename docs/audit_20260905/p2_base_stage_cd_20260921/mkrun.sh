#!/bin/bash
# Build one isolated run directory from the atomic checkpoint.
# usage: mkrun.sh <name>
set -e
S=/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/p2_base_cd
R=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
C=$R/LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7
G=$C/states/g0002_20260919T004843Z_f7485b14
D=$S/LHS1140b/models/case/$1
rm -rf $D
mkdir -p $D/output
\cp -f $C/input.inp $D/input.inp
\cp -f $G/Hydro_ioniz.txt $D/output/Hydro_ioniz_IC.txt
\cp -f $G/Ion_species.txt $D/output/Ion_species_IC.txt
\cp -f $S/EXHALE.x $D/EXHALE.x
echo $D
