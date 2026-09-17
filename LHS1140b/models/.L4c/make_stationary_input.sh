#!/bin/bash
# Write the stationary input of a molecular case from its marching input:
# the same physics, read back as a state (Load IC? True), rebuilt under PLM,
# no marching threshold, and entered into the partitioned stationary solve.
#   make_stationary_input.sh <marching input.inp> <output input.inp> [extra key]...
src=$1; dst=$2; shift 2
sed -e 's/^Reconstruction scheme:.*/Reconstruction scheme: PLM/' \
    -e 's/^Load IC?.*/Load IC? True/' \
    -e '/^du_th /d' -e '/^Restart intent:/d' -e '/^Secondary_ionization:/d' \
    "$src" > "$dst"
{ echo 'Restart intent: stationary'; echo 'Secondary_ionization: Immediate'; } >> "$dst"
for k in "$@"; do echo "$k" >> "$dst"; done
