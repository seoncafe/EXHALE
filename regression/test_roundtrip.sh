#!/bin/bash
# Restart round-trip test (Phase-2 gate).
#
# Verifies that load_IC restores the species state written by write_output:
#  [A] schema-2 (headered) IC  -> ALL species incl. metal ions restored
#  [B] legacy (headerless) IC  -> H/He restored, metals fall back to
#                                 neutral-from-abundance (historical behavior)
#
# Uses the ATES_DUMP_IC=1 hook: ATES loads the IC, writes the state as
# output/*.txt BEFORE any solver step, and stops. The comparison is then a
# pure loader test (the per-step ionization-equilibrium solve would otherwise
# re-equilibrate metals and mask a reset-to-neutral bug).

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
SRC="$HERE/wasp_full"
RT="$HERE/roundtrip"

rm -rf "$RT"; mkdir -p "$RT/output"
sed 's/^Load IC? False/Load IC? True/' "$SRC/input.inp" > "$RT/input.inp"
cp "$SRC/metals.inp" "$RT/metals.inp"
cp "$ROOT/ATES.x"    "$RT/ATES.x"

echo "[A] schema-2 (headered) IC round-trip"
cp "$SRC/output/Hydro_ioniz.txt" "$RT/output/Hydro_ioniz_IC.txt"
cp "$SRC/output/Ion_species.txt" "$RT/output/Ion_species_IC.txt"
( cd "$RT" && ATES_DUMP_IC=1 OMP_NUM_THREADS=1 ./ATES.x > dump_A.log 2>&1 )
python3 "$HERE/check_roundtrip.py" "$RT/output/Ion_species_IC.txt" \
                                   "$RT/output/Ion_species.txt" schema2

echo "[B] legacy (headerless) IC fallback"
grep -v '^ *#' "$SRC/output/Hydro_ioniz.txt" > "$RT/output/Hydro_ioniz_IC.txt"
grep -v '^ *#' "$SRC/output/Ion_species.txt" | awk '{print $1,$2,$3,$4,$5,$6,$7}' \
    > "$RT/output/Ion_species_IC.txt"
( cd "$RT" && ATES_DUMP_IC=1 OMP_NUM_THREADS=1 ./ATES.x > dump_B.log 2>&1 )
python3 "$HERE/check_roundtrip.py" "$SRC/output/Ion_species.txt" \
                                   "$RT/output/Ion_species.txt" legacy

echo "==> ROUNDTRIP PASS"
