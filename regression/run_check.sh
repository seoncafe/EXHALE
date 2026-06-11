#!/bin/bash
# Regression check for no-physics-change refactor steps.
#
# Usage:
#   ./run_check.sh golden [case...]   # snapshot current <case>/output as reference
#   ./run_check.sh check  [case...]   # rebuild, re-run each case single-thread,
#                                     #   bitwise-compare vs its golden
#
# Single-threaded (OMP_NUM_THREADS=1) => fully deterministic => bitwise compare.
# A no-physics-change refactor must produce byte-identical Hydro_ioniz.txt and
# Ion_species.txt for EVERY case in the feature matrix.
#
# Matrix cases (add more dirs as the matrix grows):
#   wasp_full      He23S on,  metals on   (exercises HeITR + metal paths)
#   wasp_he23off   He23S off, metals on   (exercises the HeITR-off branch)

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
FILES="output/Hydro_ioniz.txt output/Ion_species.txt"
DEFAULT_CASES="wasp_full wasp_he23off"

mode="$1"; shift || true
CASES="${*:-$DEFAULT_CASES}"

case "$mode" in
  golden)
    for c in $CASES; do
      mkdir -p "$HERE/golden/$c"
      for f in $FILES; do cp "$HERE/$c/$f" "$HERE/golden/$c/$(basename $f)"; done
      cp "$HERE/$c/run.log" "$HERE/golden/$c/run.log" 2>/dev/null || true
      echo "golden[$c]: $(grep -E 'final:' "$HERE/$c/run.log" | tail -n1)"
    done
    ;;
  check)
    echo "[build] rebuilding ATES.x ..."
    ( cd "$ROOT" && make >/dev/null 2>&1 ) && echo "  build OK"
    rc=0
    for c in $CASES; do
      echo "[$c] running (OMP_NUM_THREADS=1) ..."
      cp "$ROOT/ATES.x" "$HERE/$c/ATES.x"
      ( cd "$HERE/$c" && rm -f output/*.txt && OMP_NUM_THREADS=1 ./ATES.x > run.log 2>&1 )
      echo "       $(grep -E 'final:' "$HERE/$c/run.log" | tail -n1)"
      for f in $FILES; do
        b=$(basename "$f")
        # Compare numeric content only: '#' schema-header lines are excluded,
        # so adding/changing headers never breaks the bitwise data guarantee.
        if cmp -s <(grep -v '^ *#' "$HERE/$c/$f") <(grep -v '^ *#' "$HERE/golden/$c/$b"); then
          echo "       PASS $b (data identical)"
        else
          echo "       FAIL $b"
          cmp <(grep -v '^ *#' "$HERE/$c/$f") <(grep -v '^ *#' "$HERE/golden/$c/$b") | head -n1
          rc=1
        fi
      done
    done
    [ $rc -eq 0 ] && echo "==> REGRESSION PASS (all cases byte-identical)" || echo "==> REGRESSION FAIL"
    exit $rc
    ;;
  *)
    echo "usage: $0 {golden|check} [case...]"; exit 2 ;;
esac
