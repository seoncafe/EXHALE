#!/usr/bin/env bash
# Reproduce the solver-upgrade validation (Update_EXHALE_solver.tex Task 2 / Task 1).
#
# Runs a post-process-only sweep (Do only PP) over the converged Case B,
# feeding the converged eq snapshot as the initial condition, with:
#   newton/  -- the analytic-Jacobian Newton (default)
#   hybrd1/  -- the MINPACK hybrd1 reference (ATES_FORCE_HYBRD1=1, same binary)
# and copies the converged Case B as baseline/. Then compare.py writes
# comparison.txt. The Newton-vs-hybrd1 diff (identical one-step state) isolates
# the solver; the Newton-vs-baseline diff includes the one-step hydro drift.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
WASP="$(cd "$HERE/.." && pwd)"
BIN="$WASP/EXHALE.x"

run_variant () {            # $1 = label, $2 = force-hybrd1 (0/1)
  local label="$1" force="$2" d="$HERE/$1"
  rm -rf "$d"; mkdir -p "$d/output"
  cp -f "$BIN" "$WASP/metals.inp" "$d/"
  cp -f "$WASP/input.inp" "$d/input.inp"
  cp -f "$WASP"/output/*.txt "$d/output/"
  # The shipped *_IC.txt are stale -> feed the converged eq snapshot as the IC.
  cp -f "$WASP/output/Hydro_ioniz.txt" "$d/output/Hydro_ioniz_IC.txt"
  cp -f "$WASP/output/Ion_species.txt" "$d/output/Ion_species_IC.txt"
  sed -i 's/^Load IC? False/Load IC? True/; s/^Do only PP: False/Do only PP: True/' "$d/input.inp"
  if [ "$force" = 1 ]; then
    ( cd "$d" && ATES_FORCE_HYBRD1=1 ./EXHALE.x > run.log 2>&1 )
  else
    ( cd "$d" && ./EXHALE.x > run.log 2>&1 )
  fi
  echo "[$label] $(grep -i 'ioniz-eq solver' "$d/run.log" || echo 'no solver report')"
}

run_variant newton 0
run_variant hybrd1 1
rm -rf "$HERE/baseline"; mkdir -p "$HERE/baseline"
cp -f "$WASP"/output/*.txt "$HERE/baseline/"

python3 "$HERE/compare.py" > "$HERE/comparison.txt"
echo "----"
cat "$HERE/comparison.txt"
