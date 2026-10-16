#!/bin/bash
# Runs the six Yan et al. (2022) Table 4 points (H/He 92/8 and 98/2 at
# 0.5, 1, 2 F0, beta_m = 0.22) with the frozen binary, at most MAX_RUNS at
# once (8 OpenMP threads each). A case directory that already holds run.log
# is taken as started and skipped. Start/finish lines go to run_table4_grid.log.
BIN=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.P2/EXHALE_f8a38866.x
HERE=$(cd "$(dirname "$0")" && pwd)
MAX_RUNS=3
CASES="h92he8_f1 h98he2_f1 h92he8_f0p5 h98he2_f0p5 h92he8_f2 h98he2_f2"

running_here() {
    # count binaries whose working directory is a case directory of this benchmark
    local p
    for p in $(pgrep -f "^$BIN"); do
        readlink /proc/$p/cwd 2>/dev/null
    done | grep -c "^$HERE/"
}

for c in $CASES; do
    [ -f "$HERE/$c/run.log" ] && continue
    while [ "$(running_here)" -ge "$MAX_RUNS" ]; do sleep 60; done
    ( cd "$HERE/$c" && OMP_NUM_THREADS=8 nohup "$BIN" > run.log 2>&1 & )
    sleep 5
    pid=""
    for p in $(pgrep -f "^$BIN"); do
        [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$HERE/$c" ] && pid=$p
    done
    echo "$(date '+%F %T %Z') start $c pid $pid" >> "$HERE/run_table4_grid.log"
done
echo "$(date '+%F %T %Z') all cases started" >> "$HERE/run_table4_grid.log"
