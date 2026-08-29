#!/bin/bash
# Every stored run behind the six memo statements of this directory, re-solved
# on the current binary.  Three at a time, one thread each.
cd "$(dirname "$0")" || exit 1
export OMP_NUM_THREADS=1
run() { ./resolve_wind.sh "$1" "$2" 4 > "log_$2.txt" 2>&1; }
set -- \
  "../heh0p55                     heh0p55" \
  "../solar                       solar" \
  "../heh0p55_diff_ctrl           kzz0" \
  "../heh0p55_diff_kzz1e8         kzz1e8" \
  "../heh0p55_diff_kzz1e9         kzz1e9" \
  "../heh0p55_diff_kzz1e10        kzz1e10" \
  "../heh0p04_gj699               gj699_heh0p04" \
  "../heh0p06_gj699               gj699_heh0p06" \
  "../solar_gj699                 gj699_solar" \
  "../heh0p25_gj699               gj699_heh0p25" \
  "../heh0p55_gj699               gj699_heh0p55" \
  "../heh1_gj699                  gj699_heh1" \
  "../heh1000_gj699               gj699_heh1000" \
  "../flux_closure/ref/k00        fc_ref" \
  "../flux_closure/lo/k05         fc_lo" \
  "../flux_closure/hi/k06         fc_hi" \
  "../flux_closure/heh11p1/k01    fc_heh11p1"
for a in "$@"; do
  set -- $a
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do wait -n; done
  echo "START $2 $(date +%H:%M:%S)"
  run "$1" "$2" &
done
wait
echo ALL_DONE
