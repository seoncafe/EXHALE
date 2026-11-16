#!/usr/bin/env bash
# Regenerate inputdata/windae_grid/manifest.csv for the Wind-AE nearest-seed
# picker (wae_continuation.f90::pick_nearest_seed). The manifest holds one
#     <path>,<Mp[g]>,<Rp[cm]>,<Ftot[erg/cm^2/s]>
# row per converged seed CSV, parsed from each file's "#plnt_prms:" line.
#
# The grid itself (inputdata/windae_grid/, copied from the Broome et al. 2025
# solution library) and this manifest are .gitignored -- bulky, re-creatable.
# Run this from the EXHALE top directory after:
#   * first copying the grid into inputdata/windae_grid/, or
#   * dumping new seeds there via "Wind-AE seed out:" (so the picker sees them).
#
# Usage:  src/utils/make_windae_manifest.sh [grid_dir]   (default inputdata/windae_grid)
set -euo pipefail

GD="${1:-inputdata/windae_grid}"
out="$GD/manifest.csv"

if [ ! -d "$GD" ]; then
   echo "make_windae_manifest: no grid directory '$GD'" >&2
   echo "  copy the Broome grid there first, e.g.:" >&2
   echo "  cp -r <wind-ae>/Notebooks/.../data/Grid/* $GD/" >&2
   exit 1
fi

tmp="$out.tmp"
# path,Mp,Rp,Ftot  (plnt_prms = Mp,Rp,Mstar,a,Ftot,Lstar -> keep fields 1,2,5)
grep -rH '#plnt_prms:' "$GD" --include='*.csv' 2>/dev/null \
 | grep -v '/manifest\.csv:' \
 | sed 's/:#plnt_prms:/,/' \
 | awk -F',' '{gsub(/ /,"",$2); gsub(/ /,"",$3); gsub(/ /,"",$6);
               print $1","$2","$3","$6}' \
 | sort > "$tmp"

mv "$tmp" "$out"
echo "make_windae_manifest: wrote $(wc -l < "$out") rows -> $out"
