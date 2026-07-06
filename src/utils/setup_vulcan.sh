#!/bin/bash
# setup_vulcan.sh -- obtain and prepare the VULCAN photochemistry code for the
# EXHALE lower-atmosphere pre-step.  VULCAN (which ships FastChem inside it) is
# NOT committed to the EXHALE repository; this script clones it into
# EXHALE/VULCAN/ and applies the two small modifications EXHALE needs.
#
# Usage:   src/utils/setup_vulcan.sh            # clone into EXHALE/VULCAN/
#          src/utils/setup_vulcan.sh <dir>      # clone into <dir> instead
#          VULCAN_URL=<url> src/utils/setup_vulcan.sh   # alternate remote
#
# Idempotent: re-running only (re)applies the patches and rebuilds FastChem.
#
# Sources / citation (please cite in any publication using the pre-step):
#   VULCAN   https://github.com/exoclime/VULCAN  (mirror shami-EEG/VULCAN)
#            Tsai et al. 2017 (ApJS 228, 20); Tsai et al. 2021 (ApJ 923, 264)
#   FastChem https://github.com/NewStrangeWorlds/FastChem  (bundled in VULCAN)
#            Stock et al. 2018 (MNRAS 479, 865); Stock et al. 2022 (MNRAS 517, 4070)
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
EXROOT="$(cd "$HERE/../.." && pwd)"
DEST="${1:-$EXROOT/VULCAN}"
URL="${VULCAN_URL:-https://github.com/exoclime/VULCAN.git}"

# 1. clone (skip if already present)
if [ ! -d "$DEST" ]; then
    echo "[setup_vulcan] cloning $URL -> $DEST"
    git clone --depth 1 "$URL" "$DEST"
else
    echo "[setup_vulcan] $DEST already exists -- applying patches only"
fi

# 2. MODIFICATION 1: Python-3 bytes-vs-str fix in the element-conservation
#    check (np.genfromtxt must be told encoding=None on numpy>=1.14).
f="$DEST/make_chem_funs.py"
if grep -q "com_file,names=True,dtype=None)" "$f" 2>/dev/null; then
    sed -i 's/com_file,names=True,dtype=None)/com_file,names=True,dtype=None,encoding=None)/' "$f"
    echo "[setup_vulcan] patched make_chem_funs.py (encoding=None)"
fi

# 3. MODIFICATION 2: config defaults for the EXHALE pre-step.  vulcan_cfg.py is
#    the driver's template; it must have photochemistry on and live-plot off.
c="$DEST/vulcan_cfg.py"
if [ -f "$c" ]; then
    sed -i 's/^use_photo *=.*/use_photo = True/'            "$c"
    sed -i 's/^use_live_plot *=.*/use_live_plot = False/'   "$c"
    echo "[setup_vulcan] set vulcan_cfg.py: use_photo=True, use_live_plot=False"
fi

# 4. build FastChem (shipped inside VULCAN as fastchem_vulcan/)
if [ -d "$DEST/fastchem_vulcan" ]; then
    echo "[setup_vulcan] building FastChem ..."
    ( cd "$DEST/fastchem_vulcan" && make >/dev/null 2>&1 ) \
        && echo "[setup_vulcan] FastChem built" \
        || echo "[setup_vulcan] WARNING: FastChem build failed -- build it by hand in $DEST/fastchem_vulcan"
fi

# 5. drop any run products so the tree is clean
rm -f "$DEST"/output/*.vul 2>/dev/null || true

echo "[setup_vulcan] done. EXHALE can now use 'Lower atmosphere: vulcan <R_1bar>'."
