#!/bin/bash
# setup_photochem.sh -- put a Photochem source tree into the state EXHALE's
# lower-atmosphere profile handoff is run from: upstream v0.9.0 plus the five
# file changes of src/utils/photochem_exhale.patch.
#
# Photochem is a third-party code.  What EXHALE needs is not the release and
# not the conda package: Equilibrate's elemental mass-balance test accepted
# solves that left a trace element orders of magnitude out of balance, two of
# Clima's root solves ran unconstrained, and one of those differentiated its
# Jacobian at a step so small the difference was its own numerical noise.  The
# patch corrects all of it, adds the elemental-closure check on the Photochem
# side, and pins Clima v0.7.5.  Why each was needed and what it measured:
#   docs/photochem_solver_modification_investigation.md
#   docs/photochem_solver_modification_implementation.md
# What each change is, and how to build the wheel and the environment:
#   README_photochem.md      README_HOWTO.md "The Photochem environment"
#
# Usage:  src/utils/setup_photochem.sh              # act on EXHALE/photochem
#         src/utils/setup_photochem.sh <dir>        # clone into / patch <dir>
#         PHOTOCHEM_URL=<url> src/utils/setup_photochem.sh <dir>
#
# Idempotent: on an already-patched tree it verifies and exits 0 without
# touching a file.  It refuses, rather than half-applying, on a tree that is
# neither upstream v0.9.0 nor already patched.
#
# This only produces the source.  Building the wheel and the environment is a
# separate step, and its two traps on this machine (the broken
# /usr/include/numpy symlink, and the user-site NumPy that shadows every conda
# environment) are in README_HOWTO.md.
#
# Photochem  https://github.com/Nicholaswogan/photochem
#            Cite, as upstream asks, either Wogan et al. (2023),
#            doi:10.3847/PSJ/aced83, or Wogan et al. (2024),
#            doi:10.3847/2041-8213/ad2616.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
EXROOT="$(cd "$HERE/../.." && pwd)"
DEST="${1:-$EXROOT/photochem}"
PATCHFILE="$HERE/photochem_exhale.patch"
URL="${PHOTOCHEM_URL:-https://github.com/Nicholaswogan/photochem.git}"

# v0.9.0, which is where origin/main pointed when this was written.  origin/dev
# is not a release and still carries the solver behavior the patch corrects.
UPSTREAM_TAG="v0.9.0"
UPSTREAM_COMMIT="e1e872528b61e8dd1db891b84869738e219b4f97"

say() { echo "[setup_photochem] $*"; }
die() { echo "[setup_photochem] ERROR: $*" >&2; exit 1; }

[ -f "$PATCHFILE" ] || die "no patch at $PATCHFILE"

# --------------------------------------------------------------------------
# 1. the source tree
# --------------------------------------------------------------------------
if [ ! -e "$DEST" ]; then
    say "cloning $URL -> $DEST"
    git clone "$URL" "$DEST" || die "clone failed (network?)"
    git -C "$DEST" checkout --quiet "$UPSTREAM_COMMIT" \
        || die "$UPSTREAM_COMMIT ($UPSTREAM_TAG) is not in $URL"
fi
[ -d "$DEST" ] || die "$DEST exists and is not a directory"
[ -f "$DEST/photochem/extensions/gasgiants.py" ] \
    || die "$DEST does not look like a Photochem source tree"

# --------------------------------------------------------------------------
# 2. is it the commit the patch was cut against?
# --------------------------------------------------------------------------
# EXHALE/photochem is not tracked, and the copy standing there carries no .git
# -- it is the patched state, kept diffable against upstream and nothing else --
# so a missing history is reported and not treated as a failure.  Step 3 decides
# either way.  It sits inside EXHALE's working tree, so the history has to be
# one whose top level IS this directory; EXHALE's would otherwise answer for it.
toplevel="$(git -C "$DEST" rev-parse --show-toplevel 2>/dev/null || true)"
if [ -n "$toplevel" ] && [ "$(cd "$toplevel" && pwd -P)" = "$(cd "$DEST" && pwd -P)" ]; then
    head_commit="$(git -C "$DEST" rev-parse HEAD)"
    if [ "$head_commit" != "$UPSTREAM_COMMIT" ]; then
        die "$DEST is at $head_commit, not $UPSTREAM_COMMIT ($UPSTREAM_TAG).
     The patch is cut against that commit and will not be applied to another.
     git -C $DEST checkout $UPSTREAM_COMMIT"
    fi
    say "$DEST is at $UPSTREAM_COMMIT ($UPSTREAM_TAG)"
else
    say "$DEST carries no git history; deciding from the files themselves"
fi

# --------------------------------------------------------------------------
# 3. apply, or recognize that it is already applied
# --------------------------------------------------------------------------
cd "$DEST"
# `git apply` resolves the patch's paths against the top of whatever repository
# it finds, and EXHALE/photochem sits inside EXHALE's.  Stop the search
# at the parent so the paths are resolved against $DEST, which is what they are
# written relative to.  A tree with a history of its own is found first and is
# unaffected.
export GIT_CEILING_DIRECTORIES="$(cd .. && pwd -P)"
if git apply --reverse --check --whitespace=nowarn "$PATCHFILE" >/dev/null 2>&1; then
    say "already patched; nothing to do"
elif git apply --check --whitespace=nowarn "$PATCHFILE" >/dev/null 2>&1; then
    git apply --whitespace=nowarn "$PATCHFILE" || die "the patch checked clean but did not apply"
    say "applied $PATCHFILE"
else
    die "$DEST is neither upstream $UPSTREAM_TAG nor already patched: the patch
     applies in neither direction.  Start from a clean clone:
     rm -rf $DEST && $0 $DEST"
fi

# --------------------------------------------------------------------------
# 4. verify what the patch was for
# --------------------------------------------------------------------------
# Each of these is a change the handoff depends on, checked in the file it
# lives in rather than inferred from the patch having exited 0.
fail=0
need() {   # need <file> <string> <what it is>
    if grep -qF -- "$2" "$1" 2>/dev/null; then
        echo "  ok      $3"
    else
        echo "  MISSING $3   ($1)"; fail=1
    fi
}
needfile() {
    if [ -f "$1" ]; then echo "  ok      $2"; else echo "  MISSING $2   ($1)"; fail=1; fi
}

echo "[setup_photochem] verifying:"
g=photochem/extensions/gasgiants.py
need "$g" "equilibrium_mass_tol"         "gas-giant equilibrium tolerance is settable"
need "$g" "closure_rtol"                 "elemental closure is checked after each solve"
need "$g" "molfracs_atoms_condensate"    "the closure check is skipped where a condensate is present"
need "$g" "molfracs_atoms_sun.copy()"    "the solar elemental vector is copied before scaling"
c=src/dependencies/CMakeLists.txt
need "$c" 'GIT_TAG "v0.7.5"'             "Clima pinned to v0.7.5"
need "$c" "equilibrate-element-relative-mass-closure.patch"  "Equilibrate patch wired into CPM"
need "$c" "clima-bounded-scaled-solvers.patch"               "Clima patch wired into CPM"
needfile src/dependencies/patches/equilibrate-element-relative-mass-closure.patch \
                                         "Equilibrate patch file present"
needfile src/dependencies/patches/clima-bounded-scaled-solvers.patch \
                                         "Clima patch file present"
need src/dependencies/patches/clima-bounded-scaled-solvers.patch \
     "real(dp), parameter :: epsfcn = 1.0e-8_dp" \
                                         "Clima surface-temperature solve differences at a set step"
need tests/test_python.py "test_gas_giant_equilibrium_mass_tolerance_validation" \
                                         "tolerance-validation test present"
[ $fail -eq 0 ] || die "the tree is not in the expected state (above)"

say "$DEST is the source EXHALE runs on"
say "next: build the wheel and the environment -- README_HOWTO.md,"
say "      \"The Photochem environment\".  This script does not build."
