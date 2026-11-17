#!/bin/bash
# setup_photochem.sh -- put a Photochem source tree into the state EXHALE's
# lower-atmosphere profile handoff runs on: upstream v0.9.3 plus the changes
# of src/utils/photochem_exhale.patch (four files).
#
# Photochem is a third-party code.  At v0.9.3 the Clima and Equilibrate
# sources are part of the Photochem tree (src/clima/, src/equilibrate/), so
# the patch edits them directly and no dependency patch files exist.  What
# EXHALE needs beyond the release:
#   - Equilibrate (src/equilibrate/equilibrate_cea.f90): the elemental
#     mass-balance test of both CEA update paths compares every positive
#     requested element with its own relative tolerance,
#     abs(b_0(i_atom))*mass_tol, instead of one threshold scaled by the most
#     abundant element with elements below 1e-6 exempt; and the intent(out)
#     flag `converged` of ec_UPDATE_ABUNDS_SHORT and ec_UPDATE_ABUNDS_LONG is
#     set .FALSE. on entry, so a failed test never returns a stale .TRUE.
#     from an optimized build.
#   - Clima (src/clima/clima_adiabat.f90): make_profile_bg_gas brackets the
#     background partial pressure on [1e-12, 1] x P_surf in log pressure,
#     skipping invalid states, and bisects (the sign test is nested so no
#     undefined residual is read); simple_solver solves for
#     log10(T_surf - T_trop) with MINPACK hybrd at epsfcn = 1e-8, divides the
#     flux and skin-temperature rows by their scales, and rejects a returned
#     point whose normalized residual exceeds 1e-7.
#   - photochem/extensions/gasgiants.py: EvoAtmosphereGasGiant and
#     GasGiantData take equilibrium_mass_tol (default 1e-12, validated), and
#     composition_at_metallicity checks the elemental closure of every
#     condensate-free level and raises if no attempt closes.
#   - tests/photochem/test_python.py: tests of the tolerance argument, the
#     trace-element closure, the closure in condensing and non-condensing
#     columns, and successive Equilibrate solves with name-property reads in
#     between.
#
# Requirements for building and running the result (measured with this
# version on Linux, Python 3.11, NumPy 1.26):
#   - a Fortran compiler accepted by upstream (gfortran >= 14), fypp, ninja;
#   - scikit-build-core (the build backend of pyproject.toml), used with
#     `python -m pip wheel <tree> --no-build-isolation`;
#   - CMake 3.x: CMake 4.x stops in the CVODE 5.7.0 subproject;
#   - with NumPy 1.26, Cython 3.0 (3.0.12 tested): Cython 3.3 generates
#     NumPy 2 C-API calls that NumPy 1.26 headers lack;
#   - photochem_clima_data 0.3.4 importable at run time (`import photochem`
#     imports it for DATA_DIR).  This version rejects mechanism files whose
#     gas thermodynamic polynomials are discontinuous at a join, which the
#     0.3.1 data are.  The upstream CMake build downloads
#     photochem_clima_data_v0.3.4.tar.gz into the tree root and unpacks its
#     data/ there for its own tests; that does not install the Python
#     package, which is installed separately (e.g. from the same tag,
#     https://github.com/Nicholaswogan/photochem_clima_data, v0.3.4).
#
# Usage:  src/utils/setup_photochem.sh [--url <url>] [<dir>]
#           <dir>  default: <EXHALE>/photochem; cloned if absent
#           --url  repository to clone from (default: upstream GitHub)
#
# Idempotent: on an already-patched tree it verifies and exits 0 without
# touching a file.  It refuses, rather than half-applying, a tree that is
# neither upstream v0.9.3 nor already patched.  It produces the source only;
# it does not build.
#
# Photochem  https://github.com/Nicholaswogan/photochem
#            Cite, as upstream asks, either Wogan et al. (2023),
#            doi:10.3847/PSJ/aced83, or Wogan et al. (2024),
#            doi:10.3847/2041-8213/ad2616.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
EXROOT="$(cd "$HERE/../.." && pwd)"
PATCHFILE="$HERE/photochem_exhale.patch"
URL="https://github.com/Nicholaswogan/photochem.git"
DEST=""
while [ $# -gt 0 ]; do
    case "$1" in
        --url) [ $# -ge 2 ] || { echo "--url needs a value" >&2; exit 2; }
               URL="$2"; shift 2 ;;
        -*)    echo "unknown option $1" >&2; exit 2 ;;
        *)     [ -z "$DEST" ] || { echo "one destination only" >&2; exit 2; }
               DEST="$1"; shift ;;
    esac
done
DEST="${DEST:-$EXROOT/photochem}"

UPSTREAM_TAG="v0.9.3"
UPSTREAM_COMMIT="ba3542cf63fc28f28dc9a0ad3a21d1dd718d9fa7"

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
    && [ -f "$DEST/src/equilibrate/equilibrate_cea.f90" ] \
    && [ -f "$DEST/src/clima/clima_adiabat.f90" ] \
    || die "$DEST does not look like a Photochem $UPSTREAM_TAG source tree"

# --------------------------------------------------------------------------
# 2. is it the commit the patch was cut against?
# --------------------------------------------------------------------------
# The tree sits inside EXHALE's working tree, so its history counts only if
# its own top level is this directory.  A tree without a history of its own
# is decided from the files in step 3.
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
# `git apply` resolves the patch's paths against the top of the repository it
# finds; stop the search at the parent so a tree without its own history is
# not resolved against EXHALE's.
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
# 4. verify each correction in the file it lives in
# --------------------------------------------------------------------------
fail=0
need() {   # need <file> <string> <what it is>
    if grep -qF -- "$2" "$1" 2>/dev/null; then
        echo "  ok      $3"
    else
        echo "  MISSING $3   ($1)"; fail=1
    fi
}
needcount() {   # needcount <file> <string> <count> <what it is>
    local n
    n=$(grep -cF -- "$2" "$1" 2>/dev/null || true)
    if [ "${n:-0}" -eq "$3" ]; then
        echo "  ok      $4"
    else
        echo "  MISSING $4   ($1: found ${n:-0}, expected $3)"; fail=1
    fi
}

echo "[setup_photochem] verifying:"
e=src/equilibrate/equilibrate_cea.f90
needcount "$e" "mval_mass_good = abs(b_0(i_atom))*self%mass_tol" 2 \
    "Equilibrate: element-relative mass-balance test in both update paths"
needcount "$e" "mval_mass_good = MAXVAL(b_0)*self%mass_tol" 0 \
    "Equilibrate: no threshold scaled by the most abundant element"
needcount "$e" "! converged is intent(out): undefined on entry whatever the caller" 2 \
    "Equilibrate: convergence flag set on entry of both abundance updates"
c=src/clima/clima_adiabat.f90
need "$c" "Background-pressure bisection did not converge." \
    "Clima make_profile_bg_gas: bracketing scan and bisection"
need "$c" "      if (previous_valid) then" \
    "Clima make_profile_bg_gas: nested bracket test"
need "$c" "real(dp), parameter :: epsfcn = 1.0e-8_dp" \
    "Clima simple_solver: hybrd at a set difference step"
need "$c" "mv%x(1) = log10(T_guess_ - self%T_trop)" \
    "Clima simple_solver: unknown is the contrast T_surf - T_trop"
need "$c" "fvec_(1) = (ISR*rad_enhancement - OLR + self%surface_heat_flow)/flux_scale" \
    "Clima simple_solver: scaled flux residual"
need "$c" "Temperature solve returned residuals above the normalized tolerance." \
    "Clima simple_solver: returned point checked"
g=photochem/extensions/gasgiants.py
need "$g" "equilibrium_mass_tol must be finite and positive" \
    "gasgiants: equilibrium_mass_tol argument, validated"
need "$g" "self.gas.mass_tol = equilibrium_mass_tol" \
    "gasgiants: tolerance set on the equilibrium solver"
need "$g" "closure_rtol = max(1.0e-8, 100.0*gas.mass_tol)" \
    "gasgiants: elemental closure checked after each solve"
need "$g" "molfracs_atoms_condensate" \
    "gasgiants: closure check skipped where a condensate is present"
t=tests/photochem/test_python.py
for name in test_gas_giant_equilibrium_mass_tolerance_validation \
            test_equilibrate_trace_element_closure \
            test_gas_giant_composition_closure_regimes \
            test_equilibrate_hot_high_pressure_single_solve_closes_elements \
            test_equilibrate_successive_solves_independent_of_name_reads; do
    need "$t" "def $name(" "test present: $name"
done
[ $fail -eq 0 ] || die "the tree is not in the expected state (above)"

say "$DEST is the source EXHALE runs on"
say "next: build and install it (requirements in the header of this script)."
say "      This script does not build."
