#!/bin/bash
# THE CHORD GEOMETRY THE TRANSIT SPECTRA AND THEIR CENSUS SHARE.
#
# A ray at impact parameter b crosses every shell of radius r >= b, at
# line-of-sight coordinate |x| = sqrt(r^2 - b^2), with no upper radial limit:
# a shell outside the stellar disk still absorbs on the rays that fall on the
# disk.  One selection (exhale_transit_lib.chord_shell_indices) therefore has
# to serve both the optical-depth integrals and any census of the material
# they sample, or the diagnostic describes a different atmosphere from the
# spectrum.  These rows pin that selection on the geometry of the review
# (shell r = 12, ray b = 2, ray grid capped at Rib = 10) and pin the
# line-center optical-depth share the census reports for the rows whose steady
# advective correction was refused, including the additivity that makes the
# share exact.
#
# Python only: no Fortran object or binary is needed.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; lines beginning with two spaces or with DIAGNOSTIC are context.  Exit
# status is nonzero if any row fails.
#
# The second script pins what a spectrum knows about the validity of the rows
# it stands on: the three `_adv` schemas (two status fields, one field, and a
# profile that states none, whose rows are UNKNOWN and not corrected), the
# measure each row was corrected under and the fraction a corrected row is
# accurate to in it, and the comment metadata every saved curve carries, which
# must leave the numerical columns of the file untouched.  A spectrum built on
# corrected rows is accurate to no better than the fraction those rows carry,
# so that fraction has to travel out of the profile header and into the curve.
#
# Usage: src/tests/transit_census/run.sh
#   TRANSIT_CENSUS_MODULE_DIR  directory to import exhale_transit_lib from
#                              (default the repository root); point it at
#                              another copy of the module to measure that copy
#   TRANSIT_CENSUS_ROOT        directory the saved spectral files are read
#                              under (default the repository root)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
MODDIR="${TRANSIT_CENSUS_MODULE_DIR:-$ROOT}"

FAIL=0
for SCRIPT in chord_shell_geometry.py adv_validity_and_metadata.py ; do
   env PYTHONPATH="$MODDIR:$ROOT/examples${PYTHONPATH:+:$PYTHONPATH}" \
       TRANSIT_CENSUS_MODULE_DIR="$MODDIR" \
       TRANSIT_CENSUS_ROOT="${TRANSIT_CENSUS_ROOT:-$ROOT}" \
       MPLBACKEND=Agg \
       python3 "$HERE/$SCRIPT"
   [ $? -ne 0 ] && FAIL=1
done

if [ "$FAIL" -ne 0 ]; then
   echo "transit_census: FAILED"
else
   echo "transit_census: PASSED"
fi
exit "$FAIL"
