#!/bin/bash
# WHICH STATE OF A RUN A TRANSIT SPECTRUM STANDS ON.
#
# A run can hold two descriptions of one column.  The solved pair,
# `Hydro_ioniz.txt` and `Ion_species.txt`, is the state the wind solver
# produced and the certification judged; with the ionization stages
# transported its composition is the partition the stage transport produced.
# The `_adv` pair is a post-process of a marching state, a steady advective
# correction at the bulk velocity, an independent discretization of the same
# column.  A stationary solve writes no `_adv` products, so a solved run has
# no `adv` state, and one is not to be manufactured for it.
#
# These rows pin the selection (EXHALE_TRANSIT_STATE), the refusal of a
# missing `_adv` pair by name instead of a silent substitution, the solution
# selection on a run that has no post-process, the identity of the files a
# curve was built from in its own header, and the reported distance between
# the two states against a direct comparison of the same columns.
#
# Python only: no Fortran object or binary is needed.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; lines beginning with two spaces or with DIAGNOSTIC are context.  Exit
# status is nonzero if any row fails.
#
# Usage: src/tests/transit_state/run.sh
#   TRANSIT_STATE_MODULE_DIR  directory to import exhale_transit_lib from
#                             (default the repository root); point it at
#                             another copy of the module to measure that copy
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
MODDIR="${TRANSIT_STATE_MODULE_DIR:-$ROOT}"

FAIL=0
for SCRIPT in state_selection_and_provenance.py ; do
   env PYTHONPATH="$MODDIR:$ROOT/examples${PYTHONPATH:+:$PYTHONPATH}" \
       MPLBACKEND=Agg \
       python3 "$HERE/$SCRIPT"
   [ $? -ne 0 ] && FAIL=1
done

if [ "$FAIL" -ne 0 ]; then
   echo "transit_state: FAILED"
else
   echo "transit_state: PASSED"
fi
exit "$FAIL"
