#!/bin/bash
# WHICH BINARY A SUITE RAN IS STATED BY THE SUITE, AND THE STATEMENT IS THE
# ONE THAT WAS ASKED FOR (PLAN_20260919_rev1 item P5b; PLAN_20260919_review
# 7.3).
#
# WHY A SENTINEL.  A suite that reports nothing about its binary can run a
# build nobody asked for and still print passing numerical rows, because a
# number is evidence about the physics and not about the file that produced
# it.  The rows below point EXHALE_EXE at a COPY of the tree's binary whose
# md5 differs from the tree's, and require each suite to report that copy.
# A suite that ignores the request reports the tree's md5 and fails here.
# The copy is made by appending bytes past the end of the ELF image, which
# changes the md5 and leaves the program itself untouched; no row here
# executes it, because each suite is asked for its identity alone.
#
# EXHALE_EXE_IDENTITY_ONLY=1 is the contract that makes this cheap: a suite
# prints its identity block and its identity row and stops without running
# anything.  The policy itself, including the accepted aliases and the
# refusal of a conflicting pair, is stated in src/tests/exhale_exe.sh.
#
# Every assertion prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; the exit status is nonzero if any of them fails.
#
# Usage: src/tests/executable_identity/run.sh
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
TESTS="$(cd "$HERE/.." && pwd)"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/executable_identity}"
mkdir -p "$OUT"
n_fail=0

row() {   # row <PASS|FAIL> <name> <measured> <reference>
   echo "$1 $2 measured=$3 reference=$4 tol=0"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

TREE="$ROOT/EXHALE.x"
if [ ! -x "$TREE" ]; then
   row FAIL executable_identity_tree_binary missing "$TREE"
   exit 1
fi
TREE_MD5="$(md5sum "$TREE" | awk '{print $1}')"

SENT="$OUT/EXHALE_sentinel.x"
\cp -f "$TREE" "$SENT"
printf 'sentinel copy of the tree binary, src/tests/executable_identity\n' >> "$SENT"
chmod +x "$SENT"
SENT_MD5="$(md5sum "$SENT" | awk '{print $1}')"
echo "  tree binary $TREE md5 $TREE_MD5"
echo "  sentinel copy $SENT md5 $SENT_MD5"
if [ "$SENT_MD5" != "$TREE_MD5" ]; then
   row PASS the_sentinel_differs_from_the_tree_binary "$SENT_MD5" "not_$TREE_MD5"
else
   row FAIL the_sentinel_differs_from_the_tree_binary "$SENT_MD5" "not_$TREE_MD5"
   exit 1
fi

# The suites that run the production binary.  A suite added later that runs
# the binary belongs in this list.
SUITES="adv_static_limit boundary_state charge_exchange_rows
        coupled_block_jacobian coupled_source_step fuv_band_ledger
        grid_and_gates molecular_seed residual_determinism run_mode
        species_face_flux spectrum_type steady_species_rows"

# reported_md5 <log-file>: the md5 the suite printed for its own binary.
reported_md5() {
   sed -n 's/^  binary md5: \(.*\)$/\1/p' "$1" | head -n 1
}

echo ""
echo "---- the identity a suite reports is the binary that was requested ----"
for s in $SUITES; do
   log="$OUT/$s.log"
   env -u EXHALE_RESID_EXE -u EXHALE_SPECIES_EXE -u EXHALE_COUPLED_EXE \
       EXHALE_EXE="$SENT" EXHALE_EXE_IDENTITY_ONLY=1 \
       EXHALE_TEST_OUT="$OUT/$s" \
       bash "$TESTS/$s/run.sh" > "$log" 2>&1
   m="$(reported_md5 "$log")"
   [ -z "$m" ] && m=not_reported
   if [ "$m" = "$SENT_MD5" ]; then
      row PASS "the_requested_binary_is_reported_$s" "$m" "$SENT_MD5"
   else
      row FAIL "the_requested_binary_is_reported_$s" "$m" "$SENT_MD5"
      echo "     see $log"
   fi
done

echo ""
echo "---- the selection policy ----"

# With no name set the binary is the tree's own build.
log="$OUT/default.log"
env -u EXHALE_EXE -u EXHALE_RESID_EXE EXHALE_EXE_IDENTITY_ONLY=1 \
    EXHALE_TEST_OUT="$OUT/default" bash "$TESTS/run_mode/run.sh" > "$log" 2>&1
m="$(reported_md5 "$log")"; [ -z "$m" ] && m=not_reported
if [ "$m" = "$TREE_MD5" ]; then
   row PASS the_default_is_the_tree_binary "$m" "$TREE_MD5"
else
   row FAIL the_default_is_the_tree_binary "$m" "$TREE_MD5"
fi

# The alias alone is honored, and is honored by every suite and not only by
# the residual-determinism scripts it was written for.
log="$OUT/alias.log"
env -u EXHALE_EXE EXHALE_RESID_EXE="$SENT" EXHALE_EXE_IDENTITY_ONLY=1 \
    EXHALE_TEST_OUT="$OUT/alias" bash "$TESTS/run_mode/run.sh" > "$log" 2>&1
m="$(reported_md5 "$log")"; [ -z "$m" ] && m=not_reported
if [ "$m" = "$SENT_MD5" ]; then
   row PASS the_alias_alone_is_honored "$m" "$SENT_MD5"
else
   row FAIL the_alias_alone_is_honored "$m" "$SENT_MD5"
fi

# Two names that agree name one binary, and the run proceeds.
log="$OUT/agree.log"
env EXHALE_EXE="$SENT" EXHALE_RESID_EXE="$SENT" EXHALE_EXE_IDENTITY_ONLY=1 \
    EXHALE_TEST_OUT="$OUT/agree" bash "$TESTS/residual_determinism/run.sh" \
    > "$log" 2>&1
m="$(reported_md5 "$log")"; [ -z "$m" ] && m=not_reported
if [ "$m" = "$SENT_MD5" ]; then
   row PASS two_agreeing_names_are_one_request "$m" "$SENT_MD5"
else
   row FAIL two_agreeing_names_are_one_request "$m" "$SENT_MD5"
fi

# Two names that disagree are refused: nothing is chosen and nothing runs.
log="$OUT/conflict.log"
env EXHALE_EXE="$SENT" EXHALE_RESID_EXE="$TREE" EXHALE_EXE_IDENTITY_ONLY=1 \
    EXHALE_TEST_OUT="$OUT/conflict" bash "$TESTS/residual_determinism/run.sh" \
    > "$log" 2>&1
crc=$?
if [ $crc -ne 0 ] && grep -q 'binary_selection measured=conflicting_request' "$log" \
   && [ -z "$(reported_md5 "$log")" ]; then
   row PASS conflicting_names_are_refused "refused_rc_$crc" refused
else
   row FAIL conflicting_names_are_refused "rc_${crc}_$(reported_md5 "$log")" refused
   echo "     see $log"
fi

# The seed-independence script takes its binary as an argument, which is a
# request of the same kind: it is reported, and a name that disagrees with it
# is refused.
log="$OUT/argument.log"
env -u EXHALE_EXE -u EXHALE_RESID_EXE EXHALE_EXE_IDENTITY_ONLY=1 \
    bash "$TESTS/steady_selfconsistent_residual/run_seed_independence.sh" \
    "$SENT" > "$log" 2>&1
m="$(reported_md5 "$log")"; [ -z "$m" ] && m=not_reported
if [ "$m" = "$SENT_MD5" ]; then
   row PASS the_binary_named_on_the_command_line_is_reported "$m" "$SENT_MD5"
else
   row FAIL the_binary_named_on_the_command_line_is_reported "$m" "$SENT_MD5"
fi

log="$OUT/argument_conflict.log"
env -u EXHALE_RESID_EXE EXHALE_EXE="$TREE" EXHALE_EXE_IDENTITY_ONLY=1 \
    bash "$TESTS/steady_selfconsistent_residual/run_seed_independence.sh" \
    "$SENT" > "$log" 2>&1
crc=$?
if [ $crc -ne 0 ] && grep -q 'binary_selection measured=conflicting_request' "$log"; then
   row PASS a_name_against_the_argument_is_refused "refused_rc_$crc" refused
else
   row FAIL a_name_against_the_argument_is_refused "rc_$crc" refused
   echo "     see $log"
fi

echo ""
if [ $n_fail -gt 0 ]; then
   echo "executable_identity: $n_fail failing assertion(s)"
   exit 1
fi
echo "executable_identity: every suite reported the binary it was given"
exit 0
