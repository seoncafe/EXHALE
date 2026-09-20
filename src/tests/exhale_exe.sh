# THE BINARY A TEST RUNS IS SELECTED IN ONE WAY, AND THE TEST SAYS WHICH
# BINARY IT RAN.  This file is the one place the policy is stated; every
# suite under src/tests/ that runs the production binary sources it and the
# other scripts point here rather than restating it.
#
# THE POLICY.
#
#   EXHALE_EXE is the way to select the binary.  Unset, the binary is
#   $ROOT/EXHALE.x, the tree's own build.
#
#   EXHALE_RESID_EXE is an accepted alias of EXHALE_EXE, kept because the
#   residual-determinism scripts were written with it.  Two suites carry a
#   further alias of their own, EXHALE_SPECIES_EXE (species_face_flux) and
#   EXHALE_COUPLED_EXE (coupled_source_step), which also switch their
#   optional whole-binary rows on.
#
#   PRECEDENCE.  When more than one of these names is set and they resolve
#   to the SAME file, EXHALE_EXE is the name reported and the run proceeds.
#   When they resolve to DIFFERENT files the suite REFUSES: it prints one
#   FAIL row naming both settings and exits 1.  A run that asked for two
#   binaries never gets one of them chosen silently.
#
#   IDENTITY.  Before anything is run, the suite prints the path of the
#   binary and its md5, and asserts that the file it is about to run has the
#   md5 that was requested.  A suite that copies the binary into a working
#   directory asserts the identity of the COPY, which is the file that
#   executes.  A numerical row that passes is not evidence of which build
#   produced it; the identity row is.
#
#   EXHALE_EXE_IDENTITY_ONLY=1 makes a suite print the identity block and
#   the identity row and stop without running anything.  It is what the
#   executable_identity suite reads to check, cheaply and on every suite,
#   that the binary a suite reports is the binary that was requested.
#
# USAGE, in a suite that runs the binary:
#
#   . "$HERE/../exhale_exe.sh"
#   exhale_select_exe "$ROOT" <suite-name> [extra-alias-name]
#   EXE="$EXHALE_RUN_EXE"
#   exhale_announce_exe
#
# and, where the binary is copied before it is run:
#
#   cp -f "$EXE" "$WORK/EXHALE.x"
#   exhale_assert_exe_identity "$WORK/EXHALE.x" || exit 1
#
# After exhale_select_exe the caller can read:
#   EXHALE_RUN_EXE            the selected path
#   EXHALE_RUN_EXE_MD5        its md5, or "missing"
#   EXHALE_RUN_EXE_SOURCE     which name selected it, in words
#   EXHALE_RUN_EXE_REQUESTED  1 if any selecting name was set, else 0

exhale_md5() {
   md5sum "$1" 2>/dev/null | awk '{print $1}'
}

exhale_realpath() {
   readlink -f "$1" 2>/dev/null || echo "$1"
}

# exhale_select_exe <root> <suite-name> [extra-alias-name]
exhale_select_exe() {
   local root="$1"
   EXHALE_SUITE_NAME="$2"
   local extra_name="${3:-}"
   local names="" chosen="" chosen_name=""
   local n v r first_r=""
   for n in EXHALE_EXE EXHALE_RESID_EXE "$extra_name"; do
      [ -z "$n" ] && continue
      eval "v=\"\${$n:-}\""
      [ -z "$v" ] && continue
      r="$(exhale_realpath "$v")"
      if [ -z "$chosen" ]; then
         chosen="$v"; chosen_name="$n"; first_r="$r"
      elif [ "$r" != "$first_r" ]; then
         echo "FAIL ${EXHALE_SUITE_NAME}_binary_selection measured=conflicting_request reference=one_binary tol=0"
         echo "     $chosen_name=$chosen"
         echo "     $n=$v"
         echo "     the two name different files, so nothing is chosen; set one"
         echo "     of them.  The policy is stated in src/tests/exhale_exe.sh"
         exit 1
      fi
      names="$names $n"
   done

   if [ -z "$chosen" ]; then
      EXHALE_RUN_EXE="$root/EXHALE.x"
      EXHALE_RUN_EXE_SOURCE="the default, the tree's own EXHALE.x, no name set"
      EXHALE_RUN_EXE_REQUESTED=0
   else
      EXHALE_RUN_EXE="$chosen"
      EXHALE_RUN_EXE_SOURCE="$(echo "$names" | sed 's/^ //; s/ /, /g'), naming one file"
      EXHALE_RUN_EXE_REQUESTED=1
   fi
   EXHALE_RUN_EXE_MD5="$(exhale_md5 "$EXHALE_RUN_EXE")"
   [ -z "$EXHALE_RUN_EXE_MD5" ] && EXHALE_RUN_EXE_MD5="missing"
   # A subtest, in the shell or in Python, reads the selection rather than
   # repeating it.
   export EXHALE_RUN_EXE EXHALE_RUN_EXE_MD5 EXHALE_RUN_EXE_SOURCE \
          EXHALE_RUN_EXE_REQUESTED
   return 0
}

# exhale_use_exe <path> <suite-name>
# The binary is named on the command line, which is a request like any
# other: a selecting name that resolves to a different file is refused.
exhale_use_exe() {
   local given="$1"
   exhale_select_exe "$(cd "$(dirname "$given")" 2>/dev/null && pwd || echo .)" "$2"
   if [ "$EXHALE_RUN_EXE_REQUESTED" = "1" ] && \
      [ "$(exhale_realpath "$EXHALE_RUN_EXE")" != "$(exhale_realpath "$given")" ]; then
      echo "FAIL ${EXHALE_SUITE_NAME}_binary_selection measured=conflicting_request reference=one_binary tol=0"
      echo "     the binary named on the command line is $given"
      echo "     $EXHALE_RUN_EXE_SOURCE gives $EXHALE_RUN_EXE"
      echo "     the policy is stated in src/tests/exhale_exe.sh"
      exit 1
   fi
   EXHALE_RUN_EXE="$given"
   EXHALE_RUN_EXE_REQUESTED=1
   EXHALE_RUN_EXE_SOURCE="the binary named on the command line"
   EXHALE_RUN_EXE_MD5="$(exhale_md5 "$given")"
   [ -z "$EXHALE_RUN_EXE_MD5" ] && EXHALE_RUN_EXE_MD5="missing"
   return 0
}

# exhale_assert_exe_identity <path-of-the-file-that-runs>
exhale_assert_exe_identity() {
   local ran="$1" m
   m="$(exhale_md5 "$ran")"
   [ -z "$m" ] && m="missing"
   if [ "$m" != "missing" ] && [ "$m" = "$EXHALE_RUN_EXE_MD5" ]; then
      echo "PASS ${EXHALE_SUITE_NAME}_binary_identity measured=$m reference=$EXHALE_RUN_EXE_MD5 tol=0"
      return 0
   fi
   echo "FAIL ${EXHALE_SUITE_NAME}_binary_identity measured=$m reference=$EXHALE_RUN_EXE_MD5 tol=0"
   echo "     the file about to run is $ran"
   echo "     the requested binary is $EXHALE_RUN_EXE"
   return 1
}

exhale_announce_exe() {
   if [ "${EXHALE_EXE_IDENTITY_STATED:-0}" = "1" ]; then
      # The calling suite already stated and asserted the identity of this
      # same file, so a subtest does not repeat the row.
      return 0
   fi
   echo "  binary: $EXHALE_RUN_EXE"
   echo "  binary md5: $EXHALE_RUN_EXE_MD5"
   echo "  binary selected by: $EXHALE_RUN_EXE_SOURCE"
   exhale_assert_exe_identity "$EXHALE_RUN_EXE" || exit 1
   EXHALE_EXE_IDENTITY_STATED=1
   export EXHALE_EXE_IDENTITY_STATED
   if [ "${EXHALE_EXE_IDENTITY_ONLY:-0}" = "1" ]; then
      echo "  EXHALE_EXE_IDENTITY_ONLY=1: the identity is stated and nothing is run"
      exit 0
   fi
   return 0
}
