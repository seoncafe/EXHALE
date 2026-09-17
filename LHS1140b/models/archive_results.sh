#!/bin/bash
# Move the results of the tree aside, keeping the inputs, so that the catalog
# can be solved again on a new binary while the old answers stay readable.
#
#     models/archive_results.sh <tag> [group-name-substring ...]
#
# writes LHS1140b/models_<tag>/, one directory per case in the same relative
# layout, holding everything the runs produced -- output/, the closure
# iterations k00/ ... , seed/, attempt_*/, the logs, the tpm_* line products
# and REPRODUCE.md -- plus a COPY of the input files, so a preserved case is
# self-contained. What stays in models/ is the inputs alone:
#
#     input.inp  base.inp  lower_atmosphere_profile.dat
#     closure.json  input_template.inp
#
# and `python3 make_models.py` then reports them unchanged. `pick_seed.py`
# reads the newest preserved tree as its tier 0 where `models/` no longer
# carries a case's own certified state, so name the new tree in its
# PRESERVED constant after running this.
#
# ARCHIVE_EXCLUDE, an extended regular expression on "<group>/HeH<value>",
# leaves cases where they are (the 2026-09-15 re-run used
# ARCHIVE_EXCLUDE='^molecular_' to keep the molecular group, which was not
# part of it, in place).
#
# A case with nothing to move is passed over. Nothing is deleted and nothing
# is overwritten: the tag directory must not already hold the case.

set -u
MODELS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LHS=$(dirname "$MODELS")

if [ $# -lt 1 ]; then
   echo "usage: $(basename "$0") <tag> [group-name-substring ...]" >&2
   exit 1
fi
TAG=$1; shift
DEST="$LHS/models_$TAG"
: "${ARCHIVE_EXCLUDE:=}"
KEEP="input.inp base.inp lower_atmosphere_profile.dat closure.json input_template.inp"

mkdir -p "$DEST" || exit 1

CASES=$(cd "$MODELS" && python3 -c '
import os, sys
sys.path.insert(0, ".")
from make_models import GROUPS
want = sys.argv[1:]
for name, ladder in GROUPS:
    if want and not any(w in name for w in want):
        continue
    for heh in ladder:
        print("%s/HeH%s" % (name, heh))
' "$@")

moved=0; passed=0; excluded=0
for c in $CASES; do
   if [ -n "$ARCHIVE_EXCLUDE" ] && echo "$c" | grep -qE "$ARCHIVE_EXCLUDE"; then
      echo "left in place (ARCHIVE_EXCLUDE): $c"; excluded=$((excluded + 1)); continue
   fi
   src="$MODELS/$c"
   [ -d "$src" ] || continue
   if [ -d "$DEST/$c" ]; then
      echo "REFUSED: $DEST/$c already exists"; exit 2
   fi
   n=0
   for e in "$src"/* "$src"/.[!.]*; do
      [ -e "$e" ] || continue
      b=$(basename "$e")
      keepit=0
      for k in $KEEP; do [ "$b" = "$k" ] && keepit=1; done
      [ $keepit = 1 ] && continue
      mkdir -p "$DEST/$c"
      \mv -f "$e" "$DEST/$c/" || exit 1
      n=$((n + 1))
   done
   if [ $n -gt 0 ]; then
      for k in $KEEP; do [ -f "$src/$k" ] && \cp -f "$src/$k" "$DEST/$c/"; done
      moved=$((moved + 1)); echo "preserved $c ($n entries)"
   else
      passed=$((passed + 1)); echo "nothing to preserve: $c"
   fi
done

\cp -f "$LHS/MODELS.md" "$DEST/MODELS.md"
\cp -f "$MODELS/README.md" "$DEST/README_models.md"

# The binaries the preserved runs were solved with, read from their own
# records, so the tree states what produced it without being told.
MD5S=$(grep -h '^| md5 |' "$DEST"/*/*/REPRODUCE.md 2>/dev/null \
       | sed 's/.*`\(.*\)`.*/\1/' | cut -c1-12 | sort | uniq -c \
       | awk '{printf "%s (%s case records)\n", $2, $1}')

cat > "$DEST/README.md" <<EOF
# The LHS 1140 b catalog, results preserved as \`$TAG\`

Moved here by \`models/archive_results.sh $TAG\` on $(date '+%Y-%m-%d'),
before the catalog was solved again. These are the results as they stood;
the inputs beside them are copies of what \`models/\` still holds, so a case
here is self-contained. It is a record and not a runnable tree: the runners
of \`models/\` act on \`models/\`.

$moved case directories preserved, $passed with nothing to preserve, $excluded
left in place.

Binaries these runs were solved with, read from each case's \`REPRODUCE.md\`:

\`\`\`
$MD5S
\`\`\`

\`MODELS.md\` and \`README_models.md\` are the catalog and the tree-level
reproduction notes as they stood. Why the catalog was re-solved is stated in
\`MODELS.md\` section 5 and in the plan item that forced it.
EOF

echo "--- $moved preserved, $passed with nothing to preserve, $excluded left in place -> $DEST ---"
