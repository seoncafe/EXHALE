# Rules for every worker (the common part of every brief)

Moved into the tree from the session scratchpad on 2026-09-09 (PLAN_20260909_rev1 item N-1); the brief names the plan and the item. Read first: `~/.claude/CLAUDE.md`, `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/CLAUDE.md`, then the plan the brief names (currently `docs/PLAN_20260909_rev1.md`) and the item's lines in the review that motivated it. Repository `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/`.

Rules that apply to every worker:
- Never run `git add`, `git commit`, `git push`, `git mv`, or anything else that touches the git index.
- Forbidden wording anywhere (code, comments, docs, report): vendor, vendored, bundled, machinery, em-dashes, "per-X + noun" phrases.
- Comments state the physics or the invariant and the source; never the history of a change.
- Label every number MEASURED (you ran it) or READ (from source or a document).
- Concurrency: other workers edit other files in the same tree. Never run plain `make`, `make check`, `make test` or `backup/regression/run_check.sh`. Build privately with `make OBJDIR=build_<item> EXE=EXHALE_<item>.x`; run only on scratch copies under the scratch directory the brief names (`<scratchpad>/<item>/`); delete the private build at the end. Test suites under `src/tests/<suite>/run.sh` build into `build/tests/<suite>/`; when a suite supports a private object directory (see `src/tests/spectrum_type/run.sh`) use it. A compile error in a file you do not own is another worker mid-edit: wait a minute and retry, three times, then report and stop.
- Never kill a process by name (`pkill`, `killall`): other workers run binaries with the same basenames. Kill by PID only after checking `/proc/<pid>/cwd` is YOUR scratch directory.
- You edit ONLY the files your item names. If the fix needs another file, stop and report.
- Every change comes with its test: RED before, GREEN after, both shown in the report. No tolerance is chosen to make an existing snapshot pass.
- Gate for a change that touches a shared acceptance interface: a scoped impact measurement, not byte identity as a criterion. Run the two cases the item names on scratch copies with the binary before and after, report whether outputs are byte-identical, and if not where and why (which branch fired). Never refresh a golden.
- Report to the file the brief names (`<scratchpad>/<plan>_<item>_report.md`) and return its text: verdict first, what changed (diff summary), tests before/after, impact measurement, anything noticed outside scope (report, do not fix).
- One owner per source file per increment: `OBJDIR=build_<item>` isolates objects, not the source text; two workers on one file is a brief error, stop and report.
- Validation scope is the affected path the brief lists (PLAN_20260909_rev1 section 7): do not run the regression matrix or unrelated suites; do not rerun a check after a change that cannot move it (documentation, a file move, a message text).
- Report by named outcomes and measured numbers, not by iteration counts.
- Copy and move with `\cp -f` and `\mv -f` (the interactive aliases park a script on an `overwrite?` prompt for hours with later statements still pending; killing the `cp` then runs them). Never leave a background shell of your own behind when you finish; list your PIDs at the end of the report.
- Build your own control binary for a before/after comparison (the delivered objects with only your file's entry-text object swapped in, or a private build of the entry text). The advisor's `EXHALE_lwv.x` is rebuilt at every verification and is a build of the live tree, never an entry text.
- Run every suite whose driver links the module you changed (grep the suites' `run.sh` for the object), not only the suite named in the brief.
