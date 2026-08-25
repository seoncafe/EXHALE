# He/H series on the tutorial planet — Phase C first attempt

This directory holds the original Phase C attempt: the He-rich audit run on
the tutorial planet, with the He/H number ratio the only parameter changed
(`heh_1`, `heh_10`, `heh_100`, `heh_1000`, driven by `run_case.sh`).

`heh_1`, `heh_10` and `heh_1000` completed. `heh_100` did not: none of the
three marching variants tried — `heh_100`, `heh_100_plm002` (looser PLM
threshold) and `heh_100_newton_early` (Newton handed the marched state
earlier) — reached a state the audit could use, and each was terminated
rather than left running (`case_heh_100*.log`).

Phase C was therefore rescoped on 2026-08-24 onto the already converged
LHS 1140 b solutions in `../exhale/`, which are the target planet and span a
wider He/H range (0.55 to 10000) on both stellar spectra. `audit_checks.py`
now takes case directories as arguments, so it runs against those solutions
unchanged; with no arguments it still audits the four cases here.

The Phase C acceptance record is `../exhale/audit_summary.md`.
