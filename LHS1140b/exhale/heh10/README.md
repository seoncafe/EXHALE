# heh10 — cold start diverged; warm start planned

The cold start at 30 R_p broke away from the other five compositions: while
they descended to ||R|| = 0.17-0.45 by step ~180k, this case rose to
||R|| ~ 24 with a flux spread near 280 and stayed there. The record is
`run_cold_diverged.log` / `output_cold_diverged/`; it is not a solution.
(The 10 R_p attempt of every case is under `attempt_rmax10/`.)

Resolution (2026-08-24): seeded from `heh1`'s JFNK-converged solution --
`finish_case.sh heh10 heh1` -- with load_IC carrying the state onto
He/H = 10. JFNK converged at once (info=0, ||R|| = 5.3e-4,
log10 Mdot = 7.76, continuous with the neighboring compositions). The IC is
the only thing that differs from the cold-started cases; no filter, no
changed thresholds.
