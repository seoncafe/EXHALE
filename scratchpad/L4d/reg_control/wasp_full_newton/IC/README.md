# Pinned certified state of `wasp_full_newton`

Copied 2026-09-09 from `backup/regression/golden/wasp_full_newton/` (the
certified Newton state of the 2026-09-08 19:05 golden refresh, `certified=T`,
`sec_ion_step=2402`), renamed for `Load IC? True`. The suite row
`restart_intent` of `src/tests/grid_and_gates/` reads THIS pair, not the live
`output/` of the case, which the regression matrix rewrites (and leaves as a
relaxation snapshot while a matrix run is in progress). Refresh this pair only
together with the golden, and say so here.
