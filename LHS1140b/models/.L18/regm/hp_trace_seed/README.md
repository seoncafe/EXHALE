# `hp_trace_seed` -- the transported proton started from a trace

## What this case is

`hp_front` (read its README first: it carries the configuration, the reason
the initial condition is a restart, and where the IC pair lives) with one
change: the H II column of `IC/Ion_species_IC.txt` is set to 1e-12 of the H
nuclei of its own row,

    HII_new = 1e-12 * n_H,   n_H = HI + HII + 2(H2 + H2+) + 3 H3+ + HeH+
    HI_new  = HI_old + HII_old - HII_new

so the hydrogen count of every row is unchanged, for the reason given in
`hp_zero_seed/README.md` (`load_IC.f90` rescales or refuses a column whose
He/H has moved). Where the row's own H I + H II is smaller than the trace the
value is clamped to it, so hydrogen is never taken from a molecular carrier;
MEASURED, no row of this file reached that clamp.

MEASURED 2026-09-05 on the written file: `x(H II) = 1.000e-12` in every row,
against 1.54e-8 (row 3, the base), 4.17e-4 (row 251) and 0.800 (row 502, the
top) in the file it was made from.

## What it is for, beside `hp_zero_seed`

Zero and 1e-12 are physically the same state and numerically are not. The
pair separates a solve that cannot start from an identically zero proton
(a rate divided by n(H+), a logarithm, a relative convergence measure) from
one that genuinely fails to build the front: a difference between the two
cases at 100 steps that is larger than the difference between 1e-12 and 0 of
the H nuclei is a property of the method, not of the initial condition.

## Step count and what was measured

`maxsteps` = 100. MEASURED 2026-09-05, single-threaded, gfortran 16.2 build
`e779298d0b482bfbc40ae2f9aa80d60a`: 5.8 s wall, `final: count=100
du= 1.2819E+00`, Mdot 10.51 (log10 g/s), no NaN and no abort. `du` agrees
with `hp_zero_seed` to the five printed digits and differs from `hp_front`
(1.2828), which is what a seed 1e-12 below the front should do. The two
states are nevertheless not identical: worst relative difference 1.0 against
`hp_zero_seed` on a column that is zero in one and 1.07e-9 in the other.

## Provenance

Built 2026-09-05 for Phase 0 item 4 of
`docs/development_plan_20260905_execution.md`. Not in `DEFAULT_CASES`;
reference output in `baseline_post170_20260905/hp_trace_seed/`.
