# `hp_zero_seed` -- the transported proton started from no protons at all

## What this case is

`hp_front` (read its README first: it carries the configuration, the reason
the initial condition is a restart, and where the IC pair lives) with one
change: the H II column of `IC/Ion_species_IC.txt` is zero in every row.

The H nuclei are conserved by the edit. The proton density of each row is
moved into H I, not deleted:

    HI_new  = HI_old + HII_old
    HII_new = 0

so H I + H II + 2(H2 + H2+) + 3 H3+ + HeH+ is bitwise what
`mol_base_handoff` wrote. That matters: `load_IC.f90` compares the loaded
He/H of every physical cell against the input's `He/H number ratio` and
rescales the whole column when they differ by more than 1e-6 -- and refuses
outright when the state carries HeH+, which this one does. An edit that
simply deleted the protons would change the hydrogen count and be rejected
before the seed reached the solver.

MEASURED 2026-09-05 on the written file: `x(H II) = 0` in every row, against
1.54e-8 (row 3, the base), 4.17e-4 (row 251) and 0.800 (row 502, the top) in
the file it was made from.

## What it is for

The transported proton has no local root to fall back on: with
`Ionization transport: True` the sweep is handed the carrier value and solves
the remaining rows against it. Starting from an identically zero proton asks
whether the transport-chemistry solve can build the ionization front out of
nothing -- whether the source term alone populates H II, and at what rate --
and it is the state in which a division by the proton density, or a rate
normalized by it, is reached with the numerator zero as well.

## Step count and what was measured

`maxsteps` = 100. MEASURED 2026-09-05, single-threaded, gfortran 16.2 build
`e779298d0b482bfbc40ae2f9aa80d60a`: 5.8 s wall, `final: count=100
du= 1.2819E+00`, Mdot 10.51 (log10 g/s), no NaN and no abort. The 100-step
state differs from `hp_front`'s and from `hp_trace_seed`'s (worst relative
difference 1.0 against each, on a column that is zero in one and not the
other), so the seed reaches the solver rather than being replaced on the way
in.

The case is a relaxation snapshot. It asserts that the run reproduces its
recorded output bitwise, nothing more; the derivative and front behavior the
plan asks about (rev 3 section 3.3, "derivative test and front test
separate") needs the dump hook named in `hp_front/README.md` and is not
tested here.

## Provenance

Built 2026-09-05 for Phase 0 item 4 of
`docs/development_plan_20260905_execution.md`. Not in `DEFAULT_CASES`;
reference output in `baseline_post170_20260905/hp_zero_seed/`.
