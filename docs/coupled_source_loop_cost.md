# The coupled source loop: what it costs and why

The measured record behind the constants of the coupled source step
(`csm_T_tol`, `csm_comp_tol`, `csm_max_pass`, `csm_err_safety` and the
`csm_geom_*`, `csm_extrap_*` and `csm_x_*` constants in
`src/EXHALE_main.f90`). The code site keeps the
conclusions; this memo keeps the tables and the history they were drawn
from. Items: B3c (the loop), COST2 (Aitken taken cell by cell, a coupled Newton
declined), COST3 (carried columns), FIELD-SELF (the cell's own optical
depth), COST4 (the ladder of modes), COST5 (the pass-count law), COST6 (the
global-mode extrapolation), COST7 (the stopping test on the estimated
error). Reports: `docs/Update_EXHALE.md`, the entries of those names. All numbers MEASURED unless marked otherwise; "p" is passes a
coupled step, wall times are one thread on the 72-core machine unless said.

## 1. What the tolerances cost (2026-09-06)

WHAT IT COSTS AND WHAT LOOSENING IT WOULD BUY, MEASURED 2026-09-06,
one private binary per value.  The state column is the largest
change of any column of Hydro_ioniz.txt or Ion_species.txt in units
of that column's own maximum, against the same run at 1e-10:
               wasp_full, 300 steps       mol_base_handoff, 1500
    1e-9    3.2e-10  12.28p  166.6 s    3.5e-09  8.73p  313.9 s
    1e-8    8.5e-08  10.10p  140.2 s    1.2e-07  7.19p  273.3 s
    1e-7    1.2e-07   9.10p  124.5 s    2.6e-07  5.92p  231.8 s
    1e-6    6.5e-07   7.36p  106.3 s    4.4e-07  4.72p  176.0 s
    1e-5    1.1e-05   5.53p   80.9 s    1.7e-06  3.45p  135.6 s
The pair contracts by 0.19 per pass (median over the 2729 passes of
the wasp_full run), so a decade of tolerance is 1.3 passes and about
15 per cent of the wall time.  THE TOLERANCE IS THEREFORE NOT THE
COST LEVER: even three decades of it, which the anchor above
forbids, would leave the coupling within a factor of two of what it
costs now.  What the cost is made of is stated at csm_max_pass.

The accuracy of the trajectory does not enter: the step's own
truncation error, from the step-doubling estimate of wasp_full in
phys mode, is 5.3e-3 relative, five decades above the loosest value
in the table.

## 2. Where the cost is, and carrying the column (COST2, COST3, FIELD-SELF)

WHERE THE COST IS, MEASURED 2026-09-07 with timers inside ioniz_eq
(wasp_full, 300 steps, one thread).  Of the marching wall time, the
attenuated photoionization field was 38 per cent, the cell chemistry
loop 25 per cent and the two traversals of the cooling 16 per cent;
the run's own timers put 84 per cent of it inside ioniz_eq and 16
per cent inside the energy solve.  The field was the largest single
term AND it set the pass count: it was built once per pass from the
composition the pass was HANDED, so information reached a cell from
its outer neighbour only on the next pass and the ionization front
walks inward one cell per pass (the cell carrying the largest
composition move goes 206, 196, 191, 187, 185, 183, 181, 180, 179,
178, 177 over the passes of one wasp_full step) although each cell
reaches its own root in one Newton step.

CARRYING THE SOLVED COLUMN FORWARD REMOVES THAT WALK AND DOES NOT
PAY.  The sweep can be given the columns of the composition it has
ALREADY solved instead (ionization_equilibrium,
xuv_field_block_cells); that cell then stops walking and stands at
the front (202, 198, 198, 197, 197, 197, 197, 197, 196, 196, 196 on
the same step), and the pass count falls from 10.10 to 7.13.
But a pass costs the same field work either way, and the field of a
whole grid is one loop over independent cells that a machine solves
all at once: giving that width up costs more at sixteen threads than
the pass count saves (wasp_full at 300 steps, 17.0 s whole-grid
against 42.5 s at the width with the best pass count).  The measured
table is at that constant's declaration, and the default is the
whole grid.

THE COUNT IS A LADDER OF NESTED MODES AND IT IS THE SLOWEST OF THEM,
NOT THEIR SUM (item COST4, MEASURED on wasp_full at 300 steps by
removing one composition coupling at a time): the production loop
takes 10.09 passes; with the stellar field's dependence on the
composition removed, 7.16; with the temperature removed as well (the
energy solve replaced by its own converged answer), 4.94; with the He
recombination rates removed too, 2.00.  Each coupling has its own
contraction rate and the loop stops when the slowest has decayed, so
removing a subordinate one changes the total by nothing.  That is why
every lever tried has failed to move the number: the He recombination
lag is the fastest mode (10.09 to 10.01, and removing the coupling
outright gives the same 10.01); the H(n=2) Balmer continuum is not a
lag at all, since excited_H_update fills its arrays ONCE before this
loop and they are constants of it (refreshing it every pass reads
10.10 for 1.5 per cent more wall time); and closing the stellar
field's own lag in full, a one-cell block and eight self-field
passes, lands on the next rung, 8.18 passes for 1.56 times the wall
time at one thread and, from COST3, 2.8 times at sixteen.  The
temperature is not the residual either: the composition sweep alone,
at the temperature the step CONVERGED to, still takes 10.94 passes.
To go below seven every rung has to come out at once.

WHAT THE COUNT IS, AS A FORMULA (item COST5, MEASURED and then
TESTED against a prediction): the loop's error decays geometrically
with a ratio stable from the second pass, so
  passes = 2 + log10(the step's own excursion / the stopping
                     tolerance) / (decades per pass) ,
and every term of it is measured.  On mol_base_handoff the
composition contracts by 0.171 a pass and the temperature by 0.072,
i.e. 0.77 decades a pass, from an excursion of 7.0e-5 to a tolerance
of 1e-8: 7.00 predicted against 6.93 measured.  On wasp_full the
ratios are 0.222 and 0.199 from an excursion seven times larger: 9.2
predicted against 10.09.  The law says the count depends on the
stopping tolerance through one logarithm and on nothing else, at
1.30 passes per decade, and moving both tolerances together gives
8.06, 6.96, 5.68, 4.04 at 1e-9, 1e-8, 1e-7, 1e-6, i.e. 1.34.  It is
not the cell solve either: six decades of hybrd1 tolerance leave the
count at 6.97.  So no coupling is slow and nothing is waiting; the
loop walks four to five decades at two thirds of a decade a pass.
The molecular path carries TWO modes of nearly equal cost, the
composition and the temperature, each about six passes on its own,
which is why removing either leaves the other standing and why its
ladder looks flat; with both out it reaches 2.01.  The one lag that
limits the composition mode's RATE is n_tot, the third body of the
three-body molecular reactions, built from the entry composition:
holding it takes the composition ratio from 0.171 to 0.073 and the
count from 6.93 to 6.11, the rest being taken up by the temperature.

## 3. The ladder of nested modes (COST4)

See section 2 above, the paragraph beginning "THE COUNT IS A LADDER"; the
full matrix is in the COST4 entry of `docs/Update_EXHALE.md`.

## 4. The pass-count law (COST5)

See section 2 above, the paragraph beginning "WHAT THE COUNT IS, AS A
FORMULA"; the tolerance scan (8.06 / 6.96 / 5.68 / 4.04 at 1e-9 / 1e-8 /
1e-7 / 1e-6) and the cell-tolerance scan (6.96 / 6.97 / 6.97 / 6.98 at
1.5e-8 / 1e-10 / 1e-12 / 1e-14) are in the COST5 entry of
`docs/Update_EXHALE.md`.

## 5. The extrapolation: why the object is the global mode (COST6)

THE OBJECT IS THE WHOLE STATE'S DOMINANT MODE, NOT A CELL'S RATIO,
and that distinction is the whole of it.  An earlier attempt (B3c,
reproduced by item COST2) extrapolated with a ratio taken component
by component and MEASURED IT WORSE, 12.00 passes against 11.05 on
wasp_full.  The reason was then measured too: within one step the
contraction ratio spreads by 122 per cent over the global sequence's
tail and by 200 per cent from cell to cell (10th to 90th percentile
0.070 to 0.260), because the ionization front walks through a cell
while the iteration runs, and the amplification theta/(1 - theta) is
unbounded as theta approaches one.  What IS stable is the ratio of
the whole increment VECTOR, a Rayleigh quotient
    theta = <d_k, d_{k-1}> / <d_{k-1}, d_{k-1}> ,
i.e. the projection of the error onto the iteration's dominant decay
mode; that quantity was measured at 0.072 (temperature) and 0.171
(composition) on mol_base_handoff and 0.199 and 0.222 on wasp_full,
stable from the second pass onward.  One such scalar is taken for
the temperature block and one for the composition block, in the
norms the stopping test itself uses (relative for the temperature,
absolute mass fraction for the composition), and the candidate is a
single point on the single ray they define.

## 6. The reach guard (COST6)

WHAT THE FACTOR IS WORTH, MEASURED on mol_base_handoff at 300 steps
against a control of 6.96 passes per coupled step:
    reach     1      2      3     10     30    100
    passes  6.61   6.15   6.36   7.81   7.74   7.75
Above 3 the model is asked for more than it has and the loop pays a
pass for each candidate its own next pass then undoes, which is the
same statement as the 200 per cent spread of the cell ratios that
defeated the earlier attempt, measured from the other side.

## 7. What the loop's stopping test measures (COST6, closed by COST7)

Until item COST7 the test stopped on the last INCREMENT, while the error of
the accepted state is |theta/(1 - theta)| times it. Decision item 8 of
`docs/To_be_determined_by_user_20260906.md` was taken on 2026-09-08, option
(b): the test is now taken on the estimated error at the same two
tolerances. Section 8 is that item's measurement, including how much of the
one-to-two decades COST6 estimated survives the norm the estimate is taken
in.

## 8. The stopping test on the estimated error (COST7)

Decision item 8, option (b), implemented 2026-09-08: the loop stops when the
ESTIMATED REMAINING ERROR of the pair is within `csm_T_tol` and
`csm_comp_tol`, and falls back to the last increment on a pass whose
sequence gives no estimate. The estimate is
`csm_err_safety * |theta/(1 - theta)|` times the pass's own increment, with
theta the block Rayleigh quotient of section 5 and its four guards, formed
in one routine (`coupled_pair_geometric_error_estimate`) that both the
stopping test and the extrapolation read.

THE ESTIMATE IS TAKEN IN A DIFFERENT NORM FROM THE TEST, AND THAT IS THE
ITEM'S FINDING. theta is an L2 Rayleigh quotient; the two tolerances are
max-norm tests that one cell and one species gate, and that cell is not the
one the L2-dominant mode lives in. MEASURED with the probe
(`EXHALE_CSM_ERR_PROBE`, which holds the accepted state, walks the same
sequence on to an increment of 1e-12 and records the distance between the
two), the max-norm distance of a state accepted near the tolerance is up to
2.4 times the projected estimate in the composition norm, so an unfactored
estimate accepts states OUTSIDE that tolerance. (On a state accepted far
inside the tolerance the ratio reaches 12, which costs nothing: what the
factor has to cover is the states accepted at the tolerance.) The factor is
measured, on mol_base_handoff at 300 steps against a control of 6.96 passes:

| csm_err_safety | 1 | 2 | 3 | increment test |
|---|---|---|---|---|
| worst distT | 4.5e-09 | 2.5e-09 | 2.0e-09 | -- |
| worst distC | 1.8e-08 | 1.0e-08 | **7.4e-09** | -- |
| over 1e-08 | 100/299 | 1/299 | **0/299** | -- |
| passes | 6.14 | 6.76 | **6.82** | 6.96 |

and on wasp_full at 300 steps against a control of 10.09:

| csm_err_safety | 1 | 2 | 3 | increment test |
|---|---|---|---|---|
| worst distT | 6.3e-09 | 3.1e-09 | 1.4e-09 | -- |
| worst distC | 7.6e-09 | 3.8e-09 | **2.1e-09** | -- |
| over 1e-08 | 0/299 | 0/299 | **0/299** | -- |
| passes | 9.13 | 9.43 | **10.01** | 10.09 |

3 is the smallest of the three that leaves no accepted state outside either
tolerance, and it is the default. The step from 2 to 3 costs 0.06 passes on
the first case, because past a factor of 2 most exits are taken on the
increment anyway: the guards refuse a pair whose increments have reached
round-off, and MEASURED 48 of 300 exits carry an estimate on
mol_base_handoff against 271 of 300 on wasp_full.

WHAT THE CHANGE IS WORTH IN PASSES IS SMALL, and that corrects COST6's
estimate of 1.0 to 1.6 passes. That estimate treated the projected L2
remainder as the max-norm error; with the norm carried honestly the admitted
increment is only as much looser than the tolerance as
`csm_err_safety*|theta/(1 - theta)|` allows, and MEASURED from the pass
trace that quantity is 0.18 (temperature) and 0.21 (composition) on
mol_base_handoff, i.e. 0.74 and 0.68 decades, but 0.74 and 0.72 on
wasp_full, i.e. 0.13 and 0.14 decades. The whole table, control against the
new test, one thread and sixteen:

| case | steps | passes | sweeps | wall 1 thr | wall 16 thr | log10 Mdot |
|---|---|---|---|---|---|---|
| `wasp_full` | 300 | 10.09 -> 10.01 | 3028 -> 3002 (-0.9%) | 139.4 -> 138.7 s | 16.2 -> 16.5 s | 7.21 -> 7.21 |
| `mol_base_handoff` | 1500 | 7.19 -> 6.96 | 10780 -> 10442 (-3.1%) | 272.1 -> 264.3 s | 35.6 -> 35.7 s | 7.80 -> 7.80 |
| `mol_metals` | 400 | 6.99 -> 6.93 | 2796 -> 2772 (-0.9%) | 143.7 -> 142.5 s | 17.2 -> 17.3 s | 7.94 -> 7.94 |
| `mol_carrier` | 300 | 6.01 -> 5.01 | 1802 -> 1503 (-16.6%) | 46.6 -> 39.4 s | 6.3 -> 5.5 s | 7.95 -> 7.95 |
| `hydrostatic_column` | 300 | 1.01 -> 1.01 | 302 -> 303 (+0.3%) | 5.3 -> 5.3 s | 0.9 -> 1.0 s | 8.27 -> 8.27 |

`mol_carrier` is the one case where the lever is real, a whole pass and a
sixth of the sweeps, and its 300 exits are all taken on the estimate.
`hydrostatic_column` gains ONE sweep: its loop takes 1.01 passes and reaches
a third pass on one step only, where the estimate says the error exceeds the
increment and the test is therefore the STRICTER of the two. Where the
estimate is available at the exit pass is what decides the size: MEASURED
783 of 1500 exits on mol_base_handoff, 271 of 300 on wasp_full, 300 of 300
on mol_carrier, 39 of 400 on mol_metals and 1 of 300 on
hydrostatic_column. The state moves by 4.3e-10 (hydrostatic_column),
8.1e-09 (mol_metals), 1.4e-08 (wasp_full), 4.9e-08 (mol_base_handoff) and
1.1e-07 (mol_carrier) in units of the largest column entry, which is what a
different stopping point inside the same ball has to do.

What the change buys instead is the MEANING of the two tolerances: the
accepted state's distance to its fixed point is now measured at 2.2e-09 to
7.4e-09 against tolerances of 1e-08, where under the increment test it was
one to two decades inside them.

THE FALLBACK IS THE OLD TEST, BIT FOR BIT. With the estimate refused
(`EXHALE_CSM_GEOM_REFUSE`) the new binary reproduces the pre-COST7 binary on
every written file of wasp_full at 300 steps (10 files) and
mol_base_handoff at 1500 (7 files), at the same pass count (10.09 and 7.19)
and with every exit taken on the increment.

THE PROBE DOES NOT MOVE THE TRAJECTORY. It restores the pair the accepted
pass was taken FROM and retakes that pass, so a probe run writes the same
files as a run without it: MEASURED on all four cases probed
(wasp_full 300, mol_base_handoff 1500, mol_metals 400, mol_carrier 300),
every written file identical and the same pass count. That is the same
statement COST6's rollback rests on, measured again: the pair (T, f_sp) is
the whole iterate of this loop. Over those four runs no accepted state is
outside either tolerance: worst distance 1.4e-09 / 2.1e-09 (wasp_full),
2.0e-09 / 7.4e-09 (mol_base_handoff), 2.0e-09 / 7.4e-09 (mol_metals) and
2.1e-09 / 2.2e-09 (mol_carrier), temperature norm first.

THE EXTRAPOLATION OF SECTION 5 NO LONGER BUYS ANYTHING and stays default
off: MEASURED on mol_base_handoff at 300 steps, 6.79 passes and 2038 sweeps
with it on against 6.82 and 2046 without, where against the increment test
it was worth 6.15 against 6.96. Its window is now the band between the
tolerance and `csm_x_reach` times it, and the passes it used to save are the
ones the error test no longer spends: 890 attempts, 851 refused by the
guards, 39 candidates taken, all 39 met by the stopping test on the very
next pass, and the admissibility record identically zero.

THE PROBE'S OWN LIMIT, for the record. On mol_base_handoff at 1500 steps
127 probes of 1500 reach the probe's 59-pass cap with the composition
increment still at 5.6e-09, so on those steps the "fixed point" is not
resolved to 1e-12 and the distance measured against it carries that
uncertainty. The headline numbers are unaffected: the worst distance comes
from the 1373 probes that DID converge (7.4e-09), and the stalled ones give
smaller distances (6.4e-09). What the stall says is that the composition
increment of that case has a floor of order 5e-09, which is where a
tolerance of 1e-08 sits.
