# `atomic_elem_newton` -- the atomic element-row coupled solve that makes no progress

## Purpose

The fixture for item N7 (B5m) of `docs/PLAN_20260909_rev1.md`: HD 209458 b,
atomic (no molecular chemistry), solar trace metals, `He_diffusion: True`,
`He_metal_diffusion: True`, `Coupled carrier solve: True`, `Include He23S? True` (the helium triplet is ON, which sends the sweep to `hybrd1(ion_system_HeH_TR_metals)` directly, never to `solve_ieq`; N32), `Solver: Newton`,
cold start (`Load IC? False`). Before N1 (2026-09-09) the stationary solve
registered fourteen unknowns per cell (three hydrodynamic, helium, and all
ten trace elements whether or not `metals.inp` carries them); since N1 it
registers eleven (the three absent elements Si, K, S are no longer rows).
On the entry text it stood still: `||R||` 1.892 to every printed digit for eleven
iterations, GMRES iterations 0, a zero dogleg step and a zero trust radius,
1020 residual samples admitted and none refused (2026-09-09, advisor build
of the B5k tree; `run.log`). Not in `DEFAULT_CASES`; the case is a
diagnostic entry point, not a golden.

## Contents

- `input.inp`, `metals.inp`: the configuration.
- `EXHALE_setup.out`, `EXHALE_resolved.out`, `run.log`: the 2026-09-09 run.
- `output/Hydro_ioniz.txt`, `output/Ion_species.txt`: the state handed back
  by that run (exit status 2, "not certified").
- `run.sh`: reruns the case in place with the binary given as its argument
  (default `../../../EXHALE.x`). `Load IC? False`, so a rerun marches from
  the cold start again; the marching phase is the larger part of the wall
  time.

## Context

`examples/14_diffusion` (2026-08-09) reached info 2 at `||R||` 5.9e-3 and
`benchmarks/hd209` (2026-08-15) info 0 at 7.4e-4 in the old norm, so the
configuration has converged before with earlier solver text; what stands
still now is the element-row route after B5j/B5k. The first two iterations
of the corrected solver (N1 to N3) are to be logged here with named linear
outcomes before anything else is changed (PLAN rev 1, section 4, N7).

## Reload entry (added 2026-09-09)

`IC/Hydro_ioniz_IC.txt` and `IC/Ion_species_IC.txt` are the handed-back
state of the 2026-09-09 run renamed for `Load IC? True`. Copy them into
`output/` of a run directory, set `Load IC? True` in its `input.inp`, and
the run skips the cold-start marching (two CFL steps, then the Newton): this
is the route for N7's first-two-iterations log. The grid guard accepts them
(same `Grid cells:`).

## Control reload, 2026-09-09 (advisor build of the B5l tree, entry text of N1 to N3)

MEASURED (`OMP_NUM_THREADS=16`, 8 s): the reload takes two CFL steps and
enters the Newton at `||R||` 1.888; it reports 1500 species unknowns on a
bound with smallest adopted value 0; the initial trust radius is 0.000 because
the first `pgmres` call took 0 products of 40 (relative residual reached
1.0); every iteration then reports "the dogleg step is zero"; `info = 2`
after 11 iterations, `||R||` unchanged. The 1500 unknowns are the three
elements `metals.inp` does not carry: Si (3 stages), K (2), S (2) are
identically zero in all 504 cells of `IC/Ion_species_IC.txt`, and 3 absent
elements x 500 physical cells = 1500. An absent element's mass-fraction
unknown sits exactly on its lower bound with an identity-like row, and the
box cut refuses a direction whose component there points below zero, which
is how the first Krylov direction is never sampled. Log:
`<scratchpad>/N7/ctrl_reload/run.log` of the 2026-09-09 session; the same
run is reproduced from `IC/` in seconds.

## After N1 to N3 (2026-09-09, Stage A binary md5 `dd233d6a`)

MEASURED on the same reload: unknowns on a bound 0 (the absent elements are
not registered), smallest adopted species unknown 2.820e-5, initial radius
1.909e-2 measured on a Krylov step that reached its tolerance (3 products,
relative residual 6.1e-2 against 0.1 asked), first step accepted with ratio
0.995. `||R||` 1.888 falls to about 1e-1 within 130 iterations and then
oscillates between 7.5e-2 and 1.1e-1 with iterations refused because the
model and the true slope differ in size (the worker's 119-iteration run and
the advisor's rerun agree). Not converged; the next obstruction is named in
`docs/ISSUES_20260909.md` section 3.1 (d).

## N7 diagnostic (2026-09-09)

The refused trials of the oscillation have a CORRECT model (predicted and
actual reductions agree to 3e-4); what refuses them is the ray test whose
probe sits 11x above its assumed floor, and what keeps `||R||` swinging is
that the merit reads the element rows (0.44 of it) while `||R||` does not.
`EXHALE_ELEM_DIAG=1` prints the full ordered diagnostic at iterations 1, 2
and 60 (`=n` moves the third). `Solver: Newton 100.0` does not cap the
iterations (the third word is the hand-off flux metric); the JFNK cap is the
literal 500 in `EXHALE_main.f90`. Report: the N7 report of the 2026-09-09
session; summary in `docs/ISSUES_20260909.md` section 3.1 (d).

## After N20 and N8a (2026-09-09, later)

The two sections above describe the element-row solve between N1 and N20. N20 (decision
20 a: the merit reads every row on its certification scale) removed the
merit/gate mismatch and, on this fixture, drove the solve to a no-descent
abort at `||R||` 1.78: the certification row scales span 11.4 decades and
at iteration 11 the carbon unknown of cell 500 sat at zero with its column
scale on the absolute floor 1e-20, so its Jacobian column was exactly empty
and the banded preconditioner singular (`dgbcon` 4e-27). N8a floored the
element column scale at 1e-4 of the reservoir and equilibrated the
factorized band (the model untouched): the reload then reaches `||R||`
3.8e-2 at the 150-iteration cap, still descending (2.5e-1 to 3.8e-2 over the
last thirteen iterations), `dgbcon` 1.7e-6 at iteration 1. The ray test is
no longer what refuses steps. Not converged; N8a report of the 2026-09-09
session; `docs/ISSUES_20260909.md` 3.1 (d).

## After N21 (2026-09-09)

From the N8a state the single solve stops on the stagnation detector at
iteration 227 at `||R||` 2.671e-2 (170 accepted, 55 refused, 54 of them
"the reduction ratio is below eta"); no reset of the trust-region state
helps (N21). The run's later solves start from MARCHED states (`info = 1`
returns the run to marching for 2000 steps): 6.06e-2 -> 5.27e-4 after the
second, 4.26e-4 after the third. Below `||s||` 1e-7 the predicted decrease
stops scaling with the step (N22).

## After N22 (2026-09-09)

The refused trials of the stall were those whose species unknown was written
onto a face of its box; with the step held at the iterate there the solve
reaches `||R||` 1.356e-2 in 113 iterations (91 accepted, 19 refused) on the
stagnation detector, and the predicted decrease scales with the step again.
`Restart intent: stationary` (N10) measures the state AS LOADED at `||R||`
1.8924 without a CFL step.

## After N23 and N24 (2026-09-10)

110 trust-region iterations (91 accepted, 19 refused), `||R||` 1.356e-2 on
the stagnation detector. After outer iteration 59 the binding row is the
sodium element row of cell 500 (r 4.06 Rp): 0.44 of the merit, the whole
judged distance, the smallest scaled row and column of the band; the Krylov
cycle reaches 0.58 of its right-hand side with the Arnoldi image 21 percent
off the operator. The image-only leg test is worse (2.30e-2). N25 works on
the Arnoldi gap, a direction on the ball, and the row's outer boundary.

## After N26 (2026-09-10)

With the upper ghosts of the transported columns following the iterate,
the reload reaches `||R||` 3.719e-4 (judged distance 6.2e4, flux spread
5.5e-7) in 252 iterations on the stagnation detector; the sodium row of
cell 500 never binds; the binding row is the helium row of cell 246
(r 1.155). Not certified (N27 diagnoses that row).

## After N27 and N26c (2026-09-10)

The binding helium row of cell 246 is held by the linear solve (Krylov 40 of
40 at 0.99). The solve is chaotic at the ulp level: a one-ulp change of the
upper ghosts moves `||R||` at iteration 40 from 3.2e-2 to 1.1 (N26c), so this
fixture's `||R||` at a fixed iteration count is NOT an acceptance quantity;
compare named outcomes (binding rows, refusal counts, flux spread) instead.

## State on 2026-09-10 (after N33/N34, tree 3414478 plus N31 to N34)

MEASURED by N34 with a control build of this tree at `EXHALE_JFNK_MAXIT=250`,
8 threads: the stagnation detector stops the solve at iteration 167 at
`||R||` 3.047e-4, NOT CERTIFIED (11 entries refuse), the binding rows the
element Fe rows of cells 132 to 137 (the N26/N27 numbers above, helium row
of cell 246 at 3.719e-4 in 252 iterations, are the earlier tree's; the solve
is chaotic at the last bit, N26c, so the numbers move with the tree and
the thread count). Cap 40 with every hook off, 8 threads: 1.888, 1.891,
1.819, 1.583, 1.861, 8.964e-1, 7.607e-2 at 1, 5, 10, 15, 20, 30, 40, handed
back 5.377e-2 (advisor, N31 to N34, identical on every tree). The base cell
width is 1.955e-4 (`Base grid` key, `dr_base = 2.0e-4`), `r^2/dV` 5.1e3:
the residual's rounding floor of N33 is proportional to one over it.

## Partitioned route, 2026-09-11 (`docs/PLAN_20260911_partitioned_solver.md` item P7)

The measurement that adopted the partitioned route. Same binary
(`EXHALE_P6.x`, built from the tree of that day's items P1 to P6), 8 threads,
500 cells, reloaded from `IC/`; the whole account is
`docs/Update_EXHALE.md` section 8.

### Recipe

The `input.inp` beside this file is the cold-start COUPLED configuration. Three lines
make it the partitioned reload, on a copy of the case directory:

```
Load IC? True                       # was False
Coupled carrier solve: False        # was True
Restart intent: stationary          # added
```

and then

```
EXHALE_OUTER_PASSES=<n> EXHALE_JFNK_MAXIT=80 OMP_NUM_THREADS=8 ./run.sh /abs/path/to/EXHALE.x
```

`EXHALE_OUTER_PASSES` is the outer pass budget of
`steady_wind_with_element_diffusion` and `EXHALE_JFNK_MAXIT` the cap of each
hydrodynamic solve inside a pass. `Restart intent: stationary` also fixes the
pseudo-time start of the reload at 1.0, which since 2026-09-11 is that
route's default (`EXHALE_PTC_DTAU0` overrides it). That start is what makes
this reload solvable at all: at the CFL start of the older text the first
hydrodynamic solve stopped on its stagnation detector at `||R||` 1.490 and
the pass was then refused because the element relaxation found no admissible
advance, where at 1.0 the same solve hands back 1.167E-08.

The 8-thread convention differs from the 16 threads of the earlier sections;
the solve is chaotic at the last bit (N26c), so compare named outcomes.

### Measured ladder, partitioned, 12 passes

The `(EXHALE_main) outer pass N:` line of each pass. All MEASURED. The
element under-relaxation omega starts at 0.500 and is halved once, to 0.250
at pass 3, by the no-fall rule of the outer contract.

| pass | hydro info | worst gated species row (cell) | hydro mass | hydro momentum | hydro energy | omega | seconds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | O 2.61e-3 (268) | 2.04e-09 | 8.84e-13 | 1.17e-08 | 0.500 | 119.37 |
| 2 | 2 | O 3.81e-3 (268) | 5.12e-10 | 4.30e-13 | 1.15e-08 | 0.500 | 115.82 |
| 4 | 2 | O 9.45e-4 (268) | 7.86e-10 | 1.77e-12 | 1.17e-08 | 0.250 | 75.57 |
| 6 | 2 | O 2.36e-4 (268) | 1.32e-09 | 6.50e-13 | 1.05e-08 | 0.250 | 74.89 |
| 8 | 2 | O 5.97e-5 (268) | 8.88e-10 | 4.36e-13 | 1.08e-08 | 0.250 | 70.78 |
| 10 | 2 | O 1.50e-5 (268) | 8.65e-10 | 1.61e-13 | 8.98e-09 | 0.250 | 94.56 |
| 11 | 2 | O 7.94e-6 (268) | 6.92e-10 | 2.79e-13 | 9.72e-09 | 0.250 | 89.47 |
| 12 | 2 | He/H 5.47e-6 (288) | 1.05e-09 | 1.67e-13 | 1.50e-08 | 0.250 | 54.21 |

Wall time 1043 s for the twelve passes (the sum of those seconds, which
agrees with the run directory's timestamps to under a second). The one rise,
pass 1 to pass 2, is the alternation's own feedback across a hydrodynamic
solve, which is why a single rise does not end the loop
(`outer_no_fall_max = 3`).

Outcome: NOT CERTIFIED on **one** entry, the hydrodynamic mass row at
1.054E-09 against 3.0E-12 at cell 16 (r = 1.0031). Every elemental wind row
is within its 1.0E-05: He/H 5.468E-06 at cell 288, O 3.720E-06, Ca
1.648E-06, N 1.638E-06, Fe 1.630E-06, C 1.619E-06, Na 1.588E-06, Mg
1.582E-06. So the alternation has done its work and what is left is the
base-layer mass row of the three-unknown solve: all twelve hydrodynamic
solves end with `||R||` between 8.980E-09 and 1.501E-08 against the run's
`Resid tol` of 1E-08, and the last one ends on "no descent direction exists
for the banded model at this state" at 1.501E-08. A shorter run of the same
The coupled solve, 5 passes and 528 s, still has 9 entries refusing (the worst gated row
4.64E-04), so the refusal count falls 9 to 1 between pass 5 and pass 12.

### The coupled solve of the same day, same binary and threads

`input.inp` with `Load IC? True` and `Restart intent: stationary` added and
`Coupled carrier solve` left at True, `EXHALE_JFNK_MAXIT=250`, 8 threads.
This coupled solve still alternates the element relaxation, because `He_diffusion` is
on; what "coupled" changes is that the element rows are carried as Newton
unknowns inside each hydrodynamic solve. Three passes completed, at 1380.59
s, 482.04 s and 528.98 s (2392 s in all), and the run was ended in the
fourth.

| pass | hydro info | worst gated species row (cell) | hydro mass | hydro momentum | hydro energy |
| --- | --- | --- | --- | --- | --- |
| 1 | 2 | He/H 2.90e-4 (288) | 1.18e-04 | 2.42e-05 | 3.15e-04 |
| 2 | 2 | He/H 1.58e-4 (292) | 1.90e-06 | 3.66e-07 | 6.32e-06 |
| 3 | 2 | He/H 8.29e-5 (295) | 1.16e-06 | 2.89e-06 | 2.83e-06 |

So on this fixture the partitioned route stands strictly below the coupled
the coupled solve on every judged row (1.054E-09 against 1.16E-06, 1.665E-13 against
2.89E-06, 1.501E-08 against 2.83E-06, worst gated species row 5.468E-06
against 8.29E-05) in 1043 s against 2392 s, and neither route certifies it.
The state AS LOADED, for reference, refuses on 11 entries at mass 1.509,
momentum 3.348E-01, energy 1.892 and a worst gated element row of 5.840E-01.

### Re-measured the same day by item P10

Item P10 made the best-iterate ledger rank iterates by
`max(mass/3e-12, momentum/1e-8, energy/1e-6)` instead of by `||R||` over the
run's `Resid tol`. The same twelve-pass recipe, re-run against a control
built from the entry text of `steady_newton.f90`. All MEASURED, READ from the
P10 report.

| | control (the entry text) | item P10 |
| --- | --- | --- |
| passes, wall time | 12, 1043 s | 12, **1064 s** |
| the one refusing entry, tolerance 3.0E-12 | mass row 1.054E-09 at cell 16 | mass row **8.893E-10 at cell 3** |
| solves reporting a stalled best iterate | 5 | 3 |
| solves ending "no descent direction exists" | 11 | 10 |
| worst gated species row of each of the twelve passes | the ladder above | identical, cell included |

The two differ by the ledger's choice of iterate and not by the elemental
relaxation, which is identical in all twelve passes; that identity is the
check that says the change stays inside the hydrodynamic solve.

What the re-measurement adds, and it is the reason this entry is an
ACCEPTANCE question and not a solver item: the mass row at this level carries
no signal. At pass 1 the ledger's best distance is 4.407E+02, that is a mass
row of 1.32E-09, while the state handed back carries **4.11E-09** after the
final evaluation at its own composition. The row IS the N33 rounding of the
base-layer flux difference, so ranking iterates by it ranks the luckiest
rounding of a quantity that moves when the same state is evaluated again,
which is why the mass rows of the ladder above scatter rather than descend in
the control and in the new run alike.

The 2 percent in wall time is inside the contention of a shared machine and
is not a claim.

## After P16 and P17 (2026-09-11)

The same twelve-pass partitioned recipe, after the user's decision of that
afternoon anchored the continuity row's tolerance on its measured rounding
floor (item P16, `max(3e-12, min(1, 10 * floor(j)))`) and the best-iterate
ledger was made to read that row at the cell's own tolerance (item P17). Each
item measured its own build against a private control built from the entry text
of the files it owns; 8 threads, `EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80
EXHALE_DIFF_OMEGA=0.5`, reloaded from `IC/`. All MEASURED, READ from the two
reports.

**THIS FIXTURE NOW CERTIFIES**, at pass 12, in P16's build and in both builds of
P17:

```
(EXHALE_main) outer pass 12: ACCEPTED -- every active equation of this
                             state is within its own tolerance.
```

P16's control is the entry text of P16's two files, so it predates the anchor;
P17's control is the entry text of P17's file on the tree that already carries
the anchor, so it certifies as well, and what P17's two builds differ in is the
ledger and not the verdict.

| | control of P16 (before the anchor) | P16 | P17 (both builds certify; the column is P17's fix) |
| --- | --- | --- | --- |
| verdict | NOT CERTIFIED, 1 entry | **CERTIFIED at pass 12** | **CERTIFIED at pass 12**, and so is its own control |
| the mass row's verdict | 7.159E-10 at cell 13, against the fixed 3.0E-12, REFUSED | cell 34, r = 1.00665: 6.542E-10 against **4.6E-09, the rounding anchor**, distance 0.141 | cell 50, r = 1.00977: 4.245E-10 against **3.2E-09**, distance 0.132 |
| largest mass row of that state | 7.159E-10 at cell 13 | 7.703E-10 at cell 2 | 9.262E-10 at cell 1 |
| momentum, energy | 3.983E-14, 1.513E-08 | 1.562E-14, 9.981E-09 | 2.439E-14, 2.173E-08 |
| worst gated species row at pass 12 | 6.84E-06 (O, cell 268) | identical to its control to three digits at every pass | 7.73E-06 (He/H, cell 289), its control 6.84E-06 (O, cell 268) |
| passes, summed pass seconds | 12, 1040.83 | 12, 1046.99 | 12, 1096.24 (its own control 1047.35) |

P16's two builds are identical in every reported number at passes 1 to 5 and
part at pass 6, the first pass whose hydrodynamic solve stops on the loop-top
test because its mass row is inside its own tolerance (`info = 0` where the
control's is `info = 2`). The anchored tolerance still REFUSES states on this
fixture: the mass row is admitted at 11 of the 12 passes and refused at pass 3
at distance 69.9 (cell 160). The eight elemental wind rows of the last pass
are identical in every digit between P16's control and P16's fix (He/H
2.793E-05, C 1.000E-05, O 1.238E-05, N 1.496E-05, Mg 4.377E-05, Ca 1.942E-05,
Na 1.941E-05, Fe 8.966E-05), which is the check that the change stays inside
the acceptance and does not touch the relaxation.

**What P17 moved is the quantity the ledger ranks on**, and it is the number
that says the item worked: the `best judged iterate: distance` of the twelve
passes runs

| | control | P17 |
| --- | --- | --- |
| best judged distance, the twelve passes | 4.407E+02, 6.310E+05, 1.828E+03, 2.802E+02, 2.340E+02, 2.245E+02, 3.102E+02, 3.080E+02, 2.602E+02, 2.350E+02, 2.716E+02, 2.568E+02 | **1.570E-01, 1.646E-01, 1.416E-01, 1.555E-01, 1.576E-01, 1.588E-01, 1.469E-01, 1.868E-01, 1.445E-01, 1.640E-01, 1.400E-01, 1.321E-01** |
| solves reporting "best iterate has not improved" | 4 | 7 |
| solves ending "no descent direction exists" | 8 | 11 |
| solves ending STAGNATED | 1 | 0 |
| outer iterations, seconds each | 861, 1.216 s | 909, 1.206 s |

So the control's ledger ranked every iterate of every solve as 2.2E+02 to
6.3E+05 times outside certification, on a row its own acceptance admits at a
distance of 0.13, while the fixed ledger reports the number the verdict is
taken on. The elemental relaxation is not damaged: the worst gated species row
descends 2.61E-03 to 7.73E-06 against the control's 2.61E-03 to 6.84E-06, both
inside 1E-05 at pass 12, and the certified state's elemental rows are SMALLER
in P17's fix on every element (Fe 5.086E-05 against 8.966E-05, C 3.70E-06
against 1.00E-05, O 3.34E-06 against 1.24E-05). The cost of the 500
rounding-floor evaluations added to each judged residual is not separable from
the machine's noise (1.206 s against 1.216 s an outer iteration).

The 2 percent and 5 percent movements in wall time above are inside the
contention of a shared machine and are not claims. The controls of P16 and P17
are the CURRENT tree's entry text and not P10's binary, so their ladders are
not P10's number for number (P10 read 8.893E-10 at cell 3 after twelve
passes, P16's control 7.159E-10 at cell 13); what every control agrees on is
the finding, that before the anchor the ONE entry refusing after twelve passes
was the mass row of a base-layer cell.
