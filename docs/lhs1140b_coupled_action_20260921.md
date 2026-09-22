# LHS 1140 b: the coupled action against the probe arc, and the standing coupled-block failure

2026-09-21. Phase 6b of `PLAN_20260920_rev9.md` section 15: the Stage B of
section 6 run ONCE MORE in the COUPLED configuration, which section 9.4 makes
an explicit prerequisite of the phase-7 trial, together with item R5 of
`PLAN_20260919_rev1.md`, the one standing failure of
`src/tests/coupled_block_jacobian/`.

## 1. The two verdicts, first

**(a) Is the coupled action resolved?** **Yes on the unknowns the coupled
system adds, and no on the base cell's mass row and on the banded model's own
proposed step.** In the terms section 8.2 sets, at
`molecular_scalar_gj1132_kzz1e9/HeH0.083`, generation
`g0004_20260920T053734Z_79b42a03`, with `Coupled carrier solve: True` and the
He element row and the H2 carrier row registered (5 unknowns per cell, 2500
entries):

- **Stabilizes in `h`?** Direction by direction, relative to the action at the
  production arc, over `h` from 1E-02 to 1: `cell_species` (the species
  unknowns of the cell whose carrier row stands furthest outside its
  certification) **2.8E-06**, and it holds that to `h = 100`; `smooth`
  **1.4E-07**; `cell_energy` **4.5E-04**; `cell_mass` **1.4E-02**;
  `banded_step` **3.6E-01**, with no plateau anywhere.
- **Moves with `k`?** No. `k = 1` reproduces the converged-closure residual
  **bitwise**, and at the production arc the action at a fixed count differs
  from the action at the run's own converged closure by at most **4.1E-06**
  over the counts `k >= 2` on the four directions other than `banded_step`
  (`cell_species` 4.1E-06, `smooth` 1.3E-09, `cell_energy` 5.5E-10,
  `cell_mass` 3.6E-10), and by **2.7E-02** on `banded_step`. The one larger
  entry, 6.5E-05 on `cell_species` at `k = 1`, is the endpoint-pass asymmetry
  of the converged reference that section 8.1 of the plan names and not a
  movement of the map: at `h >= 1` the endpoints of the converged row
  completed two passes and the base one, so `k = 2` and not `k = 1` is the
  count that matches them there, and at `k = 2` the distance is exactly zero.
- **Is the inner map resolved?** Yes: the elimination exits after ONE pass at
  an increment of 1.461E-13 against the 1.0E-08 it is asked for, and **all 315
  actions of the grid are admissible**: `ok=TTT` at every point, no cell left
  without a chemical root, no refusal, no blocked component, no backward
  sample, and the arc requested equal to the arc used everywhere.

  What limits the two directions that do not resolve is **the base cell's mass
  row**, the object sections 13.2b and 11.3 of the plan and section 7.2 of
  `docs/lhs1140b_closure_and_probe_20260921.md` already name, reached here from
  a fourth side. The residue is a fixed ABSOLUTE size, about 6E+04 of the row
  scales, independent of the arc over three decades and sitting in the mass
  rows; it is not rounding (it does not grow as `1/h`) and not a closure effect
  (it does not move with `k`).

**(b) Does the carrier action really disagree with a central difference?**
**No. The test is measuring its own probe arc.** The statistic the failing row
gates on is the distance between the production FORWARD difference and a
CENTRAL difference at the SAME arc, which is the first-order truncation of the
forward difference, and it is exactly proportional to the arc:

| probe arc | 1E-04 | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|---|---|
| energy direction, mean over the rows it moves | 1.2390E-06 | 1.2388E-05 | 1.2396E-04 | 1.2381E-03 | 3.7034E-03 | **1.2379E-02** | 3.4546E-02 | 1.1897E-01 |
| largest response in the support | 2.4963E-03 | 2.4964E-03 | 2.4968E-03 | 2.5010E-03 | 2.5104E-03 | 2.5428E-03 | 2.6340E-03 | 2.9339E-03 |

Four decades of exact first power, with no floor reached; the reference the
ratio is taken against moves by 0.19 per cent over the four decades from
1E-04 to 1E-01 and by 1.9 per cent from 1E-04 to the production arc, so the
DERIVATIVE is resolved and it is the forward difference's distance from it
that the row reads. **The disagreement does not survive inside the resolved
interval**: at arc 1E-01 the same statistic is **1.2381E-03**, ten times
inside the stated 1.0E-02, with the central difference unchanged to 0.16 per
cent. Section 10 says what the test should assert instead.

**Consequence for phase 7.** (a) does not fail, so the gate of section 15 is
not closed against the trial, but it is passed with a named limit: the coupled
action is resolved to 3E-06 on the species unknowns and to 1E-07 on a smooth
direction, and to no better than 1E-02 on the base cell's mass row and 3E-01
along the step the banded model proposes. A phase-7 verdict read off the
species rows and the outer column is supported by this measurement; a
phase-7 verdict that turns on the mass row of cells 1 to 3 at a relative
accuracy better than 1E-02 of the action is not.

## 2. Classification, and what was inherited

Every new number here is **Mode R** (rev9 section 3): the current operator
evaluating an immutable stored state. No physical step was taken, no state was
published, nothing under `LHS1140b/models/` was written and nothing under
`backup/regression/` was read for writing or written. Each invocation ran on
its own copy of the stored state pair in its own scratch directory.

Inherited and not re-measured:

- **The residual is a state function at these checkpoints**
  (`docs/lhs1140b_repeatability_20260921.md`), so a directional difference is a
  statement about the operator and the independently measured evaluation
  uncertainty section 6 item 2 asks for is zero.
- **The uncoupled Stage B of the same checkpoint**
  (`docs/lhs1140b_closure_and_probe_20260921.md` sections 8.3 and 8.4). It is
  not merely cited: the same grid was re-run here with the new binary and its
  whole diagnostic block is **byte-identical** to the one that document
  reports (section 4 below), which is what makes the coupled and the uncoupled
  numbers comparable.
- **`EXHALE_RESID_QUAD` cannot fire on a diagnostic route**, INSPECTED at
  `steady_residual.f90:246-248` and MEASURED again by every invocation here:
  `ieq_sweep_state_kind 1, is the marching kind: T`.

## 3. Identities

### 3.1 The binaries

| | |
|---|---|
| measured build | `<scratch>/coupled_stageB/EXHALE_test.x`, md5 MEASURED `d4065fc49b9351f6c762d6ffccfff608` |
| built from | the working tree at repository HEAD `93eed86` (dirty), plus the extension of section 4 |
| build command | `make OBJDIR=<scratch>/obj EXE=<scratch>/EXHALE_test.x -j8` |
| source digest | `1a377a84ccd34fd7d7e7f1741956369b`, the md5 of the sorted md5 list of the 245 `.f90`/`.inc` files of `src/` at the time of the build (MEASURED). A live tree with other workers in it, so this records what stood beside the binary |
| tree binary | `EXHALE.x` md5 MEASURED `89ec67149aa452eea3deeb19052f83d1`, unchanged and not rebuilt. It produced the cross-check of section 10.3, which needs no new entry point; every table that needs the extension names the measured build |
| compiler strings | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`; `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` |
| BLAS and LAPACK | `/opt/miniconda3/lib/libopenblas.so.0` |
| OpenMP runtime | `/opt/miniconda3/lib/libgomp.so.1` |
| Fortran runtime | `/opt/miniconda3/lib/libgfortran.so.5`, with `libquadmath.so.0` |

### 3.2 The two checkpoints

| | `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `molecular_scalar_gj1132_wellmixed/HeH0.55` |
|---|---|---|
| why | the phase-7 family, whose Stage A0 identity record is complete (`docs/lhs1140b_identity_records_20260921.md` section 3.4) | the fixture `src/tests/coupled_block_jacobian/run.sh` uses, so that (b) is read against the failing row itself |
| state | `g0004_20260920T053734Z_79b42a03` | the suite's own restart pair, `output/*_IC.txt`, which is state `g0002_20260915T221824Z_d7347336` |
| `Hydro_ioniz` md5 MEASURED | `fd882d8918052069f8ed229dd376f8c4` | `5a9de6ae680246c6badbec424d60b8dd` |
| `Ion_species` md5 MEASURED | `87288b03ad88a92b7b8b8b6c2b1bdd66` | `b0635219ad30ce5551a066e3b890c90e` |
| coupled system (MEASURED) | 5 unknowns per cell over 500 cells, 2500 entries, 2 species rows: the He element row (`He_diffusion: True`) and the H2 carrier row | 4 unknowns per cell, 2000 entries, 1 species row: the H2 carrier |
| species unknown space (MEASURED) | a carrier is ln n, an element is its mass fraction | the same |
| seed composition digest (MEASURED) | `34C30AC3DDA15380` | `732802108C90B266` |
| state digest, coupled packing (MEASURED) | `1857722DF4D03008` | `0420A853A4E34B83` |
| admissibility of the configuration | `Coupled carrier solve: True` is accepted: `Molecular carrier transport: True` is set and `Ionization transport` is absent, the two conditions `input_read.f90:2311` and `:2371` refuse on (INSPECTED; that file is being edited by other work, so read the two `ERROR: "Coupled carrier solve"` screens by their text) | the same |

Both are Mode R and neither is a reproduction of its historical certificate.

### 3.3 The effective configuration at every evaluation

MEASURED, printed by every invocation:

| | `kzz1e9/HeH0.083` | `wellmixed/HeH0.55` |
|---|---|---|
| reconstruction | WENO3, PLM to WENO3 continuation off | the same |
| judged row scaling (`EXHALE_JUDGED_ROWS`) | F | F |
| `ieq_sweep_state_kind` | 1, the marching kind | the same |
| active species bounds held | 0 | 0 |
| closure policy of the run | sweep to the fixed point, at most 25 passes, cap -1, increment asked for 1.0E-08 | the same |
| carrier floor (code density units) | smallest 1.694E-31, largest 8.154E-22 | smallest 9.788E-32, largest 1.119E-21 |

The reconstruction MODE is frozen by `select_stationary_reconstruction`; the
WENO smoothness weights, the contact-mode dissipation and the positivity
repairs are RECOMPUTED at every evaluation, because each evaluation runs the
whole flux assembly. The active set is empty and the projection displacement
is zero at every point of both grids (section 9), so every difference below is
a difference of the unconstrained map.

**The row scaling is not the same array as on the uncoupled route, and that is
the code's own rule, not a choice made here.**
`row_scaling_of_the_linear_model` (`steady_newton.f90:2785-2810`) returns
`cell_row_scales` wherever a species row is registered and the state column
scales otherwise, so `Dr` below is the certification row scale on the coupled
grid and the column scale on the uncoupled one. Numbers divided by `Dr` are
therefore comparable ACROSS ARCS AND COUNTS within one grid, and not across
the two grids; the judged scale `Dj`, which is the tolerance times the
certification scale on every row, is given beside it where that matters.

**The second factor the brief names.** Registering species rows also switches
the globalization to the scaled trust region (`EXHALE_main.f90:7066` then
`steady_newton.f90:17451`). Nothing in this document runs a solve, so the
globalization never enters; it is stated here because every comparison below
with an uncoupled number is a comparison of two operators only, and would not
be if a solve were involved.

## 4. What was added to the control, and what it cost

`EXHALE_CLOSURE_PROBE=1` already existed (phase 4) and its reader already
registers the species rows from `carrier_in_newton`
(`EXHALE_main.f90:1571-1582`, INSPECTED), so the coupled configuration was
reachable without a new key. Four things were added inside the diagnostic, all
off the production route, and one of them was a defect the coupled
configuration exposed.

1. **The carrier unknown space and the species box are now formed, in the order
   `solve_steady_jfnk` has.** INSPECTED: the solve packs the carriers as
   densities, evaluates the residual once (which is what makes the carrier
   operator freeze the element budget), then calls `form_carrier_unknown_space`
   (`:17502`), which rewrites the carrier unknowns as `ln n`, and only then
   `freeze_species_unknown_box` and the scales. The probe did none of that, so
   it measured a DIFFERENT unknown space: MEASURED before the repair, the
   column scale of the H2 carrier where H2 is absent was its own vanishing
   density, the banded model's proposal came out at `||D^-1 s|| = 1.1947E+29`
   instead of `5.4880E+01`, and the failing row of (b) read 2.457E-04 at the
   production arc instead of the 1.2379E-02 the suite reports. After the
   repair the two agree to four digits at every arc (section 10.2). The base
   residual is taken after the space is armed, because `exp(ln n)` is not `n`
   again in its last bits and the base cell's mass row amplifies that by about
   3.5E+01 into a judged row of order one.
2. **The judged size of a species row.**
   `judged_size_of_the_hydrodynamic_rows` writes only the three conserved
   slots by construction and leaves the rest to its caller, so on a
   species-row system the probe was reading uninitialized memory for the
   species slots of `Dj`. They are now filled with the same product the
   certification forms, `cert_tol_carrier_at(r(j))` or
   `cert_tol_element_at(r(j))` times the row's certification scale.
3. **Two further directions.** `cell_species`, the species unknowns of the cell
   whose species row stands furthest outside its certification, which exists
   only where a species row does; and `named_block`, one unknown kind over a
   range of cells, given by `EXHALE_CLOSURE_PROBE_DIRECTION` as
   `<kind>,<first cell>,<last cell>`, built exactly as
   `jacobian_action_of_the_coupled_block` builds the direction it reports on,
   so that the arc ladder can be run along the failing direction itself. Absent
   the variable, no such direction is built.
4. **The statistic the coupled-block measurement reports, at every point of the
   grid**: the relative distance between the production forward difference and
   the central difference of the same map, row by row, over the cells the
   direction sits on, averaged over every row of that support and over the rows
   the direction MOVES at the suite's own `moved_floor = 1.0E-03`.

**The three-unknown route is untouched, and that is measured and not asserted.**
The evaluation is moved across the scales only where a species row exists; on
the three-unknown route the text and the order stand character for character.
MEASURED: the whole 2329-line diagnostic block of the uncoupled grid at
`kzz1e9/HeH0.083` is **byte-identical** (md5 `2001609a92667bbc8f28dff20b9927bd`)
to the block the same case produced before any of this was written; its state
digest is `A52CCB800ACD9772`, the one
`docs/lhs1140b_closure_and_probe_20260921.md` section 3.2 records for the
molecular checkpoint; and its numbers are the ones sections 8.3 and 8.4 of that
document publish, digit for digit (`banded_step` 2.890E-01, 3.038E-02,
3.768E-03, 1.068E-03, 0, 4.331E-04, 3.765E-04, 3.743E-04, 5.243E-04 in `h`). The intermediate build that moved the evaluation on BOTH
routes is what showed why the branch is needed: it moved all 252 digests and
the fifth digit of the largest forward-to-central distances, because the
column scales are formed from module data an evaluation refreshes.

**The regression with the keys unset.** Two fast cases, on a scratch copy of
the harness and of its references, so that nothing under `backup/regression/`
was touched (a full matrix was running there under the repository harness at
the time):

```
REGRESSION_EXE=<scratch>/EXHALE_test.x ./run_check.sh check hp_front roundtrip
[build] skipped: REGRESSION_EXE=... (b3c669cefd25)
[hp_front]  PASS Hydro_ioniz.txt / Ion_species.txt / Hydro_ioniz_adv.txt / Ion_species_adv.txt (data identical)
[roundtrip] PASS Hydro_ioniz.txt / Ion_species.txt / Hydro_ioniz_adv.txt / Ion_species_adv.txt (data identical)
==> REGRESSION PASS (all cases byte-identical)
```

MEASURED under `env -i` carrying no `EXHALE_*` variable at all. The log is
`docs/audit_20260905/coupled_action_20260921/regression/regression_final.log`;
it names the build that first carried the branch, and the final binary
`d4065fc4` differs from it only in where the base residual of the diagnostic is
taken, which no production route reaches. No golden was written and `golden`
mode was never invoked.

## 5. What the grid measures

The grid of phase 4, unchanged in its two dimensions and its criteria:

- **`k`, the eliminated-closure sweep count**: 1, 2, 3, 4, 6, 8, and `k = 0`
  standing for the run's own policy, taken FIRST so that it is the reference.
- **`h`, the probe arc**: 1, 1E-03, 1E-02, 1E-01, 3E-01, 3, 10, 100, 1000
  times the production arc, the production arc taken first so that the plateau
  is read against it. `EXHALE_JV_COLUMN_SCALE` was left off, so the arc is the
  unscaled rule `sqrt(epsilon)(1 + ||Y||)/||v||`.

63 pairs `(k, h)` times the directions the state carries: five at
`kzz1e9/HeH0.083` (315 actions) and six at `wellmixed/HeH0.55` (378). Each
action costs three residual evaluations, the production `jv_product` and the
two endpoints at the SAME `k` from the SAME seed.

**The directions, built once and held over the whole grid**, every one a unit
vector in the column-scaled coordinates:

| name | what it is | `kzz1e9/HeH0.083` | `wellmixed/HeH0.55` |
|---|---|---|---|
| `banded_step` | the step the solve proposes at this state: the Newton direction of the banded model of the FULL residual in the scaled coordinates the linear system is written in. It is not the Krylov step, which is itself a function of `h` and `k` | `||D^-1 s|| = 5.4880E+01`, 0 colors unresolved, LAPACK info 0 | `2.1276E+02`, 0, 0 |
| `cell_mass` | the three conserved unknowns of the cell whose MASS row stands furthest outside its certification | cell 1 | cell 2 |
| `cell_energy` | the same at the cell whose ENERGY row stands furthest outside | cell 3 | cell 246 |
| `smooth` | `sin(0.1 j)` on every unknown, scaled by the column scales | 2500 entries | 2000 entries |
| `cell_species` | every species unknown of the cell whose species row stands furthest outside | cell 499, the H2 carrier row at 81.5 judged units | cell 306, the H2 carrier row at 2.363E+03 judged units |
| `named_block` | the energy unknown of cells 300 to 312, the direction `src/tests/coupled_block_jacobian` reports its failing row on | not built | 13 entries |

**A cross-check the first direction supplies, and it does not pass here.** On
the uncoupled route the action along `banded_step` returns `||Dr^-1 F||`
divided by `||D^-1 s||` to three digits (phase 4). On the coupled grid at
`kzz1e9/HeH0.083`, MEASURED: `||Dr^-1 J v||_2 = 7.02E-02` against
`9.7197E-02 / 5.4880E+01 = 1.771E-03`, a factor of 40. The colored
finite-difference band is a much poorer model of the full action on the
coupled route at this state than on the three-unknown one, and that is
recorded here because it is the preconditioner phase 7 would solve with.

## 6. The gate of section 4.2

Every invocation ran in a directory of its own, created fresh, holding its own
copy of the input file, its own copy of the stored state pair as the restart
pair, its own copy of the executable and its own empty `output/`. The complete
environment was imposed with `env -i`. The driver is
`docs/audit_20260905/coupled_action_20260921/run_coupled_probe.sh`; the
sidecars, the diagnostic blocks, the reduced tables, the arc ladders of
section 10, the regression log, the reduction `summarize.py` and `MD5SUMS.txt`
over all of them are beside it.

| invocation | case | coupled | threads (OMP, BLAS) | wall s | exit | stage reached | block lines |
|---|---|---|---|---|---|---|---|
| `C0_kzz1e9_uncoupled` | `kzz1e9/HeH0.083` | False | 8, 1 | 222.6 | 0 | the grid completed | 2329 |
| `C1_kzz1e9_coupled` | `kzz1e9/HeH0.083` | True | 8, 1 | 281.8 | 0 | the grid completed | 2898 |
| `C2_kzz1e9_coupled_p2` | `kzz1e9/HeH0.083` | True | 8, 1 | 283.2 | 0 | the grid completed | 2898 |
| `B1_wellmixed_coupled` | `wellmixed/HeH0.55` | True | 8, 1 | 333.5 | 0 | the grid completed | 3467 |
| `B2_wellmixed_coupled_p2` | `wellmixed/HeH0.55` | True | 8, 1 | 333.1 | 0 | the grid completed | 3467 |
| `C4_kzz1e9_coupled_t1` | `kzz1e9/HeH0.083` | True | 1, 1 | 2005.4 | 0 | the grid completed | 2898 |

There is no "no output", no timeout and no unresolved measurement among them.
Each left a nonempty standard output ending in `---- end of the grid ----`,
and the restart pair each was given is byte-identical to the stored state
after the run.

**Reproducible across processes AND across thread counts.** MEASURED: the whole
2898-line block of `C1`, `C2` and `C4` has one md5,
`33ff769e24b4f493e58871ae0f187211`, although `C4` ran on one OpenMP thread and
took 2005 s against the 282 s of the eight-thread runs, which is the evidence
that the threads were used; the 3467-line block of `B1` and `B2` has one md5,
`00eae01bae1d9e1580ff44000b8618a4`. Every number below is therefore the same
number in two independent processes, and on the coupled grid of the phase-7
checkpoint in two thread configurations as well.

## 7. The closure count at the coupled state

`|| Dr^-1 (F - F_jac) ||` at `kzz1e9/HeH0.083`, `F` the residual at the run's
own closure policy and `F_jac` at the fixed count. `||Dr^-1 F||_2` by class,
MEASURED: mass 1.1655E-08, momentum 7.6622E-13, energy 2.3277E-08, species
9.7197E-02, all 9.7197E-02. `||Dj^-1 F||_2`: mass 1.9459E+00, momentum
7.6622E-05, energy 2.3277E-02, species 3.7137E+02, all 3.7137E+02. The largest
judged species row is the H2 carrier of cell 499 at **81.5**.

| k | passes completed | admissible | cells without a chemical root | refusal | closure increment against the converged closure | channel | `||Dr^-1 (F-F_jac)||_2` all | `||Dj^-1 (F-F_jac)||_2` all |
|---|---|---|---|---|---|---|---|---|
| converged | 1 | T | 0 | none | 0 | (none moved) | 0 | 0 |
| 1 | 1 | T | 0 | none | 0 | (none moved) | 0 | 0 |
| 2 | 2 | T | 0 | none | 7.147E-13 | heat - cool | 1.1679E-08 | 9.2207E-01 |
| 3 | 3 | T | 0 | none | 5.931E-14 | heat - cool | 1.3720E-08 | 9.6668E-01 |
| 4 | 4 | T | 0 | none | 7.145E-13 | heat - cool | 1.3669E-08 | 1.0872E+00 |
| 6 | 6 | T | 0 | none | 7.144E-13 | heat - cool | 1.6327E-08 | 1.1503E+00 |
| 8 | 8 | T | 0 | none | 7.144E-13 | heat - cool | 1.6984E-08 | 1.2803E+00 |

`k = 1` reproduces the converged-closure residual **bitwise**, exactly as on
the uncoupled route. Beyond one pass the elimination moves the radiative pair
of cell 499 by 7.1E-13 and that reaches the rows at 0.92 to 1.28 in the units
the certification reads them in, of which 0.92 to 1.28 sits in the MASS rows;
the SPECIES class carries 3E-07 of it. That is the finding phase 4 reported on
the uncoupled grid of this state, unchanged by the coupled registration and
now shown not to be a species-row effect at all.

## 8. The two-dimensional table, `kzz1e9/HeH0.083`, coupled

### 8.1 The plateau in `h`

Distance from the action at the production arc, relative, at `k = 1`
(MEASURED; the other `k` agree to the digits shown except where section 8.2
says otherwise):

| direction | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| `banded_step` | 3.94E+01 | 3.03E+00 | 3.59E-01 | 1.24E-01 | 0 | 3.18E-02 | 3.37E-02 | 3.24E-02 | 3.24E-02 |
| `cell_mass` | 1.38E-02 | 1.38E-02 | 1.18E-02 | 5.26E-03 | 0 | 8.68E-03 | 1.34E-02 | 1.78E-01 | 2.08E-01 |
| `cell_energy` | 4.49E-04 | 4.49E-04 | 4.39E-04 | 2.65E-04 | 0 | 3.28E-04 | 5.23E-03 | 2.60E-01 | 4.32E-01 |
| `smooth` | 1.34E-06 | 1.43E-07 | 1.80E-08 | 1.21E-08 | 0 | 6.70E-02 | 6.70E-02 | 6.79E-02 | 3.08E-02 |
| `cell_species` | 2.06E-05 | 2.76E-06 | 6.34E-07 | 6.21E-07 | 0 | 2.08E-07 | 2.68E-07 | 2.75E-07 | 4.03E-06 |

Read direction by direction, because they do not say one thing:

- **`cell_species` is the best-resolved direction of either grid.** The action
  is `3.8490E+01` to five digits at EVERY arc from 1E-03 to 1E+03, and the
  relative distance from the production arc is at most 2.8E-06 over
  `h` from 1E-02 to 100. The unknowns the coupled system adds, sampled at the
  cell where its own row refuses, carry a derivative that is resolved over five
  decades of arc.
- **`smooth` has a clean plateau `h` from 1E-02 to 1** at 1.4E-07, with a
  rounding floor visible below it (1.34E-06 at 1E-03, falling as `1/h`), and
  then a STEP at `h = 3` to 6.70E-02 which does not grow with the arc: 6.70E-02,
  6.70E-02, 6.79E-02 at `h` = 3, 10, 100. That is not curvature, which would
  grow as the square of the arc.
- **`cell_energy`** holds 4.5E-04 from 1E-03 to 3 and then leaves it.
- **`cell_mass`, at the base cell, never resolves better than 5E-03**, and its
  residue does not fall as the arc is shortened: 1.383E-02 at both 1E-03 and
  1E-02.
- **`banded_step` has no plateau**: `1/h` below the production arc and a floor
  of 3.2E-02 above it.

**Where the residue lives, and it is one object.** MEASURED, the class
decomposition of `||Dr^-1 (J(h) - J(1))||` at `k = 1`: for `cell_mass` at
`h = 1E-03` it is mass 5.63E+04, momentum 0.15, energy 3.77E+03, species
2.32E+04; for `smooth` at `h = 3` it is mass 6.01E+04, momentum 0.06, energy
1.72E+03, species 3.11E+04. **The same fixed absolute size, about 6E+04 of the
row scales, and the same distribution**, on two directions that have nothing
else in common. `cell_species`, which sits at cell 499 and is the only
direction that never touches the base cells, does not carry it: its residue is
1E-04 of the row scales.

This is the phenomenon section 10 of
`docs/lhs1140b_closure_and_probe_20260921.md` records on the FORWARD difference
of the uncoupled route, a symmetric offset of fixed absolute size switching on
above a direction-dependent arc. Here it is in the CENTRAL difference of the
coupled route, and phase 4 found the central difference immune. Its cause is
not established here either. What is established is that it is the base cells'
mass rows, that it is not rounding and not the closure count, and that a
direction which avoids those cells is unaffected.

### 8.2 The movement with `k`

Distance from the action at the converged closure, same arc, relative
(MEASURED):

| direction \ h | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| `banded_step`, k=1 | 0 | 0 | 0 | 0 | 1.13E-01 | 1.09E-01 | 1.09E-01 | 1.09E-01 | 1.09E-01 |
| `banded_step`, k=8 | 6.29E-01 | 1.12E+00 | 3.22E-01 | 1.96E-01 | 2.72E-02 | 1.09E-02 | 3.00E-03 | 2.65E-04 | 2.55E-05 |
| `cell_mass`, k=8 | 4.40E-07 | 2.71E-08 | 5.92E-09 | 4.15E-10 | 3.43E-10 | 3.30E-10 | 3.14E-10 | 3.66E-10 | 3.73E-10 |
| `cell_energy`, k=8 | 2.02E-08 | 1.41E-08 | 3.21E-09 | 1.87E-09 | 5.35E-10 | 5.16E-10 | 5.12E-10 | 5.94E-10 | 5.59E-10 |
| `smooth`, k=8 | 1.18E-06 | 1.19E-07 | 1.95E-08 | 4.00E-09 | 9.50E-10 | 3.50E-10 | 1.08E-10 | 9.44E-12 | 9.44E-13 |
| `cell_species`, k=8 | 6.21E-05 | 6.88E-05 | 6.13E-05 | 6.13E-05 | 3.89E-06 | 3.48E-07 | 3.42E-07 | 5.83E-08 | 1.48E-08 |

At the production arc, over the fixed counts `k >= 2`, the movement with `k` is
at most **4.1E-06** on the four directions other than `banded_step`, and the
three directions that touch the base cells move by less than 1.3E-09. **The
closure count is not what limits the coupled action.**

The `k = 1` entries at `h >= 1` are the endpoint-pass asymmetry section 8.1 of
the plan warns of, and the record names it: at `k = 0` the endpoints of those
differences completed 2 passes while the base completed 1, so the two sides are
different maps, and the `k = 2` row is exactly zero there for the same reason.
At every fixed `k` both endpoints completed exactly `k` passes.

## 9. Admissibility, and the endpoint record

Section 8.1 asks whether `F` is admissible for derivative use at every count,
and section 8.2 rules out any count at which the inner map leaves roots
unresolved.

**Every evaluation of the coupled grid at the phase-7 checkpoint was
admissible.** MEASURED over all 315 actions and their 7 base evaluations:

- the probe, the forward endpoint and the backward endpoint all returned
  `admissible` at every point, **315 of 315**;
- no cell was left without a chemical root at any endpoint, at any `k`, at any
  arc;
- no refusal screen fired;
- no component of any direction was blocked and no probe was taken on the
  backward side, which is what an empty active set means;
- the arc requested and the arc used are the same number at every point: no
  step was shortened by the species box or by the shared element rows.

**There is therefore no `k` and no arc at which the coupled map was
inadmissible for a derivative at this checkpoint.**

At `wellmixed/HeH0.55`, 21 of 378 points are not clean, all of them at
`h = 100` and `h = 1000`: the backward endpoint of `named_block` and of
`cell_energy` is inadmissible there and 13 cells are left without a chemical
root. They are outside every plateau and no statement here rests on them.

## 10. (b) The standing coupled-block failure, against the arc

### 10.1 What the failing row measures, INSPECTED

`the_assembled_action_reproduces_the_central_difference_energy` gates the mean,
over the rows the direction moves, of

```
e = |Jv - cd1| / max(|Jv|, |cd1|)
```

with `Jv` the action the solve uses, a FORWARD difference at the probe arc
(`jacobian_action_of_direction` then `jv_product`), and `cd1` a CENTRAL
difference of the same residual at the SAME arc
(`steady_newton.f90:7533`, `:7543` and `:7552`). A row counts as moved when its
response reaches 1E-03 of the largest response in the support
(`moved_floor`, `:7431`). `e` is therefore the distance between a forward and
a central difference of one map at one arc, which to leading order is
`(h/2) |F''| / |F'|`: **the first-order truncation of the forward
difference**, not a property of the assembly.

MEASURED today with the measured build, the same as the brief reports:
carrier direction 4.6313904E-04 mean and 8.9552239E-03 worst over 30 moved
rows; energy direction **1.2378866E-02** and 4.1032325E-02 over 21 moved rows,
against a stated tolerance of 1.0E-02. The largest response in the support is
2.5428167E-03 (energy) and 2.8031307E-08 (carrier); the 4.79E-12 the brief
quotes is today's value of neither.

### 10.2 The ladder, on the suite's own fixture and with the suite's own statistic

Two independent routes, and they agree.

**Through the suite's own hook.** `EXHALE_JV_PROBE_ARC` is read by
`read_species_unknown_space_controls` and multiplies
`probe_length_of_the_jacobian_action` (`steady_newton.f90:11640-11641`), which
is what `jacobian_action_of_the_coupled_block` takes its probe length from, so
the suite's measurement can be made at a named arc with NO code change at all.
MEASURED, on the suite's fixture, cells 300 to 312, energy direction:

| arc | 1E-04 | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|---|---|
| probe length | 6.868E-06 | 6.868E-05 | 6.868E-04 | 6.868E-03 | 2.060E-02 | 6.868E-02 | 2.060E-01 | 6.868E-01 |
| mean over the rows it moves | 1.2390E-06 | 1.2388E-05 | 1.2396E-04 | 1.2381E-03 | 3.7034E-03 | **1.2379E-02** | 3.4546E-02 | 1.1897E-01 |
| worst moved row | 4.4404E-06 | 4.4406E-05 | 4.2221E-04 | 4.1291E-03 | 1.2353E-02 | 4.1032E-02 | 1.1440E-01 | 3.7855E-01 |
| largest response in the support | 2.4963E-03 | 2.4964E-03 | 2.4968E-03 | 2.5010E-03 | 2.5104E-03 | 2.5428E-03 | 2.6340E-03 | 2.9339E-03 |

Exactly the first power of the arc over four decades, and no floor is reached
even at 1E-04. The reference the ratio is taken against changes by 0.19 per
cent from 1E-04 to 1E-01 and by 1.9 per cent from 1E-04 to 1.

And the carrier direction, which PASSES the same row, has the other half of the
V:

| arc | 1E-03 | 1E-02 | 1E-01 | 1 | 10 |
|---|---|---|---|---|---|
| probe length | 1.804E-09 | 1.804E-08 | 1.804E-07 | 1.804E-06 | 1.804E-05 |
| mean over the rows it moves | 2.9313E-01 | 1.8570E-02 | 2.3836E-03 | **4.6314E-04** | 5.0222E-05 |
| rows counted as moved | 35 | 30 | 30 | 30 | 30 |
| largest response in the support | 2.9300E-08 | 2.8012E-08 | 2.8042E-08 | 2.8031E-08 | 2.8031E-08 |

Falling with the arc from 10 down to 1E-02 and then turning up hard at 1E-03,
where it reads 0.293 with a worst row of exactly 2.0, the maximum this
relative measure can take: that is the ROUNDING FLOOR of the difference, and it
is reached three decades below the production arc. The carrier direction's
production arc sits on the useful part of the V; the energy direction's sits
four decades above its own optimum.

**Through the closure probe, independently.** The same fixture, the same
direction built the same way, through the diagnostic of section 4 at the
loaded restart, with the grid's own `MOVED` statistic. MEASURED, `k = 1`:

| h | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|---|
| mean over the rows it moves | 1.239E-05 | 1.240E-04 | 1.238E-03 | 3.703E-03 | 1.238E-02 | 3.455E-02 | 1.190E-01 |
| rows counted as moved | 21 | 21 | 21 | 21 | 21 | 21 | 21 |
| `||J(h) - J(1)|| / ||J(1)||`, central | 4.766E-03 | 4.731E-03 | 4.316E-03 | 3.384E-03 | 0 | 1.058E-02 | 6.562E-02 |
| movement with `k` at that arc, `k = 8` | 3.99E-08 | 6.19E-09 | 4.69E-10 | 1.93E-10 | 3.49E-11 | 1.22E-11 | 5.02E-12 |

The first row is the suite's row to four digits at every arc, from a different
entry point, a different direction normalization and a different reference
residual. The third row is the Stage B statement the brief asks for: **the
CENTRAL difference along the failing direction is flat to 4.8E-03 over
`h` from 1E-03 to 1 and the production arc is inside that interval**, so the
derivative along the very direction the row fails on is resolved there. The
fourth row says the closure count moves it by at most 1.2E-05.

### 10.3 The cross-check with the tree binary

The arc ladder needs no new entry point, so it was repeated with the tree's own
`EXHALE.x` (md5 `89ec67149aa452eea3deeb19052f83d1`). MEASURED, energy
direction: 1.2381487E-03 at arc 1E-01 and 1.2378866E-02 at arc 1; carrier
direction: 2.3835770E-03 and 4.6313904E-04. **Digit for digit the measured
build's numbers**, so nothing in section 4 touches this result and it is a
property of the tree as it stands.

### 10.4 The verdict, and what the test should assert instead

**The test is measuring its own probe arc.** The quantity it gates is the
forward difference's first-order truncation at the production arc; it is
proportional to the arc over four decades on the energy direction and reaches
its rounding floor three decades below the production arc on the carrier one,
while the central difference it is measured against holds to a fifth of a per
cent across the same range. Inside the resolved interval the disagreement does
not survive: **1.2381E-03 at arc 1E-01, ten times inside the tolerance**. It
is not an implementation error, and the row as written cannot distinguish one
from a probe arc that is four decades longer than that direction's optimum.

This item does not edit the test. What it should say, in order of preference:

1. **Assert the first power.** Take the statistic at two arcs, the production
   one and a tenth of it, and require the ratio to be 10 within a stated
   factor, on each direction. That is a statement about the assembled action
   and is free of the arc: a forward difference that does not converge to the
   central difference as the arc is shortened is an implementation error, and
   one that does is not. It costs one extra invocation per direction, 17 s
   each MEASURED.
2. **Read the gate at an arc inside the resolved interval**, named in the
   test, with `EXHALE_JV_PROBE_ARC` set to it, and keep 1.0E-02. This must be
   done per direction and not globally: arc 1E-01 puts the energy direction at
   1.2381E-03 but the carrier direction at 2.3836E-03, still passing, while arc
   1E-03 puts the carrier direction at 0.293, at its floor, and would turn a
   passing row into a failing one for no physical reason.
3. **Gate the truncation constant**, `mean / arc`, which is 1.2379E-02 per unit
   arc for the energy direction and 4.6314E-04 for the carrier one. It is the
   second derivative of the residual along the direction divided by the first,
   and it is the quantity that actually drifted (see below).

Whichever is chosen, **the test should print the probe length it used**. It
already writes it into `jac.txt`; the README does not carry it, and that is
why the drift the brief names cannot be attributed.

**On that drift.** The README of the suite quotes 3.792E-03 for the energy mean
as of 2026-09-17 and 1.238E-02 for 2026-09-21, a factor 3.3, and 3.692E-04 to
4.631E-04 for the carrier one, a factor 1.25. Since the statistic is exactly
the first power of the arc, and the arc is `sqrt(epsilon)(1 + ||Y||)` times the
scale factor, a factor 3.3 is a factor 3.3 in either the arc or the curvature
along the direction. Today's arc is 6.868E-02 (energy) and 1.804E-06
(carrier); no earlier arc is recorded anywhere in the tree, so **the two cannot
be separated from what is stored**, and this document does not claim to have
explained the drift. Recording the arc in the README, or gating the constant of
item 3, would make the next such drift attributable.

## 11. What the coupled registration changes, against the same state uncoupled

Both grids were run on the same stored state with one binary, so the two are a
controlled pair. MEASURED:

| | uncoupled (`C0`) | coupled (`C1`) |
|---|---|---|
| system | 3 unknowns per cell, 1500 entries, 0 species rows | 5 unknowns per cell, 2500 entries, 2 species rows |
| `Dr`, the row scaling of the linear model | the state column scales | the certification row scales |
| `||Dj^-1 F||_2`, mass class | 1.9018E+00 | 1.9459E+00 |
| `||Dj^-1 F||_2`, energy class | 2.4573E-02 | 2.3277E-02 |
| largest judged mass row | 5.3329E-01 at cell 1 | 5.0734E-01 at cell 1 |
| largest judged species row | none | 8.1496E+01, the H2 carrier of cell 499 |
| `||D^-1 s||` of the banded model's proposal | 9.2545E-10 | 5.4880E+01 |
| the banded model against the full action along its own step | agrees to 3 digits | a factor of 40 apart |
| plateau of the localized mass direction | 1.27E-05 at `h = 1E-03` | 1.38E-02 at `h = 1E-03` |
| plateau of `smooth` | 2.63E-05 at `h = 1E-03`, 4.6E-08 at `h = 3` | 1.34E-06 at `h = 1E-03`, 6.70E-02 at `h = 3` |

Two things follow and they are separate.

- **The hydrodynamic residual of the SAME state is not the same map on the two
  routes**: the judged mass class moves from 1.9018 to 1.9459, 2.3 per cent,
  and the energy class from 2.4573E-02 to 2.3277E-02. That is expected, since
  the carrier content is a Newton unknown on one route and whatever the
  elimination returns on the other, and it is recorded so that no phase-7
  comparison reads a hydrodynamic row of the coupled run against an uncoupled
  number as though they were the same quantity.
- **The relative plateau measures are not comparable between the two columns**,
  because `Dr` is a different array; what IS comparable is the shape, and the
  shape says the base cells' mass rows limit the coupled grid where they did
  not limit the uncoupled one, at this state, with this row scaling.

## 12. The diagnostics of section 7 that become reachable

The brief notes that the Ritz values, the band difference, the Krylov size scan
and the binding row sit inside `trust_region_step` behind `elem_diag_here` and
become reachable once species rows are registered. **They were NOT taken**, and
the reason is not cost: they are produced inside a SOLVE, at outer iterations
1, 2 and `elem_diag_third` (`steady_newton.f90:17701`), and every measurement
in this document evaluates a frozen state and stops before any solve. Taking
them would mean running the coupled solve, which is phase 7 itself and is
Mode C. They belong to that phase, where they will be produced by the trial
run at no extra cost, and this document states that rather than leaving the
brief's item silently unanswered.

## 13. What was not excluded

- **Two states, and both have their elimination at its fixed point after one
  pass.** The `k` dimension was again exercised on a map that is already
  converged at `k = 1`. This remains the clearest limit of the closure half of
  the result, as it was in phase 4.
- **The cause of the fixed-size residue in the base cells' mass rows.** It is
  localized, it is not rounding, it is not the closure count, and a direction
  that avoids those cells does not carry it. A limiter, a positivity repair, a
  smoothness-weight switch and a mass row whose terms cancel to that level
  would all produce it, and nothing here distinguishes them. It is the same
  open item as section 10 of `docs/lhs1140b_closure_and_probe_20260921.md`.
- **The factor of 40 between the banded model and the full action** along the
  banded model's own proposed step on the coupled route. Measured, not
  explained. It concerns the preconditioner, not the action.
- **Thread independence of the `wellmixed/HeH0.55` grid**, which was run only at
  eight threads; the coupled grid of the phase-7 checkpoint was measured at one
  and at eight and is byte-identical (section 6).
- **`EXHALE_JV_COLUMN_SCALE` was left off**, so only the unscaled arc rule was
  swept. On the coupled route this is a live question and not a formality: the
  column-scaled rule displaces every unknown by the same fraction of its own
  scale, and it is the natural candidate for the `banded_step` direction, whose
  components span the conserved unknowns and the logarithm of a carrier. It
  would need its own grid.
- **The drift of the coupled-block statistic since 2026-09-17** (section 10.4).
- Processor affinity, other hosts, other library versions.

## 14. The closing verdict

**(a)** In the words section 8.2 sets, at the phase-7 checkpoint with the
carrier and the diffused element rows in the unknown set:

- **the coupled action stabilizes in `h`** on the species direction (2.8E-06
  over five decades of arc), on a smooth direction (1.4E-07 over `h` from
  1E-02 to 1) and on the energy row of the third cell (4.5E-04); it does not
  stabilize better than 5E-03 on the base cell's mass row, nor at all along the
  banded model's proposal;
- **it does not move with `k`**: at most 4.1E-06 at the production arc over the
  counts `k >= 2`, and the residual at `k = 1`, the count the solve forms at
  this iterate, is bitwise the converged-closure residual;
- **the inner map is resolved**, one pass to an increment of 1.461E-13 against
  the 1.0E-08 asked for, with every one of the 315 actions admissible, no
  missing chemical root, no refusal, no blocked component and no shortened arc.

So **the coupled action is resolved where the coupled system's own unknowns
live, and the limit is the base cells' mass rows, which were already the open
item of phases 4 and 5.** Phase 7 is not blocked by this measurement, and the
limit is stated so that its verdict is not read where the operator is only
accurate to a per cent.

**(b)** **The test is measuring its own floor, or rather its own arc.** The
disagreement is the first-order truncation of the forward difference at a probe
arc four decades above that direction's optimum; it falls exactly as the arc
and is 1.2381E-03 at arc 1E-01, ten times inside the tolerance, with the
central difference it is compared against unchanged to 0.16 per cent. It is not
an implementation error, and the row should assert either that the forward
difference converges to the central one with the first power of the arc, or the
truncation constant, rather than one number at one arc that is not recorded.
The test was not edited here.
