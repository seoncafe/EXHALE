# L4g: the pseudo-time ramp of the partitioned hydrodynamic solve

Item L4g of `docs/PLAN_20260913_lhs_stationary.md`, the item
`docs/lhs1140b_stationary_L4e_20260914.md` sections 7 and 10 end on. Every
number is MEASURED on this tree unless it is marked READ. The file changed is
`src/modules/time_step/steady_newton.f90` and nothing else.

## 1. Verdict

**The pseudo-time start `dtau0 = 1.0` that `Restart intent: stationary` fixes
is below the scale these states are solved at, and the ramp cannot leave it
because it is read from a functional that is flat on them and floored at the
start itself. Both are repaired, and the LHS 1140 b states that L4b, L4d, L4e
and L4f measured as insoluble now certify from that same start with no
continuation.**

`dtau` is in the code's time unit `R0/v0`, `v0 = sqrt(k_B T0/m_H)`, which for
this case is 8.07e+03 s (section 2). On the state the LHS 1140 b C/N/O case
enters its eighth pass at:

| | code units | seconds |
|---|---|---|
| the start `Restart intent: stationary` fixes | 1.0 | 8.07e+03 |
| the explicit-stable interval `eval_dt` returns (CFL 0.6) | 6.93e-05 | 0.56 |
| the largest cell crossing time in the column | 1.90e-01 | 1.54e+03 |
| the flow time `r/v` at the outer cell | 3.58e+01 | 2.89e+05 |

So `dtau0 = 1.0` is 1.4e+04 explicit-stable intervals -- fully implicit in
the acoustic modes -- and ONE THIRTY-SIXTH of the flow time through the
column: the shift `I/dtau` stands a factor 36 above the rate of the advective
modes, and those are the modes a stationary wind is made of. The iteration at
that pseudo-time freezes exactly what it has to move.

**That the ramp could not leave it is now measurable directly.** With the
pseudo-time printed beside every outer iteration (new, section 3), the entry
text on the L4d reproduction carries `dtau` from **1.00 to 1.42 over twelve
iterations**.

**The change.** The pseudo-time DOUBLES on an accepted step, which is the
growth the published practice allows a favorably converging step, and the cut
of a step the line search cannot admit stops at the explicit-stable interval
of the state instead of at `dtau0`.

**The ratio the item proposed -- the judged distance's fall -- was
implemented in five forms and none of them ramps** (section 5). The reason is
one measured fact: the first Newton step of these states raises the judged
distance at EVERY pseudo-time (1.03e+03 to 2.94e+06 at `dtau0 = 1.0`, and at
1e+08 to a continuity row of 3.3e-02 that the step after it removes and
certifies), so its rise is not a statement about `dtau`; read as an unclipped
ratio it takes `dtau` to 1.8e-04 on that one step, and read clipped it is
1.00 as soon as the distance is flat, which is what these states are.

**GREEN, at the delivered default `dtau0 = 1.0` and with no continuation:**

| state | entry text | this item |
|---|---|---|
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` pass-8 entry | STAGNATED, 263 iterations, `info = 2`, the state returned unchanged | **CERTIFIED at outer pass 1**, 18 Newton iterations |
| L4d's pass-2 entry of the same case | STAGNATED at 40; 98 iterations after L4d's detector (READ) | `info = 0` at 20 iterations |
| the RAW archived seed of `atomic_scalar_gj1132_kzz1e9/HeH2.13` | `info = 0` at 16 iterations, `\|\|R\|\|` 7.03e-08 | `info = 0` at 23 iterations, 1.36e-08, the same element row |

and the new case the item ends on, **`atomic_scalar_gj1132_kzz1e9/HeH4.0`
seeded from the certified HeH2.13 as L11 does, CERTIFIED at outer pass 5 in
one run**, `omega` never halved, where L11 records that case REFUSED at pass
4 with `omega` halved twice (section 8).

The recipe L4e had to write into `REPRODUCE.md` -- solve at `dtau0 = 1`, then
resume the written state at `dtau0 = 1e8` -- is no longer needed: one start
serves the raw seed and the near-certified state, which is what this item was
for.

## 2. The pseudo-time in the units the code keeps it in

`dtau` is in the code's time unit, `t_s = R0/v0` with
`v0 = sqrt(k_B T0/m_H)` (`input_read.f90`; the symbol `mu` there is the
hydrogen atom mass and not a mean molecular weight, READ). For
`atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`: `R0 = 0.157692 R_J = 1.1024e+09`
cm, `T0 = 226` K, so `v0 = 1.3655e+05` cm/s and **`t_s = 8.073e+03` s**.

The intervals of the column, on the state the eighth pass is entered at
(`output/Hydro_ioniz_IC.txt` of that state, cell widths from the same file):

| | cell | r [R_p] | seconds | code units |
|---|---|---|---|---|
| smallest `dr_j/(\|v_j\|+c_j)` | 50 | 1.0097 | 0.933 | 1.156e-04 |
| the same at the base | 1 | 1.0000 | 1.430 | 1.771e-04 |
| the same at the outer cell | 500 | 29.031 | 1.535e+03 | 1.902e-01 |
| flow time `r/v`, outer cell | 500 | 29.031 | 2.893e+05 | 3.584e+01 |

`eval_dt` returns `CFL min_j dr_j/(|v_j|+c_j)` with `CFL = 0.6` (the default;
this case states no `CFL:` key), which is **6.93e-05** in code units, and
that is the number the new floor takes -- confirmed against the value the cut
stops at in the runs of section 5.

This is also why `dtau0 = 1.0` was right where it was measured. It replaced
the CFL start for the HD 209458 b element reload and the hot-Uranus carrier
reload (READ, `docs/Update_EXHALE_stage2.md` section 8: those two stall from
the CFL interval and reach 1.2e-08 and 1.5e-09 from 1.0). Those reloads enter
far from their root, where a march is what is wanted. What is new here is a
route that enters NEAR the root, where the same number is four decades too
small -- and the repair has to serve both, which is section 7.3.

## 3. What changed in the source

All of it in `src/modules/time_step/steady_newton.f90`.

1. **The floor of the pseudo-time is the explicit-stable interval of the
   state the solve is entered at**, `eval_dt` of that state, and no longer
   `dtau0`. It is formed once, before the first iteration, because the floor
   is what the cut is allowed to reach and a floor that follows the iterate
   is a bound the iterate can move; it is never taken above the start, so a
   route that already starts inside that interval -- the direct steady route,
   which starts from this same number -- keeps the arithmetic it has always
   had; and a state whose ghosts the caller has not filled gives no usable
   interval, in which case the start is used. The comment at the old cut said
   the floor was "the explicit-stable dtau0", which for the stationary
   restart it was not; that statement is now true because the quantity is
   computed. The module gains one `use eval_time_step, only: eval_dt`.

2. **The three places that cut or gate on the floor read it**: the cut after
   a failed line search, the cut after a singular banded factorization, and
   the gate of the damped Gauss-Newton escape, which fires when the
   pseudo-time can no longer be reduced and therefore has to ask the floor
   and not the start.

3. **The ramp of an accepted step is a doubling.** The line-search merit's
   ratio and the `max(lam, 0.1)` factor no longer enter it. Section 5 is the
   measurement that chose this over the ratio the item proposed.

4. **The pseudo-time is printed** beside `||R||`, the merit and `lam` at
   every outer iteration of the JFNK route, as the direct pseudo-transient
   route has always printed it, and like that route it is the value the NEXT
   step is taken at. The field is a write of a variable the routine already
   holds; nothing is evaluated for it.

The file's diff against `43bc28c` moves from L4e's +589/-108 to +703/-115,
so this item's increment is +114/-7, of which the greater part is the two
comment blocks that state the measurements above.

One option, read once per solve by `read_species_unknown_space_controls`:
`pseudo_time_doubles_on_an_accepted_step`, **on by default**,
`EXHALE_PTC_RAMP_DOUBLE=0` restores the merit ratio and the floor at `dtau0`,
which is the arithmetic of every result before this item and is bit-identical
to it (section 7.1). The default is argued in section 7.5.

`solve_steady_ptc`, the direct pseudo-transient route, is NOT changed. Its
ramp has the same form, but it is entered from the CFL interval and is the
continuation it was designed as; nothing in this plan runs through it (the
marching hand-off and the stationary restart both enter `solve_steady_jfnk`),
and it is reported here rather than changed without a measurement of its own.

## 4. RED: the three reproductions, on the entry text

The control build is `43bc28c` plus the working tree's `steady_newton.f90`,
which is L4d's and L4e's change and is not in HEAD; the measured build
differs from it in that file alone. Both were built in an isolated copy of
the tree, because other workers were editing source files while this item ran
(section 9). 8 threads, `EXHALE_OUTER_PASSES=1`.

### (a) The pass-8 entry state of `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` at `dtau0 = 1.0`

No iteration cap:

```
(JFNK) STAGNATED: neither the judged distance nor the merit improved in 40
       iterations, monotone search included -- stopping at a judged distance
       1.993E+06, best 1.032E+03, merit over the window 4.772E-05,
       ||R||=  1.032E-03
(JFNK) done info=2 ||R||=  1.032E-03  non-monotone accepts=120
(EXHALE_main) outer pass 1: hydro info=2, worst gated species row 6.37E-06
       of 1.0E-05 at cell 283, mass 9.22E-09, momentum 1.51E-13,
       energy 1.03E-03, omega 0.500, 953.12 s
```

263 printed iterations; `||R||` 1.032e-03 in and out. Every number of that
pass line is L4e section 7's to every printed digit (its wall clock there was
783.61 s on a less contended machine, READ).

### (b) The same state at other starts, control build

| `dtau0` | outcome | Newton iterations | pass 1 seconds |
|---|---|---|---|
| 1.0 | STAGNATED, `info = 2`, the entry state returned | 263 | 953.1 |
| 1e1 | `info = 0`, `\|\|R\|\|` 2.824e-08, **outer pass 1 ACCEPTED** | 19 | 90.5 |
| 1e2 | `info = 0`, 3.654e-08, **ACCEPTED** | 8 | 37.6 |
| 1e4 | `info = 0`, 2.665e-08, **ACCEPTED** | 4 | 17.7 |
| 1e8 | `info = 0`, 6.972e-08 (READ, L4e section 5) | 5 | 22.8 |

One decade is the whole difference between a solve that returns its entry
state after 263 iterations and one that certifies in 19.

### (c) The raw archived seed of `atomic_scalar_gj1132_kzz1e9/HeH2.13`

`models/.stopped/HeH2.13_archive_seed_output`, the state
`map_state_to_grid.py` wrote from `archive_20260830`, with that case's
`input.inp`. Its judged distance at entry is 3.472e+11, eight decades above
the pass-8 state's 1.032e+03.

| `dtau0` | outcome |
|---|---|
| 1.0 | `info = 0` at **16** iterations, `\|\|R\|\|` 7.03e-08, 107.9 s; the outer loop goes on because the element row is 2.37e-04 at cell 319 |
| 1e8 | `info = 1` at the 25-iteration cap, `\|\|R\|\|` 1.653, mass row 1.65, judged distance 1.47e+08, 374.6 s, the damped Gauss-Newton escape firing at EVERY iteration because `dtau` sits on its floor, and its floor is `dtau0` |

### (d) The entry text's ramp, measured

With the pseudo-time printed (the option off, so the arithmetic is the entry
text's), the L4d pass-2 reproduction at `dtau0 = 1.0`, twelve iterations:

```
dtau = 1.06  1.08  1.12  1.17  1.22  1.27  1.31  1.34  1.36  1.38  1.40  1.42
```

That is the obstruction in one line.

## 5. The ratio the item proposed, in five forms, and why none of them ramps

The item's design was "`dtau` grows with the judged distance's fall and
shrinks with its rise". Five forms of that were implemented and MEASURED on
reproduction (a) at the delivered start `dtau0 = 1.0`, each with the floor
already lowered to the explicit-stable interval. `jref` below is
`maxval(jhist)`, the same five-iterate window the line search accepted the
step against, and `dj_try` the judged distance of the accepted trial.

| form | what `dtau` does | outcome |
|---|---|---|
| `max(lam,0.1) * jref/dj_try`, unclipped | 1.0 -> **1.76e-04 on the first step** | the solve stands still at 1.8e-04 |
| the same clipped to [0.1, 2], the published bounds | 5.0e-02, 2.9e-02, 6.8e-03, ... **6.93e-05 at iteration 8, 7.18e-05 at iteration 26** | the ratio is 1.00 once the distance is flat; `dtau` sits on the floor for the rest of 80 iterations |
| the window's verdict, x2 inside it and x0.5 outside, times `max(lam,0.1)` | 2.5e-01, 6.3e-02, 1.6e-02, 7.8e-03, 3.9e-03, then the floor | the first five accepted steps take the distance above the window (1.03e+03, 2.94e+06, 4.41e+06, 5.14e+06, 5.50e+06, 5.68e+06), so every one of them halves `dtau`; the line search then finds nothing and the cut takes it the rest of the way |
| the same verdict without the `lam` factor | 5.0e-01, 2.5e-01, 1.3e-01, 6.3e-02, 3.1e-02, then the floor | the same |
| the entry text while the merit falls, `max(2, jref/dj_try)` when it does not | 5.6e-01, 1.1, 2.2, ... **5.09 at iteration 26** | once `dtau` is a few the merit falls by a few percent on nearly every iteration, so the entry-text branch governs and the ramp creeps |

**The reason is one measured fact: the first Newton step of these states
raises the judged distance at every pseudo-time.** At `dtau0 = 1.0` it goes
1.032e+03 -> 2.937e+06; at 1e+08 the same step raises the continuity row to
3.29e-02 (READ, L4e section 5) and the step after it removes the excursion
and certifies. A rise of the judged distance is therefore not a statement
about `dtau`, and a rule that reads it as one cuts the pseudo-time on exactly
the step it should be raised over. Once `dtau` is cut, the state stops moving
and the distance goes flat, and a ratio of a flat functional is one.

A sixth form was measured and is the one NOT delivered for a different
reason: `max(2, jref/dj_try)` on every accepted step, no `lam`. It certifies
reproduction (a) in 12 iterations and (c) in 15, better than what is
delivered -- and it drives the HD 209458 b element reload, which needs the
transient, to `||R||` 1.53 at iteration 74 of a first pass that reaches
2.2e-09 with the plain doubling, and its second pass to a mass row of 2.2e-01
against 7.5e-10. The judged distance of that fixture falls fast early, the
ratio then reaches 254 and 3.6e+04 in single steps, and the pseudo-time is
four decades past the state in six iterations. **The published clip exists
for that**, and taking it as the rule is what serves both routes.

## 6. GREEN: the same states, one start, the delivered defaults

Measured build, `EXHALE_PTC_DTAU0=1.0` (the value `Restart intent:
stationary` fixes; it is stated in the command only because the RED runs
state it), `EXHALE_OUTER_PASSES=1`, 8 threads.

| state | entry text | this item | the pseudo-time it ran at |
|---|---|---|---|
| (a) pass-8 entry, `dtau0 = 1.0` | STAGNATED, 263 iterations, `info = 2`, 953.1 s | **CERTIFIED at outer pass 1**, `info = 0` at **18** iterations, `\|\|R\|\|` 4.110e-08 | 2, 4, 8, ... 2.62e+05 |
| (b) the same at `dtau0 = 1e8` | `info = 0` at 5 iterations, 22.8 s (READ) | **CERTIFIED at outer pass 1**, `info = 0` at **5**, 4.386e-08, 23.8 s | 2e+08 ... 3.2e+09 |
| (b') L4d's pass-2 entry, `dtau0 = 1.0` | STAGNATED at 40; 98 iterations after L4d's detector (READ) | `info = 0` at **20** iterations, 4.964e-08, 167.8 s | 2 ... 1.05e+06 |
| (c) the raw archived seed, `dtau0 = 1.0` | `info = 0` at 16, `\|\|R\|\|` 7.03e-08, 107.9 s | `info = 0` at **23**, **1.356e-08**, 243.9 s | 2 ... 8.39e+06 |

(a)'s certified pass reads `worst gated species row 6.20E-06 of 1.0E-05 at
cell 217 (elemental transport He/H partition), mass 1.40E-09, momentum
5.43E-13, energy 4.11E-08`, then `ACCEPTED -- every active equation of this
state is within its own tolerance`. It is the same solution L4e reached by
resuming that state at `dtau0 = 1e8`, now from the start the restart contract
fixes.

(c) is the check that the ramp does not replace one wrong start with another:
the raw seed, eight decades further from certification, reaches the SAME
place -- the element row 2.37e-04 at cell 319 to three digits -- at a smaller
`||R||` than the entry text's, in 23 iterations against 16. That is the price
of the ramp on a state the old start was already right for, and it is paid in
iterations and not in the answer.

**What is still NOT a recipe: a raw seed at `dtau0 = 1e8`.** With the floor
lowered the cut now works there (`dtau` 1.28e+10 -> 3.20e+09 at the first
iteration that finds no admissible step) but 1e8 is fourteen cuts above the
scale that state wants, and the solve ends `info = 1` at the 25-iteration cap
at `||R||` 1.940, as it does on the entry text (1.653). The delivered start
does not need it.

## 7. Controls

### 7.1 With the option off, the binary is the entry text, bitwise

The L4d pass-2 reproduction, `EXHALE_PTC_RAMP_DOUBLE=0`, twelve iterations,
measured build against control build: `output/Hydro_ioniz.txt` and
`output/Ion_species.txt` **IDENTICAL byte for byte** apart from the
provenance line's timestamp, and every printed iteration of the log the same
to every digit apart from the new `dtau` field.

### 7.2 `backup/regression/carrier_model_a_newton`, the four-unknown route

`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40 OMP_NUM_THREADS=8 ./run.sh`, on
scratch copies, control build against measured build with the delivered
defaults.

**Same verdict and the same two refusing entries to every digit**: `carrier
balance H2: gated row measure 6.996E-04 above 1.0E-05 at cell 363` and
`carrier balance H+: 7.425E-04 at cell 323`; the worst gated species row of
passes 11 and 12 identical (1.68E-03 at cell 248, 7.42E-04 at cell 323);
every pass `hydro info = 0` in both.

| file | max relative movement |
|---|---|
| `output/Hydro_ioniz.txt` | 1.768e-10 (row 5, column 2) |
| `output/Ion_species.txt` | 1.635e-10 (row 210, column 1) |

Seven decades inside the harness's 1e-3. No golden of this case moves and
none was refreshed.

### 7.3 `backup/regression/atomic_elem_newton`, the element reload

HD 209458 b, atomic with solar trace metals and binary element diffusion,
reloaded from `IC/` with `Load IC? True`, `Coupled carrier solve: False`,
`Restart intent: stationary`; `EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80
EXHALE_DIFF_OMEGA=0.5`, 8 threads, the README's own recipe. This is the
fixture the LHS 1140 b repair had to be measured against, because `dtau0 =
1.0` was adopted FOR it.

**CERTIFIED in both**: the control at outer pass 11, the measured build at
outer pass 12.

| pass | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| worst gated species row, control | 2.91e-3 | 4.21e-3 | 1.77e-3 | 9.42e-4 | 4.66e-4 | 2.42e-4 | 1.23e-4 | 6.30e-5 | 3.22e-5 | 1.64e-5 | 8.42e-6 | -- |
| the same, measured | 2.61e-3 | 3.81e-3 | 1.76e-3 | 9.31e-4 | 4.61e-4 | 2.31e-4 | 1.19e-4 | 5.76e-5 | 3.11e-5 | 1.41e-5 | 1.04e-5 | 7.73e-6 |

The elemental relaxation is the same relaxation: the rows agree to two digits
at every pass. What the item buys this fixture is its FIRST pass, the
hydrodynamic solve the old ramp could not finish there either:

| pass 1 | mass | momentum | energy |
|---|---|---|---|
| control | 2.91e-02 | 4.73e-05 | 4.11e-01 |
| measured | **2.16e-09** | **4.24e-13** | **1.55e-08** |

and what it costs is one outer pass, because the element relaxation halves
`omega` to 0.250 at pass 2 where the control holds 0.500. The two very nearly
cancel: 11 passes against 12, and the same verdict.

**Annotation of 2026-09-16 (item L16).** The pass-1 rows of the two tables
above were measured at 8 threads on the builds of 2026-09-14, whose object
directories were deleted and for which no manifest exists (the first binary
manifest is `db87b88d1ce5` of 2026-09-15), so they cannot be re-measured on
the build that produced them. The same recipe on the tree of 2026-09-16
(binary `c2e9c9990b9f`, item L15 stage 2, single-threaded and 8-threaded
alike) ends pass 1 at mass 4.57e-10, momentum 3.48e-14, energy 2.72e-08, and
the control of item L4h measured 3.08e-02 / 1.90e-05 / 8.34e-01 on its own
build of 2026-09-15. The three sets are three builds, and the difference
between them is attributed to nothing here: the source moved between them
(items L7e, L13, L14, L18, L19, L21 and the R4 Jacobian corrections), and
item L15 stage 2 found no thread dependence on the 2026-09-16 tree, so the
8-thread measurement is not by itself the cause. The rows stay as the
record of what those builds gave.

### 7.4 The test suites

Run from the isolated tree against its private object directory, measured
build.

| suite | outcome |
|---|---|
| `krylov_and_dogleg` | all rows passed |
| `steady_species_rows` | 195 PASS, **FAIL 0**, the reload rows included (`EXHALE_SPECIES_EXE` set to the measured binary) |
| `certification` | every assertion passed |
| `attempted_step` | every assertion passed |
| `carrier_retry` | every assertion passed |

`steady_species_rows` and `steady_selfconsistent_residual` are the two suites
whose drivers link `steady_newton.o`; the other three were run because the
brief names them. `steady_selfconsistent_residual` needs a run log
(`EXHALE_STEADY_RUNLOGS`) and was not run.

### 7.5 Why the default is ON

The option reaches every stationary restart and the marching hand-off, so the
question is not whether it helps LHS 1140 b but whether it damages anything
that works. What was measured:

- the four-unknown route, which the ramp reaches on every solve: same
  verdict, same refusing rows, movement 1.8e-10 (7.2);
- the element reload, the fixture `dtau0 = 1.0` was adopted for: the same
  certification, one pass later, with a first pass seven decades better
  (7.3);
- the five suites: no row moves (7.4);
- with the option off, bit-identical to the entry text on both routes (7.1).

Against that, the route it repairs is the one the whole of
`docs/PLAN_20260913_lhs_stationary.md` is about, and the repair is not a
speed-up: three states that could not be solved at all are solved, and the
case ladder of L11 closes (section 8). An option left off would have to be
set by hand in every campaign run, and the campaign is what this plan
delivers. **It is on, and the movement it causes is reported above and in
7.6.**

## 8. The case the item ends on: `atomic_scalar_gj1132_kzz1e9/HeH4.0`

The L11 ladder step: the case seeded from the CERTIFIED `HeH2.13` solution of
the same group, mapped onto this grid by `map_state_to_grid.py
--reservoir-HeH` (the pair `models/atomic_scalar_gj1132_kzz1e9/HeH4.0/output/
*_IC.txt` that case's `seed.log` records writing), run from that case's own
`input.inp` in one go with no continuation, `EXHALE_PTC_DTAU0=1.0`, 8
threads, the outer pass budget left at its default.

**CERTIFIED at outer pass 5**, `omega` never halved:

| pass | hydro `info` | worst gated species row | mass | momentum | energy | `omega` | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.93e-05 at cell 292 | 8.21e-10 | 1.47e-14 | 1.98e-08 | 0.500 | 179.2 |
| 2 | 0 | 4.71e-05 at cell 294 | 9.13e-10 | 7.97e-14 | 1.37e-08 | 0.500 | 97.3 |
| 3 | 0 | 2.48e-05 at cell 296 | 7.11e-10 | 5.93e-14 | 1.59e-08 | 0.500 | 81.1 |
| 4 | 0 | 1.30e-05 at cell 297 | 2.02e-09 | 2.18e-12 | 3.51e-08 | 0.500 | 76.5 |
| 5 | 0 | 6.79e-06 at cell 299 | 8.02e-10 | 9.51e-13 | 2.54e-08 | 0.500 | 75.9 |

`ACCEPTED -- every active equation of this state is within its own
tolerance`, `CERTIFIED`, flux gate accepted 4.7044e-11, 509 s in all.

L11 records the same case from the same seed REFUSED at pass 4, "the joint
distance of the state has not fallen in 3 consecutive passes", with `omega`
halved twice and hydro `info = 1` from pass 2 on (READ, plan row L11 and
`docs/lhs1140b_stationary_L4e_20260914.md` section 3 (c), which reproduces
it: passes 1 to 4 at 8.93e-05, 4.80e-05, 3.68e-05, 3.25e-05 and `omega` 0.500
to 0.125). Pass 1 is identical in the two runs -- 8.93e-05, the seed is in
the basin either way -- and what differs from pass 2 on is that the
hydrodynamic solve now reaches its tolerance, so the joint distance falls and
the relaxation keeps its step.

This is the first case of the L11 ladder to certify from its seed in one run.

## 9. Reproduction

The measurements were made from an ISOLATED copy of the tree, `git archive
HEAD` (`43bc28c`) plus the working tree's `steady_newton.f90`, for the reason
L4e gives: other workers were editing source files while this item ran, and a
binary built from the live tree would carry unfinished physics into every
number above. The control build is that tree with L4e's `steady_newton.f90`;
the measured build is the same tree with this item's.

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
T=$EX/scratchpad/L4g/tree
mkdir -p $T && cd $EX && git archive HEAD | tar -x -C $T
\cp -f $EX/src/modules/time_step/steady_newton.f90 $T/src/modules/time_step/
cd $T && PATH=/usr/bin:$PATH make OBJDIR=build_L4g EXE=EXHALE_L4g.x -j8
# the delivered source text builds to md5 d9178c554ebdb1999269e54a2a6ea5c4
# here, and reproduces section 6 (a): info = 0 at 18 iterations, ||R||
# 4.110E-08, outer pass 1 ACCEPTED

# (a) and the GREEN of section 6: the pass-8 entry state of
#     atomic_scalarCNO_gj1132_kzz1e9/HeH2.13 (input.inp, base.inp and
#     output/*_IC.txt are $EX/scratchpad/L4e/red_A's)
cd <case> && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 \
             EXHALE_JFNK_MAXIT=80 <binary>

# (b') the pass-2 entry state: $EX/scratchpad/L4e/red_B
# (c)  the raw seed: models/atomic_scalar_gj1132_kzz1e9/HeH2.13/input.inp
#      with models/.stopped/HeH2.13_archive_seed_output/*_IC.txt in output/
#      (the Spectrum file: line of that input.inp is relative and has to be
#      made absolute when the case is run outside models/)

#   EXHALE_PTC_RAMP_DOUBLE=0   the merit ratio and the floor at dtau0, i.e.
#                              the entry text, bitwise

# the fixtures, on copies of the case directories
#   atomic_elem_newton: Load IC? True, Coupled carrier solve: False,
#   Restart intent: stationary, IC/*_IC.txt into output/, then
#   EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80 EXHALE_DIFF_OMEGA=0.5
#   OMP_NUM_THREADS=8 ./run.sh <binary>
#   carrier_model_a_newton: EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40
#   OMP_NUM_THREADS=8 ./run.sh <binary>
# the suites, against the private object directory
#   EXHALE_OBJDIR=$T/build_L4g EXHALE_TEST_OUT=$T/build/tests/<suite> \
#     [EXHALE_SPECIES_EXE=$T/EXHALE_L4g.x] $T/src/tests/<suite>/run.sh
# section 8: models/atomic_scalar_gj1132_kzz1e9/HeH4.0's input.inp and its
#   output/*_IC.txt pair, copied to a scratch directory, OMP_NUM_THREADS=8
```

> Note added 2026-09-14 15:20: the private builds of this memo were made with `PATH=/usr/bin:$PATH make`, which selects the Ubuntu gfortran 13.1 and its LAPACK, not the conda-forge gfortran 16.2 + OpenBLAS the tree binary `EXHALE.x` is built with (bare `make`). The control and the measured build of the memo share one toolchain, so every comparison here stands; reproduce with bare `make` (`OBJDIR=... EXE=...` as written) to obtain the tree's toolchain.

## 10. Scope, and what was not measured

- **What was verified is the path the change touches**: the stationary solve
  of the three-unknown (partitioned) route on the three reproductions and the
  new case, the four-unknown route on the carrier reload, the element reload,
  the two regression cases that reach `solve_steady_jfnk` through the
  marching hand-off (7.6), and the five suites. The trust-region branch takes
  its step from `tr_delta` and not from `dtau`, and the only line of it this
  change touches is the restart that sets `dtau = dtau0`, which is unchanged.
- **No golden was refreshed.** The movement measured on the two cases that
  carry goldens is in 7.2 and 7.6.
- **`solve_steady_ptc` keeps its ramp and its floor** (section 3). It is
  reached only by the direct steady route with `EXHALE_STEADY_JFNK` unset, no
  fixture runs it, and it is entered from the CFL interval where the ratio's
  premise holds.
- **The wall times are contended**: other workers held the machine
  throughout, and the same run took 131.8 s and 83.8 s on two occasions. Only
  iteration counts and named outcomes are claims.
- **`EXHALE_PTC_DTAU0=1e8` on a raw seed is still not a route** (section 6).
  The floor makes the cut able to descend, but not in 25 iterations from 1e8.

## 11. Noticed outside this item, reported and not acted on

- **`backup/regression/` carries eleven case directories whose names use the
  noun "arm"** (`armA_LW`, `armA_noLW`, `armD_D2`, `armD_D2_LW`,
  `armD_D2_LW_newton`, `armD_D2_newton`, `armD_D2_newton_bigstack`,
  `armHeH_0p3`, `armHeH_3`, `armHeH_10`, `armHeH_30`, `arm_heh1_x2matched`).
  The 2026-09-11 sweep took that noun out of every `.md`, `.tex` and `.f90`;
  these are directory names a rename would have to follow through
  `run_check.sh`'s case lists and every memo that cites them, so they are
  reported here rather than moved.
  **DONE 2026-09-16** (PLAN_20260916_rev3 section 10): all sixteen were
  renamed; the old-to-new mapping is `docs/named_case_audit.md` section 6.
- **`solve_steady_ptc`'s cut comment** says the floor is "the explicit-stable
  dtau0", the statement this item corrected on the JFNK route. On that route
  the start IS the CFL interval, so the comment is true there; it is noted
  because the two routes now say the same thing for different reasons.
### 7.6 The two regression cases that reach this code path through the marching hand-off

`wasp_full_newton` and `mol_base_handoff`, run by the harness on scratch
copies with `REGRESSION_EXE` pointing at the measured build and
`REGRESSION_GOLDEN_DIR` at the tree's `golden/`, single-threaded as the
harness fixes it. Both reach `solve_steady_jfnk` through the marching
hand-off, so the ramp acts on their Newton finish.

| case | file | verdict |
|---|---|---|
| `wasp_full_newton` | `Hydro_ioniz.txt` | PASS within 1e-3, max rel 1.149e-08 |
| | `Ion_species.txt` | PASS within 1e-3, max rel 2.171e-08 |
| | `Hydro_ioniz_adv.txt` | **FAIL**, max rel 1.000e+00 at row 62 column 10, 0 against 3.724e-14 |
| | `Ion_species_adv.txt` | PASS within 1e-3, max rel 2.171e-08 |
| `mol_base_handoff` | `Hydro_ioniz.txt` | PASS within 1e-3, max rel 8.703e-09 |
| | `Ion_species.txt` | PASS within 1e-3, max rel 7.231e-08 |
| | `Hydro_ioniz_adv.txt` | PASS within 1e-3, max rel 6.005e-07 |
| | `Ion_species_adv.txt` | PASS within 1e-3, max rel 7.231e-08 |

`wasp_full_newton`'s Newton finish is a real solve and the ramp runs in it:
`||R||` 6.194e-09, `info = 0` at ten iterations with `dtau` 2, 4, ... 1.02e+03.
Its two solved files stay within 1e-08 of the golden, and L4e measured the
CONTROL build of this same tree at 1.321e-08 and 2.512e-08 against the same
golden (READ), so the movement this item adds is of the size the tree already
carries and no golden moves at 1e-3.

**The one FAIL is the one L4e reported on the control build**: a golden zero
compared with a quantity at the rounding floor, which is a relative measure
of one by construction. L4e read it at row 259 column 10 (0 against 4.001e-14)
and it is here at row 62 column 10 (0 against 3.724e-14); the cell it lands in
moves with the last bits of the state, the nature of it does not. It is the
state of the tree and not something this item did, and it is reported and not
acted on.

`mol_base_handoff` is a relaxation snapshot pinned at 12000 steps, so its
hand-off is a fixed-step stop and not a converged solve; its movement is the
tree's own against the golden and is the same figure L4e measured on the
control build (8.7e-09 to 7.2e-08, READ).

