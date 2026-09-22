# LHS 1140 b, P2 at the base: Stage C and Stage D

Phase 9 of `docs/PLAN_20260920_rev9.md` section 15, restricted as section 17.3
restricts it: the directions probed are the base ones, and the rows attributed
are the mass row at cells 1 and 2 and the energy row at cell 1. Written
2026-09-21.

## 1. The result, first

**At the base rows of this checkpoint, under the current binary, none of the
four readings of the decision table in section 6 holds, because there is
nothing left there to classify: the base rows are solved.** From the stored
generation `g0002_20260919T004843Z_f7485b14` the default three-unknown route
takes the mass row of cell 2 from 9.806E-01 to 2.439E-08, the energy row of
cell 1 from 1.071E+00 to 2.559E-07 and the momentum row of cell 1 from
8.724E-06 to 1.017E-14 in 25 accepted Newton steps of one outer pass, ending
`info = 0` with every certified hydrodynamic row inside its own tolerance
(MEASURED, observation group, `obs_run.log:395`). The line search never cut a
step in that solve (`lam = 1.00E+00` at all 25), no unknown was ever held at a
bound (`nact = 0` at every iteration of the trace, and `[diag 5] unknowns held
at a bound 0`), and the requested linear tolerance was reached in every one of
the 25 cycles, in one to four operator products.

**The base face mass flux, which is the quantity section 17.1 built its
argument on, comes to the wind value in that one pass.** MEASURED from the
conservation export of section 13.2, evaluated on the entry state and on the
state the pass returned:

| face, in units of the wind window's own r^2 (rho v) | entry `g0002` | after outer pass 1 |
|---|---|---|
| base face (lower face of cell 1) | 82.1824 | 1.0000 |
| upper face of cell 1 | 51.5850 | 1.0000 |
| upper face of cell 2 | 1.0000 | 1.0000 |

So the 82.18 of section 17.1 and the 82.2 of the boundary memo of 2026-09-19
are properties of that stored state, not of the base law the current binary
solves. **This contradicts the framing of section 17.1** ("the base state
these cases are given is not compatible with a stationary wind at that base"),
and it decides between the two readings section 17.2 left open in favor of the
second: at this checkpoint the base is far from its own fixed point and the
path to it is not hard.

**What does refuse, after the base is solved, is not at the base.** The pass
that reaches `info = 0` leaves the elemental transport He/H partition at
7.218E-05 against 1.0E-05 at cell 295, a wind cell; the composition update of
the alternation then raises the hydrodynamic energy row from 2.559E-07 back to
6.512E-06, and the second hydrodynamic solve stagnates with its worst row at
cell 266 to 268 (r = 1.47 to 1.49) and its worst certified entry at cell 196
(r = 1.1388). Both groups end there, in the same place, from different step
controls.

## 2. Classification

| activity | label | what was performed |
|---|---|---|
| frozen evaluation of `g0002` under the current binary | Mode R | `EXHALE_RESIDUAL=1`, one assembly, no step |
| frozen evaluation of the state outer pass 1 returned | Mode R | the same, on the state the observation group wrote |
| bounded solve, default route, `EXHALE_LINEAR_ROWS=1` | Mode C | observation group, budget in section 4 |
| bounded solve, `EXHALE_TRUST_REGION=1` and the element diagnostic | Mode C | named numerical variant, its own step control AND its own merit |

Nothing was published into `LHS1140b/models/`. Every run stood in its own
directory under the session scratch, with a copy of the case `input.inp`, the
stored state pair under the names `load_IC` reads, a copy of the executable
and an empty `output/`. The two Mode C runs wrote `output/` inside their own
run directories only.

The checkpoint was written under boundary model
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2` and the binary
solves `..._ghost_fixed_point_seed_reservoir_row_v3`, so every run here is the
CURRENT OPERATOR ON AN OLD STATE; `load_IC` states it and rebuilds the
boundary (`obs_run.log:31-34`).

## 3. Identity

### 3.1 The binary and the tree

| | |
|---|---|
| path | `EXHALE.x` of the tree, the binary of record; copied into each run directory |
| md5 MEASURED | `89ec67149aa452eea3deeb19052f83d1` |
| repository HEAD | `93eed8667483e11dffe51b46008d40081292e104`, working tree dirty (MEASURED, 59 entries) |
| source digest | `a2f0da209aec8243f70caf6d7e339c5d`, the md5 of the sorted md5 list of `Makefile` and the 246 `.f90`/`.inc` files of `src/` (MEASURED 2026-09-21) |
| compiler strings | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`; `GCC: (conda-forge gcc 16.2.0-4) 16.2.0`; `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` (MEASURED, `strings -a`) |
| linear algebra | OpenBLAS (`libopenblas.so.0`), thread pool pinned to 1 by `OPENBLAS_NUM_THREADS=1` in every run and echoed by the run |
| threads | `OMP_NUM_THREADS=8` in the two Mode C runs, `1` in the Mode R evaluations; each run states its own count |

### 3.2 The checkpoint

| | |
|---|---|
| case | `LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` |
| generation | `g0002_20260919T004843Z_f7485b14` |
| Stage A0 record | `docs/lhs1140b_identity_records_20260921.md` section 3.1 |
| `Hydro_ioniz.txt` md5 MEASURED | `25e11562ec3cd1f6d204bf36f8e503eb` |
| `Ion_species.txt` md5 MEASURED | `78ac49d7e1ce90a3a4f271785400efda` |
| `input.inp` md5 MEASURED | `e43c6697786711fb916d3b8eb3c113aa` |
| `base.inp` | the case uses none |

All three md5 values equal the ones the generation manifest records and the
ones the conservation audit of 2026-09-21 MEASURED.

### 3.3 The effective parsed configuration

READ from each run's own `EXHALE_setup.out` and `EXHALE_resolved.out` and from
the run log: `Grid type: Mixed`, 500 cells, `Base grid [dr,cells]:
1.9999999494757503e-4 50`, `Numerical flux: ROE`, `Reconstruction scheme:
PLM`, `Solver: Newton`, `Restart intent: stationary`, `Well balanced: True`,
`He_diffusion: True` with `He_Kzz: 1.0e9`, `Secondary_ionization: Immediate`,
`Resid tol: 5.0e-5`, base BC pressure-anchored at 1.0 microbar with derived
`n0 = 3.2049E+13 cm^-3`. No `metals.inp`, no `opacity.inp`. `EXHALE_RESID_QUAD`
was left unset in every run, so the assembly selector of section 13.6 stands at
its default; the export header of the Mode R runs states the assembly actually
used, `rows_kind=production_double`, `momentum_branch=weno3_well_balanced`,
`recon=WENO3`.

The reconstruction key says PLM and the assembly says WENO3: the run is inside
the reconstruction continuation, which the export header also states
(`rec_method WENO3 (the flag pair does not name the branch while the
continuation is armed; the branch above does)`). This is READ, not a defect
found here; it is the same statement the conservation audit recorded.

### 3.4 `Energy solver` and the marching map

Section 17.3 requires `Energy solver: Explicit` for any comparison of the
marching update with the stationary residual at the base. **No comparison of
that kind is made in this document.** Every run here carries `Restart intent:
stationary`, which measures the loaded state and then solves; the run states
it (`the loaded state is measured as it stands; no CFL step is taken`), so the
marching operator and the semi-implicit energy step are never entered and the
requirement does not bind. This is recorded rather than assumed.

## 4. The declared budgets and the recorded endings

Estimated, as section 15 asks, from one MEASURED frozen evaluation on the
actual checkpoint: 0.71 s of one thread for one residual assembly of the
loaded state (`/usr/bin/time -v`, exit 0).

| group | declared budget | ending, RECORDED |
|---|---|---|
| observation | 60 Newton iterations per inner solve, 2 outer passes | pass 1 inner solve returned `info = 0` at iteration 26 of 60 (loop-top stop, all certified rows within tolerance); pass 2 inner solve returned `info = 2`, STAGNATED at iteration 60 of 60 after the monotone restart; the outer loop "spent its budget of 2 passes without an accepted state: 2 active equations refuse it" |
| named variant | the same 60 and 2 | pass 1 `info = 0` at iteration 25 of 60; pass 2 `info = 1`, cap reached; the same outer ending |

Spent, MEASURED: observation group 2 m 25 s wall at 8 threads (pass 1 53.29 s,
pass 2 91.52 s), 600 residual evaluations in pass 1 and 2332 in pass 2, 25 and
60 Krylov cycles. Named variant 8 m 24 s wall (pass 1 159.87 s, pass 2
347.25 s), 1913 and 8508 residual evaluations, 24 and 60 cycles. No run was
stopped from outside and no timeout was converted into an ending.

Two supporting runs: one inner iteration with the same keys and one outer
pass, to place the array dump at the entry state of pass 1 (3 s); and one
solve of 60 iterations with `EXHALE_OUTER_PASSES=1` and the trace off, which
reproduced the observation group's pass 1 iteration for iteration (same
`||R||`, `||Fs||2`, `lam`, `dtau` and cycle length at all 26 lines). That is
the control that the trace changes nothing.

## 5. Stage D, the ledger

### 5.1 What each field is, and where it comes from

| field of section 7.2 | observation group | named variant |
|---|---|---|
| state identity | the hash record of `EXHALE_L15_TRACE`: `Y`, the composition, `F`, `F_jac`, `D`, `Drow`, the bound mask and `dtau`, once per iteration | the same key was set; the same record |
| active unknowns | 1500, the three hydrodynamic unknowns of 500 cells; `set_transported_species_rows` is gated on `carrier_in_newton`, false here | the same |
| unscaled row measures | in the run's output, `[diag] residual, both norms`, at the entry and at the state handed back only | the same |
| tolerance-normalized distances with their cells | in the run's output at every iteration: `judged rows: distance ... hydrodynamic ... element ... carrier` and `cells outside the tolerance of their row` | the same |
| actual scaled displacement | NOT IN THE RUN'S OUTPUT on this route; needs a diagnostic capture | in the output: `||s||` on every `[TR] it` line |
| step length or trust radius | `lam` (line search) and `dtau` (the pseudo-transient shift) on every iteration line | `delta` on every `[TR] it` line |
| predicted and observed reduction under fixed weights | NOT DIRECTLY IN THE OUTPUT; see section 6.2 for what is available and what it is not | in the output: `pred`, `actual`, `ratio` on every `[TR] it` line |
| reduced Krylov residual | `the cycle reached <rel> in <n> product(s) of <m>` | `[diag 4] Krylov leg: ... relative residual returned ... TRUE relative residual of the returned step` |
| independently evaluated residual of the returned step | `linear residual left, over the row of F`, from a fresh action of the Jacobian on the returned direction | the same, plus the TRUE relative residual of `[diag 4]` |
| bound hits with their projection displacement | NOT APPLICABLE: there are no species unknowns, so no box and no bound. `nact = 0` in every trace record, and `[diag 5] unknowns held at a bound 0` | NOT APPLICABLE, same evidence |
| limiter and contact events | the positivity-limited face count is not printed at each iteration; NOT IN THE RUN'S OUTPUT | the same |
| change of composition or radiation closure | in the output at the end of each solve: `composition elimination: N sweeps over M residual evaluations, k per evaluation; the Newton model carried p Picard term(s)` | the same |
| termination status | distinguished in the output: inner `info` and its named reason; outer `spent its budget`; no external interruption in any run | the same |

Where a field is not in the output it is named so above; none of them is
recorded as a zero.

### 5.2 The base state of the scalings and the arrays, Stage C's precondition

MEASURED from the binary dump the trace writes at the entry of pass 1
(`dump_entry_summary.txt`), 1500 unknowns, `nvar_jac = 3`, `dtau = 1.0`:

| array | md5 of the array bytes | min | max | 2-norm |
|---|---|---|---|---|
| `Y` | `a2f5b5e4872b7d91f4ae3203f6da07fe` | -1.657008E-05 | 3.489290E+00 | 9.167074E+00 |
| `F` | `8fc4647d5f009a7bd4faa887d5efed80` | -1.345809E-03 | 3.231763E-03 | 3.821128E-03 |
| `F_jac` | `8fc4647d5f009a7bd4faa887d5efed80` | the same | the same | the same |
| `D` | `eee3b830b4e0fa85cef6454d7ad64379` | 5.951745E-10 | 3.489290E+00 | 1.100684E+01 |
| `Drow` | `eee3b830b4e0fa85cef6454d7ad64379` | the same | the same | the same |

Three facts this settles, all MEASURED:

1. **`F` and `F_jac` are bitwise the same array**, so the model base point is
   the state itself and the closure count of the model equals the closure
   count of the residual at this iterate. That is the condition section 8.1
   asks for, confirmed here at the actual iterate.
2. **`Drow` and `D` are bitwise the same array** on this route, as the code
   comment at `steady_newton.f90:1478` states for the judged row scaling off.
   So the merit the line search minimizes, the norm the Krylov cycle's
   relative tolerance is stated in, and the fixed scale the stagnation counter
   reads are one norm, and the predicted and observed reductions of section
   6.2 are comparable without a change of weights.
3. **The bound mask is empty**: `nact 0` in every record of both traces.

The entry values of the two arrays at the base, per row:

| row | `Drow` cell 1 | `Drow` cell 2 | `F` cell 1 | `F` cell 2 |
|---|---|---|---|---|
| mass | 3.489290E+00 | 3.229696E+00 | -8.143560E-04 | -1.345809E-03 |
| momentum | 2.321185E+00 | 2.148537E+00 | 3.231763E-03 | -9.187373E-04 |
| energy | 1.389698E+00 | 1.286352E+00 | -6.551809E-05 | -9.132226E-04 |

Every one of those six `F` entries equals, digit for digit, the `R_mass`,
`R_momentum`, `R_energy` column of the conservation export at the same cell,
so the vector the Newton system carries is the assembled residual itself.

**Where the merit lives at the entry state.** MEASURED from the same dump:
`||F/D||_2 = 1.6898687476565501E-03`, which the run prints as
`||Fs||2 = 1.69E-03`, and cells 1 and 2 carry 1.000000 of its square, cells 3
to 500 carrying 0 to the last bit. The six largest contributions are the
momentum row of cell 1 (0.6788 of the merit), the energy row of cell 2
(0.1765), the momentum row of cell 2 (0.0640), the mass row of cell 2
(0.0608), the mass row of cell 1 (0.0191) and the energy row of cell 1
(0.0008). **The merit and the certification measure localize at the same two
cells at this state**, although they are different functionals; that is not
true in general and is a property of this iterate.

### 5.3 The ledger of the observation group

The full table is `docs/audit_20260905/p2_base_stage_cd_20260921/ledger_obs.txt`,
one line per inner iteration, with the outer pass named in the first column.
The columns are: the forcing term asked of the cycle, the relative residual
the cycle reached and the number of operator products, the three hydrodynamic
rows entering the step, the linear residual left over the row of `F`, the
three rows of the state the iteration ended at, the judged distance of the
state and of the trial, `lam`, `dtau`, `||R||`, `||Fs||2`, the worst row with
its cell, and how many cells of each row stand outside their tolerance.

The part of it that answers the question, outer pass 1:

| it | eta | cycle reached | products | mass in | energy in | r_lin/F mass | r_lin/F energy | mass out | energy out | judged state to trial | lam | dtau |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.00E-01 | 4.689E-04 | 1 | 9.806E-01 | 1.071E+00 | 5.100E-04 | 8.073E-06 | 1.572E+00 | 8.632E-01 | 5.049E+07 to 2.069E+07 | 1.00 | 2.00E+00 |
| 2 | 2.85E-02 | 1.998E-04 | 1 | 1.572E+00 | 8.626E-01 | 1.736E-05 | 6.570E-04 | 1.967E+00 | 7.959E-01 | 2.069E+07 to 1.928E+07 | 1.00 | 4.00E+00 |
| 5 | 1.00E-01 | 8.750E-04 | 1 | 1.712E+00 | 1.735E+00 | 4.537E-05 | 2.561E-03 | 1.997E+00 | 1.778E+00 | 1.032E+07 to 5.511E+06 | 1.00 | 3.20E+01 |
| 9 | 1.00E-01 | 9.939E-03 | 1 | 7.201E-01 | 7.443E-01 | 1.335E-03 | 2.780E-03 | 2.065E-01 | 4.595E-01 | 7.443E+05 to 4.595E+05 | 1.00 | 5.12E+02 |
| 13 | 1.00E-01 | 1.994E-03 | 2 | 2.397E-02 | 9.238E-02 | 1.643E-02 | 1.683E-02 | 1.630E-02 | 8.061E-02 | 9.238E+04 to 8.061E+04 | 1.00 | 8.19E+03 |
| 20 | 1.00E-01 | 6.429E-02 | 3 | 3.886E-03 | 1.348E-02 | 4.276E-01 | 4.933E-01 | 7.067E-04 | 3.147E-03 | 1.348E+04 to 3.147E+03 | 1.00 | 1.05E+06 |
| 25 | 4.10E-03 | 2.869E-03 | 4 | 2.324E-07 | 1.218E-06 | 7.280E-03 | 6.145E-03 | 2.439E-08 | 2.559E-07 | 1.218E+00 to 2.559E-01 | 1.00 | 3.36E+07 |

`lam = 1.00E+00` at every one of the 25 iterations, the judged distance falls
at every one of them, and the cycle reaches the tolerance it was asked for at
every one of them. The row measures are not monotone: the mass row rises from
9.806E-01 to 1.967E+00 over iterations 1 to 3 and the energy row from
1.071E+00 to 1.845E+00 over iterations 4 to 6 before both fall. That is the
transient the code's own stagnation counter is written around, and it is the
reason the ranking is taken on the judged distance and not on `||R||`.

Outer pass 2, the same table, in short: the composition update between the
passes raises the energy row from 2.559E-07 to 6.512E-06 and moves the worst
row to cell 204 (r = 1.160); the first three iterations are cut to
`lam = 9.77E-04, 7.81E-03, 9.77E-04` with `dtau` falling from 1.00E-01 to
1.00E-03, and the solve then circles between an accepted full step and a
collapsed one at `||Fs||2 = 3.07E-09`, firing the damped Gauss-Newton escape
20 times, until the stagnation test stops it at the 60th cycle, after 40
iterations in which neither the judged distance nor the merit improved. The
best iterate is handed back, `||R|| = 6.512E-06`.

### 5.4 The ledger of the named numerical variant

`EXHALE_TRUST_REGION=1` replaces the line search by the scaled trust region
AND moves the merit onto the model row scaling (`steady_newton.f90:7734`).
**Its acceptance record is therefore not comparable with the observation
group's and is not compared with it here.** What the two share is the entry
state, which is bitwise the same: the first record of the two traces agrees on
every one of the eight hashes it carries, `Y 1B75B5AB91B05F7D`, `comp
37C52467DF120579`, `F` and `Fjac 23989D9F4D429AF1`, `D` and `Drow
93E32CCA029F9940`, the bound mask and `dtau`, and both runs print `start ||R||
= 1.071E+00, ||Fs||2 = 1.69E-03`. They also share the outcome of the pass.

Its own ledger is `ledger_tr_steps.txt` (the step control) beside
`ledger_tr_lin.txt` (the same linear-row record as the observation group).
Outer pass 1:

| it | delta | pred | actual | ratio | \|\|s\|\| | cuts | model_ok | accepted | worst row after |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 7.746E-05 | 3.461E-08 | 3.461E-08 | 1.000 | 3.873E-05 | 0 | T | T | momentum of cell 1 |
| 2 | 1.549E-04 | 9.009E-08 | 9.009E-08 | 1.000 | 7.746E-05 | 0 | T | T | momentum of cell 1 |
| 4 | 6.197E-04 | 5.927E-07 | 5.926E-07 | 1.000 | 3.098E-04 | 0 | T | T | momentum of cell 1 |
| 8 | 2.479E-03 | 3.014E-11 | 3.013E-11 | 1.000 | 1.844E-03 | 0 | T | T | mass of cell 3 |
| 11 | 9.915E-03 | 2.033E-14 | 9.158E-15 | 0.451 | 4.372E-03 | 0 | T | T | energy of cell 61 |
| 13 | 1.726E-03 | 7.421E-15 | 9.147E-17 | 0.012 | 6.904E-03 | 0 | T | F | energy of cell 37 |
| 22 | 3.452E-03 | 7.116E-17 | 7.113E-17 | 1.000 | 1.136E-04 | 0 | T | T | momentum of cell 37 |
| 24 | 3.452E-03 | 4.833E-21 | 4.823E-21 | 0.998 | 2.521E-08 | 0 | T | T | momentum of cell 78 |

`cuts = 0` at every iteration of the pass, one step rejected (iteration 13),
`info = 0` at iteration 25. **While the worst row is a base row, that is
iterations 1 to 8, the ratio of the observed to the predicted reduction is
1.000 to three decimals.** The model degrades only after the refusal has
left the base: 0.451 at iteration 11 (worst row cell 61), 0.012 at iteration
13 (cell 37).

In outer pass 2 the trust radius collapses from 2.507E-06 to 1.2E-15 with
`model_ok = F` on 8 of the first 20 iterations and ratios of -34, -2.7, -11,
+7.6 and +12.7; the predicted and observed reductions there are 1E-21 to
1E-24 against a merit of 3.58E-09, that is, at the arithmetic floor of the
merit, so a ratio taken on them classifies nothing.

## 6. Stage C, the three quantities kept apart

### 6.1 The three, named

```
r_model = -F - J(dY) - dY/dtau      what EXHALE_LINEAR_ROWS prints, per row
r_trial = F(Y + s)                  the judged slots of the trial the search measured
r_state = F(Y_adopted)              the rows of the state the iteration ended at
```

They are kept apart in the ledger: `r_model` is the column "linear residual
left, over the row of F", `r_trial` is the "of the trial" half of the judged
slots line, and `r_state` is the column "rows now" together with the `||R||`
of the iteration line. At outer pass 1 iteration 1 they read, in that order:
`5.100E-04, 9.469E-07, 8.073E-06` relative to their own rows of `F`; judged
trial distance 2.069E+07; and rows `1.572E+00, 2.104E-08, 8.632E-01`, that is,
`r_model` small, `r_state` larger than the state it started from in the mass
row and smaller in the energy row, and the judged distance down by a factor
2.44. A small `r_model` beside a raised mass row classifies nothing by itself,
which is why the trust-region group's `pred` and `actual` are the number that
decides, and they say 1.000.

**What `r_model` as printed is, exactly, and what it is not.** It is the
residual of the SHIFTED system `(J + I/dtau) dY = -F`, so it contains the
pseudo-transient shift `dY/dtau` and is not the residual of the Newton system.
It is also a maximum over cells of a row, divided by the maximum over cells of
the same row of `F`, and not the norm the cycle's relative tolerance is stated
in. The clean number for the cycle is `gm_resid_rel`, the relative residual
the cycle reports, and the clean number for what the returned direction leaves
against the operator is the TRUE relative residual of `[diag 4]` in the
variant group. **The shift is not reported separately by either hook**, and
section 7's order of inference asks for it; that is a gap in the instruments,
recorded here rather than worked around. It shows plainly in outer pass 2 of
the observation group, where `dtau` falls to 1.0E-01 and the printed
`r_model/F` reaches 2.989E+05 in the mass row while the cycle reports
4.678E-02 relative and "the requested linear tolerance was reached": the large
number is the shift, not a linear residual.

### 6.2 Predicted against observed reduction

Because `Drow` and `D` are bitwise one array on this route (section 5.2), the
merit `||F/D||_2` the run prints as `||Fs||2` and the norm the cycle's relative
tolerance is stated in are the same norm, and the two reductions can be put
beside each other without a change of weights. Outer pass 1 of the observation
group:

| it | cycle's relative residual | merit before | merit after | observed relative merit |
|---|---|---|---|---|
| 1 | 4.689E-04 | 1.69E-03 | 2.00E-04 | 0.118 |
| 2 | 1.998E-04 | 2.00E-04 | 1.82E-04 | 0.910 |
| 3 | 3.276E-04 | 1.82E-04 | 1.45E-04 | 0.797 |
| 9 | 9.939E-03 | 3.47E-06 | 1.57E-06 | 0.452 |
| 23 | 3.910E-02 | 3.61E-10 | 4.76E-11 | 0.132 |

The observed reduction is one to three decades short of the linear model's,
which is expected and is not by itself a classification: the step solved is
the SHIFTED one, and with `dtau` between 2 and 3.4E+07 the shift deliberately
shortens it. The quantity that removes the ambiguity is the variant group's
`pred` and `actual`, which are formed on one merit at one iterate, and while
the worst row is a base row they agree to 1.000.

### 6.3 What holds the cycle on this system

All three diagnostics of section 7 were reached, at outer pass 1 iteration 1
of the variant group, and all three PRODUCED OUTPUT.

**The Krylov size scan** (`[diag 16]`, `||b|| = 1.690E-03`): subspace sizes
40, 80, 160 and 320, with and without a second orthogonalization, all return
`reduced 1.27889E-02, TRUE 1.27907E-02, ||x|| 3.176E-03` after 3 products, and
all reach the requested tolerance. Enlarging the subspace buys nothing at this
iterate, and the reduced and the true residual agree to five digits, so the
reduced quantity is a faithful statement about the system.

**The Ritz values of the preconditioned operator** (`[diag 15]`, 200
products): largest magnitude 9.99385E-01, smallest 2.52128E-04, ratio
3.96381E+03; 6 of 200 below 1E-2 and none below 1E-4; 194 of 200 in complex
pairs. The smallest Ritz vector carries 0.7850 in the mass rows and 0.2150 in
the energy rows and is spread over 231 cells, its largest at 116, 31, 169 and
113, that is, in the wind and not at the base. **There is no near-null mode of
the preconditioned operator, and the slow modes that exist do not live at the
base.**

**The band difference** (`[diag 17]`), on the smallest Ritz direction:
`||A v|| = 6.91613E-01`, `||(A - A_band) v|| = 1.82248E-03`, ratio
2.63512E-03. What the band omits is 0.26 per cent of the action, and **nine
tenths of it sits in 2 cells, the largest of them 1 and 2**, with the energy
rows carrying 0.6205 of the difference against 0.2241 for the mass rows and
0.1553 for the momentum rows. Column by column the omission is entirely the
frozen-radiation part: for the mass column of cell 116 the part outside the
band is 8.441E-07 of a column of norm 3.702E+01, and the two parts sum to the
whole difference exactly.

So the base IS where the banded preconditioner is least faithful, and the
amount it is unfaithful by, 0.26 per cent of the action along the slowest
direction, is far too small to obstruct a cycle that reaches its tolerance in
one to four products.

The band itself at that iterate (`[diag 3]`): `dgbtrf info 0`,
`|U(j,j)|` between 9.140E-01 and 1.132E+04 with the smallest at the energy of
cell 500, 1-norm of the shifted band 4.011E+04, `dgbcon rcond` 1.732E-06,
pivot ratio 8.074E-05.

**The binding row, term by term** (`[diag 10]`). The row the merit is largest
on at the entry state is **the momentum row of cell 1**, scaled residual
1.392E-03 on a row scale of 2.321E+00, and its derivatives are:

| unknown | d(row)/d(unknown) | scaled by the row |
|---|---|---|
| mass of cell 1 | 1.40232E+02 | 6.04138E+01 |
| momentum of cell 1 | 8.75900E+02 | 3.77350E+02 |
| energy of cell 1 | 1.20013E+03 | 5.17032E+02 |
| mass of cell 2 | 3.69106E+01 | 1.59016E+01 |
| momentum of cell 2 | -1.12605E+03 | -4.85119E+02 |
| energy of cell 2 | 1.27426E+03 | 5.48969E+02 |
| mass of cell 3 | -5.53370E+00 | -2.38400E+00 |
| momentum of cell 3 | 2.80543E+02 | 1.20862E+02 |
| energy of cell 3 | -1.47758E+02 | -6.36565E+01 |
| mass of cell 4 | -6.88024E-06 | -2.96411E-06 |

The row is coupled to its own cell and its two neighbors at comparable
strength, with nothing beyond cell 3, and no entry is small enough for the row
to be singular in any of its own unknowns. There is no row standing against
every scaling here.

`EXHALE_FRONT_ROW` covers species rows only and was not set; it could not have
spoken on this system, which carries no species row.

## 7. The base rows, attributed term by term

From the conservation export of section 13.2, evaluated on the two states.
Terms are contributions already divided by the cell volume, in code units.

### 7.1 The mass row at cells 1 and 2, at the checkpoint

The mass row has no source, so the row IS the flux divergence and the row
measure is `|Phi_lo - Phi_hi| / max(|Phi_lo|, |Phi_hi|)` with `Phi` the area
weighted face mass flux.

| | cell 1 | cell 2 |
|---|---|---|
| `r^2 Phi` through the lower face | 4.230098E-07 | 2.655186E-07 |
| `r^2 Phi` through the upper face | 2.655186E-07 | 5.147205E-09 |
| the same, in wind-window units | 82.1824 to 51.5850 | 51.5850 to 1.0000 |
| `dF_mass` = `R_mass` | -8.143560E-04 | -1.345809E-03 |
| row scale, DERIVED from the export as `max(|r^2 Phi|)/volume` | 2.18730E-03 | 1.372425E-03 |
| row measure, DERIVED | 3.7229E-01 | 9.8061E-01 (the run's own 9.8061E-01, which is how the derivation is checked) |

**The whole content of the mass row at the base is that the base face carries
82.18 times the wind's mass flux and the column sheds it within two cells.**
The measure at cell 2 is, identically, `1 - 1/51.585 = 0.98061`. From cell 3
upward `r^2 Phi` is 5.147205E-09 with a row residual of 1E-13, that is, at the
rounding floor.

### 7.2 The energy row at cell 1, at the checkpoint

| term | value |
|---|---|
| `dF_energy` (the flux divergence, gravitational work already inside it with a positive sign) | -6.116629E-05 |
| of which the gravitational work | +3.653793E-05 |
| `S_energy` | 0 in every branch |
| `heat` | 4.366900E-06 |
| `cool` | 1.510249E-08 |
| `heat - cool` | 4.351798E-06 |
| `R_energy` | -6.551809E-05 |
| row scale, `max(|dF|, |S|, heat, cool)` | 6.116629E-05, which is `|dF|` |
| row measure | 1.0711 |

So at cell 1 the transport term is 14.1 times the radiative net, the scale IS
the transport term, and the measure is identically
`1 + (heat - cool)/|dF_energy| = 1.0711`. The energy row at cell 1 refuses
because energy is being carried into the cell through the base face 14 times
faster than the radiation can account for, which is the same statement the
mass row makes, in energy.

### 7.3 The same two rows after outer pass 1

| | cell 1 | cell 2 |
|---|---|---|
| `r^2 Phi` lower face, in wind units | 1.0000 | 1.0000 |
| `R_mass` | -2.458778E-13 | -2.708739E-13 |
| `dF_energy` | 4.369037E-06 | 4.533413E-06 |
| `heat - cool` | 4.369037E-06 | 4.533414E-06 |
| `R_energy` | -9.141415E-14 | -1.763792E-13 |

The energy row of cell 1 has become the radiative balance it should be, to
eleven digits, and the mass row is at its rounding floor. The re-evaluation of
that state as a frozen state gives mass 2.438614E-08, momentum 1.017326E-14,
energy 2.583606E-07, against the 2.439E-08, 1.017E-14, 2.559E-07 the run
itself reported, that is, inside the known round-trip band of the energy row.

## 8. The verdict, in the terms of section 6's decision table

**At the base rows: none of (i) to (iv).**

| reading | verdict at the base, and its evidence |
|---|---|
| (i) resolved action, small linear residual, nonlinear prediction fails | NO. The action is resolved (phase 4), the cycle reaches its tolerance in 1 to 4 products, and where the worst row is a base row the variant group's observed reduction equals its prediction to 1.000 at 8 consecutive iterations. The row measures rise transiently, the judged distance and the merit do not. |
| (ii) the linear residual itself stalls | NO. The relative residual reached is 1.998E-04 to 9.455E-02 against forcing terms between 4.10E-03 and 1.00E-01, reached at every iteration of pass 1; the size scan shows 40, 80, 160 and 320 give the identical reduced and true residual after 3 products; the Ritz spectrum has no value below 1E-4 and a ratio of 3.96E+03; the band omits 0.26 per cent of the action. |
| (iii) one row at one cell standing against every scaling | NO. The binding row of the merit at entry is the momentum row of cell 1, and its derivatives with respect to its own and its neighbors' unknowns are O(10) to O(500) scaled, with no small entry. All three hydrodynamic rows fall below their certification tolerances in the same solve. |
| (iv) steps cut at the same bound every pass | NO. There are no species unknowns and therefore no box. `nact = 0` in every record of both traces, `[diag 5] unknowns held at a bound 0`, `cuts = 0` at every trust-region iteration, and `lam = 1.00E+00` at all 25 line-search iterations of pass 1. |

**What the decision table does classify, and it is not at the base.** After
the base is solved the case refuses at two places, both above r = 1.1:

- the elemental transport He/H partition, 7.218E-05 against 1.0E-05 at cell
  295, a wind cell, in both groups;
- the hydrodynamic rows after the composition update, whose worst certified
  entry is the energy row at cell 196 (r = 1.1388) at 6.47E-06 to 6.48E-06 in
  both groups, and whose worst row during the second solve is the mass row of
  cell 266 to 268 (r = 1.47 to 1.49).

The second solve's behavior is reading (i) in the wind, not at the base: the
cycle reaches its tolerance and the trust region's model ratio is -34, -2.7,
-11, +7.6, +12.7, that is, the model has no predictive value there. The
qualification is that the reductions it is taken on, 1E-21 to 1E-24 against a
merit of 3.58E-09, are at the arithmetic floor of the merit, so what that
second solve measures may be the floor rather than the model. Deciding that
needs a measurement of the merit's own reproducibility at those iterates
against the size of the reductions asked for; the variant group prints the
first (`merit reproducibility at the iterate 0.000E+00` at some iterates,
2.118E-22 at the entry of pass 1) but the two were not put together here.
**That is the honest limit of this phase: it decides the base and it does not
decide the wind.**

## 9. What this does not establish, and what contradicts the plan

1. **It is one checkpoint.** `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` and
   `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` were not run. The second
   of those refuses with `mass 1.163E+00 at cell 1` in its 2026-09-20 log, so
   the base question is not closed for it by anything here.
2. **It is not a statement about the historical runs.** The campaign runs of
   2026-09-20 were seeded from the certified x0.10 case and interpolated onto
   this grid; this work starts from `g0002`. They are different initial states
   and their outcomes are not comparable. What is comparable is that both were
   asked of the same case at the same XUV level and that this one reaches
   `info = 0` at the base.
3. **Section 17.1 is contradicted in its inference, not in its measurement.**
   The 82.18 is reproduced here exactly on the same state. What does not
   follow from it is "the base state these cases are given is not compatible
   with a stationary wind at that base": the current operator moves that face
   to 1.0000 in 25 Newton steps. Section 17.2's reading 2 is the one supported
   at this checkpoint.
4. **Section 15's row for phase 9 asks for a mechanism, not a correlation.**
   The mechanism found at the base is that the stored state carries a base
   face mass flux 82 times the wind's, that the two refusing rows are exactly
   the statement of that excess and of the energy it carries, and that the
   route removes it without a cut, a bound or a stalled cycle. The mechanism
   of what remains, in the wind, is not established here.
5. **Two fields of section 7.2 are not in the run's output on the default
   route**: the actual scaled displacement and the limiter or contact event
   count. They are named NOT IN THE RUN'S OUTPUT in section 5.1 and were not
   substituted by anything.
6. **The pseudo-transient shift is not reported apart from the linear
   residual** by either hook, which section 7's order of inference asks for.
   Recorded as a gap in the instruments.
7. `EXHALE_FRONT_ROW`, `EXHALE_KRYLOV_ON_THE_BALL`, `EXHALE_KRYLOV_TOL_ABS`,
   `EXHALE_GM_ORTHO` and `EXHALE_RESID_QUAD` were NOT SET in any run here, so
   nothing they would have measured is claimed. `EXHALE_CARRIER_ROW_TERMS`
   would have written nothing: this system carries no carrier row.

## 10. Artifacts

Archived under `docs/audit_20260905/p2_base_stage_cd_20260921/`, with
`MD5SUMS.txt`:

| file | what it is |
|---|---|
| `R_frozen_run.log` | Mode R, the frozen evaluation of `g0002` |
| `R_afterpass1_run.log` | Mode R, the frozen evaluation of the state outer pass 1 returned |
| `dump_entry_run.log`, `dump_entry_l15_trace.txt`, `dump_entry_summary.txt` | the array dump at the entry of pass 1 and its summary |
| `obs_run.log`, `obs_l15_trace.txt`, `ledger_obs.txt` | the observation group |
| `obs1pass_run.log` | the control with the trace off, one outer pass |
| `tr_run.log`, `ledger_tr_steps.txt`, `ledger_tr_lin.txt` | the named numerical variant |
| `base_terms_entry.txt`, `base_terms_after_pass1.txt` | the base face fluxes and the cell terms of both states |
| `mkrun.sh`, `base_terms.py`, `read_dump.py`, `ledger.py` | the commands that made and read them |

The run directories themselves, with `output/`, the conservation exports and
the binary array dumps, stand in the session scratch at
`<scratchpad>/p2_base_cd/LHS1140b/models/case/{R_frozen, R_afterpass1,
dump_entry, obs, obs1pass, tr}/`.

Every invocation ran under `env -i` with the full environment set written
beside it, in its own directory, and exited 0.
