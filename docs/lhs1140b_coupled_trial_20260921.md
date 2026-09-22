# The M1 coupled trial on `molecular_scalar_gj1132_kzz1e9/HeH0.083`

Phase 7 of `docs/PLAN_20260920_rev9.md` section 15, the coupled trial of
section 9.4. Mode C (rev9 section 3): a retained trajectory advanced under the
current solver from a stored generation, with a budget declared before the run
and a recorded ending. Nothing is published: the trial ran on a scratch copy
outside the repository and no file of the case directory was written, replaced
or removed.

Written 2026-09-21. Every number carries a label: MEASURED (this trial
produced it), READ (from a source file, a stored product or an earlier
document) or INSPECTED (read from the source text).

---

## 1. The verdict, first

**The coupled route removes the movement bound and does not remove the drift.
It replaces a composition half that reached its own fixed point under a
saturated bound by a solve that does not move the state at all.** On this case,
from generation `g0004_20260920T053734Z_79b42a03`, with
`Coupled carrier solve: True` and therefore also the scaled trust region:

- the carrier relaxation does not run: MEASURED, the pass-1 block reads
  `carrier residual of the returned composition 0.00E+00 ... [the run carries
  no relaxed carrier]`, so the movement bound that saturated at cells 418 to
  426 in every pass of the control is gone, and with it every one of the 25 to
  30 refused carrier trials a control pass carried;
- the first hydrodynamic-and-species solve spent **41 trust-region iterations
  and 6632 residual evaluations in 7171.84 s**, took the merit from 9.720E-02,
  the value its own line reports at the entry iterate, to a best of 9.718E-02
  over its acceptance window, then **STAGNATED** (`info = 2`, `neither the judged distance
  nor the merit improved in 40 iterations, monotone search included`) and
  handed back the iterate it had started from;
- the certificate of pass 1 therefore reads **the entry values to four
  digits**: carrier balance H2 8.150E-04 at cell 499, elemental transport He/H
  2.503E-04 at cell 280. The control reached exactly those two numbers at the
  same pass with `info = 0` in 8.6 s.

**No pass certifies**, and the two entries that refuse are the two the entry
state already carried. The one completed pass of the trial is a localized
recorded failure of the kind the gate of rev9 section 15 admits, not a
certification.

**The budget of six outer passes was NOT spent.** One pass completed; the
second solve was still running at iteration 59 when the 18000 s safety ceiling
stopped the process (`exit = 124`). That is an EXTERNAL TERMINATION and is
never read here as an ending. Section 4 gives what it cost and why six passes
were out of reach.

**The result is a statement about the combined configuration.** Registering the
carrier and element rows changes the unknown set AND switches the globalization
to the scaled trust region in one step, and it removes the carrier half of the
alternation as well (section 3.1). Nothing here attributes the outcome to one
of the three alone.

---

## 2. Identity

### 2.1 The operator

| field | value |
|---|---|
| binary | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x`, md5 `89ec67149aa452eea3deeb19052f83d1` (MEASURED), the tree's binary of record and the binary phase 6 used |
| binary mtime | 2026-09-21 13:52:38 (MEASURED) |
| repository | `93eed86` with a dirty working tree (READ, `git log -1`, `git status --porcelain` counts 66 entries) |
| the source text and the binary | NOT the same text. MEASURED: six files under `src/` are newer than the binary (`EXHALE_main.f90`, `util_ion_eq.f90`, `steady_newton.f90`, `molecular_seed_from_atomic_state.f90`, `input_read.f90`, `load_IC.f90`), and the md5 of the sorted `md5sum` of every `src/**/*.f90` and `*.inc` is `dcca114508a9ac4b459a60df6cdc3e60` now against `fb3b5dd2c3410a4764d6d24f00190872` when phase 6 measured it (READ from `docs/lhs1140b_alternation_audit_20260921.md` section 3.1). Two other workers are editing this tree. Every INSPECTED citation below is therefore of the text as it stands today and is quoted by its wording, not by a line number; every behaviour is the binary's and was MEASURED in the run |
| threads | `OMP_NUM_THREADS=8`; the run pins the BLAS pool to one thread itself (READ, the run's own line) |
| `EXHALE_RESID_QUAD` | UNSET, so the assembly selector of rev9 section 13.6 took its default. It is in no configuration record and is stated here because only the setter can state it |

### 2.2 The checkpoint, and that the trial starts where the control started

| field | value |
|---|---|
| case | `LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083` |
| generation | `g0004_20260920T053734Z_79b42a03`, the `latest_complete` of `state_index.json` (READ); `latest_certified` is `none` (READ) |
| identity record | `docs/lhs1140b_identity_records_20260921.md` section 3.4 |
| `Hydro_ioniz.txt` md5 | `fd882d8918052069f8ed229dd376f8c4` (MEASURED on the trial's `_IC` pair, equal to the manifest value READ in the identity record) |
| `Ion_species.txt` md5 | `87288b03ad88a92b7b8b8b6c2b1bdd66` (MEASURED, equal to the READ manifest value) |
| boundary model of the loaded state | `characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3`, this run's own (READ, the run's `load_IC` line) |
| the binary that wrote the state | `git=3c73905ca8a7`, another build (READ, `load_IC`), so this is the CURRENT operator on an OLD state |

**The composition refresh stage is the same in the two trajectories, and the
entry certificate proves it.** Both runs state `Restart intent: stationary`,
which measures the state as loaded and takes no CFL step, and neither states
`stationary equilibrate`. The trial's own entry certification, before any
relaxation and before the first solve, MEASURED:

| row | control entry (READ, phase 6) | trial entry (MEASURED) |
|---|---|---|
| hydrodynamic mass | 1.393E-08 at cell 1 | 1.393E-08 at cell 1 |
| hydrodynamic energy | 7.651E-09 at cell 3 | 7.651E-09 at cell 3 |
| carrier balance H2, gated | 8.150E-04 at cell 499 | 8.150E-04 at cell 499 |
| elemental transport He/H, gated | 2.503E-04 at cell 280 | 2.503E-04 at cell 280 |
| `||R||` of the work state | 1.3929E-08 | 1.3929E-08 |

So the two trajectories start from the same state at the same stage, which is
what rev9 section 9.4 requires of the comparison.

**And the trial reproduces phase 6b's coupled reading of the same state.**
Phase 6b measured, on this generation with the same two rows registered, a
largest judged species row of 8.1496E+01, the H2 carrier of cell 499 (READ,
`docs/lhs1140b_coupled_action_20260921.md` section 11). The trial's own
`judged rows` line at the entry iterate reads `distance 8.150E+01 window
8.150E+01 ||R|| 1.393E-08 hydrodynamic 7.175E-01 element 2.503E-04 carrier
8.150E-04` (MEASURED). The frozen probe of phase 6b and the solve of phase 7
are looking at the same object.

### 2.3 The one line that differs

MEASURED, `diff` of the two input files: exactly one line,
`Coupled carrier solve: False` against `Coupled carrier solve: True`. Both
files carry the two edits phase 6 recorded and justified (an absolute
`Spectrum file:` and `Restart intent: stationary` in place of the case's
`stationary equilibrate`), so the trial inherits the control's configuration
unchanged in every other respect. The trial input md5 is
`f510b39051b6bdc8e0bbfe9bd660c1c5` and `base.inp` is
`1e33d197ead9269661dc6bb0d88c0e26`, the same `base.inp` the control used
(MEASURED).

---

## 3. Admissibility, checked before the case was selected

`Coupled carrier solve: True` is refused together with
`Ionization transport: True`. INSPECTED, `input_read.f90`, inside the block
that opens `if (ionization_transport) then`: the reader writes
`"Ionization transport: True" with "Coupled carrier solve: True" is refused.
The coupled row registry carries every carrier unknown as the species mass
fraction of its f_sp column, and a stage unknown is a fraction per element
nucleus` and stops the run with `error stop 1`. A second requirement in the
same file refuses `Coupled carrier solve: True` on a molecular configuration
whose carriers are not transported.

The case satisfies both:

| key | value | where |
|---|---|---|
| `Molecular chemistry` | True | READ, `input.inp` |
| `Molecular carrier transport` | True | READ, `input.inp` |
| `Ionization transport` | not stated, so False | READ, `input.inp`; and `ionization_transport F` in the trial's own `EXHALE_resolved.out` (MEASURED) |
| `carrier_in_newton` | T | MEASURED, the trial's `EXHALE_resolved.out` |

MEASURED and not only read: the run launched, parsed, loaded the restart and
entered the solver without a refusal, and its `load_IC` block states
`the restart changes the ROUTE and not the equations: carrier_newton=F ->
carrier_newton=T`. **The combination is admissible in the binary of record.**

### 3.1 What registering the species rows also changes

Two things move together and the trial cannot separate them.

1. **The unknown set.** INSPECTED, `EXHALE_main.f90` calls
   `set_transported_species_rows(block_now .or.
   coupled_block_jacobian_action_requested())` immediately before
   `solve_steady_jfnk`, with `block_now = carrier_in_newton`. MEASURED in the
   trial: `cost: 5 unknowns per cell, 29 colors` against the three of the
   control's route.
2. **The step control.** INSPECTED, `steady_newton.f90`:
   `use_tr = (nspec_row .gt. 0) .or. trust_region_on_the_hydrodynamic_rows`,
   with the comment `the scaled trust region wherever the system carries a
   species row`. `EXHALE_TRUST_REGION` was left UNSET in this trial, so the
   trust region is the route here and not a variant; the diagnostics of rev9
   section 7 that live inside `trust_region_step` are therefore OBSERVATION on
   this route.

3. **The carrier half of the alternation stops running.** INSPECTED,
   `EXHALE_main.f90` guards `relax_photochemical_composition` with
   `transported_rows_exist() .and. .not. block_now`. MEASURED, the trial's
   pass-1 block: `carrier residual of the returned composition 0.00E+00 ...
   [the run carries no relaxed carrier]`. The element half still runs, because
   `species_alternated = he_diffusion .or. (transported_rows_exist() .and.
   .not. block_now)` and `he_diffusion` is T, so the pass cap stays at the
   six the run asked for.

So the trial's result concerns the combined configuration: five unknowns a
cell, the scaled trust region, and one half of the alternation removed.

---

## 4. The budget, declared before the run, and what it bought

Predeclared in the units rev9 section 15 asks for, and estimated from the
control's own measured cost (READ, `docs/lhs1140b_alternation_audit_20260921.md`
section 2.3: 6 passes, 6 solves, 150 JFNK iterations, 133 Krylov products,
1759.1 s).

| limit | value | why |
|---|---|---|
| outer passes | 6, `EXHALE_OUTER_PASSES=6` | the control's budget, so that the two trajectories are comparable pass for pass |
| hydrodynamic-and-species solves | at most 6, one per pass | the route takes one solve per outer pass |
| inner iterations | no cap imposed | `EXHALE_JFNK_MAXIT` caps the RUN and not the solve (READ, phase 6 section 2.2), so setting it would have ended the trajectory after one pass |
| wall ceiling | 18000 s (5 h), `timeout 18000` | a SAFETY ceiling only. A stop at it is recorded as an external termination and is never converted into an ending (rev9 section 12) |
| threads | `OMP_NUM_THREADS=8`, one trajectory at a time | two other workers hold the machine |

**The environment actually in force** (rev9 section 4.2 item 1), MEASURED from
the launch script that ran it: `OMP_NUM_THREADS=8`, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=6`, `EXHALE_DIFFUSION_CHECK=1`,
`EXHALE_CARRIER_ROW_TERMS=1`, `EXHALE_CARRIER_DEBUG=1`, `EXHALE_ELEM_DIAG=1`,
`EXHALE_PRECOND_SPECTRUM=1`, `EXHALE_BAND_DIFFERENCE=1`,
`EXHALE_KRYLOV_SIZE_SCAN=1`; and explicitly unset: `EXHALE_RESID_QUAD`,
`EXHALE_JFNK_MAXIT`, `EXHALE_CARRIER_TRUST`, `EXHALE_CARRIER_TRUST_HOLD`,
`EXHALE_DIFF_OMEGA`, `EXHALE_TRUST_REGION`, `EXHALE_LINEAR_ROWS`,
`EXHALE_FRONT_ROW`. The control's set was the same minus the last four of the
set list (READ, phase 6 section 5).

### 4.1 What was spent, and how the trajectory ended

MEASURED. Launched 2026-09-21T18:24:50+09:00, stopped
2026-09-21T23:24:50+09:00 by the safety ceiling, `exit = 124`.

| quantity | spent |
|---|---|
| outer passes completed | **1 of 6** |
| hydrodynamic-and-species solves | 1 completed (`info = 2`), 1 in progress at the termination |
| trust-region outer iterations | 41 in solve 1, 59 in solve 2, **100 in all** |
| Krylov products inside the solves | 40 at every one of the 100 iterations, the whole subspace each time (section 6.2), so **4000** |
| residual evaluations | 6632 in solve 1 (MEASURED, the solve's own `cost:` line, which counts the diagnostic products separately); solve 2 did not print its own count because it did not end, recorded as NOT AVAILABLE |
| composition elimination sweeps | 120018 over the 6632 evaluations of solve 1, 18.10 per evaluation (MEASURED) |
| carrier trials | 0, the carrier relaxation does not run on this route |
| wall of the completed pass | 7171.84 s (MEASURED, the pass line) |
| wall of the whole run | 18000 s, ended by the ceiling |

**Why six passes were out of reach.** One completed pass of the trial cost
7171.84 s against the control's 8.6 s for the same pass and 1759.1 s for all
six of its passes (READ, phase 6 section 2.3). Two things carry that factor and
this record does NOT separate them, because the run instruments neither
separately:

- the route itself: 5 unknowns a cell and 29 colors against 3 unknowns, a
  solve that has real work to do where the control's first solve had none (the
  merit at the entry iterate is 9.720E-02 on the coupled route against
  `||Fs||2 = 8.45E-11` on the control's, READ from the control's own start
  line), and a trust-region iteration that evaluates the residual several times
  where the line search evaluated it once;
- the four diagnostics of rev9 section 7. What IS measured of their cost is the
  operator work: each firing runs 1080 products in the subspace scan and 600 in
  the Ritz measurement, plus the band difference's own, and there were four
  firings, against the 4000 products the 100 iterations of the two solves took.
  So the diagnostics at least doubled the operator work of the run. Their wall
  share is not separately measured and is not claimed.

A repeat of this trial WITHOUT those four keys is what would spend the six
passes, and it is named here as the next measurement rather than taken as
done.

---

## 5. The two trajectories, pass for pass

Both start from generation `g0004_20260920T053734Z_79b42a03` at the same
composition refresh stage (section 2.2), under the same binary, with
`OMP_NUM_THREADS=8`, and differ in the one input line of section 2.3. The
control column is READ from `docs/lhs1140b_alternation_audit_20260921.md`
section 6.1 and, for `||R||`, from the six `(JFNK) done info=` lines of the
control's own `run.log`. The trial column is MEASURED.

| pass | route | hydro ending | `||R||` at that ending | carrier balance H2, gated | elemental transport He/H, gated | joint residual | certified hydrodynamic rows | wall [s] |
|---|---|---|---|---|---|---|---|---|
| 1 | control, `False` | `info = 0` | 1.836E-08 | 8.150E-04 at cell 499 | 2.503E-04 at cell 280 | 8.150E+01 | mass 1.836E-08 cell 1, energy 7.651E-09 cell 3, both within | 8.6 |
| 1 | **trial, `True`** | **`info = 2`, STAGNATED** | 1.107E-08 | **8.150E-04 at cell 499** | **2.503E-04 at cell 280** | 8.150E+01 | mass 1.107E-08 cell 1, momentum 5.248E-13 cell 1, energy 5.803E-09 cell 3, all within | **7171.8** |
| 2 | control, `False` | `info = 0` | 5.817E-09 | 6.290E-03 at cell 500 | 4.532E-04 at cell 280 | 6.290E+02 | mass 3.681E-09, energy 5.817E-09, within | 413.7 |
| 2 | trial, `True` | NOT REACHED: the safety ceiling stopped the process inside outer iteration 60 of this solve, whose last completed iteration line is 59 | -- | -- | -- | -- | -- | > 10828 |
| 3 to 6 | control, `False` | `info = 0` every pass | 5.571E-09, 7.687E-09, 7.573E-09, 1.354E-08 | 6.258E-03, 6.388E-03, 2.950E-03, 2.893E-03 (cells 500) | 6.069E-04, 6.961E-04, 6.968E-04, 7.063E-04 (cell 280) | 6.26E+02, 6.39E+02, 2.95E+02, 2.89E+02 | within at every pass | 346.8, 344.3, 354.1, 291.5 |
| 3 to 6 | trial, `True` | NOT REACHED | -- | -- | -- | -- | -- | -- |

**The one row that is a comparison is pass 1, and it is the same two numbers
twice.** The control solved the hydrodynamic subsystem and then moved the
composition by 2.740E-02, of which the carrier half was 2.740E-02 and the
element half 1.340E-03. The trial's solve, given the carrier as a Newton
unknown, made no progress at all and returned its entry iterate; its element
relaxation then applied exactly the same 1.340E-03 the control's did, and the
carrier half applied nothing.

### 5.1 The iterates the trial reached inside its two solves

These are ITERATES, not certifications, and are labeled so. They are the
`judged rows` line the solve writes at every accepted iterate, whose element
and carrier entries are the same gated measures the certification reports
(they agree with the entry certificate to four digits at iteration 0).

| point | judged hydrodynamic distance | `||R||` | element, gated | carrier, gated |
|---|---|---|---|---|
| entry, before any step | 7.175E-01 | 1.393E-08 | 2.503E-04 | 8.150E-04 |
| solve 1, after the first accepted step | 5.874E+03 | 4.349E-07 | 2.504E-04 | 8.151E-04 |
| solve 1, the worst iterate before the monotone restart | 1.281E+05 | 5.947E-05 | 2.493E-04 | 8.036E-04 |
| solve 1, after the restart to the best iterate | 5.884E-01 | 1.142E-08 | 2.503E-04 | 8.150E-04 |
| solve 2, entry (after the element relaxation of pass 1) | 4.018E+06 | 5.652E-02 | 1.282E-04 | 8.150E-04 |
| solve 2, iteration 59, the last completed iterate before the ceiling | 3.870E+06 | 5.442E-02 | 1.063E-04 | 7.637E-04 |

**What this says and what it does not.** Over 100 trust-region iterations the
gated carrier row moved from 8.150E-04 to 7.637E-04, 6.3 per cent, and the
gated element row from 2.503E-04 to 1.063E-04; but the element fall is the
pass-1 element RELAXATION's, since 2.503E-04 becomes 1.282E-04 between the two
solves and not inside either of them, and the whole of it is bought at a hydrodynamic
judged distance that rises from 7.175E-01 to 3.870E+06. The control, over its
six passes, moved the carrier the other way, 8.150E-04 to 2.893E-03, and the
element from 2.503E-04 to 7.063E-04, with the hydrodynamic rows inside their
tolerances at every one of the six certifications.

So the honest comparison on the species rows is: **the control drifts away from
the entry values with a wind in solution; the trial does not drift, because its
solve does not move, and where it does move it moves a state whose wind is six
decades outside certification.** Neither route approaches a joint fixed point on
this case.

### 5.2 The compact state norms: NOT AVAILABLE, and why

Rev9 section 9.2 asks for `Y(k+1) - Y(k)` and `Y(k+1) - Y(k-1)` in the blocks
phase 6 used. They are read from `output/diffusion_pass_profiles.txt`, which
`write_diffusion_pass_profile` appends one 500-row block to at the END of each
outer pass. MEASURED: the file exists, 67187 bytes, 2 header lines and 500 data
rows, and its first column carries the single value 1. **One pass completed, so
the file holds one block and neither difference exists.** This is recorded as
NOT AVAILABLE with its reason and is not substituted by a difference between
iterates, which is a different quantity.

### 5.3 The gate of rev9 section 4.2, for every diagnostic file this trial used

| requirement | how it is met |
|---|---|
| the complete environment in force, and the mode | section 4; Mode C |
| the file, its existence, size, header and row count | `output/carrier_row_terms.txt`: present, 288177 bytes, 9 header lines, 500 cell rows, 500 face rows and 500 channel rows, one carrier (`H2`), the count taken from the file's own header (MEASURED). `output/diffusion_pass_profiles.txt`: present, 67187 bytes, 2 header lines, 500 data rows, one pass block (MEASURED). `output/element_flux_profile.txt`: present, 96310 bytes (MEASURED) |
| which state the rows belong to | the writer emits no identifier, so the invocation identifies it. All three files carry the modification time 20:24:22, which is the certification of outer pass 1, and no later certification ran, so they are the rows of the pass-1 state. The entry-state file was overwritten by that certification and is NOT separately available from this run; the entry state's own carrier row terms are in the phase 6 archive (`c1_kzz1e9_HeH0.083_eval`), measured from the same generation |
| the exit class | `exit = 124`, the 18000 s safety ceiling. An external termination, recorded as one |

---

## 6. The diagnostics of rev9 section 7, taken on this route

These four keys fire inside `trust_region_step` behind `elem_diag_here`, which
is true at outer iterations 1, 2 and `elem_diag_third` (INSPECTED,
`steady_newton.f90`: `elem_diag_here = elem_diag_on .and. (iter .le. 2 .or.
...)`, with `elem_diag_third` defaulting to 60). Phase 6b left them to this
phase because every measurement there evaluated a frozen state and stopped
before a solve (READ, `docs/lhs1140b_coupled_action_20260921.md` section 12).
On this route the trust region is the route and not a variant, so they are
observation.

**Which produced output** (rev9 section 4.2). MEASURED, counted in the trial
log by their own tags:

| key | tag | fired | verdict |
|---|---|---|---|
| `EXHALE_ELEM_DIAG=1` | `[diag 10]`, `[diag 11]`, `[diag 12]` | yes, at iterations 1 and 2 of each solve | PERFORMED |
| `EXHALE_PRECOND_SPECTRUM=1` | `[diag 15]` | yes | PERFORMED |
| `EXHALE_KRYLOV_SIZE_SCAN=1` | `[diag 16]` | yes | PERFORMED |
| `EXHALE_BAND_DIFFERENCE=1` | `[diag 17]` | yes | PERFORMED |
| `EXHALE_FRONT_ROW` | `[diag 13]` | NOT SET in this trial | NOT PERFORMED, by configuration |

**One field rev9 section 7 asks for is not in the instrument.** INSPECTED,
`ritz_values_of_the_preconditioned_operator` writes the largest and smallest
magnitudes, their ratio, the counts below 1E-2 and 1E-4, the number in complex
pairs, the six smallest values and, for three of them, the row-class shares and
the cells of the Ritz vector. It computes and prints **no Ritz residual**
`||A v - theta v||`. The review that asks for them
(`docs/PLAN_20260920_review.md`, the paragraph that also states that the
species and hydrodynamic compressions are not Schur complements) is therefore
NOT satisfied on that item, and this is recorded as an instrument gap and not
substituted by another number. Phase 9 recorded the same absence on the atomic
state.

**`elem_diag_third` was never reached.** It defaults to 60 and
`EXHALE_ELEM_DIAG=1` does not move it (INSPECTED: the reader sets
`elem_diag_third = 60` and only an integer greater than 1 replaces it). Solve 1
ended at iteration 41 and solve 2 was at iteration 60 when the ceiling fired,
so the third scheduled firing of each solve produced the element-row blocks
(`[diag 10]`, `[diag 11]`) of solve 2 and was cut before `[diag 15]`,
`[diag 16]` and `[diag 17]`. **Four complete firings** are on record: outer
iterations 1 and 2 of each of the two solves.

### 6.1 The Ritz values of the preconditioned operator

MEASURED, 200 products each, on the whole operator and on its two compressions.

| firing | whole operator: largest, smallest, ratio | below 1E-2 / below 1E-4 | species compression: largest, smallest, ratio | below 1E-2 / 1E-4 | hydrodynamic compression: largest, smallest, ratio |
|---|---|---|---|---|---|
| solve 1, iteration 1 | 1.51623E+01, 5.81277E-04, **2.60844E+04** | 6 / 0 | 9.91482E-01, 3.42540E-05, 2.89450E+04 | 14 / **1** | 8.48446E+01, 6.04786E-03, 1.40289E+04 |
| solve 1, iteration 2 | 4.07957E+01, 2.98764E-04, **1.36548E+05** | 6 / 0 | 1.00037E+00, 1.24987E-04, 8.00378E+03 | 14 / 0 | 1.14849E+00, 1.02325E-02, 1.12239E+02 |
| solve 2, iteration 1 | 1.57497E+01, 4.77846E-04, **3.29597E+04** | 6 / 0 | 9.91480E-01, 4.30001E-05, 2.30576E+04 | 14 / **1** | 6.38282E+01, 6.32510E-03, 1.00913E+04 |
| solve 2, iteration 2 | 2.80642E+01, 3.95268E-04, **7.10004E+04** | 6 / 0 | 1.00036E+00, 1.27767E-04, 7.82958E+03 | 14 / 0 | 1.14793E+00, 1.02170E-02, 1.12355E+02 |

**Where the near-null directions live.** The routine also prints the row-class
shares of the three smallest Ritz vectors and the cells they concentrate in.
MEASURED at solve 1 iteration 1: the smallest Ritz vector (magnitude
5.81277E-04) carries **0.9618 in the carrier rows**, 0.0305 in the mass rows
and 0.0076 in the energy rows, and nine tenths of it sits in **2 cells, cells 1
and 2**; the third smallest has the same shares and the same two cells. The
second smallest is the one direction that spreads, over 195 cells with cells
431, 435, 430, 490 the largest.

**Two things this does not license.** First, the review that asks for this
measurement states that the species and hydrodynamic compressions **are not
Schur complements**: they omit the cross-block feedback, so a well-conditioned
compression is no statement about the coupled system
(`docs/PLAN_20260920_review.md`, the paragraph beginning "Report Ritz
residuals"). Here it does not have to be invoked as a caveat, because neither
compression is well conditioned at the two iterations where the whole operator
is worst: the species one carries a value at 3.4E-05 and the hydrodynamic one a
ratio of 1.4E+04. Second, **the Ritz RESIDUALS the same paragraph asks for are
not produced by the instrument** (section 6, the last paragraph), so how well
each printed value approximates an eigenvalue is unmeasured here.

### 6.2 The same linear system over a ladder of subspace sizes

MEASURED. `reduced` is the residual of the reduced problem, `TRUE` the
independently evaluated relative residual of the step that cycle returned.

| firing | 40 | 80 | 160 | 320 | 160 orthogonalized twice | 320 orthogonalized twice |
|---|---|---|---|---|---|---|
| solve 1, it 1, reduced | 8.70118E-01 | 8.56542E-01 | 7.62250E-01 | 6.67633E-01 | 7.62251E-01 | 6.68078E-01 |
| solve 1, it 1, TRUE | 8.72736E-01 | 1.00662E+00 | 3.71352E+01 | **7.98372E+01** | 3.71311E+01 | 7.99112E+01 |
| solve 1, it 2, reduced | 8.61796E-01 | 7.96900E-01 | 7.66691E-01 | 7.62010E-01 | 7.66617E-01 | 7.61930E-01 |
| solve 1, it 2, TRUE | 9.44588E-01 | 4.26306E+00 | 4.66223E+00 | 4.83629E+00 | 4.65809E+00 | 4.84352E+00 |
| solve 2, it 1, reduced | 6.94024E-01 | 6.72065E-01 | 5.26881E-01 | 4.24483E-01 | 5.26883E-01 | 4.22292E-01 |
| solve 2, it 1, TRUE | 6.97058E-01 | 7.93023E-01 | 9.65518E+01 | **2.15531E+02** | 9.65776E+01 | 2.20254E+02 |
| solve 2, it 2, reduced | 6.88963E-01 | 6.56777E-01 | 6.27632E-01 | 5.59022E-01 | 6.27614E-01 | 5.61292E-01 |
| solve 2, it 2, TRUE | 7.02381E-01 | 2.28770E+00 | 4.05383E+01 | 4.49141E+01 | 2.73270E+00 | 4.05383E+01 |

**Two statements, and they are separate.**

1. **The cycle stagnates at every size.** Eight times the subspace and the
   reduced residual falls from 0.87 to 0.67, from 0.86 to 0.76, from 0.69 to
   0.42 and from 0.69 to 0.56. In the words of the routine's own header, that
   is "an operator with a part the preconditioner does not touch, and no amount
   of subspace removes it". The solve's own cycle says the same at every one of
   the 100 iterations: 40 products of 40, the subspace exhausted, the relative
   residual reached between 5.836E-01 and 8.782E-01 against the 1.00E-01 it was asked
   for. **Not once in 100 iterations did the linear solve reach its tolerance.**
2. **The step the cycle returns gets WORSE as the subspace grows.** The true
   residual of the returned step rises from 0.873 to 79.8, from 0.697 to 215.5.
   Reorthogonalizing the basis changes nothing (7.62250E-01 against 7.62251E-01
   at 160, 5.26881E-01 against 5.26883E-01), so this is not a lost
   orthogonality: the reduced problem and the true operator disagree, and the
   disagreement grows with the basis. The solve's own cycle reports the same
   object as a "relative gap of the Arnoldi image" of 1.469E-01 at 40 products.

### 6.3 What the band omits of the Jacobian

MEASURED, on the three smallest Ritz directions and on the first two Arnoldi
directions of the solve's own right-hand side. The entry is
`||(A - A_band) v|| / ||A v||`.

| firing | Ritz 1 | Ritz 2 | Ritz 3 | Arnoldi 1 | Arnoldi 2 |
|---|---|---|---|---|---|
| solve 1, iteration 1 | **2.38158E+00** | 1.51188E+00 | 2.39190E+00 | 6.77423E-04 | 1.69946E-03 |
| solve 1, iteration 2 | **2.33035E+00** | 1.05792E+00 | 6.46478E-04 | 6.77417E-04 | 1.63412E-03 |
| solve 2, iteration 1 | **2.38195E+00** | 1.25048E-03 | 6.58109E-04 | 6.61626E-03 | 6.65371E-03 |
| solve 2, iteration 2 | **2.36684E+00** | 9.70412E-01 | 2.34845E+00 | 6.61200E-03 | 6.89809E-03 |

On the direction the spectrum names as the worst, the band misses **more than
twice the whole action**; on the directions the right-hand side actually
excites it misses 0.07 to 0.7 per cent. For comparison, on the atomic
checkpoint with the trust region as a variant the band omitted 0.26 per cent on
the smallest Ritz direction (READ, rev9 section 7.1b).

**Where the omission sits, and the limit that applies to it.** MEASURED at
solve 1 iteration 1 on the smallest Ritz direction: 0.8078 of the difference is
in the **mass rows**, 0.1921 in the carrier rows and 0.0001 in the energy rows,
and nine tenths of it sits in **one cell, cell 1**, although the direction
itself carries only 0.0305 in the mass rows. The columns the routine then reads
say the same from the other side: for the carrier H2 column of cell 1 and the
mass column of cell 1 the part OUTSIDE the band is exactly 0.00000 and the
whole difference is inside it, against the band's own entries, at 0.00071 of
the column.

That attribution lands on the object phase 6b named as the limit of the coupled
action: an arc-independent residue of about 6E+04 of the row scales in the base
cells' mass rows, which leaves the action on a localized base mass direction
resolved to no better than 1.4E-02 (READ,
`docs/lhs1140b_coupled_action_20260921.md` section 1). **This document therefore
does not turn any conclusion on that attribution.** What survives the limit is
the part that does not depend on the base mass rows: the size of the band
difference on the Ritz directions is 1.5 to 2.4, two decades above the 1.4E-02
the action is resolved to on any base mass direction, and the same measurement
on the Arnoldi directions is three decades smaller, so the band's failure is
direction-selective and is not the action's own noise floor.

---

## 7. Stage D, the ledger (rev9 section 7.2)

A recording protocol, not an interpretation. Where a field does not apply or is
not in the run's output it says so; no field carries a zero that could be read
as convergence.

### 7.1 The fields, and where each comes from

| field (rev9 section 7.2) | on this route |
|---|---|
| the state identity and the active unknowns | section 2; 5 unknowns a cell, 2500 entries, the He element row and the H2 carrier row, 29 colors (MEASURED, the solve's `cost:` line) |
| the unscaled row measures AND the tolerance-normalized distances, with the cell of each | in the run's own output: the certification block per pass (section 5) and the `judged rows` line per iterate (section 5.1) |
| the actual scaled state displacement, the step length or trust radius | in the output on this route: `||s||` and `delta` of the `[TR] it` line, section 7.2 below. On the three-unknown route phase 9 had to record this as NOT APPLICABLE |
| the predicted and the observed reduction under fixed comparison weights | `pred` and `actual` of the same line, with the solver's own note that the merit is compared frozen-mode against trial-mode and their ratio printed (1.000 at every iteration measured here) |
| the reduced Krylov residual and the independently evaluated residual of the step it returned | the reduced one is in the `Krylov leg` line of every iteration; the independently evaluated one is NOT produced per iteration, only inside the subspace scan at the four scheduled firings (section 6.2) |
| bound hits with their projection displacement | MEASURED: `unknowns held at a bound 0` at every `[diag 5]` line, and the solve's own summary `most species unknowns held on an active bound by one step: 0`, `trials written onto the faces of the species box: 0`, `adopted carriers at or below that floor: 0`, `adopted carriers above that budget: 0`. **No bound was active at any iteration of either solve** |
| limiter and contact events | `restoration phase: trials restored 0, cell-steps 0`; `carrier cells whose element budget face was raised to the iterate: 0`; `shared element rows: cell-rows whose iterate stood outside its own budget 0, rows dropped as already spanned 0` (MEASURED, solve 1's summary). Solve 2 did not print its summary and those counts are NOT AVAILABLE for it |
| any change of the composition or radiation closure | the closure map is reported at the entry of each solve: solve 1 `exited at increment 1.458E-13 after 1 passes; the reference ran 4 passes`, base amplification 2.065E+04; solve 2 `1.789E-13 after 4 passes; the reference ran 7`, base amplification 6.457E+04 (MEASURED) |
| the termination status, distinguishing an inner return, a completed outer route and an external interruption | solve 1: an inner return, `info = 2`, STAGNATED at a judged distance 5.957E+03 with the best at 8.150E+01. Outer pass 1: completed, `gate NOT met: carrier row 2.983E-03 >= 1.000E-03`. Solve 2 and the run: an EXTERNAL INTERRUPTION, the 18000 s safety ceiling, `exit = 124`. No counter was reset by anything but the solver's own restarts, which are named in section 7.3 |

### 7.2 Every trust-region iteration of the trial

MEASURED, one row per outer iteration, read from the single line the solver
writes for each, accepted or not. `delta` is the trust radius, `||s||` the
scaled length of the step actually taken, `cuts` the number of radius
reductions inside the iteration, and the Krylov column the leg the dogleg was
built from. The full table with the Krylov column is in the artifacts
(`trial_tr_ledger.md`); reproduced here without that column, which reads
`40 of 40, the subspace was exhausted with a finite approximate step` at every
one of the 100 rows.

| solve | it | delta | pred | actual | ratio | ||s|| | cuts | model_ok | accepted | relative residual of the Krylov leg | why no step |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 1.000e-04 | 5.960e-07 | 5.882e-07 | 0.987 | 5.000e-05 | 0 | T | T | 8.701E-01 | -- |
| 1 | 2 | 2.000e-04 | 4.608e-07 | 5.170e-07 | 1.122 | 1.000e-04 | 0 | T | T | 8.618E-01 | -- |
| 1 | 3 | 4.000e-04 | 7.282e-06 | 7.563e-06 | 1.039 | 2.000e-04 | 0 | T | T | 8.338E-01 | -- |
| 1 | 4 | 4.000e-04 | 1.445e-06 | 9.274e-07 | 0.642 | 4.000e-04 | 0 | T | T | 7.756E-01 | -- |
| 1 | 5 | 8.000e-04 | 1.906e-06 | 1.857e-06 | 0.974 | 4.000e-04 | 0 | T | T | 7.720E-01 | -- |
| 1 | 6 | 2.000e-04 | 7.828e-06 | -3.630e-05 | -4.637 | 8.000e-04 | 0 | F | F | 7.189E-01 | the reduction ratio is below eta |
| 1 | 7 | 2.000e-04 | 7.248e-07 | 4.838e-07 | 0.667 | 2.000e-04 | 0 | T | T | 7.755E-01 | -- |
| 1 | 8 | 4.000e-04 | 9.473e-07 | 9.630e-07 | 1.017 | 2.000e-04 | 0 | T | T | 7.728E-01 | -- |
| 1 | 9 | 1.000e-04 | 4.040e-06 | -1.540e-05 | -3.812 | 4.000e-04 | 0 | F | F | 7.191E-01 | the reduction ratio is below eta |
| 1 | 10 | 1.000e-04 | 3.634e-07 | 2.482e-07 | 0.683 | 1.000e-04 | 0 | T | T | 7.755E-01 | -- |
| 1 | 11 | 2.000e-04 | 4.812e-07 | 4.958e-07 | 1.030 | 1.000e-04 | 0 | T | T | 7.727E-01 | -- |
| 1 | 12 | 5.000e-05 | 1.974e-06 | -6.243e-06 | -3.163 | 2.000e-04 | 0 | F | F | 7.192E-01 | the reduction ratio is below eta |
| 1 | 13 | 1.250e-05 | 1.833e-07 | 1.265e-07 | 0.690 | 5.000e-05 | 0 | F | F | 7.758E-01 | the model and the true slope differ in size |
| 1 | 14 | 2.500e-05 | 6.169e-08 | 6.915e-08 | 1.121 | 1.250e-05 | 0 | T | T | 8.631E-01 | -- |
| 1 | 15 | 5.000e-05 | 5.320e-07 | 5.857e-07 | 1.101 | 2.500e-05 | 0 | T | T | 8.340E-01 | -- |
| 1 | 16 | 5.000e-05 | 1.820e-07 | 1.239e-07 | 0.681 | 5.000e-05 | 0 | T | T | 7.760E-01 | -- |
| 1 | 17 | 1.000e-04 | 2.326e-07 | 2.449e-07 | 1.053 | 5.000e-05 | 0 | T | T | 7.734E-01 | -- |
| 1 | 18 | 2.500e-05 | 9.887e-07 | -3.036e-06 | -3.071 | 1.000e-04 | 0 | F | F | 7.200E-01 | the reduction ratio is below eta |
| 1 | 19 | 2.500e-05 | 9.156e-08 | 6.282e-08 | 0.686 | 2.500e-05 | 0 | T | T | 7.763E-01 | -- |
| 1 | 20 | 5.000e-05 | 1.173e-07 | 1.231e-07 | 1.049 | 2.500e-05 | 0 | T | T | 7.731E-01 | -- |
| 1 | 21 | 1.250e-05 | 4.689e-07 | -1.336e-06 | -2.850 | 5.000e-05 | 0 | F | F | 7.194E-01 | the reduction ratio is below eta |
| 1 | 22 | 2.500e-05 | 4.886e-08 | 4.798e-08 | 0.982 | 1.250e-05 | 0 | T | T | 7.750E-01 | -- |
| 1 | 23 | 5.000e-05 | 1.174e-07 | 1.319e-07 | 1.123 | 2.500e-05 | 0 | T | T | 7.725E-01 | -- |
| 1 | 24 | 1.250e-05 | 4.614e-07 | -1.233e-06 | -2.672 | 5.000e-05 | 0 | F | F | 7.191E-01 | the reduction ratio is below eta |
| 1 | 25 | 3.125e-06 | 4.892e-08 | 4.816e-08 | 0.984 | 1.250e-05 | 0 | F | F | 7.753E-01 | the model and the true slope differ in sign |
| 1 | 26 | 7.812e-07 | 1.445e-08 | 1.644e-08 | 1.138 | 3.125e-06 | 0 | F | F | 8.619E-01 | the model and the true slope differ in sign |
| 1 | 27 | 1.562e-06 | 3.318e-07 | 3.729e-07 | 1.124 | 7.812e-07 | 0 | T | T | 8.561E-01 | -- |
| 1 | 28 | 3.125e-06 | 1.845e-08 | 1.828e-08 | 0.991 | 1.562e-06 | 0 | T | T | 8.699E-01 | -- |
| 1 | 29 | 7.812e-07 | 1.441e-08 | 1.629e-08 | 1.130 | 3.125e-06 | 0 | F | F | 8.618E-01 | the model and the true slope differ in sign |
| 1 | 30 | 1.562e-06 | 3.399e-07 | 3.817e-07 | 1.123 | 7.812e-07 | 0 | T | T | 8.561E-01 | -- |
| 1 | 31 | 3.125e-06 | 1.825e-08 | 1.821e-08 | 0.998 | 1.562e-06 | 0 | T | T | 8.698E-01 | -- |
| 1 | 32 | 6.250e-06 | 1.433e-08 | 1.635e-08 | 1.141 | 3.125e-06 | 0 | T | T | 8.616E-01 | -- |
| 1 | 33 | 1.250e-05 | 1.595e-07 | 1.693e-07 | 1.062 | 6.250e-06 | 0 | T | T | 8.350E-01 | -- |
| 1 | 34 | 3.125e-06 | 4.953e-08 | 4.824e-08 | 0.974 | 1.250e-05 | 0 | F | F | 7.756E-01 | the ray difference is below its cancellation floor |
| 1 | 35 | 7.812e-07 | 1.443e-08 | 1.650e-08 | 1.143 | 3.125e-06 | 0 | F | F | 8.617E-01 | the ray difference is below its cancellation floor |
| 1 | 36 | 1.562e-06 | 5.289e-08 | 5.331e-08 | 1.008 | 7.813e-07 | 0 | T | T | 8.781E-01 | -- |
| 1 | 37 | 3.125e-06 | 1.818e-08 | 1.804e-08 | 0.992 | 1.562e-06 | 0 | T | T | 8.698E-01 | -- |
| 1 | 38 | 7.812e-07 | 1.444e-08 | 1.663e-08 | 1.152 | 3.125e-06 | 0 | F | F | 8.617E-01 | the model and the true slope differ in sign |
| 1 | 39 | 1.562e-06 | 5.298e-08 | 5.346e-08 | 1.009 | 7.812e-07 | 0 | T | T | 8.782E-01 | -- |
| 1 | 40 | 3.125e-06 | 1.815e-08 | 1.793e-08 | 0.988 | 1.563e-06 | 0 | T | T | 8.698E-01 | -- |
| 1 | 41 | 7.813e-07 | 1.444e-08 | 1.633e-08 | 1.131 | 3.125e-06 | 0 | F | F | 8.617E-01 | the model and the true slope differ in sign |
| 2 | 1 | 1.000e-04 | 1.919e-05 | 1.952e-05 | 1.017 | 5.000e-05 | 0 | T | T | 6.940E-01 | -- |
| 2 | 2 | 2.000e-04 | 2.827e-06 | 3.120e-06 | 1.103 | 1.000e-04 | 0 | T | T | 6.890E-01 | -- |
| 2 | 3 | 4.000e-04 | 2.879e-06 | 3.206e-06 | 1.114 | 2.000e-04 | 0 | T | T | 6.790E-01 | -- |
| 2 | 4 | 8.000e-04 | 3.862e-06 | 1.481e-05 | 3.835 | 4.000e-04 | 0 | T | T | 6.308E-01 | -- |
| 2 | 5 | 1.600e-03 | 6.471e-06 | 2.095e-05 | 3.238 | 8.000e-04 | 0 | T | T | 6.202E-01 | -- |
| 2 | 6 | 4.000e-04 | 2.509e-05 | -2.831e-04 | -11.284 | 1.600e-03 | 0 | F | F | 5.836E-01 | the reduction ratio is below eta |
| 2 | 7 | 8.000e-04 | 3.820e-06 | 1.454e-05 | 3.805 | 4.000e-04 | 0 | T | T | 6.323E-01 | -- |
| 2 | 8 | 1.600e-03 | 6.452e-06 | 2.001e-05 | 3.102 | 8.000e-04 | 0 | T | T | 6.220E-01 | -- |
| 2 | 9 | 4.000e-04 | 2.444e-05 | -2.702e-04 | -11.054 | 1.600e-03 | 0 | F | F | 5.843E-01 | the reduction ratio is below eta |
| 2 | 10 | 8.000e-04 | 3.772e-06 | 1.423e-05 | 3.774 | 4.000e-04 | 0 | T | T | 6.338E-01 | -- |
| 2 | 11 | 1.600e-03 | 6.414e-06 | 1.910e-05 | 2.979 | 8.000e-04 | 0 | T | T | 6.232E-01 | -- |
| 2 | 12 | 4.000e-04 | 2.394e-05 | -2.643e-04 | -11.042 | 1.600e-03 | 0 | F | F | 5.860E-01 | the reduction ratio is below eta |
| 2 | 13 | 8.000e-04 | 3.733e-06 | 1.405e-05 | 3.764 | 4.000e-04 | 0 | T | T | 6.351E-01 | -- |
| 2 | 14 | 1.600e-03 | 6.360e-06 | 1.825e-05 | 2.870 | 8.000e-04 | 0 | T | T | 6.245E-01 | -- |
| 2 | 15 | 4.000e-04 | 2.395e-05 | -2.448e-04 | -10.223 | 1.600e-03 | 0 | F | F | 5.869E-01 | the reduction ratio is below eta |
| 2 | 16 | 8.000e-04 | 3.701e-06 | 1.375e-05 | 3.715 | 4.000e-04 | 0 | T | T | 6.364E-01 | -- |
| 2 | 17 | 1.600e-03 | 6.375e-06 | 1.740e-05 | 2.730 | 8.000e-04 | 0 | T | T | 6.261E-01 | -- |
| 2 | 18 | 4.000e-04 | 2.408e-05 | -2.354e-04 | -9.776 | 1.600e-03 | 0 | F | F | 5.885E-01 | the reduction ratio is below eta |
| 2 | 19 | 8.000e-04 | 3.665e-06 | 1.352e-05 | 3.689 | 4.000e-04 | 0 | T | T | 6.377E-01 | -- |
| 2 | 20 | 1.600e-03 | 6.293e-06 | 1.635e-05 | 2.597 | 8.000e-04 | 0 | T | T | 6.269E-01 | -- |
| 2 | 21 | 4.000e-04 | 2.313e-05 | -2.119e-04 | -9.164 | 1.600e-03 | 0 | F | F | 5.908E-01 | the reduction ratio is below eta |
| 2 | 22 | 8.000e-04 | 3.639e-06 | 1.329e-05 | 3.651 | 4.000e-04 | 0 | T | T | 6.389E-01 | -- |
| 2 | 23 | 1.600e-03 | 6.287e-06 | 1.558e-05 | 2.479 | 8.000e-04 | 0 | T | T | 6.285E-01 | -- |
| 2 | 24 | 4.000e-04 | 2.438e-05 | -2.359e-04 | -9.677 | 1.600e-03 | 0 | F | F | 5.900E-01 | the reduction ratio is below eta |
| 2 | 25 | 8.000e-04 | 3.625e-06 | 1.311e-05 | 3.618 | 4.000e-04 | 0 | T | T | 6.400E-01 | -- |
| 2 | 26 | 1.600e-03 | 6.319e-06 | 1.504e-05 | 2.380 | 8.000e-04 | 0 | T | T | 6.303E-01 | -- |
| 2 | 27 | 4.000e-04 | 2.428e-05 | -2.111e-04 | -8.694 | 1.600e-03 | 0 | F | F | 5.915E-01 | the reduction ratio is below eta |
| 2 | 28 | 8.000e-04 | 3.594e-06 | 1.285e-05 | 3.574 | 4.000e-04 | 0 | T | T | 6.411E-01 | -- |
| 2 | 29 | 1.600e-03 | 6.266e-06 | 1.408e-05 | 2.247 | 8.000e-04 | 0 | T | T | 6.314E-01 | -- |
| 2 | 30 | 4.000e-04 | 2.449e-05 | -2.345e-04 | -9.579 | 1.600e-03 | 0 | F | F | 5.922E-01 | the reduction ratio is below eta |
| 2 | 31 | 8.000e-04 | 3.584e-06 | 1.268e-05 | 3.537 | 4.000e-04 | 0 | T | T | 6.421E-01 | -- |
| 2 | 32 | 1.600e-03 | 6.133e-06 | 1.292e-05 | 2.106 | 8.000e-04 | 0 | T | T | 6.314E-01 | -- |
| 2 | 33 | 4.000e-04 | 2.402e-05 | -2.129e-04 | -8.861 | 1.600e-03 | 0 | F | F | 5.936E-01 | the reduction ratio is below eta |
| 2 | 34 | 8.000e-04 | 3.568e-06 | 1.243e-05 | 3.483 | 4.000e-04 | 0 | T | T | 6.432E-01 | -- |
| 2 | 35 | 1.600e-03 | 6.137e-06 | 1.216e-05 | 1.982 | 8.000e-04 | 0 | T | T | 6.329E-01 | -- |
| 2 | 36 | 4.000e-04 | 2.340e-05 | -2.032e-04 | -8.683 | 1.600e-03 | 0 | F | F | 5.952E-01 | the reduction ratio is below eta |
| 2 | 37 | 8.000e-04 | 3.564e-06 | 1.221e-05 | 3.427 | 4.000e-04 | 0 | T | T | 6.440E-01 | -- |
| 2 | 38 | 1.600e-03 | 6.117e-06 | 1.175e-05 | 1.921 | 8.000e-04 | 0 | T | T | 6.340E-01 | -- |
| 2 | 39 | 4.000e-04 | 2.394e-05 | -1.794e-04 | -7.495 | 1.600e-03 | 0 | F | F | 5.956E-01 | the reduction ratio is below eta |
| 2 | 40 | 8.000e-04 | 3.552e-06 | 1.195e-05 | 3.364 | 4.000e-04 | 0 | T | T | 6.449E-01 | -- |
| 2 | 41 | 1.600e-03 | 6.079e-06 | 1.071e-05 | 1.761 | 8.000e-04 | 0 | T | T | 6.349E-01 | -- |
| 2 | 42 | 4.000e-04 | 2.266e-05 | -1.643e-04 | -7.252 | 1.600e-03 | 0 | F | F | 5.977E-01 | the reduction ratio is below eta |
| 2 | 43 | 8.000e-04 | 3.539e-06 | 1.168e-05 | 3.301 | 4.000e-04 | 0 | T | T | 6.458E-01 | -- |
| 2 | 44 | 1.600e-03 | 6.112e-06 | 1.032e-05 | 1.689 | 8.000e-04 | 0 | T | T | 6.364E-01 | -- |
| 2 | 45 | 4.000e-04 | 2.391e-05 | -1.784e-04 | -7.462 | 1.600e-03 | 0 | F | F | 5.975E-01 | the reduction ratio is below eta |
| 2 | 46 | 8.000e-04 | 3.542e-06 | 1.145e-05 | 3.232 | 4.000e-04 | 0 | T | T | 6.465E-01 | -- |
| 2 | 47 | 1.600e-03 | 6.048e-06 | 9.425e-06 | 1.558 | 8.000e-04 | 0 | T | T | 6.368E-01 | -- |
| 2 | 48 | 4.000e-04 | 2.292e-05 | -1.644e-04 | -7.172 | 1.600e-03 | 0 | F | F | 5.980E-01 | the reduction ratio is below eta |
| 2 | 49 | 8.000e-04 | 3.542e-06 | 1.117e-05 | 3.155 | 4.000e-04 | 0 | T | T | 6.474E-01 | -- |
| 2 | 50 | 1.600e-03 | 5.958e-06 | 8.066e-06 | 1.354 | 8.000e-04 | 0 | T | T | 6.372E-01 | -- |
| 2 | 51 | 4.000e-04 | 2.310e-05 | -1.493e-04 | -6.462 | 1.600e-03 | 0 | F | F | 5.986E-01 | the reduction ratio is below eta |
| 2 | 52 | 8.000e-04 | 3.542e-06 | 1.088e-05 | 3.072 | 4.000e-04 | 0 | T | T | 6.482E-01 | -- |
| 2 | 53 | 1.600e-03 | 5.977e-06 | 7.642e-06 | 1.279 | 8.000e-04 | 0 | T | T | 6.384E-01 | -- |
| 2 | 54 | 4.000e-04 | 2.307e-05 | -1.344e-04 | -5.826 | 1.600e-03 | 0 | F | F | 6.000E-01 | the reduction ratio is below eta |
| 2 | 55 | 8.000e-04 | 3.532e-06 | 1.051e-05 | 2.977 | 4.000e-04 | 0 | T | T | 6.490E-01 | -- |
| 2 | 56 | 1.600e-03 | 5.984e-06 | 7.058e-06 | 1.179 | 8.000e-04 | 0 | T | T | 6.394E-01 | -- |
| 2 | 57 | 4.000e-04 | 2.348e-05 | -1.348e-04 | -5.740 | 1.600e-03 | 0 | F | F | 6.007E-01 | the reduction ratio is below eta |
| 2 | 58 | 8.000e-04 | 3.520e-06 | 1.017e-05 | 2.890 | 4.000e-04 | 0 | T | T | 6.497E-01 | -- |
| 2 | 59 | 1.600e-03 | 5.967e-06 | 6.754e-06 | 1.132 | 8.000e-04 | 0 | T | T | 6.400E-01 | -- |

### 7.3 What the ledger says, counted

MEASURED, from the 100 rows above.

| | solve 1 | solve 2, to the ceiling |
|---|---|---|
| outer iterations | 41 | 59 |
| steps accepted | 27 | 41 |
| steps refused | **14** | **18** |
| the reduction ratio is below eta | 6 | 18 |
| the model and the true slope differ in sign | 5 | 0 |
| the model and the true slope differ in size | 1 | 0 |
| the ray difference is below its cancellation floor | 2 | 0 |
| iterations with `model_ok = F` | 14 | 18 |
| iterations with a radius cut inside them | 0 | 0 |
| trust radius, first / largest / smallest / last | 1.000E-04 / 8.000E-04 / **7.812E-07** / 7.813E-07 | 1.000E-04 / 1.600E-03 / 1.000E-04 / 1.600E-03 |
| accepted `||s||`, largest / smallest | 4.000E-04 / **7.812E-07** | 8.000E-04 / 5.000E-05 |
| reduction ratio, smallest / largest | **-4.637** / 1.152 | **-11.284** / 3.835 |
| iterations with `|ratio - 1| > 0.2` | 12 of 41 | **54 of 59** |
| unknowns held on a bound | 0 at every iteration | 0 at every iteration |

Two events of the solver's own are part of the record and are named rather
than smoothed:

- **the monotone restart.** After outer iteration 20 of solve 1 the solver
  wrote `neither the judged distance nor the merit has improved in 20
  iterations; restarting from the best iterate (judged distance 1.280E+05,
  best 8.150E+01, merit over the window 9.704E-02, ||R||= 1.872E-08) with a
  monotone line search`. The best iterate is the ENTRY state: its judged
  distance 8.150E+01 and its `||R||` 1.872E-08 are the entry's own. So the
  first twenty iterations were discarded by the solver itself;
- **the stagnation.** After iteration 41, `STAGNATED: neither the judged
  distance nor the merit improved in 40 iterations, monotone search included
  -- stopping at a judged distance 5.957E+03, best 8.150E+01, merit over the
  window 9.718E-02, ||R||= 1.872E-08`, and the state handed back is again the
  entry state.

**Where the step came from.** MEASURED, solve 1's summary:
`[TR] steps that were the Cauchy point alone: 0` and
`[TR] dogleg trials that left the approximate-gradient leg out of the step: 41
(shorter than cauchy_leg_min_fraction of the radius)`. At outer iteration 1 the
Cauchy step was 5.090E-06 of the radius. **Every step of solve 1 was therefore
the truncated Krylov direction scaled to the radius, that is, the direction
whose own linear residual was 5.836E-01 to 8.782E-01 against the 1.00E-01 asked
for.**

---

## 8. The four questions the phase asks

### 8.1 Does removing the movement bound remove the drift?

**No, and the reason is not that the bound was helping.** The bound is gone:
the carrier relaxation does not run at all on this route, so the 25 to 30
refused trials and the saturation at cells 418 to 426 that every control pass
carried have no counterpart here. What the drift is replaced by is a solve that
returns its own entry iterate.

Pass for pass, on the gated rows, with the trial's one completed pass beside
the control's six (READ for the control, MEASURED for the trial):

| pass | carrier gated, control | carrier gated, trial | elemental gated, control | elemental gated, trial |
|---|---|---|---|---|
| entry | 8.150E-04 | 8.150E-04 | 2.503E-04 | 2.503E-04 |
| 1 | 8.150E-04 | **8.150E-04** | 2.503E-04 | **2.503E-04** |
| 2 | 6.290E-03 | NOT REACHED | 4.532E-04 | NOT REACHED |
| 3 | 6.258E-03 | NOT REACHED | 6.069E-04 | NOT REACHED |
| 4 | 6.388E-03 | NOT REACHED | 6.961E-04 | NOT REACHED |
| 5 | 2.950E-03 | NOT REACHED | 6.968E-04 | NOT REACHED |
| 6 | 2.893E-03 | NOT REACHED | 7.063E-04 | NOT REACHED |

At the only pass both trajectories have, the two routes stand at the same two
numbers, the control having spent 8.6 s and the trial 7171.8 s. The last
ITERATE the trial reached inside its unfinished second solve carries carrier
7.637E-04 and element 1.063E-04 (section 5.1), which is below the control's
pass-2 values, but it is an iterate of a state whose hydrodynamic judged
distance is 3.870E+06 and it is not a certificate.

**So the question is answered for the one pass measured and left open beyond
it**, and the answer for that pass is that the coupled route does not drift
because it does not move. The comparison is of the COMBINED configuration and
not of the bound alone.

**One asymmetry of the comparison has to be said, because it is the reason the
two solves are not the same problem.** The merit the solve minimizes carries the
registered rows. MEASURED at the same entry state: on the control's route
`||Fs||2 = 8.45E-11`, so its first solve had nothing to do and returned
immediately; on the coupled route `||Fs||2 = 9.720E-02`, of which the carrier
H2 row holds a share of 0.7701 and the element He row 0.2299, with the mass,
momentum and energy rows holding 0.0000 between them. The coupled solve is
therefore asked to solve the species balances that the control hands to a
bounded relaxation, and the measurement of this trial is that it does not solve
them: the worst row of the merit is `element He of cell 2` at every one of the
40 iteration lines of solve 1, and at the two iterations where the diagnostic
reads its value it stands at 4.656E-02 and 4.655E-02.

### 8.2 Does the refreshed full certificate of any pass certify the state?

**No.** The one refreshed full certification the trial produced, at the end of
outer pass 1, reads (MEASURED):

```
NOT CERTIFIED: 2 entry/entries of the inventory refuse it
  carrier balance H2: gated row measure 8.150E-04 above 1.0E-05 at cell 499 (a wind cell)
  elemental transport He/H partition: gated row measure 2.503E-04 above 1.0E-05 at cell 280 (a wind cell)
```

with the hydrodynamic mass row 1.107E-08 at cell 1, the momentum row 5.248E-13
at cell 1 and the energy row 5.803E-09 at cell 3, all within their tolerances,
the He 2 3S level row and the eliminated-species closure within, and 0 cells
without a chemical root. The column maxima of the two refusing rows, which the
certification reports without gating, are 3.458E-02 at cell 1 for the carrier
and 4.859E-02 at cell 2 for the element. The state was NOT published: the run
was stopped by the ceiling before it wrote one, and its `output/` directory
holds only the diagnostic files and the `_IC` pair it was given.

### 8.3 What the Ritz spectrum, the band difference and the subspace scan say

Section 6 carries the tables. The three statements they support, in the order
of rev9 section 7's order of inference, are:

1. **The linear model is never solved on this route.** At every one of the 100
   outer iterations the Krylov cycle used its whole subspace, 40 products of
   40, and returned a relative residual between 5.836E-01 and 8.782E-01 against the
   1.00E-01 it was asked for. The subspace scan says more subspace does not fix
   it: eight times the size buys a fall from 0.87 to 0.67, and the true
   residual of the step the larger cycles return RISES, to 79.8 at 320 products
   at one firing and to 215.5 at another, with reorthogonalization changing
   nothing to five digits.
2. **The spectrum says where the obstruction sits and the band says what it
   is.** The preconditioned operator carries a ratio of 2.6E+04 to 1.4E+05 with
   its smallest value at 3.0E-04 to 5.8E-04, and the smallest Ritz vector
   carries 0.96 of itself in the carrier rows concentrated in cells 1 and 2. On
   exactly those directions the banded preconditioner misses more than twice
   the whole action (2.33 to 2.39 at all four firings), while on the two
   Arnoldi directions of the right-hand side it misses 0.07 to 0.7 per cent.
   The routine's own column reading says the omission is not outside the band
   at all: for the carrier and mass columns of cell 1 the part outside the band
   is exactly zero and the whole difference is the coloring contaminating
   entries inside it.
3. **Where that localization may and may not be pressed.** The band difference
   attributes 0.81 of itself to the mass rows of cell 1, and phase 6b's named
   limit is that the coupled action on a localized base mass direction is
   resolved to no better than 1.4E-02. This document therefore reads the
   attribution as a coincidence of location and draws no conclusion from it.
   What stands without it is the SIZE, 1.5 to 2.4 against a resolution of
   1.4E-02, and the contrast with the Arnoldi directions three decades below.
4. **Against the atomic state.** On the atomic checkpoint with the trust region
   as a named variant, phase 9 recorded a Ritz ratio of 3.96E+03 with nothing
   below 1E-4, the band omitting 0.26 per cent, the subspace scan identical at
   40, 80, 160 and 320 after three products, and observed reduction equal to
   the prediction to 1.000 at eight consecutive iterations (READ, rev9 section
   7.1b). Every one of those four is qualitatively different here. The two
   states differ in the case, in the physics and in the unknown set at once, so
   this is a contrast and not an attribution.

### 8.4 Does the factor 40 of phase 6b show up as rejected steps or shrinking radii?

**Yes, as both, and the ledger of section 7.2 is the record.** Phase 6b
measured, on this very state, that the banded model proposes a step with
`||D^-1 s|| = 5.488E+01` against 9.25E-10 uncoupled, and that the full action
along that proposed step is a factor 40 from what the model predicts, where the
two agree to three digits on the uncoupled route (READ,
`docs/lhs1140b_coupled_action_20260921.md` sections 1 and 11).

In the solve:

- **32 of the 100 outer iterations refused their step**, 14 in solve 1 and 18
  in solve 2. The dominant reason is the prediction test itself: `the reduction
  ratio is below eta` 24 times, with observed-to-predicted ratios of -4.637,
  -3.812, -3.163, -2.850 in solve 1 and down to -11.284 in solve 2, that is,
  the merit ROSE by three to eleven times the fall the model predicted. Five
  further refusals in solve 1 are `the model and the true slope differ in sign`
  and one `differ in size`, which is the same disagreement read on the
  directional derivative instead of on the step;
- **the trust radius collapses by three decades in solve 1**, from 8.000E-04 to
  7.812E-07, and the accepted step with it, from 4.000E-04 to 7.812E-07. In
  solve 2 it does not collapse but it does not grow either: it cycles between
  1.000E-04 and 1.600E-03 for 59 iterations while `|ratio - 1|` exceeds 0.2 at
  54 of them;
- **no radius cut occurred inside any iteration** (`cuts = 0` at all 100 rows)
  and **no unknown was ever held on a bound**, so neither the inner cut logic
  nor the species box is what limits the step. What limits it is the model.

The chain is therefore closed on the measurement side: the Krylov cycle cannot
solve the banded-preconditioned system (section 8.3 item 1), the step it
returns is the one the dogleg takes (section 7.3), and the model's prediction
for that step is wrong by factors of 3 to 11 often enough to refuse a third of
the iterations. What phase 6b measured on a frozen state is what the solve then
spends its budget on.

---

## 9. What this does not establish

- **Five of the six passes.** The trajectory holds one completed outer pass.
  Everything said about passes 2 to 6 of the trial is "not reached", and the
  ending of the run is an external termination at a safety ceiling.
- **Which of the three changes produced the outcome.** The unknown set, the
  globalization and the removal of the carrier relaxation move together
  (section 3.1). A trust-region control on the three-unknown route, or a
  coupled run with `EXHALE_TRUST_REGION=0`, would separate the first two; rev9
  section 9.4 says such a control is worth running only if the attribution
  becomes necessary, and this document does not claim it is.
- **Whether the same happens without the four diagnostics.** They fire at
  outer iterations 1 and 2 of every solve and are the larger part of the wall
  clock (section 4.1). They adopt nothing and their residual samples are held
  out of the solve's counts (INSPECTED, the subspace-scan routine holds and
  restores the residual evaluation products and the refusal statistics), so
  they are not expected to change the trajectory; that expectation is NOT
  measured here.
- **The cause of the coloring contamination inside the band.** It is measured
  and localized and not explained, and its localization falls where the coupled
  action is only resolved to 1.4E-02 (section 8.3 item 3).
- **The Ritz residuals.** The instrument does not produce them (section 6), so
  how well each printed value approximates an eigenvalue is unmeasured.
- **Anything about the other three M1 cases.** One case was run, as rev9
  section 9.4 asks.

---

## 10. Artifacts

All raw artifacts are in the scratch directory of this trial, outside the
repository:

```
/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/m1_coupled/
```

| path | content |
|---|---|
| `trial/input.inp`, `trial/base.inp` | the configuration actually run, md5 `f510b39051b6bdc8e0bbfe9bd660c1c5` and `1e33d197ead9269661dc6bb0d88c0e26` |
| `trial/output/Hydro_ioniz_IC.txt`, `Ion_species_IC.txt` | the loaded generation, md5 `fd882d8918052069f8ed229dd376f8c4` and `87288b03ad88a92b7b8b8b6c2b1bdd66` |
| `trial/run.log`, `trial/run.err` | the Mode C trajectory, 2009 lines |
| `trial/EXHALE_resolved.out`, `trial/EXHALE_setup.out` | the effective configuration |
| `trial/output/carrier_row_terms.txt`, `diffusion_pass_profiles.txt`, `element_flux_profile.txt` | the pass-1 certification's diagnostic files |
| `launch_trial.sh` | the launch, carrying the complete environment |
| `trial_start.txt`, `trial_exit.txt` | the start and the recorded ending, `exit=124` |
| `artifacts/SNAPSHOT.txt` | the md5 of every input and state file as set up, and of the binary |
| `artifacts/trial_tr_ledger.md` | the 100-row trust-region ledger with the Krylov column |
| `artifacts/trial_linear_model.txt` | the four firings of the section 7 diagnostics, labeled by solve and outer iteration |
| `artifacts/trial_gated_by_pass.txt` | the gated rows of the pass certification |
| `artifacts/control_case1_pass_table_READ.md` | the control's table as READ from phase 6, with the `||R||` column taken from the control's own log |
| `analysis/` | the readers, including `tr_ledger.py` and `linear_model.py` written for this phase; the rest are phase 6's, copied unchanged |

The control's own raw artifacts stay where phase 6 left them,
`.../scratchpad/m1_alternation/c1_kzz1e9_HeH0.083/`.
