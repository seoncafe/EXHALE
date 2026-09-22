# LHS 1140 b, M1: the alternation audit of the four molecular models whose chemistry rows refuse

Phase 6 of `docs/PLAN_20260920_rev9.md` section 15, the audit of that
document's sections 9.2 and 9.3, carried out on 2026-09-21.

**Mode C throughout for the trajectories.** Every pass table below advances a
retained trajectory under the current solver from a named stored generation,
with a declared budget and a recorded ending, and is labeled Mode C in its own
table (rev9 section 3). Two kinds of number beside them are NOT Mode C and say
so where they stand: the identity probe of section 3, which is Mode R (one
evaluation of the current operator on an immutable checkpoint, no step taken),
and the values read out of stored manifests and certificates, which are a
historical audit.

Every number is labeled MEASURED (this document ran it), READ (from a named
file) or INSPECTED (source, with `file:line`). Nothing was published, no
golden was refreshed, and nothing under `LHS1140b/models/` or
`backup/regression/` was written: each case was copied into a scratch
directory and run there.

---

## 1. Verdict, first

**The alternation is approaching a joint fixed point in ONE of the four cases
and in none of the other three, and the three failures are three different
mechanisms.** All four spent their full six-pass budget; none was stopped.

| case | six-pass verdict | the gated rows over the six passes (MEASURED) | what the record of each pass says holds it |
|---|---|---|---|
| 1, `kzz1e9/HeH0.083` | **not approaching a fixed point: a drift under the movement bound** | carrier 8.150E-04 (cell 499) to 2.893E-03 (cell 500), a rise of 3.5; elemental 2.503E-04 to 7.063E-04 (cell 280), a rise of 2.8 | the carrier relaxation ended on the movement bound in every pass and the bound was SATURATED (its measure equals the bound to five digits at cells 418 to 426, r = 8.0 to 8.7); the hydrodynamic rows stayed within tolerance throughout, so nothing else obstructs |
| 2, `wellmixed/HeH0.083` | **not approaching a fixed point: the composition update takes the wind out of solution** | carrier 7.371E-02 to 7.668E-02 at cell 500 throughout, flat within 30 per cent, no trend | the entry state is a converged wind; ONE bounded composition update moves the mass flux by 9.85E-01 of its own column scale, and from pass 2 on the hydrodynamic solve returns `info = 2` with its mass and energy rows seven decades above tolerance and never recovers |
| 3, `kzz1e9/HeH0.55` | **approaching a fixed point, and the budget is what ended it** | joint residual 4.63E+05, 1.39E+03, 3.39E+02, 7.20E+02, 2.43E+02, 9.72E+01: a fall of 3.7 decades in five passes, with one excursion | nothing throttles it: the refusing energy row of the entry state (4.634E-01 at cell 257) is cured by the first composition update, and from pass 2 the carrier relaxation reaches its OWN fixed point at the fixed wind (residual 2E-11, 65 to 72 accepted trials, **zero** refusals, the bound never attained) |
| 4, `wellmixed/HeH0.55` | **not approaching a fixed point: the gated row rises monotonically** | carrier 2.363E-02, 3.300E-02, 5.029E-02, 7.919E-02, 9.505E-02, 1.066E-01, with its cell walking 306, 319, 327, 336, 345, 349 (r = 1.95 to 3.02) | the same two-factor failure as case 2 (the wind is lost at pass 2 and the mass row stands at 1.2E-01) together with a front that keeps moving outward, so the row and its cell both run away |

**A flat maximum was not taken as evidence of anything, and in the one case
where the maximum is flat (case 2) the record of each pass says what it is.** The
four discriminators rev9 section 9.2 asks for were each measured and each
decided something:

- **an inexact inner solve** is excluded in cases 1 and 3, where the element
  relaxation reached the fixed point of its own operator in every pass and, in
  case 3 from pass 2 on, the carrier relaxation did too;
- **a movement bound** is what holds case 1: the trial ledger shows the bound
  measure pinned at 1.0000E-02 and later 5.0000E-03 to five digits, with 25 to
  30 accepted trials and 25 to 29 refusals in every pass;
- **rejected updates** are counted in every pass and are reported in the tables
  as accepted/refused trials;
- **a two-cycle** is excluded in all four: the two-pass state difference
  `||Z(k+1) - Z(k-1)||` is never small against the one-pass differences it
  spans. Its ratio to their sum is 0.50 at worst and 1.00 at best over the
  sixteen measured steps (section 7.2), where aligned steps give 1 and a pair
  of states that undo each other gives 0;
- **a moving worst cell** is excluded as an illusion in cases 1, 2 and 3, where
  the refusing cell is fixed at 500 from pass 2 on, and is REAL in case 4,
  where it walks 43 cells outward while the row rises by a factor 4.5.

**On where the refusal sits (rev9 section 9.3).** In cases 1, 2 and 3 the
refusing carrier row follows the OUTER BOUNDARY, and the term that carries it
is named: at cell 500 the outer diffusive face is exactly zero,
`fdif_out = 0.0000000E+00` in every one of the eight files measured (four entry states and four ending states, MEASURED), which is
`carrier_outflow_ghost` continuing the last physical value into the ghost
(INSPECTED, `diffusive_photochemistry.f90:4167-4169` and `:4181-4188`), so the
inner diffusive face delivers H2 that neither the outer face nor the chemistry
of a 29 R_p cell can remove. In case 4 the refusal follows the WIND: its
largest gated measure is a localized interior peak at r = 3.015 whose loss is
dominated by the He+ and H2 channel, 4.373E-03 of a total loss of 4.425E-03.

---

## 2. The classification, the budget and the ending

### 2.1 Classification

| activity of this document | label |
|---|---|
| the four six-pass trajectories of section 6 | **Mode C**, a retained trajectory on the current branch |
| the four `stationary evaluate` probes of section 3 | **Mode R**, one evaluation of the current operator on an immutable checkpoint, no step taken |
| the certificates, manifests and `not_solved.md` notes quoted for context | historical audit, carrying the HISTORICAL executable identity |

No table mixes them.

### 2.2 The budget, declared before the runs

Predeclared in outer passes and residual evaluations, as rev9 section 15 asks,
and not only in wall clock:

| limit | value | why |
|---|---|---|
| outer passes per case | 6 | `EXHALE_OUTER_PASSES=6`; four cases, so **24 outer passes** is the whole budget |
| hydrodynamic solves per case | at most 6, one per pass | the route takes one solve per outer pass (`EXHALE_main.f90:7102-7126`) |
| residual evaluations | counted from the logs and reported in section 5; no cap was imposed on the inner solve, because `EXHALE_JFNK_MAXIT` caps the RUN and not the solve (INSPECTED, `EXHALE_main.f90:6983-6997`: after the first solve under a named cap the run refuses to enter the solver again), so setting it would have ended each trajectory after one pass |
| wall ceiling per case | 6000 s (100 min), `timeout 6000` | a stop at the ceiling is recorded as a stop and is never converted into a measured ending (rev9 section 12) |
| threads | `OMP_NUM_THREADS=8`, four cases concurrently on 72 cores | |

### 2.3 What was spent, and how each trajectory ended

MEASURED. Every run ended by spending its six passes and exiting 0; none
reached the 6000 s ceiling, so no ending in this document is a timeout.

| case | outer passes spent | hydrodynamic solves | JFNK iterations | Krylov products | carrier trials | sum of the pass wall clocks [s] | ending |
|---|---|---|---|---|---|---|---|
| 1, `kzz1e9/HeH0.083` | 6 of 6 | 6 | 150 | 133 | 264 | 1759.1 | the budget was spent; `Restart intent: stationary`, the stationary solve returned `info = 1`, the state was written `certified=F`, reason `no_stationary_claim`, exit 0 |
| 2, `wellmixed/HeH0.083` | 6 of 6 | 6 | 130 | 360 | 188 | 2111.0 | the same, and the run also printed `WARNING: the written state is NOT the state the gates accepted` |
| 3, `kzz1e9/HeH0.55` | 6 of 6 | 6 | 336 | 321 | 330 | 3293.4 | the same, no write-back warning |
| 4, `wellmixed/HeH0.55` | 6 of 6 | 6 | 194 | 580 | 192 | 2741.3 | the same, with the write-back warning |

Total spent: **24 outer passes, 24 hydrodynamic solves, 810 JFNK iterations,
1394 Krylov products, 974 carrier trials**, 9905 s of pass wall clock over
four runs made concurrently between 16:24 and 17:19 KST on `lart4`. The four
`stationary evaluate` probes of section 3.2 cost one residual evaluation each
and about 20 s each.

**One launch was discarded and is not evidence.** The four trajectories were
started once without `EXHALE_CARRIER_DEBUG=1` and stopped about ten minutes
later, before any of them had left its second pass, because without that key
the log carries no trial-by-trial record and rev9 section 9.2 asks for the
reason of every reduction and rejection. Their logs are kept under
`discarded_first_launch/` in the artifact directory so the discarded work is
visible; no number in this document comes from them.

---

## 3. Identity

### 3.1 The operator

| field | value |
|---|---|
| binary | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x`, md5 `89ec67149aa452eea3deeb19052f83d1` (MEASURED), the tree's binary of record |
| source text behind it | md5 of the sorted `md5sum` of every `src/**/*.f90` and `*.inc`: `fb3b5dd2c3410a4764d6d24f00190872` (MEASURED) |
| repository | `93eed86` with a dirty working tree (READ, `git log -1` and `git status`) |
| compiler identities in the image | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`, `GCC: (conda-forge gcc 16.2.0-4) 16.2.0`, `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` (MEASURED, `strings`) |
| linked libraries | `libopenblas.so.0`, `libgfortran.so.5`, `libgomp.so.1`, `libquadmath.so.0`, all from `/opt/miniconda3/lib`; `libm`, `libmvec`, `libpthread`, `libdl`, `librt`, `libc` from the system (MEASURED, `ldd`) |
| threads | `OMP_NUM_THREADS=8`; the run pins the BLAS pool to one thread itself |
| `EXHALE_RESID_QUAD` | UNSET in every run of this document, so the assembly selector of rev9 section 13.6 took its default. It is in no configuration record and is stated here because only the setter can state it |

The binaries that wrote the four stored states are older ones
(`EXHALE_75d55d9d.x` and `EXHALE_3146d11b.x`, READ from the manifests), so
every evaluation here is the CURRENT operator on an OLD state.

### 3.2 The four snapshots, and that the current operator reproduces their refusal

Each case was copied into its own scratch directory; the `latest_complete`
generation of its `state_index.json` was placed as the `_IC` pair. The state
md5 of each was MEASURED equal to the md5 its own manifest records.

| case | generation (READ, `state_index.json` `latest_complete`) | `Hydro_ioniz.txt` md5 (MEASURED = READ) | `Ion_species.txt` md5 (MEASURED = READ) |
|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `g0004_20260920T053734Z_79b42a03` | `fd882d8918052069f8ed229dd376f8c4` | `87288b03ad88a92b7b8b8b6c2b1bdd66` |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `g0004_20260920T053815Z_1ab3de51` | `74260bbcc103a26006f90616f7ffe835` | `840e46512e0ac6bba610e3e51d6d5336` |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | `g0005_20260920T053754Z_2c1363df` | `fd9619de632e274073ff3e68aee3e966` | `fe69ce07009fc964229abbde7d7aaf20` |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `g0004_20260920T053837Z_46ab9e80` | `2e342e79c8ec5c6f52f4da3a7bb7daee` | `90a23d0598a080765a77a4dfef592bce` |

No generation of any of the four is certified: `latest_certified` is `none` in
all four indexes (READ).

**The Mode R probe.** Each snapshot was evaluated once with
`Restart intent: stationary evaluate`, which measures the state as loaded and
takes no step. Every gated species row and every refusing energy row comes
back at the value the stored certificate carries, to all four printed digits:

| case | row | stored certificate (READ) | Mode R under `89ec6714` (MEASURED) |
|---|---|---|---|
| 1 | carrier balance H2, gated | 8.150E-04 at cell 499 | 8.150E-04 at cell 499 |
| 1 | elemental transport He/H, gated | 2.503E-04 at cell 280 | 2.503E-04 at cell 280 |
| 2 | carrier balance H2, gated | 7.371E-02 at cell 500 | 7.371E-02 at cell 500 |
| 3 | hydrodynamic energy row | 4.634E-01 at cell 257 | 4.634E-01 at cell 257 |
| 3 | carrier balance H2, gated | 8.922E-03 at cell 261 | 8.922E-03 at cell 261 |
| 4 | hydrodynamic energy row | 1.000E+00 at cell 246 | 1.000E+00 at cell 246 |
| 4 | carrier balance H2, gated | 2.363E-02 at cell 306 | 2.363E-02 at cell 306 |

What does NOT come back identically is the hydrodynamic mass and energy rows
of the base cells and the He 2 3S level row, which move by tens of per cent of
themselves at values four to eleven decades below their tolerances: case 1
mass 1.393E-08 both times but energy 6.461E-09 at cell 31 against 7.651E-09 at
cell 3; case 2 mass 2.339E-08 at cell 1 against 2.164E-08; case 3 mass
3.805E-09 against 3.765E-09; case 4 mass 1.672E-08 both times. That is the
re-entry behaviour `REPRODUCE.md` of these cases already states: the first
equilibrium sweep of a re-entry moves the composition by 1e-14 to 1e-11 and
the base flux assembly stands at its rounding floor. None of it moves a
verdict, and every refusing row is reproduced exactly.

### 3.3 The effective configuration

READ from `EXHALE_resolved.out` of each scratch run (the derived route keys,
`write_setup_report.f90:796` onward) and from `EXHALE_setup.out` (the prose
keys, which the resolved record does not carry):

| key | cases 1 and 3 (`kzz1e9`) | cases 2 and 4 (`wellmixed`) |
|---|---|---|
| `he_diffusion` | T | **F** |
| `carrier_transport` | T | T |
| `carrier_in_newton` | F | F |
| `carrier_newton_on_stall` | F | F |
| `ionization_transport` | F | F |
| `oxygen_chemistry` | F | F |
| `well_balanced` | T | T |
| numerical flux | HLLC | HLLC |
| reconstruction method | PLM | PLM |
| boundary model of the loaded state | `characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3` | the same |

So the alternation of cases 2 and 4 has ONE composition half, the carrier
relaxation; only cases 1 and 3 run the element relaxation as well. Their
certificates say the same thing from the other side: 7 active equations in
cases 1 and 3, 6 in cases 2 and 4 (READ).

### 3.4 The one edit made to a case input, and why

The four case directories are used unchanged except for two lines, both
recorded here:

1. `Spectrum file:` was made absolute. The case states it relative to the case
   directory and the run's working directory is the scratch copy.
2. `Restart intent:` was set to `stationary` in all four. Cases 2, 3 and 4
   already state exactly that; case 1 alone still states
   `stationary equilibrate`, which puts the loaded composition on its own
   fixed point BEFORE the evaluation and so changes the state the trajectory
   starts from (INSPECTED, `input_read.f90:76-83`). An audit of the
   alternation has to start at the snapshot it names, so the equilibration was
   not taken, and case 1 is then configured exactly as the other three.

`Coupled carrier solve: False` is stated in all four (cases 2, 3 and 4 state
it themselves; it was appended to case 1, which leaves the key out and gets
the same value by default). It is the route rev9 section 9.1 describes.

---

## 4. What the route is, and where each recorded quantity comes from

INSPECTED, `EXHALE_main.f90`, the outer pass loop from `:7100`. One pass is:

1. the hydrodynamic solve, `solve_steady_jfnk` at `:7122`;
2. the refresh: `U_to_W`, `ioniz_eq`, and `T = p/(n_tot + n_e)` of the
   composition beside it (`:7128-7146`);
3. **the joint test**: `assemble_residual` then `certification_evaluate` at
   `:7172-7178`, on that refreshed state. This is the acceptance criterion,
   and it is what the pass line reports as its worst gated species row and as
   the coupled residual of the whole state;
4. the progress control (`:7345` onward), which shortens the movement bound or
   the under-relaxation factor when neither the coupled nor the composition
   residual fell;
5. `relax_element_composition` at the face mass flux of this state
   (`:7475-7481`), under `he_diffusion`;
6. `relax_photochemical_composition` (`:7540-7544`), under
   `transported_rows_exist() .and. .not. block_now`;
7. the log line and the diagnostics.

So the composition the pass hands over is consumed by the NEXT pass's solve,
and the certification of pass k is of the state the solve of pass k produced
from the composition of pass k-1. That is the residual of the coupled map, as
rev9 section 9.1 says, and not of the last relaxation alone.

**Two state sequences, kept apart.** `Y(k)` is the state at the joint test of
pass k, which is the one the certificate of pass k describes. `Z(k)` is the
state at the END of pass k, after both relaxations, which is what
`write_diffusion_pass_profile` exports (`EXHALE_main.f90:8258`, under
`EXHALE_DIFFUSION_CHECK=1`) and what the next solve is entered with. The
compact norms of section 7 are norms of the `Z` sequence, because that is the
sequence the export carries.

**On the "map defect".** No second map was evaluated. For the implemented
deterministic update the accepted displacement IS `G(Y) - Y`, and rev9 section
9.2 is followed: what is reported is the change REQUESTED, the change
ACCEPTED, the reason for every reduction and rejection, the two state
differences, and the refreshed joint residual.

**Where each recorded quantity is written.** All INSPECTED in
`EXHALE_main.f90` unless stated:

| quantity | the line that writes it |
|---|---|
| joint residual, worst gated row and cell | the pass line, `:7805-7812` |
| element relaxation ending, inner map distance, displacement applied, steps | `:7813-7819` |
| carrier relaxation ending, displacement kept, transport steps | `:7821-7828` |
| carrier residual of the returned composition, with `|res|` and its own row terms | `:7844-7848` |
| element residual of the returned composition, with cell and element | `:7849-7864` |
| accepted composition displacement, element and carrier separately | `:7866-7870` |
| coupled residual of the whole state, and which measure decided progress | `:7871-7875` |
| the cell the movement bound was refused on, and the measure it was refused at | `:7896-7915` |
| every carrier trial, its growth multiplier, its verdict and its measure | `diffusive_photochemistry.f90:8020-8034` (refused) and `:8060` onward (accepted), under `EXHALE_CARRIER_DEBUG=1` |
| the H2 front positions `x2 = 0.5` and `x2 = 1e-2` | `:7989-7992` |

**A window difference that has to be said once.** The carrier and element
residuals each pass reports are COLUMN MAXIMA over all 500 cells. The certification
gates the same rows only in the wind, at `r >= 1.20`, and reports the cells
below that radius without gating. So on case 1 the pass line's carrier
residual 1.70E-02 is the base cell 1 and the gated verdict 8.15E-04 is cell
499: the two numbers are the same row measured over two different windows, and
they are never compared with each other below.

---

## 5. The gate of section 4.2, for every diagnostic file this audit used

Every diagnostic invocation ran in its own output directory, and no existing
product was touched.

| requirement (rev9 section 4.2) | how it is met |
|---|---|
| the complete environment actually in force | Mode R probes: `OMP_NUM_THREADS=8`, `EXHALE_CARRIER_ROW_TERMS=1`; nothing else set. Mode C runs: `OMP_NUM_THREADS=8`, `EXHALE_PTC_DTAU0=1.0`, `EXHALE_OUTER_PASSES=6`, `EXHALE_DIFFUSION_CHECK=1`, `EXHALE_CARRIER_ROW_TERMS=1`, `EXHALE_CARRIER_DEBUG=1`; `EXHALE_RESID_QUAD`, `EXHALE_JFNK_MAXIT`, `EXHALE_CARRIER_TRUST`, `EXHALE_CARRIER_TRUST_HOLD` and `EXHALE_DIFF_OMEGA` all UNSET |
| the file, its existence, nonzero size, header and row count | `output/carrier_row_terms.txt`: present in all eight invocations, 288177 bytes each, 9 header lines, **500 cell rows, 500 face rows and 500 channel rows**, one carrier (`H2`), and **15 values on every channel row**, the count taken from the file's own `# chan cell ...` header line and not from a constant (MEASURED). `output/diffusion_pass_profiles.txt`: present in all four Mode C runs, 402187 bytes each, 2 header lines and 3000 data rows, which is 500 rows in each of the six blocks the file labels 1 to 6 in its first column (MEASURED) |
| which state the rows belong to | the writer emits no identifier (INSPECTED, `diffusive_photochemistry.f90:2784-2790`; the path is opened `status='replace'` and the certification does not remove it), so the invocation is what identifies it: in the Mode R probes the file can only be the loaded state, because no step is taken; in the Mode C runs the file is the LAST certification of the run, which is the `final state, as written` certification of `EXHALE_main.f90:4380`, so it belongs to the state the run's `output/Hydro_ioniz.txt` carries. Intermediate passes overwrite the same path and their content is not used |
| the exit class | all eight exited 0; the four Mode C runs ended on their own pass budget and none reached the wall ceiling |

The four Mode R files and the four Mode C files are archived separately with
this audit's artifacts, so the entry state and the ending state of each case
are two files and not one path measured twice.

---

## 6. The four trajectories, pass by pass

Each table is **Mode C**: a retained trajectory advanced under
`EXHALE.x` md5 `89ec67149aa452eea3deeb19052f83d1` from the generation named in
section 3.2, budget six outer passes, ending as recorded in section 2.3. Every
value is MEASURED, read from the line of `run.log` named in section 4.

Conventions used in every table:

- **joint residual** is `prog_worst`, the largest certification entry over its
  own tolerance, evaluated on the refreshed state of that pass. It is the
  acceptance criterion.
- **worst gated row, cell** is the same certification's worst GATED species
  row, which is gated at `r >= 1.20`.
- **element residual returned** and **carrier residual returned** are COLUMN
  MAXIMA, over all 500 cells, of each relaxation's own operator at the
  composition the pass hands back. They are not the gated numbers.
- **displacement ADOPTED** is what the pass actually moved the composition by;
  **bound REQUESTED** is the largest bound measure any trial of that pass asked
  for, with the cell that refused it. A requested value above the bound was
  rejected and the growth multiplier halved.
- **trials acc/ref** counts the carrier relaxation's accepted and refused
  trials of that pass.
- pass 6 of every case takes no composition update: it is the last pass and no
  solve is left to consume one (`EXHALE_main.f90`, the update is guarded by
  `it_diff .lt. pass_cap`).

### 6.1 Case 1, `molecular_scalar_gj1132_kzz1e9/HeH0.083`

| pass | hydro info | joint residual (worst gated entry over its tolerance) | worst gated row, cell | element relaxation ending | element residual returned | carrier relaxation ending | carrier residual returned | displacement ADOPTED (element, carrier) | bound: REQUESTED / cell | trials acc/ref | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.150e+01 | 8.150e-04 carrier balance H2, cell 499 | the fixed point of the element operator | 5.670e-02 | the movement bound, with the last step inside it | 1.700e-02 | 2.740e-02 (1.340e-03, 2.740e-02) | 1.1990e-02 at cell 421 | 25/26 | 8.6 |
| 2 | 0 | 6.290e+02 | 6.290e-03 carrier balance H2, cell 500 | the fixed point of the element operator | 7.430e-02 | the movement bound, with the last step inside it | 1.680e-02 | 2.800e-02 (2.270e-03, 2.800e-02) | 1.1420e-02 at cell 424 | 26/27 | 413.7 |
| 3 | 0 | 6.260e+02 | 6.260e-03 carrier balance H2, cell 500 | the fixed point of the element operator | 8.380e-02 | the movement bound, with the last step inside it | 1.660e-02 | 2.820e-02 (2.940e-03, 2.820e-02) | 1.1070e-02 at cell 423 | 30/29 | 346.8 |
| 4 | 0 | 6.390e+02 | 6.390e-03 carrier balance H2, cell 500 | the fixed point of the element operator | 9.600e-02 | the movement bound, with the last step inside it | 1.670e-02 | 1.430e-02 (1.620e-03, 1.430e-02) | 5.6890e-03 at cell 425 | 26/27 | 344.3 |
| 5 | 0 | 2.950e+02 | 2.950e-03 carrier balance H2, cell 500 | the fixed point of the element operator | 9.650e-02 | the movement bound, with the last step inside it | 1.660e-02 | 1.410e-02 (1.610e-03, 1.410e-02) | 6.9490e-03 at cell 418 | 23/25 | 354.1 |
| 6 | 0 | 2.890e+02 | 2.890e-03 carrier balance H2, cell 500 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | 0.000e+00 (0.000e+00, 0.000e+00) | -- | 0/0 | 291.5 |

front and worst-cell record:
| pass | H2 front x2=0.5 | x2=1e-2 | carrier steady residual, worst cell | displacement worst cell | bound held at cell | trust |
|---|---|---|---|---|---|---|
| 1 | 1.3837 | > r(N) | 1.700e-02 at cell 1 (r 1.000) | 427 (r 8.853) | 422 (r 8.197) | 1.000e-02 |
| 2 | 1.3837 | > r(N) | 1.680e-02 at cell 1 (r 1.000) | 430 (r 9.274) | 426 (r 8.717) | 1.000e-02 |
| 3 | 1.3905 | > r(N) | 1.660e-02 at cell 1 (r 1.000) | 429 (r 9.131) | 424 (r 8.453) | 1.000e-02 |
| 4 | 1.3905 | > r(N) | 1.670e-02 at cell 1 (r 1.000) | 430 (r 9.274) | 426 (r 8.717) | 5.000e-03 |
| 5 | 1.3905 | > r(N) | 1.660e-02 at cell 1 (r 1.000) | 425 (r 8.584) | 420 (r 7.950) | 5.000e-03 |
| 6 | 1.3905 | > r(N) | 1.620e-02 at cell 1 (r 1.000) | 425 (r 8.584) | -- (r --) | 5.000e-03 |

The gated rows of each pass, and the column maxima beside them:

```
pass 1: hydrodynamic mass row 1.8360e-08 at cell 1 (within) | hydrodynamic energy row 7.6510e-09 at cell 3 (within) | carrier balance H2 gated 8.1500e-04 at cell 499 [column max 1.7440e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 2.5030e-04 at cell 280 [column max 4.6560e-02 at cell 2, ABOVE]
pass 2: hydrodynamic mass row 3.6810e-09 at cell 1 (within) | hydrodynamic energy row 5.8170e-09 at cell 3 (within) | carrier balance H2 gated 6.2900e-03 at cell 500 [column max 1.6030e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 4.5320e-04 at cell 280 [column max 7.9960e-02 at cell 2, ABOVE]
pass 3: hydrodynamic mass row 4.7500e-09 at cell 1 (within) | hydrodynamic energy row 5.5710e-09 at cell 3 (within) | carrier balance H2 gated 6.2580e-03 at cell 500 [column max 1.5640e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 6.0690e-04 at cell 280 [column max 9.6870e-02 at cell 2, ABOVE]
pass 4: hydrodynamic mass row 5.2130e-09 at cell 1 (within) | hydrodynamic energy row 7.6870e-09 at cell 3 (within) | carrier balance H2 gated 6.3880e-03 at cell 500 [column max 1.5490e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 6.9610e-04 at cell 280 [column max 1.0670e-01 at cell 2, ABOVE]
pass 5: hydrodynamic mass row 7.5730e-09 at cell 1 (within) | hydrodynamic energy row 6.3630e-09 at cell 3 (within) | carrier balance H2 gated 2.9500e-03 at cell 500 [column max 1.6240e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 6.9680e-04 at cell 280 [column max 1.0580e-01 at cell 2, ABOVE]
pass 6: hydrodynamic mass row 1.3540e-08 at cell 1 (within) | hydrodynamic energy row 6.1530e-09 at cell 59 (within) | carrier balance H2 gated 2.8930e-03 at cell 500 [column max 1.6160e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 7.0630e-04 at cell 280 [column max 1.0750e-01 at cell 2, ABOVE]
```

The carrier relaxation trial ledger of each pass:

```
pass 1: 51 trials, 25 accepted, 26 refused {'the movement bound': 26}
    largest REQUESTED bound measure 1.1990e-02 at cell 421 (trial 21, grow 1.108e+03), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 422 (trial 40, grow 5.137e-01)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 2: 53 trials, 26 accepted, 27 refused {'the movement bound': 27}
    largest REQUESTED bound measure 1.1420e-02 at cell 424 (trial 21, grow 1.108e+03), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 426 (trial 51, grow 9.766e-04)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 3: 59 trials, 30 accepted, 29 refused {'the movement bound': 29}
    largest REQUESTED bound measure 1.1070e-02 at cell 423 (trial 21, grow 1.108e+03), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 424 (trial 55, grow 3.810e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 4: 53 trials, 26 accepted, 27 refused {'the movement bound': 27}
    largest REQUESTED bound measure 5.6890e-03 at cell 425 (trial 19, grow 4.926e+02), refused
    last ACCEPTED bound measure   5.0000e-03 at cell 426 (trial 51, grow 9.766e-04)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 5: 48 trials, 23 accepted, 25 refused {'the movement bound': 25}
    largest REQUESTED bound measure 6.9490e-03 at cell 418 (trial 18, grow 9.853e+02), refused
    last ACCEPTED bound measure   5.0000e-03 at cell 420 (trial 44, grow 3.568e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 6: no carrier trial recorded
```

The compact state norms (section 7 states the blocks and the scales):

```
# fixed block scales S_b, from pass 1 (the entry state):
#   HeH_ratio  1.000334e+00   (He/H)/HeH, the elemental partition
#   mass_flux  7.271635e-06   r^2 rho v, the wind
#   T          2.383802e+03   temperature [K]
#   x_H2       7.888761e-01   2 n(H2)/n_H, the carrier fraction
#   n_H3p      5.199933e+04   n(H3+) [cm^-3]
#   n_HeII     2.556521e+02   n(He II) [cm^-3]
#   n_e        1.583673e+07   n_e [cm^-3]
# k   ||Y(k+1)-Y(k)||  ||Y(k+1)-Y(k-1)||    HeH_ratio   mass_flux           T        x_H2       n_H3p      n_HeII         n_e
  1->2   2.1078e-02         n/a           2.189e-03   2.124e-03   2.752e-03   2.108e-02   1.178e-02   5.416e-03   3.115e-03
  2->3   2.1262e-02      4.2329e-02       2.820e-03   2.076e-03   2.617e-03   2.126e-02   1.183e-02   7.085e-03   3.001e-03
  3->4   1.0824e-02      3.2083e-02       1.557e-03   2.146e-03   2.634e-03   1.082e-02   8.416e-03   6.063e-03   4.798e-03
  4->5   1.0713e-02      2.1480e-02       1.543e-03   9.208e-04   1.262e-03   1.071e-02   5.493e-03   4.211e-03   1.341e-03
  5->6   4.0254e-03      1.0713e-02       0.000e+00   1.046e-03   1.245e-03   0.000e+00   2.754e-03   4.025e-03   3.082e-03
```

**Verdict: the alternation is NOT approaching a joint fixed point, and what
holds it is the movement bound on the carrier half.**

The hydrodynamic half is never the obstruction here: `info = 0` in all six
passes and the mass, momentum and energy rows are inside their tolerances at
every one of them. Both gated species rows nevertheless move AWAY from the
entry state over the six passes, the carrier from 8.150E-04 to 2.893E-03 and
the element from 2.503E-04 to 7.063E-04.

The trial ledger says what the carrier relaxation is doing. In each of the
five update passes the bound measure is driven to the bound and pinned there:
the last accepted trial stands at 1.0000E-02 to five digits in passes 1 to 3
and at 5.0000E-03 in passes 4 and 5, against a bound of exactly those values,
and 25 to 29 trials per pass are refused on it while the growth multiplier is
halved from as much as 1.1E+03 down to the floor 9.77E-04. The cell that holds
the bound is 418 to 426, r = 7.95 to 8.72, which is neither the refusing cell
(499, then 500) nor the H2 front (r = 1.38, which moves by 0.007 R_p in six
passes). One cell of the far wind is rationing the movement of the whole
column.

The element relaxation, by contrast, reaches the fixed point of its own
operator in every pass, so no inner inexactness is available as an
explanation, and its returned residual rises monotonically 5.67E-02, 7.43E-02,
8.38E-02, 9.60E-02, 9.65E-02 while it does so: the element half is converging
to a fixed point that the hydrodynamic refresh then moves.

The state norms rule out a cycle. `||Z(k+1) - Z(k)||` is 2.11E-02, 2.13E-02,
1.08E-02, 1.07E-02, 4.03E-03, falling in step with the halving of the bound at
pass 4 and not with any residual; `||Z(k+1) - Z(k-1)||` is 4.23E-02, 3.21E-02,
2.15E-02, 1.07E-02, which is the SUM of the two neighbouring one-pass norms to
three digits at the first three and 0.73 of it at the last. Consecutive passes
therefore move the state in the same direction, by an amount the bound sets,
and the joint residual does not fall. The dominant block is `x_H2` itself at
the first four steps.

### 6.2 Case 2, `molecular_scalar_gj1132_wellmixed/HeH0.083`

| pass | hydro info | joint residual (worst gated entry over its tolerance) | worst gated row, cell | element relaxation ending | element residual returned | carrier relaxation ending | carrier residual returned | displacement ADOPTED (element, carrier) | bound: REQUESTED / cell | trials acc/ref | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.370e+03 | 7.370e-02 carrier balance H2, cell 500 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 8.000e-02 | 1.730e-02 (0.000e+00, 1.730e-02) | 1.4320e-02 at cell 272 | 17/21 | 5.4 |
| 2 | 2 | 4.510e+07 | 6.000e-02 carrier balance H2, cell 500 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 8.120e-02 | 1.730e-02 (0.000e+00, 1.730e-02) | 1.3620e-02 at cell 275 | 16/21 | 472.6 |
| 3 | 2 | 4.790e+07 | 7.630e-02 carrier balance H2, cell 500 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 5.890e-02 | 8.570e-03 (0.000e+00, 8.570e-03) | 6.3170e-03 at cell 276 | 16/21 | 462.9 |
| 4 | 2 | 6.790e+07 | 6.800e-02 carrier balance H2, cell 500 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 5.890e-02 | 8.550e-03 (0.000e+00, 8.550e-03) | 6.1180e-03 at cell 278 | 18/22 | 411.3 |
| 5 | 2 | 3.840e+07 | 7.730e-02 carrier balance H2, cell 500 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 8.900e-02 | 8.530e-03 (0.000e+00, 8.530e-03) | 5.9840e-03 at cell 280 | 15/21 | 419.8 |
| 6 | 2 | 3.130e+07 | 7.670e-02 carrier balance H2, cell 500 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | 0.000e+00 (0.000e+00, 0.000e+00) | -- | 0/0 | 339.0 |

front and worst-cell record:
| pass | H2 front x2=0.5 | x2=1e-2 | carrier steady residual, worst cell | displacement worst cell | bound held at cell | trust |
|---|---|---|---|---|---|---|
| 1 | 1.4044 | 18.8174 | 8.000e-02 at cell 1 (r 1.000) | 273 (r 1.535) | 271 (r 1.516) | 1.000e-02 |
| 2 | 1.4261 | 17.6173 | 8.120e-02 at cell 1 (r 1.000) | 276 (r 1.564) | 274 (r 1.544) | 1.000e-02 |
| 3 | 1.4337 | 17.0480 | 5.890e-02 at cell 500 (r 29.031) | 278 (r 1.584) | 276 (r 1.564) | 5.000e-03 |
| 4 | 1.4413 | 16.4981 | 5.890e-02 at cell 500 (r 29.031) | 280 (r 1.604) | 278 (r 1.584) | 5.000e-03 |
| 5 | 1.4491 | 15.9671 | 8.900e-02 at cell 1 (r 1.000) | 281 (r 1.615) | 280 (r 1.604) | 5.000e-03 |
| 6 | 1.4491 | 15.9671 | 3.640e-01 at cell 1 (r 1.000) | 281 (r 1.615) | -- (r --) | 5.000e-03 |

The gated rows of each pass, and the column maxima beside them:

```
pass 1: hydrodynamic mass row 3.2070e-08 at cell 1 (within) | hydrodynamic energy row 2.8430e-08 at cell 2 (within) | carrier balance H2 gated 7.3710e-02 at cell 500 [column max 9.2640e-02 at cell 1, ABOVE]
pass 2: hydrodynamic mass row 3.4650e-01 at cell 1 (ABOVE) | hydrodynamic energy row 1.1310e-01 at cell 1 (ABOVE) | carrier balance H2 gated 6.0010e-02 at cell 500 [column max 2.8310e-01 at cell 1, ABOVE]
pass 3: hydrodynamic mass row 2.5160e-01 at cell 2 (ABOVE) | hydrodynamic energy row 6.9920e-01 at cell 1 (ABOVE) | carrier balance H2 gated 7.6350e-02 at cell 500 [column max 2.1690e-01 at cell 2, ABOVE]
pass 4: hydrodynamic mass row 9.3420e-01 at cell 2 (ABOVE) | hydrodynamic energy row 1.1350e+00 at cell 2 (ABOVE) | carrier balance H2 gated 6.7960e-02 at cell 500 [column max 7.2680e-01 at cell 2, ABOVE]
pass 5: hydrodynamic mass row 6.1860e-01 at cell 1 (ABOVE) | hydrodynamic energy row 8.9630e-01 at cell 1 (ABOVE) | carrier balance H2 gated 7.7260e-02 at cell 500 [column max 4.9750e-01 at cell 2, ABOVE]
pass 6: hydrodynamic mass row 3.9000e-01 at cell 1 (ABOVE) | hydrodynamic energy row 1.0540e+00 at cell 1 (ABOVE) | carrier balance H2 gated 7.6680e-02 at cell 500 [column max 3.6450e-01 at cell 1, ABOVE]
```

The carrier relaxation trial ledger of each pass:

```
pass 1: 38 trials, 17 accepted, 21 refused {'the movement bound': 21}
    largest REQUESTED bound measure 1.4320e-02 at cell 272 (trial 12, grow 8.650e+01), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 271 (trial 35, grow 2.506e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 2: 37 trials, 16 accepted, 21 refused {'the movement bound': 21}
    largest REQUESTED bound measure 1.3620e-02 at cell 275 (trial 12, grow 8.650e+01), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 274 (trial 33, grow 3.341e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 3: 37 trials, 16 accepted, 21 refused {'the movement bound': 21}
    largest REQUESTED bound measure 6.3170e-03 at cell 276 (trial 12, grow 2.883e+01), refused
    last ACCEPTED bound measure   5.0000e-03 at cell 276 (trial 32, grow 6.682e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 4: 40 trials, 18 accepted, 22 refused {'the movement bound': 22}
    largest REQUESTED bound measure 6.1180e-03 at cell 278 (trial 12, grow 2.883e+01), refused
    last ACCEPTED bound measure   5.0000e-03 at cell 278 (trial 38, grow 9.766e-04)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 5: 36 trials, 15 accepted, 21 refused {'the movement bound': 21}
    largest REQUESTED bound measure 5.9840e-03 at cell 280 (trial 12, grow 2.883e+01), refused
    last ACCEPTED bound measure   5.0000e-03 at cell 280 (trial 34, grow 9.766e-04)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 6: no carrier trial recorded
```

The compact state norms (section 7 states the blocks and the scales):

```
# fixed block scales S_b, from pass 1 (the entry state):
#   HeH_ratio  1.000000e+00   (He/H)/HeH, the elemental partition
#   mass_flux  9.543801e-06   r^2 rho v, the wind
#   T          2.835758e+03   temperature [K]
#   x_H2       9.723396e-01   2 n(H2)/n_H, the carrier fraction
#   n_H3p      1.677188e+05   n(H3+) [cm^-3]
#   n_HeII     2.265708e+04   n(He II) [cm^-3]
#   n_e        1.367694e+07   n_e [cm^-3]
# k   ||Y(k+1)-Y(k)||  ||Y(k+1)-Y(k-1)||    HeH_ratio   mass_flux           T        x_H2       n_H3p      n_HeII         n_e
  1->2   9.8483e-01         n/a           0.000e+00   9.848e-01   1.254e-02   1.726e-02   3.255e-02   4.592e-02   5.692e-03
  2->3   3.5458e-02      9.8695e-01       0.000e+00   3.577e-03   1.450e-02   8.567e-03   3.546e-02   3.276e-02   3.500e-03
  3->4   2.1421e-02      5.6798e-02       0.000e+00   1.959e-03   9.160e-03   8.546e-03   2.134e-02   2.142e-02   2.957e-03
  4->5   1.8096e-02      3.7500e-02       0.000e+00   3.183e-03   5.459e-03   8.524e-03   1.315e-02   1.810e-02   2.030e-03
  5->6   1.7467e-02      3.2307e-02       0.000e+00   3.139e-03   7.461e-03   0.000e+00   1.747e-02   1.424e-02   3.602e-03
```

**Verdict: the alternation is NOT approaching a joint fixed point, and the
obstruction is the splitting: one bounded composition update takes the wind out
of solution and the hydrodynamic solve never gets it back.**

The entry state is a converged wind: at pass 1 the mass row is 3.207E-08 and
the energy row 2.843E-08, both within, and only the carrier row refuses. The
carrier relaxation of pass 1 ends on the movement bound with 17 accepted and 21
refused trials and an accepted displacement of 1.73E-02. What the next solve is
handed then differs from the state the hydrodynamic rows were last solved at by
`||Z(2) - Z(1)|| = 9.85E-01`, carried ENTIRELY by the mass flux block: the wind
moved by 98.5 per cent of its own column maximum. From pass 2 the solve returns
`info = 2` at every pass and the hydrodynamic rows stand at 3.5E-01 to 1.1E+00
against tolerances of 2.8E-08 and 1.0E-06, and they are still there at pass 6.

The gated carrier row is flat over the six passes, 7.371E-02, 6.001E-02,
7.635E-02, 6.796E-02, 7.726E-02, 7.668E-02, always at cell 500. **That flat
maximum is not a fixed point of anything.** It is a row evaluated on a state
whose hydrodynamics the alternation has destroyed, and the two-pass norms
(9.87E-01, 5.68E-02, 3.75E-02, 3.23E-02) show the state still moving by one to
four per cent of its own scale at the last pass. This case is the clearest
instance in the four of what rev9 section 9.4 states about the splitting: there
is no step between a composition good enough for the carrier row and one inside
the range in which holding the wind is admissible.

This run is also one of the two that printed
`WARNING: the written state is NOT the state the gates accepted`.

### 6.3 Case 3, `molecular_scalar_gj1132_kzz1e9/HeH0.55`

| pass | hydro info | joint residual (worst gated entry over its tolerance) | worst gated row, cell | element relaxation ending | element residual returned | carrier relaxation ending | carrier residual returned | displacement ADOPTED (element, carrier) | bound: REQUESTED / cell | trials acc/ref | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 4.630e+05 | 8.920e-03 carrier balance H2, cell 261 | the fixed point of the element operator | 2.610e-03 | the movement bound, with the last step inside it | 4.760e-03 | 4.820e-02 (3.690e-05, 4.820e-02) | 1.0380e-02 at cell 266 | 24/26 | 520.8 |
| 2 | 0 | 1.390e+03 | 1.390e-02 carrier balance H2, cell 500 | the fixed point of the element operator | 4.400e-02 | the carriers stopped moving at the fixed wind | 2.080e-11 | 1.060e-01 (1.060e-01, 4.780e-02) | 7.4980e-03 at cell 277 | 72/0 | 711.0 |
| 3 | 0 | 3.390e+02 | 3.390e-03 carrier balance H2, cell 500 | the fixed point of the element operator | 1.070e-02 | the carriers stopped moving at the fixed wind | 1.650e-11 | 3.890e-02 (2.500e-02, 3.890e-02) | 5.3300e-03 at cell 259 | 65/0 | 774.4 |
| 4 | 0 | 7.200e+02 | 7.200e-03 carrier balance H2, cell 500 | the fixed point of the element operator | 8.980e-03 | the carriers stopped moving at the fixed wind | 1.850e-11 | 2.440e-02 (8.420e-03, 2.440e-02) | 2.8030e-03 at cell 266 | 72/0 | 436.1 |
| 5 | 0 | 2.430e+02 | 2.430e-03 elemental transport He/H partition, cell 335 | the fixed point of the element operator | 1.140e-02 | the carriers stopped moving at the fixed wind | 2.140e-11 | 1.380e-02 (1.380e-02, 4.820e-03) | 5.1200e-04 at cell 394 | 71/0 | 467.6 |
| 6 | 0 | 9.720e+01 | 9.720e-04 carrier balance H2, cell 500 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | 0.000e+00 (0.000e+00, 0.000e+00) | -- | 0/0 | 383.5 |

front and worst-cell record:
| pass | H2 front x2=0.5 | x2=1e-2 | carrier steady residual, worst cell | displacement worst cell | bound held at cell | trust |
|---|---|---|---|---|---|---|
| 1 | 1.1271 | 2.1332 | 4.760e-03 at cell 302 (r 1.888) | 262 (r 1.441) | 265 (r 1.465) | 1.000e-02 |
| 2 | 1.1271 | 13.5728 | 2.080e-11 at cell 3 (r 1.001) | 340 (r 2.723) | -- (r --) | 1.000e-02 |
| 3 | 1.1271 | 9.7185 | 1.650e-11 at cell 25 (r 1.005) | 385 (r 4.776) | -- (r --) | 1.000e-02 |
| 4 | 1.1271 | 2.0943 | 1.850e-11 at cell 298 (r 1.828) | 384 (r 4.710) | -- (r --) | 1.000e-02 |
| 5 | 1.1271 | 2.2151 | 2.140e-11 at cell 7 (r 1.001) | 393 (r 5.341) | -- (r --) | 1.000e-02 |
| 6 | 1.1271 | 2.2151 | 9.720e-04 at cell 500 (r 29.031) | 393 (r 5.341) | -- (r --) | 1.000e-02 |

The gated rows of each pass, and the column maxima beside them:

```
pass 1: hydrodynamic mass row 3.8050e-09 at cell 2 (within) | hydrodynamic energy row 4.6340e-01 at cell 257 (ABOVE) | carrier balance H2 gated 8.9220e-03 at cell 261 [column max 8.9220e-03 at cell 261, ABOVE] | elemental transport He/H partition gated 2.1390e-06 at cell 339 [column max 4.4600e-05 at cell 2, within]
pass 2: hydrodynamic mass row 4.0580e-09 at cell 1 (within) | hydrodynamic energy row 8.0460e-09 at cell 2 (within) | carrier balance H2 gated 1.3860e-02 at cell 500 [column max 2.4190e-02 at cell 1, ABOVE] | elemental transport He/H partition gated 2.6600e-03 at cell 307 [column max 1.1750e-01 at cell 2, ABOVE]
pass 3: hydrodynamic mass row 5.0120e-09 at cell 2 (within) | hydrodynamic energy row 1.1410e-08 at cell 2 (within) | carrier balance H2 gated 3.3900e-03 at cell 500 [column max 3.3900e-03 at cell 500, ABOVE] | elemental transport He/H partition gated 3.1190e-03 at cell 333 [column max 2.3540e-02 at cell 2, ABOVE]
pass 4: hydrodynamic mass row 7.4050e-09 at cell 1 (within) | hydrodynamic energy row 1.2870e-08 at cell 2 (within) | carrier balance H2 gated 7.2030e-03 at cell 500 [column max 7.2030e-03 at cell 500, ABOVE] | elemental transport He/H partition gated 1.0900e-03 at cell 314 [column max 7.8830e-03 at cell 2, ABOVE]
pass 5: hydrodynamic mass row 3.9290e-09 at cell 2 (within) | hydrodynamic energy row 8.9320e-09 at cell 2 (within) | carrier balance H2 gated 1.7510e-03 at cell 500 [column max 1.7510e-03 at cell 500, ABOVE] | elemental transport He/H partition gated 2.4280e-03 at cell 335 [column max 1.6300e-02 at cell 2, ABOVE]
pass 6: hydrodynamic mass row 5.3970e-09 at cell 1 (within) | hydrodynamic energy row 9.4120e-09 at cell 2 (within) | carrier balance H2 gated 9.7220e-04 at cell 500 [column max 9.7220e-04 at cell 500, ABOVE] | elemental transport He/H partition gated 7.3530e-04 at cell 336 [column max 1.2000e-02 at cell 2, ABOVE]
```

The carrier relaxation trial ledger of each pass:

```
pass 1: 50 trials, 24 accepted, 26 refused {'the movement bound': 26}
    largest REQUESTED bound measure 1.0380e-02 at cell 266 (trial 19, grow 4.926e+02), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 265 (trial 47, grow 1.338e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 2: 72 trials, 72 accepted, 0 refused {}
    largest REQUESTED bound measure 7.4980e-03 at cell 277 (trial 45, grow 5.598e+07)
    last ACCEPTED bound measure   7.4980e-03 at cell 277 (trial 72, grow 1.414e+12)
    growth multiplier from 1.000e+00 to 1.414e+12
pass 3: 65 trials, 65 accepted, 0 refused {}
    largest REQUESTED bound measure 5.3300e-03 at cell 259 (trial 33, grow 4.314e+05)
    last ACCEPTED bound measure   5.3300e-03 at cell 259 (trial 65, grow 1.861e+11)
    growth multiplier from 1.000e+00 to 1.861e+11
pass 4: 72 trials, 72 accepted, 0 refused {}
    largest REQUESTED bound measure 2.8030e-03 at cell 266 (trial 27, grow 3.788e+04)
    last ACCEPTED bound measure   2.7560e-03 at cell 267 (trial 72, grow 1.414e+12)
    growth multiplier from 1.000e+00 to 1.414e+12
pass 5: 71 trials, 71 accepted, 0 refused {}
    largest REQUESTED bound measure 5.1200e-04 at cell 394 (trial 50, grow 4.251e+08)
    last ACCEPTED bound measure   5.1200e-04 at cell 394 (trial 71, grow 1.414e+12)
    growth multiplier from 1.000e+00 to 1.414e+12
pass 6: no carrier trial recorded
```

The compact state norms (section 7 states the blocks and the scales):

```
# fixed block scales S_b, from pass 1 (the entry state):
#   HeH_ratio  1.000528e+00   (He/H)/HeH, the elemental partition
#   mass_flux  1.763819e-05   r^2 rho v, the wind
#   T          3.699513e+03   temperature [K]
#   x_H2       6.664166e-01   2 n(H2)/n_H, the carrier fraction
#   n_H3p      1.965481e+04   n(H3+) [cm^-3]
#   n_HeII     4.115996e+04   n(He II) [cm^-3]
#   n_e        1.959909e+07   n_e [cm^-3]
# k   ||Y(k+1)-Y(k)||  ||Y(k+1)-Y(k-1)||    HeH_ratio   mass_flux           T        x_H2       n_H3p      n_HeII         n_e
  1->2   6.0268e-01         n/a           8.415e-02   1.362e-02   1.061e-01   2.177e-02   4.962e-02   6.027e-01   4.308e-02
  2->3   1.7047e-01      4.4241e-01       1.054e-02   2.813e-03   1.526e-02   1.676e-02   7.951e-03   1.705e-01   3.597e-02
  3->4   2.2602e-01      3.9149e-01       9.626e-03   2.844e-03   1.940e-02   8.358e-03   9.590e-03   2.260e-01   2.439e-02
  4->5   9.3873e-02      1.6004e-01       7.840e-03   1.122e-03   7.264e-03   1.693e-03   4.300e-03   9.387e-02   1.565e-02
  5->6   1.9757e-03      9.2220e-02       0.000e+00   8.225e-05   7.482e-04   0.000e+00   1.568e-03   1.976e-03   1.088e-03
```

**Verdict: the alternation IS approaching a joint fixed point here, and six
passes were not enough to reach it. This is a budget statement and not a
failure.**

The entry state refuses on the HYDRODYNAMIC ENERGY row, 4.634E-01 at cell 257,
and the solve of pass 1 returns `info = 2`. One composition update cures it:
from pass 2 on the solve returns `info = 0` at every pass and the energy row
sits at 8.05E-09, 1.14E-08, 1.29E-08, 8.93E-09, 9.41E-09 against 1.0E-06, with
the mass row at 4E-09 to 7E-09.

With the wind solvable, nothing throttles the composition. From pass 2 the
carrier relaxation ends on **the fixed point of its own operator**, not on the
bound: 72, 65, 72 and 71 trials, **all accepted, none refused**, the growth
multiplier running up to 1.4E+12 and the returned carrier residual falling to
2E-11. The bound measure stays at 7.50E-03, 5.33E-03, 2.80E-03 and 5.12E-04,
below the bound of 1.0E-02, which is why the progress control never halved it.

The joint residual falls 4.63E+05, 1.39E+03, 3.39E+02, 7.20E+02, 2.43E+02,
9.72E+01: 3.7 decades in five passes, with one excursion at pass 4 which is the
elemental row overtaking the carrier row (at pass 5 the worst gated entry is
the elemental transport row at cell 335). The accepted state displacement
collapses with it, 6.03E-01, 1.70E-01, 2.26E-01, 9.39E-02, 1.98E-03, and the two-pass
norms give a ratio of 0.57, 0.99, 0.50 and 0.96 to the sums they span, so the
state turns as it descends but never undoes a pass: this is a descending
trajectory and not a cycle.

At the ending the state still refuses, carrier 9.722E-04 at cell 500 and
elemental 7.353E-04 at cell 336, which is 97 and 74 times their tolerances. A
geometric extrapolation of the last three passes would put the carrier row
inside 1.0E-05 in a handful of further passes, but that is an extrapolation and
is not measured here.

### 6.4 Case 4, `molecular_scalar_gj1132_wellmixed/HeH0.55`

| pass | hydro info | joint residual (worst gated entry over its tolerance) | worst gated row, cell | element relaxation ending | element residual returned | carrier relaxation ending | carrier residual returned | displacement ADOPTED (element, carrier) | bound: REQUESTED / cell | trials acc/ref | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 1.000e+06 | 2.360e-02 carrier balance H2, cell 306 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 2.620e-01 | 2.830e-02 (0.000e+00, 2.830e-02) | 1.2050e-02 at cell 207 | 19/23 | 128.6 |
| 2 | 2 | 7.870e+07 | 3.300e-02 carrier balance H2, cell 319 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 2.230e-01 | 2.760e-02 (0.000e+00, 2.760e-02) | 1.1220e-02 at cell 246 | 17/21 | 667.8 |
| 3 | 2 | 5.630e+07 | 5.030e-02 carrier balance H2, cell 327 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 1.800e-01 | 2.900e-02 (0.000e+00, 2.900e-02) | 1.3070e-02 at cell 228 | 18/22 | 794.5 |
| 4 | 2 | 5.030e+07 | 7.920e-02 carrier balance H2, cell 336 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 2.490e-01 | 2.830e-02 (0.000e+00, 2.830e-02) | 1.3300e-02 at cell 282 | 16/21 | 441.5 |
| 5 | 2 | 5.190e+07 | 9.500e-02 carrier balance H2, cell 345 | no element relaxation in this run | 0.000e+00 | the movement bound, with the last step inside it | 1.390e-01 | 1.330e-02 (0.000e+00, 1.330e-02) | 5.7110e-03 at cell 284 | 15/20 | 322.8 |
| 6 | 2 | 3.780e+07 | 1.070e-01 carrier balance H2, cell 349 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | NO UPDATE: this is the last pass and no solve is left to consume one | 0.000e+00 | 0.000e+00 (0.000e+00, 0.000e+00) | -- | 0/0 | 386.1 |

front and worst-cell record:
| pass | H2 front x2=0.5 | x2=1e-2 | carrier steady residual, worst cell | displacement worst cell | bound held at cell | trust |
|---|---|---|---|---|---|---|
| 1 | 1.1626 | 1.4987 | 2.620e-01 at cell 1 (r 1.000) | 211 (r 1.181) | 208 (r 1.171) | 1.000e-02 |
| 2 | 1.1684 | 1.5636 | 2.230e-01 at cell 1 (r 1.000) | 241 (r 1.306) | 247 (r 1.339) | 1.000e-02 |
| 3 | 1.1744 | 1.6370 | 1.800e-01 at cell 1 (r 1.000) | 235 (r 1.275) | 239 (r 1.295) | 1.000e-02 |
| 4 | 1.1806 | 1.7325 | 2.490e-01 at cell 1 (r 1.000) | 238 (r 1.290) | 281 (r 1.615) | 1.000e-02 |
| 5 | 1.1806 | 1.7854 | 1.390e-01 at cell 1 (r 1.000) | 243 (r 1.317) | 284 (r 1.648) | 5.000e-03 |
| 6 | 1.1806 | 1.7854 | 1.970e-01 at cell 1 (r 1.000) | 243 (r 1.317) | -- (r --) | 5.000e-03 |

The gated rows of each pass, and the column maxima beside them:

```
pass 1: hydrodynamic mass row 1.9580e-08 at cell 1 (within) | hydrodynamic energy row 1.0000e+00 at cell 246 (ABOVE) | carrier balance H2 gated 2.3630e-02 at cell 306 [column max 2.9850e-01 at cell 1, ABOVE]
pass 2: hydrodynamic mass row 1.2050e-01 at cell 1 (ABOVE) | hydrodynamic energy row 4.7240e-01 at cell 1 (ABOVE) | carrier balance H2 gated 3.3000e-02 at cell 319 [column max 3.2630e-01 at cell 1, ABOVE]
pass 3: hydrodynamic mass row 1.4350e-01 at cell 2 (ABOVE) | hydrodynamic energy row 7.5450e-01 at cell 1 (ABOVE) | carrier balance H2 gated 5.0290e-02 at cell 327 [column max 1.9700e-01 at cell 1, ABOVE]
pass 4: hydrodynamic mass row 6.4100e-01 at cell 1 (ABOVE) | hydrodynamic energy row 9.7160e-01 at cell 1 (ABOVE) | carrier balance H2 gated 7.9190e-02 at cell 336 [column max 6.5320e-01 at cell 1, ABOVE]
pass 5: hydrodynamic mass row 1.1150e-01 at cell 3 (ABOVE) | hydrodynamic energy row 6.1300e-01 at cell 2 (ABOVE) | carrier balance H2 gated 9.5050e-02 at cell 345 [column max 2.2800e-01 at cell 3, ABOVE]
pass 6: hydrodynamic mass row 7.1330e-02 at cell 1 (ABOVE) | hydrodynamic energy row 4.3480e-01 at cell 1 (ABOVE) | carrier balance H2 gated 1.0660e-01 at cell 349 [column max 2.1160e-01 at cell 1, ABOVE]
```

The carrier relaxation trial ledger of each pass:

```
pass 1: 42 trials, 19 accepted, 23 refused {'the movement bound': 23}
    largest REQUESTED bound measure 1.2050e-02 at cell 207 (trial 12, grow 8.650e+01), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 208 (trial 38, grow 2.819e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 2: 38 trials, 17 accepted, 21 refused {'the movement bound': 21}
    largest REQUESTED bound measure 1.1220e-02 at cell 246 (trial 14, grow 6.487e+01), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 247 (trial 36, grow 1.253e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 3: 40 trials, 18 accepted, 22 refused {'the movement bound': 22}
    largest REQUESTED bound measure 1.3070e-02 at cell 228 (trial 13, grow 1.297e+02), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 239 (trial 37, grow 1.879e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 4: 37 trials, 16 accepted, 21 refused {'the movement bound': 21}
    largest REQUESTED bound measure 1.3300e-02 at cell 282 (trial 13, grow 1.297e+02), refused
    last ACCEPTED bound measure   1.0000e-02 at cell 281 (trial 32, grow 6.682e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 5: 35 trials, 15 accepted, 20 refused {'the movement bound': 20}
    largest REQUESTED bound measure 5.7110e-03 at cell 284 (trial 12, grow 2.883e+01), refused
    last ACCEPTED bound measure   5.0000e-03 at cell 284 (trial 31, grow 4.454e-03)
    growth multiplier from 1.000e+00 to 9.766e-04
pass 6: no carrier trial recorded
```

The compact state norms (section 7 states the blocks and the scales):

```
# fixed block scales S_b, from pass 1 (the entry state):
#   HeH_ratio  1.000000e+00   (He/H)/HeH, the elemental partition
#   mass_flux  2.204999e-05   r^2 rho v, the wind
#   T          5.033268e+03   temperature [K]
#   x_H2       9.327891e-01   2 n(H2)/n_H, the carrier fraction
#   n_H3p      7.803788e+04   n(H3+) [cm^-3]
#   n_HeII     1.604197e+06   n(He II) [cm^-3]
#   n_e        1.700169e+07   n_e [cm^-3]
# k   ||Y(k+1)-Y(k)||  ||Y(k+1)-Y(k-1)||    HeH_ratio   mass_flux           T        x_H2       n_H3p      n_HeII         n_e
  1->2   9.9079e-01         n/a           0.000e+00   9.908e-01   7.600e-02   2.756e-02   1.653e-01   4.299e-01   4.447e-02
  2->3   3.9470e-01      9.9278e-01       0.000e+00   2.297e-03   1.979e-02   2.901e-02   8.883e-02   3.947e-01   4.537e-02
  3->4   2.7357e-01      5.9600e-01       0.000e+00   4.439e-03   2.648e-02   2.827e-02   1.038e-01   2.736e-01   3.577e-02
  4->5   1.0383e-01      3.5900e-01       0.000e+00   2.793e-03   2.997e-02   1.326e-02   4.387e-02   1.038e-01   5.197e-02
  5->6   3.2038e-02      1.0746e-01       0.000e+00   8.053e-04   1.167e-02   0.000e+00   3.204e-02   9.628e-03   2.991e-02
```

**Verdict: the alternation is NOT approaching a joint fixed point; the gated
carrier row rises monotonically and its refusing cell runs outward with the H2
front.**

The entry state refuses on the hydrodynamic energy row, 1.000E+00 at cell 246,
with the mass row within at 1.958E-08. From pass 2 the composition update has
taken the wind as well: `info = 2` at every pass, the mass row at 1.2E-01 to
6.4E-01 and the energy row at 4.3E-01 to 9.7E-01, and the energy row's worst
cell moves from 246 to cell 1 and stays there. `||Z(2) - Z(1)|| = 9.91E-01`,
again carried entirely by the mass flux block.

The gated carrier row rises at every pass without exception, 2.363E-02,
3.300E-02, 5.029E-02, 7.919E-02, 9.505E-02, 1.066E-01, and its cell walks
outward 306, 319, 327, 336, 345, 349, which is r = 1.95 to 3.02. The x2 = 1E-02
front moves outward with it, 1.4987, 1.5636, 1.6370, 1.7325, 1.7854. The
carrier relaxation ends on the movement bound in every update pass, with 15 to
19 accepted and 20 to 23 refused trials, and the bound is held at cells 208 to
284.

So this case has neither a flat maximum nor a moving worst cell masking a fall:
the row rises by a factor 4.5 over six passes while the wind is out of
solution. It is the same two-factor failure as case 2 with a front that has not
finished moving.

This run also printed
`WARNING: the written state is NOT the state the gates accepted`.

---

## 7. The compact state norms: what they are and what they exclude

### 7.1 The vector, the blocks and the scales

`Z(k)` is the state at the END of outer pass `k`, after both relaxations, which
is the state the next hydrodynamic solve is entered with. It is read from
`output/diffusion_pass_profiles.txt`, written by
`write_diffusion_pass_profile` (`EXHALE_main.f90:8258`, under
`EXHALE_DIFFUSION_CHECK=1`), which appends one block of 500 rows per pass.

The seven exported blocks:

| block | quantity | units |
|---|---|---|
| `HeH_ratio` | `(He/H)/HeH`, the elemental partition against the run's reservoir | dimensionless |
| `mass_flux` | `r^2 rho v`, the wind | code units |
| `T` | temperature | K |
| `x_H2` | `2 n(H2)/n_H`, the carrier fraction | dimensionless |
| `n_H3p` | `n(H3+)` | cm^-3 |
| `n_HeII` | `n(He II)` | cm^-3 |
| `n_e` | electron density | cm^-3 |

**What the export does NOT carry, said plainly.** The density and the velocity
appear only through their product `r^2 rho v`, and there is no pressure column.
The norms below are therefore norms of these seven blocks and not of the full
conserved vector. That is the whole of what the current code exports per pass,
and adding to it would be a source change this audit did not make.

**The block scales are FIXED for a case**: `S_b = max_j |Z_b(1, j)|`, the entry
state's own column maximum, printed at the head of every norm block in section
6 so that the divisor is on the page. The compact norm is

```
||Z(a) - Z(b)||  =  max over blocks b of  ( max over cells j |Z_b(a,j) - Z_b(b,j)| / S_b )
```

and the value of each block is printed beside it, so a reader sees which block
carries the number.

### 7.2 What the two differences decided

MEASURED. `r` is the ratio of the two-pass difference to the sum of the two
one-pass differences that span it. Aligned steps give `r` near 1; a two-cycle,
in which pass `k+1` undoes pass `k`, gives `r` near 0.

| case | `\|\|Z(k+1)-Z(k)\|\|` over the five steps | `\|\|Z(k+1)-Z(k-1)\|\|` | `r` | the block that carries the one-pass norm |
|---|---|---|---|---|
| 1 | 2.108E-02, 2.126E-02, 1.082E-02, 1.071E-02, 4.025E-03 | 4.233E-02, 3.208E-02, 2.148E-02, 1.071E-02 | 1.00, 1.00, 1.00, 0.73 | `x_H2` at the first four, `n_HeII` at the last (where `x_H2` is exactly zero because pass 6 takes no composition update) |
| 2 | 9.848E-01, 3.546E-02, 2.142E-02, 1.810E-02, 1.747E-02 | 9.869E-01, 5.680E-02, 3.750E-02, 3.231E-02 | 0.97, 1.00, 0.95, 0.91 | `mass_flux` at the first step, then `n_H3p` and `n_HeII` |
| 3 | 6.027E-01, 1.705E-01, 2.260E-01, 9.387E-02, 1.976E-03 | 4.424E-01, 3.915E-01, 1.600E-01, 9.222E-02 | 0.57, 0.99, 0.50, 0.96 | `n_HeII` at every step |
| 4 | 9.908E-01, 3.947E-01, 2.736E-01, 1.038E-01, 3.204E-02 | 9.928E-01, 5.960E-01, 3.590E-01, 1.075E-01 | 0.72, 0.89, 0.95, 0.79 | `mass_flux` at the first step, then `n_HeII` and `n_H3p` |

**No case shows a two-cycle.** The smallest `r` anywhere is 0.50, at case 3
steps 3 and 4, where the state is genuinely turning while it descends; every
other value is 0.57 or above and most are within a few per cent of 1. Nothing
here is near the `r` close to zero that an oscillating pair of states would
give, so the flat scalars of cases 2 and 4 cannot be read as a hidden cycle.

The three statements the block breakdown adds:

- in cases 2 and 4 the FIRST step is carried by `mass_flux` at 0.99 of its own
  column maximum, which is the composition update moving the wind, and after
  it the wind is out of solution for good;
- in case 1 the one-pass norm is carried by `x_H2` itself and its value falls
  in step with the halving of the movement bound at pass 4, not with any
  residual;
- in case 3 it is carried by `n_HeII` throughout and falls by two and a half
  decades, which is the only sequence of the four that looks like a
  contraction.

Case 4 is the measurement that earns rev9 section 9.2's warning its place. Its
accepted displacement falls by a factor 31 over five passes while its gated
carrier row rises by a factor 4.5. A displacement is not a residual and was
never read as one here.

---

## 8. The outer cells: the carrier row of cells 499 and 500, term by term

### 8.1 What the file's four face entries are, and how they were used

INSPECTED, `diffusive_photochemistry.f90:5754-5756` and `:5860-5865`: with
`Kj = 1/(cv(j) R0cb)`, `sL = fa(j-1) R0sq` and `sR = fa(j) R0sq`, the file holds
`-Kj sL J(j-1)` and `Kj sR J(j)`. The area and the cell volume are ALREADY
applied. This audit therefore applied no geometry to them and summed them only
within one cell, where `fdif_in + fdif_out` reproduces the `diffusive` column
and `fadv_in + fadv_out` the `advective` column. MEASURED: they do, at every
cell inspected. No volume-weighted sum over cells was formed from them; that is
the export of rev9 section 13.2 and is not this file.

The H2 channel record is read with its names and its count from the file's own
header: **15 values**, `R15_3body R9_H2p_H R6_H3p_e R11_H3p_H P_H2_photo
LW_photdis R10R13_Hp R12_thermal R14_edis R8_H2p_H2 R17R23_Hep R18_HeHp
HeI23S_H2 oxy_prod oxy_loss`, the last two zero in all four cases because
`oxygen_chemistry` is F.

### 8.2 Case 1 at the entry state: the row at cells 498 to 500

MEASURED, from the Mode R probe of the generation `g0004_...79b42a03`:

| cell | r [R_p] | measure | residual | scale | diffusive | advective | net source | `fdif_in` | `fdif_out` |
|---|---|---|---|---|---|---|---|---|---|
| 498 | 28.073 | 7.5684E-04 | +4.1035E-05 | 5.4218E-02 | +1.4600E-03 | -1.4460E-03 | -2.7000E-05 | +3.7190E-03 | -2.2590E-03 |
| 499 | 28.549 | **8.1496E-04** | +3.7676E-05 | 4.6231E-02 | +1.3245E-03 | -1.3114E-03 | -2.4552E-05 | +2.1514E-03 | -8.2683E-04 |
| 500 | 29.031 | 9.8108E-05 | -3.9425E-06 | 4.0185E-02 | +7.9275E-04 | -8.1905E-04 | -2.2361E-05 | +7.9275E-04 | **0.0000E+00** |

The row is `transport - net source` and it closes: at cell 499,
`(1.3245E-03 - 1.3114E-03) - (-2.4552E-05) = 3.765E-05` against the residual
3.7676E-05 the file reports. The net source is photodissociation and nothing
else: of a total loss 2.4654E-05 the two photolysis channels carry
`LW_photdis 1.3473E-05` and `P_H2_photo 1.1169E-05`, and the largest production
channel `R9_H2p_H` is 9.87E-08, three decades below.

So the refusal at cell 499 is the failure of a transport divergence that is
itself the difference of two numbers a hundred times larger to cancel against
the photodissociation of H2 at 28.5 R_p.

### 8.3 Case 2 at the entry state: cell 500

MEASURED: measure 7.3714E-02, residual +2.1287E-04, scale 2.8878E-03,
`fdif_in = +2.4696E-04`, `fdif_out = 0.0000E+00`, advective sum -3.6691E-05,
net source -2.5961E-06. The residual is, to two figures, the inner diffusive
inflow that the outer face does not carry away and the advective divergence
does not absorb: `2.4696E-04 - 3.669E-05 + 2.60E-06 = 2.13E-04`.

In case 1's cell 500 the same zero outer face is present and the advective
divergence happens to over-compensate it, which is why case 1's maximum is at
cell 499 and case 2's at cell 500.

### 8.4 Does the refusal follow the outer boundary or the wind?

**Both, and the two can be told apart.** MEASURED over the whole gated window
`r >= 1.20`:

| case, state | gated maximum | its cell, r | the shape of the gated profile |
|---|---|---|---|
| 1, entry | 8.1496E-04 | 499, 28.55 | rises MONOTONICALLY outward over the last ~150 cells: 2.97E-06 at r = 1.86, 1.36E-05 at 3.05, 1.20E-04 at 8.58, 3.71E-04 at 19.1, 8.15E-04 at 28.5; cell 500 drops to 9.81E-05 |
| 1, ending | 2.8931E-03 | 500, 29.03 | the same rise, and now cell 500 is the maximum, a factor 3.5 above cell 499 |
| 2, entry | 7.3714E-02 | 500, 29.03 | the same monotone outward rise, 2.32E-03 at r = 1.23 to 4.66E-02 at cell 499, with cell 500 a further factor 1.6 above it |
| 2, ending | 7.6679E-02 | 500, 29.03 | unchanged in shape |
| 3, entry | 8.9224E-03 | 261, 1.43 | an INTERIOR plateau at the H2 front; the outer cells are at 1E-06, four decades below |
| 3, ending | 9.7219E-04 | 500, 29.03 | the interior plateau is gone and the maximum is the last cell, **ten times its own neighbour** (cell 499 is 8.97E-05): a one-cell boundary feature |
| 4, entry | 2.3629E-02 | 306, 1.95 | an interior peak, with an outer rise reaching 1.99E-02 at cell 500, within 16 per cent of it |
| 4, ending | 1.0660E-01 | 349, 3.02 | a localized interior peak falling off on both sides; cell 490 is 1.45E-02, a factor 7 below |

**The mechanism of the outward rise is measured and it is a normalization, not
a growing imbalance.** In case 1 the absolute residual FALLS monotonically
outward, 1.36E+00 at cell 225 to 3.77E-05 at cell 499, while the row's own
physical scale falls faster, 9.78E+04 to 4.62E-02. The ratio the certification
gates on therefore grows outward although the equation is being violated less
and less in absolute terms. Case 2 is the same, residual 1.26E+02 to 1.67E-04
against a scale 5.43E+04 to 3.59E-03.

**The last cell is separated from that rise by the closure.** The outer face
rule continues the last physical value into the ghost faces
(`diffusive_photochemistry.f90:4167-4169`) and `carrier_outflow_ghost` does the
same for the state (`:4181-4188`), and the arithmetic consequence is visible in
every one of the eight files: `fdif_out = 0.0000000E+00` exactly at cell 500.
Whether that makes cell 500 the maximum (cases 2 and 3 at their ending, case 1
at its ending) or a local minimum (case 1 at entry) depends on whether the
advective divergence of that cell happens to cancel the truncated inflow.

**Case 4 is the one that follows the wind.** Its ending maximum at cell 349,
r = 3.015, has a diffusive divergence -9.3279E-03 against a net source
-4.0066E-03, and the loss is dominated by the He+ and H2 channel
`R17R23_Hep = 4.3731E-03` of a total 4.4252E-03. That is a chemistry and
transport imbalance in the interior of a helium-rich wind, it moves outward with
the H2 front, and it has nothing to do with the outer face.

**What this does NOT say.** A residual that follows a boundary localizes the
question; it does not say which boundary law is right, and it is not a reason to
drop the row from certification. A radius experiment, holding the inner mesh and
extending the domain outward, would separate the closure from the physics. It
was not run here and is not in this phase's budget.

---

## 9. Answers to the questions the phase asks

**1. Is the alternation approaching a joint fixed point?** In case 3 yes, and
the ending is a budget statement: the joint residual falls 3.7 decades in five
passes, nothing throttles the composition, both relaxations sit at their own
fixed points from pass 2, and the accepted state displacement is a factor 305
smaller at the last step than at the first.
In cases 1, 2 and 4 no, and the record of each pass says which of the plan's five
explanations applies to each: the movement bound in case 1, the splitting
itself in cases 2 and 4, with a moving front on top of it in case 4. In no case
was a flat maximum taken as the answer.

**2. Is a movement bound, an inner accuracy, an outer closure or a domain
problem present?** All four, in different cases, and they are now separated:

| | case 1 | case 2 | case 3 | case 4 |
|---|---|---|---|---|
| movement bound attained | every pass, SATURATED | every pass | pass 1 only | every pass |
| inner relaxation reached its own fixed point | element yes every pass, carrier never | no element relaxation in this case; carrier never | element yes every pass, carrier yes from pass 2 | no element relaxation in this case; carrier never |
| the wind stayed in solution | yes, all six passes | no, lost at pass 2 | no at pass 1, yes from pass 2 | no, lost at pass 2 |
| the refusal sits at the outer closure | yes | yes | yes at the ending | no, in the wind |

**3. What would decide what this measurement does not.** Three things are
named and none of them was run here:

- whether the outer refusal of cases 1, 2 and 3 is the closure or the physics:
  the radius experiment of rev9 section 9.3, inner mesh held, domain extended
  outward, fluxes and inner profiles compared;
- whether case 3 reaches its joint fixed point: more passes on the same
  recipe. Its last two steps fall by factors 3.0 and 2.5, so a handful more
  passes would settle it, but that is an extrapolation and not a measurement;
- whether the splitting failure of cases 2 and 4 is the splitting or the bound:
  only the coupled route of rev9 section 9.4 separates them, and that trial is
  phase 7 and has its own prerequisite (the coupled Stage B, phase 6b).

**4. The base and the wind: one phenomenon or two?** The earlier phases put the
ATOMIC imbalance at cells 1 and 2 (`docs/lhs1140b_conservation_audit_20260921.md`)
and the M2 refusal at the molecular base reservoir
(`docs/lhs1140b_m2_conversion_audit_20260921.md`); the four cases here refuse at
cells 246 to 500. **They are most likely the same KIND of defect at the two ends
of the domain, and not one phenomenon carried from one end to the other.** The
measured ground for saying so is narrow and is stated as such: at cell 1 the
carrier row's INNER diffusive face is exactly zero in all four cases
(`fdif_in = -0.0000000E+00`, MEASURED in all eight files), just as the OUTER diffusive face is
exactly zero at cell 500, so both end cells carry a one-sided diffusive
divergence that the interior cells do not. What separates them is that the base
cell has a large advective inflow to balance against (-1.57E+06 at cell 1 of
case 1, against a scale 3.38E+06) while the outer cell has a scale seven
decades smaller, so the same one-sided truncation shows up at the base as a
1.7E-02 row and at the edge as a 8.2E-04 row that the gate still refuses.
Case 4's refusal is neither: it is an interior front. This is a statement about
where the rows sit and about one shared arithmetic feature; it is NOT a
demonstration that the same physical defect produces both, and that
demonstration is the radius experiment together with the base work already open
in rev9 section 16.

---

## 10. Artifacts

All raw artifacts are in the scratch directory of this audit,

```
/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/m1_alternation/
```

with `artifacts/README_artifacts.md` describing the layout. Per case:
`<case>/run.log` and `run.err` (the Mode C trajectory),
`<case>/output/diffusion_pass_profiles.txt` and `output/carrier_row_terms.txt`,
`<case>_eval/run_eval.log` and `<case>_eval/output/carrier_row_terms.txt` (the
Mode R probe), and under `artifacts/` the derived tables of sections 6 to 8.
`SNAPSHOTS.txt` carries the generation and the md5 of every input and state file
as set up; `modec_exit.txt` the exit status and end time of each run. The
readers are in `analysis/` and take the channel count and names from the file's
own header.
