# LHS 1140 b, the three atomic P2 cases: a longer alternation from their own stored states

The experiment `docs/PLAN_20260920_rev9.md` section 17.3 item 2 names as the
next one, carried out on 2026-09-21. It is not one of the numbered phases of
that document's section 15 table, whose entry 10 is Stage F; it follows
phase 9 (`docs/lhs1140b_p2_base_stage_cd_20260921.md`) and reads the same
three cases phase 6 read for the molecular side
(`docs/lhs1140b_alternation_audit_20260921.md`).

**Mode C throughout for the trajectories** (rev9 section 3): every pass table
below advances a retained trajectory under the current solver from a named
stored generation, with a budget declared before the run and an ending read
out of the run's own words. Two kinds of number stand beside them and say so
where they stand: the entry certification each run prints before it takes a
pass, which is one evaluation of the current operator on an immutable
checkpoint and is therefore Mode R, and the values read out of stored
manifests, indexes and `not_solved.md` notes, which are a historical audit.

Every number is labeled MEASURED (this document ran it), READ (from a named
file) or INSPECTED (source, with `file:line`). Nothing was published, no
golden was refreshed, and nothing under `LHS1140b/models/` or
`backup/regression/` was written: each case was copied into a scratch
directory and run there.

**All three stored states were written under the superseded boundary model
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2`, and the binary
solves `..._ghost_fixed_point_seed_reservoir_row_v3`, so every run here is the
CURRENT OPERATOR CONTINUING AN OLD STATE.** The run says so itself and rebuilds
the boundary from the physical column and this run's reservoir before it takes
anything (READ, the `load_IC` note in each `run.log`). The statement is
repeated in every table.

---

## 1. Verdict, first

**None of the three approaches a certified state, and the three fail in two
different ways.** In the two scalar cases the alternation is CONVERGING, and
it is converging to a state that is not a solution: the accepted state
displacement halves every pass, the gated elemental row ends BELOW its
tolerance, and the hydrodynamic energy row rises monotonically to 29 and 64
times its own tolerance and settles there. In the profile case the
hydrodynamic solve never leaves its entry state at all, and the run ends
itself on its own no-progress rule after 5 passes.

| case | verdict after the budget | the evidence |
|---|---|---|
| 1, `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | **converging, to a state that is not a solution** | pass 1 certifies every hydrodynamic row (`info = 0`); from pass 2 on every solve returns `info = 2` and gives back, to four printed digits, the state it was entered with. The gated elemental row falls 7.218E-05, 4.440E-05, 2.219E-05, 1.038E-05, turns up to 2.273E-05 and then falls again to **8.190E-06, below its 1.0E-05 tolerance**, while the energy row rises monotonically 2.559E-07 to 2.911E-05 with lifts that decay geometrically at 0.71 per pass. Two refusing entries at the ending, both hydrodynamic |
| 2, `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | **the same, more cleanly** | the same shape: pass 1 `info = 0`, `info = 2` after; the gated elemental row 1.200E-04 to **5.210E-06**, below tolerance, with one excursion at passes 4 and 5; the energy row 8.946E-07 to 6.387E-05 with lifts decaying at 0.65 per pass. Two refusing entries at the ending, both hydrodynamic |
| 3, `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | **not approaching: the hydrodynamic solve cannot move the state** | every one of the five solves starts at `||R|| = 1.949E+00`, ends `info = 2` at the same number, and hands back a state whose mass row is 6.925E-01, momentum 7.77E-07 and energy 1.949E+00, identical to all printed digits at every pass. The composition does move, and shrinks, and changes nothing. The run ends ITSELF at pass 5: `the joint distance ... has not fallen in 3 consecutive passes` |

**The question the phase asked is answered: the elemental row does keep
falling, and reaching it is not enough.** On the two scalar cases the row that
phase 9 saw refusing at 7.2E-05 after one pass goes under its tolerance by
pass 11 or 12. What stands in its place is the hydrodynamic energy row, which
the composition update lifts by an amount that decays geometrically and which
the hydrodynamic solve, from pass 2 onward, hardly takes back: the recovered
part is at most 1.6 per cent of the lift at every pass of case 1 and at most
11.2 per cent at every pass of case 2, and is exactly zero to the four digits
the solver prints on 11 of the 22 passes. The alternation therefore has a joint
fixed point on these two cases and that fixed point is about 31 and 65 times
the energy tolerance.

---

## 2. The classification, the budget and the ending

### 2.1 Classification

| activity of this document | label |
|---|---|
| the three trajectories of section 6 | **Mode C**, a retained trajectory on the current branch, with a declared budget and a recorded ending |
| the entry certification each run prints before it takes a pass (`Restart intent: stationary`: the loaded state is measured as it stands, no CFL step) | **Mode R**, one evaluation of the current operator on an immutable checkpoint |
| the manifests, indexes, `ENDING` and `not_solved.md` notes quoted for context | historical audit, carrying the HISTORICAL executable identity |

No table mixes them.

### 2.2 The budget, declared before the runs

| limit | value | why |
|---|---|---|
| outer passes per case | **12**, `EXHALE_OUTER_PASSES=12`; three cases, so 36 outer passes is the whole budget | twice the six of phase 6, because the question is whether the elemental row keeps falling past the two passes phase 9 took |
| hydrodynamic solves per case | at most 12, one per pass (INSPECTED, `EXHALE_main.f90:7102-7126`) | |
| Newton iterations per inner solve | **3000, the case's own setting.** No case names a cap, and `EXHALE_JFNK_MAXIT` was UNSET in every run, so each solve entered at `jfnk_outer_iterations_default = 3000` (INSPECTED, `EXHALE_main.f90:347`; READ, each run's `(JFNK) outer iteration cap 3000, from the caller`). The brief forbids bounding a solve with `EXHALE_JFNK_MAXIT`, and INSPECTED at `EXHALE_main.f90:6993-7010` it is the RUN's cap in any case: after the first solve under a named cap the run refuses to enter the solver again |
| wall ceiling per case | 3600 s, `timeout 3600`, a SAFETY ceiling only | a stop at the ceiling is an external termination, recorded as one and never converted into an ending of the alternation (rev9 section 12). **It never fired**: the largest wall of the three was 1703 s |
| threads | `OMP_NUM_THREADS=8`, `OPENBLAS_NUM_THREADS=1`; the three cases ran one after another and never concurrently | |

The estimate the budget was set from: phase 9 MEASURED 53.3 s for outer pass 1
and 91.5 s for pass 2 of case 1 at 8 threads
(READ, `docs/lhs1140b_p2_base_stage_cd_20260921.md` section 4). Twelve passes
at that rate is about 1200 s, which is what set the 3600 s ceiling at three
times the estimate.

### 2.3 What was spent, and how each trajectory ended

MEASURED. All three exited 0. None reached the wall ceiling, so **no ending in
this document is a timeout**.

| case | outer passes spent | hydrodynamic solves | JFNK iterations | Krylov cycles | Krylov products | residual evaluations | element relaxation steps | wall [s] | ending, in the run's own words |
|---|---|---|---|---|---|---|---|---|---|
| 1, `kzz1e9/HeH9.7` | **12 of 12** | 12 | 691 | 701 | 813 | 26443 | 336 | 1280.2 | `the stationary outer iteration spent its budget of 12 passes without an accepted state: 2 active equations refuse it`; then `the stationary solve returned info = 1`, the state written `certified=F`, reason `no stationary claim`, exit 0 |
| 2, `kzz1e9/HeH2.13` | **12 of 12** | 12 | 502 | 512 | 714 | 13087 | 310 | 572.5 | the same |
| 3, `kzzprofile/HeH9.7` | **5 of 12** | 5 | 300 | 305 | 547 | 10951 | 146 | 1702.7 | **the alternation's own rule, not the budget**: `outer pass 5: REFUSED -- the joint distance of the state (its largest entry over its own tolerance) has not fallen in 3 consecutive passes; the alternation is not approaching a joint fixed point at these step lengths`; then the same `info = 1` write, exit 0 |

Total spent: **29 outer passes, 29 hydrodynamic solves, 1493 JFNK iterations,
2074 Krylov products, 50481 residual evaluations, 792 element relaxation
steps**, 3555 s of pass wall clock, between 18:04 and 19:04 KST on `lart4`.
No run printed `WARNING: the written state is NOT the state the gates
accepted` (MEASURED, zero occurrences in all three logs), which the two
`wellmixed` runs of phase 6 did print.

`EXHALE_CARRIER_DEBUG=1` was set in every run from the start, as the brief
requires. It produced **zero trials** in all three, and that is correct rather
than a loss: `carrier_transport` is `F` in all three cases (READ, each
`EXHALE_resolved.out`), so there is no carrier relaxation to trace and each
pass says `the run carries no relaxed carrier` where the carrier measure would
stand.

---

## 3. Identity

### 3.1 The operator

| field | value |
|---|---|
| binary | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x`, md5 `89ec67149aa452eea3deeb19052f83d1` (MEASURED), the tree's binary of record; copied into each run directory, each copy MEASURED equal |
| repository | `93eed8667483e11dffe51b46008d40081292e104`, working tree dirty with 63 entries (MEASURED, `git log -1` and `git status --porcelain`) |
| source digest | **not quoted, and the reason is recorded.** Other workers are editing this tree while these runs were made: `src/` holds 245 `.f90`/`.inc` files now against the 246 the phase 9 record counted this morning, and four source files carry mtimes later than the phase 9 runs. A digest taken now would not describe the sources this binary was built from. The binary md5 is the anchor; every `file:line` below was read from the tree as it stands at the time of writing |
| compiler identities in the image | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`, `GCC: (conda-forge gcc 16.2.0-4) 16.2.0`, `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` (MEASURED, `strings -a`) |
| linked libraries | `libopenblas.so.0`, `libgfortran.so.5`, `libgomp.so.1`, `libgcc_s.so.1`, `libquadmath.so.0`, all from `/opt/miniconda3/lib`; `libm`, `libmvec`, `libpthread`, `libdl`, `librt`, `libc` from the system (MEASURED, `ldd`) |
| host, threads | `lart4`, 72 cores; `OMP_NUM_THREADS=8`, `OPENBLAS_NUM_THREADS=1`, both echoed by each run |
| `EXHALE_RESID_QUAD` | UNSET in every run, so the assembly selector of rev9 section 13.6 stands at its default. It is in no configuration record and is stated here because only the setter can state it |
| the complete environment of every run | `PATH`, `HOME`, `OMP_NUM_THREADS=8`, `OPENBLAS_NUM_THREADS=1`, `EXHALE_OUTER_PASSES=12`, `EXHALE_DIFFUSION_CHECK=1`, `EXHALE_CARRIER_DEBUG=1`, and nothing else; every run was launched under `env -i` and wrote its own `ENV.txt` beside it. `EXHALE_JFNK_MAXIT`, `EXHALE_PTC_DTAU0`, `EXHALE_CARRIER_TRUST`, `EXHALE_CARRIER_TRUST_HOLD`, `EXHALE_DIFF_OMEGA`, `EXHALE_TRUST_REGION`, `EXHALE_ELEM_DIAG` and `EXHALE_LINEAR_ROWS` all UNSET |

`EXHALE_PTC_DTAU0` being unset, every run entered its continuation at the
default `dtau0 = 1.0` and printed it (READ,
`(EXHALE_main) stationary restart: dtau0 = 1.00E+00`; INSPECTED,
`EXHALE_main.f90:5931-5934`).

### 3.2 The three snapshots

Each case was copied into its own run directory and the generation its index
names `latest_complete` was placed as the `_IC` pair. Every state md5 was
MEASURED equal to the md5 its own manifest records.

| case | generation (READ, `state_index.json` `latest_complete`) | `Hydro_ioniz.txt` md5 (MEASURED = READ) | `Ion_species.txt` md5 (MEASURED = READ) | the binary that WROTE it (READ, manifest) |
|---|---|---|---|---|
| 1, `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `g0002_20260919T004843Z_f7485b14` | `25e11562ec3cd1f6d204bf36f8e503eb` | `78ac49d7e1ce90a3a4f271785400efda` | `EXHALE_7670f310.x`, md5 `7670f31031fb4db91d27b44cb0da6f70` |
| 2, `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | `g0002_20260919T004834Z_d59ec9ff` | `741dc1415efcebfb4f2e13e818ac6b5a` | `60ab4afeff02fcb97e5ab37c39edb754` | the same |
| 3, `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | `g0004_20260919T094145Z_259fe9c3` | `c19b3df722fb19c91662c633bc5f1bc1` | `9f04c20d6c365f87067c081ec2054b53` | `EXHALE_2c3b0acc.x`, md5 `2c3b0acc9983aec03bb4f844294fed18` |

`latest_certified` is `none` in all three indexes (READ): no generation of any
of the three carries a certificate covering its own bytes. The boundary model
of all three state headers is
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2` (READ, the
`# boundary_model` line, which is part of the hashed bytes).

**The Mode R identity check, and what it reproduces.** Each run measures its
loaded state before it takes a pass, and on case 1 that measurement returns
phase 9's numbers to every printed digit: mass row 9.806E-01 at cell 2, energy
row 1.071E+00 at cell 1, momentum row 8.724E-06 at cell 1, the base face at
82.180718 wind means. The FIRST pass then reproduces phase 9's pass 1 to every
printed digit as well, `info = 0` with mass 2.44E-08, momentum 1.02E-14,
energy 2.56E-07 and the gated elemental row 7.22E-05 at cell 295, and the
SECOND pass reproduces phase 9's pass 2, `info = 2` with 4.44E-05 at cell 284
and energy 6.48E-06, **although phase 9 capped that solve at 60 iterations and
this one was uncapped and ran 59**. The stagnation is therefore the solve's own
and not the cap: both runs print
`STAGNATED: neither the judged distance nor the merit improved in 40
iterations, monotone search included`.

### 3.3 The effective configuration

READ from each run's own `EXHALE_resolved.out` (the derived route keys, which
`write_setup_report.f90:994` onward writes) and `EXHALE_setup.out` (the prose
keys the resolved record does not carry).

| key | case 1 | case 2 | case 3 |
|---|---|---|---|
| `he_diffusion` | T | T | T |
| `carrier_transport` | **F** | **F** | **F** |
| `carrier_in_newton`, `carrier_newton_on_stall`, `ionization_transport`, `oxygen_chemistry` | F | F | F |
| `well_balanced` | T | T | T |
| `lower_profile_present` | F | F | **T** |
| `base_inp_present` | F | F | F |
| `HeH_number_ratio` | 9.7 | 2.13 | 9.711073282076663 |
| grid cells | 500, Mixed, base grid 1.9999999494757503e-4 x 50 | the same | the same |
| numerical flux, reconstruction | ROE, PLM | ROE, PLM | ROE, PLM |
| secondary ionization | Immediate, on from step 0 | the same | the same |
| `Resid tol` (the SOLVER's target, not a certification tolerance) | 5.0e-5 | 5.0e-5 | **4.0e-6** |
| `K_zz` | scalar `He_Kzz = 1.0e9` | the same | from the profile, `1.000E+09 to 1.000E+09 cm^2/s` on the grid |
| metals | none (`No metals.inp found`) | the same | **the profile's elemental reservoirs**, `C/H 2.778E-04 O/H 9.028E-07 N/H 8.191E-05` in the state header, and the closure row is `eliminated-species closure System_HeH_TR_metals` |
| boundary model of the loaded state | `..._contact_upwind_v2` | the same | the same |
| boundary model this run solves | `..._ghost_fixed_point_seed_reservoir_row_v3` | the same | the same |

**No case input was edited.** All three state `Restart intent: stationary`
themselves and all three name their spectrum relative to the case directory,
so a symbolic link `LHS1140b/sed` inside the scratch tree was enough and the
`input.inp` of each run directory is byte-identical to the case's own (MEASURED
md5: case 1 `e43c6697786711fb916d3b8eb3c113aa`, case 2
`c355d9827b63e92b3d2f711bf73e26a7`, case 3 `6fc953139ef22e7e88397b465d997b64`).

---

## 4. The route, and where each recorded quantity comes from

INSPECTED, `EXHALE_main.f90`, the outer pass loop from `:7100`. One pass is:

1. the hydrodynamic solve, `solve_steady_jfnk` at `:7122`;
2. the refresh: `U_to_W`, `ioniz_eq`, and `T = p/(n_tot + n_e)` of the
   composition beside it (`:7128-7146`);
3. the joint test: `assemble_residual` then `certification_evaluate` at
   `:7172-7178`, on that refreshed state. This is the acceptance criterion and
   it is what the pass line reports;
4. the progress control (`:7345` onward), which halves the under-relaxation
   factor `omega` when neither the coupled nor the composition residual fell;
5. `relax_element_composition` at the face mass flux of this state
   (`:7475-7481`), under `he_diffusion`;
6. `relax_photochemical_composition` (`:7540-7544`), which **does not run in
   any of these three cases**, because `carrier_transport` is `F`;
7. the log line and the diagnostics.

So the alternation of these three cases has ONE composition half, the element
relaxation. That is the mirror image of the two `wellmixed` cases of phase 6,
whose one half was the carrier relaxation.

**REQUESTED and ADOPTED.** `relax_element_composition` returns two numbers
(INSPECTED, `binary_element_diffusion.f90:1199` onward and the header at
`:1265-1280`): the `map_distance`, the UNDAMPED distance
`max |X_relaxed - X_old| / X_base` over every element the pass relaxed, and
the `displacement`, the damped step actually applied,
`X <- X_old + omega (X_relaxed - X_old)`. The header states in as many words
why the two are reported separately, so that `omega` cannot make the distance
look small by moving less. In the tables below the inner map distance is
therefore the REQUESTED displacement and the applied displacement is the
ADOPTED one, and their ratio is `omega` at every pass of all three cases
(MEASURED: 2.620E-05 and 1.310E-05 at `omega = 0.500`, 2.120E-03 and
5.310E-04 at `omega = 0.250`, 1.730E-03 and 2.170E-04 at `omega = 0.125`).

**The only reason a reduction is ever taken here is `omega`**, and the run
names it where it takes one:
`-> neither the joint distance nor the composition distance fell;
under-relaxation omega = 0.250` (`EXHALE_main.f90:7386-7389`). There is no
movement bound in this route: the bound of phase 6 is `carrier_trust`, which
belongs to the carrier relaxation these cases do not run, so the `trust` field
of the pass line is the unused setting and never acted. No trial was refused
in any pass of any of the three, and no closure refusal was printed.

**Two state sequences, kept apart.** `Y(k)` is the state at the joint test of
pass `k`, which the certificate of pass `k` describes. `Z(k)` is the state at
the END of pass `k`, after the element relaxation, which
`write_diffusion_pass_profile` exports (`EXHALE_main.f90:8258`, under
`EXHALE_DIFFUSION_CHECK=1`) and which the next solve is entered with. The
compact norms of section 7 are norms of the `Z` sequence, because that is the
sequence the export carries.

**Where each recorded quantity is written.** All INSPECTED in
`EXHALE_main.f90` unless stated:

| quantity | the line that writes it |
|---|---|
| joint residual, worst gated row and cell, the three hydrodynamic rows, `omega`, `trust` | the pass line, `:7805-7812` |
| element relaxation ending, inner map distance, displacement applied, steps | `:7813-7819` |
| element residual of the returned composition, with its absolute value, its row terms and its cell | `:7849-7864` |
| accepted composition displacement | `:7866-7870` |
| coupled residual of the whole state, and which measure decided progress | `:7871-7875` |
| the reduction of `omega`, and its reason | `:7386-7389` |
| the rows of the state each solve is ENTERED with | `(JFNK) start ||R||` and `(JFNK) below r_esc [..] |R|: mass= .. mom/grav= .. energy=`, `steady_newton.f90` |
| the solve's ending, its merit `||Fs||2` per iteration, its stagnation statement | `(JFNK) it ..`, `(JFNK) STAGNATED ..`, `(JFNK) done info=..`, `steady_newton.f90:18786` onward |
| the mass-flux spread of the state each solve hands back | `(JFNK) done .. flux spread=` and `(JFNK) flux spread by window` |

**A window difference that has to be said once.** The element residual each
pass reports is a COLUMN MAXIMUM over all 500 cells. The certification gates
the same row only in the wind, at `r >= 1.20`, and reports the cells below that
radius without gating. The two numbers are the same row over two windows and
are never compared with each other below.

---

## 5. The gate of rev9 section 4.2, for every diagnostic file used

| requirement | how it is met |
|---|---|
| the complete environment actually in force | section 3.1, and `ENV.txt` beside every run |
| the file, its existence, nonzero size, header and row count | `output/diffusion_pass_profiles.txt`: present in all three runs, 2 header lines, and 500 rows in each of the blocks its first column labels (12, 12 and 5 blocks; MEASURED). `output/element_flux_profile.txt`: present in all three |
| which state the rows belong to | `write_diffusion_pass_profile` appends one block per pass and labels each block with its pass number in column 1, so the `Z` sequence is identified by the file itself |
| the exit class | all three exited 0; two ended on their own pass budget, one on the alternation's own no-progress rule, and none reached the wall ceiling |

---

## 6. The three trajectories, pass by pass

Each table is **Mode C**: a retained trajectory advanced under `EXHALE.x` md5
`89ec67149aa452eea3deeb19052f83d1` from the generation named in section 3.2,
under the current operator on a state written by an older boundary model.
Every value is MEASURED, read from the line of `run.log` named in section 4.

Conventions:

- **solve entry rows** are the mass, momentum and energy rows of the state the
  solve of that pass is ENTERED with, which is the state the composition
  update of the previous pass returned.
- **ending `||R||`** is the solve's own `done info=.. ||R||=..`, printed to
  four digits; the pass line prints the same rows to three.
- **joint residual** is the largest certification entry over its own tolerance,
  evaluated on the refreshed state of that pass. It is the acceptance
  criterion.
- **gated elemental row** is the certification's elemental transport He/H
  partition, gated at `r >= 1.20`.
- **element residual returned** is a COLUMN MAXIMUM over all 500 cells of the
  element operator's own row at the composition the pass hands back.
- the last pass of each run takes no composition update and says so.

### 6.1 Case 1, `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7`

| pass | solve entry rows (mass, momentum, energy) | info | ending \|\|R\|\| | rows after the solve (mass, momentum, energy) | gated elemental row, cell | joint residual | element ending | REQUESTED | ADOPTED | omega | element residual returned, cell | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 9.806e-01, 8.724e-06, 1.071e+00 | 0 | 2.559e-07 | 2.44e-08, 1.02e-14, 2.56e-07 | 7.220e-05, 295 | 7.22e+00 | the fixed point of the element operator | 2.620e-05 | 1.310e-05 | 0.500 | 1.790e-03 at 4 | 51.9 |
| 2 | 6.196e-08, 1.017e-14, 6.512e-06 | 2 | 6.475e-06 | 3.82e-08, 1.02e-14, 6.48e-06 | 4.440e-05, 284 | 6.48e+00 | the same | 1.790e-05 | 8.970e-06 | 0.500 | 8.930e-04 at 4 | 90.6 |
| 3 | 2.439e-08, 1.017e-14, 7.942e-06 | 2 | 7.919e-06 | 1.58e-06, 1.13e-13, 7.92e-06 | 2.219e-05, 288 | 7.92e+00 | the same | 1.950e-05 | 9.770e-06 | 0.500 | 4.460e-04 at 4 | 94.7 |
| 4 | 1.589e-06, 1.143e-13, 1.111e-05 | 2 | 1.110e-05 | 1.60e-06, 1.16e-13, 1.11e-05 | 1.038e-05, 263 | 1.11e+01 | the same | 1.920e-05 | 9.610e-06 | 0.500 | 2.230e-04 at 4 | 93.4 |
| 5 | 1.597e-06, 1.147e-13, 1.463e-05 | 2 | 1.460e-05 | 1.52e-06, 5.67e-12, 1.46e-05 | 1.847e-05, 269 | 1.46e+01 | the same | 1.700e-05 | 8.490e-06 | 0.500 | 1.120e-04 at 4 | 95.0 |
| 6 | 1.496e-06, 5.674e-12, 1.809e-05 | 2 | 1.809e-05 | 1.46e-06, 5.68e-12, 1.81e-05 | 2.246e-05, 272 | 1.81e+01 | the same | 1.400e-05 | 6.980e-06 | 0.500 | 5.580e-05 at 4 | 125.3 |
| 7 | 1.438e-06, 5.677e-12, 2.123e-05 | 2 | 2.123e-05 | 1.44e-06, 5.68e-12, 2.12e-05 | 2.273e-05, 273 | 2.12e+01 | the same | 1.090e-05 | 5.450e-06 | 0.500 | 2.790e-05 at 4 | 125.1 |
| 8 | 1.445e-06, 5.680e-12, 2.382e-05 | 2 | 2.382e-05 | 1.44e-06, 5.68e-12, 2.38e-05 | 2.060e-05, 274 | 2.38e+01 | the same | 8.200e-06 | 4.100e-06 | 0.500 | 1.750e-05 at 275 | 125.9 |
| 9 | 1.442e-06, 5.680e-12, 2.582e-05 | 2 | 2.582e-05 | 1.91e-06, 5.70e-12, 2.58e-05 | 1.750e-05, 275 | 2.58e+01 | the same | 5.990e-06 | 2.990e-06 | 0.500 | 1.410e-05 at 275 | 125.3 |
| 10 | 1.910e-06, 5.696e-12, 2.731e-05 | 2 | 2.731e-05 | 2.14e-06, 5.70e-12, 2.73e-05 | 1.410e-05, 275 | 2.73e+01 | the same | 4.270e-06 | 2.140e-06 | 0.500 | 1.090e-05 at 276 | 125.3 |
| 11 | 2.140e-06, 5.704e-12, 2.837e-05 | 2 | 2.837e-05 | 2.30e-06, 5.71e-12, 2.84e-05 | 1.090e-05, 276 | 2.84e+01 | the same | 2.990e-06 | 1.500e-06 | 0.500 | 8.200e-06 at 276 | 104.5 |
| 12 | 2.298e-06, 5.706e-12, 2.912e-05 | 2 | 2.911e-05 | 2.53e-06, 5.72e-12, 2.91e-05 | **8.190e-06**, 276 | 2.91e+01 | NO UPDATE: the last pass, no solve is left to consume one | -- | -- | 0.500 | -- | 123.2 |

Radii of the cells that carry the gated row: 263 is `r = 1.4491`, 269
`1.4987`, 272 `1.5256`, 273 `1.5349`, 274 `1.5443`, 275 `1.5539`, 276
`1.5636`, 284 `1.6482`, 288 `1.6951`, 295 `1.7854` (MEASURED from the state's
own radius column). The refusing cell walks 33 cells INWARD over the
trajectory, from 295 to 263, then 13 cells outward and stops at 276.

**The two entries that refuse at the ending** (READ, the last certification):
`hydrodynamic mass row: row measure 1.067E-09 above 3.9E-10 at cell 236`
(`r = 1.2800`) and
`hydrodynamic energy row: row measure 2.911E-05 above 1.0E-06 at cell 260`
(`r = 1.4261`). The elemental row does not refuse.

The mass-flux spread of the state each solve hands back, which is the gate the
wind's stationarity is judged on: 1.146E-11, 1.146E-11, 1.513E-10, 1.534E-10,
1.462E-08, 1.463E-08, 1.463E-08, 1.464E-08, 1.472E-08, 1.476E-08, 1.478E-08,
1.482E-08. It rises three decades over the twelve passes and stays nine
decades below the 6.5E-03 of case 3.

### 6.2 Case 2, `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`

| pass | solve entry rows (mass, momentum, energy) | info | ending \|\|R\|\| | rows after the solve (mass, momentum, energy) | gated elemental row, cell | joint residual | REQUESTED | ADOPTED | omega | element residual returned, cell | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 9.589e-01, 5.539e-06, 1.112e+00 | 0 | 8.946e-07 | 4.06e-08, 2.92e-14, 8.95e-07 | 1.200e-04, 248 | 1.20e+01 | 8.310e-05 | 4.160e-05 | 0.500 | 4.140e-03 at 4 | 50.3 |
| 2 | 6.951e-08, 2.916e-14, 1.744e-05 | 2 | 1.745e-05 | 6.95e-08, 2.92e-14, 1.75e-05 | 5.910e-05, 254 | 1.75e+01 | 5.660e-05 | 2.830e-05 | 0.500 | 2.050e-03 at 4 | 96.9 |
| 3 | 1.117e-07, 2.916e-14, 3.170e-05 | 2 | 3.081e-05 | 5.45e-06, 1.07e-11, 3.08e-05 | 2.650e-05, 263 | 3.08e+01 | 4.230e-05 | 2.110e-05 | 0.500 | 1.020e-03 at 4 | 58.8 |
| 4 | 5.452e-06, 1.073e-11, 4.153e-05 | 2 | 4.153e-05 | 5.45e-06, 1.07e-11, 4.15e-05 | 3.070e-05, 223 | 4.15e+01 | 3.150e-05 | 1.580e-05 | 0.500 | 5.110e-04 at 4 | 40.0 |
| 5 | 5.452e-06, 1.073e-11, 4.936e-05 | 2 | 4.929e-05 | 8.63e-06, 1.13e-11, 4.93e-05 | 3.290e-05, 227 | 4.93e+01 | 2.290e-05 | 1.140e-05 | 0.500 | 2.550e-04 at 4 | 40.8 |
| 6 | 8.691e-06, 1.135e-11, 5.480e-05 | 2 | 5.425e-05 | 9.63e-06, 1.54e-11, 5.43e-05 | 3.000e-05, 229 | 5.43e+01 | 1.610e-05 | 8.070e-06 | 0.500 | 1.270e-04 at 4 | 43.8 |
| 7 | 9.522e-06, 1.542e-11, 5.796e-05 | 2 | 5.795e-05 | 8.21e-06, 1.54e-11, 5.80e-05 | 2.500e-05, 231 | 5.80e+01 | 1.110e-05 | 5.560e-06 | 0.500 | 6.380e-05 at 4 | 41.1 |
| 8 | 8.235e-06, 1.541e-11, 6.044e-05 | 2 | 6.016e-05 | 1.01e-05, 1.73e-11, 6.02e-05 | 1.960e-05, 232 | 6.02e+01 | 7.520e-06 | 3.760e-06 | 0.500 | 3.180e-05 at 4 | 41.6 |
| 9 | 1.007e-05, 1.730e-11, 6.177e-05 | 2 | 6.177e-05 | 1.01e-05, 1.73e-11, 6.18e-05 | 1.470e-05, 233 | 6.18e+01 | 5.000e-06 | 2.500e-06 | 0.500 | 1.590e-05 at 4 | 39.9 |
| 10 | 1.007e-05, 1.730e-11, 6.281e-05 | 2 | 6.281e-05 | 1.01e-05, 1.73e-11, 6.28e-05 | 1.060e-05, 234 | 6.28e+01 | 3.280e-06 | 1.640e-06 | 0.500 | 7.950e-06 at 4 | 39.9 |
| 11 | 1.007e-05, 1.730e-11, 6.346e-05 | 2 | 6.346e-05 | 1.01e-05, 1.73e-11, 6.35e-05 | 7.520e-06, 234 | 6.35e+01 | 2.120e-06 | 1.060e-06 | 0.500 | 5.210e-06 at 235 | 39.8 |
| 12 | 1.007e-05, 1.730e-11, 6.387e-05 | 2 | 6.387e-05 | 1.01e-05, 1.73e-11, 6.39e-05 | **5.210e-06**, 235 | 6.39e+01 | -- | -- | 0.500 | -- | 39.7 |

The element relaxation ended on the fixed point of the element operator at
every pass. Radii: 223 is `r = 1.2230`, 227 `1.2392`, 229 `1.2477`, 231
`1.2565`, 232 `1.2611`, 233 `1.2657`, 234 `1.2704`, 235 `1.2751`, 248
`1.3455`, 254 `1.3837`, 263 `1.4491`. The refusing cell walks 25 cells inward
and then 12 outward.

**The two entries that refuse at the ending**:
`hydrodynamic mass row: row measure 2.811E-08 above 1.6E-09 at cell 194`
(`r = 1.1340`) and
`hydrodynamic energy row: row measure 6.387E-05 above 1.0E-06 at cell 157`
(`r = 1.0697`). Again the elemental row does not refuse.

Mass-flux spread per solve: 1.172E-09, 1.172E-09, 7.562E-08, 7.562E-08,
8.075E-08, 1.173E-07, 1.199E-07, 1.377E-07, and 1.377E-07 for the last four.
Two decades over the twelve passes, then flat.

### 6.3 Case 3, `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`

| pass | solve entry rows (mass, momentum, energy) | info | ending \|\|R\|\| | rows after the solve (mass, momentum, energy) | gated elemental row, cell | joint residual | REQUESTED | ADOPTED | omega | element residual returned | wall [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6.925e-01, 2.926e-07, 1.949e+00 | 2 | 1.949e+00 | 6.92e-01, 7.77e-07, 1.95e+00 | 3.540e-02, 389 | 5.07e+06 | 5.260e-03 | 2.630e-03 | 0.500 | 3.540e-02 at 389 | 330.4 |
| 2 | 6.925e-01, 2.926e-07, 1.949e+00 | 2 | 1.949e+00 | the same | 3.540e-02, 389 | 5.07e+06 | 3.350e-03 | 1.680e-03 | 0.500 | 3.540e-02 at 389 | 348.1 |
| 3 | 6.925e-01, 2.926e-07, 1.949e+00 | 2 | 1.949e+00 | the same | 3.540e-02, 389 | 5.07e+06 | 2.120e-03 | 5.310e-04 | **0.250** | 3.540e-02 at 389 | 348.1 |
| 4 | 6.925e-01, 2.926e-07, 1.949e+00 | 2 | 1.949e+00 | the same | 3.540e-02, 389 | 5.07e+06 | 1.730e-03 | 2.170e-04 | **0.125** | 3.540e-02 at 389 | 347.9 |
| 5 | 6.925e-01, 2.926e-07, 1.949e+00 | 2 | 1.949e+00 | the same | 3.540e-02, 389 | 5.07e+06 | NO UPDATE: this pass ends the iteration | -- | 0.125 | -- | 328.3 |

`omega` was halved at pass 3 and again at pass 4, each time with the run's own
reason, `neither the joint distance nor the composition distance fell`, and
the loop then refused at pass 5. **Every solve is bitwise the same
experiment**: the same entry `||R|| = 1.949E+00`, the same 19 non-monotone
accepts, the same `flux spread = 6.539E-03`, the same ending. The composition
does move and the movement does shrink, from 2.630E-03 to 2.170E-04 accepted,
and it changes no row of the certification at the digits printed.

**The gated elemental row of this case is a ratio of two vanishingly small
numbers.** At every pass the run prints
`element residual of the returned composition 3.54E-02 = |res| 1.25E-27 over
its own row terms 3.55E-26 at cell 389` (pass 1; `6.27E-28` over `1.77E-26` at
pass 2). Cell 389 is `r = 5.0484` (MEASURED from the state's own radius
column). This is the collapsing-denominator mechanism phase 6 identified in
the outer cells (`docs/lhs1140b_alternation_audit_20260921.md` section 8.4),
here at 5 planetary radii rather than at 28: the row's own physical scale has
fallen 26 decades, so a balance that is absolutely negligible leaves a large
relative measure.

**What refuses at the ending**: four entries, the three hydrodynamic rows and
the elemental row, with `elemental transport He/H partition: gated row measure
3.538E-02 above 1.0E-05 at cell 389 (a wind cell)`. The flux gate reads
`accepted 6.5390E-03 as written 6.5390E-03`, that is, the wind of this stored
state has a mass-flux spread of 6.5E-03 and the current operator does not
reduce it.

---

## 7. The compact state norms

`Z(k)` is the state at the END of pass `k`, read from
`output/diffusion_pass_profiles.txt`. The seven exported blocks are the
He/H ratio, `r^2 rho v`, `T`, `x_H2`, `n(H3+)`, `n(He II)` and `n_e`; in these
atomic runs `x_H2` and `n(H3+)` are identically zero. **The export does not
carry the density and the velocity separately, only their product, and carries
no pressure column**, so these are norms of those blocks and not of the full
conserved vector; that is the whole of what the current code exports per pass.
The block scales are FIXED per case at the entry state's own column maximum,
and the compact norm is the largest scaled block difference.

| case | one-pass norms `\|\|Z(k+1)-Z(k)\|\|` | ratio of consecutive one-pass norms | two-pass `r` | the block that carries it |
|---|---|---|---|---|
| 1 | 1.694E-04, 8.724E-05, 4.471E-05, 2.277E-05, 1.149E-05, 5.724E-06, 3.802E-06, 2.616E-06, 1.788E-06, 1.208E-06, 9.878E-07 | 0.515, 0.512, 0.509, 0.505, 0.498, 0.664, 0.688, 0.684, 0.675, 0.818 | 0.81 to 1.00 | `n_HeII` at every step but the last, where `n_e` carries it |
| 2 | 1.043E-04, 5.450E-05, 3.125E-05, 2.105E-05, 1.437E-05, 9.710E-06, 6.474E-06, 4.260E-06, 2.780E-06, 1.784E-06, 1.070E-06 | 0.522, 0.573, 0.674, 0.683, 0.676, 0.667, 0.658, 0.653, 0.642, 0.600 | 0.69 to 1.00 | the same |
| 3 | 1.385E-03, 5.511E-04, 1.968E-04, 8.831E-05 | 0.398, 0.357, 0.449 | 0.87 to 0.98 | the same |

`r` is the two-pass difference over the sum of the two one-pass differences it
spans; aligned steps give 1 and a pair of states that undo each other gives 0.
**A two-cycle is excluded in all three cases**: the smallest `r` measured
anywhere is 0.69.

**The `mass_flux` block never carries the norm in any of the three.** Its
largest single-pass movement is 5.98E-09 of its own column maximum in case 1,
1.94E-08 in case 2, and exactly zero at every pass of case 3. Phase 6's `wellmixed` failure, in
which ONE bounded composition update moved the mass flux by 0.985 of its own
column scale and the wind was never recovered, **does not happen here**: all
three carry element diffusion, and none of them loses its wind to a
composition update. Case 3's wind was never in solution to begin with.

---

## 8. Answers to the four questions the brief asks

### 8.1 Does the gated elemental row fall geometrically, stall, or rise?

**Case 1: it falls, turns, and falls again, and ends BELOW its tolerance.**
The measure over the twelve passes is 7.220E-05, 4.440E-05, 2.219E-05,
1.038E-05, 1.847E-05, 2.246E-05, 2.273E-05, 2.060E-05, 1.750E-05, 1.410E-05,
1.090E-05, 8.190E-06, with pass-to-pass ratios (previous over current)
1.63, 2.00, 2.13, 0.56, 0.82, 0.99, 1.10, 1.18, 1.24, 1.29, 1.33. The first
three passes are a clean halving; passes 5 to 7 reverse it; the last five fall
at a ratio that is itself rising toward 4/3.

**Case 2: the same shape with a shallower excursion, and it also ends below
tolerance.** 1.200E-04, 5.910E-05, 2.650E-05, 3.070E-05, 3.290E-05,
3.000E-05, 2.500E-05, 1.960E-05, 1.470E-05, 1.060E-05, 7.520E-06, 5.210E-06;
ratios 2.03, 2.23, 0.86, 0.93, 1.10, 1.20, 1.28, 1.33, 1.39, 1.41, 1.44.

**Case 3: it stalls exactly.** 3.540E-02 at cell 389 at every one of the five
passes, at the printed resolution.

**The element operator is not what limits any of them.** In all three the
element relaxation ended on the fixed point of its own operator at every pass,
with no refused trial and no closure refusal, and in cases 1 and 2 the element
residual's COLUMN MAXIMUM halves exactly every pass while `omega = 0.5`
(case 1: 1.790E-03, 8.930E-04, 4.460E-04, 2.230E-04, 1.120E-04, 5.580E-05,
2.790E-05, then it leaves the base and stands at the wind cell). So the
turn of the GATED row at pass 5 of case 1 is not the element relaxation
failing to converge; it is the wind moving under it, and the worst cell moves
with it.

### 8.2 Does the energy row oscillate, lifted by each update and pulled back by each solve?

**There is no oscillation. There is a ratchet.** The composition update lifts
the energy row at every pass, and the state the hydrodynamic solve hands back
carries at most 1.6 per cent of that lift away in case 1 and at most 11.2 per
cent in case 2.

Case 1, in units of 1E-06, `lift` being the entry row of pass `k+1` minus the
ending `||R||` of pass `k`, and `recovered` the entry row of pass `k` minus its
own ending `||R||`:

| pass | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lifted by the update | +6.256 | +1.467 | +3.191 | +3.530 | +3.490 | +3.140 | +2.590 | +2.000 | +1.490 | +1.060 | +0.750 |
| recovered by the solve | -0.037 | -0.023 | -0.010 | -0.030 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | -0.010 |

Case 2, the same, in units of 1E-06:

| pass | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lifted by the update | +16.55 | +14.25 | +10.72 | +7.83 | +5.51 | +3.71 | +2.49 | +1.61 | +1.04 | +0.65 | +0.41 |
| recovered by the solve | +0.01 | -0.89 | 0.000 | -0.07 | -0.55 | -0.01 | -0.28 | 0.000 | 0.000 | 0.000 | 0.000 |

Pass 1 is the exception and it is the whole of the solve's work in both cases:
it takes the energy row from 1.071E+00 to 2.559E-07 in case 1 and from
1.112E+00 to 8.946E-07 in case 2, and certifies every hydrodynamic row.

**The lifts decay geometrically.** Their ratios settle at 0.71 in case 1
(0.90, 0.82, 0.77, 0.74, 0.71, 0.71 over the last six) and at 0.65 in case 2
(0.67, 0.67, 0.65, 0.65, 0.62, 0.63). Summing the remaining geometric tail
puts the limit of the energy row at about 3.1E-05 in case 1 and about 6.5E-05
in case 2, that is, 31 and 65 times its tolerance. **That extrapolation is an
extrapolation and is labeled one**; what is MEASURED is that the row rose
monotonically to 2.911E-05 and 6.387E-05 and that its increments fell by a
factor 8 and 40 over the run.

In case 3 there is neither a lift nor a recovery: the entry and the ending
rows are the same numbers at every pass.

### 8.3 Does any pass certify?

**No pass of any of the three certifies the state.** Every pass of all three
ends with `NOT CERTIFIED`, and the number of refusing entries is 2, 2 and 4 at
the endings.

What DOES certify, and is worth naming because it is what phase 9 measured:
**pass 1 of cases 1 and 2 certifies every HYDRODYNAMIC row** (`info = 0`, all
three rows within tolerance), and from pass 11 or 12 the gated elemental row
is within tolerance as well. The two never happen on the same state. The
certificate of the last pass of case 1 reads

```
   NOT CERTIFIED: 2 entry/entries of the inventory refuse it
     hydrodynamic mass row: row measure  1.067E-09 above  3.9E-10 at cell 236
     hydrodynamic energy row: row measure  2.911E-05 above  1.0E-06 at cell 260
```

and of case 2

```
   NOT CERTIFIED: 2 entry/entries of the inventory refuse it
     hydrodynamic mass row: row measure  2.811E-08 above  1.6E-09 at cell 194
     hydrodynamic energy row: row measure  6.387E-05 above  1.0E-06 at cell 157
```

No state was published.

### 8.4 Is the stagnation the arithmetic floor every time?

**On the default route the quantity the question names cannot be formed, and
that is said rather than substituted.** Phase 9's "predicted and observed
reductions of 1E-21 to 1E-24 against a merit of 3.58E-09" are `pred` and
`actual` of the TRUST-REGION variant
(`docs/lhs1140b_p2_base_stage_cd_20260921.md` section 5.4, the table of
`ledger_tr_steps.txt`), not of the observation group. `EXHALE_TRUST_REGION`
changes the step control AND the merit, so running it would not have continued
these trajectories. It was not run.

What the default route does print is the merit `||Fs||2` at every iteration,
the `||R||` of every iterate, and the judged distance of the best iterate the
solve keeps. **Those three separate the three cases, and on one of them they
say the opposite of "the arithmetic floor": the solve is not out of arithmetic,
it finds much better iterates and its own acceptance functional throws them
away.** MEASURED, over every solve after the first, as the ratio of the entry
value to the smallest value the solve reached:

| case | `\|\|R\|\|` entry over the smallest reached | merit entry over the smallest reached | `\|\|R\|\|` of the state handed back |
|---|---|---|---|
| 1, `kzz1e9/HeH9.7` | **3.73 to 4.18 at every one of passes 2 to 12** | 1.92 to 2.09 | the entry value, to four digits, at every pass |
| 2, `kzz1e9/HeH2.13` | 1.32 at pass 2, then **1.00 to 1.03 at passes 3 to 12** | 1.89 at pass 2, then 1.04 to 1.06 | the entry value at every pass |
| 3, `kzzprofile/HeH9.7` | **1.42 at all five passes**, identically | **4.43 at all five**, identically | the entry value at every pass |

So the answer is case by case:

- **case 2 is the one the phrase fits.** From pass 3 its solve cannot improve
  `||R||` at all, to two decimal places, and moves the merit by 6 per cent in
  39 to 43 iterations. There is nothing left for it to do.
- **case 1 is the opposite.** Every solve from pass 2 reaches an iterate whose
  `||R||` is about four times smaller than the one it was entered with, halves
  its merit, and hands back the entry state, because the measure it accepts on
  is the JUDGED distance, the largest certified row over its own tolerance, and
  that never improves. `best judged iterate: distance 2.911E+01` at pass 12 is
  exactly the entry energy row 2.912E-05 over its 1.0E-06, and the run says so:
  `stopping at a judged distance 9.290E+03, best 2.911E+01, merit over the
  window 1.331E-08`, with `outer iterations that met the acceptance gate while
  a certified row refused: 44` and `accepted steps that raised the judged rows
  above their five-iterate window: 5`.
- **case 3 is the same shape as case 1 and more extreme**: the merit falls by
  a factor 4.43 and `||R||` by 1.42 in every one of its five solves, to every
  printed digit, and the judged distance stands at 5.075E+06 at the entry and
  at the end of each.

Every solve after the first ends on the SAME statement, `STAGNATED: neither
the judged distance nor the merit improved in 40 iterations, monotone search
included`, which is a self-declared stop and not a cap (INSPECTED,
`steady_newton.f90:18786`); case 1's pass 2 is the proof, capped at 60 by
phase 9 and uncapped here, ending at iteration 59 with the same verdict.

**What this changes in rev9 section 17.2.** Reading the second pass of case 1
as "the hydrodynamic solve has nothing left to do" is not supported over
twelve passes: on that case the solve has a factor four in `||R||` left and
takes it, and what refuses the result is the acceptance functional, not the
arithmetic. The `pred`/`actual` pair phase 9 quoted belongs to the
trust-region variant, which changes the step control AND the merit; it is not
formed here and is not substituted for.

One finding that phase 9 could not make: **the stagnation is not the 60
iteration cap.** Pass 2 of case 1 was capped at 60 there and ran uncapped
here, and it ended at iteration 59 with the same `info = 2`, the same
`||R|| = 6.475E-06`, the same gated row 4.44E-05 at cell 284 and the same
energy row 6.48E-06.

---

## 9. A defect found in the progress control, and what it did to these runs

**The test that asks whether a hydrodynamic row carries the state's distance
can never be true.** INSPECTED, `EXHALE_main.f90:7235`:

```fortran
            if (cert_now%e(icert)%name(1:14) .eq. 'hydrodynamic ')        &
               hydro_worst_dist = max(hydro_worst_dist, dist_here)
```

The literal `'hydrodynamic '` is 13 characters. Fortran pads the shorter
operand of a character comparison with blanks, so the right-hand side is
compared as `'hydrodynamic  '` with two trailing blanks. The names the
certification gives these rows are `'hydrodynamic mass row'`,
`'hydrodynamic momentum row'` and `'hydrodynamic energy row'` (INSPECTED,
`certification.f90:1299-1301`), whose first fourteen characters are
`'hydrodynamic m'`, `'hydrodynamic m'` and `'hydrodynamic e'`. None equals
`'hydrodynamic  '`, so the branch never fires and `hydro_worst_dist` is
identically zero at every pass of every run. Its one consumer is the progress
control (`:7358`):

```fortran
            composition_closing = (hydro_worst_dist .lt. prog_worst) .and. &
                                  ((n_comp_updates .lt. 2) .or.           &
                                   (comp_residual_last .lt. comp_residual_before))
```

so the first factor is always true and the control cannot tell a state whose
distance is carried by a hydrodynamic row from one whose distance is carried
by the composition.

**MEASURED proof from this experiment's own logs, not from reading alone.** On
every pass from 3 to 12 of case 1 and from 3 to 12 of case 2 the run prints

```
    coupled residual of the whole state  2.12E+01 ...; progress decided by:
    the composition residual fell while the state's distance is carried by a
    row that is not hydrodynamic
```

while the certification of that same pass reads
`hydrodynamic energy row ... max= 2.123E-05 ... tol= 1.0E-06 ABOVE` and the
worst gated species row stands at only 2.27 times its own tolerance. The
largest certification entry, 21.2, IS the hydrodynamic energy row, and the
sentence says it is not.

**What it did to this measurement.** With the test working,
`hydro_worst_dist` would equal `prog_worst` on those passes, `composition_closing`
would be false, and `n_no_fall` would have counted every pass on which the
joint residual rose. The joint residual of case 1 rises at every pass from 3
to 12 (7.92, 11.1, 14.6, 18.1, 21.2, 23.8, 25.8, 27.3, 28.4, 29.1) and of case
2 from 2 to 12 (17.5, 30.8, 41.5, 49.3, 54.3, 58.0, 60.2, 61.8, 62.8, 63.5,
63.9). With `outer_no_fall_max = 3` (INSPECTED, `EXHALE_main.f90:6885`) case 1
would have reached its third consecutive pass without a fall at pass 5 and case
2 at pass 4, and each would have refused there on the same rule case 3 ended
on. The trajectory would not have been otherwise unchanged, because `omega`
would have been halved on each of those passes, so this says which rule would
have fired and when, not what the states would then have been. Case 3 shows
the rule CAN fire: there the composition residual did not fall either, so the
second factor of `composition_closing` was false on its own and the loop
refused at pass 5.

**The repair is one digit**, `name(1:13)` in place of `name(1:14)`, or
better `index(name, 'hydrodynamic ') .eq. 1`. **It was not applied here** and
the reason is recorded: `src/EXHALE_main.f90` is modified in the working tree
by another worker, the brief names no source file for this item, and every
number in this document was MEASURED with the binary of record, which a repair
would invalidate. This is reported for the owner of that file to fix.

**Nothing in section 8 depends on the defect.** The rows, the ratios, the
lifts and the recoveries are all measurements of states; the defect changes
only how many passes the loop would have taken before stopping, and a longer
trajectory is what this phase wanted.

---

## 10. What this decides about the three P2 cases, and what it does not

### 10.1 What it decides

1. **The question of rev9 section 17.3 item 2 is answered for the two scalar
   cases. The elemental row keeps falling; it goes under its tolerance; and
   the state still does not certify.** After twelve passes the gated elemental
   transport He/H row stands at 0.82 and 0.52 of its own tolerance, and the
   refusing entries are the two hydrodynamic rows. The obstruction named in
   rev9 section 17.2, "the elemental transport He/H partition row in the wind",
   is therefore NOT the standing obstruction of these two cases. It is a
   transient of the first passes.
2. **The standing obstruction of the two scalar cases is the alternation's own
   fixed point.** The composition update lifts the energy row by an amount that
   decays geometrically, the hydrodynamic solve after pass 1 returns the state
   it was given, and the sequence converges. The limit is a state whose energy
   row is tens of times its tolerance. The accepted state displacement halves
   every pass and a two-cycle is excluded, so this is a convergence and not a
   wandering.
3. **The hydrodynamic solve hands back its entry state after pass 1 on all
   three cases, and on two of them not for want of arithmetic.** Every solve
   from pass 2 stagnates by its own criterion and not by a cap, and phase 9's
   pass-2 observation is reproduced exactly and generalizes to all eleven later
   passes. But on case 1 every such solve reaches an iterate whose `||R||` is
   about four times smaller than its entry, and on case 3 one whose merit is
   4.43 times smaller, and discards it because the judged distance, the
   certified rows over their tolerances, never improves. Only on case 2 is the
   solve genuinely out of work. Section 8.4 carries the numbers, and rev9's
   reading of the second pass as "the arithmetic floor" does not survive them.
4. **The profile case is a different problem and should not be grouped with
   the two scalar ones.** Its stored state has a mass-flux spread of 6.5E-03,
   its first solve does not certify anything, no solve moves any row, and the
   run refuses itself after five passes. Its gated elemental row is a ratio of
   1.25E-27 to 3.55E-26 at `r = 5.05`, which is the collapsing-denominator
   mechanism and not a physical imbalance.
5. **Element diffusion does not reproduce phase 6's `wellmixed` failure.** In
   no pass of any of the three does the composition update take the wind out of
   solution: the `mass_flux` block moves by at most 2E-08 of its own column
   scale in the two scalar cases and by exactly zero in the profile one,
   against the 0.985 of phase 6's cases 2 and 4, and no run printed the
   write-back warning.
6. A defect in the progress control, section 9, which made the no-progress rule
   unable to fire on the two scalar cases.

### 10.2 What it does not decide

- **Why the acceptance functional refuses the better iterates of case 1.** It
  is measured that the solve reaches `||R||` four times smaller and that the
  judged distance does not improve, so the certified rows and `||R||` move
  apart. Which row does that, and whether it is the mass row's cell-dependent
  tolerance or the energy row's own scale, is a question for the Stage C and D
  instruments on the state of a LATE pass, and those instruments were not run
  here.
- **Whether the limit state of the two scalar cases is a fixed point of the
  alternation or only looks like one over twelve passes.** The evidence for it
  is the geometric decay of the lifts and the halving of the state
  displacement; the extrapolated limits (3.1E-05 and 6.5E-05) are
  extrapolations.
- **Whether a smaller `omega` would give a different limit.** `omega` stayed at
  0.500 for all twelve passes of both scalar cases, because the progress
  control could not fire; that is the defect of section 9 and a run with the
  repair is the measurement that would answer it.
- **Whether the profile case's refusal is the profile handoff, the metals it
  carries, its tighter `Resid tol` or its stored state.** Its entry state
  already refuses by six decades on the mass row and this phase changed
  nothing about that.
- **Anything about the molecular cases**, which are not in this phase.

---

## 11. Artifacts

All raw artifacts are in the scratch directory of this phase,

```
/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/p2_atomic_alternation/
```

with `artifacts/README_artifacts.md` describing the layout. Per case,
`LHS1140b/models/<case>/<HeH>/` holds `input.inp` (the case's own, unedited),
`output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt` (the stored
generation), `ENV.txt`, `run.log`, `run.err`, `TIMING.txt`,
`EXHALE_setup.out`, `EXHALE_resolved.out` and `output/`. `SNAPSHOTS.txt`
carries the generation, byte count and md5 of every state file as set up
beside the md5 its manifest records; `modec_exit.txt` the exit status of each
run; `mkrun.sh`, `runcase.sh` and `runall.sh` the commands that made and ran
them. Under `artifacts/<case tag>/`: `pass_table.md`, `gated_by_pass.txt`,
`solve_entry_exit.txt`, `solve_work.txt`, `best_judged.txt`,
`energy_alternation.md`, `state_norms.txt`, `two_cycle.txt`,
`flux_spread.txt`, `base_face.txt`, `trial_ledger.txt` (empty, see section
2.3) and `pass_lines.txt`. The readers are in `analysis/`;
`pass_record.py`, `gated_by_pass.py`, `trial_ledger.py` and `state_norms.py`
are the phase 6 readers of `docs/lhs1140b_alternation_audit_20260921.md`,
copied so that the two documents read the same lines the same way, with
`gated_by_pass.py` widened to the momentum row, the He 2^3S level row and the
eliminated-species closure; `solve_entry_exit.py`, `atomic_pass_table.py`,
`energy_row_alternation.py`, `flux_spread.py`, `base_face_by_pass.py`,
`solve_work.py` and `two_cycle_ratio.py` are new here.
