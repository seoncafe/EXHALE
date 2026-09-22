# LHS 1140 b: Stage A repeatability of the stationary residual

2026-09-21. Phase 3 of `PLAN_20260920_rev9.md` section 15, the measurement
section 5.2 asks for, on two of the checkpoints whose Stage A0 identity
records are `docs/lhs1140b_identity_records_20260921.md` sections 3.1 and 3.4.

## 1. The result, first

**Every evaluation of the stationary residual performed here was bitwise
reproducible.** At five thread configurations, in two independent processes
each, on both checkpoints, the instrument reported `0 of N entries differ`
for both of its comparisons, and printed the same verdict line every time.
The full standard output of the diagnostic is byte-identical across
processes and across thread counts except for the two lines that state the
thread configuration itself. The serialized residual vector of the loaded
state is byte-identical across processes and across thread counts.

There is therefore no measured evaluation noise to subtract at either
checkpoint, and Stage B can be read as measuring the operator. The gate of
section 15 phase 3, "a reproducible operator or a quantified uncertainty",
is met by the first of the two.

Three qualifications belong with that result and are stated in full below:

1. The vector the instrument replays is 1500 entries, three unknowns per
   cell, for BOTH checkpoints. Neither case registers a species row, because
   both carry `carrier_newton=F` and the registry is gated on that single
   flag. The molecular checkpoint does NOT give a longer vector in its own
   configuration. A supplementary route change was run to reach 2000 and
   2500 entries and is reported separately in section 12.
2. The molecular checkpoint was written under the SAME boundary model the
   current binary solves, `..._ghost_fixed_point_seed_reservoir_row_v3`. Only
   the atomic checkpoint is an evaluation across a boundary model change.
   Both are Mode R; only the atomic one carries the rebuilt boundary.
3. `EXHALE_RESID_QUAD` had no effect on either route. INSPECTED and MEASURED
   in section 9: the assembly selector of section 13.6 never fires on the
   evaluate or the replay route, because both leave `ieq_sweep_state_kind` at
   `ieq_state_marching`.

## 2. Classification

Every new number in this document is **Mode R** (rev9 section 3): the current
binary evaluating an immutable stored generation. No physical step was taken,
no state was published, and no file of a case directory was written. Each run
was made on a copy of the stored generation in a scratch directory of its own.

Values taken from the stored manifests, certificates and state headers are
READ and carry the historical identity. Values this work produced are
MEASURED. Statements about the source are INSPECTED with `file:line`.

## 3. Identities

### 3.1 The binary of record

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 MEASURED | `b51164071cc21f89f1c4412a8838b5a3` |
| built | 2026-09-21 07:59 (mtime MEASURED) |
| repository HEAD | `93eed8667483`, working tree dirty (MEASURED, `git rev-parse`, `git status --porcelain`) |
| sources newer than the binary | none (MEASURED, `find src Makefile -type f -newer EXHALE.x` returned nothing) |
| source digest | `972f1e6714ed455881ec388555771e77`, the md5 of the sorted md5 list of `Makefile` and the 244 `.f90`/`.inc` files of `src/` (MEASURED 2026-09-21) |

The source digest is a snapshot. This is a live tree with concurrent workers,
so it records what stood beside the binary at the time of these runs, not a
permanent property of the binary.

### 3.2 Compiler and linked libraries

Section 12 of the plan asks for these, because they set the floating-point
and reduction behavior a repeatability measurement is about. MEASURED with
`strings -a EXHALE.x | grep GCC:` and `ldd EXHALE.x`:

| | |
|---|---|
| compiler strings in the binary | `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)`; `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`; `GCC: (conda-forge gcc 16.2.0-4) 16.2.0` |
| BLAS and LAPACK | `/opt/miniconda3/lib/libopenblas.so.0` -> `libopenblasp-r0.3.34.so`, md5 `d77a762c672e6e8c39a458506f73cfd1`; the library's own string is `OpenBLAS 0.3.34 DYNAMIC_ARCH NO_AFFINITY` |
| OpenMP runtime | `/opt/miniconda3/lib/libgomp.so.1`, md5 `230180d19a34912535fea4e32e869d06` |
| Fortran runtime | `/opt/miniconda3/lib/libgfortran.so.5`, md5 `d0b302fc7706b0387920430bde531a17` |
| also linked | `libquadmath.so.0`, `libgcc_s.so.1`, `libm`, `libmvec`, `libdl`, `libpthread`, `librt`, `libc` |

`libquadmath` is linked because the kind-generic assembly of section 13.6 has
a quadruple instantiation. Section 9 shows that it was not reached here.

### 3.3 The two checkpoints

Both state pairs were copied unchanged out of `states/<generation>/` into a
scratch run directory as `output/Hydro_ioniz_IC.txt` and
`output/Ion_species_IC.txt`. The md5 of every copy MEASURED equal to the md5
the generation manifest records, which is the same value the Stage A0 record
MEASURED.

| | `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `molecular_scalar_gj1132_kzz1e9/HeH0.083` |
|---|---|---|
| generation | `g0002_20260919T004843Z_f7485b14` | `g0004_20260920T053734Z_79b42a03` |
| Stage A0 record | section 3.1 | section 3.4 |
| historical binary (READ) | `7670f310` | `75d55d9d` |
| `Hydro_ioniz.txt` md5 | `25e11562ec3cd1f6d204bf36f8e503eb` | `fd882d8918052069f8ed229dd376f8c4` |
| `Ion_species.txt` md5 | `78ac49d7e1ce90a3a4f271785400efda` | `87288b03ad88a92b7b8b8b6c2b1bdd66` |
| `input.inp` md5 | `e43c6697786711fb916d3b8eb3c113aa` | `c3916719a26a2f5a7d1d7025a7cd71fa` |
| `base.inp` md5 | (the case uses none) | `1e33d197ead9269661dc6bb0d88c0e26` |
| boundary model of the STORED state (READ, state header) | `characteristic_face_ps_reservoir_C_minus_contact_upwind_v2` | `characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3` |
| boundary model the CURRENT binary solves (MEASURED, `load_IC` report) | `..._ghost_fixed_point_seed_reservoir_row_v3` | `..._ghost_fixed_point_seed_reservoir_row_v3` |
| relation | **the current operator on an old state, boundary rebuilt** | **the current operator on a state of its own boundary model** |
| model options (READ, state header) | `mol=F carrier=F carrier_newton=F he_diff=T sec_ion=T wellbal=T` | `mol=T molbase=T carrier=T carrier_newton=F he_diff=T sec_ion=T wellbal=T` |
| replayed system (MEASURED) | 3 unknowns per cell over 500 cells, 1500 entries, species rows 0 | 3 unknowns per cell over 500 cells, 1500 entries, species rows 0 |

**Both rows are Mode R and neither is a reproduction of its historical
certificate.** For the atomic checkpoint the current binary states the
difference itself, MEASURED verbatim from the run:

```
 (load_IC) NOTE: the restart was produced under boundary model characteristic_face_ps_reservoir_C_minus_contact_upwind_v2
   and this run solves characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3. The boundary is rebuilt from the physical column
   and this run's reservoir, and its ghost composition is solved to this model's own seed and fixed point,
   so the state is loaded; its residual is not the residual it was written with.
```

For the molecular checkpoint the same reporting line reads, MEASURED:

```
 (load_IC) boundary model of the restart: characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3 (this run's).
```

This corrects the expectation the item was written with, which had both
states under `..._contact_upwind_v2`. The stored header of
`g0004_20260920T053734Z_79b42a03` says `_v3`, and so does its manifest and the
`_v3` column of the Stage A0 table in section 2.1 of the identity records.

## 4. What the instrument does

`EXHALE_RESID_DETERMINISM=1` is read at `src/EXHALE_main.f90:1395`, and the
block that follows runs to `:1533`, the `stop` at `:1532`, and then ends the run. INSPECTED, in
order:

1. **It installs the stationary reconstruction** (`:1397`, and
   `src/modules/states/stationary_operator.f90:96-99`): WENO3, `use_plm`
   false, and the PLM to WENO3 continuation disabled. That satisfies the
   endpoint restriction section 13.6 imposes, and it agrees with the `recon=WENO3`
   the two stored state headers carry.
2. **It registers the system a solve of this configuration would carry**
   (`:1400`, `call set_transported_species_rows(carrier_in_newton)`). The
   registry, `steady_newton.f90:2314-2394`, adds an element row for helium
   when `he_diffusion` and `thereis_He`, the trace element rows under
   `he_metal_diffusion`, and one row for every solved carrier, but ONLY when
   its single argument is true. With `carrier_in_newton` false the registry
   stays empty and the replayed vector is the three hydrodynamic unknowns per
   cell.
3. **It packs the state into Y** (`:1408-1412`), freezes the species box, and
   computes the state and row scales the Newton step control divides by.
4. **It builds three trial states** (`:1417-1425`): two hydrodynamic, of the
   size a line search takes (a `1e-6` move on slot 3 and a `1e-4` move on slot
   2 of every cell), and one that moves every registered species slot by
   `1e-5`.
5. **The control** (`:1429-1431`): `newton_residual` is called twice on the
   same Y with nothing in between, into two different output arrays, and
   `replay_distance` compares them. Any difference here is state the
   evaluation left behind for itself.
6. **The replay** (`:1445-1449`): F at Y, then at each of the three trials,
   then at Y again, all reusing ONE composition output workspace, so that each
   evaluation's composition is carried into the next call. The seed is named
   in every call and is the same array every time.
7. **A deliberately discarded evaluation** (`:1454-1455`):
   `eval_residual(..., state_is_discarded=.true.)` on the species trial. That
   is what a point on a dogleg ray is; the routine puts back the products it
   overwrote (`steady_newton.f90:3980-3982` and `:4507`, `put_products_back`).
8. **A Jacobian-vector probe** (`:1460-1465`): `jv_product` on a deterministic
   direction `sin(0.1 j)` normalized and scaled by the unknowns' own
   magnitudes, with the components on an active bound held.
9. **F at Y once more** (`:1466`), and two comparisons against the FIRST
   evaluation (`:1467-1470`): `Fa` against `Fb`, after the three trials, and
   `Fa` against `Fc`, after the discarded trial and the probe.
10. **The report** (`:1433-1505`): the number of entries of the full vector
    that differ at all, and the largest difference in row-scale units with its
    cell and slot, separately for the hydrodynamic rows and the species rows,
    for the control and for both replays; then the first twelve cells and
    slots that moved, with radius and scaled size, if any did.
11. **The verdict** (`:1516-1529`): the invariant is against the control and
    not against zero. `det_rs = max(control hydrodynamic, control species,
    1e-9)`, and both replay distances must be at or below it. Inside that,
    a zero entry count prints "bitwise identical" and a nonzero one prints
    "within the control".

`replay_distance` (`src/modules/time_step/steady_newton.f90:8625-8662`) is
where "every component" is enforced: it loops over all `nvar_jac*N` entries,
counts `Fa(is) .ne. Fb(is)` exactly, and tracks the maximum of
`|Fa - Fb| / max(Drow, 1e-300)` with its cell and slot in two classes, slots 1
to 3 and slots above 3.

What it does NOT do: it does not print the largest ABSOLUTE difference, only
the row-scaled one; it does not serialize the vector, so two processes cannot
be compared entry by entry through it; and it fixes no thread count of its
own. Sections 6, 7 and 8 address the first two; section 7 addresses the third.

## 5. The gate of section 4.2

Every invocation ran in a directory of its own, created fresh, holding its own
copy of the input files, its own copy of the restart pair, its own copy of the
executable and its own empty `output/`. Nothing in `LHS1140b/models/` was read
except the stored generations and the input files, and nothing there was
written. The complete environment was imposed with `env -i`, so the set
recorded IS the set in force.

Each invocation left a sidecar naming the invocation, the mode, the case, the
stage, the run directory, the executable and its md5, the complete environment
including the variables set to `0`, the md5 of every input file and of the
restart pair, the start and end time, the wall seconds, the exit status, the
size, row count and md5 of every file the invocation left, and the verdict
line it printed. The driver is
`docs/audit_20260905/repeatability_20260921/run_stage.sh`; the sidecars and
the diagnostic blocks are beside it, with `MD5SUMS.txt` over all of them.

Every invocation in this document exited with status 0 and left a nonempty
standard output containing the block the stage is expected to print. There is
no "no output", no timeout and no unresolved measurement among them.

## 6. The verdicts, verbatim

The verdict line, MEASURED, byte for byte, at **every one of the twenty
determinism invocations** of section 7, on both checkpoints:

```
   VERDICT: the residual is a state function (bitwise identical).
```

The whole diagnostic block, MEASURED, for
`atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` at `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, and byte for byte the same at the other four thread
configurations and in the second process of each:

```
 (resid_determinism) system replayed: 3 unknowns per cell over 500 cells = 1500 entries, of which species rows 0
 (resid_determinism) control, F(Y0) twice with nothing between: 0 of 1500 entries differ; row-scale max (hydrodynamic, species) with cell and slot   0.000E+00    0  0  0.000E+00    0  0
 (resid_determinism) F(Y0) re-evaluated, one workspace reused throughout:
   after three trials: 0 of 1500 entries differ
   after a discarded trust-region trial and a Jacobian-vector probe: 0 of 1500 entries differ; the probe sampled: T
   hydrodynamic rows, row-scale distance (control, after trials, after probes) with cell and slot   0.000E+00    0  0  0.000E+00    0  0  0.000E+00    0  0
   species rows, row-scale distance (control, after trials, after probes) with cell and slot   0.000E+00    0  0  0.000E+00    0  0  0.000E+00    0  0
   VERDICT: the residual is a state function (bitwise identical).
 (EXHALE_main) EXHALE_RESID_DETERMINISM=1: stopping.
```

The block for `molecular_scalar_gj1132_kzz1e9/HeH0.083` is the same text,
including the entry count: MEASURED, the ten molecular blocks and the ten
atomic blocks all have md5 `b29fd718eca3ce51e83fee3767303256`. That equality
is NOT evidence that the two checkpoints have the same residual. The block
carries no state-dependent number when nothing moved: every printed value is
either the system size, which happens to coincide, or a zero. The
state-dependent cross-process comparison is section 8.

## 7. The full-vector comparison, and the thread configurations

### 7.1 What was run

Five thread configurations, both variables always set explicitly, two
independent processes at each, on each checkpoint. Twenty determinism
invocations in all. The setup report's own statement of what the linear
algebra library was asked for is quoted from each run's standard output
(MEASURED).

| OMP_NUM_THREADS | OPENBLAS_NUM_THREADS | the setup report's two lines, MEASURED | atomic wall s (p1, p2) | molecular wall s (p1, p2) |
|---|---|---|---|---|
| 1 | 1 | `- Using  1 OMP threads (from OMP_NUM_THREADS)` / `- LAPACK library: OpenBLAS; BLAS threads OPENBLAS_NUM_THREADS stated:1 thread(s), left as stated` | 3.47, 3.49 | 10.29, 10.31 |
| 1 | 8 | `- Using  1 OMP threads (from OMP_NUM_THREADS)` / `- LAPACK library: OpenBLAS; BLAS threads OPENBLAS_NUM_THREADS stated:8 thread(s), left as stated` | 3.54, 3.55 | 10.31, 10.24 |
| 4 | 1 | `- Using  4 OMP threads (from OMP_NUM_THREADS)` / `- LAPACK library: OpenBLAS; BLAS threads OPENBLAS_NUM_THREADS stated:1 thread(s), left as stated` | 0.94, 0.94 | 2.87, 2.86 |
| 8 | 1 | `- Using  8 OMP threads (from OMP_NUM_THREADS)` / `- LAPACK library: OpenBLAS; BLAS threads OPENBLAS_NUM_THREADS stated:1 thread(s), left as stated` | 0.54, 0.53 | 1.60, 1.64 |
| 8 | 8 | `- Using  8 OMP threads (from OMP_NUM_THREADS)` / `- LAPACK library: OpenBLAS; BLAS threads OPENBLAS_NUM_THREADS stated:8 thread(s), left as stated` | 0.78, 0.56 | 1.59, 1.63 |

The complete environment of every one of them, MEASURED from the sidecars:
`PATH=/usr/bin:/bin`, `HOME`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`OMP_DYNAMIC=FALSE`, `OMP_MAX_ACTIVE_LEVELS=1`, `EXHALE_RESID_DETERMINISM=1`,
`EXHALE_RESID_QUAD=0`. Nothing else was set. The host has 72 processors
(MEASURED, `nproc`).

The wall times are the evidence that the thread count took effect and that
this is not a comparison of five runs that all used one thread: the atomic
case falls from 3.5 s to 0.54 s from 1 to 8 OpenMP threads, the molecular case
from 10.3 s to 1.6 s. The BLAS thread count does not move the wall time,
which is expected on this route: no solve runs, so the banded factorizations
that would use BLAS are never reached. It is set and recorded anyway, because
a count that is not stated is a count the library chooses, and the setup
report shows it would otherwise have been 8.

### 7.2 Every component, at every configuration

| checkpoint | configuration | entries in the vector | entries differing, control | entries differing, after three trials | entries differing, after the discarded trial and the probe | largest scaled difference, cell, slot | largest absolute difference |
|---|---|---|---|---|---|---|---|
| atomic | all five, both processes | 1500 | 0 | 0 | 0 | 0.000E+00, cell 0, slot 0 | 0 exactly |
| molecular | all five, both processes | 1500 | 0 | 0 | 0 | 0.000E+00, cell 0, slot 0 | 0 exactly |

The count is over every component: `replay_distance`
(`steady_newton.f90:8646-8656`) tests `Fa(is) .ne. Fb(is)` for every one of
the `nvar_jac*N` entries and increments on any difference whatsoever. A count
of zero therefore means that every entry compared equal bit for bit, and the
largest ABSOLUTE difference, which the instrument does not print, is exactly
zero by the same test. Cell and slot are zero because no maximum was ever
taken: the initialization `dcl = 0; jcl = 0; kcl = 0` (`:8648`) was never
overwritten.

The instrument's own comment (`EXHALE_main.f90:1506-1515`) anticipates a
control that is not clean, "cells 1-2, row-scale 1e-12", from the tables the
first evaluation of a run builds, and sets the acceptance floor at `1e-9` for
that reason. MEASURED here, on both checkpoints and at every thread
configuration, the control was exactly zero, so the floor was never used and
the comparison ran against a bitwise ideal. That makes the comment
conservative on these two states rather than wrong.

## 8. A fresh process against the in-process control

The instrument's control is an in-process comparison: `F(Y0)` twice inside one
run. Item 2 of section 5.2 asks for a fresh process, and for evidence that the
fresh process read exactly the same serialized state. Three separate
measurements, because the instrument serializes nothing of its own.

### 8.1 The state each process read

`EXHALE_DUMP_IC=1` (`EXHALE_main.f90:1212-1248`) writes the state exactly as
`load_IC` restored it, before the first equilibrium sweep, and stops. Six
independent processes per checkpoint were run this way, at `(1,1)`, `(8,1)`
and `(8,8)`, two processes each. The state is written by list-directed output
(`write_output.f90:283-350`), which round trips a binary64 value exactly, so
this is a 17-digit comparison and not a truncated one.

MEASURED: with the one provenance line that carries the run's wall-clock
timestamp excluded, all six re-serializations of `Hydro_ioniz.txt` are
byte-identical to each other, and all six of `Ion_species.txt` are byte-
identical to each other, for each checkpoint. The atomic digests are
`c958c73c03145af7a35c9c5bc7cb2647` and `ee1e83d75ce5781bdead117e898337a2`; the
molecular ones are `98c5f8a38df79bffa2d4838a802ef6de` and
`15f04c1d48effd3c0aaa3fe336fecdc6`. `Ion_species.txt` carries no timestamp
line, so its six copies are byte-identical without any exclusion at all.

Against the stored generation itself, the same re-serialization of the atomic
checkpoint compares column by column over all 504 rows (MEASURED):

| column | rows differing | largest relative difference | what it is |
|---|---|---|---|
| r | 0 of 504 | 0 | the grid |
| rho | 2 of 504 | 9.68E-02 | the two lower ghost rows, which the current boundary model rebuilds |
| v | 2 of 504 | 9.57E-01 | the same two rows |
| p | 2 of 504 | 1.40E-01 | the same two rows |
| T | 60 of 504 | 4.82E-02 | 2 rows are those same ghost rows, at 4.8E-02 and 1.9E-02; the other 58 differ at 3.9E-16 or less, one unit in the last place of the pressure to temperature map |
| heat | 504 of 504 | 1.0 | this path sets heat and cool to zero before writing; it measures nothing |
| cool | 504 of 504 | 1.0 | the same |

That is the boundary model change of section 3.3 made visible: the interior
column of the stored state is read back unchanged to the last bit, and only
the two ghost rows the current model derives for itself differ. The rebuild
is confined to the ghosts.

### 8.2 The residual each fresh process computed, entry by entry

`EXHALE_RESIDUAL=1` (`EXHALE_main.f90:1262-1342`, the `stop` at `:1341`) evaluates the stationary
residual of the loaded state and writes `output/residual_profile.txt`: one
row per physical cell with the radius, the number density, the velocity, the
temperature, the mass, momentum and energy rows, and the two face mass fluxes
the cell differences, at nine significant digits. It reaches the residual
through `assemble_residual`, the same dispatcher the Newton route's
`eval_residual` calls (rev9 section 5.1).

Twenty invocations, the same five thread configurations and two processes each
on both checkpoints. MEASURED: all ten atomic files have md5
`6a94f1dac83778329fec462072f29919` and all ten molecular files have md5
`2ab3f0f4c58f4a0908228fe54ad5dee7`. Every entry of the serialized residual
vector, 500 cells by 3 rows, together with the state and the two face fluxes,
is identical across processes and across thread counts to the nine digits the
file carries. The two checkpoints give different digests, so the comparison is
discriminating.

### 8.3 The Newton route itself, across processes

`EXHALE_NEWTON_TEST=1` (`EXHALE_main.f90:1535-1563`) calls `newton_residual`,
the exact routine the replay uses, and prints the relative residual per row.
Twenty invocations over the same matrix. MEASURED, identical in all ten runs
of each checkpoint:

| checkpoint | k=1 (mass) | k=2 (momentum) | k=3 (energy) |
|---|---|---|---|
| atomic | 9.806145E-01 | 8.723527E-06 | 1.071147E+00 |
| molecular | 9.745141E-09 | 5.136539E-13 | 7.651254E-09 |

Beside the atomic row, the stored certificate of
`g0002_20260919T004843Z_f7485b14` READS `9.806E-01`, `8.724E-06` and
`1.071E+00` for the same three rows. The current operator on the old state
reproduces the historical row maxima to the digits the certificate prints,
although its boundary was rebuilt. That is an observation, not a
reproduction of the certificate: the certificate is a measurement of the
state under the historical binary and boundary model, and these numbers are
a measurement under the current ones. It says the rebuild did not move these
three row maxima at the printed precision, and nothing more.

### 8.4 The whole standard output

The strongest available cross-process statement from the instrument itself.
MEASURED: for each checkpoint and each stage, the complete standard output of
all ten invocations is byte-identical once the two lines that state the thread
configuration are removed. A direct `diff` of the full logs of
`atomic, OMP 1, BLAS 1, process 1` against `atomic, OMP 8, BLAS 8, process 2`
returns exactly those two lines and nothing else; the two processes at one
configuration differ in nothing at all. The same holds for the molecular
checkpoint and for the `EXHALE_RESIDUAL`, `EXHALE_NEWTON_TEST` and coupled
stages.

That covers the `load_IC` report, which prints the state-dependent diffused
helium profile at the top cell, the weight of the species columns against the
density column and the cell that carries its worst value, the whole setup
report, and the diagnostic block.

## 9. The assembly selector and the sweep state kind

Section 13.6 asks that the actual `rows_kind` and `ieq_sweep_state_kind` of
every evaluation be recorded, because the selector is read at every call.

**INSPECTED**, `steady_residual.f90:244-246`: `assemble_residual` sets
`rows_kind = ROWS_PRODUCTION` and replaces it by
`generic_precision_rows_selected()` ONLY when
`ieq_sweep_state_kind .ne. ieq_state_marching`.

**INSPECTED**, `ionization_equilibrium.f90:549-552`: `ieq_sweep_state_kind` is
a saved module variable initialized to `ieq_state_marching`, which is 1, and
it is changed only by `set_ioniz_eq_sweep_state_kind`
(`ionization_equilibrium.f90:4472-4484`).

**INSPECTED**: `set_ioniz_eq_sweep_state_kind` is called from thirteen places
(MEASURED with `grep -c`), all of them in `steady_newton.f90`, and none of
them on the path of `EXHALE_RESID_DETERMINISM`, `EXHALE_RESIDUAL` or
`EXHALE_NEWTON_TEST`: `EXHALE_main.f90` contains no call to it at all, and
neither `eval_residual` (`steady_newton.f90:3826-4510`) nor `jv_product`
(`steady_newton.f90:10937-11128`) contains one. Three further places touch the
variable itself, and none is on this path either: the checkpoint of the
marching loop saves and restores it (`attempted_step.f90:846`, `:963`), its
consistency check compares it (`:1059-1060`), and one test routine,
`attempted_step_perturb_checkpointed_state_for_test` (`:1156` onward),
deliberately increments it at `:1199` so that the checkpoint comparison has
something to catch. The only reader outside the owning module is
`steady_residual.f90:245`, the selector itself.

**Therefore, at every evaluation of this work:**

| | |
|---|---|
| `ieq_sweep_state_kind` | `ieq_state_marching` (= 1) |
| `rows_kind` | `ROWS_PRODUCTION` |
| effective reconstruction | WENO3, continuation disabled, `use_plm` false, `use_weno3` true, flags consistent with the assembled endpoint (`stationary_operator.f90:96-99`) |
| transport active | atomic: none; molecular: the H2 carrier is transported (`carrier=transported` READ from the state header), but no carrier is a Newton unknown |

**MEASURED confirmation, not inference alone.** The determinism and evaluate
routes were repeated with `EXHALE_RESID_QUAD` set to `1` (quadruple rows) and
`2` (the kind-generic text at the double kind), on both checkpoints. All three
settings produced the same `residual_profile.txt`: atomic
`6a94f1dac83778329fec462072f29919` at `EXHALE_RESID_QUAD` 0, 1 and 2;
molecular `2ab3f0f4c58f4a0908228fe54ad5dee7` at the same three. If the
selector had fired, mode 1 would have assembled the flux in quadruple
precision and mode 1 and mode 0 could not agree at nine digits on a state
whose mass row stands at 9.8E-01. The selector did not fire.

This is a Mode R statement of what happened on THIS route. It says nothing
about the generic assembly itself. In particular, the mode 2 bitwise
equivalence gate of section 13.6 is **NOT PERFORMED** here, exactly as that
section says to record it: a mode 2 request evaluated in the marching state
kind exercises the production assembly, and comparing production against
production would look like a passing equivalence test.

## 10. Thread counts, and what a difference would have meant

No difference between thread counts was observed, so nothing has to be
attributed to anything. Item 4 of section 5.2 is nonetheless answered, because
the design of the measurement is what would have let it be answered.

Had a difference appeared, the candidates are reduction order, shared mutable
state, a race, initialization, and a chemical iteration that takes a different
branch. The measurement separates them as follows, and what it excludes is
stated only where it actually excludes it:

- **Reduction order and races inside the parallel region** would move the
  residual between 1, 4 and 8 OpenMP threads. MEASURED: they did not. Every
  entry of the serialized residual vector, the whole standard output and all
  three Newton row norms are identical at 1, 4 and 8 threads, and the wall
  times show the threads were used.
- **BLAS reduction order** would move it between 1 and 8 BLAS threads.
  MEASURED: it did not. This is a weak exclusion on this route, and it is
  labeled as one: no linear solve runs before the instrument stops, so the
  banded factorizations are not reached and the BLAS thread count has nothing
  to reorder. It is fixed and recorded so that the Stage B and Stage C
  measurements, which do reach them, inherit a stated value rather than a
  library default of 8.
- **Shared mutable state carried between evaluations** is what the instrument
  is built to catch, and it is the one candidate this measurement addresses
  directly: F is evaluated after three trials, after a discarded evaluation
  and after a Jacobian-vector probe, all reusing one composition workspace,
  and it came back bit for bit. MEASURED: 0 of 1500 entries moved.
- **Initialization**, the tables a first evaluation builds, is what the
  control measures. MEASURED: 0 of 1500 entries between the first and second
  evaluation, so no table building shows up in the residual at all.
- **A branch-sensitive chemical iteration** would show as a small number of
  cells moving rather than a spread: the instrument prints the first twelve
  cells and slots that moved for exactly that reason. MEASURED: none moved, so
  the branch data of every cell was identical in every evaluation and in
  every process.

Two things this does NOT exclude, stated so that they are not read into the
result: process placement and processor affinity were not varied (`env -i` left
`OMP_PROC_BIND` and `GOMP_CPU_AFFINITY` unset, and OpenBLAS reports
`NO_AFFINITY`); and every run was on one host, `lart4`, with one library set,
so nothing here speaks to reproducibility across machines or library versions.

## 11. What this adds to the suite `make test` already runs

`src/tests/residual_determinism/run.sh` runs the same instrument, so the
question is what these runs add.

READ from that driver and from `run_residual_determinism.sh`: the suite runs
two cases that live in the tree, `backup/regression/mol_base_handoff` and
`backup/regression/atomic_elem_newton`, at `OMP_NUM_THREADS=1`, with the BLAS
thread count left unset; it checks that the verdict string is present and that
the unknown count is above three on the element reload; and it then runs the
seed family of contracts 2 and 3. It does not read the vector, does not run a
second process, does not fix BLAS threading and does not vary the thread count.

MEASURED, with the binary of record and `EXHALE_TEST_OUT` pointed at a scratch
directory so that nothing of `build/tests/` was touched:

```
PASS replay_is_a_state_function_mol_base_handoff measured=state_function reference=state_function tol=0
     mol_base_handoff replays 3 unknowns per cell
PASS replay_is_a_state_function_atomic_elem_newton measured=state_function reference=state_function tol=0
     atomic_elem_newton replays 11 unknowns per cell
PASS replay_carries_the_species_rows measured=11 reference=>3 tol=0
FAIL closure_spread_within_the_row_tolerance_atomic_elem_newton measured=1.715E+04 reference=<1 tol=0
```

The failing row is contract 2, the closure spread, not contract 1, and it is
the known standing value: `docs/ISSUES_20260909.md:1026` and
`docs/lhs1140b_p5b_20260919.md:127` both READ `1.715E+04` for it. It is
unchanged and is not a finding of this work.

What this document adds, none of which the suite covers:

1. the two LHS 1140 b checkpoints the plan works on, neither of which is in
   the tree as a test case;
2. five thread configurations with the BLAS thread count fixed and recorded at
   each, instead of one configuration with it unset;
3. two independent processes at every configuration, and a cross-process
   comparison of the serialized state, the serialized residual vector and the
   whole standard output, which the verdict string alone cannot give;
4. the measurement that the assembly selector never fires on this route.

## 12. The supplementary coupled route

The item expected the molecular checkpoint to carry carrier and element rows
and therefore a longer vector. It does not, in its own configuration: both
checkpoints replay 1500 entries, because the registry
(`steady_newton.f90:2314-2394`) is gated on `carrier_in_newton`, and the
resolved configuration of both cases has it false.

To exercise the longer vector at all, both checkpoints were rerun with one
line appended to `input.inp`, `Coupled carrier solve: True`. **This is a route
change and a different experiment**, on the same immutable state, and it is
reported separately for that reason. It is not the configuration either
checkpoint was written under and none of section 7 depends on it.

MEASURED, with the same complete environment, at `(1,1)` and `(8,1)`, two
processes each:

```
 (resid_determinism) system replayed: 4 unknowns per cell over 500 cells = 2000 entries, of which species rows 1
 (resid_determinism) control, F(Y0) twice with nothing between: 0 of 2000 entries differ; row-scale max (hydrodynamic, species) with cell and slot   0.000E+00    0  0  0.000E+00    0  0
 (resid_determinism) F(Y0) re-evaluated, one workspace reused throughout:
   after three trials: 0 of 2000 entries differ
   after a discarded trust-region trial and a Jacobian-vector probe: 0 of 2000 entries differ; the probe sampled: T
   hydrodynamic rows, row-scale distance (control, after trials, after probes) with cell and slot   0.000E+00    0  0  0.000E+00    0  0  0.000E+00    0  0
   species rows, row-scale distance (control, after trials, after probes) with cell and slot   0.000E+00    0  0  0.000E+00    0  0  0.000E+00    0  0
   VERDICT: the residual is a state function (bitwise identical).
```

for the atomic checkpoint, where the one species row is the helium element row
the registry adds under `he_diffusion`, and

```
 (resid_determinism) system replayed: 5 unknowns per cell over 500 cells = 2500 entries, of which species rows 2
 (resid_determinism) control, F(Y0) twice with nothing between: 0 of 2500 entries differ; row-scale max (hydrodynamic, species) with cell and slot   0.000E+00    0  0  0.000E+00    0  0
 (resid_determinism) F(Y0) re-evaluated, one workspace reused throughout:
   after three trials: 0 of 2500 entries differ
   after a discarded trust-region trial and a Jacobian-vector probe: 0 of 2500 entries differ; the probe sampled: T
   hydrodynamic rows, row-scale distance (control, after trials, after probes) with cell and slot   0.000E+00    0  0  0.000E+00    0  0  0.000E+00    0  0
   species rows, row-scale distance (control, after trials, after probes) with cell and slot   0.000E+00    0  0  0.000E+00    0  0  0.000E+00    0  0
   VERDICT: the residual is a state function (bitwise identical).
```

for the molecular one, where the two species rows are the helium element row
and the H2 carrier row. All four atomic blocks share md5
`58585ad4ff7f08d05541165a187f6173` and all four molecular ones share
`82c1a55e09e2c1b4afaefe55e5fda9a4`; the full standard output is again
byte-identical across processes and thread counts except the two lines that
state the thread configuration.

So the species rows, when they are in the vector at all, are as reproducible
as the hydrodynamic ones. That is the part of phase 6b, the coupled Stage B,
that can be settled before the ladder runs.

## 13. Where the evidence is

`docs/audit_20260905/repeatability_20260921/`:

- `run_stage.sh`, the driver every invocation was made with;
- `sidecars/<invocation>.txt`, 88 files, one per invocation, each carrying
  what section 5 lists;
- `determinism_blocks/<invocation>.txt`, the diagnostic block of each of the
  32 determinism invocations;
- `products/*.residual_profile.txt`, the serialized residual vector of four of
  the evaluate invocations, the ones the cross-process digests of section 8.2
  are taken from;
- `MD5SUMS.txt` over all of the above.

The run directories themselves were made under the session scratch directory
and are not kept: every invocation's sidecar records its inputs by md5, and
those inputs are the stored generations, which the tree holds.

## 14. What is still unresolved

1. **The instrument cannot compare two processes entry by entry.** It prints
   counts and row-scaled maxima and serializes nothing. The cross-process
   statements of section 8 are therefore taken on three other products: the
   re-serialized state at 17 digits, the `EXHALE_RESIDUAL` vector at nine
   digits, and the whole standard output. The Newton vector `Fvec` itself is
   never written to a file by any of these routes, so "the fresh process gave
   the same F bit for bit" is supported at nine digits through a route that
   shares `assemble_residual`, and not at seventeen through the replay route.
   Closing that would need a writer at the instrument, which this item did not
   add.
2. **Only the hydrodynamic rows of these two states were replayed.** The
   element and carrier balances that the certificates of both checkpoints
   actually refuse, He/H at cell 389 and 280, H2 at cell 499, are evaluated
   inside `eval_residual` but are not packed into `Fvec` when
   `carrier_in_newton` is false, so they are outside the vector the
   instrument compares. Section 12 reaches them, but only under a route
   change. Whether the rows a certificate gates on are reproducible in the
   configuration the certificate was written under is therefore not settled
   by this work for the gated species rows.
3. **One host, one library set.** Everything here was measured on `lart4` with
   OpenBLAS 0.3.34 and the conda-forge GCC 16.2 runtime. Nothing is claimed
   about another machine, another library version, or another processor
   dispatch of the same `DYNAMIC_ARCH` library.
4. **Affinity was not varied.** `OMP_PROC_BIND` and `GOMP_CPU_AFFINITY` were
   left unset in every run and OpenBLAS reports `NO_AFFINITY`, so thread to
   processor placement was free to differ between runs. It did not produce a
   difference, which is a result; it was not controlled, which is a limit.
5. **The mode 2 equivalence gate of section 13.6 remains unperformed**, for
   the reason section 9 gives: the selector cannot fire on the routes this
   phase uses. Whoever performs it has to reach a route that sets
   `ieq_sweep_state_kind` away from the marching value first.
6. **Reproducibility is not accuracy.** A residual that is a state function
   says the derivative tests of Stage B measure the operator rather than
   noise. It says nothing about whether the operator is right, and nothing
   about whether either state has a stationary solution.
