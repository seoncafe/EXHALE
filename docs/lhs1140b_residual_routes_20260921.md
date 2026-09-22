# LHS 1140 b: residual-route equivalence, term by term

2026-09-21. Phase 2 of `PLAN_20260920_rev9.md` section 15, the measurement
section 5.1 asks for.

**Classification, rev9 section 3: every evaluation below is Mode R, the
CURRENT operator on a STORED generation.** No state was published, no
certificate was refreshed, and no golden was touched. Eleven of the
invocations perform ONE explicitly identified disposable update (one
marching step) from a restored snapshot and publish nothing; the rest take
no step at all.

**The warning that belongs in every table below.** The atomic checkpoint
was written under boundary model
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2` and the binary
of record solves
`characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3`.
Every number measured on it is the current operator on an old state and is
NOT the residual the state was written with. The molecular checkpoint was
written under `..._v3`, which IS the model the binary solves, so it carries
no such caveat; the difference between the two is used below.

Every number is labeled MEASURED (this document ran it), READ (from a named
file) or INSPECTED (executable statements, with `file:line`).

---

## 1. Identity

### 1.1 The binary of record

| item | value |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `b51164071cc21f89f1c4412a8838b5a3` (MEASURED; matches the brief) |
| bytes, mtime | 4044168, 2026-09-21 07:59 (MEASURED) |
| compiler | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0` (MEASURED, `strings`) |
| linked libraries | `libopenblas.so.0`, `libgfortran.so.5`, `libgomp.so.1`, `libquadmath.so.0`, all from `/opt/miniconda3/lib` (MEASURED, `ldd`) |
| source manifest | none exists for this md5 under `LHS1140b/models/` (MEASURED). Working tree at `93eed86` with modifications; md5 of the sorted md5 list of `src/**/*.{f90,inc}` is `30a17ea80ed173ae56e1d4367f5f85f2` (MEASURED), recorded as the source-text fingerprint in place of a manifest |
| threads | `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` in every invocation (MEASURED, recorded in each sidecar) |

### 1.2 The two checkpoints

| item | checkpoint A | checkpoint M |
|---|---|---|
| case | `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `molecular_scalar_gj1132_kzz1e9/HeH0.083` |
| generation | `g0002_20260919T004843Z_f7485b14` | `g0004_20260920T053734Z_79b42a03` |
| `Hydro_ioniz.txt` md5 | `25e11562ec3cd1f6d204bf36f8e503eb` | `fd882d8918052069f8ed229dd376f8c4` |
| `Ion_species.txt` md5 | `78ac49d7e1ce90a3a4f271785400efda` | `87288b03ad88a92b7b8b8b6c2b1bdd66` |
| both md5 verified | MEASURED before the first run and again after the last, both unchanged, both equal to the md5 the manifest READs | same |
| boundary model of the state | `..._contact_upwind_v2` (READ, state header) | `..._ghost_fixed_point_seed_reservoir_row_v3` (READ, state header) |
| boundary model of this binary | `..._ghost_fixed_point_seed_reservoir_row_v3` (READ, run log) | the same |
| writing binary | `7670f310` (READ, manifest) | `75d55d9d` (READ, manifest) |
| options (READ, state header) | `mol=F carrier=F he_diff=T sec_ion=T wellbal=T visc=F cond=F` | `mol=T molbase=T carrier=T he_diff=T sec_ion=T wellbal=T visc=F cond=F` |
| flux, reconstruction (READ, `input.inp`) | ROE, PLM | HLLC, PLM |
| full Stage A0 record | `docs/lhs1140b_identity_records_20260921.md` section 3.1 | the same document, section 3.4 |

---

## 2. The four routes, and what each of them is

INSPECTED, with the line numbers of the executable statements.

| label | route | entry | what it does |
|---|---|---|---|
| 1 | loaded-state diagnostic | `EXHALE_main.f90:1262-1340`, armed by `EXHALE_RESIDUAL=1` | optional `equilibrate_loaded_composition`; `Apply_BC`; `U_to_W`; densities; `comp_T_from_p`; `ioniz_eq`; densities; `Apply_BC`; `assemble_residual` (`:1300`); writes `output/residual_profile.txt`; stops |
| 2 | Newton | `steady_newton.f90:3826` `eval_residual`, reached from `EXHALE_main.f90:1548` `newton_residual` under `EXHALE_NEWTON_TEST=1` | `unpack_U`; densities; `Apply_BC` (`:4139`); `U_to_W`; densities; `comp_T_from_p`; `ioniz_eq`; densities; `Apply_BC` (`:4247`); **`U_to_W` + `get_species_densities` again (`:4299-4302`)**; `assemble_residual` (`:4315`); `pack_R` |
| 3 | evaluate | `EXHALE_main.f90:5543` `stationary_state_of_the_loaded_restart`, entered by `Restart intent: stationary evaluate` | densities; `Apply_BC`; `U_to_W`; densities; `comp_T_from_p`; `ioniz_eq`; densities; caloric pressure; `Apply_BC` (`:5751`); `assemble_residual` (`:5759`); certification; products |
| 4 | marching operator | `reconstruction_continuation_rhs` at `EXHALE_main.f90:2284`, `:2329`, `:2365`, bracketed by the update-map snapshots `u_umA` (`:2479`), `u_umB` (`:3063`), `u_umC` (`:3073`), `u_umD` (`:3095`), armed by `EXHALE_UPDATE_MAP` | ADVANCES a state. `update_map_begin_step` (`:6205`, its `assemble_residual` at `:6271`) assembles the stationary residual of the frozen state by the same sequence as route 1, then the production step runs and `update_map_end_step` (`:6460`) writes `R`, `G_dt = (u_after - u_before)/dt` and the operator-split stage rates per cell to `output/update_map.txt` |

**Route 4 cannot be evaluated without advancing time**, so, as the brief
requires, every route-4 comparison below is labeled a RIGHT-HAND SIDE
comparison and never a stationary residual comparison. The right-hand side
is recovered as the dt -> 0 limit of `G_dt`; the dt ladder is measured, not
assumed.

**One fact that bounds the whole exercise** (INSPECTED). `assemble_residual`
(`steady_residual.f90:239-247`) replaces `ROWS_PRODUCTION` by
`generic_precision_rows_selected()` only when
`ieq_sweep_state_kind .ne. ieq_state_marching`. That variable initializes to
`ieq_state_marching` (`ionization_equilibrium.f90:552`) and the only setter,
`set_ioniz_eq_sweep_state_kind` (`:4472`), is called from
`steady_newton.f90` (lines 6337, 6447, 6495, 6560, 6721, 6963, 10661, 10718,
16701, 16956, 17006) and from `src/tests/krylov_and_dogleg/`, and from
NOWHERE ELSE. None of routes 1, 2, 3 or 4 as invoked here passes through
those call sites. Section 7 below reports the measurement that confirms it.

---

## 3. The protocol actually run

### 3.1 Non-mutation, and how it was obtained

The plan asks that the pre-route state be saved and restored between routes.
This was done in the stronger form available to a fixed binary: **each route
ran in its own operating-system process, in its own directory, with its own
freshly copied restart pair.** No module workspace, cache, ghost state,
composition, reconstruction mode or branch flag can survive from one route
to the next, because no route shares an address space with another.

- Each invocation has its own directory `<scratch>/p2/models/grp/<name>/`
  with its own `input.inp`, its own `output/` and its own
  `Hydro_ioniz_IC.txt` / `Ion_species_IC.txt` copied from the stored
  generation.
- The restart pair's md5 was MEASURED before and after every invocation.
  **All 24 invocations: unchanged** (recorded in each `SIDECAR.txt`).
- The two stored generations' md5s were MEASURED before the first run and
  after the last: unchanged, and equal to the md5 the manifest READs.
- **Nothing was written into any `LHS1140b/models/<case>/` directory.**
  The case directories were opened read-only.

The within-process question, which the fresh process cannot answer, is
measured separately: the evaluate route assembles the residual TWICE on the
same conserved array, once before and once after the face-flux report
installs a boundary of its own (`EXHALE_main.f90:5759` and `:5848`,
labels `C_residual_rows` and `F_residual_rows_after_report`). Section 6.3
reports it.

### 3.2 The gate of section 4.2

Every invocation carries `SIDECAR.txt` naming: the invocation, the working
directory, the executable and its md5, the case and generation, the mode,
the restart intent, the input md5, the complete environment set in force
(including the variables set to `0`), the restart-pair md5 before and after,
the exit status and its class, and, for every file the run wrote, its bytes,
line count, md5 and first header line. All 24 invocations exited 0 with a
normal stop; there is no timeout and no unresolved measurement in this
document.

### 3.3 The commands

The run directory is built by `<scratch>/p2/mkrun.sh` and the invocation by
`<scratch>/p2/drive.sh`; both are kept with the artifacts. What they do:

```bash
# one run directory = one invocation
rm -rf   $D && mkdir -p $D/output
\cp -f   $CASE/input.inp $D/input.inp
sed -i   "s|^Restart intent:.*|Restart intent: $INTENT|" $D/input.inp
\cp -f   $CASE/base.inp $D/base.inp                   # molecular case only
\cp -f   $CASE/states/$GEN/Hydro_ioniz.txt $D/output/Hydro_ioniz_IC.txt
\cp -f   $CASE/states/$GEN/Ion_species.txt $D/output/Ion_species_IC.txt
cd $D
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export <the invocation's own variables>
timeout 1800 /nfs/.../EXHALE_v1.00/EXHALE.x > run.log 2> run.err
```

The run directories sit three levels below `<scratch>/p2/`, and
`<scratch>/p2/sed` is a symbolic link to `LHS1140b/sed/`, so the case's own
`Spectrum file: ../../../sed/...` resolves without editing that line. The
spectrum file is therefore the case's own, unmodified.

The invocations:

| name | checkpoint | intent | variables |
|---|---|---|---|
| `A1_residual` | A | stationary evaluate | `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1` |
| `A1e_residual_eqrel` | A | stationary evaluate | `EXHALE_RESIDUAL=1 EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1` (no `EXHALE_RELOAD_EQ`, so the loaded composition IS equilibrated first) |
| `A2_newton` | A | stationary evaluate | `EXHALE_NEWTON_TEST=1 EXHALE_RELOAD_EQ=0 EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1` |
| `A3_evaluate` | A | stationary evaluate | `EXHALE_RELOAD_EQ=0 EXHALE_EVAL_STATE_DUMP=1 EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1` |
| `A4_updatemap`, `A4_um_10e1`, `A4_um_10e3`, `A4_um_10e5` | A | relaxation | `EXHALE_UPDATE_MAP=1.0e-2 / 1.0e-1 / 1.0e-3 / 1.0e-5`, `EXHALE_RELOAD_EQ=0` |
| `A5_um_explicit` | A | relaxation | `EXHALE_UPDATE_MAP=1.0e-1 EXHALE_RELOAD_EQ=0`, input with `Energy solver: Explicit` appended |
| `A6_cfl1_si`, `A6_cfl1_ex` | A | relaxation | `EXHALE_UPDATE_MAP=1.0 EXHALE_RELOAD_EQ=0`, the second with `Energy solver: Explicit` |
| `A7_quad2`, `A8_quad1` | A | stationary evaluate | `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 EXHALE_RESID_QUAD=2 / 1` |
| `A9_newton_quad1` | A | stationary evaluate | `EXHALE_NEWTON_TEST=1 EXHALE_RELOAD_EQ=0 EXHALE_RESID_QUAD=1` |
| `A10_determinism` | A | stationary evaluate | `EXHALE_RESID_DETERMINISM=1 EXHALE_RELOAD_EQ=0` |
| `M1_residual` | M | stationary evaluate | `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 EXHALE_BOUNDARY_TRACE=1` |
| `M1b_entrybnd` | M | stationary evaluate | the same plus `EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition` |
| `M2_newton` | M | stationary evaluate | `EXHALE_NEWTON_TEST=1 EXHALE_RELOAD_EQ=0` |
| `M3_evaluate` | M | stationary evaluate | `EXHALE_RELOAD_EQ=0 EXHALE_EVAL_STATE_DUMP=1 EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1` |
| `M4_um1e1`, `M4_um1e3` | M | relaxation | `EXHALE_UPDATE_MAP=1.0e-1 / 1.0e-3`, `EXHALE_RELOAD_EQ=0` |
| `M5_ex_10e1`, `M5_ex_10e3` | M | relaxation | the same with `Energy solver: Explicit` |
| `M7_quad2` | M | stationary evaluate | `EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 EXHALE_RESID_QUAD=2` |

`EXHALE_RELOAD_EQ=0` is set wherever the route would otherwise put the
loaded composition on its own fixed point first, so that routes 1, 2, 3 and
4 all measure the state the file carries. `A1e_residual_eqrel` is the one
invocation that does not, and it exists to size that policy.

---

## 4. Term-by-term table, checkpoint A (atomic)

All four routes evaluated the same stored state. Boundary caveat of the
header applies to every row of this table.

| term | route 1 | route 2 | route 3 | route 4 | verdict |
|---|---|---|---|---|---|
| row measure, mass `max_j \|R_1\|/s_1` | 9.806145E-01 | 9.806145E-01 | 9.806145E-01 | (route 4's own stationary residual, same value) | AGREE, every printed digit |
| row measure, momentum | 8.723527E-06 | 8.723527E-06 | 8.723527E-06 | " | AGREE |
| row measure, energy | 1.071147E+00 | 1.071147E+00 | 1.071147E+00 | " | AGREE |
| mass row vector `R_1(j)`, j = 1..500 | `residual_profile.txt` | not emitted | not emitted | `update_map.txt` column `R` | route 1 vs route 4: **0 of 500 entries differ** |
| momentum row vector `R_2(j)` | " | not emitted | not emitted | " | **0 of 500 differ** |
| energy row vector `R_3(j)` | " | not emitted | not emitted | " | **0 of 500 differ** |
| cellwise maximum and its cell | mass 9.8061E-01 at cell 2 (r 1.00039), momentum 8.7235E-06 at cell 1, energy 1.0711E+00 at cell 1 | -- | identical, cell for cell | -- | AGREE |
| integrated norm, 3 rows | 1.2851E-01 / 8.1806E-07 / 2.4773E-01 | -- | identical | -- | AGREE |
| window contributions, 4 windows x 3 rows (24 numbers) | printed | -- | identical in all 24 | -- | AGREE |
| worst cell signed terms: `R`, flux divergence, source, heat, cool, heat-cool, scale, rho, rho v | -6.551809E-05, -6.116629E-05, 0.0, 4.366900E-06, 1.510249E-08, 4.351798E-06, 6.116629E-05, 3.489290E+00, -1.180386E-05 | -- | identical in all nine | -- | AGREE |
| base face mass fluxes, cell 1 | `Phi_lo` 4.23009766E-07, `Phi_hi` 2.65518592E-07 (these carry the `r_edg^2` factor) | -- | `F_lower` 1.11839E-07-class column: 4.2292800228252566E-07, `F_upper` 2.6536466983879691E-07 (no `r_edg^2`) | -- | AGREE after the geometric factor: ratios 1.00019335 and 1.00058 reproduce `r_edg(0)^2` and `r_edg(1)^2` of this grid |
| base continuity row and its scale and floor, cells 1..12 | -- | -- | `boundary_trace.txt`, two independent assemblies | -- | see section 6.3 |
| base face state gap (cached vs the boundary of the installed composition) | 0.0000000000000000E+000 | not emitted | 0.0000000000000000E+000 twice | 0.0000000000000000E+000 | AGREE, exactly zero everywhere |
| prescribed reservoir (version 1) p, T, particles per unit mass, level radius | -- | -- | 1.0000000001000000E+00, 1.0, 2.7072294627320281E-01, 1.0 | identical | AGREE |
| solved ghost counts | not emitted | not emitted | heavy 1.0458870006292493E+00, electrons 1.1039102313710122E-07 | not emitted | route 3 only; see section 8 |
| particle count, temperature, composition, per cell | not emitted | not emitted | `eval_state_dump.txt`, 2008 lines, 3 blocks at 17 digits | not emitted | route 3 only; see section 8 |
| heat, cool arrays | worst cell only | not emitted | full column in `Hydro_ioniz.txt` | full column | route 3 vs route 4 agree to 1E-9 relative in cells 1..4; route 4's state has taken one step, which is the difference |
| MARCHING RIGHT-HAND SIDE, mass | -- | -- | -- | `\|G_dt + R\|`: 4.927E-05 (dt 1e-1 CFL), 4.991E-06 (1e-2), 4.988E-07 (1e-3), 2.080E-07 (1e-5) | first order in dt, **converges to zero**; the 1e-5 value is the rounding floor of the difference quotient |
| MARCHING RIGHT-HAND SIDE, momentum | -- | -- | -- | 1.652E-05, 1.637E-06, 1.635E-07, 1.636E-09 | first order, **converges to zero** |
| MARCHING RIGHT-HAND SIDE, energy | -- | -- | -- | 3.643E-05, 7.603E-06, 4.807E-06, 4.569E-06 | **does NOT converge to zero**; plateau 4.57E-06. Attributed in section 6.4 |
| composition policy, route 1 without vs with `equilibrate_loaded_composition` | mass 0/500 differ, momentum 0/500, energy 136/500 differing by at most 2.0E-21 absolute; n, v, T, `Phi_lo`, `Phi_hi` 0/500 | | | | the policy costs 2E-21 in the energy row at this checkpoint |
| assembly selector `EXHALE_RESID_QUAD` = 1, = 2 | 0 of 500 entries differ from the unset run, in every row and in both face-flux columns | 0 change in the printed norms | | | the selector is INERT on these routes; section 7 |

Two facts about precision belong with this table. `residual_profile.txt` and
`update_map.txt` both serialize at `ES16.8`, so "0 of 500 entries differ"
means the two routes' row vectors are identical in all nine significant
digits that are written, not that they were compared as raw doubles.
`boundary_trace.txt` and `eval_state_dump.txt` carry 17 digits and are
compared as written.

---

## 5. Term-by-term table, checkpoint M (molecular)

This checkpoint's stored boundary model IS the one the binary solves, so the
header caveat does not apply to it. `mol=T carrier=T`, so it carries a
carrier row the atomic checkpoint does not.

| term | route 1 | route 2 | route 3 | route 4 | verdict |
|---|---|---|---|---|---|
| row measure, mass | 1.392918E-08 | **9.745141E-09** | 1.392918E-08 | -- | **routes 1 and 3 agree; route 2 DIFFERS by a factor 0.6996** |
| row measure, momentum | 4.875714E-13 | **5.136539E-13** | 4.875714E-13 | -- | **route 2 DIFFERS by a factor 1.0535** |
| row measure, energy | 7.651254E-09 | 7.651254E-09 | 7.651254E-09 | -- | AGREE, every printed digit |
| mass, momentum, energy row vectors | `residual_profile.txt` | not emitted | not emitted | `update_map.txt` column `R`, at dt 1e-1 and 1e-3 | route 1 vs route 4: **0 of 500 entries differ in every row, at both dt** |
| assembly selector `EXHALE_RESID_QUAD=2` | 0 of 500 differ in every row | -- | -- | -- | INERT; section 7 |
| boundary derived from the ENTRY composition instead of the post-sweep one (`EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition`) | mass 4.269520E-02, momentum 3.625712E-06, energy 2.844228E-02 | -- | -- | -- | **six decades** in the mass row against the 1.39E-08 of the same route with the rebuild in place |
| base face state gap | 0.0 | not emitted | 0.0 twice | 0.0 at the update-map assembly | AGREE |
| "cached base face state belonged to another composition at a read and was recomputed" | absent; the evaluate route's own counter READs 0 | **present** | absent; counter 0 | present once, at the certification AFTER the step, not at the update-map assembly | this is the mechanism of the route-2 difference; section 6.2 |
| refusing certification entries | -- | -- | carrier balance H2 8.150E-04 at cell 499; elemental transport He/H partition 2.503E-04 at cell 280 | -- | reproduces the historical certificate READ in `lhs1140b_identity_records_20260921.md` section 2.1, digit for digit |
| closure defect of the work state | -- | -- | composition against the sweep's own root 1.459E-13; largest accepted reaction residual 7.430E-17; cells without a chemical root 0; thermal 1.454E-13 | -- | route 3 only |
| CARRIER row, stationary residual `max_j \|R_c\|` | -- | -- | -- | 5.9005E+04 (identical at both dt) | -- |
| CARRIER row, marching rate `max_j \|G_c\|` | -- | -- | -- | 1.0346E+05 (dt 1e-1), 5.9462E+04 (dt 1e-3) | -- |
| CARRIER row, `\|G_c + R_c\|` | -- | -- | -- | 4.4457E+04 (dt 1e-1), 4.5690E+02 (dt 1e-3) | first order in dt, **converges to zero**: the carrier transport operator and `carrier_steady_residual` are the same map |
| MARCHING RIGHT-HAND SIDE, mass | -- | -- | -- | `\|G+R\|` 8.107E-05 (dt 1e-1), 8.246E-07 (1e-3) | first order in dt; the limit is not resolved below 8.2E-07, five decades above `\|R\|` = 8.06E-12. **No disagreement above that floor** |
| MARCHING RIGHT-HAND SIDE, momentum | -- | -- | -- | 1.810E-05, 1.885E-07 | the same statement, against `\|R\|` = 1.12E-11 |
| MARCHING RIGHT-HAND SIDE, energy | -- | -- | -- | `\|G+R\|` 1.409E-04 (dt 1e-1), **3.906E-03 (dt 1e-3)**; the stage rates are hydro 3.748E-03 and chem 3.888E-03 at dt 1e-1, hydro 3.887E-03 and chem 1.576E-04 at dt 1e-3 | **grows as dt falls**, against `\|R\|` = 2.385E-11. Attributed in section 6.4 |
| the same with `Energy solver: Explicit` | -- | -- | -- | `\|G+R\|` energy 1.409E-04 (dt 1e-1), **1.433E-06 (dt 1e-3)**; chem 3.888E-03 at both dt | with the explicit source update the energy right-hand side **converges to zero**, first order in dt |

---

## 6. Every difference, and which cause it belongs to

The plan names three candidate causes: the composition and closure policy,
the assembly selector of section 13.6, and the module state a route leaves
behind. Each difference below is assigned to one of them, or explicitly left
open.

### 6.1 The composition policy of the loaded state, checkpoint A

MEASURED: running route 1 with and without `equilibrate_loaded_composition`
leaves the mass and momentum row vectors identical (0 of 500), and moves 136
of 500 entries of the energy row by at most 2.0E-21 absolute; the state
columns `n`, `v`, `T` and both face mass fluxes are identical. The evaluate
route READs the reason: the loaded composition stands 3.295E-15 from the
sweep's own root at this checkpoint (1.459E-13 at checkpoint M), so the
sweep has nothing to move.

**Cause: composition and closure policy. Size: negligible at both
checkpoints.** It would not be negligible at a state whose stored
composition is far from its own root, and this measurement does not bound
that case.

### 6.2 Route 2's mass and momentum rows at checkpoint M

MEASURED: 9.745141E-09 against 1.392918E-08 (mass) and 5.136539E-13 against
4.875714E-13 (momentum), with the energy row identical to seven digits.
Route 2's log carries, and routes 1, 3 and the update-map assembly of route
4 do not, the line

> `[base boundary] the cached base face state belonged to another composition at a read and was recomputed`
> `(a caller refreshed the composition without deriving the boundary again; counted for the run)`

INSPECTED, the mechanism: `eval_residual` derives the boundary at
`steady_newton.f90:4247` and then, at `:4299-4302`, calls `U_to_W` and
`get_species_densities` AGAIN on the refreshed ghosts before
`assemble_residual` at `:4315`. That second refresh moves the
composition-derived module state, the caloric arrays among it, after the
boundary was derived, so the base face state the reconstruction reads is
recomputed at the read. The loaded-state route (`EXHALE_main.f90:1297-1300`)
and the evaluate route (`:5751-5759`) call `Apply_BC` and go straight to
`assemble_residual` with no refresh in between. The comment at
`steady_newton.f90:4288-4298` states the purpose of that refresh: it makes
the particle count the residual uses the one formed from the refreshed
ghosts, and it says that "without conduction or viscosity the count enters
no row and nothing changes". That last clause is where the measurement
disagrees: this case has `visc=F cond=F`, and the row DOES change, because
what moved is not only the particle count but the base face state the
reconstruction reads.

The same difference is absent at checkpoint A, where
`caloric_mixture_active` is false on every cell (READ, route 3: "the mixture
is atomic on every cell, where the two rules are the same map") and the
composition refresh therefore does not move the caloric map.

The size is bounded from above by the same route's boundary-order
experiment: deriving the boundary from the entry composition instead of the
post-sweep one moves the mass row from 1.392918E-08 to 4.269520E-02, six
decades. The route-2 difference is a factor 1.43 on the same row, so it is a
small second-order effect of the same variable and not an independent one.

**Cause: composition and closure policy, specifically WHICH composition the
base face state belongs to at the read.** Attributed.

### 6.3 Module state left behind

Two measurements.

Between routes: none is possible, because every route ran in its own
process. That is stated, not measured.

Within a process: the evaluate route assembles the residual twice on the
same conserved array, before and after the face-flux report installs a
boundary of its own, and `EXHALE_BOUNDARY_TRACE=1` records the base
continuity row, both face mass fluxes, the row scale and the row's rounding
floor for cells 1..12 at each.

- Checkpoint A: **all 12 rows bit-identical in all five columns**, and both
  face-state gaps exactly 0.0.
- Checkpoint M: **11 of 12 rows bit-identical**; cell 1 differs in the
  ROUNDING-FLOOR column alone, 1.9412824659667027E-09 against
  1.9412824659667043E-09, one unit in the last place, while its signed row,
  both face fluxes and its scale are bit-identical. Both face-state gaps
  exactly 0.0.

**Cause: module state left behind. Size: zero in the rows; one ulp in a
diagnostic quantity at one cell.** The source comment at
`EXHALE_main.f90:5845` says the second assembly "need not reproduce"
the first; at these two checkpoints it does.

### 6.4 The marching operator's energy row loses the radiative source

This is the one large disagreement, it appears at BOTH checkpoints, and it
is fully attributed.

MEASURED, checkpoint A, `|G_dt + R|` of the energy row against dt:

| dt / CFL | 1e-1 | 1e-2 | 1e-3 | 1e-5 |
|---|---|---|---|---|
| mass | 4.927E-05 | 4.991E-06 | 4.988E-07 | 2.080E-07 |
| momentum | 1.652E-05 | 1.637E-06 | 1.635E-07 | 1.636E-09 |
| energy | 3.643E-05 | 7.603E-06 | 4.807E-06 | **4.569E-06** |

The mass and momentum rows fall linearly in dt and reach their rounding
floor. The energy row stops at 4.57E-06, which is `max_j |heat - cool|` of
the same state (4.351798E-06 at cell 1, READ from the route-1 and route-3
worst-cell breakdown). Cell by cell at dt = 1e-5 CFL, `G_dt(3,1)` = 6.09E-05
against `-(dF - S)` = 6.116629E-05, so the marching right-hand side of the
energy row converges to `-(dF - S)` and NOT to
`-(dF - S - (heat - cool)) = -R`.

MEASURED, the same comparison with `Energy solver: Explicit` at dt = 1e-1
CFL, checkpoint A: the "chem" stage rate of the energy row at cell 1 is
4.35180523E-06, which is the state's own `heat - cool` = 4.351798E-06; with
the default semi-implicit solver it is -1.118E-10. At a FULL CFL step the
max over cells is 4.527E-06 explicit against 8.921E-07 semi-implicit, and
**128 of 500 cells, every cell below r = 1.042, keep less than one part in a
thousand of the source** while 372 cells keep more than 90 per cent of it.

MEASURED, checkpoint M, where the state is close to local energy balance so
the effect is large: with the default solver `|G+R|` of the energy row is
1.409E-04 at dt = 1e-1 CFL and 3.906E-03 at dt = 1e-3 CFL, growing as dt
falls; with `Energy solver: Explicit` it is 1.409E-04 and 1.433E-06,
falling first order in dt. The "chem" stage rate is 3.888E-03 at both dt
under the explicit update and drops to 1.576E-04 under the semi-implicit
one.

INSPECTED, the mechanism. `energy_balance_init`
(`energy_semi_implicit.f90:386-398`) forms

```
R0 = internal_energy_per_particle(T_old) - u_old + a*(cool0 - heat)
sc = |u_old| + a*(|heat| + |cool0|)
if (abs(R0) .le. energy_res_tol*sc) then    ! energy_res_tol = 1.0d-9, :128
   T_new = T_old ;  cycle                   ! the cell is declared already balanced
```

with `a = dt/(n_tot + n_e)` (`:801`). Where the composition did not move,
the first two terms cancel and `R0 = a*(cool - heat)`, which is proportional
to dt, while `sc` is dominated by the cell's own energy per particle, which
is not. So the test `|R0| <= 1e-9 * sc` becomes true for every cell once dt
is small enough, and the semi-implicit update then returns `T_new = T_old`
and the radiative source never reaches `u(3,:)`. At checkpoint A, cell 1,
dt = 1 x CFL: `|R0|/sc` is about 2.5E-10, below the 1E-9 threshold, so the
source is dropped at the production step size and not only in the
diagnostic's shortened ones.

**Cause: composition and closure policy, specifically the convergence test
of the semi-implicit energy update.** It is NOT the assembly selector and
NOT leftover module state; the stationary residual assembly is identical in
routes 1, 3 and 4, as section 4 and section 5 show.

**What it means for the plan.** The marching operator and the stationary
residual do not share the energy row wherever this test fires, so they
cannot share a fixed point there. The affected region at checkpoint A is
r < 1.042, which contains the two cells at which the hydrodynamic rows
refuse. The update-map comparison of the energy row therefore cannot be made
to converge by reducing dt, and any future use of `EXHALE_UPDATE_MAP` on the
energy row has to state which energy solver was in force. This document does
not repair anything: it names the location, the constant and the two
measurements that isolate it, and leaves the decision to the plan.

### 6.5 Route 4's stage labels do not match the code

INSPECTED and MEASURED, noticed while attributing section 6.4 and reported
here because a future reader of `update_map.txt` will be misled by it.

`update_map_end_step` (`EXHALE_main.f90:6460`, its column list at `:6466-6472` and the `energy` label at `:6470`) documents its five stage
columns as `hydro`, `chem`, `energy`, `transport`, `filter`, and says that
`energy` is "the semi-implicit heating-cooling step". The snapshots are
taken at `u_umA` (`:2479`, after the RK3), `u_umB` (`:3063`), `u_umC`
(`:3073`) and `u_umD` (`:3095`), while `solve_energy_semi_implicit` is
called at `:2776`, that is BETWEEN `u_umA` and `u_umB`. The only statement
between `u_umB` and `u_umC` is `call Apply_BC(u)` (`:3067`), which writes the
ghosts and returns the interior unchanged.

MEASURED, consistent with that: the `energy` column is exactly 0.000E+00 in
every row of every one of the eleven update-map invocations, on both
checkpoints and at every dt, and the radiative source appears in the `chem`
column instead. The `chem` column therefore carries the composition sweep,
the element diffusion, the photochemical transport AND the heating and
cooling together, which is not what the header says and is not a separation
a reader can undo.

**Reported, not fixed:** the file is `EXHALE_main.f90`, which this item does
not own.

---

## 7. The assembly selector was never exercised, and that is a correction

Section 5.1 of the plan says that two callers can share `assemble_residual`
and still evaluate different assemblies, and section 13.6 makes the selector
the third gate in front of everything. At these checkpoints, on these four
routes, **the selector cannot fire at all.**

INSPECTED: `assemble_residual` (`steady_residual.f90:244-247`) consults
`generic_precision_rows_selected()` only when `ieq_sweep_state_kind` is not
`ieq_state_marching`. That variable is initialized to `ieq_state_marching`
(`ionization_equilibrium.f90:552`); its only setter,
`set_ioniz_eq_sweep_state_kind` (`:4472`), is called from eleven sites, all
inside `steady_newton.f90`, and from `src/tests/krylov_and_dogleg/`. Routes
1, 2, 3 and 4 as invoked here reach none of them: route 2 enters
`eval_residual` through `newton_residual` from `EXHALE_main.f90:1548`,
outside any solver pass, and there is no setter inside `eval_residual`
(MEASURED by grep over `steady_newton.f90:3826-4510`).

MEASURED, the confirmation:

| invocation | `EXHALE_RESID_QUAD` | result |
|---|---|---|
| `A7_quad2` | 2 | `residual_profile.txt` differs from the unset run in 0 of 500 entries, in all three rows and both face-flux columns |
| `A8_quad1` | 1 | the same, 0 of 500 |
| `A9_newton_quad1` | 1 | the three printed norms identical to `A2_newton` |
| `M7_quad2` | 2 | 0 of 500 in all three rows |

So, per section 13.6's own rule: **the generic comparison is NOT PERFORMED.**
Every row in every table of this document was assembled by the PRODUCTION
routines, `rows_kind = ROWS_PRODUCTION`, `ieq_sweep_state_kind =
ieq_state_marching`. The actual values were established by inspection of the
setters plus the four null experiments above, because the binary of record
prints neither variable; section 8 records that as a gap.

Two consequences for the plan. First, phase 2's ambiguity about the selector
is closed for these routes: it is not a cause of anything measured here, and
it cannot be, whatever the environment says. Second, a future phase that
wants the mode-2 bitwise gate of section 13.6 has to arrange for the state
kind to be a steady one, which today means entering a solver pass, not
setting an environment variable.

---

## 8. What could not be measured, and why

1. **Route 2's residual as a vector.** `EXHALE_NEWTON_TEST`
   (`EXHALE_main.f90:1535-1562`) packs `Fvec` into `Rres` and prints
   `residual_norms` alone; no writer of the vector cell by cell is reachable from
   `newton_residual`, and `output/residual_profile.txt` is written only
   inside the `EXHALE_RESIDUAL` branch. So routes 1, 3 and 4 were compared
   entry by entry over 500 cells and route 2 only through three row measures
   and the boundary-recompute line. Where route 2 differs (section 6.2), the
   cell at which it differs is therefore NOT known. Closing this needs one
   writer call in the `EXHALE_NEWTON_TEST` block, which this item does not
   own.

2. **The unknown space route 2 carries.** MEASURED, `A10_determinism`:
   "3 unknowns per cell over 500 cells = 1500 entries, of which species rows
   0". The `EXHALE_NEWTON_TEST` and `EXHALE_RESID_DETERMINISM` hooks
   therefore replay the three hydrodynamic rows only, while checkpoint A's
   certification carries an elemental transport He/H partition row and
   checkpoint M carries that plus a carrier balance H2 row. The carrier and
   element rows were compared between routes 3 and 4 (certification and
   update map) and NOT against route 2.

3. **`rows_kind` and `ieq_sweep_state_kind` as printed values.** Section
   13.6 asks that the ACTUAL values be recorded at every evaluation. The
   binary prints neither. They are established here by inspection of every
   setter plus four null experiments (section 7), which is weaker than a
   printed value and is labeled as such.

4. **The composition sweep count and the seed identity per route.** No
   route prints the number of `ioniz_eq` sweeps it ran or a hash of the seed
   it handed them. What IS available and was recorded is the evaluate
   route's closure block: the distance of the loaded composition from the
   sweep's own root, the largest accepted reaction residual, the count of
   cells without a chemical root, and the thermal defect. Routes 1, 2 and 4
   print none of those.

5. **The ghost state as a vector, for routes 1, 2 and 4.**
   `ghost_record_of_this_evaluation` is called from the evaluate route only
   (`EXHALE_main.f90:5761`), so `output/ghost_record.txt` exists for route 3
   alone. The cross-route ghost comparison rests instead on the base face
   mass fluxes (route 1 and route 3, agreeing after the `r_edg^2` factor),
   on the face-state gap (exactly zero on routes 1, 3 and 4) and on the
   prescribed reservoir block (identical on routes 3 and 4).

6. **`residual_profile.txt` and `update_map.txt` are 9-digit
   serializations.** "0 of 500 entries differ" is a statement at `ES16.8`,
   not a raw-double comparison. The 17-digit comparisons in this document are
   the ones taken from `boundary_trace.txt`.

7. **Repeatability across processes and thread counts** is phase 3 and was
   not attempted. The one repeatability number taken here is
   `A10_determinism` at one thread: at checkpoint A the route-2 residual is
   bitwise identical to itself over the control, after three trials, and
   after a discarded trust-region trial and a Jacobian-vector probe, 0 of
   1500 entries differing in each. That bounds the noise floor of the
   route-2 comparison at zero, so the differences of section 6.2 are not
   noise.

8. **Generality.** Two checkpoints of one planet. The atomic one carries the
   `_v2` boundary caveat, the molecular one does not. Nothing here bounds a
   state whose stored composition is far from its own root, or a
   configuration with viscosity or conduction on, where the particle-count
   refresh of `steady_newton.f90:4299-4302` enters the rows directly.

---

## 9. Two observations that qualify the plan's text

**The current operator reproduces the historical certificates.** Rev9
section 3 states that a run which rebuilds the boundary is the current
operator on an old state and that "its residual is not the residual the
restart was written with". MEASURED, at checkpoint A, under the rebuilt
`_v3` boundary and a different binary: the refusing rows are
mass 9.806E-01 at cell 2, momentum 8.724E-06 at cell 1, energy 1.071E+00 at
cell 1, which reproduce the historical certificate READ in
`lhs1140b_identity_records_20260921.md` section 2.1 in every printed digit,
including the cell indices. The same holds at checkpoint M for its carrier
(8.150E-04 at cell 499) and element (2.503E-04 at cell 280) rows. The plan's
statement is correct as a statement of principle and should not be read as a
prediction that the numbers move: at these two checkpoints they do not, at
the four digits a certificate prints.

**The boundary order, not the boundary model, is what moves these rows.**
MEASURED at checkpoint M: leaving the boundary of the entry composition
standing instead of deriving it from the composition the residual is
assembled with moves the mass row from 1.392918E-08 to 4.269520E-02 and the
energy row from 7.651254E-09 to 2.844228E-02, six decades in both. Any
future comparison of two routes on a molecular state has to record which
composition the boundary belongs to before anything else, because nothing
else measured here is within six decades of that.

---

## 10. Verdict

**Routes 1, 3 and 4 assemble the same stationary residual, term for term.**
At both checkpoints the mass, momentum and energy row vectors of route 1 and
route 4 are identical in all 500 cells; route 3 reproduces route 1's three
row measures, its three cellwise maxima and their cell indices, its three
integrated norms, all 24 window contributions and all nine worst-cell signed
terms; the base face mass fluxes agree after the recorded geometric factor;
the face-state gap is exactly zero on all three. That closes the ambiguity
section 5.1 opened for these three routes.

**Route 2 joins them at the atomic checkpoint and departs at the molecular
one**, by a factor 1.43 in the mass row and 1.0535 in the momentum row, with
the energy row identical. The cause is identified and inspected: the second
composition refresh at `steady_newton.f90:4299-4302`, after the boundary was
derived at `:4247`, which makes the base face state be recomputed at the
read. It is a composition and closure policy difference, it is absent where
the mixture is atomic, and it is three decades smaller than the full
boundary-order effect on the same row.

**The marching operator's right-hand side equals minus the stationary
residual in the mass, momentum and carrier rows**, to the dt -> 0 limit,
first order in dt, at both checkpoints. **It does not in the energy row**,
and the cause is the semi-implicit energy update's convergence test
(`energy_semi_implicit.f90:128, 386-398`), which drops the radiative source
in any cell where the source's energy change over the step is below 1E-9 of
the cell's energy scale. At a full CFL step that is 128 of 500 cells at
checkpoint A, all of them below r = 1.042, and it is most of the energy row
at checkpoint M.

**The assembly selector is not a cause of anything here and cannot be**, on
any of the four routes, because `ieq_sweep_state_kind` never leaves
`ieq_state_marching` outside a solver pass.

---

## 11. Artifacts

All raw products are in the worker scratch and nothing was written into the
repository except this file.

```
/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/
  59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/p2/
    mkrun.sh                 builds one isolated run directory from a stored generation
    drive.sh                 one invocation, its sidecar, its exit class, its file inventory
    sed -> LHS1140b/sed      so the case's own relative spectrum path resolves unedited
    models/grp/<name>/       one directory per invocation, 24 of them
      SIDECAR.txt            the gate of rev9 section 4.2, one per invocation
      input.inp  run.log  run.err
      output/residual_profile.txt     routes 1 (r, n, v, T, R_mass, R_mom, R_energy, Phi_lo, Phi_hi)
      output/update_map.txt           route 4 (dt, j, k, r, R, G_dt, hydro, chem, energy, transport, filter;
                                      k = 4 rows are the carrier row and carry 9 columns)
      output/boundary_trace.txt       face-state gaps and the base mass rows of cells 1..12
      output/ghost_record.txt         route 3 only
      output/eval_state_dump.txt      route 3 only, 17 digits, blocks loaded / work_before_products /
                                      work_after_products / loaded_kept
      output/Hydro_ioniz.txt, Ion_species.txt, *_adv.txt, Heating_breakdown.txt,
      Cooling_breakdown.txt           routes 3 and 4
```

The 24 invocation names are those of the table in section 3.3. No background
process of this work is left running.
