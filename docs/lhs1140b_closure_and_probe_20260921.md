# LHS 1140 b: the eliminated-closure count against the probe arc

2026-09-21. Phase 4 of `PLAN_20260920_rev9.md` section 15, the experiment
section 8 asks for, on the two checkpoints whose Stage A0 identity records are
`docs/lhs1140b_identity_records_20260921.md` sections 3.1 and 3.4.

## 1. The verdict, first

**The directional derivative stabilizes in the probe arc and it does not move
with the closure count, at both checkpoints and along all four directions.**
In the terms section 8.2 sets:

- **Does it stabilize in `h`?** Yes. Over the arcs the plateau covers, the
  central difference of one and the same map changes by at most **1.2E-03**
  of itself on the atomic checkpoint and **1.9E-03** on the molecular one.
  The plateau is bounded below by rounding and above by curvature, both of
  which the grid shows on at least one direction of each checkpoint, and the
  production arc `h = 1` sits inside it on every direction of both.
- **Does it move with `k`?** No. At the production arc the action at a fixed
  closure count differs from the action at the run's own converged closure by
  at most **1.1E-04** (atomic) and **8.7E-04** (molecular) of itself, and on
  three of the four directions by less than **1E-06**. At `k = 1` the two are
  **bitwise the same vector** at both checkpoints, because the elimination
  reaches its fixed point in one pass at both.
- **Is the inner map resolved?** Yes, at these two states. The elimination
  exits after ONE pass, at an increment of 3.3E-15 (atomic) and 1.5E-13
  (molecular) against the 1.0E-08 it is asked for, with no cell left without a
  chemical root and no refusal at any of the 7 counts tried.

**So the closure-bias question is closed for these two states and Stage C can
proceed on the linear model alone.** Two qualifications belong with that, both
measured here and neither a closure effect:

1. `|| Dr^-1 (F - F_jac) ||` is exactly zero at `k = 1`, which is the count
   the solve forms at these iterates. At `k >= 2` it is 9E-17 of the row
   scales on the atomic checkpoint and 2E-11 to 4E-11 on the molecular one.
   In the units the CERTIFICATION reads the rows in, the molecular number is
   **0.86 to 1.29 against a residual of 1.90**: at a state this close to its
   gate, letting the elimination take one pass more moves the residual by half
   of what the residual is. It is the amplification of a 7E-13 composition
   change, it lands in the MASS row of cell 1, and it is a statement about how
   little that state's mass row survives its own cancellation, not about the
   Jacobian.
2. The FORWARD difference, which is what the production action takes, acquires
   a symmetric offset of fixed absolute size once the arc exceeds about three
   times its production value along some directions. The CENTRAL difference is
   immune to it, and at the production arc the offset is 2.3E-04 of the action
   (atomic, the smooth direction) and 1.5E-05 (molecular, the mass-row cell).
   This is a Stage B property of the residual along a direction, not a closure
   count effect, and its cause was not identified here.

## 2. Classification, and what was inherited

Every new number in this document is **Mode R** (rev9 section 3): the current
binary evaluating an immutable stored generation. No physical step was taken,
no state was published, nothing under `LHS1140b/models/` was written, and
nothing under `backup/regression/` was written. Each run was made on its own
copy of a stored generation in a scratch directory of its own.

Three results are INHERITED and were not re-measured; the experiment rests on
them and says so:

- **The residual is a state function at both checkpoints**, bitwise identical
  over 28 invocations, five thread configurations and two processes
  (`docs/lhs1140b_repeatability_20260921.md`). A directional difference
  measured here is therefore a statement about the operator and not about its
  noise, and the evaluation uncertainty section 6 item 2 asks to be measured
  independently is exactly zero.
- **The loaded-state route, the evaluate route and the stationary assembly of
  the marching route agree term for term**; the Newton route departs on a
  caloric mixture; the marching ENERGY row is a different map
  (`docs/lhs1140b_residual_routes_20260921.md`). This experiment runs on the
  Newton route's own `eval_residual`, which is the route whose derivative is
  in question.
- **`EXHALE_RESID_QUAD` cannot fire on a diagnostic route.** INSPECTED,
  `steady_residual.f90:246-248`: the assembly selector replaces
  `rows_kind = ROWS_PRODUCTION` only when the sweep state kind is not the
  marching one, and this route never calls
  `set_ioniz_eq_sweep_state_kind`. MEASURED here, printed by every
  invocation: `ieq_sweep_state_kind 1, is the marching kind: T`. So
  `rows_kind = ROWS_PRODUCTION` throughout, the mode 2 comparison of section
  13.6 is again NOT PERFORMED, and this was not chased further.

## 3. Identities

### 3.1 The binary

| | |
|---|---|
| path | `<scratch>/closure_probe/EXHALE_test.x` (private build; the tree's `EXHALE.x` was not rebuilt and not replaced) |
| md5 MEASURED | `55e3bf7bb89166072ff867fc7b2a24db` |
| built from | the working tree at repository HEAD `93eed8667483` (dirty), plus the diagnostic of section 4 |
| build command | `make OBJDIR=<scratch>/obj EXE=<scratch>/EXHALE_test.x -j8` |
| source digest | `41972a7390f86e5c658465d27c84d54d`, the md5 of the sorted md5 list of the 245 `.f90`/`.inc` files of `src/` at the time of the build (MEASURED). This is a live tree with concurrent workers, so it records what stood beside the binary, not a permanent property of it |
| the tree binary of the brief | `EXHALE.x` md5 `a55d90086af7a18699ed86d2d39b31b2` (MEASURED), unchanged and unused for these runs, since the grid needs the new entry point |

Compiler and linked libraries, MEASURED with `strings -a` and `ldd`, because
they set the floating-point and reduction behavior a difference quotient is
measured against:

| | |
|---|---|
| compiler strings | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`; `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` |
| BLAS and LAPACK | `/opt/miniconda3/lib/libopenblas.so.0` |
| OpenMP runtime | `/opt/miniconda3/lib/libgomp.so.1` |
| Fortran runtime | `/opt/miniconda3/lib/libgfortran.so.5`; `libquadmath.so.0` also linked |

### 3.2 The two checkpoints

| | `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `molecular_scalar_gj1132_kzz1e9/HeH0.083` |
|---|---|---|
| generation | `g0002_20260919T004843Z_f7485b14` | `g0004_20260920T053734Z_79b42a03` |
| Stage A0 record | section 3.1 | section 3.4 |
| `Hydro_ioniz.txt` md5 MEASURED | `25e11562ec3cd1f6d204bf36f8e503eb` | `fd882d8918052069f8ed229dd376f8c4` |
| `Ion_species.txt` md5 MEASURED | `78ac49d7e1ce90a3a4f271785400efda` | `87288b03ad88a92b7b8b8b6c2b1bdd66` |
| boundary model of the STORED state (READ) | `..._contact_upwind_v2` | `..._ghost_fixed_point_seed_reservoir_row_v3` |
| boundary model the CURRENT binary solves | `..._ghost_fixed_point_seed_reservoir_row_v3` | the same |
| relation | **Mode R, the current operator on an old state, boundary rebuilt** | **Mode R, the current operator on a state of its own boundary model** |
| model options (READ, state header) | `mol=F carrier=F carrier_newton=F he_diff=T sec_ion=T wellbal=T` | `mol=T molbase=T carrier=T carrier_newton=F he_diff=T sec_ion=T wellbal=T` |
| system the grid measures (MEASURED) | 3 unknowns per cell over 500 cells, 1500 entries, species rows 0 | the same |
| seed composition digest (MEASURED) | `288F904E7B774270` | `34C30AC3DDA15380` |
| state digest (MEASURED) | `CF430FA7D6AF677C` | `A52CCB800ACD9772` |

Every table in this document is Mode R and neither row is a reproduction of
its historical certificate.

**A cross-check the grid supplies for free.** On the atomic checkpoint the
largest mass row in the units the certification judges it in is MEASURED
`5.0488E+07` at cell 2. The Stage A0 table of the plan READS, for the same
generation, `mass 9.806E-01 against 1.9E-08 at cell 2 (distance 5.049E+07)`.
The two agree, so the scale this document divides by is the scale the
certificate was written on.

### 3.3 The effective configuration at every evaluation

MEASURED, printed by every invocation:

| | atomic | molecular |
|---|---|---|
| reconstruction | WENO3 | WENO3 |
| PLM to WENO3 continuation armed | F | F |
| judged row scaling (`EXHALE_JUDGED_ROWS`) | F | F |
| `ieq_sweep_state_kind` | 1, the marching kind | the same |
| active species bounds held | 0 | 0 |
| closure policy of the run | sweep to the fixed point, at most 25 passes, cap -1, increment asked for 1.0E-08 | the same |

The reconstruction MODE is frozen: `select_stationary_reconstruction`
(`stationary_operator.f90:80-101`) installs WENO3 with `use_plm` false and the
blend disabled, which is the endpoint condition section 13.6 requires and
which agrees with the `recon=WENO3` both stored headers carry. The WENO
smoothness weights, the contact-mode dissipation and the positivity repairs
are RECOMPUTED at every evaluation, because each evaluation runs the whole
flux assembly; nothing here freezes them. The active-set mask is empty and the
projection displacement is zero at both checkpoints, because the species-row
registry is empty (`carrier_newton=F`), so every difference below is a
difference of the unconstrained map and not of a projected one.

## 4. The control that was added

Section 8.1 asks for the closure count to be CONTROLLED, and states that
`EXHALE_RESID_SC_MAX` is a cap and not the control. Every statement rev9
section 8.1 makes about the control holds; its line numbers do not, because
`steady_newton.f90` has grown since rev9 was written: by other work above
line 6300, and by the 543 lines the diagnostic of this section adds at
`:6298`, which shift everything below it. The current ones, INSPECTED with
that diagnostic in the file, and the rev9 number beside each:

- `eval_residual` runs the elimination to a convergence test or, when the
  caller names `n_eq_sweeps_fixed`, for exactly that many passes and takes no
  test at all (`steady_newton.f90:4177-4186`, rev9 `:4175-4183`; and `:4251`).
- The Newton loop forms the base point its model is differenced about at
  `n_eq_sweeps_fixed = n_eq_sweeps_model` (`:17605`, rev9 `:17038`), and the
  matrix-free action samples its endpoint at the same count (`:11647`, inside
  `jv_product` at `:11491`, rev9 `:10964`).
- `n_eq_sweeps_model` starts at 1, is REASSIGNED from `n_eq_sweeps_last`
  (`:11221`, `:17291`; rev9 `:10653`, `:16723`) and raised with `max`
  (`:11320`, `:18235`; rev9 `:10752`, `:17667`). It is fixed within one
  Jacobian evaluation and varies along a run.
- The probe-arc controls are read at `:8227-8240` (rev9 `:7659-7674`) and the
  assembly selector stands at `steady_residual.f90:246-248` (rev9
  `:239-247`).
- **And the right-hand side is not the base point.** `pgmres` is handed
  `F_jac` as the point the action is differenced about and `-F/Drow` as the
  right-hand side (`:17978`), `F` being the residual at the run's own closure
  policy. `|| Dr^-1 (F - F_jac) ||` is therefore exactly the mismatch between
  the operator of the linear model and its right-hand side, which is why
  section 8.2 asks for it.

**What was added: one diagnostic entry point, off by default, that evaluates
and stops.** `EXHALE_CLOSURE_PROBE=1` (`EXHALE_main.f90:1536-1572`) installs
the stationary reconstruction, registers the system a solve of this
configuration would carry, reads the composition-elimination controls, calls
`jacobian_action_against_closure_count_and_probe_arc`
(`steady_newton.f90:6352-6836`), with two routines of its own beside it,
`digest_of_a_real_array` (`:6300-6315`), which says whether two evaluations
were handed the same numbers, and `scaled_row_class_norms` (`:6319-6348`),
which is the row distribution every norm in this document is reported as, and
stops the run. It reaches no solve, adopts
no state and writes no file. The count it imposes is the existing argument
`n_eq_sweeps_fixed` of `eval_residual`, not a new control; `n_eq_sweeps_model`
is set inside the routine only so that `jv_product` samples the same map, and
the run ends before any solve could read it. The key is documented at its
reader, in the comment style the neighboring residual hooks use.

**The regression that says nothing moves with the key unset.** Two fast cases,
run on scratch copies of the case directories and their references, so that
nothing under `backup/regression/` was read for writing or written at all:

```
REGRESSION_EXE=<scratch>/EXHALE_test.x ./run_check.sh check hp_front roundtrip
[build] skipped: REGRESSION_EXE=... (55e3bf7bb891)
[hp_front] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=100) ...
            final: count=100  du= 6.6225E-01  dtu= 3.9604E-03
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
[roundtrip] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=40) ...
            final: count=40  du= 2.5794E+00  dtu= 6.4339E-04
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
==> REGRESSION PASS (all cases byte-identical)
```

MEASURED, with `EXHALE_CLOSURE_PROBE` absent from the environment (`env`
carried no `EXHALE_*` variable at all). Eight files, all byte-identical to the
stored references. The logs are in
`docs/audit_20260905/closure_and_probe_20260921/regression/`. No golden was
written, and `golden` mode was never invoked.

## 5. What the grid measures

At one frozen state, two dimensions:

- **`k`, the eliminated-closure sweep count**: 1, 2, 3, 4, 6, 8, and `k = 0`
  standing for the run's own policy, sweep to the convergence test. The
  converged row is taken FIRST so that it is the reference the fixed counts
  are read against.
- **`h`, the probe arc**: 1, 1E-03, 1E-02, 1E-01, 3E-01, 3, 10, 100, 1000
  times the production arc, imposed on the module variable
  `EXHALE_JV_PROBE_ARC` sets (`steady_newton.f90:8230-8240`), so that the
  whole sweep runs in ONE process on one loaded state with one identity.
  `EXHALE_JV_COLUMN_SCALE` was left off, so the arc is the unscaled rule
  `sqrt(epsilon)(1 + ||Y||)/||v||`.

63 pairs `(k, h)`, four directions each, so 252 actions per checkpoint. Each
action costs three residual evaluations, the production `jv_product` and the
two endpoints `Y +/- eps v` at the SAME `k` and from the SAME seed, and each
`k` costs one base residual besides. 763 evaluations per checkpoint, plus the
18 the banded model of the first direction takes.

**The four directions, built once and held over the whole grid.** All four are
unit vectors in the column-scaled coordinates, so their actions are read on
one scale.

| name | what it is | atomic | molecular |
|---|---|---|---|
| `banded_step` | the step the solve proposes at this state: the Newton direction of the banded model of the FULL residual in the scaled coordinates the linear system is written in, `(J D / Dr) z = -F/Dr`, `s = D z`. It is the step the outer iteration takes when the Krylov cycle accepts its own preconditioner's proposal with no pseudo-transient term. It is NOT the Krylov step, and it cannot be: the Krylov step is itself a function of `h` and `k`, so it is not a direction a sweep in `h` and `k` can be run along | `||D^-1 s|| = 5.5078E-02`, 0 colors unresolved, LAPACK info 0 | `9.2545E-10`, 0 colors unresolved, info 0 |
| `cell_mass` | supported on the three unknowns of the cell whose MASS row stands furthest outside what the certification allows | cell 2 | cell 1 |
| `cell_energy` | the same at the cell whose ENERGY row stands furthest outside | cell 1 | cell 3 |
| `smooth` | `sin(0.1 j)` on every unknown, scaled by the column scales: the deterministic direction the run-level self-tests use, which reaches every cell | 1500 entries | 1500 entries |

**A self-consistency check the first direction supplies.** The action along
`banded_step` should return `-F/Dr` divided by `||D^-1 s||`, since that
direction solves the banded model. MEASURED on the atomic checkpoint,
`||Dr^-1 J v||_2 = 3.0661E-02` against `||Dr^-1 F||_2 / ||D^-1 s|| =
1.6899E-03 / 5.5078E-02 = 3.068E-02`: they agree to the three digits the
banded model is accurate to.

**And the forward difference reported here IS the production action.**
MEASURED: at all 252 grid points of the atomic checkpoint the forward
difference computed here and the vector `jv_product` returned print the same
five class norms, digit for digit.

## 6. The gate of section 4.2

Every invocation ran in a directory of its own, created fresh, holding its own
copy of the input file, its own copy of the stored state pair as the restart
pair, its own copy of the executable and its own empty `output/`. The complete
environment was imposed with `env -i`, so the set recorded IS the set in
force. The driver is
`docs/audit_20260905/closure_and_probe_20260921/run_probe.sh`; the sidecars,
the diagnostic blocks, the reduced tables, the two regression logs and
`MD5SUMS.txt` over all of them are beside it, and the reduction is
`summarize.py` in the same directory.

| invocation | case | threads (OMP, BLAS) | wall s | exit | evaluation stage reached | block lines |
|---|---|---|---|---|---|---|
| `A01_atomic_p1` | atomic | 8, 1 | 107.91 | 0 | the grid completed | 2077 |
| `A02_atomic_p2` | atomic | 8, 1 | 106.21 | 0 | the grid completed | 2077 |
| `A03_atomic_t1` | atomic | 1, 1 | 793.34 | 0 | the grid completed | 2077 |
| `M01_molecular_p1` | molecular | 8, 1 | 223.04 | 0 | the grid completed | 2077 |
| `M02_molecular_p2` | molecular | 8, 1 | 221.60 | 0 | the grid completed | 2077 |

There is no "no output", no timeout and no unresolved measurement among them.
Each left a nonempty standard output containing the whole block, ending in
`---- end of the grid ----`, and the restart pair each was given is
byte-identical to the stored generation after the run. Each `stderr` holds one
line, the gfortran note `IEEE_UNDERFLOW_FLAG IEEE_DENORMAL`, which the runtime
prints at every `STOP`.

**The diagnostic is itself reproducible across processes AND across thread
counts.** MEASURED: the whole 2077-line block of `A01`, `A02` and `A03` has
one md5, `80fcee034a9b0be60f69733c5da78ecf`, although `A03` ran on one OpenMP
thread and took 793 s against the 108 s of the eight-thread runs, which is the
evidence that the threads were used; of `M01` and `M02`, one md5,
`5b48da9932d2e145ead921d7597ed0b7`. Every number in every table below is
therefore the same number in two independent processes, and on the atomic
checkpoint in two thread configurations as well.

## 7. The closure count at the state

`|| Dr^-1 (F - F_jac) ||`, `F` the residual at the run's own closure policy
and `F_jac` the residual at the fixed count, with its row distribution. `Dr`
is the scaling the linear model is divided by, which on this three-unknown
route is the state column scales (`row_scaling_of_the_linear_model`,
`steady_newton.f90:2783-2808`, `judged_row_scaling_on` false). `Dj` is the
size each hydrodynamic row is JUDGED at, which is the scale the certification
reads, and it is given beside it because a row-scale norm and a judged-scale
norm answer different questions.

### 7.1 `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7`

`||Dr^-1 F||_2` by class, MEASURED: mass 4.7761E-04, momentum 1.4565E-03,
energy 7.1150E-04, all 1.6899E-03. `||Dj^-1 F||_2`: mass 5.8014E+07, momentum
9.1260E+02, energy 1.4688E+06, all 5.8032E+07.

| k | passes completed | increment the sweep exited on | admissible | cells without a chemical root | refusal screen | closure increment against the converged closure | channel | `||Dr^-1 (F-F_jac)||_2` mass / momentum / energy / all | `||Dj^-1 (F-F_jac)||_2` all |
|---|---|---|---|---|---|---|---|---|---|
| converged | 1 | 3.295E-15 | T | 0 | none | 0 | (none moved) | 0 / 0 / 0 / 0 | 0 |
| 1 | 1 | not taken at a fixed count | T | 0 | none | 0 | (none moved) | 0 / 0 / 0 / 0 | 0 |
| 2 | 2 | not taken at a fixed count | T | 0 | none | 9.545E-14 | heat - cool | 0 / 0 / 8.9635E-17 / 8.9635E-17 | 1.0295E-07 |
| 3 | 3 | not taken at a fixed count | T | 0 | none | 6.482E-16 | heat - cool | 0 / 0 / 5.5796E-19 / 5.5796E-19 | 3.6339E-09 |
| 4 | 4 | not taken at a fixed count | T | 0 | none | 9.545E-14 | heat - cool | 0 / 0 / 8.9635E-17 / 8.9635E-17 | 1.0296E-07 |
| 6 | 6 | not taken at a fixed count | T | 0 | none | 9.545E-14 | heat - cool | 0 / 0 / 8.9635E-17 / 8.9635E-17 | 1.0296E-07 |
| 8 | 8 | not taken at a fixed count | T | 0 | none | 9.545E-14 | heat - cool | 0 / 0 / 8.9635E-17 / 8.9635E-17 | 1.0296E-07 |

`k = 1` reproduces the converged-closure residual **bitwise**: the mass,
momentum and energy row differences are exactly zero, not small. Beyond one
pass the elimination moves only the radiative pair of cell 486 by 9.5E-14, and
that reaches the energy rows at 1.0E-07 of the tolerance they are judged at.
**The closure bias at this state is zero in the count the solve uses and
1.0E-07 of a tolerance at any other count.**

### 7.2 `molecular_scalar_gj1132_kzz1e9/HeH0.083`

`||Dr^-1 F||_2` by class, MEASURED: mass 3.7956E-11, momentum 3.7985E-11,
energy 4.2935E-11, all 6.8753E-11. `||Dj^-1 F||_2`: mass 1.9018E+00, momentum
7.1791E-05, energy 2.4573E-02, all 1.9019E+00. The largest judged row is the
mass row of cell 1 at 5.3329E-01.

| k | passes completed | increment the sweep exited on | admissible | cells without a chemical root | refusal screen | closure increment against the converged closure | channel | `||Dr^-1 (F-F_jac)||_2` mass / momentum / energy / all | `||Dj^-1 (F-F_jac)||_2` all |
|---|---|---|---|---|---|---|---|---|---|
| converged | 1 | 1.459E-13 | T | 0 | none | 0 | (none moved) | 0 / 0 / 0 / 0 | 0 |
| 1 | 1 | not taken at a fixed count | T | 0 | none | 0 | (none moved) | 0 / 0 / 0 / 0 | 0 |
| 2 | 2 | not taken at a fixed count | T | 0 | none | 7.150E-13 | heat - cool | 2.7841E-11 / 5.4386E-12 / 1.9840E-11 / 3.4617E-11 | 9.3196E-01 |
| 3 | 3 | not taken at a fixed count | T | 0 | none | 5.931E-14 | heat - cool | 1.3209E-11 / 6.8192E-12 / 1.8075E-11 / 2.3403E-11 | 8.5753E-01 |
| 4 | 4 | not taken at a fixed count | T | 0 | none | 7.150E-13 | heat - cool | 2.1629E-11 / 8.2686E-12 / 1.9128E-11 / 3.0034E-11 | 9.4982E-01 |
| 6 | 6 | not taken at a fixed count | T | 0 | none | 7.150E-13 | heat - cool | 3.1928E-11 / 1.3438E-11 / 2.3111E-11 / 4.1642E-11 | 1.2105E+00 |
| 8 | 8 | not taken at a fixed count | T | 0 | none | 7.150E-13 | heat - cool | 3.1654E-11 / 1.4445E-11 / 2.4307E-11 / 4.2443E-11 | 1.2890E+00 |

Again `k = 1` is bitwise the converged residual. Beyond one pass the story is
different from the atomic one and it is worth stating plainly:

- The composition change is 7.2E-13, in the radiative pair of cell 499.
- It reaches the rows at 3E-11 of the column scales, which is **half of the
  whole residual**, whose own norm is 6.9E-11.
- In the units the certification judges the rows in it is **0.86 to 1.29
  against a residual of 1.90**, and 0.93 to 1.29 of it sits in the MASS rows,
  worst at cell 1.
- The sequence over `k` is 2.78, 1.32, 2.16, 3.19, 3.17 (units 1E-11 of the
  row scales at `k` = 2, 3, 4, 6, 8). It does not converge; it oscillates at
  that level. The elimination does not settle to one composition beyond its
  first pass, it cycles at 1E-13.

**How to read that.** The amplification from the composition increment to the
row is about 35 (7.2E-13 of the radiative pair to 2.5E-11 of the mass row
scale), which is the map `measure_the_closure_map_at_the_base` exists to
measure and which is not available on this route because that routine returns
at once when no species row is registered (`steady_newton.f90:6232`). It is
NOT evidence that the Newton model at this state is built on the wrong
derivative: the count the solve actually forms here is 1, at which the two
residuals are the same bits. It IS evidence that the mass row of the base cell
of this state carries so little above its own cancellation that a 1E-13
composition wobble is half of it, which is the finding to carry into the
conservation budgets of phase 5.

## 8. The two-dimensional table

Each entry of the first table of a pair is `||Dr^-1 J v||_2` over the whole
vector, the CENTRAL difference at that `(k, h)`; the second is the distance
from the same `k` at the production arc, relative; the third is the distance
from the converged closure at the same arc, relative. The full tables, with
the row distribution by class and the endpoint record of every point, are
`docs/audit_20260905/closure_and_probe_20260921/tables/`.

### 8.1 Atomic, the plateau in `h`

Distance from the action at the production arc, relative, at `k = 1`
(MEASURED; the other `k` agree to the digits shown):

| direction | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| `banded_step` | 3.63E-01 | 3.50E-02 | 3.91E-03 | 1.19E-03 | 0 | 3.96E-04 | 3.79E-04 | 3.84E-04 | 3.85E-04 |
| `cell_mass` | 4.94E-07 | 1.15E-07 | 2.39E-09 | 1.82E-09 | 0 | 1.55E-08 | 1.92E-07 | 1.86E-05 | 1.69E-03 |
| `cell_energy` | 6.66E-07 | 1.46E-07 | 2.78E-08 | 1.11E-08 | 0 | 1.46E-08 | 1.67E-07 | 1.68E-05 | 1.47E-03 |
| `smooth` | 2.28E-04 | 2.28E-04 | 2.28E-04 | 1.04E-04 | 0 | 1.61E-04 | 2.08E-04 | 6.95E-04 | 8.28E-04 |

The two localized directions show the textbook shape with the rounding floor
below the resolution of the table: flat at 1E-09 to 1E-07 from `h = 1E-03` to
`h = 10`, then curvature taking over, 1.7E-03 at a thousand times the
production arc. `banded_step` is the direction along which the action is
smallest (3.07E-02, against 4.48E+03 for `cell_mass`), so its relative measure
is the most sensitive one; it is rounding-limited below `h = 0.1`, falling
exactly as `1/h`, and flat at 4E-04 from `h = 0.3` upward. `smooth` is flat at
1E-04 to 8E-04 throughout.

**The plateau the four directions share is `h = 0.3` to `h = 10`: over it the
action of one and the same map changes by at most 1.2E-03 of itself on the
worst direction, `banded_step`, and by at most 1.9E-07 on the two localized
ones.** The production arc is inside it.

### 8.2 Atomic, the movement with `k`

Distance from the action at the converged closure, same arc, relative
(MEASURED):

| direction \ h | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| `banded_step`, k=1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.11E-04 | 1.11E-04 |
| `banded_step`, k=8 | 1.17E-04 | 1.11E-04 | 2.14E-04 | 1.27E-04 | 1.11E-04 | 1.11E-04 | 1.11E-04 | 3.62E-08 | 3.62E-08 |
| `cell_mass`, k=8 | 2.93E-13 | 1.62E-13 | 1.56E-13 | 1.57E-13 | 6.38E-10 | 2.13E-10 | 1.57E-13 | 1.57E-13 | 1.57E-13 |
| `cell_energy`, k=8 | 1.44E-13 | 1.55E-07 | 1.47E-13 | 2.58E-09 | 7.75E-10 | 1.47E-13 | 7.76E-11 | 7.86E-12 | 5.60E-12 |
| `smooth`, k=8 | 4.87E-09 | 4.78E-09 | 4.78E-09 | 4.90E-09 | 4.78E-09 | 4.78E-09 | 4.78E-09 | 1.95E-12 | 1.95E-12 |

The largest movement with `k` anywhere in the atomic grid is 1.7E-02, at
`banded_step`, `k = 2`, `h = 1E-03`, which is the point at which the quotient
is rounding and not a derivative at all; at every `h` in the plateau the
movement is at most **1.1E-04**, and on the three directions other than
`banded_step` at most **5E-09**.

### 8.3 Molecular, the plateau in `h`

Distance from the action at the production arc, relative, at `k = 1`
(MEASURED):

| direction | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| `banded_step` | 2.89E-01 | 3.04E-02 | 3.77E-03 | 1.07E-03 | 0 | 4.33E-04 | 3.77E-04 | 3.74E-04 | 5.24E-04 |
| `cell_mass` | 1.27E-05 | 1.66E-07 | 1.40E-07 | 4.67E-08 | 0 | 1.77E-01 | 2.64E-01 | 2.97E-01 | 3.10E-01 |
| `cell_energy` | 1.75E-06 | 1.63E-07 | 5.98E-08 | 3.44E-08 | 0 | 2.85E-07 | 3.61E-04 | 9.51E-04 | 7.63E-02 |
| `smooth` | 2.63E-05 | 2.90E-06 | 2.46E-07 | 9.25E-08 | 0 | 4.58E-08 | 4.93E-08 | 9.32E-03 | 8.83E-03 |

Here the plateau is bounded on BOTH sides and the measurement says where. The
`cell_mass` direction is flat to 1.7E-07 from `h = 1E-02` to `h = 1` and then
leaves it hard: the action falls from 1.315E+04 to 1.154E+04 between `h = 1`
and `h = 1000`, a 12 per cent change, which is real curvature of the residual
along that direction and not a closure or a rounding effect. `banded_step` is
rounding-limited below `h = 0.1`, exactly as on the atomic checkpoint, and
flat at 4E-04 above `h = 0.3`.

**The plateau the four directions share is `h = 0.3` to `h = 1`: over it the
action of one and the same map changes by at most 1.9E-03 of itself on the
worst direction, `banded_step`, and by at most 1E-07 on the two localized
ones.** `cell_mass` and `cell_energy` hold their plateau down to `h = 1E-02`
as well; `banded_step` does not, and that is what bounds the shared range from
below. The production arc is at the upper edge of the shared plateau, which is
where the arc rule is designed to put it.

### 8.4 Molecular, the movement with `k`

Distance from the action at the converged closure, same arc, relative
(MEASURED, at `k = 8`; the other counts agree to the digits shown):

| direction \ h | 1E-03 | 1E-02 | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| `banded_step` | 4.08E-01 | 3.22E-02 | 3.46E-03 | 1.36E-03 | 7.54E-04 | 5.83E-04 | 5.60E-04 | 5.58E-04 | 2.70E-07 |
| `cell_mass` | 1.20E-05 | 1.82E-07 | 8.00E-08 | 2.36E-08 | 1.97E-08 | 4.65E-10 | 5.22E-10 | 8.92E-11 | 7.06E-11 |
| `cell_energy` | 2.76E-06 | 1.68E-07 | 6.82E-08 | 1.27E-08 | 8.98E-09 | 2.57E-09 | 6.94E-10 | 1.17E-10 | 6.27E-11 |
| `smooth` | 2.44E-05 | 2.56E-06 | 7.26E-07 | 6.83E-07 | 6.76E-07 | 6.76E-07 | 3.49E-09 | 1.43E-09 | 8.99E-11 |

Inside the shared plateau the movement with `k` is at most **8.7E-04**
(`banded_step`, `k = 6`, `h = 1`) and at most **7.1E-07** on the other three
directions. Outside it, at `h = 1E-03`, `banded_step` moves by up to 0.75, and
that number says nothing about the closure: at that arc the quotient is
rounding, and the `k` rows and the `h` rows carry the same 1/h growth.

## 9. Admissibility, and the endpoint record

Section 8.1 asks whether `F` is admissible for derivative use at every count,
and section 8.2 rules out any count at which the inner map leaves roots
unresolved or changes branch.

**Every evaluation of both grids was admissible.** MEASURED over all
2 x 252 actions and their 2 x 7 base evaluations:

- the probe, the forward endpoint and the backward endpoint all returned
  `admissible` at every point (`ok=TTT`, 252 of 252 on each checkpoint);
- no cell was left without a chemical root at any endpoint, at any `k`, at any
  arc (`roots 0 0` everywhere);
- no refusal screen fired (`refusal 0` at every base evaluation);
- no component of any direction was blocked and no probe had to be taken on
  the backward side (`blocked 0`, `backward F` everywhere), which is what an
  empty active set means;
- the arc requested and the arc used are the same number at every point: no
  step was shortened by the species box or by the element rows.

**There is therefore no `k` at which the map was inadmissible for a
derivative, at either checkpoint.** That is a statement about these two
states; it is not a statement about the states of the plan's section 9 whose
chemistry rows refuse.

**The one endpoint asymmetry, and it is in the converged row.** At a fixed
`k` both endpoints completed exactly `k` passes at every point of both grids,
so both sides of every fixed-count difference are the same map. At `k = 0`,
the run's own policy, the base completed 1 pass while the endpoints completed
3 or 4 at the longest arcs:

| checkpoint | where | passes at the base | passes at the two endpoints |
|---|---|---|---|
| atomic | `banded_step` and `smooth`, `h = 100` and `h = 1000` | 1 | 3, 3 |
| molecular | `smooth`, `h = 10` and `h = 100`; `banded_step`, `h = 1000` | 1 | 3, 3 |
| molecular | `smooth`, `h = 1000` | 1 | 4, 4 |

At those eight points the converged-closure forward difference is a difference
of two maps with different pass counts, which is exactly what section 8.1
warns of. They are all outside the plateau, the central difference at them
stands at most 9.3E-03 from the same direction's action at the production arc,
and no conclusion here rests on them. Every point of the plateau has 1 pass at
the base and 1 at both endpoints.

## 10. One thing the grid found that is not a closure effect

The FORWARD difference, which is the production action, acquires a symmetric
offset of FIXED ABSOLUTE SIZE once the arc exceeds a direction-dependent
threshold; the CENTRAL difference is exactly immune to it. MEASURED on the
atomic checkpoint along `smooth` at `k = 1`, `||Dr^-1 (J_forward -
J_central)||_2` over the whole vector:

| h | 1E-01 | 3E-01 | 1 | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|
| `||Dr^-1 (Jf - Jc)||_2` | 6.92E-05 | 1.27E-01 | 2.30E-01 | 2.91E+03 | 8.73E+02 | 8.73E+01 | 8.77E+00 |
| times `2 eps`, the second difference | 7.8E-12 | 4.3E-08 | 2.6E-07 | 9.87E-03 | 9.87E-03 | 9.87E-03 | 9.91E-03 |
| `||Dr^-1 J_central||_2` | 1.0122E+03 | 1.0122E+03 | 1.0122E+03 | 1.0122E+03 | 1.0122E+03 | 1.0122E+03 | 1.0122E+03 |

The second difference `F(Y + eps v) - 2 F(Y) + F(Y - eps v)` is **constant at
9.87E-03 of the row scales over the three decades from `h = 3` to
`h = 1000`**, which is the signature
of an offset that both endpoints acquire and the base does not, and not of
curvature (which would fall as the square of the arc) nor of rounding (which
would fluctuate). At the production arc it sits almost entirely in the MASS
rows (2.2999E-01 of a total 2.3002E-01); at `h = 3` and above the MOMENTUM
rows carry most of it (2.5140E+03 of 2.9109E+03 at `h = 3`). The central
difference cancels it exactly, and stands at 1.0122E+03 to five digits at
every arc from 1E-03 to 1000.

At the production arc the offset is 2.30E-01 against an action of 1.0122E+03,
that is **2.3E-04**. On the molecular checkpoint the same `1/h` shape appears
along `smooth` at `h >= 100` (3.947 at 100, 0.3947 at 1000: a factor ten down
for a factor ten of arc); along `banded_step` there the forward-to-central
difference GROWS in proportion to the arc instead, 1.16E-02 at `h = 10` to
1.16E+00 at `h = 1000`, which is the ordinary first-order truncation of a
forward difference and not this offset. At the production arc the
forward-to-central difference is 1.5E-05 (`cell_mass`), 1.0E-06
(`cell_energy`), 8.1E-08 (`smooth`) and 1.4E-03 (`banded_step`).

**What it is was not established here.** It is symmetric in the sign of the
step, so it is a property of the size of the displacement and not of its
direction; it is concentrated in the mass rows; and it appears between arcs of
5.6E-08 and 1.7E-07 in the column-scaled coordinates. A limiter, a positivity
repair or a smoothness-weight switch that fires on any displacement of that
size would produce it, and so would a mass row whose terms cancel to that
level; nothing here distinguishes them. It belongs to Stage B, section 6 item
7, and it is recorded rather than chased because it does not touch the
question this phase was asked.

## 11. What was not excluded

- **Two states, not a class.** Both checkpoints have their elimination
  converging in ONE pass, so the `k` dimension was exercised on a map that is
  already at its fixed point at `k = 1`. A state whose elimination needs
  several passes would be a different and stronger test, and neither of these
  two checkpoints is one. This is the clearest limit of the result.
- **The species rows were never in the vector.** Both checkpoints carry
  `carrier_newton=F`, so the registry is empty, the vector is 1500 entries of
  three hydrodynamic unknowns, and the active set is empty. The closure
  question on a COUPLED route, where a probe can be blocked by the species box
  and the endpoints can be projected differently, is untouched; that is
  phase 6b.
- **The direction the solve takes is the banded model's proposal, not the
  Krylov step.** The reason is stated in section 5 and it is a property of the
  experiment and not a shortcut: the Krylov step depends on `h` and `k`, so
  holding it fixed across a sweep in `h` and `k` is not possible.
- **`EXHALE_JV_COLUMN_SCALE` was left off**, so only the unscaled arc rule was
  swept. The column-scaled rule is a different probe and would need its own
  grid.
- **Two thread configurations on the atomic checkpoint, one on the molecular
  one.** The whole atomic grid is byte-identical at 1 and at 8 OpenMP
  threads; the molecular grid was run only at 8, and its thread independence
  is inherited from the repeatability record rather than re-measured over the
  grid.
- Processor affinity, other hosts, other library versions.

## 12. The closing verdict

In the words section 8.2 sets:

**The derivative stabilizes in `h`**: over the plateau each action changes by
at most 1.2E-03 (atomic) and 1.9E-03 (molecular) of itself, and by less than
2E-07 on the directions localized at the refusing cells. The plateau is
bounded below by rounding and above by curvature, both of which the grid shows
on at least one direction of each checkpoint, and the production arc is inside
it at both checkpoints and along every direction.

**The derivative does not move with `k`**: at the production arc, at most
1.1E-04 (atomic) and 8.7E-04 (molecular) of itself, and at `k = 1`, the count
the solve forms at these iterates, it is the same vector bit for bit.

**The inner map is resolved at these two states**: one pass reaches the fixed
point, at an increment five (molecular) to seven (atomic) decades below what
it is asked for, with
no unresolved chemical root and no refusal at any count tried.

**Therefore the closure-bias question is closed for these two checkpoints and
Stage C can proceed on the linear model alone.** What remains open is not the
derivative but two things the grid measured on the way: on the molecular
checkpoint the elimination does not settle past its first pass but cycles at
1E-13, and that cycle is amplified about 35 times into the mass row of the base
cell, where it is half of the residual in the units the gate reads; and the
production forward difference carries a symmetric offset of fixed size above a
direction-dependent arc, 2.3E-04 of the action at the production arc, which
the central difference does not.
