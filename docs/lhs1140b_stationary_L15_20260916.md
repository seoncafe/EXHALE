# L15: where the eight-thread JFNK solve stops being reproducible

Item L15 of `docs/PLAN_20260913_lhs_stationary.md`. First step of the
investigation, 2026-09-16.

## 1. Verdict

**The residual evaluation is not where the thread dependence enters.** The
stationary residual, its certification rows and the state files it writes are
**bitwise identical** over five repeats at one thread and five at eight, on a
frozen state -- no solve, no update, no accepted step. What L4h measured (two
eight-thread runs of the same fixture parting after about 21 JFNK iterations)
therefore enters at the KRYLOV stage or later: the matrix-vector products, the
norms and inner products of the linear cycle, the banded preconditioner solve,
or the line search that reads them.

**This is not called harmless rounding, and it must not be until the shared
writes are excluded.** A difference that appears only at eight threads is a
race until proved otherwise; the order of a reduction is one explanation of
such a difference and a scratch array written by more than one thread is
another, and the two are told apart by reading the parallel regions, not by
the size of the difference. The next step of this item is to enumerate the
parallel regions of the JFNK path and, for each, the arrays written inside it,
the variables carried across iterations and the reductions formed -- and only
then to attribute.

## 2. What was measured

The frozen state is the certified fiducial of item L21,
`atomic_scalar_gj1132_kzz1e9/HeH2.13`
(`LHS1140b/models/.L21/resolve_new/output`), handed back to the binary as its
own IC pair with `Restart intent: stationary evaluate`, which evaluates the
residual and certifies without advancing anything. Ten runs, five at
`OMP_NUM_THREADS=1` and five at `OMP_NUM_THREADS=8`, each in its own
directory: `LHS1140b/models/.L21/L15/t{1,8}_r{1..5}`. Binary
`EXHALE_L21b.x`, md5 `6463f0fdf181bd1859d59da3226e7819` (the delivered
build adds one comment block and is `142ed07857266ea832f9389ecc210d59`; the
fiducial re-evaluates bit for bit identically under it).

| what | result |
|---|---|
| `output/Hydro_ioniz.txt`, data lines, md5 | `8f15a045baa3` in all ten |
| `output/Ion_species.txt`, data lines, md5 | `8a427fb6e00f` in all ten |
| hydrodynamic mass row | 2.175E-09 (cell 1) in all ten |
| hydrodynamic momentum row | 5.417E-13 (cell 500) in all ten |
| hydrodynamic energy row | 2.126E-08 (cell 1) in all ten |
| `\|\|R\|\| = max_k` | 2.1257E-08 in all ten |
| the whole `run.log`, timing and thread lines removed | one md5 within each thread count |

`Ion_species.txt` carries the composition and `Hydro_ioniz.txt` the heating
and the cooling, so the ionization sweep, the radiative rates and the
hydrodynamic row assembly are all inside this comparison.

**The one difference between the two thread counts is in neither.** The
1-thread and 8-thread logs differ in exactly one line, the floating-point
exception summary the runtime prints at exit:

```
< Note: The following floating-point exceptions are signalling: IEEE_DENORMAL
> Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
```

That summary is the union of the flags raised on the threads that ran, so more
threads can only add to it; it says that some cell handled by a worker thread
underflowed, which single-threaded the master thread reached too but without
setting the flag on the path the summary reads. No physical array differs.

## 3. What this does NOT show

- It does not show that the residual assembly has no shared write. A race that
  is not exercised by this state, this grid and this thread count is not
  excluded by ten runs; it shows only that on the path this evaluation takes,
  the answer does not depend on the thread count.
- It does not show that the JFNK path's difference is a reduction order. That
  is the hypothesis L4h recorded, and it is still a hypothesis.
- It says nothing about the molecular carrier path, which this state does not
  run.

## 4. Reproduce

```bash
cd /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L21/L15/t8_r1
OMP_NUM_THREADS=8 ../../../../../EXHALE_L21b.x > run.log 2>&1
grep -v '^#' output/Hydro_ioniz.txt | md5sum
```

The `input.inp` of each directory is the fiducial's with
`Restart intent: stationary evaluate`, and `output/` starts holding only the
`_IC` pair copied from `LHS1140b/models/.L21/resolve_new/output`.

## 5. Next

Enumerate the parallel regions on the JFNK path and what each writes:
`src/modules/time_step/steady_newton.f90` (owned by items L4h and L17 --
coordinate before editing), the residual assembly it calls, the
preconditioner's banded factor and solve, and the norms and inner products of
the Krylov cycle. For each region: arrays written, variables carried between
iterations, reductions formed and their order. Attribute only after that list
exists.

---

# Stage 2: the solve path, traced along the failing trajectory

Section 2 of `docs/PLAN_20260916_rev3.md`, 2026-09-16. Stage 1 above stands
unchanged; this section is the step it named, taken with instrumentation that
records the solve path itself rather than the state files it ends on. The
record of every identity and every raw comparison is
`docs/audit_20260905/L15_stage2_20260916/`.

## 1. Verdict

**The thread dependence is not localized, because on this tree it does not
occur.** The whole JFNK trajectory of the element reload is bitwise identical
at 1, 8 and 16 OpenMP threads: 105 outer iterations, 152 products of the
Jacobian action, and for each product the preconditioner input, the
preconditioner output, the action of the operator and the Arnoldi column, all
recorded as hashes of the bit patterns and all equal. The four trace files of
the three thread counts have ONE md5,
`bd12e7c8ad917beda538d789080e10e3` (MEASURED).

Three 8-thread runs of one binary, and an 8-thread run of a second binary
built from the entry text of `steady_newton.f90`, produce run logs that differ
from each other in exactly one field, the wall-clock seconds of the pass line.
The same holds in the configuration of the L4h memo section 5.0
(`EXHALE_PTC_RAMP_GUARD=0`, `EXHALE_JFNK_MAXIT=80`, 8 threads, two builds),
which is the comparison that recorded the divergence: the two builds now agree
in every printed digit of all 80 iterations.

What this does not say is that the L4h observation was wrong. The trajectory
itself has moved: L4h ended pass 1 of this fixture at a mass row of 3.08e-02
(guard off) and 1.65e-02 (guard on), where the same recipe on this tree ends
it at 6.11e-10 and 4.57e-10 (MEASURED, below). A divergence measured on a
trajectory that ran along the stalled ramp is not excluded by agreement on a
trajectory that no longer goes there. What IS excluded, for this fixture and
this path, is a shared scratch array written by more than one thread, a
mutable cache read across threads, and a reduction whose result follows the
thread count: each of the three would have to show at one of the 923 recorded
quantities, and none does.

## 2. What was instrumented

`src/modules/time_step/steady_newton.f90`, three environment keys, all off
unless set, and with none of them set no statement of the block is entered:

- `EXHALE_L15_TRACE=<path>`: one record per outer iteration carrying FNV-1a 64
  hashes of the bit patterns of the unknowns, the composition, the residual,
  the model base point, the column scales, the row scales, the active bound
  set and the pseudo-time, and one record per Krylov product carrying the same
  hash of the right-hand side, of the preconditioner input and output, of the
  action of the operator and of the Arnoldi column. The record is a hash of
  the bits and never a floating sum, so it cannot itself carry an order of
  addition, and every write happens after the quantity is complete and outside
  every loop and reduction of the solve.
- `EXHALE_L15_DUMP_AT=<k>`: at outer iteration k, those arrays themselves in
  stream access, with the faces of the species box.
- `EXHALE_L15_REPEAT=<k>`: at outer iteration k, the model base residual and
  one action of the Jacobian evaluated again, twice, inside the same process
  at the same state, with the evaluation products held and put back around
  each repeat.

## 3. What was measured

`backup/regression/atomic_elem_newton` reloaded from `IC/` per that case's
README (`Load IC? True`, `Coupled carrier solve: False`, `Restart intent:
stationary`, which fixes the pseudo-time start at 1.0), `EXHALE_OUTER_PASSES=1
EXHALE_DIFF_OMEGA=0.5`, on scratch copies, with the private build
`EXHALE_L15.x` md5 `f23da8d84fafc2f0587b610add56e714` and its control
`EXHALE_L15ctl.x` md5 `7fece24a7561d6413fb30d07ba9e90c8`, which is the same
tree with the entry text of the one file this item edits. The delivered source
differs from the built one in three lines of one comment and its build, md5
`1e8afdcdf6b2a84b723562ecdce4b5a6`, reproduces the repeat run's trace and its
checkpoint exactly. All MEASURED.

| comparison | outcome |
|---|---|
| two 8-thread runs of `EXHALE_L15.x`, cap 200, traces | one md5; run logs differ in the wall seconds alone |
| a third 8-thread run | the same |
| `EXHALE_L15ctl.x` at 8 threads against them | the same |
| 1 thread against 8 threads, traces | one md5; the logs differ further only in the two thread-count lines and in the exit-time floating-point exception summary |
| 16 threads against 8 threads, traces | one md5 |
| `EXHALE_PTC_RAMP_GUARD=0`, cap 80, 8 threads, the two builds | every printed digit equal; wall seconds 308.37 and 307.52 |
| the checkpoint of outer iteration 21, 8 threads against 1 thread | byte-identical, md5 `3830aab11ec1018e740db6c5ab4fb11f` |
| the residual repeated inside one process at outer iteration 21 | bitwise equal, at 1 and at 8 threads, and equal across them |
| the action of the Jacobian repeated the same way | bitwise equal, at 1 and at 8 threads, and equal across them |

The solve itself, for the record: 105 outer iterations under a cap of 200,
ending `info = 2` at `||R||` 2.720E-08 with 38 non-monotone accepts and 152
Krylov products; pass 1 hands back the worst gated species row 2.84E-03 of
1.0E-05 at cell 301 (elemental transport O), mass 4.57E-10, momentum
3.48E-14, energy 2.72E-08. With the ramp guard off and the cap at 80 the same
fixture ends at `||R||` 2.664E-08 and mass 6.11E-10.

The trace of one run holds 105 iteration records, 105 merit records and 713
records of the linear cycle.

## 4. Limitations

- Finite repetition proves no absence. Three 8-thread runs, one 16-thread run
  and one 1-thread run of one fixture exercise one schedule distribution on
  one state sequence; a race not exercised by this grid, this composition and
  this load is not excluded by them. The machine carried other single-threaded
  work throughout, which makes the schedule less regular, not more, so the
  agreement is taken under load and not on an idle host.
- The comparison is of what the record covers. The unknowns, the composition,
  the residual, the model base point, both scalings, the active bound set, the
  pseudo-time, the right-hand side, the preconditioner input and output, every
  action of the operator and every Arnoldi column are in it. The Givens
  rotations, the least-squares coefficients and the line-search trials are not
  recorded separately; they are functions of what is, and the merit and the
  judged rows of every iteration are recorded, so a difference in them would
  have to cancel exactly to stay hidden.
- The checkpoint dump does not carry the module caches of the residual
  assembly that live outside `steady_newton` (`ionization_equilibrium`,
  `util_ion_eq`, `constrained_chemical_equilibrium`,
  `diffusive_photochemistry`), which have no accessor this module may read. A
  fresh process therefore cannot be put back into the checkpoint state, and
  the comparison used in its place is the repeated evaluation inside one
  process, which reads exactly those caches as the solve left them.
- The molecular configuration of step 5 was not run. The atomic cause is not
  understood in the sense the plan meant, because there is nothing to
  understand on this tree; whether the molecular carrier path behaves the same
  way is a separate measurement.
- Both binaries carry the uncommitted edits the other items of 2026-09-16 had
  in about twenty source files at the time of the build. That is why the
  control is a private build of the entry text of the one file this item
  edits, and not the tree's `EXHALE.x`, which is a build of an earlier state of
  those same files.

## 5. What follows

The recommendation of section 2 of the plan, to annotate the 8-thread
REPRODUCE.md records, is now a question of which statement to annotate with.
The measurement supports "the solve path of the element reload was bitwise
reproducible at 1, 8 and 16 threads on 2026-09-16", and it does not support a
general claim about every configuration; the `jfnk_thread_reproducibility`
suite the plan asks for would be this fixture's trace comparison, run at two
thread counts, with the trace kept in the suite directory. The fix candidates
of section 2 have nothing to act on: no located defect, so no correction and
no golden movement.

## 6. Step 5, the molecular configurations

Run after the sections above, on the path that carries the OpenMP regions the
atomic fixture never enters: the carrier relaxation and the H2 transport of
`diffusive_photochemistry.f90`, the constrained chemical equilibrium and
`System_HeH_mol.f90`. Private build `EXHALE_L15.x` md5
`21ba370bf1ae1dd0816c299fc2b7abc0`, rebuilt for this step because the tree's
other files had moved; every comparison below is inside that one binary. All
MEASURED.

**Identical, in both configurations.**

`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13`, the certified case,
reloaded from its own `output/*_IC.txt` with `Load IC? True` and `Restart
intent: stationary`, keeping the case's `Secondary_ionization: Immediate` and
`Molecular carrier transport: True`, `EXHALE_PTC_DTAU0=1.0
EXHALE_OUTER_PASSES=3`: two runs at 8 threads and one at 1 thread give ONE
trace md5, `d5768c1c3c906a65a60fe0e62043e125`, over 158 outer iterations, 158
merit records and 1098 records of the linear cycle. The run logs differ in the
three wall-second fields of the pass lines and, for the 1-thread run, in the
two lines that report the thread count, and in nothing else: the carrier
relaxation's ending, its kept displacement and its transport-step count are
equal line for line (pass 2 "the carriers stopped moving at the fixed wind;
displacement kept 1.17E-06 in 45 transport steps", pass 3 "7.01E-02 in 69
transport steps"), and so are the judged rows of all three passes (pass 3
`info = 0`, worst gated species row 5.61E-01 of 1.0E-05 at cell 223, carrier
balance H2, mass 5.34E-09, momentum 9.69E-13, energy 1.58E-08). The written
states agree in every data line (`Hydro_ioniz.txt`
`4b7ce7f1d9b86f34f3f6ea189943cac7`, `Ion_species.txt`
`0c29e3e53d11b4d4a5f7c8b6297a7de7`, all three runs). The 1-thread run took
2190, 2297 and 2870 s per pass against 302, 341 and 418 s at 8 threads.

`LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.083`, the STALLED
state, the trajectory the campaign's 8-thread hours were actually spent on,
same recipe with its own `input.inp` and `base.inp`: two runs at 8 threads
give ONE trace md5, `e465602dae6f7cfdae3dac4320891a55`, over 145 outer
iterations and 921 records of the linear cycle; the run logs differ in the
three wall-second fields alone, the carrier relaxation ending on the movement
bound with the same kept displacement and the same step count in both (pass 2
"1.67E-02 in 12 transport steps", pass 3 "1.65E-02 in 11 transport steps");
the written states agree in every data line (`Hydro_ioniz.txt`
`806f3e3a7fba26fb89bf6fe6c9173538`, `Ion_species.txt`
`6fa32cde8fb3c735edb72224c6716379`).

So the molecular path is reproducible in the same sense as the atomic one, on
a certified state and on a stalled one, and the limitation of section 4 stands
unchanged: this is five runs of two fixtures, not a proof that no schedule
anywhere exercises a shared write.

**The catalog binary itself**, `EXHALE.x` md5
`c2e9c9990b9f14f1be8cd77abca68945`, the one the 120 certified states of
`LHS1140b/models/` were produced with (it carries no trace key, so this is a
comparison of the run logs and the written states). MEASURED 2026-09-17, two
runs at 8 threads and one at 1 thread of each fixture. On the atomic element
reload (cap 200, `EXHALE_OUTER_PASSES=1 EXHALE_DIFF_OMEGA=0.5`) the two
8-thread logs differ in the wall seconds of the pass line alone and the
1-thread log further in the two thread-count lines and the exit-time
floating-point exception summary; every `(JFNK)` line of all three is equal,
and the written states agree in every data line (`Hydro_ioniz.txt`
`c2512194a8934a1d527903e46afea9c4`, `Ion_species.txt`
`9ce8f44c24e4b23da2eff923943a9602`, the same pair the private build wrote).
On `molecular_scalar_gj1132_kzz1e9/HeH2.13` (`EXHALE_PTC_DTAU0=1.0
EXHALE_OUTER_PASSES=3`) the three logs differ only in the thread-count lines
and the wall seconds, both passes reporting `info = 0` with worst gated
species rows 5.64E-07 at cell 244 and 1.68E-07 at cell 227 (carrier balance
H2), and the written states agree in every data line (`Hydro_ioniz.txt`
`59caf2c63097df1473b8bff304f20b5f`, `Ion_species.txt`
`dcb31893b8d72f3d8c08d19cde716bd5`; this pair differs from the private
build's because the day's other items moved the tree's sources between the two
binaries, which is a difference of physics text and not of threads).

**The item closes here**: no thread dependence was found on the catalog binary
or on the 2026-09-16 tree, atomic or molecular, at 1, 8 or 16 threads. The
2026-09-15 divergence is not reproducible and is attributed to nothing.
