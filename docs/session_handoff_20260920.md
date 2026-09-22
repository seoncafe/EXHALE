# Session handoff, 2026-09-20: the series on `docs/PLAN_20260919_rev1.md`

Written 2026-09-20 at the close of the series. The previous handoff is
`docs/session_handoff_20260919.md` (the D1 to D9 series of
`docs/PLAN_20260918_rev2.md`); its outline is the outline of this one. The plan
this series executed is `docs/PLAN_20260919_rev1.md`, the revision of
`docs/PLAN_20260919.md` after the review `docs/PLAN_20260919_review.md`; its
sections 7 and following carry an execution paragraph for every item, and the
update log records the series in `docs/Update_EXHALE_stage2.md` section 14.

**Every number in this document is READ from the item memo, the item report,
the plan or the file named beside it.** Where a worker labeled a number
MEASURED, it is quoted here as READ and the memo is named; the measurement is
that worker's. Four things were checked for this document and are said to be:
the md5 of the tree's `EXHALE.x`, whether `make -q` finds it current, which
source files differ between the catalog binary's source list and the tree, and
the movement table of the golden refresh, read from the strict check that
preceded it.

**The binary of record** is the tree's `EXHALE.x`, md5
`f89569405648f9ef8ff8dc78b1ab6810`, `make -q` clean (MEASURED here). The
catalog ran on the copy `LHS1140b/models/EXHALE_75d55d9d.x`, md5
`75d55d9d4fd0e748cd01d6e35713e34c`, which differs from it by one file
(section 2.3).

**GOLDEN REFRESH: one, on 2026-09-20**, on the binary of record, over eighteen
cases, after a strict check of all eighteen on it. The movement table is in the
update log's section 14 and is summarized in section 2.3 below.

---

## 0. What the series did, in one page

**Where it started.** The previous series left the boundary as a contact
(D2b), and with it three molecular states whose certificates the new boundary
refused. A review of the first plan text (`docs/PLAN_20260919_review.md`)
raised ten findings; the plan of record checked each against the source, and
against an execution where the finding was about behavior, and kept nine as
confirmed and one (the solved ghost electron count) with a narrower reading.
The order of work was changed so that measurement and record correctness came
before the boundary investigation, because the investigation's evidence is
read through them.

**What landed.** Every item of the plan but P2, which the user has not opened.

- **P0**: the catalog counts are now OF a published list, with the selection
  rule in words and as the function that is the rule.
- **P4a**: how a run ended is four facts kept apart, read from the route's own
  terminal event, with one rule for the runner, the publisher and the status
  tool.
- **P5a**: a certification report that is full stops the run instead of
  overwriting its last entry, and its capacity is the size of the inventory
  and not a literal. The capacity was exactly consumed in production.
- **P4b, P4b-2, P4c, P4d**: the seed of a solve is resolved once and travels
  to the manifest and to `REPRODUCE.md`; a flux-closure rung publishes at all,
  with its seed; the provenance of the 398 published generations is recovered
  into attachments beside them and nothing is rewritten; a record is bound to
  its log by bytes.
- **P5b**: one name selects the binary a test runs, a conflicting pair is
  refused, and every suite states and asserts the md5 of the file it runs.
- **R2 and P1 steps 1, 1b, 2 and 2b**: the refusal of a certified molecular
  state on reload was a seed dependence of the ghost composition solve, not a
  failure of the physical column. There is one root and the sweep did not
  reach it. The boundary now applies its own composition map again until the
  ghost is a fixed point of it, to an accuracy stated in the composition's own
  units and bracketed by what the base row can resolve and what the
  composition solve delivers.
- **P6b**: the time step is the face restriction of the cells that are
  evolved, on the spherical geometry the divergence uses, with the boundary
  Riemann problems kept.
- **P6c**: there is no 4.7 per cent inconsistency; the number compared a count
  prescribed at the base level with a count solved one cell below it, and it
  follows the grid.
- **P6d, P6e, P6f**: the charge-exchange turnover states its cell; the
  elemental base flux is a record of one evaluation, with both faces and both
  halves; the transit tool reads a solved state, with three defects of its own
  fixed.
- **The user's instruction of 2026-09-20**: a state that passes the
  certificate is certified whatever the run that produced it ended as. Two
  molecular cases are certified for the first time.
- **The catalog** was re-evaluated twice under the new boundary and the
  goldens were refreshed once.

**What is left.** P2, the four models whose hydrodynamic Newton makes no
progress, which is the one decision for the user (section 4), and the
technical items of section 4.2.

---

## 1. The items, one by one

### P0: the evidence scope

READ, `docs/lhs1140b_p0_p4a_20260919.md`. `status.py inventory [--write]`
writes `LHS1140b/models/CASE_INVENTORY.{md,json}`: the rule, the instant, the
binary with its source list, the criterion for certified, the totals by kind
and one row per state index. 137 state indexes, 83 certified at the time (84
ladder cases, 46 closure rungs, 4 `diffusion_check`, 2 archives, 1 study
directory). The three earlier counts differ by the tree and by how narrowly
the table was read, not by the rule.

### P4a: the termination classification

READ, the same memo. The route is read from the run's own `Restart intent:`
echo and only that route's terminal event is a verdict; an inner
`(JFNK|PTC) done info=` line is evidence. `state_missing`, `nonfinite_state`,
`stale_state`, `conflicting_evidence`, `marching_stop`, `no_certificate` and
(after the advisor's correction of the marching route's six endings)
`no_integration` are the classes. One rule for `run_case.sh`,
`publish_state.py` and `status.py`. `tests/termination_classification.sh` 14
of 14, RED 1 of 14. Read only over 398 records: 44 of the 306 bound records
would classify differently, all `no_verdict` becoming `no_integration`.

### P5a: the certification report's capacity

READ, `docs/lhs1140b_p5a_20260919.md`. `add_entry` returns the index it
created, every writer assigns through it, and an entry that does not fit stops
the run naming the equation, the capacity and the last entry. `cert_max_entries`
is the sum of the inventory and stays 32, which a transported-stage
configuration consumes exactly (`active equations 11 of 32 in the inventory`).
`certification` 123 of 123; `wasp_full` and `mol_carrier` byte-identical.

### P4b: the seed of a new solve

READ, `docs/lhs1140b_p4b_20260919.md`, contract `MODELS.md` section 9.12. The
runner resolves the seed once into `runs/<run id>/seed_identity.json`
(`exhale_seed_identity/1`) and hands that one file to the publisher and to
`write_reproduce.py`. Who chose it, what was done to it (`as_is`, `mapped`,
`converted`) and whether it is a generation of this case (then the parent) or
of another (then not) are all recorded. `tests/seed_provenance.sh` 33 of 33,
RED 5 of 33.

### P5b: the executable a test selects

READ, `docs/lhs1140b_p5b_20260919.md`. One policy in `src/tests/exhale_exe.sh`:
`EXHALE_EXE` selects, three names are aliases, a conflicting pair is refused,
every suite prints path, md5 and which name selected it and asserts the md5 of
the file that runs. New sentinel suite `src/tests/executable_identity`, 21
rows, none failing.

### R2 and P1 step 1, experiments 1 and 2

READ, `docs/lhs1140b_p1_step1_20260919.md`. The solved ghost electron count is
now `calc_ne`'s own count of the post-sweep composition (it moves by 1.022e-2
relative and nothing numerical moves). Three default-off diagnostics were
added (`EXHALE_GHOST_COMPOSITION_SEED`, `EXHALE_GHOST_SEED_WRITE`,
`EXHALE_GHOST_RECORD`). At one frozen column the ghost's seed decides the
certificate: 8.3353e-08 refused from the reservoir row against 2.9229e-09
certified from the state file's ghost rows, the two ghosts 4.7e-3 to 5.5e-2
apart in the trace ions and both accepted at 3e-26.

### P1 step 1, experiments 3 to 6

READ, `docs/lhs1140b_p1_step1b_20260919.md`. One root, which the sweep does
not reach: six decades of the inner stopping accuracy move the two states not
at all; the refused ghost is not a fixed point of the solve seeded with it;
and the seed-to-ghost map has a measured gain of 4.9e-3. The interior refresh
is not a term in the base row, the temperature closure removes most of the
difference, and the reload ladder settles with no drift and no cycle.

### P1 step 2

READ, `docs/lhs1140b_p1_step2_20260920.md`. `ioniz_eq` applies its own map
again at the ghost it returned, interior restored, until the ghost is a fixed
point; twenty applications stop the run by name. The seed is stated by the
boundary model and the model identity moves to v3, with a v2 state reported
and loaded rather than refused. The refused generation certifies.

### P1 step 2b

READ, `docs/lhs1140b_p1_step2b_20260920.md`. The contract's accuracy was in
the wrong units and stopped a run whose boundary was converged: the ghost
settles into an exact two-cycle whose moves are 6.5e-12 of the cell and
5.1e-3 of a stage's own value. The move is now absolute, in the composition's
own units, at 1e-11, with the count leg at 1e-6 against its own value.
`grid_and_gates` 305 of 305.

### P6d and P6e

READ, `docs/lhs1140b_p6d_p6e_20260919.md`. `cx_add_to_turnover` takes the cell
and goes through `cx_require_cell`; instrumented runs show the entry point is
not reached at all on five production cases, so the defect was a missing
contract. `element_base_flux` replaces a record that was identically zero at
face 0, and carries both faces and both halves with the evaluation that wrote
them; the advective half at the base is not steady on `lower_profile`,
changing sign between evaluations.

### P6f

READ, `docs/lhs1140b_p6f_20260919.md`. `EXHALE_TRANSIT_STATE=solution` is the
way a stationary run is read; three defects of the tool were fixed (the bare
`FileNotFoundError`, curves that named their inputs by path only, and O I
levels read from the `_adv` file whatever the selection). What the choice
costs on `wasp_full`: He I 10830 -0.64 per cent, H-alpha -1.28, H-beta -3.48,
Ly-alpha unchanged at its saturated ceiling.

### P6b and P6c

READ, `docs/lhs1140b_p6b_p6c_20260920.md`. The time step of an evolved cell is
`CFL V_j / max(A S)` over its two faces, minimized over cells 1 to N; the
entry text was too long by `dr/r` (1.96e-4) and let a ghost dictate the step
(89.7 per cent on a constructed case). The solved steady state does not move
(5.7e-12, same log10 Mdot 13.33). P6c: the four relations hold at their own
locations and the 4.7 per cent follows the grid (4.87e-2 against 5.50e-3 on
two grids).

### P4b-2, P4c and P4d

READ, `docs/lhs1140b_p4cd_20260920.md`, contracts `MODELS.md` sections 9.13
and 9.14. The closure runner publishes its rungs through the one publisher and
resolves each rung's seed first; `recover_provenance.py` attaches what can be
recovered beside each of 398 generations (176 established, 10 inferred, 212
not established) and rewrites nothing; records are bound to their logs by
bytes, and of the 311 that name a log on disk, 226 still carry the block and
85 do not.

### The certification rule of 2026-09-20

READ, the fourth pass of `docs/lhs1140b_catalog_refresh_20260919.md` and
`MODELS.md` section 9.15. A certificate is a property of the state, measured;
the records never confer it. `--evaluate` measures a named generation, or
`latest_certified`, or `latest_complete` where the index names none, and says
which; a certified evaluate child becomes `latest_certified`. Four bookkeeping
defects went with it (section 4 of the update log's subsection).
`tests/certified_evaluation.sh` 18 of 18, RED 3 of 18.

### The catalog and the du comparison

READ, the third and fourth passes of
`docs/lhs1140b_catalog_refresh_20260919.md`, and
`docs/du_stop_vs_stationary_20260920.md`. Section 3 below carries the catalog.
The du comparison answers what a state that meets the literature's flux-spread
criterion and fails this code's certification costs: at most 1.9 per cent in
any profile on the case measured, 1.3 per cent in the He I 10830 red-pair
equivalent width, and nothing at the two decimals of log10 Mdot.

---

## 2. The state of the tree

### 2.1 What this series changed, and by which item

The repository carries the uncommitted work of several series, so
`git status --short` does not by itself say what this series touched. What it
touched is read from the two binary manifests and from the memos.

**Production sources.** MEASURED here: of the 169 files of the catalog
binary's source list, nine differ between `BINARY_MANIFEST_2c3b0acc9983.txt`
(the binary this series started on) and `BINARY_MANIFEST_75d55d9d4fd0.txt`
(the binary the catalog ran on), and one more differs between the second and
the tree now.

| file | item |
|---|---|
| `src/modules/radiation/ionization_equilibrium.f90` | P1 step 2 and step 2b (the ghost fixed point), R2 (the electron count), the three diagnostics |
| `src/modules/states/base_boundary.f90` | P1 step 2 (the model identity v3) and step 2b, and the wording P6c corrected |
| `src/modules/files_IO/load_IC.f90` | P1 step 2: a v2 state is reported, not refused |
| `src/modules/time_step/eval_dt.f90` | P6b |
| `src/modules/radiation/charge_exchange.f90` | P6d |
| `src/modules/functions/binary_element_diffusion.f90` | P6e |
| `src/modules/time_step/certification.f90` | P5a |
| `src/modules/states/boundary_state_trace.f90` | named by the P1 step 1b and P0/P4a memos |
| `src/EXHALE_main.f90` | no memo of the series names it; what changed in it was not established |
| `src/modules/time_step/viscous_conduction.f90` | the P6c wording at `conduction_base_level_T`, corrected by the advisor after that item closed (comment only) |

**Catalog tools and contract**: `LHS1140b/models/run_case.sh`,
`publish_state.py`, `status.py`, `pick_seed.py`, `write_reproduce.py`,
`run_closure.sh`, the new `recover_provenance.py`, and `LHS1140b/MODELS.md`
sections 7, 8, 9.1, 9.2, 9.3, 9.5, 9.8 and the new 9.12 to 9.15 (P4a, P4b,
P4b-2, P4c, P4d and the certification rule).

**Post-processing**: `EXHALE_transit.py` and `exhale_transit_lib.py` (P6f).

**Tests**: new suites `src/tests/executable_identity` (P5b) and
`src/tests/transit_state` (P6f); new fixtures under `src/tests/boundary_state`
(`ghost_seed_and_ghost_record.sh`, `interior_composition_held.sh`,
`ghost_composition_fixed_point.sh`); new rows in `src/tests/certification`,
`src/tests/charge_exchange_rows`, `src/tests/element_operator` and
`src/tests/grid_and_gates`; new catalog fixtures
`LHS1140b/models/tests/{termination_classification,seed_provenance,closure_seed_provenance,provenance_recovery,certified_evaluation}.sh`.

**Documents**: `docs/Update_EXHALE_stage2.{md,tex,pdf}` section 14, this
handoff, the item memos of section 5, `docs/EXHALE_user_manual`,
`docs/EXHALE_physics_and_algorithms`, `docs/EXHALE_BC_and_IC` and the
`docs/physics_overview/*.tex` sections behind them,
`docs/lhs1140b_catalog_state_20260922.tex` with its three figures.

### 2.2 New keys and their defaults

**No input key of `input.inp` was added by any item of this series** (READ,
the memos and the header of `BINARY_MANIFEST_75d55d9d4fd0.txt`). Every
diagnostic this series added is default off: unset, no file is opened, no
branch is entered and no number moves.

| key | default | item |
|---|---|---|
| `EXHALE_GHOST_COMPOSITION_SEED` | unset: the sweep is entered at the rows the caller installed. Set, it names the source (`the_state_file_ghost_rows`, `the_previous_sweep_ghost`, or a path) | P1 step 1 |
| `EXHALE_GHOST_SEED_WRITE` | unset; set, writes the installed ghost rows in the format the seed key reads | P1 step 1 |
| `EXHALE_GHOST_RECORD` | unset; set, writes one record per evaluation to `output/ghost_record.txt` | P1 step 1 |
| `EXHALE_INTERIOR_COMPOSITION_HELD` | unset (off): the declared diagnostic policy that holds the interior's eliminated species | P1 step 1b |
| `EXHALE_EXE` | unset: `$ROOT/EXHALE.x`. `EXHALE_RESID_EXE`, `EXHALE_SPECIES_EXE`, `EXHALE_COUPLED_EXE` are aliases; names resolving to different files are refused | P5b |
| `EXHALE_EXE_IDENTITY_ONLY` | unset; `1` states the identity and runs nothing | P5b |
| `EXHALE_BIN` | unset: `status.py` takes the newest installed catalog copy; set, it names the binary of record | the catalog passes |
| `RUN_CLOSURE_PUBLISH_ONLY` | unset; `1` defines the closure publication and runs nothing | P4b-2 |
| `pick_seed.py --seed-provenance` | off: every existing caller sees the same lines | P4c |
| `run_case.sh --generation <id>` | absent: `latest_certified` where the index names one, else `latest_complete`, with the run saying which | the certification rule |
| `ghost_composition_fixed_point_move`, `ghost_count_fixed_point_move`, `ghost_composition_fixed_point_passes` | 1e-11 (the composition's own units), 1e-6 (against its own value), 20 (constants) | P1 step 2 and step 2b |

What changed in the contract without a key: the boundary model identity is
`..._ghost_fixed_point_seed_reservoir_row_v3` and travels in every state file;
a v2 state is loaded with a note instead of being refused; the transit tool's
`solution` selection is how a stationary run is read, and no `_adv` file is
manufactured for one; a certification report that fills stops the run.

### 2.3 The binaries and the goldens

**The binary of record** is the tree's `EXHALE.x`, md5
`f89569405648f9ef8ff8dc78b1ab6810`, built 2026-09-20 11:32, `make -q` clean
(both MEASURED here). It carries one file beyond the catalog copy
`LHS1140b/models/EXHALE_75d55d9d.x` (md5
`75d55d9d4fd0e748cd01d6e35713e34c`, built 2026-09-20 06:12, manifest
`BINARY_MANIFEST_75d55d9d4fd0.txt`): the comment correction in
`viscous_conduction.f90` (MEASURED here: the md5 of every other file of that
manifest is the file's md5 in the tree now). The catalog's 83 re-evaluations,
its five re-solves and the nine molecular measurements all ran on the catalog
copy; the regression matrix and the golden refresh ran on the tree's binary,
which `run_check.sh` copies into each case directory (MEASURED: the copy in
`wasp_full_newton/` and in `mol_carrier/` carries the tree's md5).

**The goldens were refreshed once, on 2026-09-20 13:53**, on the binary of
record, over eighteen cases, and the previous set was replaced in place; the
newest archived set in `backup/regression/` is `golden_20260917/`. The
movement against the 2026-09-19 set, read from the strict check that ran
before the refresh, is tabulated in `docs/Update_EXHALE_stage2.md` section 14.
In summary: `iontrans_metals`, the certified stationary case, is identical;
`wasp_full_newton`, the converged marching case, moves by 1.4e-11 and 1.1e-10;
the two `du`-stop cases move by 1.9e-1 in the velocity of the base rows near
`r` = 1.075, where it alternates in sign cell to cell and the absolute
movement is 18 to 26 cm/s, and by 6.6e-4 to 9.0e-4 elsewhere; every other case
is a fixed-step snapshot of a run that is not stationary and moves with its
trajectory, by 9.0e-5 to 1.9e0. The shorter time step of P6b moves every
trajectory, and the molecular cases carry the boundary composition contract as
well.

**The strict re-check after the refresh** finished with `REGRESSION PASS (all
cases byte-identical)` over all eighteen cases (MEASURED,
`REGRESSION_REL_TOL=0`, `OMP_NUM_THREADS=1`), in both state files and in the
two `_adv` products where the case writes them. `parse_golden/` was not
touched by this series.

---

## 3. The catalog

READ, `LHS1140b/models/CASE_INVENTORY.md` (taken 2026-09-20 14:42:33 KST), the
third and fourth passes of `docs/lhs1140b_catalog_refresh_20260919.md`, and
`LHS1140b/MODELS.md` sections 7, 8 and 9.

**137 state indexes, 85 certified**, by kind 84 ladder cases (76 certified),
46 closure rungs (9), 4 `diffusion_check` (0), 2 archives (0), 1 study
directory (0). The criterion is unchanged: the index names a
`latest_certified` generation whose entry carries no `stale_under`. A case
with no `latest_certified` is one no measurement has been recorded for under
the binary of record, which says nothing about whether its state solves the
equations.

**The 83 atomic cases that carried a certificate all certify again** under the
new boundary, with log10 Mdot unchanged in all 83 and the largest He 10830
equivalent-width change 2.6e-14 per cent A. The boundary composition contract
applies to none of them, its condition being a molecular network with helium.

**The molecular states now carry measurements where they carried none.** Nine
states hold a molecular index: seven case directories and the two `pre_D9s3`
archival copies, whose states are by md5 the generations of two of those
cases. Measured through the runner on the binary of record, three certify and
six are refused, and every one of the nine reaches the ghost composition fixed
point in two or three applications at 6e-15 to 5e-12 against 1e-11, so the
contract is not what separates them.

**The two cases that certify**, the first certified molecular cases of the
catalog:

| case | generation | log10 Mdot | He 10830 EW [%A] |
|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | `g0004_20260919T103629Z_9b8394a3` (and the `pre_D9s3` state `g0002_...552403ab`) | 7.94 | 2.4061 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | the `pre_D9s3` state `g0002_...19d7ee4d` | 7.90 | 1.5603 |

Their new `latest_certified` are the evaluate children
`g0008_20260920T052941Z_c3c466a4` and `g0008_20260920T053924Z_43df1674`.
`HeH2.13`'s own newer `latest_complete`, the 2026-09-20 solve
`g0004_...7033d8a4`, is NOT certified.

**The five that do not, and by which row:**

| case | the rows that refuse it |
|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` at `g0004_20260920T025822Z_7033d8a4` | mass 3.868E-05 above 3.0E-12 at cell 285; momentum 3.848E-04 above 1.0E-08 at cell 500; energy 7.977E-01 above 1.0E-06 at cell 2; carrier balance H2 7.747E-02 above 1.0E-05 at cell 224 (a wind cell); elemental transport He/H 1.731E-04 above 1.0E-05 at cell 217 |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | mass 1.041E+00 above 3.4E-10 at cell 2; momentum 1.450E-03 above 1.0E-08 at cell 1; energy 1.021E+00 above 1.0E-06 at cell 1; carrier balance H2 1.000E+00 above 1.0E-05 at cell 432; elemental transport He/H 1.494E-04 above 1.0E-05 at cell 217 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | carrier balance H2 8.150E-04 above 1.0E-05 at cell 499; elemental transport He/H 2.503E-04 above 1.0E-05 at cell 280 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | energy 4.634E-01 above 1.0E-06 at cell 257; carrier balance H2 8.922E-03 above 1.0E-05 at cell 261 |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | carrier balance H2 7.371E-02 above 1.0E-05 at cell 500 |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | energy 1.000E+00 above 1.0E-06 at cell 246; carrier balance H2 2.363E-02 above 1.0E-05 at cell 306 |

The table has six rows because `kzz1e9/HeH2.13` appears twice: it is certified
at the state named above and refused at its own newer `latest_complete`. The
five cases that carry no certificate at all are the other five rows. A carrier
balance H2 in the wind refuses all of them, an interior energy row three, and
in the photochemical HeH9 the base rows at order one.

**Nothing was deleted.** Every case-level file an evaluation replaced is in
`runs/<run>/superseded_case_products/`, the two `pre_D9s3` copies were
measured without being run or repaired, and `molecular_photochem_.../HeH2.09`
and `molecular_scalar_..._wellmixed/HeH2.13` publish no index and hold no
state.

---

## 4. Open items and the one user decision

### 4.1 What the user decides: P2

**P2, the hydrodynamic Newton makes no progress on four models.** The plan
removed `.L22/i3_alt` from the first text's list of five, because it was
published as a generation and never re-solved. What was measured on the four,
READ from the plan's P2 and from the third pass of the catalog memo:

| model | what was measured |
|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | mass row frozen at 1.52e-1 from pass 12, `info = 2` at every pass, stopped by the 3 h ceiling in pass 33. At the class budget of the third pass: 5 outer passes, `info = 2` at every one, mass 2.46e-1, momentum 3.29e-7, energy 1.22e-1, unmoved from pass 4 to 5, gated species row 4.38e-2 at cell 331 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | energy row 0.998 at cell 1 from pass 1. At the class budget: 4 passes, `info = 2` at every one, energy row 9.98e-1 at the base, mass 7.33e-8, momentum 8.25e-14, gated species row falling 1.3e-2 to 9.77e-4 at cell 354 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | mass row 0.578 at cell 1, stopped by the stall rule twice. At the class budget: 4 passes, mass 1.16 and energy 1.10 unmoved from pass 3 to 4, gated species row 3.29e-2 at cell 387 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | 40 passes then the 6 h ceiling; at the third pass the gated carrier balance H2 bands between 5e-2 and 2e-1 from pass 20 to the ceiling with no trend and the hydrodynamic rows stay at 1e-2 to 1 |

Six times the wall clock and 2.25 times the passes left the stuck rows where
they were, so budget is not the answer. The three x0.01 cases are open because
their old certificates were earned on the pre-D2b blended base face, and under
the accepted contact law their base face flux is 38 to 178 times the window
mean: the refusal is physical and a solve is what is missing. If the user
opens the item, its design is the five steps of the plan's P2, ending at the
measurement and the report, with no solver control change and no reopening of
the held movement bound. The two `coupled_block_jacobian` failures and the
`residual_determinism` closure spread below belong inside it.

`molecular_scalar_gj1132_kzz1e9/HeH9.7` is not in this class and is worth
saying so: it converges and does not arrive, its gated carrier row falling
monotonically to 5.13e-4 against 1e-5 when the 40 passes are spent, a factor
1.2 from where the previous pass left it.

### 4.2 Technical items noticed and not fixed

Each READ from the memo named.

- `coupled_block_jacobian` fails 2 of 7 rows on the control build as well: the
  carrier action at 8.954e-1 against 1e-2, whose central difference has a
  largest response in support of 4.79e-12, and the species-row scale at the
  last digits of 3 of 13 (R5; they need a noise-floor ladder and an arithmetic
  contract, not a bitwise demand).
- `residual_determinism`'s
  `closure_spread_within_the_row_tolerance_atomic_elem_newton` measures
  1.715E+04 against a bound of 1, identically before and after P5b (P5b).
- `coupled_source_step` fails 2 of its 44 rows on the entry text as well
  (`mean_passes_never_exceeds_the_worst` measured `absent`,
  `extrapolated_state_is_admissible` `0.000e+00_over_0`), which reads as a
  suite that no longer reaches the events it measures (P1 step 1b, step 2).
- The three call sites of `cx_add_to_turnover` still state no cell; the one in
  `constrained_chemical_equilibrium.f90` sits on a path none of the five
  instrumented runs entered (P6d).
- `element_base_flux` has no reader, and adding the base face to
  `element_flux_profile.txt` would move a reported product of `lower_profile`
  (P6e).
- The ghost's `p` and `rho` are built with the reservoir's particles per unit
  mass while the caloric map uses the composition's own count, 2.7e-08 and
  1.3e-07 apart on the two states measured (P6c).
- The base continuity row still rests on a near cancellation, and the atomic
  ghost is outside the composition contract (P1 step 2, step 2b).
- The contract's accuracy 1e-11 sits a factor 1.5 above the lowest plateau
  measured, so a configuration whose composition solve plateaus higher stops
  the run by name; the answer is to measure that plateau (P1 step 2b).
- `src/tests/grid_and_gates/mass_flux_spread_functional.py` hardcodes its work
  directory to the tree's `build/tests/` and ignores `EXHALE_TEST_OUT`, which
  every other target of that suite honors (P1 step 2b).
- The summary line of the hydrodynamic mass row prints `tol=` of the verdict
  cell beside `max=` and `cell=` of the cell carrying the maximum, which reads
  as a contradiction until the verdict line below it is read (P1 step 1).
- A metals-off run writes its metal columns at a floor of order 1e-301 that
  does not survive a load and a write exactly, so an all-column comparison of
  two state files must floor them (P1 step 1b).
- `status.py`'s wall-clock and outer-pass cells are still unlabeled
  inferences (P4c); the two `pre_D9s3` copies still refuse
  `publish_state.py verify`, by the standing decision that an archival copy is
  not repaired (the certification rule report).
- `molecular_scalar_gj1132_kzz1e9/HeH9.7` `g0004` carries both a
  `stale_under` from the previous binary's refusal and a `certified_under`
  from this binary's measurement; both are kept (the certification rule
  report).
- The ten provenance facts inferred from a `REPRODUCE.md` window stay
  inferred and no reader uses them (P4c).

---

## 5. Paths

**The plan and its review**: `docs/PLAN_20260919_rev1.md` (the plan of record,
with the execution sections), `docs/PLAN_20260919.md` (the superseded first
text), `docs/PLAN_20260919_review.md`.

**Memos of this series** (all under `docs/`):
`lhs1140b_p0_p4a_20260919.md`, `lhs1140b_p5a_20260919.md`,
`lhs1140b_p4b_20260919.md`, `lhs1140b_p5b_20260919.md`,
`lhs1140b_p1_step1_20260919.md`, `lhs1140b_p1_step1b_20260919.md`,
`lhs1140b_p1_step2_20260920.md`, `lhs1140b_p1_step2b_20260920.md`,
`lhs1140b_p6d_p6e_20260919.md`, `lhs1140b_p6f_20260919.md`,
`lhs1140b_p6b_p6c_20260920.md`, `lhs1140b_p4cd_20260920.md`,
`lhs1140b_stationary_D1b_20260919.md` (D1b, executed at the opening of this
series), `lhs1140b_catalog_refresh_20260919.md` (third and fourth passes) and
`du_stop_vs_stationary_20260920.md`.

**The catalog summary document**: `docs/lhs1140b_catalog_state_20260922.tex`
and its `.pdf` (7 pages), with `docs/figures/du_stop_vs_stationary_profiles.pdf`,
`docs/figures/du_stop_vs_stationary_differences.pdf` and
`docs/figures/du_stop_vs_stationary_he10830.pdf`.

**Item reports** (the session scratch directory
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/`):
`P0P4a_report.md`, `P5a_report.md`, `P4b_report.md`, `P5b_report.md`,
`P1step1_report.md`, `P1step1b_report.md`, `P1step2_report.md`,
`P1step2b_report.md`, `P6de_report.md`, `P6f_report.md`, `P6bc_report.md`,
`P4cd_report.md`, `certify_solution_report.md`, `catalog_contract_report.md`,
`du_vs_stationary_report.md`, with each item's run directories and private
builds under `<item>/`, and the golden-refresh logs `gr1.log` (the strict
check before the refresh, which is the movement table), `gr2.log` (the
refresh) and `gr3.log` (the strict re-check).

**The update log**: `docs/Update_EXHALE_stage2.md` section 14 and its typeset
twin `docs/Update_EXHALE_stage2.pdf` (290 pages).

**The catalog contract**: `LHS1140b/MODELS.md` sections 9.12 to 9.15;
`LHS1140b/models/{run_case.sh,run_closure.sh,publish_state.py,status.py,pick_seed.py,write_reproduce.py,recover_provenance.py}`;
the fixtures under `LHS1140b/models/tests/`;
`LHS1140b/models/CASE_INVENTORY.{md,json}`.

**The test suites of this series**: `src/tests/executable_identity/`,
`src/tests/transit_state/`, the three new fixtures of
`src/tests/boundary_state/`, and the new rows of `src/tests/certification/`,
`src/tests/charge_exchange_rows/`, `src/tests/element_operator/` and
`src/tests/grid_and_gates/`.
