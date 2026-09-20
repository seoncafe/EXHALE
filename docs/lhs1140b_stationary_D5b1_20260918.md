# D5b-1: the evaluation routes derive the boundary after the composition sweep

Item D5b-1 of `docs/PLAN_20260918_rev2.md` (section D5, "the item"), the
bounded first increment of D5b: the largest of the causes D5a isolated, and
the only one that needs no new interface. Opened by
`docs/lhs1140b_stationary_D5a_20260918.md` section 9 item 1 and its table 5.
2026-09-18.

Every number below is labelled MEASURED (this item ran it) or READ (from a
file, a log or a source line).

---

## 1. Verdict

**The evaluate route, the steady trial of the update map and the
`EXHALE_RESIDUAL` diagnostic assembled their residual on a boundary built
from the composition the state was LOADED with, while the sources, the
pressure map and the sound speeds of the same residual were those of the
composition one sweep later.** Under the caloric equation of state those are
different gases, and Rec_BC puts the cached `base_face_W` straight into the
left slot of the base face, so the assembled row is the row of neither. The
three routes now derive the boundary once more, from the composition the
residual is assembled with, which is the call the JFNK residual already makes
for the stated reason (READ, `steady_newton.f90`, the second `Apply_BC` of
`eval_residual`) and the order the marching loop has always had.

MEASURED, the cell-1 continuity row of one no-step evaluation, in floors of
that row's own rounding floor (gate ten):

| state | entry text | this item | D5a's reading |
|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | 15.98 | **9.06** | 15.98 -> 9.06 |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | 31.24 | **3.09** | 31.24 -> 3.09 |
| `.L22/i3_alt` | 19.04 | **0.08** | 20.93 -> 0.08 |
| `.L26/fid_resolve` (atomic) | 1.52 | **1.52** | no movement |

All three molecular states pass from NOT CERTIFIED, on one refusing entry
each and that entry the hydrodynamic mass row, to CERTIFIED. The atomic state
does not move at all, and the reason is measurable: the boundary of an atomic
mixture carries no composition, so the two orders give the same boundary
bitwise (section 4).

**The boundary the residual reads is now the boundary of the state installed
at that instant, exactly.** MEASURED, the largest relative difference between
the cached base face state, ghost cell averages and lower face state and the
boundary `base_boundary_states` gives for the installed state, at
`assemble_residual`: `0.0000000000000000E+000` on all four states, against
4.3504e-11, 5.9572e-12 and 1.4504e-11 on the three molecular states in the
entry text's order (and 0 on the atomic one). The 4.3504e-11 is D5a's
table-1 entry "C to D, `base_face_W` v 4.350e-11" met again.

**Nothing moves on a route that already rebuilt the boundary.** MEASURED,
byte identity of the data lines of `Hydro_ioniz.txt` and `Ion_species.txt`
against a build of the entry text: `carrier_model_a_newton`
(`Restart intent: stationary`, three outer passes, JFNK cap five) identical
including every header line; `wasp_full` at 400 marching steps identical in
every data line, the only difference in any file being the run timestamp of
the provenance header.

**What this item does NOT do.** It does not change the boundary rule, the
ghost composition prescription, the reservoir counts or the loader. The ghost
composition is still READ from the restart pair and still closed by
alternation across sweeps, the prescribed reservoir and the solved ghost
counts still share one name, and the cached face state still carries no
validity. Those are D5b-2, and the second and third of them are the remaining
9.06 floors of the molecular reference state (READ, D5a section 8).

---

## 2. What changed

| file | change |
|---|---|
| `src/EXHALE_main.f90` | `Apply_BC(u)` after the composition sweep and the pressure refresh and before `assemble_residual`, at three sites: the evaluate route of `stationary_state_of_the_loaded_restart`, the steady residual of `update_map_begin_step`, and the `EXHALE_RESIDUAL` diagnostic block. The invariant is stated once, at the first of them; the other two name it. Also `base_mass_rows_of_this_evaluation`, a default-off writer of the base continuity rows, and one default-off re-assembly after the face-flux report. |
| `src/modules/states/boundary_state_trace.f90` (new) | the measurement, all of it behind `EXHALE_BOUNDARY_TRACE`: the gap between the cached boundary and the boundary of the installed state, and the signed base mass rows with their scales and floors. `EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition` leaves the boundary of the entry composition standing, so that one binary measures both orders. |
| `Makefile` | one `SRC` line. |
| `src/tests/boundary_state/` (new) | `run.sh` and `boundary_of_the_installed_composition.sh`. |

The gap measure is the largest relative difference over `base_face_W`,
`base_ghost_W(:,1-Ng:0)` and `base_face_lower_W` between the values cached in
`BC_Apply` and the values `base_boundary_states` returns for the state
installed at the moment of the call, derived through `U_to_W_interior`
exactly as `Apply_BC` derives them. The branch diagnostics
`base_boundary_states` keeps of its last evaluation are saved and put back,
so an armed measurement does not enter the counts a run reports.

---

## 3. Binaries, snapshot and provenance

The working tree is under edit by other items, so the measurement works on a
frozen copy of `Makefile` and `src/`, and the tree's own `EXHALE.x` and
`build/` were never touched.

| | value |
|---|---|
| snapshot taken | 2026-09-18, from the working tree at git `3c73905` plus its uncommitted files |
| entry-text manifest (md5 of the md5 list of every `.f90`/`.inc`) | `f0aeaf8f3bc0ea147f2554a0fae7d5ca` |
| measured manifest, the same text plus this item's change | `1d1195c8ca5ad96ffde34b752282768c` |
| control build (entry text) | `0d4ff2af447ac931fe630eb88809fd7e` |
| measured build | `4e02c6807587f0d88cdacb944a65ed25` |
| compiler | conda-forge gfortran 16.2.0, `-O3 -fopenmp`, OpenBLAS of the same prefix |
| every run | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, `Restart intent: stationary evaluate` unless stated |

**The diagnostics change nothing when they are off**, and nothing when they
are on. MEASURED on all four states: the measured build with every key unset
and the same build with `EXHALE_BOUNDARY_TRACE=1` produce
`Hydro_ioniz.txt`, `Ion_species.txt`, `Hydro_ioniz_adv.txt`,
`Ion_species_adv.txt`, `Cooling_breakdown.txt` and `Heating_breakdown.txt`
with identical data lines, the only difference in any file being the run
timestamp of the provenance header.

**The experiment key reproduces the entry text exactly.** MEASURED: the
measured build with `EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition`
gives 15.9842, 31.2415, 19.0439 and 1.5202 floors on the four states, which
are the control build's readings to every printed digit, and every product's
data lines are identical between the two.

**One reading of D5a is not reproduced, and the reason is the tree, not
this item.** `.L22/i3_alt` reads 1.4863E-08 (19.04 floors) on the entry text
of this tree where D5a's snapshot binary read 1.6332E-08 (20.93); 1.4863E-08
is what L37 recorded for that state's evaluate route (READ, D5a section 0),
so the tree has moved back to L37's value since D5a's snapshot was taken at
08:33 on the same day. The corrected reading, 0.08 floors, is D5a's.

---

## 4. The four states, before and after

MEASURED with `EXHALE_MASS_FLOOR_SCAN=1`. The row measure and the floor are
that block's own cell-1 columns, the scaled continuity row and its rounding
floor.

| state | row measure, entry text | row measure, this item | floor | floors before | floors after |
|---|---|---|---|---|---|
| `HeH2.13` molecular | 1.24750E-08 | 7.07240E-09 | 7.80460E-10 | 15.9842 | 9.0618 |
| `HeH9.7` molecular | 2.33980E-08 | 2.31320E-09 | 7.48940E-10 | 31.2415 | 3.0886 |
| `.L22/i3_alt` molecular | 1.48630E-08 | 6.37280E-11 | 7.80460E-10 | 19.0439 | 0.0817 |
| `.L26/fid_resolve` atomic | 1.20160E-09 | 1.20160E-09 | 7.90440E-10 | 1.5202 | 1.5202 |

The certification verdicts of the same runs, READ from the run logs:

| state | entry text | this item |
|---|---|---|
| `HeH2.13` molecular | NOT CERTIFIED, 1 refusing entry: hydrodynamic mass row, 1.248E-08 above 7.8E-09 at cell 1; exit 2 | CERTIFIED in the wind, r >= 1.20; mass row within; exit 0 |
| `HeH9.7` molecular | NOT CERTIFIED, 1 refusing entry: hydrodynamic mass row, 2.340E-08 above 7.5E-09 at cell 1; exit 2 | CERTIFIED in the wind; exit 0 |
| `.L22/i3_alt` molecular | NOT CERTIFIED, 1 refusing entry: hydrodynamic mass row, 1.486E-08 above 7.8E-09 at cell 1; exit 2 | CERTIFIED in the wind; exit 0 |
| `.L26/fid_resolve` atomic | CERTIFIED in the wind; exit 0 | unchanged, every printed measure identical |

The three molecular states' files all claimed `certified=T
cert_reason=certified_in_wind`; on the entry text each printed "original
claim: NOT REPRODUCED", and each now reproduces it. This is the certification
reading of a molecular state on the evaluate route moving, which is the
correction; no tolerance was touched.

**Why the atomic state does not move.** MEASURED, the gap between the cached
boundary and the boundary of the installed state at `assemble_residual`, in
the entry text's order: 4.3504466172418882E-011 (`HeH2.13`),
5.9571936364851129E-012 (`HeH9.7`), 1.4503878633850953E-011 (`i3_alt`) and
exactly zero on the atomic state. The ghost continuation of an atomic mixture
reads `adiabatic_index_from_state`, which returns `gamma_ad` with no
composition in it (READ, D5a section 4), so one sweep cannot move that
boundary and there is nothing for a second derivation to correct.

---

## 5. The face cache, recorded

The report that follows the residual (`stationary_face_mass_flux`) installs a
boundary of its own on a copy of the conserved array, and the cached face
state carries no validity; making it carry one is D5b-2. The signature D5a
item 5 measured is therefore recorded here and not gated: the continuity row
of cell 1 re-assembled on the SAME conserved array after the report, against
the row the residual gave, in floors.

| state | entry-text order, at the residual / after the report | this item, at the residual / after the report |
|---|---|---|
| `HeH2.13` molecular | 15.9846 / 12.7990 (**-3.19 floors**, base face mass flux 2.4862E-09 relative) | 9.0619 / 9.0619 (**0**, flux 0) |
| `HeH9.7` molecular | 31.2421 / 33.0239 (**+1.78**, 1.3345E-09) | 3.0886 / 3.0886 (**0**, 0) |
| `.L22/i3_alt` molecular | 19.0434 / 18.0675 (**-0.98**, 7.6168E-10) | 0.0817 / 0.0817 (**0**, 0) |
| `.L26/fid_resolve` atomic | 1.5202 / 1.5202 (**0**, 0) | 1.5202 / 1.5202 (**0**, 0) |
| `carrier_model_a_newton` IC, a relaxation snapshot | 5.94127E+11 / 5.94127E+11 (**-1.455E+05 floors**, 2.2997E-07) | 5.94130E+11 / 5.94130E+11 (**-0.82**, 1.2963E-12) |

MEASURED. The entry-text column reproduces D5a's -3.18 and +1.78 floors on
the two states they share; `.L22/i3_alt` differs for the tree reason of
section 3.

**On the four converged states the signature vanishes, and the reason is not
a validity.** With the boundary derived from the installed composition, the
boundary the report rebuilds on its copy is the same boundary, because it
refreshes the composition-derived state from the same `f_sp` and the same
conserved array before its own `Apply_BC` (READ, `stationary_operator.f90`
lines 158 to 164). The two caches then coincide by construction, not by a
check. On a state far from stationary the coincidence fails: the
`carrier_model_a_newton` snapshot still moves by 0.82 floors after the
report, because the composition-derived state installed at the end of the
route is the post-process's and not the work state's. The cache still needs
its validity.

---

## 6. Byte identity on the routes that already rebuilt the boundary

MEASURED, control build against measured build, single-threaded, on scratch
copies; no regression directory was written to.

| case | how it was run | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|
| `backup/regression/carrier_model_a_newton` | `Restart intent: stationary`, `EXHALE_OUTER_PASSES=3 EXHALE_JFNK_MAXIT=5` | identical, every line including the headers | identical, every line |
| `backup/regression/wasp_full` | marching, `EXHALE_MAXSTEPS=400` | every data line identical; two header lines differ, both the provenance run timestamp | identical, every line |

`carrier_model_a_newton` enters the evaluate route's residual before its
solve, so the new derivation runs there too; it changes the ghosts of a state
whose next reader (the JFNK residual) derives the boundary itself before
using them, and the solve is bit for bit what it was. MEASURED, its `run.log`
differs in sixteen lines: the pre-solve residual norms of the state, which is
the reading this item corrects (mass 9.390325E-01 -> 9.390271E-01, 5.7e-6
relative; momentum 4.964298E-01 -> 4.964295E-01), the last digits of the
pre-solve face mass flux report (base face 9.425012038911570E+00 ->
9.425012038899352E+00, 1.3e-12 relative), and the wall-clock time of each of
the three outer passes. Every solver line, every `info`, every gated row and
every written column is identical.

---

## 7. The suite

`src/tests/boundary_state/boundary_of_the_installed_composition.sh`, on a
copy of `backup/regression/carrier_model_a_newton/IC/` evaluated with
`Restart intent: stationary evaluate`. A molecular base is required: on an
atomic mixture the measure is zero in both orders and the test would assert
nothing.

| row | measured | reference |
|---|---|---|
| `face_state_gap_at_assemble_residual` | zero (`0.0000000000000000E+000`) | zero |
| `face_state_gap_of_the_entry_composition_boundary` | nonzero (`4.2856189743416476E-007`) | nonzero |

The second row exists so that a row which passes by measuring nothing is
caught. RED, MEASURED: with
`EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition` set for the
whole suite, the first row reads `4.2856189743416476E-007` against a
reference of zero and fails; on the control build of the entry text the trace
does not exist at all and the suite fails on its first row. GREEN, MEASURED:
both rows pass on the delivered build.

The suite also prints, ungated, the cell-1 continuity row at the residual and
after the report, in floors. The fixture is a relaxation snapshot and its row
stands at 5.9E+11 floors in both orders; nothing about its stationarity is
asserted, and the measure under test is a property of the evaluation.

Suites run whole, MEASURED with `EXHALE_OBJDIR`/`EXHALE_EXE` pointing at the
private build: `boundary_state` (2 assertions, 0 failing), `steady_species_rows`
(201 assertions, 0 failing, with `EXHALE_SPECIES_EXE` set so the atomic and
carrier reload rows run), `certification` (84 assertions, 0 failing) and
`grid_and_gates` (260 assertions, 0 failing, every test program passed).
Those are the suites whose drivers link the production objects and therefore
the new module; no suite driver links `EXHALE_main.o`.

---

## 8. Noticed outside the scope, reported and not changed

- **Two further sites have the same order, and both are a solve's own
  certification.** READ, `EXHALE_main.f90`: inside
  `steady_wind_with_element_diffusion` the joint test's residual is assembled
  after that pass's `ioniz_eq` sweep with no `Apply_BC` between, and the
  final certification of a stationary restart is assembled after that routine
  returns. Both were left alone: this item's scope is the three evaluation
  routes, and a change at either of those two moves the certification of
  every solved state, which is a measurement of its own. The
  `carrier_model_a_newton` identity of section 6 is evidence that the
  boundary a solve leaves standing and the boundary of its returned
  composition are close on that case, not that they are the same.
  `report_marching_stop` is NOT such a site: the marching loop ends its step
  with `Apply_BC` after the composition update, so the state it reports on
  already carries the boundary of its own composition.
- **`stationary_face_mass_flux` leaves the composition-derived state of its
  own refresh installed** (READ, `stationary_operator.f90` lines 158 to 164),
  which is what makes the coincidence of section 5 hold on converged states
  and fail on a snapshot. It is D5b-2's to make explicit.
- The three items D5a reported outside its scope stand unchanged: the comment
  at `ionization_equilibrium.f90` lines 3176 to 3190 defends a base pressure
  boundary condition `base_boundary` replaced; `write_output.f90` writes the
  state files with list-directed output and no declared round-trip format;
  and `.L26/fid_resolve/output` and
  `atomic_scalar_gj1132_kzz1e9/HeH2.13/output` are the same state pair.

---

## 9. How to repeat any of it

- one reading: copy `input.inp`, `base.inp` and the `_IC` pair into a scratch
  directory, set `Restart intent: stationary evaluate`, run with
  `EXHALE_MASS_FLOOR_SCAN=1`; the cell-1 line of the `[mass floor scan]`
  block carries `|R_1|/s_1` and the floor.
- the boundary measure and the base mass rows: add `EXHALE_BOUNDARY_TRACE=1`;
  the blocks are written to `output/boundary_trace.txt` and the gap is
  printed to the log.
- the entry text's order from the same binary: add
  `EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition`.
- the suite: `EXHALE_EXE=<binary> EXHALE_TEST_OUT=<scratch>
  src/tests/boundary_state/run.sh`.
