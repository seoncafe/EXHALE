# D4b: the returned-state contract, and what a refused trial puts back

Item D4b of `docs/PLAN_20260918_rev2.md`, taken on the measurement of
D4a (`docs/lhs1140b_stationary_D4a_20260918.md`). 2026-09-18 (KST). Every
number is MEASURED with the suites and the builds named in sections 2
and 3 unless it is marked READ. No tolerance was moved: `ieq_res_tol`, `chem_cycle_tol` and the
fixture's 1e-20 presence rule stand at the values D4a read from source.

D4a settled which of the two outcomes D4 left open holds: production's
closure did not exit before its promise was met, and the failing test row
compared an abundance movement with a residual tolerance that bounds it
only twelve decades away. Three findings sat beside that verdict. D4b acts
on all three.

---

## 1. What changed

### 1.1 A refused trial puts back everything the next operation reads

`relax_photochemical_composition` takes each trial on the live state and
undoes it when the trial is refused. Its refusal branch put back the
composition, the frozen cell state and the primitive pair, and left the
rate state the closure evaluator reads, the caloric maps, the base ghost
and the carrier subsystem as the refused trial wrote them. The
consequence D4a measured: asking for the chemical residual of the state
the pass hands back gives 5.33e-4 where the state's own is 1.61e-16.

The enumeration is now one type with one save and one restore:

- `ieq_rate_state` with `save_ieq_rate_state` and `restore_ieq_rate_state`
  in `ionization_equilibrium.f90`, beside the declarations of the arrays
  they hold (`ieq_rate_cell`, `ieq_ne_cell`, `ieq_TK_cell`,
  `ieq_ntot_cell`, `ieq_met_coef`, `ieq_nonroot_streak`,
  `ieq_rates_ready` and the three stored sizes), so that an array added to
  the rate state is added to the snapshot in the same place.
- `thermochemical_state` with `save_thermochemical_state` and
  `restore_thermochemical_state` in `diffusive_photochemistry.f90`: the
  composition, `p`, `T`, `heat`, `cool`, `eta`, `bg_cell`, the rate state
  above, `dp_bc`, `n_part_cell1` and the carrier checkpoint
  (`carrier_checkpoint_take` / `_restore`, which already enumerate the
  carrier module arrays).

The caloric maps are not copied, they are rebuilt: the arrays are private
to `caloric_eos` and its own header states the rule (a path that installs
a different composition must refresh them BEFORE it maps an energy density
to a pressure), so the restore reinstates the composition and then calls
`get_species_densities`, the one routine that turns `(rho, f_sp)` into
number densities and refreshes them. The base ghost pair is written back
after that call, because `get_species_densities` writes `n_part_cell1`
itself.

### 1.2 The test asserts the returned state, not the movement of a sweep

In `src/tests/carrier_retry/carrier_retry.f90` the row
`one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
is replaced by
`every_species_row_of_the_returned_state_is_a_root_at_the_returned_temperature`:
the normalized reaction residual of every physical cell, evaluated with
`ionization_closure_residual_profile` at the returned temperature and
electron density with the transported fractions held by their ownership,
against `ieq_res_tol`. A companion row asserts that the profile is
available at all (the evaluator refuses a state whose element totals have
moved away from the rates').

The abundance movement of one further sweep stays, as a line that gates
nothing, printed with the binding species and cell AND with the
conditioning of that row: the sensitivity of its normalized residual to
its own unknown, `|R|/eps`, and the relative abundance a residual of
`ieq_res_tol` therefore admits.

### 1.3 The closure's header states what the code does

The header of `equilibrate_chemistry_at_fixed_conserved_state` named
`ieq_res_tol, 1e-6` as the tolerance the cycle stops on; the code tests
`chem_cycle_tol = 1e-5`. It also claimed that a temperature change below
that tolerance cannot move a composition beyond it. The header now states
the exit test, states that the composition is the last sweep's at the
temperature of the cycle before the exit while the temperature is the one
computed after it, and replaces the claim by what D4a measured, for
abundance and for residual separately.

---

## 2. The tests

### 2.1 The restore, item by item

`src/tests/carrier_retry/closure_returned_state.f90` (the D4a driver) now
gates its table 4. A refusal is forced through the closure's own cycle cap
(`chem_cycles_cap_for_test = 1`), which makes the closure spend its budget
on every trial; the pass ends on `carrier_relax_chemistry_refused` with no
step kept. The pass is run twice: the first is entered with the pressure
the column was seeded with through `comp_p_from_T`, so `p` and `T` come
back one caloric round trip (1e-14) away and that difference is not a
restore; the second is entered with exactly the pair a refused pass
returns, so there `p` and `T` are gated too.

Control: the ENTRY TEXT itself, compiled with this driver, which needs no
routine the entry text lacks.

| item | entry text | corrected |
|---|---|---|
| composition `f_sp` | 0 PASS | 0 PASS |
| frozen cell state `bg_cell` | 0 PASS | 0 PASS |
| rate state `ieq_rate_cell` | 4.723783293125725e+7 FAIL | 0 PASS |
| `ieq_ne_cell` [cm^-3] | 2.031037144688021e+7 FAIL | 0 PASS |
| `ieq_TK_cell` [K] | 1.487024919043597e+1 FAIL | 0 PASS |
| `ieq_ntot_cell` [cm^-3] | 4.723783293125725e+7 FAIL | 0 PASS |
| `ieq_met_coef` | 5.661042215266668e-7 FAIL | 0 PASS |
| stored layout and `ieq_rates_ready` | 0 PASS | 0 PASS |
| non-root persistence counter | 0 PASS | 0 PASS |
| caloric map, as `p` from the unchanged `rho e` | 1.071138215601453e-3 FAIL | 0 PASS |
| base ghost electron term `dp_bc` | 1.052540458382342e-6 FAIL | 0 PASS |
| base ghost particle count `n_part_cell1` | 1.804717559635849e-4 FAIL | 0 PASS |
| carrier subsystem matches | false FAIL | true PASS |
| pressure (second pass) | 0 PASS | 0 PASS |
| temperature (second pass) | 0 PASS | 0 PASS |
| heating, cooling, efficiency | 0 PASS | 0 PASS |

Nine rows RED on the entry text, all GREEN after. The entry-text numbers
are D4a's table 4 to every digit it printed. On the entry text the SECOND
pass measures zero throughout: the first pass left the module at a fixed
point of the defect, and the state the second pass is entered with is the
stale one it then reproduces. The row that fails on such a text is the
first pass's.

### 2.2 The returned-state residual, and that the row has teeth

`src/tests/carrier_retry/run.sh`, section 11, on the same column. This
driver calls the new save and restore pair, so it cannot be compiled
against the entry text; its control is the corrected text with
`restore_thermochemical_state` reduced to exactly what the entry text's
refusal branch put back (composition, frozen cell state, primitive pair),
built in the session scratch directory. That the reduction IS the entry
text is checked by the first row, which reproduces D4a's table 2 reading.

| row | entry text | corrected |
|---|---|---|
| `every_species_row_of_the_returned_state_is_a_root_at_the_returned_temperature` | 5.332326643642249e-4 FAIL | 4.440892098500626e-16 PASS |
| `the_diagnostic_sweeps_left_the_returned_state_as_they_found_it` | 5.290941096322854e-4 FAIL | 0 PASS |
| `and_that_trial_is_undone_bit_for_bit` | 5.291399172167668e-4 FAIL | 0 PASS |
| `the_returned_state_residual_row_fails_on_the_rate_state_of_another_trial` | 7.797539205884153e-5 PASS | 7.797539205884153e-5 PASS |

The RED reading 5.3323e-4 is D4a's table 2 number for the same state
(5.332326629e-4 in the D4a driver, whose sequence enters the pass with
another pressure; D4a section 2 measures that history dependence).

The last row is the sensitivity of the gate, built from production
routines: one transport interval is taken from the returned composition
and its chemistry is closed at the fixed conserved state, which is what a
trial of the relaxation does, and the residual of the RETURNED state read
against the rate state that trial leaves is 7.80e-5, 78 times
`ieq_res_tol`. It is a row that passes by failing to be small, and it is
what makes the gate above a statement and not an accident. The trial is
then undone with `restore_thermochemical_state` and the profile re-read:
the same bits.

### 2.3 The movement, and why it is not a gate

Printed, gating nothing (MEASURED on the corrected text):

```
one further sweep moves f_sp column 37 at cell 1 by 1.1005e-06;
the sensitivity of that row to its own unknown is |R|/eps 5.6749e-13,
so the sweep tolerance admits a relative abundance of 1.7621e+06
```

Column 37 is HeH+, 1.3e-12 of the gas in that cell. D4a's table 5 gives
5.456e-13 and 1.8328e+6 for the same row on its own sequence; the two
differ by the history dependence of section 2 of that memo and not in
their reading. A movement of 1.1e-6 against a tolerance that bounds the
same movement at 1.8e+6 is twelve decades from being a statement.

---

## 3. Impact

Control: a build of the entry text. Measured: the same text with the two
changed production files swapped in, so that the comparison carries this
item's change and no other worker's. Both single-threaded, on scratch
copies of the case directories; nothing under `backup/regression/` was
written and the tree's `EXHALE.x` was not touched.

### 3.1 The builds

| build | source | md5 of `EXHALE_D4b.x` |
|---|---|---|
| control | the working tree's `Makefile` + `src/` copied 2026-09-18 15:44 KST (repository HEAD `3c73905`) | `772c0fee25a54048197cda03f2af7de3` |
| measured, `wasp_full` and `carrier_model_a_newton` | the same, with `ionization_equilibrium.f90` and `diffusive_photochemistry.f90` replaced by this item's | `13f8c3c7ea203e26a8d74d90f6f23eea` |
| measured, the two 12000-step molecular cases | the same two files at an earlier point of the item: the delivered text minus one comment block later added to the type header, no executable statement differing | `559503ca008f3fce05ba41bd84c2502d` |

Compiler: conda-forge gfortran 16.2.0 (`/opt/miniconda3/bin/gfortran`),
bare `make`, `-O3 -fopenmp`, OpenBLAS of the compiler's own prefix. Every
run `OMP_NUM_THREADS=1`.

### 3.2 Key-off: the three regression cases

`relax_photochemical_composition` is called from ONE place, the outer
Picard loop of the stationary solver (`src/EXHALE_main.f90`, the
`outer_ending .eq. outer_running` branch), READ from the source. The three
cases below are step-capped marching runs, and none of their logs carries
a `carrier relaxation` or an `outer pass` line, so the branch this item
changes is not entered in any of them; the comparison measures the rest of
the binary. The molecular pair additionally writes the advection-corrected
profiles and the two breakdown files, and those are identical too.

| case | steps | Hydro_ioniz | Ion_species | `*_adv` |
|---|---|---|---|---|
| `wasp_full` | 400 | identical | identical | identical |
| `mol_diffusion` | 12000 | identical | identical | identical |
| `mol_carrier` | 12000 | identical | identical | identical |

"identical" is every data line; the `# provenance` comment carries the
wall-clock time of the run and differs.

### 3.3 Key-on: `carrier_model_a_newton`

`EXHALE_OUTER_PASSES=3 EXHALE_JFNK_MAXIT=5 EXHALE_CARRIER_DEBUG=1`, from
the pinned IC of `backup/regression/carrier_model_a_newton`.

**It does take refused trials**: 31 `REFUSED on the movement bound` lines
in the log of each build, and the pass of the third outer iteration ENDS
on one (`carrier relaxation ended on the movement bound, with the last
step inside it`, the last trial refused at `grow = 9.766e-4`), so the
branch this item changes runs and runs last.

**Movement: none.** `Hydro_ioniz.txt` and `Ion_species.txt` are identical
line for line apart from the provenance time, and the two run logs differ
only in their wall-clock columns (6 lines of 927). The outer-pass
residuals agree digit for digit (pass 1 mass 2.33e-2, momentum 7.86e-3,
energy 8.05e-1; pass 3 mass 9.42e-3, momentum 3.24e-3, energy 2.60e-1),
and so does the JFNK sequence.

The reading: what the refusal branch used to leave behind is rebuilt by
the caller before anything reads it on this configuration. The caller
re-derives the primitive state from the conserved variables and refreshes
the composition (`U_to_W` then `get_species_densities` then
`comp_T_from_p`, `src/EXHALE_main.f90`), the next trial's transport step
rebuilds the carrier background from `f_sp`, and the next chemistry
closure rebuilds the rate state and the caloric maps from the composition
it is handed. The state that IS read stale is the one a certification or
a closure-residual evaluation asks for straight after a pass that ended on
a refusal, which is what section 2.2 measures and what the fixture in
`src/tests/carrier_retry/` exercises.

That the leftover is inert on these two states is a measurement of the
states, not a reason to leave the branch as it was: a composition handed
back beside another state's rates, heat capacity and boundary is wrong
whether or not a run reads it before it is rebuilt.

### 3.4 Every suite whose driver links the two changed modules

Selected by running `src/tests/physics_probe/source_closure.py` on every
test program and keeping those whose closure contains
`diffusive_photochemistry.f90` or `ionization_equilibrium.f90`.

PASS, no failures: `acceptance_classes` (22), `adv_static_limit` (53),
`attempted_step` (70), `carrier_boundary_jacobian` (12),
`carrier_constraint_attribution` (10), `carrier_helium_inventory` (15),
`carrier_outer_boundary` (16), `carrier_reference_scales` (16),
`carrier_retry` (156), `carrier_returned_state_acceptance` (36),
`certification` (84), `element_operator` (40), `grid_and_gates` (279),
`ionization_imposed_fractions` (48), `ionization_stage_flux` (40),
`krylov_and_dogleg` (334), `species_masses` (9), `steady_species_rows`
(197), and the D4a driver `closure_returned_state` (34).

`coupled_block_jacobian` was not run: it is driven by a binary and reads
the repository's LHS 1140 b case, and its two failing rows are the
standing L22 coupled-block item, unchanged by anything here.


---

## 4. Noticed outside the item

1. **`grid_and_gates` fails one row when it is run against the tree's own
   `EXHALE.x`.** `base_grid_resolved_record_present` reads
   `missing_rows` for `base_cell_width_Rp`, `base_uniform_cells` and
   `base_cell_width_source`, because the suite's default executable is
   `$ROOT/EXHALE.x`, which in this shared working copy is another item's
   build and is behind the sources. With `EXHALE_EXE` pointing at a build
   of the sources under test the case passes and the suite is 279 of 279.
   Not fixed: the tree's binary is not this item's to rebuild.

2. **`docs/session_handoff_20260918.md` lists the replaced row as an open
   item** ("`src/tests/carrier_retry/`'s knife-edge row"). D4b closes it.
   Not edited: the handoff is the session's record and not this item's
   file.

3. **`docs/lhs1140b_stationary_D4a_20260918.md` describes the refusal
   branch as it stood when D4a measured it.** It is left as it is on
   purpose: it is the control this item's RED numbers are measured
   against, and its provenance section pins the snapshot it was taken on.

4. **`ionization_closure_residual_profile` refuses a state whose element
   totals have moved more than 1e-8 relative from the rates'** ("a
   different gas, whose closure these rates cannot state"). That refusal
   is why the new gate carries a companion row asserting the profile is
   available: a returned-state residual that is silently unavailable would
   otherwise read zero and pass.

5. **Restoring `ieq_nonroot_streak` is a choice, and it is stated at the
   type.** The counter exists to stop a run RESTING on a non-root, so the
   sweeps of a trial the caller discarded are not sweeps of the state the
   run carries and do not count towards that persistence. On this fixture
   the difference is zero either way, so no test separates the two
   readings.
