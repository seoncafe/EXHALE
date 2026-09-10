# Phase A: elemental invariants, the oxygen-carrier double count, and the
# deterministic tests

Date: 2026-09-03.  Working copy only; nothing under
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00` was modified.
The tree was copied to the scratch directory at 2026-09-03 18:16:22 KST
(`src/` md5 of the sorted per-file md5 list: `5c69f0fbf040ab7220a81f50a21cbb30`;
`Makefile` md5 `6ef7ef75c6c8e51e10d963181ecc7b45`).

This covers item A1 (centralized elemental accounting with debug assertions),
item A3 (focused deterministic tests) and item 12.8 (the oxygen-carrier double
count) of `docs/open_defects_20260903_review.md` section 6, plus what could be
established for items 12.4 and 12.6.

---

## 1. A1: one elemental accounting routine

New module `src/modules/functions/element_census.f90`.

| routine | what it states |
|---|---|
| `element_nuclei_and_charge(rho, f_sp, n_nuc, n_chg, rho_comp)` | nucleus density of H, He and the ten metal elements, the positive charge density, and the mass density the composition implies |
| `element_census_take(label, rho, f_sp, snap)` | record that state under a label |
| `element_census_verify(snap, rho, f_sp [, rho_is_fixed])` | compare a later state with the record |
| `element_census_reservoir(label, rho, f_sp)` | the absolute statement against the reservoirs the run resolved (`melem_ab`, `HeH`, `rho/mass_per_H`) -- the live-state form of `src/utils/element_budget.py` |

**No multiplicity is written in the new module.**  It reads `bsp_nH`,
`bsp_nHe`, `bsp_nO`, `bsp_nC`, `bsp_charge`, `bsp_mass`,
`bsp_is_excited_level`, `mion_elem`, `mion_stage` and `melem_A` from
`species_table`, so H2 and H2+ carry two H nuclei, H3+ three, HeH+ one H and
one He, OH one H and one O, H2O two H and one O, CO one C and one O, each
metal ion stage one nucleus of its element, and the He 2^3S column is skipped
because it is a level of He I.  Adding a carrier is adding a row to
`species_table` and this module follows.

### What it consolidated, and what it did not

`utilities.f90::hydrogen_helium_nuclei_density` is the existing single
definition of the H and He nucleus totals and it is **left as it is**: it is
called from the equilibrium sweep, the heating dump and the advection
post-process, its statement order is the floating-point add order the goldens
were snapshotted with, and it has no metals.  The new routine is the C/N/O and
metal counterpart it never had; the two are checked against each other by test
E1 through `calc_ne` (the charge the census returns is the electron density
`calc_ne` builds from the same state).  Merging them would move the add order
inside `hydrogen_helium_nuclei_density` and break byte-identity for no gain,
so it was not done -- this is the one place where two routines count the same
nuclei, and the test is what keeps them from drifting.

`src/utils/element_budget.py` keeps its own species list because it reads the
finished output by column label; `element_census_reservoir` is the in-code
statement of the same three closures.

### Where the assertions sit

Off unless `EXHALE_ELEMENT_ASSERT` is `1` (report and continue) or `2` (report
and stop); when off, `take` and `verify` return before allocating anything.
That follows the code's existing convention (`EXHALE_CARRIER_DEBUG`,
`EXHALE_DIFFUSION_CHECK`, `EXHALE_MAXSTEPS`); the tree has no preprocessor use
in any `.f90`, so a `-D` flag would have been a new mechanism.

| # | site | file | statement |
|---|---|---|---|
| 1 | around `carrier_write_back` | `diffusive_photochemistry.f90`, `photochemical_transport_step` | absolute nucleus density of every element, per cell |
| 2 | around the whole `photochemical_transport_step` | same | same (the write-back restores the entry element totals, so the step is element-neutral cell by cell) |
| 3 | around every `ioniz_eq` sweep, tagged marching / steady iterate / steady candidate | `ionization_equilibrium.f90` | every ratio `n_El/n_H`; the mass closure is reported |
| 4 | around `relax_photochemical_composition` | `diffusive_photochemistry.f90` | absolute nucleus density of every element |
| 5 | around the restart equilibration | `EXHALE_main.f90`, `equilibrate_loaded_composition` | every ratio `n_El/n_H` |
| 6 | the accepted state of each steady outer pass | `EXHALE_main.f90`, `steady_wind_with_element_diffusion` | against the resolved reservoirs |
| 7 | immediately before each final `write_output` (marching route and direct steady route) | `EXHALE_main.f90` | against the resolved reservoirs |

Site 3 covers "the `ioniz_eq` of the steady outer loop" and every other caller
in one hunk, which also keeps the diff out of `steady_newton.f90` and
`steady_residual.f90`.

The ratio, not the absolute density, is what is gated across `ioniz_eq`:
`n_io` is `intent(inout)` and the sweep rewrites it from its own `calc_rho`, so
all densities may move by one common mass-closure factor.  That factor is
reported separately, because "every element off by the same amount" is the
signature of a moved `n_H` denominator and "one or two elements off" is the
signature of a repartition defect, and the report tells them apart by printing
the per-element worst departure for all twelve elements.

### Tolerance

Default `1.0e-9`, overridable with `EXHALE_ELEMENT_ASSERT_TOL`.  These are
exact bookkeeping identities -- the operators rescale stage populations onto
element totals -- so the only error is floating point.  Measured accumulated
round-off on the marching states is `<= 1e-12` (this work measured 3.5e-16 per
operator on `mol_carrier`); the cell solver's own `xtol = sqrt(eps) = 1.5e-8`
is the largest departure a converged cell solve can leave.  `1e-9` sits between
them.  It is not a percentage of abundance.

### On failure

The element, the cell index, its radius, the before/after element and hydrogen
totals, the before/after `rho*n0` and `rho(composition)`, and every species
that carries that element in that cell with its multiplicity -- then the
per-element worst departure over the whole grid.

---

## 2. Item 12.8: the oxygen-carrier double count

### Verdict

**Not reproduced in the current working tree.**  The assertions were run
through the A2 oxygen configuration on the steady/Picard path that
`TO_BE_DONE.md` item (T) names, including 20 outer passes of
`relax_photochemical_composition`, and **no element budget broke**: every
`ioniz_eq` sweep (marching, steady-iterate and steady-candidate), every carrier
write-back, every transport step and the restart equilibration closed inside
`1e-9`, and the finished output closes all five elements at `1.77e-10` through
`src/utils/element_budget.py`.

This is a *not reproduced on the equivalent path*, not a *proved fixed*.  What
could not be done is stated in section 5.

### What was run

Configuration: `examples/18_oxygen_chemistry/input.inp` + `metals.inp`
(HD 209458 b, `Molecular chemistry: True`, `Oxygen chemistry: True`, and hence
`Molecular carrier transport` on by the `thereis_oxychem` default), restarted
from that example's own step-12000 output as the `_IC` pair, driven straight
into the steady solver with `EXHALE_PTC=1 EXHALE_PTC_JFNK=1`,
`OMP_NUM_THREADS=4`, `EXHALE_ELEMENT_ASSERT=1`.

* Run 1 (`$S/run_a2`): the JFNK returned `info = 2` and the outer loop left
  before the Picard pass, so 42 steady-iterate sweeps and 1230 probe sweeps
  were exercised but no carrier relaxation.  No assertion fired; output budgets
  `2.6e-13`.
* Run 2 (`$S/run_a2b`): the same, with a scratch-only diagnostic
  (`EXHALE_FORCE_CARRIER_PICARD=1`, **not** part of any delivered patch and
  removed from the tree afterwards) that lets the outer loop take its carrier
  pass from an unconverged solve, so the refresh under suspicion actually runs.
  20 outer passes, each with one `relax_photochemical_composition` step and its
  three `ioniz_eq` sweeps.  No assertion fired.  Output budgets: H, He, C, N, O
  all `1.768e-10` at `r = 1.0129` -- the same number for every element, which is
  the moved-`n_H`-denominator signature at an amplitude four orders below the
  reported `5.18e-2` and eight below the reported carbon excess.

### Ranked reading of why

The suspects were followed in code and none of them creates a nucleus:

* `carrier_state` (`diffusive_photochemistry.f90:612`) assembles `nH_free`,
  `nO_free`, `nC_free` as the full element totals from `f_sp` through the
  `bsp_n*` table plus the metal ion stages; `carrier_write_back`
  (`:1866`) rescales the free stages onto exactly those totals less what the
  carriers took.  The pair is element-neutral by construction and test E2
  measures it at `1.9e-16`, including the CO-collapse configuration section 111
  reports.
* `ioniz_eq` (`ionization_equilibrium.f90:583-633`) builds `nm_el` (the
  conserved element) as the stage sum plus the carriers and `nm_tot` as the free
  family after CO is removed, and the write-out at `:1808-1827` puts the stages
  back on `nm_tot` and subtracts OH/H2O from the O I column.  Entering and
  leaving, the convention is the same object, so the sweep is a fixed point in
  the element sense.
* The restart oxygen seeding (`load_IC.f90:466-512`) takes the carriers out of
  the ion stages of the same element.

The likeliest explanation is that the defect was removed by the carrier rework
of `docs/Update_EXHALE_stage1.md` section 128 (2026-09-02), which is *after* item (T)
was opened (2026-08-31) and which explicitly split the two disagreeing
definitions of the hydrogen the carriers may hold
(`hydrogen_available_to_carriers` versus the element total).  Section 128's own
note says `carrier_source` and `limit_to_element_budget` "had two different
definitions of the hydrogen the carriers may hold, and the looser one let them
take" nuclei the frozen stages were holding.  That is the same class of defect
as (T), in the same routines.

### What is left over, and is real

The steady/Picard path closes at `1.77e-10` where the marching path closes at
`2.6e-13`, and the departure is identical for all five elements -- so it is the
hydrogen denominator, not the elements.  It is three decades inside the `1e-8`
gate and eight decades inside the reported defect, and it is reported here
rather than chased.

---

## 3. A3: the deterministic tests

### `src/tests/element_census_tests.f90` -- `make element_census_tests`

Follows the `diffusion_tests.x` convention (standalone executable, explicit
source list, synthetic columns set by the driver, `stop 1` on failure).  It
links every module the production binary does except the main program, the way
`cce_probe.x` does, because the carrier write-back reads the ionization sweep's
cell state.  **20 tests, 20 passing**:

| id | statement | measured |
|---|---|---|
| E1a-d | `element_nuclei_and_charge` reproduces the species-table stoichiometry on a hand-computed state | exact (H 13, He 2, O 4, C 2 nuclei per unit density) |
| E1e-f | the positive charge it returns is the electron density `calc_ne` builds from the same state | exact |
| E2a | element closure through one carrier write-back, generic perturbation | 0.0 |
| E2b | element closure when the transported CO collapses by 2.5 decades -- the configuration item 12.8 reports | 1.9e-16 |
| E3a | H nuclei per H2 photoionization, over 15-200 eV of the Chung et al. (1993) branching | 2 to 1e-14 |
| E3b | net charge per H2 photoionization including the released electron | 0 to 1e-14 |
| E3c | the Lyman-Werner channel H2 -> H + H | closes by construction |
| E3d-e | the branching stays in [0,1] and is held at its 124 eV value above the measured range | 0.20477 |
| E4a-c | with no handoff the base H2 fraction comes from the chemical-equilibrium fit at `(p_base, T0)` and the ghost H2 nucleus fraction is its `mu_mixture` conversion | q = 0.838401, x2 = 0.984428 at He/H = 0.0793 |
| E4d-e | with a handoff the handoff value is used and the base composition is imposed | 0.84 |
| E4f-g | the fit is attainable at He/H = 0.0793 and exceeds the ceiling `0.5/(0.5+He/H)` at He/H = 0.2, which `input_read` refuses | 0.8384 against 0.7143 |

E3 also **reports** the two channels the branching does not resolve, rather
than assuming them small: double ionization is inside `f_di` and is about 4.3%
of the H2 ionizations above 80 eV with its co-product mis-assigned, and neutral
dissociation over 33-41 eV (up to 7% of `sigma_H2` at 37.5 eV) is handed to the
H2+ branch by `1 - f_di`, which creates one H2+ and one electron the event does
not make.  These are item 12.1's known approximation; the gate for them belongs
with the explicit-channel cross sections of phase E1, not here.

The test needed two routines made public in `diffusive_photochemistry`
(`carrier_state`, `carrier_write_back`), the same way `carrier_slope` and
`n_carrier_max` are already public for the order-of-accuracy test.

### Restart round trip -- `backup/regression/roundtrip`, rebuilt

The fixture item 12.4 asks for, built from a **short deterministic run** rather
than from an expensive converged product:

1. run the case cold for `EXHALE_MAXSTEPS` steps at `OMP_NUM_THREADS=1`, keep
   `output/Hydro_ioniz.txt` and `output/Ion_species.txt` as the reference;
2. hand those two back as the `_IC` pair and reload with `EXHALE_DUMP_IC=1`,
   which writes the loaded state before any equilibrium sweep and stops;
3. `check_roundtrip_state.py`: every species column of `Ion_species` matched
   **by label** to `1e-12`, `rho`/`v`/`p` of `Hydro_ioniz` likewise, and the
   `# coupling:` header preserved.

It replaces `backup/regression/check_roundtrip.py`, which compares by column
position and whose case `backup/regression/roundtrip` cannot run at all
(`output/Hydro_ioniz_IC.txt` is zero bytes, there is no `Ion_species_IC.txt`,
its logs still say `ATES_DUMP_IC` and `(ATES_main)`, and it is in no
`DEFAULT_CASES` list) -- item 12.4 confirmed.

**The case is rebuilt in place** (coordinator's direction, 2026-09-03), in two
halves that need no coordination: stage A is an ordinary named case -- 40
deterministic steps of the hot-Uranus molecular gate, `input.inp` + `base.inp`
+ `maxsteps`, so `run_check.sh check roundtrip` bitwise-compares it like any
other -- and stages B/C are `roundtrip/roundtrip_check.sh`, run after it, which
restores `input.inp` and the stage-A outputs on exit so `golden roundtrip`
still snapshots the right state. It stays OUT of `DEFAULT_CASES` and is
described in the `run_check.sh` header, because its second half is not a
bitwise comparison and needs its own invocation. The retired contents,
`check_roundtrip.py` among them, go to `roundtrip/.superseded_20260903/` rather
than being deleted, since `backup/` is not in the git remote.

Exercised on a scratch replica of `backup/regression/` (`$S/fakeroot`): the
installer runs `--dry-run` and for real, `run_check.sh` still parses,
`check roundtrip` reports `PASS ... (data identical)` against its own golden,
`roundtrip_check.sh` passes, and the case directory is left with no `_IC`
leftovers and `input.inp` bit-identical. Run on a binary WITHOUT the 12.4
patch, the same case fails on all four molecular columns and on the coupling
header -- so it guards exactly the two defects.

**Run on `backup/regression/mol_base_handoff` at 40 steps, the fixture found
two real defects and now passes.**  Both are in the delivered patch
`12.4_restart_roundtrip.patch`:

1. **`EXHALE_DUMP_IC` reported a molecular restart as atomic.**
   `write_output` takes the H2 / H2+ / H3+ / HeH+ / OH / H2O / CO columns from
   `nmol_eq` and `nox_eq`, which only the equilibrium sweep fills; the dump hook
   runs before the first sweep, so all four molecular columns came out exactly
   zero against the 5.34e12 cm^-3 of H2 the restart file carried.  Fixed with
   `molecular_carrier_densities_from_state(rho, f_sp)` in
   `ionization_equilibrium.f90`, called from the dump hook.
2. **`EXHALE_DUMP_IC` wrote the wrong coupling header.**  The rule that adopts
   the file's `# coupling:` state (`sec_ion_active = ...` in `EXHALE_main.f90`)
   runs about 280 lines *after* the dump hook, so the dump wrote
   `sec_ion=F sec_ion_step=-1` over a file that stated `sec_ion=T
   sec_ion_step=40`.  Fixed by applying the same rule inside the dump block.

Both are defects of the diagnostic, not of `load_IC`: the species fractions and
the parsed coupling fields were restored correctly all along.  They are the
reason the round trip could not previously be tested for a molecular case.

### The two halves that are run-level, not module-level

* The ionization-sweep half of A3(a) and the "repeated `F(Y)` is independent of
  the preceding trial state" test of A3 are exercised by the runtime assertions
  with `EXHALE_ELEMENT_ASSERT=2`, not by the standalone driver: a synthetic
  column cannot reach `ioniz_eq` without the whole radiation setup.
* The ghost-composition half of item 12.6 (does the chemical-equilibrium branch
  give the ghost the fitted composition) is likewise run-level.  E4 covers the
  scalar half.

---

## 4. Verification

* **Byte-identity.**  A binary built from the untouched tree and a binary built
  from all three patches produce bitwise-identical `output/Hydro_ioniz.txt` and
  `output/Ion_species.txt` on `mol_base_handoff`, `mol_carrier` and `wasp_full`,
  300 steps each at `OMP_NUM_THREADS=1`.  This is the scope the changes reach:
  the new module is inert unless the environment variable is set, the assertion
  calls return before allocating anything when it is not, and the two restart
  fixes are inside the `EXHALE_DUMP_IC` block, which ends in `stop`.
  `make check` itself was not run -- the harness lives in `backup/`, which is
  outside the scratch copy.
* **Cost.**  `mol_carrier`, 300 steps, `OMP_NUM_THREADS=1`: 11.96 s with the
  assertions off, 11.77 s with `EXHALE_ELEMENT_ASSERT=2` (898 checks).  Within
  noise.
* **Assertions actually fire on the real path.**  The same run reports every
  operator closing at 3.5e-16 to 6.2e-16.
* **Patches apply, including to the tree as it moved under this work.**  All
  three apply with `patch -p1` in the order A1 -> 12.4 -> A3 to a fresh copy of
  the 18:16 tree, and the result rebuilds and reproduces the staged tree
  exactly.  The tree was edited by another worker at 19:00 -- `EXHALE_main.f90`,
  `input_read.f90`, `write_output.f90`, `write_setup_report.f90`,
  `utilities.f90`, `parameters.f90`, `steady_newton.f90`, `steady_residual.f90`
  and `Makefile` -- and the three patches still apply to that revision with line
  offsets only, no rejected hunk (`$S/rebase`).  That build passes
  `element_census_tests.x` 20/20 and reproduces the unmodified 19:00 build's
  `Hydro_ioniz.txt` and `Ion_species.txt` on `mol_carrier` and `wasp_full` at
  200 steps, once the new `# provenance:` header (which carries the build
  timestamp and so differs between any two builds) is stripped the way
  `run_check.sh` strips comments.  Neither `ionization_equilibrium.f90` nor
  `diffusive_photochemistry.f90` -- where most of the A1 diff sits -- was
  touched by that revision.
* **Tests.**  `element_census_tests.x` 20/20; `roundtrip_fixture.sh` on
  `mol_base_handoff` at 40 steps passes all three statements.

---

## 5. What was not established

* **The section-111 failure was not reproduced.**  That state is a converged A2
  wind reached from a ~25000-step march plus a JFNK finish that returns
  `info = 0`, and it does not exist anywhere in the tree
  (`grep -rl "Oxygen chemistry: True" --include=input.inp` matches only
  `examples/18_oxygen_chemistry` and `examples/19_molecular_ir_bands`, and the
  example-18 output is the diverged step-12000 state, `du = 1939`).  Rebuilding
  it is hours of wall time.  What was exercised instead is the same *code path*
  -- JFNK probes, steady-iterate sweeps, the Picard carrier relaxation, the
  write-back and the restart equilibration -- from a different state.  A state
  that reaches the accepted A2 root is still what would settle item (T)
  outright, and the instrument to settle it in one run now exists:
  `EXHALE_ELEMENT_ASSERT=2` stops at the first broken refresh and names it.
* **The exact offending refresh was therefore not isolated**, and no fix for
  12.8 is delivered.  `patches/` contains no 12.8 patch.
* **`make check` was not run** (the harness is under `backup/`, deliberately
  outside the scratch copy).  Byte-identity was measured directly against a
  pristine build on three cases instead.
* **Item 12.6's ghost composition** was not measured; only the scalar branch is
  under test.  Adding a matrix case that omits the handoff -- the review's
  recommendation -- is still open.
* The `1.77e-10` common departure of the steady/Picard path against the
  marching path's `2.6e-13` is reported, not diagnosed.

---

## 6. Files

```
patches/A1_element_census.patch            A1: the census module, its call sites, SRC entry
patches/12.4_restart_roundtrip.patch       the two EXHALE_DUMP_IC defects the case found
patches/A3_tests.patch                     the test driver, its make target, two public::
patches/12.4_roundtrip_case_install.sh     rebuilds backup/regression/roundtrip (--dry-run first)
patches/run_check_header.snippet           the header block that installer inserts
roundtrip_case/                            what the installer drops into the case dir
docs/phaseA_report.md                      this file
docs/Update_EXHALE_section149_draft.md     section 149 draft
docs/TO_BE_DONE_T_draft.md                 replacement wording for item (T)
stage1/ stage2/ stage3/                    the trees the three patches were cut between
rebase/ rebase_ref/                        the same patches on the 19:00 tree, and its own build
run_a2/ run_a2b/                           the two A2 reproduction attempts
bi_*/ rb_*/                                the byte-identity pairs
fakeroot/                                  scratch replica used to exercise the installer
rt_test/ rt_test_unpatched/                the rebuilt case, passing and failing
```

Apply order, on the coordinator's signal (after section 147 and block H), each
re-based on the tree as it stands then and dry-run first:

1. `patch -p1 --dry-run < patches/A1_element_census.patch`, then apply
2. the same for `patches/12.4_restart_roundtrip.patch`
3. the same for `patches/A3_tests.patch`
4. `patches/12.4_roundtrip_case_install.sh <repo root> --dry-run`, then for real
5. `backup/regression/run_check.sh check roundtrip` and then
   `golden roundtrip` (first snapshot), then `check roundtrip` again
6. `( cd backup/regression/roundtrip && ./roundtrip_check.sh )`

Nothing has been written into the tree.
