# Session handoff, 2026-09-19: the series on `docs/PLAN_20260918_rev2.md`, items D1 to D9

Written 2026-09-19 at the close of the series. The previous handoff is
`docs/session_handoff_20260918.md` (the L27 to L37 series of
`docs/PLAN_20260917.md`); its outline is the outline of this one. The plan this
series executed is `docs/PLAN_20260918_rev2.md`, the revision of
`docs/PLAN_20260918_rev1.md` and of the first text `docs/PLAN_20260918.md` after
the two reviews `docs/PLAN_20260918_review.md` and `docs/PLAN_20260918_review1.md`.
Its table "Decisions made on 2026-09-18" carries a DONE or RESULT paragraph for
every decided item, and the update log records the series in
`docs/Update_EXHALE_stage2.md` section 13.

**Every number in this document is READ from the item memo, the item report,
the plan or the source line named beside it.** Where a worker labeled a number
MEASURED, it is quoted here as READ and the memo is named; the measurement is
that worker's. Two things were checked for this document and are said to be:
the md5 of the tree's `EXHALE.x` and whether `make -q` finds it current.

**The binaries.** The tree the series worked on is repository HEAD `3c73905`
plus the uncommitted L27 to L37 series; the closing binary of that series is
`74b96cdcf887dcfee301ec547f012d88`
(`LHS1140b/models/BINARY_MANIFEST_74b96cdcf887.txt`). Every item of this series
measured on private builds (`OBJDIR=build_<item>`, `EXE=EXHALE_<item>.x`) of a
frozen snapshot or of the live tree with its own files at their entry text, so
every control is that item's own entry-text build. **No closing binary was
built for this series**: the tree's `EXHALE.x` is stale against the sources
(section 2.3).

**GOLDEN REFRESH: none in this series**, by the user's decision of 2026-09-18
(the corrections continue on the same tree; the goldens stay the 2026-09-17
set; one refresh when the corrections have landed, after the reference states
are re-solved under the current boundary model). No worker wrote under
`backup/regression/golden/`.

---

## 0. What the series did, in one page

**Where it started.** The previous series ended with the plan's nine items:
five of the first text, corrected twice (D1 the default base cell width and its
key; D2 the base boundary at zero and reversed flow; D3 the transported He II row
and its acceptance norm; D4 the thermochemical closure's returned state; D5 one
boundary-state operation and the restart that reproduces it), one from the first
review (D6 the helium held in HeH+ not reserved from the stage simplex), one
from L36e extended by the second review's confirmed sign defect (D7 metal charge
exchange), and two from the L34c and L36e reports (D8 run generations and
continuation; D9 pass budgets, seed compatibility and the reference refresh).

**What the user decided on 2026-09-18.** No strict matrix on the closing binary;
no running table of movement against the goldens, each item measuring against
its own entry-text control; no golden refresh while the corrections continue.
On `docs/DECISION_D2a_D8_review.md`: D2a is **A** at the pressure-balanced
stationary contact (the reservoir owns the thermodynamic state of the level,
including at rest), the rest of the Codex D2a text kept; D8 is the Codex
generation structure as reviewed, with legacy import and D8a first. D3b was
approved as recommended (the He II gate and denominator stay, the stage-sum
entries gate at their executed bound, one bounded experiment on the movement
bound). The D5b-2 shape was approved as D5a proposed it. D4b was decided by
D4a's measurement.

**What landed.** All items but D1b and D9 step 3.

- **D1a**: the default base cell width is one named constant, its provenance a
  flag, and `EXHALE_resolved.out` states it at round-trip precision; no value
  moved.
- **D3a / D3b**: the He II row of the hot Uranus is stiff by three decades and
  arithmetically clean (net 3.3e11 above its cancellation floor); its refusal is
  an abundance error of 0.16 per cent; the gate and its denominator stay, and the
  two stage-sum entries now gate at `2 (5 nk + 2) eps g` (3.109e-15 H, 5.329e-15
  He at the measured g = 1).
- **D4a / D4b**: production's closure keeps its promise (the returned state is a
  root at 2.6e-14); the failing test compared an abundance movement with a
  residual tolerance twelve decades looser and is replaced by a returned-state
  residual gate; a refused trial now puts back the rate state, the caloric maps,
  the base ghost pair and the carrier checkpoint.
- **D5a / D5b-1 / D5b-2 / D5b-3**: the lower boundary is one operation whose
  inputs are the physical column, the prescribed reservoir, the radiation context
  and the model options; the ghost rows of a file are not read; the ghost's H2
  partition is solved against its own ionization; every evaluation site derives
  the boundary after the composition sweep. Consequence: the three molecular
  states whose certificates depend on the old ghost are refused (stale).
- **D6**: the helium ion stages are bounded by the helium not held in molecules
  or in the metastable level, in every place the partition is formed; the proton
  ledger closes.
- **D7a / D7b / D7c / D7d**: the generic charge-exchange helium row has an
  explicit orientation; the transported stages carry the metal charge exchange
  through one reaction-source interface; every ground-state helium reaction of
  Table 4 reacts from the ground singlet.
- **D8a / D8b**: the runner classifies an ending before it continues; a case
  publishes immutable generations through one atomic index, `latest_complete` and
  `latest_certified` apart; the catalog was imported without recomputation.
- **The bounded experiment**: no item on the movement bound.
- **D9 steps 1 and 2**: a pass and wall ceiling by configuration class, read off
  recorded solves; seed compatibility before distance in every tier.
- **D2b**: the lower boundary is a contact, acoustically matched, then upwinded
  with no width; the reservoir owns the level at inflow and at rest; conduction,
  element diffusion and carrier transport each carry their own base condition.
  One departure from the decision text is for the user (section 4.1).

**What is left.** The D2b direction rule (for the user), D1b, the retention
policy of the 84 legacy generations that carry no verdict, D9 step 3 (re-solve the
reference states under boundary model v2, then the one golden refresh), and the
technical items of section 4.2.

---

## 1. The items, one by one

### D1a: the base cell width disclosed

READ, `docs/lhs1140b_stationary_D1a_20260918.md`. `dr_base_default` in
`parameters.f90` (the exact decimal `1.9999999494757503d-4`, bits
`3F2A36E2E0000000`), `dr_base_from_key` set by the reader,
`round_trip_decimal` (`ES26.17E3`) in `write_setup_report.f90`, seven new rows
of `EXHALE_resolved.out`, a `## The grid` section in `REPRODUCE.md`, row K33b of
`docs/input_schema.md`. The key stating the resolved value rebuilds the grid
bit for bit; `2.0e-4` in the key moves the centers by 3.818119e-09 on the
`roundtrip` grid. `wasp_full` and `hydrostatic_column` byte-identical.

### D3a: the transported He II row, channel by channel

READ, `docs/lhs1140b_stationary_D3a_20260918.md`. 26 channels (32 after D7b)
from both kernels behind `EXHALE_STAGE_CHANNELS`; the ledger closes to 3.0e-16.
At cell 248 of `carrier_model_a_newton`: P and L 7.6e2, net 1.2153, floor
3.717e-12, the residual an abundance correction of 1.58e-3 of x(He II) (2.35e-5
on the locally equilibrated state). The present scale admits 6.0e-13 to 2.1e-2
of the stage fraction over one column, the turnover scale 2.0e-5 at every stiff
cell. The recorded normalized order 2.886 is a dimensional 1.923. New suite
`stage_row_balance`.

### D3b: the stage-sum entries gate

READ, `docs/lhs1140b_stationary_D3b_20260918.md`. `ionization_stage_sum_entry`
in `certification.f90` gates at `ionization_stage_sum_rounding_bound(nk, g)`
with `g` measured at the face that carries the measure; anchor (8). The three
production states read 0.057 to 0.079 of the bound; broken constructions ten
decades above. `certification` 103 of 103.

### D4a: the returned state measured

READ, `docs/lhs1140b_stationary_D4a_20260918.md`. Residual of the returned state
at the returned temperature 2.593e-14 over the column (5.99e-19 on HeH+); the
HeH+ row admits a relative abundance of 1.8328e+6 at `ieq_res_tol`. Three
findings: the refusal branch left module state behind (5.33e-4 where 1.61e-16 is
right), the header named the wrong tolerance, and its temperature-to-composition
claim is false for abundance.

### D4b: the returned-state contract

READ, `docs/lhs1140b_stationary_D4b_20260918.md`. `ieq_rate_state` and
`thermochemical_state` with save and restore; nine restore rows RED then 0;
`every_species_row_of_the_returned_state_is_a_root_at_the_returned_temperature`
5.3323e-4 on the reduced control and 4.44e-16 delivered; the header corrected.
No tolerance moved; key-off identical; `carrier_model_a_newton` identical.

### D5a: where the boundary comes from

READ, `docs/lhs1140b_stationary_D5a_20260918.md`. Diagnosis only. Closure-stage
mixing worth 6.92 of the 15.98 floors of the molecular reference state; the
ghost composition determining the boundary (9.06 to 7.4e8 floors over four
admissible ghosts); `ntot_bc + dp_bc` two quantities under one name; the face
cache worth up to 3.2 floors; serialization innocent at 6.1e-16; a 2-cycle of
1.9 floors that no boundary operation removes.

### D5b-1: the boundary after the sweep

READ, `docs/lhs1140b_stationary_D5b1_20260918.md`. Three evaluation routes of
`EXHALE_main.f90` apply `Apply_BC` after the sweep; `boundary_state_trace.f90`
and `src/tests/boundary_state/` new. Molecular rows 15.98, 31.24, 19.04 floors
to 9.06, 3.09, 0.08, all three then CERTIFIED; atomic unchanged.

### D5b-2: one boundary-state operation

READ, `docs/lhs1140b_stationary_D5b2_20260918.md`. The operation's contract in
`Apply_BC.f90`, the reservoir version and model identity in `base_boundary.f90`,
the ghost H2 closure (tolerance 1.0e-10, 30 passes, `error stop` on failure) in
`ionization_equilibrium.f90`, the lower ghost rows not read by `load_IC.f90`, the
boundary headers in `write_output.f90`. Acceptance all bitwise. The three
molecular states now read 19.97, 110.0 and 27.13 floors: NOT CERTIFIED, the
ghost's trace ions moving by 3e-3 to 5e-2 inside the sweep's acceptance band.

### D5b-3: the four rows and the fourth site

READ, the D5b-2 memo section 11 and the D5b-2 report's D5b-3 part. The periodic
certification derives its own boundary (recomputes 19 and 21 to 0); three
round-trip rows compare the physical cells; the stagnation budget 100 passes;
`EXHALE_resolved.out` labels the prescribed counts. `grid_and_gates` 285 of 285.

### D6: helium reserved from the stage simplex

READ, `docs/lhs1140b_stationary_D6_20260918.md`. `n(He II) + n(He III) <=
n_He,nuc - n_He(molecules) - n_He(2^3S)` in the projection, the limiter, every
trial, the difference step, the source and the write-back; the singlet not
clipped; the proton ledger closed. New suite `carrier_helium_inventory`: the
counterexample's inventory error 5.000e-2 to 1.110e-16. The constraint does not
bind on `hp_front` or `carrier_model_a_newton`.

### D7a: the helium row's equation basis

READ, `PLAN_20260918_D7a_report.md` in the scratch directory (its path is in section 5
below; D7a wrote no memo under `docs/`). A mandatory `he_row_sign` on `cx_add_to_fvec` and
`cx_add_to_jac` (-1 in the triplet system, +1 elsewhere), the cell-ownership
guard `cx_require_cell`, the interim refusal of transported ionization with
metals, new suite `charge_exchange_rows`. The review's counterexample RED to every
digit, GREEN at 0. `cx_full 1` on `wasp_full` moves He II by 1.817e-3.

### D7b: metal charge exchange in the stage sources

READ, `docs/lhs1140b_stationary_D7b_20260918.md`. `cx_reaction_rates` and
`charge_exchange_stage_sources` in `charge_exchange.f90`; the missing term was
1.813e-3 of the H II source; the transported and the row-implied sources now
agree to 2.2e-16 and exactly on the solved column. The refusal lifted. The key-on
`wasp_full` solve does not certify, with or without metals.

### D7c: the He <-> H pair from the ground singlet

READ, `docs/lhs1140b_stationary_D7c_20260918.md`. The two triplet systems passed
the summed He I; the barrier of the Table 4 rate is the ground state's.
`wasp_full` moves by 3.98e-6 at most (He I where it is 3.4e-7 cm^-3).

### D7d: group C and the advected pair from the ground singlet

READ, `docs/lhs1140b_stationary_D7d_20260918.md`. The advection system and the
three group C sites; `cx_full` off byte-identical; the certified atomic
fiducial's `_adv` profiles move by 2.43e-13.

### D8a: endings classified before a continuation

READ, `docs/lhs1140b_stationary_D8a_20260918.md`. `classify_ending` with seven
classes, the continuation for `hydrodynamic_refusal` alone, `solve_<n>/` and
`ENDING`, the evaluate pass in `eval/`, the campaign's `pipefail` and status
file. `tests/run_case_policy.sh` 16 of 16.

### D8b: generations and one index

READ, `docs/lhs1140b_stationary_D8b_20260918.md`. Section 3 below.

### The bounded experiment of D8 step 8

READ, `docs/lhs1140b_stationary_D8bound_20260918.md`. `EXHALE_CARRIER_TRUST_HOLD`
(default unset). `carrier_model_a_newton`: the bound never binds, held and halved
bitwise equal, the He II distance from the root falls 7.52 per cent a pass and
crosses zero. `wellmixed/HeH0.083`: the row rises at the held bound.
`kzz1e9/HeH0.083`: 8.15e-04, 6.29e-03, 1.93e-03 and `info = 2`, stopped at three
passes by a line-search collapse. No item.

### D9 steps 1 and 2: budgets and seed compatibility

READ, `docs/lhs1140b_stationary_D9_20260918.md`. Four configuration classes
(40 p / 30 m, 20 p for each solve / 45 m, 40 p / 6 h, 90 p / 6 h), authorized
overrides, UNCERTIFIED INCOMPLETE at a ceiling; `pick_seed.py` compatibility in
every tier; `element_flux_closure.py` sourcing the runner's policy. Fixtures 12 of 12,
10 of 10 and 10 of 10.

### D2b: the lower boundary as a contact

READ, `docs/lhs1140b_stationary_D2b_20260918.md`. Matching by C-, direction
after it, `w_rev` in {0,1}; model identity
`characteristic_face_ps_reservoir_C_minus_contact_upwind_v2`; conduction reads
the level's temperature (`conduction_base_level_T`, budget entry
`conduction_base_heat_flux`); element diffusion records
`element_base_diffusive_flux`. Movement below 1e-3 everywhere measured
(`lower_profile` 3.312e-04 the largest). `boundary_state` 35 of 35,
`grid_and_gates` 290 of 290 (the advisor's verification, READ from the plan).

---

## 2. The state of the tree

### 2.1 Source files changed in this series, and by which item

READ from the memos and reports; the file list itself is `git status --short`
and `git diff --stat`, run read-only for this document. **The tree also carries
the uncommitted L27 to L37 changes of the previous series** (their owners are in
`docs/session_handoff_20260918.md` section 2.1), so the `git diff` of a file below
is the sum of both series; a file changed by the L series alone is not listed.

| file | items of this series |
|---|---|
| `src/EXHALE_main.f90` | D5b-1 (three evaluation routes), D5b-2 (the joint test, the stationary restart's final certification, each marching attempt, the reports), D5b-3 (the periodic certification), the bounded experiment (`EXHALE_CARRIER_TRUST_HOLD`) |
| `src/modules/states/Apply_BC.f90` | D5b-2 (the operation's contract, the input-set tag, the recompute in `Rec_BC`) |
| `src/modules/states/base_boundary.f90` | D5b-2 (reservoir version, model identity v1, solved ghost counts), D2b (the contact law, model identity v2, `base_reservoir_temperature_at`, `EXHALE_BASE_BRANCH_ON_CELL1` removed) |
| `src/modules/states/boundary_state_trace.f90` | D5b-1 (NEW, 188 lines, diagnostic) |
| `src/modules/radiation/ionization_equilibrium.f90` | D4b (`ieq_rate_state`), D5b-2 (the ghost H2 closure, the write-back removed) |
| `src/modules/radiation/charge_exchange.f90` | D7a (orientation, guard), D7b (`cx_reaction_rates`, `charge_exchange_stage_sources`), D7d (header) |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | D6 (the helium reservation, the proton ledger), D4b (`thermochemical_state`, the closure's header), D3b (`g` at the measured face), D7b (the frozen metal stages, the stage sources, the atomic pair's singlet), D2b (comment: the carrier base condition) |
| `src/modules/functions/ionization_stage_transport.f90` | D6 (`xsum_max`), D3b (the rounding bound and its ceiling, `Ek_out`) |
| `src/modules/nonlinear_system_solver/ion_residual_core.f90` | D3a (26 channels), D7b (32 channels, `transport_operator`) |
| `src/modules/nonlinear_system_solver/System_HeH_mol.f90` | D3a (the channels of `mol_heh_rows`) |
| `src/modules/nonlinear_system_solver/System_HeH_TR.f90` | D7c |
| `src/modules/nonlinear_system_solver/System_HeH_TR_metals.f90` | D7a (-1), D7c, D7d |
| `src/modules/nonlinear_system_solver/System_HeH_metals.f90` | D7a (+1) |
| `src/modules/nonlinear_system_solver/System_HeH_mol_metals.f90` | D7a (+1), D7d |
| `src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90` | D7a (+1), D7d (group C singlet, the dead assignment, `constrained_network_balance_rows_of_cell`) |
| `src/modules/nonlinear_system_solver/System_implicit_adv_HeH_TR.f90` | D7d |
| `src/modules/time_step/certification.f90` | D3b (`ionization_stage_sum_entry`), the advisor's wording fix |
| `src/modules/time_step/viscous_conduction.f90` | D2b |
| `src/modules/functions/binary_element_diffusion.f90` | D2b (the Dirichlet statement, `element_base_diffusive_flux`) |
| `src/modules/files_IO/load_IC.f90` | D5b-2 (ghost rows not read, `base_reservoir_composition_row`) |
| `src/modules/files_IO/write_output.f90` | D5b-2 (the boundary headers) |
| `src/modules/files_IO/input_read.f90` | D1a, D7a (the refusal), D7b (lifted) |
| `src/modules/files_IO/write_setup_report.f90` | D1a (`round_trip_decimal`, the grid rows), D5b-3 (the reservoir labels), the advisor's label fix |
| `src/modules/init/parameters.f90` | D1a |
| `Makefile` | D5b-1 (one `SRC` line) |
| `src/utils/map_state_to_grid.py` | D8b (index resolution) |
| `src/utils/element_flux_closure.py` | D9 (the continuation sourced), D8b (reader resolution) |

Test files: `src/tests/boundary_state/` (NEW: D5b-1, D5b-2, D2b's
`contact_law_at_the_base.f90`); `src/tests/charge_exchange_rows/` (NEW: D7a, D7b,
D7c, D7d); `src/tests/carrier_helium_inventory/` (NEW, D6);
`src/tests/stage_row_balance/` (NEW: D3a, D7b, the advisor's comment fix);
`src/tests/carrier_retry/` (D4a's `closure_returned_state.f90` and runner, D4b's
`carrier_retry.f90`); `src/tests/certification/certification_contexts.f90` (D3b);
`src/tests/ionization_stage_flux/` (D6 case 4b, D3a);
`src/tests/ionization_imposed_fractions/imposed_ionization_rows.f90` (D6, one
`cx_set_cell`); `src/tests/grid_and_gates/`: `base_cell_width_provenance.f90` and
`base_grid_key_reproduces_default.sh` (NEW, D1a), `ionization_transport_atomic.sh`
(D7a), `restart_round_trip.sh` (D5b-3, the advisor's band), 
`output_state_consistency.sh` and `run.sh` (D5b-3), `carrier_movement_bound_hold.sh`
(NEW, the bounded experiment), `base_branch_discriminant.f90` and
`base_boundary_continuity_probe.f90` (D2b). Outside `src/`, and outside the git
remote: `LHS1140b/models/run_case.sh`, `run_campaign.sh` (D8a, D9, D8b),
`pick_seed.py` (D9, D8b), `status.py`, `write_reproduce.py` (D1a, D8b),
`publish_state.py`, `import_legacy_states.py`, `legacy_state_map.txt` (NEW, D8b),
`tests/run_case_policy.sh` (D8a), `tests/campaign_budget.sh`,
`tests/seed_compatibility.sh`, `tests/closure_continuation.sh` (D9),
`tests/state_generations.sh` (D8b); tracked: `LHS1140b/MODELS.md` sections 3, 4, 6
and 9 (D9, D2b, D8a, D8b). Documents: `docs/input_schema.md` (D1a, D5b-2, D2b),
`docs/certification_tolerance_anchoring_20260910.md` anchor (8) (D3b).

The D4a driver wrote its saved return, `D4a_returned_state.txt`, into whatever
directory the caller started from, which left copies at the repository root and in
`src/tests/carrier_retry/`. `run_closure_returned_state.sh` now runs the driver in
the suite's own object directory (34/0). The new output is byte-identical to both
copies, and both copies were removed on 2026-09-19.

**Added after this handoff was written (2026-09-19, the user's decisions on D2b, D1b and retention):**
the Mixed-grid default base cell width is now exactly `2.0d-4` R_p (item D1b,
`docs/lhs1140b_stationary_D1b_20260919.md`); every live input with stored
results carries `Base grid [dr,cells]: 1.9999999494757503e-4 50`, the frozen
trees carry `GRID_DEFAULT_NOTE.md`, and `make_models.py` writes the pinned line
for new catalog cases. A scratch copy of an input taken before the pinning
builds the new grid and will refuse catalog states: copy the input again. The
certification report now states whether the base contact's direction, where
the wind window decided it, agrees with the base face mass flux (D2b memo
section 6b). All generations are kept until D9 step 3 has re-solved their
replacements.

### 2.2 New keys and their defaults

**No input key of `input.inp` was added.** The environment variables and
constants the series added or changed, READ from the memos:

| key | default | item |
|---|---|---|
| `EXHALE_STAGE_CHANNELS` | unset: nothing reached, no file opened; a list of up to eight cells writes `output/stage_channels.txt` | D3a (`ion_residual_core.f90`) |
| `EXHALE_BOUNDARY_TRACE` | unset; set, writes `output/boundary_trace.txt` and prints the gap | D5b-1 (`boundary_state_trace.f90`) |
| `EXHALE_TRACE_EXPERIMENT` | unset; `boundary_from_the_entry_composition` reproduces the pre-D5b-1 order from the same binary | D5b-1 |
| `EXHALE_CARRIER_TRUST_HOLD` | unset (off); a positive value holds the composition movement bound | the bounded experiment (`EXHALE_main.f90`) |
| `EXHALE_BOUND_HOLD_PASSES` | the fixture's pass count (14) | the bounded experiment (`carrier_movement_bound_hold.sh`) |
| `base_face_mach_blend` | 1e-8, now the supersonic-outflow regularization at `M_i = -1` only (overridable by the existing `EXHALE_BASE_MACH_BLEND`) | D2b (`base_boundary.f90`) |
| `EXHALE_BASE_BRANCH_ON_CELL1` | REMOVED | D2b |
| `base_ghost_closure_tol`, `base_ghost_closure_passes` | 1.0e-10, 30 (constants) | D5b-2 (`ionization_equilibrium.f90`) |
| `base_reservoir_prescription_version` | 1 | D5b-2 |
| `RUN_CASE_POLICY_ONLY`, `RUN_CASE_POLICY` | unset; `RUN_CASE_POLICY_ONLY=1` defines the runner's policy and returns | D8a, D9 |
| `CAMPAIGN_WALL`, `models/budget_overrides.txt` | the class table; an override without `authorized:` is refused | D9 |

What changed in the contract without a key: `Ionization transport: True` with
metals was refused by D7a and is accepted again since D7b; state files carry
`# boundary_model` and `# boundary_reservoir` lines, and the lower ghost rows of a
restart file are not read (D5b-2); `EXHALE_resolved.out` carries the grid rows
(D1a) and the labeled reservoir counts (D5b-3); the base face state follows the
contact law of model v2 (D2b).

### 2.3 The binaries and the goldens

**The tree's `EXHALE.x` is stale against the sources.** Checked 2026-09-19 for
this document: md5 `7bcd4f59ded4c6e069f6ce897ad652e3`, modified 2026-09-18 16:48,
and `make -q` returns nonzero. That md5 is D3a's final measured build (READ, the
D3a memo) and D7b's control; it carries neither D3b, D7b to D7d, D4b, D5b-2,
D5b-3, D8's key nor D2b. The suites that default to `$ROOT/EXHALE.x` read it
unless `EXHALE_EXE` names a current build (D4b noticed
`base_grid_resolved_record_present` failing for that reason). The
`EXHALE.x` in `backup/regression/wasp_full/` (and in every other case
directory there) is not a stray: `run_check.sh` copies the binary under test
into each case directory before it runs it (line 287), so the file is the
binary of that directory's last check (corrected 2026-09-19; the first text
of this handoff called it stray).

Binaries of record for the catalog: `LHS1140b/models/EXHALE_3146d11b.x` (L34, the
md5 every `REPRODUCE.md` of 2026-09-18 names). The D8a and D9 fixtures run on it.

**No golden was refreshed.** `backup/regression/golden/` is the 2026-09-17 set and
no longer matches the tree on the `hp_*` cases (stale before this series,
`hp_front` 1.6e-1, `hp_zero_seed` 3.9e-3, `hp_trace_seed` 3.3e-3,
`carrier_model_a_newton` 1.2e-4) and on every case D7c, D5b-2, D5b-3 and D2b
moved (`wasp_full` by D7c; `mol_base_handoff`, `mol_carrier` by D5b-2, D5b-3 and
D2b; `lower_profile` by D2b). The table of every movement the memos measured is
in `docs/Update_EXHALE_stage2.md` section 13, "The golden refresh: still deferred
by decision". `wasp_full_newton`, `mol_lyman_werner`, `mol_ir_bands`,
`mol_sec_ion` and `oxygen_chemistry` were run by no item of this series.

---

## 3. The catalog

READ, `docs/lhs1140b_stationary_D8b_20260918.md` and `LHS1140b/MODELS.md`
section 9.

**The generation contract.** A state directory publishes each state as an
immutable `states/<generation_id>/` (the two state files, `certification.txt`,
`manifest.json`); `state_index.json` carries `latest_complete` and
`latest_certified` as separate fields; one publisher,
`LHS1140b/models/publish_state.py`, writes both under a lock and replaces the
index by an atomic rename on one filesystem; `runs/<run_id>/` holds the inputs of
one attempt; every reader (`run_case.sh`, `run_campaign.sh`, `pick_seed.py`,
`status.py`, `write_reproduce.py`, `map_state_to_grid.py`,
`element_flux_closure.py`) resolves the index once and reads both halves from one
generation. Generations are copies, because the binary truncates `output/` in
place. The restart is a warm restart (`MODELS.md` section 9.6): the primitives and
species are authoritative, the lower ghost rows are re-derived, every solver
control restarts. No crash durability across a server failure, no exact
continuation and no `best` reference are claimed.

**The import.** 136 state directories (90 cases and 46 flux-closure iterates), 222
generations, **134 indexes, 86 with a `latest_certified`**, 48 with none, +132 MB
(the `models/` tree 1992 MB to 2124 MB); a second run publishes nothing. The
certificate attachment is a surrogate (the header's second plus agreement of
verdicts) and every manifest says so: 135 of 222 attached, 87 not, 84 of those
being `output_pre_L34/` states whose own log was overwritten. 81 of the 134
published states were written by the evaluate pass over a solve that is not on
disk.

**The two defects exposed.** `molecular_scalar_gj1132_kzz1e9/HeH0.083` published
the seed of L33's forty passes (md5 `1229b961323d7ac710146650da72454a`
= `.L22/i4_kz0083/output/Hydro_ioniz_IC.txt`): that seed is now its first
generation and the state the passes left (md5
`1ee5430a3c8560c40f017bc5eb618a8b`) is `latest_complete`, neither certified.
`molecular_scalar_gj1132_wellmixed/HeH0.083` has no `ENDING` file; its ending is
classified from the log (no case of the catalog carries an `ENDING` yet, none
having run since D8a).

**The stale certificates.** Since D5b-2 the three molecular states
`molecular_scalar_gj1132_kzz1e9/HeH2.13`, `molecular_scalar_gj1132_kzz1e9/HeH9.7`
and `.L22/i3_alt` are refused on the cell-1 continuity row (19.97, 110.0 and 27.13
floors) with their physical cells unchanged. The first two carry
`stale_under: boundary_model_v1` in their manifests and are no longer offered as
certified seeds; their `REPRODUCE.md` still states the certification of the solve.
`.L22/i3_alt` is not an indexed directory (it holds no `states/`, READ from a
listing), so its stale certificate is recorded in the D5b-2 memo only. Under D2b
the two indexed states read the same refused numbers (1.558e-08 and 8.240e-08
against 7.8e-09 and 7.5e-09); the atomic `kzz1e9/HeH2.13` stays CERTIFIED.

A further oddity: `atomic_scalar_gj1132_kzz1e9/HeH2.13`'s `output/Hydro_ioniz_IC.txt`
is byte for byte the seed of the earlier run, not the state it wrote, so which
state the L34 evaluate pass was handed is worth checking in L34's record.

---

## 4. Open items and user decisions

### 4.1 What the user decides

- **The direction rule of D2b.** The decision text said the direction is decided
  from the sign of the matched face velocity and the window must not override a
  local reverse flow. D2b reads it from the wind window's mass flux where the
  window has standing (`|M_wind| >= 1e-12`, `d_window <= 2e-3`) and from `v_b`
  elsewhere, with a rounding floor near 1e-15 of the sound speed at which the
  reservoir owns the contact. The measured reason (READ, the plan's D2a row and
  the D2b memo section 2): on the three certified catalog states `v_b` is -7.3,
  -2.1 and -3.4 cm/s while the base face carries +0.9995 to +0.99998 of the
  window flux; `v_b` there is set by the 1e-4 pressure residual a converged state
  leaves at the base, a hundred times the flow's own velocity, so a bare sign
  would put every certified state on the reverse branch. No state has been
  measured where the window has standing and the local face flux reverses.
- **D1b**: move the default to `2.0d-4` for new configurations, only after every
  historical omitted-key input is pinned from the D1a record.
- **The retention policy** for the 84 `output_pre_L34/` generations that carry no
  verdict. Nothing deletes them.
- **D9 step 3**: re-solve the reference states under boundary model v2, including
  the three molecular states whose certificates are stale since D5b-2, then the
  one golden refresh, with its movement against the 2026-09-17 set reported.

### 4.2 Technical items noticed and not fixed

Each READ from the memo or report named.

- `eval_dt` takes the CFL minimum over the ghost cells too, so a boundary refresh
  at a monitor step moves `dt` and the marching trajectory (D5b-3).
- `cx_add_to_turnover` reads `cx_kc` without a cell guard; it bounds a
  normalization, not a source (D7a, D7b).
- `element_base_diffusive_flux` is recorded and not printed; face 0 in
  `element_flux_profile.txt` would change a gated product of `lower_profile`
  (D2b).
- The ghost's `p/n_part` and the level's temperature differ by 4.7 per cent on a
  molecular base; conduction reads the level's, advection the ghost's (D2b, D5a).
- The ghost closure has no temperature leg: the molecular partition is closed
  against the ionization at a fixed temperature, and the ghost's ionization is
  defined only to the sweep's acceptance band, which the base continuity row
  amplifies to tens of floors (D5b-2).
- Fourteen invocations in `src/tests/grid_and_gates/run.sh` ignore
  `EXHALE_TEST_OUT` and run in the tree's `build/tests/grid_and_gates/` (D5b-3).
- `add_entry` returns silently when the certification report is full
  (`cert_max_entries = 32`) and every caller then writes into the previous entry;
  not reachable on any present configuration (D3b report).
- The key-on `wasp_full` stationary solve does not certify with or without
  metals, so no certified metal-bearing transported reference exists (D7b).
- `kzz1e9/HeH0.083` at the held bound: a JFNK line-search collapse (`lam` near
  1e-6) costs the whole iteration cap before the outer loop can judge the pass
  (the bounded experiment).
- `steady_newton.f90` near line 14313 builds a diagnostic's carrier densities
  with the mass-row density for every carrier, wrong for a stage row;
  `element_inventory_report` skips helium; the hydrogen and oxygen closures of
  `carrier_source` still clip a negative remainder (D6).
- The stage operator's boundary rows are first order (0.93 against 1.92 in the
  interior); `row_terms_phys` carries no advective term for a stage row (D3a).
- `write_output.f90` writes the state files list-directed with no declared
  round-trip format (D5a); `parse_dump.txt` still prints `ntot_bc` and `dp_bc`
  unlabeled because it is compared against `parse_golden` (D5b-3).
- An entry's `reason` text prints only when it refuses (D3b).
- `docs/audit_20260905/plan_20260918_rev1_source_probe.f90` no longer compiles
  against the current `cx_add_to_fvec`; it is the dated record of that review
  (D7a, D7c, D7d).
- `src/tests/grid_and_gates/ionization_transport_atomic.sh` drops `metals.inp`
  from case A since D7a; the metal-bearing atomic branch is supported again and
  the fixture could carry it (D7b report).
- `coupled_block_jacobian`'s two standing rows (the L22 coupled block) are
  unchanged by every item.

---

## 5. Paths

**The plan and its reviews**: `docs/PLAN_20260918_rev2.md` (the plan of record,
with the DONE and RESULT paragraphs), `docs/PLAN_20260918_rev1.md`,
`docs/PLAN_20260918.md`, `docs/PLAN_20260918_review.md`,
`docs/PLAN_20260918_review1.md`, their retained checks under
`docs/audit_20260905/`; `docs/DECISION_D2a_D8_codex.md` and
`docs/DECISION_D2a_D8_review.md`.

**Memos of this series** (all under `docs/`):
`lhs1140b_stationary_D1a_20260918.md`, `_D3a_`, `_D3b_`, `_D4a_`, `_D4b_`,
`_D5a_`, `_D5b1_`, `_D5b2_` (which also records D5b-3 in its section 11),
`_D6_`, `_D7b_`, `_D7c_`, `_D7d_`, `_D8a_`, `_D8b_`, `_D8bound_`, `_D9_`, all
dated 20260918, and `lhs1140b_stationary_D2b_20260918.md` (landed 2026-09-19).
D7a wrote no memo under `docs/`.

**Item reports** (the session scratch directory):
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/PLAN_20260918_<item>_report.md`
for D1a, D2b, D3a, D3b, D4a, D5a, D5b1, D5b2 (with D5b-3 appended), D6, D7a, D7b,
D7c, D7d, D8a, D8b, D8bound and D9; the run directories and builds of each item
under the same directory, `<item>/`.

**The update log**: `docs/Update_EXHALE_stage2.md` section 13 and its typeset
twin `docs/Update_EXHALE_stage2.pdf` (273 pages).

**The catalog contract**: `LHS1140b/MODELS.md` section 9;
`LHS1140b/models/publish_state.py`, `import_legacy_states.py`,
`legacy_state_map.txt`; the fixtures under `LHS1140b/models/tests/`.

**The test suites of this series**: `src/tests/boundary_state/`,
`src/tests/charge_exchange_rows/`, `src/tests/carrier_helium_inventory/`,
`src/tests/stage_row_balance/`, and the new rows of `grid_and_gates`,
`certification`, `carrier_retry` and `ionization_stage_flux`.
