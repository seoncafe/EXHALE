# Session handoff, 2026-09-11

The state to restart from. Read `docs/code_status_20260910.md` first (its
section 3 now carries the 2026-09-11 reopening paragraph), then
`docs/ISSUES_20260909.md` (sections 3.1, 3.7 and 5),
`docs/PLAN_20260911_partitioned_solver.md` (the plan of record, with its
status section at the close) and `docs/Update_EXHALE.md` **section 8** (the
dated record of items P1 to P9, of P10 to P12 which the user's further
instruction of the same day added, and of P14 to P18, the afternoon's Codex
review findings and the user's two decisions, with the two subsections "The
Codex adversarial review (2026-09-11)" and "User decisions of 2026-09-11
(afternoon)"). The previous handoffs
(`session_handoff_20260910.md`, `_20260908.md`, `_20260907.md`,
`_20260905.md`) are records and are not updated; the 2026-09-10 one is the
entry state of this session.

## Verified state

- Tree: branch `v1.00`, commit `a0f4292` (2026-09-10 23:30 KST,
  "Well-balanced flux option and stationary-solver diagnostics"), which is
  the entry text every worker of this session measured its control against.
  On top of it, UNCOMMITTED: the 40 modified and 9 untracked paths listed
  below, of which `src/tests/test_columns.f90` is a new file. Item P18 was
  still running when this list was taken, so `steady_newton.f90` and
  `steady_species_rows_tests.f90` are mid-edit in it.
- Goldens: unchanged since the 2026-09-10 refreshes (Stage D gate at 00:40,
  `lower_profile` again at 11:10 for N29; the previous sets are
  `backup/regression/golden_pre_stageD_20260910/` and
  `golden_pre_n29_20260910/`). Nothing this session refreshed a golden.
  Every marching impact measurement of items P1, P2, P3 and P4 came out
  byte-identical apart from the provenance `run=` timestamp of
  `Hydro_ioniz.txt`, each against a control built from the entry text of
  the one file that item owned.
- The regression gate: `make check` on the shared `build/` rebuilt from
  this tree (advisor, 2026-09-11 09:40 KST): **REGRESSION PASS**, every
  DEFAULT_CASES output and `wasp_full_newton` byte-identical to the
  goldens. No golden moved, because no default case enters the stationary
  solve or the fixed-wind relaxations. Second `make check` after P10 to P12
  (advisor, 2026-09-11 15:15 KST, the shared build rebuilt from the merged
  tree): **REGRESSION PASS**, every DEFAULT_CASES output and the cold-start
  `wasp_full_newton` byte-identical to the goldens; no golden moved.
  **Third `make check` after P14 to P18** (advisor, 2026-09-11 21:10 KST,
  the shared build rebuilt from the final tree): **REGRESSION PASS**, every
  DEFAULT_CASES output and the cold-start `wasp_full_newton` byte-identical
  to the goldens; no golden moved.
- Suites, advisor's snapshot build at the close of P1 to P9 (READ from the
  closing brief): `element_operator` 24/0, `carrier_retry` 122/0,
  `krylov_and_dogleg` 290/0, `certification` 57/0, `steady_species_rows`
  194/0 and 195/0 after P8, `grid_and_gates` 158 PASS with 3 FAIL that are
  missing files of the snapshot copy (`inputdata`, `examples`) and not
  assertions. P4 measured the whole `grid_and_gates` suite at 178/0 in the
  live tree. `residual_determinism` still has the row that FAILs by design.
- Suites after P10 to P12, each MEASURED by the item that owns the file
  (READ from the three reports): `krylov_and_dogleg` **309/0** (P10, 19 new
  rows), `grid_and_gates` **198/0** with `hydrostatic_residual` **60/0**
  (P11, 20 new rows), `element_operator` **24/0** and `attempted_step`
  **70/0**, `carrier_returned_state_acceptance` **36/0**, `carrier_retry`
  **122/0**, `carrier_reference_scales` **14/0**,
  `carrier_constraint_attribution` **10/0**, `constrained_network_layout`
  **8/0**, `species_masses` **9/0** (P12, in a tree with no `build/` and a
  work directory outside it as well as in the repository),
  `certification` 57/0 and `adv_static_limit` 53/0 unchanged.
  `steady_species_rows` stood at 193/2 on P10's and P11's objects, the two
  rows being the superseded semantics the advisor then replaced; on the
  merged tree (advisor's snapshot build after P10 to P12, 2026-09-11) it
  is **195/0**, and the thirteen suites that link the changed modules
  all PASS on that build (`element_operator` 24/0, `carrier_retry` 122/0,
  `krylov_and_dogleg` 309/0, `certification` 57/0, `grid_and_gates`
  187/1 with the one FAIL a WASP-52b SED file absent from the snapshot
  copy, `attempted_step` 70/0, `carrier_returned_state_acceptance` 36/0,
  `species_masses` 9/0, `carrier_reference_scales` 14/0,
  `carrier_constraint_attribution` 10/0, `constrained_network_layout`
  8/0, `steady_completion_flag` 3/0; MEASURED).
- The bitwise guard: the `wasp_full_newton` reload ends `info = 0`, `||R||`
  1.667e-9, CERTIFIED, `log10 Mdot` 13.30, unchanged, on the merged tree
  (advisor) and on each worker's own pair.
- Fixtures: `backup/regression/atomic_elem_newton` and
  `backup/regression/carrier_elem_newton`, both with a new dated section
  "Partitioned route, 2026-09-11" in their README carrying the recipe and
  the measured ladder beside the same day's coupled control, and a second
  dated section "After P16 and P17 (2026-09-11)" carrying the afternoon's
  outcomes (the atomic fixture CERTIFIED at pass 12, the carrier fixture's
  mass verdict at its own anchor and its refusing entries 2 to 1, and the
  110-pass long continuation in the carrier README);
  `backup/regression/wasp_full_newton/IC` unchanged.

## What this session settled

The user's instruction of 2026-09-11 reopened the stage-2 solver program as
the partitioned route and asked for the errors of
`docs/solver_approach_analysis_20260910.md` to be corrected, the experiments
to be run, and the route to be adopted if it measured better. All three were
done.

1. **Three defects of the partitioned route were real and are repaired**, each
   with its reproduction: a failed element composition solve was handed back
   as an update (P1; 17.7 percent mass-closure error and a nonfinite
   chemistry become the entry composition restored to the last digit and an
   explicit failure), the carrier movement bound was tested after the step
   was applied (P2; drift 2.81e-2 at `trust = 0.01` becomes 9.996e-3, and
   2.98e-3 at `trust = 0.003`), and the stationary solve's stop test was
   weaker than its own acceptance test (P3; a loop-top `info = 0` refused at
   return on a mass row of 5.954e-12 becomes an iteration that continues and
   an honest `info = 1` on its budget).
2. **The momentum row's reference scale under `Well balanced: True` is
   fixed** (P4). It was exactly 1 at every cell by construction; it is now
   read against the pressure force of the cell's own hydrostatic
   equilibrium, and with the key on `wasp_full_newton` **certifies**
   (`info = 0`, `||R||` 9.719e-10, momentum row 9.719e-10 against 1e-8)
   where it was refused on that one entry at exactly 1.000. The key is still
   default off and the default path is byte-identical.
3. **The outer iteration is a stated contract** (P6): a state is accepted
   only by the certification of the full set of active equations on the
   refreshed state, a hydrodynamic `info = 0` with an alternated species row
   refusing is not an accepted state, a nonzero hydrodynamic flag does not
   end the loop, the endings are named, the last pass takes no composition
   update, and a pass that fails to move the joint measure shortens both
   updates.
4. **The pseudo-time start of `Restart intent: stationary` was wrong** and is
   now 1.0, the marching hand-off's value (`EXHALE_PTC_DTAU0` overrides; the
   direct steady route keeps its CFL start). At the CFL start the carrier
   reload's hydrodynamic solve was handed back at `||R||` 2.685e-1 after 40
   iterations and the atomic one stagnated at 1.490; at 1.0 they hand back
   1.498e-9 and 1.167e-8. No regression case sets the key, so no golden
   moves.
5. **The measurement (P7) adopts the partitioned route.** Neither fixture
   certified at the time of that measurement, so the plan's first adoption
   branch was not met; the second is met on both, a strictly lower joint
   residual in no more wall time. (The atomic fixture CERTIFIES after the
   afternoon's item 13 below; the carrier one still does not.) The
   table and every number are in `docs/Update_EXHALE.md` section 8. Short
   form: atomic, 12 passes in 1043 s, every elemental wind row within 1e-5
   and only the hydrodynamic mass row refusing, against the coupled
   control's three passes in 2392 s with the worst gated row at 8.29e-5;
   carrier, 8 passes in 173 s strictly below the coupled control on every
   judged row in 321 s, and 40 passes in 831 s bringing the H2 wind row from
   7.38e-2 to 2.23e-2. `Coupled carrier solve` already defaults to False,
   so nothing in the default configuration changes; `True` remains the
   coupled measurement arm.
6. **The documents that carried the withdrawn inference are corrected**
   (P5, P8): the spatial truncation error is not a lower bound on the
   algebraic residual (the frozen-background H2 wind row measures
   2.449021e-13 on the same 500-cell grid), the Ritz ratios are not
   condition numbers, and decision 22 (a) stands with a different
   justification in all three documents.

Then the user asked for the open items the plan had left to the user to be
fixed as well, the pressure term among them. Three more items ran (P10, P11,
P12 with P12b) and the advisor made two further changes.

7. **The best-iterate ledger reads the hydrodynamic rows against their own
   certification tolerances** (P10). It divided the slot by the run's
   `Resid tol`, so the ledger could not see a mass row above 3e-12 while
   `||R||` was below 1e-8; P10 made the slot
   `max(mass/3e-12, momentum/1e-8, energy/1e-6)` (P17 then made it cell by
   cell and P18 left one expression of it, `..._by_cell`), and
   `d < 1` is exactly the condition the loop-top stop and the acceptance at
   return impose. MEASURED: `wasp_full_newton` reloaded is unmoved, all
   eleven output files byte-identical, `info = 0`, `||R||` 1.667e-9,
   `log10 Mdot` 13.30; the carrier reload hands back a lower mass row at a
   higher `||R||` from pass 2 (4.85e-12 at 2.459e-9 against the control's
   6.98e-12 at 1.593e-9), which is the trade the certification asks for and
   the superseded ledger could not see; the HD 209458 b element reload ends
   its twelve passes on the same single refusing entry, mass 8.893e-10 at
   cell 3 against the control's 1.054e-09 at cell 16, in 1064 s. Tests
   `krylov_and_dogleg` 290/0 to 309/0, with 300/9 on a RED build.
8. **The momentum row's reference scale is the largest physical term of the
   momentum equation** (P11), the user's "pressure term". It was built from
   the pieces of the discretization; it is now `max` of the ram divergence,
   the pressure gradient and the weight (and `S_visc`), gathered by one
   routine from the stored face data on both reconstructions and under both
   values of the well-balanced key, with the residual itself untouched.
   MEASURED: the three terms add up to the row they scale (0 under WENO3,
   1.7e-16 default and 5.3e-15 key-on under PLM); the zero-gravity PLM
   artifact, a scale of exactly 1.000 of `2 p/r` where every physical term
   is zero, now reads the `tiny` floor; on supersonic uniform flow the scale
   agrees with the ram divergence to 1.9e-14 where it stood 6.7e-2 above it;
   on the discrete equilibrium it was 28 percent below the weight it
   balances and is now never below it. `wasp_full_newton` stays CERTIFIED
   with the key off (momentum row 8.469e-10 to 1.165e-9, states apart by at
   most 2.6e-9 relative) and with it on (9.719e-10 to 9.621e-11, outputs
   byte-identical); `mol_carrier` marching is byte-identical. Tests
   `hydrostatic_residual` 50/10 to 60/0, `grid_and_gates` 188/10 to 198/0.
9. **Two test-hygiene defects** (P12, P12b). The `build_stamp` fallback was
   written where the compile cannot find it, which fires for any work
   directory outside `$root/build`, that is for every concurrent worker:
   seven `run.sh` now take `${EXHALE_OBJDIR:-$root/build}/build_stamp.f90`
   and compile whichever stamp is in force as the first file of the closure,
   and four gained the `EXHALE_TEST_OBJDIR` override without which the
   no-build arm cannot be run at all. MEASURED RED (0/1, exit 1, `Cannot
   open module file 'build_stamp.mod'`) in six suites and GREEN at their
   full counts after. The synthetic mass-closed column, written twice, is
   now the one module `src/tests/test_columns.f90`, with `element_operator`
   24/0 and `steady_species_rows` 195/0 before and after and the two full
   logs byte-identical.
10. **Two advisor changes beside them**: the two rows of
   `steady_species_rows_tests.f90` that stated the superseded ledger
   semantics are replaced by
   `hydrodynamic_distance_reads_each_row_against_its_own_tolerance`
   (`1e-9/3e-12` from a mass row of 1e-9 beside a momentum row of 1e-13 and
   an energy row of 1e-8), and the `Coupled carrier solve` key comment in
   `input_read.f90` no longer says the alternation cannot converge: it names
   the default, what the alternation is judged by, and the two reloads of
   section 8 for what it reaches.

Then a Codex adversarial review was run on the working tree and the user took
two decisions on the plan's open items. Four more items ran (P14, P15, P16,
P17) and P18 followed; three further changes were made by the advisor.

11. **Every ending of the stationary outer iteration hands back one consistent
   state** (P14, the Codex review's first finding, verified in the source
   before the item ran). On a stagnation ending the composition had been
   updated while the particle densities and the temperature written beside it
   had not. MEASURED: re-evaluating the written state moved its temperature by
   **1.317089e-8** relative, where the same fixture ending on its pass budget,
   which takes no update, moves by 5.921e-16, and after the fix by 6.499e-14.
   The progress control now stands between the certification of the pass and
   the update, so the pass that ends the iteration takes no update; the ending
   is announced with the other endings; and the element path refreshes the
   particle densities and the temperature between `relax_element_composition`
   and its equilibrium sweep, which before equilibrated the new composition at
   the old temperature. `grid_and_gates` 199/1 on the entry text to **200/0**,
   the new arm asserting the round trip against 1e-12. `wasp_full_newton`
   reloaded unmoved (CERTIFIED, `||R||` 1.495e-9, `log10 Mdot` 13.30); the
   carrier reload identical pass line by pass line; the element reload parts
   from pass 2 and both arms end the same way.
12. **The kind-generic residual arm hands back the face departures its
   momentum row was built from** (P15, the Codex review's second finding, and
   the third entry of the plan's "what remains"). `store_the_interface_fluxes`
   wrote `face_flux`/`face_p` alone, so under `Well balanced: True` the
   momentum row's pressure-gradient term, the row scale and every measure read
   against it belonged to the last `RK_rhs` evaluation; they are now the
   assembled state's. The reconstruction continuation, which CAN reach the arm,
   is refused rather than weighted, no single pair of departures satisfying the
   blended row identity. MEASURED RED 2.098877e+1 (PLM) and 2.032204e+1
   (WENO3) on the pressure gradient, 3.262392 and 3.428893 on the row scale,
   1.190363e+1 and 1.181676e+1 on the row rebuilt from the stored face data;
   GREEN 0 bitwise for the generic-double instantiation and 2.6e-16 to 2.7e-15
   for the quadruple one. `krylov_and_dogleg` 309/0 to **321/0**. The one
   quantity that moves in a run is the measure inside the loop, iteration-1
   `||R||` 1.756e-3 against 1.753e-3, overstated by 0.17 percent; the
   certification at return was never on the arm (the sweep kind is reset to
   marching before it), so `wasp_full_newton` stays CERTIFIED at the same
   momentum row 9.621e-11 with the key on, and `mol_carrier` marching is
   byte-identical.
13. **The mass row's tolerance is anchored on its measured rounding floor**
   (P16, the user's first decision). `tol(j) = max(3e-12, min(1, 10 *
   floor(j)))` with the floor `eps rho(|v| + c_s) A / dV` over the row's own
   scale, the verdict taken cell by cell. The brief's first floor expression
   was inert, MEASURED at 4.4409e-16 on all 500 cells of all four states
   (between one and two ulps, the row's scale being the larger addend and the
   `dV` cancelling); the adopted estimate bounds the measured one-ulp step
   within 1.5931, 1.8857, 2.6806 and 1.7733 on states whose base Mach numbers
   span three decades, and `c_round = 10` is 3.7 times above the largest of
   those. **The HD 209458 b element reload CERTIFIES, at pass 12**, where the
   control ends the same twelve passes NOT CERTIFIED on the mass row
   (7.159e-10 against 3e-12 at cell 13); `wasp_full_newton` is unmoved, its
   floor below 3e-13 at every cell. `certification` 57/0 to **70/0**.
14. **The best-iterate ledger reads the mass row against the cell's own
   tolerance** (P17, P16's own outside-scope finding). The ledger's
   hydrodynamic slot is now `max over cells and rows of |R|/scale/tol_k(j)`
   with the mass row's tolerance taken through P16's `mass_row_cell_verdict`,
   which is the functional the acceptance applies. MEASURED on the HD 209458 b
   reload: the best judged distance over the twelve passes runs **2.2e+2 to
   6.3e+5 in the control and 1.321e-1 to 1.868e-1 after**, while the acceptance
   certifies the pass-12 state in both arms at a mass-row distance of 0.13 to
   0.14. On the carrier reload the hydrodynamic solve returns `info = 0` at all
   three passes against the control's 0, 2, 2, and the refusing entries of the
   last pass fall **2 to 1**, the H2 carrier row alone. `krylov_and_dogleg`
   321/0 to **334/0**, with 327/7 on a semantic RED build. Cost not measurable
   (1.206 s against 1.216 s an outer iteration).
15. **One entry point for the hydrodynamic distance** (P18, DONE): the
   fixed-tolerance form of the distance, which had no production reader
   left after P17, is deleted; every test row that read it restates its fact
   through `hydrodynamic_distance_from_certification_by_cell` on a one-cell
   column whose floor puts it under the fixed tolerance. MEASURED (P18):
   `krylov_and_dogleg` 334/0, `steady_species_rows` 195/0, `certification`
   70/0, no row lost (two renamed); `wasp_full_newton` reload byte-identical,
   run logs differing in the execution time only; `grep` finds no reader of
   the deleted name in `src/`.
16. **Three advisor changes beside them**:
   `assert_written_state_is_the_accepted_one` carries an absolute floor of
   1e-12 beside its 1e-10 relative test (P14 measured the false alarm,
   3.0982e-14 as written against 2.7487e-14 accepted on a state flat to
   fourteen digits); the three `grid_and_gates` shell rows that wrote to a
   fixed directory now honor `EXHALE_TEST_OUT`
   (`momentum_row_from_fluxes_only.sh`, `base_level_single_statement.sh`,
   `sed_coverage_stop.sh`), so two concurrent runs of the suite no longer
   collide; and the paragraph at `assemble_residual` in `steady_residual.f90`
   that stated what the arm does not hand back is rewritten to what P15 made
   true.

## What is open, with its evidence and its next item

1. **The H2 front of the hot-Uranus carrier reload, and the wind floor the
   long continuation exposed.** The user's second decision of 2026-09-11 ran
   the continuation to its stop: the pass budget was set to 400 and the
   progress control ended the run at pass **110**, the carrier trust having
   been halved 1.0e-2 to 5.0e-3 to 2.5e-3 at passes 109 and 110. MEASURED (the
   advisor's snapshot build of the merged tree after P1 to P12 and before P14
   to P17, 8 threads; the table is in `docs/Update_EXHALE.md` section 8, "User
   decisions of 2026-09-11 (afternoon)", and in
   `backup/regression/carrier_elem_newton/README.md`): the H2 front does NOT
   stop, `x2 = 0.5` moving **1.0726 to 1.1230 R_p** over 108 passes at about
   5e-4 R_p in a pass with no sign of slowing; and from pass 40 the worst gated
   row is no longer at the front but stands at **2.2e-2, nearly uniform in r**,
   its cell walking outward through the wind 443 to 500 and reaching the outer
   boundary at pass 108 while the value does not fall. The reading: the bounded
   carrier pass, trust 1e-2 on the H2 mixing-ratio maximum and met by one
   transport step, advances the front but leaves the wind cells, whose H2
   mixing ratio is orders of magnitude below the maximum, at one implicit step
   in a pass, so their own H2 balance never relaxes, and the drift measure,
   absolute in the grid maximum, does not see them. Two questions stand, and
   they are the whole of what is open here. (a) **Does the front stop inside
   2 R_p?** That is physical, not a solver question, and its route is the P23
   comparison against Koskinen et al. (2022) Model A. (b) **A relaxation step
   for the wind cells' H2 balance, distinct from the bounded front advance.**
   A bounded front advance and a wind-cell relaxation are two different
   operations; the second needs its own step, a cell-relative drift measure or
   a wind-window relaxation after the bounded pass, before the H2 wind row can
   be judged at 1e-5 at all. Recorded and NOT adopted, and no further
   continuation of the present operator answers either question.
2. Whether the two fixture certifications should be re-measured with
   `Well balanced: True`, which is now meaningful on both paths (P4 made the
   key-on momentum scale a term of the equation and P11 the same for the
   default path), is open as a measurement and not as a defect.
3. Everything the 2026-09-10 handoff listed as open items 2 to 5 stands
   unchanged: the layer's element flux conservation, the element operator's
   discretization order (whose acceptance P8 restated as a measured order
   and the accuracy it buys, not a comparison against a gate), the hygiene
   list, and the `_adv` product on states that are not stationary.

## Paths, and where the numbers of this session live

- **In the tree, permanent**: `docs/Update_EXHALE.md` section 8 (every number
  of items P1 to P8 and of P10 to P12, and the P7 comparison table);
  `docs/PLAN_20260911_partitioned_solver.md` (the plan and its status
  section); `backup/regression/atomic_elem_newton/README.md` and
  `backup/regression/carrier_elem_newton/README.md` (the recipe and the
  ladder of the partitioned route beside the same day's coupled control);
  `docs/solver_approach_analysis_20260910.md` (corrected, with section 9
  indexing the corrections); `docs/ISSUES_20260909.md` sections 2, 3.1, 3.7
  and 5; `docs/code_status_20260910.md` section 3;
  `docs/certification_tolerance_anchoring_20260910.md`; `TO_BE_DONE.md`.
- **In the tree, the reproductions**:
  `docs/audit_20260905/partition_experiment_20260911/` (the standalone
  driver `transport_wind_experiment.f90` and its build and run scripts,
  which items P1, P2 and P3 each linked against their own object directory
  by a private copy of those scripts) and
  `docs/solver_partition_experiment_20260911.md`, the experiment report the
  three defects were named from.
- **TEMPORARY, and will be gone**: everything under the session scratch
  root `/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/a37fb7bf-6fd2-41fe-a495-7983c09938ff/scratchpad/`,
  including the worker reports `P1/P1_report.md` to `P19/P19_report.md`, the
  110-pass long continuation `P7/carrier_long400/` (its `run.log` is the
  source of every number of decision 2, and its `SUMMARY.md` the advisor's
  reading of it; both are recorded in section 8 and in the carrier fixture
  README), the
  suite logs before, after and RED of each item, and the run directories
  `P6/runs/*` and `P7/*` whose `run.log` files are the source of every P7
  number quoted above. Nothing in the tree depends on them; what had to
  survive was copied into the documents named above. The reload fixtures
  reproduce the runs: the recipes are in the two fixture READMEs.
- The experiment's own private build is retained outside the tree at
  `/tmp/exhale_partition_20260911.ZxXJEB` per its section 9; it can be
  rebuilt from the archived scripts.

## Uncommitted files (`git status --short` at the close)

```
 M docs/ISSUES_20260909.md
 M docs/Update_EXHALE.md
 M docs/Update_EXHALE.pdf
 M docs/Update_EXHALE.tex
 M docs/certification_tolerance_anchoring_20260910.md
 M docs/code_status_20260910.md
 M docs/solver_approach_analysis_20260910.md
 M src/EXHALE_main.f90
 M src/modules/files_IO/input_read.f90
 M src/modules/functions/binary_element_diffusion.f90
 M src/modules/lower_atmosphere/diffusive_photochemistry.f90
 M src/modules/time_step/RK_rhs.f90
 M src/modules/time_step/certification.f90
 M src/modules/time_step/hydrodynamic_rows.f90
 M src/modules/time_step/hydrodynamic_rows_body.inc
 M src/modules/time_step/steady_newton.f90
 M src/modules/time_step/steady_residual.f90
 M src/tests/attempted_step/run.sh
 M src/tests/carrier_constraint_attribution/run.sh
 M src/tests/carrier_reference_scales/run.sh
 M src/tests/carrier_retry/carrier_retry.f90
 M src/tests/carrier_retry/run.sh
 M src/tests/carrier_returned_state_acceptance/run.sh
 M src/tests/certification/certification_contexts.f90
 M src/tests/constrained_network_layout/run.sh
 M src/tests/element_operator/element_operator_tests.f90
 M src/tests/element_operator/run.sh
 M src/tests/grid_and_gates/README.md
 M src/tests/grid_and_gates/base_level_single_statement.sh
 M src/tests/grid_and_gates/hydrostatic_residual.f90
 M src/tests/grid_and_gates/momentum_row_from_fluxes_only.sh
 M src/tests/grid_and_gates/output_state_consistency.sh
 M src/tests/grid_and_gates/sed_coverage_stop.sh
 M src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90
 M src/tests/species_masses/run.sh
 M src/tests/steady_completion_flag/run.sh
 M src/tests/steady_selfconsistent_residual/run.sh
 M src/tests/steady_species_rows/run.sh
 M src/tests/steady_species_rows/steady_species_rows_tests.f90
 M src/utils/update_log_to_tex.py
?? docs/PLAN_20260911_partitioned_solver.md
?? docs/audit_20260905/partition_experiment_20260911/
?? docs/audit_20260905/solver_approach_checks_20260910.py
?? docs/session_handoff_20260911.md
?? docs/solver_approach_analysis_20260910_review.md
?? docs/solver_partition_experiment_20260911.md
?? poster/WASP-52b_poster_v1.nav
?? poster/WASP-52b_poster_v1.snm
?? src/tests/test_columns.f90
```

`docs/session_handoff_20260911.md` is this file; the two `poster/` files are
LaTeX byproducts of an earlier session and belong to nothing here.
`src/tests/test_columns.f90` is the new module of P12.
`docs/Update_EXHALE.tex` and `.pdf` are GENERATED from the Markdown
(`python3 src/utils/update_log_to_tex.py` plus `latexmk`); they were
regenerated once after P9 and are behind the P10 to P12 and P14 to P18 record
added to the Markdown afterwards, so they are to be regenerated at the gate,
and `src/utils/update_log_to_tex.py` itself is modified in this list. Files
`TO_BE_DONE.md` and everything under `backup/` do not appear because
`.gitignore` excludes them; both were edited this session (the fixture
READMEs and the TO_BE_DONE entries).

## User-gated (unchanged)

Physical reruns of the planet folders and `benchmarks/`; `Load IC` case
regeneration on the current grid; a converged run reaching 3000 K for the CO
domain; paper and poster re-convergence; the LHS 1140 b write-up; the GitHub
public switch. None is picked up autonomously.
