# Session handoff, 2026-09-12

The state to restart from. Read `docs/code_status_20260910.md` first, then
`docs/ISSUES_20260909.md` (section 2 now carries the rows of the stage-2
review's findings), `docs/PLAN_20260912_review_fixes.md` (the plan of record
of this day, with its status section at the close) and
`docs/Update_EXHALE_stage2.md` **section 9** (the dated record of items Q1 to
Q4). The previous handoff, `session_handoff_20260911.md`, is the entry state
of this session and is not updated.

## Verified state

- Tree: branch `v1.00`, HEAD `ffbefad` ("P21, docs", 2026-09-12 early
  morning), the entry text every measurement of this day is controlled
  against. On top of it, UNCOMMITTED: the paths listed at the close.
- The user's instruction of the day: analyze
  `docs/Update_EXHALE_stage2_review_20260912.md` against the code and, where
  right, plan and start the corrections. All seven findings R1 to R7 were
  confirmed in the source and corrected (Q1 to Q4); the confirming lines are
  at the head of the plan, the corrections and their measurements in
  section 9 of the update log.
- The state contract at a fixed hydrodynamic state (decided for R2/R3, the
  one rule every composition update at fixed hydro follows now): the
  conserved variables `u` are held; pressure and temperature are recomputed
  from the unchanged thermal energy through the caloric EOS of the new
  composition. `relax_photochemical_composition` takes `u` and returns `p`
  and `T`; the archived experiment drivers under
  `docs/audit_20260905/partition_experiment_20260911/` and the review's own
  probes under `docs/audit_20260905/stage2_review_20260912/` call the older
  signatures and are records, not built by anything.
- Goldens: unchanged. `mol_carrier` (12000 steps, one thread),
  `wasp_full_newton` (reload), `mol_diffusion` and `lower_profile` (300 steps)
  measured identical to their controls in the data rows; the regression gate
  of the merged tree: `make check` REGRESSION PASS, all cases byte-identical
  (2026-09-12 10:35 KST, section 9 "Gates").
- Suites on the merged tree (private build, 2026-09-12): `carrier_retry`
  138/0, `certification` 84/0, `element_operator` 28/0, `attempted_step`
  70/0, `steady_species_rows` 195/0, `krylov_and_dogleg` 334/0,
  `diffusion_tests` 35/0, `carrier_returned_state_acceptance` 36/0,
  `carrier_reference_scales` 14/0, `carrier_constraint_attribution` 10/0,
  `species_masses` 9/0, `steady_completion_flag` 3/0; `grid_and_gates` 198/1
  with the one FAIL the stagnation row on its OLD configuration (the suite
  was started before that script was re-anchored), and the re-anchored
  `output_state_consistency.sh` 4/4 on its own afterwards.

## What is open

- The joint progress measure changes when the outer loop halves omega and
  the carrier movement bound: on the 12-pass hot-Uranus production loop no
  pass triggered it (H2 row 7.38e-2 to 5.49e-2, as the 110-pass ladder); on
  a starved hydro solve (`EXHALE_JFNK_MAXIT=5`) the loop now keeps going
  while the hydro rows fall, where it used to stop on the flat species row.
  Whether that is the better behavior on a real production run has not been
  measured beyond these two configurations.
- `newton_drop` of the helium mass-fraction Newton moved from 1e-6 to 2e-5
  (Q2): on the HD 209458 b element reload the judged rows at pass 2 are
  unchanged and omega differs (0.5 against 0.25). A wider look at whether
  other columns sit in the 1e-6 to 2e-5 band has not been taken.
- The margin of 10 behind the mass-row tolerance was measured on three
  fixtures only, none a near-zero-flow or molecular-EOS column (addendum of
  `docs/certification_tolerance_anchoring_20260910.md`).
- Everything section 3 of `docs/ISSUES_20260909.md` lists as open stands:
  no stationary solve with a species row certifies yet; the hot-Uranus H2
  front walks outward under a bounded pass and leaves the domain under an
  unbounded one (P21, re-measured on the new contract: 61 steps, front at
  3.079 R_p).

## Paths, and where the numbers of this session live

Scratch (session-local, `/tmp/claude-1000/.../scratchpad/`): `Q1/` (the
production runs `m1_unbounded`, `m2_12pass`, `atomic3_new` and `_ctrl`,
`wasp_new` and `_ctrl`, `mol_carrier`, `q4_phys`, `q4_init`, `q4_phys_ctrl`,
the stagnation configurations `stag_*`, `Q1_measurements.md`), `Q1probe/`
(the review's 500-cell probe adapted and its log), `Q1red/` (the RED tree),
`Q2/Q2_report.md`, `Q3/Q3_report.md`, the suite logs `q1_suite_*.log` and
`q1m_suite_*.log`. Nothing of it is in the repository beyond what section 9
quotes.

## Uncommitted files (`git status --short` at the close)

```
 M docs/ISSUES_20260909.md
 M docs/Update_EXHALE_stage2.md
 M docs/certification_tolerance_anchoring_20260910.md
 M docs/code_status_20260910.md
M  examples/README.md
 M src/EXHALE_main.f90
 M src/modules/functions/binary_element_diffusion.f90
 M src/modules/lower_atmosphere/diffusive_photochemistry.f90
 M src/modules/time_step/attempted_step.f90
 M src/modules/time_step/certification.f90
 M src/tests/carrier_retry/carrier_retry.f90
 M src/tests/carrier_retry/run.sh
 M src/tests/certification/certification_contexts.f90
 M src/tests/element_operator/element_operator_tests.f90
 M src/tests/grid_and_gates/README.md
 M src/tests/grid_and_gates/output_state_consistency.sh
```

## User-gated (unchanged)

Physical reruns of the planet folders and `benchmarks/`; `Load IC` case
regeneration on the current grid; a converged run reaching 3000 K for the CO
domain; paper and poster re-convergence; the LHS 1140 b write-up; the GitHub
public switch. None is picked up autonomously.
