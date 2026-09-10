# Session handoff, 2026-09-08

## State on 2026-09-10 (read this first; the sections below are the 2026-09-08 record)

- Plan of record `docs/PLAN_20260909_rev1.md`; items N-1 to N30 DONE and
  logged in `docs/Update_EXHALE.md` section 7; decisions 13 to 22 taken
  (`docs/To_be_determined_by_user_20260906.md`). The state of the code is
  `docs/code_status_20260910.md`; the problems are `docs/ISSUES_20260909.md`
  (sections 3.1 to 3.8 and 5 current as of N30).
- Goldens: refreshed 2026-09-10 00:40 (`golden_pre_stageD_20260910/` is the
  previous set), `lower_profile` again at 11:10 (`golden_pre_n29_20260910/`);
  `make check` PASS on the shared build at that gate.
- Fixtures: `backup/regression/atomic_elem_newton/IC` (atomic element arm,
  README with the control of every arm), `backup/regression/wasp_full_newton/IC`
  (the certified atomic state, bitwise guard); the carrier reload is made by
  the scratch recipe of the 2026-09-09 session (B5k `mkrun.sh carrier`) and
  is NOT pinned in the tree yet.
- Documentation brought to the stage-2 state 2026-09-10 (DOC1 to DOC3):
  `EXHALE_BC_and_IC`, `composition_restart_and_base_handoff`,
  `transmission_spectrum`, `steady_solver_memo`, `molecular_hydrogen_treatment`,
  `lower_atmosphere_coupling`, `EXHALE_user_manual` (tex and pdf), `README.md`,
  `README_HOWTO.md`, `TO_BE_DONE.md`, `steady_solver_design.md` section 22.
- Next: N31, the column preconditioner (brief: why the finite-difference
  Jacobian action is not additive, then a column-scale or central-difference
  arm, accepted by named outcomes; `<scratchpad>/planrev1_brief_N31.md`).


Continues `docs/session_handoff_20260907.md` (the PLAN rev2 series gate).
This file covers the two batches that followed it: B5, B5b, B5c,
FIELD-SELF, DIFT-LINK, THREAD-DET, HYG-PY, DOCS-LINES (gate 1), then B5d,
B5e, COST4, COST5, COST6, HYG-T9, HYG-4 (gate 2), with COST7 (decision 8)
in progress. The record is
`docs/Update_EXHALE.md` section 7; `docs/Update_EXHALE.{tex,pdf}` is its
generated twin (`python3 src/utils/update_log_to_tex.py`, then latexmk).

## Verified state

Gate 4 (2026-09-08, 15:39 to 18:45), tree after COST7, B5e, B5g, B5h,
ADV-STATIC, ADV-ENERGY, ADV-REFUSE (decisions 8b, 9a; B5i, HYG-SETUP and
HYG-CHECKED verified on the advisor build alongside):
- shared `build/` remade, `make test` all suites green, `EXHALE.x` current.
- `make check` FAIL against 1e-3 as expected: marching files carry exactly
  gate 3's (COST7's) movement, log10 Mdot unchanged on all 14 cases; the
  `_adv` files move by design (eighth column `adv_status`; refused base
  cells at equilibrium).
- **Goldens refreshed 2026-09-08 19:05** by the advisor on the user's
  instruction (previous set in `golden_pre_cost7_20260908/`); **gate 5
  (`make check`, 19:03 to 22:02) PASS, 64 files byte-identical.** Since
  then B5j landed (`steady_newton.f90`: bound-aware region, the new
  elimination increment measure; `wasp_full_newton` moves 5.9e-9 / 9.6e-9
  against the refreshed golden, no refresh), B5k (decision 11a, the carrier
  unknown as ln n; `wasp_full_newton` bitwise against B5j) and B5l (the
  element budget as a face of the box; `wasp_full_newton` bitwise against
  B5k, marching bitwise on `mol_carrier`) landed with no matrix movement.
  The shared `build/` and `EXHALE.x` are BEHIND B5j/B5k/B5l (only the
  advisor's `build_lwv`/`EXHALE_lwv.x` is current, 2026-09-09 01:04); the
  `make test` suites will fail their staleness precheck until `make` is
  run. **Code work paused by the user's instruction of 2026-09-09.**

Gate 3 (2026-09-08, 09:29 to 12:00), tree after COST7 (decision 8b):
- shared `build/` remade, `make test` all suites green, `EXHALE.x` current.
- `make check` FAIL against 1e-3 as expected (the stopping test changed);
  every movement examined and at the tolerance level or a pre-existing base
  artifact (Update log, "Gate 3"); log10 Mdot unchanged on all 14 cases.
- COST7 suites on the advisor build: `coupled_source_step` 44/0,
  `physics_probe` 1491/0.

Gate 2 (2026-09-08, 05:04 to 08:08), tree as of B5e + HYG-4:
- `make check` REGRESSION PASS, all 65 files of the 16 cases byte-identical
  to the 2026-09-07 goldens; `wasp_full_newton` re-checked separately with
  the B5e binary (identical, info 0, CERTIFIED, log10 Mdot 13.30).
- every `make test` suite green on the advisor's private build of the same
  tree (the batch script's `make test` failed 7 suites on the build-staleness
  precheck alone, because the shared `build/` was remade only afterwards).
- COST6 (default off) verified on top: `coupled_source_step` 39/0,
  `physics_probe` 1491/0, `mol_base_handoff` 12000 and `wasp_full` at matrix
  length IDENTICAL to the goldens.
- `run_fcheck.sh` not re-run this batch.

Gate 1 (2026-09-07 20:54 to 2026-09-08 00:05), for the record:

- `EXHALE.x` md5 `eba2d7f2260f`, gfortran 16.2 (conda) + OpenBLAS,
  `make -q` up to date.
- `make test` green (the two steady suites need
  `EXHALE_STEADY_RUNLOGS=backup/regression/wasp_full_newton/run.log`; the
  `attempted_step` whole-binary rows need `EXHALE_ATTEMPTED_STEP_EXE`).
  `diffusion_tests` links and passes again (35/35) after DIFT-LINK.
- `make check`: REGRESSION PASS, 16 cases byte-identical, 2026-09-07 20:54
  to 2026-09-08 00:05. Goldens unchanged since the 2026-09-07 refresh
  (`golden_pre_series_20260907/` is the set before that).
- `run_fcheck.sh`: CLEAN as of the 2026-09-07 gate; not re-run for this
  batch (no new bounds-sensitive path landed; B5c and DIFT-LINK were run
  under `-fcheck` privately by their workers).

## What this batch settled

- One discretization of the material advective term exists in the code: the
  RK stages, the stationary carrier and element rows, the carrier relaxation
  and the element relaxation all form it as the divergence of the face
  species fluxes (B4-1, B4-1c, B5b, B5c).
- The stationary Newton system carries a row and an unknown per transported
  balance (B5, B5b); `Ionization transport` with `Solver: Newton` is
  supported by derivation.
- The thread count is no longer an input to any result (THREAD-DET: fixed
  even blocks, not a block count from `omp_get_max_threads`).
- The coupled loop's pass count is a LADDER OF NESTED MODES and equals the
  slowest of them (COST4, after COST3 and FIELD-SELF): 10.09 passes on
  `wasp_full`, 7.16 with the stellar field's composition dependence removed,
  4.94 with the temperature removed too, 2.00 with the He recombination
  rates as well. Every lever tried removes a subordinate rung and leaves the
  next standing, so nothing changed; the decision is to change nothing.
- The count is arithmetic (COST5): `2 + log10(excursion / tolerance) /
  (decades per pass)`, every term measured and the law tested against a
  tolerance scan; the loop is one to two decades stricter than its
  tolerances say because the test reads the increment and the error is
  |theta/(1 - theta)| times it. The user decided (decision 8, option b) to
  stop on the estimated error: item COST7. Memo:
  `docs/coupled_source_loop_cost.md`.
- The global-mode extrapolation (COST6) exists, default off
  (`EXHALE_CSM_EXTRAP=1`), worth 6.3 percent of the sweeps on
  `mol_base_handoff` under the present test.
- The first stationary solve with a species row meets its target (B5e,
  `mol_diffusion` element row, 4.3e-9 against 1e-8); the trust region is the
  step control of the coupled route; `cert_tol_element` has one measurement
  under it (4.6e-11) and is not moved.

## Open items, in order

(The problem-centered account of the whole of stage 2, with evidence, is
`docs/ISSUES_20260909.md`; this list is its short form. Its review,
`docs/ISSUES_20260909_review.md`, and the plan derived from that review,
`docs/PLAN_20260909.md`, its review `docs/PLAN_20260909_review.md` and the
revision `docs/PLAN_20260909_rev1.md` (items N-1 to N17, stages A to E,
decisions 13 to 19; rev 1 is the plan of record), supersede the ordering
below wherever they differ. Rev 1's preparation N-1 and N17 are DONE
2026-09-09: `backup/regression/atomic_elem_newton/`, `docs/worker_rules.md`,
`docs/audit_20260905/run_issues_review_20260909.out`, ISSUES section 3.8;
Stage A (N0, N3, N1, N2, N16a) is DONE and its gate PASSED 2026-09-09
(advisor build `build_lwv`, md5 `dd233d6a`; see the update log). B5m's
root cause was the absent elements Si, K, S registered as unknowns on their
bound; the atomic element solve now moves to `||R||` 7e-2 and oscillates.
Since then (all 2026-09-09, logged in the update log): N13 (transit census
mask), N12 (upstream caloric energy + the stale `x_h2_cell` defect; molecular
`_adv` goldens will move, refresh at the Stage D gate), N9 (restart He/H was
a message defect; He/H is moved by the coupled-carrier Newton), N7 (the
atomic stall is the ray test + a merit/gate mismatch, decision item 20),
N4a (the element inventory map, reporting only), N7b (sound ratio opens the
cycle; the limiter is now the linear solve). N5 done (replay GREEN on four
configurations; conditional closure exceeds the tolerances, two suite rows
FAIL by design; B5h's two roots do not reproduce with H2 transported).
Decisions 14 (route i) and 20 (one functional) taken 2026-09-09. N20 done (one functional: carrier
0.446, atomic arm collapses on conditioning), N16b/N16d/N16e done (species
rows in the estimate; the checkpoint restore complete; p = 3 for the RK3
rows, p = 1 for the split rows). N8a done (one empty Jacobian
column was the collapse; column floor + equilibrated preconditioner; atomic
reload 3.8e-2 at 150 and 5.3e-4 after a second 150, carrier 0.151). N4b done (element-conserving
write-back, He/H fixed at 1.4e-11; carrier reload 0.090; the projection
never fires on the fixtures). Advisor measured the atomic reload at cap 600:
one solve stagnates at 2.7e-2 where 150 + 150 reaches 5.3e-4, so a fresh
trust-region state helps. N21 done: no reset helps, the linear
solve is not the limiter, and the second solve's gain came from 2000
marching steps between the solves (the advisor's earlier reading was wrong);
the limiter is a predicted decrease that stops scaling below `||s||` 1e-7.
Decisions 15, 16, 17 taken (a).
N13b done (transit metadata). N11 done (mass-operator stationarity, two
fields); N11b running (conditional tolerance 1e-2 + the measure column,
advisor decision). N10 brief ready (restart contract), launches after N11b
frees `write_output.f90`. N22 done (face writes were the cause;
atomic reload 1.36e-2 in 113 iterations), N10 done (restart contract; decision
item 21 on option-changing restarts). N10b done (decision 21 a: `Restart
option change` key; the arm directories written before the RJ unification
are configurations, not restart sources). N23 and N24 done (the atomic arm's
obstruction is the sodium row of cell 500 and the linear model there).
**Stage D gate PASSED 2026-09-10 00:27 and the goldens were refreshed at
00:40** (previous set `golden_pre_stageD_20260910/`; marching files identical
on all 16 cases; `_adv` files moved by design). N25 done: the Arnoldi gap is the
nonlinearity of the finite-difference action (not orthogonality), the ball
truncation is a measured arm (default off), and the sodium row's upper ghost
is FROZEN at the pre-solve composition during a stationary solve (a defect
of `eval_residual`, not of the operator). N26 done: the upper ghost follows
the iterate; atomic arm 1.36e-2 -> 3.7e-4, binding row now the interior He
row of cell 246; carrier `||R||` 0.25 (same binding cell). N27 done: the helium row of
cell 246 is held by the linear solve (Krylov 40/40 at 0.99, Arnoldi gap
0.147 not tracking the step), not by physics, bound, closure or tolerance.
**N26b's ghost refill in the element operator BREAKS the atomic arm** (1.199
at iteration 40 where N26's arm reached 3.2e-2; advisor confirmed): N26c reverted it: the arm is chaotic
at the ulp level (a one-ulp ghost change breaks it as badly), so its
`||R||` at a fixed iteration is not an acceptance quantity; the linear solve
is the limiter (N21, N24, N27). The user chose the advisor's
recommendation: N8b done (decision item 22 written: 1e-8 is refused
by discretization and element flux conservation, not by any floor; the
species of the candidate do not carry its density, 6.7e-3, a defect of
`project_elements`); N29 done (mass created by
`project_elements` from the reservoir metal/H; fixed to round-off;
`lower_profile` moves, golden to refresh at the next gate); N28 and N28b done
(cleanup; three uninitialized derived-type arrays, a garbage logical
`x_h2_fixed`). Gate after N29 PASSED (15 identical, `lower_profile` moved by
design and its golden refreshed, previous in `golden_pre_n29_20260910/`). N8c
re-measured the anchors (reading unchanged). **Decision 22 taken: (a) and
implemented (N30)**: 1e-5 gating in the wind, reported below; no state
certifies. No item running. Next: the column preconditioner question (a),
the layer's element flux conservation, N6 steps 1-4, the manual PDF. The user manual PDF is rebuilt. Next: N8b (tolerances) once a
candidate certifies, N6 steps 1-4, the cleanup items of ISSUES 3.7, and the
PDF of the user manual (two keys behind its source).)

1. **The carrier arm of the stationary solve** (B5g, B5h, B5j, B5k;
   decisions 11a and 12c taken): the carrier unknown is `ln n` by default,
   `mol_carrier` `||R||` 1.415 -> 0.9957, no bound held, no Cauchy-only
   step; still not converged. **Next obstruction (B5k, MEASURED)**: the
   element budget (the carrier CEILING) shortens every step (157 refused
   samples, cut-backs of 3-5 halvings, radius at its ceiling), and it is
   deliberately not a face of the unknown box. Item B5l (DONE 2026-09-09):
   the budget is a face of the box (`freeze_species_unknown_box`), budget
   refusals 157 -> 0 and long cut-backs 38 -> 0, `||R||` 0.9957 -> 1.174 (a
   stall either way), `Mdot` unchanged; `gm_m` 40/80/160 measured equal
   (the cycle stagnates at 0.6 of its right-hand side), so the literal 40
   stays. **Next obstruction (B5l, MEASURED)**: the Krylov leg of the dogleg
   (banded preconditioner or operator), not bounds and not subspace. Code
   work is paused by the user's instruction of 2026-09-09; B5m stays
   planned.
   Superseded text follows for the record. the H2 carrier of cell 500 sits on its lower
   bound (-7e-19), the dogleg cuts every step by 2^-25..2^-49 against it,
   `||R||` stays at 1.87-1.99. Backward differences removed the sample
   refusals (49/0/0) but not the obstruction. Either a bound-aware step
   (projected trust region) or a reformulated unknown (log density, or the
   carrier fraction clamped as the element fraction already is in
   `write_species_rows_into_composition`). A decision on the unknown space,
   then an item. `cert_tol_carrier` and `cert_tol_element` UNCONFIRMED;
   B5e/B5f's element rung stands as a measurement of one state only.
2. **The ionization front has two exact roots** of the network at the same
   incident field on `mol_diffusion`'s state (x(H I) 19 percent apart at
   r = 1.086, both class 1 at 1e-25), reached from seeds far apart (B5h).
   Which root a stationary residual may land on, or a statement that the
   front is unresolved on this grid, is a physics rule the code does not
   have. Study item before any decision.
3. DONE (B5j): the elimination stops on the change of the three quantities
   the residual reads (particle count, H2 caloric share, radiative pair);
   `eq_sweep_reltol` 1e-12 reached in 33 passes on `wasp_full_newton`;
   `EXHALE_EQ_MEASURE=species` restores the old measure.
4. The O(dt) kick on a reloaded stationary state (B5f): two CFL steps
   raise `||R||` from 8e-9 to 1.2 before the Newton starts. Known
   mechanism (P54 AD), first seen from a state stationary to 1e-8.
5. Documentation debts found by the sweeps: the 63 stale "F<n>" line
   tokens of `docs/input_schema.md` (delete the convention); no
   `Ionization transport` entry in the user manual; two documents describe
   the `_adv` validity conditions (manual, HOWTO); four `# columns` parsers
   in python; `Coupled carrier solve: True` silently inert in an atomic run
   without `He_diffusion` (report belongs at the row registry).
6. `n_tot` lag on the molecular path (COST5): worth 0.8 passes at most;
   revisit only with a coupled (T, c) reformulation.

Done this session (2026-09-08): B5d/B5e verified, COST4, COST5, COST6 (off),
COST7 (decision 8b), HYG-4, HYG-T9, B5f, B5g, B5h, B5i (decision 10a),
ADV-STATIC, ADV-ENERGY, ADV-REFUSE (decision 9a), HYG-SETUP, HYG-CHECKED.
Decisions 8, 9, 10 taken by the user; the decision document has no open
item.

## User-gated (unchanged from 2026-09-07)

D-physical reruns of the planet folders and `benchmarks/` (every H-alpha and
H-beta transit number from a `Jlya escape-prob: True` run is superseded by
LYA-BETA, and Mg II by the metals work); `Load IC` case regeneration on the
current grid; a converged run reaching 3000 K for the CO domain; paper and
poster re-convergence, LHS 1140 b, the GitHub public switch.

## Not in git

As listed in the 2026-09-07 handoff, plus this file, `src/tests/steady_species_rows/`,
and the FIELD-SELF and B5 probes. `git status --short | grep '^??'` lists them.
