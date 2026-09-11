# Execution order for the development plan (rev 3 with section 10)

Companion to `docs/development_plan_20260905_rev3.md`. That document says what
and why; this one says in which order the work is done, what closes each
step, and who acts (user decision, advisor, worker). Item labels (A1, D0, C1,
10.1-1 ...) refer to rev 3.

Working method for every code item: the advisor writes a brief (files,
lines, the gate, the known traps), a worker implements it with its test in a
worktree, the advisor reads the diff, runs the item's gate and the bounded
regression (`run_check.sh check <cases>`), measures the movement on one atomic
and one molecular case, and records a dated section in `docs/Update_EXHALE_stage1.md`.
Nothing is accepted on the worker's report alone. Goldens are refreshed once
at the end of a phase, after the movement of every item in it has been
measured separately.

## Step 0. Decisions that gate the first code change (user, before Step 2)

| Decision (rev 3 number) | Needed by | Default the plan proceeds on if undecided |
|---|---|---|
| 1. Golden refresh: deferred to the end of Phase 1 (Step 2c); a scratch baseline, not `golden/`, is the control until then | Step 2c | as stated |
| 13. Sub-Lyman field source | Step 2c | DECIDED 2026-09-05: one spectrum type builds every band. Power-law input: power law everywhere (XUV and below 13.6 eV); Planck input (new `Spectrum type: Planck`): Planck everywhere; SED table: table everywhere. The Balmer continuum follows the run's type. Setup report states the sub-13.6 eV source and the integrated grid flux. |
| 14. He mass and Jupiter radius | Step 2c | DECIDED 2026-09-05: m_He/m_H = 3.9715; RJ = 7.1492e9 cm in code and bridge (R0 of every run +2.26%); MJ, MSUN to the IAU nominal values |
| 15. `du` functional | Step 2c | DECIDED 2026-09-05: signed spread over the mean, as `flux_spread` |
| 16. Every planet example and benchmark runs on the SED its reference literature used or assumed | Step 1.6, Step 2d | DECIDED 2026-09-05: inventory + file per planet under `inputdata/sed/`, inputs switched to `Spectrum type: Load`; reruns stay instruction-gated |
| 17. SED ending above the band an active absorber needs (He 2^3S 2583 A, low-IP metals) | Step 2a | DECIDED 2026-09-05: always STOP; the message offers the two remedies (turn the triplet off or remove the metal in the input, or supply a covering SED); no continue key; replaces the warning-only branch of `sed_read` |
| 3. A2 interface: status return and vacuum branch in `Num_flux` | Step 2a | DECIDED 2026-09-05: the default (status + edge speeds; HLLE vacuum branch); `docs/a2_roe_interface.md` APPROVED 2026-09-05, implement as written |
| 4. A3 retry design: carrier substep retry or outer step controller | Step 2b | increment (a) only until decided |

Decisions 2, 5, 6 (energy convention, charged closure, CO model) are taken
inside D0 (Step 3); 7 (H3+ emission model) before B2 (Step 4); 10 (atomic
ionization transport) in Step 6; 8 (author request) any time; 9 and 12 as
stated in rev 3.

## Step 1. Phase 0: evidence and gates (advisor + workers in parallel, about one week)

**DONE 2026-09-05** (`Update_EXHALE_stage2.md` sections 2 and 3).

1. Baseline without touching `golden/`: the post-170 outputs already in the
   case directories are copied once to
   `backup/regression/baseline_post170_<date>/` (a scratch reference, not a
   golden of record), and `run_check.sh` gains `REGRESSION_GOLDEN_DIR` so a
   `check` can compare against that directory. The pre-170 goldens stay in
   place, untouched, until Step 2c. Reason: a golden refresh before the
   Phase 1 fixes would be redone after them; the only thing Phase 1 needs is a
   byte-identity control, and the scratch baseline is that control.
2. Worker W1: `src/tests/physics_probe/` derived from `docs/audit_20260905/`
   (which stays): LLF/ROE Riemann tests including the equal-pressure and
   vacuum-producing expansions; Chung detector inversion; H3+ surface checks;
   photoevent energy sums on the accepted equation; H2 thermodynamic identity.
   Each with a tolerance and a nonzero exit on failure.
3. Worker W2: the new tests of section 10.4: grid-width identity `dr_j(j) =
   r_edg(j) - r_edg(j-1)`; photon quadrature of P_HI, P_HeI, P_HeII against a
   fine reference with and without the sub-Lyman grid; `du` functional on a
   sign-changing window; output header against a re-sweep; base level stated
   once; restart onto a different `Grid type` refused.
4. Worker W3: regression matrix additions: `hp_zero_seed`, `hp_trace_seed`,
   `hp_front` (from verified restart files, state asserted before
   `carrier_source`); `hydrostatic_column` (mechanical); one oxygen-chemistry
   case from `examples/18_oxygen_chemistry`; `_adv` files added to `FILES` in
   `run_check.sh`; `make test` target aggregating everything; `CASES=` on
   `make check`; the `run_check.sh` build line made fatal with its log kept.
5. Advisor: the A2 interface (decision 3) written as a signature and shown to
   the user; the Phase 1 briefs.
6. Worker W4: the SED inventory of decision 16: for each planet of record the
   SED its reference paper used or assumed, read in the PDF (title, proxy
   star, scaling, coverage), the file obtained or built into
   `inputdata/sed/` with a README, coverage extended to the 4.8 eV floor
   where He 2^3S is on (the Koskinen solar file stops at 1995 A); the
   unreferenced `inputdata/scaled_solar_hd189.sed` either adopted with its
   provenance or removed.

Gate: `make test` exists and is RED on exactly the known defects (LLF, ROE,
detector ratio, H3+ joins, `dr_j`, quadrature, `du`, header); `check`
against the scratch baseline is byte-identical (nothing has changed yet);
the new cases have baseline outputs in the same scratch directory.
**Gate met 2026-09-05.** `src/tests/physics_probe/` 41 assertions, 15 pass,
26 fail; `src/tests/grid_and_gates/` 6 programs, all RED; `make test` runs five
suites; `check wasp_full mol_base_handoff` against the scratch baseline
byte-identical on all four files per case; the five new cases reproduce their
own baseline copies. `element_census_tests` did not compile at 35d9dd5 and was
repaired the same day (20 assertions pass).

## Step 2. Phase 1: bounded corrections (workers in parallel, advisor verifies; two to three weeks)

Ordered by golden movement so that each movement is attributable.

2a. **DONE 2026-09-05** (`Update_EXHALE_stage2.md` section 4; seven items, seven
    workers on disjoint files).
    No golden movement expected (opt-in paths, dead code, guards):
    A2 ROE with the caller contract; opacity model P applied to the
    cross sections (10.2-7); Lyman-alpha photolysis attenuation through the
    existing `lya_rt` beam (10.2-8); Wind-AE bridge constants and He mass
    conversion (10.2-9, 10.2-10); `x_hp_fixed` honored in the constrained
    continuation (10.2-12); the P2 guards and dead code of 10.3 (minloc,
    `j_min` lower clamp, `get_word`, `report_unmet_steady_gates`, Mdot index,
    conduction floor counter, `okres`, `Tcode`, `carrier_mass_amu` by name,
    `adv_corr` zeroed on allocation, `oxygen_rates` header, the FUV shell
    comment, dead constants); decision 17 in `sed_read` (stop by default when
    the file ends above the band an active absorber needs; the message
    offers the two remedies, turn the absorber off in the input or supply a
    covering SED; no continue key). Gate: `REGRESSION_GOLDEN_DIR=<baseline>
    run_check.sh check` byte-identical on the eleven original cases; the
    item's own test green.
    **Gate met 2026-09-05, with one expected exception.** `make test` on the
    advisor's single rebuild runs seven suites: every assertion added in this
    batch is green and only the known defects owned by later items remain red.
    Regression, sixteen cases against the scratch baseline, single-threaded:
    fifteen byte-identical on all four files (the eleven of `DEFAULT_CASES`,
    `hydrostatic_column`, the three `hp_*`); `oxygen_chemistry` moved, as the
    Lyman-alpha item predicts, and its scratch baseline was re-snapshotted with
    the previous copy kept.

2b. **DONE 2026-09-05** (`Update_EXHALE_stage2.md` section 5; five code items and
    one documentation pass, disjoint files, one advisor rebuild and gate).
    Bounded, moves some goldens:
    A1 LLF speed with the four counters (positivity path reached or not);
    A3 increment (a): element reference table, residual on the returned
    state, status, terminating failure; A8 detector inversion with the `q_D`
    admissibility rule, the rename and the integration check; `k_lw` cell
    mean (10.3) as a first D1 increment if the worker can do it in the same
    routine; He I singlet density in `eval_cool`; He II recombination-cooling
    coefficient of the triplet path. Gate: item test green; movement measured
    on `wasp_full` and `mol_base_handoff`, recorded.
    **Gate met 2026-09-05.** Sixteen cases against the scratch baseline:
    `wasp_he23off` and `hydrostatic_column` byte-identical, fourteen moved, Mdot unchanged to two decimals in all sixteen; the baseline was
    re-snapshotted with the post-2a copy kept (section 5 of the log has the
    largest relative change of every case).

2c. **Code DONE 2026-09-05** (`Update_EXHALE_stage2.md` section 6; fourteen items,
    every one measured alone, one advisor rebuild, `make test`, sixteen-case
    gate against the post-2b baseline, `hp_*` IC regenerated on the new grid,
    baseline re-snapshotted, `golden/` archived to `golden_pre170_20260905/`
    with its README). **Golden refreshed 2026-09-06 on the user's
    instruction** (sixteen cases from the 89-bin binary, then
    `oxygen_chemistry` after the carrier Newton budget was raised; the
    loaded-SED grid (D0 item C5) and the D0 hygiene items landed the same
    day; `make check` and `run_fcheck.sh` close the phase).
    Moves every golden (each measured alone first, then all applied):
    10.1-1 `dr_j` shift; 10.1-2 threshold bins on edges (or physical-fraction
    weights); 10.1-3 as decision 13 (`Spectrum type: Planck` added, the
    Balmer continuum on the run's type, setup-report lines; the power-law
    default is unchanged in the numbers, so this item alone moves no golden
    of a power-law case; **the `Spectrum type: Planck` branch and the
    setup-report lines landed early, in 2a, so what is left here is the Balmer
    continuum following the run's type**); 10.2-1 momentum kick
    removed (denominator floored); 10.2-2 one state per output file (re-sweep
    before the final write, header describes the relaxed state); 10.2-3
    restart keeps the run's grid or refuses; 10.2-4 base level from
    `p_match_bar` on a profile run; 10.2-5 signed `du` per decision 15;
    10.2-6 C I Voronov row; 10.2-11 He mass per decision 14.
    Procedure for each: branch with only that change, run `wasp_full` and
    `mol_base_handoff`, record Mdot and the T, v, n_e, x(H+) profiles against
    the scratch baseline (plotted), then merge. After all: copy `golden/`
    to `golden_pre170_<date>/` with a README (provenance from the headers,
    unknowns marked unknown), `run_check.sh check` against the baseline to
    record the total movement, `run_check.sh golden <cases>`, `make check`
    again, `run_fcheck.sh` once, one dated section listing the movement of
    each item. This is the one golden refresh of Phase 1.

2d. Record, without re-running: a note in the READMEs of `HD209458b/`,
    `HD189733b/`, `WASP-121b/` that their He 10830 depths were computed under
    the power-law spectrum type, whose sub-13.6 eV field is the power law
    (decision 13), and with the shifted `dr_j` (10.1-1). Decision 16: every
    planet folder and benchmark switches to `Spectrum type: Load` with the
    SED its reference literature used or assumed, from the inventory of
    Step 1.6; the switch is a configuration change recorded in the folder
    README, and the rerun stays instruction-gated.

Gate for the phase: `make test` green except the items owned by Phases 3-5
(H3+ surface, closed cell, photoevent sums); goldens refreshed once with
provenance; the movement table in `Update_EXHALE_stage1.md`.

**From 2026-09-06 the order of Steps 3 to 7 is superseded by
`docs/PLAN_20260906_rev2.md`** (rev 0, rev 1 and their reviews are kept). Phase 2 gate met 2026-09-06: B1 accepted (`docs/b1_target_system_20260906.md`) (acceptance contracts first, then the target
specification separate from the D0 audit, then the coupled source step; the
spectra and configuration switches run independently). The text below is kept
as the original estimate.

## Step 3. Phase 2: D0, the governing system on paper (advisor writes, user accepts; two weeks)

Contents as rev 3 section 4.2 items 1-10 with the amendments of section 10
(cell-mean photon field, one field for 4.8-13.6 eV named with its source,
independent species space for C1, rollback contract). Decisions 2, 5, 6
resolved in the text; the two-cell radiation ownership rule and the
thermal/chemical flux split written as equations; the verification matrix
with units and tolerances. Gate: the eight acceptance questions of review 3
section 5 answered, and the user's acceptance. No Phase 3 code before it.

## Step 4. Phase 3: sources and thermodynamics (sequential increments; six to ten weeks)

1. A4+B3 increment (i): scalar T solve with residual test, rates at the
   returned T (replaces `energy_semi_implicit`). Gate: energy-solve test;
   movement measured.
2. A4+B3 increment (ii): equilibrium composition solved with the energy
   identity of D0 item 3 (the fixed-T projection of N1 removed). Gates:
   isolated photoionization at threshold and 20 eV; recombination with
   escaping energy; photon-driven cycle; closed reacting cell; two-cell
   exchange; one-cell photon balance; the four Newton-finished planets before
   and after, plotted.
3. A5 photoevent ledger (double channel sum, threshold region, excited
   fragments, stopping lengths) and A6 reaction heats (associative term,
   oxygen ledger or label, R5/R16 product state). Gate: E1 energy sums green.
4. B1 single H2 model with standard-state functions and population
   consistency; B4 rate ranges. Gate: identity and cycle tests.
5. B2 H3+ after decision 7: emission model, temperature coverage, linear
   low-density asymptote, high-density end, continuity. Gate: the four
   independent checks of section 6.
6. A4+B3 increment (iii): finite-rate integration of the non-eliminated
   species with the rollback contract (if A3 increment (b) has not landed,
   it lands here). Gate: source-step order measured; rejected trial never
   written back.
Goldens refreshed once at the end; movement recorded.

## Step 5. Phase 4: the spatial operator (four to six weeks)

D3 boundary count documented and adjusted; D1 face-flux species advection
with the thermal enthalpy flux, the chemical face flux, the charged closure
of D0 item 5, and the cell-mean photon field of every absorber (the
inner-face column of `calc_column_dens` and `k_lw`); D2 discrete hydrostatic
residual first on the mechanical column, then the base cycle on
`mol_sec_ion`; validity diagnostics written per cell. Gates: passive scalar
order and budget; diffusion-only column closing mass, elements, current and
thermal plus chemical energy; hydrostatic residual decreasing with
resolution. Goldens refreshed once.

## Step 6. Phase 5: the steady solver (six to ten weeks)

C1 species rows on the independent space of D0, unsplit stationary residual,
refusals lifted one at a time (the `x_hp_fixed` fix of Step 2a is a
prerequisite); C2 relaxation-time elimination and He+ transport; C3
composition convergence gate and reload test; C4 finite-difference action
verified over `eps`, inner closure, preconditioner omissions, the
`wasp_full_newton` base energy row. Decision 10 taken here. Gates: Newton
solution against a time-converged reference; `hp_*` cases; option-off path
byte-identical to the preceding build.

## Step 7. Phase 6: validation and extension (open-ended)

Koskinen Model A term by term (the local paper suffices to start; the
spectrum request is decision 8); lower/upper matching sweep; `_adv`
consistency; band photon conservation; CO kinetics replacing the ceiling (done 2026-09-06,
items B3b-CO, B3b-CO2 and CEILING-DEL); molecular heat capacities; net-exchange infrared form; the Koskinen 2009
H3+ table as a keyed option. Only after these do the published profiles
become acceptance targets. The instruction-gated tail (paper and poster
re-convergence with the corrected sub-Lyman field, LHS 1140 b, public
switch) follows only on explicit instruction.

## The first three actions

1. Decisions 3, 13, 14, 15 are taken (2026-09-05); decision 4 (A3 retry
   design) remains open and only blocks A3 increment (b). No golden action is
   needed now.
2. Advisor: briefs for W1, W2, W3 and the A2 signature; launch W1-W3.
3. Workers: Step 2a items can start the same day as W1-W3, because they move
   no golden and need no D0.
