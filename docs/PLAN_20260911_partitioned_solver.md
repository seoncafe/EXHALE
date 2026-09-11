# Plan of 2026-09-11: repair the partitioned stationary route and measure it

User instruction (2026-09-11): read `docs/solver_approach_analysis_20260910.md`,
its review `docs/solver_approach_analysis_20260910_review.md` and
`docs/solver_partition_experiment_20260911.md`; correct the errors they
establish; run the solver experiments; adopt the partitioned route in the
production solver if it measures better than the coupled species-row solve.
This supersedes the 2026-09-10 instruction that closed the species-row
program after N37/N38.

Rules for every worker: `docs/worker_rules.md`. The plan as written names
items P1 to P7; P8 to P13 were added as it ran and after it, and P14 to P19
in the afternoon of the same day, on the two findings of a Codex adversarial
review and on the user's two decisions (see the status section). Scratch root
for this plan: the session scratchpad, `<scratchpad>/P<n>/`.

## Phase 1: the confirmed defects (parallel, one file owner each)

| item | defect | file owned | tests owned |
| --- | --- | --- | --- |
| P1 | D1: a failed element composition solve is handed back unchecked (`solve_mass_fraction` accepts its last halved trial; `element_diffusion_step` clips and continues; `relax_element_composition` exits on movement, never on a solve result). Reproduced: 17.7 percent mass-closure error handed to the chemistry, which went nonfinite. | `src/modules/functions/binary_element_diffusion.f90` | `src/tests/element_operator/`, `src/tests/diffusion_tests.f90` |
| P2 | D2: `relax_photochemical_composition` tests `drift > trust` AFTER the step is applied and exits without restoring it; `trust = 0.01` returned drifts of 0.028. | `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | `src/tests/carrier_retry/`, `src/tests/carrier_returned_state_acceptance/` |
| P3 | D3: the loop-top stop of `solve_steady_jfnk` is `steady_gates_met` (`\|\|R\|\| < Resid tol`, flux, chemistry, carrier gate) while the final acceptance additionally requires every certified row within its own tolerance (mass 3e-12, momentum 1e-8, energy 1e-6, species 1e-5 in the wind); a solve stops early on the first and is refused by the second (`info = 2`, mass 5.95e-12 against 3e-12 in the molecular partitioned run). | `src/modules/time_step/steady_newton.f90` | `src/tests/krylov_and_dogleg/`, `src/tests/steady_completion_flag/` |
| P4 | Under `Well balanced: True` `store_row_terms` scales the momentum row by `max(\|dF(2)\|, \|S(2)\|, \|Smom\|)` with `S(2) = 0` by construction, so every scaled momentum residual is exactly 1 (N37, ISSUES 3.7). The reference must be the pressure force the arm cancels analytically, `\|rho_j (A+ (phi_i(j) - phi_c(j)) + A- (phi_c(j) - phi_i(j-1)))\|/dV`, beside the dynamic and other source terms. | `src/modules/time_step/steady_residual.f90`, `src/modules/time_step/RK_rhs.f90` | `src/tests/grid_and_gates/` (`hydrostatic_residual`, `momentum_row`) |
| P5 | The analysis document's errors the review established (truncation error read as a residual floor; Ritz ratios read as conditioning; "removes the coupling"; the outer loop said to be new; the hydro finish attributed to CETIMB; well-balanced said to remove the floor; Koskinen 2013 lower boundary is section 2.1.1; the three-unknown path does not use the trust region). | `docs/solver_approach_analysis_20260910.md` | none |

## Phase 2: the outer iteration as a stated contract (P6)

`steady_wind_with_element_diffusion` (`src/EXHALE_main.f90`): consume the
P1 status; a controlled composition update; accept an outer state only by
the certification of the FULL equations on the refreshed state; explicit
failure at the pass budget; no exit on a hydro `info = 2` that P3 has made
unreachable except through the chemistry at return. Default paths (no
element diffusion, no transported carrier) unchanged bit for bit.

## Phase 3: the measurement (P7)

Both fixtures (`backup/regression/atomic_elem_newton`,
`backup/regression/carrier_elem_newton`), same binary, same thread count,
same initial state: the coupled solve (`Coupled carrier solve: True`)
against the partitioned route (`False`), judged by the certification of
the returned state on every active equation, the worst wind species row,
and wall time. The decision rule: the route is adopted as the recommended
configuration only if it certifies where the coupled solve does not, or
reaches a strictly lower joint residual in no more wall time.

## Status at the close (2026-09-11)

Two items were added to the plan while it ran: **P8**, propagating the
established corrections into the documents and tests that carried the same
errors, and **P9**, writing the record. Four more were added AFTER it,
on the user's further instruction of the same day to fix the open items the
plan had left to the user ("fix all of them", and "the pressure term too"):
**P10**, **P11**, **P12** with **P12b**, and **P13**, the record of those
three. Six more were added in the afternoon: **P14** and **P15** on the two
findings of a Codex adversarial review of the working tree, **P16** and
**P17** on the user's first decision (the mass tolerance re-anchored above the
fixture's measured rounding floor) and its consequence for the best-iterate
ledger, **P18** on the one entry point P17 left, and **P19**, the record of
those five. The user's second decision, the long continuation of the carrier
reload, was run by the advisor and is recorded with them. The full dated
account, with every number and its source log, is `docs/Update_EXHALE.md`
section 8; the state to restart from is `docs/session_handoff_20260911.md`.

| item | state | measured outcome |
| --- | --- | --- |
| P1 (D1) | DONE | The reproduction goes from a returned composition missing rho by 1.768571764708571E-01, a nonfinite chemical refresh and exit status 1, to the entry composition restored at 1.046888214610824E-15 and an explicit named failure. `element_operator` 14/0 to 24/0, 22/2 on the RED build. `mol_diffusion` and `lower_profile` byte-identical apart from the provenance timestamp. The first inadmissible operation is now located: the composition solve leaves `[0,1]` and the range clip absorbed it. |
| P2 (D2) | DONE | The advertised bound now holds on the returned state and scales: drift 2.8126757844770485e-02 at `trust = 0.01` becomes 9.9956436739254070e-03, and 2.9845242387380568e-03 at 0.003. `carrier_retry` 111/10 RED to 122/0. `mol_carrier` byte-identical apart from the provenance timestamp. |
| P3 (D3) | DONE | The loop-top stop now reads the certified rows by the same evaluation as the acceptance at return. On the reproduction the solve no longer claims `info = 0` at iteration 14 to be refused at return on a mass row of 5.954e-12: it continues, 27 of 40 iterations meet the gate with a row refusing, and it ends `info = 1` on its budget at 5.945e-12. `krylov_and_dogleg` 279/0 to 290/0; `steady_completion_flag` RED to GREEN on both real logs. `wasp_full_newton` unchanged, `info = 0`, `\|\|R\|\|` 1.667e-09. Cost 0.44 percent of one outer iteration, paid only behind the gate. |
| P4 | DONE | The momentum reference scale under `Well balanced: True` is the pressure force of the cell's own hydrostatic equilibrium beside the dynamic and remaining source terms. With the key on `wasp_full_newton` CERTIFIES (`info = 0`, `\|\|R\|\|` 9.719E-10, momentum row 9.719E-10 against 1E-08) where it stood at exactly 1.000E+00. `hydrostatic_residual` 28/12 RED to 40/0; the suite 162/0 to 178/0. Key off byte-identical. Negative result recorded: the atomic reload with the key on stops earlier (`info = 2` at iteration 24, 1.702) than the control (1.415 at cap 40), while the carrier reload improves from 1.721 to 1.026E-01. |
| P5 | DONE | All nine corrections in place with inline marks and their sources, +463/-115 (114 to 592 lines). Three numbers of the original disagreed with their sources and are corrected; one was wrong (the N4b row read 1.4 to 0.09 where the measurement is 0.151 to 0.090). |
| P6 | DONE | `steady_wind_with_element_diffusion` accepts a state only by the certification of the full set of active equations on the refreshed state; a hydrodynamic `info = 0` with an alternated species row refusing is not accepted; a nonzero hydrodynamic flag does not end the loop; five endings are named; the last pass takes no composition update; a pass that fails to move the joint measure halves both updates and three such passes end the loop. Default path byte-identical (`wasp_full_newton` `info = 0`, `\|\|R\|\|` 1.667E-09, `log10 Mdot` 13.30, `Ion_species.txt` and `Hydro_ioniz.txt` identical between the control and the new build). |
| P7 | DONE, and the route is adopted | Both fixtures, one binary, `Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`, 8 threads, 500 cells. Neither route certifies either fixture. The partitioned route reaches a strictly lower joint residual in no more wall time on both, which is the plan's second adoption branch: atomic 12 passes in 1043 s with every elemental wind row within 1e-5 and only the hydrodynamic mass row refusing (1.054E-09 against 3.0E-12), against the coupled control's three passes in 2392 s at a worst gated row of 8.29E-05; carrier 8 passes in 173 s strictly below the coupled control on every judged row in 321 s, and 40 passes in 831 s taking the H2 wind row from 7.38E-02 to 2.23E-02. |
| P8 (added) | DONE | The withdrawn truncation-as-floor inference corrected in the three documents that still carried it, with decision 22 (a) kept and only its justification changed; `steady_selfconsistent_residual` assertion 2 corrected (RED to GREEN on the real log); the synthetic column of `relaxation_drift_covers_every_element` made mass-closed (1.66042E-02 to 6.58484E-16), the suite 194/0 to 195/0. |
| P9 (added) | DONE | This section, `docs/Update_EXHALE.md` section 8, `docs/session_handoff_20260911.md`, the dated sections of the two fixture READMEs, the reopening paragraphs of `docs/code_status_20260910.md` section 3 and `docs/ISSUES_20260909.md` section 5, and the new `TO_BE_DONE.md` entries. |

| P10 (added after) | DONE | `distance_from_certification` now reads the hydrodynamic slot as a distance, `max(mass/3e-12, momentum/1e-8, energy/1e-6)` formed by the one new public `hydrodynamic_distance_from_certification`, so `d < 1` is the row condition the loop-top stop and the acceptance at return impose, and the ledger ranks on it on every route. `krylov_and_dogleg` 290/0 to 309/0, with 300/9 on a RED build reading `maxval(rc)/1e-8`. `wasp_full_newton` reloaded is unmoved, all eleven outputs byte-identical, `info = 0`, `||R||` 1.667E-09, `log10 Mdot` 13.30. The re-ranking is visible on the carrier reload, which from pass 2 hands back a lower mass row at a higher `||R||` (4.85E-12 at 2.459E-09 against the control's 6.98E-12 at 1.593E-09). The HD 209458 b element reload still refuses on the mass row alone after twelve passes, 8.893E-10 at cell 3 in 1064 s against the control's 1.054E-09 at cell 16 in 1043 s, and the measurement says why no ledger can move it: the ledger's best distance of pass 1 is a mass row of 1.32E-09 where the same state re-evaluated carries 4.11E-09, so that row is the N33 rounding of the base-layer flux difference and ranking by it ranks rounding noise. |
| P11 (added after) | DONE | The momentum row's reference scale is the largest PHYSICAL term of `d(rho v)/dt + div(rho v v) + dp/dr + rho dphi/dr = S_visc`, the three terms kept as module arrays of `RK_integration`, filled on both reconstructions, under both key values and on the positivity-repair path, and gathered by one routine; the residual is untouched and the three terms add up to the row they scale (MEASURED 0 under WENO3, 1.7e-16 default and 5.3e-15 key-on under PLM). The PLM zero-gravity artifact, a scale of exactly 1.000 of `2 p/r` where every physical term is zero, now reads the `tiny` floor; the supersonic scale agrees with the ram divergence to 1.9e-14 where it stood 6.67E-02 above it; on the discrete equilibrium it was 2.79E-01 below the weight it balances and is now never below it. `hydrostatic_residual` 50/10 to 60/0, `grid_and_gates` 188/10 to 198/0. `wasp_full_newton` stays CERTIFIED with the key off, the states apart by at most 2.6E-09 relative and the momentum row 8.469E-10 to 1.165E-09, and with the key on, momentum row 9.719E-10 to 9.621E-11 and all eleven outputs byte-identical; `mol_carrier` marching byte-identical. |
| P12 with P12b (added after) | DONE | The `build_stamp` fallback was written where the compile cannot find it, which fires for any work directory outside `$root/build`, that is for every concurrent worker. Seven `run.sh` now take `${EXHALE_OBJDIR:-$root/build}/build_stamp.f90` and compile whichever stamp is in force as the first file of the closure, and four gained the `EXHALE_TEST_OBJDIR` override without which the no-build arm cannot be run on them at all. MEASURED RED 0/1 with exit 1 and `Cannot open module file 'build_stamp.mod'` in six suites, GREEN at 70/0, 36/0, 14/0, 10/0, 122/0, 9/0 and 8/0 after, in a tree with no `build/` and in the repository, with the PASS/FAIL lines unchanged. The mass-closed synthetic column, written twice, is now the one module `src/tests/test_columns.f90` (125 lines): `element_operator` 24/0 and `steady_species_rows` 195/0 before and after with both full logs byte-identical. |
| P13 (added after) | DONE | This status section, the P10 to P12 bullets of `docs/Update_EXHALE.md` section 8 with the "what did not certify" and "reported and not fixed" statements brought up to date, `docs/session_handoff_20260911.md`, the closed entries of `TO_BE_DONE.md`, the P10 and P11 re-measurements in the two fixture READMEs, and `docs/ISSUES_20260909.md` section 3.7. |

| P14 (added, Codex finding 1) | DONE | Every ending of the outer iteration hands back one consistent state. On a stagnation ending the composition was one update ahead of the particle densities and the temperature written beside it: re-evaluating the written state moved its temperature by 1.317089E-08 relative, against 5.921E-16 for the same fixture ending on its pass budget and 6.499E-14 after the fix. The progress control now stands between the certification of the pass and the update, the ending is announced with the other endings, and the element path refreshes the particle densities and the temperature before its equilibrium sweep. `grid_and_gates` 199/1 RED to **200/0** with the new round-trip arm at 1e-12. `wasp_full_newton` reloaded unmoved (`info = 0`, CERTIFIED, `\|\|R\|\|` 1.495E-09, `log10 Mdot` 13.30); the carrier reload identical pass line by pass line; the element reload parts from pass 2 and both arms end the same way. |
| P15 (added, Codex finding 2) | DONE | The kind-generic rows hand back the face departures their momentum row was built from, and the reconstruction continuation is refused on that arm. `krylov_and_dogleg` 309/0 to **321/0**; RED 2.098877E+01 (PLM) and 2.032204E+01 (WENO3) on the pressure gradient against the production term, 3.262392 and 3.428893 on the row scale, 1.190363E+01 and 1.181676E+01 on the row rebuilt from the stored face data, GREEN 0 bitwise for the generic-double instantiation and 2.6E-16 to 2.7E-15 for the quadruple one. `wasp_full_newton` unmoved with the key off and CERTIFIED at the same momentum row 9.621E-11 with the key on; the one quantity that moves is the iteration-1 measure inside the loop, 1.756E-03 to 1.753E-03; `mol_carrier` marching byte-identical. Third entry of the earlier "what remains" list, closed. |
| P16 (added, user decision 1) | DONE | The continuity row's tolerance is `max(3e-12, min(1, 10 * floor(j)))` with the floor `eps rho(\|v\| + c_s) A / dV` over the row's own scale, the verdict taken cell by cell. The brief's first floor expression was inert, MEASURED 4.4409E-16 on all 500 cells of all four states; the adopted estimate bounds the measured one-ulp step within 1.5931, 1.8857, 2.6806 and 1.7733, and `c_round = 10` is 3.7 times above the largest. **The HD 209458 b element reload CERTIFIES at pass 12**, where the control ends the same twelve passes NOT CERTIFIED on the mass row (7.159E-10 against 3.0E-12 at cell 13); `wasp_full_newton` unmoved, its floor below 3E-13 at every cell. `certification` 57/0 to **70/0**, RED 0/1. First entry of the earlier "what remains" list, closed. |
| P17 (added) | DONE | The best-iterate ledger reads the mass row against the cell's own tolerance, through P16's `mass_row_cell_verdict`, over the whole physical column: the ledger's best judged distance over the twelve HD 209458 b passes runs 2.2E+02 to 6.3E+05 in the control and **1.321E-01 to 1.868E-01** after, while the acceptance certifies the pass-12 state in both arms at a mass-row distance of 0.13 to 0.14. On the carrier reload the hydrodynamic solve returns `info = 0` at all three passes against the control's 0, 2, 2, and the refusing entries fall **2 to 1**. `krylov_and_dogleg` 321/0 to **334/0**, 327/7 on a semantic RED build. Cost 1.206 s against 1.216 s an outer iteration. A separate fix in the same file: the mixed-cell refusal message now names the binding cell throughout. |
| P18 (added) | RUNNING at the time of this record | One entry point for the hydrodynamic distance: the row of `src/tests/steady_species_rows/steady_species_rows_tests.f90` that keeps the bound form `hydrodynamic_distance_from_certification` alive moves onto `hydrodynamic_distance_from_certification_by_cell`, and the bound function and the two rows that read it are deleted. Nothing in the solver changes. |
| P19 (added) | DONE | This status section, the P14 to P18 bullets of `docs/Update_EXHALE.md` section 8 with the two new subsections (the Codex review, the user's two decisions) and the "what did not certify" and "reported and not fixed" statements brought up to date, `docs/session_handoff_20260911.md`, the `TO_BE_DONE.md` entries, the "After P16 and P17" sections of the two fixture READMEs, and `docs/ISSUES_20260909.md` sections 2 and 5. |

Two changes outside the item list were made by the advisor while reading
P10 to P12: the two rows of
`src/tests/steady_species_rows/steady_species_rows_tests.f90` that stated the
superseded ledger semantics are replaced by
`hydrodynamic_distance_reads_each_row_against_its_own_tolerance` (a mass row
of 1e-9 beside a momentum row of 1e-13 and an energy row of 1e-8 gives
`1e-9/3e-12`), and the `Coupled carrier solve` key comment of
`input_read.f90` no longer states that the alternation cannot converge.

One change outside the item list, made by the advisor while P6's runs were
being read: the `Restart intent: stationary` route now starts its
pseudo-time at 1.0, the marching hand-off's value, instead of the CFL
interval of the reloaded state, with `EXHALE_PTC_DTAU0` as the override and
the direct steady route keeping its CFL start. It was the hidden obstruction
of every partitioned run before it: at the CFL start the carrier reload's
hydrodynamic solve is handed back at `||R||` 2.685E-01 after 40 iterations
and the atomic one stagnates at 1.490, while at 1.0 the same solves hand
back 1.498E-09 and 1.167E-08. No regression case sets the key.

Three changes outside the item list were made by the advisor in the afternoon,
while reading P14 to P17: `assert_written_state_is_the_accepted_one`
(`src/EXHALE_main.f90`) carries an absolute floor of 1e-12 beside its 1e-10
relative test, which P14 had reported as a false alarm on a state whose mass
flux is flat to fourteen digits (3.0982E-14 as written against 2.7487E-14
accepted); the three `grid_and_gates` shell rows that wrote to a fixed
directory now honor `EXHALE_TEST_OUT` (`momentum_row_from_fluxes_only.sh`,
`base_level_single_statement.sh`, `sed_coverage_stop.sh`), so two concurrent
runs of the suite no longer collide; and the paragraph at `assemble_residual`
(`src/modules/time_step/steady_residual.f90`) stating what the arm does not
hand back is rewritten to what P15 made true.

### What remains

One item, and it is two questions about the hot-Uranus carrier reload's H2
front. Everything else of the earlier list is closed: the HD 209458 b mass
row by P16 (the fixture CERTIFIES at pass 12 under anchor (6)), the
kind-generic residual arm by P15, the `distance_from_certification` ledger by
P10, the default path's PLM geometric pressure term by P11, the `build_stamp`
fallback and the single home for `mass_closed_column` by P12 with P12b, and
the `Coupled carrier solve` key comment by the advisor, all on 2026-09-11.

The user's second decision ran the continuation to its stop. MEASURED (the
advisor's snapshot build of the merged tree after P1 to P12 and before P14 to
P17, `EXHALE_OUTER_PASSES=400 EXHALE_JFNK_MAXIT=40`, 8 threads): the progress
control ended the run at pass **110**, the H2 front had NOT stopped (`x2 =
0.5` moving 1.0726 to 1.1230 R_p over 108 passes at about 5e-4 R_p in a pass),
and from pass 40 the worst gated row was no longer at the front but stood at
**2.2e-2, nearly uniform in r**, its cell walking outward through the wind 443
to 500 and reaching the outer boundary at pass 108 while the value did not
fall. The table and the two-regime reading are in `docs/Update_EXHALE.md`
section 8, "User decisions of 2026-09-11 (afternoon)", and in
`backup/regression/carrier_elem_newton/README.md`.

1. **Does the front stop at all inside 2 R_p?** Physical, not a solver
   question. Its route is the P23 comparison against Koskinen et al. (2022)
   Model A, not a longer continuation of the present operator.
2. **A relaxation step for the wind cells' H2 balance, distinct from the
   bounded front advance.** The bounded carrier pass, trust 1e-2 on the H2
   mixing-ratio maximum and met by one transport step, advances the front but
   leaves the wind cells, whose H2 mixing ratio is orders of magnitude below
   the maximum, at one implicit step in a pass, so their own H2 balance never
   relaxes and the drift measure, absolute in the grid maximum, does not see
   them. The second operation needs its own step, a cell-relative drift measure
   or a wind-window relaxation after the bounded pass, before the H2 wind row
   can be judged at 1e-5 at all. Recorded and NOT adopted: it changes the
   acceptance of the carrier route and is the user's to decide.
