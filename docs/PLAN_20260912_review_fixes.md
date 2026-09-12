# Plan of 2026-09-12: the seven findings of `Update_EXHALE_stage2_review_20260912.md`

User instruction (2026-09-12): analyze the review against the code and, where it
is right, plan and start the corrections. Advisor's verification: all seven
findings are confirmed in the source (R1 `diffusive_photochemistry.f90:1847`
sets the history flag on any uncovered interval outside physical mode, so a
rejected stationary trial marks it; R2 the carrier branch of
`steady_wind_with_element_diffusion` keeps `u(3,:)` while the relaxation held
`p`; R3 `equilibrate_chemistry_at_fixed_pressure` hands `ioniz_eq` the T of
the entry particle count; R4 `element_diffusion_step` judges only when
`status` is present and the marching call at `EXHALE_main.f90:2297` passes
none; R5 `mass_row_cell_verdict` clips the tolerance at 1 and tests `q/tol < 1`
where the comment beside `cert_tol_mass_ceiling` says such a row cannot be
judged; R6 `build_molecular_hydrogen_column` is not mass-closed; R7 the outer
summary divides `row_max` by the binding cell's tolerance).

Rules: `docs/worker_rules.md`. Items Q1 to Q4. Scratch `<scratchpad>/Q<n>/`.

## The state contract decided for R2 and R3 (advisor)

At a fixed hydrodynamic state the quantities the solve holds are the conserved
variables `u = (rho, rho v, E)`. A composition update at fixed hydro therefore
keeps `u` and recomputes the pressure and the temperature from the unchanged
thermal energy `rho e = E - (rho v)^2/(2 rho)` with the caloric EOS of the NEW
composition (`pressure_from_energy_density`, then `comp_T_from_p`). The
chemistry refresh iterates (densities, p from E, T, sweep) until T stops
moving, so the state handed on is a joint fixed point of transport, chemistry
and the EOS at fixed `u`; a held pressure is not a state the solve knows.

| item | findings | files | worker |
| --- | --- | --- | --- |
| Q1 | R1, R2, R3, R6, R7 and the outer progress metric of review section 10 | `diffusive_photochemistry.f90`, the carrier branch and the summary of `steady_wind_with_element_diffusion` in `EXHALE_main.f90`, `carrier_retry.f90` | one worker |
| Q2 | R4 in the operator | `binary_element_diffusion.f90`, `element_operator_tests.f90`, `diffusion_tests.f90` | one worker |
| Q3 | R5 | `certification.f90`, `src/tests/certification/` | one worker |
| Q4 | R4 in the marching caller (after Q1 releases `EXHALE_main.f90`) | `EXHALE_main.f90` marching loop, `attempted_step.f90` | after Q1 |

## Status at the close of 2026-09-12 (advisor)

| item | status | record |
| --- | --- | --- |
| Q1 | DONE (advisor; the Opus workers assigned to it stopped on the weekly limit and the advisor implemented it) | `docs/Update_EXHALE_stage2.md` section 9: the state contract, the R1 trial argument, the conserved-state relaxation with the closed chemistry refresh, the mass-closed carrier retry column with absolute rows (138/0; RED on the entry text 5 rows), the binding-cell attribution and the joint progress measure; the review's 500-cell probe re-measured on the new contract (0.0 / 0.0 / 6.9e-16 against 5.6e-9 / 3.3e-3 / 3.3e-3). |
| Q2 | DONE (Sonnet worker, verified by the advisor from the diff and the suites) | Unconditional judgment in `element_diffusion_step`, `element_step_last_status`, `element_step_outcome_text`; the Newton-floor defect of `solve_mass_fraction` found and fixed beside it (`newton_drop` 2e-5). `element_operator` 28/0, `diffusion_tests` 35/0, `steady_species_rows` 195/0. |
| Q3 | DONE (Sonnet worker, verified by the advisor from the diff and the suites) | RESOLVED / UNRESOLVED outcome of `mass_row_cell_verdict`, `mass_row_column_verdict`, the unresolved cells counted and named in the report; addendum in `docs/certification_tolerance_anchoring_20260910.md`. `certification` 84/0, `krylov_and_dogleg` 334/0. |
| Q4 | DONE (advisor) | The marching loop refuses the attempted step on a refused element step in phys mode (`as_reject_element`, retaken at half dt) and counts the refusals in init mode; test knob `EXHALE_ELEMENT_REFUSE_STEPS`; the reason census bound of `attempted_step_note_step` corrected. The stagnation row of `grid_and_gates` re-anchored to `EXHALE_JFNK_MAXIT=40`, `EXHALE_CARRIER_TRUST=1e-6` (4/4 on the merged build). |

Gates of the merged tree are in section 9 of `docs/Update_EXHALE_stage2.md`
("Gates"). The handoff of record is `docs/session_handoff_20260912.md`.
