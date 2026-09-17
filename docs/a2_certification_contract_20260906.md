# A2: one evaluator, four acceptance contexts (PLAN_20260906_rev2, Step A2)

Status: design contract written by the advisor before A2 is coded (review 2,
F1 and R5). It consumes `docs/b1a_active_equation_inventory_20260906.md`
(the equations, closures, constraints and validity states) and
`docs/a0_run_mode_contract_20260906.md` (the run states). Every statement
about the present code is READ from the source named.

## 1. What exists today, and why it is not enough

- The stationary gate `steady_gates_met` (`steady_residual.f90`)
  tests a hydrodynamic residual norm, the mass-flux spread, an optional
  carrier residual, and since A1 an optional count of cells without a
  chemical root. The marching caller (`EXHALE_main.f90`, the residual and
  gate block of the marching stop test) passes
  `.false., 0.0d0` for the carrier condition and cannot pass the chemical
  count because the marching loop never receives the sweep ledger (A1 class
  report). The steady solver passes the count at the three points that
  declare a state solved. (Read again 2026-09-07: the marching caller now
  takes the sweep ledger and passes `n_cells_without_chemical_root` of it;
  only the carrier condition is still `.false., 0.0d0` there.)
- `carrier_steady_residual` (`diffusive_photochemistry.f90`)
  zeroes its norms and returns when `thereis_mol`, `carrier_transport` or
  `bg_ready` is false, so an active model with no background reads as a
  converged zero; despite `intent(in)` arguments it refreshes base state,
  photolysis, advection correction, diffusivities, `row_terms`, the scales
  and (under `weno_mode .eq. 1`) the H2 headroom. Its `rvol` is a ratio of
  sums over all carriers of a region and can hide an unbalanced species
  (review R7).
- Since A1, `carrier_returned_state_verdict` decides a carrier update on the
  returned state, row by row, with the same row scale `row_terms`
  (`|res|/max(terms, 1e-300)` against `carrier_accept_tol = 1e-8`); the
  steady residual uses the same expression with `dt = 1e30`. One scale, two
  equations (time-discrete and stationary).
- No state-consistent evaluator exists for the elemental transport balances,
  the eliminated-species closures (evaluated only inside the sweep that
  rewrites `f_sp` and `rho`), the H(n=2) system, or the marching energy
  update (two fixed Newton iterations, no residual test). The hydro rows have
  one, but it needs `heat`/`cool` from a preceding mutating `ioniz_eq`
  (B1a section 2).
- Validity states: the CO destruction model has a cumulative domain record
  (the thermodynamic CO ceiling that this line first described was deleted on
  2026-09-06, item CEILING-DEL); the energy floor has counters; the
  composition projection, the `rho` rewrite and the remaining out-of-domain
  closures are recorded nowhere (B1a section 4).

## 2. The evaluator

One routine, `certification_evaluate(state, field, context, report)`, with:

- **Inputs**: the state to be judged (conserved variables, species
  fractions, carrier fractions, level populations, temperature), the
  radiation field and columns OF THAT STATE (recomputed inside an isolated
  workspace, never taken from a cache that a later operation may have
  refreshed), the active-equation inventory (B1a, derived from the
  configuration flags at startup), and the context (section 3).
- **No side effects**: the routine works on copies; module-level saved
  arrays that today are refreshed by the residual routines (`row_terms`,
  scales, headroom, backgrounds, photolysis arrays, advection corrections,
  `bg_cell_adopted`, streak counters) are either not touched or restored
  before return. Test: repeated evaluation before and after an optional
  output call leaves every module array bit-identical.
- **Never re-equilibrates a transported partition** to measure it: the
  eliminated-species closures are evaluated as residuals of the closure
  equations at the given composition, not by re-solving and comparing.
- **Output, per equation of the inventory**: an evaluation status
  (`not_applicable`, `evaluated`, `unavailable`), the row values, the row
  scale with its units and floor, and finiteness. An active equation whose
  inputs are not available (`bg_ready` false, missing columns) is
  `unavailable`, never zero. Plus the constraint residuals (element totals,
  charge, mass), the energy identity of the last accepted step where a
  time-discrete context asks for it, and the validity state (B6) with the
  five categories of B1a section 4.

## 3. The four contexts

| Context | Caller | Required from the report | Not required | Decision |
|---|---|---|---|---|
| **Probe** (residual / Jacobian evaluation for a Newton direction or a finite difference) | `eval_residual` and its Jacobian-action callers in `steady_newton.f90` | every requested row `evaluated` and finite; domain checks of the closures (a class-4 chemistry is `unavailable`, not a value); no persistent mutation | any convergence | usable / not usable |
| **Stationary Newton trial** | the line search and trust-region acceptance of `steady_newton.f90` (`trust_region_step` and the line search of `solve_steady_jfnk`) | physical bounds (positivity, simplex), valid local closures (`n_cells_without_chemical_root = 0` on the trial), then the solver's own merit-decrease or trust-region rule on the residual it minimizes | a zero global residual; any stationary tolerance | admissible / rejected trial |
| **Physical time step** (A0 physical mode) | the adoption boundary of the marching step (B3a controller; before it, the point where the step is accepted today) | the time-discrete balances of every active equation within their tolerances (the carrier verdict of A1 for carriers; the analogous returned-state conditions for the energy update after B2 and for the composition step after B3c), admissibility, element and charge invariants to round-off, the energy identity of the step, the integration-error requirement (B3) | a zero stationary residual | accepted step / rejected step (restore) |
| **Stationary certification** | the point where a state is declared converged and written: `steady_gates_met` callers in `EXHALE_main.f90` and `steady_newton.f90` | every active stationary equation `evaluated` and within its tolerance (hydro rows, each transported carrier balance as transport minus reaction, each elemental transport balance, the level equations), every constraint, `n_cells_without_chemical_root = 0`, validity state clear (or the run excluded and the state written as uncertified), exact output-state consistency (B3c) | nothing may be skipped or replaced by the last time-discrete residual | certified / written uncertified with the failing rows named |

Class-4 chemistry: `unavailable` in the probe, rejects a trial, rejects a
physical step, refuses certification. It is admissible only in A0's
initialization mode, where the evaluator still runs and reports, and its
report is labeled initialization.

## 4. Norms and tolerances

- One condition per active independent balance; the combined region ratio
  `rvol` is reported, never decisive.
- Row scale: the sum of the row's own term magnitudes plus an absolute floor
  in the row's units (the carrier floor of A1: 1e-20 of the free element
  density in rate units; the hydro rows: the existing `scale` of
  `steady_residual`; elemental transport: to be defined with its evaluator,
  same construction). Units and floors stated in the evaluator's header and
  printed with every report.
- Local safeguard: the maximum over cells of the row measure, not a volume
  average, so a front or a thin layer is not averaged away; the volume
  average is reported beside it.
- Tolerances per context are parameters of the evaluator with their own
  names (`cert_tol_hydro`, `cert_tol_carrier`, `cert_tol_element`,
  `step_tol_*`), set in the implementation brief from measured convergence
  studies, never to make an existing snapshot pass.

## 5. What replaces what

| Today | After A2 |
|---|---|
| `steady_gates_met(resid_max, u, resid_th, fspread, .false., 0.0d0)` in the marching stop test of `EXHALE_main.f90` | the certification context on the exact state about to be written, with the sweep ledger of that state (the marching loop asks `ioniz_eq` for it, or the evaluator recomputes the closure residuals) |
| the trial count of cells without a chemical root in `steady_newton.f90` | unchanged as a local-closure validity check in the trial context (A1 class); no stationary condition added |
| `carrier_steady_residual` called by the steady driver with its side effects | the same equation evaluated through the evaluator's isolated workspace; `rvol` reported only |
| no elemental transport residual | a stationary residual of the He/H partition and the trace-element balances, evaluated on the state (new routine in `binary_element_diffusion.f90`, state-consistent, side-effect-free) |
| no H(n=2) or He 2^3S residual | the level balance residuals at the state's field |
| no record of the composition projection, the `rho` rewrite, out-of-domain closures | validity-state producers (B6): counters and a header field; until B3c removes the projection, its activity is a recorded unbudgeted correction and refuses certification of a physical step |

## 6. Tests (from the plan's matrix, A2 rows)

| Test | Required |
|---|---|
| active carriers with `bg_ready` false | status `unavailable`; certification refused; never a zero |
| inactive model | `not_applicable`, no row |
| stale background after a state change | the evaluator's own field, not the cache; detected by a constructed mismatch |
| repeated evaluation around output calls | every module array bit-identical before and after |
| small-rate species with a large fractional imbalance beside a balanced dominant one | caught by the row measure, invisible to `rvol` |
| oxygen carriers, H+ transport, active element diffusion | all evaluated, none silently absent |
| positive, merit-decreasing Newton trial with a large global residual | admissible trial, not certified |
| valid transient with a nonzero stationary balance | accepted physical step, no stationarity claim |
| trial with invalid eliminated-species chemistry | `unavailable` closure, trial rejected |
| a class-4 cell at the final state | written uncertified, cell named |

## 7. Order of implementation

1. The evaluator skeleton with the inventory, statuses, the isolated
   workspace and the hydro and carrier rows (both already have residual
   code; the carrier one is wrapped so its side effects stay inside the
   workspace).
2. The certification context at the two declaration points, replacing the
   `.false., 0.0d0` call and consuming the class count (closes the A1 open
   item).
3. The probe and trial contexts wired into `steady_newton.f90` without
   changing its merit rule.
4. The elemental transport and level residuals (new evaluators).
5. The validity-state producers that exist today (CO record, floors) read
   into the report; the missing producers land with B2, B3c, B4.
6. The physical-step context lands with B3a (the adoption boundary).

## 8. Decided by the user (2026-09-06)

- A run whose final state is written uncertified exits with a NONZERO status
  (the state and its headers are still written; the summary names the
  failing entries), so scripts can tell.
- The tolerance values: proposed to be fixed by a convergence study on
  `wasp_full_newton` and `mol_carrier` in the implementation brief, and
  reported with the study.

## 9. Tolerance study

Section 4 says the tolerances of the evaluator are set from measured
convergence studies and never to make an existing snapshot pass. This section
is that study, and it is a measurement only: nothing in the production source
was changed for it. Every run below used the current binary on a scratch copy
of a regression case with `OMP_NUM_THREADS=1`, one control parameter varied at
a time.

Every number in a table is MEASURED (a run made for this study). Every number
attributed to a source file is READ. `max` is the quantity the evaluator
decides on, the maximum over the physical cells of `|res|/max(scale, 1e-300)`;
`vol` is the volume-weighted companion the report prints beside it. An
"observed order" p is defined by `measure` proportional to `control^p` over the
whole ladder, and is stated so that a p near zero can be read for what it is:
the control and the measure are not related.

### 9.1 The marching relaxation cannot reach the hydrodynamic rows

`wasp_full` (WASP-121b, He 2^3S and metals), stopped by its own `du`
threshold, with the WENO3 entry of `du_th [PLM,WENO3]` tightened by a decade
in two steps. `du` is the fractional radial spread of `rho v r^2`.

| `du_th` (WENO3) | steps | `du` at the stop | mass max (vol) | momentum max (vol) | energy max (vol) | log10 Mdot |
|---|---|---|---|---|---|---|
| 1.0e-3 | 16230 | 9.995e-4 | 1.326e-02 (2.051e-03) | 1.785 (7.367e-01) | 7.245e-01 (2.825e-02) | 13.29 |
| 3.0e-4 | 22645 | 2.991e-4 | 1.234e-02 (1.823e-03) | 1.860 (6.095e-01) | 4.014e-01 (8.440e-03) | 13.30 |
| 1.0e-4 | 27962 | 9.921e-5 | 8.356e-03 (2.156e-03) | 1.932 (5.608e-01) | 5.920e-01 (1.121e-02) | 13.31 |

Observed order over the decade of `du`: mass p = 0.20, momentum p = -0.034
(the row RISES), energy p = 0.087 and not monotone. The worst cell stays where
it was: the mass row at cells 11-16 (the base), the momentum row at cells
495-496 (the outer edge), which is exactly where the `du` functional, measured
over the escape window, has nothing to say.

**The `du` stop and the hydrodynamic row measures are not measures of the same
thing.** A decade of `du` buys a factor 1.6 on one row, nothing on the second
and nothing reliable on the third, while the mass-loss rate moves by 0.02 dex.
No `du` threshold fixes a tolerance, because no `du` threshold moves the rows.

The same ladder in the step budget, on `mol_carrier` (the hot-Uranus molecular
gate with H2 carrier transport), which stops at a step cap rather than at `du`:

| steps | mass max (vol) | momentum max (vol) | energy max (vol) | H2 carrier max (vol) |
|---|---|---|---|---|
| 12000 | 1.807e-01 (1.393e-02) | 1.633 (9.742e-01) | 1.851 (7.783e-01) | 2.957e-01 (3.301e-02) |
| 24000 | 2.299e-02 (2.187e-03) | 8.663e-01 (3.416e-02) | 6.792e-01 (6.709e-02) | 2.841e-01 (5.257e-03) |
| 48000 | 2.258e-02 (1.801e-03) | 5.402e-01 (7.182e-03) | 6.628e-01 (5.076e-02) | 2.637e-01 (5.195e-03) |

Observed order over the factor 4 in steps: mass p = 1.50, momentum p = 0.80,
energy p = 0.74, **H2 carrier p = 0.082**. The hydrodynamic rows fall between
the first two rungs and then stop; the carrier row is flat to 12 per cent over
a four-fold step budget, seven decades above `cert_tol_carrier`. A relaxation
of this configuration does not approach a stationary carrier balance, and
`mol_carrier` run to a `Solver: Newton` hand-off does not reach it either: in
83000 steps `du` fell only to 6.3e-2 against the 1.0e-2 hand-off threshold, so
the JFNK never started.

### 9.2 The stationary solve moves the rows, and lands at a fixed multiple of its own target

`wasp_full_newton`, the same model with `Solver: Newton`, with `Resid tol`
added and tightened by two decades. The run stops in the same place every time
(JFNK finish at step 4659), so the ladder isolates the solve.

| `Resid tol` | `||R||` the JFNK printed | mass max (vol) | momentum max (vol) | energy max (vol) | energy / target | log10 Mdot |
|---|---|---|---|---|---|---|
| 1.0e-5 (the default) | 4.088e-04 | 1.729e-11 (3.049e-12) | 1.343e-07 (5.082e-09) | 4.088e-04 (1.233e-04) | 40.9 | 13.26 |
| 1.0e-6 | 4.486e-05 | 9.604e-13 (1.753e-13) | 2.140e-09 (1.831e-10) | 4.486e-05 (1.419e-05) | 44.9 | 13.26 |
| 1.0e-7 | 4.700e-06 | 2.630e-13 (3.968e-14) | 1.280e-09 (5.689e-11) | 4.700e-06 (1.378e-06) | 47.0 | 13.26 |

Observed order over the two decades: energy p = 0.97, mass p = 0.91, momentum
p = 1.01 for the first decade and 0.23 for the second (it is flattening). The
energy row is the binding one at every rung and equals the `||R||` the solve
prints, as it must, since `||R||` is the maximum over the three rows.

Two things follow, and they are the reason a tolerance cannot simply be tied
to `Resid tol`.

- **The energy row follows the tolerance, one for one.** Tightening the target
  by a decade tightens the achieved row by a factor 9.1 and then 9.5. So the
  solve is not at a floor at 1e-7 and the row is limited by the target, not by
  the arithmetic (section 9.4 measures how far the floor is: five decades
  further down).
- **The achieved row is 41 to 47 times the target at every rung.** That ratio
  is not noise, it is structural, and section 5 of the A2 report names the
  cause: `steady_newton.f90` sets `info = 1` once at the start of the solve,
  the loop-top gate sets `info = 0` and leaves the loop, and the post-loop
  re-evaluation at the state's own composition -- the honest number -- can
  only confirm `info = 0` and never take it back. So a solve declared
  converged returns a state whose measured energy row is about 45 times the
  number it was asked for. **`cert_tol_hydro = resid_th` is therefore a
  condition no JFNK finish can meet, for a reason that has nothing to do with
  the state.**

### 9.3 The diffusive runs: the marching route cannot measure this operator

`mol_diffusion` (He/H binary diffusion on the hot-Uranus gate) and
`lower_profile` (a profile handoff with trace metals and `K_zz(p)`), both
stopped by a step cap.

| case | steps | He/H partition max (vol) | trace elements max (vol) |
|---|---|---|---|
| `mol_diffusion` | 12000 | 6.422e-05 (5.583e-06) | not carried |
| `mol_diffusion` | 24000 | 1.091e-04 (7.893e-06) | not carried |
| `mol_diffusion` | 48000 | 7.438e-05 (4.707e-06) | not carried |
| `lower_profile` | 12000 | 4.810e-03 (1.657e-04) | 1.234e-02 (3.247e-03) |
| `lower_profile` | 24000 | 4.402e-04 (6.673e-05) | 3.462e-03 (3.216e-04) |
| `lower_profile` | 48000 | 6.589e-04 (8.523e-05) | 4.434e-03 (6.581e-04) |
| `lower_profile`, `Solver: Newton` | 60000 | 6.961e-04 (9.292e-05) | 4.394e-03 (8.755e-04) |

Neither falls. `mol_diffusion` has observed order p = -0.11 over the factor 4
in steps; `lower_profile` falls by a decade between the first two rungs and
then comes back up and settles at 5e-4 to 7e-4 and 3e-3 to 4e-3 through 60000
steps. **These numbers are a transient, not the operator's floor, and they
cannot fix a tolerance.**

The reason is READ, and it is decisive. On the marching route the composition
is moved by `element_diffusion_step(rho, v, T, f_sp, dt_loc)`
(the attempted step of `EXHALE_main.f90`): one backward-Euler step at the
hydrodynamic CFL step. The header of `relax_element_composition`
(`binary_element_diffusion.f90`) states what that costs: at the base
of HD 209458 b the cell diffusion time is about 1e6 s against a CFL step of
about 1 s, so a relaxation driven at the CFL step moves the composition by
about 1e-6 of the way per step and no practical number of steps reaches the
steady state. The measured plateaus are that statement.

The route that does reach it is the steady one, where
`relax_element_composition` steps on the composition time scale and grows the
step geometrically to 1e12 times its start, so its late passes are direct
steady solves. **That route was not reached by any diffusive configuration in
this study**: `mol_diffusion` with `Solver: Newton` ran 60000 steps and
`lower_profile` with `Solver: Newton` 60000 steps without `du` ever crossing
the 1.0e-2 JFNK hand-off. Five further runs -- those two continued to 90000
steps, `mol_carrier` on the steady route, and restarts of all three cases from
their own 48000-step states -- were still marching at `du` between 5.4e-2 and
6.1e-2 when this study closed, falling by about 2e-3 per 4000 steps. Nothing
came within reach of the hand-off.

So the reachable value of the elemental transport rows is READ from the
operator's own controls rather than measured here, and the operator states it
on exactly the measure the certification uses -- `composition_residual` returns
the residual "measured RELATIVE to the size of the terms of its own row", the
same construction as `row_terms`:

- `newton_tol = 1.0d-12`, the target the inner Newton is asked for;
- `newton_floor = 1.0d-8`, the residual below which a pass that no longer
  halves is taken as converged, declared to be "the arithmetic floor";
- the reachable floor stated with its range: "1e-12 in a wind, ~1e-5 in the
  K_zz = 2e12 homopause column of test T7, where the time term is negligible
  against the eddy term";
- `he_relax_tol = 1.0d-12` for the outer relaxation, on `max|dX|` per step
  over `X_base`.

### 9.4 The round-off floor of the row measures

The row measure is formed from a state that is itself known only to the
arithmetic of the machine, so it has a floor. Measured on the certified-looking
state of section 9.2, `wasp_full_newton` at `Resid tol` 1.0e-7: the residual
was evaluated at the final state, then at the same state with the temperature
of every cell moved by ONE ULP at fixed density and velocity (one ulp of the
cell's internal energy added to its conserved energy), then at the first state
again. The measurement instrument was a scratch copy of the tree outside the
repository; no production file was touched.

| row | row measure at the state | repeat control | one ulp in T | worst cell |
|---|---|---|---|---|
| mass | 2.667e-13 | 0.000 | 1.018e-13 | 77 |
| momentum | 2.657e-09 | 0.000 | 5.127e-12 | 492 |
| energy | 4.852e-06 | 0.000 | 1.643e-11 | 4 |

The repeat control is exactly zero in every row: the residual is bitwise
reproducible, so the floor measured here is the state's own ulp sensitivity
and nothing else. (The row measures in the first column differ from section
9.2's by a few per cent because this instrument re-sweeps the composition
through `newton_residual` while the certification uses the residual the run
assembled; the orders are the same.)

What it says, row by row:

- **the mass row is already an arithmetic statement.** At the tightest rung it
  reads 2.6 times its own one-ulp floor. There is no room below 1e-12 in this
  row and no purpose in asking for it;
- **the momentum row has about two decades left** (2.7e-09 against a floor of
  5.1e-12), which is consistent with its ladder flattening between the 1e-6
  and 1e-7 rungs;
- **the energy row is five decades above its floor** (4.9e-06 against
  1.6e-11), so what stops it is the solver, not the machine. Every decade the
  ladder bought was real and more are available.

The rows that carry their own scale were not perturbed this way, because the
case that certifies carries no carrier and no transport. Their floors are READ:
`carrier_row_roundoff = 64*epsilon = 1.42e-14` of the row's full terms, with
the note that the carrier Newton is MEASURED (in
`src/tests/carrier_retry/`) to reach 2.6e-16, about one eps. The closure and
level rows are dimensionless against their own turnover, and the values this
study measured for them (below) sit at the same place.

### 9.5 Recommended tolerances

Each recommendation is anchored on what a converged solve was measured to
reach, with a factor of about ten of margin, and is checked to stand at least
a decade above the arithmetic floor of its row. None was chosen to make a
present snapshot pass; section 9.6 lists what does and does not pass, and the
answer is that most of the matrix does not.

| entry | today | recommended | anchored on | margin over the floor |
|---|---|---|---|---|
| `cert_tol_hydro` (mass) | `resid_th` | **3e-12** | 2.630e-13 reached at `Resid tol` 1e-7 | 29x (floor 1.018e-13) |
| `cert_tol_hydro` (momentum) | `resid_th` | **1e-8** | 1.280e-09 reached at the same rung | 2000x (floor 5.127e-12) |
| `cert_tol_hydro` (energy) | `resid_th` | **5e-5** | 4.700e-06 reached at the same rung | 3e6x (floor 1.643e-11) |
| `cert_tol_carrier` | 1e-8 | **1e-8, kept** | `newton_floor` of the carrier Newton on the same measure (READ) | 7e5x (floor 1.42e-14, READ) |
| `cert_tol_element` | 1e-6 | **1e-8** | `newton_floor` of the composition Newton on the same measure (READ) | see the caveat below |
| `cert_tol_closure` | 1e-6 = `ieq_res_tol` | **1e-6, kept** | the sweep's own root tolerance | 6 to 12 decades measured |
| `cert_tol_level` | 1e-6 = `ieq_res_tol` | **1e-6, kept** | the same | 8 to 15 decades measured |

**One number per hydrodynamic row, not one for the three.** The three rows
reach values seven decades apart on the same converged state (2.6e-13,
1.3e-09, 4.7e-06) and their floors are two decades apart. A single
`cert_tol_hydro` either certifies a mass row that is eight decades above its
own floor, or refuses an energy row that no solve in this code has ever
reached. Splitting it is the change this study asks for, and the names should
say which row refused: `cert_tol_mass`, `cert_tol_momentum`,
`cert_tol_energy`.

**`cert_tol_hydro` must stop inheriting `resid_th`.** Section 9.2 measures a
structural factor of 41 to 47 between a JFNK target and the row the solve
hands back, so the inherited value is a condition the solver cannot meet
whatever the state is. The fixed numbers above replace it. The right repair is
in the solver -- `info = 0` is decided before the state's own residual is
re-measured -- and when that is made, this study should be repeated, because
the anchor values will move down.

**`cert_tol_carrier` is kept at 1e-8 and is UNCONFIRMED.** No configuration in
the matrix produces a stationary carrier state: `mol_carrier` is a step-capped
relaxation whose H2 row is flat at 0.26 to 0.30 over a four-fold step budget,
and the same model on the steady route did not reach its JFNK hand-off in
83000 steps. So this study cannot say what a converged carrier balance
reaches. What it can say is that 1e-8 is not arbitrary: it is the operator's
own `newton_floor` on the operator's own measure, it stands six decades above
the assembly's round-off bound `64*eps`, and the carrier Newton is recorded as
reaching one eps on a molecular column. The measurement that would confirm it
is a JFNK finish with `Molecular carrier transport: True`, and it is not
available today.

**`cert_tol_element` should be 1e-8, not 1e-6, and it too is unconfirmed by a
run.** 1e-6 was inherited from the carriers with no measurement behind it;
1e-8 is the number the composition operator itself declares to be what its
discretization resolves, on a residual it forms relative to its own row terms,
which is the identical construction. The caveat is stated at the same place in
the operator: the reachable floor of that residual is 1e-12 in a wind but
about 1e-5 in a K_zz-dominated homopause column where the eddy term dominates
and the time term is negligible. A cell in that regime cannot meet 1e-8, and
the certification should say so -- report the regime rather than refuse
silently -- rather than have the tolerance loosened by two decades everywhere
to accommodate it. That is an open item for whoever wires the elemental
transport evaluator's report.

**`cert_tol_closure` and `cert_tol_level` are kept, with the margin measured.**
Across the fifteen final states of the ladders the eliminated-species closure reads
8.182e-18 to 1.561e-12 and the level rows 5.423e-21 to 7.482e-14 (He 2^3S)
and 1.973e-16 to 2.022e-15 (H(n=2)) -- six to fifteen decades below 1e-6 --
and both fall as a run relaxes (the `wasp_full` closure reads 1.561e-12 at the
`du` 1e-3 stop and 2.538e-14 at the 3e-4 stop). No rung of any ladder brings
either row within two decades of its tolerance. There is nothing here to
change.

### 9.6 What certifies at the recommended values

MEASURED, on the final state of every run of this study, against the
recommended numbers of section 9.5:

| run | mass vs 3e-12 | momentum vs 1e-8 | energy vs 5e-5 | other rows | certified? |
|---|---|---|---|---|---|
| `wasp_full_newton` + `Resid tol: 1.0e-7` | 2.630e-13 within | 1.280e-09 within | 4.700e-06 within | closure and both levels within by 9 to 15 decades | **yes** |
| `wasp_full_newton` + `Resid tol: 1.0e-6` | 9.604e-13 within | 2.140e-09 within | 4.486e-05 within (by 11 per cent) | the same | **yes, marginally** |
| `wasp_full_newton` as shipped (target 1e-5) | 1.729e-11 ABOVE | 1.343e-07 ABOVE | 4.088e-04 ABOVE | the same | no |
| `wasp_full` at any of the three `du` stops | 8e-03 to 1.3e-02 ABOVE | 1.8 to 1.9 ABOVE | 4e-01 to 7e-01 ABOVE | the same | no |
| `mol_carrier` at 12000, 24000, 48000 | ABOVE | ABOVE | ABOVE | H2 carrier 2.6e-01 to 3.0e-01 ABOVE 1e-8 | no |
| `mol_diffusion` at 12000, 24000, 48000 | ABOVE | ABOVE | ABOVE | He/H 6e-05 to 1.1e-04 ABOVE 1e-8 | no |
| `lower_profile` at 12000, 24000, 48000 | ABOVE | ABOVE | ABOVE | He/H and trace ABOVE 1e-8 | no |

**No case in the regression matrix as shipped certifies at these values, and
one added input line makes one of them certify.** That is the outcome the
study was for: the tolerances are not what stands between the code and a
certified state on the marching cases (the rows there are of order one, seven
decades away, and no relaxation control moves them), while on the stationary
route the gap is one to two decades of `Resid tol` plus the solver's own
acceptance defect.

Two recommendations follow that are not tolerances:

1. `wasp_full_newton` should state a `Resid tol` rather than take the default
   1e-5, because the default is one to two decades short of what the solve
   reaches for a cost of nothing (the run stops at step 4659 in all three
   cases, and log10 Mdot is 13.26 at every rung);
2. the JFNK's `info` handling should be repaired before this study is used to
   tighten anything further, since the 41-to-47 factor it introduces sits
   between every target and every achieved row.

## 10. The tolerance study repeated on the self-consistent residual

Section 9 closed with two statements that have since been acted on: the JFNK
returned a state whose rows were 41 to 47 times the target it was given, and
the anchors of section 9.5 would move once that was repaired. B5a made the
steady residual self-consistent -- the composition is eliminated at the state
being evaluated, at every residual evaluation and not only at the hand-back --
and B3a is applying section 9.5's provisional split, `cert_tol_mass` 3e-12,
`cert_tol_momentum` 1e-8, `cert_tol_energy` 5e-5. This section re-measures
the study on that code.

It is a measurement only: no production source was changed for it. Every run
used one privately built binary of the current tree, `OMP_NUM_THREADS=1`, on
scratch copies of the regression cases outside the repository. MEASURED means
a run made for this section; READ means taken from a source file or a run's
own printed record.

**The finding, first. The 41-to-47 factor is gone: the row the solve hands
back is now 0.36 to 0.93 times the target it was given, at four rungs spanning
three decades, and the state handed back scores exactly what its accepted
iterate scored (ratio 1.00 at every rung, printed by the solver).** What
remains is a simple arithmetic relation and it settles the configuration
question: since the maximum row IS the returned `||R||` and the binding row is
now momentum, a run cannot certify unless its `Resid tol` is at least as tight
as `cert_tol_momentum`. `wasp_full_newton` at `Resid tol: 1.0e-8` is the first
state in this study to be certified by the solver itself, `info = 0`, every
row within tolerance.

### 10.1 The `Resid tol` ladder on `wasp_full_newton`

The same configuration as section 9.2 with one added input line, on the
current binary. The marching phase is identical at every rung (JFNK hand-off
at step 4656; section 9.2 measured step 4659 on the older code), so the ladder
isolates the solve. The three "certification rows" are the rows of the
certification of the state as written.

| `Resid tol` | outer its | wall clock | `\|\|R\|\|` handed back | mass max (vol) | momentum max (vol) | energy max (vol) | max row / target | `info` | certified |
|---|---|---|---|---|---|---|---|---|---|
| 1.0e-5 (the default) | 11 | 590.0 s | 7.390e-06 | 1.102e-09 (1.131e-10) | 7.390e-06 (9.787e-07) | 8.375e-08 (8.015e-09) | 0.74 | 2 | no |
| 1.0e-6 | 12 | 608.3 s | 9.330e-07 | 1.392e-10 (1.429e-11) | 9.330e-07 (1.236e-07) | 1.057e-08 (1.012e-09) | 0.93 | 2 | no |
| 1.0e-7 | 13 | 627.8 s | 7.281e-08 | 1.086e-11 (1.115e-12) | 7.281e-08 (9.646e-09) | 1.398e-09 (4.429e-10) | 0.73 | 2 | no |
| 1.0e-8 | 14 | 658.1 s | 3.574e-09 | 5.333e-13 (5.467e-14) | 3.574e-09 (4.730e-10) | 6.275e-10 (1.634e-10) | 0.36 | **0** | **yes** |

`log10 Mdot = 13.27` at all four rungs. Observed order over the three decades:
mass p = 1.10, momentum p = 1.10, energy p = 0.71. One outer iteration buys one
decade, exactly, at every rung; the iteration lines agree digit for digit
across the rungs up to the point where the looser target stops.

Four things follow, and each replaces a statement of section 9.2.

- **The relation between the target and the returned rows is now 0.36 to 0.93,
  not 41 to 47.** The last Newton step overshoots its target, as it should. The
  solver's own printed check, "residual of the state handed back at its own
  composition ... ratio to the accepted iterate's own `||R||`", reads
  `1.00E+00` at all four rungs; section 9.2's 7.95 for the same quantity was
  the lag.
- **The binding row moved from energy to momentum.** At every rung the maximum
  row is the momentum row at cell 499, the outer edge, and it equals the
  returned `||R||` to all printed digits. The energy row, which section 9.2
  measured at 4.088e-04 at the 1e-5 rung, now reads 8.375e-08 at the same rung:
  a factor 4.9e3 lower, and it is the lag that was removed. The mass and
  momentum rows moved the other way -- 1.729e-11 to 1.102e-09 and 1.343e-07 to
  7.390e-06 at that rung -- because the solve now stops when the true maximum
  meets the target instead of being driven further by an energy row that was
  measuring the previous evaluation's composition.
- **The rows are limited by the target, not by the arithmetic.** The ladder is
  clean over three decades and section 10.4 measures the floors: the mass row
  at the certified rung is 10 times its own one-ulp floor, the energy row 40
  times, the momentum row three decades above.
- **`info` is now decided by the certification.** The refusal line names the
  row: at the three loose rungs, "hydrodynamic momentum row measure ... above
  1.0E-08 at cell 499". `cert_tol_hydro = resid_th` is gone with the split, and
  what stands between a solve and `info = 0` is whether its own `Resid tol` is
  tight enough to reach the fixed certification tolerance.

Cost: 68 s, 11.5 per cent of the whole run, buys the three decades from the
shipped default to the certified rung. The marching phase (4656 steps) is the
bulk of every run and is common to all four; the JFNK adds about 25 s per
decade. The composition elimination costs 14.05 to 14.10 sweeps per residual
evaluation and the Newton model carries 15 Picard terms, at every rung.

### 10.2 The PTC route

MEASURED. `EXHALE_PTC=1` on the `wasp_full_newton` configuration (the
environment variable is READ from `EXHALE_main.f90`, which also takes
`EXHALE_PTC_DTAU0`, `EXHALE_PTC_NFIX` and `EXHALE_PTC_JFNK`; the route is
refused outright when the ionization state is transported).

**It runs and it does not certify.** The route is entered before the time
loop, takes its pseudo-time step from the CFL step (`PTC dtau0 = 1.01E-04`)
and solves from the default isothermal initial condition. It starts at
`||R|| = 1.952` and stalls at 1.958: the line-search parameter collapses from
iteration 6 onward (`lam = 9.54E-07` at almost every iteration thereafter, with
`||F||2` frozen at 3.55e+02), and 41 iterations produced no digit of progress.
The run was stopped at a 2700 s wall-clock cap without reaching its 3000
iterations, so it has no certification record at all.

Cost, for the comparison the item asks for: about 66 s per PTC iteration
against about 5 s per JFNK iteration on the marched state. The two are not
solving the same problem -- the PTC route starts from a cold isothermal column
where the JFNK starts from a state marched to `du < 1e-2` -- but the arithmetic
per iteration is also heavier, because the PTC route builds a banded Jacobian
whose every colour is now a 15-pass composition elimination.

`ptc_warm`, the regression case that exercises this route from a warm restart,
could not be run (section 10.3).

### 10.3 The other `Solver: Newton` cases

Scratch copies, current binary, `OMP_NUM_THREADS=1`, each case's own
`maxsteps` exported where it has one.

| case | outcome | its | wall | `\|\|R\|\|` / target | mass | momentum | energy | `info` | certified |
|---|---|---|---|---|---|---|---|---|---|---|
| `wasp_he23off_newton` | JFNK finish, step 4656 | 12 | 539.6 s | 1.988e-06 / 1e-5 | 2.927e-10 | 1.988e-06 | 2.292e-08 | 2 | no |
| `newton_rsw01` (du switch 0.1) | JFNK finish, step 3684 | 13 | 559.0 s | 1.009e-06 / 1e-5 | 1.506e-10 | 1.009e-06 | 1.144e-08 | 2 | no |
| `newton_rsw05` (du switch 0.5) | JFNK finish | 9 | 421.9 s | 2.110e-06 / 1e-5 | 3.149e-10 | 2.110e-06 | 1.212e-07 | 2 | no |
| `roundtrip` | 40-step cap, no JFNK | 0 | 1.8 s | -- | 1.180 | 1.588 | 1.979 | -- | no |
| `oxygen_chemistry` | 1000-step cap, no JFNK | 0 | 225.5 s | -- | 1.974 | 1.972 | 1.966 | -- | no |

`log10 Mdot = 13.27` for the three WASP-121b cases that reach the solve, the
same value the ladder gives. All three refuse on the same two rows and in the
same place -- mass and momentum, cell 499 -- with the momentum row 1.0e-06 to
2.1e-06 against a 1e-5 target, that is 0.10 to 0.21 of the target, the same
relation as section 10.1. None of the three states a `Resid tol`.

`oxygen_chemistry` is the only run in this study that reaches the carrier rows:
H2, OH, H2O and CO all read exactly 1.000 (`vol` 2.8e-01, 1.0, 1.0, 6.8e-03),
and the record adds "CARRIER HISTORY NOT CERTIFIABLE: an interval was left
uncovered and the march went on from the entry carriers". It is a step-capped
relaxation, so this is a transient and not the operator's floor;
`cert_tol_carrier` is still unconfirmed by any run, exactly as section 9.5 said.

**Four cases are BLOCKED and could not be started**: `jfnk_hd189`,
`jfnk_hd189_tight`, `ptc_warm` and `heh_1_newton` (renamed 2026-09-16 from `armD_D2_newton`; mapping in `docs/named_case_audit.md` section 6), all four `Load IC? True`.
`load_IC` refuses each with the same message: the stored
`output/Hydro_ioniz_IC.txt` holds cell centers of a different grid from the one
the run builds, worst relative difference 2.22e-02 to 2.25e-02 against a
1e-10 tolerance. That is the planetary radius decision of 2026-09-05 (`RJ =
7.1492e9` everywhere, +2.26 per cent on `R0`): every stored restart predates it
and every radius in it is short by that factor.

What regenerating them needs, stated so it can be scheduled: for each case, a
run of the SAME configuration with `Load IC? False` on the current grid, marched
to the state the stored `_IC` represented, whose `_IC` files then replace the
stored ones. For `jfnk_hd189` / `jfnk_hd189_tight` (HD 189733 b) and `ptc_warm`
(WASP-121b) that is a full marching relaxation each; `heh_1_newton` restarts
from the 12000-step state of `heh_1_12k` (renamed 2026-09-16 from `armD_D2`), so its regeneration is that case
re-run under its own step cap. The alternative the error message names is a
state-mapping step that interpolates a stored state onto the grid the run
builds, which does not exist today
(`docs/development_plan_20260905_rev3.md` section 10.2 item 3). Note also that
these four restarts predate the `# coupling:` header of Update_EXHALE_stage2.md
section 144, so even regenerated they will not reproduce their stored logs.

The remaining cases carrying `Solver: Newton` -- `mol_base_handoff`,
`mol_carrier`, `mol_diffusion`, `mol_ir_bands`, `mol_lyman_werner`,
`mol_metals`, `mol_sec_ion`, `lower_profile`, `heh_1_40k`/`heh_1_lw_40k`, `heh_0p3`/`heh_3`/`heh_10`/`heh_30`, `heh_1_12k`,
`heh_1_lw_12k` (renamed 2026-09-16 from `armD_D2_LW`) -- stop at their step cap long before `du` reaches the 1e-2
hand-off, so they never enter the steady route. Sections 9.1 and 9.3 measured
their rows (order one on the hydrodynamic rows, 6e-05 to 1.1e-04 on the He/H
partition), and B5a measured `mol_base_handoff` byte-identical across the
residual change because `eval_residual` is reached only from the steady solver.
They were not re-run here, and `roundtrip` and `oxygen_chemistry` above are two
members of that group that were, confirming the pattern.

### 10.4 The round-off floor on the certified state

Same instrument and same construction as section 9.4, on the state of the
1.0e-8 rung: the residual at the final state, then at the same state with the
temperature of every cell moved by ONE ULP at fixed density and velocity, then
at the first state again. The instrument is one added block in a scratch copy
of the tree outside the repository; no production file was touched, and the
copy is deleted.

| row | row measure at the state | repeat control | one ulp in T | worst cell | measure / floor |
|---|---|---|---|---|---|
| mass | 5.3327e-13 | 0.0000 | 5.1629e-14 | 108 | 10.3 |
| momentum | 3.5744e-09 | 0.0000 | 1.9683e-12 | 494 | 1816 |
| energy | 7.0220e-10 | 0.0000 | 1.7647e-11 | 36 | 39.8 |

The repeat control is exactly zero in every row again: the residual is bitwise
reproducible. The mass and momentum measures agree with the certification's to
every printed digit (5.333e-13, 3.574e-09), which they did not in section 9.4 --
the instrument and the certification now form the row from the same
self-consistent composition. The energy row differs by 12 per cent (7.022e-10
against 6.275e-10) and is at a different cell, which is what a row four decades
into its own noise looks like.

The floors themselves moved a little from section 9.4, all in the same
direction: mass 1.018e-13 to 5.163e-14, momentum 5.127e-12 to 1.968e-12,
energy 1.643e-11 to 1.765e-11.

**The state that certifies stands 10 times its own mass floor and 40 times its
own energy floor.** Those two rows have less than two decades of room left; the
momentum row, which is the one that binds, has three.

### 10.5 Re-anchored tolerances

| entry | value B3a applies | recommended | anchored on (MEASURED at the 1.0e-8 rung unless said) | margin of the tolerance over the floor |
|---|---|---|---|---|
| `cert_tol_mass` | 3e-12 | **3e-12, kept** | 5.333e-13 reached, 5.6x under the tolerance | 58x (floor 5.163e-14) |
| `cert_tol_momentum` | 1e-8 | **1e-8, kept** | 3.574e-09 reached, 2.8x under the tolerance | 5.1e3x (floor 1.968e-12) |
| `cert_tol_energy` | 5e-5 | **1e-6** | 6.275e-10 reached; 1.212e-07 is the largest energy row of any state in this study whose solve met its own target | 5.7e4x (floor 1.765e-11) |
| `cert_tol_carrier` | 1e-8 | **1e-8, kept, still UNCONFIRMED** | `newton_floor` of the carrier Newton (READ); no stationary carrier state exists in the matrix | 7e5x (floor 1.42e-14, READ) |
| `cert_tol_element` | 1e-8 | **1e-8, kept, still unconfirmed by a run** | `newton_floor` of the composition Newton on the identical measure (READ) | the homopause caveat of section 9.5 stands |
| `cert_tol_closure` | 1e-6 | **1e-6, kept** | 4.7e-14 to 9.0e-16 measured across this section's states | 12 decades |
| `cert_tol_level` | 1e-6 | **1e-6, kept** | 1.9e-21 to 1.0e-20 (He 2^3S), 1.0e-15 to 2.7e-15 (H(n=2)) | 9 to 15 decades |

**Only one number should move, and it should move DOWN.** `cert_tol_energy` =
5e-5 was anchored in section 9.5 on 4.700e-06, the energy row of the tightest
rung of the LAGGED ladder, and that number was the lag rather than the state:
the same configuration at the same target now reads 1.398e-09, three decades
lower, and at the certified rung 6.275e-10. Left at 5e-5 the energy row cannot
refuse anything -- it would take a state five decades worse than any measured
here -- so the entry would certify by construction rather than by measurement.
1e-6 is a decade above 1.212e-07, the largest energy row measured on any state
whose solve met its own target (`newton_rsw05`), and 5.7e4 above the row's
one-ulp floor. It refuses nothing in this study and it restores the entry's
ability to refuse.

The mass and momentum anchors, chosen this morning against the lagged ladder,
survive the re-measurement, but for a different reason than they were chosen
for: 3e-12 and 1e-8 are now the numbers a run must be given a tight enough
`Resid tol` to REACH, and at the 1.0e-8 rung the state clears both with a
factor of 5.6 and 2.8. Neither has a decade of margin, which is honest: the
momentum row is 1816 times its own arithmetic floor, so the room to tighten
`cert_tol_momentum` further exists, but it would have to be bought with a
`Resid tol` of 1e-9 and that rung was not measured.

**`wasp_full_newton` should state `Resid tol: 1.0e-8` as its configuration of
record.** It is the only setting in this study under which the case certifies,
it costs 68 s (11.5 per cent) and three Newton iterations over the shipped
default, and `log10 Mdot` is 13.27 at every rung of the ladder, so nothing
physical is bought or lost by the choice -- what is bought is a state the code
can certify. Section 9.6's recommendation of 1e-7 is superseded: on the
self-consistent residual 1e-7 leaves the mass row at 1.086e-11 and the momentum
row at 7.281e-08, both above their tolerances.

The same line is what `wasp_he23off_newton`, `newton_rsw01` and `newton_rsw05`
need, and by the same arithmetic (their momentum rows are 0.10 to 0.21 of a
1e-5 target, so a 1e-8 target should land them at 1e-9 to 2e-9), but that is an
inference from the relation of section 10.1 and was not measured on those three
cases.

### 10.6 What certifies now

MEASURED, on the final state of every run of this section, against the
recommended values of section 10.5:

| run | mass vs 3e-12 | momentum vs 1e-8 | energy vs 1e-6 | certified? |
|---|---|---|---|---|
| `wasp_full_newton` + `Resid tol: 1.0e-8` | 5.333e-13 within | 3.574e-09 within | 6.275e-10 within | **yes** |
| `wasp_full_newton` + `Resid tol: 1.0e-7` | 1.086e-11 ABOVE | 7.281e-08 ABOVE | 1.398e-09 within | no |
| `wasp_full_newton` + `Resid tol: 1.0e-6` | 1.392e-10 ABOVE | 9.330e-07 ABOVE | 1.057e-08 within | no |
| `wasp_full_newton` as shipped (target 1e-5) | 1.102e-09 ABOVE | 7.390e-06 ABOVE | 8.375e-08 within | no |
| `wasp_he23off_newton`, `newton_rsw01`, `newton_rsw05` (target 1e-5) | 1.5e-10 to 3.1e-10 ABOVE | 1.0e-06 to 2.1e-06 ABOVE | 1.1e-08 to 1.2e-07 within | no |
| `roundtrip`, `oxygen_chemistry` (step-capped, no JFNK) | order 1 ABOVE | order 1 ABOVE | order 1 ABOVE | no |
| `EXHALE_PTC=1` on the `wasp_full_newton` configuration | stalled at `\|\|R\|\| = 1.958`, no certification record | | | no |
| `jfnk_hd189`, `jfnk_hd189_tight`, `ptc_warm`, `heh_1_newton` | blocked: stored `_IC` on the pre-`RJ` grid | | | not started |

One case certifies, by one input line, and the two recommendations section 9.6
closed with are both now answered: the `info` handling is repaired (section
10.1), and the `Resid tol` the case should state is 1.0e-8 rather than 1.0e-7.

---

## 11. The species rows of the stationary system, and their anchors (B5, 2026-09-07)

Sections 9 and 10 could not anchor `cert_tol_carrier` or `cert_tol_element`
because no configuration of the matrix produced a stationary carrier state.
B5 made the transported balances rows of the stationary system, which is the
first of the two things the anchor needs; the second, a solve of such a
configuration that meets its own target, was still not reached. **Both
constants therefore stand where sections 9 and 10 left them and neither is
moved by this item.**

### 11.1 What changed in the acceptance

The completion flag of the stationary solve reads the species rows of the
system it solved, and not only the three hydrodynamic ones
(`stationary_rows_of_the_returned_state`, `steady_newton.f90`). A solve whose
unknown vector held `n(H2)` or `n(H+)` and whose residual drove those balances
to zero cannot end `info = 0` while one of them is above `cert_tol_carrier`:
reading only the hydrodynamic rows called such a state converged on three of
its four or five equations. A balance the system did **not** carry is still
measured and reported by the certification but does not decide the flag; it is
an equation the solve was never given.

### 11.2 Why the anchor is still not measured

The anchor's method (sections 9 and 10) is the row of a state whose solve met
its own target, with a decade of margin and at least a decade above the row's
own floor. MEASURED, 2026-09-07, on the two configurations that were brought
to a stationary solve for this item (hot Uranus, `Load IC`, `Solver: Newton
100.0`, `Coupled carrier solve`, secondary ionization off, WENO3 throughout):

| configuration | species rows | carrier row of the returned state | flag |
|---|---|---|---|
| H2 alone | 1 | 7.464e-03 volume-weighted, 1.49e-01 at cell 208 | `info = 2`, refused on the carrier row |
| H2 and H+ (`Ionization transport`) | 2 | 5.822e-03 volume-weighted, 8.43e-02 at cell 214 | `info = 2`, refused on the carrier row |

Neither solve met its own target, so neither row is an anchor: both are five
to six decades above `cert_tol_carrier = 1e-8` and what they measure is a solve
that stopped, not a stationary state. The elemental rows were not exercised at
all, because the elemental balances are not yet rows of the stationary system
(B5 target 1 for the element rows is not implemented; the report of the item
says what is left).

The observation of section 9 therefore stands unchanged and is now stated with
a stationary system behind it rather than only with a diagnostic: the carrier
tolerance is **UNCONFIRMED**, and the first configuration whose stationary
solve converges with a carrier row will be its anchor.

### 11.3 What B5b changed, and the anchors after it (2026-09-07)

Two things moved, and neither of them moves a constant.

**The species rows now measure the terms the marching operator carries.** The
advective term of a carrier row and of an elemental row was a cell-velocity
upwind difference while the marching stages advected the same species as the
divergence of the face species fluxes; item B5b made both the same expression
(`docs/steady_solver_design.md` section 14). A row's VALUE therefore changed
on every state that is not stationary in mass, because
`div(F_rho Y) = Y div(F_rho) + F_rho grad Y` and the cell-velocity form was
the second term alone. MEASURED on the 300-step relaxation snapshots of the
matrix, one thread:

| case | row | before | after |
|---|---|---|---|
| `mol_carrier` | carrier balance H2 | 3.912e-01 at cell 210 | 8.358e-01 at cell 211 |
| `mol_diffusion` | elemental transport He/H | 5.637e-04 at cell 212 | 9.987e-01 at cell 212 |
| `lower_profile` | elemental transport He/H | 4.799e-03 at cell 246 | 9.832e-01 at cell 246 |

Those states carry mass rows of 1.09 to 1.91 of their own terms, so a species
row that now reads O(1) is reporting the mass row's own imbalance, which is
the honest statement: a state that is not stationary in mass is not stationary
in a species it carries either. None of it touches a marching profile; the
three cases are byte-identical in `Hydro_ioniz.txt` and `Ion_species.txt`.

**The elemental balances are reported one entry per element.** The trace entry
used to be the worst element at each cell, which could not say which element
was out and could not be a row of a stationary system. There are now
`elemental transport <name>` entries, and an element the run holds none of at
the base is `not applicable` rather than a satisfied balance. MEASURED on
`lower_profile`, which carries C, O and N and no other metal: C 9.774e-01,
O 9.212e-01, N 9.683e-01, and Mg, Si, Ca, Na, K, S, Fe not applicable.

**The anchors are unchanged and still UNCONFIRMED.** `cert_tol_carrier` and
`cert_tol_element` stand at 1e-8. No stationary solve of a configuration with
a species row met its own target in this item either, so the method of
sections 9 and 10 still has nothing to anchor on. **No number was chosen to
pass.**

### 11.4 What B5d measured, and why the anchors still have nothing to stand on (2026-09-08)

B5d asked whether the element solve's failure to converge is the equation or the
Newton model, and the answer is the model
(`docs/steady_solver_design.md` section 16). Two defects of the model were
found and repaired, and both are repairs to the SOLVE and not to the
acceptance:

* the base element row was structurally empty, because cell 1 is the Dirichlet
  reservoir of the element operator and carries no transport balance, so the
  banded model was singular by construction. Its equation is now the boundary
  condition the operator states;
* the element row scale, which is the sum of the row's own terms floored at
  1e-20 of the flow-time rate, collapses on a damaged iterate and turns an
  ordinary derivative into a scaled band entry of 4.890e+22. The Newton's row
  scale is now floored at the flow-time rate of the transported quantity.

**Neither touches a certification measure.** The certification calls
`element_transport_residual` itself, over cells 2..N and on the operator's own
floor, so the row measures, the tolerances and the completion flag are the
same expressions they were.

**And no anchor was produced.** MEASURED on the hot Uranus reload, the
configuration B5b used, one thread:

| model | how the solve ended | `||R||` returned | elemental He/H row of the returned state |
|---|---|---|---|
| as B5b left it | `no descent direction exists`, iteration 8, `\|\|grad merit\|\| = 1.104e+44` | 1.974 | 9.987e-01 at cell 212 |
| `EXHALE_TRUST_REGION=1`, same text | radius collapsed to 3.8e-16, iteration 24 | 1.528 | 2.560e-01 at cell 195 |
| with both repairs | `no descent direction exists`, iteration 18, `\|\|grad merit\|\| = 4.813e+03` | 1.723 | 9.994e-01 at cell 192 |

None of the three met its own target, so none of their rows measures a
stationary state and none is an anchor by the method of sections 9 and 10.
**`cert_tol_carrier` and `cert_tol_element` stand at 1e-8 and remain
UNCONFIRMED. No number was chosen to pass.** What changed is that the third
row's refusal is now a statement about the state: the gradient it reports is a
number the model's own entries can produce, where the first row's was not.

### 11.5 The first stationary element solve, and what it anchors (2026-09-08, B5e)

B5d closed with "no stationary solve with an element row met its own target,
so the anchors have nothing to stand on". That is no longer true.

MEASURED on the same case B5d and B5b used -- `mol_diffusion` (hot Uranus,
`He_diffusion: True`, `He_Kzz: 1.0e9`, molecular base) reloaded from its own
300-step marching snapshot with `Coupled carrier solve: True`,
`Resid tol: 1.0e-8`, `Secondary_ionization: False`, one thread -- under the
scaled trust region, which is now the step control of the coupled route
(`docs/steady_solver_design.md` section 17):

```
(JFNK) it  66  ||R||=  5.485E-09
(JFNK) [TR] steps accepted 54, rejected 12
(JFNK) residual of the state handed back at its own composition: ||R|| = 4.345E-09
```

and the certification of the state handed back:

| row | measure | tolerance | verdict |
|---|---|---|---|
| hydrodynamic mass | 6.211e-11 at cell 1 | 3.0e-12 | ABOVE |
| hydrodynamic momentum | 3.966e-09 at cell 500 | 1.0e-08 | within |
| hydrodynamic energy | 4.345e-09 at cell 6 | 1.0e-06 | within |
| **elemental transport He/H partition** | **4.633e-11 at cell 500** (volume-weighted 2.797e-12) | **1.0e-08** | **within** |
| level balance He 2^3S | 4.448e-19 | 1.0e-06 | within |
| eliminated-species closure | 8.444e-18 | 1.0e-06 | within |

**`cert_tol_element` = 1e-8 now has a measurement under it, and the
measurement supports it.** A converged stationary element solve reaches
4.633e-11 on the maximum over cells, 215 times inside the tolerance, and
2.797e-12 on the volume-weighted companion. That is the same relation the
hydrodynamic anchors of section 9.5 were set by: the tolerance stands about
two decades above what a converged solve was measured to reach.

**The number is NOT moved, and this is why.** By the method of sections 9 and
10 an anchor wants a ladder -- the achieved row measured against two or three
values of the solve's own target -- so that a row limited by the target can be
told from a row at its floor. This is one rung. What it settles is the
question B5d could not answer: whether any state reaches the tolerance at all.
It does, with two decades to spare, so 1e-8 is not a condition no solve can
meet (which is what section 9.2 found `cert_tol_hydro = resid_th` to be).
Tightening it wants the ladder, and the ladder wants two more solves of this
case at `Resid tol` 1e-9 and 1e-10.

**`cert_tol_carrier` stays 1e-8 and stays UNCONFIRMED.** MEASURED on
`mol_carrier` reloaded from its own 12000-step snapshot with
`Coupled carrier solve: True`, `Solver: Newton 100.0`, `Resid tol: 1.0e-8`,
one thread, neither step control reaches a stationary state:

| step control | how it ended | `||R||` | H2 carrier row of the returned state |
|---|---|---|---|
| line search | `no descent direction exists`, iteration 16 | 1.855 | 2.081e-01 max, 1.27e-02 volume-weighted, at cell 205 |
| trust region | STAGNATED at iteration 48 | 1.873 | 1.617e-01 max, 1.40e-02 volume-weighted, at cell 205 |

The obstruction is named by the new exit ledger and it is not the step
length: **18 of the trust region's 49 iterations ended on "no Krylov
direction could be sampled"**, and 500 trials were refused for the element
budget (115 on the line-search route). Seven decades separate both carrier
rows from 1e-8, which is where section 9.5 left them, so this study still
cannot say what a converged carrier balance reaches.

**One caveat on the state above.** It is not certified: the hydrodynamic mass
row of cell 1 reads 6.211e-11 against 3.0e-12. That is the base cell of the
molecular gate and a hydrodynamic matter; it is not an element-row failure,
and it does not affect what the element row measures.

### 11.6 The `Resid tol` ladder on the element case, and the floor it finds (2026-09-08, B5f)

Section 11.5 closed with "tightening it wants the ladder, and the ladder
wants two more solves of this case at `Resid tol` 1e-9 and 1e-10". Those two
solves were made, together with the 1e-8 rung re-measured on the current
tree, and the ladder does not exist: **the three targets give one state.**

MEASURED. The same case section 11.5 used, one privately built binary of the
current tree, `OMP_NUM_THREADS=1`, one run at a time, on scratch copies
outside the repository, with only the `Resid tol` line differing and bitwise
the same restart files:

| `Resid tol` | outer its | accepted / rejected / model refused | final radius | how it ended | `\|\|R\|\|` handed back | elemental He/H row |
|---|---|---|---|---|---|---|
| 1.0e-8 | 66 | 54 / 12 / 8 | 3.203e-12 | the target was met at the loop top | 4.345e-09 | 4.633e-11 (vol 2.797e-12) |
| 1.0e-9 | 88 | 58 / 31 / 27 | 9.322e-23 | ABORTED, no descent step in 12 consecutive iterations; best iterate returned | 4.345e-09 | 4.633e-11 (vol 2.797e-12) |
| 1.0e-10 | 88 | 58 / 31 / 27 | 9.322e-23 | the same abort, the same ledger | 4.345e-09 | 4.633e-11 (vol 2.797e-12) |

Every certification row of the state handed back is identical at all three
targets, and so is `log10 Mdot = 10.60` and the refusing row (hydrodynamic
mass, cell 1, 6.211e-11 against 3.0e-12, `info = 2`). The 1e-9 and 1e-10 runs
write bitwise identical output files, and against the 1e-8 run they differ in
four of the 504 rows of `Hydro_ioniz.txt` and two of `Ion_species.txt`: the
two base ghost rows in the twelfth digit of the velocity and the two outer
ghost rows in the sixteenth of pressure and temperature. **Every physical
cell is byte-identical.** The 1e-8 rung reproduces section 11.5 digit for
digit, so COST7's change to the coupled loop's stopping test leaves this
stationary solve untouched.

Where the two tighter runs stall is in their exit ledger, and it is the ray
test's own arithmetic: of 31 refusals, 9 ended on "the ray difference is
below its cancellation floor". That difference is two merit evaluations
differenced, and below `eq_sweep_reltol` times their size it is the noise of
the composition elimination and carries no verdict; with no verdict the
region shrinks, which makes the next difference smaller still, and the radius
runs to 9.3e-23 in twelve iterations.

**The row is FLOOR-limited, and the ladder that shows it is the one the run
prints at its own iterates.** The solver reports the certification's row
measures at every accepted iterate, which is a ladder in the achieved
residual instead of in the target. MEASURED, from the 1e-9 run:

| `\|\|R\|\|` at the iterate | elemental transport He/H row |
|---|---|
| 1.122e-02 | 4.213e-05 |
| 6.094e-04 | 1.919e-05 |
| 7.026e-05 | 8.273e-06 |
| 9.695e-06 | 4.906e-07 |
| 1.834e-06 | 1.838e-08 |
| 2.949e-07 | 1.266e-09 |
| 1.056e-07 | 4.700e-11 |
| 1.518e-08 | 4.678e-11 |
| 1.734e-08 | 4.638e-11 |
| 4.054e-09 | 4.633e-11 |
| 9.099e-09 | 4.631e-11 |

Two decades of `||R||` above the knee buy 2.5 decades of the element row, so
there the row is limited by the solve. Below `||R||` about 1e-07 it stops:
the last factor 26 in `||R||` moves it by 1.5 per cent and the last factor 3
not at all. **4.63e-11 is a floor of the row, reached three decades before
the solve stops.**

**The homopause caveat of section 9.5 did not materialize on this
configuration.** The element row's maximum is at cell 500, r = 4.62, the
outer wind edge, in all three runs; the `K_zz = 1.0e9` base cells, which the
caveat says cannot meet 1e-8, are not where this row is worst. That is one
configuration and not a refutation of the caveat, whose own reference case is
`K_zz = 2e12`.

**The proposed value, the rule, and why it is not applied.** By the rule of
sections 9.5 and 10.5 -- a factor of about ten of margin over what a
converged solve was measured to reach, standing at least a decade above the
row's arithmetic floor -- an anchor on the measured floor gives 4.633e-11
times ten, that is **5e-10**. `cert_tol_element` is left at **1.0e-8**, and
the three reasons are these.

1. There is one measurement and not a ladder. Two of the three rungs are not
   solves that met their target and they return the first rung's state.
2. The second half of the rule fails. A value a decade above ONE
   configuration's floor is not a decade above the row's floor in general,
   and the operator states its own reachable floor with a range (READ, the
   header of the element Newton in `binary_element_diffusion.f90`): 1e-12 in
   a wind, about 1e-5 in a `K_zz`-dominated homopause column. 5e-10 would
   refuse every configuration whose floor sits above 5e-11.
3. No second configuration exists to measure. `lower_profile` is the only
   other case of the matrix carrying an elemental transport row and it stops
   at its step cap without entering the steady route (section 10.3).

What changes is the evidence under the number: from `newton_floor` of the
composition Newton, READ, to a floor MEASURED on one configuration at
4.633e-11, 215 times inside the tolerance. **No number was chosen to pass and
no golden was refreshed.**

**The mechanism of the floor: two candidates eliminated, one open.** At
`dt_stationary = 1e30` the time term is identically zero and the row is
`K (s_R J_R - s_L J_L) + div(F_rho X)` measured over the sum of the
magnitudes of those same terms.

* The row scale floor is RULED OUT by that construction: a cell whose scale
  sits at `1e-20 rho X_base/(R0/v0)` is a cell whose transport terms all
  vanish, so its measure cannot be the maximum over the column.
* Round-off of the operator's own terms is RULED OUT by size: pure round-off
  of terms whose magnitudes form the denominator would read a few times
  2.2e-16, and 4.633e-11 is 2.1e5 times that.
* Whether the residue is the interpolation of the face fluxes at the outer
  boundary or a real stationary imbalance the model cannot descend is NOT
  SETTLED. Both are consistent with what was measured and not with each
  other.

Two circumstantial facts, recorded as such. The row's maximum is at the last
physical cell, whose outer face is the domain boundary, and the hydrodynamic
momentum row peaks at the same cell. And the row's advective half is
`div(F_rho X)`, which for a nearly uniform `X` is `X div(F_rho)`, the mass
row's own imbalance of the same state, and that mass row reads 6.211e-11, the
same order as the element row's floor. Neither is a measurement of the
element row cell by cell.

**What the next attempt needs.** The one-ulp instrument of sections 9.4 and
10.4, extended to the element row, which means `res_he`, `sc_he` and the two
term halves printed at the row's maximum cell. Nothing carries them today:
`EXHALE_RESIDUAL=1` writes `output/residual_profile.txt` with the three
hydrodynamic rows only, and the certification reports the species rows as a
maximum, a volume-weighted companion and a cell index. Any anchoring of
`cert_tol_element` or `cert_tol_carrier` on a floor will want the terms at
that cell.

**One measurement that did not work, recorded so it is not repeated.** The
1e-9 target re-run with `EXHALE_RESID_EQ_TOL=1.0e-10` does not reach deeper:
the solve leaves the loaded state and collapses at iteration 17 at
`||R|| = 1.958`, 5 trials accepted of 18, element row 9.993e-01, with the
elimination costing 15.7 sweeps per residual evaluation against 14.1 at the
default. So the coupled route's convergence on this case depends on
`eq_sweep_reltol` and 1.0e-8 is the value it works at; the run says nothing
about whether the element row's floor is the elimination's accuracy.

> **Superseded 2026-09-10 (decision 22 a, item N30):** the element and carrier
> row tolerances are no longer the single 1e-8 quoted here; they are 1e-5,
> gating for cells at r >= 1.20 (`cert_regime_wind_r`), and the rows of
> cells below that radius are reported and do not gate. Anchoring in
> `docs/certification_tolerance_anchoring_20260910.md`.
