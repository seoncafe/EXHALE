# B3a design: the attempted-step controller

Status: **design for approval**. Nothing in this document is implemented. It
is the document `PLAN_20260906_rev2.md` step B3a is briefed from, and it is
written so the advisor can approve or change the boundary, the checkpoint
list, the acceptance set and the rejection policy before any code is written.

Provenance. `src/EXHALE_main.f90` was READ at **11:12 KST on 2026-09-06**
(file mtime 10:23). A0-impl owns that file and is editing it, so every line
number below carries drift; the anchors that do not drift are the routine
names and the operator-split seams (`u_umA`, `u_umB`, `u_umC`, `u_umD`,
`update_map_end_step`). Module line numbers quoted from
`docs/b1_target_system_20260906.md` section 7 were READ there at 10:29 and
some had already moved by 11:12: `ionization_equilibrium.f90`'s
`bg_cell_adopted` had shifted by about sixty lines from where B1's table put
it, because A0-impl turned `ieq_marching_ledger` into an array of two. The
identity of the item, not its line, is what the checkpoint is written against,
and on 2026-09-07 (item DOCS-LINES) the citations of this document were
re-anchored to routine and block names accordingly.

Inputs READ: `docs/b1_target_system_20260906.md` section 7 (T7.1 to T7.4,
AT-7) and section 1.3/1.5 (T1.5, AT-1a to AT-1c);
`docs/a0_run_mode_contract_20260906.md`;
`docs/a2_certification_contract_20260906.md` section 3;
`docs/PLAN_20260906_rev2.md` (B2, B3a, B3b, B3c, B4);
`docs/PLAN_20260906_review2.md` sections 5.1 to 5.3;
`docs/PLAN_20260906_review.md` R10; the worker reports
`planrev2_A3_report.md`, `planrev2_B2_report.md`, `planrev2_A2_report.md`,
`planrev2_A0impl_report.md` and the brief `planrev2_A0impl_brief.md`.

---

## 1. Where the controller sits

### 1.1 The step in program order, as it is at 11:12 KST

READ from `src/EXHALE_main.f90`. The numbering is B1 section 7.1's, kept so
the two documents index the same operations.

| # | operation | lines (11:12) | has its own rejection today? |
|---|---|---|---|
| 0 | `eval_dt(W,dt,dt_loc)`, update-map scaling of `dt`, `dt_loc` | 1081-1091 | no (it sets the step) |
| 1 | `u_old = u` | 1104 | the anchor of the inner loop |
| 2 | `retry_step`: three RK stages, `Apply_BC`, `positivity_limited_fluxes`, the `EXHALE_REJECT_STEP` probe | 1145-1285 | **yes**: positivity test, `dt` halving, `n_dt_halve_max = 20` |
| -- | A0-impl's acceptance point: `n_steps_accepted`, `t_phys += dt*R0/v0` | 1302-1314 | (moves to the new boundary, section 5) |
| -- | `u_umA = u` | 1323 | |
| 3 | `U_to_W`, `get_species_densities`, `comp_T_from_p` | 1331-1341 | no |
| 4 | `element_diffusion_step` (`he_diffusion`) | 1348 | no |
| 5 | `photochemical_transport_step` | 1358 | **yes**: A3's carrier controller, internal to the operator |
| 6 | `excited_H_update` (`use_excited_H`) | 1363 | no |
| 7 | `ioniz_eq` | 1373 | partly: acceptance classes by cell, no step-level rejection |
| 8 | `get_species_densities`, `comp_p_from_T`, `W` assembly, `W_to_U` (the `C1` projection) | 1376-1388 | no |
| -- | `u_umB = u` | 1391 | |
| 9 | `solve_energy_semi_implicit`, or the explicit `u(3,:)` update | 1403 / 1406 | **yes as a rejection**, in both run modes: the routine reports its status, assembles nothing, and the controller refuses the attempted step (`as_reject_energy`) |
| 10 | `Apply_BC(u)` | 1409 | no |
| -- | `u_umC = u` | 1411 | |
| 11 | `viscous_conduction_step` with its `U_to_W`, `comp_T_from_p`, `Apply_BC` | 1422-1426 | no (it clamps at a floor) |
| -- | `u_umD = u` | 1429 | |
| 12 | `shapiro_filter(u)`, `Apply_BC(u)` | 1431-1437 | no |
| -- | **adoption boundary** `update_map_end_step` | 1441 | |
| 13 | post-adoption reads, control flags, staged secondary-ionization flips | 1444 onward | not part of the step |

So three of the fourteen operations can refuse something today, and the three
refusals are of three different kinds:

- **Row 2** refuses an *attempt* and retakes it at half `dt`. This is a real
  rejection with a restore (`u = u_old`) and a retry cap, but its criterion is
  admissibility (`rho > 0`, `rho e > 0`) and nothing else, and its restore
  covers `u` only.
- **Row 5** refuses a *substep of one operator* and retakes it, entirely
  inside `diffusive_photochemistry`. Its criterion is A1's returned-state
  verdict, its checkpoint is `carrier_checkpoint`, its cap is
  `carrier_retry_max = 8`, and on exhaustion it restores the interval entry
  state and stops. Nothing outside the operator learns that a substep was
  refused.
- **Row 9** does not refuse; it *stops*. B2's header already records that B3a
  replaces the stop by a rejection of the whole attempted step.

### 1.2 The boundary

**The single adoption boundary of the attempted step is placed immediately
before the `update_map_end_step` call of the marching loop
(`EXHALE_main.f90`), and everything from `u_old = u`
(row 1) through the Shapiro filter (row 12) becomes one trial.** Three reasons,
in order of weight:

1. It is the first point at which the physical state of the step is complete.
   Row 12 is the last operation that writes `u`; row 13 only reads it.
2. It is where B1 section 7.1 places it, and where the update-map
   instrumentation already places its last seam, so the boundary and the
   diagnostic that measures the operator split coincide.
3. A0-impl's clock is currently one boundary too early (section 5), and the
   comment at 1296-1301 says so in the source.

The controller is a loop around rows 1 to 12, not a wrapper on each row. Rows
2, 5 and 9 keep their own local behavior but stop owning the outcome of the
step (section 4.4).

**What stays outside the trial.** Row 0 sets the step and is re-executed for
each outer attempt with the reduced `dt`, so it is inside the retry loop but
before the checkpoint restore point: the controller owns `dt` and `dt_loc`
directly and does not restore them from the checkpoint (it sets them). Row 13
is outside because it changes control state, not physical state, and because
the staged secondary-ionization flip there is mode state that belongs to A0.

---

## 2. The checkpoint

### 2.1 What is saved, by module

The rule that decides membership: **an item is in the checkpoint if some
operation of rows 1 to 12 writes it and some later read of it, in this step or
a later one, precedes a write.** Derived quantities that row 3 rebuilds from
`u` and `f_sp` before anything reads them are not saved; they are hashed in
the tests (section 7) so the rule itself is checked and not assumed.

| group | items | shape at N = 500 | reuse |
|---|---|---|---|
| **main, physical** | `u` | (3, 504) | `u_old` exists and is exactly this; the controller renames its role from "inner retry anchor" to "the step's checkpoint of `u`" |
| **main, composition** | `f_sp` | (504, 40) | new; nothing saves it today |
| **main, source fields** | `heat`, `cool`, `eta` | 3 x 504 | new. Written by row 7, read by row 9 |
| **main, step size** | `dt`, `dt_loc` | scalar + 504 | owned by the controller, not restored from the checkpoint |
| **carriers** (`diffusive_photochemistry`) | the whole of A3's `carrier_checkpoint`: `fc`, the `carrier_module_state` pair of A2 (11 `cbg_*`, `cph_klw`, `cph_jh2o`, `cph_joh`, `row_terms`, `row_terms_phys`, `col_scale_car`, `row_scale_car`, `headroom_car`, `headroom_set`, `pct_worst_j`, `pct_worst_ic`), `pct_Dco`, `pct_cell_constrained`, the `pct_*` solve diagnostics and `pct_verdict`. The CO ceiling record left this list when the ceiling was deleted on 2026-09-06 (item CEILING-DEL), and the CO destruction domain record that replaced it is deliberately **not** checkpointed | see 2.3 | **reused whole**: `carrier_checkpoint_take` / `_restore` / `_matches` already exist and already assert their own completeness. B3a calls them; it does not re-enumerate the module |
| **element diffusion** (`binary_element_diffusion`) | `he_fraction_over_one`, `he_fraction_under_zero`, `he_fraction_newton_steps`, `he_fraction_newton_resid`, `trace_ratio_under_zero` | 5 scalars | new accessor pair needed; all five are `protected`, so the module must expose take/restore of its own |
| **excited hydrogen** (`parameters` + `excited_hydrogen` + `lya_rt`) | `gph_balmer_HI`, `heat_balmer`, `Jlya_arr`, `n2s_arr`, `n2p_arr`, `Sproton_arr`, `Hpe_arr`, `Hdx_arr`; `jlya_rt_loaded`, `jlya_rt_grid`, `Tdiag`, `nhidiag`, `nediag`, `nhiidiag`, `taulya`; `jint_arr`, `jstar_arr` | 16 x 504 + 1 logical | new |
| **ionization sweep** (`ionization_equilibrium`) | `nmol_eq` (504,4); `NH2_col_lw`, `f_shield_lw`, `k_lw_diss`, `tr_lines_lw`, `p_lw_single`, `p_lw_absorbed`, `P_H2_eq`; `nox_eq` (504,3), `n_o1d_eq`, `NH2O_col`, `NOH_col`, `heat_fuv`, `heat_chem`; `j_h2o_fuv`, `j_oh_fuv`, `tau_fuv`, each (504,5); `bg_cell`, `bg_cell_adopted`, `bg_cell_best`, each 504 of `ion_rates`; `ieq_nonroot_streak`; `ieq_sweep_state_kind`; the three ledgers `ieq_marching_ledger(2)`, `ieq_steady_iterate_ledger`, `ieq_steady_candidate_ledger`; `ieq_acc_nprint` | see 2.3 | new; follows A2's `save_carrier_module_state` pattern exactly, including the `_matches` round-trip assertion |
| **energy update** (`energy_semi_implicit`) | `n_energy_floor_hits`, `n_energy_floor_hits_family(2)`, `energy_floor_first_step`, `energy_floor_last_step`, `energy_floor_cell_hits(:)`, and B2's `energy_update_last_*` status fields | 504 integers + scalars | new accessor pair |
| **conduction** (`viscous_conduction`) | `n_conduction_floor_hits`, `n_conduction_floor_hits_family(2)`, `conduction_floor_first_step`, `conduction_floor_last_step`, `conduction_floor_cell_hits(:)` | 504 integers + scalars | new accessor pair |
| **hydro attempt counters** | `n_faces_flux_positivity_limited`, `n_faces_flux_positivity_limited_accepted`, `n_dt_halve`, `n_steps_dt_halved`, `n_dt_halvings` | scalars | attempt statistics, class 2 below |
| **Shapiro** | nothing: the filter holds no module state (READ, `Apply_BC.f90`) | | |

A design decision worth stating: the checkpoint is **one derived type per
owning module, taken and restored by that module**, and one aggregate in the
controller that holds them. The alternative, a single flat list in
`EXHALE_main`, was rejected because it puts the enumeration of a module's
private state in a file that cannot see it, which is exactly the drift A2 and
A3 avoided by having the module own its own save routine.

### 2.2 The three counter classes and their restore policy (T7.2)

| class | members | on a rejected attempt |
|---|---|---|
| **physical accumulations** | the reaction and heating budgets, `t_phys`, `n_steps_accepted`, the physical ledger family. The CO ceiling record stood here until the ceiling was deleted on 2026-09-06 (item CEILING-DEL); the CO destruction domain record that replaced it is **not** in this class, because whether a state lay inside the destruction model's domain has an answer whether or not the step that read it was accepted | **restored**, so a rejected trial leaves no contribution (review 2 section 5.3, A0 section 5) |
| **attempt statistics** | `n_steps_attempted`, `n_dt_halve`, `n_steps_dt_halved`, `n_dt_halvings`, the new outer `n_step_rejections`, A3's `carrier_substeps_*` and `carrier_intervals_retried` | **not restored**: they count attempts |
| **diagnostic extrema** | `he_fraction_over_one`, `he_fraction_under_zero`, `trace_ratio_under_zero`, `he_fraction_newton_resid`, the floor-hit counters of rows 9 and 11, `pct_worst_limit` | **kept, not restored** (B1 decision 10, 2026-09-06), labeled at the declaration as an attempt statistic, and never reported as a property of the adopted state |

The consequence for the checkpoint type is that class 2 and class 3 items are
*in* the saved set only so the tests can hash them and assert they did **not**
move back; they are excluded from the restore. Stating this in the type, one
field group each, keeps the AT-7 (iii) and (iv) assertions on the same list.

### 2.3 Size and time cost, ESTIMATED from the declarations

All shapes READ from the declarations and allocation statements named above,
at `N = 500` (the `N` declaration in `parameters.f90`), `Ng = 2` so 504 cells,
`n_species = 40`, `n_carrier = 5` (H2, OH, H2O, CO, H+, the maximum),
`n_fuv_band = 5`, `real*8` = 8 bytes. `ion_rates` is 34 `real*8` plus 3
`logical`, ESTIMATED at 288 bytes with padding.

| group | bytes | note |
|---|---|---|
| `u` | 12,096 | already paid: this is `u_old` |
| `f_sp` | 161,280 | the single largest main-program item |
| `heat`, `cool`, `eta` | 12,096 | |
| carriers (A3's `carrier_checkpoint`) | about 203,000 | 12 grid vectors, 2 (504,5), 4 (504,5) arrays, 3 vectors, `fc`, `pct_Dco` |
| excited hydrogen and Ly-alpha | 64,512 | 16 grid vectors |
| ionization sweep, arrays | about 139,000 | `nmol_eq`, `nox_eq`, the 3 FUV band arrays, 12 grid vectors |
| ionization sweep, `ion_rates` x 3 | about 435,000 | `bg_cell`, `bg_cell_adopted`, `bg_cell_best` |
| floor counters, ledgers, scalars | about 6,000 | |
| **total** | **about 1.03 MB** | about 2.1 kB per cell, everything active |

Time: one take is one copy of about 1 MB, one restore the same. At a nominal
10 GB/s of memory bandwidth that is about 0.1 ms for a take plus a restore.
The comparison figure is MEASURED in the B2 report: `wasp_full` at 300 steps
runs 16.9 s, i.e. **56 ms per step**. So the checkpoint is ESTIMATED at about
**0.2% of a step** with every optional physics on, and less in a run that
carries neither carriers nor oxygen chemistry (the two largest groups are then
unallocated and are not copied). The estimate is a bandwidth argument, not a
measurement, and B3a's report must replace it with a MEASURED number from
`EXHALE_PROFILE`.

A cheaper variant exists and is **not** proposed: taking the checkpoint lazily,
each module saving its own state on first write within the step. It halves the
copy in configurations that skip operators, and it makes the completeness
assertion of AT-7 (i) impossible to state as one hash over one list. The 0.2%
does not buy that.

---

## 3. The evaluation at the boundary

A2 section 3, row "Physical time step", fixes what the controller must ask
for: the time-discrete balances of every active equation within tolerance,
admissibility, element and charge invariants to round-off, the energy identity
of the step, and the integration-error requirement. What follows is that list
against what exists at 11:12 KST.

### 3.1 What exists as a returned-state verdict today

| what | who returns it | measure and tolerance |
|---|---|---|
| **positivity of `u`** | `positive_density_and_internal_energy` (row 2) | `rho > 0`, `rho e > 0`; already a rejection |
| **carrier balances** | A1's returned-state acceptance, read through `pct_verdict` (row 5) | worst unconstrained row, `carrier_accept_tol` = the Newton floor, 1e-8 |
| **energy update** | B2's `energy_update_last_status` / `_reason` / `_cell` / `_residual` (row 9) | `abs(R)/scale <= energy_res_tol = 1e-9` at the returned temperature with the cooling evaluated there |
| **chemical root class per cell** | `ioniz_eq`'s acceptance classes and `ieq_marching_ledger` (row 7) | class 4 is the relaxation amnesty and is not a root; A0-impl's handoff check already refuses a class-4 cell at step 0 |
| **element and charge invariants, and the simplex** | inside `ioniz_eq`'s acceptance of each cell and in the carrier `limit_to_element_budget` | present as a test of each cell; not assembled into one step-level statement |

### 3.2 What does not exist today

**These are the operations of the fourteen with no returned-state verdict.**

| # | operation | what is missing | who owes it |
|---|---|---|---|
| 2 | RK stages and fluxes | the **time-discrete residual of the hydrodynamic rows**: `(u - u_old)/dt + (dF - S)` measured against a scale. A2 measures the *stationary* hydro rows; nothing measures the time-discrete ones | B3a defines it; the row scale can reuse `residual_row_scale` and its `row_terms_describe_state(u)` state guard |
| 3 | `U_to_W`, densities, `comp_T_from_p` | nothing to verify: it is a change of variables. An inadmissible result here is a rejection, not a verdict | -- |
| 4 | `element_diffusion_step` | no returned-state test at all. It clips the helium fraction to a range, counts the clip, and returns. Its balance (the He/H partition transport equation) is `unavailable` in A2 | B4 for the equation, B3a for the time-discrete form once B4 gives it |
| 6 | `excited_H_update` | no verdict. The H(n=2) level equations return `unavailable` in A2, and the closure is a lag: it is evaluated one outer pass behind the state it closes, and the lag error is not measured | B3b / B4 |
| 7 | `ioniz_eq` | the class of a cell is a **local root test of an equilibrium**, not a time-discrete residual of a composition *step*. There is no statement that the composition changed by the amount the rates over `dt` imply | B3c |
| 8 | the `C1` projection | no verdict, and it is itself the defect: it rewrites `p` and `u(3,:)` from the new composition with no source behind the change | B3c removes the producer; until then section 6 applies |
| 10 | `Apply_BC(u)` | no verdict. `check_base_inflow_is_subsonic` runs after the boundary and is diagnostic only | B4 (the characteristic face boundary) |
| 11 | `viscous_conduction_step` | no returned-state test; the temperature floor is an adopted clamp, counted, with no residual. B2's report names this as owed to the same treatment it gave row 9 | B2b |
| 12 | `shapiro_filter` | no verdict, and by construction none is possible while the energy it removes is unbudgeted. Default off (`shapiro_eps = -1`) | B3b (budget it) or leave it as a validity state of every step it runs on |

So **eight of the fourteen operations have no returned-state verdict today**
(rows 2, 4, 6, 7, 8, 10, 11, 12), of which rows 2, 4, 6, 10 and 11 have an
equation that could be measured and rows 8 and 12 have a change with no
equation behind it at all.

### 3.3 The acceptance predicate the controller evaluates

Proposed, in the order it should short-circuit (cheapest and most decisive
first):

1. **finiteness and positivity** of `u`, `T`, `f_sp` over the physical cells;
2. **element and charge invariants to round-off**: hydrogen, helium and every
   active trace element nucleus count, and `n_e` formed from the stage charges
   of the same composition, compared against the values row 1 saved plus the
   transport across the two faces of the domain. This is the one invariant
   that is available today, from the same simplex test `ioniz_eq` already
   performs per cell, assembled once;
3. **every active returned-state verdict of 3.1**: carriers, energy, and the
   chemical-root count being zero (a class-4 cell rejects a physical step, A2
   section 3);
4. **the time-discrete hydrodynamic residual** of 3.2 row 2, the one new
   measure B3a itself introduces;
5. **the thermal identity of the step** (T15, 2026-09-07; section 9 decision 7
   states it in full): the thermal energy the step changed equals the thermal
   source it was given plus the MEASURED transport contribution, to the energy
   row's own tolerance. `Delta u_form` is reported beside it as the reservoir's
   change and is not added to the row, because this code's `heat` and `cool`
   are the net thermal sources and already carry every reservoir transfer that
   passes through the thermal pool (B3c);
6. **the integration-error estimate** (3.4).

Entries 1 to 4 are checkable now. Entry 5 was not closable before B3c: the
`C1` projection of row 8 changed `u(3,:)` by an amount with no source, so the
identity failed by that amount by construction, and B1's AT-1a names exactly
this as its RED reference. B3c removed the projection and T15 restated the
identity as the thermal one above; it now closes at the energy row's tolerance
and gates. Entry 6 is a new mechanism, proposed below.

The evaluator itself should be **A2's `certification_evaluate` called in the
physical-step context**, not a second implementation. A2 already has the
context argument, already measures the hydro rows and the carrier balances on
a passed-in state with one row measure, and already reports `evaluated` /
`unavailable` / `not_applicable` per entry. What B3a adds to it is the
time-discrete row set; what B3a must not do is grow a parallel evaluator whose
tolerances can drift from the certification's.

### 3.4 The integration-error estimate (R10)

R10 requires an error policy distinct from the nonlinear convergence
tolerance, and names step doubling or a justified fixed-step refinement study
as either being sufficient.

**Proposed: step doubling at the whole-step level, on a duty cycle.** The
controller takes the step once at `dt` and, every `n_err_every` accepted steps,
takes the same step as two of `dt/2` from the same checkpoint, and forms

```
e = max over cells and rows of |u_(dt) - u_(dt/2, twice)| / (atol + rtol |u|)
```

with `e < 1` required and `dt` adjusted by the usual controller factor. The
two half steps use the same checkpoint restore, so the estimate is of the
whole operator split and not of one operator, which is the point R10 makes
about splitting error that an estimate of one operator misses.

Cost: one extra pair of steps every `n_err_every` steps, so a factor
`1 + 2/n_err_every` on run time. At `n_err_every = 20` that is +10%; at 50,
+4%. The duty cycle is the reason to prefer this over an embedded estimate.

**Why not an embedded estimate.** The hydro stage is SSP-RK3 and an embedded
pair could be read off it almost free, but it would measure the hydro
truncation error alone, and the error this code has reason to fear is the
splitting error across rows 4 to 12 and the backward-Euler error of the
implicit source updates. A3 MEASURED the carrier operator to be first order in
its substep (ratio 2.001), so an estimator that ignores the source operators
would be reporting the accuracy of the one stage that is third order. Step
doubling over the whole step is the estimate that matches the quantity.

**What the estimate must not become.** A3's report contains the finding that
governs this: the carrier acceptance measure carries an explicit `1/dt` time
term in its row scale, while the imbalance a clamped neighbor leaves is a
flux term that does not scale with `dt`, so **halving the step can bring any
bounded defect under a fixed relative floor without the defect changing**
(MEASURED there: a refusal at 1.042e-7 becomes an acceptance at 7.888e-9 after
10 attempts, same cell). An outer controller that reduces `dt` on rejection
therefore has a route to accepting a state by shrinking the measure rather
than by fixing the state. This is an open choice for the advisor (section 8,
choice 1) and it should be settled before the controller's acceptance is read
as evidence of anything.

---

## 4. Rejection policy

### 4.1 What is restored

Every group of section 2.1 except the class-2 and class-3 counters of 2.2. The
restore is by module, each module restoring what it saved, and the controller
then asserts the round trip with the `_matches` routines A2 and A3 established
(a mismatch is loud, not silent).

Row 3 rebuilds `W`, `rho`, `v`, `p`, `T`, the species densities, `n_e` and
`n_tot` from the restored `u` and `f_sp` at the top of the next attempt, so
those are not restored. The AT-7 (i) hash covers them anyway, so that this is
a checked property and not an assumption.

### 4.2 How `dt` is reduced

- **Factor**: 0.5 on an admissibility or verdict rejection, matching the
  inner positivity path so the two cannot disagree about what a halving is.
  On an integration-error rejection, the standard controller factor
  `max(0.2, 0.9 e^(-1/(p+1)))` with `p = 1` (the order the source operators
  set, MEASURED for the carrier operator by A3), floored at 0.2 so one
  rejection cannot collapse the step.
- **Retry cap**: `n_step_retry_max = 8`, i.e. a floor of `dt/2**8` for the
  outer controller. Chosen equal to A3's `carrier_retry_max` so the two caps
  are one number in two places and the failure of a carrier interval and the
  failure of a step are reached at the same depth. The hydro `n_dt_halve_max`
  is 20 today and stays what it is, for the reason in 4.4.
- **Floor**: the cap and the factor together define it; no second absolute
  floor is introduced. A rejection arriving with the step already at the floor
  is the exhaustion.
- **Growth**: `dt` returns to `eval_dt`'s value at the top of the next step,
  because row 0 is re-evaluated. The controller does not carry a shrunken
  step forward; if it should, that is choice 3 of section 8.

### 4.3 Counters, and the exhaustion action

On a rejection: `n_steps_attempted` advances (A0-impl already advances it once
per pass through `retry_step`; the outer controller advances it once per outer
attempt, and section 4.4 says how the two are kept from double counting),
`n_step_rejections` advances, the reason and the failing operation are
recorded, and nothing else moves. `t_phys`, `n_steps_accepted` and the
physical ledger family do not move, which is A0 section 5 and AT-7 (ii).

On exhaustion: restore the **last accepted state**, that is the checkpoint of
this step, print the reason, the operation that failed, the failing cell, the
measure and its tolerance, the number of attempts, the step reached and the
retry history (the `dt` of each attempt with the operation and reason that
refused it), then stop **in both run modes** with exit status 2, the status a
refused stationary claim uses (`certification_stop_uncertified`), so that one
status means one thing to a caller: the code could not produce an acceptable
state. The state written before the stop is the last accepted one and carries
its own `mode` and `t_phys`, so a run that ends this way still writes a
consistent state. This replaces B2's `error stop 1` and A3's stop for the
interval it owns.

### 4.4 How the existing rejection paths nest

| existing path | after B3a |
|---|---|
| **hydro `retry_step`** (row 2) | **kept as an inner loop, with its outcome demoted.** Its positivity criterion is local to one RK stage and cannot be replaced by an evaluation at the outer boundary, because the state it refuses is one the next stage would take a square root of. It keeps its own halving and its own cap. What it loses is the acceptance: `step_accepted` no longer means the step is accepted, it means the hydro stage produced an admissible state, and the outer controller decides. Its `n_steps_attempted` increment stays where it is; the outer controller increments a **separate** `n_outer_attempts`, and the end-of-run report prints both with their definitions, so the two counts are never added |
| **carrier controller** (row 5) | **kept, unchanged, absorbed as a subordinate.** Its refusals are substeps of one operator over one interval, at a frozen background, and they are the right granularity for that operator: forcing them up to the outer boundary would discard rows 1 to 4 for a failure the operator can fix at a tenth of the cost. What changes is its exhaustion: instead of `error stop 1`, `carrier_transport_interval` reports `completed = .false.` to the outer controller, which turns it into an outer rejection with the operation named. A3 already returns that flag and already restores the entry state, so this is a change at the call site in `EXHALE_main`, not in the operator |
| **energy failure** (row 9) | **stop became rejection** (2026-09-06, item 1 of `To_be_determined_by_user_20260906.md`). The stop and the flag that suppressed it are gone: the routine reports one compact line pair and returns its status, the controller reads it and refuses the attempted step at `as_op_energy` with reason `as_reject_energy`, in BOTH run modes. The initialization leniency is about an iterate that exists; a failed update produced no temperature at all |
| **conduction floor** (row 11) | unchanged until B2b, and its activity remains an unbudgeted accepted correction that the validity state reports |

The double-counting rule, stated once: **each loop counts its own attempts and
nothing else.** An outer attempt that contains three inner halvings and two
carrier substep refusals contributes 1 to `n_outer_attempts`, 3 to
`n_steps_attempted` (which becomes, in name and in the report, the hydro stage
attempt count), and 2 to `carrier_substeps_rejected`. Accepted steps are
counted at one place, the outer boundary.

---

## 5. Interaction with A0-impl

A0-impl has landed (report READ). The hooks it provides, by name, all in
`global_parameters` unless stated:

| name | what it is | what B3a does with it |
|---|---|---|
| `run_mode`, `run_mode_init`, `run_mode_phys` | the mode of the run | the controller's acceptance is required in `phys`; in `init` it runs as a diagnostic (section 6) |
| `t_phys` | the clock, seconds | **moves from its increment in the marching loop of `EXHALE_main.f90` to the outer adoption boundary.** The increment stays `dt*R0/v0` |
| `n_steps_accepted` | accepted steps | moves with the clock, to the same boundary |
| `n_steps_attempted` | attempts | **stays where it is**, at its increment inside the Runge-Kutta retry loop, and is renamed in the report to say it counts hydro stage attempts; the outer count is new |
| `ledger_family`, `ledger_family_init`, `ledger_family_phys` | the ledger-family switch | **the switch to `ledger_family_phys` moves to the adoption boundary**, so that no producer inside the trial can write to the physical family. A0-impl sets it at the accepted handoff and holds it at `init` around the stationary solve; B3a narrows the `phys` window to the adopted step only |
| `physical_handoff_check` (`EXHALE_main.f90`) | the step-0 admissibility and source-consistency test | unchanged; it is the entry condition, the controller is the condition of each step. Its three tests are the same ones entries 1 to 3 of section 3.3 evaluate, so the two should call one routine rather than two |
| `reject_step_probe` (`EXHALE_REJECT_STEP`) | forces a rejection at a named step, through the positivity path | **reused and extended**: B3a needs the same probe at each of the fourteen operations for its injected-failure tests (section 7), so the probe gains an operation index |
| `trace_step_clock` (`EXHALE_STEP_CLOCK`) | prints count, attempted, accepted, `dt`, `t_phys` per accepted step | extended with the outer attempt count and the rejection reason |
| `ic_run_mode`, `ic_run_mode_present`, `ic_t_phys` (`load_IC`) | restart mode and clock | untouched |
| `n_energy_floor_hits_family(2)`, `n_conduction_floor_hits_family(2)`, `ieq_marching_ledger(2)` | the family-indexed ledgers | their producers keep writing through `ledger_family`; the narrowed window of the row above is what makes a rejected trial write to neither family |

One consequence to state plainly: **A0-impl's own test
`a_rejected_trial_is_an_attempt_and_not_a_step` remains valid and its meaning
sharpens.** Today it exercises a rejection at the hydro boundary, which is the
only rejection there is. After B3a the same probe exercises a rejection at the
inner boundary, and a new probe exercises one at the outer boundary; both must
show no clock advance.

---

## 6. The interim contract

Review 2 section 5.1 and B1's T7.3, restated as the contract B3a certifies on
the day it lands, before B3c removes the `C1` composition reset.

**What B3a certifies in the interim.**

1. Rows 1 to 12 are restored on a rejected attempt, bit for bit, verified by
   hashing the saved set before and after (AT-7 (i)).
2. Physical accumulations are restored; attempt statistics are not; diagnostic
   extrema are kept and labeled (T7.2, decision 10).
3. `t_phys` and `n_steps_accepted` advance only after an adopted step
   (AT-7 (ii)).
4. No physical accumulation, and in particular no reaction or heating budget,
   carries a contribution from a rejected attempt (AT-7 (iii)); the attempt
   statistics do (AT-7 (iv)).
5. The acceptance predicate of section 3.3 entries 1 to 4 is evaluated at the
   boundary and a step failing any of them is rejected.
6. A failed energy update rejects the attempted step in both run modes, and
   the exhaustion of the retry budget ends the run with status 2 in both
   (2026-09-06, item 1 of `To_be_determined_by_user_20260906.md`). The interim
   stop inside `energy_semi_implicit` is gone, and so is the flag that
   suppressed it.

**What B3a does not certify, stated so it cannot be misread.**

- It does **not** certify that the step it adopts is a physically valid time
  step. Row 12's Shapiro filter is still in the trial, and A2 section 5
  refuses certification of any state it was active in. A step that the
  controller adopts while the filter is active is a step that was **restorable
  and admissible**, not a step whose every correction is budgeted; its thermal
  sources do now balance, which is what section 3.3 entry 5 gates.
- The thermal identity of the step (section 3.3 entry 5) was **evaluated and
  reported, and its failure was not a rejection** while the projection was in
  the trial. It gates now (T15), under `energy_res_tol` unchanged at 1e-9 and
  the carrier tolerance unchanged at 1e-8: **no tolerance was changed to make
  anything pass.** What is still carried without a budget is row 12, the
  Shapiro filter, so no state certifies while the filter is active.
- A demonstration that the controller restores state is a demonstration about
  memory, not about physics. The report of B3a must say this in its verdict,
  in those terms, and must not describe a passing rollback matrix as a
  validated step.

The interim contract ends when B3c lands: entry 5 becomes a gate, entries for
rows 7 and 8 gain their verdicts, and the sentence above about what is not
certified is deleted rather than weakened.

---

## 7. Tests

The suite is `src/tests/step_controller/`, on the pattern of
`src/tests/carrier_retry/` (a Fortran driver plus `run.sh` honoring a private
object directory), with the whole-binary rows as shell tests on the pattern of
`src/tests/run_mode/run.sh`.

### 7.1 The rollback matrix

| test | assertion | RED reference |
|---|---|---|
| **mutated trial fully restored** | perturb every item of the checkpoint list (one distinct value per item, allocating what is unallocated), restore, and require `_matches` true and the hash equal, item by item and allocation status included | today there is no checkpoint outside `u_old` and the carrier module: the test cannot compile |
| **injected failure at each of the fourteen operations** | fourteen rows. For each of rows 1 to 12, force a rejection immediately after it through the extended `reject_step_probe`, and require: the full restore of the row above, no clock advance, no accepted-step advance, no physical accumulation moved, the attempt statistics advanced by exactly one, and the reported failing operation equal to the injected one. Rows 0 and 13 are refusals of the test itself: injecting there must be reported as out of the trial | rows 4 to 12 leave module state mutated today; a retry re-enters with a partially updated `bg_cell_adopted`, `ieq_nonroot_streak` and the floor counters (T7.1) |
| **elapsed time correct** | over a run with `k` injected rejections, `t_phys` equals the sum of the accepted global `dt`, to the last bit, and the accepted count equals the attempted count minus the rejections. A0-impl's `EXHALE_STEP_CLOCK` already MEASURED this identity to the last bit for the inner path; the same measurement over the outer path | |
| **no rejected contribution in any accumulated record** | AT-7 (iii): with a rejection injected in a step whose carrier operator deposits reaction heat, the accepted heating and reaction budgets are bit-identical to their values before the step, and the attempted ledger carries the deposit. The CO destruction domain record is outside this row by construction: it is deliberately not restored (item CEILING-DEL, 2026-09-06), and `carrier_retry` asserts instead that it comes through a restore unchanged | A3 already asserts this for the carrier interval; the new row is the whole step |
| **exhaustion** | with every attempt refused, the run restores the last accepted state, names the operation and the reason, stops, and the state it wrote is the last accepted one with its own `t_phys` | |

### 7.2 The integration-error test

Two rows, both MEASURED and neither asserting a size:

1. **the estimate falls with the step**: on a fixed state, `e` at `dt`, at
   `dt/2` and at `dt/4`, with the ratio reported. For a first-order splitting
   the ratio is 2; the assertion is that halving reduces it, which is the form
   A3 used and the only form that is a property of the method rather than of
   the case.
2. **the estimate is of the split, not of the hydro**: the same measurement on
   a configuration with the source operators switched off must show the
   third-order ratio of the RK stage, and with them on the first-order one.
   This is what distinguishes the proposed estimator from an embedded one, and
   it should be measured rather than argued.

### 7.3 Byte identity

**Expected: a run that rejects nothing is byte-identical to the run before
B3a**, in the numeric content of every output file. The controller adds a copy
and an evaluation; it changes no arithmetic on the accepted path. Two
qualifications, both to be stated in the report rather than discovered:

- The evaluation itself must not mutate state. A2 met the same requirement
  with `save_carrier_module_state` around `carrier_steady_residual` and with
  `row_terms_describe_state(u)`; the controller's evaluation reuses both.
- The integration-error duty cycle **is not byte-identical** when it is on,
  because it takes extra steps from the checkpoint. It must therefore default
  to off (`n_err_every = 0`), and the plan's rule that a new physics or
  numerics option defaults to off applies to it unchanged.

The scoped impact measurement for the item is `wasp_full` and
`mol_base_handoff` on scratch copies, before and after, with the controller
active and rejecting nothing.

---

## 8. Open choices for the advisor

1. **The dt-scaled acceptance measure.** A3 MEASURED that the carrier
   acceptance row scale carries `1/dt`, so an outer controller that halves
   `dt` can convert a refusal into an acceptance without the state improving;
   decide whether the measure is fixed before the controller is relied on, or
   whether the controller ships with the finding recorded at its acceptance
   site.
2. **The retry cap and the floor.** `n_step_retry_max = 8` is proposed to match
   A3; the hydro path uses 20 today, and the two can be unified at one number
   or deliberately kept different.
3. **Step-size memory.** Whether a `dt` reduced by a rejection is carried into
   the next step (a real controller) or discarded because row 0 re-evaluates
   `eval_dt` each step (the proposal above).
4. **The integration-error duty cycle.** Whether `n_err_every` defaults to off
   with a recommended value, or is required in `phys` mode, given the +10% at
   20 steps.
5. **The evaluator's home.** Whether the physical-step context extends A2's
   `certification_evaluate` (proposed) or is a separate routine in the
   controller, given that A2's file is shared and the two tolerances must not
   drift.
6. **Ownership of `EXHALE_main.f90`.** A0-impl owns it now and B3a rewrites its
   marching loop; the two must be sequenced, and B3a should start from a state
   of that file the advisor names.
7. **The interim energy identity.** Whether T1.5 is reported at every step in
   the interim (a number at each step, useful as the AT-1a RED trace and costing an
   evaluation) or only when a probe asks.
8. **Rows 8 and 12.** Whether the `C1` projection and the Shapiro filter are
   simply carried inside the trial until B3c and B3b reach them, or whether the
   controller refuses to run in `phys` mode while either is active.

## 9. Advisor decisions on section 8 (2026-09-06)

1. The measure is fixed first: item A1scale (in progress) removes the `1/dt`
   time term from the returned-state and stationary carrier measures; B3a is
   briefed only after it lands, and the controller's acceptance reads the
   dt-independent measure.
2. Outer cap `n_step_retry_max = 8` and floor `dt/2^8` as A3; the nested hydro
   `retry_step` keeps its 20 for now, deliberately, and the two are
   reported separately; unification is measured later, not assumed.
3. Step-size memory: a `dt` reduced by a rejection bounds the next step
   (`dt_next = min(dt_CFL, 2 dt_accepted)`), a real controller; a run that
   rejects nothing is unchanged.
4. Integration error: in `phys` mode the step-doubling estimate runs every
   `n_err_every = 20` steps by default (the measured +10 percent is the
   price of a physical-mode claim); in `init` mode it is off unless asked.
5. The physical-step context extends A2's `certification_evaluate`; one
   evaluator, one set of tolerances.
6. B3a starts from the state of `EXHALE_main.f90` that A0-impl hands back,
   and owns the marching loop from then on.
7. The energy identity T1.5 is evaluated and printed every step in `phys`
   mode (a sum over cells, no new physics evaluation) and only on request
   in `init`.

   **What that line means now (T15, 2026-09-07).** The identity the
   instrument evaluates is the THERMAL identity of the step, not
   `Delta u_th + Delta u_form = integral Q_ext`. The reason is the B3c
   decision: in this code `heat` and `cool` are the NET THERMAL sources,
   so photoionization deposits `h nu - E_th` and never `h nu`, collisional
   ionization removes `E_th` as the `coio` cooling term, and every
   collisional reaction heat is a difference of the one species
   formation-energy table. The formation reservoir is therefore already
   accounted for by those two arrays, and adding `Delta u_form` to the row
   would count every collisional transfer twice and charge the gas for the
   ionization energy the photons paid. What the instrument reports and
   gates is

   ```text
   sum_j dV_j [u_th^{n+1} - u_th^n]_j
       =  sum_j dV_j dt (heat - cool)_j  +  (the transport contribution)
   ```

   with the transport contribution MEASURED, not modelled: the marching
   loop marks the thermal energy on both sides of the coupled source step,
   the difference of the marks is what the sources did, and the rest of the
   step's thermal change is by definition what the hydrodynamic stages, the
   diffusion and carrier operators, the boundary condition, the conduction
   stage and the filter did. The residual is then the closure of the source
   step alone. `Delta u_form` stays on the printed line as the reservoir's
   own change, labeled as such and never added to the row.

   **It is now a GATE in `phys` mode**, under the energy row's own
   tolerance and its own scale (`energy_res_tol = 1e-9`,
   `scale = sum_j dV_j [|u_th^n| + dt(|heat|+|cool|)]`, both read from
   `energy_semi_implicit`). No tolerance was chosen to make a run pass: the
   summed residual is bounded by the sum of the cell residuals, each of
   which the cell-local row already holds to `energy_res_tol` of its own
   scale, so the summed test cannot be stricter than the row test it
   aggregates. The gate is inert while the two marks are absent, because
   without them there is no measured transport contribution to subtract.

   MEASURED on 300 `Run mode: phys` steps, one tree with and without this
   change, worst `|residual/scale|` over the run: `wasp_full` **7.44e-1 ->
   5.40e-10**, `mol_base_handoff` **8.63e-1 -> 1.88e-10**. No step is
   refused by the gate in either case, and the state written
   (`Hydro_ioniz.txt`, `Ion_species.txt`) is byte-identical to the run
   without it apart from the provenance timestamp. The worst value sits
   within a factor of two of the tolerance, which is what an aggregate of
   cell rows each held to that same tolerance should give: the gate is a
   real bound, not a slack one.

8. Rows 8 (the composition projection) and 12 (the Shapiro filter) are
   carried inside the trial until B3c and B3b reach them, each recorded as
   an unbudgeted accepted correction (B6 category 4) when it fires, so the
   controller runs in `phys` mode but no state can certify while either is
   active; the interim contract of section 6 stays in force.

   *Closed for row 8, 2026-09-07.* B3c removed the producer of the `C1`
   composition projection, so row 8 can no longer contribute an unbudgeted
   accepted correction. Its counter `n_projection_applied`, the report block
   that printed it and the certification entry that summed it are deleted;
   only the Shapiro filter of row 12 remains under this decision, and it is
   still what refuses certification while it is active.
