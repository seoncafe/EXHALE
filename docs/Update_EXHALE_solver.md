# Update_EXHALE_solver: solver / numerics changes vs. the original ATES

The single, complete record of the **solver / numerics** changes in EXHALE
relative to the original ATES (`ATES/ATES-Code-main`): time integration of the
stiff source terms, the convergence algorithm, the radial-grid convergence window,
the advection post-processor's solve, and the nonlinear solvers (ionization
equilibrium and the energy equation in each cell). It **excludes** the physics changes
(metals, cooling, opacity, charge exchange, excited hydrogen, Lyα), which are in
`Update_EXHALE_stage1` and `Update_EXHALE_stage0`. It absorbs in full the
former `energy_semi_implicit_solver.tex` (§2) and the solver portion of
`convergence_fix_and_validation.tex` (§1); those standalone memos are superseded by
this file. Companion: `Update_EXHALE_solver.tex`/`.pdf`. Sections are chronological.

---

## 1. Convergence-algorithm fix (2026-06-02)

*Absorbed in full from the former `convergence_fix_and_validation.tex` (its Phase
1/2 **physics** validation (opacity dispatcher, tabulated opacity, trace-metal
smoke test) is intentionally not reproduced here; see `Update_EXHALE_stage0` /
`Update_EXHALE_stage1`).*

**Overview.** Fixed a non-terminating (effectively infinite) time-integration loop,
found on the first physical test case (HD 209458b; power-law SED, solar He/H,
default numerics).

**Symptom.** A full run never terminated; `du` settled on a plateau:

| iteration | du |
|---:|:--|
| 26,905 | 0.62 |
| 52,239 | 0.077 |
| 667,569 | 0.0109 (still running) |

**Diagnosis.** The loop in `EXHALE_main.f90` terminated on a single criterion,
`do while( .not.is_mom_const .or. force_start)` with `is_mom_const = (du < du_th)`,
`du_th = 1e-3`. The steady-state criterion `is_zero_dt = (dtu < dtu_th)` was
computed but **never used in the loop condition** (only in a commented-out print).
For this setup `du` floors at ≈0.011 (a transonic-wind numerical limit), ~10× above
`du_th`, so with `is_zero_dt` disconnected nothing could terminate. Measuring `dtu`
directly (near-converged IC) gave ≈5.5e-8: the solution *had* reached steady state
but the strict thresholds were physically unreachable.

**Fix** (in `EXHALE_main.f90` + `parameters.f90`):
1. **Restore `is_zero_dt`** as an OR-branch of the termination test.
2. **Loosen `du_th`** 1e-3 → 2e-2 (a deliberate, global user choice).
3. **Stall detection:** if the relative change in `du` stays below
   `stall_tol = 1e-6` for `N_stall = 2000` consecutive iterations, declare a plateau
   (steady state) and stop. `stall_tol = 1e-6` is chosen so a slowly-*decreasing*
   `du` (relative change ~3e-6 in the descent) is not mistaken for a plateau.
4. **Hard cap** `count_max = 1e6` + a termination report.

Resulting loop:
```fortran
do while( ( .not.is_mom_const .and. .not.is_zero_dt .and.
            .not.is_stalled   .and. count < count_max ) .or.
          force_start )
```
*Touched files:* `parameters.f90` (new flags/params, `du_th`); `EXHALE_main.f90`
(loop condition, stall logic, report). Stdout now prints `count, du, dtu`.

**Verification.** Plateau IC + loosened `du_th`: stops in 2 iterations
(`du=0.011<0.02`). Stall path (`du_th=1e-3`, before loosening): a previously
non-terminating run stopped at iteration 554,485 via stall detection at the
`du=0.011` floor (confirming `stall_tol` doesn't trip during the descent).
Profiles unchanged to 3 sig figs.

**Caveat.** `du_th=2e-2` is loose and global; tighten it (and rely on stall
detection) if a study needs a stricter steady state.
*2026-08-15 note: reverted to `du_th = 1e-3` (`parameters.f90`), the original
ATES-Code-main value; the two-stage `du_th [PLM,WENO3]` key now carries the
loose threshold on the PLM stage only.*

---

## 2. Semi-implicit energy solver (2026-06-03)

*Absorbed in full from the former `energy_semi_implicit_solver.tex`.* Replacing the
explicit forward-Euler energy update with a constant-derivative 2-iteration
Newton-Raphson scheme resolves the stiffness of radiative heating/cooling:
unconditional stability and ~2× speedup on the hydro loop.

**Introduction / problem.** Radiative `tau_cool = p/((g-1)C)` is far shorter than
the CFL step `tau_CFL = dr/(|v|+c_s)`. The original ATES integrated the source
explicitly, `u3^{n+1} = u3^n + dt*(H - C(T^n))`, with `dt` set by the hydro CFL,
unstable in the stiff regime: high-frequency T/p oscillations, NaN crashes when a
big step over-cools to `p<0`, and stalled convergence (>2e5 iterations).

**Formulation.** Operator-split; ρ, v fixed during the source step, so the energy
equation reduces to `n_all/(g-1) dT/dt = H - C(T)` (`n_all = n_tot + n_e`,
adimensional). Backward Euler: `T^{n+1} - T^n = dt*(g-1)/n_all*(H - C(T^{n+1}))`.

**Original explicit solver (before).** `u3^{n+1} = u3^n + dt*(H(T^n) - C(T^n))`,
then `p = (g-1)(u3 - 0.5 rho v^2)`, `T = p/(n_all k_B)`. Stability needs
`dt <= tau_cool`, but `tau_cool` can be 1e-3-1e-6 × `dt_CFL`, giving the
oscillations / NaNs / stalls above.

**Numerical algorithm.** Residual `F(T) = T - T^n - dt*(g-1)/n_all*(H - C(T))`,
`dF/dT = 1 + dt*(g-1)/n_all * dC/dT`. Two enhancements:
- **Non-negative derivative safeguard:** `dC/dT <- max(0, dC/dT)`, so `dF/dT >= 1`
  always → unconditional stability, clamps unphysical excursions.
- **Constant-derivative 2-iteration scheme:** `C(T)` is expensive (`eval_cool`:
  Gaunt factors, recombination, collisional ionization, metal lines). Compute
  `dC/dT` *once* per step from `T^n` and `dT = max(1e-5, 1e-5 T^n)`, then take 2
  Newton steps with that fixed derivative: T0=T^n; eval C(T0) (call 1) and
  C(T0+dT) (call 2) → derivative; iter 1 → T1; eval C(T1) (call 3); iter 2 → T2;
  set T^{n+1}=T2 with floor T>=0.01 (adimensional). 3 cooling evals/step vs 10 for a
  5-iteration Newton → 70% less cooling overhead, ~2× faster overall.

**Implementation.** New module `src/modules/time_step/energy_semi_implicit.f90`,
subroutine `solve_energy_semi_implicit(u, W, dt, heat, cool, f_sp)`. In
`EXHALE_main.f90` the explicit `u(:,3) = u(:,3) + dt*(heat - cool)` was replaced by
`call solve_energy_semi_implicit(u,W,dt,heat,cool,f_sp)`; Makefile updated.

**Verification (HD 209458b, metals on/off comparison).** Metals-ON converged in
**71,993** steps vs **236,344** explicit (~33% less wall-clock, 14 min 7 s). Same
state: log10 Ṁ 9.49061 vs 9.49064; T_max 6525.4 vs 6525.6 K; T_min 1009.3 K
(exact). Metals-OFF: explicit stalled at du≈0.033; semi-implicit reached du=0.02.

---

## 3. Escape-radius / empty-window grid guard (Phase 4)

*Context: `Update_EXHALE_stage1` Phase 4.* In `src/modules/init/define_grid.f90`, the
convergence diagnostics (`du`, `dtu`) are computed over `[j_min:N]`, where `j_min`
is the first cell with `r >= r_esc` (`Escape radius`). When the escape radius is
outside the (L1-truncated Roche) domain (`r_esc > r_max`, e.g. WASP-121b Case D:
r_max~1.35 < r_esc=1.5), no cell has `r >= r_esc`, the DO loop leaves `j = N+Ng+1`,
and `[j_min:N]` is **empty**. `maxval`/`minval` over the empty slice return ∓HUGE,
so `du=Inf` while `dtu=-HUGE < dtu_th` → a spurious "steady state" exit at step 0.

**Guard.** If `j_min > N`, re-anchor `j_min` to the first cell with
`r >= 1 + 0.5*(r_max-1)` (mid-domain) and print a loud WARNING to set a smaller
`Escape radius [R_p]:`. Converts a silent false-convergence into a visible,
recoverable warning.

---

## 4. Breathing-base advection guard (post-process, option c)

*Full account and validation: `Update_EXHALE_stage1`, "Post-process advection guard at a
breathing base".* In `src/modules/post_process/post_process_adv.f90`: for strongly
Roche-filling planets the 1D wind is subsonic at L1 and the dense base recirculates
(small *negative* inflow velocities, stagnation point v=0 near the wind base). The
advection post-processor re-solves, in each cell, the energy/ionization balance assuming
an outflow (upwinds from the next-inner cell); in the dense base the non-monotone
metal cooling gives the energy equation a spurious *hot* root, and the upwind
coupling cascades it (sawtooth + spike).

**Guard.** Where `v <= 0` the advection correction is skipped and the converged
equilibrium ionization/temperature are kept (as first written, gated on
`pp_metal_on`, so a metals-off run was byte-identical), breaking the upwind cascade.
§5 makes the same energy solve robust by construction (Brent).

*Since generalized: the validity test in `post_process_adv.f90` is now three
conditions and none of them is gated on `pp_metal_on` -- (i) inflow `v <= 0` on
either face, (ii) a Damkohler number `Da = (dr/v)*nu_relax` above
`Da_local_equilibrium`, with `nu_relax` the SLOWEST relaxation rate among the
species the advection system solves (the He 2^3S row, `A31 = 1.27e-4` 1/s, usually
sets it; gating on the hydrogen rate alone froze the metastable at equilibrium in
cells where it is in fact advected), and (iii) an equilibrium ion fraction below
`xHII_adv_min`, where the solver's absolute resolution on `x_HI` makes the
extracted ion density meaningless. Details:
`docs/postprocess_advection_validity.md`.*

---

## 5. Brent for the scalar energy equation (Task 1, 2026-06-07)

`T_equation` (the post-process energy balance in each cell) is a single nonlinear equation
in `x = T/T0`. With metal cooling its residual is non-monotone and has a spurious
*hot* root that the general Newton/Powell solver (`hybrd1`) can land on (the root
cause of the §4 spike). A bracketing solver that selects the lowest (physical) root
removes the failure mode structurally. The scalar T solve runs only in the single
post-process pass, so the scan cost is negligible.

**Solver** (in `T_equation.f90`, module `equation_T`): `Tres(xx, params)` (scalar
residual wrapper); `solve_T_brent(params, x_guess, x_out, ok)`: log-spaced upward
scan from `max(0.05*x_guess, 1 K)` to `4*x_guess`, take the **first sign change**
(= lowest = physical root), polish with Brent; `ok=.false.` → caller falls back to
eq T; `brent_root(...)` (standard Brent).

**Wiring** (`post_process_adv.f90`): when `pp_metal_on` (and `use_brent_tsolve`, §7)
call `solve_T_brent`; else legacy `hybrd1` + 2×-band reject; metals-off keeps the
original `hybrd1`. The §4 `v<=0` guard is retained.

**Validation.** PP-only on the converged Case B: wind cells (v>0, 244) match the
`hybrd1` baseline to **6.5e-12** (machine precision; same physical root); base
[1.02,1.10]Rp has 0 T sign-changes (smooth); log10 Ṁ 13.38 unchanged. Brent
reproduces the validated result while *structurally* selecting the physical root.

---

## 6. Analytic-Jacobian Newton for the ionization equilibrium (Task 2, 2026-06-07)

All ionization systems were solved with MINPACK `hybrd1` (finite-difference
Jacobian). Task 2 supplies an analytic Jacobian + a self-contained damped-Newton
solver, keeping `hybrd1` as a fallback (zero regression risk).

**Newton solver** (new `newton_solver.f90`): `newton_dense` (damped Newton +
backtracking Armijo line search on `||f||_2`; `J dx = -f` via Gaussian elimination
with partial pivoting, `gauss_solve`; `info=1` converged, `0` = singular/no-progress/
NaN/maxit). `solve_ieq` tries Newton, on any failure restores x and calls `hybrd1`;
counters `nt_calls`/`nt_fallback` reported by `EXHALE_main`.

**Analytic Jacobians.** `jac_system_H` (1×1), `jac_system_HeH` (3×3) direct.
`jac_system_HeH_metals` built structurally as
`J(i,k) = [photoion + dC_i/dx_local * n_e]_local + C_i * dn_e/dx_k` (block-diagonal
"direct" part + rank-1 n_e coupling). Charge exchange added by `cx_add_to_jac` (new
in `charge_exchange.f90`), the exact derivative mirror of `cx_add_to_fvec` (each
rate kc*D*A is bilinear in two reactant densities, both linear in the unknowns).
Absent elements and the unused upper stage of two-stage elements are pinned to
identity rows. `System_HeH_TR` (He triplet) stays on `hybrd1` (no Jacobian; at the
time of this change mutually exclusive with metals and unused here, the merged
`System_HeH_TR_metals` came later).

**Wiring/build.** `ionization_equilibrium.f90`: the H-only, H/He, and H/He+metals
`hybrd1` calls go through `solve_ieq(..., jac_*, ...)`; triplet unchanged.
`newton_solver.f90` added to Makefile. Builds clean (gfortran).

**Validation.** PP-only sweep over converged Case B (504 cells):
- **Newton 504/504 solves (100%, zero fallback)**: analytic Jacobian (H/He + 10
  metals + charge exchange) converges on every cell.
- **A/B at identical state** (Newton vs forced-`hybrd1` via `EXHALE_FORCE_HYBRD1=1`,
  so the one-step hydro drift cancels): `Ion_species` agrees to **7.2e-7** over all
  stages, T to 3.5e-4 → Jacobian correct.
- log10 Ṁ 13.38 unchanged.
- A naive Newton-vs-baseline `Ion_species` compare showed 63% on trace neutral C I:
  the 0.26% one-step T drift amplified, not the solver (confirmed by the 7.2e-7
  A/B agreement).

---

## 7. Input switches (new default, legacy hybrd1 optional)

Both upgrades are the default and can be turned off per run via optional `input.inp`
keywords (anywhere in the optional block):

| keyword | default | effect when `False` |
| :-- | :-- | :-- |
| `Newton solver: True/False` | `True` | ionization equilibrium uses legacy MINPACK `hybrd1` |
| `Brent solver: True/False`  | `True` | energy solve uses legacy `hybrd1` + 2×-band reject |

Absent keywords keep the new solvers. Flags `use_newton_ieq` / `use_brent_tsolve`
(global_parameters), parsed in `input_read.f90`. The `EXHALE_FORCE_HYBRD1=1`
environment variable also forces the `hybrd1` ionization path for validation
without editing the input.

---

## 8. Speed

Full-run wall-clock (hydro + radiation + ionization) over the converged Case B with
`Load IC` + `Force start` (~2025 steps, identical trajectory), `OMP_NUM_THREADS=4`:

| solver | steps | wall (run 1 / 2) | ms/step |
| :-- | :-: | :-- | :-: |
| Newton | 2023 | 54.7 / 54.3 s | 26.9 |
| hybrd1 | 2027 | 67.8 / 68.9 s | 33.7 |

The analytic-Jacobian Newton is **~25% faster per step (×1.25)** than legacy
`hybrd1`; step counts match to 0.2% (fair step-for-step comparison). The `ioniz_eq` cell
loop is **serial**, so the relative speedup grows with thread count (hydro/radiation
parallelize, the solve does not).
*2026-08-15 note: no longer serial, the `ioniz_eq` cell sweep was parallelized in
2026-06 (`!$omp parallel do` over the cell loop in
`src/modules/radiation/ionization_equilibrium.f90`; see
`docs/openmp_parallelization.md`), so the "grows with thread count" argument no
longer holds. The measured ×1.25 speedup per step above is unaffected: it was taken
step-for-step at fixed `OMP_NUM_THREADS=4`.*
This is on top of the ~2× hydro-loop speedup from
the semi-implicit energy solver (§2). Raw numbers in
`WASP-121b/solver_validation/timing.txt`.

---

## 9. Validation data (persisted) and reproduction

The Task 1/2 comparison data lives in `WASP-121b/solver_validation/` (kept):
- `comparison.txt`: results table (sections A/B/C).
- `newton/`, `hybrd1/`: the two PP-only runs (Newton vs the `hybrd1` reference)
  over the converged Case B, each with full `output/` and `run.log`.
- `baseline/`: a copy of the converged Case B (`WASP-121b/output/`).
- `run_validation.sh`, `compare.py`, regenerate everything:
  `cd WASP-121b/solver_validation && bash run_validation.sh`.
- `solver_comparison.ipynb` (built by `build_solver_nb.py`), plots: (1) T/n/v and
  key ion densities for baseline/hybrd1/Newton overlaid (coincide); (2) the pure
  solver difference Newton−hybrd1 vs radius (~1e-6-1e-7 over [1.05,2]Rp); (3) the
  Brent `_adv` base temperature (smooth; machine-precision match in the wind).
- `timing.txt`: the §8 speed table.

Latest `comparison.txt`:
```
A. Newton vs hybrd1  (identical state -> pure SOLVER difference)
     eq T          max rel diff : 3.499e-04
     Ion_species   max rel diff : 7.186e-07   (worst: FeI)
B. Newton vs baseline (converged; includes one-step PP hydro drift)
     eq T          max rel diff : 2.578e-03
     Ion_species   max rel diff : 6.266e-01   (worst: CI, drift-amplified)
C. Brent _adv T vs baseline _adv (Task 1, scalar energy equation)
     wind (v>0)    max rel diff : 6.504e-12
     all cells     max rel diff : 3.499e-04
```
Section A is the meaningful solver test (identical state, so the one-step drift
cancels); the 63% on the trace neutral C I in section B is that 0.26% T drift
amplified, not the solver.

---

## 10. Steady-state Newton--Krylov solver and convergence overhaul (2026-06-10/11)

The largest solver change to date; summarized here for completeness of this
log and **documented in full in `steady_solver_memo.pdf`** (development
record, validation, hardening), with user-facing usage in the user manual
Sects. 2.4--2.5.

**Why.** The marching `du` stop (relative spread of rho*v*r^2) is provably
blind to the operator-split energy imbalance: states satisfying `du < 1e-3`
can carry steady residuals ||R||_E ~ 30 ("premature dip"), shifting Mdot by
~20% (WASP-121b: log10 Mdot 13.631 premature vs. 13.711 converged). The fix
is to measure convergence by the *steady residual* R of the discrete
equations and, near the fixed point, to solve F(Y)=0 directly.

**What was added** (all opt-in via `input.inp` keywords; every legacy
default is byte-identical, regression-gated):

- *Two-stage marching* (enabled by `Reconstruction scheme: PLM+WENO3`, with
  thresholds from `du_th [PLM,WENO3]:`): PLM warm-up, switch to WENO3 at the
  first threshold, stop at the second. `Reconstruction scheme: PLM` or `WENO3`
  alone is single-stage and uses only the first threshold.
- *Residual monitor and residual-based stop* (`Resid tol:`, `Level tol:`):
  periodic max ||R||_inf over mass/momentum/energy, independent of `du`.
- *Preconditioned JFNK steady solver*
  (`src/modules/time_step/steady_newton.f90`): matrix-free Newton--Krylov
  with right-preconditioned GMRES, colored-FD banded preconditioner
  (LAPACK dgbtrf/dgbtrs, kl=ku=8), PTC (I/dtau + J) with SER ramp, and the
  five ingredients each found via a localized stall: smooth ||D^-1 F||_2
  merit, diagonal scaling, smooth base valve (`Valve eps:`), frozen WENO
  weights, non-monotone (Grippo) line search.
  *2026-08-15 note: two of the five were replaced on 2026-08-10/11. (i) The
  diagonal scaling is no longer the global `build_scaling` with its
  `1e-6*max_j|rho v|` momentum floor but `cell_state_scales`
  (`steady_newton.f90`), which builds every row scale from the cell's own
  state (rho, rho(|v|+c_s), E). (ii) The WENO weights are still frozen for
  the Newton model (`weno_mode = 1/2` around the Jacobian build), but the
  line-search trials are evaluated with `weno_mode = 0`, i.e. against the
  true residual the solve is driving to zero; see
  `docs/newton_scaling_and_base_wall.md` §§3 and 10.*
- *Production wiring* (`Solver: Newton`): marching warm-up with du-stops
  suspended until the monitored ||R|| < 5e-2, then the JFNK finish to
  `Resid tol` (1e-3), then standard outputs and post-processing. Validated
  end-to-end from a cold IC (WASP-121b: identical Mdot to the warm-started
  reference; HD189733b: 33 s where marching spent >3 h).
  *2026-08-15 note: the hand-off criterion is the flux metric `du`, not
  ||R||. `EXHALE_main.f90` hands over once the warm-up has flattened the wind
  to `du < newton_du_switch`, default `1.0e-2` (`parameters.f90`), with the
  trigger armed on a descending crossing and a plateau escape at
  `5*newton_du_switch` after `N_stall` steps.*
- *Failure containment* (2026-06-11): best-iterate tracking, fail-fast
  after 15 outer iterations without a new best residual (`info=2`); on
  failure `EXHALE_main` keeps the best iterate, disables the Newton mode, and
  resumes plain marching --- an unconverged Newton state is never accepted
  as the final answer.
  *2026-08-15 note: the watchdog was replaced on 2026-08-10. It now fires on
  `n_no_descent_max = 12` consecutive failed line searches
  (`steady_newton.f90`), not on 15 outer iterations without a new best
  residual; see `docs/newton_scaling_and_base_wall.md` §4. The rest of the
  containment (best-iterate tracking, fall back to marching) is unchanged.*

Ready-to-run configurations for every solver combination (legacy marching,
two-stage, Newton, Newton-from-state, warm-seed IC, ...) live under
`examples/01_legacy_marching/` ... `examples/16_molecular_metals/`
(user manual, Table 4).
