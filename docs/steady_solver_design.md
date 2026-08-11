# Steady-state solver design — PTC + banded Newton (Path B, step (ii))

Date: 2026-06-10. Goal: solve the finite-volume steady residual `F(U) = 0`
directly, so the true wind is reached without the tens-of-thousands of
marching steps the residual study showed are needed (and which `du`/`dtu`
stops miss entirely — the premature WASP golden had `R_energy ~ 30`).

## 1. What we already have (building blocks)

- `eval` of the steady residual is proven (the `EXHALE_RESIDUAL=1` hook and the
  in-loop monitor): `R(:,1)=dF-S`, `R(:,2)=dF-S`, `R(:,3)=dF_E-S_E-(heat-cool)`,
  assembled from `Reconstruct` + `RK_rhs` (hydro) and `ioniz_eq` (heat/cool via
  local ionization-equilibrium elimination). HLLC ignores `alpha`.
- LAPACK is available (`-llapack`, ATLAS 3.7.1; `dgbtrf`/`dgbtrs` link-tested).
- Deterministic single-thread regression harness + the criterion-independent
  `||R||` measure to validate the result against.

## 2. Unknowns and residual

- Unknown vector `Y` = hydro conserved variables on physical cells:
  `Y = [u(j,1), u(j,2), u(j,3)]` for `j = 1..N` → `neq = 3*N` (~1500 for N=500).
  Ionization/temperature are NOT unknowns — they are eliminated locally inside
  each residual evaluation (`ioniz_eq`), the "method A / local elimination"
  choice from `ATES_sundials_solver_plan.md` (smaller system, closest to the
  current code; revisit monolithic coupling only if Newton struggles).
- `F(Y)`: unpack `Y → u` (physical cells); `Apply_BC` fills ghosts; then the
  existing evaluation pipeline produces `R`; pack `R` over physical cells.

## 3. Pseudo-transient continuation (PTC)

Pure Newton from a cold isothermal IC will not converge. Wrap it:

    G(Y) = (Y - Y_old)/dtau + F(Y) = 0

- Small `dtau` ≈ explicit marching (robust startup); large `dtau` → Newton.
- Ramp `dtau` on success (e.g. SER: `dtau <- dtau * ||F_prev|| / ||F||`),
  cut on a failed/again-rising Newton step.
- Start `dtau` near the CFL step; `Y_old` = previous PTC iterate.

## 4. Jacobian (banded, LAPACK)

1-D finite volume → `J = dG/dY` is block-tridiagonal: cell `j` couples to
`j-1, j, j+1` through the reconstruction/flux stencil (WENO3 is 2-cell each
side → pentadiagonal in cells; with 3 components that is bandwidth
`ku = kl = 3*2 = 6`, store as `(2*kl+ku+1, neq)` for `dgbtrf`).

- **First cut: finite-difference banded Jacobian with graph coloring.** A
  banded matrix of half-bandwidth `b` needs only `b+1` residual evaluations
  (perturb every `(b+1)`-th unknown simultaneously), i.e. ~13–19 `F` evals per
  Jacobian — cheap at this size. Factor with `dgbtrf`, solve with `dgbtrs`.
- The radiation column density makes `F` weakly non-local (one-sided
  cumulative). Ignore that coupling in the banded `J` (treat as a
  preconditioner); the PTC `1/dtau` diagonal + line search absorb the error.
  If convergence stalls, switch to a matrix-free Newton–Krylov (GMRES) using
  the banded `J` as preconditioner.

## 5. Positivity / robustness

- Line search / step damping so `rho, p > floor`; reject steps that produce
  NaN or negative pressure and cut `dtau`.
- Scale residuals/unknowns by reference values (`rho_bc`, `p0`-equivalents)
  so the Newton system is well-conditioned.
- Lower-boundary: keep the existing `Apply_BC` (ghost fill) inside `F`; treat
  only interior cells as unknowns initially.

## 6. Integration & validation

- New module `src/modules/time_step/steady_newton.f90`; runtime opt-in
  `Solver: Newton` (default `Marching`, byte-identical preserved).
- Reuse the `||R||` monitor for the stop test and to compare the Newton
  solution against the long-marching `cold-35k` reference (must agree, and
  reach a much smaller `||R||`, in far fewer residual evaluations).
- Gate: on WASP-121b, Newton (from the cold IC, or warm-started from a cheap
  marching pass) reaches `||R|| << 1e-2` and `Mdot ≈ 13.71` in
  O(10²) residual evals, vs >35k marching steps; then HD189733b.

## 7. Incremental plan (each step validated)

1. **[DONE]** Factor the residual evaluation into a shared module
   `steady_residual.f90` (`assemble_residual`, `residual_norms`), used by the
   diagnostic hook, the in-loop monitor, and Newton. Byte-identical re-confirmed;
   `EXHALE_RESIDUAL` numbers unchanged.
2. **[DONE]** `steady_newton.f90`: `pack_U`/`unpack_U` (unknowns = hydro on
   physical cells `1..N`, `neq=3*N`; ghosts set by `Apply_BC`) + `newton_residual
   (Y, f_sp, F)` (local ionization elimination). Validated via `EXHALE_NEWTON_TEST=1`:
   pack/unpack are exact inverses, and `F(Y)` reproduces the `EXHALE_RESIDUAL`
   diagnostic to all digits. KEY FIX: the residual must `Apply_BC` BEFORE the
   ioniz_eq/heat-cool step so `F(Y)` depends only on the interior `Y` (ghosts =
   BC); the diagnostic hook was reordered to match (this shifted the consistent
   cold-35k energy residual from 1.12e-2 to 1.52e-2; the golden >> cold ranking
   is unchanged).
3. **[NEXT]** FD banded Jacobian builder (`dgbtrf`/`dgbtrs`; half-bw 6 for
   WENO3 3-var; graph coloring → ~13-19 `F` evals per Jacobian); unit-check
   against directional finite differences of `F`.
4. **[PARTIAL]** PTC driver (`solve_steady_ptc`): `(I/dtau + J) dY = -F` via
   LAPACK `dgbtrf`/`dgbtrs`, backtracking line search, positivity, SER `dtau`
   ramp; `EXHALE_PTC=1` hook. Status on WASP-121b:
   - Cold IC: stalls (expected — `du`/thermal transient too far for Newton).
   - **KEY FIX**: the line-search merit must be the smooth `||F||_2`, not the
     component-wise max-relative `rnorm` (non-smooth → rejected all steps).
   - With `||F||_2` merit + full-residual banded Jacobian
     (`build_banded_jac_full`, includes local `d(heat-cool)/dE`): warm-started
     from the cold-35k state it now DESCENDS (`||F||_2` 0.28 -> 0.11) at
     `lam=1`, but only LINEARLY (~2%/iter), `dtau` stays ~1e-4 (won't grow),
     and the max-relative `||R||` drifts UP (2.1e-2 -> 5.7e-2). I.e. the banded
     Jacobian is too approximate: it omits the non-local radiation
     (column-density) coupling and suffers colored-FD contamination, so the
     Newton step is only marginally better than relaxation and biases the
     worst cell.
   - Conclusion: a banded preconditioner alone is insufficient here. NEXT is
     step 5 (JFNK), not more PTC tuning.
5. **[NEXT]** Matrix-free Newton-Krylov (JFNK): exact J*v via a directional
   FD of the FULL residual (captures the non-local radiation), GMRES with the
   banded `build_banded_jac_full` as right-preconditioner, inside the same PTC
   `dtau` continuation. Add variable scaling (density/momentum/energy and
   residual) — the un-scaled system is part of why `dtau` won't grow. Validate
   vs the marching reference (~13.71), then HD189733b.

## 8. JFNK (step 5) — built, and the blocker localized (2026-06-11)

Implemented matrix-free Newton-Krylov (`jv_product` = directional FD of the
FULL residual; `pgmres` = right-preconditioned GMRES(m) with the banded
`build_banded_jac_full` factorization as preconditioner; `solve_steady_jfnk`
= PTC + GMRES + `||F||_2` line search). Hook: `EXHALE_PTC=1 EXHALE_PTC_JFNK=1`
(`EXHALE_PTC_DTAU0=<v>` overrides dtau0). LAPACK linked (`LDLIBS=-llapack`).

Results on WASP-121b (warm start from the cold-35k state, dtau0 = 1):
- GMRES works: inner solve converges in ~9 Krylov vectors.
- JFNK clearly beats banded-only PTC: `||F||_2` drops 0.281 -> 0.106 -> 0.062
  in TWO Newton steps (banded-only managed ~2%/step), dtau grows to ~2.7.
- BUT both methods then STALL (line search `lam -> 0`, `||F||_2` frozen).
  **The stall is localized**: the worst relative residual is always the
  ENERGY equation at cell `j = j_min`, i.e. `r = 1.50 = r_esc` (the escape
  radius / inner edge of the supersonic wind window). It is NOT the base
  (base-breathing is not the issue here).

~~Initial interpretation: sonic-point non-smoothness at `r_esc`.~~
**RETRACTED (2026-06-11)** — that localization was an artifact of searching
the worst cell only inside `[j_min:N]`: with the search widened to the whole
domain the worst residual is the MOMENTUM at the FIRST cells (`j=1`,
`r=1.000`), i.e. the lower-boundary region, and the actual sonic point of the
state is at `r=1.61`, not `1.50`. Also `r_esc` enters no physics (only the
convergence window), and the `Rate/2` 2D factor is a global constant.

## 9. Reduction tests (2026-06-11) — base-freeze does NOT unblock; scaling is
the prime suspect

- Frozen-base experiment (`EXHALE_PTC_NFIX=n`: anchor rows `F_j = Y_j - Y^fix_j`
  for the first n cells): the worst cell simply FOLLOWS the slab edge
  (n=0 -> j=1; n=2 -> j=5; n=60 -> j=61, r=1.012, momentum), and with any
  freezing the very FIRST Newton step already fails its line search (whereas
  un-frozen JFNK took two good steps first). The anchor interface itself
  introduces a kink; base-freezing is not the fix.
- Convergent evidence for a SCALING problem: the unknowns span ~7 decades in
  code units (base rho ~1 vs wind ~1e-7). The FD steps use absolute floors
  (`max(|Y_j|,1)*sqrt(eps)` per column; global eps in `jv_product`), so wind
  cells receive 15-60% RELATIVE perturbations -> the Jacobian columns and
  J*v products for wind rows are inaccurate. This matches the observed
  behavior precisely: the first 1-2 Newton steps (base-scale corrections,
  where FD is accurate) are excellent, then the residual moves to small-scale
  cells and no computed direction descends.

Next concrete step: diagonal scaling D applied to each unknown (e.g.
`D_i = max(|Y_i|, c_k * max_j|u(:,k)|)` with a floor c_k ~ 1e-6 for each component):
solve the scaled system `D^-1 J D z = -D^-1 F`, take FD steps relative to D,
and use `||D^-1 F||_2` as the line-search merit. Then re-run the JFNK probes
(and only afterwards revisit base-BC smoothness if a residual blocker remains).

## 10. Scaled JFNK + smooth-valve probes (2026-06-11)

Implemented: `build_scaling` (D_i = max(|Y_i|, 1e-6*max|Y_k|)); the banded
Jacobian's FD column steps are now D-relative; `pgmres` operates on the scaled
system `D^-1 (I/dtau + J) D` with the scaled banded factorization as
preconditioner; the line-search merit is `||D^-1 F||_2`. Also added the
opt-in smooth base valve (`Valve eps: <v_eps>`, code units): ghost velocity
`0.5*(v + sqrt(v^2 + eps^2))` instead of `max(v,0)` (differentiable; -> 0 for
v << -eps, -> v for v >> eps; default off = exact legacy valve).

Probe results (warm start from cold-35k, dtau0 = 1):
- **Scaling is a major win**: GMRES now needs 2-6 Krylov vectors (was 9);
  `||Fs||_2` drops 29 -> 1.2 (24x) before stalling — vs 2.6x unscaled.
- The stall ALWAYS lands on the base cells (j = 1-2, momentum/energy), with
  the state hugging the valve kink: v(1) = +4.5e-3 c_s at the warm state.
- Smooth valve, eps = 1e-2 (too large, ~v(1) itself): perturbs the start
  state (||Fs|| 29 -> 878) but the solver then descends ~10x with healthy
  line searches THROUGH the base — evidence that differentiability there is
  what Newton needs. eps = 1e-4 (minimal perturbation): 4.6x descent
  (29 -> 6.3), then stalls at base momentum AGAIN.
- Conclusion so far: the base-cell residual contains FURTHER non-smoothness
  beyond the valve. Prime remaining suspect: the **WENO3 nonlinear weights**
  (smoothness indicators) — a classic Newton blocker for FV steady solves —
  plus possibly the HLLC wave-speed switches and ioniz_eq solver noise.

Prioritized next experiments:
1. **Frozen-weights Newton**: evaluate the residual for the Jacobian/GMRES
   with the WENO3 nonlinear weights FROZEN at the outer iterate (or solve
   with PLM/linear reconstruction inside Newton, exact WENO3 outside) — the
   standard remedy in steady FV solvers.
2. If base cells still block: solve with the base pinned by a PHYSICAL
   boundary condition rather than anchor rows (e.g. prescribe the cell-1
   state from hydrostatic balance, removing those rows' kinks).
3. Wire the converged path into a `Solver: Newton` runtime option once a
   full-residual convergence is demonstrated.

## 11. CONVERGED (2026-06-11) — frozen weights + non-monotone line search

Implemented (a) `weno_mode` (global_parameters; Reconstruction stores the
ESWENO3 smoothness factors S0/S1 on a mode-1 pass and reuses them on mode-2
passes; default 0 = byte-identical) — `solve_steady_jfnk` freezes the weights
at each outer iterate and runs all inner evaluations (banded FD, GMRES J*v,
line-search trials) with them frozen; and (b) a **non-monotone (Grippo,
memory 5) line search** — the monotone Armijo test was rejecting valid steps
once the required decrease fell below the iterative-chemistry noise floor of
the residual.

**Result (WASP-121b He23S+metals, warm start from cold-35k, dtau0=1, smooth
valve eps=1e-4, diagonal scaling):**
- `Resid tol 1e-3`: **converged in 8 Newton iterations** (all lam=1, GMRES
  2-3 vectors) — roughly ~200 residual evaluations, i.e. ~200 marching-steps
  of cost, vs >35,000 marching steps that never got this far.
- `Resid tol 1e-4`: converged in 22 iterations to `||R|| = 6.0e-5` —
  territory marching cannot reach at all. Convergence is linear (Picard
  refresh of the frozen weights + chemistry noise floor), but every step is
  a full Newton step.
- **Physical validation**: raw `log10 4*pi*rho*v*r^2(r_max) = 13.7111`, the
  flattest state produced to date (`rho*v*r^2` spread for r >= 1.5 of
  **5.0e-4**, vs 1.2e-3 at marching's 35k-step state and 1.8e-3 at the old
  premature golden), `T_max = 10,747 K` (0.2% from the marching reference).
  The Newton answer sits exactly on the marching asymptote (~13.71).

The full ingredient list that was necessary (each verified by its own probe):
smooth `||F||_2`-type merit -> **diagonal scaling** -> **smooth base valve**
-> **frozen WENO weights** -> **non-monotone line search**. Removing any one
of these reproduces a stall that was observed and documented along the way.

Remaining hardening (next):
- Cold-start strategy (current solver needs a warm start; a short marching
  warm-up -> Newton finish is the natural hybrid, or a dtau ramp from small).
- `Solver: Newton` runtime wiring + valve-eps sensitivity study (does the
  answer depend on eps = 1e-4 vs 1e-5?).
- HD189733b application — the real prize (its marching took >1e6 steps).

Everything remains opt-in; the marching default and all regressions are
byte-identical (re-confirmed after each increment, including weno_mode=0 and
the valve-off default).

## 12. Scaling replaced; the "base blocker" of §9-§10 reassessed (2026-08-10)

`build_scaling` (`D_i = max(|Y_i|, 1e-6 max_j |Y_k|)`), introduced in §10, is
replaced by `cell_state_scales`, which builds every scale from the cell's own
state: `rho`, `rho(|v| + c_s)`, `E`. Everything else in §10-§11 stands — the
scaled system, the D-relative FD steps, the merit `||D^-1 F||_2`, the smooth
valve, the frozen WENO weights, the non-monotone line search.

The old momentum floor `1e-6 max_j |rho v|` was set by the base cell, which
holds the global maximum of `|rho v|` while carrying no wind, so cells where
the flow reverses were scaled by a number unrelated to their own state. The
repeated observation in §9-§10 that "the stall ALWAYS lands on the base cells
(j = 1-2)" was in large part that floor plus a worst-cell diagnostic that also
normalized by `max_j |u(k,j)|`; the base momentum row itself turns out to be
satisfiable to within a 2.9 ppm ghost-pressure change. The stagnation watchdog
was also changed, from "no new best `||R||` in 15 iterations" to a count of
consecutive failed line searches. Measurements and the cases this unblocks:
`docs/newton_scaling_and_base_wall.md`.
