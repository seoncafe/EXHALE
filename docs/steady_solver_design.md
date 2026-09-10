# Steady-state solver design: PTC + banded Newton (Path B, step (ii))

Date: 2026-06-10. Goal: solve the finite-volume steady residual `F(U) = 0`
directly, so the true wind is reached without the tens-of-thousands of
marching steps the residual study showed are needed (and which `du`/`dtu`
stops miss entirely, the premature WASP golden had `R_energy ~ 30`).

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
  Ionization/temperature are NOT unknowns: they are eliminated locally inside
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
  (perturb every `(b+1)`-th unknown simultaneously), i.e. ~13-19 `F` evals per
  Jacobian: cheap at this size. Factor with `dgbtrf`, solve with `dgbtrs`.
- The radiation column density makes `F` weakly non-local (one-sided
  cumulative). Ignore that coupling in the banded `J` (treat as a
  preconditioner); the PTC `1/dtau` diagonal + line search absorb the error.
  If convergence stalls, switch to a matrix-free Newton-Krylov (GMRES) using
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
   - Cold IC: stalls (expected, `du`/thermal transient too far for Newton).
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
   residual): the un-scaled system is part of why `dtau` won't grow. Validate
   vs the marching reference (~13.71), then HD189733b.

## 8. JFNK (step 5): built, and the blocker localized (2026-06-11)

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
**RETRACTED (2026-06-11)**, that localization was an artifact of searching
the worst cell only inside `[j_min:N]`: with the search widened to the whole
domain the worst residual is the MOMENTUM at the FIRST cells (`j=1`,
`r=1.000`), i.e. the lower-boundary region, and the actual sonic point of the
state is at `r=1.61`, not `1.50`. Also `r_esc` enters no physics (only the
convergence window), and the `Rate/2` 2D factor is a global constant.

## 9. Reduction tests (2026-06-11): base-freeze does NOT unblock; scaling is
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
  `||Fs||_2` drops 29 -> 1.2 (24x) before stalling: vs 2.6x unscaled.
- The stall ALWAYS lands on the base cells (j = 1-2, momentum/energy), with
  the state hugging the valve kink: v(1) = +4.5e-3 c_s at the warm state.
- Smooth valve, eps = 1e-2 (too large, ~v(1) itself): perturbs the start
  state (||Fs|| 29 -> 878) but the solver then descends ~10x with healthy
  line searches THROUGH the base, evidence that differentiability there is
  what Newton needs. eps = 1e-4 (minimal perturbation): 4.6x descent
  (29 -> 6.3), then stalls at base momentum AGAIN.
- Conclusion so far: the base-cell residual contains FURTHER non-smoothness
  beyond the valve. Prime remaining suspect: the **WENO3 nonlinear weights**
  (smoothness indicators), a classic Newton blocker for FV steady solves,
  plus possibly the HLLC wave-speed switches and ioniz_eq solver noise.

Prioritized next experiments:
1. **Frozen-weights Newton**: evaluate the residual for the Jacobian/GMRES
   with the WENO3 nonlinear weights FROZEN at the outer iterate (or solve
   with PLM/linear reconstruction inside Newton, exact WENO3 outside), the
   standard remedy in steady FV solvers.
2. If base cells still block: solve with the base pinned by a PHYSICAL
   boundary condition rather than anchor rows (e.g. prescribe the cell-1
   state from hydrostatic balance, removing those rows' kinks).
3. Wire the converged path into a `Solver: Newton` runtime option once a
   full-residual convergence is demonstrated.

## 11. CONVERGED (2026-06-11): frozen weights + non-monotone line search

Implemented (a) `weno_mode` (global_parameters; Reconstruction stores the
ESWENO3 smoothness factors S0/S1 on a mode-1 pass and reuses them on mode-2
passes; default 0 = byte-identical): `solve_steady_jfnk` freezes the weights
at each outer iterate and runs all inner evaluations (banded FD, GMRES J*v)
with them frozen; and (b) a **non-monotone (Grippo, memory 5) line search**:
the monotone Armijo test was rejecting valid steps once the required decrease
fell below the iterative-chemistry noise floor of the residual.

*[2026-08-15 note: two corrections to the paragraph above. (i) The line-search
trials were dropped from the frozen-evaluation list, since 2026-08-11 the
trials run with `weno_mode = 0`, i.e. the smoothness weights are recomputed at
the trial state, so acceptance is decided on the true residual the solve is
driving to zero (`steady_newton.f90`, the `weno_mode = 0` set just before the
`do ls = 1, 20` search). (ii) The "iterative-chemistry noise floor"
justification for the non-monotone window was measured and refuted: four repeat
residual evaluations at fixed `Y`, resetting `f_sp` as the line search does,
reproduced the residual bitwise, so there is no such floor. The window is kept,
but for a different measured reason. Both in
`docs/newton_scaling_and_base_wall.md` §10.]*

**Result (WASP-121b He23S+metals, warm start from cold-35k, dtau0=1, smooth
valve eps=1e-4, diagonal scaling):**
- `Resid tol 1e-3`: **converged in 8 Newton iterations** (all lam=1, GMRES
  2-3 vectors): roughly ~200 residual evaluations, i.e. ~200 marching-steps
  of cost, vs >35,000 marching steps that never got this far.
- `Resid tol 1e-4`: converged in 22 iterations to `||R|| = 6.0e-5`,
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
- HD189733b application: the real prize (its marching took >1e6 steps).

Everything remains opt-in; the marching default and all regressions are
byte-identical (re-confirmed after each increment, including weno_mode=0 and
the valve-off default).

## 12. Scaling replaced; the "base blocker" of §9-§10 reassessed (2026-08-10)

`build_scaling` (`D_i = max(|Y_i|, 1e-6 max_j |Y_k|)`), introduced in §10, is
replaced by `cell_state_scales`, which builds every scale from the cell's own
state: `rho`, `rho(|v| + c_s)`, `E`. Everything else in §10-§11 stands, the
scaled system, the D-relative FD steps, the merit `||D^-1 F||_2`, the smooth
valve, the frozen WENO weights, the non-monotone line search.
*[2026-08-15 note: "the frozen WENO weights stands" no longer covers the line
search, since 2026-08-11 the trials are evaluated with the weights recomputed
(`weno_mode = 0`); the freezing applies to the Jacobian/GMRES evaluations only.
The non-monotone window also stands on a different measured reason than the one
recorded in §11. See `docs/newton_scaling_and_base_wall.md` §10.]*

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

## 13. The unknown space is the transported set, not the hydrodynamic triple (2026-09-07, B5)

Section 2 states the unknown vector as the three hydrodynamic conserved
variables per physical cell, everything else eliminated locally. That is right
for a configuration whose only equations are the hydrodynamic ones, and it
stopped being the whole story when the code gained transported balances: a
carrier or an element that is moved by its own operator is not a local
equilibrium of its cell, and eliminating it locally solves a different system
from the one the marching path integrates.

**The unknown vector is now three hydrodynamic variables plus one unknown for
every transported balance the configuration activates.** The set is the one
`carrier_set_init` fixes once the keys are parsed (H2 always; OH, H2O and CO
under `Oxygen chemistry`; H+ under `Ionization transport`), read through
`carrier_solved`, so the stationary system and the transport operator cannot
disagree about which balances exist. The mapping row -> physics lives in one
registry, `set_transported_species_rows` in `steady_newton.f90`, and a species
slot carries the transported density in the code's own units, `n_i/n0`, for the
reason the H2 slot already gave: the Krylov step size is
`sqrt(eps)(1 + ||Y||)/||v||`, so a component thirteen decades from the
hydrodynamic ones would set that step by itself.

The species unknown enters the residual as the composition the equilibrium
sweep is handed: `ioniz_eq` already imposes a transported partition instead of
solving it (`x_h2_fix`, `x_oh_fix`, `x_h2o_fix`, `x_hp_fix`), and the only
change is that the number comes from the Newton unknown rather than from the
last transport step. The row is the stationary balance
`div(F_s + Phi_s) = P_s - L_s`, assembled by the module that owns it
(`carrier_steady_residual`) and converted by the code's own time scale so the
merit compares one clock.

### 13.1 The band geometry is derived

WENO3 reaches two cells and within a cell every variable couples, so in the
flat ordering `Y(nvar*(j-1)+k)` a column couples to rows within
`|row - col| <= 2 nvar + (nvar - 1) = 3 nvar - 1`, and the coloring stride is
`kl + ku + 1 = 6 nvar - 1`. The formula reproduces the two geometries the code
carried as constants, `nvar = 3 -> (8, 17)` and `nvar = 4 -> (11, 23)`, which
is the check that it is the geometry already in use. By configuration:
`nvar = 4` for the H2 balance, `5` with the proton, `7` with the oxygen cycle,
`8` with both. A species row itself reaches only one cell either way, so the
hydrodynamic rows still set the bandwidth and each extra variable only widens
the stride.

### 13.2 The Jacobian action is asserted on the system that is solved

The run-level self-test (`EXHALE_JAC_TEST=1`) builds its Jacobian before any
solve is entered, so the registry is empty and it measures the three-unknown
system in every configuration. `EXHALE_SPECIES_JAC_TEST=1` makes the same
comparison inside the solve, on `build_banded_jac_full` and `eval_residual`,
which is the Jacobian the solve actually uses, and it reports the species rows
separately so a mismatch confined to them cannot be averaged away by the
hydrodynamic rows. MEASURED on the hot Uranus with the H2 and H+ rows:
`||J r - dFD||_inf / ||J r||_inf` = 7.19e-07 over all rows and 4.60e-05 over
the species rows, with no unresolved color.

### 13.3 The frozen residual states no species balance

`frozen_residual` holds the radiation fixed and re-evaluates no chemistry, so
it carries no species row and writes those slots to zero rather than leaving
them as it found them. Before 2026-09-07 it left them unwritten, and a caller
that differenced two of its outputs divided the difference of two pieces of
stale memory by `eps`; MEASURED, the directional check of a two-species system
then reported a relative error of 1.0. The banded Jacobian built from it is a
hydrodynamic preconditioner whose species block is empty, which is what the
PTC route uses and why that route refuses a system with species rows.

### 13.4 `Ionization transport` with `Solver: Newton`

Supported by derivation when the proton is an unknown, refused otherwise. With
`Coupled carrier solve` the sweep is handed the proton fraction and solves the
other stages against it, while the proton's own equation is the stationary
balance the Newton drives to zero, so nothing is undone. Without that key the
proton is outside the unknown space and the last sweep does put back the local
ionization state the option exists to leave; that is refused at
`input_read.f90` with the reason at the site. `EXHALE_PTC=1` keeps its own
refusal of the combination, and the PTC route refuses a species row in any
case.

## 14. The species rows carry the marching operator's own terms (2026-09-07, B5b)

Section 13 made the transported balances rows of the stationary system. What
those rows BALANCED was still not what the marching path integrates, and item
B5b closed that.

### 14.1 The advective term is the divergence of the face species fluxes

Increment B4-1c moved the material advection of every transported species into
the Runge-Kutta stages as `F_s = F_rho Y_s^face` and removed the cell-velocity
term from the marching rows. The stationary balance and the fixed-wind
relaxation kept the cell-velocity upwind difference `n_tot v df/dr` plus a
deferred van Leer correction, so one term of one equation had two
discretizations and therefore two fixed points. The acceptance contract
requires the state a Newton converges on and the state the marching converges
on to be one object, so that is a defect and not a choice.

Both now form the term through the same three routines
(`species_face_fraction`, `species_face_flux`, `species_flux_divergence`) on
the same faces, areas and volumes as the mass row, and
`species_advective_update` calls the last of them too, so the expression
exists once. The carrier row multiplies the mass divergence by `msum/m_c` to
return it to the `cm^-3 s^-1` the row is written in; the helium row by
`n0 mu msum` to return it to `g cm^-3 s^-1`; a trace element row by
`n0 msum/A_X` and by `1/n_H` to return it to the mixing ratio per second the
row is written in.

MEASURED (`src/tests/steady_species_rows/`, a synthetic column whose face mass
flux changes sign inside it, so one cell has an inflow at both faces and no
velocity of its own): the row's advective term now agrees with the rate of
change a stage-1 species update produces on the same state to a relative
7.5e-14, which is the cancellation floor of forming a rate by differencing a
state; the cell-velocity form it replaced differs from that same update by a
relative 76.7 at worst.

### 14.2 What that changes in a row's value, and why

`div(F_rho Y) = Y div(F_rho) + F_rho grad Y`. The cell-velocity form was the
second term alone, i.e. the element or carrier equation MINUS the transported
quantity times the MASS equation. The two agree exactly where the mass row is
zero and nowhere else, so a species row now inherits the mass row's own
imbalance. That is the conservative statement and it is the one the marching
operator integrates: a state that is not stationary in mass is not stationary
in an element it carries either.

MEASURED on the 300-step relaxation snapshots of the regression cases, whose
mass rows are far from stationary (`mol_diffusion` cell 212: mass row 1.913 of
its own terms): the He/H stationary row of `mol_diffusion` moves from 5.637e-04
to 9.987e-01, and the H2 row of `mol_carrier` from 3.912e-01 to 8.358e-01. The
marching path is untouched by all of this and the profiles are byte-identical.

### 14.3 The base face needs no special case

The face composition is reconstructed from the inner ghosts wherever the face
mass flux flows inward, and those ghosts hold the handoff partition where a
handoff states one (`q_H2_base`) and the base cell's own where none does. The
Dirichlet record the cell-velocity form needed -- which carrier is pinned, at
what value, and whether the base face is inflowing -- is gone with that form:
`carrier_base_state`, `fc_base`, `base_dirichlet`, `base_inflow`,
`carrier_advection_correction`, `adv_corr`, `carrier_slope` and the advective
coefficient `ntv` no longer exist in `diffusive_photochemistry.f90`.

### 14.4 The elemental balances are rows

`element_transport_residual` returns one residual per element instead of the
worst element at each cell, and the registry carries `srow_element_he` (the
unknown is the helium mass fraction `X = rho_He/rho`) under `He_diffusion` and
one `srow_element_trace` per metal element under `He_metal_diffusion`. The
element rows are registered BEFORE the carriers and written back before them,
because an element unknown is an element TOTAL and its write-back
(`project_element_mass_fractions`, the same map the advective stages use) is a
projection of the whole species vector: a carrier written first would be
rescaled by it.

The geometry follows the same formula: `nvar = 4` for helium alone,
`4 + n_melem` with the trace elements, and the carriers add to that. The
certification reports one entry per trace element
(`elemental transport <name>`), and the completion flag of the stationary
solve reads every element row the registry carried.

### 14.5 Which key asks for the coupled system

`EXHALE_main.f90` now calls `set_transported_species_rows(carrier_in_newton)`,
the run-level `Coupled carrier solve` key alone. It used to pass
`thereis_mol .and. carrier_transport .and. carrier_in_newton`, so an
element-only configuration (`He_diffusion` with no molecular carriers)
registered no rows however the key was set. WHICH balances become rows is the
registry's own question, answered from the configuration; what the key asks is
only whether the transported balances are solved WITH the wind instead of
alternated with it, which is one choice for all of them. The key's NAME is
still carrier-specific and renaming it is a user-visible change that was not
taken here.

## 15. The elemental relaxation is that same operator (2026-09-07, B5c)

Section 14 left one advective term outside the single discretization. The
Runge-Kutta stages, the stationary carrier rows, the stationary element rows
and the fixed-wind CARRIER relaxation all formed the term as the divergence of
the face species fluxes; `element_diffusion_step` and
`relax_element_composition` still advected the elements as `rho v dX/dr` on a
SMOOTHED steady mass flux `mdot/(4 pi r^2)`, with `mdot` the median of
`r^2 rho v` over the escape window. The elemental Picard alternation of the
direct-steady outer pass -- solve the wind at a fixed composition, relax the
composition at a fixed wind -- therefore ran one operator while the row that
judged its result was another, and a converged alternation could only land on
a state the row does not read as stationary.

### 15.1 The recorded rationale for the smoothed flux does not survive

What the smoothing was for is written at the code site and it is true as far
as it goes: `rho (dX/dt + v dX/dr) = -div(r^2 J)/r^2` is the conservative
equation only where `rho` and `v` satisfy continuity, and below ~1.02 R_p the
converged states do not (the spread of `r^2 rho v` there is 10^2 to 10^4 times
its own median, because the base carries a standing sound wave whose sign
alternates from cell to cell). Relaxed on that field the NON-CONSERVATIVE form
converges to the composition of a flow that neither conserves mass nor exists.

The repair was needed because the form was non-conservative, and the
conservative form does not need it. `div(F_rho X) = X div(F_rho)
+ F_rho grad X`: the cell-velocity form was the second term alone, i.e. the
element equation minus `X` times the mass equation, so it had the mass row's
own error subtracted out of it and a base that does not conserve mass had to
be hidden from it. The divergence form carries that error instead of removing
it, and a cell can lose only the fraction of its mass the mass row loses. The
face mass flux of the current state is available from the same Riemann
evaluation the hydrodynamic rows use (`face_mass_flux_of_state`), which is
what section 14 already gave the carriers and the stationary element rows.

`element_diffusion_step` therefore takes `Frho_in`, the face mass flux, in
place of the advecting momentum density `rhov_in`, and forms the term through
`element_advective_divergence` -- `species_face_fraction`,
`species_face_flux`, `species_flux_divergence`, the same three routines on the
same faces, areas and volumes as the mass row. `relax_element_composition`
is GIVEN the flux, once, for the whole relaxation, exactly as the wind it
belongs to is read once: the caller obtains it with
`face_mass_flux_of_state` and passes it in, so the element module carries no
dependence on the steady residual (item DIFT-LINK; a `functions/` module
reaching up into `time_step/` dragged the whole tree plus LAPACK into every
program linking it, and `diffusion_tests.x` stopped linking). The trace-element arm takes the divergence of the same face
element mass fluxes, converted by `n0 msum/A_X` and `1/n_H` into the mixing
ratio per second its row is written in, and its solve became a deferred
correction: the matrix carries the diffusive half and the donor-cell part of
the advective term, the right-hand side carries what the reconstruction and
the bounding add on top of the donor cell at the current iterate, so the fixed
point of the iteration is the FULL operator. With no advective term nothing is
deferred, the system is linear and the loop takes one pass, which is the
direct solve it was before.

`element_transport_residual` now takes both halves from
`composition_residual` and the trace half from `trace_composition_residual`,
the same two routines the step solves, at a step long enough that the time
term is absent. **The relaxation's fixed point is the row's zero by
construction and not by agreement.**

### 15.2 What is asserted, and what was measured

`src/tests/steady_species_rows/` gained
`elemental_relaxation_fixed_point_is_the_row`: a synthetic 60-cell column with
a settling helium profile is relaxed at a GIVEN face mass flux until the
helium mass fraction stops moving, and the stationary elemental row of that
same flux is then measured on the relaxed state. MEASURED: the row reads
**8.65e-13** of its own terms, and the advective term is **0.993** of those
terms (the second row of the pair, which is what says the balance is not the
diffusive one on a column the wind never touched).

The flux is one with `r^2 F_rho` constant, and that is not a free choice: on a
flux that does not conserve mass, `X div(F_rho)` is a source of the element
with no sink and the conservative balance has no bounded steady composition at
all. MEASURED on the same column with a sign-changing flux: the relaxation
drives X onto 0 and 1 and stays there. That is a property of the equation, and
it is the reason the object this row is measured on is the wind the outer pass
hands over.

The marching path does not move: `element_diffusion_step` is called there
without a flux and carries the diffusive half alone, and `mol_diffusion` and
`lower_profile` (the latter with `He_metal_diffusion`, so the trace arm runs
too) are byte-identical at 300 steps before and after.

## 16. The banded model an element row can be solved on (2026-09-08, B5d)

Sections 14 and 15 made the element balances rows of the stationary system and
put one operator behind them. What they did not produce was a solve: on
`mol_diffusion` reloaded from its own 300-step snapshot with
`Coupled carrier solve: True`, the four-unknown solve left `||R||` at 1.974
and aborted at outer iteration 8 with `no descent direction exists for the
banded model at this state`, while the THREE-unknown solve of the same state
falls to 1.352e-06 in 13 iterations. That pair is the experiment this section
is about: same state, same physics, same solver, one row and one column
different.

### 16.1 The model is the defect, not the equation

MEASURED on that reload, one thread, the same start:

| route | `\|\|F/Drow\|\|_2` start -> end | `\|\|R\|\|` start -> end | how it ended |
|---|---|---|---|
| default (Krylov + the damped Gauss-Newton escape) | 1.44e+02 -> 1.28e+02 | 1.986 -> 1.974 | `no descent direction exists for the banded model`, iteration 8 |
| `EXHALE_TRUST_REGION=1` (the escape skipped, dogleg globalization) | 1.44e+02 -> 5.87e-01 | 1.986 -> 1.528 | radius collapsed to 3.8e-16, 12 model refusals by the ray test, iteration 24 |

The true merit falls by 2.4 decades under the second route on the state the
first calls stationary-or-nothing, so the state is not the obstacle: **what
could not produce a descent direction was the banded model the escape
factors.**

### 16.2 The reach of the band is NOT the reason

The expected reason was the column integral: an element unknown at cell `j`
changes the composition, hence the opacity, hence the field of every cell
inside `j`, so its true column would carry a dense block below the diagonal
that `|row - col| <= 3 nvar - 1` cannot hold. **Measured, that is not what
happens.** `jacobian_column_reach_beyond_the_band`
(`EXHALE_JAC_COLUMN_REACH=1`) forms a scaled Jacobian column by one full
residual difference per unknown, with no coloring, and reports the fraction of
its 2-norm outside the band. MEASURED on the same state:

| unknown | cell 1 | cell 100 | cell 200 | cell 300 | cell 400 |
|---|---|---|---|---|---|
| mass | 0.0000 | 0.0003 | 0.0014 | 0.0018 | 0.0000 |
| momentum | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| energy | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| element He | 0.0000 | 0.0006 | 0.0009 | 0.0634 | 0.0005 |

The element column is as local as a hydrodynamic one to within a factor of a
few: at worst 6.3 percent of its own norm lies outside the band, on the
column with the smallest norm of the five. A preconditioner correction
carrying the column integral -- a diagonal correction, a rank-one update, or
an identity block for the element rows -- would therefore be a remedy for a
defect that is not there, and none was written.

### 16.3 Two defects that ARE there

**(1) The base element row was structurally empty.** Cell 1 is the Dirichlet
reservoir of the element operator (`solve_mass_fraction` holds
`Xhe(1-Ng:1) = X_base`; the trace solve holds the base mixing ratio at its own
entry value), so `element_transport_residual` fills its rows over `j = 2..N`
and reports cell 1 as zero against the floor. That is the right statement for
a relative MEASURE of the balance and the wrong one for a ROW: the residual
slot of the base element unknown was identically zero for every `Y`, the
unknown was unconstrained, and the banded model was singular by construction.
The equation of that unknown is the boundary condition the operator states,
and it is written as such in `eval_residual`: `F = X(1) - X_reservoir`, with
the reservoir recorded from the state the unknown vector was packed from
(`srow_base_value`), and with the unknown's own scale as its row scale rather
than the transport floor. MEASURED, before and after: structurally empty rows
1 -> 0, and the largest pivot of the banded LU 8.413e+19 -> 1.433e+06 against
a largest band entry of 5.229e+03 either way, i.e. a growth factor of 1.6e+16
-> 274.

**(2) The element row scale collapses, and it is the row scale that broke the
model.** `element_transport_residual` scales a row by the sum of its own terms
with a floor of `1e-20 rho X_base / (R0/v0)`. As a relative measure that is
right: a cell whose terms all vanish reads zero instead of dividing by zero.
As the Newton's ROW SCALE it divides the row's DERIVATIVES as well, and a
floor twenty decades under the row's own rate manufactures entries that are
not sensitivities. MEASURED at the iterate the default route then abandoned:
the helium row scales spanned **1.277e-29 to 2.360** over the column, 29
decades against the 9 decades of a healthy iterate (1.671e-07 to 2.714e+02),
an ordinary derivative `dF/dY = 46.7` became a scaled band entry of
**4.890e+22**, `J^T J` of that band was that one row, and the escape reported
`||grad merit|| = 1.104e+44` and no descent.

This is the defect `docs/p54_base_layer_mass_flux.md` section 10.4 measured for
the HYDRODYNAMIC rows, which is why those rows are not scaled by their own
terms at all: scaling them so left the molecular hot Uranus with no descent
direction after 179 iterations against 9 without. The element rows were
scaled by their own terms and reproduced it. `cell_row_scales` therefore
floors an element row scale at the flow-time rate of the quantity that row
transports -- `rho X` per unit code time for the helium row, and the element's
own column scale for a trace row, which is its mixing ratio up to the smooth
factor `msum/(A_X n_H)`. That is the hydrodynamic convention
(`residual_row_scale` divides the mass row by `rho`, the state itself over one
code time) applied to the transported quantity, and it carries no chosen
constant. **The floor is inactive on a healthy state** -- MEASURED, the row
scales of the first six iterates are unchanged to every digit printed and the
starting merit is 1.443e+02 either way -- and it bites only where the terms
have collapsed below what the state supports.

### 16.4 The certification is not touched

The certification measures the elemental balances by calling
`element_transport_residual` itself, over cells `2..N` and on the operator's
own floor, so neither the base anchor nor the Newton row-scale floor moves a
row measure, a tolerance or a completion flag. What they change is the model
the solve descends on.

## 17. What the solve minimizes, what it is judged by, and why the region gave up (2026-09-08, B5e)

Section 16 left a model an element row can be solved on and a solve that
still did not converge. This section is about the two things that stood
between them: the acceptance and the step-length control.

### 17.1 Three functionals, and the acceptance bounded one of them

READ, from `solve_steady_jfnk`, `resid_relnorm` and
`stationary_rows_of_the_returned_state`:

| what | expression | kind |
|---|---|---|
| the merit the step control descends on | `\|\| F/Drow \|\|_2` over all `nvar*N` entries, `Drow` the Newton's own scales | a 2-NORM |
| the acceptance | `f2_try < (1 - 1e-4 lam) max(f2 of the last five iterates)` | a window on that 2-norm |
| the gate `info = 0` rests on | `max_k max_j \|F_kj\| / residual_row_scale(k,j)`, layer and wind by the larger | a MAXIMUM, on the row's own largest term |
| the certification's rows | the same three, PLUS one entry per registered species row on the operator's own scale | maxima, one tolerance each |

They are not one functional, and the acceptance bounds neither of the last
two. MEASURED on the `mol_diffusion` reload with an element row and the
line-search route: 15 non-monotone acceptances, and **7 accepted steps left
the judged rows above the worst of the previous five iterates**, with `||R||`
running 1.978, 12.53, 1.702, 2.467, 1.907, 5.682 while the merit fell.

**Making it a condition of the STEP was tried and refused by measurement.**
A trial refused when `distance_from_certification` rises above the same
five-iterate window: on the element solve the trust region shrank on 13 such
refusals and stalled at `||R|| = 1.92`, against the 4.345e-09 the same text
reaches without it. The judged measure is a maximum over cells and over
rows, far less smooth than the merit, and a solve legitimately trades one row
against the others on its way down; a step control bounds a smooth merit.

**Where the judgement does belong is the ledger that chooses the state to
hand back**, because that is where the solve makes its claim. It ranked on
`||R||`, which covers the three hydrodynamic rows, while the state handed
back is judged on those AND on every species row the registry carried. It
now ranks on

```
d(state) = max( ||R||/resid_tol , element row/cert_tol_element , carrier row/cert_tol_carrier )
```

the same rows the certification reads, each over the tolerance it is judged
against, so `d < 1` is the certification condition itself. With no species
row `d` is `||R||/resid_tol` and the ranking is the one it always was; the
three-unknown branch keeps the original comparison character for character,
because dividing both sides by `resid_tol` can round two distinct values onto
one.

### 17.2 The eleven exits of the trust region, and which one gave up

`trust_region_step` had seven exits that all end with `model_ok = .false.`
and a shrunken radius, and the summary line called every one of them "model
refused by the ray test". Each exit now carries a code and a sentence, and
under `EXHALE_TR_TRACE=1` a refusal prints the radius in and out, the step
length, the theta cuts, the predicted and actual merit change, the two slopes
the ray test compared with the cancellation floor of their difference, and
the worst scaled row and cell of the refusing trial.

MEASURED, that names the collapse at once. On the element solve the region
fell to 1.109e-16 through **eleven consecutive outer iterations that ended on
"the model promises no reduction"**, with `pred` negative and falling by
exactly a factor 4 at each factor-4 shrink of the radius:

```
delta 2.977e-08  pred -3.610e-04      delta 2.908e-11  pred -3.526e-07
delta 7.269e-12  pred -8.816e-08      delta 1.817e-12  pred -2.204e-08
delta 4.543e-13  pred -5.509e-09      delta 1.136e-13  pred -1.377e-09
```

A predicted reduction exactly linear in the step length with a NEGATIVE slope
is a direction that ascends the model, and shortening a step along it cannot
help: every radius the region tries predicts an increase, so the region is
shrunk again, down to the floor.

### 17.3 The direction that ascends, and the repair

The dogleg's first leg is the Cauchy direction `-g` with
`g = A^T r0` taken from the BANDED `A`, because a transpose is not available
matrix-free. Nothing makes that a descent direction of the model it is used
in: the slope of `1/2 ||r0 - t A g||^2` at `t = 0` is `-(r0 . A g)`, which is
`-||A^T r0||^2` for the exact gradient and has no sign at all for the banded
one. MEASURED at the collapsing iterates, `r0 . A g = -2.088e+13` and
`-2.218e+19`.

So the leg is DROPPED where `r0 . A g <= 0` rather than shortened. The dogleg
then reduces to the Krylov leg cut to the ball, and that leg always predicts
a reduction: GMRES returns `||r0 + A sN|| < ||r0||`, hence
`r0 . A sN < -||A sN||^2/2`, and the model change along `t sN` is positive for
every `t` in `(0,1]`. The radius update itself is unchanged and is the
published one: shrink to `0.25 ||s||` below a ratio of 1/4, double above 3/4
on the boundary, accept above eta (Nocedal and Wright, *Numerical
Optimization*, 2nd ed., Algorithm 4.1, with the classical Powell dogleg).

`pred` is also no longer computed as `1/2(||r0||^2 - ||r0 + As||^2)`. That
differences two sums of the same size, so on a short step its value is the
rounding of `||r0||^2` and its SIGN -- the thing that decides whether the
region shrinks again -- is noise. The algebraically identical
`-(r0 . As) - 1/2 ||As||^2` carries no cancellation.

The Cauchy step LENGTH is left as it was, `||g||^2/||A g||^2`. Replacing it
by the line minimizer along the direction actually used,
`(r0 . A g)/||A g||^2`, is the more defensible expression and it was written
and measured: it moves the Cauchy point at every iterate that uses the leg,
and the element solve then took a path that was still crawling at
`||R|| = 1.8e-04` after 81 outer iterations. Nothing measured says the
original length is wrong, so only the ascent case is changed.

### 17.4 The trust region is the step control of the coupled route

MEASURED on the `mol_diffusion` reload (one element row, `nvar = 4`), from
the same 300-step marching snapshot, one thread, everything else equal:

| step control | how it ended | `\|\|R\|\|` | elemental He/H row of the returned state |
|---|---|---|---|
| the pseudo-transient line search | `no descent direction exists for the banded model`, iteration 18 | 1.702 | above tolerance |
| the scaled trust region | `\|\|R\|\| < Resid tol` at iteration 66, 54 accepted steps | 4.345e-09 | 4.633e-11 against `cert_tol_element` = 1e-08 |

So the trust region is now the step control wherever the system carries a
species row, and `EXHALE_TRUST_REGION=0` restores the line search there for
measurement. **No case of the regression matrix carries a species row**, so
this default moves no golden; what it changes is the route a coupled solve
takes, and it changes it from one that does not converge to one that does.

The three-unknown route never reaches one statement of the trust region.
MEASURED, `wasp_full_newton` reloaded from its own converged state is
byte-identical before and after, `info = 0`, CERTIFIED, at the same 6 outer
iterations and `log10 Mdot = 13.30`.

### 17.5 The carrier arm is not rescued by this, and the ledger says why

MEASURED on `mol_carrier` reloaded with `Coupled carrier solve: True`,
`Solver: Newton 100.0`, `Resid tol: 1.0e-8`, one thread: the line search
aborts at iteration 16 on `no descent direction exists` at `||R|| = 1.855`
and the trust region STAGNATES at iteration 48 at 1.873, the H2 carrier row
of the returned state 2.081e-01 and 1.617e-01 respectively, both seven
decades above `cert_tol_carrier`.

**The obstruction there is not the step length**, and the exit ledger is what
says so: 18 of the trust region's 49 iterations ended on `no Krylov direction
could be sampled`, and 500 trial states were refused for the element budget
(115 on the line-search route). A carrier solve that has no Newton model to
step on is a different problem from an element solve that had one and could
not descend it.

## 18. A residual evaluation is a function of its arguments (2026-09-08, B5g)

Section 17 left one measurement it could not explain: `EXHALE_SPECIES_JAC_TEST=1`
moved the solve it was measuring. This section is about that, and about what
it uncovered.

### 18.1 The by-products of an evaluation, and who reads them back

READ, from `eval_residual` and the modules it calls. `F(Y)` is written into
the caller's argument, but the evaluation also leaves this behind:

| what | where it lives | who reads it back |
|---|---|---|
| `erow_he`, `escale_he`, `erow_tr`, `escale_tr`, `erow_tr_carried` | `steady_newton` | `cell_row_scales` (the Newton ROW SCALE) and `certified_row_measures` |
| `carrier_relnorm_last`, `carrier_cellmax_last`, `carrier_cell_worst` | `steady_newton` | the acceptance gate, the best-iterate ledger, `certified_row_measures` |
| `col_scale_car`, `row_scale_car`, `row_terms`, `headroom_car`, the frozen backgrounds | `diffusive_photochemistry` | `cell_state_scales`, `cell_row_scales`, the element headroom |
| `bg_cell` | `ionization_equilibrium` | the carrier and element operators of the NEXT evaluation |
| `n_headroom_iterate` | `steady_newton` | the element-budget comparison of every trial |

The solve evaluates the residual at three kinds of point: the state it holds,
a trial it may adopt, and a PROBE -- a Jacobian column, a Krylov product, a
point on a trust-region ray, a self-test. A probe is not a state anyone
keeps, and until this item it left the table above holding its own numbers.

**MEASURED, the size of that.** `mol_diffusion` reloaded with
`Coupled carrier solve: True` and one element row, with and without
`EXHALE_SPECIES_JAC_TEST=1`: iterations 1 to 5 agree to every printed digit,
iteration 6 reads `||R|| = 1.958` against 1.956 and merit 4.34e+01 against
4.30e+01, and iteration 7 takes 32 Krylov iterations against 3. The
self-test's probes at `Y + eps r` had left that state's element row scales
for the first outer iteration to scale itself by.

`eval_residual` now takes `state_is_discarded`, and an evaluation so marked
holds the table above aside on entry and puts it back on exit
(`residual_evaluation_products`). It is the statement
`save_carrier_module_state` already makes for the certification, raised to
the residual: measuring a balance is not an operation on the state. The two
self-tests additionally tag their sweeps as candidate sweeps, so a
diagnostic no longer writes the run's own acceptance ledger or its non-root
streak.

### 18.2 The row scale is taken after the iterate is evaluated

`cell_row_scales` claims to be the row scaling of `Y`. It is not a function
of `Y`: it reads `escale_he` and `carrier_row_term_scale`, which the LAST
residual evaluation formed. At the top of an outer iteration it used to be
called BEFORE the iterate's own evaluation, so the scaling of iteration
`k + 1` was that of whatever state iteration `k` evaluated last -- the
accepted trial where a step was taken, and a point on a trust-region ray
where it was not. It is now called after that evaluation, so the scaling and
the residual it divides come from one evaluation of one state.

### 18.3 The acceptance gate reads the state the solve holds

`carrier_relnorm_last` is the carrier measure of the last EVALUATION.
`steady_gates_met`, the best-iterate ledger and the choice of the state to
hand back were all given it, at points where the last evaluation was a
refused trial or a Krylov probe. They are now given
`carrier_relnorm_state`, refreshed wherever `n_no_chem_root_state` is: at
the iterate, at an accepted trial, and at the best-iterate restore.

### 18.4 What this moved, and the defect underneath it

MEASURED, `mol_diffusion` with the element row, before and after: the solve
that reached `||R|| = 4.345e-09` in 66 outer iterations now STAGNATES at
`||R|| = 1.957` after 17. The two runs are identical for five iterations and
part at iteration 6 in the trust region's `actual` reduction -- 4.389e+03
against 4.334e+03 for the same `pred` and the same step -- which is a
percent-level difference in `Drow`, not a round-off one.

**Two evaluations of the SAME state give row scales that differ by percent,
because the eliminated residual of that state is not single-valued.**
MEASURED with `EXHALE_RESID_SC_PROBE=1` on the hand-off state of this case:
the residual evaluated from the run's own composition and from a displaced
seed differ by **9.577e-01 in the units of `Resid tol`**, against
`||R|| = 1.986` -- the seed dependence is half the residual itself. The
elimination stopped after 2 passes from one seed and 19 from the other, both
satisfying the same `1e-8` test on the composition increment, at two
different compositions.

**Tightening the elimination does not remove it**, which is what says the
stopping test is not the cause. MEASURED, the same probe at three settings
of (`EXHALE_RESID_EQ_TOL`, `EXHALE_RESID_SC_MAX`):

| eq tol | pass cap | passes taken (own seed / displaced) | seed dependence |
|---|---|---|---|
| 1e-8 | 25 | 2 / 19 | 9.577e-01 |
| 1e-12 | 200 | 200 / 200 | 4.186e+00 |
| 1e-14 | 500 | 500 / 500 | 9.693e-01 |

At the two tighter settings the increment test is never satisfied from
either seed: the Picard iteration runs the whole cap and the composition is
still moving. Five hundred passes from each of two seeds leave the residual
of one state differing by about one unit of `Resid tol`. The elimination
therefore does not have a fixed point it reaches from either seed on this
state, and `F(Y)` is not a function of `Y` at the level of the residual
itself.

So the 66-iteration convergence of section 17.1 is a property of one
particular sequence of evaluations of this case and not of the state: two
independent perturbations of order 1e-8 in the row scale -- the self-test's,
and the change of which evaluation the scale is taken from -- both send it to
the same stagnation. **What is recorded here is that the element arm of the
coupled solve rests on a seed-dependent residual**, and that no step control
can be judged on this case until the elimination has one fixed point. The
element row measure of the stagnated state is 9.992e-01 against
`cert_tol_element` = 1e-08, so nothing about that tolerance follows from it
either.

The three-unknown route is untouched by all of it: with no species row
`erow_*` and the carrier arrays are never written, `carrier_relnorm_state`
is zero and the gate ignores it, and the scale reordering has no reader
between the two positions.

### 18.5 The carrier arm: what "no Krylov direction could be sampled" was

Section 17.5 left the carrier arm stagnating with 18 of its 49 iterations
ending on `no Krylov direction could be sampled` and a count of "500 trials
refused for the element budget". Both statements are now measured, and the
second was mislabeled.

`eval_residual` records WHICH screen refused a sample, at which cell and
which unknown (`eval_refusal`), and the solve prints the census of the
screens and, at every exit that ends on "no sample", the last refusal.
MEASURED on `mol_carrier` reloaded with `Coupled carrier solve: True`:

```
(JFNK) [TR]      no step: no Krylov direction could be sampled
(JFNK) [TR]      refused: carrier H2 of cell 500 is negative,-7.444E-19;
                 1 species unknown(s) of this state sit on their lower bound
(JFNK) residual samples: admitted 1411, refused for a negative species unknown 514,
       an element fraction above one 0, cells outside the element budget 0,
       a non-finite sweep 0, a non-finite row 0, no chemical root 0
```

All 18 exits carry that same line, and every one of the 514 refusals of the
solve is the positivity screen; not one is the element budget. The counter
the older line read, `n_trial_headroom`, is the SUM over the screens and its
label named one of them, so it is corrected and the census above is what a
reader should use.

The mechanism: `n(H2)` of cell 500, the outermost physical cell of a wind
whose H2 is long gone, is zero. `jv_product` samples `F(Y + eps v)`;
wherever the Krylov direction's component at that unknown is negative the
sample puts that density at about `-7e-19` of the code density unit and the
screen refuses it, at every step length, because the unknown is AT the
bound. GMRES then returns `gm_iters = 0`, the dogleg has no Newton leg, and
the trust region takes no step. The banded model is NOT singular there (the
empty-row statement never fires on this registry) and the Jacobian COLUMNS
are unaffected: the solve reports `0 Jacobian color(s) zeroed, 18 GMRES
cycle(s) truncated`, because a column probe steps one unknown in the
positive direction and never leaves the lower bound.

`jv_product` therefore takes the other side of the same derivative,
`[F(Y) - F(Y - eps v)]/eps`, when no forward step lands on a describable
state. It is first-order accurate in `eps` exactly as the forward quotient
is, it steps away from the bound the forward step crossed, and it is tried
only after every forward halving has been refused, so no case in which the
forward sample exists can move. MEASURED, before against after on the same
state:

| | without | with |
|---|---|---|
| steps accepted / rejected / model refused, of 49 | 31 / 18 / 18 | **49 / 0 / 0** |
| iterations ending "no Krylov direction" | 18 | **0** |
| final trust radius | 8.006e-03 | 2.050e+00 |
| GMRES cycles truncated / colors zeroed | 18 / 0 | 0 / 0 |
| `\|\|R\|\|` returned, H2 row of the returned state | 1.873, 1.62e-01 | 1.873, 1.62e-01 |

**The arm is not converged by this and the repair does not claim to be.**
What it removes is an iteration with no Newton model; what remains is the
same bound seen from the trial side, the dogleg cutting its step back by
`2^-25` to `2^-49` to keep it describable and the merit not moving. A trust
region cannot express that: a radius bounds a norm and the obstruction is a
face of the feasible set. Handling it needs a decision about the unknown
space -- a probe context that admits a sample outside the bound, a positive
parametrization of a carrier density, or a projection onto the free
variables -- and none of the three is taken here.

## 19. A coupled steady solve may not eliminate H2 (2026-09-08, B5h, B5i)

The coupled route of section 14 registers a row and an unknown for every
transported balance the configuration activates. Section 18 asked whether one
residual evaluation is a function of its arguments; this section records the
one configuration in which it is not, and the refusal that now keeps it out.

**The configuration.** `Coupled carrier solve: True` with
`Molecular chemistry: True` and `Molecular carrier transport: False`. The
carriers are then not registered, so `n(H2)` is not an unknown of the
stationary system: it is whatever the local-equilibrium sweep returns for the
composition it was seeded with, and `R(U) = L(U) + S(U, c^*(U))` is a function
of `U` only if that sweep reaches one composition `c^*(U)`.

**In the shielded layer it does not.** The fast chemistry there cycles
`H2 -> H2+ -> H3+ -> H2` without changing the number of H2 nuclei, so the local
balance rows fix only the partition among the molecular species and leave their
sum free; the content is set by the slow formation and dissociation and by
transport, neither of which a local equilibrium sees. MEASURED on the hot
Uranus element state (`mol_diffusion` reloaded, `Resid tol` 1e-8, one element
row, one thread), by seeding the sweep with the state's own composition and the
H2 partition scaled by `(1 + delta)` and running the elimination to its fixed
point:

| delta | H2 eliminated | H2 a Newton unknown |
|---|---|---|
| 1e-08 | 3.165e-05 | 5.685e-10 |
| 1e-06 | 3.155e-03 | 9.736e-10 |
| 1e-04 | 2.405e-01 | 1.485e-06 |
| 1e-02 | 9.714e-01 | 1.450e-02 |

The column is `||R_d - R_0||`, the difference of the residual of the SAME state
from two seeds. With H2 eliminated it is linear in `delta` over four decades
with a gain of 3.2e3, so it is a derivative and not a threshold: about 0.77 of
any seed perturbation of the layer's H2 content survives every pass of the
sweep, for ever (the same 0.77 at 1e-6 and at 1e-2, so the map's eigenvalue in
that direction is 1), and the base cell's energy row, being a
near-cancellation of the fluxes its continuous-temperature ghost produces,
amplifies the surviving part by ~4e3. One `Resid tol` is reached at a seed
change of 3e-12 relative, and two evaluations inside a Newton differ in their
seed by vastly more than that: the iterate's own composition against the
composition the previous trial's sweep left. Every convergence claim on that
element solve rests on one sequence of evaluations rather than on a root.

**With the carriers transported the residual is a function of its unknowns
again.** `Molecular carrier transport: True` sets `ieq_cell%x_h2_fixed` at
every cell the transport operator owns, so the H2 row of the network becomes
`x - x_fix` with `x_fix` the Newton unknown and the content is solved from its
transport balance. The same measurement then reads 9.7e-10 at the seed
perturbations a Newton generates, below `Resid tol`, and `x(H2)` of the front
cell is bit-identical from every seed of the ladder.

**The refusal (item B5i, decision (a) of
`docs/To_be_determined_by_user_20260906.md` section 10).** `input_read.f90`
refuses `carrier_in_newton .and. thereis_mol .and. .not. carrier_transport`
after every key is resolved, so that the default `carrier_transport` takes from
the oxygen chemistry is final, and the message names
`Molecular carrier transport: True` as the remedy. The rule this follows is the
one the SED-coverage and retired-key refusals follow: a configuration the
method cannot solve is refused with the remedy named, rather than repaired
behind the user's back. Turning the transport on automatically was the
alternative and was declined, because it would change the unknown set of the
solve from a key the user did not set. The scope is narrow: a run whose only
transported balance is an element (`He_diffusion` in an atomic gas) still uses
the coupled route, and marching is untouched, never eliminating a quantity it
also has to determine. `src/tests/grid_and_gates/coupled_carrier_h2_row.sh`
states the three outcomes.

**What this does not settle.** With H2 carried, a seed displaced by five cells
still reaches a SECOND fixed point of the elimination, 19 percent apart in
`x(H I)` at r = 1.086 and 12 percent in `x(He I)` at r = 1.100, both exact
roots of the network at the same incident field. Which root a stationary
residual may land on, or a statement that the ionization front is unresolved on
this grid, is a physics rule the code does not have.

## 20. The carrier unknown of the coupled solve is `ln n` (2026-09-08, B5j, B5k)

Section 17 gave the coupled route a trust region and section 19 fixed which
configurations may use it. This section records what the carrier unknown IS,
which is a different question from how a step in it is controlled, and the
decision that settled it (`docs/To_be_determined_by_user_20260906.md`
section 11, option (a), user, 2026-09-08).

**The problem.** A carrier density is non-negative and the stationary system
has no way of knowing it. Carried as `n` itself, the unknown has a face of a
box under it: the Newton proposes negative densities, and an unknown that has
REACHED zero cannot be sampled at all, because a two-sided finite difference
along a direction with an outward component leaves the admissible set on both
sides. MEASURED on the coupled `mol_carrier` reload (item B5j, section 1.2 of
its report): at outer iterate 10 the carrier row of cell 215 has its root at
`n = 2.75e-06` of the code density unit, `x(H2) = 8.1e-05` of that cell's own
mass, while the Newton holds the unknown at exactly zero and the row reads
`-1.00000` of its own scale. Twenty-six to thirty-two species unknowns of one
state sit on their bounds at once.

**The two arms, MEASURED on `mol_carrier`** with `Coupled carrier solve: True`,
`Molecular carrier transport: True`, `Solver: Newton 100.0`, `Resid tol` 1e-8,
one thread, ONE binary, the two runs an environment variable apart, each run
to its own exit (item B5k):

| | the density unknown | `ln n` |
|---|---|---|
| outer iterations, how it ended | 63, stall detector | 66, stall detector |
| `\|\|R\|\|` of the state handed back | 1.415 | **0.9957** |
| best merit `\|\|F/Drow\|\|_2` | 2.52 | **0.960** |
| best judged distance | 1.415e+08 | **1.004e+08** |
| carrier row of that state, worst cell | 1.000 at j = 215 | **0.4223 at j = 205** |
| carrier row, volume-weighted | **0.0151** | 0.0215 |
| steps accepted / rejected | 64 / 0 | 63 / 4 |
| steps that were the Cauchy point alone | 28 | **0** |
| most unknowns held on an active bound | 10 | **0** |
| products taken on the backward side | 539 | **0** |
| samples refused for a negative carrier | 0 | 0 |
| samples refused outside the element budget | 52 | 157 |
| smallest adopted carrier density | **0** | 2.209e-11 |
| that density over the floor of its cell | **0** | 2.53e+15 |
| adopted carriers at or below the floor | **716** | **0** |
| `log10 Mdot` | 10.62 | 10.63 |

The one entry in which the density arm reads better is the volume-weighted
carrier row, and it is the same fact as the entry above it: pinning the
carrier of one cell at zero makes THAT cell's row read 1 and leaves the rest
of the column alone, while the logarithm spreads a smaller error over more
cells. The judged distance, which is what the best-iterate ledger ranks on,
prefers the logarithm by 29 percent.

`ln n` is better on every measure, and it is better for a reason rather than by
luck: the parametrization removes the bound instead of negotiating with it. It
takes no projection, holds no active bound, needs no Cauchy-only step and no
backward sample. The a-priori worry that `dR/d ln n = n dR/dn` would empty the
banded model's columns where the carrier vanishes is refuted by measurement:
structurally empty rows 0, empty columns 0, smallest `|U(j,j)|` of the banded
LU 1.83.

**It is the default.** `EXHALE_CARRIER_LOG_UNKNOWN=0` restores the density
unknown, which exists to be measured against. The element fractions are
unchanged: they are bounded on both sides, they are linear in the composition,
and the bound-aware region of section 17 keeps them.

### 20.1 What the logarithm costs, and where each cost is stated

Three things, each of them written at the code site that carries it
(`steady_newton.f90`).

* **The Jacobian column carries the factor `n`.** The column scale of the
  unknown is 1 -- a logarithm is its own scale -- so the banded model, which
  differentiates the residual column by column, and the matrix-free action,
  which differentiates it along one direction, both step `ln n` by a multiple
  of `sqrt(eps_mach)` and therefore both probe a RELATIVE change of the
  density. Stated at `cell_state_scales` and at the declaration of
  `carrier_unknown_is_logarithmic`.
* **The trust region is asymmetric in the density.** Its radius bounds a norm
  of the scaled step, so a radius `Delta` permits `n -> n exp(+/-Delta)`: the
  same radius multiplies and divides the carrier by the same factor and can
  never carry it to zero. The region is symmetric in the ratio and not in the
  density, which is the geometry a quantity spanning decades asks for.
* **`ln 0` is not a number**, and zero is what the outermost cell of an ionized
  wind holds, so the space needs a smallest representable density.

### 20.2 The floor, and why it is the element budget

The floor under `ln n` is numerical, because a carrier density has no physical
lower bound other than zero. What it may not be is a constant chosen in the
solver: the quantity it bounds spans the column, and a fixed number would be a
hard bound in the base and unreachable noise in the wind.

It is **1e-20 of the element budget of that carrier in that cell**, cell by
cell and carrier by carrier -- the largest density the carrier can reach there,
the free density of its element divided by the nuclei one carrier holds
(`carrier_headroom`, `diffusive_photochemistry.f90`), times the fraction the
operator that owns the carrier row already uses for the same statement.
That fraction is not new here: `carrier_residual` adds `1e-20` of the same
element's density to its row scale as an absolute floor, with the reason
written there, "a density that small cannot change any observable, so a row
below it IS converged". So the unknown space and the row measure agree about
which densities the discretization resolves: the space cannot name a density
the row would call converged, and it can name every density the row can
resolve. The two references differ only in the direction that is safe: the
headroom is the element's free density less the nuclei sitting in the stages
the carrier step holds frozen, divided by the nuclei one carrier holds, so
this floor is at or below the row's own and the space never forbids a density
the row measure can still resolve.

MEASURED on `mol_carrier`, N = 500: the floor runs from 5.599e-28 to 4.608e-21
in code density units, seven decades across the column, which is the spread a
constant would have had to stand in for.

**When the space is armed.** The budget is the one the carrier operator freezes
at the iterate, and it does not exist until that operator has assembled once,
which first happens inside the solve's own opening residual evaluation. So the
opening state is packed and evaluated as a density and the carrier slots are
rewritten as logarithms immediately after, with the floor that evaluation
produced (`arm_carrier_log_unknown`). The residual is a function of the state
and not of the coordinates the state is named in, so that opening residual
stands for the rewritten vector as it stood for the packed one. A cell whose
budget is not a positive finite density leaves the solve in the density
unknown and says so in one line, rather than substituting a number of its own.

**Whether the floor is ever reached is measured, not asserted.** The solve
reports the smallest ratio of an adopted carrier density to the floor of its
own cell and the count of adopted carriers at or below it. MEASURED on
`mol_carrier`, one binary and one environment variable apart: with `ln n` the
ratio is 2.53e+15 and the count 0; with the density unknown the ratio is 0 and
the count 716, which is the same fact as the exact zero of the table above. On
`mol_diffusion` with H2 carried the log arm reads 4.02e+08 and 0.
`src/tests/steady_species_rows/run.sh` carries both as rows
(`carrier_is_above_its_own_floor_*`, `carrier_never_reached_its_floor_*`),
GREEN on the log runs and RED on the density one.

### 20.2a What limits the log arm now, MEASURED

Neither arm converges. What stops the log arm is no longer a lower bound: on
`mol_carrier` it holds no unknown on a bound, takes no Cauchy-only step, needs
no backward sample and has no sample refused for a negative carrier. Its steps
are model-accurate -- the reduction ratio of the last ten accepted steps runs
0.39 to 1.09, eight of the ten within 12 percent of one -- and SHORT: `||s||`
is 0.089 to 0.17 against a region radius of 2.535, because the dogleg is halved
three to five times before the trial is describable (`cuts` at least 3 in 38 of
68 iterations), and the screen doing the
refusing is the element budget: **157 residual samples refused for cells
outside it, against 0 for a negative carrier**. The region itself is not the
limit; its radius spent 30 of 68 iterations at its own ceiling.

So the completion the density arm was missing -- Bertsekas' two-metric
projection, Newton on the free variables and a scaled gradient on the bound
ones -- is **not** what the log arm needs: it has no bound variables. The
constraint now shortening its steps is the carrier CEILING, the element budget,
which `species_unknowns_outside_their_bounds` deliberately does not treat as a
face of the box ("the element budget above it is a comparison and not a bound")
and which is enforced only by refusing an evaluated state. Giving the model
that constraint is a different repair from the one the density arm asked for.

### 20.3 What the certification is not asked to do

`docs/To_be_determined_by_user_20260906.md` section 12, option (c): no
certification change. A carrier at an active bound is an error to the
certification and a constrained optimum to the solver, and under `ln n` it
cannot occur, so a Karush-Kuhn-Tucker branch for it would be a rule nobody
exercises. The statement the certification keeps making is the physical one: a
wind carries a positive density of every carrier.

## 21. The element budget is a face of the unknown box (2026-09-08, B5l)

The step of the coupled carrier solve was short because it kept walking into a
constraint the model could not see. B5k measured it: with the carrier in
`ln n` the steps are model-accurate (reduction ratio 0.39 to 1.09) and
`||s||` runs 0.089 to 0.17 against a trust radius of 2.535 that sits at its own
ceiling for 30 of 68 iterations, each trial cut by three to five halvings
before it is describable, and the screen doing the refusing is the ELEMENT
BUDGET: 157 residual samples refused for cells outside it against 0 for a
negative carrier.

The budget is the statement that a carrier cannot hold more nuclei of its
element than the cell has free (`carrier_element_headroom`, and
`hydrogen_available_to_carriers` for what the frozen stages must be left).
That is a bound on the unknown, and it was enforced only by refusing an
already evaluated state, so the region had no way to size a step against it:
a radius bounds a norm and this obstruction is a face.

### 21.1 One box, three readers

`freeze_species_unknown_box` is now the only place the admissible set is
written down, as one lower and one upper bound per unknown of the flat vector,
in the coordinates the unknown is carried in. Three routines read it and none
of them forms a bound of its own any more:

| reader | what it does with the box |
|---|---|
| `species_unknowns_outside_their_bounds` | writes a trial onto the nearest face |
| `fix_active_species_bounds` | the epsilon-active set: which unknowns the step may not move |
| `largest_step_inside_the_species_box` | the fraction-to-the-boundary rule for a finite-difference probe |

The three had the bounds written out three times, and they had **drifted**:
the projection read the ceiling from the TRIAL's own density slot while the
other two read it from the ITERATE. So the set a trial was written onto and
the set the Krylov cycle was sampled under were not the same set. They are one
set now, and it is frozen for the step, which is what the projection needs to
be a projection and what keeps the sampled operator linear in its argument.

### 21.2 The faces of a carrier unknown

* Below, the smallest density the space can name: zero in the density
  unknown, `carrier_log_floor` in `ln n` (section 20.2).
* Above, the cell's own density: a carrier cannot be more of the cell than
  the cell is.
* Above, the element budget, which is the tighter of the two wherever it
  exists. For H2 it is half the hydrogen the carriers may take.

**The budget is FROZEN at the outer iterate, not evaluated at the trial**, and
that is the carrier operator's own choice, not a convenience of the step
control: `headroom_car` is assigned in `carrier_steady_residual` under
`weno_mode == 1`, the one evaluation per outer iteration that IS the iterate,
because the budget depends on the very unknown it constrains (`n_H_free`
counts H2's own two nuclei) and evaluating it at the trial would make it move
with the step. Refreshing it at every non-probe evaluation was tried and is
worse: the last line-search trial then leaves its own budget behind and the
escape that follows finds the ITERATE outside it. The box therefore reads one
fixed number per cell and carrier for the whole of a step, and the box and the
screen of `eval_residual` describe one state.

**The budget face is raised to the iterate where the iterate stands outside
it.** Because the budget is frozen at a state that depends on it, the iterate
itself can sit a cell or two outside -- which is exactly why the screen in
`eval_residual` is a comparison with the iterate's own count and not a
threshold ("a rule that rejects the neighborhood of the point it stands on
cannot be used to leave that point"). A face cutting through the point the
step starts from would leave that unknown no room on either side and the
fraction-to-the-boundary rule would return zero for the whole direction.
Raising the face keeps the iterate feasible and still forbids every cell
inside the budget from leaving it, which is the screen's statement made cell
by cell instead of as a count. How often it is raised is counted and reported
(`n_budget_face_at_the_iterate`). The cell-density face is not raised: it is a
hard statement about the state and a carrier above it is to be pulled back in.

The face is stored one floating-point step inside the budget
(`budget_face_inside = 1 - 8 eps`). That is a representability guard and not a
margin: in `ln n` the face is stored as `log(budget)`, and `exp` of it can
land one unit in the last place above the budget the screen compares against,
which would count a cell written exactly onto the face as outside it.

**The predicted reduction is recomputed for the projected step**, as B5j's
projection already did: the projected step leaves the plane the dogleg was
built in, so `A` applied to it is not a combination of `A s_U` and `A s_N`,
and one matrix-free product buys the model of the step actually taken.

### 21.3 A box side here, and a hyperplane the code does not state

The face made here is a side of a box, and that is exactly what the
constraint the code states is. `carrier_element_headroom` bounds each carrier
by the FULL free density of its element, independently of what the other
carriers of that element hold: `0.5 nH_avail` for H2, `nH_avail` for H+,
`nO_free` for OH and for H2O, `min(nO_free, nC_free)` for CO. One bound per
carrier, so the projection treats it exactly, and the shipped set that the two
measured cases carry is `{H2}` alone -- a single carrier, a single element, a
single side.

**Carriers do share elements, and the joint constraint is a hyperplane that
`carrier_element_headroom` does not state.** With `ionization_transport` the
set is `{H2, H+}` (or `{H2, OH, H2O, CO, H+}`) and both draw hydrogen; with
`thereis_oxychem` the set is `{H2, OH, H2O, CO}`, of which OH, H2O and CO all
draw oxygen. The write-back names the physical statement: its hydrogen
remainder is

    rest = nH_free - 2 n(H2) - n(OH) - 2 n(H2O) [- n(H+)]

so the hydrogen-bearing carriers are jointly bounded by

    2 n(H2) + n(OH) + 2 n(H2O) + n(H+) <= nH_avail,

and the oxygen-bearing ones by `n(OH) + n(H2O) + n(CO) <= nO_free`. Two
consequences, both READ from the source and neither of them repaired here
(`carrier_element_headroom` is the carrier operator's, not the solver's):

* the joint bound is a hyperplane, not a box side, so the exact projection
  written here does not cover it and a different projection would be needed;
* OH and H2O hold hydrogen nuclei (one and two) that their own headroom does
  not mention at all, so a state can satisfy every headroom and still leave
  the write-back's hydrogen remainder negative, which is the branch that sets
  H I, H2+ and H3+ to zero.

Neither reaches the two cases measured here, whose carrier set is `{H2}`.

### 21.4 An unknown whose box is a point

Where a carrier's element has no free density in the cell, the budget is zero
and the only admissible carrier density is zero. Such an unknown cannot move
whichever way the model wants to move it, so `fix_active_species_bounds` holds
it without asking the gradient. Left free, its component of any direction
would give the fraction-to-the-boundary rule no room and cost the whole
Krylov cycle rather than only that component.

### 21.5 What the face changed, MEASURED

`mol_carrier` reloaded with `Coupled carrier solve: True`,
`Molecular carrier transport: True`, `Solver: Newton`, `Resid tol: 1.0e-8`,
one thread, the carrier in `ln n`; one binary, `EXHALE_SPECIES_BUDGET_FACE=0` (retired by N4b, 2026-09-09; the measurement arm is now `EXHALE_ELEMENT_CONSTRAINT_ROWS=0`)
apart. The switched-off arm reproduces the entry text iteration for iteration
and byte for byte in both output files, so the face is the only behaviour this
item changed on this route.

| `mol_carrier` | budget not a face | **budget a face** |
|---|---|---|
| outer iterations, how it ended | 66, stall detector | 52, stall detector |
| `\|\|R\|\|` handed back | **0.9957** | 1.174 |
| best merit `\|\|F/Drow\|\|_2` | **0.960** | 1.10 |
| steps accepted / rejected | 63 / 4 | 52 / 1 |
| trials halved 3 or more times | **38 of 68** | **0 of 53** |
| residual samples refused, element budget | **157** | **0** |
| trials written onto a face | 0 | 29 |
| largest adopted carrier over its own budget | 0.9886 | 0.9878 |
| adopted carriers above their own budget | 0 | 0 |
| iterations ending "the model promises no reduction" | 0 | **19** |
| worst Krylov relative residual reached (asked 0.1) | 0.938 | 0.569 |
| `log10 Mdot` | 10.63 | 10.63 |

So the face does what it was built to do: the halvings against the budget and
the refusals for leaving it are gone, exactly. It does NOT converge the arm,
and it does not improve `||R||` on this case; it moves the limit to a
different, named place. **The Krylov leg is now what refuses the step**: with
the trial on a face, 19 of 53 iterations end because the model of the projected
step promises no reduction, and the cycle reaches 0.53 to 0.57 of its
right-hand side where it was asked for 0.1.

**The iterate never stood outside its own budget on either coupled case**
(`n_budget_face_at_the_iterate` = 0 in both), so the raise of section 21.2 is
a rule that held nothing here. It is not idle: with
`EXHALE_SPECIES_BOUND_HOLD=all` on `mol_diffusion` it fired 36 times.

On `mol_diffusion` with H2 carried the face does the same to the step control
and the solve stalls at the same place it always has (the H2 front at cell
212, B5h's obstruction; `||R||` 1.982 against 1.943, `log10 Mdot` 7.97 either
way): trials halved 6 to 11 times in 50 of 50 iterations against 0 of 41,
budget refusals 388 against 0. What appears there instead is the probe: 419
Jacobian-vector products taken on the backward side of their direction, 26
unknowns held on an active bound, 12 steps the Cauchy point alone and 18
truncated cycles, against none of any of those without the face. The element
row is 1.000 of its own scale at cell 212 in both, which is what that case
stalls on.

### 21.6 The Krylov budget is now measurable

`gm_m = 40` is a literal at the `solve_steady_jfnk` call site in
`src/EXHALE_main.f90`, the only number in the coupled route's step control a
user cannot reach. `EXHALE_GM_M` overrides it for measurement (the caller's
value is the default), and `pgmres` now reports the relative residual the
cycle REACHED against the tolerance it was asked for, which `gm_iters` alone
cannot say: a cycle that ran `m` products may have converged at exactly `m` or
may have been cut off there. The dogleg falls monotonically along its path
only when its Newton leg is the model's own minimizer (section 159.3 of the
update log), so the subspace size is part of the step control and not only of
its cost.

**MEASURED, and the answer is that the subspace is not the limit.** Same case
and protocol, the face on, `EXHALE_GM_M` = 40, 80, 160:

| `gm_m` | 40 | 80 | 160 |
|---|---|---|---|
| outer iterations | 52 | 75 | 52 |
| `\|\|R\|\|` handed back | 1.174 | 1.103 | 1.271 |
| best merit | 1.10 | 0.943 | 0.916 |
| worst relative residual reached (asked 0.1) | 0.569 | 0.630 | 0.614 |
| `\|\|s\|\|` of the last ten accepted steps | 0.076 to 0.21 | 0.039 to 0.24 | 0.051 to 0.20 |
| `log10 Mdot` | 10.63 | 10.64 | 10.64 |

Four times the subspace buys nothing: the cycle spends every product it is
given and stops at the same 0.6 of its right-hand side, the step keeps the
same length, and `||R||` moves inside the spread of a stall. **The default
stays 40.** The cycle is not short of subspace, it is stagnating on this
operator, so what to change is the preconditioner or the operator and not the
budget.

> **Superseded 2026-09-10 (decision 22 a, item N30):** every `cert_tol_element`
> / `cert_tol_carrier` = 1e-8 quoted above is historical; the values are
> 1e-5 gating at r >= 1.20 with the rows below reported, through
> `cert_tol_element_at(r)` / `cert_tol_carrier_at(r)`.

## 22. Stage 2 outcome: what the design became (2026-09-10)

Sections 1 to 21 are the record of how the stationary solver was built and
what each measurement refuted along the way; they are left as they were
written. This section states what the design IS at the end of stage 2, after
items N1 to N30, and what it does not yet do. Numbers are LOGGED from
`docs/Update_EXHALE.md` section 7 and `docs/ISSUES_20260909.md` unless marked
otherwise.

### 22.1 The system

The unknown space is the transported set (section 13): the three
hydrodynamic rows of every cell, plus one row and one unknown for every
transported balance the configuration activates. An element `metals.inp`
does not carry is no longer registered at all: three absent elements
(Si, K, S) had been registered as mass-fraction unknowns with an identity
row sitting exactly on their lower bound, and the fraction-to-the-boundary
rule then returned zero for any direction with a component there, so the
first Krylov cycle took 0 products of 40 and the trust radius was
initialized to zero, an absorbing state (N1). The carrier unknown is
`ln n` (decision 11 a, section 20), which removes the bound rather than
projecting against it.

A coupled steady solve may not eliminate H2 (section 19, decision 10 a): a
molecular configuration without `Molecular carrier transport` is refused at
startup, because the local-equilibrium closure does not determine the
shielded layer's H2 content and the residual would not be a function of its
unknowns. `Ionization transport: True` with `Solver: Newton` is supported
exactly when `Coupled carrier solve: True` is set, and refused without it.

### 22.2 What each iterate must return

Two constraints were added to the write-back after section 21.

- **The element totals are returned as handed** (N4b, decision 14 route i).
  Every coupled-carrier state of B5j and B5k had departed from the input's
  He/H by 4 to 8 percent (density unknown) or 41 to 43 percent (`ln n`
  unknown), while 12000 marching steps with the carrier transported keep it
  to 1e-11: the Newton moved n(H2) with no element row holding the hydrogen
  total, and the sweep then read the nuclei from the state. A state that
  fails this is not a state of the configured atmosphere whatever its
  residual. The shared constraint row now replaces the five corners as the
  step limit; He/H holds at 1.4e-11 over 68 iterates where it had drifted
  to 0.20.
- **The projection conserves the cell's own mass** (N29). `project_elements`
  took the mixture mass from the RESERVOIR metal/hydrogen ratio instead of
  the cell's and created mass at every call (6.4e-3 on the atomic
  candidate); with the cell's own mass the closure is round-off (1.9e-9,
  the file's precision) and the writer-to-loader round trip returns the
  state, which it did not before.

### 22.3 The step control

The trust region is entered only when a species row exists. Its shape at the
end of the series:

- The merit reads every row on its CERTIFICATION scale, and the merit's
  largest row is `||R||`'s row (decision 20 a, N20). Before it the merit and
  the gate were different functionals: the element rows held 0.44 of the
  merit and nothing of `||R||`, accepted steps swung `||R||` by 40 percent
  while the merit moved 1e-4 of itself.
- The element budget is a face of the unknown box (section 21), written once
  per outer iteration and read by all three bound readers.
- A trial whose species unknown lands on a face of its box is held at the
  iterate there instead of being refused with a fixed ascent, which removed
  a period-two cycle of an ascending long step and an accepted short one
  (N22).
- A reduction ratio within 1e-2 of unity is no longer overruled by the
  finite-difference slope test (N7b), which had refused 82 of 91 steps whose
  reduction ratio was 0.999.
- The upper ghost of a transported element or carrier column follows the
  iterate with a zero gradient, the marching path's own rule (N26). Frozen
  at the pre-solve composition it left the outermost element row
  unconstrained from outward.
- The approximate-gradient leg may be admitted by its share of the model
  image (N23); it fires zero times on the present trajectories and the
  image-only variant is worse (N24).
- The preconditioner is a two-sided equilibration of the factorized band,
  with a reservoir-tied column floor (N8a). One empty Jacobian column, a
  carbon unknown at zero on the 1e-20 floor whose probe step underflowed,
  had made the band singular at `dgbcon` 4e-27.
- The Krylov subspace is 40 and stays 40: 80 and 160 reach the same 0.6 of
  the right-hand side (section 21.6); a restarted second cycle is worse on
  both fixtures, and none of the five trust-region resets helps (N21).

### 22.4 What the design achieves, and what it does not

The atomic three-unknown solve converges and certifies
(`wasp_full_newton`: info 0, `||R||` 4.3e-9, log Mdot 13.30).

No solve carrying a species row certifies. The two candidates:

| arm | `\|\|R\|\|` handed back | binding row in the gated window | verdict |
|---|---|---|---|
| atomic element (eight element rows, HD 209458 b) | 2.3e-4 | helium, 2.9e-4 of its scale against 1e-5 | NOT certified |
| molecular carrier (hot Uranus, H2 transported) | 0.25 | H2, 7.3e-2 against 1e-5 | NOT certified |

The obstruction is named and agreed on by three independent items (N21,
N24, N27): **the linear solve**. It reaches 0.58 to 0.99 of its right-hand
side where 0.1 was asked, spending all 40 products, and its Arnoldi image is
15 to 21 percent off the operator. N25 showed the gap is the nonlinearity of
the finite-difference action rather than a loss of orthogonality. The row
scaling is fixed by decision 20 a, so the untried freedom is the COLUMN
preconditioner; the other route is a Newton whose linear tolerance is
actually reached.

### 22.5 Two facts about this solver that change how it is judged

- **The atomic arm is chaotic at the ulp level** (N26c). A one-ulp change of
  the upper ghosts, with no boundary posed, moves `||R||` at iteration 40
  from 3.2e-2 to 1.1, because the Krylov leg's one-digit tolerance decides
  which vector crosses it at the last bits of the residual. A physically
  correct ghost refill (N26b) broke the arm for this reason and was
  reverted; the element operator's caller dependence stays, measured, in
  `src/tests/element_operator`. **Consequence: the arm's `||R||` at a fixed
  iteration is not an acceptance quantity.** What is quoted instead are
  named outcomes: which row binds, which screen refused, the flux spread.
- **Marching relaxes what the Newton cannot**, and much of the reason is now
  known: the marching path refreshes the ghosts every step (N26), and it
  already enforces the shared element sums through
  `limit_to_element_budget` (N4a). The 150 + 150 iteration route that
  reached 5.3e-4 did so because the run marched 2000 steps between the two
  solves, not because the trust region was fresh (N21).

### 22.6 What the solver is judged by

Not by `Resid tol`, which is what the solver was asked for, but by the
certification: one evaluator, one condition per active balance, each with a
tolerance anchored by measurement. The two species-row tolerances were
anchored on 2026-09-10 (decision 22 a, N30,
`docs/certification_tolerance_anchoring_20260910.md`): `1e-5` gating for
`r >= cert_regime_wind_r` = 1.20, the rows below measured, reported and not
gating, because the same balance is a cancellation of advective terms in the
wind and of eddy and settling terms in the layer and reads 2 to 11 times
smaller in the wind on the same state. The band 1.10 to 1.20 is reported and
does not gate: it holds the candidates' binding cell at r = 1.153, where the
element operator's own discretization error is 4.6e-4, so a 1e-5 gate there
would ask the residual to fall 1.5 decades below the error of its own
discrete equation. One reduction, `certification_species_row_gate`, decides a
species row over the gated cells; no gated cell means refused, never a
satisfied row.

What is left of the anchoring: the element flux is conserved to 2.6e-2 in
the wind and not at all in the layer (a factor 39), and the operator's
discretization is order 1.6, so N = 1000 halves the binding row.
