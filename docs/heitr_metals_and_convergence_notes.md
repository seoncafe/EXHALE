# He 2$^3$S + metals merge, and convergence working notes

*Working notes (2026-06-08). The code-change sections describe what is in the
tree and verifiable. The convergence sections are **provisional observations and
hypotheses, not settled conclusions** — several interpretations made during this
work were later revised, so the wording below is deliberately tentative.*

---

## 1. Code changes (in the tree, compiled)

These are concrete, verifiable changes to `EXHALE`.

### 1.1 Merged He 2$^3$S (HeITR) + metals ionization solver
Previously the He 2$^3$S triplet and trace metals could not be solved together:
both placed an extra unknown at `x(4)`, and `input_read.f90` allocated the metal
unknowns only `if (thereis_metals .and. .not. thereis_HeITR)`.

- New module `src/modules/nonlinear_system_solver/System_HeH_TR_metals.f90`
  (`ion_system_HeH_TR_metals`): the He/H/triplet rows are taken verbatim from
  `System_HeH_TR` (`x1=HII, x2=HeII, x3=HeIII, x4=HeITR`) and the metal rows from
  `System_HeH_metals`, shifted to `x(5+2*(e-1))`, sharing one electron density
  `ne = HII + HeII + 2HeIII + Σ(X+ + 2X++)` (the neutral triplet is excluded).
  Solved with `hybrd1`.
- `charge_exchange.f90`: a settable module variable `cx_metal_base` (default 4;
  the merged dispatch sets 5) shifts the metal fvec rows for charge exchange.
- `input_read.f90`: `N_eq = 4 + 2*n_melem` when both are active.
- `ionization_equilibrium.f90`: a generalized metal-unknown base index
  (`mbase = 4`, or 5 with HeITR) in the initial guess / solution scatter, plus a
  both-on dispatch branch.
- `Makefile`: the new source added to `SRC`.

The radiation grid, cross sections, photo-heating, IC setup and post-processing
already handled both species; only the equilibrium **solver** needed merging.

### 1.2 Lyman-alpha escape-prob input guard
`input_read.f90`: if `Jlya escape-prob: True` (`jlya_mode=2`) but
`Stellar Lya flux [erg/cm2/s]:` is missing or $\le 0$, the run now stops with an
explanatory message. Rationale: in mode 2, `J_lya = J_int + J_star` with
`J_star ∝ F_Lya_star`; with `F_Lya_star = 0` the stellar beam vanishes and the
`Lya stellar halfwidth/boost` settings become silent no-ops. (Separately,
`use_excited_H` only activates when both `Stellar Teff` and `Stellar radius` are
set.)

### 1.3 Convergence threshold and two-stage reconstruction
- `parameters.f90`: `du_th` is now a runtime variable (was a compile-time
  `parameter`) with default `1.0e-3` — the original ATES-Code-main value.
  (EXHALE had it at `2.0e-2`; see §2.1.) `dtu_th = 1e-8` unchanged.
- The single-stage vs. two-stage choice is set by the `Reconstruction scheme:`
  line. `Reconstruction scheme: PLM+WENO3` selects a **two-stage** run: the
  input line `du_th [PLM,WENO3]: <du_plm> <du_final>` supplies both thresholds
  (`du_th_plm` and `du_th`), and `EXHALE_main.f90` starts in PLM and switches
  `rec_method` to WENO3 once `du < du_plm` (or PLM stalls), then converges at
  `du < du_final`. `Reconstruction scheme: PLM` or `WENO3` is single-stage and
  uses only the first `du_th` value (the second is ignored).
- `CFL` is now runtime + settable via the input line `CFL: <value>`.

### 1.4 Other source-term changes (kept, with caveats)
- `energy_semi_implicit.f90`: `dF_dT = 1 + c_factor*|dC_dT|` (was
  `max(0,dC_dT)`). This adds damping on the falling cooling branch; it is a
  no-op on the rising branch. It removed an oscillation in one test but did
  **not** by itself fix the convergence issues discussed below.
- `Cool_coeff.f90`: `interp_cool_table` now uses monotone PCHIP (C1) rather than
  linear (C0) interpolation. This was implemented to test the hypothesis that
  interpolation kinks drove non-convergence; that hypothesis was **not**
  supported (§2.3). It is retained as a smoothness improvement; a keep-vs-revert
  decision is open. [2026-08-15: decided — PCHIP is kept.]

---

## 2. Convergence — provisional observations

*Everything in this section is tentative. It is recorded as a research log, not
as established fact.*

### 2.1 The threshold had been loosened (observation)
`du` is the relative spatial spread of the mass flux $\rho v r^2$ over the wind
region `[j_min:N]`. EXHALE had `du_th = 2.0e-2`, whereas the upstream
ATES-Code-main uses `1.0e-3`. At `2e-2`, a run flagged "converged" can still show
~2% spread in the supersonic mass flux. We restored `1.0e-3`. This *appears* to
explain why several earlier "converged" profiles did not look flat.

### 2.2 Two-stage PLM$\to$WENO3 (per the ATES docs)
The ATES README recommends PLM from general ICs until $\Delta\dot M/\dot M
\lesssim 0.5\text{--}1$, then a restart with `Load IC` + WENO3. We had been
running single-stage PLM, which *appears* to stall near `du ~ 0.02`. With the
two-stage workflow (now enabled via `Reconstruction scheme: PLM+WENO3`, §1.3)
and strict `du_th`, WASP-121b
(He 2$^3$S + metals) reached `du ~ 1e-3` with a supersonic mass-flux spread of
~0.16% (vs ~2% for single-stage PLM). We interpret this as the two-stage being
necessary for a flat solution, though we have only tested it on a few cases.

### 2.3 A methodological caution (recorded so it is not repeated)
Two issues invalidated several intermediate "non-convergence" readings:
1. An over-short step cap: some weak-wind cases were judged "floored" at ~25k
   steps while `du` was still descending. One case (HD209458b, He+metals) later
   reached `du < 2e-2` near ~50k steps.
2. A run-script bug: a "cap" killed the backgrounding subshell rather than the
   `EXHALE.x` child, so runs continued orphaned and the reported `du` was a
   premature snapshot.

Because of these, the earlier framing that "metal line cooling destabilizes the
wind / causes non-convergence" should be treated as **not established**. It is
*possible* the metals matter, but the evidence we gathered for it was confounded.
(Use `pkill -x EXHALE.x` to stop runs; the subshell kill leaks.)

### 2.4 What we currently see across planets (provisional)
With the two-stage strict-`du_th` pipeline:
- **WASP-121b** (inflated, strongly irradiated): converges quickly and flat
  (~0.16% spread). This case looks solid.
- **HD209458b**: reaches a small `du` (~0.02) that then oscillates rather than
  dropping to `1e-3`. *Tentatively* usable as a quasi-steady state.
- **HD189733b** (strongly bound, $\beta_0 \approx 192$): `du` oscillates around
  ~2 in a sustained limit cycle and does not reach the threshold. This is the
  problematic case.

[2026-08-27: the marching `du` reading of this subsection is superseded. The
route to a converged HD 189733 b is the JFNK finish, not a lower marching
`du`: the base-wall hypothesis of section 2.5 was **refuted** and the actual
causes were the Newton scaling floor and the watchdog, then the line-search
merit, which accepted steps measured on a residual with the WENO3 weights
frozen at the previous iterate. Both are fixed --
`docs/newton_scaling_and_base_wall.md` sections 1-4 and 10,
`docs/Update_EXHALE.md` sections 40-41 -- together with the metal-line
escape probability and the coronal-fit validity floor at a cold base
(section 42). Anything read off this subsection should be re-measured
against those.]

A plausible (but unproven) reading is that convergence difficulty tracks how
weakly driven / strongly bound the wind is, rather than the metal content per se.

### 2.5 Base-breathing hypothesis and the fixes that did NOT work
The HD189733b oscillation *resembles* the known ATES base-breathing behavior,
and the lower boundary (`Apply_BC.f90`, `BC_component_constrho`) hard-pins the
ghost density and pressure to cold reservoir values with a one-way velocity
valve `max(v1,0)` — a configuration that can reflect acoustic waves. This makes
a reflecting-BC limit cycle a reasonable *hypothesis*. However, the fixes we
tried did **not** confirm it and did **not** help:

| Test | Outcome |
| :-- | :-- |
| Lower CFL (0.6 → 0.2) | `du` floor barely moved (~2.8 → ~1.9) |
| Zero-gradient base velocity (drop the valve) | `du` still ~2.3 |
| Velocity from outgoing Riemann invariant $J^- = v-2c/(\gamma-1)$ | did not help HD189733b **and** broke WASP-121b's convergence (drove a large spurious base inflow when cell 1 is hot) |

All three were reverted; the base BC is back to its original form. So the
base-breathing reading is *consistent with* the observations but is **not
demonstrated**, and the simple BC variants we tried are not solutions.

### 2.6 Open directions (not recommendations to do, just options)
If this is pursued further, the candidate routes — both substantial and of
uncertain payoff for a wind this marginal — appear to be:
- a carefully constructed characteristic / non-reflecting (NSCBC-style) base BC;
- a steady-state Newton / BVP solver (which, if a steady solution exists, would
  not orbit a time-marching limit cycle);
  [2026-08-15: implemented — the JFNK steady solve in
  `src/modules/time_step/steady_newton.f90`, selected by `Solver: Newton`.]
- accepting HD189733b as a marginal case and reporting a time-averaged
  quasi-steady state.
- a cheaper, separate idea worth trying for *speed* (not the oscillation):
  local (cell-by-cell) time-stepping, since the global `dt` is currently set by the
  smallest base cell (`eval_dt.f90`: `dt = CFL·min(dr/(|v|+cs))`).
  [2026-08-15: implemented — the input key `Time stepping: Local`.]

See `docs/numerical_methods.md` for a fuller discussion of solver options.

---

## 3. Summary of state (factual)

- He 2$^3$S + metals can now be solved together; this was exercised end-to-end
  on WASP-121b with the full physics stack (triplet + metals + excited-H + Lyα)
  and reached a flat, converged solution.
- The convergence pipeline (`du_th=1e-3`, `Reconstruction scheme: PLM+WENO3`
  two-stage, input-settable `CFL`) is in place.
- HD189733b's marginal-wind oscillation is unresolved and is recorded above as
  an open problem with tentative hypotheses, not a closed result.
  [2026-08-15: resolved since. HD 189733 b reaches a Newton-grade steady state
  (`info = 0`) after the JFNK line-search fix and the beta(tau)/CHIANTI-guarded
  cooling of 2026-08-11; see `docs/hd189_base_checkerboard.md` §10 and
  `docs/newton_scaling_and_base_wall.md`.]
