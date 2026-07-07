# Stage I Extensibility Refactor — Change Log

Date: 2026-06-10. Companion to `ATES/ATES_refactor_plan_from_metal.md` (Path B).
Every step below is a no-physics-change refactor gated by the regression
harness; "PASS" means the gate ran green immediately after the step.

## Regression harness (`regression/`)

- `run_check.sh {golden|check} [case...]` — rebuilds, re-runs each matrix case
  single-threaded (`OMP_NUM_THREADS=1`, fully deterministic), and compares
  `Hydro_ioniz.txt` / `Ion_species.txt` byte-wise against the stored golden,
  excluding `#` header lines.
- Matrix cases (WASP-121b, converges to `du < 1e-3` in ~5 min/case):
  `wasp_full` (He 2^3S on + metals on; exercises the HeITR and metal paths) and
  `wasp_he23off` (He 2^3S off + metals on; exercises the HeITR-off branches).
- `test_roundtrip.sh` + `check_roundtrip.py` — restart loader test using the
  `EXHALE_DUMP_IC=1` hook in `EXHALE_main.f90`, which writes the state exactly as
  loaded and stops. The dump happens *before* the first ionization-equilibrium
  solve: the per-step equilibrium re-solve would otherwise re-derive the metal
  fractions and mask a loader that resets metals to neutral.

## Phase 1 — named species constants + composition module (PASS, bitwise)

- `species_table.f90`: `isp_HI..isp_HeTR` named constants for the fixed
  `f_sp(:,1:6)` H/He/HeITR layout; literal column indices 1–6 eliminated
  across `EXHALE_main`, `energy_semi_implicit`, `set_IC`, `load_IC`.
- New `src/modules/functions/composition.f90`: `get_species_densities`
  (rho, f_sp → all number densities + ne + n_tot), `comp_T_from_p`,
  `comp_p_from_T`. The three near-identical extraction blocks in `EXHALE_main`
  collapse into calls; the electron/total-density policy now lives in one
  place. `eos_include_metals` flag added (default `.false.`, byte-identical;
  hook for a future fully-coupled-metal EOS experiment).

## Phase 2 — output schema + restart preservation (PASS, data-identical + roundtrip)

- `write_output.f90`: schema-2 `#` headers on `Hydro_ioniz*.txt` and
  `Ion_species*.txt`; the species-label line is generated from
  `species_table`, so it stays correct when species are added. NumPy readers
  (`EXHALE_plots.py`, `EXHALE_transit.py`, `examples/exhale_io.py`) are unaffected
  (`np.loadtxt` skips `#` by default).
- `load_IC.f90` rewritten: schema-2 files are read by label mapping
  (order-free), restoring **all** species including the metal ions — metal
  restarts now preserve the ionization state (verified to rtol 1e-12).
  Headerless legacy files keep the historical behavior exactly (H/He read,
  metals to neutral-from-abundance). Per-element rule: an element is restored
  only if all of its ion stages are present in the file.

## Phase 3 — abundance unification (PASS, bitwise)

- The hard-coded per-element metal blocks in `set_IC` (and the fallback in
  `load_IC`) are element loops over `melem_ab(:)` / `melem_i0` / `melem_top` /
  `mion_fsp`. The `X_C..X_Fe` scalars survive only as the input-parsing
  targets that fill `melem_ab`; no downstream physics code reads them.

## Phase 4a — rate dispatch by ion index (PASS, bitwise + roundtrip)

- `species_table.f90`: per-ion index constants `im_CI..im_FeIII` (canonical
  mion order). Note: Fortran identifiers are case-insensitive, so sulfur is
  `im_S_I`/`im_S_II` — `im_SII` would collide with `im_SiI` (neutral Si).
- `Cool_coeff.f90`: `rec_coeff_by_ion`, `ion_coeff_by_ion`,
  `cool_coeff_by_ion`, `cool_coeff_by_ion_scalar` select on the ion index;
  the case lists exist in one place only. The old name-keyed entry points
  remain as thin wrappers (`ion_index_of` + delegate) for compatibility.
- Call sites switched to index dispatch: `util_ion_eq` (rec/ion/cool loops
  and the Fe II density-dependent override) and `T_equation` (post-process
  metal coolant sum). No per-call `trim(mion_name(i))` string comparisons
  remain on the rate path.

## Phase 5 — ionization-solver context (deferred by design)

Step 1 of the v3 plan ("wrap the per-cell coefficients in an explicit
setter") already exists as `set_metal_coeffs` in `System_HeH_metals.f90`.
The remaining steps (local context for the analytic Newton, thread-local or
serial-only MINPACK fallback, active-only unknown mapping) are prerequisites
for cell-parallelizing the ionization loop, not for adding species; they are
deferred until that parallelization is actually scheduled.

## What "add an element" requires after Stage I

1. Rows in `species_table.f90` (indices, element metadata, thresholds,
   flags, `im_*` constant).
2. Rate data: a recombination/collisional-ionization/cooling routine (or
   table) per new ion, plus one `case` line in each `*_by_ion` dispatcher.
3. An abundance entry (`metals.inp` label; `melem_ab` slot via input_read).
4. Nothing in `EXHALE_main`, `set_IC`, `load_IC`, `write_output`, or the
   composition/EOS path — those are all table-driven now.

## Appendix — convergence criterion + local-time-stepping study (2026-06-10)

Two runtime options were added (both opt-in; absent => global-dt path is
byte-identical, regression-gated):
- `Time stepping: Local` — per-cell `dt_j = CFL*dr_j/(|v|+cs)`.
- `Level tol: <val>` — mass-flux LEVEL-stability gate: a converged/stalled
  stop additionally requires the mean `|rho v r^2|` over `[j_min:N]` to be
  unchanged (relative `< lev_th`) across the last `N_stall` steps. `du` is the
  spatial SPREAD of the flux and is blind to a uniform level drift. Default
  `lev_th <= 0` (disabled) until a steady-state reference calibrates it.

**Headline (WASP-121b, He 2^3S + metals) — corrects an earlier wrong reading:**
the `du < 1e-3` stop is unreliable here. `du` oscillates 1e-3..1e-2 and dips
below 1e-3 transiently long before the wind settles, so the original golden
(stopped at step 7,288, raw `log10 4*pi*rho*v*r^2(r_max) = 13.631`) was a
**premature-dip snapshot, not the steady state.** Marching the SAME cold start
much further (to step ~35,085, with the du stop suppressed) the level keeps
rising monotonically — `lev_rel` falls 0.23 -> 2.7e-3 — toward **raw
log10 ~13.71** (`Mdot 13.707` at 35k, still creeping up), with `T_max` 11021 ->
10767 K. Independent paths agree on this deeper state: warm-restart from the
local-dt result settles at 13.716, and naive local-dt itself reaches 13.716 in
**1,702 steps / 69 s**.

Revised conclusions (superseding the first draft of this appendix):
- The true WASP-121b converged state appears to be `raw log10 ~13.71`, reached
  by cold global-dt only after >35k steps. **Naive LTS is therefore NOT
  "answer-shifting"** — it reaches essentially the correct deeper state ~20x
  faster. The earlier "LTS shifts the answer" claim was an artifact of
  comparing LTS against the premature-dip golden; with no trustworthy
  reference, the comparison was meaningless. (Tentative: LTS's exact agreement
  with the global asymptote needs one more confirmation at tighter level tol.)
- The level gate behaves correctly: it refuses the premature du-dip stop and
  keeps cold marching past step 7,288; cold (heading to ~13.71) and warm
  (already ~13.716) agree under it. It is necessary but, on its own, an
  expensive way to find the steady state (tens of thousands of steps).
- Net implication unchanged and reinforced: the principled fix is a direct
  steady-state residual solve (PTC + banded Newton), which yields the true
  fixed point cheaply and gives the "truth" needed to calibrate `lev_th` and
  to re-baseline the regression goldens (whose current du<1e-3 snapshots are
  fine as deterministic refactor references but are NOT the converged wind).

## Steady-residual diagnostic (2026-06-10) — criterion-independent convergence

`EXHALE_RESIDUAL=1` (env hook in EXHALE_main, like EXHALE_DUMP_IC) loads a state,
evaluates the finite-volume steady residual R = du/dt once (reusing
Reconstruct + RK_rhs for dF, S and ioniz_eq for heat,cool; WENO3), and stops.
R(:,1)=dF-S (mass), R(:,2)=dF-S (mom), R(:,3)=dF_E-S_E-(heat-cool). Reported as
max_j|R|/max|u| over [j_min:N]. This is independent of du/dtu/stall.

WASP-121b He23S+metals, three "converged" candidates:

| state                     | R mass  | R mom   | R energy |
|---------------------------|---------|---------|----------|
| golden 13.631 (du<1e-3)   | 9.9e-2  | 4.7e-1  | **3.07e+1** |
| warm 13.716 (LTS+cont.)   | 2.4     | 2.8     | 7.4e-1   |
| cold-35k 13.707 (global)  | **2.1e-2** | **1.6e-2** | **1.1e-2** |

Conclusions (tentative but strong):
- The golden's ENERGY residual is ~30 (heat != cool grossly): `du < 1e-3` is
  blind to energy imbalance because it only measures the mass-flux *spread*.
  The premature-dip reading is confirmed by an independent measure.
- cold-35k is by far the best converged; true Mdot ~13.71. Even it has
  R_E ~1e-2 (not fully steady; lev_rel was still 2.7e-3 climbing) -> a direct
  steady solve is needed to reach R=0 cheaply.
- warm has the right LEVEL (13.716) but large hydro residual: naive LTS reaches
  the level fast but leaves per-cell profile artifacts needing global cleanup.
  (So LTS is "fast but needs polishing", not simply vindicated.)
- ||R|| is the criterion-independent convergence measure the project needed,
  and the building block for PTC + Newton (which drives ||R|| -> 0 directly).
  Next: (i) periodic ||R|| monitor / residual-based stop in the marching loop;
  (ii) PTC + banded-Newton steady solve.
