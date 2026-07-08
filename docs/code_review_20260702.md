# EXHALE full-code review (2026-07-02)

**Method.** (1) Full warning build (`gfortran -Wall -Wextra -fcheck=bounds,do,mem`), (2) an
automated scan for local declarations that case-insensitively shadow `global_parameters`
names (the `t0`/`T0` bug class), (3) four parallel subsystem reviews (hydro core;
radiation/ionization; nonlinear/steady solvers; IO/init/post-process), with every substantive
claim re-verified by hand before acting. 61 in-tree modules covered (`wind_ae/` port skimmed
only).

---

## Bugs found and FIXED

### 1. NaN-unsafe cooling-table interpolation → silent out-of-bounds read (LATENT, REAL)
`Cool_coeff.f90`, vector routines `interp_cool_table` / `interp_cool_table_2d`.
A transient NaN temperature during relaxation (observed on WASP-121b at step ~584 under
`-fcheck=bounds`) makes `pos = NaN`; both range tests are false for NaN, so
`k = int(NaN) = -2147483648` indexes `logL` far out of bounds. In production `-O3` builds
this reads arbitrary memory silently (undefined behavior). The **scalar** mirrors already had
the NaN-safe `.not.(pos > 1)` guard (with an explanatory comment) — the vector routines had
been left with the unsafe `(pos <= 1)` form. **Fix:** applied the same inversion to the
vector routines (no-op for finite `pos`). Found only because the review rebuilt with runtime
checks; recommend a periodic `-fcheck` regression run.

### 2. Post-process advection systems used the global `HeH` for the electron density
`System_implicit_adv_HeH{,_TR}.f90`: `xe = xhii + HeH*(xheii + 2 xheiii)` assumed the He/H
ratio equals the input constant in every cell. Correct for legacy runs, but with
`He_diffusion` the local He/H differs from `HeH` by up to ~6× aloft, so the `_adv`
re-ionization (which feeds `EXHALE_transit.py`) misstated n_e. **Fix:** new params slot (15 non-TR / 23
TR) carries the effective He/H — packed as `HeH` when diffusion is off (byte-identical
legacy) and as the local `nhe(j)/nh(j)` when on. Effect on the version comparison: HD 209458b
He 10830 27.2% → 27.1% (trivial); docs updated. Also fixed two wrong params comments
(`xhi_old = nhi/nh`, `xhei_old = nhei/nhe`).

### 3. Robustness guards (no-op in healthy runs)
- `util_ion_eq.f90` (2 sites): `q = Hea_1/q_abs` → `/max(q_abs, 1d-99)` — a fully
  transparent/unilluminated cell gave 0/0 → NaN heating efficiency.
- `sed_read.f90`: `log10(LX_int)`/`log10(LEUV_int)` → guarded — a loaded SED with no bins in
  a band gave `log10(0) = -Inf`.

**Regression after fixes:** production rebuild; HD 209458b and WASP-121b diffusion runs
reproduce all headline numbers (10830: 70.5/27.1, 53.7/53.6; peaks 102/2383; Ṁ 9.04/13.32);
diffusion-off path is arithmetically identical by construction.

---

## Claims investigated and REFUTED (no action)

- **"`Solver` keyword steals `Newton solver:`/`Brent solver:` lines"** — refuted: Fortran
  `index()` is case-sensitive; `'Solver'` (capital S) does not match `...solver:`.
- **`do-subscript` warnings in `ionization_equilibrium.f90:355-367`** — false positives: the
  `j = N+Ng` iteration takes the dedicated first branch, so `j+1` accesses are always in
  bounds.
- **gfortran `maybe-uninitialized` at `EXHALE_main.f90:496/813`** — false positives
  (assigned by `ioniz_eq` before use).
- **`cx_metal_base` threadprivate hazard** — safe: statically initialized per OpenMP rules;
  the runtime `= 5` assignment happens inside the cell-by-cell OMP loop (private to each thread).
- **Missing analytic Jacobians for the TR systems** — by design (documented); MINPACK
  `hybrd1` retained for the triplet/merged systems.

## Shadowing scan results (t0-class)

No active shadowing bug found (beyond the already-fixed `t0`). Hazardous-but-benign namings
recorded: `b0` (WENO smoothness β₀ in `Reconstruction.f90`; recombination array dummy in
`System_HeH_metals`), `mu` (viscosity array in `Apply_BC` vs global hydrogen mass), `g`/`r`
(Brent locals, loop indices), `ng=400` (quadrature size in `excited_hydrogen` vs ghost `Ng`).
The code already defends in one place (`steady_residual` names its residual `Res` explicitly
because of the `r` collision). **Recommendation:** avoid new locals named
`t0,n0,r0,v0,b0,q0,r,n,ng,mu,g,heh` — a rename pass is low-value churn, but new code should
not add to the list.

## Accepted/deferred observations (documented, not fixed)

- **Riemann wave-speed estimate** (`speed_estimate_*`): `Q = p_max/p_min` with no p_min
  floor — safe for healthy states (positivity maintained upstream); add
  `p_min = max(p_min, tiny*p_max)` if pathological states ever appear. Not touched to keep
  the golden-validated hydro core byte-identical.
- **OMP CRITICAL in `PH_heat_HHe`** serializes the result writes in each cell — a scaling
  (not correctness) limitation; measured 2.07×@16 threads previously. Optimization candidate.
- **Cooling/Fe II table clamps** extrapolate flat outside [1e3,1e5] K × [1,1e14] cm⁻³;
  fine for current regimes, revisit for denser/colder bases.
- **Negative intermediate fractions** (`xheii = 1-x2-x3`) during solver exploration —
  standard practice; converged states are physical.
- **Diffusion minor approximations** — see `design_hehe_diffusion.md` §7e (metal mass not
  returned to H, outflow-only top advection, one-step T lag, `He_metal_diffusion` inert
  without `He_diffusion`).
- **Unit-number reuse** (units 1/2/3 across IO modules) — currently safe (sequential
  open/close); prefer `newunit=` in new code.
- 2512 `-Wtabs` warnings (cosmetic), ~40 unused variables (cosmetic).

## Overall assessment

The hydro core, boundary conditions, equilibrium ionization systems, and params
packing/unpacking all verified clean. The two real defects found (NaN-unsafe vector
interpolation; global-HeH in the adv systems) are fixed with regression-verified,
legacy-identical changes. The dominant residual risks are (i) transient NaN temperatures on
hard cases (now harmless at the cooling tables, but worth instrumenting at the source), and
(ii) the documented Phase-1/2 diffusion approximations.
