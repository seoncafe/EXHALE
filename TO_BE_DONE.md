# TO BE DONE

Running list of known limitations and planned improvements.

## OPEN — Wind-AE IC continuation stage C-2 (self-consistent base BCs) — deferred, non-blocking

**What it is.** The in-process Wind-AE IC generator (`IC mode: windae`,
`src/modules/wind_ae/`) currently ships continuation **stage C-1**: it ramps
the system parameters (Ftot -> gravity -> star) from the bundled seed while
holding the seed's *base boundary conditions* fixed. C-2 is the deferred
second stage that **re-converges the base BCs at each ramp step** instead of
holding them fixed, so the lower boundary stays physically consistent with
the target planet as the parameters move far from the seed.

**Why C-1 alone is not enough.** With static BCs, ramping from the mc09
fiducial seed (F ~ 1095) to a strongly-bound, high-flux planet such as
HD 189733 b (F ~ 24000, a ~22x Ftot ramp) leaves the base density/pressure
inconsistent with the new gravity and flux. The relaxation then fails every
step, the adaptive `ramp_var` keeps halving its step, and the ramp stalls
before it finishes even the first variable. This is exactly why
`examples/11_windae_ic/` is shipped as a **deliberate non-converging**
demonstration (the manual points the user to `IC mode: auto` for that case),
while the seed-adjacent `examples/12_windae_ic_hd209/` (HD 209458 b) works.

**The four pieces to port** (from the original Wind-AE Python wrapper; see
`md_backup/wind_ae_fortran_port_plan.md` and `docs/wind_ae_solver.pdf`):

1. **`base_bcs`** — closed-form base state for the target planet: skin /
   effective temperature (`T_skin`, `T_eff` from `F_opt = L*/(4 pi a^2)`),
   the slant-path tau=1 density, and a microbar reference pressure, giving
   `Rmin`, `rho_rmin`, `T_rmin`. ~30 lines, no solver iterations — easiest
   piece.
2. **`ramp_base_bcs`** — try to jump straight to the `base_bcs` target; on
   failure, ramp `Rmin` / `T_rmin` / `rho_rmin` individually.
3. **`converge_Ncol_sp` / `self_consistent_Ncol`** — the sonic-point column
   density BC as a fixed-point loop: `Ncol_sp` = integral of neutral number
   density from the sonic point outward; set -> re-solve -> recompute until
   the mean change is < 8%.
4. **`erf_velocity`** (optional) — the `erf_drop` bolometric heating/cooling
   transition BC; small effect for metals-off warm starts, can be skipped
   initially.

**Companion: HX composition ramp.** Wind-AE solves a fixed H/He mixture; to
hand EXHALE an IC at the run's *exact* He/H ratio, the wrapper also ramps the
composition (`HX`) to the target. Deferred together with C-2 (the IC writer's
density rescaling already absorbs most of the mismatch for warm-start use).

**Cost / difficulty.** Moderate. The risky numerical core — `solvde` plus the
outward `odeint`/`bsstep` integration — is already ported and validated
(windsoln rtol <= 1.9e-11 on three planets; C-1 ramp Mdot matches Python to
5 digits). C-2 adds only fixed-point loops and closed-form formulas on top.
The genuinely tricky part is the **orchestration heuristics** (when, during
the ramp, to re-converge the BCs; the wrapper re-converges after several
consecutive failures) — that is wrapper know-how rather than new physics.

**Validation gate if implemented.** `examples/11_windae_ic/` (HD 189733 b)
converges in-process -> retire its "deliberate non-converging" caveat in
`examples/README.md`, the user manual, and the root `README.md`.

**Suggested order.** `base_bcs` (formula, easiest) -> `converge_Ncol_sp` ->
`ramp_base_bcs` -> (optional) `erf_velocity` -> `HX` ramp.

**Why it is deferred (not blocking).** C-1 plus the IC writer's density
rescaling already covers seed-adjacent hot Jupiters and the warm-start use
case. For the strongly-bound / far-from-seed planets C-2 would unlock,
EXHALE's own `IC mode: auto` already produces a working initial condition. So
C-2's payoff is **code uniformity / completeness** (Wind-AE ICs in-process for
*every* planet), not a capability EXHALE otherwise lacks. An alternative that
avoids C-2 entirely is to **bundle multiple seeds** (one per regime) so every
target is "seed-adjacent" and C-1 suffices.

## DONE 2026-06-13 — Metals in the mass/charge budget (`eos_metals`, default ON)

*(Original entry 2026-06-12, implemented 2026-06-13. Logged as
`docs/Update_EXHALE` §17; the same-day `IC mode: auto` selector is §18.)*

Metals used to be strictly passive trace species: they affected
heating/cooling and the ionization equilibrium, but were excluded from the
mass density, the base density, and the electron/particle densities used
outside the ionization solver. (The only self-consistent place was the
charge balance **inside** the MINPACK systems, `System_HeH_metals.f90`.)

**Implemented** (runtime key `eos_metals 0|1` in `metals.inp`; `1` = new
default, `0` = legacy trace approximation; metals-off runs identical
either way):

- **Base density / composition constants** — `rho_bc` now includes the
  metal mass: `mass_per_H = 1 + 4*HeH + sum(melem_ab*melem_A)` and
  `rho_bc = mass_per_H/(1 + HeH)` (`input_read.f90`); atomic weights in
  `species_table.f90:melem_A`. `n0` keeps its H/He-nuclei meaning, so
  `n_H = n0/(1+HeH)` is unchanged. (Correction to the original entry:
  `v0`, `p0`, `b0` are m_H-based normalizations and never contained the
  mean molecular weight — they were never affected.)
- **Mass density** — `calc_rho` takes an optional `nm` argument and adds
  `melem_A(elem)*nm` (all stages); passed from `ionization_equilibrium`.
- **Electron density** — `calc_ne` takes an optional `nm` and adds
  `stage*nm`; passed everywhere nm is available: `composition`
  (main-loop EOS), `ionization_equilibrium` (opacity/cooling),
  `util_ion_eq` (eval_cool + cooling dump), `energy_semi_implicit`,
  `excited_hydrogen`, `post_process_adv`. The scalar `T_equation` already
  carried metal electrons via `pp_nm_cell`.
- **Particle density** — `calc_ntot` likewise adds the metal nuclei
  (EOS `T = p/(n_tot + n_e)`).
- **Ghost pressure / dp_bc** — the lower-BC ghost pressure is now
  `ntot_bc + dp_bc` with `ntot_bc = (1 + HeH + sum(melem_ab))/(1 + HeH)`
  (so T = T0 holds exactly at the base), and `dp_bc` includes the metal
  electrons. `set_IC`/`load_IC` use `mass_per_H`/`ntot_bc` consistently.

**Validation (WASP-121b regression matrix, 2026-06-13):**
- `eos_metals 0` is **byte-identical** to the pre-change goldens for both
  cases (wasp_full count=7296, wasp_he23off count=7300) — the legacy path
  is untouched.
- `eos_metals 1` (new default): converges cleanly (count 7145/7142, ~2%
  fewer steps); profile shifts are at the expected ~1% mass scale —
  median |dn|~2.4%, |dp|~2.3%, |dT|~0.9% (max |dT| = 132 K at the
  breathing-base region r~1.09, peak rel. 3.8% at r~1.007);
  number flux n*v*r^2 -1.1%, i.e. log10 Mdot 13.31 -> 13.30 g/s — the
  -1.1% number flux offsets the +1.2% mass per particle, so the MASS
  loss rate is essentially unchanged.
- Goldens re-snapshotted under the new default; a repeat run reproduces
  them byte-identically (single-thread determinism intact).

**Still open (small):** the per-cell scalar Brent/Newton `T_equation`
energy balance uses `pp_nm_cell` only in the post-process path; the
time-marching per-cell solve receives ne via `energy_semi_implicit`
(now metal-aware). The `use_2lev_cool` legacy branch remains
AIOLOS-hard-coded (default off).
