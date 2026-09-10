# Base-breathing convergence: progress and current state

_Working status note (2026-06-14). Companion to the formal write-up
`docs/base_breathing_investigation.tex`; captures everything since, including
the CETIMB reference findings and the new knobs._

> **Note added 2026-08-10.** Two claims below have been overtaken.
> (i) The Newton/PTC stalls attributed here to the base cells (`worst cell
> j=1/2, momentum`) were a solver artifact, not a property of the base: the
> diagonal scaling floored every momentum scale at `1e-6` of the *base*
> `|rho v|`, and the worst-cell print normalized by the same quantity. With a
> local scale and a watchdog on consecutive failed line searches the affected
> configurations converge: `docs/newton_scaling_and_base_wall.md`. The
> *marching-time* breathing described here is a separate observation and
> stands. (ii) The selectable base BC "planned" at the end of this note is
> implemented (`Base BC: density|pressure [<p_ubar>]`, default density).

## Summary

HD189733b does not fully converge in EXHALE: the time-marching settles into a
bounded limit cycle (the "base breathing"), and the steady (Newton/PTC) solvers
stall. Reading the reference code: **CETIMB** (Koskinen et al. 2013a, 2022),
the model behind **Huang et al. 2023** (WASP-121b, our upgrade target):
revealed both the mechanism and the fixes. Several are now implemented in
EXHALE (Shapiro filter, mass-flux base velocity → new defaults; a gated
viscosity foundation). **No single knob fully converges HD189733b yet**; the
remaining pieces are explicit viscosity/conduction (Phase-2) and the residual
normalization.

> [2026-08-15: all three items of that last sentence have since closed.
> HD 189733 b converges to Newton grade (2026-08-11) after the JFNK line-search
> fix and the beta(tau)/CHIANTI-guarded cooling; explicit viscosity and
> conduction exist as `src/modules/time_step/viscous_conduction.f90` behind the
> `Viscosity:` and `Conduction:` keys; and the residual normalization was
> settled: the volume-weighted norm is the default (`Resid norm: vol`). See
> `docs/hd189_base_checkerboard.md` §10.]

## The problem (localized)

- The steady residual is **purely momentum** (mass ~2e-3, energy ~5e-2, momentum
  ~4.6). A cell-by-cell dump (`EXHALE_RESIDUAL=1` → `output/residual_profile.txt`,
  added this session) shows it **peaks at the dense base** (r≈1.0002,
  R_mom≈180, driven by the gravity source ρg at n0=1e14) and decays outward.
- The convergence metric uses `max|R(j_min:N)|` (excludes the base cells), and
  divides by `max|ρv|` over the escape region, which is **tiny**, so small
  absolute residuals are **inflated** to ~4.6 (a normalization artifact).
- JFNK/PTC, started from a (filtered, smooth) state, still **stall at the base**
  (worst cell j=1/2, r=1.000, momentum; line-search λ→0 = no descent).
- **HD209458b breathes too** (baseline spikes ||R||=216): this is _general_,
  not unique to HD189733b. So the new defaults do not break a clean case.

## CETIMB (reference): how Huang/Koskinen converges

Time-dependent relaxation (like EXHALE/ATES), but with stabilizers EXHALE lacks:

- **Numerics:** operator-split, advection (van Leer, flux-conservative) then
  the Lagrangian terms (viscosity, conduction) via **semi-implicit
  Crank-Nicholson**. Δt=1 s. _Not_ Godunov/HLLC.
- **Explicit viscosity + heat conduction** (Koskinen 2022 B2/B3/B5/B6).
- **A periodic Shapiro (1970) filter.** Koskinen 2013a §2.1.3 states the
  instability explicitly: _"the primary source of the instabilities is pressure
  fluctuations (sound waves) that are not balanced by gravity. We used a
  two-step Shapiro filter periodically to remove numerical instabilities."_
  → this is exactly EXHALE's base breathing.
- **Lower BC at p=1 μbar** (≈3.6e12 cm⁻³, ~30× less dense than EXHALE's n0=1e14):
  fix T, p, composition (from a coupled photochemical model); **base velocity
  from mass-flux continuity ρ₀v₀r₀² = F_c**, not a valve.

## What is implemented in EXHALE (all gated; see input reference below)

| Capability | input.inp key | default | status |
|---|---|---|---|
| **Shapiro low-pass filter** | `Shapiro filter: <eps> <every>` | **off** (`shapiro_eps = -1`) | works for the sound-wave component; made default-on briefly in 2026-06, reverted per the recommendation below |
| **Mass-flux base velocity** | `Base velocity: massflux\|valve` | **off** (`base_v_massflux = .false.`) | ≈valve here, EMA-stabilized; also reverted |
| **Base BC anchor** | `Base BC: density\|pressure [<p_ubar>]` | density (legacy) | the pressure mode planned at the end of this note is implemented |
| Momentum-consistent base p | `Hydrostatic base: True` | off | null for this case |
| Roche IC base blend | (auto, `tidalforce`-gated) | on for tidal | correct improvement |
| Cell-by-cell residual dump | env `EXHALE_RESIDUAL=1` | - | the key diagnostic |
| **Viscosity + heat conduction** | `Viscosity: True`, `Conduction: True` | **off** | complete and validated; ~1e-4 of the base momentum residual, so it does not help (see below) |

Files touched: `src/modules/init/parameters.f90`, `.../files_IO/input_read.f90`,
`.../states/Apply_BC.f90` (shapiro_filter, base-v branch),
`.../time_step/viscous_conduction.f90` (viscosity + heat conduction),
`EXHALE_main.f90` (loop wiring + residual dump). All default-gated; both
`EXHALE.x` and `wind_ae_ic.x` build clean.

## Test results: what each knob does (HD189733b / HD209458b)

- **Shapiro filter**: breaks the violent limit cycle into a smooth descent
  (HD189: 0.15↔2.4 eruptions → smooth ~0.45; HD209: spike 216→65). Confirms the
  breathing _is_ gravity-unbalanced sound waves, and the filter is the right
  damping. _Does not_ converge alone (over-diffusion floor + the base residual
  remain).
- **Mass-flux base velocity**: ≈the valve (the IC already conserves mass flux at
  the base, so v0_massflux≈v0_valve; the base blocker is _momentum/gravity_, not
  mass). The instantaneous F_c fed transient noise into v0 (du blew up to 1353);
  fixed with a **slow EMA** `F_c ← 0.99 F_c + 0.01·mean(ρvr²|[j_min:N])` (du →
  319). Kept default (user's call) for the rho·v·r² overshoot-near-Rp issue.
- **Lower n0 (1 μbar-like, 14→12.5)**: much milder transient (HD189 spike
  233→14.6), confirming the **dense base** drives the violent eruptions; but the
  limit cycle persists and the normalized metric inflates → still no convergence.
- **PTC/JFNK from the IC**: stall (λ→0), worst at the base cell, with/without the
  filter, frozen-base, or smooth valve.

**Conclusion as written (2026-06):** the breathing = (1) sound waves [Shapiro
fixes] + (2) the dense-base momentum residual + (3) the small-ρv normalization.
CETIMB converges via the _combination_; EXHALE's biggest gap is **explicit
viscosity/conduction**.

> [2026-08-15: the last clause did not hold. Viscosity and conduction were
> implemented and measured, and the viscous force is ~1e-4 of the base momentum
> residual: see "(b) viscosity" below. Item (2) was the solver's scaling, not
> the base; item (3) closed with the volume-weighted norm.]

### Knob-combination sweep on HD209458b (cold IC, tidal): `backup/HD209458b_test/`

A controlled 4-knob sweep (S = Shapiro, V = mass-flux base velocity, W = warm
windae IC, N = Newton finish), each model in its own sub-folder, plotted by
`backup/HD209458b_test/plot_models.ipynb`:

The full **2×2 of S × V** on the same cold IC (each run carried to a natural
finish, _adv written), plus the two windae warm-starts:

| model | S | V | W | N | v_out [km/s] | infall% | logṀ | result |
|---|---|---|---|---|---|---|---|---|
| `S0_V0_W0_N0_cold_tidal` | off | off | cold | off | **+6.77** | 44%* | 9.72 | **clean outflow** |
| `S0_V1_W0_N0_cold_tidal` | off | **on** | cold | off | **+8.58** | **1.0%** | 10.50 | **clean outflow (cleanest base)** |
| `S1_V0_W0_N0_cold_tidal` | **on** | off | cold | off | **−13.34** | 99.6% | NaN | **infall** (converged to infall root) |
| `S1_V1_W0_N0_cold_tidal` | **on** | **on** | cold | off | **−13.32** | 99.6% | NaN | **infall → NaN** (~314k steps) |
| `S0_V0_W1_N0_windae_spherical` | off | off | windae | off | +3.68 | 31% | 10.08 | du-stop @count=2 (≈IC) |
| `S1_V1_W1_N0_windae_spherical_force` | on | on | windae | off(Force) | −20.03 | 86% | NaN | infall → NaN |

\* the 44% in `S0_V0` are ~−2 m/s near-hydrostatic base noise, not real infall.

**Key finding (the culprit is the Shapiro filter S, NOT the velocity BC V):**
the 2×2 shows a clean **main effect of S**. With **S off** the cold IC relaxes to
a clean outflow for _both_ V settings (+6.77 / +8.58 km/s); with **S on** it
drifts to **infall** for _both_ V settings (−13.3 km/s, 99.6% of cells negative;
S1_V1 additionally NaNs as V's mass-flux feedback amplifies the infall). **V
(mass-flux base velocity) is benign**, on its own (`S0_V1`) it gives the
_cleanest_ base of all (1.0% vs 44% negative), because setting v₀ from the mass
flux suppresses the ±m/s base oscillation. (An earlier note here blamed V from a
single mid-run snapshot of a run I accidentally killed; the completed 2×2
**refutes** that: V alone never causes infall.)

So the periodic 1-2-1 **Shapiro filter, added (default-ON) to damp the base
breathing, appears to inject a systematic diffusive bias into the momentum that
pushes the clean cold-IC solution off the transonic-wind saddle into the infall
attractor.** The early-transient du _improved_ (216→65), exactly the
"trust-the-number" trap: the residual looked better while the physics was
sliding to infall. Cost: V on also slows convergence ~5× (554k vs 107k steps).

**Implication / action (recommended):** the two new defaults were wrong calls.
- `shapiro_eps` should be **default OFF** (opt-in, breathing cases only): it
  _breaks_ clean cases. This is the important one.
- `base_v_massflux` is physically harmless and cleans the base, but slows
  convergence ~5×; safest as **default OFF** too (opt-in), keeping the fast
  legacy behavior as the default.
Both capabilities stay available; only the defaults revert.

**Pressure-base result (`Pbase_1ubar_S0_V0_W0_N0_cold_tidal`):** the new
`Base BC: pressure 1.0` (n0 derived = 4.9975e12 cm⁻³, ~20× below 1e14) run with
everything else all-off relaxes to a **clean outflow** (+8.36 km/s) with the
**cleanest base of all models, 1/504 negative** (vs 223 for the density all-off),
logṀ 10.40, no NaN, 586k steps. So the less-dense pressure base is a genuine
improvement: it removes the dense-base ±m/s noise and confirms the n0=1e14 base
contributes to the base imbalance, _without_ changing the Shapiro conclusion
(this run has S off). (Not yet tested: S on + pressure base, i.e. whether a light
base lets the Shapiro filter run without driving infall.)

## (b) viscosity: DONE 2026-08-10, and it does not work

The full Navier-Stokes viscous force, its dissipation, and heat conduction are
implemented (`Viscosity: True` / `Conduction: True`, Crank-Nicolson, consistent
between the marching update and the steady residual;
`docs/viscosity_conduction.md`). Calibrated coefficients: Watson et al. (1981)
`κ(T) = 4.45e4 (T/1000 K)^0.7` with `μ(T)` tied to it by the monatomic
Chapman-Enskog relation. Koskinen (2022) never prints a transport coefficient,
and its B5/B6 as printed are not the Navier-Stokes expressions, so the correct
ones are used instead.

**The measurement kills the idea.** On the converged WASP-121 b state the
viscous force is 7.8e-4 of the volume-weighted momentum residual and 1.8e-4 of
it in the worst base cell; `μ` would have to be ~5.7e3 times the physical value
to cancel the base imbalance, and even then the imbalance just moves to the next
cell. The base momentum residual is not a viscous phenomenon. *(The follow-on
reading (that what is left is the lower boundary condition itself) was
withdrawn on 2026-08-10; the base momentum row is satisfiable to within a
2.9 ppm ghost-pressure change, and the JFNK floor was the solver's scaling and
watchdog. See `docs/newton_scaling_and_base_wall.md`.)*

Still open from the original list: revisit the **residual normalization** (don't
divide by tiny ρv).

> [2026-08-15: closed. The volume-weighted residual norm is the default
> (`resid_vol = .true.`); `Resid norm: Linf` selects the old max-over-cells
> norm.]

## Diagnostics / how to reproduce

- Cell-by-cell residual: `EXHALE_RESIDUAL=1 ./EXHALE.x` → `output/residual_profile.txt`
  (columns r, n, v, T, R_mass, R_mom, R_energy).
- Direct steady solve from the IC: `EXHALE_PTC=1` (+ `EXHALE_PTC_NFIX`, `EXHALE_PTC_DTAU0`).
- Disable the new defaults to recover legacy: `Shapiro filter: -1`,
  `Base velocity: valve`.

## HD209458b: the dense base drives _infall_, not just breathing

Re-examining HD209458b (the cleaner case, with a real Koskinen solution to
compare to) exposed a sharper symptom: the relaxed EXHALE solution is **infall
(accretion), not an outflow wind**, v < 0 for r > 1.017 (down to −25 km/s, 459
of 504 cells negative). Diagnosis:

- The **Wind-AE IC is a clean outflow** (0/504 negative, base v ≈ +0.1, outer v
  ≈ +3.7e5 cm/s).
- EXHALE time-marching introduces **infall at the first interior cell**
  (r ≈ 1.0002) **within 2 steps** (154 cells negative), growing with marching
  (355 negative by step ~7000). This is the **dense-base momentum imbalance**
  (ρg at n0=1e14) manifesting as negative velocity: the _same_ root cause as
  the breathing.
- **Newton makes it worse**: it converges to the global **infall root** (459
  negative), a valid steady solution of the equations, but the wrong branch.
  The transonic wind is a _saddle_ (passes through the sonic point); the
  nearby attractors are infall / breeze, so both time-marching and a
  root-finder drift off it.
- The base velocity BC (valve or mass-flux, both give v0 ≥ 0) is **not** the
  cause: the infall is generated in the interior, not the ghost.
- "Wind-AE IC + no Newton" du-stops at count=2 (the flux-flat IC trips the
  du-stop) and so _returns ~the outflow IC_: usable as an outflow solution,
  but not a real relaxation.

### Confirmed: cold IC + no Newton (legacy workflow) → clean outflow

Running `HD209458b/input.inp` (cold IC = no `IC mode`; no Newton = no `Solver`
line; tidal; legacy base, with the new Shapiro/mass-flux defaults turned off)
reproduces the old **clean outflow**: ~106,805 steps (~45 min), Mdot = 10^9.72
g/s, no NaN, a proper transonic wind, v = 0 at the base, +0.07 km/s @1.3 Rp,
+0.95 @2 Rp, **+6.77 km/s @4.2 Rp** (the "223 negative cells" are ~−2 m/s
near-hydrostatic base noise, _not_ real infall); T 1450 → 1167 K. Versus the
windae+Newton infall (−25 km/s, 459 cells).

**The discriminator is the IC's base.** The cold IC uses EXHALE's _own_
hydrostatic base (consistent → it sits in the wind basin → time-marching climbs
to the wind saddle and stays). The Wind-AE IC uses Wind-AE's base (mismatched
with EXHALE's base BC → perturbed off the saddle to the infall attractor →
Newton converges to the infall root). So the **Wind-AE warm-start's base
mismatch** is the issue: the same dense-base / base-BC root as the breathing.
Cold + time-march works today; the windae path needs the base-BC reconciliation.

_Gotcha: cold-IC relaxation is slow, tens of minutes to many hours (even tens
of hours). Run detached; do not poll frequently._

## Selectable base BC: density (legacy) vs pressure / 1 μbar (CETIMB)

_Planned here, implemented since: `Base BC: density | pressure [<p_ubar>]`,
`base_bc_mode` in `parameters.f90`, default density._

The original ATES/EXHALE base BC fixes the **density** (n0, `Log10 lower
boundary number density`) and T (=Teq) at r=Rp. n0=1e14 is very dense → a large
gravity source ρg → the base momentum imbalance that drives both the breathing
and the infall. CETIMB instead anchors at a **pressure** (p = 1 μbar ≈ 3.6e12
cm⁻³, ~30× less dense), deriving ρ from the ideal gas. Lowering n0 to ~12.5 (a
1 μbar-like density) already gave a much milder transient (spike 233 → 14.6),
confirming the direction.

**Plan as written (since carried out):** make the base BC a selectable mode,

- `Base BC: density`, fix ρ=rho_bc(n0) + T=Teq (legacy, default, backward-
  compatible), and
- `Base BC: pressure`, fix p (=1 μbar, CETIMB-style) + T, derive ρ.

The pressure mode should substantially reduce the dense-base momentum imbalance
(hence breathing _and_ infall). Validate against Koskinen's HD209458b
temperature/velocity profiles (`references/Koskinen_2013Icarus_226_1678.pdf`).
The new `examples/12_windae_ic_hd209/plot_output.ipynb` plots the result.
