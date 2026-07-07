# Design sketch: automatic initial-condition selection (`IC mode: auto`)

**Status: IMPLEMENTED (v2 tree) and validated, 2026-06-13.** Phases A and
B of §6 are done; results in §6.1 below. Phase C (retry ladder) remains
open by design.

**Motivation:** Kubyshkina et al. (2018) build a 7000-model grid by
*automatically selecting an initial atmospheric profile for each planet*.
This is what lets them launch winds robustly across the full escape-
parameter range, including the boil-off (low-gravity, nearly
hydrostatic) cases that do **not** launch from a cold static start — the
same regime where EXHALE's breathing base is hardest to converge
(e.g. HD 189733 b). EXHALE already has the building blocks; this sketch
automates the *choice* among them. See `code_comparison.tex` for the
full motivation.

---

## 1. What already exists

`set_IC.f90` implements three IC families, selected by **manual** input
keys (`input_read.f90`):

| Family | Key | What it builds |
|---|---|---|
| cold hydrostatic (default) | (none) | isothermal hydrostatic density + small linear velocity seed |
| transonic | `Transonic IC: True` | steady isothermal-wind Bernoulli profile (`wind_profile`); **falls back to hydrostatic if no interior sonic point** |
| hot-Parker warm seed | `Hot Parker IC: <T_wind>` | cold hydrostatic density + warm T/ionization/velocity overlay (smoothstep `hp_base_rtr`) |

Two facts make automation cheap:

1. **The decision variable is already computed.** `input_read.f90:486`
   sets
   ```
   b0 = (Gc*Mp*mu)/(kb_erg*T0*R0)
   ```
   which is *exactly* Kubyshkina's escape parameter
   Λ = G M_pl m_H / (k_B T_eq R_pl) (here `mu` = m_H, `T0` ≈ T_eq).
   So `b0` is the natural regime selector: **large `b0` = strongly
   bound** (cold hydrostatic works), **small `b0` = boil-off/blow-off**
   (needs a flowing or warm seed).

2. **A sonic-point probe already exists.** `find_sonic(c2, rc, have_rc)`
   bisects the isothermal critical condition φ'(r_c) = 2c²/r_c and
   returns `have_rc` = whether an interior sonic point exists in the
   (possibly Roche-truncated) domain. `transonic_ic` already calls it
   and falls back gracefully when absent.

So the auto-selector is essentially: *run the cheap probes that already
exist, then set `transonic_ic` / `hot_parker_ic` accordingly* — instead
of asking the user to set them by hand.

## 2. The physical regimes

Ordered by escape parameter `b0` ≡ Λ (at the cold base sound speed
c₀² = (ntot_bc + dp_bc)/rho_bc):

| Regime | b0 | Behavior | Right IC |
|---|---|---|---|
| **Strongly bound** | large (≳ 15–25) | modest wind, sonic point far out or outside domain; cold static atmosphere evaporates fine | cold hydrostatic (default) |
| **Transonic-launch** | intermediate | interior sonic point present (often Roche-driven, near L1); cold static start breathes / won't launch | transonic |
| **Boil-off** | small (≲ a few) | weak gravity, nearly hydrostatic but rapid escape; cold static density is too extended and drains onto the base | warm seed (hot-Parker), or transonic at a warm sound speed |

The thresholds are **not** sharp and **must be calibrated** (Section 5);
treat the numbers above as placeholders. The robust *structural* signal
is `have_rc` from `find_sonic`, which is exact (it is the actual sonic
topology of the chosen potential), not a guessed cutoff.

## 3. Proposed decision logic

A new IC mode `auto` that runs at the end of `set_IC` setup, before the
density is built:

```
select_IC_auto:
    ! cold-base sound speed (same c2 the transonic IC uses)
    c2_cold = (ntot_bc + dp_bc) / rho_bc
    call find_sonic(c2_cold, rc, have_rc)

    if (b0 <= b0_boiloff) then
        ! boil-off: cold static density is globally inconsistent ->
        ! warm seed (stable density + warm T/ioniz/velocity overlay).
        hot_parker_ic = .true.
        T_wind_ic     = max(T_wind_ic, f_boil * T0)   ! warm target
        reason = 'boil-off (low b0): hot-Parker warm seed'

    else if (have_rc) then
        ! interior sonic point exists -> a flowing transonic IC both
        ! launches the wind and matches the eventual topology.
        transonic_ic = .true.
        reason = 'interior sonic point: transonic IC'

    else
        ! strongly bound, no interior sonic point -> default cold
        ! hydrostatic; the wind launches from radiative heating.
        reason = 'strongly bound: cold hydrostatic IC'
    endif

    force_start = (transonic_ic .or. hot_parker_ic) .and. .not. do_only_pp
    write(*,*) '(select_IC_auto) ', trim(reason), '  (b0=', b0, ')'
```

Key properties:

- **Conservative / non-regressing.** The `else` branch is the current
  default, so high-`b0` planets (the bulk of validated cases: WASP-121b
  Case A/B, tutorial, HD 209458 b) are unchanged. Auto is **opt-in** via
  `IC mode: auto`; the existing manual keys keep working and override.
- **Uses only exact probes** (`b0`, `find_sonic`) — no new physics, no
  new tables.
- **Self-healing.** If `transonic_ic` is chosen but `wind_profile`
  finds no usable root at run time, the existing fallback to hydrostatic
  still fires — auto never makes things worse than manual.
- **Logged.** The chosen branch and `b0` go to stdout (and should go to
  `write_setup_report.f90`) so the choice is auditable in `run.log`.

### 3.1 Refined decision tree (v2 — recommended)

Working through the numbers shows the `b0 <= b0_boiloff` branch of §3 is
**unnecessary**: the cold sonic-point probe alone separates the regimes,
because for a point-mass potential the cold sonic radius is
`r_c = b0/(2 c0²)` (dimensionless), so *low* `b0` automatically pulls the
sonic point inside the domain, and deep-RLOF cases are caught by the
Roche topology (`dφ/dr → 0` at L1) in the *actual* potential that
`find_sonic` bisects. Representative point-mass estimates
(c0² ≈ 0.75 for H/He):

| Case | b0 | r_c(cold) | r_max | interior? | v2 choice |
|---|---|---|---|---|---|
| HD 209458 b (classic HJ, 10 Rp) | ~81 | ~54 Rp | 10 | no | cold hydrostatic |
| HD 189733 b (breathing case) | ~179 | ~119 Rp | 10 | no | cold hydrostatic |
| WASP-121b (RLOF, Roche) | ~48 | (Roche: interior via L1) | 1.57 | yes | transonic |
| generic boil-off sub-Neptune | ~9 | ~6 Rp | 10 | yes | transonic |

So the v2 tree is a **single exact probe, no tunable threshold**:

```
if (have_rc at cold c2)  ->  transonic IC
                             (catches deep-RLOF via Roche topology,
                              boil-off via low b0, and supersonic-from-
                              base extremes where find_sonic returns
                              rc = base)
else                     ->  cold hydrostatic (classic EUV-heated wind;
                              launches fine from a cold start)
```

with `b0` **demoted to a logged diagnostic** (still worth printing — it
labels the regime for the user) and the **hot-Parker warm seed removed
from the auto path**: the IC benchmark showed it gives no speedup when
the cold start works, and a "warm sonic point exists" criterion would
fire on essentially every hot Jupiter. The warm seed remains available
manually, and as the escalation step of the optional retry ladder
(§Phase C below). Note this also means auto does *not* claim to fix the
HD 189733 b breathing stall — no IC does (benchmarked); that case needs
a different mechanism entirely.

## 4. Implementation steps

1. **Input key.** Parse `IC mode: <cold|transonic|hot_parker|auto>` in
   `input_read.f90` (alongside the existing `Transonic IC` /
   `Hot Parker IC` keys, which remain as explicit overrides). Default
   `cold` to preserve current behavior. Add an `ic_mode` enum/string to
   `parameters.f90`.
2. **Selector.** Add `select_IC_auto` to `set_IC.f90` (it needs `b0`,
   `find_sonic`, `rho_bc`, `ntot_bc`, `dp_bc`, all already in scope).
   Call it from `set_IC` when `ic_mode == auto`, *before* the
   density-building branch, so it just sets the existing
   `transonic_ic` / `hot_parker_ic` flags and the rest of `set_IC` runs
   unchanged.
3. **force_start.** Move the `force_start` assignment (currently in
   `input_read.f90:415`) so it is re-evaluated after `select_IC_auto`
   (or recompute it inside the selector), since auto may flip the flags
   after input parsing.
4. **Report.** Emit the chosen branch + `b0` in `write_setup_report.f90`.
5. **Docs.** Document `IC mode: auto` in `EXHALE_user_manual.tex` (§input)
   and cross-reference from `EXHALE_BC_and_IC.tex` §5.

No change to the physics, the BC, or the golden outputs for any planet
that does not opt into `auto` (or that auto routes to the cold default).

## 5. Calibration and validation plan

*(Written for the v1 tree of §3. Under the recommended v2 tree (§3.1)
there are no thresholds to calibrate — `find_sonic` is exact — and the
validation reduces to the regime checks below, restated concretely as
the phase gates of §6. Kept for the record.)*

The v1 thresholds `b0_boiloff` and `f_boil` would need calibration;
`find_sonic` does not. Proposed validation:

1. **Regression invariance.** With `IC mode: cold` (default) and the
   existing manual keys, every golden is byte-identical. With
   `IC mode: auto` on the high-`b0` cases (WASP-121b A/B, HD 209458 b),
   auto must route to cold hydrostatic and reproduce the goldens.
2. **Transonic regime.** WASP-121b Case D (sonic point ≈ L1) and similar
   deep-RLOF cases: auto should select transonic and launch, matching
   the manual `Transonic IC: True` result.
3. **Boil-off regime.** Pick a low-`b0` planet (low gravity, high T_eq;
   a sub-Neptune near its host). Confirm auto selects the warm seed,
   that it launches without the breathing-base divergence, and that the
   steady-state `Mdot` is IC-independent (rerun from a different IC and
   compare — Salz's and Kubyshkina's stated property).
4. **Threshold scan.** Vary `b0_boiloff` around the transition and check
   that the final steady state (not the path) is insensitive — i.e. the
   selector affects convergence robustness/speed, not the answer.

## 6. Phased implementation plan (concrete)

Follows the project's phase-gate workflow: each phase has a hard
validation gate; do not proceed past a failed gate. Call-order fact
underpinning the placement: `input_read` -> `init` (`define_grid` ->
`set_gravity_grid` -> `set_IC`/`load_IC`), and `force_start` is only
read later in the main loop -- so a selector inside `set_IC` has the
grid, the potential, `find_sonic`, and all composition globals
available, and may still set `force_start`.

### Phase A — plumbing + selector (no behavior change by default)

1. `parameters.f90`: add `integer :: ic_mode = 0`
   (0 = cold, 1 = transonic, 2 = hot_parker, 3 = auto) next to the
   existing `transonic_ic`/`hot_parker_ic` declarations.
2. `input_read.f90` (optional-key loop, alongside the
   `'Transonic IC'` / `'Hot Parker IC'` branches): parse
   `IC mode: <cold|transonic|hot_parker|auto>` via
   `index(line,'IC mode')` + `get_word(line, 3)`. Mapping: the explicit
   words set the corresponding legacy flags (synonyms for the existing
   keys); `auto` sets `ic_mode = 3` only. The legacy keys keep working
   and take precedence over `auto` (explicit beats automatic).
3. `set_IC.f90`: new module subroutine `select_IC_auto()` implementing
   the v2 tree (§3.1):
   - `c2_cold = (ntot_bc + dp_bc)/rho_bc`
   - `call find_sonic(c2_cold, rc, have_rc)`
   - `have_rc` -> `transonic_ic = .true.`; else leave flags off.
   - set `force_start = (transonic_ic .or. hot_parker_ic) .and. .not. do_only_pp`
     (mirrors input_read.f90:415, which cannot see auto's decision).
   - log: chosen branch, `b0`, `rc`/`have_rc`, one line, prefixed
     `(select_IC_auto)`.
   Call it at the very top of `set_IC` when `ic_mode == 3`, before the
   `hot_parker_ic / transonic_ic` branch at set_IC.f90:43.
4. `write_setup_report.f90`: report the IC family actually used (and
   `auto` provenance if `ic_mode == 3`).
5. Doc: one paragraph in `EXHALE_user_manual.tex` §input-file (optional
   keys), cross-ref in `EXHALE_BC_and_IC.tex` §5.

**Gate A (regression invariance):**
- Without any new key: full regression matrix byte-identical
  (`run_check.sh check`).
- With `IC mode: cold`: byte-identical.
- With `IC mode: auto` on a no-interior-sonic-point case
  (HD 209458 b-like spherical tutorial): selector must route to cold
  hydrostatic and the run must be byte-identical to the no-key run.

### Phase B — behavior validation on the three regimes

1. **RLOF / transonic**: WASP-121b Case D (`input.inp` in `WASP-121b/`,
   r_esc = 1.20): `IC mode: auto` must select transonic
   (have_rc via L1) and reproduce the manual `Transonic IC: True`
   launch behavior (wind launches; subsonic-at-L1 geometry caveat
   unchanged).
2. **Boil-off**: construct a low-`b0` test planet (target `b0` ~ 5-10,
   e.g. ~6 M_earth / ~5 R_earth / T0 ~ 1100 K, spherical 10 Rp;
   verify the printed `b0` lands in range). Auto must select transonic
   (interior cold sonic point). Check: launches and converges where the
   cold start breathes/stalls, and steady `Mdot` is IC-independent
   (rerun the same planet with `IC mode: cold` long enough, or with the
   manual warm seed, and compare `Mdot` to <~1%).
3. **Classic**: WASP-121b Case A/B and tutorial under `auto`:
   Case A/B are Roche-mode -> auto may legitimately pick transonic;
   verify the *steady state* matches the cold-start golden `Mdot`
   (path differs, answer must not). Tutorial (spherical) routes to
   cold -> byte-identical.

**Gate B:** all three regimes behave as designed; any case where auto
picks a *different* family than manual best practice must match the
manual result in `Mdot` to <~1%.

### Phase C (optional, separate decision) — retry ladder

Escalation on failure instead of prediction: if a run NaN-crashes or
hits the stall detector with no wind launched, restart once with the
warm seed (`hot_parker_ic`, T_wind = 1e4 K), then once from the last
good state. This is where the hot-Parker family re-enters, and where
Murray-Clay-style continuation (§7) would slot in if ever needed.
Implementation cost is higher (main-loop restart logic); decide after
Phase B data shows whether any real case still fails.

### Edge cases to handle explicitly

- `Load IC? True`: `set_IC` is bypassed (init.f90 calls `load_IC`), so
  the selector never runs -- correct by construction; document it.
- `do_only_pp`: selector may run but `force_start` stays false
  (already in the formula); IC is irrelevant for PP-only runs.
- Spherical mode with user `r_max`: `find_sonic` brackets over the
  actual domain, so a sonic point just outside `r_max` correctly routes
  to cold -- no special case needed.
- `find_sonic` returning `rc = base` (supersonic from base): transonic
  branch handles it (`wind_profile` takes the supersonic root).

### Estimated scope

~60-80 new Fortran lines (selector + key parsing + report), no new
files, no physics change, goldens untouched by default. Phase A+B is a
single working session plus the boil-off calibration run.

### 6.1 Validation results (2026-06-13)

Implemented exactly as planned (`ic_mode` in `parameters.f90`, `IC mode:`
key in `input_read.f90` with the legacy keys taking precedence,
`select_IC_auto` at the top of `set_IC`, IC family in the setup report).
All gates run:

| Gate | Setup | Result |
|---|---|---|
| A-1 | no new key, full regression matrix | **PASS** — byte-identical |
| A-2 | `IC mode: cold` on wasp_full | **PASS** — byte-identical, same step count (7145) |
| A-3 | `IC mode: auto` on spherical tutorial (HD 209458 b-like) | **PASS** — selector logs b0 = 83.2, no interior cold sonic point, routes to cold hydrostatic; `IC_dump.txt` byte-identical to the no-key run |
| B-1 | `IC mode: auto` on WASP-121b Case D | **PASS** — interior cold sonic point at r_c = 1.297 Rp (= the L1 crossing), b0 = 35.2 -> transonic, identical to the manual key |
| B-2 | constructed boil-off planet (0.0189 MJ / 0.446 RJ / 1100 K, spherical 10 Rp) | **routing + launch PASS, speed verdict mixed** — b0 = 8.5, r_c = 5.2 Rp (predicted 8.3 / 5.5) -> transonic; the IC starts at residual ~0.19 vs ~14 from cold (70x closer) and the wind launches (13 km/s, T ~ 5400 K). However the COLD control also launched and converged first (138k steps, level-stable, log Mdot = 11.20), while the auto run relaxed monotonically toward the same flux level from above (5.2 -> 3.2 -> toward 1.5e15) without finishing in the 50-min budget: the cold-c2 isothermal Parker profile *overestimates* the flux of this weakly-heated wind, so the transonic start must shed mass first. Tentative reading: for this mild boil-off case auto is not harmful (same attractor) but not faster either; the launch-failure regime that motivates auto is likely more extreme (lower b0 / Roche-truncated). |
| B-3 | `IC mode: auto` on wasp_full (Roche) | **PASS at the fixed point** — auto picks transonic (r_c = 1.548 ~ L1, b0 = 49.9) and reaches the du-stop in 5876 steps vs 7145 cold (-18%). The du-stop states differ by +4.8% in flux (both with rho*v*r^2 flat only to 1-2% — the known false-convergence spread), but Newton-finishing BOTH states lands on the SAME fixed point: flux ratio 0.9996 (-0.04%), flatness 0.005% each, median profile difference 0.003%. |

Two byproduct findings worth recording:

1. **The du = 1e-3 stop carries a path-dependent spread of up to ~5%**
   in the mass flux (B-3): two legitimate runs of the same physics from
   different ICs both "converge" by the du criterion at states 4.8%
   apart, each with percent-level flux non-flatness. The Newton finish
   collapses both to the same fixed point (0.005% flat). For
   quantitative Mdot work, finish with `Solver: Newton` (or `Resid
   tol`); do not trust bare du-stops.
2. **The transonic IC can speed up Roche-mode runs** (-18% steps on
   wasp_full to the du-stop) — unlike the warm seed on spherical runs
   (no speedup, see `initial_condition_benchmark`). Plausibly because
   in Roche mode the transonic IC has the correct outflow topology from
   step 0. The effect does NOT carry over to the spherical boil-off
   test (B-2), where the cold-c2 Parker flux overshoots the true wind
   and relaxation is slower than from cold — the speedup appears to be
   specific to cases whose steady wind is close to the cold isothermal
   transonic solution (deep-RLOF).

## 7. Alternative / complementary strategies from the literature

The regime-switch selector above is the cheapest robust improvement, but
two other reference codes suggest complementary ideas worth recording:

- **Continuation / homotopy (Murray-Clay et al. 2009).** Rather than
  picking one IC, they reach a transonic solution by *solving
  successively more complicated problems*: start from an isothermal wind
  with no photoionization, then add photoionization, heating, cooling,
  and tidal gravity "one by one ... in a diluted form first and
  gradually strengthened to full amplitude." For EXHALE this maps to a
  **physics ramp**: launch the wind with a simplified source (e.g. a
  fraction of the heating, no metal cooling), then ramp the full
  microphysics in over the first phase of the run. This is a proven way
  to reach the transonic branch for cases that breathe from a cold
  start, and it composes with the SER/PTC ramp the JFNK solver already
  has. It is more invasive than the IC selector but attacks the same
  non-launching cases from the source-term side.

- **Sonic-point anchoring (Murray-Clay et al. 2009).** They stress that
  "for every transonic wind there are an infinite number of breeze
  solutions," and codes that do not enforce the Parker critical-point
  conditions "found ... breezes instead." EXHALE's transonic IC already
  builds the solution *through* the sonic point (`find_sonic` +
  `wind_profile`), so the selector inherits this; the deeper option is to
  have the JFNK steady-state solver itself impose the sonic-point
  conditions, which is a larger change (noted in `code_comparison.tex`
  §"A different solver philosophy").

- **Per-grid-point auto IC (Kubyshkina et al. 2018).** The original
  motivation: their automatic per-planet initial-profile selection is
  exactly the `select_IC_auto` logic of §3, generalized to a parameter
  grid. If EXHALE is ever used to mass-produce a model grid, the selector
  should be driven directly by the grid coordinates (`b0`, T_eq, orbital
  separation) rather than re-probed per run.

- **Wind-AE warm-start IC (tested 2026-06-13).** The MC09 successor
  Wind-AE (Broome et al. 2025; in-tree at `wind-ae-main/`, builds and
  runs on this machine, ~14 s per solve) provides converged steady BVP
  solutions that can be interpolated onto the EXHALE grid and loaded via
  `Load IC` (converter: `src/utils/windae_to_exhale_ic.py`). Hands-on
  result on
  HD209458b (H/He, spherical, 10 Rp): the loaded state is flux-flat *by
  construction*, so plain marching trips the `du` stop **in 3 steps — a
  false convergence** (the most extreme instance of the du-stop trap of
  §6.1); with `Solver: Newton` the state instead relaxes smoothly toward
  EXHALE's own attractor (flux 10.78 → 10.68 dex, spread 2.5e-4 →
  6.6e-2 → 2.1e-2 at 146k steps, no NaN), while the **cold control
  NaN-crashed at the breathing base** (step 172k, du still 3.7). So the
  Wind-AE IC is a genuine robustness tool for the hardest
  (weakly-driven) cases, but such runs must always be Newton-finished.
  Wind-AE also showed a steady HD189733b solution exists (its BVP
  converges; the EXHALE limit cycle is likely numerics), and its
  continuation stalls for the WASP-121b near-RLOF corner — it
  complements, not replaces, the Roche machinery. Details:
  `code_comparison.tex` §"Wind-AE ... verified in-tree".

## 8. Open questions / caveats

- **Boil-off needs more than an IC.** Kubyshkina's robustness in
  boil-off also comes from the **molecular-H₂ lower boundary** (a deeper,
  nearly-hydrostatic base), not just the IC. EXHALE's atomic H/He
  base may still breathe in the most extreme low-gravity cases even with
  the right IC. A warm-seed IC is necessary but possibly not sufficient;
  a deeper/molecular base is a larger, separate change.
- **`b0` uses the base T₀, not a wind temperature.** For boil-off the
  relevant sound speed is the warm/ionized one; the selector could
  optionally probe `find_sonic` at the warm `c2` (as the hot-Parker IC
  already does) to decide transonic-vs-warm more sharply.
- **Convergence is momentum-gated.** The IC benchmark showed warm seeds
  give no step-count speedup for already-launching cases — so the payoff
  of auto is *robustness* (launching the hard cases at all), not speed
  for the easy ones. Set expectations accordingly.

## 9. Summary

The selector is small and low-risk: it reuses `b0` (= Λ, already
computed) and `find_sonic` (already called) to set the existing IC flags
automatically, defaults to the current cold-hydrostatic behavior, is
opt-in, self-healing, and logged. The main payoff is robust launching in
the transonic and boil-off regimes — EXHALE's weakest area — closing
the one clear capability gap vs. Kubyshkina+2018 identified in
`code_comparison.tex`. The deeper boil-off fix (molecular/deeper base)
is noted as a larger, separate follow-up.
