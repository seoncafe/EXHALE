# Design: He/H diffusive separation in EXHALE (Phase 1 + Phase 2)

**Status:** Phase 1 and Phase 2 implemented and validated (2026-07-02). Flag `He_diffusion`
default OFF (opt-in). See §7d/§7e for results.
**Scope:** Phase 1 = He element relative to H (§2–§7). Phase 2 = ambipolar-corrected
settling (P2b), thermal diffusion α_T (P2c), and trace-metal diffusion element by element (P2d);
self-consistent n_H (P2a) attempted but abandoned (§7e).

---

## 1. Why this is non-trivial in EXHALE

EXHALE is single-fluid: hydro evolves total `rho`, momentum, energy with one bulk velocity
`v`. Composition is stored as ion mass-fractions `f_sp` and re-solved each step by **local**
ionization equilibrium, which conserves the element ratio. There is currently **no transport
of the element ratio at all** — not even advection — so He/H is globally fixed at the input
`HeH`. To get diffusive separation we must add a genuine **He-element transport equation**
(advection at bulk `v` **plus** a diffusive drift), then feed the resulting radius-dependent
`nhe(j)/nh(j)` into the ionization solver (which already works cell by cell).

---

## 2. Governing equation (He element)

Track the He **mass fraction** of the gas, `Y ≡ rho_He/rho` (equivalently the element ratio
`f ≡ n_He/n_H`). Evolve

```
∂(rho Y)/∂t + (1/r²) ∂/∂r[ r² (rho Y v + F_diff) ] = 0
```

- Advective term uses the **same mass flux** `rho v` the hydro already computes → consistent
  and bounded; if `Y` is uniform it reduces to continuity × `Y` (no spurious separation).
- No chemical source: ionization moves He between stages but conserves the He **element**.
- `F_diff` is the diffusive mass flux of He relative to the bulk (below).

### Diffusive drift (binary He-in-H, diffusion approximation)

Relative drift of He w.r.t. H (Banks & Kockarts form), with `f = n_He/n_H`:

```
w_d = -D · [ (1/f) ∂f/∂r  +  (Δm · g)/(k T)  +  α_T · ∂lnT/∂r ]
```

so the diffusive **number** flux is `Φ_He = n_He · w_d` and `F_diff = m_He · Φ_He`, i.e.

```
Φ_He = -D · n_H · [ ∂f/∂r  +  f (Δm g)/(kT)  +  f α_T ∂lnT/∂r ]
```

- **Δm = m_He − m_H = 3 amu.** Justification for Phase 1: in the highly-ionized wind the
  ambipolar field lifts ions; for He⁺⁺ in an H⁺ background the net settling force works out to
  ≈ 3 m_H g — the **same** as the neutral-atom value — so the neutral Δm=3 is a good first
  approximation and lets us defer the explicit ambipolar field to Phase 2.
- **Gravity term** drives He settling (He/H falls with altitude); **gradient term** opposes;
  **advection** (bulk `v`) drags He up. Their competition sets the profile — He/H declining
  with altitude, as in Taylor et al. (2025: 8%→2.5%) and Xing et al. (2023).
- **α_T (thermal diffusion): 0 for Phase 1** (small; a documented refinement).

### Diffusion coefficient D (He in H)

Binary molecular diffusion of a minor species in a background, **Banks & Kockarts (1973,
*Aeronomy*)** form:

```
D(He-H) = 1.52e18 · (1/m_H + 1/m_He)^{1/2} · T^{1/2} / n_tot        [cm² s⁻¹]
        = 1.70e18 · T^{1/2} / n_tot        (m in amu, T in K, n_tot in cm⁻³)
```

Phase-1 default uses this **neutral** binary value; the ionized-plasma (resonant/Coulomb)
correction is a documented caveat and a Phase-2 item. The coefficient sets **where**
separation turns on via the timescale ratio `τ_D/τ_v = v H / D` (Koskinen 2013): near the
base `v` is small ⇒ diffusion dominates ⇒ He settles (He/H drops with altitude); aloft
`v` is large ⇒ advection freezes the composition (He/H flattens) — reproducing the
Taylor/Xing profile shape.

---

## 3. Numerics

- **Conservative finite volume** on the existing grid (`r`, `r_edg`, `dr_j`); fluxes at faces
  `r_edg` with `r²` areas, update cell centers.
- **Advection of `Y`** with the hydro face mass-flux (reuse the hydro's `rho v` at faces;
  upwind on `Y`) — consistent multi-species advection, keeps `Y ∈ [0,1]`.
- **Diffusion** term explicit, **sub-cycled** under its own stability limit
  `dt_diff ≤ 0.5 dr²/D` (D can be large where `n` is low → a handful of sub-steps).
- **Operator split**, placed in the relaxation loop **after** ionization equilibrium /
  before the next hydro step (EXHALE_main, near the energy update). One diffusion+advection
  update of `Y` per relaxation iteration, sub-cycled internally.

## 4. Coupling to ionization / metals / rho

Diffusion conserves total mass locally (the implied H flux is `F_H = -F_He`, i.e. rho owned
by hydro is preserved). Each step, after updating `Y(j)`:

```
n_He(j) = rho(j) Y(j) / m_He
n_H(j)  = rho(j) (1 − Y(j)) / (m_H + Σ_metal A_metal m_metal)   ! metals slaved to H
```

then rescale `f_sp` He-stages to sum to `n_He`, H-stages to `n_H` (preserving the within-
element ionization split), metals `= A_metal n_H`; ionization equilibrium re-solves with the
new `nh(j)`, `nhe(j)` in each cell (interface already accepts these independently).

## 5. Boundary conditions

- **Base (lower):** fixed composition `Y = Y(HeH)` (Dirichlet) — the homopause-anchored
  reservoir, as in Koskinen (fixed base composition) and Taylor (μbar base).
- **Top (outer):** zero **diffusive** flux (`Φ_He = 0`); advection carries He out with the
  wind (outflow). This lets He/H reach its wind value aloft without an artificial pile-up.

## 6. Integration points (files)

- `input_read.f90` / `parameters.f90`: new flag `He_diffusion` (+ optional `Kzz`/`alpha_T`
  hooks, unused in Phase 1); **default OFF** (goldens unchanged; metals-off byte-identical).
- new module `species_diffusion.f90` (generic: given element mass, abundance array, D(T),
  return updated abundance) under `src/modules/` (radiation/ or a new `transport/`).
  [2026-08-15: it landed as `src/modules/functions/species_diffusion.f90`.]
- `EXHALE_main.f90`: call the diffusion substep in the relaxation loop; carry a persistent
  `Y(1-Ng:N+Ng)` (He mass fraction) array across iterations.
- `composition.f90` / `ionization_equilibrium.f90`: derive `nh(j)`,`nhe(j)` from `rho` and
  `Y(j)` instead of the global `HeH`; rescale `f_sp` accordingly.
- `set_IC.f90`: initialize `Y(:)` from `HeH` (uniform).
- `write_output.f90` / `write_setup_report.f90`: already write species in each cell; optionally
  add a He/H(r) column and report the flag.

## 7. Validation gate (must pass before default-on)

1. **Flag OFF ⇒ byte-identical** to current goldens (no path change when disabled).
2. **He/H declines with altitude** on HD 209458b (target Taylor: ~8%→2.5%) and WASP-121b;
   qualitatively matches Xing (2023) fractionation.
3. **Timescale check:** separation appears where `τ_D/τ_v = vH/D ≳ 1`, negligible where ≪1.
4. **Conservation:** total He number flux `r²(n_He v + Φ_He)` constant at steady state;
   `rho` unchanged by the diffusion operator (mass conserved).
5. **He 10830** recomputed (EXHALE_transit.py) — expect a weaker line consistent with reduced upper-
   atmosphere He (the physical point of the whole exercise).

## 7b. Validation status (2026-07-01) — GATE NOT PASSED

The implementation is **numerically stable** (implicit tridiagonal solve,
no NaN on HD 209458b; flag OFF is byte-identical) and **produces diffusive separation**
(He/H falls with altitude). **But it over-separates**: on HD 209458b (He23S on, metals on,
15k steps) He/H holds ~8% to 1.2 R_p then collapses to ~3e-5 by 1.4 R_p and →0 above —
versus Taylor's mild 8%→2.5% plateau. The separation also **deepens monotonically with
step count** (6k: 0.06 at 1.4 R_p; 15k: 3e-5), i.e. it is not converging to a mild balance
but draining He.

Diagnostics: at 1.4 R_p the settling velocity w_s = D·G ≈ 1.2e3 cm/s is ~160× smaller than
the wind v ≈ 1.9e5 cm/s (τ_D/τ_v = vH/D ≈ 460), so **physically advection should dominate
and keep He nearly well-mixed** — the simulated total depletion is therefore too strong.
Adding eddy diffusion (K_zz = 1e9 cm²/s) did not help: the homopause (D = K_zz) sits below
1.2 R_p, so molecular settling still drains everything above it.

Leading hypotheses for the over-drain (to resolve before enabling by default):
1. **Discretization inconsistency:** advection is applied in non-conservative
   (material-derivative, upwind) form while diffusion is a conservative face-flux; the two
   are not flux-consistent, which can act as a spurious He sink. A fully conservative He
   transport (evolve n_He with r²(n_He v + Φ_He) in one divergence, derive n_H from ρ) is
   the most likely fix.
2. **Near-base stiffness + boundary layer:** the cold dense base has a tiny He scale height;
   settling there may drain the column faster than advection refills it. A physically higher
   homopause (larger/height-dependent K_zz) and/or ambipolar relief in the ionized region
   (deferred Δm approximation) may be required.

Until this is resolved the flag stays **default OFF**; normal runs are unaffected.

## 7c. Conservative reformulation (option A) — status 2026-07-01

The transport was rewritten to be **fully conservative**: the He number density
`n_He` is evolved with advection (bulk `v`, upwind) and the diffusive flux combined in a
**single conservative divergence** of the total face flux `J = n_He v − D n_H(∂f/∂r + fG)`
(`n_H` lagged), solved implicitly (tridiagonal). This replaced the earlier
non-conservative material-advection form.

Findings:
- **Settling OFF (G=0):** He/H stays ~0.8×HeH out to ~2 R_p (no runaway drain) — the
  advection+diffusion transport is essentially well-behaved. So the conservative form fixed
  the gross transport inconsistency.
- **Settling ON:** still **over-separates** — He/H collapses across a thin transition
  (~1.25 R_p) to ~0 above. A near-base numerical pile-up (He/H → 8, then NaN at r≈1.24)
  appears when `n_He` overshoots against the lagged `n_H`; a physical cap `fHe ≤ HeH`
  removes the NaN but the aloft depletion remains.
- The earlier apparent "8%→2.5% plateau" was an artifact of a NaN-corrupted run, not a real
  result.

Unresolved: the settling term drains He far more than the coefficients imply. Hand estimates
give a settling velocity `w_s = D·G ≈ 1.2e3 cm/s` vs wind `v ≈ 1.9e5` (w_s/v ≈ 0.006), and
in the assembled flux the advection coefficient (~1.9e5) dominates diffusion (~4.4e3) and
settling (~6e2) by 40–300×, so the discrete steady state *should* be advection-dominated
(nearly well-mixed). The observed total settling implies a ~100× discrepancy not yet located.

Leading remaining suspects: (i) the lagged-`n_H` linearization biases the steady state
when He redistributes strongly (needs an inner Newton/Picard iteration on `n_H` or a
nonlinear solve); (ii) the strong molecular-diffusion term drives toward the *settled*
equilibrium `∂f/∂r = −fG` and the operator-split steady-state selection with the top
advective-outflow BC may not be picking the physical (advection-dominated) branch; (iii) a
full Koskinen-style multicomponent flux (ambipolar + proper `n_H` coupling) may be needed.

Next concrete step: instrument the *converged* flux balance per cell (advective vs diffusive
vs settling fluxes, and the total He flux `r²(n_He v + Φ)` which must be constant at steady
state) to localize where He is being lost, and/or validate against a hand-integrated 1-D
steady advection–diffusion ODE with the same coefficients.

## 7d. Root cause found and fixed — GATE PASSED (2026-07-01)

The flux-balance instrumentation showed the settling flux `F_set` exceeding the advective
flux `F_adv` by up to ~4× where physics demands the opposite (`w_s/v ≈ 0.006`). Printing the
internal coefficients revealed the settling coefficient `G` was ~5800× too large because the
**temperature `TK` had been floored to 1 K**: the routine declared a local time scale named
`t0` (`t0 = R0/v0`), which — Fortran being **case-insensitive** — silently shadowed the
**global temperature normalization `T0`** used one line earlier in `TK = Tcode*T0`. At that
point the local `t0` was still unassigned (≈0), so `TK = Tcode*0 → 0 → floored to 1`, and
`G ∝ 1/TK` blew up. The two constants coincidentally have nearly equal magnitude
(`R0/v0 ≈ 2.83e4 s`, `T0 ≈ 2.83e4 K`), which hid the bug in spot checks.

**Fix:** rename the local time scale `t0 → tscale`. With that one change the He/H profile
becomes physical:

| r [R_p] | 1.05 | 1.1 | 1.4 | 2.0 | 3.0 | 4.0 |
|---|---|---|---|---|---|---|
| (He/H)/HeH | 0.63 | 0.61 | 0.83 | 0.84 | 0.17 | 0.24 |

i.e. mild diffusive separation — He/H ~0.6–0.9× the base value in the inner thermosphere,
falling to ~0.17–0.24× aloft — consistent in **magnitude** with Taylor et al. (2025)
(8%→2.5%, i.e. →0.3×) and Xing et al. (2023). Stable (no NaN), flag OFF byte-identical.

Remaining Phase-1 limitations (documented, non-blocking):
- A near-base numerical **pile-up** from the **lagged-`n_H`** nonlinearity (as He
  concentrates, `n_H` should drop at fixed `rho`, capping `f=n_He/n_H`; lagging `n_H` breaks
  that feedback). Held in check by a physical limiter `fHe ≤ HeH`; the aloft separation (the
  observable) sits below HeH and is unaffected. Proper fix = a self-consistent `n_H`
  (Picard/Newton inner iteration) — Phase 2.
- The inner-region profile is not perfectly monotonic (mild wiggles), and quantitative
  mass-loss needs a Newton-finish (step-capped `Ṁ` is not the converged value).
- Ambipolar field, thermal diffusion, and metal diffusion remain Phase-2 items.

Flag stays **default OFF** (validated, opt-in) pending the Phase-2 refinements above.

## 7e. Phase 2 — implemented and validated (2026-07-02)

Extends the Phase-1 He/H kernel (now factored into `solve_1elem`) with:

- **P2b — ambipolar-corrected settling (default ON, `He_ambipolar`).** In the ionized wind
  the polarization field `eE = -(1/n_e)dp_e/dr` lifts ions; for an H⁺ background it supports
  ~½ a proton weight per unit charge, so a species of charge Z has effective mass
  `m − Z·m_H/2`. The He-vs-H relative settling mass becomes
  `Δm_eff = 3 − 0.5(Z̄_He − Z̄_H)` (mean charges from the local ionization state) → 3 at the
  neutral base, 2.5 in the fully-ionized He⁺⁺/H⁺ wind. Effect on HD 209458b is modest
  (~2% less depletion aloft: (He/H)/HeH at 3 R_p 0.165→0.168), as expected.

- **P2c — thermal diffusion (default `He_alphaT = 0`, no-op).** Adds `α_T ∂lnT/∂r` to the
  settling coefficient. Off by default; activatable via `He_alphaT`.

- **P2d — metal diffusion element by element (default OFF, `He_metal_diffusion`).** Each trace metal
  element is diffused **independently** against the (post-He-rescale) background n_H with its
  own mass `melem_A`, binary D (`1.52e18(1/m_H+1/m_X)^½ T^½/n`), and mean-charge ambipolar
  correction, via the shared `solve_1elem` kernel; the element's ion stages are rescaled to
  the diffused total (cap at the reservoir metal/H). Metals are trace (~1e-3 of the mass) so
  their diffusion does not feed back on n_H (no mass-budget coupling). **Validated
  (HD 209458b, C/N/O):** heavier elements deplete more — at 3 R_p, element/base ≈ He 0.171,
  C 0.158, N 0.136, O 0.105 (monotonic in mass), reproducing the Xing et al. (2023)
  "heavier fractionates more strongly" result. Regression: with `He_metal_diffusion` OFF,
  He/H is unchanged and C/H stays 1.0 (frozen).

- **P2a — self-consistent n_H (NOT adopted).** A Picard iteration coupling n_He↔n_H by mass
  conservation was tried to remove the fHe≤HeH cap. It **diverges**: the near-base settling
  genuinely concentrates He, driving n_H→0 (all mass He), so `f=n_He/n_H` and the 1/n_H flux
  coefficients blow up. The concentration is real, not a lag artifact; the physical cap
  `fHe≤HeH` (He/H cannot exceed the reservoir value) remains the correct, stable limiter.

**Settling discretization (robustness for heavy metals).** The settling drift uses a
**Peclet-based central/upwind hybrid**: central differencing where the drift is well-resolved
(`|D·G|·dr < 2(D+K_zz)`, e.g. light He — less numerical diffusion), first-order upwind where
not (heavy metals like Fe have large `G` that would otherwise break central differencing and
NaN). These were needed to run **WASP-121b** (Fe, mass 56) stably.

**Metal-rescale ratchet bug (found in code review, 2026-07-02, FIXED).** The first version of
the metal-stage rescale had two bugs that together made metal depletion irreversible: (i) an
`rX ≤ 1` clamp — redundant for the metal/H ≤ reservoir cap (already enforced by
`min(nX, fXbase·nHl)`) but forbidding any *replenishment* of a previously depleted cell (a
one-way ratchet); and (ii) a skip of cells where the element was negligible, which made
exhausted cells permanent holes (once ~0, never refilled). During early relaxation (wind
undeveloped) settling transiently depletes metals; the ratchet locked that in. **Fix:** rX is
no longer clamped above (the cap alone bounds it), and an exhausted cell that the solve
replenishes is re-seeded through the neutral stage (the ionization equilibrium re-partitions
next step). Consequence: the previously reported WASP-121b "metal homopause" (Fe → 0 by
~1.25 R_p) was an **artifact and is retracted** — with the fix, WASP metals track He (even Fe
is advection-dominated there, w_s/v ~ 1e-3), and the HD 209458b results return to the
pre-hybrid values (10830: 70.5%→27.1%).

**Two-planet validation (`docs/version_compare.{md,tex,pdf}`).** v1.0 (no diffusion) vs
current on HD 209458b and WASP-121b shows two regimes: on the gentler HD 209458b helium
itself separates (He I 10830 70.5%→27.1%, −2.6×) and heavier C/N/O deplete more than He aloft
(He > C > N > O at 3 R_p); on the furious WASP-121b (log Ṁ ≈ 13.3) advection dominates for
every species — He *and* metals are dragged out essentially unfractionated (10830 unchanged;
metal/H tracks He/H) — the (settling)/(wind) scaling of Koskinen (2013) / Xing (2023). Both
planets stable (no NaN); Ṁ essentially unchanged (metals trace).

**Known minor approximations (from the 2026-07-02 code review; accepted, documented):**
- *Metal mass is not returned to H when metals deplete*: the depleted metal mass simply
  leaves the `Σ m_s f_sp = 1` budget (≤ ~1% where all metals settle out; negligible now that
  metals barely deplete in practice). A correct return would need the `eos_metals` mass-budget
  flag; not worth the complexity at trace abundances.
- *Outer-boundary inflow carries no He/metal*: the top face uses `max(v,0)` (outflow-only
  advection + zero diffusive flux), so transient inflow at the outer boundary brings in
  element-free gas. Safe (no NaN), correct at the steady outflow.
- *One-step temperature lag*: the diffusion step uses the pre-diffusion `T`; the ionization
  and pressure updates that follow use the new composition, so the inconsistency is one
  relaxation step and vanishes at the fixed point.
- `He_metal_diffusion` is silently inert without `He_diffusion` (the metal loop lives inside
  `he_diffusion_step`); pair the flags.

**Remaining (genuinely future):** the near-base pile-up limiter (a more physical base BC or
stronger near-base K_zz would remove it).

*Update.* The second item listed here — that a Newton finish refines the hydro at frozen
composition, because the diffusion operator is outside the Newton path — no longer describes
the code. `Solver: Newton` together with `He_diffusion` now runs an outer co-convergence
(`EXHALE_main.f90`, the `it_diff` loop): JFNK solve, then 500 diffusion relaxation steps of
the He/H field at the converged wind, repeated until the He/H field moves by less than 1e-3
between passes, at most 5 passes. The diffusion operator is still *not* part of the JFNK
residual; the outer loop is what makes the two consistent.

## 8. Original Phase-2 plan (now implemented above; historical)

- Metals diffuse with their own `D_i`, `m_i` (loop the generic operator) — heavier ⇒ stronger
  depletion; needed for metal-line transmission (Mg II/Ca II/Na I in EXHALE_transit.py).
- Explicit **ambipolar** field `eE = -(1/n_e) dp_e/dr` (replaces the Δm≈3 approximation).
- **Thermal diffusion** α_T ≠ 0.
- Ion-stage-resolved diffusion (vs element-level) and Coulomb/resonant D corrections.

---

## Decisions requested before coding

- **D1. Diffusion coefficient:** neutral binary Mason&Marrero `D=b/n`, `b≈1.04e18 T^0.732`,
  ionized correction deferred. (recommend: yes)
  [2026-08-15: not what was implemented. `species_diffusion.f90` uses the
  Banks & Kockarts (1973) binary form
  `D = 1.52e18 (1/m_H + 1/m_He)^(1/2) T^(1/2) / n_tot` cm^2/s, i.e. a
  `T^(1/2)` law, not the Mason & Marrero `T^0.732` fit.]
- **D2. Ambipolar field:** defer to Phase 2 (Δm=3 approximation). (recommend: defer)
- **D3. Thermal diffusion α_T:** 0 in Phase 1. (recommend: 0)
- **D4. Activation:** new `He_diffusion` flag, **default OFF**; enable explicitly, validate,
  then decide on default. (recommend: default OFF)
