# Connecting EXHALE to the lower atmosphere: survey and proposal

**Sources surveyed** (references/): Taylor et al. 2025, 2026 (μbar coupling to Lavvas &
Arfaux); Koskinen et al. 2013a, 2022 (CETIMB lower boundary; analytic lower column; H₃⁺);
Lavvas et al. 2014 / Lavvas & Koskinen 2017-era / Lavvas & Arfaux 2021 (the lower/middle
atmosphere model family); Salz et al. 2016 (TPCI = PLUTO+CLOUDY full-column coupling);
Huang et al. 2017 (hydrostatic atomic layer on an analytic molecular spacer).

**EXHALE today:** base at the input "Planet radius" (≈1 μbar level, user-guessed), fixed
T₀ = T_eq, fixed n₀ (or the pressure-anchored `base_p_ubar` mode), purely atomic H/He +
trace metals at prescribed abundances. No molecules. The Wind-AE port carries a *simplified*
molecular layer (μ-blend `molec_adjust` + bolometric heating/cooling with automatic
shutoff when the atomic transition reaches the base).

---

## 1. What the literature actually does at the interface

| Approach | Lower atmosphere | Hands to the wind model at ~1 μbar | Cost / availability |
|---|---|---|---|
| **Salz 2016 (TPCI)** | none — CLOUDY solves the whole column (molecules disabled) | n/a (single domain, base at n=10¹⁴ cm⁻³, T_eq) | ~3×10⁵ CPU-h for 18 planets; interface not public |
| **Koskinen 2022** | *analytic*: isothermal-T_eq hypsometric column from 1 bar, μ(p) from chemical-equilibrium H₂/H/He fits (Visscher) | base **radius** r₀(1 μbar), T₀=T_eq, q_H₂/q_H/q_He; base v from mass conservation | trivial (closed-form; all equations in the paper) |
| **Huang 2017** | *analytic spacer*: isothermal μ=2.3 column, 1 bar → 10 μbar (their Eq. 15) | base radius of the atomic layer; Lyα **absorbing** bottom boundary (H₂ accidental resonances, τ=1 at N_H₂≈10¹⁴ cm⁻²) | trivial |
| **Taylor 2025/26** | full Lavvas & Arfaux model (RC temperature + ~150-species kinetics + haze microphysics, 10³→10⁻⁶ bar) | T(1 μbar), altitude of 1 μbar, composition (with a *tolerated discontinuity*), K_zz | in-house, **not public** ("reasonable request"); coupling is a manual splice, iterated ~3× |

Key quantitative facts:

- **Base insensitivity for strong winds (Salz §3.5–3.6):** T_base ±40% → Ṁ ±5–10%
  (linear); base density 10¹⁴→10¹² cm⁻³ → Ṁ ≲40% lower. "We do not need to know the
  temperature at the lower boundary with high precision." Adopted error budget: base
  density ±50%, base T ±10% — both dwarfed by the XUV uncertainty (×3).
  ⇒ **EXHALE's fixed-T_eq base is already defensible for hot Jupiters.**
- **Where the fixed base breaks:** (i) planets **near radiative equilibrium** (high
  gravitational potential — stable thermospheres), where everything happens in the lower
  region; (ii) **cool/sub-Neptune planets** where H₂ survives into the wind — Salz's
  GJ 1214 b molecular test: H₂ 10–25% throughout, H⁻ a significant coolant, Ṁ −15%;
  Taylor 2026 GJ 1214 b: adding H₂ **halves the He 10830 signal** (14%→7%); Koskinen
  2022 hot Uranus: base q_H₂ ≈ 0.84 and **H₃⁺ IR emission is the dominant radiative
  coolant below ~2.3 R_p**.
- **What a photochemical lower model uniquely provides** (Lavvas 2014/2021): the
  H₂→H transition can sit **exactly at 1 μbar** and is temperature-sensitive there;
  photochemistry (OH-catalysed H₂ destruction) dissociates H₂ even at 500–1000 K where
  thermal equilibrium would not; **atomic-metal release profiles** (Na, K, Mg, Ca, Fe,
  Al, Si return to atomic form below ~10⁻³ bar) — i.e. the physically-motivated *base
  abundances for EXHALE's trace metals*; photoionization-dominated electron density
  (~10⁸ cm⁻³, 100× Saha) between 10 μbar and 1 bar; radical/haze heating of ±100–400 K
  right at the 1 μbar handoff level.
- **The published coupling is loose, not monolithic:** Lavvas 2014 takes T(p>1 μbar…top)
  *from* the Koskinen thermosphere and uses "species < 3 amu escape at the wind
  velocity" as its upper BC; Taylor tolerates a composition discontinuity at 1 μbar.
  Nobody runs a tight per-timestep two-way coupling except TPCI (and it disables
  molecules).

---

## 2. Proposal: three tiers (phase-by-phase, each with a validation gate)

### Tier 1 — analytic lower column (Koskinen-2022 style). *Recommended first; ~days.*

Implement the closed-form lower/middle atmosphere of Koskinen 2022 §2.2 as a small
module (Fortran `lower_column.f90` or a Python pre-processor):

- Integrate the hypsometric equation from p = 1 bar (at the observed transit radius) to
  the 1 μbar base, isothermal at T_eq, with μ(p) from the Visscher chemical-equilibrium
  H₂/H/He partition (their Eqs. 11–13, fully printed in the paper; validated against
  NASA CEA at 1000–2500 K). Roche-potential gravity along the column (EXHALE's
  `grav_field` already provides it).
- Deliverables to EXHALE: **base radius r₀(1 μbar)** — replacing the user-guessed
  "Planet radius" with a derived quantity (or a consistency check on it), **base
  H₂/H/He fractions**, and the correct base mean molecular weight.
- Same machinery gives the TPM/Balmer "molecular spacer" of Huang 2017 (their Eq. 15 is
  the μ=const degenerate case) — useful for the transmission-continuum radius.
- **Validation gate:** reproduce Koskinen 2022 Model A (hot Uranus, 0.05 au):
  r₀ = 1.34 R_p, T₀ = 1140 K, q_H₂ ≈ 0.84 at 1 μbar. And for a hot Jupiter confirm
  q_H₂(1 μbar) ≈ 0 (thermal dissociation), i.e. EXHALE's atomic assumption is verified,
  not assumed.
- Caveat to document: chemical equilibrium *underestimates* H dissociation at
  T_eff ≈ 1000–2000 K (photochemistry ignored) — the Tier-3 hook.

### Tier 2 — molecular extension of EXHALE (H₂, H₂⁺, H₃⁺, HeH⁺ + H₃⁺ cooling). *~weeks; required only for warm Neptunes / sub-Neptunes.*

The complete reaction network **with all rate coefficients is printed in Koskinen 2022
Table 1** (P1–P5, R1–R23): H₂ photoionization + dissociative photoionization, thermal
(R12) and electron-impact (R14) dissociation, H₂⁺+H₂→H₃⁺+H (R8), H₃⁺ dissociative
recombination (R6/R7, Larsson 2008), HeH⁺ chain (R16–R20), three-body H₂ formation
(R15). Work items:

- extend the coupled ionization system by 4 species (a `System_HeH_mol*` sibling; the
  merged-system pattern from `System_HeH_TR_metals` applies directly);
- EOS: μ and n_tot with H₂ (already centralized in `composition.f90` / `calc_ntot` —
  one policy point), γ/dof for a diatomic, and the 4.48 eV dissociation energy sink in
  the energy equation;
- **H₃⁺ infrared cooling**, optically thin: start with the Miller et al. (2013) LTE
  per-molecule emission fits (analytic log-polynomial in T); defer the non-LTE
  correction factor (Koskinen et al. 2009) with a documented caveat — same
  fit-first/refine-later pattern as our CHIANTI cooling work;
- optional interim step (**Tier 2a, ~days**): promote the Wind-AE port's *passive*
  molecular-base treatment (μ-blend `molec_adjust` + bolometric heating/cooling +
  automatic shutoff) into EXHALE's base cells — no chemistry, but the correct μ and
  energy budget where a thin H₂ layer overlaps the domain base.
- **Validation gates:** (i) hot Jupiter with molecules on ⇒ H₂ vanishes at the base,
  results reproduce the atomic code (regression); (ii) GJ 1214 b-like case ⇒ He 10830
  drops by ~×2 with H₂ on (Taylor 2026), Ṁ shifts ≲40% (Koskinen 2022 photodissociation
  sensitivity was ×1.4; Salz −15%).

### Tier 3 — coupling to a real photochemical/RC lower-atmosphere model. *~months; the metal-abundance payoff.*

The Lavvas & Arfaux code is not public, so the practical paths are:

- **(a) Open-source stack (recommended):** `VULCAN` (public 1-D photochemical kinetics,
  C/H/N/O/S networks, K_zz, photolysis; the closest public analogue of the Lavvas 2014
  kinetics) for composition, plus a temperature model — either `HELIOS` (public RC) or,
  cheaper, an analytic Guillot/picket-fence T(p) as the first iteration. A thin Python
  driver (`lower_atmosphere/run_lower.py`) runs the stack for the planet and writes a
  new optional EXHALE input file, e.g. `base.inp`:

  ```
  # base.inp — written by the lower-atmosphere driver (VULCAN + T(p))
  T_base   [K]      1140.0
  r_base   [RJ]     1.402      # radius of the 1 ubar level
  q_H2              0.03
  q_H               0.88
  q_He              0.09
  # atomic-metal release fractions at 1 ubar (replace solar totals in metals.inp)
  Na  1.6e-6
  Mg  3.1e-5
  ...
  Kzz_base [cm2/s]  1.0e9
  ```

  EXHALE-side work is small: extend `input_read.f90` to consume `base.inp`
  (override T₀, base radius, He/H, metal abundances; `He_Kzz` default), exactly the
  keyword-block pattern already used for `metals.inp`/`opacity.inp`.
- **(b) Collaboration path:** request Lavvas & Arfaux outputs for target planets
  ("shared on reasonable request") — the fastest route to a Taylor-grade coupled model
  for a specific system, at the cost of dependency.
- **(c) TPCI-style full coupling: rejected.** ~3×10⁵ CPU-h class, interface not public,
  and Salz's own sensitivity results show it is unnecessary for our strong-wind targets.

**Coupling protocol (one-way first, then iterate):**
1. lower model → `base.inp` → EXHALE run (Dirichlet base; composition discontinuity
   tolerated, as in Taylor);
2. optionally feed back: EXHALE's transmitted XUV spectrum at the base (we already
   compute the attenuated flux — add a dump of F_ν(r_base)) as the photochemical
   model's top irradiation, and Ṁ (or per-species escape velocity) as its upper BC —
   precisely the Lavvas 2014 upper BC ("species < 3 amu escape at the model wind
   velocity"). Taylor found ~3 iterations suffice for the analogous Lyα loop.

**Validation gate:** HD 209458 b with the full stack vs. published Taylor 2025 base
values (T ≈ 1270 K at 10⁻⁶ bar, their Model C/best-fit family), and metal base
abundances vs. Lavvas 2014 Fig. 9 (Mg/Fe/Si ionic above 10⁻⁶ bar, Na/K atomic).

### What each science case actually needs

| Science case | Needed tier |
|---|---|
| Hot Jupiters (HD 209458 b, HD 189733 b, WASP-121 b) — current program | Tier 1 only (base radius + verified atomic base); fixed-T_eq base already within ±10% on Ṁ (Salz) |
| He 10830 / Hα population work | Tier 1 (+ the Huang absorbing-bottom check on `lya_rt`'s lower BC — audit item) |
| Metal transmission lines (Phase 5, WASP-121 b) | Tier 3(a): photochemical atomic-metal release at the base instead of assumed solar totals (condensation/molecule sequestration) |
| Warm Neptunes / sub-Neptunes (GJ 1214 b-class, future) | Tier 2 mandatory (H₂ in the wind, H₃⁺ cooling, ×2 He-signal effect) + Tier 1 |

### Effort summary

| Tier | New code | Effort | Risk |
|---|---|---|---|
| 1 analytic column | `lower_column.f90` (~150 SLOC) or Python pre-processor | days | low (closed-form, gate = published numbers) |
| 2a passive molecular base | reuse `wind_ae` molec machinery | days | low |
| 2 full molecular chemistry | `System_HeH_mol*`, EOS, H₃⁺ cooling | weeks | medium (solver stiffness; γ/EOS consistency) |
| 3(a) open-source lower stack | Python driver + `base.inp` reader | weeks–months (mostly VULCAN/HELIOS learning + per-planet setup) | medium (network/opacity choices); EXHALE-side change is small |

**Recommended order:** Tier 1 → (science-driven fork) hot-Jupiter metals ⇒ Tier 3(a);
sub-Neptunes ⇒ Tier 2 (with 2a as the quick bridge). Tier 1 is worth doing
unconditionally: it removes the last *ad hoc* number in EXHALE's setup (the base
radius) at negligible cost, and its H₂-fraction output is the automatic switch that
tells us when Tier 2 physics is actually required for a given planet.


---

## Implementation status (2026-07-02)

All Fortran lives in the new folder `src/modules/lower_atmosphere/`; the Python driver in
`src/utils/`. Everything is **opt-in / default off** — with no keys and no `base.inp`,
standard runs are unchanged (regression: HD 209458 b He2³S+metals 3000-step Ṁ = 8.97,
identical to EXHALE_v1.0).

| Tier | Deliverable | Status | Validation |
|---|---|---|---|
| 1 | `lower_column.f90` — hypsometric column, Visscher chem-eq H₂/H/He (+ fully-atomic bracket); key `Lower column: <R_1bar RJ>` reports r₀(1 μbar), base q_H₂/q_H/q_He/μ vs the input "Planet radius" | **done** | Koskinen 2022 Model A gate: r₀/R₁ᵦₐᵣ = 1.343 (paper 1.34), q_H₂ = 0.838 (0.84), q_H = 0.027 (0.026) |
| 2 (foundation) | `h3p_cooling.f90` — Miller+2013 Table-5 LTE emission fits + Table-6 non-LTE factor s(T,n_H₂), bilinear; `mol_rates.f90` — Koskinen 2022 Table-1 rates R1–R23 (verified against the PDF, saved as `references/Miller_2013_JPCA_117_9770.pdf` / `Koskinen_2022...`) | **done (standalone)** | fit reproduces Miller Table-4 anchors at 500–5000 K to <0.5%; s(1000 K, 10¹⁰ cm⁻³)=0.4955 exact |
| 2 (core) | **`System_HeH_mol.f90`** — coupled H⁺/He⁺/He⁺⁺(+2³S) + H₂/H₂⁺/H₃⁺/HeH⁺ equilibrium (7–8 unknowns, hybrd1; atomic rows use EXHALE's own rates so the molecule-free limit reproduces the atomic systems); **σ_H₂** (Yan+1998 Eqs. 17–19, `cross_sec.f90`) wired into opacity/photoionization/heating (`PH_heat_HHe`); H₃⁺ cooling in the `cool` array; EOS (`calc_ne/ntot/rho`, `composition`) molecule-aware; output/IC columns H2/H2p/H3p/HeHp; key `Molecular chemistry: True` (v1: metals excluded) | **done (core)** | σ_H₂ reproduces Yan Table 7 (0.04761/0.006169/0.001739 Mb at 100/200/300 eV); Gate 0: mol-off byte-equivalent (HD209 3000-step Ṁ=8.97); Gate 1 (HD209 mol-on): fully molecular base + sharp H₂→H front at r=1.019 R_p, wind above front ≈ atomic (Ṁ 8.91 vs 8.97, He2³S peak +1.4%); Gate 2 (hot-Uranus-like): front rises to r=1.148, H₃⁺ active in the molecular layer |
| 2 (remaining) | local-equilibrium caveat (no molecular *advection* — Koskinen's high-altitude H₂ replenishment not reproduced); Lyman–Werner photodissociation; dissociative/double photoionization channels (P4/P5); 4.48 eV dissociation energy sink; diatomic γ; metals+molecules merge; `_adv` post-process; GJ 1214 b He-halving validation (full M-dwarf setup) | open | gates defined in §Tier-2 |
| 2a | `Molecular base: True` — EOS-only base correction: removes the H₂-bound particles from `ntot_bc` via the equilibrium fit (lower base pressure / heavier base μ; chemistry stays atomic — crude, documented) | **done** | HD 209458 b: q_H₂(1 μbar,1450 K)=0.831 → ntot_bc 1.0→0.546 |
| 3 | `base.inp` reader in `input_read` (T_base / r_base / HeH_base / Kzz_base, echo + no-op when absent) + `src/utils/run_lower.py` driver (isothermal or Guillot 2010 semi-grey T(p); writes base.inp with a molecular-base warning) | **done (analytic stack)** | end-to-end: driver → base.inp → EXHALE consumes and echoes; iso vs Guillot: r₀ 1.4723 vs 1.4660 R_J, T_base 1450 vs 1313 K (HD 209458 b) |
| 3 (upgrade) | **VULCAN end-to-end**: public VULCAN cloned (`../VULCAN`, FastChem compiled), HD 189733 b SNCHO photochemical run converged (2356 steps); converter `src/utils/vulcan_to_base.py` (.vul → base.inp with photochemical q_H2/q_H, VULCAN-μ/T hypsometric r_base; molecular mixing ratios as comments; **no metal release** — outside VULCAN's scope, stays Lavvas-only) | **done (H/C/N/O composition)** | HD 189733 b at 1 μbar: **q_H2=0.63, q_H=0.23 — photochemistry dissociates ~11× more H than the equilibrium column (q_H=0.020)**, directly quantifying the Tier-1 caveat; r_base 1.168 vs analytic 1.174 R_J; T_base 863 K (Moses11 T(p)); EXHALE consumes the file (overrides echoed) |

Notable physics finding from the Tier-1 gate work: the Visscher **equilibrium** fit keeps the
1 μbar base strongly molecular up to T ≈ 2000 K (fully atomic only above ~2400 K), so for
Teq ~ 1000–2000 K hot Jupiters the *equilibrium* base contradicts the standard atomic
assumption — the dissociation there is photochemical (Moses 2011; Koskinen 2013a), which is
exactly what the Tier-3 VULCAN upgrade supplies. The Tier-1 report therefore always prints
the atomic-bracket radius alongside the equilibrium one.

---

## Example applications (2026-07-02): HD 209458 b, HD 189733 b, WASP-121 b, WASP-52 b

Full tables and interpretation: the tex/pdf version (§Examples); run configs:
`examples/13_lower_atmosphere/`. Isothermal-Teq column (r in R_J):

| planet | input R0 | r0(1 μbar) chem.eq | atomic bracket | q_H2(base) | verdict |
|---|---|---|---|---|---|
| HD 209458 b | 1.401 | 1.472 | 1.582 | 0.83 | input ~5% below bracket |
| HD 189733 b | 1.193 | 1.174 | 1.205 | 0.84 | input INSIDE bracket ✓ |
| WASP-121 b | 2.2075 | 2.010 | 2.133 | 0.005 | **atomic base VERIFIED** (UHJ; even in equilibrium) |
| WASP-52 b | 1.270 | 1.437 | 1.615 | 0.84 | input = transit radius → base 0.17–0.35 R_J too deep; **recomputed** (below) |

Guillot T(p) variant: r0 = 1.466 / 1.171 / 1.966 / 1.428 R_J, T_base = 1313 / 1072 /
2136 / 1181 K respectively. Common pattern: every Teq ≲ 1500 K planet has a strongly
molecular *equilibrium* base (q_H2 ≈ 0.83) — the atomic base rests on photochemical
dissociation, the quantified motivation for the VULCAN tier. WASP-121 b's production base
radius (2.2075, from the Huang 2023 log-g anchoring) sits ~4–10% above the
isothermal-Teq bracket, consistent with a dayside-hot lower atmosphere.

**WASP-52 b recomputation (2026-07-02).** The `fxuv1p0_he98` configuration was rerun
(Newton-converged, current binary) with the corrected equilibrium-column base
r0 = 1.437 R_J (vs the old 1.27 = transit-radius shortcut), all else identical:

| | old (1.270) | new (1.437) | change |
|---|---|---|---|
| log₁₀ Ṁ [g/s] | 11.97 | 12.14 | ×1.5 |
| He 10830 line-center absorption | 27.11% | 27.56% | ×1.02 |
| He 10830 equivalent width | 0.713 Å | 0.722 Å | ×1.01 |
| peak n(He 2³S) [cm⁻³] | 150.5 | 147.6 | ~unchanged |

The line-forming effective radius is set by the τ=1 surface of the extended wind, which
stays nearly fixed in *physical* units, so the He 10830 observable is essentially
insensitive to the base misplacement — the bias is in **Ṁ (×1.5)**, not in the line. The
He 10830 scan conclusions therefore appear robust to the base-radius shortcut.