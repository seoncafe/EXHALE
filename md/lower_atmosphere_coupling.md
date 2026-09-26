# Connecting EXHALE to the lower atmosphere: survey and proposal

> **Figures:** the tex/pdf version now carries four result figures (Tier-1 4-planet
> brackets, Tier-2 dissociation fronts/molecular species, VULCAN-vs-equilibrium,
> WASP-52b He 10830 base-radius comparison): sources and data snapshots in
> `docs/lower_atmosphere_figs/` (`make_figures.py`).

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

> [2026-08-15: this paragraph describes the tree as of 2026-07 and is no longer
> current. The proposal below has been built: the coupled molecular network
> (H₂/H₂⁺/H₃⁺/HeH⁺, `Molecular chemistry: True`), H₃⁺ cooling, Lyman-Werner
> photodissociation, the analytic lower column and the `base.inp` handoff all
> exist in `src/modules/lower_atmosphere/`. All of them default off, so an
> ordinary run is still the atomic base described here. The per-tier status
> table further down is the authoritative account of what is done.]

---

## 1. What the literature actually does at the interface

| Approach | Lower atmosphere | Hands to the wind model at ~1 μbar | Cost / availability |
|---|---|---|---|
| **Salz 2016 (TPCI)** | none: CLOUDY solves the whole column (molecules disabled) | n/a (single domain, base at n=10¹⁴ cm⁻³, T_eq) | ~3×10⁵ CPU-h for 18 planets; interface not public |
| **Koskinen 2022** | *analytic*: isothermal-T_eq hypsometric column from 1 bar, μ(p) from chemical-equilibrium H₂/H/He fits (Visscher) | base **radius** r₀(1 μbar), T₀=T_eq, q_H₂/q_H/q_He; base v from mass conservation | trivial (closed-form; all equations in the paper) |
| **Huang 2017** | *analytic spacer*: isothermal μ=2.3 column, 1 bar → 10 μbar (their Eq. 15) | base radius of the atomic layer; Lyα **absorbing** bottom boundary (H₂ accidental resonances, τ=1 at N_H₂≈10¹⁴ cm⁻²) | trivial |
| **Taylor 2025/26** | full Lavvas & Arfaux model (RC temperature + ~150-species kinetics + haze microphysics, 10³→10⁻⁶ bar) | T(1 μbar), altitude of 1 μbar, composition (with a *tolerated discontinuity*), K_zz | in-house, **not public** ("reasonable request"); coupling is a manual splice, iterated ~3× |

Key quantitative facts:

- **Base insensitivity for strong winds (Salz §3.5-3.6):** T_base ±40% → Ṁ ±5-10%
  (linear); base density 10¹⁴→10¹² cm⁻³ → Ṁ ≲40% lower. "We do not need to know the
  temperature at the lower boundary with high precision." Adopted error budget: base
  density ±50%, base T ±10%, both dwarfed by the XUV uncertainty (×3).
  ⇒ **EXHALE's fixed-T_eq base is already defensible for hot Jupiters.**
- **Where the fixed base breaks:** (i) planets **near radiative equilibrium** (high
  gravitational potential, stable thermospheres), where everything happens in the lower
  region; (ii) **cool/sub-Neptune planets** where H₂ survives into the wind, Salz's
  GJ 1214 b molecular test: H₂ 10-25% throughout, H⁻ a significant coolant, Ṁ −15%;
  Taylor 2026 GJ 1214 b: adding H₂ **halves the He 10830 signal** (14%→7%); Koskinen
  2022 hot Uranus: base q_H₂ ≈ 0.84 and **H₃⁺ IR emission is the dominant radiative
  coolant below ~2.3 R_p**.
- **What a photochemical lower model uniquely provides** (Lavvas 2014/2021): the
  H₂→H transition can sit **exactly at 1 μbar** and is temperature-sensitive there;
  photochemistry (OH-catalyzed H₂ destruction) dissociates H₂ even at 500-1000 K where
  thermal equilibrium would not; **atomic-metal release profiles** (Na, K, Mg, Ca, Fe,
  Al, Si return to atomic form below ~10⁻³ bar), i.e. the physically-motivated *base
  abundances for EXHALE's trace metals*; photoionization-dominated electron density
  (~10⁸ cm⁻³, 100× Saha) between 10 μbar and 1 bar; radical/haze heating of ±100-400 K
  right at the 1 μbar handoff level. *(2026-08-19: EXHALE still has no
  O/OH/H₂O species of its own, so the OH-catalyzed destruction reaches the code
  only through the handoff.)*
- **The published coupling is loose, not monolithic:** Lavvas 2014 takes T(p>1 μbar…top)
  *from* the Koskinen thermosphere and uses "species < 3 amu escape at the wind
  velocity" as its upper BC; Taylor tolerates a composition discontinuity at 1 μbar.
  Nobody runs a tight two-way coupling at every timestep except TPCI (and it disables
  molecules).

---

## 2. Proposal: three tiers (phase-by-phase, each with a validation gate)

### Tier 1: analytic lower column (Koskinen-2022 style). *Recommended first; ~days.*

Implement the closed-form lower/middle atmosphere of Koskinen 2022 §2.2 as a small
module (Fortran `lower_column.f90` or a Python pre-processor):

- Integrate the hypsometric equation from p = 1 bar (at the observed transit radius) to
  the 1 μbar base, isothermal at T_eq, with μ(p) from the Visscher chemical-equilibrium
  H₂/H/He partition (their Eqs. 11-13, fully printed in the paper; validated against
  NASA CEA at 1000-2500 K). Roche-potential gravity along the column (EXHALE's
  `grav_field` already provides it).
- Deliverables to EXHALE: **base radius r₀(1 μbar)**, replacing the user-guessed
  "Planet radius" with a derived quantity (or a consistency check on it), **base
  H₂/H/He fractions**, and the correct base mean molecular weight.
- The same handling gives the TPM/Balmer "molecular spacer" of Huang 2017 (their Eq. 15 is
  the μ=const degenerate case): useful for the transmission-continuum radius.
- **Validation gate:** reproduce Koskinen 2022 Model A (hot Uranus, 0.05 au):
  r₀ = 1.34 R_p, T₀ = 1140 K, q_H₂ ≈ 0.84 at 1 μbar. And for a hot Jupiter confirm
  q_H₂(1 μbar) ≈ 0 (thermal dissociation), i.e. EXHALE's atomic assumption is verified,
  not assumed.
- Caveat to document: chemical equilibrium *underestimates* H dissociation at
  T_eff ≈ 1000-2000 K (photochemistry ignored): the Tier-3 hook.

### Tier 2: molecular extension of EXHALE (H₂, H₂⁺, H₃⁺, HeH⁺ + H₃⁺ cooling). *~weeks; required only for warm Neptunes / sub-Neptunes.*

The complete reaction network **with all rate coefficients is printed in Koskinen 2022
Table 1** (P1-P5, R1-R23): H₂ photoionization + dissociative photoionization, thermal
(R12) and electron-impact (R14) dissociation, H₂⁺+H₂→H₃⁺+H (R8), H₃⁺ dissociative
recombination (R6/R7, Larsson 2008), HeH⁺ chain (R16-R20), three-body H₂ formation
(R15). Work items:

- extend the coupled ionization system by 4 species (a `System_HeH_mol*` sibling; the
  merged-system pattern from `System_HeH_TR_metals` applies directly): done, and the
  same pattern then carried the metals in (`System_HeH_mol_metals`);
- EOS: μ and n_tot with H₂ (already centralized in `composition.f90` / `calc_ntot`,
  one policy point), γ/dof for a diatomic, and the 4.48 eV dissociation energy sink in
  the energy equation;
- **H₃⁺ infrared cooling**, optically thin: start with the Miller et al. (2013) LTE
  emission fits for each molecule (analytic log-polynomial in T); defer the non-LTE
  correction factor (Koskinen et al. 2009) with a documented caveat: same
  fit-first/refine-later pattern as our CHIANTI cooling work;
  [2026-08-15: the non-LTE factor was not deferred. `h3p_nonlte_factor(T, nH2)`
  in `src/modules/lower_atmosphere/h3p_cooling.f90` bilinearly interpolates the
  Miller et al. (2013) Table 6 departure factor, which the LTE emission is
  multiplied by.]
- optional interim step (**Tier 2a, ~days**): promote the Wind-AE port's *passive*
  molecular-base treatment (μ-blend `molec_adjust` + bolometric heating/cooling +
  automatic shutoff) into EXHALE's base cells, no chemistry, but the correct μ and
  energy budget where a thin H₂ layer overlaps the domain base.
- **Validation gates:** (i) hot Jupiter with molecules on ⇒ H₂ vanishes at the base,
  results reproduce the atomic code (regression); (ii) GJ 1214 b-like case ⇒ He 10830
  drops by ~×2 with H₂ on (Taylor 2026), Ṁ shifts ≲40% (Koskinen 2022 photodissociation
  sensitivity was ×1.4; Salz −15%).

When `Include He23S?` is also on, the metastable couples to H₂ through Penning
ionization He(2³S) + H₂ → He(1¹S) + H₂⁺ + e⁻: the dominant He(2³S) loss toward
an H₂-dominated base (García Muñoz 2025, A&A 698, A199, Fig. 4). Its rate
coefficient `penning_HeI23S_H2` (Cool_coeff.f90) is an analytic fit to García
Muñoz Table A.5 (Cohen & Lane 1977 cross sections), k(T) = 5.3791e-12
T^0.6760 exp(−695.21/T) cm³ s⁻¹, reproducing the tabulated 500-10000 K points to
≤0.13%. It enters `System_HeH_mol` as an H₂ loss / H₂⁺ source / He(2³S) sink and
adds (E[2³S] − IP[H₂]) ≈ 4.4 eV of electron heating; the minor associative
H + HeH⁺ branch (~10%) is folded into the Penning channel.

**The `_adv` post-process is atomic-only (documented design), so the Penning-H₂
term does not appear there.** The advection-corrected reconstruction in
`post_process_adv.f90` (header, lines 4-13) treats the gas as H/He + trace
metals and excludes the molecular species entirely: the molecular densities are
never passed in or re-solved, and its `calc_ne`/`calc_ntot` calls omit the
`nmol` argument. Consequences for a molecular run: (i) the `_adv` n_tot/n_e
omit the neutral-H₂ particle count and the molecular-ion electrons; (ii) the
He(2³S)+H₂ Penning sink cannot act in the `_adv` solve (no n_H₂ there), so the
`_adv` outputs (and the transit module that reads them) carry the
*equilibrium* molecular-base triplet suppression (39× at the HD 209458 b base,
430× on the hot-Uranus case) without an advection correction. This is
acceptable where `_adv` is physically meaningful (the atomic/ionized escape
flow above the H₂→H front); below the front the atomic `_adv` approximation is
outside its validity domain anyway. A molecular-aware `_adv` is a separate
extension, naturally paired with the open "molecular advection" item in the
Tier-2 remaining list below.

### Tier 3: coupling to a real photochemical/RC lower-atmosphere model. *~months; the metal-abundance payoff.*

The Lavvas & Arfaux code is not public, so the practical paths are:

- **(a) Open-source stack (recommended):** `VULCAN` (public 1-D photochemical kinetics,
  C/H/N/O/S networks, K_zz, photolysis; the closest public analogue of the Lavvas 2014
  kinetics) for composition, plus a temperature model, either `HELIOS` (public RC) or,
  cheaper, an analytic Guillot/picket-fence T(p) as the first iteration. A thin Python
  driver (`src/utils/run_lower.py`) runs the stack for the planet and writes a
  new optional EXHALE input file, e.g. `base.inp`:

  ```
  # base.inp: written by the lower-atmosphere driver (VULCAN + T(p))
  T_base    1140.0            # K
  r_base    1.402             # R_J, radius of the 1 ubar level
  HeH_base  0.0793
  Kzz_base  1.0e9             # cm2/s
  q_H2_base 0.03              # photochemical H2 mixing ratio at the base
  p_base    1.0e-6            # bar, the level all of the above refer to
  # q_H = 0.88, q_He = 0.09 (comments; implied by q_H2_base and HeH_base)
  # atomic-metal release fractions at 1 ubar (replace solar totals in metals.inp)
  Na  1.6e-6
  Mg  3.1e-5
  ...
  ```

  EXHALE-side work is small: extend `input_read.f90` to consume `base.inp`
  (override T₀, base radius, He/H, metal abundances; `He_Kzz` default), exactly the
  keyword-block pattern already used for `metals.inp`/`opacity.inp`. As of 2026-08-10
  `read_base_inp` consumes T_base / r_base / HeH_base / Kzz_base / q_H2_base / p_base;
  the metal release fractions are still not read (no public photochemical source
  produces them: Tier 3(a) remains the Lavvas collaboration path).
- **(b) Collaboration path:** request Lavvas & Arfaux outputs for target planets
  ("shared on reasonable request"): the fastest route to a Taylor-grade coupled model
  for a specific system, at the cost of dependency.
- **(c) TPCI-style full coupling: rejected.** ~3×10⁵ CPU-h class, interface not public,
  and Salz's own sensitivity results show it is unnecessary for our strong-wind targets.

**Coupling protocol (one-way first, then iterate):**
1. lower model → `base.inp` → EXHALE run (Dirichlet base; composition discontinuity
   tolerated, as in Taylor);
2. optionally feed back: EXHALE's transmitted XUV spectrum at the base (we already
   compute the attenuated flux, add a dump of F_ν(r_base)) as the photochemical
   model's top irradiation, and Ṁ (or the escape velocity for each species) as its upper BC,
   precisely the Lavvas 2014 upper BC ("species < 3 amu escape at the model wind
   velocity"). Taylor found ~3 iterations suffice for the analogous Lyα loop.

   [2026-08-27: step 2 is half implemented, in the *elemental-flux* form rather
   than the Ṁ form. `src/utils/element_flux_closure.py` drives the pair of codes
   to the fixed point Phi_El = F_El(Phi_El): the trial flux is imposed on the
   chemistry's upper boundary (`--trial-flux-H` / `--trial-flux-He` of
   `photochem_to_lower_profile.py`, g/s outward positive), the wind's own
   elemental flux is measured face by face by `write_element_flux_profile`
   (`binary_element_diffusion.f90`) and reduced over a radial window, and the
   update is under-relaxed with omega = 0.5 (halved, floor 0.125, on any
   iteration whose residual failed to fall). The convergence test reads the
   undamped residual, and a residual smaller than the spread of its own window
   is reported UNRESOLVED rather than converged. The radiation half (dumping
   F_ν(r_base) and handing it to the photochemical model) is still not
   implemented: the adapters take a stellar flux file and a dilution instead.]

**Validation gate:** HD 209458 b with the full stack vs. published Taylor 2025 base
values (T ≈ 1270 K at 10⁻⁶ bar, their Model C/best-fit family), and metal base
abundances vs. Lavvas 2014 Fig. 9 (Mg/Fe/Si ionic above 10⁻⁶ bar, Na/K atomic).

### What each science case actually needs

| Science case | Needed tier |
|---|---|
| Hot Jupiters (HD 209458 b, HD 189733 b, WASP-121 b): current program | Tier 1 only (base radius + verified atomic base); fixed-T_eq base already within ±10% on Ṁ (Salz) |
| He 10830 / Hα population work | Tier 1 (+ the Huang absorbing-bottom check on `lya_rt`'s lower BC, audit item) [2026-08-15: implemented as the key `Lya absorbing bottom: True`, default off; `src/modules/radiation/lya_rt.f90`] |
| Metal transmission lines (Phase 5, WASP-121 b) | Tier 3(a): photochemical atomic-metal release at the base instead of assumed solar totals (condensation/molecule sequestration) |
| Warm Neptunes / sub-Neptunes (GJ 1214 b-class, future) | Tier 2 mandatory (H₂ in the wind, H₃⁺ cooling, ×2 He-signal effect) + Tier 1 |

### Effort summary

| Tier | New code | Effort | Risk |
|---|---|---|---|
| 1 analytic column | `lower_column.f90` (~150 SLOC) or Python pre-processor | days | low (closed-form, gate = published numbers) |
| 2a passive molecular base | reuse the `wind_ae` molecular routines | days | low |
| 2 full molecular chemistry | `System_HeH_mol*`, EOS, H₃⁺ cooling | weeks | medium (solver stiffness; γ/EOS consistency) |
| 3(a) open-source lower stack | Python driver + `base.inp` reader | weeks-months (mostly VULCAN/HELIOS learning + setup for each planet) | medium (network/opacity choices); EXHALE-side change is small |

**Recommended order:** Tier 1 → (science-driven fork) hot-Jupiter metals ⇒ Tier 3(a);
sub-Neptunes ⇒ Tier 2 (with 2a as the quick bridge). Tier 1 is worth doing
unconditionally: it removes the last *ad hoc* number in EXHALE's setup (the base
radius) at negligible cost, and its H₂-fraction output is the automatic switch that
tells us when Tier 2 physics is actually required for a given planet.


---

## Implementation status (2026-07-02)

All Fortran lives in the new folder `src/modules/lower_atmosphere/`; the Python driver in
`src/utils/`. Everything is **opt-in / default off**, with no keys and no `base.inp`,
standard runs are unchanged by the lower-atmosphere code (the molecular chemistry
off-switch reproduces the atomic systems byte-for-byte; the absolute production
baseline itself moved with the 2026-07-23 defaults: staged secondary ionization,
He_rec_coupling, He-H charge exchange).

| Tier | Deliverable | Status | Validation |
|---|---|---|---|
| 1 | `lower_column.f90`: hypsometric column, Visscher chem-eq H₂/H/He (+ fully-atomic bracket); key `Lower column: <R_1bar RJ>` reports r₀(1 μbar), base q_H₂/q_H/q_He/μ vs the input "Planet radius" | **done** | Koskinen 2022 Model A gate: r₀/R₁ᵦₐᵣ = 1.343 (paper 1.34), q_H₂ = 0.838 (0.84), q_H = 0.027 (0.026) |
| 2 (foundation) | `h3p_cooling.f90`: Miller+2013 Table-5 LTE emission fits + Table-6 non-LTE factor s(T,n_H₂), bilinear; `mol_rates.f90`, Koskinen 2022 Table-1 rates R1-R23 (verified against the PDF, saved as `references/Miller_2013_JPCA_117_9770.pdf` / `Koskinen_2022...`) | **done (standalone)** | fit reproduces Miller Table-4 anchors at 500-5000 K to <0.5%; s(1000 K, 10¹⁰ cm⁻³)=0.4955 exact |
| 2 (core) | **`System_HeH_mol.f90`**: coupled H⁺/He⁺/He⁺⁺(+2³S) + H₂/H₂⁺/H₃⁺/HeH⁺ equilibrium (7-8 unknowns, hybrd1; atomic rows use EXHALE's own rates so the molecule-free limit reproduces the atomic systems); **σ_H₂** (Yan+1998 Eqs. 17-19, `cross_sec.f90`) wired into opacity/photoionization/heating (`PH_heat_HHe`); H₃⁺ cooling **inside `eval_cool`**, with its own `Cooling_breakdown` column, so the marching temperature update and the steady residual balance the same cooling function (2026-08-13); EOS (`calc_ne/ntot/rho`, `composition`) molecule-aware, and the cooling (`eval_cool`) reads the same `calc_ne` electron density as the equilibrium solver, molecular ions included (2026-08-13); output/IC columns H2/H2p/H3p/HeHp; key `Molecular chemistry: True`; trace metals may be solved in the same system (`System_HeH_mol_metals`, 2026-08-13) | **done (core)** | σ_H₂ reproduces Yan Table 7 (0.04761/0.006169/0.001739 Mb at 100/200/300 eV); Gate 0: mol-off byte-equivalent (molecular chemistry off reproduces the atomic systems); Gate 1 (HD209 mol-on): fully molecular base (x_H₂≈0.996) + sharp H₂→H front at r=1.020 R_p, wind above the front ≈ atomic; below the front the molecular base suppresses He2³S by orders of magnitude (He2³S+H₂ Penning destruction, ~10⁴× at the front, rising back to the atomic value above ~1.3 R_p); Gate 2 (hot-Uranus-like: 0.0457 M_J, R_p=0.49 R_J, T_eq=1140 K, HD209 orbit/spectrum): front at r=1.156, H₃⁺ active in the molecular layer (peak ~7×10⁴ cm⁻³ at r≈1.04). The H₂→H fronts and molecular base composition are essentially unchanged from the 2026-07-16 gates (front 1.019→1.020, 1.166→1.156), confirming the dissociation-front result is robust to the 2026-07-23 production defaults (staged secondary ionization, He_rec_coupling, He-H charge exchange). Gate numbers refreshed 2026-07-23 at a shared 12000-step relaxation-snapshot convention: both HD209 runs read Ṁ = log₁₀ 10.58 (relaxation snapshots, not flux-flat converged, under these defaults the HD209 atomic gate plateaus near 4% mass-flux spread). Gate inputs pinned at `lower_atmosphere_figs/data_g*/input.inp` (Update_EXHALE_stage1 §35). 2026-08-13 (Update_EXHALE_stage1 §54): the H₃⁺ cooling moved inside `eval_cool`, so the **marching temperature update feels it for the first time**, H₃⁺ carries >99% of the radiative cooling from the base to the front (7.3× the local photoheating at r≈1.04) on the metals-off gate. At the 12000-step convention the gate observables are unchanged to <0.1% (front, H₃⁺ peak, Ṁ) because the layer's H₃⁺ cooling time is 200-2000 t_s; T has only begun to fall (−0.5% at r≈1.08). The converged molecular thermal structure is therefore **not** pinned by this gate |
| 2 (Lyman-Werner) | **`Stellar LW flux [erg/cm2/s]:`**, H₂ + hν(912-1110 Å) → H + H in the coupled network, with the temperature-dependent self-shielding of the star-ward H₂ column from Richings, Schaye & Oppenheimer (2014) eqs. (3.12)-(3.15) and 0.4 eV of heating per dissociation (Black & Dalgarno 1977). New module `src/modules/lower_atmosphere/lyman_werner.f90`; diagnostic `output/Lyman_Werner.txt`. Default off (band flux 0) | **done** | section "H₂ photodissociation in the Lyman-Werner bands" below |
| 2 (remaining) | local-equilibrium caveat (no molecular *advection*: Koskinen's high-altitude H₂ replenishment not reproduced); dissociative/double photoionization channels (P4/P5); 4.48 eV dissociation energy sink; diatomic γ; `_adv` post-process; GJ 1214 b He-halving validation (full M-dwarf setup) | open | gates defined in §Tier-2 |
| 2a | `Molecular base: True`, EOS-only base correction (`Molecular chemistry: True` turns this on by itself, since an atomic `ntot_bc` under a molecular base is inconsistent): removes the H₂-bound particles from `ntot_bc` (lower base pressure / heavier base μ; chemistry stays atomic, crude, documented). q_H₂ comes from the photochemical handoff when `base.inp` carries `q_H2_base`, and from the equilibrium fit otherwise (2026-08-10; `composition.f90`) | **done** | HD 209458 b: q_H₂(1 μbar,1450 K)=0.831 → ntot_bc 1.0→0.546 (equilibrium fit) |
| 3 | `base.inp` reader in `input_read` (T_base / r_base / HeH_base / Kzz_base, echo + no-op when absent) + `src/utils/run_lower.py` driver (isothermal or Guillot 2010 semi-grey T(p); writes base.inp with a molecular-base warning). [2026-08-26: the reader also takes `q_H2_base`, `p_base` and one `<El>_H_base` key per element (any of the ten `species_table` symbols), the El/H **nuclei** ratio at the handoff level, which overrides `metals.inp` for that element, turns the metal system on by itself when it is the only nonzero abundance, and on a restart renormalizes that element's whole loaded column onto the stated ratio (`load_IC.f90`). Every key now carries a *category* (provenance / EOS boundary / elemental reservoir / initial guess / boundary constraint), which states what it is allowed to do to the wind; table in `md/input_schema.md` §2c.] | **done (analytic stack)** | end-to-end: driver → base.inp → EXHALE consumes and echoes; iso vs Guillot: r₀ 1.4723 vs 1.4660 R_J, T_base 1450 vs 1313 K (HD 209458 b) |
| 3 (upgrade) | **VULCAN end-to-end**: public VULCAN cloned (`EXHALE/VULCAN/`, FastChem compiled), HD 189733 b SNCHO photochemical run converged (2356 steps); converter `src/utils/vulcan_to_base.py` (.vul → base.inp with photochemical q_H2/q_H, VULCAN-μ/T hypsometric r_base; molecular mixing ratios as comments; **no metal release**, outside VULCAN's scope, stays Lavvas-only). Since 2026-08-10 the converter also writes the read keys `q_H2_base` and `p_base`, so the photochemical H₂ partition **replaces the chemical-equilibrium fit** in the molecular-base particle count instead of being recorded as a comment | **done (H/C/N/O composition)** | HD 189733 b at 1 μbar: **q_H2=0.63, q_H=0.23, photochemistry dissociates ~11× more H than the equilibrium column (q_H=0.020)**, directly quantifying the Tier-1 caveat; r_base 1.168 vs analytic 1.174 R_J; T_base 863 K (Moses11 T(p)); EXHALE consumes the file (overrides echoed) |
| 3 (profile) | [2026-08-27] **`Lower atmosphere profile: <file>`**, the lower atmosphere handed over as a *table over an interval of pressure* instead of the single-level scalars of `base.inp`. Reader `src/modules/files_IO/lower_atmosphere_profile.f90`, schema module `src/utils/lower_profile_schema.py`, producers `src/utils/photochem_to_lower_profile.py` (production path; `--tp-file` or a `--climate` radiative-convective solution) and `src/utils/vulcan_to_lower_profile.py` (cross-check, VULCAN carries no climate model). At its matching pressure the profile sets T₀, R₀, p_base, q_H2, He/H and each elemental reservoir it carries; `K_zz` is taken as a *profile* interpolated onto the grid, which is the point, the homopause is where K_zz crosses the molecular diffusion coefficient and a scalar `He_Kzz` cannot locate it. `n_tot` and `rho` are carried, not imposed. A `base.inp` beside a profile may no longer state those quantities: the EOS-boundary, elemental-reservoir and boundary-constraint keys are refused by name, and the file must carry a `# solution_id` matching the profile's | **done** | `lower_profile` regression case (HD 209458 b: reader, matching-level base state, elemental reservoirs, K_zz(p) in place of the scalar, the accepting branch of the solution_id pairing, and the elemental-flux window statistics); example `examples/17_lower_profile/`; schema `md/input_schema.md` §2d |

Notable physics finding from the Tier-1 gate work: the Visscher **equilibrium** fit keeps the
1 μbar base strongly molecular up to T ≈ 2000 K (fully atomic only above ~2400 K), so for
Teq ~ 1000-2000 K hot Jupiters the *equilibrium* base contradicts the standard atomic
assumption, the dissociation there is photochemical (Moses 2011; Koskinen 2013a), which is
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
| WASP-52 b | 1.270 | 1.437 | 1.615 | 0.84 | input = transit radius → base 0.17-0.35 R_J too deep; **recomputed** (below) |

Guillot T(p) variant: r0 = 1.466 / 1.171 / 1.966 / 1.428 R_J, T_base = 1313 / 1072 /
2136 / 1181 K respectively. Common pattern: every Teq ≲ 1500 K planet has a strongly
molecular *equilibrium* base (q_H2 ≈ 0.83), the atomic base rests on photochemical
dissociation, the quantified motivation for the VULCAN tier. WASP-121 b's production base
radius (2.2075, from the Huang 2023 log-g anchoring) sits ~4-10% above the
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
insensitive to the base misplacement: the bias is in **Ṁ (×1.5)**, not in the line. The
He 10830 scan conclusions therefore appear robust to the base-radius shortcut.

---

## Converged Tier-2 solution (2026-08-13)

Every Tier-2 number quoted above comes from a 12000-step relaxation snapshot: until
now no molecular configuration in the tree had ever reached a steady state
(`examples/15_molecular` ran 135694 steps to `du = 2.9e-2`; it is converged as of
2026-08-13, last section). The first
Newton-converged molecular solutions were obtained on the hot-Uranus gate
(`docs/lower_atmosphere_figs/data_g2/input.inp`, the configuration the
`mol_base_handoff` regression case runs, with `q_H2_base = 0.75`) once the H₃⁺
infrared cooling was inside `eval_cool` (Update_EXHALE_stage1 §54).

**Recipe.** Restart the 12000-step snapshot with `Load IC? True` (which arms the
`du` triggers immediately), a raised hand-off threshold `Solver: Newton 5.0e-2`, and
`Max steps: 150000`. Marching then takes `du` from 2.63 down to ~5e-2 over some
20000 steps (PLM→WENO3 at step 5858, secondary ionization flipped in at step
20556) and the JFNK finish engages at step 22556 and returns `info = 0`. The
descent is not monotonic: `du` reaches 8.8e-2 near step 11500 and bounces back to
7e-1 before resuming. The raised hand-off threshold is what let the run hand over
inside the step budget; whether `du` would have reached the 1e-2 default on its
own was not tested. Metals off: 15 Newton
iterations, `||R||` 1.186e-1 → 8.972e-4, 8 min on 8 threads. Metals on: 8
iterations, 1.156e-1 → 8.804e-4, 9 min. The mass flux `4πr²ρv` comes out flat to
4e-3 (metals off) and 9e-3 (metals on) over `r > 1.2`, against a spread of
2.3e2-3.8e3 in the snapshots it started from.

**The default residual target is too loose for a molecular run.** `||R||` is the
maximum relative residual over the wind, and in a molecular run it is set by the
first two or three cells above the base, where the energy residual is ~10³ times
the value anywhere else. Stopping at the default `||R|| < 1e-3` therefore leaves the
molecular layer still cooling: its cell-by-cell energy residual is then still 0.7×
the local H₃⁺ cooling rate, with the same sign in 173 of 180 cells. Re-solving the
same state with `Resid tol: 1.0e-5` moves the layer again (by 60-70% in T) and
only there does the solution become a fixed point of the procedure: a third solve
from it reproduces T to ≤0.14% and the H₂→H front to five digits while `||R||`
falls to 6.7e-7 (metals off) / 1.4e-6 (metals on), and the mass flux is flat to
2.4e-3. **A molecular run needs `Resid tol: 1.0e-5`; the converged numbers below
are from those solutions.**

| | metals off | metals on (solar C/N/O/Mg/Ca/Na/Fe) |
|---|---|---|
| residual norm at exit | 6.7e-7 | 1.4e-6 |
| H₂→H front (x_H₂ = 0.5) | 1.16037 → **1.03581** | 1.16244 → **1.07974** |
| base T (cell 1, pinned by the BC) | 1213.42 → 1212.86 K | 1213.46 → 1213.45 K |
| H₃⁺ peak | 7.03e4 @ 1.0419 → **3.75e5 cm⁻³ @ 1.0148** | 52.1 @ 1.1355 → **1.01e3 @ 1.0789** |
| T at r = 1.02 | 1304.9 → **299.0 K** | 1312.8 → **928.5 K** |
| T at r = 1.04 | 1305.9 → **192.9 K** | 1316.3 → **701.7 K** |
| T at r = 1.08 | 1302.6 → 1075.5 K | 1319.8 → **208.9 K** |
| log₁₀ Ṁ [g/s], spherical | (not flat) → 10.542 | (not flat) → 10.599 |
| log₁₀ Ṁ [g/s], as reported (Ṁ/2) | 10.58 → 10.24 | 10.58 → 10.30 |
| metal share of nₑ at the base | - | 1.000 → 1.000 (unchanged) |

(The snapshot Ṁ is quoted only because the code prints it; its mass flux varies by
a factor 2e2-4e3 across the wind, so no single number represents it.)

**What the converged solution says physically.** The molecular layer collapses to
190-300 K (far below T_eq = 1140 K) and stops there because the coolant switches
itself off: the H₃⁺ cooling time at r = 1.02 rises from 1.3e7 s in the snapshot to
1.3e10 s in the converged state. With metals the collapse is driven by the
saturated ground-term fine-structure lines instead (C I 609/370 μm at the base,
O I 63/145 μm above it: 89% and 91% of the local cooling at r = 1.002 and 1.040),
and H₃⁺ never exceeds 2.3%. Either way the layer is being cooled by lines treated as
**optically thin**, with no thermal-infrared heating from below and no radiative
equilibrium: the base temperature is pinned in one cell by the boundary condition
and nothing holds the column above it. A real H₂ atmosphere at these column
densities is optically thick in exactly these lines and sits near the
radiative-convective profile. The converged Tier-2 thermal structure should
therefore be read as **the steady state of the model as written, and as a
quantitative statement that the model is missing an infrared escape probability
(or a radiative-equilibrium floor) below the H₂→H front**, not as a prediction of
the temperature of a warm Neptune's lower thermosphere. The wind above the front is much less affected: Ṁ moves by
0.085 dex (metals off) and 0.043 dex (metals on) between the `1e-3` and `1e-5`
solutions, because it is launched above the collapsed layer.

**Why marching alone can never do this.** The H₃⁺ cooling time in the layer is
1.3e7 s at r = 1.02 while the CFL-limited marching step of this configuration is
0.30 s: 4.3e7 steps per cooling time. The layer is out of reach of any marching
budget; it is the Newton solve that reaches it.

**HD 209458 b (2026-08-13): the same recipe on a real planet.** The second
Newton-grade molecular solution is `examples/15_molecular`, the metals-off
HD 209458 b configuration of Gate 1, the one that had run 135694 steps to
`du = 2.9e-2` without ever converging. The three keys above (`Solver: Newton
5.0e-2`, `Resid tol: 1.0e-5`, `Max steps: 150000`) are the whole change: no warm
restart was needed, the run went from the default cold IC to `info = 0` in one
pass, 10 min 25 s on 8 threads. PLM → WENO3 at step 13299, secondary ionization
flipped in at step 23837 when `du` first crossed 5e-2, JFNK engaged at step 25837
and took `||R||` from 4.031e-2 to **9.542e-6 in 272 iterations**. Here too the
`du` descent is not monotonic (it reached 2.4e-1 near step 16000 and bounced
back above 1 before resuming), so the raised hand-off threshold matters as much
as it did on the hot-Uranus gate. The residual is dominated by the same place:
the worst cell is `j = 1-3`, `r ≈ 1.000-1.001` (the two or three cells above the
base) for 267 of the 272 iterations, split about evenly between the momentum and
the energy row.

| | 12000-step snapshot (`data_g1m`) | converged (residual norm 9.5e-6) |
|---|---|---|
| H₂→H front (x_H₂ = 0.5) | 1.01975 | **1.00891** |
| base T (cell 1, pinned by the BC) | 1466.55 K | 1466.36 K |
| T at r = 1.005 / 1.02 / 1.08 | 1909 / 1916 / 2855 K | **401** / 2206 / 7574 K |
| coldest cell below r = 1.3 | 1463.6 K @ 1.0002 | **397.3 K @ 1.0056** |
| H₃⁺ peak | 9.57e4 cm⁻³ @ 1.0004 | **5.18e5 cm⁻³ @ 1.0017** |
| x_H₂ at the base | 0.9964 | 0.9963 |
| log₁₀ Ṁ [g/s], spherical | 10.873 | 10.532 |
| log₁₀ Ṁ [g/s], as reported (Ṁ/2) | 10.572 | **10.231** |
| mass-flux spread over r > 1.05 | 1.23 | **4.2e-3** |

The solution is a fixed point of the procedure: re-solving from it with the same
tolerance returns `info = 0` at `||R|| = 9.0e-6` and reproduces the front to
1.00879 (0.01%), Ṁ to 10.229 (0.002 dex) and the H₃⁺ peak to 0.4%; only the
coldest cell of the collapsed layer still moves (397 → 353 K). Note the reload
does not preserve the residual: `Load IC? True` re-derives the ionization state
from the file, so the second solve starts at `||R|| = 8.7e-4`, not at the 9.5e-6
the first one exited on. The physical
reading is the hot-Uranus reading, on a hot Jupiter: the molecular layer between
the base and the front collapses to ~400 K, H₃⁺ carries **100%** of the radiative
cooling from the base out to r ≈ 1.005 (10⁴ times the local photoheating at the
base, 14 times it at r = 1.002), and the same optically-thin treatment and
missing radiative-equilibrium floor apply. The
front moves inward by 0.011 R_p and Ṁ falls by 0.34 dex relative to the snapshot,
so for HD 209458 b the relaxation-snapshot Ṁ is not a substitute for the
converged one.

**With metals the same planet does not converge.** Running `examples/16_molecular_metals`
(the identical configuration plus solar C/N/O) under the same three keys reaches the
hand-off (PLM → WENO3 at step 15263, secondary ionization at 61158) and the JFNK
then stalls: `||R||` falls from 2.338e-2 to 3.6e-4 over 281 iterations and the line
search finds no descent step for 12 consecutive iterations (`info = 2`, worst cell
`j = 93`, `r = 1.021`, momentum), after which the run reverts to marching. This is
the reverse of the hot-Uranus gate, where the metals-on case converged as readily as
the metals-off one, and it is the same base-adjacent momentum stall, displaced outward to the front. `16` is therefore left with
the plain `Solver: Newton` line and no converged solution.

---

## The infrared field of the lower atmosphere (`Base IR field`, 2026-08-13)

Why the converged Tier-2 molecular layer sits at 190-300 K was left open, with
the optically thin line cooling named as the suspect. The
measurement below says the suspect is the right one for the metals-off case and
the wrong one for the metals-on case, and the fix follows the measurement.

### What the lines actually see (measured on the converged metals-on solution)

Line-center optical depths of the eight ground-term fine-structure lines, from the
converged `Resid tol: 1.0e-5` metals-on solution of the hot-Uranus gate, with the
level populations and Doppler widths the cooling module itself uses. τ_up is the
column from the cell to the top of the domain (the one the escape probability was
built from), τ_dn the column to the bottom of the domain.

| r/R_p | T [K] | τ_up(C I 609) | τ_dn | β | τ_up(O I 63) | τ_dn | β |
|---|---|---|---|---|---|---|---|
| 0.9998 | 1213 | 0.088 | 1.7e-4 | 0.904 | 2.01 | 4.0e-3 | 0.211 |
| 1.0102 | 1033 | 0.069 | 0.019 | 0.924 | 1.57 | 0.44 | 0.265 |
| 1.0302 | 820 | 0.040 | 0.048 | 0.955 | 0.89 | 1.1 | 0.421 |
| 1.0501 | 573 | 0.020 | 0.068 | 0.978 | 0.43 | 1.6 | 0.633 |
| 1.0749 | 289 | 0.005 | 0.083 | 0.995 | 0.093 | 1.9 | 0.898 |
| 1.0832 | 188 | 0.001 | 0.086 | 0.998 | 0.028 | 2.0 | 0.968 |

So **escape was never the problem**: outward the collapsed layer is thin
(β ≥ 0.90 for C I, ≥ 0.63 for O I above r = 1.05), and trapping can only reduce
the cooling by a factor of a few at the base. What the treatment left out is the
other direction. The gas *below* the base is not vacuum. Taking the base cell
opacity and one pressure scale height of an isothermal exponential reservoir
(base pressure 9.0 μbar):

| line | τ per scale height at the base | pressure where τ = 1 | τ at 1 bar |
|---|---|---|---|
| O I 63 μm | 1.08 | 17 μbar | 1.2e5 |
| O I 145 μm | 0.32 | 37 μbar | 3.6e4 |
| C I 370 μm | 0.092 | 106 μbar | 1.0e4 |
| C I 609 μm | 0.046 | 205 μbar | 5.1e3 |
| C II 158, N II 205/122, O I 44 μm | ≤ 1.2e-10 | - | ≤ 3e-2 |

Every line that carries the cooling is black downward within one to three scale
heights below the model base. The layer therefore faces a blackbody, not a void,
and the ratio of the incident mean intensity to the line source function,
J̄/S = (1/2)B_ν(T_base)/B_ν(T), crosses 1 wherever T drops below ~610 K: at the
coldest cell it is 3.4 (C I 609 μm) to 5.7 (O I 63 μm). Those lines were being
made to cool a gas that they should have been heating.

### But that is not what makes the metals-on layer cold

Two measurements say the fine-structure lines are not the cause of the 190-300 K
in the metals-on solution.

- **The specific entropy rises monotonically outward through the whole collapsed
  layer**, p/ρ^γ going from 1.00 at the base to 1.67 at r = 1.05 and 5.69 at
  r = 1.083 (metals-on, normalized to the base cell). The gas is being net
  heated everywhere; it is cold because ρ has fallen by a factor 10². The layer
  is an expansion, not a radiative collapse.
- **The metals-off solution collapses further with no fine-structure cooling at
  all**: 115 K at r = 1.033, where its total cooling is 3.8e-12 erg cm⁻³ s⁻¹.

The energy budget makes the same point directly. The metals-on cooling is 36% of the local
photoheating at r = 1.01 and 2-16% above r = 1.02 (3.5e-9 against 3.5e-8 at
r = 1.05),
and heat − cool is balanced by the expansion term ρv[dε/dr + p d(1/ρ)/dr] to
within a few percent. The one place radiation dominates is the metals-off base:
there H₃⁺ carries 100% of the cooling, the entropy *falls* by a factor 3 between
the base and r = 1.02, and the layer really is radiating itself down.

### The closure

`Base IR field: True` (default `False`) gives both families of infrared coolant
the field they sit in. The lower atmosphere is taken to be black at their
wavelengths and to radiate B_ν(T₀) over the sky fraction
f = 1 − sqrt(1 − (R_p/r)²), i.e. the whole lower hemisphere at the base and the
usual dilution far away.

For the eight fine-structure lines the incident field enters *inside* the
ground-term statistical equilibrium, as the photon occupation number
n̄ = β₁(τ_dn) f / (exp(E_ul/T₀) − 1) with radiative rates βA(1 + n̄) down and
βA(g_u/g_l)n̄ up, so the returned power is the net one, emission minus
absorption. The escape probability becomes two-sided, β = β₁(τ_up) + β₁(τ_dn)
with β₁ the single-face Hollenbach & McKee (1979) / de Jong et al. (1980) form;
with the field off it reduces to 2β₁(τ_up), the previous value, bit for bit. Each
line then stops cooling at its own radiative equilibrium temperature,
575.8 K (C I 609 μm) to 641.7 K (O I 44 μm) for T₀ = 1140 K and half-sky
coverage, and heats below it.

For H₃⁺ there is no line list: Miller et al. (2013) fit the total emission, so
the absorption is taken from that fit at the radiating temperature, which is
what Kirchhoff's law makes of an emission integral:
Λ_net = n(H₃⁺) 4π s(T,n_H₂) [E(T) − W E(T₀)]. This is the same net-exchange form
the H₂, H₂O and CO channels use (§10), and it is bounded at every temperature.
Its radiative equilibrium temperature is 975.4 K at the base for T₀ = 1140 K.
**Between 2026-08-13 and 2026-08-31 this exchange was instead closed on ONE
effective band, the ν₂ fundamental at 2521.3 cm⁻¹ (E/k = 3627.5 K), as
Λ_net = Λ_emit(T)[1 + n̄ − n̄ exp(E/T)], with equilibrium temperatures of 936.1 K
and then 941.2 K.** Pairing a total emission fit with the Boltzmann factor of one
of its transitions diverges below about 150 K and is wrong by a factor 1.5-1.9
over 400-1000 K. The approximations and where they break are
written at `h3p_net_cooling_rate` in `src/modules/lower_atmosphere/h3p_cooling.f90`
and at `fine_structure_line_transfer` in `src/modules/radiation/Cool_coeff.f90`.
Fe II (a precomputed statistical-equilibrium table), the coronal remainders and
every permitted line keep the optically thin, no-incident-field limit; the
estimated heating the remainders would add is ~0.4% of the losses of this layer.

### What it does

Both runs restart the converged `Resid tol: 1.0e-5` solutions of the section
above with the switch on, and both return `info = 0`.

| | metals off | | metals on | |
|---|---|---|---|---|
| | field off | **field on** | field off | **field on** |
| residual norm at exit | 6.7e-7 | 1.3e-6 | 1.4e-6 | 9.9e-6 |
| T at r = 1.005 | 626.8 K | **910.0 K** | 1099.0 K | 1098.9 K |
| T at r = 1.02 | 299.6 K | **859.4 K** | 928.9 K | 928.9 K |
| T at r = 1.03 | 140.2 K | **840.2 K** | 820.0 K | 820.1 K |
| T at r = 1.05 | 476.8 K | **753.0 K** | 573.0 K | 573.4 K |
| coldest cell | 115.1 K @ 1.0334 | **229.7 K @ 1.0947** | 187.9 K @ 1.0831 | 189.1 K @ 1.0831 |
| H₂→H front (x_H₂ = 0.5) | 1.03582 | **1.08833** | 1.07976 | 1.07979 |
| H₃⁺ peak [cm⁻³] | 3.75e5 @ 1.0148 | 3.19e5 @ 1.0083 | 1.01e3 @ 1.0789 | 1.00e3 @ 1.0789 |
| log₁₀ Ṁ [g/s], spherical | 10.539 | 10.617 | 10.596 | 10.596 |
| mass-flux spread, r > 1.2 | 2.4e-3 | 2.4e-3 | 2.4e-3 | 2.4e-3 |

**Stale, and not re-measurable from here: the field-on columns and the predicted
floors below were computed with the single-band H₃⁺ closure, which was replaced
on 2026-08-31.** The converged solutions these
restart from no longer exist as run directories, so the table is left as the
record of what that closure gave. The same four-way comparison re-measured with
the current closure, at the 12000-step protocol of the `mol_ir_bands` regression
case, is in §112; re-running §55's converged runs needs those solutions rebuilt.

The metals-off layer stops collapsing and settles **on** the H₃⁺ radiative
equilibrium curve: the predicted floor with the local sky fraction is 911 K at
r = 1.005, 886 K at 1.02 and 874 K at 1.03, and the solution sits at 910, 859 and
840 K, just below, by the margin the expansion takes out. Its H₃⁺ channel is a
net heating term (−1.1e-7 erg cm⁻³ s⁻¹ at r = 1.02) instead of the −1.6e-10 of
cooling it had. The front moves out by 0.053 R_p and Ṁ by +0.078 dex. It is a
fixed point of the procedure like the solutions above: re-solving from it returns
`info = 0` at `||R|| = 8.8e-6` and reproduces T below r = 1.3 to 1e-4.

The metals-on layer does **not** move (0.1-1.2 K), for the reason measured above:
its temperature is set by the expansion, and the fine-structure lines it does
carry are 5-20% of the local heating. They now correctly turn into a heating term
above r ≈ 1.06, which is physically right and numerically almost invisible.

### What is still open

The metals-on molecular layer, and the outer part of the metals-off one past the
H₂→H front (229.7 K at r = 1.095), remain far below any radiative equilibrium
temperature, and this closure cannot reach them: it only gives the *existing*
infrared coolants their incident field, and above the front there are no
molecular coolants left. Holding that gas requires a continuum infrared coupling
to the deep atmosphere (H₂ collision-induced absorption and the H₂O/CH₄/CO
bands), which the model does not have: the whole radiative budget between the
base and the front is line channels with a radiative time of 3e8 s against a flow
time of 5e6 s. The missing molecular infrared cooling of the layer below the
front is therefore narrowed, not closed.

## H2 photodissociation in the Lyman-Werner bands (`Stellar LW flux`, 2026-08-13)

The Tier-2 network inherited from Koskinen et al. (2022) Table 1 has no
photodissociation of neutral H2: its photo-rates start at the 15.4 eV
photoionization edge, so below that the only H2 losses are thermal (R12),
electron impact (R14) and ion chemistry. Their own note says the omission is
deliberate and that adding it (Backx et al. 1976 cross section, dissociation
probability 0.125) moved their Ṁ by ≤ 1.4×. The consequence for us is that the network wants a base more molecular than
any photochemical code gives, and pinning `q_H2` at the base (Route 1) patches
over the missing physics rather than supplying it. This section supplies it.

### The band is not in the code's own radiation field

Checked in the source, not assumed. `set_energy_vectors` builds the photon grid
from 13.6 eV up; it extends below the H I edge only when the He 2³S metastable
(4.8 eV) or a low-IP metal is active, and then only to that species' threshold.
`read_sed` selects SED rows by the same `e_low`. So for a numerical SED the
11.2-13.6 eV band is simply not read, and for a power-law SED (every Tier-2 gate
case) anything below 13.6 eV is the XUV power law extrapolated downward, which
has no relation to a star's FUV. Even where the grid does reach into the band,
its opacity is continuum photoionization, whereas Lyman-Werner absorption is a
forest of saturated lines whose attenuation is nothing like exp(−τ_continuum).
**The band flux is therefore a separate user input, and the transfer is done
separately.**

```
Stellar LW flux [erg/cm2/s]: 343.0    # 912-1110 A, integrated, at the planet
```

Default 0 = off. With `Molecular chemistry` off the key is inert and says so.

### Rate

Draine & Bertoldi (1996), ApJ 468, 269 (DB96), published version. Their band is
912-1110 Å: 912 Å is the H Lyman edge, and longward of 1110 Å the H2 absorptions
out of v = 0 are negligibly weak (their footnote 4). They characterize a field
by the photon flux in that band, F ≡ c n_phot (their eq. 21), and tabulate it
with the unshielded dissociation rate. For the flat-F_λ spectrum (u_ν ∝ ν^−2,
their eq. 24) at χ = 1: F = 1.208e7 photons cm⁻² s⁻¹ (Table 1), ζ_pump =
3.09e-10 s⁻¹ with ⟨p_diss⟩ = 0.135 (Table 2), hence ζ_diss(0) = 4.17e-11 s⁻¹
(Fig. 7 caption). The dissociation rate per band photon is then an effective
cross section

    sigma_LW = 4.17e-11 / 1.208e7 = 3.452e-18 cm^2 ,

and with the mean photon energy of a flat-F_λ band, ⟨hν⟩ = 2hc/(912 + 1110 Å) =
12.2635 eV,

    k_LW,thin = 1.757e-7 * F_LW    [s^-1, F_LW in erg cm^-2 s^-1] .

Reproducing DB96's own Table 1 numerically from their eqs. (22)-(24) confirms
the flux convention (Habing 1.2220e7 against their 1.222e7; ν^−2 1.2084e7
against 1.208e7; Draine 1978 1.2313e7 against 1.232e7) and the closure of the
calibration (the formula returns 4.17e-11 s⁻¹ for the ν^−2 field at χ = 1, their
value to 0.04%).

**The approximation and its range.** σ_LW depends on the shape of the spectrum
*within* the band, because the Lyman/Werner lines sample it unevenly. The same
arithmetic for the much softer Draine 1978 field (ζ_pump = 2.78e-10,
⟨p_diss⟩ = 0.119, F = 1.232e7 at χ = 1) gives 2.685e-18 cm², 22% lower, and the
formula above would overpredict that field's rate by 26%. The two spectra
bracket color temperatures 1.3e4-2.9e4 K, and DB96 note that PDR properties are
insensitive to the spectrum for T_color ≳ 1e4 K. **Take the adopted value as
good to ±25% for a stellar FUV band that is not strongly tilted**; a band
dominated by a single emission line at one end is outside it.

### Self-shielding

DB96 eq. (37), their fit to the full multiline calculation *including line
overlap*:

    f_shield(N_H2) = 0.965/(1 + x/b5)^2
                   + 0.035/(1+x)^0.5 * exp[-8.5e-4 (1+x)^0.5] ,
    x = N_H2 / 5e14 cm^-2 ,   b5 = b / 1e5 cm s^-1 ,

and their eq. (40), ζ_diss = f_shield e^{−τ_dust} ζ_diss(0).

**The rate no longer uses that `f_shield`.** DB96's fit is calibrated for cold
gas in which only the lowest rotational states of H2 are populated; Richings,
Schaye & Oppenheimer (2014, MNRAS 442, 2780) compared it against a
level-resolved CLOUDY calculation and found it overestimates the shielding
factor by a factor ~3. Our molecular layer runs at 945-1307 K, where the
overestimate measured with our own columns is 2.8-5.4×. The photodissociation
rate therefore uses their temperature-dependent eqs. (3.12)-(3.15), with the
same thermal Doppler parameter (their fit was made for purely thermal
broadening, so no turbulent term is added). DB96 eq. (37) is kept, and is still
used, for the *band equivalent width* of the shared 912-1110 Å beam: their
eq. (39) is the closed-form integral of eq. (37) and of nothing else, and
Richings publishes no corresponding equivalent width. `lyman_werner.f90` §2 sets
out why the two coexist and which quantity each owns.

`b` is the H2 Doppler parameter, taken thermal, b = (2kT/m_H2)^{1/2}: 3.07 km/s
at 1140 K, next to the 3 km/s of DB96's own figures.

`N_H2` is the star-ward column, built by the same `calc_column_dens_one` radial
integration and `opa_pf` weighting as every other absorber column, from the
incoming (pre-solve) H2 density: the same lagging the photoionization columns
use.

Three things are deliberately *not* attenuating the band, all noted at the code:

- **dust**: EXHALE's metals are atomic and trace, so there are no grains and the
  e^{−τ_dust} of eq. (40) is identically 1;
- **trace-metal continuum**: the neutral low-IP metals do photoionize inside the
  band, but at solar abundance and σ ~ 1e-18 cm² their optical depth is ~4e-23
  N_H, i.e. ≲ 0.05 at a base column where f_shield is already below 1e-4;
- **H Lyman-series lines** (Lyβ 1025.7, Lyγ 972.5, … all lie in the band): DB96
  include them in the equivalent width their fit was built on and state, §4.2,
  "We will see below that absorption by the H Lyman lines has only a small
  effect on the H2 pumping rates". They are not treated separately here. This is
  the least controlled of the three, because a wind has a much larger N_HI/N_H2
  than the PDRs DB96 fitted; it can only reduce the rate.

### Chemistry and heating

H2 + hν → H + H enters the H2 balance row of the coupled system next to the H2
photoionization (`mol_heh_rows`, `System_HeH_mol`), so it is solved together
with everything else rather than in a side loop. Both products are neutral H,
which the H-nucleus closure supplies automatically; no other row changes. The
system is solved by `hybrd1` with a numerical Jacobian, so there is no analytic
derivative to update.

The fragments carry kinetic energy. Black & Dalgarno (1977), ApJS 34, 405,
p. 418: "Fluorescent dissociation of H2 gives rise to a pair of energetic
hydrogen atoms (Milgrom, Panagia, and Salpeter 1973; Stephens and Dalgarno
1973); for a typical ultraviolet radiation field, the yield is about 0.4 eV per
atom pair, but it varies slightly with depth." We adopt 0.4 eV per dissociation,
held constant with depth, added to `heat` in `ioniz_eq` next to the Penning
heating terms. **The 4.48 eV bond energy is paid by the absorbed photon, not by
the gas, and is not a sink of this channel** (a thermal dissociation-energy sink
for R12/R14 remains a separate open item).

`output/Lyman_Werner.txt` is written whenever a molecular run carries a band
flux: r, T, x_H2, n_H2, N_H2, f_shield, k_LW and the photodissociation heating.

### The band flux for a planet

The code's own SED cannot supply it, so it is integrated externally. The VULCAN
stellar spectra shipped with the tree (`EXHALE/VULCAN/atm/stellar_flux/`, flux
at the *stellar surface*, erg cm⁻² s⁻¹ nm⁻¹) integrated over 91.2-111.0 nm and
diluted by (R_star/a)²:

| star / planet | surface band flux | R_star, a | F_LW at the planet | k_LW,thin |
|---|---|---|---|---|
| Sun, `Gueymard_solar.txt`, at 1 AU | 2.74e4 | 1 R_sun, 1 AU | 0.592 | 1.04e-7 s⁻¹ |
| HD 209458 b (solar spectrum as proxy) | 2.74e4 | 1.155 R_sun, 0.048 AU | **343** | **6.03e-5 s⁻¹** |
| HD 189733 b, `sflux-HD189_Moses11.txt` | 4.23e4 | 0.805 R_sun, 0.03142 AU | 600 | 1.05e-4 s⁻¹ |
| HD 189733 b, `sflux-HD189_B2020.txt` | 1.29e5 | 0.805 R_sun, 0.03142 AU | 1.83e3 | 3.21e-4 s⁻¹ |

(erg cm⁻² s⁻¹ throughout.) The `VULCAN_run_*/atm/stellar_flux/` directories of
the workspace are empty (casualties of the 2026-08-09 loss), so no HD 209458
spectrum exists in the tree and the solar spectrum stands in for it; HD 209458
is a G0 star of comparable activity, but **343 erg cm⁻² s⁻¹ is a proxy, not a
measurement**. The calibration check on the same file is the solar Lyα
irradiance it returns at 1 AU, 6.12 erg cm⁻² s⁻¹ against an observed 6-8. The
factor 3 between the two HD 189733 spectra is a fair measure of how uncertain a
real stellar FUV band flux is, and it dwarfs the ±25% of the band-shape
approximation above.

### What it does: hot-Uranus gate, A/B on the key alone (metals off)

Both runs restart the `mol_base_handoff` configuration (`q_H2_base = 0.75`,
`p_base = 1e-5` bar) with `Solver: Newton 5.0e-2`, `Resid tol: 1.0e-5`,
`Max steps: 150000` and `Base IR field: True`, and differ only in
`Stellar LW flux`. Both return `info = 0`.

| | LW off | **LW on (343)** |
|---|---|---|
| residual norm at exit | 8.15e-6 | 4.86e-6 |
| H₂→H front (x_H2 = 0.5) | 1.08806 | **1.08622** |
| x_H2 at the base | 0.99921 | 0.99883 |
| q_H2 at the base | 0.8618 | 0.8612 |
| x_H2 at r = 1.05 | 0.9777 | 0.9629 |
| x_H2 at r = 1.10 | 2.54e-2 | **1.51e-4** |
| H₃⁺ peak | 3.19e5 @ 1.0079 | 3.20e5 @ 1.0025 |
| T at r = 1.005 / 1.05 / 1.10 | 901.5 / 752.3 / 255.6 K | 886.6 / 751.8 / 249.9 K |
| log₁₀ Ṁ [g/s], spherical | 10.6196 | 10.6248 |
| mass-flux spread, r > 1.2 | 2.38e-3 | 2.39e-3 |

The self-shielding works as designed and is visible in the dump:

| r | N_H2 [cm⁻²] | f_shield | k_LW [s⁻¹] |
|---|---|---|---|
| 1.000 (base) | 4.32e21 | 9.79e-7 | 5.90e-11 |
| 1.020 | 1.59e21 | 4.31e-6 | 2.60e-10 |
| 1.050 | 3.07e20 | 2.29e-5 | 1.38e-9 |
| 1.087 (front) | 6.23e18 | 2.85e-4 | 1.72e-8 |
| 1.100 | 3.43e14 | 4.67e-1 | 2.82e-5 |
| ≥ 1.15 | ≤ 7e9 | 1.000 | 6.026e-5 |

i.e. the band is unattenuated in the wind, where it recovers exactly the
unshielded 6.03e-5 s⁻¹ computed above, and suppressed by 10⁶ at the base.

### Route 2 verdict: the physics was missing, but it does not close the gap

The Route 2 criterion was whether the
network's own base composition moves toward the photochemical partition once the
missing dissociation is supplied, without pinning it. **It does not**, and the
measurement says why.

The base q_H2 moves from 0.8618 to 0.8612 against a photochemical target of 0.75
(the value that case's `base.inp` carries): 0.5% of a gap of 0.11. Yet
Lyman-Werner *is* the largest single H2 loss term at the base after the change:
5.90e-11 s⁻¹ against 2.52e-11 for H⁺ + H2 (R10 + R13), 5.56e-12 for He⁺ + H2 and
3.87e-13 for thermal dissociation. The two facts are consistent because the base
partition is a formation-destruction balance in which the three-body reaction
R15 (H + H + M) sets n_H2 ∝ n_H², so

*(Note added 2026-09-17: R15's third body is no longer the total
heavy-particle density with the M = H2 coefficient but the collider sum over
H2, H and He with the published coefficients of Cohen & Westberg 1983, which at
this base lowers the association rate by 22 per cent; the balance argument
below is unchanged in form.)*


    n_H / n_H2  ∝  (total H2 destruction rate)^(1/2) ,

and the square root is brutal: the measured 2.3× rise in the destruction rate
raises the atomic fraction by only 1.5×, from 7.9e-4 to 1.2e-3. Reaching
q_H2 = 0.75 means x_H2 = 0.925, i.e. an atomic fraction 64× larger, i.e. a
destruction rate **4100× larger**: about 3.7e-7 s⁻¹. The unshielded band would
supply that many times over (6.0e-5 s⁻¹); it is stopped by the column. The base
of this planet sits under N_H2 = 4.3e21 cm⁻², where f_shield = 9.8e-7, and
f_shield reaches the required ~6e-3 only near N_H2 ~ 1e17-1e18 cm⁻², four orders
of magnitude shallower.

**So Lyman-Werner photodissociation cannot be the reason photochemical codes
find more atomic H at 1 μbar.** That reading is not new physics: this document
already records it in §1, "photochemistry (OH-catalyzed H₂ destruction)
dissociates H₂ even at 500-1000 K where thermal equilibrium would not". The
catalytic cycles run on O, OH and H₂O, and EXHALE's molecular network is
H2/H2+/H3+/HeH+ with atomic metals; it has nowhere to put them. Closing the base
gap needs those species, not a stronger radiation field. Route 1 (pinning
`q_H2_base`) therefore remains the only way EXHALE can carry the photochemical
base partition today, and it is now a documented modeling choice rather than a
patch over a missing term.

What Lyman-Werner *does* change is the top of the molecular layer, where the
column has fallen enough for the band to bite: at r = 1.10 the H2 fraction drops
by a factor 168 (2.5e-2 → 1.5e-4), the front moves inward by 0.0018 R_p, and
Ṁ rises by 0.005 dex. Small, and in the direction and of the order Koskinen
et al. (2022) reported for their own sensitivity test.

### The same A/B on HD 209458 b

The `examples/15_molecular` configuration (metals off, `Base IR field: True`
added, both runs restarted from its converged solution) reproduces the pattern on
a real planet. Here the JFNK finish would not reach `Resid tol: 1.0e-5` with the
band on (it stalled three times at ‖R‖ = 1.5e-5 to 6.8e-5 on the base-adjacent
momentum row) so **both** members of the pair were
run at `Resid tol: 2.0e-5`, where both return `info = 0` and both mass fluxes are
flat to 4.1e-3 over r > 1.2. The pair still differs only in `Stellar LW flux`.

| | LW off | **LW on (343)** |
|---|---|---|
| residual norm at exit | 1.77e-5 | 1.96e-5 |
| H₂→H front (x_H2 = 0.5) | 1.01578 | **1.01290** |
| x_H2 at the base | 0.99628 | 0.99605 |
| q_H2 at the base | 0.8512 | 0.8509 |
| x_H2 at r = 1.02 | 2.37e-4 | **4.22e-6** |
| H₃⁺ peak | 4.81e5 @ 1.0010 | 5.46e5 @ 1.0010 |
| T at r = 1.005 / 1.02 | 937.7 / 986.3 K | 792.7 / 1275.8 K |
| log₁₀ Ṁ [g/s], spherical | 10.5711 | 10.5600 |

The base is again untouched (N_H2 = 2.6e21 cm⁻², f_shield = 2.1e-6), the front
moves in by 0.0029 R_p, the H2 fraction just above it falls by 56×, and Ṁ moves
by −0.011 dex. The chemical-equilibrium fit gives q_H2 = 0.831 at this base and
photochemistry would give less still; the network sits at 0.851 with the band on
as without it.

## The molecular infrared bands (`Molecular IR bands`, 2026-08-31)

The `Base IR field` closure of the section above gave the infrared coolants the
code already had (the eight ground-term fine-structure lines and the H₃⁺ bands)
the field they sit in, and it stopped there. Its own closing paragraph named
what was left: past the H₂→H front there are no molecular coolants in the model
at all, and holding that gas needs the coolants a real H₂ atmosphere carries.
This section supplies three of them and measures what they do.

### Which coolants, and why these three

| coolant | where it comes from | why it is here |
|---|---|---|
| H₂ quadrupole + magnetic dipole lines | Roueff et al. (2019), A&A 630, A58, table 2 (VizieR J/A+A/630/A58) | H₂ is the bulk gas. Its transition probabilities are 10⁻¹⁰-10⁻⁶ s⁻¹, but it outnumbers every other molecule by 10³-10⁴, and the two cancel: at n(H₂) = 5×10¹³ cm⁻³ and 1000 K it radiates 1.4×10⁻⁷ erg cm⁻³ s⁻¹, four times the local photoheating (the pure rotational part alone is 4.8×10⁻⁸; above ~600 K the vibrational bands carry most of it). It is the only one of the three that exists in a Tier-2 run without the oxygen option. |
| H₂O vibration-rotation bands | HITEMP 2010 through the correlated-*k* coefficients distributed with Photochem | The strongest infrared coolant of a solar-composition H₂ atmosphere below ~2000 K, and the one the oxygen option switches the [O I] fine-structure lines off in favour of. |
| CO vibration-rotation bands | HITEMP 2019, same route | Holds essentially all the carbon at the base; the 4.7 μm fundamental carries four to five orders of magnitude more than the pure rotational lines above 300 K. |

H₂-H₂ and H₂-He collision-induced absorption is deliberately **not** a channel,
and that is a measurement rather than a scope decision. It scales as n², and on
the CIA tables shipped with Photochem its optical depth over a pressure scale
height at the 1 μbar base is 8×10⁻⁹ (hot Uranus) to 2×10⁻⁹ (HD 189733 b), with
unit optical depth reached only near 0.2-0.4 bar: consistent with Lavvas &
Arfaux (2021), who put CIA's onset at p > 1 bar. It is what makes the reservoir
*below* the base black in the windows between the bands, which is the assumption
`Base IR field` already encodes, not a local coolant of the modelled domain.

### The closure

For each species the net radiative loss per unit volume is the optically thin,
LTE exchange with the diluted field of the lower atmosphere:

```
Lambda_net = n_X [ E_X(T) - W E_X^abs(T; T0) ]
```

with, for a band absorber,

```
E_X(T)          = 4 pi Sum_b sigma_b(T) Int_b B_nu(T)  dnu
E_X^abs(T; T0)  = 4 pi Sum_b sigma_b(T) Int_b B_nu(T0) dnu
```

and, for the H₂ line sum, the same quantity assembled transition by transition
with the full two-level net factor `1 + nbar - nbar exp(dE/T)`. `W` is the same
dilution the H₃⁺ closure uses, `½ [1 − sqrt(1 − (R_p/r)²)]`, and it is zero when
`Base IR field` is off.

**Keeping the stimulated-emission term is not decoration.** It is what makes the
bracket vanish identically at `T = T0, W = 1`: gas buried in a blackbody at its
own temperature neither cools nor heats. For the band absorbers that property is
automatic, because the HITRAN line intensities the cross sections are built from
already carry the `1 − exp(−h nu/kT)` factor, so `sigma_b` is the coefficient
that pairs with `B_nu`. For a line sum it has to be written down. The H₃⁺
closure of the previous section had dropped it; it was put back the same day,
and because `Ek(nu2)/T0` is large that cost only 2.2% of the emission and 5 K of
equilibrium temperature there (936 → 941 K at `T0 = 1140 K`). On a far-infrared
band the same omission would be an order-unity error, which is why the new
module carries it explicitly. **The single-band H₃⁺ form the 941 K belongs to
was itself replaced the same day** (a total emission fit cannot be paired with
the Boltzmann factor of one transition at all temperatures) and the H₃⁺ column
below is the replacement's.

Each channel therefore has its own radiative equilibrium temperature, and it is
the fixed point rather than the magnitude that answers item (G):

| `T0` | H₂O | CO | H₂ | (H₃⁺, for comparison) |
|---|---|---|---|---|
| 900 K | 705 K | 752 K | 783 K | 780.9 K |
| 1140 K | 889 K | 918 K | 984 K | 975.4 K |
| 1183 K | 922 K | 947 K | 1019 K | 1010.7 K |
| 1450 K | 1120 K | 1120 K | 1226 K | 1231.6 K |
| 1800 K | 1369 K | 1338 K | 1486 K | 1517.6 K |
| 2358 K | 1749 K | 1668 K | 1883 K | 1938.6 K |

(half-sky coverage, `W = ½`). How much rests on that assumption is measurable:
for H₂O at `T0 = 1183 K` the equilibrium temperature is 726.5 / 922.1 / 1065.0 /
1183.0 K at `W = 0.25 / 0.5 / 0.75 / 1`, and CO 784.7 / 946.5 / 1073.1 /
1183.0 K, at `W = 1` both are exactly `T0`, the closure's fixed point recovered
numerically. Because `W` multiplies only the absorption, these
temperatures are insensitive to first order to any grey escape probability
applied to both terms, which is why the optically thin limit is enough to fix
the *temperature* even where a band is marginally thick.

### Optical depth, and why the thin limit is used

At a 1 μbar base with solar oxygen the band-mean optical depth of the H₂O
rotational bands over a scale height is a few times 10⁻²; the strongest
*g*-points inside a band reach order unity. The H₂ lines are thinner still: the
line-centre opacity of 0-0 S(1) at 17 μm gives τ ≈ 3×10⁻⁴ over a scale height at
n(H₂) = 5×10¹³ cm⁻³ and 1000 K. No escape probability is applied, and the
Planck-mean optical depths of the H₂O and CO columns to the top of the domain
are written to `output/Cooling_breakdown.txt` so the assumption is measured
rather than asserted. `output/FUV_bands.txt` carries the same statement as a
conservation check: the column-integrated absorbed power of the three bands is
compared against the incident infrared flux `W sigma T0⁴`, and the run warns if
it exceeds it.

### Where the numbers come from, and three independent checks

The band-mean cross sections are `Sum_g w_g k_g` of the correlated-*k* tables,
which is exactly the bin-mean cross section. Line strengths do not depend on the
broadening, only the line shape does, so this mean is pressure independent to
better than 2% from 10⁻⁶ to 1 bar (measured) and the whole EXHALE domain uses
the 10⁻⁶ bar slice, which is also the lowest tabulated pressure. The tables run
50-2000 K and are clamped outside; the H₂ line sum has no such limit.

Three checks, against sources that share none of that line processing:

- The Planck-mean H₂O cross section agrees with the independent table
  `aiolos/inputdata/1H2-16O_T400.aiopa` to within 11% over 600-2000 K (ratio
  0.889-1.111 on the 14 shared grid temperatures). Below 600 K that table is
  not monotonic in T and the comparison is not meaningful; the low-temperature
  end is covered by the Neufeld & Kaufman comparison instead.
- Restricted to λ > 9 μm, i.e. to the pure rotational band, the H₂O emission
  agrees with the optically thin LTE rate of Neufeld & Kaufman (1993), ApJ 418,
  263, table 2 to a factor 0.70 (100 K), 0.81 (200 K), 0.90 (400 K), 1.02
  (1000 K). The same cut on CO against their table 3 gives 0.81-0.93 over
  300-2000 K. At 100 K the CO ratio falls to 0.10; the cause is most likely
  the eight-point *g*-quadrature resolving the linear bin mean poorly when the
  spectrum inside a bin is a few isolated narrow lines, which CO's rotational
  ladder is at that temperature: it is not the coarse long-wavelength binning,
  since 97% of the 100 K emission lands in the eighteen 9-250 μm bins. It does
  not matter: CO radiates 1.6×10⁻²⁰ erg s⁻¹ per molecule there against
  1.8×10⁻¹⁵ for H₂O.
- The H₂ line sum agrees with the LTE rates of Hollenbach & McKee (1979),
  ApJS 41, 555, eqs. (6.37) and (6.38) to 1-11% over 200-2000 K.

Two numbers from those comparisons are the reason the whole band system is
integrated instead of a published *rotational* cooling function being adopted:
above 300 K CO's total emission is four to five orders of magnitude above its
rotational part, and H₂O's is 3.2× its rotational part at 1000 K. A
rotation-only coolant would have missed almost all of it. Neufeld & Kaufman's
tables also stop at an optical depth parameter of 10¹⁹ cm⁻² (km/s)⁻¹, a ceiling
set by hydrostatic equilibrium in a self-gravitating cloud, which does not apply
to a planetary atmosphere and which a base column exceeds.

The generator, the provenance and the validity statements are
`cooling_data/molecular_infrared_bands.py`; the tables it writes are
`src/modules/lower_atmosphere/molecular_infrared_data.f90`; the physics is
`src/modules/lower_atmosphere/molecular_infrared_cooling.f90`.

### What it does

Four converged solutions of the hot-Uranus gate (`data_g2` restarted with
`Load IC? True`, `Solver: Newton 5.0e-2`, `Resid tol: 1.0e-5`,
`Max steps: 150000`), all with `Base IR field: True`, so each pair differs only
in `Molecular IR bands`. This configuration has no oxygen chemistry, so the only
new channel is the H₂ line sum. The bands-off members reproduce the field-on
solutions of the section above to 1-3 K (metals on) and 5-30 K (metals off, the
difference being the H₃⁺ stimulated-emission fix).

| | metals off, off | metals off, **on** | metals on, off | metals on, **on** |
|---|---|---|---|---|
| JFNK iterations | 89 | **21** | 25 | 26 |
| residual norm at exit | 8.77e-6 | 8.39e-6 | 2.00e-6 | 2.82e-6 |
| T at r = 1.005 | 940.5 K | 912.1 K | 1097.0 K | 1044.1 K |
| T at r = 1.02 | 878.6 K | 857.6 K | 926.1 K | 887.0 K |
| T at r = 1.05 | 750.3 K | 721.9 K | 572.2 K | 588.7 K |
| T at r = 1.08 | 420.9 K | 345.4 K | 215.1 K | 236.1 K |
| coldest cell below 1.3 | 245.6 K @ 1.0972 | 225.4 K @ 1.0923 | 190.4 K @ 1.0842 | 198.2 K @ 1.0853 |
| H₂→H front | 1.08817 | 1.08532 | 1.07916 | 1.07970 |
| log₁₀ Ṁ, as reported | 10.32 | 10.31 | 10.29 | 10.30 |
| mass-flux spread, r > 1.2 | 2.39e-3 | 2.38e-3 | 2.38e-3 | 2.38e-3 |

**The solution's own zero crossing lands on the predicted equilibrium
temperature.** In the metals-on solution the H₂ channel changes sign at
r = 1.0122, where the gas is at 951.7 K and the sky dilution is W = 0.423; the
closure evaluated at that dilution gives 951.5 K. The second crossing (r =
1.1239, T = 850.6 K, W = 0.272) is predicted at 872.7 K. Nothing was fitted:
the tables come from HITEMP and Roueff et al., the dilution from the geometry.

The H₂ channel is **first order and does not repair the collapse**, and those
are not in tension. In the metals-on solution H₂ carries 56% of the whole
radiative rate at r = 1.002-1.005 (cooling there, since T = 1044-1124 K is
above its radiative equilibrium temperature) and essentially all of the net at
r = 1.02-1.05, heating. The layer's radiative time falls with it:
`t_rad/t_flow` at r = 1.05 goes from **8.22 to 0.263**. But the coldest cell
sits *at or beyond* the H₂→H front in every run, where there is no H₂ to absorb
anything and `t_rad/t_flow` is 0.7-1.9; that gas is expansion-dominated and no
molecular coolant can reach it.

With the oxygen chemistry on, where H₂O and CO exist, the magnitudes are
different by five orders: on HD 189733 b the two bands carry 4.8e-5 and 2.7e-5
erg cm⁻³ s⁻¹ at r = 1.005 against a total of 7.5e-5, and a 12000-step A/B moves
the base from 760.2 to 829.1 K and `x_H2` from 0.330 to 0.483, both toward the
photochemical reference (864 K, 0.910). Neither of those runs is converged.

### What is assumed

- **LTE level populations.** The rotational levels of all three species
  thermalize far below the densities of this layer. The vibrational levels are
  the constraint (n_crit(H₂ v = 1) is 10⁹-10¹¹ cm⁻³ against 10¹³-10¹⁴ cm⁻³ at
  the base), but the vibrational bands are a small share of the H₂ emission
  where the density is low: 0.7% at 400 K, 5.7% at 500 K, 18% at 600 K, 66% at
  1000 K. The LTE overestimate is therefore confined to the warm, dense part of
  the layer, where LTE is also best justified.
- **The incident field is not attenuated by the intervening molecular column**
  of the domain itself, the same approximation `Base IR field` makes for H₃⁺.
- **OH is left out.** It has no cross-section table here and carries a few
  percent of the oxygen where the water carries most of it.
- The atomic advection post-process (`T_equation`) does not carry these
  channels, exactly as it does not carry H₃⁺.
