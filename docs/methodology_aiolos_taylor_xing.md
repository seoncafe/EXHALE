# Methodology Comparison: AIOLOS, Taylor et al. (2025, 2026), and Xing et al. (2023)

**Purpose.** This note compares the methodology of three independent modeling efforts for
XUV-driven atmospheric escape from close-in exoplanets, and assesses which elements could
be adopted into **EXHALE** (our extended ATES fork). The three efforts are:

- **AIOLOS** — the general-purpose 1-D multi-species radiation-hydrodynamics code of
  Schulik & Booth (2023, "SB23"), distributed on GitHub (`github.com/Schulik/aiolos`).
- **Taylor et al. (2025)**, *ApJ* 989:68 — "A Multispecies Atmospheric Escape Model with
  Excited Hydrogen and Helium: Application to HD209458b," built on the Koskinen et al.
  (2013a,b; 2022) thermosphere–ionosphere code.
- **Taylor et al. (2026)**, *ApJ* 999:214 — "Helium Escape in Context: Comparative
  Signatures of Four Close-in Exoplanets," a direct application/extension of the Taylor
  (2025) model to four planets (HD 209458b, HD 189733b, HD 149026b, GJ 1214b).
- **Xing et al. (2023)**, *ApJ* 953:166 — "The Mass Fractionation of Helium in the Escaping
  Atmosphere of HD 209458b," a **multi-fluid** model built inside the PLUTO code.

Throughout, "Taylor" refers to the shared Koskinen-based framework used in both Taylor
papers unless a specific paper is named. EXHALE's own methodology is summarized in
Section 2 so the adoption analysis (Section 6) has a concrete baseline.

---

## 1. One-paragraph characterization of each code

**AIOLOS** is a *general-purpose* radiation-hydrodynamics code, not a dedicated escape
code. It solves the time-dependent Euler equations (mass, momentum, energy) for an
arbitrary number of species — each a **separate fluid** with its own density, momentum,
and internal energy — coupled by an inter-species **friction/drag** solver. It uses
finite-volume **HLLC** fluxes with **PLM (MC-limited)** reconstruction and second-order
TVD Runge–Kutta time-stepping, on Cartesian/cylindrical/**spherical** log-spaced grids.
Radiation is treated with **flux-limited diffusion (FLD)** over arbitrary in/out spectral
bands, solved implicitly via a block-tridiagonal system; chemistry ranges from a fast
C2Ray ionization scheme to a full implicit multi-reaction network. Escape is one of many
problems it can run (alongside dust formation, mini-Neptune interiors, tidal heating).

**Taylor (2025/2026)** is a *dedicated* 1-D spherical thermosphere–ionosphere +
hydrodynamic-escape model (Koskinen et al. 2013/2022), integrated to steady state, with a
**single bulk momentum equation** but **multi-species molecular + eddy diffusion** that
redistributes species (so He and H can separate diffusively). Its distinguishing strengths
are physical completeness of the *microphysics*: a full He(2³S) and H(n=2) non-LTE network
solved by the KPP kinetic preprocessor, a **Lyα Monte Carlo** radiative-transfer model
(Huang et al. 2017) iterated with the hydro solution for Hα, a new high-resolution He(2³S)
photoionization cross section (B-spline K-matrix), temperature-dependent Penning ionization,
and self-consistent coupling to a photochemical lower/middle atmosphere at the μbar level.

**Xing (2023)** is a *dedicated* 1-D spherical escape model whose distinguishing feature is
that it is genuinely **multi-fluid**: H, H⁺, He, He⁺, and e⁻ each carry their **own
velocity**, coupled by resonant/non-resonant collisional drag. This is the only one of the
three that can capture **He/H mass fractionation** (heavy He under-dragged by escaping H)
self-consistently in the dynamics rather than through a diffusion coefficient. It is built
on PLUTO (HLL Riemann solver, 3rd-order TVD RK) with a novel Riemann split to carry the
electron/electric-field pressure term. The He(2³S) 10830 Å line is post-processed, not
solved in the hydro loop.

**EXHALE (baseline for comparison)** is a 1-D spherical **single-fluid** steady-state
escape code (ATES heritage), with a two-stage PLM→WENO3 reconstruction and an approximate
Riemann solver, integrated to steady state on a mass-flux (ρvr²) convergence criterion.
Ionization of H, He, and trace metals (C/N/O) is solved as one coupled nonlinear system
(analytic-Jacobian Newton with a MINPACK `hybrd1` fallback). It already includes a He(2³S)
metastable network, an H(n=2)/Balmer non-LTE population (Christie et al. 2013), a Lyα
radiative-transfer suppression via a Neufeld escape-probability treatment, Badnell RR+DR
recombination, Voronov collisional ionization, Kingdon & Ferland charge exchange, CHIANTI
metal-line cooling (Mg/Ca/Na/Fe II, C/N/O), a Roche/tidal potential (RLOF), and a
transmission-spectrum post-processor (`EXHALE_transit.py`).

---

## 2. Side-by-side methodology table

| Axis | AIOLOS (SB23) | Taylor 2025/2026 | Xing 2023 | **EXHALE (current)** |
|---|---|---|---|---|
| Purpose | General multi-species RHD | Dedicated escape + line diagnostics | Dedicated escape, He/H fractionation | Dedicated escape + line diagnostics |
| Dimensionality | 1-D (sph/cyl/cart) | 1-D spherical | 1-D spherical | 1-D spherical |
| Time | Time-dependent → steady | Time-dependent → steady | Time-dependent → steady | Steady-state (relaxation) |
| Fluid model | **Multi-fluid** (v, ρ, E for each species) + friction | **Single bulk** momentum + multispecies diffusion | **Multi-fluid** (v for each species) | **Single fluid** (one bulk v) |
| He/H separation | Via friction/drag | Via molecular + eddy diffusion | Via multi-fluid dynamics (fractionation) | **None** (fixed He/H) |
| Reconstruction | PLM + MC limiter | (finite-difference, Koskinen 2013) | PLM (PLUTO) | **PLM → WENO3** (two-stage) |
| Riemann solver | HLLC | — | HLL (+ electron-pressure split) | Approximate (ATES/PWN) |
| Time integrator | 2nd-order TVD RK | forward to steady | 3rd-order TVD RK | RK relaxation |
| Grid | Log-spaced | radial | 1024 pts, stretched (1→10 Rp) | Stretched radial |
| Radiation transport | **FLD**, multi-band, implicit block-tridiag | XUV attenuation + Lyα Monte Carlo | Radial XUV attenuation (53 bins) | Radial XUV attenuation |
| Heating efficiency | From C2Ray energy split | Photoelectron efficiency (tuned 20–40%) | Frequency-averaged η (tuned 0.1–0.5) | Tuned efficiency |
| Ionization solve | C2Ray (Brent) or implicit network | KPP kinetic preprocessor | Semi-implicit source terms | **MINPACK / analytic-Jacobian Newton** |
| Species (fiducial) | Arbitrary (gas/dust/ions) | H, H(n=2), He, H⁺, He⁺, He²⁺, He(2³S), e⁻ (+metals, +H₂) | H, H⁺, He, He⁺, e⁻ | H, He (all stages), C/N/O, e⁻, He(2³S), H(n=2) |
| Metals | Via general chemistry | Optional (solar abundances) | None | **C/N/O ionization + Mg/Ca/Na/Fe cooling** |
| He(2³S) 10830 | Not a design focus | In-loop non-LTE network | Post-processed (Yan 2022) | In-loop metastable network + `EXHALE_transit.py` |
| H(n=2) / Hα | Not a focus | Non-LTE, Lyα Monte Carlo, iterated | Post-processed | Non-LTE (Christie 2013) + Lyα escape prob. |
| Cooling | Recomb/line/free-free (Black 1981), dust, H₃⁺ | Recomb + H I line (Huang 2023) + H₃⁺ | **Lyα only** (+ adiabatic, frictional) | Lyα + recomb + free-free + CHIANTI metals |
| Tidal / Roche | `USE_TIDES` quadrupole | Optional (Taylor 2025); off in 2026 | Stellar tidal term in a_ext | Roche potential (RLOF) |
| Lower boundary | Configurable | **μbar, coupled to photochemistry** | 1 Rp, fixed n, T=1500 K | ~μbar base, fixed T/n |
| Transmission spectrum | — | Ray-traced Voigt (He 10830, Hα) | Post-processed 10830 | `EXHALE_transit.py` (10830, Lyα, Hα/β, metal doublets) |

---

## 3. What they have in common

1. **1-D spherical, XUV-driven escape** with photoionization heating balanced by
   Lyα/recombination cooling and adiabatic expansion is the shared physical backbone of
   all four (AIOLOS when configured for escape, Taylor, Xing, EXHALE).
2. **Steady-state target.** All reach a stationary transonic wind — the three external
   codes by evolving the time-dependent equations forward; EXHALE by direct relaxation.
   The physical end state (a transonic outflow through a sonic point) is the same.
3. **Tunable heating efficiency.** None derives the photoelectron heating efficiency from
   first principles; all treat it as the primary free knob (Taylor 20–40%, Xing 0.1–0.5,
   EXHALE similar), calibrated against He 10830 / mass-loss constraints.
4. **Same core atomic physics for H/He.** Voronov collisional ionization, radiative +
   dielectronic recombination, and H↔He charge exchange appear in Xing, Taylor, and EXHALE
   with near-identical rate coefficients (differing mainly in the source table: Storey &
   Hummer vs Badnell vs Verner).
5. **He(2³S) as the key observable.** Taylor, Xing, and EXHALE all target the He 10830 Å
   metastable line (Taylor and EXHALE additionally target Hα); all three build a non-LTE
   2³S population from recombination, collisions, radiative decay, Penning ionization, and
   photoionization. (AIOLOS does not target this line specifically.)
6. **Finite-volume, TVD-limited, high-resolution shock-capturing hydrodynamics** underlies
   AIOLOS, Xing (PLUTO), and EXHALE (WENO3). Taylor's Koskinen code is the outlier
   (finite-difference thermosphere solver).

---

## 4. Where they genuinely differ

### 4.1 Fluid model — the single most important axis
- **Xing is multi-fluid** (each species its own velocity). This is what lets it reproduce
  **He/H mass fractionation**: because He is 4× heavier, it is under-dragged by the H wind,
  so the *escaping* He/H (0.039) is roughly half the bulk abundance (0.086) — reproducing
  the observed weak 10830 without invoking a subsolar He abundance.
- **Taylor is single bulk momentum but adds molecular + eddy diffusion**, which produces
  the *same* qualitative diffusive He/H separation (8% → 2.5% with altitude) at far lower
  cost than a full multi-fluid solve. Taylor explicitly argues this reproduces Xing's and
  Schulik & Owen (2025)'s multi-fluid results.
- **AIOLOS is multi-fluid but escape-agnostic**: its friction solver would in principle
  yield fractionation, but the code is not tuned for the transonic-wind escape problem.
- **EXHALE is strictly single-fluid** with a *fixed* He/H at all radii — it **cannot**
  represent fractionation or diffusive separation at all. This is the largest physical gap
  relative to the modern He 10830 literature.

### 4.2 Radiation transport
- **AIOLOS** is the only one with a full **flux-limited-diffusion** thermal-radiation
  solver (block-tridiagonal implicit, multi-band) — necessary for its optically-thick
  interior/dust problems but overkill for the optically-thin escaping thermosphere.
- **Taylor** is the only one that does **resonant-line (Lyα) Monte Carlo** radiative
  transfer, iterated with the hydro solution, to set the H(n=2)/Hα population accurately.
- **Xing and EXHALE** use simple **radial XUV optical-depth attenuation** for the driving
  radiation. EXHALE additionally treats Lyα trapping through a Neufeld escape-probability
  factor — a middle ground between Xing's neglect and Taylor's full Monte Carlo.

### 4.3 Ionization/chemistry solver
- **Taylor**: KPP kinetic preprocessor (large stiff ODE network, general).
- **AIOLOS**: C2Ray (fast, quasi-steady) or a full implicit reaction matrix.
- **Xing**: semi-implicit source terms folded into the PLUTO RK stages.
- **EXHALE**: a single coupled **nonlinear algebraic** ionization system (Newton with
  analytic Jacobian, MINPACK fallback) — appropriate for a steady-state code, and already
  including trace metals inside the coupled system.

### 4.4 Metals and cooling
- **EXHALE** has the richest *closed-form* metal-cooling treatment (CHIANTI fits for
  Mg/Ca/Na/Fe II and C/N/O), and solves C/N/O ionization in the coupled system.
- **Taylor** can include solar-abundance metals but treats the escape structure as H/He
  dominated; metals mainly lower the required heating efficiency.
- **Xing** has **no metals** and **only Lyα cooling** — deliberately minimal.
- **AIOLOS** carries cooling through its general chemistry (Black 1981 line/recomb/free-free,
  H₃⁺), but again not specialized to the metal-line escape problem.

### 4.5 Line diagnostics
- **Taylor** self-consistently computes both He 10830 and Hα in-loop, with the most
  complete non-LTE microphysics (temperature-dependent Penning ionization, B-spline He 2³S
  cross section, proton de-excitation, photoelectron impact excitation).
- **Xing** post-processes 10830 only (Yan et al. 2022 NLTE), decoupled from the hydro.
- **EXHALE** sits between: it builds the 2³S and H(n=2) populations *in the run* and
  produces the full transmission spectrum (10830, Lyα, Hα, Hβ, metal doublets) in `EXHALE_transit.py`.
- **AIOLOS** does not target these lines.

### 4.6 Lower boundary
- **Taylor** couples to a full photochemical lower/middle atmosphere and places the escape
  base at ~1 μbar — the most physically motivated boundary of the four.
- **Xing** and **EXHALE** use a fixed-density, fixed-temperature base near 1–2 μbar
  (1 Rp / T ≈ 1500 K for Xing).
- **AIOLOS** leaves the boundary fully configurable.

---

## 5. Adoption assessment for EXHALE

The three efforts offer very different things to EXHALE. EXHALE is already a *dedicated,
line-diagnostic-capable* escape code, so it is closest in spirit to **Taylor** and **Xing**,
and furthest from **AIOLOS** (whose value is a general RHD architecture EXHALE does not
need). Below, each candidate element is rated by **physics payoff** and **implementation
cost**, given EXHALE's single-fluid steady-state architecture.

### 5.1 High payoff

**(A) He/H diffusive separation / mass fractionation — the top priority.**
This is EXHALE's biggest physical gap and directly controls the He 10830 depth EXHALE is
built to predict. Two adoption routes:
- **Taylor route (recommended):** add a **molecular + eddy diffusion** term that lets the
  He/H ratio vary with radius while keeping the single bulk-momentum solve. This is a
  bolt-on to EXHALE's existing single-fluid framework (a diffusion flux in the species
  continuity equations), moderate cost, and Taylor shows it reproduces the multi-fluid
  answer. It also naturally supplies the altitude-dependent He/H that `EXHALE_transit.py` needs.
- **Xing route (high cost):** convert EXHALE to genuinely multi-fluid (velocities
  for each species + collisional drag + electron-pressure Riemann split). This is the physically
  purest treatment of fractionation but is effectively a rewrite of the hydro core and its
  steady-state relaxation, and is hard to reconcile with EXHALE's coupled-nonlinear
  ionization solve. **Not recommended** unless fractionation in the *transonic/decoupled*
  regime becomes a primary research goal.

**(B) Temperature-dependent Penning ionization rate (Taylor 2025).**
Taylor replaces the long-used temperature-independent He(2³S)+H Penning rate
(~5×10⁻¹⁰ cm³ s⁻¹) with a Maxwell–Boltzmann-averaged, temperature-dependent fit (from
Morgner & Niehaus 1979 and Cohen & Lane 1971). Penning ionization is often the **dominant
2³S loss near the base**, so this directly affects EXHALE's 10830 prediction. **Low cost,
high value** — a drop-in replacement in the He(2³S) network. Worth checking whether EXHALE
currently uses the constant rate.

**(C) High-resolution He(2³S) photoionization cross section (Taylor 2025, B-spline).**
Taylor's B-spline K-matrix 2³S cross section resolves sharp EUV autoionization resonances
that overlap strong stellar coronal lines. EXHALE has **already** moved in this direction
(the recent TOPbase-based high-energy extension + two-VFKY96 + bridge in `cross_sec.f90`);
Taylor's machine-readable cross section is a natural cross-check / possible upgrade for the
resonance region. **Low cost** (data swap), moderate value (resonance overlap is
star-dependent).

### 5.2 Medium payoff

**(D) Whole-atmosphere lower-boundary coupling (Taylor).**
Placing the escape base at ~1 μbar with composition/temperature from a photochemical
lower/middle atmosphere (rather than a fixed base) removes an ad hoc boundary assumption.
For EXHALE this would mean coupling to (or importing profiles from) a separate
photochemistry model. **Moderate-to-high cost**; valuable mainly for sub-Neptune / molecular
regimes (Taylor 2026's GJ 1214b), less so for the hot-Jupiter cases EXHALE currently targets.

**(E) Molecular chemistry (H₂, H₂⁺, H₃⁺, HeH⁺) and H₃⁺ cooling (Taylor 2026).**
Needed only if EXHALE is extended toward **cooler / higher-μ sub-Neptunes** (GJ 1214b-like),
where H₂ survives into the wind and H₃⁺ becomes a significant coolant. **Moderate cost**;
defer unless the science scope broadens beyond hot Jupiters.

**(F) Lyα Monte Carlo (Taylor / Huang 2017) vs EXHALE's Neufeld escape probability.**
Taylor's iterated Monte Carlo is the gold standard for the Hα-driving Lyα field, but EXHALE
already approximates Lyα trapping with a Neufeld escape probability (chosen precisely because
the WASP-121b τ₀ ≈ 5×10⁷ makes brute-force Monte Carlo infeasible). Adopting full Monte Carlo
is **high cost for incremental gain** given EXHALE's targets — keep the escape-probability
approach, but Taylor's setup is a useful validation benchmark.

### 5.3 Low payoff for EXHALE

**(G) AIOLOS's FLD radiation, HLLC solver, and general chemistry architecture.**
These are excellent engineering for a *general* RHD code but solve problems EXHALE does not
have: FLD matters for optically-thick interiors/dust, not the optically-thin escaping
thermosphere; HLLC vs EXHALE's approximate Riemann solver is not a limiting factor for the
smooth transonic wind; and EXHALE's coupled-nonlinear ionization solve is already
well-suited to steady state. **Not recommended** as targeted adoptions. AIOLOS remains
valuable as an **independent cross-check** (its multi-species friction solver is an
alternative reference for fractionation physics) rather than a source of components.

**(H) Xing's electron-pressure Riemann split.**
Elegant, but only relevant if EXHALE goes multi-fluid (route A/Xing), which is not
recommended. **Skip** unless (B)-multi-fluid is pursued.

### 5.4 Recommended adoption order

1. **Temperature-dependent Penning ionization rate** (B) — cheapest, directly improves 10830.
2. **He/H diffusive separation via a diffusion term** (A, Taylor route) — largest physics
   gain for EXHALE's core observable; keeps the single-fluid architecture.
3. **Cross-check / optionally upgrade the He(2³S) cross section** against Taylor's B-spline
   data (C) — EXHALE is already most of the way there.
4. Consider **whole-atmosphere boundary coupling** (D) and **molecular chemistry** (E)
   *only if* the science scope extends to sub-Neptunes.
5. Treat **AIOLOS and Xing multi-fluid** as **validation references**, not component donors.

---

## 6. Summary

- **AIOLOS** is a general-purpose, multi-fluid RHD code with the most sophisticated
  *numerics* (HLLC, FLD, implicit chemistry) but the least specialization to the He/Hα
  escape-diagnostic problem EXHALE cares about. Best used as an independent cross-check.
- **Taylor (2025/2026)** is the closest match to EXHALE's goals and the richest source of
  *microphysics* upgrades: temperature-dependent Penning ionization, a high-resolution
  He(2³S) cross section, single-fluid diffusive He/H separation, non-LTE H(n=2)/Hα, and
  whole-atmosphere coupling. Most of EXHALE's near-term gains come from here.
- **Xing (2023)** provides the physically purest treatment of **He/H mass fractionation**
  via genuine multi-fluid dynamics, but at an architectural cost EXHALE should avoid;
  Taylor's diffusion approximation captures most of the same effect far more cheaply.
- **EXHALE's single largest gap** is the absence of any radius-dependent He/H separation.
  Adopting a **Taylor-style diffusion term** (plus the **temperature-dependent Penning
  rate**) is the highest-value, lowest-risk path to bringing EXHALE's He 10830 predictions
  in line with the current state of the art, without abandoning its single-fluid,
  steady-state, metal-aware design.
