# Extending the Cherubim et al. (2026) LHS 1140b analysis to the lower atmosphere

**Status: superseded, kept as a record.  The plan of record is
`lhs1140b_lower_atmosphere_plan_new.md` (2026-08-22), written after the
external review `lhs1140b_lower_atmosphere_plan_review.md`; that file's
section 1 lists the claims of this draft that were corrected against the
source (normalization scale, diffusion generalization, crossover-mass
claim, SED-path history, transit metrics, Photochem adapter status,
molecular reaction set, fluid-validity wording).  Do not quote this draft.**

Source papers: `references/Cherubim_2026Science.pdf` and
`references/Cherubim_2026Science_Supplement.pdf` (C. Cherubim et al., Science,
10.1126/science.aea9708, first release 2026-07-16).  All numbers quoted below
are from those PDFs (read with `pdftotext -layout`), not re-derived.

Related plan documents: `oxygen_chemistry_new_plan.md` (phases P1/P2 are
prerequisites for WP2 below), `base_composition_handoff_plan.md`,
`lower_atmosphere_coupling.md`, `design_hehe_diffusion.md`.

---

## 1. What the paper computed and what it only estimated

The observation: metastable He 10830 absorption from LHS 1140b (5.60 +/- 0.19
M_Earth, 1.730 +/- 0.025 R_Earth, P = 24.7 d, Teq = 226 +/- 4 K, inactive
M dwarf host at 14.96 pc) — detected in 2024, not detected in 2025
(detection limit ~0.6%).

**Computed (p-winds, upper atmosphere only).**  A 1-D isothermal Parker wind
forward model fit by MCMC to the 2024 transmission spectrum:

- excess absorption 1.24 +0.22/-0.23 % (blended red pair),
  0.25 +0.14/-0.12 % (blue line); FWHM 0.86 +0.15/-0.27 A;
  equivalent opaque radius 1.52 R_p;
- mass-loss rate 2.03 +0.58/-0.67 x 10^8 g/s, T_wind ~ 5160 K (grid to
  6400 K), fiducial F_XUV = 0.033 W/m^2 (= 33 erg/cm^2/s at the orbit);
- **H:He number ratio <~ 1e-3.**  The physical driver: at this low XUV flux,
  H at >~ 1% would absorb the <= 911 A photons that ionize ground-state He,
  suppressing metastable He production.  The wind is helium-dominated.
- Numerical stability of p-winds itself was checked down to H:He ~ 1e-8
  (Supplement, "Retrieving atmospheric properties").

**Estimated analytically (no model was run below the wind).**

- Fractionation: crossover mass < 9 amu at any T < 10,000 K
  (Hunten, Pepin & Walker 1987 closed form); minimum escape rate to drag
  atomic O is ~2 x 10^9 g/s — so O, C, N stay behind and accumulate.
- Cold trap: T_skin = 2^(-1/4) Teq = 194 K, tropopause pressure assumed
  0.1 bar, giving f_H2O ~ 7 ppm at the tropopause — the water (and hence
  hydrogen) supply to the upper atmosphere.
- Origin of the He-dominated envelope: cited from a magma-ocean +
  escape-fractionation evolution model (their ref. 39), not computed.
- Variability: steady-state p-winds models at 1-33% of the fiducial XUV
  bracket the 2025 non-detection.

The missing link is exactly the lower atmosphere: **where the trace hydrogen
comes from, how much of it reaches the wind, and how diffusive separation
sets the H:He ratio that p-winds had to leave as a free parameter.**

## 2. Goal

Replace the paper's free parameters and closed-form estimates with one
self-consistent EXHALE chain:

```
radiative-convective lower atmosphere (cold trap -> f_H2O)
   -> photochemical H production (Photochem, base.inp handoff)
   -> diffusive separation through the homopause (H trace in He background)
   -> He-dominated EXHALE wind (energy balance, He ionization, He 2^3S)
   -> He 10830 forward model (EXHALE_transit.py)
   -> compare: depth 1.24%, FWHM 0.86 A, Mdot ~ 2e8 g/s, H:He <~ 1e-3
```

The retrieval targets H:He and Mdot become **outputs**; the inputs become the
stellar SED and the bulk-atmosphere composition.

## 3. Why EXHALE cannot run this today — the three gaps

### Gap 1: the He-dominated regime is untested

EXHALE normalizes everything to hydrogen nuclei (`He/H number ratio` key,
`input_read.f90:133`; `mass_per_H = 1 + 4*HeH`).  Entering HeH ~ 1e3 is
formally possible but the regime is unexercised:

- solver paths that divide by n_H must survive the trace-H limit (p-winds
  needed an explicit stability check there; expect the same here);
- the cooling budget is H-centric (Ly-alpha is the dominant coolant in every
  validated run).  In an H-poor gas the energy equation must close on He
  channels — He recombination, free-free, He 2^3S lines.  The HeITR module
  exists; whether it closes the budget alone has never been checked;
- the electron budget flips from H-ionization- to He-ionization-dominated;
  `System_HeH_*` must be audited in that limit.

### Gap 2: the diffusion formulation is written the wrong way around

`src/modules/functions/species_diffusion.f90` diffuses He (and trace metals)
**relative to a hydrogen background** — the kernel, the base normalization
`fbase = nX/nH`, and the settling term all reference n_H.  For LHS 1140b the
roles invert: H is the trace species, He the background.  The binary
diffusion coefficient (Banks & Kockarts) is symmetric, so the physics
carries over, but the bookkeeping (background density, element normalization,
settling sign for a *lighter*-than-background species: H floats up, not
settles) must be generalized.

Payoff: crossover mass and O/C/N drag become code results instead of the
Hunten 1987 closed form — this is the same diffusion generalization already
listed for oxygen in `oxygen_chemistry_new_plan.md` phase P3.

### Gap 3: the lower boundary assumes a hot H2 atmosphere

Tier 1 (`lower_atmosphere/lower_column.f90`) is the Koskinen et al. (2022)
analytic column: isothermal at Teq with a chemical-equilibrium H2/H/He
partition — built for hot Jupiters.  At Teq = 226 K the relevant lower
atmosphere physics is different:

- a radiative-convective structure with a cold trap sets f_H2O (~7 ppm),
  and photolysis of that water is the hydrogen source for the wind — this
  is the physical input that should reproduce H:He <~ 1e-3;
- the natural route is the existing Tier-3 handoff: run **Photochem**
  (Wogan et al. 2025, `references/Wogan_2025_Planet._Sci._J._6_256.pdf`;
  working clone at `../photochem/`) for a temperate He-dominated
  atmosphere and pass the result through `base.inp`.  Photochem, unlike
  VULCAN, is built for temperate planets and does its own climate/cold-trap
  handling.  This makes phases P1 (matched Photochem comparison) and P2
  (elemental handoff contract) of `oxygen_chemistry_new_plan.md`
  prerequisites, not parallel work;
- in a He-dominated gas the Tier-2 molecular network's HeH+ channels
  (`mol_rates.f90` R16-R21) change from side reactions to main chemistry
  and need re-examination in that limit.

## 4. Work packages

### WP0 — He-dominated validity audit (before any new physics)

Stress-test the existing code at HeH >> 1: ionization equilibrium, energy
closure, electron budget, and NaN/overflow behavior of every n_H-normalized
path (run under `-fcheck`).  Deliverable: a memo stating which routines are
valid in the trace-H limit and which needed guards.  Acceptance: a converged
pure-attempt run at HeH ~ 1e3 with a closed energy budget, or a precise list
of what breaks.

### WP1 — Generalize diffusive separation (H trace in He background)

Rewrite the element-diffusion kernel in `species_diffusion.f90` so the
background is the dominant species (or, cleaner, the total gas) rather than
n_H.  Verify against the analytic diffusive-equilibrium profile in an
isothermal test column, both orientations (He-in-H and H-in-He).
This is shared work with oxygen plan P3.

### WP2 — Temperate lower boundary via Photochem handoff

Run Photochem for LHS 1140b (He-dominated bulk, Teq 226 K, cold-trapped
water); extend the `base.inp` handoff contract to carry what this case
needs (elemental H inventory at the base, base temperature/pressure from
the radiative-convective solution).  Depends on oxygen plan P1/P2.
Acceptance: the handoff reproduces f_H2O ~ ppm at the model base and the
resulting EXHALE wind returns H:He in the retrieved range without tuning.

### WP3 — Stellar inputs

Construct the XUV SED as the paper did: scaled GJ 1132 and GJ 699 spectra
plus the XMM-Newton X-ray constraint (Supplement, "X-ray observations and
analysis"); feed it through the existing `Spectrum file` / `sed_read` path.
Note: every current regression case uses a power-law spectrum, so this run
is the first real exercise of the SED path in production — check the
resolved configuration in the setup report.

### WP4 — Validation ladder and science runs

1. **p-winds oracle** (clone at `../p-winds/`): reproduce the paper's
   retrieval point (Mdot 2.03e8 g/s, T 5160 K, H:He 1e-3) locally — same
   oracle pattern as Wind-AE.
2. EXHALE steady-state run with the same inputs; compare density/velocity/
   metastable-He profiles against the p-winds Parker solution.
3. He 10830 forward model with `EXHALE_transit.py` (He line already
   implemented); compare depth, FWHM, blue/red ratio.
4. XUV grid at 1-33% of fiducial to bracket the 2025 non-detection, as in
   the paper's Fig. 4.
5. **Fluid-validity check (a posteriori, mandatory):** F_XUV is ~3 orders
   of magnitude below the hot-Jupiter runs EXHALE is validated on.  Compute
   the Knudsen number / exobase location on the converged solution; if the
   sonic point sits above the exobase the hydrodynamic Mdot is an
   overestimate and the result must say so.

Suggested order: WP0 -> WP1 and WP3 in parallel -> WP4 steps 1-3 with a
*prescribed* H:He first (decouples the wind work from WP2) -> WP2 ->
WP4 with the computed H:He -> WP4 steps 4-5.

## 5. Explicitly out of scope (1-D)

- The leading tail (and tentative trailing tail) asymmetry — 3-D geometry.
- Stellar-wind interaction and true time-dependent response of the outflow.
- The Gyr-scale fractionation history that produced the He-dominated
  envelope (the paper cites an evolution model for this; we take the
  present-day composition as a boundary condition).

The paper itself treated all three qualitatively; a 1-D steady-state chain
reproducing the 2024 line and bracketing the 2025 non-detection is the
deliverable.

## 6. Relation to the oxygen chemistry plan

This target is a driver, not a competitor, for `oxygen_chemistry_new_plan.md`:
WP2 *is* phases P1/P2 exercised on a real science case, and WP1 is the P3
diffusion item.  The genuinely new physics opened here is the He-dominated
limit (WP0) — everything else reuses planned work.
