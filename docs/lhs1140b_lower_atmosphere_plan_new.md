# LHS 1140b lower-atmosphere extension — plan of record (revised)

**Status: plan of record, 2026-08-22.  Supersedes
`lhs1140b_lower_atmosphere_plan.md` (the first draft) after the external
implementation review `lhs1140b_lower_atmosphere_plan_review.md`; both stay
as records.  Every disputed claim below was re-verified against the source
tree, not taken from either document.**

Source papers: `references/Cherubim_2026Science.pdf` and
`references/Cherubim_2026Science_Supplement.pdf` (C. Cherubim et al., Science,
10.1126/science.aea9708, first release 2026-07-16), read with
`pdftotext -layout`.  Related plans: `oxygen_chemistry_new_plan.md` (phases
P1/P2 are prerequisites for Phase E below),
`base_composition_handoff_plan.md`, `lower_atmosphere_coupling.md`,
`design_hehe_diffusion.md`.

---

## 1. Verified baseline

Where the first draft and the review disagreed, the code was read.  The
review was right on every checked point; one of its counter-claims needed a
nuance (row 6).

| # | First-draft claim | What the code actually does | Verified at |
|---|---|---|---|
| 1 | "EXHALE normalizes everything to hydrogen nuclei" | The hydrodynamic density scale is **H+He nuclei**: `n0` is an H+He scale with `n_H = n0/(1+HeH)`, and `comp_rho_bc = comp_mass_per_H()/(1+HeH)` — so `HeH = 1000` does **not** blow up the density normalization (rho/(m_H n0) -> ~4). What remains H-referenced is the species algebra: Newton initial guesses, advection-correction ratios, and H-relative diagnostics. | `input_read.f90:961-975`, `composition.f90:188-192`, `ionization_equilibrium.f90:596-632`, `post_process_adv.f90:438-491` |
| 2 | Cooling is "H-centric"; He channels "may need" to close the budget | The He channels already exist: He II/III recombination, He I/II/2^3S collisional ionization and excitation, He free-free with the charge factor, He 10830 and singlet conversion. The open question is not existence but whether their **sum closes the energy equation** in a helium-rich low-XUV solution, and whether the photon-escape assumptions hold there. | `util_ion_eq.f90:604-705` |
| 3 | Diffusion kernel "diffuses trace species against an n_H background" — generalize by swapping the background | Confirmed trace-form: conservative face flux `J = n_He v - D n_H (df/dr + f G)` with lagged n_H, post-solve rescale, and a hard cap `n_He/n_H <= HeH`; metals reuse the same kernel independently, no feedback on the bulk. **But swapping in "the total gas" is not a valid fix**: it does not enforce `sum_i rho_i w_i = 0` (equal-and-opposite binary fluxes) and can corrupt the bulk mass flux. A conservative binary (Maxwell-Stefan-consistent) formulation with a defined bulk velocity is required. | `species_diffusion.f90:136-166,197-255` |
| 4 | (not in first draft) | **Molecular chemistry and `He_diffusion` are mutually exclusive** — the parser `error stop`s on the combination. The planned chain (molecular lower boundary -> diffusive separation -> wind) cannot run until this is resolved by design, not by flipping the gate. | `input_read.f90:902-906` |
| 5 | WP1 makes crossover mass and O/C/N drag "code results" | Too strong. The metal loop is independent trace diffusion against H with no momentum feedback on the H/He bulk; reversing the H/He orientation does not create a multicomponent drag solution. Crossover mass stays an **analytic diagnostic** until a multicomponent transport package exists (not planned here). | `species_diffusion.f90:197-255` |
| 6 | "Every current regression case uses a power-law spectrum, so this is the first production use of the SED path" | Half right, half wrong. The `backup/regression/*` cases are all `Spectrum type: Power-law` (checked), **but** `benchmarks/wasp52` runs `Spectrum type: Load from file..` with `../../WASP-52b/eps_eri_sed_fxuv1p0.txt` and has a logged converged solution — the external-SED path is already exercised in production. What LHS 1140b adds is a *new spectrum*, not a new code path. | `benchmarks/wasp52/input.inp:11-12`, `backup/regression/*/input.inp` |
| 7 | (not in first draft) | When He 2^3S is on, `sed_read` extends the selected band down to the triplet threshold **4.8 eV** (~2590 A). An "XUV" SED that stops at the H I edge silently omits the photons that photoionize metastable He and directly biases the 10830 A population. The LHS 1140b spectrum must cover 4.8 eV - X-ray. | `sed_read.f90:37-38` |
| 8 | WP4 compares transit depth/FWHM/blue-red ratio | The transit script reads planet parameters from `input.inp` **only** — never `base.inp` — so a handoff that moves the base radius changes the wind but not the transit geometry. Measured on `benchmarks/wasp52`: input.inp 1.27 R_J / 1304 K / HeH 0.0204 vs. base.inp 1.40492 R_J / 1181.2 K / 0.0959. Also, `FWHM_HeTR` is the **instrument** LSF width, and no astrophysical FWHM, blue-component depth, windowed red depth, or blue/red ratio extractor exists. These metrics require new code. | `EXHALE_transit.py:80-93,398`, `exhale_transit_lib.py:307-331`, `benchmarks/wasp52/base.inp` |
| 9 | Lower boundary via "a Photochem `base.inp` handoff" | No Photochem adapter exists. The automatic route recognizes `Lower atmosphere: analytic | vulcan | none` only, and `read_base_inp` carries exactly six scalars (`T_base, r_base, HeH_base, Kzz_base, q_H2_base, p_base`), ignoring unknown keys. The contract cannot carry a water abundance, an H production rate, an elemental flux, or any profile. | `input_read.f90:1183-1201,1228-1280` |
| 10 | Molecular re-examination "R16-R21" | The molecular residual stores R16-R20 and **R23**; **R21/R22 (H-He charge exchange) are deliberately excluded** and come in through the shared `charge_exchange` path to preserve the atomic limit. The audit set is R16-R20, R23, the shared CX path, and He 2^3S + H2 Penning. | `System_HeH_mol.f90:21-25,84-89` |
| 11 | "If the sonic point sits above the exobase the hydrodynamic Mdot is an overestimate" | Overclaimed, and no Knudsen/exobase diagnostic exists in the source (grepped). The defensible statement: the continuum solution is then **unvalidated through its critical point**; sign and size of the error need a kinetic or transitional-flow comparison. The diagnostic itself (collision model, species-resolved mean free paths, exobase and critical-point definitions) is new code. | `src/` grep for knudsen/exobase: empty |

Also verified: the `q_H2_base`/`HeH_base` same-solution requirement is only
*printed*, not enforced (`input_read.f90:1206-1217`), and the flux-closure
argument of the review stands on physical grounds: a steady lower-atmosphere
chemistry solution cannot determine the hydrogen supply independently of its
upper-boundary flux, so without a flux-continuity condition at a matching
surface, `HeH_base` remains a prescribed input and "the wind returns H:He"
is circular.

## 2. Target

Cherubim et al. (2026) detect metastable He 10830 from LHS 1140b
(5.60 +/- 0.19 M_Earth, 1.730 +/- 0.025 R_Earth, P = 24.7 d,
Teq = 226 +/- 4 K, inactive M dwarf at 14.96 pc) in 2024, not in 2025
(limit ~0.6%).  Their p-winds retrieval (isothermal Parker wind, upper
atmosphere only): depth 1.24 +0.22/-0.23 % (blended red), 0.25 +0.14/-0.12 %
(blue), FWHM 0.86 +0.15/-0.27 A, Mdot 2.03 +0.58/-0.67 x 10^8 g/s,
T_wind ~ 5160 K, **H:He <~ 1e-3**, fiducial F_XUV = 33 erg/cm^2/s.
Everything below the wind — crossover mass < 9 amu, cold-trap water
f_H2O ~ 7 ppm at a 0.1 bar tropopause (T_skin = 194 K), the He-dominated
origin — is closed-form estimation.

Goal: make the chain

```
lower atmosphere (climate + cold trap -> f_H2O)
  -> photochemical H production            (Photochem, new adapter)
  -> elemental-H flux matched at a defined surface   (new closure)
  -> H trace diffusing through the He background     (new formulation)
  -> He-dominated EXHALE wind                        (audited regime)
  -> He 10830 metrics                                (new extractor)
```

self-consistent, so H:He and Mdot become outputs and the inputs are the
stellar spectrum and the bulk-atmosphere composition.  The closure condition
(elemental-H flux continuity at the matching surface, iterated to
agreement) is what makes this non-circular — see Phase E.

## 3. Work structure

The first draft's WP0-WP4 let later stages start before their inputs and
metrics existed.  This ordering (following the review) fixes that.  Phases
A-B need no new physics and de-risk everything after them.

### Phase A — freeze the target and the reference

- Re-verify the paper values against the published PDFs (the numbers in
  section 2 are quoted, not yet independently re-extracted for this plan).
- Build a **reproducible p-winds oracle**: a frozen local driver
  (`../p-winds/` clone) with pinned revision, radius grid, boundary ion
  fractions, stellar spectrum, and instrument convolution — the three
  scalars (Mdot, T, H:He) are not sufficient to reproduce a profile.
- Define the observational line metrics (windowed red depth, blue depth,
  blue/red ratio, astrophysical FWHM) exactly as the paper measured them.
- Construct the stellar spectrum: scaled GJ 1132 + GJ 699 with the
  XMM-Newton X-ray constraint, covering **4.8 eV through the X-ray band**
  (baseline row 7), with documented stitching, normalization at the orbit,
  and sampling across the H I, He I, He II, He 2^3S edges; keep the resolved
  file in the run record; sensitivity test on the 4.8-13.6 eV band.

**Acceptance:** the frozen p-winds driver reproduces the published line
profile and metrics within stated tolerances; the SED file passes the
coverage checks in the setup report.

### Phase B — one resolved configuration for wind and transit

- Have EXHALE write a structured resolved-configuration product (the values
  it actually used after `base.inp`), and make `EXHALE_transit.py` consume
  the resolved base radius/temperature from it instead of re-parsing
  `input.inp` (baseline row 8; the WASP-52b mismatch is the standing bug
  this fixes, independent of LHS 1140b).
- Implement the astrophysical metric extractor (FWHM, blue depth, windowed
  red depth, blue/red ratio) with documented windows and continuum
  treatment; validate on synthetic profiles with known answers.

**Acceptance:** wind and transit report the same resolved radius; the
extractor recovers known synthetic metrics; a test where `base.inp` moves
the radius shows the transit following it.

### Phase C — audit the helium-rich atomic wind (revised WP0)

Not a normalization audit (the hydro scale is already H+He, baseline
row 1) but an audit of the remaining **H-referenced algebra**:

- run the atomic H/He model at `HeH = 1, 10, 100, 1000` under `-fcheck`;
- audit Newton initial guesses and the advection-correction divisions
  (`ionization_equilibrium.f90:596-632`, `post_process_adv.f90:438-491`)
  for conditioning as n_H approaches its floor; record every active floor
  and show the converged result is insensitive to reasonable floor changes;
- verify charge neutrality from the actual H/He/electron densities and
  that the channel-resolved heating/cooling sums reconstruct the energy
  source term; **report** the dominant channels rather than assuming them
  (the He channels exist — baseline row 2 — the question is closure);
- check the triplet and line-cooling photon-escape assumptions in the
  helium-rich column.

**Acceptance:** a converged `HeH = 1000` atomic run with finite state
everywhere, closed charge and energy budgets, documented floors, and
stability under resolution and floor changes.

Executed 2026-08-24 on the converged LHS 1140 b solutions rather than on a
separate tutorial-planet series; the acceptance record is
`LHS1140b/exhale/audit_summary.md`. Not covered by it: resolution
sensitivity, and the photon-escape bullet above.

### Phase D — conservative H/He separation (revised WP1)

Replace, not extend, the trace-orientation assumption in
`species_diffusion.f90`: evolve a conserved element variable with a binary
flux that enforces `sum_i rho_i w_i = 0` about a defined bulk velocity;
gravity, the ambipolar term, and (new) eddy diffusion defined consistently
with that bulk velocity.  Required tests: zero-flux diffusive equilibrium;
recovery of **both** dilute limits (He-in-H reproduces the current kernel;
H-in-He); elemental conservation in a closed column; zero net diffusive
mass flux at every face; uniform-mixture preservation; grid/timestep
convergence; a homopause (molecular + eddy) test; a moving-wind elemental
face-flux test; neutral/partial/full ionization limits.

Also in this phase: **decide the design that removes the molecular-chemistry
/ `He_diffusion` mutual exclusion** (baseline row 4) — one transport
formulation valid through the matching region, or a defined transition above
the molecular layer.  The chain of Phase E cannot execute without it.

Crossover mass and O/C/N drag remain **analytic diagnostics** on the
converged wind (baseline row 5); a multicomponent momentum treatment is out
of scope.

**Acceptance:** all listed tests pass; the He-in-H limit reproduces the
current validated behavior; the incompatibility has a decided, implemented
resolution.

**Done 2026-08-26 (milestones M1-M4).**  The design is
`binary_diffusion_design.md`; the operator is
`src/modules/functions/binary_element_diffusion.f90` and the trace kernel
`src/modules/functions/species_diffusion.f90` named above is deleted, so the
row-3 and row-5 code references of section 1 are historical.  The tests are
`src/tests/diffusion_tests.f90` (`make diffusion_tests`).  The molecular /
`He_diffusion` exclusion of baseline row 4 is resolved and the `error stop`
is gone from `input_read.f90`; the changelog record is `Update_EXHALE.md`
sections 67-69 and 73, and the route is pinned by
`backup/regression/mol_diffusion`.

### Phase E — Photochem coupling with flux closure (revised WP2)

Everything here is new code (baseline row 9).  Two parts:

1. **Adapter and contract.**  A Photochem driver (configuration generator,
   result reader, failure policy, version/config fingerprint) alongside the
   existing VULCAN route, and an extended handoff carrying what the six
   scalars cannot: matching pressure/radius, T and total density, elemental
   H and He abundances **and upward elemental fluxes**, molecular/atomic H
   partition, Kzz at the match, selected molecular abundances the EXHALE
   network needs, and the source fingerprint.  Whether this stays a scalar
   file or becomes a profile file is an interface decision to confirm with
   the user before implementation.  Enforce (not print) the same-solution
   consistency of the handoff pair.
2. **Climate and closure.**  Photochem's chemistry API takes T(P) and Kzz
   as *inputs* (`EvoAtmosphere` signature; the obsolete `evolve-climate`
   setting is rejected, H2O condensation goes through the mechanism's H2O
   particle, rainout needs a set tropopause) — so the climate/cold-trap
   solution must be produced explicitly (Photochem's climate part or an
   external radiative-convective model) and documented.  Then iterate:
   trial elemental-H flux at the Photochem upper boundary -> chemistry to
   the matching pressure -> EXHALE wind -> measured lower-boundary elemental
   fluxes -> update until they agree.

**Acceptance:** lower and upper models agree on elemental H and He fluxes
at the match; the cold-trap water abundance is documented; the converged
H:He is insensitive to the initial trial flux.  Depends on
`oxygen_chemistry_new_plan.md` phases P1 (matched Photochem comparison) and
P2 (handoff contract) — this phase is those two exercised on a real target,
and as of 2026-08-27 it has exercised them.

**Done 2026-08-27 (milestones E1-E5).**  The interface decision of part 1 was
taken in favor of a profile file, and the handoff pair is enforced rather
than printed.  The implementation design is
`phase_e_flux_closure_design.md`; the changelog record is
`Update_EXHALE.md` sections 76 (the profile and the reader), 77 (the
Photochem and VULCAN adapters), 78 (the climate step and the cold trap) and
79 (the closure driver and the LHS 1140 b result).  The three acceptance
criteria, in the order they are stated above:

1. **The two models agree on the elemental fluxes at the match** — the
   closure residual `|F_measured - Phi_trial| / |F_measured|` converges to
   0.015-0.031 across three starts, inside the 0.05 criterion and above
   neither the flux window's own spread (0.004-0.009) nor the 0.025 dex
   reproducibility floor.  Test T-E6, `Update_EXHALE.md` section 79.  The
   flux is flat to 0.44 per cent (H) and 0.46 per cent (He) against the mass
   flux's own 0.45 per cent, but only in the **steady-flux** window
   `r >= r_esc`; the overlap window the design first proposed is 0.003 R_p
   thick and the flux across it spreads by 13-35, so it is measured,
   rejected on physical grounds, and the substitution is written into
   `closure.log` and the history table rather than made silently.  That is
   test T-E5, and the rejection is the result rather than a defect.
2. **The cold-trap water abundance is documented** — `f_H2O = 1.3e-7` at a
   1.031 bar, 185.0 K tropopause, 0.13-0.33 ppm over the deep-boundary scan,
   against the ~7 ppm Cherubim et al. estimate at an assumed 0.1 bar.  The
   departure is stated and not tuned.  Gas-phase `O/H` falls 6.06e-4 ->
   4.96e-7 across the trap; C, N and He are untouched.
   `Update_EXHALE.md` section 78, milestone E3.
3. **The converged H:He is insensitive to the initial trial flux** — three
   starts spanning a factor of 10 in the trial fluxes converge to
   `He/H = 2.09235` within **5.5e-6** and `log10 Mdot = 7.500` within the
   0.005 the log prints.  Test T-E7, `Update_EXHALE.md` section 79.

Where the numbers live: the closure runs are
`LHS1140b/exhale/flux_closure/{ref,lo,hi}/` (one `k**` directory per
iteration, with `closure_history.txt` and `closure.log` per arm), the
profile and its provenance-only `base.inp` are `LHS1140b/lower_profile/`,
and the execution board row is `LHS1140b/WORKPLAN.md` step E.  The route is
pinned in the regression matrix by `backup/regression/lower_profile`
(milestone E5, test T-E10).

**What the phase does not close.**  The composition at the match is now an
output, but on this planet the output is the well-mixed value: at
`K_zz = 1e9 cm^2/s` eddy mixing holds the matching-level `He/H` to within
2e-4 of it while the wind fractionates strongly above (`He/H` 2.09 at the
base, 0.183 at 30 R_p).  The answer therefore rests on the adopted eddy
coefficient (`LHS1140b/kzz_decision.md` section 0), not on the escape flux.
And the layer's *energy* is still unconstrained — item (G), the missing
continuum IR coupling, is the other half of the same weakness and is
untouched here.

### Phase F — science runs and applicability (revised WP4)

1. Prescribed-composition comparison first (H:He = 1e-3 fixed): EXHALE vs.
   the Phase-A p-winds oracle on density/velocity/temperature/ionization/
   metastable-He profiles — decoupled from Phase E.
2. Coupled-composition run (Phase E output) — the headline result.
3. 2024 line metrics through the Phase-B extractor; XUV grid at 1-33% of
   fiducial to bracket the 2025 non-detection.
4. **Collisional-validity diagnostic (new code, mandatory):** define the
   collision model (He-He, H-He neutral, ion-neutral, Coulomb), compute
   species-resolved Knudsen numbers and coupling times across the heating/
   acceleration region, locate exobase and critical point with stated
   definitions.  If the critical region is not collisional, the hydrodynamic
   result is reported as **unvalidated** (not "overestimated" — baseline
   row 11) and a kinetic estimate bounds it.

**Acceptance:** every reported observable traces to one resolved
configuration, and the applicability statement follows from computed
collision and flow scales.

Molecular chemistry audit along the way: R16-R20, R23, the shared H-He
charge-exchange path, He 2^3S + H2 Penning, and the electron/H-nucleus
closures in the He-dominated limit (baseline row 10).

[2026-08-27: item 4 is implemented as `src/utils/collisional_validity.py`,
written up in `docs/collisional_validity.md` and recorded in
`Update_EXHALE.md` section 81.  The molecular-chemistry audit is
`docs/molecular_chemistry_audit_he_rich.md`, recorded in section 82; it
found the three-body third body M was being passed as `rho/m_H` rather than
the particle density and replaced it with `calc_ntot`.  Items 1-3, the
science runs themselves, are tracked in `LHS1140b/WORKPLAN.md` and are not
asserted here.]

## 4. Explicitly out of scope

- The leading/trailing tail asymmetry, stellar-wind interaction, true
  time dependence (3-D / kinetic; the paper treats them qualitatively too).
- Multicomponent (beyond binary H/He) momentum-coupled transport — heavy
  element drag stays analytic (baseline row 5).
- The Gyr fractionation history that produced the He-dominated envelope;
  present-day composition is a boundary condition.

## 5. Relation to the oxygen chemistry plan

Unchanged from the first draft in substance: Phase E *is* oxygen-plan
P1/P2 driven by a real science case — and since 2026-08-27 it has been run
as such: P1's Photochem arm is the production adapter
(`src/utils/photochem_to_lower_profile.py`) with VULCAN as the cross-check
arm, and P2's handoff contract is what the profile schema and its refusal
rules replace at the profile level.  The diffusion work overlaps the
P3 diffusion item — with the correction that neither plan's diffusion step
yields crossover masses without a multicomponent treatment neither plan
contains.  The genuinely new physics opened here is the He-dominated limit
(Phase C) and the flux-closure coupling (Phase E part 2).

## 6. Document ownership

- **This file** — the plan of record.
- `lhs1140b_lower_atmosphere_plan.md` — first draft, kept as a record; its
  claims corrected in section 1 should not be quoted from it.
- `lhs1140b_lower_atmosphere_plan_review.md` — the external review that
  prompted this revision; its findings were re-verified against the source
  (section 1) and all checked ones confirmed.
- `oxygen_chemistry_new_plan.md` — owns P1/P2/P3 that Phase D/E build on.
