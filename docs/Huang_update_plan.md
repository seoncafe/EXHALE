# Upgrade Plan: EXHALE → Huang et al. (2023) WASP-121b Model

This plan describes how to upgrade `EXHALE` to reproduce the ultrahot-Jupiter
escape model of **Huang, Koskinen, Lavvas & Fossati (2023, ApJ 951, 123)** — the
WASP-121b model with trace-metal chemistry, excited hydrogen H(n=2) / Lyα
radiative transfer, Roche-lobe-overflow (RLOF) tidal dynamics, and a velocity-
broadened NUV/optical transmission spectrum. It supersedes the earlier draft;
the physics here is taken directly from the paper (equation/figure/table numbers
refer to Huang+2023 unless noted).

> [!IMPORTANT]
> **Read this first — why the previous attempt was set aside.**
> The earlier attempt (now in `ATES/EXHALE_something_wrong/`) implemented
> *all* phases in one pass and verified only that the code **compiled** and the
> LaTeX built. It was never validated against the paper's numerical results, and
> the physics came out wrong. **The governing rule of this plan is: implement one
> phase at a time and pass a quantitative physical-validation gate (against a
> specific Case/figure/table in the paper) before starting the next.** Building
> successfully is *not* a validation.

---

## 0. What EXHALE already has (do **not** reinvent these)

The current tree is already well past a bare ATES v2.0:

- **Ionization network:** a MINPACK steady-state system for **H, He, C, N, O**
  (He/C/N/O include 2nd ionization → 9 ratio-equations), in
  `src/modules/nonlinear_system_solver/System_*.f90`, assembled in
  `radiation/ionization_equilibrium.f90`.
- **Rates:** Badnell RR+DR recombination, **Voronov (1997)** collisional
  ionization, **Kingdon & Ferland (1996)** charge transfer with H, Verner
  photoionization cross sections via the A/C/P/T opacity dispatcher
  (`radiation/opacity_models.f90`, `functions/cross_sec.f90`).
- **Cooling:** AIOLOS/Black metal-line cooling (+ optional two-level
  fine-structure via `use_2lev_cool`) in `radiation/Cool_coeff.f90`.
- **Excited H (post-processing only):** `EXHALE_transit.py` already builds the non-LTE
  H(n=2) population (Christie+2013) and Hα/Hβ/Lyα/He-10830 transmission spectra
  with a **decoupled, parameterized** `J_Lyα` (Huang+2017 Eq. 6) and a
  Balmer-continuum n=2 photoionization estimate. This is *post-processing*, not
  the in-code coupled scheme Huang+2023 uses.
- **Runtime config:** `metals.inp` / `opacity.inp` (abundances and opacity model
  selected without recompiling); incremental `Makefile` build.

So the gap to Huang+2023 is: **(a)** more metal species solved in-code, **(b)**
their specific rate sources, **(c)** in-code coupled excited-H/Lyα-RT heating &
ionization, **(d)** RLOF tidal dynamics, and **(e)** a velocity-broadened 2D
transmission spectrum.

---

## 1. Physics specification (corrected to match the paper)

> [!WARNING]
> The earlier draft had three concrete errors. Use the specification below.

### 1.1 Species set (Section 2.4) — **12 elements**
H, He, **Mg, Fe, Si**, O, C, N, **S**, **Ca, Na, K** (solar abundances,
Asplund+2009). *(The old draft omitted **S**.)*

Ionization states actually solved by Huang+2023:
- **2nd ionization** (neutral, +, ²⁺): **Mg, Fe, Si, Ca** only.
- **neutral + 1st ion** (neutral, +): H, **O, C, N, S, Na, K**.
  *(Note: Huang treats C/N/O only to 1st ionization; EXHALE currently
  carries C/N/O to 2nd. See Open Question Q1.)*

Counting ratio-unknowns (Huang's choice): H(1) + He(2) + Mg/Fe/Si/Ca(2×4=8) +
O/C/N/S/Na/K(1×6=6) = **17**. If EXHALE keeps C/N/O at 2nd ionization, it is
**20**. The "19-equation" figure in the old draft was imprecise — size the
solver from this table, not from a remembered number.

### 1.2 Photoionization & collisional ionization (Section 2.4.1)
- Cross sections: **Verner+1996** (outer shell) + **Verner & Yakovlev 1995**
  (inner shell); **Opacity Project / TOP** high-resolution tables near threshold
  for H I, He I, C I, N I, O I, Na I, Mg I, Mg II, Si I, Si II, Ca I, Ca II, and
  **H I(2s), H I(2p)**; **Fe I from Zatsarinny+2019**. EXHALE already has Verner —
  extend the dispatcher and tables to the new species.
- Collisional ionization: **Voronov (1997)** for all species (already present).
- Photoelectron heating efficiency ≈ **0.93** (50 eV, electron mixing ratio 0.1;
  Cecchi-Pestellini+2009).

### 1.3 Recombination (Section 2.4.3) — **element-specific, NOT all Badnell**
- **H:** case B, Eq. (3): `α_B = 1.00e-11 (T/300)^-1.02`; recomb. energy Eq. (4).
- **He:** Yelle (2004).
- **Mg, Si, S, Ca:** Shull & van Steenberg (1982).
- **O, C, N:** UMIST / Woodall+2007.
- **Na, K:** Verner & Ferland (1996); Landini & Fossi (1991).
- **Fe I, Fe II:** explicit radiative+dielectronic fits, Eqs. (5) and (6).

  *(The old draft said "Badnell RR+DR for Mg, Fe, Si, Ca, Na, K" as the rates to
  port. That is **not** what Huang+2023 uses — so for the **validation** phase use
  the element-specific sources above, or the comparison to Cases A–D is not
  apples-to-apples.)*

> [!IMPORTANT]
> **Recombination endgame — Huang's rates *for validation*, then Badnell *for production*.**
> The Badnell RR+DR fits (Badnell 2006; the OPEN-ADAS / Badnell tabulations) are
> considered the more accurate recombination rates, and are already what
> EXHALE uses for C/N/O. **The intended final state is therefore to replace
> all metal recombination rates (H excepted — keep case B) with Badnell RR+DR.**
> Plan accordingly:
> 1. Build the recombination code so the per-ion rate source is **swappable**
>    (a dispatcher or table keyed by ion), with Huang's §1.3 sources as one option
>    and Badnell as the other.
> 2. Do the full Huang reproduction (Phase 1 → Case D) with **Huang's** rates so
>    the validation is faithful.
> 3. **Then** switch the production default to **Badnell** and re-run the Phase-1
>    and Phase-5 gates, **documenting the shift** in ionization balance, T(r), Ṁ,
>    and transit depths relative to Huang's published numbers. The recombination
>    coefficient sets the ion/neutral ratio directly, so expect a modest but
>    non-zero change — that change is the point (Badnell is the better physics),
>    not a bug. See Phase 6.

### 1.4 Charge exchange (Section 2.4.2, **Table 4**)
**65 reactions** with H and He. Forward rates from Kingdon & Ferland (1996),
Glover & Jappsen (2007), Stancil et al., Zhao et al., Neufeld & Dalgarno (1987),
Rutherford & Vroom (1972), etc.; **reverse rates via microscopic balance**,
Eq. (1): `K(T)=k_f/k_r=exp(-ΔG/RT)` with `ln K = ln a + b ln T − c/T`
(Gordon & McBride 1999 thermodynamic fits). Port the full Table 4 verbatim into
`Cool_coeff.f90` (or a dedicated `charge_exchange.f90`). On WASP-121b the
dominant relevant exchange is **Fe + H⁺ ↔ Fe⁺ + H** near the molecule→atom layer
(Fig. 11/27), but it is a *cursory theoretical estimate* — flag it as uncertain.

### 1.5 Radiative cooling (Section 2.5) — update to CHIANTI-based rates
- Two-level approximation, Eqs. (7)–(8).
- **Mg I λ2853; Mg II λ2796, λ2804** (as Huang+2017; cf. Fig. 4).
- **Ca II λ3934, λ3968** — CHIANTI electron-collisional rates.
- **Fe II** — CHIANTI: 1105 of 4339 transitions with A or oscillator strength;
  **Van Regemorter** Eq. (10) for ~1025 permitted lines (A>1e6); Hollenbach &
  McKee (1989) fine-structure/forbidden (Fig. 6).
- **Fe I** — levels ≥ 64th depleted, Boltzmann below (Fig. 5).
- **Na I**, and **free–free for ions**, Eq. (9): `Λ_ff = 1.9095e-25 Z² T4^0.55`.
- **Lyα cooling:** replace Black (1981) with **CHIANTI collisional-excitation**
  rates (Fig. 7 — Black is ~10× too high). This alone changes the thermal
  balance; do it early.

### 1.6 Excited hydrogen H(n=2) & Lyα radiative transfer (Section 2.7) — **hardest**
- **H(2p):** Eq. (11), `n_2p = n_1s (g_2p/g_1s)(c²/2hν³) J̄_Lyα`.
- **H(2s):** rate-equilibrium Eq. (14) with electron/proton ℓ-mixing 2s↔2p,
  Eqs. (12)–(13) (Seaton 1955); recombination to 2s, Eq. (15)
  (Storey & Sochi 2015); CHIANTI 2s↔1s; two-photon `A_2s=8.26 s⁻¹`;
  H(2s) photoionization `Γ_2s=576 s⁻¹`.
- **Lyα mean intensity `J̄_Lyα`:** Huang+2017 plane-parallel **Monte Carlo RT**
  with the **Hummer-IIB partial redistribution function** (not complete
  redistribution); three source channels — stellar (incident flux ÷2 for dayside
  redistribution), electron-impact excitation, and recombination cascade (case B
  → 1 Lyα/recomb). Smooth `J̄_Lyα(r)` with a B-spline; **iterate MC ↔ hydro to
  self-consistency**.
- Feed back **Balmer-continuum photoionization of H(n=2)** into the H ionization
  balance and **photoelectric + collisional-deexcitation heating** into the
  energy equation (Figs. 10, 26, 27).
- H₂ Lyα sink (`σ_H₂≈4.0e-19 cm²` at 2500 K) is **negligible** on WASP-121b —
  do not bother (the old draft's emphasis on it was misplaced).

### 1.7 Lower/middle atmosphere & boundary (Section 2.6)
Bottom boundary at **P = 1 μbar** (all molecules dissociated above). Below it,
Huang couples a hydrostatic **photochemical model (Lavvas+2014)** for
1e-6 < P < 100 bar providing radius, T, and atomic mixing ratios at the
boundary. Transit/optical continuum radius **R_p set at P = 4 mbar**. For
WASP-121b the 1 μbar level sits at **1.46 R_p**.

### 1.8 Tidal potential / RLOF (Section 5) — **the key to matching observations**
Spherical Case A **underestimates** the observed transit depths and line widths;
RLOF is required. Add the **Roche potential** Eq. (17)
`φ(x,y,z) = −GM_p/r − GM_*/|r−a| − ½Ω²[(M_*a/(M_*+M_p) − x)² + y²]`,
solve along the **substellar streamline**, locate **L1** Eq. (18), and map the
substellar profile to the terminator (`y–z`) plane assuming uniform properties on
equipotential surfaces, scaling radii by `R_RL/R_L1 ≈ 2/3` at the RL boundary.
The spherical-equivalent mass-loss is reduced by `(2/3)² = 4/9` vs the 1D
substellar rate. Triaxial radii `R_px,R_py,R_pz` with `R_py R_pz = R_p²`
(Table 2). Energy-limited cross-check: Eq. (19) with Erkaev K(η) correction.

---

## 2. Phased implementation (each phase ends at a validation gate)

> [!CAUTION]
> Do not proceed past a gate until the listed quantity matches the paper to the
> stated tolerance. If it does not, stop and debug *that* physics — do not pile
> on the next phase (that is exactly how the previous attempt failed).

### Phase 0 — WASP-121b setup + reproduce the H/He-only and Case A references
- Add WASP-121b to `src/utils/params_table.txt` (Table 1: `R_p=1.766 R_J`,
  `M_p=1.1824 M_J`, `M_*=1.3521 M_⊙`, `R_*=1.4572 R_⊙`, `a=0.02545 au`),
  `F_XUV=1.6e6 erg cm⁻² s⁻¹`, Lyα flux `1.0e5 erg cm⁻² s⁻¹`, bottom boundary at
  1 μbar.
- Run with metals **off** (≈ Case A′, H/He) and with the current C/N/O
  (≈ Case A).
- **Gate:** temperature profile peaks ≈ 12,000 K with the 50%-ionization,
  τ=1, and sonic-point radii roughly as in **Fig. 9**, and the spherical
  mass-loss rate **Ṁ ≈ 0.052 M_p/Gyr (≈ 3.7e12 g/s)** within ~50% (**Case A**,
  Table 3). *No new physics until this baseline reproduces.*

### Phase 1 — Expand the ionization network (add Mg, Fe, Si, Ca, S, Na, K)

> [!IMPORTANT]
> **Revised 2026-06-05 — two changes to *how* Phase 1 is executed.**
> 1. **Refactor the per-species interfaces *before* adding more elements
>    (Phase 1a).** Each element added so far (C/N/O, then Mg) widened the argument
>    lists of `PH_heat_HHe`, `eval_cool`, `calc_column_dens_metals`,
>    `write_output`, the MINPACK system (`ion_system_HeH_metals`), and the driver
>    (`ionization_equilibrium.f90`) by ~3–6 named per-ion arguments, and the
>    MINPACK coefficients are passed through a hand-indexed flat `params(60)`
>    array. Threading six more elements this way is unmaintainable and
>    error-prone. Convert to a **species-metadata table + array-indexed data
>    flow** first.
> 2. **Batch the elements; do not go strictly one-by-one.** Once the refactor is
>    in, the elements that share Mg's template are added *together* in a single
>    batch (Phase 1b); only the elements needing special data or physics are
>    handled separately afterward (Phase 1c–d). (This supersedes the previous
>    "add ONE element at a time, order Mg → Fe → Si → Ca → Na → K → S" rule. Mg
>    was done singly under the old scheme — see the progress note — and is the
>    template the refactor is generalized from.)

#### Phase 1a — Interface / data-structure refactor (NO new physics)

The goal is that adding an element becomes "add rows to a table," not "thread N
new arguments through M subroutines." Target structure:

- **Species-metadata module** (new `init/species_table.f90`, populated at init):
  - Per-ion arrays of length `n_ion`: `ion_elem(i)` (element index),
    `ion_stage(i)` (0/1/2), `ion_charge(i)` (for bremsstrahlung `Z²`),
    `ion_ethr(i)` (photoionization threshold [eV]), `ion_is_phot(i)` (is this ion
    photo-ionizable, i.e. not the top stage).
  - Per-element arrays of length `n_elem`: `el_abund(e)` (Asplund+2009 relative to
    H), `el_nstage(e)` (1 or 2 ionization steps solved), `el_Z(e)` (nuclear charge
    for the Gaunt factor).
  - `sigma_tab(Nl, n_phot)`: photoionization cross sections per energy bin per
    photo-ionizable ion, filled once at init by a per-ion `sigma()` dispatcher.
    Replaces the loose globals `s_hi, s_hei, …, s_mgi, s_mgii`.
- **2D density flow.** Pass the existing `f_sp(:, :)` / a derived `n_ion(:, :)`
  array directly into `PH_heat`, `eval_cool`, `write_output`,
  `calc_column_dens`, instead of unpacking into named scalars at every call.
  Each routine **loops over ions** using the metadata table:
  - `PH_heat`: `tauE = Σᵢ sigma_tab(:,i)·N_col(:,i)`; photoheating
    `Σᵢ (1−ethr(i)/e_v)·sigma_tab(:,i)·n(:,i)`; per-ion `P(:,i)`. Returns a 2D
    photoionization-rate array, no per-ion `P_*` arguments.
  - `eval_cool`: brem `= Σᵢ charge(i)²·GF(elem)·n(:,i)`; recomb/ioniz/line-cool
    coefficients via **dispatchers keyed by ion index** returning 2D arrays.
- **Generalized MINPACK system.** `ion_system` must keep hybrd1's fixed
  `(N_eq,x,fvec,iflag,params)` signature, so feed per-cell coefficients via a
  **module-level current-cell block** (a derived type the driver sets before each
  `hybrd1` call), and assemble `fvec` by **looping over elements** (1 or 2 balance
  equations each, force-zeroing absent elements as now). *Recommended over
  keeping `params`* — it removes the magic index bookkeeping and scales to the
  full ~20-unknown system and to charge exchange without renumbering. (If
  `params` is kept instead, at least compute its offsets from the metadata, not
  by hand.)
- **Grid-floor generalization.** `set_energy_vectors.f90` currently hardwires the
  below-13.6 eV floor to Mg's 7.646 eV. Make it `min(ion_ethr)` over the *active*
  low-IP metals (with Na/K this drops to ~4.34 eV, K I). Verify the stellar SED
  input is defined down to that floor; document any extrapolation.
- **Gate (refactor validation):** with **C/N/O + Mg** active, the refactored code
  reproduces the current `mg_validation/cno_mg` run (T(r), ionization fractions,
  Ṁ) **to round-off**. Same physics, restructured plumbing — bit-for-bit-ish, no
  element added until this passes.

#### Phase 1b — Uniform-template metals added *together*: Si, Ca, Na, K, S

These five share Mg's template exactly and differ only in tabulated data, so once
Phase 1a is in they go into the metadata table in **one batch**:
- 2nd ionization (neutral/+/²⁺): **Si, Ca** (like Mg).
- 1st ionization (neutral/+): **Na, K, S**.
Per element, add only table rows: Asplund+2009 abundance; **Verner+1996** (+TOP
near threshold for Na/K) cross-section coefficients; **Voronov (1997)** ionization
fits; **Badnell** recombination via the swappable dispatcher (same deliberate
§1.3 deviation already accepted for Mg — the Phase-6 production rate applied
early). **No line cooling** (Phase 2); **no charge exchange** (Phase 1d).
- **Gate:** the now ~20-unknown MINPACK system converges without NaNs across the
  full T range; ionization structure matches **Fig. 12** — Si/Ca mostly **doubly**
  ionized in the upper thermosphere, Na/K mostly **singly** ionized; per-element
  X/H abundance conserved; the batch leaves Ṁ within a factor ~2 of Case A. Bring
  numbers + figures in a notebook (no pass/fail), as for Mg. Watch the same
  C I-type **grid-extension confound** now extended to every low-IP metal — isolate
  cleanly with batch-vs-no-metals (and per-element if a single species looks off).

#### Phase 1c — Iron (handled separately — needs non-uniform data/physics)

Fe does not fit the uniform template, so it is added after the batch:
- **Cross sections:** Fe I from **Zatsarinny+2019** (not Verner) — different source
  and format; Fe II from Verner if covered.
- **Recombination:** Fe is dielectronic-dominated and complex; use Badnell RR+DR
  via the dispatcher, falling back to Huang's explicit Fe fits Eqs (5)–(6) where
  Badnell coverage is thin.
- **Stages:** Fe / Fe⁺ / Fe²⁺ (2nd ionization).
- **Fe I level depletion** (levels ≥ 64th depleted, Boltzmann below; Fig. 5) is
  primarily a *cooling* concern (Phase 2), not the ionization balance — note it,
  don't implement yet.
- **Gate:** Fe ionization structure matches **Fig. 12** (mostly Fe²⁺ in the upper
  thermosphere); convergence holds. **Fe II line cooling — the single most
  important metal coolant in Huang — is deferred to Phase 2;** flag explicitly
  that the thermal structure will not match until then.

#### Phase 1d — Charge exchange with H (Table 4)

Charge exchange is currently implemented for **C/N/O only** (Kingdon & Ferland,
`params 40–48` in `ion_system_HeH_metals`); Mg and the Phase-1b/1c metals are added
without it. Add the remaining **Table 4** rows here as a **table-driven** step
(forward rates + microscopic-balance reverse, §1.4), most importantly
**Fe + H⁺ ↔ Fe⁺ + H** (the dominant exchange on WASP-121b near the molecule→atom
layer).
- **Gate:** in the lower thermosphere charge exchange becomes the dominant H⁺
  sink, returning H⁺ to neutral H faster than radiative recombination
  (Figs. 11/27).

> [!NOTE]
> **Progress (2026-06-05) — Mg done (under the old one-at-a-time scheme).**
> Mg I/II/III is solved in-code (`n_species=18`, 11 ratio-unknowns). Verner+1996
> cross sections, Voronov ionization, **Badnell** recombination (deliberate §1.3
> deviation — the Phase-6 production rate applied early, flagged in the log). On
> WASP-121b Mg is Mg II in the lower thermosphere and Mg III in the upper (matches
> Fig. 12 ordering); MINPACK converges. Clean isolation (Mg-only vs no-metals)
> shifts Ṁ by +0.00 dex and T_peak by −40 K — negligible, as expected for a 4e-5
> trace species with **no Mg line cooling yet** (`cool_MgI=cool_MgII=0`; Phase 2)
> and **no Mg+H⁺ charge exchange** (Phase 1d). Details + notebook:
> `Update_EXHALE.{md,tex}` §2, `mg_validation/Phase1_Mg_validation.ipynb`.

> [!NOTE]
> **Progress (2026-06-05) — Phase 1a done (interface/data-structure refactor).**
> The per-element named plumbing is now table-driven (Steps A–F): a single
> `init/species_table.f90` metadata table (per-ion `mion_*`, per-element
> `melem_*`); `write_output`, `PH_heat_HHe` + column density, and `eval_cool`
> pass 2D `nm(:,1:n_mion)` / rate arrays and loop over ions; the MINPACK system
> `ion_system_HeH_metals` keeps hybrd1's fixed signature but takes per-cell metal
> coefficients via a module-level block (`set_metal_coeffs`) and assembles `fvec`
> by looping over elements, so `params` carries only H/He + C/N/O charge transfer
> and `N_eq = 3 + 2*n_melem`; the sub-13.6 eV grid floor is now `min` over active
> low-IP metals' neutral thresholds (`set_energy_vectors.f90`, driven by a new
> `melem_ab` abundance array). **No physics change:** with C/N/O + Mg the
> refactored code reproduces `mg_validation/cno_mg` bit-for-bit (log10 Ṁ = 12.76;
> three of four output files byte-identical; `Ion_species_adv.txt` differs only at
> the one uninitialized-memory HeITR ghost-cell slot, line 1 col 7, unused in this
> run). Phase 1b should now need only metadata rows + atomic/abundance data, with
> no edits to `System_HeH_metals.f90` or `ionization_equilibrium.f90`. Charge
> exchange stays C/N/O-only (Phase 1d). Details: `Update_EXHALE.{md,tex}` §3.
> **Next: the Phase-1b batch (Si, Ca, Na, K, S).**

> [!NOTE]
> **Progress (2026-06-05) — Phase 1b done (C/N/O cross-section fix + Si/Ca/Na/K/S).**
> Two bundled changes on top of Phase 1a: **(1)** the C/N/O photoionization cross
> sections were corrected from the simpler Verner & Yakovlev 1995 form (no
> `y_0`/`y_1`) to the **full VFKY96 Eq. (1)** form (`sigma_VFKY96` in
> `cross_sec.f90`); **(2)** the five uniform-template metals **Si, Ca, Na, K, S**
> were added as one batch, completing the §1.1 set except Fe. Now `n_species=30`,
> `n_mion=24`, `n_melem=9`; the MINPACK system is `N_eq = 3 + 2*n_melem = 21` with
> a per-element `melem_top` (2nd ionization for C/O/N/Mg/Si/Ca; 1st only for
> Na/K/S, whose unused upper unknown is pinned to keep the Jacobian non-singular).
> Voronov ionization + Badnell recombination (the §1.3 deviation already accepted
> for Mg) for the new ions; **no line cooling** (`mion_iscool=.false.`, Phase 2)
> and **no charge exchange** (Phase 1d) for them. The metal solver was **renamed**
> `System_HeHCO`→`System_HeH_metals` / `set_hco_metals`→`set_metal_coeffs` /
> `hco_*`→`met_*` (the old name encoded AIOLOS's C/O-only metals). WASP-121b, three
> runs (power-law SED, Roche): correcting C/N/O moves `log10 Mdot` 12.76→**12.85**
> (+0.09 dex) and `T_peak` 9828→**12197 K** (near Huang's ~12000 K), with `ref` and
> the new C/N/O+Mg run sharing the 7.646 eV grid floor so the comparison is
> confound-free; adding Si/Ca/Na/K/S then shifts Mdot by **+0.00 dex** (passive at
> solar abundance, no line cooling). Zero NaNs; the 21-unknown system converges;
> all nine elements conserve X/H to ~1e-13 including the two-stage Na/K/S; Si/Ca
> mostly doubly and Na/K/S singly ionized in the upper thermosphere (Fig. 12
> ordering). The **thermal-balance** part of the Phase-1b gate (and Huang Fig. 12
> overlay) stays open until **Phase 2** metal line cooling. Frozen regression
> reference for Phase 1c: `phase1b_validation/.gate_ref_cno_mg_vfky96/output/`
> (round-off match on shared columns, since adding Fe will change the column
> count). Details + notebook: `Update_EXHALE.{md,tex}` §4,
> `phase1b_validation/Phase1b_validation.ipynb`.
> **Next: Phase 1c — iron (Fe I/II/III).**

> [!NOTE]
> **Progress (2026-06-05) — Phase 1c done (iron Fe I/II/III added).**
> Iron is added as the third three-stage metal, completing the §1.1 element set
> (C/N/O/Mg/Si/Ca/Na/K/S/**Fe**). Because the variable-stage solver and every metal
> loop are generic over `n_melem`, this was purely metadata rows + four rate
> functions — no solver/cooling/heating/output/post-process code changed. Now
> `n_species=33`, `n_mion=27`, `n_melem=10`, `n_mphot=17`; `N_eq = 3 + 2*n_melem = 23`.
> **Two deliberate deviations from this plan's Phase-1c text, both verified against
> primary sources first:** (1) **cross sections** use the single-shell **VFKY96**
> form for Fe I/II (Verner+1996 `phfit2.f` `PH2`; thresholds Fe I 7.902 eV NIST IP,
> Fe II 16.199 eV), **not** Zatsarinny+2019 — the multi-shell upgrade is deferred
> and flagged in `cross_sec.f90` as a high-XUV lower bound; (2) **recombination**
> uses **Huang's explicit Fe fits Eqs. (5)–(6) directly** (DR Arrhenius + RR
> power-law), **not** the Badnell dispatcher, since the Huang rates are the
> validation target for this gate (Badnell-swap deferred to Phase 6). Voronov 1997
> `cfit.dat` for collisional ionization (the `P`/`X`/`K` columns were cross-checked
> — an earlier abandoned attempt had transposed them). **No Fe line cooling**
> (`mion_iscool=.false.`, Phase 2) and **no Fe charge exchange** (Phase 1d).
> WASP-121b validation (power-law SED, Roche, 4× solar Fe = `X_Fe=1.26e-4`,
> C/N/O/Mg solar): the **regression** with Fe forced off reproduces the frozen
> Phase-1b baseline to ~**1e-13** on all shared columns (the three new Fe columns
> are identically zero, no NaNs — the zero-abundance guard makes the 23-unknown
> system return the 21-unknown answer). With Fe **on**, the **Phase-1c gate**
> (Fe ionization matches Fig. 12) is met qualitatively: **Fe II dominates the base**
> (~0.89 at r=1, peak ~0.96 near 1.02 R_p), **Fe III dominates the outflow** (0.97
> at the outer edge, crossover ~1.25 R_p), Fe I minor throughout. Iron supplies up
> to **~26%** of the free electrons at 4× solar; switching it on moves Mdot only
> +0.02 dex and T_peak +26 K (heating direction — no Fe cooling yet). X/H conserved
> to 8.3e-13 across all five elements, zero NaNs. The **thermal profile is not yet
> comparable to Huang** because Fe II is their dominant low-altitude coolant
> (Phase 2). Frozen iron-on baseline for Phase 1d/2:
> `phase1c_validation/.gate_ref_cno_mg_fe/`. Details + notebook:
> `Update_EXHALE.{md,tex}` §5, `phase1c_validation/Phase1c_validation.ipynb`.
> **Next: Phase 1d — charge exchange with H for all metals (Table 4), incl.
> Fe⁺+H⁺↔Fe²⁺+H.**

> [!NOTE]
> **Progress (2026-06-05) — Phase 1d done (charge exchange, Huang Table 4).**
> A new `radiation/charge_exchange.f90` module carries all **63 transcribed
> Table-4 rows** (23 Group A metal+H/H⁺, 2 B He+H, 6 C metal+He, 32 D
> metal+metal; verbatim from the image-verified `charge_exchange_table4.md`).
> The **default** active set is Group A (metal+H), coupling every metal's
> ionization to the H⁺ fraction; a runtime `cx_full 1` in `metals.inp` enables
> Groups B/C/D. The residual is assembled **generically** (`cx_add_to_fvec`,
> element/stage→row map), so the old hard-coded C/N/O+H charge-exchange block in
> `System_HeH_metals.f90` was **deleted** and replaced by one call; rates are set
> per cell by `cx_set_cell(T)` and the active list is built by `cx_init`. Spot
> checks at 10⁴ K match the table (Fe+H⁺=4.0e-9, Fe²⁺+H=1.26e-9, O+H⁺=2.09e-9).
> **One deliberate physics change:** Table 4 only tabulates the
> neutral↔singly-ionized CT for C/N/O, so re-sourcing to Huang **drops the
> C²⁺/N²⁺/O²⁺+H recombination terms** the old KF96 set carried — intentional, to
> match the paper, and documented. WASP-121b validation (same planet/SED/grid as
> Phase 1c; `ref` = frozen Phase-1c iron-on run): in the **ionized bulk** the net
> metal CX is ≤ **3.7e-3** of radiative recombination (median 1e-3) — negligible
> H⁺ sink, matching Huang §4.2 ("Fe–H CX < H photoionization & recombination").
> O+H⁺↔O⁺+H sits near-resonant/near-equilibrium; the base H⁺ sink is the fast
> one-sided Fe+H⁺→Fe⁺+H. Dropping the 2+ CT **raises** double-ion fractions
> (OIII 0.871→0.946, NIII 0.897→0.944, CIII 0.851→0.856) and, combined with the
> new metal+H CT, **raises T_peak 12223→12955 K and moves it inward
> 1.915→1.778 R_p**. Charge exchange adds **no cooling channel**, so metals-on
> still runs hotter than Huang (Fe II line cooling = Phase 2). Group A vs full
> Table 4: bulk T differs by ≤1.5% (max |ΔT|=556 K at r=1.008), Mdot −0.01 dex,
> outflow nearly unchanged. All three runs zero NaNs. New frozen baseline for
> **Phase 2** (the C/N/O re-sourcing shifted it off the Phase-1c snapshot):
> `phase1d_validation/.gate_ref_cno_mg_fe_cx/`. Details + notebook:
> `Update_EXHALE.{md,tex}` §6, `phase1d_validation/Phase1d_validation.ipynb`.
> **With this the ionization-network half of Phase 1 is physically complete.
> Next: Phase 2 — metal line cooling (incl. Fe II).**

### Phase 2 — Updated radiative cooling (§1.5)
- Implement Mg I/II, Ca II, Fe II, Fe I, Na I two-level/CHIANTI cooling and
  free–free in `Cool_coeff.f90`; replace Black (1981) Lyα cooling with CHIANTI.
- **Gate:** per-species cooling-rate curves reproduce **Figs. 4–7** (e.g.
  Mg II is the dominant coolant; Fe II free–free matters at T<3000 K); a metal-
  driven cooling feature appears at **1.15–1.4 R_p** (Fig. 10), though on
  WASP-121b adiabatic cooling still dominates globally.

> **Progress (2026-06-05) — Phase 2 done; gate closed (user-validated).**
> The CHIANTI-based line coolants (Mg I/II, Ca II, Na I, Fe II on top of the
> C/N/O coolants) are active in the coupled energy balance; per-species physics
> and atomic-data provenance live in `Update_EXHALE_early_phase` (Parts II–III),
> and the integrated WASP-121b gate is written up
> in `Update_EXHALE.{md,tex}` §7. The gate was closed with a new **exact
> per-channel cooling diagnostic**: `eval_cool` gained an optional `cool_chan`
> out-array (an exact split of the total `cool` into H/He recombination,
> collisional ionization, collisional excitation, bremsstrahlung, and one column
> per metal ion), and `write_cool_breakdown_eq` (utils_ion_eq, called once from
> `EXHALE_main.f90`) dumps `output/Cooling_breakdown.txt`. Consistency:
> `max |Σchannels/cool − 1| = 3.7e-16` (machine eps). **Result (WASP-121b eq):**
> in **1.15 ≲ r/Rp ≲ 1.4** the dominant radiative coolants are **Fe II (~30%)
> and Mg II (~25%)**, jointly ~55%; metals are **57–85%** of radiative cooling
> through the band. Fe II leads the inner edge (top coolant at 1.1–1.2 R_p, 53%
> at 1.1), Mg II takes over by 1.3–1.4 R_p, and H recombination only overtakes
> above ~1.5 R_p — the Fe II → Mg II handoff of Fig. 10. Mdot unchanged (12.51).
> **Deferred (not gate blockers):** Fe II IR fine-structure / density-dependent
> two-level (current Fe II is coronal, validated vs Fig. 6 within ~1.5×); Fe I
> line cooling (Fig. 5); the CHIANTI Lyα swap (intentionally a Phase-3 RT effect,
> the 0.35× suppression, not a cooling-rate change). **Next: Phase 3 — excited
> H(n=2) + Lyα RT.**

### Phase 3 — In-code excited H(n=2) + Lyα RT (§1.6)
- New `lya_rt.f90`: H(2s)/H(2p) populations, Balmer-continuum photoionization
  into the H balance, photoelectric + deexcitation heating into the energy eqn.
- **Stage it:** first drive H(n=2) with a **parameterized `J̄_Lyα`** (the
  EXHALE_transit.py / Huang+2017 Eq. 6 estimate) to get the coupling and signs right; only
  then add the **Monte Carlo RT + B-spline + outer iteration**. Keep the RT
  decoupled from the hydro sub-step (iterate between converged hydro states).
- **Gate:** photoionization of H(n=2) is a significant proton source below
  ~2 R_p (Figs. 11, 27); Balmer photoelectric heating appears in the heating
  budget (Figs. 10, 26) without destabilizing the energy solver.

### Phase 4 — RLOF / tidal potential (§1.8)
- Add the Roche potential to `grav_field.f90` and the momentum/energy source
  terms; solve the substellar streamline; compute L1; implement the
  substellar→terminator mapping and the `4/9` mass-loss conversion.
- **Gate:** with the tidal potential on (canonical params) reproduce **Case B,
  Ṁ ≈ 0.32 M_p/Gyr**; with the Case D parameter set (below) reach
  **Ṁ ≈ 1.03 M_p/Gyr** (Table 3). Outflow velocity ≈ 7× the Case A value at R_*
  (Fig. 20).

### Phase 5 — Velocity-broadened transmission spectrum (extend `EXHALE_transit.py`)
- Continuum: H⁻ (John 1988), Rayleigh-H (Lee & Kim 2004), He, H₂.
- Lines: H Balmer (Hα, Hβ, Hγ from H(2s)/H(2p)); **Mg II λ2796/λ2804;
  Ca II λ3934/λ3968; Na I λ5892/λ5898 (+ λ3303); K I λ4045/λ4048/λ7667/λ7701;
  Mg I λ2853; Ca I; Fe I (492)/Fe II (829) NIST transitions** (assume levels
  above the depletion threshold are empty, Boltzmann below).
- **Velocity broadening** (Section 4.3, Eq. 16): line-of-sight outflow
  `v(r_eff)·x/r_eff` + tidally-locked rotation `Ωy`; 2D annulus split into
  **20 angular sectors**, ~2000 radial samples; **triaxial Roche geometry**
  (R_px,R_py,R_pz). Build the transit depth via Eq. (16).
- **Gate (the headline result):** line-center transit depths approach the
  **Case D / observed** values below.

### Phase 6 (final, production) — swap recombination to Badnell RR+DR
Once Cases A–D are reproduced with Huang's recombination rates (§1.3), switch the
production default to **Badnell RR+DR** for all metals (H stays case B). Because
the recombination code is built with a swappable per-ion rate source (§1.3), this
is a configuration change, not a rewrite.
- **Gate:** re-run the Phase-1 (ionization structure) and Phase-5 (transit depth)
  checks with Badnell rates and **record the deltas** vs the Huang-rate results
  (ionization fractions, T(r), Ṁ, Mg II / Ca II / Hα / Hβ depths). Keep both rate
  sets selectable so the Huang-rate run remains reproducible for comparison.

---

## 3. Validation targets (from Table 3 / Figs. 17–24; see also `Huang_update_expected_results.md`)

| Quantity | Case A (spherical) | **Case D (preferred, RLOF)** | Observed |
| :--- | :---: | :---: | :---: |
| Ṁ | 0.052 M_p/Gyr | **1.03 M_p/Gyr** | — |
| Mg II λ2796 (4 Å bin) | 0.182 | **0.302** | 0.309 ± 0.036 |
| Ca II K λ3935 | 0.199 | **0.278** | 0.281 ± 0.009 |
| Hα λ6563 | 0.201 | **0.185** | 0.186 ± 0.003 |
| Hβ λ4861 | 0.174 | **0.135** | 0.143 ± 0.005 |
| Na D2 λ5890 | 0.152 | **0.147** | 0.147 ± 0.002 |
| χ²/N (NUV) | 2.08 | **1.20** | — |

**Case D parameter set** (relative to canonical): tidal potential **on**,
`M_p` lowered by 1σ (→1.1204 M_J), lower-atmosphere `dT = +350 K`, stellar
**XUV × 0.74**, **Lyα × 0.35**, **Fe abundance × 4** (other metals solar).
Intermediate Cases B, C bracket the path (Table 3); use them as milestones.

**Numerical gates:** the 17–20-equation MINPACK system converges without NaNs in
both hot (ionized) and cool (boundary) regimes; the semi-implicit energy solver
stays stable with the new cooling and Balmer heating; energy/mass conservation
holds across the RK loop.

---

## 4. Open questions / decisions (with recommendations)

> [!CAUTION]
> **Q1 — C/N/O ionization depth.** EXHALE carries C/N/O to 2nd ionization;
> Huang+2023 stops at 1st. *Recommend keeping EXHALE's 2nd-ion C/N/O* (more
> complete, already validated) and noting the difference, rather than reducing
> the network to match the paper exactly.

> [!CAUTION]
> **Q2 — 1 μbar boundary.** Hardcode a WASP-121b profile, or build a file-reader
> for an external photochemical/Lavvas profile? *Recommend a labelled file-reader*
> (`boundary_profile.txt` with `HI`, `MgI`, `FeII`, … rows), reusing the concept
> from the prior attempt's `load_boundary_profile`, so other planets are easy.

> [!CAUTION]
> **Q3 — Lyα RT fidelity.** Full Monte Carlo + Hummer-IIB PRF + outer iteration
> is expensive and the riskiest single piece. *Recommend the staged approach in
> Phase 3*: parameterized `J̄_Lyα` first, MC RT only once the coupling is proven.

> [!CAUTION]
> **Q4 — RLOF in 1D.** The substellar-streamline + terminator-mapping
> approximation (Section 5) is inherently multidimensional physics squeezed into
> 1D. Validate the geometry against the spherical Case A result *before* trusting
> Case B/D numbers.

---

## 5. Risk register (likely causes of the previous failure)

- **Argument-list explosion** → each element previously widened ~6 subroutine
  interfaces plus a hand-indexed `params(60)`. Do the **Phase 1a** refactor
  (species-metadata table + array-indexed flow + module-level MINPACK cell block)
  *before* adding elements, and gate it by reproducing the current C/N/O+Mg result
  to round-off — the refactor must change no physics.
- **Coupled 17–20-eq solver non-convergence** → after the refactor, add the
  uniform metals as a batch but keep good initial guesses (Saha or previous-cell
  values) and the force-zero-if-absent behavior; if a batch fails to converge,
  bisect by disabling elements via the metadata table (not by reverting code).
- **Wrong rate sources** → use the §1.3/§1.4 element-specific sources exactly;
  do not substitute Badnell for Huang's choices.
- **Lyα RT instability / wrong sign** → stage it (parameterized → MC), iterate
  between converged hydro states, B-spline-smooth `J̄_Lyα`.
- **RLOF geometry errors in the spectrum** → reproduce spherical Case A transit
  depths first; only then enable triaxial geometry + broadening.
- **"It compiles" ≠ correct** → every phase ends at a *numerical* gate tied to a
  specific Case/Figure/Table above.

---

## 6. Source-file map (where each change lands)

| Concern | File(s) |
| :--- | :--- |
| **Species metadata table (Phase 1a refactor)** | new `init/species_table.f90` (per-ion/per-element metadata, `sigma_tab`, rate dispatchers) |
| Ionization network / MINPACK | `nonlinear_system_solver/System_*.f90`, `radiation/ionization_equilibrium.f90` (generalize `params`→module-level cell block; loop over elements) |
| Per-species interfaces (de-argument) | `radiation/util_ion_eq.f90` (`PH_heat`, `eval_cool`), `functions/utilities.f90` (`calc_column_dens`), `files_IO/write_output.f90` — pass 2D `n_ion(:,:)`, loop over the metadata table |
| Recombination, Voronov, cooling, charge exchange | `radiation/Cool_coeff.f90` (consider splitting `charge_exchange.f90`) |
| Photoionization cross sections | `functions/cross_sec.f90`, `radiation/opacity_models.f90` |
| Species arrays, abundances, RLOF params, boundary | `init/parameters.f90`, `files_IO/metals_input_read.f90`, `files_IO/opacity_input_read.f90` |
| Allocation / initial & boundary conditions | `init/set_energy_vectors.f90`, `init/set_IC.f90`, `files_IO/load_IC.f90` |
| Column densities | `functions/utilities.f90` |
| Excited-H / Lyα RT | new `radiation/lya_rt.f90` |
| Roche potential / tidal source terms | `functions/grav_field.f90`, `states/Source.f90` |
| Output of new species | `files_IO/write_output.f90` |
| Transmission spectrum (broadening, geometry) | `EXHALE_transit.py` |
| Add new sources to the build | `SRC` list in `Makefile` (deps auto-generated) |
