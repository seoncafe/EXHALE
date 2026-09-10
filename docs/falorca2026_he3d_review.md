# Review: Falorca & Vidotto (2026) 3D He model, what EXHALE should adopt

**Paper:** A. Falorca & A. A. Vidotto, *A Self-Consistent 3D Hydrodynamic Model for
Helium Transit Signatures in Evaporating Hot Jupiters*, MNRAS accepted
(arXiv:2607.18193, `references/2607.18193v1.pdf`). BATS-R-US, 10 hydrodynamic
variables (6 species: H0, H+, He(1^1S), He(2^3S), He+, He++), 4-channel
radiation (mUV 7 eV / sEUV 20 eV / hEUV 40 eV / X-ray 248 eV), stellar-wind
interaction, synthetic He I 10830 transits. The population scheme closely
follows Allan et al. (2024, MNRAS 527, 4657).

> **Reading note.** Everything in §§1-6 describes EXHALE as it stood when this
> review was written. All six gaps it identifies were implemented on
> 2026-07-23; §7 is the record. Where the body says a term is missing or a flag
> defaults off: §3.1 (Penning products), §3.2/§3.3 (triplet collisional
> ionization and cooling), §3.4 (`use_he_rec_coupling`), §3.5 (He<->H charge
> exchange, and its "metals-only" scope), §4 (He(2^3S)+H2 Penning): the code
> now says otherwise. Read §7 first.

**Companion references pulled for this review** (all in `references/`):
Allan_2024_MNRAS_527_4657.pdf (population scheme source; the file has the
published erratum appended: see §2.1),
Biassoni_2024_A+A_682_A115.pdf (ATES-family He I(2^3S) extension; sigma fits),
GarciaMunoz_2025_A+A_698_A199.pdf (revised He network, photoelectrons, H2),
Oklopcic_2018_ApJL_855_L11.pdf (original triplet network),
Schulik_Owen_2025_MNRAS_542_927.pdf (arXiv accepted version; OUP blocked
automated download).

**Verdict (summary).** The paper's headline content (stellar-wind confinement,
comet-tail morphology, pre/post-transit asymmetry) is intrinsically 3D and not
adoptable in a 1D code. In the microphysics, EXHALE is already equal or ahead
on most channels (Taylor 2025 Penning rate, VFKY96-fit triplet cross section
integrated over the full SED, Voronov collisional ionization, Badnell
recombination, SvS85 secondary ionization, advected triplet in the `_adv`
post-process). However, the line-by-line comparison against their Tables A1/A2
exposed **four concrete physics gaps on the EXHALE side** (§3) and one
molecular-mode gap grounded in García Muñoz (2025) (§4). All are small,
in-place fixes; none require adopting the paper's infrastructure.

---

## 1. What the paper does (relevant scope)

- Solves the 6-species population network *coupled to* the 3D hydrodynamics
  (Eq. 5), with photoionization (4 bins + a weighting factor ω_sp,λ for each
  species, their Eq. 14), recombination, collisional ionization, collisional
  (de-)excitation, charge exchange, radiative decay (Table A1, largely
  Black 1981 / Cen 1992 / Benjamin 1999 / Bray 2000 / Koskinen 2013 rates).
- 12 cooling channels (Table A2). Notably they re-express the Black (1981)
  He-triplet cooling terms with the **explicitly computed n_He(2^3S)** and warn
  that Black's implicit steady-state triplet density (∝ T^-0.6687 n_e n_He+)
  "would have been incorrect in our case."
- He+ recombination photons (24.59 eV) act locally with rate ω·α_(A-B),
  added to the photoionization of H0, He(1^1S), He(2^3S) plus their excess-energy
  heating (following Allan et al. 2024).
- Findings: stronger stellar winds reduce the He(2^3S) volume and EW (S0.1 →
  S200: EW 26.7 → 7.1 mÅ); a young-star XUV spectrum (S40F) gives a 3.3x
  deeper transit than the old-star one at the same wind strength; advection
  matters for the triplet (their Appendix D: neglecting it over-estimates
  n_23S inside ~2 Rp and under-estimates it beyond).

## 2. Where EXHALE already matches or exceeds the paper

| Channel | Falorca & Vidotto 2026 | EXHALE | Status |
|---|---|---|---|
| α_B[He+ → 2^3S] | 2.10e-13 (T/1e4)^-0.778 | identical (`rec_HeII_23S`) | same |
| q13, q31a, q31b (e- (de-)excitation) | 2.10e-8 sqrt(13.6/kT) exp(-ΔE/kT) Υ(T), Caldiroli GitHub fits | identical forms and Υ fits (`Cool_coeff.f90`) | same (shared ATES ancestry) |
| A31 | 1.272e-4 s^-1 (Drake 1971) | identical | same |
| He(2^3S)+H loss | Ξ = 5.0e-10 n_H0 (de-excitation only, Roberge & Dalgarno-era) | Penning ionization, Taylor et al. (2025) two-branch T-dependent fit | **EXHALE ahead** (rate); but see §3.1 (products) |
| σ_ph[He(2^3S)] | Biassoni (2024) piecewise power law, one 7 eV mUV bin | VFKY96 two-wing fit, integrated over the full SED from 4.8 eV with 20 sub-13.6 eV grid points | **EXHALE ahead** |
| Collisional ionization rates (atomic mode) | Black 1981 / Cen 1992 | Voronov 1997 | EXHALE ahead (atomic mode only; see §3.2) |
| Recombination H+, He++ | Storey&Hummer / Hui&Gnedin | Badnell RR+DR − Mao&Kaastra correction | EXHALE ahead |
| Secondary ionization by fast photoelectrons | not included | SvS85 partition, `use_sec_ion` default on | EXHALE ahead |
| He++ (HeIII) | solved | solved | same |
| He+↔H charge exchange | always on (Koskinen 2013: 1.25e-15 (T/300)^0.25; 1.75e-11 (T/300)^-0.75 exp(-128000/T)) | same rates in `charge_exchange.f90` Group B (Huang 2023 Table 4) but **default off** (`cx_full=.false.`), and wired only into the metals systems | see §3.5 |
| Triplet advection | full 3D continuity | local equilibrium in the main solve; upwind continuity in `post_process_adv` (`adv_implicit_HeH_TR`), which feeds the transit module | consistent with their Appendix D conclusion |
| 4-bin radiation + ω_sp,λ weighting (their Eq. 13-14) | needed because bins are broad | not needed: EXHALE integrates σ_sp(E) F(E) e^-τ(E) over the tabulated SED, which handles species competition exactly | not adoptable (already superseded) |
| exp(-τ) normalization trick (their App. B) | 3D cost optimization | 1D integrals are cheap | not needed |

### 2.1 Allan et al. (2024) erratum: what it changes

The correction (Allan et al. 2025, MNRAS 539, 910; appended to
`Allan_2024_MNRAS_527_4657.pdf`, pages 21-23) must be consulted whenever
Allan-2024 material is used:

- **Real error:** the He(2^1S) two-photon decay rate was implemented as
  A[He(2^1S)→1^1S] = 51.3e-4 f_2^1S instead of **51.3 f_2^1S s^-1**. The fix
  lowers n_He(2^1S) by ~2 orders of magnitude and (although He(2^3S) is only
  indirectly affected) reduces the predicted He(2^3S) equivalent widths by a
  consistent **~20%** across their models. Hydrodynamic escape predictions are
  unaffected.
- **Typos only (not in the model):** α_B[H0] exponent 0.9 → −0.9; the two
  charge-exchange rate labels swapped; Ψ[H0] → Ψ[He+] and f_He++ → f_He+ in
  the last collisional-ionization row.

Relevance to EXHALE: there is **no code counterpart to the 2^1S bug**, EXHALE
does not track He(2^1S) as a state (2^3S → 2^1S/2^1P conversions are treated
as instantaneous decays, which the corrected A = 51.3 s^-1 amply justifies
against collisional rates ~1e-2 s^-1). The rate coefficients compared in this
review (α_B[2^3S], α_A−B[1^1S] = 1.54e-13, α_B*[1^1S] = 6.23e-14,
Ψ[He(2^3S)], the χ pair) are all confirmed unchanged by the erratum's
reproduced Table 2. Falorca & Vidotto (2026) postdates the correction and
omits 2^1S/2^1P states entirely, so the numbers quoted from it are unaffected.
The one practical lesson: a ~20% He 10830 EW sensitivity to the handling of
the singlet-excited channels is real, which raises the priority of validating
the §3.4 `he_rec_coupling` channel fractions (584 Å vs two-photon splits)
before defaulting it on. The erratum's corrected Table 2 also lists the
explicit singlet-excited capture coefficients
(α_B[He(2^1S)] = 5.55e-15 (T/1e4)^-0.451,
α_A−B[He(2^1P)] = 1.26e-14 (T/1e4)^-0.695), which could replace the
0.25 α_B lump in `he_rec_coupling` if a finer split is wanted.

## 3. Gaps on the EXHALE side exposed by the comparison

Ordered by likely physical relevance. Verdicts follow the physical-correctness
criterion: each item is a real process that is currently misrepresented or
absent; the impact estimates are secondary, for planning only.

### 3.1 Penning ionization products are not counted in the H balance

`ion_residual_core.f90` `tr_triplet_row` removes the triplet at the Penning
rate (−n_23S n_HI Q31), but the TR hydrogen row is
`fvec(1) = n_hi*g_hi − a_hii*n_hii*n_e`: the H+ and electron produced by
He(2^3S) + H0 → He(1^1S) + H+ + e− appear nowhere. The same event is treated
as a triplet loss but not as an H ionization. The ~6.2 eV electron energy
release is also not in the heating budget. (The He bookkeeping is consistent:
the summed He I row is unaffected because the He product is ground-state He I.)

- **Physically wrong (inconsistent ledger for one reaction).** Fix: add
  −n_23S n_HI Q31 to the H row (HI→HII source) with its analytic-Jacobian
  pieces, and optionally the 6.2 eV heating term.
- Impact estimate: Penning H-ionization per H0 is n_23S Q31 ~
  (1e-7 n_H)(2e-9) ~ 2e-16 n_H s^-1, i.e. ~1e-3 of typical photoionization
  rates: small, but the term belongs in the equation.

### 3.2 TR mode silently drops all collisional ionization

The standard rows (`heh_rows`) include collisional ionization (b_hi, b_hei,
b_heii); the TR rows (`heh_tr_rows`, "Oklopčić form") include none. Turning on
`thereis_HeITR` therefore changes the H/He ionization physics beyond adding the
triplet: e-impact ionization of H0, He(1^1S), He+ disappears, and e-impact
ionization *of* He(2^3S) (present in Falorca & Vidotto as
Ψ[He(2^3S)] ∝ exp(-55338/T), threshold 4.8 eV) was never present.

- **Physically wrong (real processes omitted, and mode-dependent physics).**
  Fix: pass the existing Voronov b coefficients into the TR rows; add a
  triplet collisional-ionization coefficient (Voronov-style or the
  Black-derived 6.41e-21 sqrt(T) exp(-55338/T) with the 4.8 eV threshold).
- Impact estimate at wind temperatures (≤1.2e4 K): H CI ≲1% of
  photoionization; He(1^1S)/He+ CI negligible (thresholds 24.6/54.4 eV);
  triplet CI ~1-2% of the q31 de-excitation loss at 1e4 K, growing steeply
  with T. Small here, but the omission is exactly the kind that becomes
  invisible-wrong if hotter regimes are ever run.

### 3.3 No triplet-channel cooling (10830 collisional excitation, triplet CI)

`eval_cool` contains no term proportional to n_He(2^3S). Falorca & Vidotto
(Table A2, following Allan 2024) include, with the explicit triplet density:

- collisional excitation cooling 2^3S → 2^3P (10830 Å escape):
  1.16e-20 sqrt(T) exp(-13179/T) n_e n_23S erg cm^-3 s^-1;
- collisional ionization cooling of 2^3S:
  6.41e-21 sqrt(T) exp(-55338/T) n_e n_23S.

The 13179 K threshold (1.14 eV) is *not* Boltzmann-suppressed at wind
temperatures, unlike Lyα (118348 K). A rough ratio at T = 1e4 K:
Λ_10830/Λ_Lyα ≈ 8.9e4 × (n_23S/n_HI); with n_23S/n_He ~ 1e-6…1e-5 and
x_HI ~ 0.05…0.3 this spans ~1e-2 to order unity, so the term may plausibly
reach the few-to-tens-of-percent level of local Lyα cooling in the
triplet-rich part of the wind. The 2^3S → 2^1S/2^1P conversions (q31a/q31b)
likewise each remove 0.8/1.4 eV of thermal energy, currently unaccounted.
Biassoni (2024) neglected triplet cooling too (their §2.1), so this is not an
ATES-heritage regression, but EXHALE explicitly tracks n_23S, so the correct
form costs nothing.

- **Physically incomplete (real coolant absent; correct explicit form is
  available).** Fix: add the two n_23S cooling terms (+ optionally the
  q31a/q31b thermal ledger) to `eval_cool`, and a `Cooling_breakdown` column
  to measure the actual share. Note the caution from Falorca & Vidotto: do
  NOT use Black's implicit steady-state form, use the computed n_23S.

### 3.4 He recombination-photon coupling is implemented but default off

In default TR mode the He+ → singlet recombination uses only
α_1 = 1.54e-13 (T/1e4)^-0.486 (ground capture, A−B), which is neither case A
(missing the 0.25 α_B capture into excited singlets, ≈6.2e-14 at 1e4 K:
the α_B[He(1^1S)] = 6.23e-14 (T/1e4)^-0.827 row of their Table A1) nor case B
(the 24.6 eV ground-capture photons are assumed to escape with no local
re-ionization). Falorca & Vidotto (and Allan 2024) keep both channels and
recycle the ground-capture photons locally. EXHALE's `he_rec_coupling`
(Draine y/z on-the-spot; `use_he_rec_coupling`) already implements a *more*
detailed version (photon energies for each channel, the 19.8 eV line, 584 Å,
two-photon fractions) and its TR branch explicitly restores the missing
0.25 α_B channel; but the flag defaults to `.false.`, so none of it acts.

- **Physically incomplete at default settings** (the singlet recombination
  coefficient is an ad hoc intermediate; He recombination radiation ionizes
  and heats nothing). The paper finds recombination heating is the *leading*
  integrated heat source in their old-star models: a 3D outer-region result
  not directly transferable, but a reminder the channel is not decorative.
- Fix path: quantify on the standard planets (He 10830 EW, Mdot, T profile)
  with `use_he_rec_coupling=.true.`, then likely make it the TR-mode default;
  document the approximation either way.

### 3.5 He+↔H charge exchange exists but is off and metals-only

Group B of `charge_exchange.f90` carries exactly the two rates the paper uses
(Koskinen 2013 via Huang 2023 Table 4), but `cx_full` defaults to `.false.`
and the assembly is called only from the metals systems
(`System_HeH_metals`, `System_HeH_TR_metals`): a metals-off TR run cannot
include He CX at all. The exothermic direction He+ + H0 → He0 + H+
(1.25e-15 (T/300)^0.25) is a real He+ loss channel where the neutral H
density is high; a crude estimate gives it ~10-20% of the He+ radiative
recombination loss near the wind base (n_H0/n_e ~ 10), i.e. it perturbs the
He+ reservoir that sources the triplet.

- **Physically real, currently unreachable in the default configuration.**
  Fix path: wire the Group B pair into the He rows of all systems (not just
  metals) and evaluate enabling it by default in TR mode; the reverse
  (endothermic, exp(-128000/T)) direction is negligible below ~3e4 K and can
  stay for completeness.

### 3.6 Minor bookkeeping notes (no action required)

- Triplet integration starts at `e_th_HeTR = 4.80` eV while the σ fit is
  nonzero from 4.78 eV: the 4.78-4.80 eV sliver is skipped. Effect ~1e-3 of
  the mUV channel; harmless, but the two constants could be unified.
- The TR-mode He I summed row keeps the Oklopčić convention (singlet+triplet
  sum), which is correct given §3.1's fix stays in the H row.

## 4. Cross-check against García Muñoz (2025): the molecular-mode gap

GM25 (references/GarciaMunoz_2025_A+A_698_A199.pdf) revises the He network for
H2-bearing atmospheres and identifies **Penning ionization of H2**,
He(2^3S) + H2 → products + e−, as the dominant triplet loss toward the lower
boundary of GJ 436 b (their Fig. 4, Table A.5; Penning:associative ≈ 0.9:0.1).
EXHALE's molecular system (`System_HeH_mol`, row 8 = `tr_triplet_row`) carries
only the atomic-H Penning term: a Tier-2 run with a molecular base loses no
triplet to H2 at exactly the depths where H2 dominates. Since `Molecular
chemistry: True` now auto-couples the molecular base, this term is the
grounded next addition for combined He-triplet + molecular runs (take the rate
from GM25 Table A.5).

Also from GM25: their Movre & Meyer-based He(2^3S)+H Penning rate agrees with
the Taylor et al. (2025) fit EXHALE uses to within ~25-65% over 5e3-1e4 K
(1.3e-9 vs 2.2e-9 at 5e3 K; 1.3e-9 vs 1.6e-9 at 1e4 K), no change needed,
but worth citing as an independent check. Their photoelectron-driven
excitation/ionization of He states changes the GJ 436 b triplet peak by ≲10%
and only in deep layers; EXHALE's SvS85 treatment covers the main (ionization)
effect, and a Monte Carlo degradation scheme is not warranted now.

## 5. Not adoptable, but citable: 3D stellar-wind effects

For the transit-spectrum interpretation (EXHALE_transit.py results and the
method paper), the paper provides quantitative 1D-caveats worth citing:

- Stellar-wind confinement alone (1/40x → 50x solar mass-loss rate, same XUV)
  reduces the He 10830 EW from 26.7 to 7.1 mÅ and removes pre-transit
  absorption; comet-tail formation gives post-transit absorption in all cases.
- At fixed wind, a young-star XUV spectrum deepens the transit 3.3x.
- Triplet-density morphology (their Fig. 5) is set by where He+ survives in
  cooled, shocked gas: physics a 1D wind cannot represent.
- Their Appendix D quantifies the advection error of local-equilibrium triplet
  populations (over-estimate ≤2 Rp, under-estimate beyond): EXHALE's `_adv`
  post-process already addresses this in 1D, and the appendix is the natural
  citation for why it matters.

## 6. Recommended actions (smallest first)

1. §3.1 Penning H+ source into the TR hydrogen row (+6.2 eV heating). Tiny
   diff, closes a ledger inconsistency.
2. §3.2 collisional ionization terms in the TR rows (reuse Voronov b's; add a
   2^3S channel).
3. §3.3 explicit-n_23S cooling terms + breakdown column; measure the 10830
   cooling share on HD 209458 b / WASP-121 b.
4. §3.4 `use_he_rec_coupling` validation runs → likely default-on in TR mode.
5. §3.5 He CX Group B available outside the metals systems; evaluate
   default-on in TR mode.
6. §4 He(2^3S)+H2 Penning term in `System_HeH_mol` row 8 (GM25 Table A.5)
   before the next Tier-2 + triplet science run.

Items 1-2 change results at the sub-percent level (expected) and need golden
re-snapshots only if adopted unconditionally; items 3-5 are potentially
percent-to-tens-of-percent on He 10830 observables and should be gated by the
usual regression + a dedicated before/after comparison.

## 7. Implementation record (2026-07-23)

All six items were implemented, validated, and adopted unconditionally on
2026-07-23 (docs/Update_EXHALE_stage1.md §39 has the full record; §38 covers the
secondary-ionization staging fix that the validation work surfaced along the
way). One-line status:

1. §3.1 Penning products: DONE. Dominant H-ionization channel in the shielded
   base (n_HII up 2-7x at r < 1.03 where g_HI is exponentially killed); wind
   and Mdot unchanged.
2. §3.2 TR collisional ionization: DONE (b terms restored; `ci_HeI23S`,
   k(1e4 K) = 3.3e-10). Also closed the pre-existing inconsistency that
   eval_cool charged CI cooling for events the TR rates never produced.
3. §3.3 triplet cooling: DONE (10830 excitation + q31a/b conversion ledger;
   q13 excluded as double-counting with the Cen He I term), and wired into the
   semi-implicit energy solver (the eval_cool terms alone were diagnostic-only).
   Measured share: WASP-121b 1-5% of total cooling (median 1.3-2.9%), ~8-11%
   of local Lyα at the 2^3S peak; HD 209458 b locally up to ~28%.
4. §3.4 `use_he_rec_coupling`: DEFAULT ON. Validated stable on wasp_full
   (base n_2^3S −10%, 10830-forming region ±2%, Mdot unchanged there); the
   HD 209458 b sensitivity is Update §37's +0.04-0.07 dex in Mdot.
5. §3.5 He<->H charge exchange: group B promoted to the default reaction set
   (`he_h_charge_exchange`, default on; `cx_full` still gates C/D), wired into
   ALL ionization systems including the advection pair, with the analytic
   Jacobian for the Newton systems. Effect confined to the neutral-H base:
   He II x0.59 at r = 1.05, −1% at r >= 1.2; Mdot unchanged.
6. §4 He(2^3S)+H2 Penning: DONE (`penning_HeI23S_H2`, GM25 Table A.5 fit
   5.3791e-12 T^0.676 exp(-695.21/T), <=0.13%). The metastable drops 39x at
   the HD209 molecular base and 430x on the hot-Uranus case exactly where
   x_H2 -> 1; upper wind unaffected; Mdot unchanged.
