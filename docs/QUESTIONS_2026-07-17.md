# Questions (2026-07-17)

Two questions raised after the 2026-07-17 rate update (commit `cb84e6e`), with
the answers as verified in the code. Item 1 is a candidate physics improvement;
item 2 is a verified non-issue, recorded here so it need not be re-checked.

## 1. He recombination radiation ionizing H I — NOT considered

**Question.** Part of the radiation emitted during He II -> He I recombination
can ionize H I. Is that channel included?

**Answer: no.** EXHALE carries no diffuse radiation field. The case-B
coefficients (`alpha_B = alpha_A - alpha_1`, Badnell minus Mao & Kaastra) only
encode the on-the-spot assumption that ground-state capture photons are
reabsorbed by the *same* species locally. Consequences:

- He I ground-capture continuum photons (>= 24.6 eV) never ionize H I. In the
  standard nebular treatment (Osterbrock & Ferland, ch. 2) a fraction
  `y = n_HI sigma_H / (n_HI sigma_H + n_HeI sigma_He)` (evaluated at ~24.6 eV)
  of those photons ionizes H instead of He.
- Cascade photons from captures to excited He I levels — the 584 A resonance
  line (21.2 eV), the 2^3S -> 1^1S line (19.8 eV), and the part of the He I
  two-photon continuum above 13.6 eV — mostly ionize H in H-dominated gas.
  None of this is modeled.
- The analogous He III -> He II channels (He II Ly-alpha at 40.8 eV, He II
  recombination continua) are likewise neglected; He III is usually a minor
  reservoir in these atmospheres.

For contrast, MoCHII transports the diffuse field explicitly by Monte Carlo
(`diffuse_mod`), so the He -> H photon coupling is automatically included
there. ATES/AIOLOS-type wind codes neglect it, as EXHALE does now.

### Osterbrock & Ferland vs. Draine (2011): same physics, different packaging

Both books treat the channel with the same on-the-spot physics (Draine cites
Osterbrock 1989 for the key numbers), but the formulations differ in a way
that matters for implementation:

- **Osterbrock & Ferland** (AGN, ch. 2) keep the decay channels explicit —
  ground-capture continuum split by the local `n_HI sigma_H : n_HeI sigma_He`
  ratio, the 584 A resonance line, the 2^1S two-photon continuum (56% of
  decays yield an H-ionizing photon), the 2^3S 19.8 eV line, and the
  density-dependent collisional 2^3S -> singlet transfer — and carry each
  through the coupled H/He transfer integrals. There is no single closed-form
  recipe to lift.
- **Draine** (2011, "Physics of the Interstellar and Intergalactic Medium")
  condenses the same bookkeeping into two parameters plus closed-form balance
  equations, which is the directly implementable version:
  - `y` (Eq. 14.16): fraction of He I *ground-capture* photons (>= 24.6 eV)
    absorbed by H rather than He,
    `y = n_HI sigma_H / (n_HI sigma_H + n_HeI sigma_He)` at 24.6 eV + kT
    (`sigma_He/sigma_H > 6` there, so `y < 0.5` once `n_HeI > 0.16 n_HI`).
    The effective He recombination becomes
    `alpha_eff(He) = alpha_B + y alpha_1s2` (Eq. 14.17) — i.e., the (1-y)
    share of ground captures is returned to He on the spot, and only that
    share may be dropped from the rate; a pure case-B `alpha_B` (what EXHALE
    uses) implicitly assumes y = 0.
  - `z` (Sec. 15.5): fraction of He *case-B* (excited-state) recombinations
    whose cascade photon ionizes H. Density-dependent through the 2^3S
    channel: `z ~= 0.96` for `n_e << n_crit(2^3S) ~= 1100 e^{1.2/T4}
    T4^{0.5} cm^-3` (the 19.8 eV line always ionizes H) dropping to
    `z ~= 0.67` for `n_e >> n_crit` (2^3S collisionally converted to the
    singlets, whose two-photon channel yields an ionizing photon only 56% of
    the time).
  - The coupled ionization balances (Eqs. 15.32-15.33) then read: the He II
    balance uses `alpha_B + y alpha_1s2` against the source `(1-y) Q1`, and
    the H II balance gains the extra ionization terms
    `y Q1 + n_HeII n_e [z alpha_B(He) + y alpha_1s2(He)]`.
  - For the He III zone (Sec. 14.3.1): He II Ly-alpha (40.8 eV, 304 A),
    the 40.8 eV two-photon continuum, and the He II Balmer continuum all
    ionize H; their production rate amounts to ~80% of the local H
    recombination rate where He is doubly ionized.

### Impact on He I 2^3S (10830 A)

The coupling never enters the 2^3S balance directly: n(2^3S) is proportional
to `alpha_3 n_HeII n_e / D`, and the triplet-capture coefficient `alpha_3 =
2.10e-13 (T/1e4)^-0.778 ~ 0.75 alpha_B` is untouched by the y/z corrections,
which redistribute only the singlet/ground channels. The effect therefore
arrives through n_HeII (and, weakly, n_e and heating).

Quantified on the converged HD 209458 b profile: `n_HeI/n_HI ~ 0.08`
throughout the wind, giving a nearly constant `y ~ 0.65-0.69` (Eq. 14.16 with
`sigma_He/sigma_H ~ 6`), while `x_HeII` rises from ~0 at the base to ~0.78 at
2-3 R_p.

- **Triplet (TR) mode — the mode that computes 10830 — is protected by an
  accidental cancellation.** The current network total is
  `alpha_1 (ground, 1.54e-13) + alpha_3 (2.10e-13) = 3.64e-13` at 1e4 K,
  treating every ground-capture photon as lost. The Draine-correct net is
  `alpha_eff = alpha_B + y alpha_1s2 ~ 2.72e-13 + 0.67 x 1.54e-13 =
  3.75e-13` — only 3-4% away, because the overcounted ground channel
  (`+(1-y) x 1.54 = +0.51`) and the missing singlet-excited channel
  (`-0.25 alpha_B = -0.68`) nearly cancel at `y ~ 0.67`. The n_HeII shift is
  `~ (1 - x_HeII) x (3-4%)`, so **n(2^3S) moves by <~ 1-4%** in the
  10830-forming region — below the stellar-EUV and Penning-rate
  uncertainties. Extra 2^3S photoionization by the diffuse 19.8-21.2 eV
  photons is negligible next to the stellar 4.8-13.6 eV flux.
- **Atomic (case-B) mode is the weaker point.** Pure `alpha_B` corresponds to
  `y = 0`, underestimating the net He recombination by the factor
  `1 + 0.57 y ~ 1.38` at `y ~ 0.67` (1e4 K) — up to ~40% too little He I,
  and 10-30% too much n_HeII in partially ionized layers. This mode does not
  compute 2^3S, so 10830 is unaffected, but the He ionization structure and
  He-related cooling carry the error.
- The physical key to the cancellation: photons below 24.6 eV (584 A at
  21.2 eV, the >13.6 eV part of the two-photon continuum, the 19.8 eV line)
  **cannot reionize ground-state He I** — they can only ionize H (or 2^3S /
  metals). Only the >= 24.6 eV ground-capture continuum can return to He,
  which is exactly why `alpha_eff = alpha_B + y alpha_1s2` touches only the
  ground channel.

**If implemented in EXHALE**, Draine's parametrization is the recipe to use:
add `n_HeII n_e [z alpha_B(He) + y alpha_1s2(He)]` (with the local `y` from
Eq. 14.16 and `z(n_e)` from the 2^3S critical density) as an extra H I
ionization rate, with the corresponding photoelectron heating (photon energies
19.8-21.2 eV, so ~6-8 eV heat each), and note that the He recombination rate
itself should then be `alpha_B + y alpha_1s2` rather than pure case B. A
diffuse-field transfer would be the full solution but is out of proportion
for a 1-D wind code.

### Implemented 2026-07-17 as `use_he_rec_coupling` (default off)

Implemented exactly this recipe as the opt-in input key `He_rec_coupling`
(flag `use_he_rec_coupling`), default off = pure case B, bit-identical to the
legacy path. Subroutine `he_rec_coupling` (`util_ion_eq.f90`), called from
`ionization_equilibrium.f90` and `post_process_adv.f90`; ground coefficient
`alpha1_HeII_mao` added to `Cool_coeff.f90`. The atomic mode uses the closed
`alpha_eff = alpha_B + y alpha_1` with the `z alpha_B + y alpha_1` extra H rate;
the 2^3S (TR) mode sums the H-ionizing photons at the real channel rates instead
of using `z` (a 2^3S destroyed by photoionization or Penning emits no 19.8 eV
photon), and its 1^1S channel coefficient gains the singlet-excited capture
channel `0.25 alpha_B` that the current network omitted. The cross-section ratio
evaluated from the code is `R = sigma_He/sigma_H(24.6 eV) = 6.004`, giving
`y ~ 0.65-0.69` across the HD 209458 b wind. See update-log section 37 for the
full write-up. Documented in `docs/input_schema.md` (K14d),
`docs/Update_EXHALE.{md,tex}`, and `docs/EXHALE_user_manual.tex`.

Example results on HD 209458 b (converged profiles loaded as IC, then relaxed
15000 steps on/off from the same IC; only the on/off difference is meaningful,
since the IC predates the 2026-07-17 rate update). Atomic mode: He I x1.04 near
the base rising to a x1.25-1.32 plateau, the He II minimum down to x0.72-0.74,
temperature within +1-3%, and log Mdot 9.58 -> 9.65 (+17%). He 2^3S mode: He I
x1.01-1.16, He II minimum x0.885, temperature within +0.3%, log Mdot 9.57 ->
9.61 (+10%). The 2^3S density drops ~11% at the base but stays within ~1-2% of
the off case through the 10830-forming region (1.5-2 R_p), so He I 10830 is
essentially unchanged. **This base drop is larger than the -4% pre-estimate
above.** The pre-estimate treated 2^3S as protected by the alpha_1/alpha_3
cancellation acting only through n_HeII; in the actual TR-mode implementation
the 1^1S channel coefficient itself changes (net ground `y alpha_1` replaces the
full `alpha_1`, and the singlet-excited `0.25 alpha_B` channel is added), and
that coefficient shift feeds directly into the coupled n_HeII the 2^3S balance
divides by -- a more direct path than the cancellation argument accounted for.
The 10830-forming region is nonetheless unaffected, as anticipated.

### Update 2026-07-23: default flipped ON

As part of the Falorca & Vidotto (2026) review series (Update_EXHALE §39,
item 4; `docs/falorca2026_he3d_review.md` §3.4/§7), `He_rec_coupling` is now
**default ON**. Rationale (physical correctness): the photons are real, and
with the coupling off the TR-mode singlet recombination used alpha_1 alone —
neither case A (missing the 0.25 alpha_B singlet-excited channel) nor case B
(no local recycling of the ground-capture photons). Additional validation on
the high-gravity wasp_full case: stable (no NaN, convergence count essentially
unchanged), base n(2^3S) -10%, the 10830-forming region within ~2%, Mdot
unchanged there — consistent with the HD 209458 b table above, whose +0.04-0.07
dex Mdot shift is the reason the paper-draft planet runs need re-convergence.
`He_rec_coupling: False` restores the legacy lost-photon path (byte-identical).
A finer split of the singlet-excited capture channel (explicit
alpha_B[He(2^1S)] = 5.55e-15 (T/1e4)^-0.451 and
alpha_A-B[He(2^1P)] = 1.26e-14 (T/1e4)^-0.695 from the Allan et al. 2025
erratum's corrected Table 2, instead of the 0.25 alpha_B lump) remains an
optional refinement; the erratum's ~20% He 10830 EW sensitivity to the
2^1S handling is the reason it may matter.

## 2. Electron contribution to the mean molecular weight — included

**Question.** Does the mean-molecular-weight (particle-count) bookkeeping
include the electrons?

**Answer: yes.** The equation of state everywhere uses
`p = (n_tot + n_e) k T`:

- `src/modules/time_step/energy_semi_implicit.f90:107`
  `T_old = p / (n_tot_ad + ne_ad)` and `:172`
  `p = (n_tot_ad + ne_ad) * T_trial`;
- `src/modules/post_process/post_process_adv.f90:653`
  `p_out = (n_tot + ne)/n0*T_out`;
- the pressure-broadening opacity factor uses
  `(nh + nhe + ne) * kb_erg * T_K`
  (`src/modules/radiation/ionization_equilibrium.f90:186`).

`calc_ne` adds the metal electrons and the molecular-ion electrons under the
`eos_metals` policy, and `calc_ntot` counts every heavy species as one gas
particle, so the effective mean molecular weight `mu = rho / (n_tot + n_e)`
carries the electron contribution to the particle count. The electron *mass*
is neglected in `rho`, which is standard and safe.
