# Oxygen chemistry: revised plan

2026-08-19. This document merges `oxygen_chemistry_options.md` (the option
enumeration) with `oxygen_chemistry_plan_code_review.md` (the implementation
review of that document) into one plan. Where the two disagreed, the
disagreement was settled by reading the source again on this date; Section 1
records each verdict with the code site. The ordering in Section 5 supersedes
Section 5 of the options document and Section 9 of the review.

Status of every item below when this was written: **proposed, nothing
implemented.** That is no longer the whole picture -- P0, P1 and P2 of the
Section 5 roadmap have since been run, and each carries its own dated result
block there; P3 has since been run and carries its own dated result block,
and P4's measurable conditions (1) and (2) are filled in. The scope decision itself (item (H)
in `TO_BE_DONE.md`) remains a user decision deferred until after the paper.

The two problems, unchanged from the options document:

- **Problem A — the neutral base.** The H2-to-H partition near 1 microbar is
  set by catalytic cycles on O, OH and H2O that EXHALE has no species for.
- **Problem B — atomic oxygen in the ionized wind.** O I/II/III are already
  solved, cooled and charge-exchanged; what is missing is the observable (an
  O I transit line) and a few completeness questions.

---

## 1. Verified baseline: where the review corrected the options document

Each row was re-checked against the source on 2026-08-19, independently of
both documents.

| Claim | Verdict | Code evidence |
|---|---|---|
| `q_H2_base` pins the base H2 abundance | **Review is right: it does not.** It is an EOS anchor only: `q_h2_base` enters `composition::h2_mixing_ratio_base` -> `h2_bound_fraction` -> `comp_ntot_bc`, i.e. the base particle count and pressure normalization. The cold-start H2 profile comes from `set_IC.f90`, which calls `q_h2_equilibrium(pbar_ic, T0)` directly and never reads `q_h2_base`; the in-cell H2 is then re-determined by the local molecular solve. | `composition.f90:150-184`, `set_IC.f90` (`qh2_ic`), `input_read.f90:978-988` |
| Oxygen charge-exchange case numbers 14/33/48 | **Review is right.** The oxygen pairs are cases 13/14 (A13/A14, `O + H+` / `O+ + H`), 30/31 (C5/C6, `O + He+` / `O+ + He`), 48/49 (D17/D18, `Fe + O+` / `Fe+ + O`). The options document's "33" is actually D2 (`C+ + Si`). Both directions of each pair are present, not just one. | `charge_exchange.f90` select block |
| Huang Table 4 has no O2+ reaction; the old Kingdon-Ferland `O2+ + H` term was removed deliberately | **Confirmed against the published table** (`references/Huang_2023_ApJ_951_123.pdf`, Table 4, pp. 25-26). The oxygen entries are exactly `O + H+`/`O+ + H` (Stancil et al. 1999), `O + He+`/`O+ + He` (Zhao et al. 2004) and `Fe + O+`/`Fe+ + O` (Rutherford & Vroom 1972); no O2+ reactant appears anywhere in the table. The published table attaches the `exp(-227/T)` barrier to `O+ + H`; the code deliberately carries it on `O + H+` (endothermic, since IP(O I) > IP(H I)), per the erratum. | Table 4 read via `pdftotext`; `charge_exchange.f90:313-315`; `HUANG2023_TABLE4_OXYGEN_ERRATUM.md` |
| O I ground-term level populations are computed but not exported | **Confirmed.** The three-level statistical equilibrium lives inside the cooling functions (`Cool_coeff.f90`, `cool_OI_ne_func` and the fine-structure block around lines 1146-1500) and nothing in `write_output.f90` exports `beta_fs`/`nbar_fs` or any level fraction. | `Cool_coeff.f90`, `write_output.f90` (no hits) |
| Transit tool ignores O I | **Confirmed.** `EXHALE_transit.py` reads exactly columns 17/23/25 (Mg II, Ca II, Na I) from `Ion_species_adv.txt`; `metal_lines` and `METAL_DOUBLETS` carry only Mg II h&k, Ca II H&K, Na I D. | `EXHALE_transit.py:816`, `:835-875` |
| `base.inp` carries six keys | **Confirmed**: `T_base`, `r_base`, `HeH_base`, `Kzz_base`, `q_H2_base`, `p_base`. | `input_read.f90:1155-1224` (`read_base_inp`) |
| Molecular chemistry excludes He diffusion | **Was true when this was written; no longer.** The `error stop` was removed on 2026-08-26 with milestone M4 of `docs/binary_diffusion_design.md` (molecular carriers in the friction, mean carrier mass and charge, mole-fraction driver). B4 (oxygen diffusion) can be run with a molecular composition. | `input_read.f90` (the check is gone), `docs/binary_diffusion_design.md` section 5 |
| The molecular solve is local steady state; `_adv` post-processing is molecule-free | **Confirmed.** `post_process_adv.f90` documents at length that molecules are excluded from the `_adv` reconstruction and that a molecular base "needs a molecular-aware post-process instead". | `post_process_adv.f90:6-17, 288, 323` |
| Photochem's Zahnle mechanism is "for rocky planets only" | **Review is right that this was too narrow.** The workspace clone carries `photochem/extensions/gasgiants.py` (`EvoAtmosphereGasGiant`), `examples/GasGiants.ipynb` and a `WASP39b/` example directory: an official gas-giant workflow generated from `zahnle_earth.yaml` with H/He/N/O/C/S. An official example is still not a validation of H/H2 at 1 microbar for a specific planet. | `photochem/extensions/`, `photochem/examples/` |
| `n_species = 37`, oxygen owned by the metal block (`im_OI..im_OIII`), total O = O I + O II + O III | **Confirmed.** | `parameters.f90:26`, `species_table.f90` |

Numbers quoted below from the earlier memos and **not** re-measured today:
the VULCAN-vs-Photochem `q_H` spread of 7.2x and its 4.0x network / 1.8x code
split (`vulcan_photochem_comparison.md`, 2026-08-09, HD 189733 b); EXHALE's
`q_H2 = 0.861` vs VULCAN's 0.75 at 1 microbar and the 0.5 %-of-gap effect of
Lyman-Werner photodissociation (`TO_BE_DONE.md` item (H)); Photochem's
15-20x runtime advantage. Photochem is not importable from any Python on this
machine (checked 2026-08-19 in the previous session).

**Superseded on 2026-08-26 by phase P1** (`vulcan_photochem_comparison.md`,
section "Phase P1"): the 7.2x is unchanged as a total but splits into
domain 1.07 / network 3.95 / code 1.70, not 4.0 / 1.8; Photochem *does* import,
from the dedicated conda environment built for the 2026-08-09
comparison — the 2026-08-19 check had used the default interpreter. The 7.2x is
also specific to HD 189733 b's 864 K base: on HD 209458 b's 2331 K base the
same three arms agree to 1.02x (network) and 1.08x (code).

---

## 2. Options, restated with the review's corrections folded in

### Problem A — the base partition

**A0 — external photochemical EOS anchor (the status quo, renamed).**
What it actually does: a photochemical H2 mixing ratio sets the base particle
count and pressure normalization. What it does not do: pin H2 at the boundary,
transport H2, or seed the H2 profile. Defensible in a paper only if stated in
these terms; the previous description ("pin `q_H2_base` from the handoff")
overstated it. A small, separable improvement: use `q_H2_base` also as the
base-cell H2 seed in `set_IC`, documented explicitly as an initial guess, not
a constraint.

**A1 — widen the handoff, in three layers with separate semantics.**
- **A1a, provenance:** record which code, which reaction network (name and
  hash), and which profile produced the `base.inp` values, as machine-readable
  keys echoed by the setup report and the output header. Today this provenance
  exists only as a comment. Given the measured 7.2x network sensitivity, an
  undocumented network is an undocumented factor of several in `q_H`.
- **A1b, elemental handoff:** `O_H_base`, `C_H_base`, `N_H_base`, `S_H_base`
  (and the metal elements Lavvas finds returning to atomic form below
  ~1e-3 bar) connected to `melem_ab`. This layer has a real consumer today:
  EXHALE already transports and redistributes elemental reservoirs, so a
  photochemically motivated base abundance replaces a scaled-solar assumption.
- **A1c, molecular-species handoff (OH/H2O/CO):** *diagnostic metadata only*
  until some part of the code actually consumes them. The parser currently
  ignores unknown keys; adding species keys without a consumer widens the
  interface without changing the model, and must not be presented as physics.

The Tier-3 `base.inp` sketch in `lower_atmosphere_coupling.md` remains the
ancestor of A1b; A1a is new, prompted by the review.

**A2 — in-code catalytic chemistry, as a single-owner C/H/O element-budget
extension.** Not "add two species". The review's inventory is adopted as the
definition of A2's scope:

- one total-oxygen closure `n_O,tot = n(O) + n(O+) + n(O++) + n(OH) + n(H2O)`,
  with the metal block's O I remainder redefined to exclude oxygen bound in
  molecules — never a second oxygen reservoir;
- H nuclei bookkeeping (one per OH, two per H2O) against the H budget;
- reactions that connect the water family back to atomic O. The conservation
  skeleton to be audited at reaction level before any coding:
  `H2O + photon -> OH + H`; `OH + H2 -> H2O + H`; `O + H2 <-> OH + H`;
  `OH + photon -> O + H` (include and branch only as the audit supports).
  Sources: Lavvas & Arfaux for the hot-Jupiter context, Moses et al. (2011)
  for the network, a kinetics compilation (Baulch 2005, KIDA/UMIST) for the
  `OH + H2` barrier, Leiden/PHIDRATES for the H2O photolysis cross sections.
  Nothing from memory;
- an FUV radiation treatment. The XUV grid starts at 13.6 eV and the numerical
  SED is extended downward only for He 2^3S and low-IP metals; H2O photolysis
  lives mostly below the H I edge. A2 therefore needs either an automatic
  extension of the SED range when water chemistry is active, or a dedicated
  band treatment analogous to the Lyman-Werner one, including OH/H2O columns
  and self-shielding;
- photolysis heating and reaction energetics; molecular mass, particle and
  electron contributions to the EOS; root physicality and element-simplex
  checks; IC/restart/output schema support;
- a timescale gate: the molecular solve is local steady state with no
  transport, so A2 is only meaningful where `tau_chem << tau_adv`. The
  chemical-to-advective time profile becomes a required output, not an
  afterthought.

Affected code areas (from the review, spot-checked against the tree):
`parameters.f90`/`species_table.f90` (metadata), `System_HeH_mol_metals.f90` +
the residual core (coupled H/O rows), `ionization_equilibrium.f90` (seeding,
extraction, validation), `mol_rates.f90` or a new oxygen-chemistry module,
`cross_sec.f90`/`sed_read.f90`/`util_ion_eq.f90` or a new FUV module,
`composition.f90` (EOS), `set_IC.f90`/`load_IC.f90`/`write_output.f90`
(schema), the element-diffusion operator (one element flux if diffusion is
ever combined; the review named `species_diffusion.f90`, which was deleted on
2026-08-25 -- it is now
`src/modules/functions/binary_element_diffusion.f90`), `EXHALE_transit.py` (atomic vs total O), and the regression matrix.
The existing merged molecular+metal residual is the right starting point; its
element bookkeeping is generalized, not replaced.

**A3 — reduced C/H/O photochemical network in-code.** Deprioritized, and for a
sharper reason than cost: EXHALE's chemistry is a local algebraic steady state,
while VULCAN/Photochem are kinetics-with-transport models. Growing the residual
to 15-25 species without transport produces a large local-equilibrium solve,
not a small VULCAN. Only worth revisiting if removing the external dependency
becomes a science requirement in itself.

**A4 — two-way offline iteration.** Possible in principle, but the feedback
variables do not exist yet: the current handoff is an EOS anchor in one
direction, with no reverse species-flux schema and no defined overlap region
in which both models are valid. The interface physics (what EXHALE returns:
upper-boundary species fluxes, escape flux, attenuated irradiation; what the
lower model accepts) must be specified before a driver is written. Does not
close item (H); it makes the imposed number self-consistent, not computed.

### Problem B — oxygen in the wind

**B1 — the O I 1302/1304/1306 triplet in the transit tool.** Still the item
with the highest scientific return (the Vidal-Madjar et al. 2004 detection on
HD 209458 b, a planet already in the paper set), but it is not a line-list
change. The three components start from the three fine-structure levels of the
O I ground term, and the transit tool only has total O I; applying the total
to all three lines would count the same atoms three times. Required sequence:

1. an O I ground-term level-population helper — preferably exporting the same
   three-level statistical equilibrium the cooling already solves (as a
   function of T, n_e, n_HI and the trapped radiation field), with a
   documented Boltzmann partition acceptable as a first stage;
2. component-resolved absorption with the correct lower-level densities;
3. a band-integrated comparison in the observers' band, with the ISM
   absorption and geocoronal treatment stated explicitly (the same
   line-center-vs-band trap already documented for Mg II);
4. the Ly-beta pumping question (O I 1025.76 A near H Ly-beta) checked in the
   literature before any agreement is quoted.

**B2 — redefined from an audit to a sensitivity question.** The Table 4 audit
the options document proposed is already done: the code carries both directions
of all three oxygen pairs, Table 4 contains no O2+ reaction, and the removal of
the Kingdon-Ferland `O2+ + H` term is recorded with its reasoning. The open
question is scientific, not clerical: does an O2+ reaction *outside* Table 4
(e.g. `O2+ + H -> O+ + H+`) materially change the O III profile at EXHALE's
temperatures and ion fractions, and if so, which modern source and which
detailed-balance policy. A bounding estimate or a rate-injection experiment
answers it.

**B3 — cooling.** Complete as far as measured (the CNO CHIANTI swap moved Mdot
by 0.004 %). The shared deliverable with B1 is the level-population helper: if
Ly-beta pumping is ever added, cooling and transit must consume the same
populations, which is exactly why the helper should be one exported field
rather than two private calculations.

**B4 — oxygen diffusion.** Runnable today with `He_diffusion` +
`He_metal_diffusion` in an atomic run only (the parser forbids the molecular
combination). One numerical experiment shows whether O separates from H enough
to matter for B1. After A2, oxygen in OH/H2O would have to join the same
element flux, so the present mechanism would not carry over unchanged.

---

## 3. Decisions this plan takes, and decisions it leaves open

Taken (they follow from the code, not from preference):

- Oxygen has exactly one owner. The metal block keeps the atomic stages; any
  molecular oxygen is closed against the same total. (Both documents agree;
  the review made it a design requirement rather than a warning.)
- `q_H2_base` is described everywhere as an EOS anchor. Documentation and any
  publication text that calls it a composition pin is wrong and gets fixed.
- A1c species keys without consumers are metadata, labeled as such.
- No O-bearing molecule enters the trace-metal ionization framework; the metal
  solver's assumptions (coronal-like rates, no three-body reactions) do not
  hold for molecules.
- No rate or cross section from memory; every reaction gets the R1-R23
  treatment (named source, validity range, written at the code site).

Left open, for the user (unchanged from the handoff's B list):

- **Scope of item (H) itself**: paper-first (P0-P3 below, A0 stated honestly)
  versus code-first (through P4/A2). This plan sequences the work so that the
  decision can be made late, after P1's measurements.
- Whether a VULCAN or Photochem dependency is acceptable in the production
  path (decides A1-vs-A2 emphasis; P1 produces the evidence).
- Whether the paper should carry O I 1302 at all, given the FUV systematics.

---

## 4. Physical stakes, and the literature already on the shelf

### 4.1 Where this ranks among everything open

Ranked by the project's acceptance criterion — physical correctness — across
all items in the 2026-08-19 handoff:

1. **Item (G), the missing continuum IR coupling, ranks first**, above this
   plan's subject. It is the only open item with a *measured* order-unity
   effect: a converged molecular layer radiates itself down to 190-400 K
   against T_eq ~ 1100-1400 K, and on the HD 209458 b molecular example the
   collapse reaches the wind at **Mdot -0.34 dex**. `TO_BE_DONE.md` item (G)
   states the consequence plainly: the converged thermal structure below the
   H2 -> H front is a property of the model, not a prediction. The
   `Base IR field` closure fixed only the existing line coolants; past the
   front there are no molecular coolants left to hand a field to.
2. **Item (H), this plan's subject, ranks second — as an unbounded
   uncertainty rather than a known error.** No term is wrong; the base
   composition is an assumption, and the assumption's input moves by 7.2x in
   `q_H` with the choice of reaction network, reaching the wind through
   `ntot_bc`. How much of that reaches Mdot is unmeasured — measuring it is
   phase P1.
3. **The two are one weakness seen from two sides.** Both live in the
   molecular layer between the 1 microbar base and the H2 -> H front: (G) is
   its energy not being constrained, (H) its composition. They feed back on
   each other (the equilibrium H2 fraction depends on T; the coolant
   inventory depends on composition), so closing either alone does not make
   the layer a prediction. A decision to go past P4 into A2 should be taken
   knowing that (G) is the other half of the same commitment — while noting
   that A2's radiation work is FUV (photolysis) and (G)'s is thermal IR:
   different spectral regions, different code, one region of the atmosphere.

Caveats: this ranking is by physics, not by paper impact — the current paper's
four planets run atomic-base configurations that carry neither the collapsed
layer nor the imposed partition, which is why both items could be deferred.
And (G)'s H2 CIA component is already decided out of the paper's scope by the
user; nothing here reopens that.

### 4.2 What `references/` already carries on this physics

Checked 2026-08-19 by `pdftotext` + targeted search (titles and the matched
passages were read, not the full papers):

| Reference | What it treats |
|---|---|
| Koskinen et al. 2013a (Icarus 226, 1678) | **Both problems, on HD 209458 b.** States "the dissociation of H2 is caused by dissociation of H2O" with the H2/H transition near 1 microbar — the Problem-A mechanism of section 2.2, already modeled — and includes H3+, CO, H2O and CH4 as "strong infrared coolants", the (G) physics. The paper-level original of what A2 plus a (G) fix would build; also the ancestor of the Tier design (Koskinen et al. 2022 supplies Tier 1 and the `mol_rates` table). |
| Lavvas et al. 2014 (ApJ 796, 15) | The catalytic-destruction and atomic-metal-release context cited by `lower_atmosphere_coupling.md` — the physical basis for A1b's metal abundances. Its transit opacities include H2 CIA and H2O. |
| Lavvas & Arfaux 2021 (MNRAS 502, 5643) | Middle-atmosphere thermal structure with CIA in the radiative transfer; notes CIA becomes significant at p > 1 bar — a literature anchor for keeping CIA out of the >= 1 microbar EXHALE domain. |
| Wogan et al. 2025 (PSJ 6, 256) | The Photochem code paper ("a general chemical and climate model") — the methods source for phase P1's external arm. |
| Robeling et al. 2026 (Kompot) | **The code-level analog of the (G) physics**: a 1-D self-consistent thermo-chemical upper-atmosphere model (Jupiter as an exoplanet analogue) carrying H2-H2 and H2-He CIA (Abel et al. 2011) plus H2O/CO/CO2/CH4 opacities in the radiative budget — what EXHALE cannot do below the front. |
| Johnstone et al. 2018 (A&A 617, A107) | Same code lineage; molecular IR cooling (CO2, NO) controlling thermospheric structure in terrestrial atmospheres — the same physics class as (G). |
| Miller et al. 2013 | The H3+ cooling function — already in EXHALE (`h3p_cooling.f90`). |
| Huang et al. 2023 | **Does not carry this physics**: an atomic domain with the 1 microbar base imposed — the shape EXHALE inherited, and the reason item (H) exists. |
| Salz et al. 2016 (TPCI) | No CIA or H2O/OH passages found by search; not confirmed to treat either problem. |

Implication for the roadmap: the methods for P1 (Wogan 2025), for the A2
reaction audit (Koskinen 2013a, Lavvas 2014), and for any future (G) work
(Robeling 2026, Johnstone 2018) are all on the shelf already — start from
these before searching further afield.

---

## 5. Phased roadmap

The B1-first order of the options document is replaced: the review showed B1
has a physics prerequisite (level populations) and B2 was already largely
answered, while the cheapest genuine progress is semantic and external.

**P0 — semantics and provenance (documentation + small code, no chemistry).**
Rename A0's description everywhere `q_H2_base` is discussed; add the A1a
provenance keys and echo them in the setup report; record Photochem's
gas-giant workflow accurately in `vulcan_photochem_comparison.md`. Gate: a
legacy six-key `base.inp` reproduces the current goldens byte-for-byte.

**P1 — external photochemistry on equal footing. GATE MET, 2026-08-26; exercised on a real target 2026-08-27.** Reinstall Photochem (it
does not import today) with its data package; rerun the HD 189733 b comparison
with matched network *and* matched vertical domain, adding Photochem's official
gas-giant H/He/N/O/C(/S) mechanism as a third arm; export the H2O/OH
production/loss budgets near 1 microbar, not only `q_H2`; and measure the
EXHALE-level effect on at least one *converged* planet (HD 189733 b itself does
not reach steady state in EXHALE, so it cannot be the only case). Gate: the
7.2x is decomposed into network / domain / code with the domain-truncation
artifact removed, and the reaction budget states whether the H2O/OH cycle
dominates H2 destruction at the handoff level.

*Result* (full record and run directories:
`docs/vulcan_photochem_comparison.md`, section "2026-08-26 - Phase P1"; script
`docs/p1_matched_comparison.py`):

- decomposition at matched model top, HD 189733 b, 1 microbar:
  **domain 1.07, network 3.95, code 1.70** (product 7.20). The truncation
  artifact was the smallest of the three, and the domain sensitivity is itself
  code-dependent - truncating VULCAN instead moves q_H by 1.70x the other way.
- budget: on HD 189733 b the H2O/OH cycle carries **99% of net H2 destruction**
  at 1 microbar in all three arms (`OH + H2 -> H2O + H` alone is 57-96% of the
  gross loss). On HD 209458 b, whose base sits at 2331 K, the same channels run
  1000x faster than the net but in near-exact balance, and `H2 + M -> H + H + M`
  carries 70% of the net: **the partition there is thermal, not photochemical.**
- `tau_chem(H2)/tau_adv ~ 1e-3` at the handoff level (P4 condition 2, met at
  that level only).
- EXHALE effect, HD 209458 b, four arms JFNK-converged (`info=0`):
  **at most 0.035 dex (8.3%) in Mdot** across the arms, 1.5% between the two
  codes on the same network, +6.6% for the sulfur mechanism. The same
  measurement at the 1e-4 bar handoff, where the arms differ by 0.1%, still
  spreads by 0.025 dex, so **0.035 dex is an upper bound, not a resolved
  signal**: the wind-level consequence of the code choice on this planet is at
  the level of the configuration's own JFNK-to-JFNK reproducibility. The
  equilibrium-fit reference at 1 microbar did not converge and is not quoted.
- the code-choice evidence table (installation, reproducibility, runtime,
  network flexibility, domain limits, q_H agreement, Mdot effect) is in the same
  section. **The decision itself is left to the user, as this plan requires.**

**P2 — the handoff contract (A1b). GATE MET, 2026-08-26; superseded at the profile level 2026-08-27.** Classify every
`base.inp` key as provenance / EOS boundary / elemental reservoir / initial
guess / boundary constraint; implement the elemental keys against `melem_ab`;
keep species keys diagnostic (A1c). Gate: element budgets (H, He, C, N, O, S)
close after the handoff; the legacy file still reproduces the goldens; every
key's category is in the manual.

*Result* (full record: `docs/Update_EXHALE.md` section 74):

- the five categories are stated in `read_base_inp`'s header and in
  `docs/input_schema.md` section 2c; provenance is still comments only (its
  keys are A1a/P0), the initial-guess category is empty by design, and the
  A1c species mixing ratios stay comments because nothing consumes them.
- `<El>_H_base` (ten elements) overrides `metals.inp` and reaches `melem_ab`.
  `thereis_metals`, `melem_ab` and `thereis_lowIP_metal` moved to the
  composition block after `read_base_inp` -- before that reordering, no
  handoff element could have been heard. `base.inp` keys are now
  label-matched like `input.inp` keys.
- `vulcan_to_base.py` writes C/N/O/S summed over every carrier, and the same
  carrier sum corrects `HeH_base` (0.09698 -> 0.096915 on HD 209458 b).
- budgets: `src/utils/element_budget.py` measures them from the output
  profiles against `EXHALE_resolved.out`. The nine elemental ratios close to
  round-off (<= 1.8e-12 over 504 cells) on a `mol_metals` copy whose
  `base.inp` doubles oxygen (`O_H_base 9.80e-4`) and introduces sulfur
  (`S_H_base 1.32e-5`). Hydrogen closes to 3.1e-13 **once the He 2^3S mass
  that `calc_rho` counts twice is removed** -- a pre-existing bookkeeping
  defect of size 2.0e-4 of `rho` in the outer wind, reported and left
  standing because fixing it moves every `He23S: True` golden. [fixed
  2026-08-27: He 2^3S is flagged `bsp_is_excited_level` in
  `src/modules/init/species_table.f90` and every budget sum skips it;
  `calc_rho` and `calc_ntot` no longer take the triplet density as an
  argument. `Update_EXHALE.md` section 75; the `He23S: True` goldens were
  re-snapshotted there.]
- legacy: `make check` 6/6 byte-identical, goldens untouched.
- the user manual `.tex` is deliberately untouched; the categories reach it
  at the next manual pass. **Still open as of 2026-08-27**: the manual pass
  of that date documented the *profile* route
  (`EXHALE_user_manual.tex`, "The lower atmosphere as a profile"), including
  which `base.inp` categories a profile refuses, but the five categories of
  the scalar file itself are still only in `input_schema.md` section 2c.

**P1 and P2 on a real target (2026-08-27).**  Phase E of
`lhs1140b_lower_atmosphere_plan_new.md` is these two phases driven by a
science case, and it has now been run end to end on LHS 1140 b
(`Update_EXHALE.md` sections 76-79; design
`phase_e_flux_closure_design.md`).  What that changes about the two results
above:

- **P1's code choice was taken.**  Photochem is the production chemistry
  (`src/utils/photochem_to_lower_profile.py`) and VULCAN the cross-check arm
  (`src/utils/vulcan_to_lower_profile.py`), on the ground P1 itself
  identified: only Photochem carries a climate model, so only it can be the
  route away from a prescribed `T(p)`.  P1's network-vs-code split reappears
  at the profile level unchanged — network 3.95, code 1.70 on an 864 K base,
  1.02/1.08 on a 2331 K base (test T-E9).
- **P2's scalar contract is now the special case.**  The five `base.inp`
  categories still hold and `vulcan_to_base.py` still writes them — four
  regression cases pin that file — but with a lower-atmosphere *profile* in
  use every scalar key of the EOS-boundary, elemental-reservoir and
  boundary-constraint categories is **refused**, and the pair must carry one
  matching `solution_id`.  The single-source problem P2 could only check is
  removed by construction for those keys.  Schema and refusals:
  `input_schema.md` section 2d.
- **The composition is an output.**  P2 left the elemental reservoirs an
  input that reached `melem_ab`; the closure driver
  (`src/utils/element_flux_closure.py`) iterates the elemental fluxes of the
  two models until they agree, and returns `He/H` at the match rather than
  taking it.  On LHS 1140 b the returned value is the well-mixed one, so the
  answer rests on `K_zz` and not on the escape flux.

**P3 — the atomic-oxygen observable.**
1. B2 as the sensitivity test of the omitted O2+ rate.
2. The O I ground-term level-population helper and output field (shared with
   cooling).
3. B1, the triplet forward model, validated against the published HD 209458 b
   band-integrated depth with the ISM/geocoronal treatment recorded in the
   result metadata. Metals-off runs skip the line, as the other metal lines
   already do.
4. B4, the atomic-O diffusion experiment.
Gates: the three lower-level populations sum to total O I; optically thin and
saturated columns both tested; comparison band-integrated, never line-center.

*Result, 2026-08-30* (full record: `Update_EXHALE.md` sections 101-105; run
directories `HD209458b/oi1302`, `..._cxO2p_1`, `..._cxO2p_3`,
`..._diffusion`, all JFNK-converged with the current executable). All four
sub-items were run; the gates are met; the validation against the published
depth does **not** reproduce it.

**Gates.** The three lower-level populations sum to the total O I to round-off
(`2.9e-16`), measured by the writer and re-measured by the transit tool, which
refuses the line above `1e-10`. Both column regimes were tested by scaling the
level densities in the real code path: at `1e-9` and `1e-8` the depth is linear
in the column to 1.3% (optically thin), and a further `1e4` in column buys only
12.8x in depth (saturated). Every comparison is band-integrated; the
line-center depth is printed but never compared.

1. **B2 (the omitted O2+ rate) is material.** `O2+ + H0 -> O+ + H+`, absent
   from Huang Table 4, is carried at the published Barragan et al. (2006) rate
   with a scale dial (`metals.inp: cx_O2p_H <scale>`), the reverse direction
   excluded by detailed balance (endothermic by 21.5 eV).
   Switching it on suppresses the O III fraction on HD 209458 b by **4.7
   decades at 1.05 Rp, 3.7 at 1.1, 2.7 at 1.2, 1.7 at 1.5**, converging on the
   Table-4-only answer only past ~3 Rp. Neutral O rises 3% below 1.2 Rp and 52%
   near 2-3 Rp; `log10 Mdot` moves 9.47 -> 9.48. It was a default-off dial when
   this was written; on 2026-08-30 it became the default (scale 1) and the
   metals-bearing goldens were refreshed (`Update_EXHALE.md` section 107).
2. **The O I ground-term populations are an output field**, from the same
   three-level statistical equilibrium the `[O I]` 63/145/44um cooling is built
   on -- one definition, two consumers, no second calculation.
   `output/OI_levels.txt` / `_adv.txt`, metal-bearing runs only.
3. **B1 is implemented and under-predicts the published depth by a factor of
   30.** Model 0.359% against Vidal-Madjar et al. (2004) 10.8 +/- 4.5% and
   Ben-Jaffel & Hosseini (2010) 8.54 +/- 4.4% (both continuum-subtracted, since
   the published values are `(R_abs/R_*)^2` and include the planetary disk).
   Setting the comparison up correctly was most of the work and changed the
   answer by a factor of 12: the band weight is the *stellar emission-line
   profile* (0.2 A FWHM, peaks 1:1.5:1.17, BJ10), not a flat 10 A average; the
   ISM is inside those measured peak ratios and must not be applied a second
   time; and the instrument LSF cancels from a band ratio, so convolving before
   weighting loses another factor of 5. All of it is written into
   `tpm_OI_band_depths.txt` beside the number. The shortfall is consistent with
   what the observational literature says about every model of this line -- BJ10
   measure stellar FWHMs "far larger than the thermal Doppler width in the
   HD209458b atmosphere" and their own scaled model reaches only 3.9% -- and
   Ly-beta pumping of O I 1025.76 A, which feeds this triplet's upper level, is
   recorded as an explicit omission whose size the literature brackets between
   ~10% of direct excitation (Earth dayglow) and ~20x collisional excitation
   (solar chromosphere), with no published calculation for an escaping
   exoplanet atmosphere.
4. **B4: oxygen fractionates by 5.55x, the observable moves 18%.** With
   `He_diffusion` + `He_metal_diffusion` + `K_zz = 1e9`, O/H falls from the
   solar 4.90e-4 to 8.82e-5 above 1.2 Rp (He/H only 1.25x), but the same
   operator raises `log10 Mdot` by 0.48 dex, and the denser wind refills most of
   what settling removes: the O I 1302 band depth goes 0.359 -> 0.294%. The two
   effects are not separable in this experiment and the fractionation scales
   with `K_zz`.

`make check` 7/7 byte-identical and `run_fcheck.sh` CLEAN throughout: nothing
in P3 changes the wind. The `Update_EXHALE` sections also record two false
labels found and fixed on the way (a run log that claimed the du convergence
criterion after a JFNK finish, and `adf48` on the dielectronic data class).

**P4 — A2 go/no-go.** Start in-code oxygen chemistry only if all four hold:
(1) P1's reaction budget shows the H2O/OH cycle dominating H2 destruction near
1 microbar; (2) `tau_chem << tau_adv` there, from the P1 profiles; (3) a
science requirement exists that an EOS-plus-elements handoff cannot meet;
(4) the validation cost of the O/H conservation refactor and the FUV treatment
is accepted. Gates if it proceeds: zero-oxygen limit reproduces the current
molecular runs exactly; the atomic limit (OH, H2O -> 0) reproduces the current
O I/II/III result; element budgets close in every cell; absorbed energy
matches when the FUV opacity is widened; IC write/read round-trips H and O;
the `_adv` limitation is either lifted or made loud; the regression matrix
passes with the new physics default off.

*Conditions (1) and (2) measured, 2026-08-30.* Both were re-measured from the
stored P1 solutions rather than quoted: the Photochem budgets from
`vulcan_work/pc_compare_p1/{hd189_toa1e-2,hd209_toa1e-2}/pc_*_budget.pkl`, the
VULCAN arm reconstructed from `VULCAN/output/HD189.vul` and
`vulcan_work/hd209_vulcan/output/HD209.vul` with the rate recipe of
`VULCAN/diagnose.py`, and every advection time computed from an EXHALE
`Hydro_ioniz.txt` with that run's own `mu` and `g`. Every number in the P1.4
table was reproduced; two numbers of `vulcan_photochem_comparison.md` were
corrected in the process (see that document).

**Condition (1) — does the H2O/OH cycle dominate H2 destruction near 1 microbar?
The two planets do not give the same answer.**

| planet, level | O family, gross loss | **O family, net loss** | largest single channel |
|---|---|---|---|
| HD 189733 b, 1 ubar, 864 K (4 arms) | 59.0-99.2% | **96.3-99.6%** | `OH + H2 -> H2O + H`, 57.2-91.9% of gross |
| HD 209458 b, 1 ubar, 2331 K (4 arms) | 95.3-99.9% | **20.3-41.7%** | `OH + H2 -> H2O + H`, 50.2-65.0% of gross |
| HD 209458 b, 1e-4 bar, 1830 K (4 arms) | 43.9-100.0% | **1.3-4.2%** | `OH + H2` 98%, or `S + H2` 48% in the sulfur arm |

- **HD 189733 b: condition (1) is met**, in all four arms and independently of
  the network.
- **HD 209458 b: condition (1) is not met.** The oxygen channels are the fastest
  ones there, but they run in near-exact balance in both directions, and the net
  is carried by thermal dissociation `H2 + M -> H + H + M` (58-96% of the net,
  the higher figure at the deeper level). **At 1e-4 bar — where the production
  configuration actually hands off — the oxygen cycle carries 2.6% of the net
  (arm C').** That is the sharpest failure of condition (1) anywhere in the
  measurement, and it is at the level the code uses.

**Condition (2) — is `tau_chem << tau_adv` there? Again the two planets differ,
and on HD 189733 b the verdict depends on which advection time is meant.**

`tau_chem(H2) = n_H2/(net H2 loss)` from the same solutions;
`tau_adv = H/|v|` with `H = p/(rho g)` at the EXHALE base cell, and `r/|v|` given
alongside because the two differ by two decades in this flow.

| planet | level | tau_chem [s] | EXHALE run | tau_adv = H/v [s] | ratio | ratio with r/v |
|---|---|---|---|---|---|---|
| HD 209458 b | 1 ubar | 1.81e5 - 2.39e5 | `pc_compare_p1/exhale_hd209/p1e-6_seed/vulcan` | 2.132e8 | **8.5e-4 - 1.2e-3** | 1.0e-5 - 1.3e-5 |
| HD 209458 b | 1e-4 bar | 1.33e6 - 2.44e6 | `.../p1e-4_seed/vulcan` | 1.965e8 | **6.8e-3 - 1.2e-2** | 4.8e-5 - 8.9e-5 |
| HD 189733 b | 1 ubar | 2.99e5 - 2.51e6 | `HD189733b/output` (JFNK info=0, 2026-08-19) | 1.498e6 | **0.20 - 1.67** | 4.0e-4 - 3.3e-3 |

Compared at the *same pressure* rather than at each model's own handoff level,
HD 209458 b gives `tau_chem/tau_adv = 1.8e-3` at the EXHALE base pressure
(22.6 dyn cm^-2) and HD 189733 b gives **16.4** at its base pressure
(7.8 dyn cm^-2).

- **HD 209458 b: condition (2) is met**, and not only at the base cell. Over the
  whole pressure range where the two grids overlap (22.6 down to 1e-2 dyn cm^-2,
  `r/R_p = 1.000` to `1.087`) the pressure-matched ratio stays at or below 0.1,
  mostly 1e-3 to 3e-2; the single exception is one breathing-base inflow cell at
  `r/R_p = 1.0004`. Further out, `tau_adv` falls below the measured `tau_chem`
  beyond `r/R_p ~ 1.17` (`H/v`) or `~1.65` (`r/v`) — i.e. a local-equilibrium
  oxygen chemistry would stop being defensible in the launch region itself, well
  inside the transit-relevant radii.
- **HD 189733 b: condition (2) is not met on `H/v`** (16.4 pressure-matched, and
  0.20-1.67 when the 1 microbar `tau_chem` is set against the base-cell
  `tau_adv`); it is met on `r/v`. The verdict on this planet therefore rests on
  which advection time the local-equilibrium assumption has to beat, which is a
  physics choice this measurement cannot make for us.

**Caveat that pushes both verdicts the same way.** Every `tau_chem` above is
evaluated on the *chemistry model's own* T(p), which is hotter than EXHALE's at
the same pressure (HD 209458 b: 2021 K vs 1429 K at 21.9 dyn cm^-2;
HD 189733 b: 883 K vs 566 K at 7.8 dyn cm^-2). In-code chemistry would run on
EXHALE's temperature, so the true `tau_chem` is *longer* than quoted: HD 189733 b
fails condition (2) by more, and HD 209458 b's factor of 1e-3 margin shrinks.
Measuring how much needs the chemistry re-run on EXHALE's T(p), which was not
done.

**What conditions (1) and (2) together say.** The two conditions are met on
*different planets*: HD 189733 b has the photochemistry (1) but not the
timescale separation (2); HD 209458 b has the timescale separation (2) but its
base partition is thermal, not photochemical, so (1) fails — and fails hardest
at the 1e-4 bar level the production configuration hands off at. **No planet
measured here meets both.** A2 built to fix HD 189733 b's partition would be a
local-equilibrium solve in a regime where the chemistry is not fast compared to
the flow; A2 built for HD 209458 b would be solving for a cycle that carries a
few percent of the net.

**Conditions (3) and (4) are the user's, and this is what they are deciding
between.** What is at stake: A2 buys a base H2/H partition that is *computed*
rather than imposed, which is the only thing that turns item (H) from an
unbounded uncertainty into a modeling result, and it is the only route to a
lower atmosphere that does not carry an external code in the production path.
What it costs is not the two species: it is one total-oxygen closure spanning
the metal block and the molecular block, H-nuclei bookkeeping against the same
budget, an FUV band treatment with OH/H2O columns and self-shielding (the XUV
grid starts at 13.6 eV and water photolysis lives below it), photolysis
energetics in the EOS and the energy equation, and the seven acceptance gates
listed above — against a measured wind-level payoff on the one converged planet
of **at most 0.035 dex in Mdot**, itself an upper bound at the level of the
configuration's own JFNK-to-JFNK reproducibility (P1.6). Item (G), the missing
molecular IR coupling, sits in the same layer with a *measured* -0.34 dex, and
section 4.1 ranks it first for that reason. The P1 decision already taken —
Photochem as the production chemistry, with a climate model — means the
external-dependency argument for A2 is weaker than when this plan was written.

If P4 says no-go, the terminal state is A0/A1 with honest semantics plus the
P3 observable — which is a publishable position, provided the network
provenance of every handoff is recorded (P0).

*P4 decision, 2026-08-30 (user).* **A2 proceeds as a runtime option, default
off.** The four conditions resolve as: (1) and (2) measured above — no planet
meets both, and where the chemistry matters (cool bases) `tau_chem > tau_adv`,
so any in-code network valid there must include vertical transport; (3) is met
by **self-containedness**: EXHALE able to stand up a molecular base without
the external Photochem stack, the same code-uniformity motivation as the
Wind-AE port; (4) accepted with the build. Two conditions attach to the
decision: (a) the scope is *local kinetics plus vertical diffusive transport*
from the start — a transport-free network would be wrong exactly where it is
needed; (b) the first gate is an A/B against the Photochem handoff on the
HD 189733 b base (`q_H2` agreement within a stated tolerance), run as a
regression so the two owners of the oxygen physics cannot drift apart
silently. The A2 design document comes first and is reviewed before any
implementation; item (H) of `TO_BE_DONE.md` tracks it.

---

## 6. Standing constraints

- Every new default stays off; old goldens reproduce until a refresh is
  deliberate and reported.
- Physical correctness is the acceptance criterion; effect size and
  byte-identity are control tools.
- Documents record what a change *does* (the `q_H2_base` lesson: the interface
  existed for nine days before its actual semantics were written down
  correctly).

## 6b. What A2 measured, and how it changes this plan (2026-08-30)

**Milestone M3 of `docs/a2_oxygen_option_design.md` is implemented, and it
contradicts one premise of section 5's roadmap: transport is not the term that
sets the base partition in EXHALE's own structure.**

The P4 condition that made vertical transport mandatory was
`tau_chem(H2)/tau_adv = 0.20-1.67` at the HD 189733 b 1-microbar level. That
ratio was `H/v` on the PHOTOCHEMICAL model's scale height and velocity at that
pressure. Measured on an EXHALE run of the same planet with the option on, the
base cell has `v = 0` exactly -- the lower boundary condition puts it there --
the flow time in the cells above it is 1e8-2e8 s against `tau_chem(H2)` of
178-552 s, and the diffusive time of the base cell is 4500 chemical times at
`K_zz = 0` (3.6 at `K_zz = 1e9`). In this code's structure the base partition
is a local quantity.

The A/B gate accordingly did not enter its band, and moved away from it:
`x_H2 = 0.769` without transport, 0.608 with it at `K_zz = 0` and 0.581 at
`K_zz = 1e9`, against 0.910 for the reference and a pass band of 0.847-0.947.
What transport does change is CO, which has no chemistry and whose mixing ratio
is therefore set by the column rather than by the local equilibrium; the oxygen
it frees raises OH and destroys more H2.

**The blocker moves to item (G).** The base of these runs sits at 1191-1403 K
against the reference's 864 K, and the run's own net H2 budget puts the thermal
channel and the oxygen cycle at +2.4 and -1.4 of the net there -- the hot-base
regime the design places outside A2's validity range. The base is hot because
the option removes the [O I] and C I/C II coolants and the code has no H2O or CO
infrared bands to replace them. Section 4.1 of this plan already ranked item (G)
above this work; the M3 measurement is the quantitative reason.

Full record: `Update_EXHALE.md` section 109, the M3 result block of
`docs/a2_oxygen_option_design.md` section 7, and item (H) of `TO_BE_DONE.md`.

## 7. Document ownership after this plan

| Document | Owns |
|---|---|
| this file | the plan of record: verified baseline, corrected options, the physical-stakes ranking and literature map (section 4), roadmap, gates |
| `oxygen_chemistry_options.md` | the original option enumeration (its Section 5 ordering superseded here) |
| `oxygen_chemistry_plan_code_review.md` | the review record (its Section 9 ordering superseded here) |
| `lower_atmosphere_coupling.md` | the tier design and the Tier-3 `base.inp` sketch |
| `base_composition_handoff_plan.md` | the implemented six-key handoff |
| `vulcan_photochem_comparison.md` | the measured code comparison and its traps |
| `HUANG2023_TABLE4_OXYGEN_ERRATUM.md` | the O/H charge-exchange barrier-direction erratum |
| `TO_BE_DONE.md` item (H) | the open blocker itself |
