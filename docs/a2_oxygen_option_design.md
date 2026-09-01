# A2: in-code oxygen chemistry with vertical transport (design)

**Status: design of record. Milestones M1, M2 and M3 executed.**
M1 (section 7) is complete and its result is `docs/a2_reaction_audit.md`, which
supersedes this document wherever the two differ on a rate coefficient or a
source. The corrections M1 forced have been folded into sections 2.2, 2.4, 2.5,
2.6 and 7 below and are listed in section 8 of the audit. The decision that
authorizes it is the P4 block of `docs/oxygen_chemistry_new_plan.md`, dated
2026-08-30: *A2 proceeds as a runtime option, default off*, with two conditions
attached — (a) the scope is local kinetics **plus vertical diffusive transport**
from the start, and (b) the first gate is an A/B against the Photochem handoff
on the HD 189733 b base, run as a regression. This document is the design that
decision asked for, and it is to be reviewed before any code is written.
`TO_BE_DONE.md` item (H) tracks the work.

Inputs read for this design, all of them read rather than recalled:

- `docs/oxygen_chemistry_new_plan.md` — the plan of record: the A2 scope
  inventory (section 2, "Problem A"), the physics ranking (section 4.1), the
  P1/P3/P4 result blocks and the seven acceptance gates of P4;
- `docs/vulcan_photochem_comparison.md` sections P1.4, P1.5 — the measured
  reaction budget and the network / code / domain decomposition;
- the P1 budget files themselves,
  `vulcan_work/pc_compare_p1/{hd189_toa1e-2,hd209_toa1e-2,hd189_toa6e-3}/pc_*_budget.pkl`
  and the matching `pc_*_solution.pkl`, re-measured for this design (section 2.2);
- `docs/binary_diffusion_design.md` sections 2, 3, 4, 7 — the transport
  formulation this design reuses;
- `docs/phase_e_flux_closure_design.md` sections 2.4, 3.2, 3.3 — the profile
  handoff, the single-source rule and its refusals, and the eddy coefficient;
- the source, read for this design: `src/modules/lower_atmosphere/mol_rates.f90`,
  `lyman_werner.f90`, `lower_column.f90`;
  `src/modules/nonlinear_system_solver/System_HeH_mol.f90`,
  `System_HeH_mol_metals.f90`; `src/modules/radiation/ionization_equilibrium.f90`,
  `util_ion_eq.f90`, `sed_read.f90`; `src/modules/functions/composition.f90`,
  `utilities.f90`, `cross_sec.f90`, `binary_element_diffusion.f90`;
  `src/modules/init/species_table.f90`, `parameters.f90`, `set_energy_vectors.f90`;
  `src/modules/files_IO/input_read.f90`, `write_output.f90`, `load_IC.f90`,
  `lower_atmosphere_profile.f90`, `write_setup_report.f90`;
  `src/modules/post_process/post_process_adv.f90`; `src/EXHALE_main.f90`.

---

## 1. Goal, acceptance, and the choices that shape everything else

### 1.1 What the option must deliver

Item (H) is the statement that EXHALE cannot compute its own base H2/H
partition: the molecular network of `mol_rates.f90` has no species to put the
catalytic cycles on, so the partition is either imposed through the handoff
(`q_H2_base`, an EOS anchor) or produced by a chemical-equilibrium fit
(`q_h2_equilibrium`, `lower_column.f90:61`) that is not photochemistry. A2 is
the option that computes it.

**Acceptance, from the P4 decision.** With the option on, EXHALE stands up a
molecular base without the external Photochem stack, and the H-nucleus
partition it produces at the HD 189733 b base agrees with the Photochem
handoff within a stated tolerance (section 6.2). With the option off, every
existing result is reproduced bit for bit.

### 1.2 The three choices that shape the rest

**C1 — species-resolved kinetics with vertical transport, not an extension of
the local equilibrium solve.** This is P4 condition (a), and the measurement
behind it is unambiguous: at the HD 189733 b base, `tau_chem(H2)/tau_adv` is
**0.20 to 1.67** on `H/v` and **16.4** when compared at the same pressure
(P4 condition (2) table). A local algebraic steady state is exactly wrong
there, and that is the one planet where the oxygen photochemistry is what sets
the partition. An A2 built as more rows of `System_HeH_mol` would inherit the
defect it exists to remove. Section 3 is therefore not an appendix to section 2;
it is half the design.

**C2 — one oxygen, one hydrogen, one owner each.** The metal block already owns
oxygen as O I / O II / O III (`species_table.f90:128`, `iel_O = 2`), and the
molecular block already owns hydrogen nuclei through `bsp_nH`
(`species_table.f90:103`). A2 adds carriers to both elements. The design
requirement, taken from the plan and not weakened here, is a single total-oxygen
closure

```
n_O,tot = n(O I) + n(O II) + n(O III) + n(OH) + n(H2O) + n(CO)
```

with the metal block's O I redefined to mean *free atomic neutral oxygen*, never
a second oxygen reservoir. Section 4 states where every one of those terms is
read and written.

**C3 — the FUV arrives as bands with named fluxes, not as an extension of the
XUV energy grid.** The photon energy grid is built from the ionization
thresholds of the species the run carries: it starts at `e_th_HI = 13.6 eV` and
is lowered only to the He 2^3S edge (4.8 eV) or to the lowest active low-IP
metal threshold, in both the power-law path
(`set_energy_vectors.f90:36-56`) and the loaded-SED path (`sed_read.f90:36-47`).
Two things follow. First, when the run's SED is a power law (`is_PL_sed`), that
power law is an XUV fit and extrapolating it into the FUV is not a spectrum.
Second, water photolysis lives mostly below the H I edge, in a wavelength range
where the absorber is a continuum, not a set of ionization edges. The existing
precedent in the code is exactly right for this: `Stellar LW flux
[erg/cm2/s]` plus a self-shielding function
(`lyman_werner.f90`, `input_read.f90:436-442`). A2 follows it. Section 2.6
states the bands.

### 1.3 What A2 is not

- It is not a photochemical model of the lower atmosphere. Photochem and VULCAN
  keep their role, and Phase E's profile route
  (`docs/phase_e_flux_closure_design.md`) keeps its. A2 is the route that lets
  EXHALE run without them, at a cost in network size that section 2.7 states.
- It is not a carbon or nitrogen or sulfur chemistry. Section 2.3 gives the
  measured reason each is excluded and the one exception (CO as a reservoir).
- It does not close item (G), the missing continuum IR coupling. Section 9.
- It does not change the wind. Every gate in section 6 is written so that a run
  without the key is byte-identical.

---

## 2. Reaction network scope

### 2.1 The network A2 extends

`src/modules/lower_atmosphere/mol_rates.f90` carries the 23 rate coefficients
of Koskinen et al. (2022), ApJ 929:52, Table 1, transcribed with their sources
at the code site (R1-R23; R21/R22 transcribed but deliberately unused, because
the H <-> He charge-exchange pair comes from the Huang et al. (2023) Table 4
rows B1/B2 that every helium-bearing system shares). The eight balance rows are
`mol_heh_rows` in `System_HeH_mol.f90`, solved per cell by `hybrd1` with the
row-turnover scaling of `set_mol_turnover_rates`, and merged with the metal
block by `System_HeH_mol_metals.f90`. The unknowns are H+, He+, He++, and the
H-nucleus fractions bound in H2, H2+, H3+ and HeH+, plus the He 2^3S metastable.

Two features of that network matter for what A2 has to add:

- It **already carries the channel that dominates a hot base.** R12
  (`H2 + M -> H + H + M`, Baulch et al. 1992) and its inverse R15 are in
  `mol_rates.f90` and in row 4 of `mol_heh_rows`. Section 2.2 shows that this is
  70-97% of the net H2 budget on HD 209458 b. A2 is therefore not needed to fix
  a hot base and must not change one.
- It has **no oxygen at all**, and the Lyman-Werner term added for exactly this
  test moves the partition by 0.5% of the gap (`TO_BE_DONE.md` item (H)). The
  gap is a missing set of *species*, not a missing rate.

### 2.2 The measured budget, and what it forces

Re-measured for this design from the P1 Photochem budget files (the pickles hold
`production` / `loss` arrays with their reaction strings on a 100-level grid;
the arms are C' = Photochem/NCHO, B' = Photochem/Zahnle, D' = Photochem/Zahnle+S).
The document's own P1.4/P1.5 tables were reproduced to the quoted digits; the
reaction-by-reaction rankings below are new and are cited to the pickles, not to the
document.

**HD 189733 b, 1 microbar, 864 K, arm C' — net H2 loss 8.617e6 cm^-3 s^-1**

| share of net | channel |
|---|---|
| 95.3% | `OH + H2 -> H2O + H` (net of its reverse) |
| 7.8% | `O(1D) + H2 -> OH + H` |
| 3.7% | `O + H2 -> OH + H` |
| 1.5% | `H2 + CN -> H + HCN` |
| -7.8% | `H2O + hv -> H2 + O(1D)` (an H2 *source*) |
| -0.5% | `H + H + M -> H2 + M` |

and the cycle that returns the water to OH: H2O loss is 54.5% `H2O + hv -> H + OH`,
38.1% `H2O + H -> OH + H2`, 5.1% `H2O + hv -> H2 + O(1D)`, 2.4%
`H2O + hv -> H + H + O`; H2O production is 100.0% `OH + H2 -> H2O + H`. OH loss
is 97.6% `OH + H2`, 2.4% `OH + H -> O + H2`, and 0.02% `OH + CO -> H + CO2`.

**HD 209458 b, 1 microbar, 2331 K, arm C' — net H2 loss 6.089e6 cm^-3 s^-1**

| share of net | channel |
|---|---|
| **70.4%** | `H2 + M -> H + H + M` |
| 20.6% | `OH + H2 -> H2O + H` |
| 7.8% | `O + H2 -> OH + H` |
| 1.7% | `O(1D) + H2 -> OH + H` |
| 0.02% | `CO + H2 -> HCO + H` |

**HD 209458 b, 1e-4 bar, 1830 K, arm C' — net 1.306e8 cm^-3 s^-1**: `H2 + M`
97.2%, `OH + H2` 1.5%, `O + H2` 1.2%, `CO + H2 -> HCO + H` 0.13%.

Four things follow, and they are the whole justification for the species list:

1. **The minimal reacting set is `OH + H2 <-> H2O + H` together with
   `O + H2 <-> OH + H` and the photolysis of H2O.** Those three carry
   95.3 + 3.7 + 7.8 (as O(1D)) percent of the HD 189733 b net and 20.6 + 7.8 + 1.7
   of the HD 209458 b net. Nothing else reaches 2%.
2. **O(1D) cannot be lumped into ground-state O.** It is 7.8% of the net H2 loss
   on HD 189733 b against 3.7% for ground O, while its mixing ratio there is
   5.35e-11 against 9.59e-8 for O — three decades less abundant and twice as
   effective. Lumping them would multiply the O + H2 channel by ~1800. It is
   carried, but it does not need transport (section 2.3).
   *M1 note:* those two shares are a correct report of a run made with the
   VULCAN NCHO rate set, whose `O(1D) + H2` coefficient is 2.6x the IUPAC
   evaluated value that M1 adopts. With the adopted rate the two channels are
   1.9% and 2.6% of the gross O1 rate, i.e. ground-state O is the larger of
   the two. The conclusion that O(1D) must be carried separately is unchanged
   and if anything stronger; the "twice as effective" ordering is not
   (`docs/a2_reaction_audit.md` section 3).
3. **CO must be carried as an oxygen reservoir even though no CO reaction
   matters.** Measured oxygen inventory, arm C':

   | case | H2O | CO | OH | O | CO2 |
   |---|---|---|---|---|---|
   | HD 189733 b, 1 microbar | 54.5% | 45.1% | ~0 | ~0 | 0.3% |
   | HD 209458 b, 1 microbar | 42.9% | 45.8% | 7.3% | 4.0% | ~0 |
   | HD 209458 b, 1e-4 bar | 54.0% | 45.8% | 0.2% | ~0 | ~0 |

   CO locks 45-46% of the oxygen at every level in every arm. A network that
   gives all of `melem_ab(iel_O)` to the water family over-supplies the OH cycle
   by about a factor of two. The largest CO channel measured anywhere is
   `CO + H2 -> HCO + H` at 0.13% of the net, so CO is required as a *sink*, not
   as a reactant.
4. **Sulfur, nitrogen and hydrocarbon channels are excluded on measurement.**
   `S + H2 -> H + HS` is 40.5% of the gross HD 189733 b loss and 0.36% of the
   net; `H2 + CN` is 1.5% of the net; `NH2 + H2 <-> NH3 + H` reaches 2.8% of the
   *gross* at the 1e-4 bar level of HD 189733 b and cancels. The one number that
   argues for revisiting sulfur is P1.5's own: it moves `q_H2` by 17% at 1e-6 bar
   on HD 209458 b. That belongs in section 8 as a decision, not in the minimal
   set.

### 2.3 Species carried

| species | `f_sp` role | transported? | why |
|---|---|---|---|
| H2, H2+, H3+, HeH+ | existing columns 34-37 | **yes** (new) | the partition being solved for; `tau_chem ~ tau_adv` at the cool base |
| H I, H II | existing columns 1-2 | **yes** (new) | closes the H-nucleus budget against H2 |
| H2O | new column | yes | 43-55% of the oxygen; the cycle's reservoir |
| OH | new column | yes | the carrier of the dominant channel |
| O (ground) | existing column 10 (O I), redefined | yes | already transported when `He_metal_diffusion` is on |
| O(1D) | new column, diagnostic only | **no** | steady state; see below |
| CO | new column | yes | inert oxygen reservoir; 45% of the element |
| He, He+, He++, He 2^3S, all metals | existing | as today | unchanged |

**O(1D) is closed by a local steady state, not transported.** Its measured
mixing ratio is 5.35e-11 (HD 189733 b, 1 microbar) against 9.59e-8 for ground O.
It has exactly one source in this set — the `H2 + O(1D)` branch of H2O
photolysis, which carries a quantum yield of 0.11 below 1201 A, 0.10 in the
Ly-alpha window and 0.00 above 1451 A (section 2.6) — and its dominant sink in
an H2-rich gas is the same reaction that makes it useful, `O(1D) + H2 -> OH + H`.
Its chemical lifetime is shorter than every other time scale in the problem by
many decades, so the steady-state approximation is the *accurate* treatment, not
a shortcut; it must be written at the code site with that justification and with
the ratio above as the evidence. Carrying it as a transported column would add a
stiff row to the transport solve for a species that cannot move.

**CO is carried but frozen chemically in the default network.** It occupies a
column, it is transported (so that its oxygen goes where the gas goes), and its
only chemical coupling is through the total-oxygen closure. Its abundance at the
base comes from the handoff or from a new `C_H`-derived partition; the choice is
a decision (section 8, D4). What must not happen is CO's oxygen silently
appearing in the OH cycle.

### 2.4 The reaction set

The rule of the plan is absolute and is repeated here: **no rate or cross
section from memory. Every reaction is transcribed from a named published
source, with its validity range, at the code site**, exactly as R1-R23 are in
`mol_rates.f90`. This design fixes *which* reactions and *which source to
transcribe from*; the numbers are copied at M1 from the publication, with the
in-workspace network files used as the cross-check that the transcription is
right.

**The sources are already on this machine, in two independent networks.**
Photochem's mechanism is `photochem/data/reaction_mechanisms/zahnle_earth.yaml`
(Kevin Zahnle's `earth_125.rx`, maintained by Wogan; provenance for each reaction in
its `ref:` keys, resolved through `photochem/data/bib.bib`). VULCAN's is
`VULCAN/thermo/NCHO_photo_network.txt` (`k = A T^B exp(-C/T)`, odd ids forward,
even ids the thermodynamic reverse). Having both means every transcription can
be checked against two independent ports of the same literature before it enters
EXHALE — which is the discipline the R16-R20 discrepancy notes in
`mol_rates.f90` were written after the fact.

Proposed set, ordered by its share of the measured budget. **The source column
is as corrected by M1**; where M1 adopted a later evaluation than the design
first proposed, the superseded source is named too, and the size of the change
is measured in `docs/a2_reaction_audit.md` section 4.

| id | reaction | why it is in | source |
|---|---|---|---|
| O1 | `OH + H2 <-> H2O + H` | 95.3% of the HD 189733 b net; 20.6% of HD 209458 b. Its reverse is 38.1% of H2O loss | **Baulch et al. (2005), JPCRD 34, 757, p. 1029**, `3.6e-16 T^1.52 exp(-1740/T)`, 250-2500 K. This supersedes the Baulch et al. (1992) recommendation that `zahnle_earth.yaml` carries under the key `Ba92`, and it is the value VULCAN already has from Oldenborg & Loge (1992), to 0.8%. The two reference networks are the two successive evaluations of one reaction, not two readings of one evaluation |
| O2 | `O + H2 <-> OH + H` | 3.7% (HD 189733 b) / 7.8% (HD 209458 b) of the net; its reverse is 2.4% of OH loss | **Baulch et al. (2005), p. 804**, `6.34e-12 exp(-4000/T) + 1.46e-9 exp(-9650/T)`, 298-3300 K. Both networks carry the superseded 1992 form and both carry it with `exp(-3160/T)` where the published value is `exp(-3163/T)`. The 2005 form is 1.84x smaller at the HD 189733 b base temperature — the largest single rate change M1 makes |
| O3 | `H2O + hv -> OH + H` | **54.5% of H2O loss** — the branch that makes the cycle catalytic rather than a one-way sink | cross section: the concatenation `photochem/data/xsections/H2O.h5` makes from Huebner & Mukherjee (2015) below 6.3 nm, Heays et al. (2017) (Leiden) to 192.056 nm and Ranjan et al. (2020) to 230.413 nm — those cut points are Photochem's, not any paper's. Branching ratios: Stief, Payne & Klemm (1975) and Slanger & Black (1982), **as summarized in** JPL Publication 19-5 entry B2 (Burkholder et al. 2020), which does **not** itself recommend H2O quantum yields and leaves the Ly-alpha branching unsettled between two disagreeing sets (`docs/a2_reaction_audit.md` section 7) |
| O4 | `H2O + hv -> H2 + O(1D)` | 5.1% of H2O loss and the second-largest net H2 *source* (-7.8%) | same |
| O5 | `H2O + hv -> O + H + H` | 2.4% of H2O loss | same |
| O6 | `O(1D) + H2 -> OH + H` | the reason O(1D) is carried at all; 1.9% of the gross O1 rate with the adopted coefficient | **Atkinson et al. (2004), Atmos. Chem. Phys. 4, 1461, data sheet I.A2.18**, `1.1e-10`, temperature-independent 200-350 K. VULCAN's `615` gives 2.87e-10 cited to Tully (1975), which is an RRKM extrapolation contradicted by the four measurements the IUPAC evaluation averages; `zahnle_earth.yaml` gives 1.5e-10 under the key `Ba92`, an attribution that does not hold because no Baulch evaluation contains O(1D) chemistry. **The factor 2.6 between these is the largest open uncertainty in the set** |
| O7 | `OH + hv -> O + H` | 0.01% of OH loss at the measured levels, but the only OH photolysis channel there is, and it is not small in the optically thin part of the column | Huebner & Mukherjee (2015) below 82.8 nm, Heays et al. (2017) to 264.9 nm; single branch, quantum yield 1 at every wavelength |
| O8 | `O + H + M <-> OH + M` | the three-body association of atomic O with H at the dense base; the oxygen analogue of R15. **Excluded** by M1: 7.8e-9 of O1 | Tsang & Hampson (1986), JPCRD 15, 1087, p. 1111, `k0 = 1.3e-29 T^-1` — an estimate with an uncertainty factor of 10 and **no stated temperature range**, the source saying there are no definitive measurements. Falloff, so the same third-body caveat `mol_rates.f90` states for R12/R13/R15 applies; the base is at `Pr ~ 1e-8`, so the high-pressure limit is irrelevant there and neither network's `k_inf` comes from the source anyway |
| O9 | `H + OH + M <-> H2O + M` | the three-body route to water that does not go through O1. **Excluded** by M1: 9.7e-8 of O1 | **Baulch et al. (1992) pp. 496-498, unchanged in Baulch et al. (2005) p. 913**, collider-resolved: `k0(N2) = 6.1e-26 T^-2.0`, `k0(Ar) = 2.3e-26`, `k0(H2O) = 3.9e-25`, all 300-3000 K. This resolves the apparent two-decade disagreement between the networks: VULCAN's `659` **is** Baulch's H2O-collider value, and Photochem's Javoy et al. (2003) number has a published validity range of 2790-3200 K. Neither is right for an H2/He bath |
| O10 | `OH + OH <-> H2O + O` | not in the measured top channels. **Excluded** by M1: 5.5e-7 of O1 | **Baulch et al. (2005), p. 1032**, `5.56e-20 T^2.42 exp(+970/T)`, 250-2400 K. `zahnle_earth.yaml` carries the superseded 1992 form and labels it `Ba92, Li91`; Lifshitz & Michael (1991) is a study of the **reverse** direction and is not in the 1992 data sheet's reference list at all. VULCAN carries the reverse direction explicitly (`O + H2O -> OH + OH`), but with no reference, and M1 found evidence that it was itself obtained by reversal |
| O11 | `O(1D) + M -> O + M` (quenching) | closes the O(1D) steady state where H2 is depleted. **Excluded** by M1, with a published reason | the IUPAC value adopted for O6 is the **total** of the reactive and the quenching channel, and the same data sheet states that the reactive channel is above 95% of it. Quenching by H2 is therefore inside O6's coefficient, under 5% of it, and has no separately evaluated rate to transcribe |
| O12 | `O(1D) + H2O -> OH + OH` | added by M1 to fill the slot section 7 names: the second O(1D) sink, and the one channel that could have upset the O(1D) steady state because H2O holds 43-55% of the oxygen. **Excluded**: 2.8e-5 of O1, and it shortens the O(1D) lifetime by 0.15% | **Atkinson et al. (2004), data sheet I.A2.19**, `2.2e-10`, 200-350 K — which is the value `zahnle_earth.yaml` carries. VULCAN's `619` is 21% below it |

O1-O7 are the closed set the measurement supports, and M1 confirmed it.
O8-O12 are the ones the plan's "conservation skeleton to be audited at reaction
level before any coding" names. **M1 excluded all five** — four of them on
measurement, at 7.8e-9 to 2.8e-5 of the dominant channel, and O11 because the
rate it would need is already inside O6's coefficient. Each is transcribed in
`oxygen_rates.f90` anyway so that the record is a value with a verdict rather
than an omission.

**Explicit reverse reactions are not written.** O1's reverse (`H2O + H -> OH + H2`)
and O2's (`OH + H -> O + H2`) do not appear as separate entries in either
network: both codes generate them thermodynamically. VULCAN went further and
*removed* its explicit `H2O + H` channel in October 2020 in favor of the reversed
`OH + H2` (its changelog line 7). Section 2.5 takes the same route, and this is
the evidence for it.

**No ion chemistry of the water family.** H3O+, OH+, H2O+ and O+ + H2 are
excluded. The plan's standing decision is that no O-bearing molecule enters the
trace-metal ionization framework, because that solver's assumptions (coronal-like
rates, no three-body reactions) do not hold for molecules; and the measured
oxygen at the levels that matter is 100% neutral. The consequence for the O
element's charge state is that A2 changes nothing about O I/II/III above the
molecular layer, which is what gate G2 (section 6) tests.

### 2.5 Reverse reactions and detailed balance

Two policies are available and they are not equivalent:

- **Explicit reverse rates**, each direction transcribed independently. Simple,
  no thermodynamic data needed, and it is what `mol_rates.f90` does today. Its
  defect is visible in the code already: the module header records that the
  R21/R22 pair does not satisfy detailed balance against itself, by about two
  decades at 1e4 K. For a pair like O1 that runs in near-exact cancellation on a
  hot base — the HD 209458 b gross rates are three decades above the net — an
  inconsistent pair does not produce a small error in the net, it produces an
  arbitrary one.
- **Thermodynamic reversal**: one direction transcribed, the other computed from
  the equilibrium constant built from tabulated enthalpies and entropies. Exact
  detailed balance by construction, and it makes the hot limit reduce to
  chemical equilibrium — the limit `q_h2_equilibrium` already encodes and the
  one the HD 209458 b budget says is the correct physics there.

**The design proposes thermodynamic reversal for O1, O2, O8, O9 and O10.** The
argument is not only the near-cancellation: **both reference codes do exactly
this, and neither carries an explicit reverse.** Photochem computes
`k_rev = k_fwd / K_c` with the standard-state factor `(k_B T / 10^6)^(dn)`
(`photochem/src/photochem_common.f90:161-186`), from **Shomate** polynomials
carried in the `species:` block of `zahnle_earth.yaml` and sourced from NIST.
VULCAN computes the same quantity from **NASA-9** polynomials in
`VULCAN/thermo/NASA9/`, sourced from Burcat's database, with the formulae
written out in `VULCAN/thermo/gibbs_text.txt`.

Cost, measured rather than estimated: five species (H, H2, O, OH, H2O) at seven
Shomate coefficients over one to four temperature ranges each — 91 numbers, or
112 with CO — plus one `gibbs_energy_shomate` function. The NASA-9 route from VULCAN
is the same size and is directly portable. Either is small; the choice between
them is D2, and it fixes which reference file the coefficients are copied from.

One consequence worth stating in advance: with thermodynamic reversal, gate G3
(the hot limit reproduces chemical equilibrium) becomes a real test of the
implementation rather than a restatement of it.

### 2.6 Photolysis and the FUV bands

Reactions O3-O5 and O7 need a photon field below the H I edge. Section 1.2 gives
the reason this cannot simply extend the XUV grid, so the field arrives as bands
with named fluxes, on the pattern of `lyman_werner.f90`.

**The band edges are not free: the branching ratios fix them.** Both reference
codes carry wavelength-resolved quantum yields for H2O photolysis built from
Stief, Payne & Klemm (1975) and Slanger & Black (1982), by way of the summary in
JPL Publication 19-5 entry B2
(`photochem/data/xsections/H2O.h5`, `photodissociation-qy`;
`VULCAN/thermo/photo_cross/H2O/H2O_branch.csv`). **The four ratio triplets are
the same numbers but the two files are not the same function of wavelength**:
VULCAN's Ly-alpha node sits at 121.0 nm, 0.567 nm short of the line, and VULCAN
interpolates linearly between nodes, so at the line center it uses
0.837 / 0.105 / 0.058 — the three-body branch is 0.058 there, not 0.12 — and it
has no flat 1231-1450 A interval at all. The table below is the reading of
`H2O.h5`, whose yields are constantly extrapolated between nodes. Reading them
onto intervals:

| interval | `OH + H` | `H2 + O(1D)` | `O + H + H` |
|---|---|---|---|
| below 1201 A | 0.89 | 0.11 | 0.00 |
| 1202-1230 A (contains Ly-alpha, 1215.67 A) | 0.78 | 0.10 | **0.12** |
| 1231-1450 A | 0.89 | 0.11 | 0.00 |
| above 1451 A | 1.00 | 0.00 | 0.00 |

The Ly-alpha interval is the only one where the three-body branch
`H2O + hv -> O + H + H` opens at all, and it is 12% there. **That is an
independent, data-driven reason to give Ly-alpha its own band** rather than
folding it into a wider FUV average — quite apart from its being a line rather
than a continuum. OH photolysis is carried as a single branch with quantum yield 1 at every
wavelength (`OH.h5`; `OH_branch.csv`), so it needs no such split — but **that
unit yield is a merge, not a measurement**: Heays et al. (2017) section 3.1
state that they do not divide the photoabsorption cross section into decay
channels except for H2O and NH3, and VULCAN's `OH_branch.csv` says in its own
header that it combines the O(1D) branch below 150 nm with the O(3P) branch
above it. O(1D) made this way feeds straight back into O6, so M2 must either
split the channel or state that it does not.

The cross-section grids, measured from the Photochem files: H2O 0.1-230.413 nm
(peak photodissociation cross section 2.75515e-17 cm^2 at 111.5 nm), OH
0.06-264.9 nm (peak 1.36397e-17 cm^2 at 100 nm). VULCAN's `H2O_cross.csv`
covers only 6.2-230.413 nm and peaks 5.7% lower, at 2.598e-17 cm^2.

Proposed channels:

| channel | interval | flux source | attenuation |
|---|---|---|---|
| B1 | 912-1201 A | the FUV band flux (section 5.1), with `Stellar LW flux` continuing to drive H2 over its own 912-1110 A interval | H2O and OH **continuum** `exp(-tau)` on their own columns; the H2 branch keeps the Draine & Bertoldi (1996) line self-shielding it has today |
| B2 | Ly-alpha, 1215.67 A | the existing key `Stellar Lya flux [erg/cm2/s]` (`F_Lya_star`); when `Jlya escape-prob` or an imported field is active, the solved `J_Lya(r)` of `lya_rt.f90` is the better field and should be used instead | H I line transfer is what `lya_rt.f90` already solves; H2O and OH absorb it as continua |
| B3 | 1231-1450 A | the FUV band flux | continuum on H2O, OH |
| B4 | 1451-2304 A | the FUV band flux | continuum on H2O (OH's cross section runs to 2649 A and is included over the same interval) |

Three constraints, each imposed by the code or the data rather than chosen:

- **H2O and OH attenuate as continua, so no self-shielding function is needed for
  them.** The `f_shield` treatment of `lyman_werner.f90` exists because H2
  absorbs in lines that saturate; a continuum absorber is `exp(-sigma N)` with
  the same column integration `calc_column_dens_one` that every other absorber
  uses (`ionization_equilibrium.f90:296`).
- **Band-averaged cross sections carry a spectral-shape sensitivity that must be
  written at the code site**, in the terms `lyman_werner.f90` uses for
  `sigma_LW` — it quotes +-25% between two bracketing color temperatures and
  states when the approximation fails (a band dominated by one emission line).
  For an M dwarf, whose FUV *is* dominated by lines, B1 and B3 are exactly that
  failure case, and the header must say so. B2 has no averaging question.
- **`Stellar LW flux` and the FUV band flux describe overlapping photons for
  different absorbers.** Over 912-1110 A the same photons are absorbed by H2 (in
  lines) and by H2O and OH (in a continuum); that is physics, not double
  counting. What would be double counting is the *energy*, so gate G4 checks each
  band's deposited heat against its own incident flux and requires the sum over
  absorbers in an overlapping interval to stay below it. A run that fails that
  check is outside the treatment's validity range and must say so rather than
  quietly over-heat.

**A run with `Stellar LW flux` set and no FUV band flux gets no oxygen
photolysis longward of 1110 A** — i.e. it loses the branch that carries 54.5% of
H2O loss. That is a warning at startup, not a silent zero.

**Photolysis energetics.** Each dissociation deposits the excess above the bond
energy as fragment kinetic energy. The precedent is `e_lw_fragment_erg`
(`lyman_werner.f90:132`), added to `heat` in `ionization_equilibrium.f90:412`
and reported in the heat breakdown (`util_ion_eq.f90:1168`). A2 adds one such
term for each photolysis channel in each band, since the excess energy depends on the
band and on the branch. Gate G4 tests that the absorbed band energy and the
deposited heat balance.

### 2.7 What is deliberately left out, and the validity range that follows

| left out | measured size | consequence |
|---|---|---|
| sulfur (`S + H2 <-> H + HS`) | 0.36% of the HD 189733 b net; 19.3% of the HD 209458 b 1-microbar net; moves `q_H2` by 17% at 1e-6 bar on HD 209458 b (P1.5) | A2 is not valid for a sulfur-rich base; the hot-base partition is thermal anyway, so the affected regime is the one where A2 is not the answer |
| nitrogen (`H2 + CN`, `NH2 + H2`) | 1.5% of the HD 189733 b net; 2.8% of the gross at 1e-4 bar | a stated few-percent term |
| hydrocarbons beyond CO (`CH3 + H2`, `CH + H2`) | <= 1.0% gross, <= 0.01% net | negligible at the measured levels |
| CO chemistry (`CO + H2 -> HCO + H`, `OH + CO -> H + CO2`) | <= 0.13% of any net budget | CO stays an inert reservoir |
| CO2 | 0.3% of the oxygen on HD 189733 b, ~0 elsewhere | folded into CO or dropped (decision D4) |
| O2, HO2, H2O2, HCO, H2CO | every one of them below `1e-8` by volume at every measured level in every arm (O2 is the largest, `1.04e-8` on HD 209458 b at 1 microbar and `1.10e-10` on HD 189733 b) | dropped outright |
| water-family ions | the measured oxygen at these levels is neutral | see section 2.4 |

**Validity range to be written into the module header.** A2 is valid where the
base is cool enough for the partition to be kinetic and oxygen-carried — the
regime P4 measured on HD 189733 b (864 K, oxygen family 96-99.6% of the net
H2 destruction). Where the base is hot enough that `H2 + M -> H + H + M` runs
the net (HD 209458 b, 2331 K, 70-97%), A2 adds a few percent and the existing
R12/R15 pair is the physics. Where the base is sulfur-rich, A2 is missing a
19% term. The header states this the way `lyman_werner.f90` states the +-25%
band-shape range: as the condition under which the approximation holds, at the
code site.

---

## 3. Transport

### 3.1 Why local kinetics fails exactly where it is needed

The P4 measurement, restated because the whole of this section follows from it:

| planet | level | `tau_chem(H2)` [s] | `tau_adv = H/v` [s] | ratio | ratio on `r/v` |
|---|---|---|---|---|---|
| HD 209458 b | 1 microbar | 1.81e5 - 2.39e5 | 2.132e8 | 8.5e-4 - 1.2e-3 | 1.0e-5 - 1.3e-5 |
| HD 209458 b | 1e-4 bar | 1.33e6 - 2.44e6 | 1.965e8 | 6.8e-3 - 1.2e-2 | 4.8e-5 - 8.9e-5 |
| HD 189733 b | 1 microbar | 2.99e5 - 2.51e6 | 1.498e6 | **0.20 - 1.67** | 4.0e-4 - 3.3e-3 |

and, pressure-matched at the EXHALE base cell, HD 189733 b gives **16.4**. The
caveat P4 attaches pushes it further: the chemistry model's `T(p)` is hotter
than EXHALE's at the same pressure (883 K against 566 K on HD 189733 b), so the
true `tau_chem` under EXHALE's temperature is *longer* than quoted.

The planet whose base partition oxygen chemistry decides is the planet whose
base partition is not in local chemical equilibrium. That is the finding, and
it is why a transport-free A2 would be a solution to the wrong problem.

### 3.2 What EXHALE transports today, and what it does not

Read from the source, because this is the load-bearing constraint on the whole
section:

- **The hydro state vector is `(rho, rho v, E)` only.** `f_sp` is not advected
  by the RK update; it appears in no flux and no reconstruction. It is
  recomputed from scratch each step by `ioniz_eq`, cell by cell, as the root of
  a local algebraic system (`ionization_equilibrium.f90:709-810`). **Species
  composition in the marching loop is local equilibrium by construction.**
- **The advection correction lives only in post-processing.** `post_process_adv.f90`
  reconstructs the advected ionization state on the converged wind, and its
  header states plainly that it **excludes the molecular species** and that
  "a molecular base needs a molecular-aware post-process instead"
  (`post_process_adv.f90:5-17`).
- **Elements move, and only elements.** `binary_element_diffusion.f90` transports
  the helium *mass fraction* `X = rho_He/rho` (its `element_nucleus_counts`
  sums `bsp_nH`/`bsp_nHe` over every carrier, molecular carriers included since
  milestone M4) and, with `He_metal_diffusion`, each metal element's mixing
  ratio `n_X/n_H`. Oxygen's nucleus count there is the bare sum over
  O I + O II + O III; there is **no `bsp_nO` and no molecular oxygen carrier
  anywhere in the tree** (`n_species = 37`, `n_bsp = 10`).
- **The operator is split and Picard-coupled.** `element_diffusion_step` is
  called after the hydro update and before `ioniz_eq`
  (`EXHALE_main.f90:696-697`); the steady route alternates
  `solve_steady_jfnk` with `relax_element_composition` under a damped Picard
  loop with `omega` starting at 0.5 and halving to a floor of 0.125
  (`EXHALE_main.f90:1385-1424`).

So A2 needs something that does not exist: **species-resolved vertical
transport**. It also has an exact template for how to build it, in the operator
that already does the element version.

### 3.3 The equation

For each transported species `i` of the A2 set, in the spherical 1-D geometry
the rest of the code uses:

```
  d n_i / dt  +  (1/r^2) d/dr [ r^2 ( n_i v + Phi_i ) ]  =  P_i - L_i        (1)
```

with `v` the hydro's mass-averaged velocity, `P_i`, `L_i` the chemical
production and loss of section 2, and the vertical flux in the standard
aeronomy form (Banks & Kockarts 1973), written on the mixing ratio
`f_i = n_i / n_tot` so that the eddy term has no preferred species — the same
argument `binary_diffusion_design.md` section 2.3 makes for `K_zz` acting on
`dX/dr` alone:

```
  Phi_i = - n_tot ( D_i + K_zz ) d f_i / dr
          - n_i D_i [ 1/H_i - 1/H_atm + (alpha_T,i) dlnT/dr
                      - (Z_i e E)/(k T) ]                                     (2)
```

`H_i = kT/(m_i g)` is the species scale height, `H_atm = kT/(mu m_H g)` the
atmospheric one, `E` the ambipolar field, and `Z_i` the charge (zero for every
neutral A2 species; the molecular ions are the exception, section 3.4). The
sign convention and the ambipolar field are **already solved** in
`binary_element_diffusion.f90`: `settling_coefficient` (lines 1388-1493) builds
`eE = -kT dln(n_e T)/dr` from the cell's solved electron pressure
(lines 1443-1452), and the design document's sign check ("helium drifts inward
relative to hydrogen") is the one every discrete form must reproduce.

**Conservation.** Summing (1) over the carriers of an element, weighted by the
nuclei each carries, must give the element equation the existing operator
solves. That is not automatic — it is a constraint on the discretization
(section 3.6) and it is gate G5.

### 3.4 What is reused, and what is new

| piece | status |
|---|---|
| `D_i` for a neutral in an H2/H/He mixture | **reuse the existing coefficient family**: `hard_sphere_pair_diffusion`, `polarization_pair_diffusion`, `ion_neutral_pair_diffusion`, `coulomb_pair_diffusion`, combined by the Blanc rule in `stage_mixture_diffusion`. All are `public` in `binary_element_diffusion`. New entries are needed in the carrier tables for the polarizabilities of OH, H2O and CO; `alpha_melem` already carries atomic O (5.3 a0^3, line 350) |
| ambipolar field `eE` | reuse `settling_coefficient`'s computation; every A2 species is neutral, so the term is zero for them, but it must not be *dropped* — the molecular ions H2+/H3+/HeH+ do carry charge, and they are in the transported set |
| `K_zz` | reuse `kzz_cell` unchanged; section 3.5 |
| Peclet-hybrid face coefficients, donor-cell switch | reuse the pattern of `drift_and_gradient_face_coefficients` (lines 1552-1567): central where `|B| dr <= 2(A+E)`, donor-cell otherwise |
| backward-Euler nonlinear step with a tridiagonal Newton | reuse the pattern of `solve_mass_fraction` (lines 1612-1774), generalized from one scalar to a block-tridiagonal system in `n_A2` species |
| the Dirichlet base and zero-gradient top | reuse the pattern (lines 485-522); section 3.7 |
| damped Picard co-convergence with the steady solver | reuse `steady_wind_with_element_diffusion` (EXHALE_main.f90:1336-1426); section 3.8 |
| the element write-back `project_elements` | **must be extended**, not reused as is: it slaves every `mion` column to `r_H` (lines 1900-1903) and re-seeds an absent element in its neutral stage. A molecule carrying both H and O is in the position `HeH+` occupies for the H/He pair, where the code already uses `rBoth = min(rH, rHe)` (line 1893) |
| `element_nucleus_counts` | **must gain an oxygen arm.** Today it sums `bsp_nH` and `bsp_nHe` only. Without a `bsp_nO` analogue, every O nucleus bound in OH, H2O or CO disappears from the metal operator's `nX` and `project_elements` scales it as if it were hydrogen |
| `m_1 = mass_per_H_nucleus_without_He()` | **must be revisited.** It is `bsp_mass(isp_HI) + sum(melem_ab*melem_A)` (`composition.f90:214-229`), i.e. the binary H/He closure assumes a *fixed* metal/H. Oxygen chemistry that moves O relative to H breaks that assumption; the mass-closure diagnostic at `binary_element_diffusion.f90:660-664` is what would report it |

The last three rows are the real integration cost of A2, and they are the
reason the P4 block calls it "one total-oxygen closure spanning the metal block
and the molecular block" rather than "add two species".

**Whether A2's transport is a new module or a generalization of
`binary_element_diffusion` is a decision** (section 8, D3). The design leans to
a new module for the species system, leaving the element operator as it is: the
two solve different variables (species number densities with chemical sources
against element mass fractions with none), and the element operator's
correctness rests on a two-component closure that a multi-species system does
not have.

### 3.5 The eddy coefficient has exactly one owner

`kzz_cell(1-Ng:N+Ng)` in `global_parameters` is the array every diffusive term
reads. It is filled once by `eddy_diffusion_on_grid`
(`lower_atmosphere_profile.f90:465-529`, called from `init.f90:79`): from the
scalar `he_kzz` when there is no profile, and by interpolation of the profile's
`Kzz` column when there is. `base.inp`'s `Kzz_base` is refused when a profile is
in use, and an `He_Kzz:` key in `input.inp` is inert and warned about
(`input_read.f90:1449-1458`).

**A2 introduces no eddy coefficient of its own and no key of its own for one.**
It reads `kzz_cell`. Any other arrangement would recreate, inside one code, the
two-owner problem section 4 exists to prevent.

One consequence must be stated rather than discovered: with A2 on and
`He_diffusion` off, `kzz_cell` is still filled (it is `he_kzz`, default 0), so
an A2 run with no `He_Kzz` and no profile has **pure molecular diffusion**. On a
lower atmosphere that is the wrong limit — eddy mixing is what holds the
composition well mixed below the homopause. A2 must therefore either require a
positive `K_zz` or state loudly in the setup report that it is running without
one. Proposed: a warning, not a refusal, so that the `K_zz = 0` limit stays
available as a test.

### 3.6 Discretization

- **Variable.** The mixing ratio `f_i = n_i/n_tot`, not `n_i`. It is the
  quantity the flux (2) is naturally written on, it is bounded, and it is what
  the eddy term mixes toward uniformity. `n_tot` is the EOS particle count
  `calc_ntot` already provides (electrons excluded), which is also the third
  body `M` the existing R12/R13/R15 use (`ionization_equilibrium.f90:267-278`).
- **Grid.** Cell-centered `f_i`, face-centered fluxes on `r_edg(0:N)`, exactly
  as the element operator does. Faces 0 and N carry zero diffusive flux.
- **Advection.** The one-sided upwind difference on the *cell* velocity
  `rho_j v_j` that the element operator uses, and for the reason recorded in its
  header (lines 44-64): the face-averaged form froze an isolated helium hole at
  the breathing base. A2 inherits that lesson rather than rediscovering it.
- **Chemistry coupling.** The chemical source `P_i - L_i` is stiff. Two options,
  and this is a decision (section 8, D5):
  - **(a) Fully coupled implicit step.** One backward-Euler step of the whole
    `n_A2 x N` system, block-tridiagonal in space and dense in species, solved
    by Newton. This is what VULCAN and Photochem do, it has no splitting error,
    and it is the only form in which the near-cancelling O1/O2 pair is safe at
    a hot base.
  - **(b) Operator split**: an implicit local chemistry substep followed by an
    implicit transport substep, iterated. Cheaper per step and closer to the
    existing code's shape, but the splitting error is largest exactly where the
    gross rates exceed the net by three decades — which the HD 209458 b budget
    says is the normal case.
  The design proposes **(a)**, on the ground the budget supplies: a scheme whose
  error scales with the gross rates cannot be used on a system whose answer is
  the difference of two numbers three decades larger than itself.
- **Positivity and the element simplex.** The existing local solve already has a
  full apparatus for this — physicality tests, three seeding attempts, and
  `clamp_fractions_to_element_budget` for a root that lands outside
  (`ionization_equilibrium.f90:738-812`). The transported system needs the same
  discipline stated up front: `f_i >= 0` enforced by the discretization where
  possible (the element operator bounds `X` in `[0,1]` by taking the two factors
  of `X(1-X)` from opposite sides of the face, lines 1588-1606), and any clamp
  counted and reported, never silent.

### 3.7 Boundary conditions

**Inner (the base).** Dirichlet, as the element operator's helium reservoir is
(`X_base` from `HeH`, lines 485-495). The values come from whichever handoff the
run is using:

- with a lower-atmosphere **profile**, from the profile at the matching level —
  it already carries the elemental reservoirs and would carry the molecular
  mixing ratios (which are today "diagnostic metadata only", A1c);
- with the scalar `base.inp`, from `q_H2_base` and `<El>_H_base`;
- with neither, from the chemical-equilibrium fit `q_h2_equilibrium` for
  hydrogen and from `melem_ab` for oxygen and carbon, with the partition among
  H2O/CO/OH/O taken at chemical equilibrium — which is a statement of what the
  base is *assumed* to be and must be printed as such.

**A Dirichlet base is a choice, not the only one.** The alternative is a
zero-flux (or specified-flux) base, which lets A2 *compute* the base partition
instead of being told it. That is the more ambitious reading of "the composition
is an output" that Phase E reached for the elements, and it is exactly what
makes the A/B gate of section 6.2 meaningful or circular. The design proposes:
**Dirichlet on the elements (H, He, O, C reservoirs) and zero-flux on the
partition among the carriers of each element.** The element totals are the
boundary data the handoff legitimately supplies; the partition among H2/H/H2O/OH
is what A2 exists to compute, so it must not be imposed at the boundary. This is
a decision (section 8, D6) because it decides whether the first gate tests
anything.

**Outer.** Zero gradient on `f_i`, as the element operator does
(`Xhe(N+1:N+Ng) = Xhe(N)`, line 522). The molecular species are negligible
there by many decades in every configuration A2 targets, so the boundary is
inert; a run in which it is not is a run outside A2's validity range and the
diagnostic must say so.

### 3.8 Time coupling and the outer loop

**Marching loop.** The A2 transport-chemistry step replaces the molecular part
of `ioniz_eq` for the species it owns, and sits in the same slot the element
diffusion occupies — after the hydro update, before the atomic/metal ionization
solve (`EXHALE_main.f90:696-705`). The atomic and metal ionization solve keeps
its local-equilibrium form: nothing in the measurement questions local
equilibrium for the ionization stages in the wind, and gate G2 requires that
limit to be reproduced.

**Steady route.** `steady_wind_with_element_diffusion` already alternates a JFNK
steady solve with a composition relaxation under damped Picard
(`omega = 0.5`, halved on failure to descend, floor 0.125, at most 20 outer
passes, drift tolerance 1e-3). A2's relaxation joins that loop as a second
composition solve, with the same damping and the same drift measure. The header
of that routine states the reason the damping exists — "nothing in a Picard
iteration of two solves ... keeps the two from chasing each other, and on the
HD 209458 b `Kzz = 0` wind they do" — and A2 adds a third participant to the
same iteration, so the risk is larger, not smaller. Milestone M2's gate is that
the loop converges on a case where it converges today.

**The steady residual does not contain the chemistry.** The existing steady
residual (`steady_newton.f90`) contains no diffusion either — that is why the
Picard loop exists. A2 keeps that structure. What must be checked, and is gate
G8, is that the fixed point of the marching loop and the zero of the steady
residual are still the same state.

---

## 4. Ownership: one oxygen, one hydrogen

### 4.1 The total-oxygen closure

The plan's requirement, stated as code:

```
n_O,tot(j) = n_OI(j) + n_OII(j) + n_OIII(j) + n_OH(j) + n_H2O(j) + n_CO(j)
           = melem_ab(iel_O) * n_H,nuclei(j)     [when nothing has moved O]
```

Today `n_O,tot` is formed in exactly one place per consumer, and each of them
has to change:

| site | today | with A2 on |
|---|---|---|
| `ionization_equilibrium.f90:236-243` (`nm_tot`) | sum of the three stages | the same sum, but it is now *free atomic* oxygen and it is no longer the element total |
| `binary_element_diffusion.f90:551-560` (`nX`) | sum of the three stages | must add the molecular carriers, via a new `bsp_nO`-style count |
| `binary_element_diffusion.f90:712-715` (`metal_hydrogen_ratio_departure`) | ratio against `melem_ab` | same, on the corrected total |
| `binary_element_diffusion.f90:1912-1918` (`project_elements`) | re-seeds an absent element in the neutral stage | needs an explicit rule for a molecule carrying two elements, like the `rBoth = min(rH, rHe)` rule HeH+ already has |
| `src/utils/element_budget.py` | closes O against `EXHALE_resolved.out` | must count the new carriers; it is the gate (G5) |

### 4.2 The metal block's O I, redefined

`System_HeH_mol_metals` solves each metal element's stage fractions against
`met_ntot(e)`, which the driver fills from `nm_tot(j,im)`. With A2 on, the value
handed to the oxygen row becomes the free-atomic-oxygen total, and the O I
column of `f_sp` (column 10, `mion_fsp(im_OI)`) means free neutral atomic O.

Three consumers read that column and each needs a verdict written into this
design rather than discovered later:

- **Cooling.** `Cool_coeff.f90` builds the [O I] 63/145/44 um fine-structure
  cooling and `cool_OI_ne_func` from the O I density. Free atomic O is the
  correct reactant, so the redefinition makes the existing cooling *more*
  correct, not less — but it changes its value wherever oxygen is molecular, and
  that is a physics change to report, not a refactor.
- **Charge exchange.** The oxygen pairs are cases 13/14, 30/31, 48/49 of
  `charge_exchange.f90`. They react with atomic O. Same verdict.
- **The transit tool.** `EXHALE_transit.py` reads the O I column and, since P3,
  `output/OI_levels.txt`. Above the molecular layer nothing changes; inside it,
  the transit tool would now see less O I. Since transit integrations are
  dominated by radii well above the base this is expected to be immaterial, and
  gate G9 measures it rather than asserting it.

### 4.3 The hydrogen budget

`bsp_nH` (`species_table.f90:103`) is the single table that says how many H
nuclei a species carries; `element_nucleus_counts` and `element_ratio_HeH` both
read it, and `binary_element_diffusion` counts nuclei through it. A2 adds rows:
OH carries 1 H, H2O carries 2, CO carries 0. Adding a species without adding its
`bsp_nH` row is the failure mode this table exists to prevent, and it is silent —
which is why gate G5 measures the H budget as well as the O budget.

### 4.4 Mass, particles, electrons

`bsp_mass` is documented as mirroring the literals the code uses today
(`H2 = 2.0`, `H3+ = 3.0`, `HeH+ = 5.0`), so that a future metadata-driven
rewrite stays byte-identical. The new rows must be consistent with **where the
oxygen mass is counted**, and the code counts a metal element's mass as
`melem_A` = 15.999 for oxygen (`species_table.f90:225`). So:

```
bsp_mass(OH)  = 1.0 + melem_A(iel_O)      = 16.999
bsp_mass(H2O) = 2.0 + melem_A(iel_O)      = 17.999
bsp_mass(CO)  = melem_A(iel_C) + melem_A(iel_O) = 28.010
```

and, critically, **the oxygen bound into OH, H2O and CO must be removed from the
metal mass sum**, or `calc_rho` counts it twice. That sum is
`comp_mass_per_H()`/`mass_per_H_nucleus_without_He()`, which today is
`sum(melem_ab*melem_A)` under `eos_include_metals` — a *fixed* metal/H
assumption. This is the same class of defect as the He 2^3S double count that
P2 found and section 75 of `Update_EXHALE.md` fixed with the
`bsp_is_excited_level` flag; the lesson from that incident is that the fix
belongs in the table, not in each summing site.

All three species are neutral: `bsp_charge = 0`, so `calc_ne` is unchanged. Each
counts as one gas particle in `calc_ntot`, as every other species does.

### 4.5 Handoff route and option route: who owns which key

`docs/phase_e_flux_closure_design.md` section 2.4 established the rule and
`input_read.f90:1371-1386` implements it: with a lower-atmosphere profile in
use, the `base.inp` keys of the EOS-boundary, elemental-reservoir and
boundary-constraint categories are **refused**, by name, with the reason
printed, and the pair must carry one matching `solution_id`. A2 is a third
producer of the same quantities and must join the same rule rather than invent
one.

| quantity | handoff route owns | A2 route owns | with both on |
|---|---|---|---|
| base `T`, `r`, `p` | `base.inp` `T_base`/`r_base`/`p_base`, or the profile | nothing | unchanged; A2 has no opinion |
| elemental reservoirs H, He, O, C | `HeH_base`, `<El>_H_base`, `metals.inp`, or the profile | nothing — A2 *consumes* `melem_ab` | unchanged |
| the H2/H partition at the base (`q_H2_base`) | `base.inp` `q_H2_base`, or the profile | **A2 computes it** | **conflict** |
| the O partition among O/OH/H2O/CO at the base | nothing today (A1c species keys are metadata) | **A2 computes it** | no conflict |
| `K_zz` | `Kzz_base` / `He_Kzz` / the profile's `Kzz` column | nothing — A2 reads `kzz_cell` | unchanged |

There is exactly **one** genuine conflict, and it is `q_H2_base`. The plan's own
verdict is that `q_H2_base` is an *EOS anchor*: it enters
`composition::h2_mixing_ratio_base` -> `h2_bound_fraction` -> `comp_ntot_bc`,
i.e. the base particle count and pressure normalization, and it does **not** pin
H2 at the boundary or seed the profile (`composition.f90:249-278`). So the
conflict is narrower than it looks: the handoff sets the base *particle count*
and A2 sets the base *partition*, and with A2 on those two would be built from
different H2 fractions.

### 4.6 Refusals

Following the `refuse_scalar_key` convention exactly — name the key, name its
category, name the other owner, name the fix, `error stop 1`:

1. **A2 on and `q_H2_base` present**: refuse, unless the run is deliberately
   pinning the EOS anchor (see D7). The message names `q_H2_base`, says A2
   computes the base partition, and offers the two resolutions: delete the key,
   or turn A2 off.
2. **A2 on and a lower-atmosphere profile in use**: *accept*, and take the
   profile's molecular mixing ratios as the base element reservoirs — the
   profile owns the region below the match and A2 owns the region above it, so
   they are not two owners of one quantity. But refuse if the profile's
   `solution_id` pairing fails, exactly as today.
3. **A2 on and `Molecular chemistry` off**: refuse. A2 has no H2 to act on.
4. **A2 on and helium absent**: refuse, inheriting the existing
   `thereis_mol .and. .not. thereis_He` refusal (`input_read.f90:948-951`).
5. **A2 on and oxygen absent** (`melem_ab(iel_O) <= 0`): refuse, naming
   `metals.inp` / `O_H_base` as the two places to set it.
6. **A2 on and `K_zz = 0` everywhere**: warn (not refuse), section 3.5.

---

## 5. Keys, files and provenance

### 5.1 Keys

One key in `input.inp`, in the keyword-extension block, matched by `lbl_match`
like every other:

```
Oxygen chemistry: True
```

Named for the physics, placed beside `Molecular chemistry` and `Stellar LW
flux`, default off. It must be tested **before** any key of which it is a
prefix, and it must be added to the label whitelist (`input_read.f90:59-62`)
or it will not be seen.

The bands B1, B3 and B4 need a flux, in the same shape as the two that exist:

```
Stellar FUV flux [erg/cm2/s]: <F>          # 912-2304 A at the planet
```

one number for the whole interval, split among B1/B3/B4 by a stated spectral
shape (flat `F_lambda` proposed, with the sensitivity written at the code site
as `lyman_werner.f90` writes its own). `Stellar LW flux` (912-1110 A, H2 only)
and `Stellar Lya flux` (B2) already exist and are reused unchanged; section 2.6
states why the overlap between the first of those and B1 is not double counting
of photons, and gate G4 is what keeps it from becoming double counting of energy.

Whether one key with an assumed shape is enough, or the three intervals want
three keys, or the whole thing wants an FUV SED file, is decision D1 — it fixes
a user-visible key and, once regression cases carry it, is expensive to rename.

**Not proposed: a file-presence switch.** `metals.inp`, `opacity.inp` and
`base.inp` switch physics by existing, but each of them carries data. A2's
switch carries no data, so it belongs in `input.inp` with the other flags. If
the reaction set later becomes user-selectable, that is when a file appears.

### 5.2 Outputs

| file | content | when |
|---|---|---|
| `output/Ion_species.txt` (`_adv`) | two new schema-2 columns, `OH H2O CO`, appended after `H2 H2p H3p HeHp` | A2 on. The header is generated from the species table (`write_output.f90:82-90`), so the columns are self-describing and `examples/exhale_io.py` picks them up without a change |
| `output/Oxygen_chemistry.txt` | `r[Rp] T[K] n_O n_OH n_H2O n_CO n_O1D x_H2 tau_chem[s] tau_adv[s] Da` per cell, plus the photolysis rates of each band | A2 on, equilibrium state only. The `tau_chem/tau_adv` profile is a **required** output, not a diagnostic: the plan makes the timescale gate part of A2's definition, and a run outside the validity range must say so in its own output |
| `output/FUV_bands.txt` | `r[Rp] N_H2O N_OH N_CO tau_B1 tau_B2 tau_B3 j_H2O j_OH heat_FUV` | A2 on with a positive band flux. Modeled on `output/Lyman_Werner.txt` (`write_output.f90:112-138`), which is the record of how deep a band penetrates |
| `output/Heating_breakdown` | one row for each photolysis channel | as `h_lw` is carried today (`util_ion_eq.f90:1168`) |

The `_adv` files are the known gap: `post_process_adv.f90` is molecule-free by
construction. Section 6, gate G7, takes the plan's wording literally — the
limitation is either lifted or **made loud** — and proposes making it loud: the
`_adv` writer prints a header line stating that the molecular species are not
advection-corrected, and the transit tool refuses to use `_adv` columns inside
the molecular layer. Lifting it is a separate piece of work and belongs in its
own design.

### 5.3 `EXHALE_resolved.out`

`write_resolved_config` (`write_setup_report.f90:430-520`) is the machine-readable
record the transit tool and the budget checkers read instead of re-parsing
`input.inp`, and it already carries `he_diffusion` for exactly the reason A2
needs: so a checker knows which quantity is a column invariant and which is a
solved profile. A2 adds:

```
oxygen_chemistry          T|F
oxygen_reaction_set       <name of the transcribed set, e.g. a2_v1>
fuv_band_LW_flux          <erg/cm2/s>
fuv_band_Lya_flux         <erg/cm2/s>
fuv_band_1110_2000_flux   <erg/cm2/s>
oxygen_base_partition     computed|handoff|equilibrium_fit
```

The last line is the provenance the P0 phase was about: which of the three
possible sources actually set the base partition of this run.

### 5.4 IC and restart

`Ion_species_IC.txt` is read by label under schema 2
(`load_IC.f90:131`), so new columns are additive: an old restart file simply
lacks them and the loader's `col_present` logic handles it. What the design must
state is what happens then — proposed: an A2 run restarted from a pre-A2 file
seeds OH/H2O/CO from the chemical-equilibrium partition of the loaded `(p, T)`
and the loaded oxygen total, prints that it did so, and does not silently start
from zero. Zero is a valid root of the water cycle and hybrd1 is already known
to be bistable from a zero molecular seed
(`ionization_equilibrium.f90:718-726`); the same trap applies here. Gate G6 is
the round trip.

---

## 6. Gates

### 6.1 The seven P4 gates, made testable

The P4 block lists them; this section turns each into something a script can
answer.

| id | P4 gate | test |
|---|---|---|
| G1 | zero-oxygen limit reproduces the current molecular runs exactly | `Oxygen chemistry: True` with `melem_ab(iel_O) = 0` is refused (section 4.6), so the limit is tested the other way: the four existing molecular regression cases (`mol_base_handoff`, `mol_metals`, `mol_lyman_werner`, `mol_diffusion`) run **without** the key and are byte-identical |
| G2 | the atomic limit (OH, H2O -> 0) reproduces the current O I/II/III result | a run with A2 on at a temperature high enough that the water family is empty reproduces the metals-only O I/II/III profile above the front to the solver tolerance; and `wasp_full` / `wasp_he23off` (atomic, metals-on) stay byte-identical |
| G3 | *(new, from section 2.5)* the hot limit is chemical equilibrium | with thermodynamic reversal, the A2 partition at high `p`, high `T` and zero photon flux reproduces `q_h2_equilibrium(p,T)` to a stated tolerance. This is the gate that makes the near-cancelling pairs trustworthy |
| G4 | absorbed energy matches when the FUV opacity is widened | band by band, `integral of 4 pi r^2 (photons absorbed) x (energy per dissociation)` equals the deposited `heat_FUV` summed over the grid, to round-off |
| G5 | element budgets close in every cell | `src/utils/element_budget.py` closes H, He, O and C with the new carriers counted, at its default tolerance, on an A2 case. The existing tolerance is 1e-8 and the measured closure on the P2 case was <= 1.8e-12 over 504 cells |
| G6 | IC write/read round-trips H and O | write `Ion_species_IC.txt` from an A2 run, restart, and reproduce the state; plus the pre-A2 restart path of section 5.4 |
| G7 | the `_adv` limitation is lifted or made loud | section 5.2: made loud, and the transit tool refuses `_adv` metal columns inside the molecular layer |
| G8 | the regression matrix passes with the new physics default off | `make check` 7/7 byte-identical throughout the series, and `run_fcheck.sh` CLEAN |

Two more that follow from this design rather than from P4:

| id | test |
|---|---|
| G9 | the O I redefinition (section 4.2) changes cooling, charge exchange and the transit depth only inside the molecular layer. Measured, not asserted: the O I 1302 band depth of the P3 forward model on an A2 run against the same run with A2 off, stated as a number |
| G10 | the species transport conserves elements: summing the transported carriers with their `bsp_nH`/`bsp_nO` weights reproduces the element flux the existing operator measures (`output/element_flux_profile.txt`), face by face, to round-off |

### 6.2 The A/B gate and its tolerance

This is P4 condition (b) and the first gate the decision names: an A/B against
the Photochem handoff on the HD 189733 b base, run as a regression so that the
two owners of the oxygen physics cannot drift apart silently.

**What is compared.** The H2 fraction at the EXHALE base level, from an A2 run
of HD 189733 b, against the Photochem arm C' solution already stored in
`vulcan_work/pc_compare_p1/hd189_toa1e-2/pc_ncho_solution.pkl`.

**On which quantity — and this matters.** `q_H2` as EXHALE defines it is a
*volume mixing ratio of the total gas*, `n_H2/(n_H2 + n_H + n_He)`
(`composition.f90:249-262`), so it depends on He/H. Two facts make that a trap
for this gate:

- the HD 189733 b input file uses `He/H number ratio: 0.083333333`, for which a
  fully molecular gas gives `q_H2 = 0.5/(0.5 + 0.08333) = 0.857`;
- the Photochem C' solution at 1 microbar has `H2 = 0.7081`, `H = 0.1401`,
  `He = 0.1509` by volume, i.e. `He/H = 0.0970` — a *different* He/H, whose
  fully molecular limit is 0.838.

The two numbers item (H) quotes side by side (EXHALE 0.861 against a
photochemical 0.75) are therefore not directly comparable, and 0.861 sits above
the fully molecular limit of its own run's He/H, which means the quoted pair
needs re-measurement before it can be a target. **The gate is therefore stated
on the H-nucleus bound fraction**

```
x_H2 = 2 n_H2 / n_H,nuclei
```

which is independent of He/H, and which the code already prints as the
`x_H2[2nH2/nH]` column of `output/Lyman_Werner.txt`
(`write_output.f90:124`). From the measured Photochem C' solution,
`x_H2 = 2 x 0.7081 / (2 x 0.7081 + 0.1401) = 0.910`.

**Tolerance, and its justification.** P1 measured, at matched network and
matched domain, the spread between two mature photochemistry codes on this exact
planet and level: **code arm 1.70x in `q_H`**, against a network arm of 3.95x
(and 1.08x / 1.02x respectively on HD 209458 b's 2331 K base). A newly written
Fortran network with a deliberately reduced reaction set cannot reasonably be
held to a tighter standard than two full codes running the same network on the
same profile. So:

- **Pass**: the atomic H-nucleus fraction agrees with the reference within a
  factor of 1.70. At `x_H2 = 0.910` the reference atomic fraction is
  `1 - x_H2 = 0.090`, so the pass band is `0.053 <= 1 - x_H2 <= 0.153`, i.e.
  **`0.847 <= x_H2 <= 0.947`**.

  P1 measured its 1.70 on `q_H`, a volume mixing ratio, not on `1 - x_H2`. The
  two are related exactly, at fixed He/H and fixed elemental hydrogen, by
  `q_H / (1 - x_H2) = 1 / [ (1 - x_H2/2) + He/H ]`, so the ratio between two arms
  differs between the two measures by at most about 6% over the pass band at
  `He/H = 0.0833`. Transferring the factor is therefore safe, and the arithmetic
  is written here so that it is checked rather than assumed.
- **Fail outright**: outside the network factor 3.95x
  (`0.023 <= 1 - x_H2 <= 0.355`). Landing outside that band means A2 is a worse
  approximation than simply choosing a different published network, which is the
  condition under which the option should not ship on by default at all.
- **The comparison pins He/H in both arms** and re-runs the Photochem arm at the
  EXHALE run's He/H if they differ, rather than comparing across it.

**Why this is not circular.** It is only a test if A2 is allowed to get the
answer wrong — which is why section 3.7 proposes a zero-flux boundary condition
on the *partition* while the element reservoirs stay Dirichlet. If the base H2
fraction were imposed from the same handoff the gate compares against, the gate
would measure nothing.

### 6.3 Regression case

One new case, added at the end of the series and snapshotted then, following the
convention of `lower_profile` and `mol_diffusion`: **`backup/regression/mol_oxygen`**,
HD 189733 b, `Molecular chemistry: True`, `Oxygen chemistry: True`, a `K_zz`
from the handoff, metals on, with a step cap in `<case>/maxsteps` if it is a
relaxation snapshot rather than a converged solution. It pins the reaction set,
the transport operator, the FUV bands and the base partition in one run. The
matrix paragraph of the project `CLAUDE.md`, the table in `README_HOWTO.md` and
`docs/input_schema.md` section 2b are updated beside it.

---

## 7. Milestones

Each milestone has a gate that can fail, and none of them is "the code
compiles".

| # | content | gate |
|---|---|---|
| **M1** | **DONE — `docs/a2_reaction_audit.md`.** Every reaction of section 2.4 traced to a published source, its validity range recorded, the conservation skeleton checked at reaction level (elements and charge balance in each reaction), and an explicit verdict on O8-O12. The thermodynamic data for section 2.5, if D2 says so. No code beyond a new `mol_rates`-style module of coefficients | the audit table exists, every row cites a publication read (not a database entry recalled), and the network conserves H, O and C reaction by reaction. A standalone driver reproduces one published rate curve per reaction against its source figure or table |
| **M2** | **DONE — local kinetics only.** The A2 species solved as a local steady state, with the FUV bands, inside the existing cell-by-cell solve. No transport. This is deliberately the state P4 says is *wrong* — it is built because it isolates the chemistry from the transport for debugging | G1, G2, G3, G4 pass. The HD 189733 b base partition is measured and reported, and it is expected to *miss* the A/B target — that miss is the measurement that motivates M3 |
| **M3** | **DONE -- vertical transport.** The species transport-chemistry solve of section 3, the Dirichlet/zero-flux boundary of 3.7, the write-back and the element bookkeeping of section 4 | G5, G10 pass; the marching loop and the steady Picard loop both converge on a case that converges today; `tau_chem/tau_adv` is an output |
| **M4** | **The A/B gate.** HD 189733 b run against the stored Photochem C' arm, He/H pinned, on `x_H2` | section 6.2's pass band. This is the milestone that decides whether A2 ships |
| **M5** | **Integration and regression.** `_adv` made loud, transit-tool guard, `EXHALE_resolved.out` provenance, IC round trip, the new regression case | G6, G7, G8, G9 pass; `make check` byte-identical with the key absent; goldens refreshed only at the end of the series and the refresh reported |
| **M6** | **Documentation.** `docs/input_schema.md` sections 2b and 4, `README.md`, `README_HOWTO.md`, the user manual's lower-atmosphere section, `docs/Update_EXHALE.md`, and the validity-range statement of section 2.7 written into the module header | the manual states what A2 computes, what it assumes, and where it is not valid, in the same terms as section 2.7 |

### M2 result, 2026-08-30

**Implemented and measured. The four gates M2 owns (G1-G4) pass, and G5 and G6
pass with them; the HD 189733 b base partition misses the A/B band, which is
what the table above says M2 is expected to do.**

What was built. `oxygen_rates.f90` is now in `SRC`; atomic carbon was added to
its thermodynamic table for the CO reservoir, and the Shomate evaluation
temperature is clamped to the tabulated range (unclamped, the quartic kept CO
fully associated at 9e5 K and held the whole carbon inventory out of the wind's
coolants). `water_photolysis.f90` is new and carries the four bands. The two
carrier rows and their coefficients live in `System_HeH_mol.f90`
(`set_oxygen_coeffs`, `oxygen_carrier_rows`, `set_oxygen_turnover_rates`,
`excited_oxygen_density`) and are wired into `System_HeH_mol_metals`, which is
the only system the option can reach -- oxygen is a metal element, so a run
without metals is refused. `oxygen_row_base()` sits between the molecular block
and the metals, and `metal_row_base()` shifts by two above it. The species table
gained `isp_OH/H2O/CO`, `bsp_nO` and `bsp_nC`; `n_species` is 40.

Two departures from what section 2 and section 6 assumed, both forced by
measurement and both stated at the code site:

- **The photolysis rate is the cell MEAN, not the face value.** Section 2.6 says
  to use the same column integration every other absorber uses. Done that way,
  the column sum of the absorbed photons fell 30% below the exact
  `N_b(1 - exp(-tau))` in the Ly-alpha band, whose optical depth crosses unity
  inside one cell -- so gate G4's "to round-off" was unreachable. Building the
  rate from the star-ward face depth and the cell's own depth reduces to the
  same rectangle rule when the cell is thin and reproduces the closed form at
  any spacing; G4 then closes at 2.3e-8 to 9.8e-6.
- **The heating charges each absorbed photon the band's own mean energy
  `<hv>`, not the cross-section-weighted `<E>`.** With one transmission for the
  whole band, the saturated limit absorbs every photon, and `<E>` would then
  deposit more energy than the band carries -- 13.5% more on B4. The price is a
  2.7-13.5% low bias in the optically thin limit, which is inside the band
  average's own spectral-shape error.

**Section 6.2's overlap check fires**, and that is a finding rather than a
failed test. With the band fluxes taken from one self-consistently split
spectrum, the Lyman-Werner and B1 treatments -- independent beams over
912-1110 A -- absorb 19.8 out of a B1 flux of 12.0 on the hot Uranus and 1480.9
out of 1248.5 on HD 189733 b. The energy of the overlap is over-counted by up to
60% of the band. The run says so; making one absorber compete with the other
there is not in M2.

**The measurement that motivates M3.** `x_H2 = 0.773` at the HD 189733 b base
against 0.910 for the stored photochemical arm and a pass band of 0.847-0.947.
It misses low, inside the outright-failure band but outside the pass band. The
run is a bounded 12000-step relaxation from cold whose base sits at 1410 K
rather than 864 K, so this is a measurement of the chemistry and not the A/B
gate of M4.

**And a second measurement M3 and (G) have to carry: the option removes
coolants it does not replace.** Where the oxygen is in H2O and CO the [O I]
fine-structure and C I/C II line cooling correctly stop, and the code has no
H2O or CO infrared bands. Measured on the hot-Uranus gate, the layer above the
base cools 15x more slowly with the option on. Every oxygen-chemistry run now
prints how much of the oxygen and carbon is molecular, where, and that those
nuclei no longer cool. Section 9 predicted this; it arrives through the
composition rather than through the temperature.

**A second defect, found by `run_fcheck.sh` and invisible at -O3.**
`x_root_best` in `ionization_equilibrium.f90` -- the thread-private array that
holds the best admissible root of a cell -- was dimensioned
`n_x_max = 8 + 2*n_melem = 28`, the largest N_eq before the oxygen carriers
existed. With every option on N_eq is 30, so the copy of the root into it ran
two elements past the end of a stack array. The optimized build did not notice;
the bounds-checked one stopped on it immediately. `n_x_max` is now
`10 + 2*n_melem`. This is the reason `run_fcheck.sh` was run against an
oxygen-chemistry case and not only against the off-path case its script uses.

**One defect found and fixed in passing, outside the new code.** The hydrogen
budget of `ionization_equilibrium.f90` counts the H nuclei bound in molecules in
three places -- the residual, `ionization_fractions_physical` and
`clamp_fractions_to_element_budget` -- and the carriers had to enter all three.
With the last two left out, a clamped root put more H into the molecules than
the cell had, the extraction floored the remainder at zero, and every element
missed its reservoir by 4.4e-4 at the worst cell. Consistent, the same run
closes at 1.0e-12.

### M3 result, 2026-08-30

**Implemented and measured. G4 now closes in its strongest form, G5 and G6 pass
with the transport on, the option is byte-identical when off, and the A/B
measurement M3 exists to make says that TRANSPORT IS NOT THE MISSING TERM on
this configuration -- with the reason, which is that EXHALE's own base cell has
neither of the two transport terms the P4 table was measured with.**

What was built. `src/modules/lower_atmosphere/diffusive_photochemistry.f90`
(D3, D10): H2, OH, H2O and CO transported by an implicit backward-Euler
diffusion-advection step solved together with the chemistry, block-tridiagonal
in space with 4x4 blocks (D5), zero diffusive flux and no imposed partition at
the base (D6), `kzz_cell` read and no eddy coefficient of its own (section 3.5).
The key is `Oxygen transport`, default True whenever the option is on; `False`
restores M2's local steady state. `P_i - L_i` comes from calling `mol_heh_rows`
row 4 and `oxygen_carrier_rows` at the trial densities rather than from
rewriting them, so the transport and the local solve cannot drift apart; the
Jacobian is a forward difference of the same calls. The local solve is then
handed the transported partition: rows 4, `iox` and `iox+1` become the value
instead of a balance, and the remaining rows are solved against it. Section 2.6's
FUV bands were also rebuilt so that the 912-1110 A interval has one incident
flux and one beam for its three absorbers; that is `Update_EXHALE.md` section
109 and it is what makes the G4 ledger close.

**H2+, H3+ and HeH+ are not transported**, and the criterion is section 2.3's
own: measured over the cells that hold H2, the three together carry at most
**1.2e-9** of the H nuclei, so leaving them local is a closure statement at that
level rather than a leak -- the H-nucleus budget closes against the element
total, not against a sum of transported species. All four transported species
are neutral, so eq. (2)'s ambipolar term is absent rather than neglected. The
ion-neutral polarization coefficient is left out of `D_i` on the same kind of
measurement: over those cells the electron fraction per H nucleus is at most
1.0e-7 (3.1e-11 at the base cell).

**Three departures from what sections 3 and 4 assumed, each forced by a
measurement.**

- **The element headroom is the element TOTAL, not the neutral stage.**
  Transport moves nuclei and the next ionization sweep re-solves which stage
  each sits in; capping the carriers at the neutral stage clamped 150-240 of 503
  cells by up to 2% every step on the HD 209458 b gate. Each cell's element
  totals are now held exactly across the step and the remainder is shared over
  the element's other species in proportion to what they already held.

- **Decision D4 does not survive transport unmodified: CO needed a thermal
  ceiling.** In the local solve "chemically inert" meant "at its own
  `CO <-> C + O` equilibrium", which dissociates CO above about 4000 K.
  Transported, inert means INDESTRUCTIBLE, and it was: on HD 189733 b the wind
  carried CO to 1.67 R_p and 2e4 K where it held **55% of the oxygen**, against
  1.2e-13 for the same run without transport. That switched off the O I and C II
  cooling of the whole wind, and it broke the carrier Newton, which reached its
  cap with a relative residual of 1 in the cell where the runaway had emptied
  the free oxygen. The transported CO is now capped at the chemical equilibrium
  of its own (n, T): inactive through the molecular layer, exact far above the
  turnover, and an over-suppression of quenched CO in the narrow interval around
  3000-4000 K. The audited set has no CO rate to do better with, and the code
  says so at the site.

- **The finite-difference step of the chemistry Jacobian needs a floor tied to
  the element.** Perturbing a carrier that sits at 1e-30 of its element by 1e-6
  of itself returns round-off, because every source term is built from densities
  of order the element total. Floored at 1e-12 of the element, the Newton
  converges in 4-16 iterations to 2e-13 where it had been hitting its cap of 30.

**The A/B measurement, on HD 189733 b with the band fluxes of the same stellar
spectrum the reference arm was run with (LW 600.1, B1 648.4, B2 1.4957e4,
B3 1374.1, B4 4.6295e4 erg cm^-2 s^-1), 12000-step relaxation snapshots from
cold, none of them a converged wind:**

| run | base T | `x_H2` | O partition at the base (O / OH / H2O / CO) |
|---|---|---|---|
| `Oxygen transport: False` (the M2 limit) | 1402.7 K | **0.769** | 0.003 / 0.006 / 0.443 / 0.549 |
| transport, `K_zz = 0` | 1214.3 K | **0.608** | 0.026 / 0.025 / 0.478 / 0.471 |
| transport, `K_zz = 1e9` | 1191.3 K | **0.581** | 0.032 / 0.028 / 0.475 / 0.465 |

against 0.910 for the stored photochemical arm, a pass band of
`0.847 <= x_H2 <= 0.947` and an outright-failure band of
`0.645 <= x_H2 <= 0.977`. The local-kinetics run reproduces M2's 0.773 (0.769
here, the 0.4% being the re-split FUV bands). **Transport moves the partition
AWAY from the target and out of the outright-failure band.**

**Why, measured rather than argued.**

1. **EXHALE's base cell has no advection at all, and its diffusion is 3-4
   decades slower than its chemistry.** `v = 0` exactly at the base cell by the
   lower boundary condition, and in the first cells above it the flow time
   `r/|v|` is 1e8-2e8 s against `tau_chem(H2) = 178-552 s`. With `K_zz = 0` the
   diffusive time of the base cell is 8.0e5 s, i.e. 4500 chemical times; with
   `K_zz = 1e9` it is 634 s, i.e. 3.6. So in EXHALE's own structure the base
   partition is a LOCAL quantity, and the section 3.1 table -- 0.20 to 1.67 --
   was `H/v` on the photochemical model's structure, not on this one. Transport
   cannot be the term that closes the gap where transport does not act.

2. **What transport does change is CO, and it changes it the wrong way.** CO has
   no chemistry, so an inert species' mixing ratio is set by the column and not
   by the local equilibrium: it falls from 0.549 of the oxygen (all of the
   carbon, the equilibrium value) to 0.471, and the freed oxygen goes to the
   water family -- free O rises by a factor 9 and OH by a factor 5 at the base.
   More OH is more H2 destruction through O1. The mechanism is decision D4's,
   not the transport operator's: an inert reservoir that is transported stops
   being an equilibrium partition of the element.

3. **The base is 330-540 K hotter than the reference's 864 K, and the run's own
   H2 budget says the partition there is not oxygen-dominated.** The net budget
   written to `output/Oxygen_chemistry.txt` has the thermal channel at +2.38 and
   the oxygen cycle at -1.38 of the net at the base cell of the `K_zz = 0` run:
   two large terms of opposite sign nearly cancelling, which is the hot-base
   regime section 2.7 puts outside A2's validity range, not the 96-99.6%
   oxygen-dominated regime the gate's reference was computed in. The base is hot
   because the option removes the [O I] and C II coolants without replacing
   them, which is item (G).

4. **Nothing here is converged, and the slow carriers are not at their
   transported steady state either.** `du` ends at 51 and `Mdot` is NaN in all
   three runs. The physical time the 12000 steps cover at the base is
   12000 x 0.589 s = 7.1e3 s, which is 0.9% of the base diffusive time -- so CO,
   which relaxes only by transport, is still near its initial condition.
   Reaching the transported steady state needs `relax_photochemical_composition`
   in the steady Picard loop, and that loop is entered only after a converged
   JFNK solve, which no oxygen-chemistry configuration in the tree reaches.

**What that leaves for M4.** The A/B gate is not answerable from a bounded
marching run: it needs a converged wind whose base sits near 864 K, and that
needs item (G) -- the H2O and CO infrared bands -- before the oxygen chemistry
is asked what partition it produces. The transport operator is in place, its
conservation and its solver are measured, and the term that is missing is
thermal, not advective.

**Gates.** G1/G8: `make check` **7/7 byte-identical** with the key absent, no
golden refreshed. G4: the band ledger closes to <= 3.4e-15 with the H2 pump
photons inside the sum, and the shared-beam block replaces the overlap warning.
G5: `src/utils/element_budget.py` closes H, He, C, N and O at **<= 2.1e-12** on
all three HD 189733 b runs and at 1.0e-12 on the HD 209458 b example, with the
transport on. G6: a restart of an A2 transport run reproduces its state (the
largest change over one marching step is 3.4e-2 and is the hydro's own, since
every species moves by the same factor), and a restart from a state written
before the option existed takes the seeding branch and closes its budgets at
8.7e-13. G10: element conservation is exact by construction rather than
measured -- the operator holds each cell's element totals across the step and
transports speciation, which is what section 3.3's constraint reduces to when
the element operator itself is not transporting that element.

**Not done here.** The steady Picard loop with the carrier relaxation is wired
(`relax_photochemical_composition`, damped with the element relaxation under the
same `omega`) but could not be exercised end to end: it runs only after
`solve_steady_jfnk` returns `info = 0`, and no oxygen-chemistry configuration in
the tree converges that far. M5's regression case, `_adv` guard and transit
guard are not added.

**M4 status update (2026-08-31, `docs/Update_EXHALE.md` section 111).** The
first converged A2 wind now exists: the du descent of the marching run turned
out to be monotonic (the "stall" of the earlier attempts was the
12000-14000-step stretch of a ~24000-step approach), and a JFNK launched from
its du-minimum state (single-stage PLM, `Solver: Newton 0.55`) converges with
`info = 0`, exercising the Picard carrier relaxation end to end for the first
time (composition drift 2.2e-3 -> 5.4e-4 over 8 passes). The gate quantity of
the converged wind is `x_H2 = 0.591` at a base ghost of 888.1 K -- below the
pass band 0.847-0.947 and below the outright-failure edge 0.645, so **the
A/B gate is not passed**, with two caveats that block a final verdict: the
default two-stage recipe destroys this solution at the PLM -> WENO3 switch
(TO_BE_DONE item (S)), and the first end-to-end Picard exercise exposed a
carrier double-count that breaks the converged state's C and O budgets by
+91%/+49% near the base (TO_BE_DONE item (T)). M4 is to be re-measured after
(T) is fixed.

---

### (G) unblocked for M4, 2026-08-31

M3's closing sentence -- **M4's A/B gate is blocked on item (G), not on M3** --
is answered on the coolant side. `Molecular IR bands` (Update_EXHALE.md section
110, `src/modules/lower_atmosphere/molecular_infrared_cooling.f90`) replaces the
[O I] and C I/C II line cooling the option switches off with the coolants that
actually carry that layer in a real H2 atmosphere: the H2O and CO
vibration-rotation bands (HITEMP, via the Photochem correlated-k coefficients)
and the H2 quadrupole plus magnetic dipole line spectrum (Roueff et al. 2019),
each as the NET exchange with the same diluted `B_nu(T0)` the `Base IR field`
closure supplies. Default off, and used together with `Base IR field`.

Measured on an HD 189733 b configuration rebuilt from the keys section 109
records -- two 12000-step relaxations from cold, `Oxygen transport: True`,
`K_zz = 1e9`, the same five band fluxes, differing only in the new key:

| | bands off | **bands on** |
|---|---|---|
| base T | 760.2 K | **829.1 K** |
| `x_H2` | 0.330 | **0.483** |
| O partition at the base (O / OH / H2O / CO) | 0.161 / 0.048 / 0.369 / 0.423 | 0.004 / 0.002 / 0.584 / 0.410 |

**Both gate quantities move toward their targets** -- the base by +69 K toward
the reference's 864 K, `x_H2` by +46% toward 0.910 -- and the mechanism is
visible in the partition: free atomic oxygen falls from 0.161 of the element to
0.004 and the water rises from 0.369 to 0.584, starving the O1 channel that
destroys H2. At r = 1.005 the H2O and CO bands carry 4.8e-5 and 2.7e-5
erg cm^-3 s^-1 of a total 7.5e-5, five orders of magnitude above the H3+ and
fine-structure channels of the same layer.

Two caveats, and they are why this is not M4 passed:

- **Neither run is converged**, and section 6.2's gate is a statement about a
  converged wind. `Mdot` is NaN in the bands-off arm and the bands-on arm's mass
  flux still varies by a factor 2.5 over the wind. `x_H2 = 0.483` is outside the
  pass band `0.847-0.947` and below the outright-failure band.
- **The bands-off arm is not the M3 table's.** Section 109 records the band
  fluxes and the protocol but not the rest of the input file, and the run
  directories no longer exist; rebuilt from the recorded keys the arm lands at
  760.2 K / 0.330 where M3's `K_zz = 1e9` row has 1191.3 K / 0.581. Only the
  difference between the two columns above is a measurement of the coolant
  change.

A converged A2 wind was attempted for this measurement and does not exist yet:
both arms were restarted with the converged-molecular recipe and marched a
further 13000 steps, and `du` plateaus at 1.10-1.60 (bands off) and 1.08-1.35
(bands on) instead of descending to the 5e-2 hand-off. **The non-convergence
survives the coolant being supplied**, which locates it in the base
hydrodynamics rather than in the thermal budget -- so M4 needs a convergence
fix of its own, not more physics.

**What M4 now needs is a converged HD 189733 b A2 wind**, which is M4's own
work: the thermal term it was waiting on exists, is validated against three
independent published sources, and moves the gate quantities the right way.
G5 still closes (`src/utils/element_budget.py`: 1.8e-14 with the bands on), and
the G4 ledger gained an infrared side -- the column-integrated emitted, absorbed
and net power of the three bands against the incident flux `W sigma T0^4` that
bounds what an optically thin layer can absorb, measured 73x below its bound on
the bands-on arm (absorbed 7.5e5 against 5.4e7 erg cm^-2 s^-1).

---

M1's deliverables are `docs/a2_reaction_audit.md`,
`src/modules/lower_atmosphere/oxygen_rates.f90` and the standalone driver
`src/tests/a2_m1/` (the coefficient module is in `SRC` since M2). M2 added the
local kinetics and the FUV bands; M3 added the transport module
`src/modules/lower_atmosphere/diffusive_photochemistry.f90`, the
`Oxygen transport` key and the shared-beam FUV field. What M3 measured moves the
remaining work: **M4's A/B gate is blocked on item (G), not on M3**, because the
base cell has no advection in EXHALE's own structure and the base sits 330-540 K
above the reference's 864 K until the H2O and CO infrared bands exist. M5's
regression case, `_adv` guard and transit guard are untouched.

---

## 8. Decisions to confirm before implementation

**All ten decided 2026-08-30 (user).** D1 departs from the proposal: **three
separate FUV keys, one per band interval** (`Stellar FUV B1/B3/B4` naming per
the table below; `Stellar Lya flux` supplies B2 as proposed) -- the
single-key flat-shape form is weakest exactly on the M dwarfs this code now
runs, so the split is adopted from the start. D2-D10 are adopted **as
proposed**: thermodynamic reversal with the NIST-Shomate data of
`zahnle_earth.yaml` (D2); a new transport module (D3); CO as an unreactive
reservoir, no CO2 (D4); fully implicit chemistry-transport coupling (D5);
elemental Dirichlet base with zero-flux carrier partition (D6); `q_H2_base`
refused while the option is on (D7); no sulfur (D8); default off stays after
the gates pass (D9); physics-based module naming (D10). Implementation
starts at M1.

Each is expensive to reverse because it fixes a user-visible key, a data file, a
physical convention, or what the first gate measures.

| # | decision | proposed | alternative, and what reversing costs |
|---|---|---|---|
| **D1** | **Key names and the FUV band split.** `Oxygen chemistry: True`; one new `Stellar FUV flux [erg/cm2/s]` covering 912-2304 A, split among B1/B3/B4 by a flat `F_lambda` shape, with `Stellar Lya flux` supplying B2 | as stated | three separate keys, one per interval, or a full FUV SED file. Renaming later breaks every stored input file and the regression case. The M-dwarf case (a line-dominated FUV) is where the single-key form is weakest |
| **D2** | **Reverse rates, and which thermodynamic table.** | thermodynamic reversal for O1, O2, O8, O9, O10 (section 2.5) — both reference codes do this and neither carries an explicit reverse | explicit reverse rates, as `mol_rates.f90` does today. Within the proposed route the sub-choice is **Shomate from `zahnle_earth.yaml`'s NIST-sourced `species:` block** against **NASA-9 from `VULCAN/thermo/NASA9/`** (Burcat). Both are about 30 numbers for the five species; the decision fixes which reference file the transcription is checked against |
| **D3** | **New transport module, or generalize `binary_element_diffusion`.** | a new module for the species system; the element operator untouched | one operator for both. Reversing is a rewrite of whichever was built, and the element operator's two-component closure is not obviously extensible |
| **D4** | **Carbon: CO as an inert reservoir, and where its abundance comes from.** | carry CO, frozen chemically, its abundance from the handoff or from a chemical-equilibrium C/O partition at the base; CO2 dropped | carry CO reacting (`CO + OH -> CO2 + H`), or drop CO and rescale the oxygen reservoir. The measurement says CO holds 45% of the oxygen, so dropping it is not free |
| **D5** | **Chemistry-transport coupling: fully coupled implicit, or operator split.** | fully coupled implicit (section 3.6) | operator split. This decides the size and shape of the solve and is expensive to swap after M3 |
| **D6** | **The base boundary condition on the partition.** | Dirichlet on the element reservoirs, **zero-flux on the partition among each element's carriers** | Dirichlet on the partition too, taken from `q_H2_base`. This decides whether the A/B gate of 6.2 measures anything |
| **D7** | **`q_H2_base` when A2 is on.** | refuse the key (section 4.6) | accept it as the EOS anchor only, with A2's computed partition used for the chemistry and the handoff value used for `comp_ntot_bc` — defensible, since they are different quantities, but it is two numbers for one gas and the plan's single-source rule argues against it |
| **D8** | **Sulfur.** | out of the minimal set | in. P1.5 measured sulfur moving `q_H2` by 17% at 1e-6 bar on HD 209458 b, which is larger than several terms that are in. It is excluded on the HD 189733 b budget (0.36% of the net), and the two planets disagree |
| **D9** | **Whether the default stays off after the gate passes.** | yes, off, indefinitely — every new default stays off until a golden refresh is deliberate | on by default once M4 passes. Turning it on moves every molecular golden and changes the base of every molecular run |
| **D10** | **Module and routine names.** proposed: `oxygen_hydrogen_rates.f90` (coefficients, mirroring `mol_rates.f90`), `water_photolysis.f90` (the FUV bands, mirroring `lyman_werner.f90`), and for the transport solve a name for the physics it computes — `diffusive_photochemistry.f90` | as stated | the naming rule forbids role names, and `species_diffusion.f90` is taken by history (it was deleted on 2026-08-25 and replaced by `binary_element_diffusion.f90`), so reusing it would be actively confusing |

---

## 9. What this does not close

**Item (G), the missing continuum IR coupling, still ranks first.** Section 4.1
of the plan of record ranks it above this work, on a measured effect: a
converged molecular layer radiates itself down to 190-400 K against
`T_eq ~ 1100-1400 K`, and on the HD 209458 b molecular example the collapse
reaches the wind at `Mdot -0.34 dex`. A2 is the composition of that layer; (G)
is its energy. They feed on each other — the equilibrium H2 fraction depends on
`T`, and the coolant inventory depends on composition — so **A2 alone does not
make the molecular layer a prediction**, and the design must not be presented as
if it does. A2's radiation work is FUV photolysis; (G)'s is thermal IR. Different
spectral regions, different code, one region of the atmosphere.

**The measured payoff on the wind is small and is an upper bound.** P1.6
measured, on the one converged planet, at most **0.035 dex in `Mdot`** across
four chemistry arms, and the same measurement at a level where the arms differ
by 0.1% still spreads by 0.025 dex — so 0.035 dex is at the level of the
configuration's own JFNK-to-JFNK reproducibility. A2 is justified by
self-containedness and by turning an unbounded uncertainty into a modeling
result, which is what the P4 decision says; it is **not** justified by a wind
effect, and no report of it should claim one.

**A2 does not remove the caveat that makes it necessary.** Section 3.1's own
table says `tau_adv` falls below the measured `tau_chem` beyond
`r/R_p ~ 1.17` on HD 209458 b, i.e. inside the transit-relevant radii. Adding
transport moves the boundary of validity outward; it does not remove it. The
`tau_chem/tau_adv` profile is an output (section 5.2) precisely so that every
A2 run states where its own assumptions stop holding.

**One P3 item remains open and is not touched here.** The published
O I 1302 depth is not reproduced (model 0.359% against 10.8 +/- 4.5%); it is
recorded on `TO_BE_DONE.md` item (H) and is not A2's subject. The second item
that stood here, whether `O2+ + H0 -> O+ + H+` should be on by default, was
closed on 2026-08-30: it is on by default now
(`Update_EXHALE.md` section 107).
