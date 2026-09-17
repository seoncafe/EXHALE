# LHS 1140 b, item L7g: is the molecular base at 1023 K physics or a heat-ledger defect?

Item L7g of `docs/PLAN_20260913_lhs_stationary.md`, opened by item L7f and
approved 2026-09-15. The question: the molecular case
`molecular_scalar_gj1132_kzz1e9/HeH2.13` carries a base at 1 microbar and
3e13 cm^-3 whose temperature is 1023 K, where the ATOMIC case of the same
planet, the same He/H and the same boundary data carries 418 K and `base.inp`
states 226 K. Is that a result, or does the code add the same energy twice?

**NOTE ADDED 2026-09-17 (items L7g inventory, L7g model, L7h).** The question
this memo answers, whether the ledger double counts, is answered here and the
answer stands; but the ledger it measured was wrong in five other ways, all of
them corrected since, so the numbers below are the OLD ledger's. (i) The
dissociative recombinations of H2+ and HeH+ deposited their whole enthalpy as
gas kinetic energy where one H atom leaves in n = 2, and now return 0.749 and
1.554 eV instead of 10.948 and 11.753. (ii) The H3+ recombination deposited its
whole 9.250 eV promptly where 2.4795 eV is internal energy of the H2 fragment.
(iii) The three-body association applied the M = H2 rate coefficient to the
TOTAL heavy-particle density and is now collider-resolved. (iv) The quench
fraction's radiative side used the v = 1 total decay rate, and the validity
statement quoted here, that the collider density is far above `n_cr` of every
level, was checked against `n_cr(v=1)` only; the all-level maximum over the
302-level ladder is 6.56 times that rate. (v) The collisional side was two
Hollenbach & McKee (1979) fits with helium absent, and is now one published
quantum calculation per collider. MEASURED effect: the molecular chemical heat
at the base of this very case falls 21.6 per cent, and the certified state this
memo measured on is no longer a root of the corrected equations. Record:
`docs/lhs1140b_stationary_L7g_model_20260916.md`,
`..._L7g_inventory_20260916.md`.

Nothing in `build/`, `EXHALE.x`, `backup/regression/` or any catalogue case
directory was written by this item: every binary is `EXHALE_L7f.x` built into
`build_L7f/`, and every run is under `LHS1140b/models/.L7f/`.

## 1. Verdict

**The reaction-heat ledger does NOT double count. Every cycle closes exactly
on the ionization potential of the ion that drives it, delivered once, and
paid for by the photon that made that ion.** The 1023 K is therefore not a
bookkeeping error; it is what this energy equation gives for this chemistry,
this K_zz and -- as item L21 established while this item was running -- the
base boundary condition as it stood.

**But the arithmetic closing is a NECESSARY condition and not a sufficient
one: it says the joules are accounted for, not who receives them.** A
radiative channel emits most of its enthalpy as light and the gas gets only
the rest, and a table difference cannot know that. Audited channel by
channel (section 3.1), exactly one of the seventeen is radiative -- R23,
H2 + He+ -> H2+ + He + hv -- and it was depositing its ~153 nm photon,
8.1035 of its 9.1615 eV, as heat. That is corrected. In this configuration
it moves the base heating by 8e-05 of itself, which is reported and is not
the reason for the change.

The construction is what makes the double count impossible. Every species
carries ONE formation-plus-excitation energy, measured from a single
reference state (`species_formation_energy`, `molecular_reaction_heat.f90`:
"every element a neutral, ground-state, free atom at rest, the free electron
at zero"), and every reaction heat is a difference of those, formed nowhere
else. A photoionization deposits `hv - IP` and no more -- the integrand is
`photoelectron_share(e_th, e) = 1 - e_th/e` (`util_ion_eq.f90`) -- and the
`IP` it did not deposit is exactly the formation energy the new ion now
carries. So an ion's binding energy is put into the reservoir once and taken
out once, by whichever reaction destroys it.

**Two published enthalpies verify the table rather than the table verifying
itself.** Schauer et al. (1989), J. Chem. Phys. 91, 4593, Sec. I print the
two He+ + H2 channels as `dH = -6.51 eV` (dissociative) and `dH = -9.16 eV`
(radiative). This code's ledger, built only from ionization potentials and
dissociation energies it holds for other reasons, gives **+6.5109 eV for R17
and +9.1615 eV for R23** (`physics_probe`, `reaction_energy_R17` and
`reaction_energy_R23`). Neither number was fitted to the paper.

## 2. The reaction heats, as the ledger forms them

MEASURED (`physics_probe/species_formation_energy_table`, which also asserts
that each is the formation-table difference to 1e-12 and that each reaction
conserves nuclei and charge). Positive is released to the gas.

| # | reaction | q [eV] | what it is |
|---|---|---|---|
| R5 | H2+ + e -> H + H | +10.948 | dissociative recombination: the H2+ binding relative to 2H + e |
| R6 | H3+ + e -> H2 + H | +9.250 | |
| R7 | H3+ + e -> 3H | +4.771 | |
| R8 | H2+ + H2 -> H3+ + H | +1.698 | |
| R9 | H2+ + H -> H+ + H2 | +1.827 | |
| R10 | H+ + H2 -> H2+ + H | **-1.827** | exact reverse of R9 |
| R11 | H3+ + H -> H2+ + H2 | **-1.698** | exact reverse of R8 |
| R12 | H2 + M -> H + H + M | **-4.478** | the bond, endothermic |
| R13 | H+ + H2 + M -> H3+ + M | +4.349 | |
| R14 | H2 + e -> H + H + e | **-4.478** | same bond |
| R15 | H + H + M -> H2 + M | **+4.478** | same bond, released |
| R16 | HeH+ + e -> He + H | +11.753 | |
| R17 | He+ + H2 -> H+ + H + He | **+6.511** | published dH = -6.51 eV |
| R18 | HeH+ + H2 -> H3+ + He | +2.504 | |
| R19 | HeH+ + H -> H2+ + He | +0.805 | |
| new | H2+ + He -> HeH+ + H | **-0.805** | exact reverse of R19 |
| R23 | H2 + He+ -> H2+ + He | **+9.161** | published dH = -9.16 eV |

Every forward/reverse pair is an exact negative (R9/R10, R8/R11, R12/R15,
R19/new), which is the signature of one table and not seventeen numbers.

## 3. The cycles, closed

The question is not whether one reaction's heat is right but whether a CLOSED
loop -- H2 destroyed, the fragments recombined, the ion neutralized, H2
reformed -- adds more than the photon supplied. Three loops carry the base.

**(A) The helium loop, 39 percent of the H2 loss at cell 1.**

| step | q [eV] |
|---|---|
| He + hv -> He+ + e | hv - 24.587 to the gas; 24.587 stored in He+ |
| He+ + H2 -> H+ + H + He (R17) | +6.511 |
| H+ + e -> H + hv' (radiative recombination) | the electron's kinetic energy only; 13.598 leaves as a photon |
| H + H + M -> H2 + M (R15) | +4.478 |
| **sum to the gas** | **hv - 24.587 + 6.511 + 4.478 = hv - 13.598** |

and 13.598 eV is exactly what the recombination photon carried away.
6.511 + 4.478 = 10.989 = IP(He) - IP(H) to five digits: the loop delivers the
charge-transfer energy once and nothing else.

**(B) The proton loop, 49 percent of the H2 loss at cell 1.**

| step | q [eV] |
|---|---|
| H + hv -> H+ + e | hv - 13.598 to the gas; 13.598 stored in H+ |
| H+ + H2 + M -> H3+ + M (R13) | +4.349 |
| H3+ + e -> H2 + H (R6) | +9.250 |
| **sum** | **hv - 13.598 + 4.349 + 9.250 = hv** |

4.349 + 9.250 = 13.599 = IP(H). The whole photon ends as heat, and that is
the physics of this loop rather than an error: the proton's binding energy is
returned to the gas by the DISSOCIATIVE RECOMBINATION of H3+ instead of being
radiated by a radiative recombination of H+. A molecular layer under an
ionized wind converts recombination energy into heat.

**(C) The Lyman-Werner loop, 3 percent.** The photon pays the 4.478 eV bond
and the fragments keep the excess (the code says so at the code site: "the
bond energy is paid by the photon, not by the gas"); R15 then returns the
4.478. Sum to the gas: the whole photon, again by construction and not twice.

**So the ledger is closed, and `heat_mol_chem` is not an extra source: it is
where the ionization energy of H+ and He+ is deposited when the molecular
channels neutralize them.**

**But closing the arithmetic is a necessary condition and not a sufficient
one, and one channel failed the second test.** The cycles above account for
every joule; they do not by themselves say WHO RECEIVES it. A reaction whose
products include a photon puts most of its enthalpy into light, and the gas
gets only the rest -- and the table difference, which is the reaction
enthalpy, knows nothing about that. Of the seventeen reactions exactly one
is radiative, and it was depositing its photon as heat.

### 3.1 The radiative-channel audit

| # | reaction | kind | photon? | what the ledger deposits |
|---|---|---|---|---|
| R5 | H2+ + e -> H + H | dissociative recombination | no | all of q, correctly: the energy is the fragments' |
| R6, R7 | H3+ + e -> H2 + H, 3H | dissociative recombination | no | all of q |
| R16 | HeH+ + e -> He + H | dissociative recombination | no | all of q |
| R8, R9, R17, R18, R19, H2+ + He | ion-neutral rearrangements | no | all of q |
| R10, R11, R12, R14 | endothermic | no | q < 0, a sink |
| R13 | H+ + H2 + M -> H3+ + M | three-body association | no | all of q, to the third body |
| R15 | H + H + M -> H2 + M | three-body association | infrared, and already held back | q times `h2_vibrational_heat_fraction` |
| **R23** | **H2 + He+ -> H2+ + He + hv** | **radiative charge transfer** | **YES, ~153 nm** | **was all of q; CORRECTED 2026-09-15** |

The correction is the published mechanism, not a guess. Boehringer & Arnold
(1986) p. 1461, on Hopper's account of this channel (READ, verbatim):

> Hopper reconsidered the reaction dynamics in more detail and suggested
> that also a reaction channel leading to H2+ should be possible. This
> process involves a radiative transition from the first formed excited
> state of the collision complex (He+.H2) to the ground state (He.H2+) which
> then decays into the products (1a). The wavelength of the emitted photon
> should be about 153 nm and H2+ should preferentially be produced in a
> vibrationally excited state (v = 2).

153 nm is **8.1035 eV** of the reaction's 9.1615 eV, and at 1530 A it is
longward of the Lyman-Werner bands (912-1110 A) and of the Lyman continuum,
so nothing where it is made absorbs it. The ledger now deposits
q(R23) - 8.1035 eV = **1.0579 eV** and lets the photon go. The reaction
ENTHALPY is untouched -- q(ir_R23) is still the formation-table difference,
which is what `physics_probe` asserts -- and what is subtracted is the part
of it that leaves as light. Stated at the code site with its range: no
re-absorption is modelled (this code carries no transfer for a 1530 A line)
and the v = 2 the ion is born in, about 0.55 eV, is deposited rather than
followed (no v-resolved H2+), so the term still OVERSTATES by at most that
0.55 eV and understates nothing.

**What it moves here: almost nothing, and that is reported rather than used
as a reason.** MEASURED on the certified state (`.L7f/g_newbc`), with
k23 = 7.2e-15:

| | R23 rate [cm^-3 s^-1] | heat with the whole enthalpy | photon removed | share of the cell's heating |
|---|---|---|---|---|
| cell 1 (808 K) | 7.087e+00 | 1.040e-10 | 9.201e-11 | **8.2e-05** |
| cell 25 (1527 K) | 1.652e-01 | 2.424e-12 | 2.144e-12 | 2.8e-05 |
| cell 100 (1535 K) | 2.060e-01 | 3.024e-12 | 2.674e-12 | 3.5e-05 |

so the base heating changes in its fifth digit. The channel is small here
because R23 is 7.2e-15 against R17's 3.8e-12 at these temperatures. It is
corrected because it is wrong, not because it matters in this case.

## 4. What that leaves

MEASURED on the state item L7f delivered (`.L7f/heatpp_noir`, a `Do only PP`
pass on `.L7e/fixed5`), the base layer's heating at cell 1 is 97.7 percent
`heat_mol_chem`, 5.205e-07 of 5.327e-07 erg cm^-3 s^-1, against a total
cooling of 1.101e-08 -- a factor 48 -- and column-integrated over the whole
domain the heating is 7.69 erg cm^-2 s^-1 (plane) with `heat_mol_chem` 2.88
of it, 37 percent. The R15 three-body association alone runs at
6.92e+04 cm^-3 s^-1 at cell 1 and releases 4.48 eV each, 4.95e-07 erg cm^-3
s^-1, essentially the whole of that term.

**What is NOT settled by section 3**: whether the code is right to thermalize
nearly all of that 4.48 eV. `h2_vibrational_heat_fraction` exists to hold
back the share that leaves in the H2 infrared quadrupole lines, and at the
base density it evaluates to very nearly 1, while the H2 infrared bands
radiate 1.17e-09 -- 0.2 percent of the deposit. Whether a real base radiates
more of the association energy than that is a question about
`h2_vibrational_heat_fraction` and the band opacity, not about the ledger,
and item L7f's IR-coolant test already showed that switching both infrared
channels on raises the cooling by only 11 percent.

## 5. Why 226 K, 418 K and 1023 K are three different things

- **226 K is an input.** It is `base.inp`'s `T_base`, which fixes the
  reservoir entropy of the boundary condition. Item L21 measured that the old
  boundary blended it almost entirely out and has corrected that; the numbers
  of section 4 are on the OLD boundary and section 6 re-measures them.
- **418 K is a solved value** -- the atomic case's own base temperature under
  the same old boundary, with no molecular chemistry and therefore no
  `heat_mol_chem` at all.
- **1023 K is a solved value with the molecular energy source switched on**,
  and section 3 says that source is the ionization energy of H+ and He+
  arriving as heat through the molecular recombination chains.

**And the photochemical column cannot disagree, because it does not solve for
temperature.** READ from the clone in this workspace: `photochem` takes T from
a prescribed pressure-temperature-eddy profile and maps it onto its grid
(`photochem_evoatmosphere_utils.f90`, `map_press_temp_edd`, and
`self%var%temperature = T_new` is the only assignment to it); there is no
energy equation and no heating rate in its state. Its 1 microbar temperature
is an input, so the `q_H2_base` it hands over was computed at whatever T that
profile states and carries no statement about the reaction heat at all.

## 6. On the corrected base boundary, and the case certifies

Item L21 corrected the base boundary while this item was running: the face
was being put on the reversal branch by the SIGN OF CELL 1's velocity, which
is the odd-even artefact of the P44 note, so 97 percent of the face state was
the interior isentrope and `T_base` was blended out. Every number of
sections 4 and 5 above was measured on the old boundary. Re-measured on the
new one, with the delivered L7f network, the same `fixed5` seed and the same
continuation recipe (`.L7f/g_newbc`, 20 outer passes requested):

| | old boundary | **new boundary** |
|---|---|---|
| ghost T at r = 1.0000 | 873.0 K | **226.02 K** -- `base.inp`'s `T_base` to four digits |
| ghost rho [mH/cm^3] | 2.903e+13 | 1.153e+14 |
| cell 1 T | 1023.0 K | **808.3 K** |
| cell 1 x2 | 0.3315 | 0.3553 |
| cell 25 T | 1544.5 K | 1527.3 K |
| cell-1 mass row | 1.68e-07 | **4.92e-09** |
| momentum row | 3.0e-11 | 1.16e-12 |
| energy row | 7.4e-07 | **9.96e-09** |
| carrier row | 2.11e-08 (40 passes) | 6.99e-06 (11 passes) |
| verdict | NOT CERTIFIED, cell-1 mass row | **ACCEPTED at outer pass 11** |

**The boundary now delivers the reservoir temperature**: the ghost at
r = 1.0000 carries 226.02 K at p = 1.0001 p0, which is `base.inp`'s stated
(226 K, 1 microbar) to four digits, where the old boundary carried 873 K
there. And the state CERTIFIES: "outer pass 11: ACCEPTED -- every active
equation of this state is within its own tolerance", the cell-1 mass row
having fallen by a factor 34 and the energy row by 74.

**What that does and does not change for L7g.** The base GHOST is now the
lower atmosphere's own 226 K gas, as it should be. The first physical CELL is
still hot -- 808 K, not 418 -- and the layer above it is unmoved (1527 K at
cell 25 against 1544). So the molecular reaction heat still raises the first
resolved cell by a factor 3.6 over the boundary it is fed from, and the
verdict of sections 1 to 3 stands unchanged: that heat is the ionization
energy of H+ and He+ arriving through the molecular recombination chains,
counted once. The open question is unchanged too and is now the only one:
whether `h2_vibrational_heat_fraction` should hold back more of the 4.478 eV
that R15 releases.

**What it changes for item L7f.** The cell-1 mass row that L7f left open --
the base-face composition step -- was in large part the old boundary's
reversal-branch artefact, not the composition step alone: the same
composition step is still there (ghost x2 = 0.99998 against cell 1 at
0.3553) and the row now certifies anyway.

## 7. Tests and regression

This item changed no code. Its measurements were made with the binary item
L7f delivered and then, for section 6, with the tree build that carries item
L21's base-boundary correction on top of it (`EXHALE.x` and `EXHALE_L7f.x`
are the same file, md5 `db87b88d1ce53facf1d61084fa535ca5`). The suites and
the regression that belong to those changes are reported in
`docs/lhs1140b_stationary_L7f_20260915.md` sections 5.1 and 5.3.

What this item adds to the record is the `Do only PP` pass that produced the
channel decomposition (`.L7f/heatpp_noir` and `.L7f/heatpp_ir`) and the two
continuations (`.L7f/HeH2.13_ir`, the infrared-coolant test of item L7f, and
`.L7f/g_newbc`, the corrected-boundary state of section 6).

## 8. What is left, named

**Whether `h2_vibrational_heat_fraction` should hold back more of the
4.478 eV that R15 releases.** It evaluates to very nearly 1 at the base
density, so the code thermalizes essentially all of the association energy
there, while the H2 infrared bands radiate 1.17e-09 erg cm^-3 s^-1 -- 0.2
percent of the 5.205e-07 deposited. A newly formed H2 is born vibrationally
hot and a real atmosphere radiates part of that bond energy in the
rovibrational lines before the molecule thermalizes; whether the fraction,
the band opacity or the association rate is what carries the difference is
not settled here. The measurement that bounds it: switching both infrared
channels on raises the total cooling by 11 percent and moves the base
temperature by 0.2 K, so the answer is not simply "the bands were off".

**The code site now says what the fraction is** (item R6 of the review of
2026-09-15, checked and corrected 2026-09-16). The association-heat line of
`molecular_reaction_heat.f90` states that `h2_vibrational_heat_fraction` is
an EFFECTIVE thermalization fraction: Hollenbach & McKee's v = 1 expression
applied to the 4.478 eV that three-body association releases, which is not a
v = 1 excitation -- the nascent molecule is formed high in the vibrational
ladder and cascades -- so it is valid only where the collider density is far
above n_cr of EVERY level of that cascade, which is what makes the outcome
independent of the level the molecule was born in. This layer satisfies it
(n(H2) = 4e9 to 5e13 cm^-3 against n_cr(v=1) = 2e6 to 6e6, and the fraction
is 1 to a part in 1e3); a shallower base or the top of a thinner layer would
not, and there the level-resolved cascade would have to be carried. The item
stays open with that condition written where the number is used.
