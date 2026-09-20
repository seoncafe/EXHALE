# L31: the molecular energy cycles, and the validity record of four approximations

Written 2026-09-17 for item L31 of `docs/PLAN_20260917.md`. Every number
below is labeled MEASURED (this item ran it) or READ (from a source file, a
published PDF or a written state). The item's brief is the plan's section 5,
which follows `docs/session_handoff_20260917_review.md` section 5.

---

## 0. Verdict

**No molecular energy share is counted twice.** Five complete cycles of the
network close to roundoff on the production assembly, including the one the
review asked about by name: the vibrational excitation a nascent H2+ is born
with at R23 is deposited once, at R23, and its dissociative recombination R5
is charged the GROUND-STATE enthalpy of the ion and adds none of it back
(section 2). A deliberate 0.55 eV double count, the size of the v = 2
excitation the R23 comment names, is caught by the cycle rows at
exactly that size (section 6).

**Four prescriptions the code called bounds are approximations, and they now
say so.** The unit n = 2 recipient of R5 and R16 is a LOW-STATE RECIPIENT
APPROXIMATION and not an upper bound on the heat; the R6 internal energy is a
MODEL ESTIMATE from the position of a published peak and not the mean the
ledger needs; the scalar quench fraction is a FIRST-EVENT BRANCHING MODEL and
not a lower bound over the ladder; and argon's three-body coefficient is an
ARGON-BASED ESTIMATE for helium and not a bound. The code text now states
each as what it is, with its range and its uncertainty (section 3).

**The scalar quench fraction stays the default, and the reduced
statistical-equilibrium model says why.** The model solves the 54
rovibrational levels the published H2-H collision data reach, with three
colliders and one explicitly uncertain group for the levels above the data.
At the three representative cells the scalar form's radiated share is 2.8 to
13.2 times the largest the ladder model gives, so the scalar keeps more energy
out of the gas than any member of the bracket does and replacing it could only
add heat. **No default is proposed for change** (section 5).

**The out-of-range extrapolation of the three-body pair is negligible, and
that is now a measurement and not an argument.** On the certified molecular
fiducial the R15 association source at the 76 cells with T outside 77-5000 K
is 2.3e-08 to 3.5e-06 of the other three H2 formation channels there; on the
hottest atomic wind of the catalog it would support an H2 abundance of at most
1.7e-13 of the hydrogen nuclei (section 4).

**The largest molecular heat uncertainty at the base is the helium third
body, not either recipient.** Cohen & Westberg's +/-0.3 in log on the argon
coefficient moves the chemical heat at the molecular depth by +41 and -21 per
cent, against 2e-07 for the recipient and 7e-07 for the quench fraction
(section 5.2). Higher up the ordering reverses: at the dilute upper column
the recipient carries 130 per cent of the nominal heat.

---

## 1. What the ledger is, in one paragraph

`molecular_chemical_heating`
(`src/modules/lower_atmosphere/molecular_reaction_heat.f90`) forms the heat of
the seventeen COLLISIONAL reactions of the H2/He network from ONE table of
species formation energies, so that a cycle returning to its starting species
releases exactly zero by construction. Photon-driven reactions, radiative
recombinations, collisional ionizations, the Penning channels and the H/He
charge exchange are excluded because their energy is already in the ledger
elsewhere. Four reactions hand part of their enthalpy to something other than
the gas: R5 and R16 leave one hydrogen atom in n = 2, R6 leaves its H2
fragment vibrationally hot, R15 puts its whole bond energy into internal
energy, and R23 emits a 153 nm photon.

READ from the production accessors on the build of this item (MEASURED by
running `qtab` against the production modules), in eV:

| quantity | value |
|---|---|
| `IP(H2)` | 15.425927 |
| `IP(He)` | 24.587389 |
| `IP(H)` | 13.598435 |
| `D0(H2)` | 4.478075 |
| `E(H n=2)` = 0.75 `IP(H)` | 10.198826 |
| `E_int(R6)`, midpoint of v = 5 and v = 6 | 2.479545 |
| `E_gamma(R23)`, the 153 nm photon | 8.103542 |
| `Q(R5)` `H2+ + e -> H + H` | 10.947852 |
| `Q(R6)` `H3+ + e -> H2 + H` | 9.249565 |
| `Q(R8)` `H2+ + H2 -> H3+ + H` | 1.698287 |
| `Q(R12)` `H2 + M -> H + H + M` | -4.478075 |
| `Q(R15)` `H + H + M -> H2 + M` | 4.478075 |
| `Q(R16)` `HeH+ + e -> He + H` | 11.753281 |
| `Q(H2+ + He -> HeH+ + H)` | -0.805429 |
| `Q(R23)` `H2 + He+ -> H2+ + He` | 9.161462 |
| `IP(H2) - D0(H2)` | 10.947852 |
| `IP(He) - D0(H2)` | 20.109314 |

---

## 2. The energy cycles, and the proof that nothing is counted twice

Each row is a CLOSED cycle: it returns every species to where it started, so
what it liberates is fixed by its input alone and by nothing inside it. The
gas column is MEASURED from the production assembly
`molecular_chemical_heating`, run on a cell that carries that step and divided
by the step's own event rate; the other columns are the production accessors
that own them. Energies in eV per event of the cycle.

### 2.1 Cycle 1, the direct route

| step | module | gas | H(n=2) level | photon | stored H2 internal |
|---|---|---|---|---|---|
| `H2 + hv -> H2+ + e` | `ionization_equilibrium` (photoelectron only) | 0 | 0 | 0 | 0 |
| R5 `H2+ + e -> H + H` | `molecular_reaction_heat` | 0.749026 | 10.198826 | 0 | 0 |
| **sum** | | **0.749026** | **10.198826** | 0 | 0 |

Input to the chemical reservoir: `IP(H2)` = 15.425927, of which `D0(H2)` =
4.478075 goes into breaking the bond. Liberated: 10.947852. Sum of the row:
10.947852. Closes.

The photon's own excess `hv - IP(H2)` is the photoelectron's and is deposited
by `PH_heat_HHe`, not by this ledger; it is outside the cycle by construction
and cannot be double counted here.

### 2.2 Cycle 2, through H3+

| step | module | gas | level | photon | stored |
|---|---|---|---|---|---|
| `H2 + hv -> H2+ + e` | `ionization_equilibrium` | 0 | 0 | 0 | 0 |
| R8 `H2+ + H2 -> H3+ + H` | `molecular_reaction_heat` | 1.698287 | 0 | 0 | 0 |
| R6 `H3+ + e -> H2 + H` | `molecular_reaction_heat` | 6.770020 | 0 | 0 | 2.479545 |
| **sum** | | **8.468307** | 0 | 0 | **2.479545** |

Sum 10.947852, equal to cycle 1's. The gas column of R6 is its prompt share
at a cell with no collider, where the thermalized fraction is exactly zero and
the whole internal share is withheld; with colliders present the internal
share moves from the last column to the first through
`h2_vibrational_heat_fraction` and the total is unchanged.

The second H2 of R8 is regenerated by R6, which is what makes this a cycle.

### 2.3 Cycle 3, through HeH+

| step | module | gas | level | photon | stored |
|---|---|---|---|---|---|
| `H2 + hv -> H2+ + e` | `ionization_equilibrium` | 0 | 0 | 0 | 0 |
| `H2+ + He -> HeH+ + H` | `molecular_reaction_heat` | -0.805429 | 0 | 0 | 0 |
| R16 `HeH+ + e -> He + H` | `molecular_reaction_heat` | 1.554455 | 10.198826 | 0 | 0 |
| **sum** | | **0.749026** | **10.198826** | 0 | 0 |

Sum 10.947852. Helium is a catalyst and returns.

### 2.4 Cycle 4, THE NASCENT ION

This is the cycle the review asked to be audited: a nascent ion's vibrational
excitation must not be deposited at formation and again at its dissociative
recombination.

| step | module | gas | level | photon | stored |
|---|---|---|---|---|---|
| `He + hv -> He+ + e` | `ionization_equilibrium` | 0 | 0 | 0 | 0 |
| R23 `H2 + He+ -> H2+ + He + hv` | `molecular_reaction_heat` | 1.057920 | 0 | 8.103542 | 0 |
| R5 `H2+ + e -> H + H` | `molecular_reaction_heat` | 0.749026 | 10.198826 | 0 | 0 |
| **sum** | | **1.806946** | **10.198826** | **8.103542** | 0 |

Input `IP(He)` = 24.587389 less `D0(H2)` = 4.478075 gives 20.109314. The row
sums to 20.109314. Closes.

**Where the nascent excitation is, and where it is not.** Boehringer & Arnold
(1986) quote Hopper's mechanism for R23 with an emitted photon "about 153 nm"
and H2+ left "preferentially ... in a vibrationally excited state (v = 2)",
about 0.55 eV (READ, quoted in the code at the R23 subtraction). The ledger
deposits that 0.55 eV ONCE, inside the 1.057920 eV R23 leaves after its
photon, because the formation table's `eps(H2+)` is the GROUND-STATE value and
the whole difference `Q(R23) - E_gamma` therefore contains it. R5 is charged
`Q(R5)`, which is again the ground-state enthalpy and adds nothing back. A
ledger that deposited the excitation at formation and again at recombination
would exceed the input by 0.55 eV, and section 6 shows that it does and that
the row catches it.

The APPROXIMATION here is the timing, not the count: the code deposits the
0.55 eV promptly rather than following a v-resolved H2+, so it overstates the
prompt deposit by at most that amount if the ion radiates the excitation
instead of relaxing collisionally. That statement is in the code at the R23
subtraction and is unchanged by this item.

### 2.5 Cycle 5, the three-body pair

| step | module | gas | level | photon | radiated in the IR lines |
|---|---|---|---|---|---|
| R15 `H + H + M -> H2 + M` | `molecular_reaction_heat` | `4.478075 f` | 0 | 0 | `4.478075 (1-f)` |
| R12 `H2 + M -> H + H + M` | `molecular_reaction_heat` | -4.478075 | 0 | 0 | 0 |
| **sum** | | | | | **0** |

A closed cycle liberates exactly zero, whatever the thermalized fraction `f`
is. This is the property the single formation-energy table guarantees and
that a reaction-by-reaction enthalpy list cannot be trusted to have.

### 2.6 The routes agree with each other

Cycles 1, 2 and 3 are the same net reaction, `H2 + hv -> H + H`, through three
different intermediates. MEASURED, they agree to 0.0 eV (the row
`the_four_routes_liberate_the_same_energy` reads exactly zero). Cycle 4 is the
same net reaction driven by a helium photoionization, and closes to
2.20e-13 eV, which is the roundoff of the two subtractions the measurement
makes.

### 2.7 What was checked in the files this item does not own

The audit had to follow the excitation out of the ledger. Nothing was edited
there; what was READ:

- `excited_hydrogen.f90` lines 237 and 485 call
  `dissociative_recombination_n2_source`, which is the SAME `rk_R5_H2p_dr` and
  `rk_R16_HeHp_dr` at the same densities the ledger charged, so one event
  cannot be counted at two rates. The n = 2 population it feeds disposes of
  the 10.199 eV through the Balmer photoionization (whose deposit is the
  photoelectron's excess over the n = 2 threshold, paid by the photon), through
  Lyman-alpha emission, and through the optional collisional de-excitation term
  `Hdx_arr`. The last of these is the only one that returns the 10.199 eV to
  the gas, and it is gated on `incl_deexc_heat`. No double count.
- `ionization_equilibrium.f90` line 1262 and `util_ion_eq.f90` lines 2113 and
  2385 all pass a ground-singlet helium density into
  `h2_vibrational_heat_fraction`. The comment block of that function said they
  could not; it is corrected (section 3.3).
- `System_HeH_mol.f90` and `diffusive_photochemistry.f90` reach the collider
  sum through one routine, verified bitwise in section 5.3.

**No double count was found in a file this item does not own.**

---

## 3. The four approximations, as the code now states them

### 3.1 The n = 2 recipient of R5 and R16: a low-state recipient approximation

`molecular_reaction_heat.f90`, at the subtraction inside
`molecular_chemical_heating`. What the sources establish (READ, and quoted in
the code): Takagi (2002) gives the n = 2 product "only for the low vibrational
molecular ion and at low collision energies"; Giusti-Suzor et al. (1983) give
the `H(1s) + H(2s)` limit for electron energies below 0.5 eV and H2+ in its
lowest three vibrational states; Guberman (1994) computed ground-vibrational
3HeH+ at electron energies 0.001 to 0.33 eV.

What the code applies it to: an H2+ whose chemical lifetime against R5, R8 and
R9 is about 1e-3 s against a radiative vibrational relaxation of about 1 s, so
it recombines in the distribution it was born in.

**The claim removed.** The comment said a vibrationally hot ion opens n >= 3,
which takes more of the enthalpy out of the kinetic channel, so the share
written is an upper bound on the heat. That does not follow: the released
translational energy is

```
E_kin = Q_ground + E_int(reactants) - E_int(products) ,
```

so the reactant's vibrational energy enters with a plus sign at the same time
as a higher product n subtracts more than 10.199 eV. A bound on the product
principal quantum number alone does not bound `E_kin` when the initial
vibrational energy is free as well.

**What the code now does.** It states the limit, states that its own
population lies outside the condition, and prints a validity line once per run
when the recipients are in force:

```
 (molecular_reaction_heat) VALIDITY: the n = 2 product of R5 and R16 is applied to every event.
 (molecular_reaction_heat) VALIDITY: H2+ lives about  1.0E-03 s against R5/R8/R9 and relaxes vibrationally in about  1.0E+00 s, so it is vibrationally hot.
 (molecular_reaction_heat) VALIDITY: the unit branching is the LOW-v, low-energy limit of Takagi (2002), Giusti-Suzor et al. (1983) and Guberman (1994); it is a low-state recipient approximation and not a bound on the heat.
```

**The excited-hydrogen-off branch** now carries its assumption explicitly, at
`dissociative_recombination_n2_source`: with `use_excited_H` false nothing
consumes the n = 2 source, and the 10.199 eV the ledger has already withheld
is assumed to LEAVE THE CELL AS RADIATION. Two disposals are then not
modelled, collisional de-excitation of the n = 2 atom and reabsorption of the
Lyman-alpha photon, and both would return part of it, so with the level off
the ledger is a LOWER bound on the gas heat of R5 and R16 by at most the full
10.199 eV of each event. With the level on, the disposal is computed.

### 3.2 The R6 internal energy: a model estimate from the peak

`molecular_reaction_heat.f90`, at the assignment of `e_int_R6`. The quantity
the ledger needs is

```
E_int = sum over (v, J) of  P(v, J) E(v, J)
```

with `P` a normalized, applicable product distribution. Neither source
tabulates one: Kokoouline, Greene & Esry (2001) quote the peak of a
DIRECT-pathway calculation, which is not a complete treatment of the indirect
pathways, and Strasser et al. (2001) measure a distribution they describe as
wide and discuss rotational excitation of both the incident ions and the
products. A peak location, and an uncertainty on that location, do not fix the
mean of a broad asymmetric distribution.

**The uncertainty that IS established is the peak's.** Strasser et al. quote
"an uncertainty of about one level"; one level near v = 5 is 0.37 to 0.40 eV,
15 per cent of the 2.479545 eV share. That is what the code now says the
+/-15 per cent is, with the distance from peak to mean stated as
unquantified.

**No figure was digitized.** Neither paper's figure was extracted, so this
memo records no digitization and no population assumptions. If one is ever
digitized, this section is where the digitization procedure and the population
assumptions belong.

### 3.3 The scalar quench fraction: a first-event branching model

`h2_vibrational_relaxation.f90`, at `h2_vibrational_heat_fraction`. Two claims
were removed and one stale paragraph corrected.

**Removed: the lower bound over the whole ladder.** `A_max` bounds every
level's radiative loss from above, but the v = 1 collisional coefficient does
not bound every level's collisional loss from below. A level with `C_u < C_1`
would have a collisional branching fraction below the adopted scalar even
though `A_u <= A_max`. Separately, a first-event branching fraction is not the
total energy fraction recovered through a multistep cascade.

**Removed: the closely-spaced-levels argument.** Lique's (2015) section 3.2
reports a modest increase with initial vibrational state for a fixed
vibrational change in the H-collision transitions he studied. That is
supporting evidence inside that data set. Selection rules, the collider's
identity, rotational redistribution, reactive destruction of the excited
molecule and upward thermal transitions all enter, and none of them is settled
by it.

**Corrected: the stale helium paragraph** (the lines the review names). It
said the three photoelectric and Lyman-Werner call sites are handed no helium
density and cannot state one. They do (READ):
`ionization_equilibrium.f90` line 1262 passes
`he_ground_singlet_density(nhei, nheiTR)` into `f_vib_quench`;
`util_ion_eq.f90` line 2113 passes `nheiS_chem` into the Lyman-Werner
fluorescence channel of the heating assembly; `util_ion_eq.f90` line 2385
passes `nheiS(j)` into the same channel of the band ledger. The optional
argument is an interface property, not a caller's limitation, and the comment
now says so.

**What replaces the bound**: the reduced statistical-equilibrium model of
`src/tests/h2_level_ladder`, section 5.

### 3.4 Argon as helium: an estimate, not a bound

`mol_rates.f90`, at `k3b_H_H_to_H2_monatomic` and in the module header. The
claim removed is `k1(He) <= k1(Ar)` from mass and polarizability. Those
qualitative properties do not order thermally averaged quantum three-body
recombination rates over 77-5000 K: Paolini, Ohlinger & Forrey (2011) treat He
and Ar with different interaction potentials and different resonance and
continuum contributions, discuss the sensitivity to the potential energy
surface, distinguish equilibrium from steady-state populations of the
intermediate complex, and find an additional exchange contribution that
matters for Ar at low temperature and not for He (READ, the review's reading
of the published paper, accepted here).

What the code now states: an ARGON-BASED ESTIMATE FOR HELIUM, range 77-5000 K,
uncertainty +/-0.3 in log `k1(Ar)` (Cohen & Westberg's own statement for
argon), with the distance from argon to helium NOT inside that uncertainty.

**One choice in both directions.** R12 and R15 multiply one collider sum built
from this coefficient, so whatever the helium estimate is, the pair stays an
exact detailed balance collider by collider and the estimate cannot displace
the equilibrium constant. That is verified bitwise in section 5.3 and by the
existing detailed-balance rows of the probe.

If a helium coefficient is ever obtained (the revised handoff's request list,
section 9 item 4), the estimate becomes a value and this section is where the
substitution is recorded.

---

## 4. The out-of-range extrapolation, measured

The argument replaced: "above the H2 front they are extrapolated and the H2
density they multiply has gone." R15 forms H2 out of ATOMIC hydrogen,
`n(H)^2 sum_M k1(M) n_M`, so the disappearance of H2 does not make it vanish.

MEASURED with the production rate coefficients of `mol_rates` on scratch
copies of two written states (`r15_out_of_range_source`, a driver built
against the production modules; both states read from their
`Hydro_ioniz_IC.txt` / `Ion_species_IC.txt` pair, which is the certified one).

### 4.1 The certified molecular fiducial

`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13`, `_IC` pair.

| quantity | value |
|---|---|
| cells with T outside 77-5000 K | 76 of 504 (all above 5000 K) |
| their temperature range | 5007.83 to 5870.73 K |
| largest R15 association source there | 1.6465e-05 cm^-3 s^-1, at r = 1.2350 R_p, T = 5007.8 K |
| the other three H2 formation channels (R9, R6, R11) at that cell | 4.6848e+00 cm^-3 s^-1 |
| `P(R15) / (P(R9) + P(R6) + P(R11))` over the 76 cells | 2.3404e-08 to 3.5145e-06 |
| largest `2 n(H2)_eq / n_H` the R12/R15 pair alone balances at | 3.9053e-11 |

The association source is at most three and a half parts in a million of the
H2 formation there, and the H2 abundance the extrapolated pair would settle at
is 4e-11 of the hydrogen nuclei. **Negligible.**

### 4.2 The hottest atomic wind of the catalog

Picked by MEASURING the maximum temperature of every physical cell of all 70
certified atomic states of `LHS1140b/models/`: the hottest is
`atomic_scalar_gj699_wellmixed/HeH1000`, maximum T = 9517.55 K. (The next
three are `atomic_scalar_gj1132_wellmixed/HeH1000` at 8662.7 K, `HeH100` at
8614.3 K and `HeH10` at 8269.4 K.)

| quantity | value |
|---|---|
| cells with T outside 77-5000 K | 184 of 504 |
| their temperature range | 5013.64 to 9517.55 K |
| largest R15 association source there | 1.3357e-10 cm^-3 s^-1, at r = 1.0789 R_p, T = 5074.9 K |
| smallest | 5.1801e-20 cm^-3 s^-1 |
| largest `2 n(H2)_eq / n_H` | 1.6745e-13 |

The run carries no molecular species at all, so the other three formation
channels are identically zero there and the ratio of section 4.1 has no
denominator. The comparison that stands is against the hydrogen budget: even
with the coefficients extrapolated to 9518 K, the H2 the pair balances at is
2e-13 of the hydrogen nuclei. **Negligible.**

### 4.3 Why the extrapolation cannot move the balance at all

R12 and R15 carry ONE collider sum. Their ratio is therefore the equilibrium
constant `K_eq(T)`, collider by collider, with no `k1` left in it. What the
extrapolation can move is the RATE at which the pair relaxes towards that
balance, never the balance itself. The `n(H2)_eq` column above is
`K_eq n(H)^2` and carries no extrapolated coefficient.

---

## 5. The uncertainty record

### 5.1 The reduced statistical-equilibrium model

`src/tests/h2_level_ladder/` (new suite; `h2_ladder_model.f90` is the model,
`h2_ladder_statistical_equilibrium.f90` the driver, `README.md` the record of
what it is built from). It evaluates the exact form the scalar fraction
approximates,

```
Q = sum over level pairs (u, l) and colliders M of
    [ n_u C_ul(M) - n_l C_lu(M) ] n_M (E_u - E_l) ,
```

on a statistical equilibrium of the levels the published data reach.

**The data** (READ, each from the file its own publication distributes):
54 rovibrational levels of X with v <= 3, the set Lique (2015) distributes
state-to-state H2-H rate coefficients for over 100 to 5000 K; helium from
Jozwiak et al. (2024), both directions on a 20 to 8000 K grid; level energies,
weights and A values from Roueff et al. (2019) table 2, the same file the
production ladder is reduced from. Upward rates are built from the downward
ones by detailed balance with `g = g_I (2J+1)` and are never read.

**The H2 collider is an assumption of the model.** Le Bourlot, Pineau des
Forets & Flower's state-to-state table is not in this tree; the helium matrix
supplies its shape and it is normalized to the published thermal v = 1 -> 0
coefficient the production module carries. At the cells this model is run on
H2 carries 6e-05 of the collider sum, so it cannot carry the answer either
way.

**The one explicitly uncertain part.** Everything above v = 3 is carried as
ONE group, which is where the nascent molecules of R15 (v = 10 to 14) and R6
(v = 5 to 6) are born. It empties into the top resolved level at a collisional
rate spanning a decade about the top level's own collisional total and at a
radiative rate spanning zero to the all-level maximum. The chemical
destruction rate of an H2 molecule spans 1e-09 to 1e-02 s^-1, which brackets
the all-level radiative rate 5.6e-06 s^-1 on both sides. The spread of the
answer over those three is the model's uncertainty.

**The tests** (MEASURED, all GREEN; 32 assertions):

| row | measured |
|---|---|
| all-level `A_max` against the production module | 5.593952e-06 against 5.5939e-06 s^-1, 9.4e-06 relative, the production value the lower because its line list is trimmed |
| Boltzmann departure in a collision-only bath, 808 K | 2.97e-14 |
| the same, 2128 K / 5083 K | 2.92e-15 / 6.12e-15 |
| net collisional heat there, against the scale of the individual terms | 2.63e-16 / 8.82e-17 / 3.33e-16 |
| steady injection: `Q + L + k_destroy U` against the injection | agrees to 3.7e-13 relative |
| one backward-Euler step from an empty ladder, the four channels with the stored term | agrees to 3.4e-16 relative |

### 5.2 The bracket on the scalar fraction

`f` is the share of the internal energy a nascent molecule is born with that
stays in the gas. ESCAPING RADIATION IS THE ONLY EXIT: the net collisional
term hands the energy to the translational pool, and the internal energy a
molecule still holds when the chemistry destroys it goes to the fragments and
is carried by the caloric equation of state, which owns the rovibrational
ladder of the bulk gas. The model's `f` is `1 - dL/P`, with `dL` the
difference in escaping radiation between injecting the molecules into the
group and injecting the SAME molecule rate into the ground level; the system
is linear, so that difference is the response to the internal energy alone and
the ladder's thermal emission cancels exactly.

MEASURED on the three cells of the certified molecular fiducial (section 5.4):

| cell | T [K] | `1 - f` scalar | `1 - f` ladder, bracket | ratio scalar / most radiative ladder |
|---|---|---|---|---|
| molecular depth | 808.32 | 6.793198e-07 | 1.20e-08 to 5.140654e-08 | 13.21 |
| H2 front | 2128.37 | 1.006850e-05 | 2.76e-07 to 2.004984e-06 | 5.02 |
| dilute upper column | 5083.44 | 5.963722e-04 | 6.99e-06 to 2.146044e-04 | 2.78 |

**The scalar form radiates at least as much as any member of the bracket at
all three cells**, so it keeps more energy out of the gas than the ladder
model does and replacing it could only add heat. It is not outside the stated
uncertainty in the direction that would matter. **The default does not
change, and no change is proposed.**

The size of the ratio, 2.8 to 13.2, is the uncertainty of the scalar form on
`1 - f` at these cells. Because `f` is 1 to a part in 1e6 at the molecular
depth and to 6e-04 at the dilute upper column, that uncertainty reaches the
heat only through the terms `E_int(R6)` and `q(R15)` multiply, which the next
section measures.

### 5.3 The uncertainty rows of the 62-assertion probe

`src/tests/physics_probe/molecular_energy_recipients.f90`, extended from 62 to
81 assertions. Five quantities are varied SEPARATELY on the three cells, each
as an exact difference against the production assembly's own output, formed
from the same production rate coefficients and reaction energies so that no
term of the ledger is restated in the driver.

MEASURED, molecular chemical heat in erg cm^-3 s^-1:

| cell | T [K] | nominal heat | recipient, 1 -> 0 | R6 branching, 1 -> 0 | R6 `E_int`, +15% | quench, `f` -> 1 | He efficiency, +0.3 dex | He efficiency, -0.3 dex |
|---|---|---|---|---|---|---|---|---|
| molecular depth | 808.32 | 8.71941e-07 | +1.79695e-13 | -1.09189e-15 | +1.63784e-16 | +5.71502e-13 | +3.61454e-07 | -1.81156e-07 |
| H2 front | 2128.37 | 8.78021e-09 | +4.00747e-10 | -3.97996e-15 | +5.96995e-16 | +3.98305e-15 | -1.75920e-11 | +8.81691e-12 |
| dilute upper column | 5083.44 | 3.15327e-12 | +4.08523e-12 | -1.56921e-21 | +2.35382e-22 | +1.59043e-21 | -2.04225e-15 | +1.02355e-15 |

As fractions of the nominal heat:

| cell | recipient | R6 branching | R6 `E_int` | quench | He efficiency |
|---|---|---|---|---|---|
| molecular depth | 2.1e-07 | 1.3e-09 | 1.9e-10 | 6.6e-07 | **+41% / -21%** |
| H2 front | **+4.6%** | 4.5e-07 | 6.8e-08 | 4.5e-07 | -0.20% / +0.10% |
| dilute upper column | **+130%** | 5.0e-10 | 7.5e-11 | 5.0e-10 | -0.065% / +0.032% |

**The ordering, and it is the useful result of this item.** At the molecular
depth the heat's uncertainty is the HELIUM THIRD BODY and nothing else: the
+/-0.3 dex of the argon coefficient moves the chemical heat by +41 and -21 per
cent, four to six orders above everything else. That is because the base heat
is dominated by the R15 and R12 terms, which the collider sum multiplies, and
because helium carries most of the collider sum there. Higher up the ordering
reverses: at the H2 front and in the dilute upper column the R15 terms have
gone and the recipient of R5 and R16 carries the heat, 4.6 per cent at the
front and more than the nominal heat itself in the dilute column, where the
ledger's molecular heat is 3e-12 erg cm^-3 s^-1 and of no consequence to the
energy budget.

**The two internal shares carry no uncertainty in this layer at all.** The R6
branching and the R6 internal energy reach the gas only through `1 - f`, which
is 6.8e-07 at the molecular depth (MEASURED), so their brackets are 1e-9 and
smaller everywhere. **This is why the R6 peak-versus-mean question, which
section 3.2 leaves open, does not block the molecular catalog re-run**: it
would take a base thin enough for `f` to fall away from one before the
distinction reached a number the run reports.

### 5.4 The three representative cells

READ from
`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13/output/{Hydro_ioniz_IC.txt,
Ion_species_IC.txt}`, the certified `_IC` pair (`certified=T
cert_reason=certified_in_wind`), rows 3, 180 and 301 of the written state,
which are physical cells 1, 178 and 299.

| | molecular depth | H2 front | dilute upper column |
|---|---|---|---|
| row of the written state | 3 | 180 | 301 |
| r [R_p] | 1.000193 | 1.101040 | 1.842260 |
| T [K] | 808.319 | 2128.367 | 5083.438 |
| `2 n(H2) / n_H` | 0.35530 | 0.17366 | 3.6e-06 |
| n(H I) [cm^-3] | 1.868e+12 | 1.564e+10 | 8.162e+07 |
| n(H2) | 5.147e+11 | 1.646e+09 | 1.534e+02 |
| n(He I) | 6.172e+12 | 3.588e+10 | 7.025e+07 |
| n_e | 9.000e+06 | 2.333e+07 | 6.901e+06 |

The dilute upper column is a CHOSEN radius, r = 1.842260 R_p (cell 299), and
not a derived rule: there `2 n(H2)/n_H` is flat and rising outward (3.43e-06
at cell 296 to 3.88e-06 at cell 303, MEASURED by L37), so no threshold on it
selects a cell. The probe's state reader (`EXHALE_PROBE_STATE`, item L37)
takes the cell nearest that radius. The H2 front is defined as the cell
where `2 n(H2)/n_H` has fallen to half its base value; it lands at r = 1.1010 R_p, inside the 1.08 to 1.10 R_p the
revised handoff's regression table quotes for the front's movement.

### 5.5 One collider sum, in three places, bitwise

The chemistry (`System_HeH_mol::h2_third_body_density`, the call
`mol_heh_rows` makes at the composition the row is evaluated at), the carrier
balance (which reaches the same call through `mol_heh_rows` and writes it into
the R15 channel of its row record) and the heat ledger
(`molecular_chemical_heating`) each form the collider sum. MEASURED at the
molecular-depth cell:

| comparison | difference |
|---|---|
| chemistry against `h2_association_collider_density` | 0.0 exactly |
| the carrier row's R15 channel against `mk15 n_third n(H)^2` | 0.0 exactly |
| the heat ledger's deposit on an R12/R15-only cell against the same third body | agrees to 1e-14 relative (8.396784550935822e-07 measured and expected) |

`write_output.f90` forms a fourth reading of the same routine for its
diagnostic; it was READ and not tested, because it enters no solution.

---

## 6. RED before, GREEN after

Every new row was shown to fail before it was shown to pass.

### 6.1 The cycle rows, on an injected double count

A deliberate fault was applied to `molecular_reaction_heat.f90`: after the
n = 2 subtraction, R5's gas share was given an extra 0.55 eV, which is the
nascent H2+ vibrational energy the R23 comment names and is exactly the double
count the audit is for. MEASURED with the fault in place (the file was then
restored from a copy taken before it):

```
FAIL cycle_photoionization_then_R5_closes      measured= 5.500000000000007E-01  tol=1e-09 absolute
FAIL cycle_R23_nascent_H2p_then_R5_closes      measured= 5.500000000002210E-01  tol=1e-09 absolute
FAIL the_four_routes_liberate_the_same_energy  measured= 5.500000000000007E-01  tol=1e-09 absolute
FAIL ledger_R5_deposits_enthalpy_less_n2       measured= 2.957362899065809E-04  reference=1.705232898701601E-04
FAIL ledger_R5_same_without_excited_H          measured= 2.957362899065809E-04  reference=1.705232898701601E-04
```

The cycle rows read the injected surplus to twelve digits. After restoring the
file, all five are GREEN and the cycle residuals are 0.0, 0.0, 0.0,
2.20e-13 and 0.0 eV.

### 6.2 The control key

`EXHALE_REACTION_HEAT_RECIPIENTS=0` restores the pre-correction ledger. With
it set, 19 of the probe's 81 assertions fail, including all five cycle rows,
the three recipient-bracket rows and all three collider-sum rows. Unset, all
81 pass.

### 6.3 The ladder suite

Two faults, each applied and then reverted from a copy:

- **Detailed balance broken** (the degeneracy ratio dropped from the
  Boltzmann factor that builds the upward rates): the Boltzmann departure
  rises from 3e-14 to 6.61, the net collisional heat of a bath from 3e-16 to
  0.55 of the term scale, and the steady partition sum goes negative. Nine
  rows FAIL.
- **The group's collisional disposal counted in both the heat and the
  radiation**: the steady partition sum reads 1.126e-07 against an injection
  of 7.143e-08, the transient sum 3.748e-10 against 2.910e-10, and the
  fraction goes negative at all three cells. Eight rows FAIL.

After reverting, 32 of 32 assertions pass.

---

## 7. What changed, and what did not

**Changed, in the three source files this item owns**: comment text only,
except for one new statement. `molecular_reaction_heat.f90` gains
`reaction_heat_recipient_validity_report`, a module procedure called once from
`molecular_chemical_heating` that writes three lines to standard output on its
first evaluation when the recipients are in force. It reads no state, changes
no number and is called from a serial context (`heating_of_composition` runs
outside every parallel region: `ionization_equilibrium.f90` closes its cell
sweep at line 2760 and calls the heating assembly at line 3045).

**Not changed**: every rate coefficient, every reaction enthalpy, every
default, and every number any run produces. No golden was refreshed and none
needs to be.

**Proposed for change**: nothing. The bracket of section 5.2 does not put the
scalar quench fraction outside the direction that would matter, so the default
stands.

**Recorded for the request list.** A published helium three-body recombination
coefficient for `H + H + He -> H2 + He` over 200 to 5000 K would turn section
3.4's estimate into a value, and section 5.3 shows it is the largest single
uncertainty in the molecular base heat: +41 / -21 per cent at the certified
base. Paolini, Ohlinger & Forrey (2011), Phys. Rev. A 83, 042713, plot a total
over 0-350 K and tabulate nothing; that is the only helium calculation in
hand. This is the revised handoff's request item 4 and this item raises its
priority.
