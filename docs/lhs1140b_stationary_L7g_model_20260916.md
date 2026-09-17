# L7g model: who receives the energy of the molecular network, and how fast the three-body channel runs

Item L7g of `docs/PLAN_20260913_lhs_stationary.md`, the MODEL half, written on
the data inventory `docs/lhs1140b_stationary_L7g_inventory_20260916.md` and on
the Cloudy c25.00 reading `docs/cloudy_h2_model_reference_20260916.md`. The
inventory established what the published record says; this item puts it in the
code.

Every number below is marked MEASURED (computed or run here) or READ (from a
source document or a source file).

---

## 0. Verdict

**Five statements the code was making about the molecular layer were
physically wrong, and all five are corrected. None of them changes a reaction
enthalpy.** What changes is the RECIPIENT of an enthalpy, the RATE of one
channel, and the COEFFICIENTS of the branching that decides whether internal
energy becomes heat.

| # | what was wrong | what it is now | source |
|---|---|---|---|
| 1 | R5 (`H2+ + e -> H + H`) and R16 (`HeH+ + e -> He + H`) deposited their whole enthalpy as gas kinetic energy | one hydrogen atom leaves in n = 2, so the gas gets 0.749 and 1.554 eV and 10.199 eV goes to the H(n=2) population, or leaves as radiation where the run carries no such population | Takagi (2002); Giusti-Suzor, Bardsley & Derkits (1983); Guberman (1994) |
| 2 | R6 (`H3+ + e -> H2 + H`) deposited its whole enthalpy promptly | 2.480 eV of its 9.250 is internal energy of the H2 fragment and branches between the infrared lines and the gas exactly as R15's does | Kokoouline, Greene & Esry (2001); Strasser et al. (2001) |
| 3 | R12 and R15 applied the M = H2 three-body coefficient to the TOTAL heavy-particle density | the collider sum `k1(H2) n(H2) + k1(H) n(H) + k1(Ar) n(He)`, with the pair kept an exact detailed balance collider by collider | Cohen & Westberg (1983), p. 559 |
| 4 | the heat fraction used the v = 1 radiative rate while claiming validity "far above n_cr of EVERY level" | the largest total spontaneous decay rate over the code's own 302-level ladder, reduced from the line list at initialization | Roueff et al. (2019) through the tree's own table, cross-checked against Wolniewicz, Simbotin & Dalgarno (1998) |
| 5 | the collisional side of the fraction was two 1979 fits, and helium was absent from it altogether | three published quantum calculations, one per collider | Lique (2015) for H; Jozwiak et al. (2024) for He; Le Bourlot, Pineau des Forets & Flower (1999) for H2 |

**The largest single effect on this layer is item 3, not the recipients.**
MEASURED on the certified base state, `heat_mol_chem` falls 21.6 per cent at
cell 1 and 9.4 per cent at cell 25, and essentially all of that is the
collider-resolved three-body rate: R15 carries 95 per cent of the channel and
its rate falls 22.2 per cent at cell 1. The recipient corrections move almost
nothing HERE, because the thermalized fraction is 1 to a part in 1e6 in this
layer and because the two carriers that leave an atom in n = 2 are traces; they
matter for a shallower or thinner base, which is what item 4 of the inventory
said they would.

**The largest correction in absolute size is item 5, and it is not a small
one**: the Hollenbach & McKee H2-H2 fit the code carried is 107 times the
scattering calculation at 808 K and 21 times it at 1527 K. In this layer the
collider sum is 99.8 per cent atomic hydrogen, so the fraction does not move;
in a cold H2 base with no atomic hydrogen that fit was claiming fifty times too
much thermalization.

---

## 1. What was implemented, and where

### 1.1 R5 and R16: the recipient is the n = 2 level, not the gas

`src/modules/lower_atmosphere/molecular_reaction_heat.f90`,
`molecular_chemical_heating`.

The heat of one event becomes `q(R) - E(n=2)`, with `E(n=2)` taken from the ONE
entry of the formation-energy table that holds it,
`species_formation_energy(isp_eps_H_n2)` = 10.19883 eV = 0.75 IP(H), the same
energy the Balmer photoionization and the Lyman-alpha channel are charged. The
enthalpy `q(R)` itself is untouched, exactly as the R23 photon correction of the
predecessor item left it: the formation table is not reachable from the site
that subtracts.

MEASURED from the code's own table: q(R5) = 10.9478 eV leaving 0.7490 eV, and
q(R16) = 11.7534 eV leaving 1.5545 eV.

The 10.199 eV is delivered to the H(n=2) population by a new source term in
`src/modules/radiation/excited_hydrogen.f90`:
`n2_rate_matrix` gains a chemical production `S_chem = k5 n(H2+) n_e +
k16 n(HeH+) n_e`, formed by `dissociative_recombination_n2_source` in the module
that owns the subtraction, so the energy that left the gas and the level that
receives it are one event counted once. `n2_populations` and
`excited_hydrogen_level_residual` both pass it, so the residual measures the
balance the closure solved.

**The l-partition is an assumption and is recorded as one.** The source is
split between 2s and 2p by statistical weight, `g2s : g2p = 1 : 3`. Takagi
states the product as "n = 2" without an l, and Guberman's HeH+ C state
likewise dissociates to "an excited n = 2 H atom". Giusti-Suzor et al. DO name
the H(1s) + H(2s) limit, but for H2+ in its lowest three vibrational states at
electron energies below 0.5 eV, and this code's H2+ is vibrationally hot (its
lifetime against R5, R8 and R9 is about 1e-3 s against a radiative vibrational
relaxation of about 1 s, READ from the code's own site note), so that condition
does not hold here. The split is not neutral: 2p decays by Lyman-alpha and 2s
only by the two-photon continuum, so an all-2s assignment would route the same
energy out of the gas by a different channel and leave a larger n = 2
population. **Open.**

**Validity, written at the site.** Takagi's n = 2 product holds "only for the
low vibrational molecular ion and at low collision energies"; outside that
condition "the halves of product atoms are distributed to the highly excited
states". The code's H2+ is outside it. A vibrationally hot ion opens n >= 3,
which takes MORE of the enthalpy out of the kinetic channel, so the share
written is an UPPER BOUND on the heat and a lower bound on the excitation.
Guberman computed 3HeH+ in its ground vibrational state at electron energies
0.001 to 0.33 eV and calls the cross sections isotopomer sensitive.

### 1.2 R6: the H2 fragment is born hot

Same routine. The heat becomes `(q(R6) - E_int) + E_int * f_quench`, with
`E_int` the midpoint of the two levels the measured product distribution peaks
between, taken from the code's own H2 ladder: `e_vib_v5_eV` = 2.292710 and
`e_vib_v6_eV` = 2.666380 eV, so `E_int` = 2.479545 eV, 26.8 per cent of q(R6).
Both constants are added to `h2_vibrational_relaxation.f90` beside the v = 1 and
v = 2 term energies it already carried from the same table.

**Uncertainty, stated at the site because it is not small.** This is the
POSITION OF THE PEAK, not the mean: both papers describe the distribution as
broad, neither tabulates it, and Strasser et al. quote "an uncertainty of about
one level". One level near v = 5 is 0.37 to 0.40 eV, i.e. +/-15 per cent of the
share. Where `f_quench` is 1 the share is deposited in full either way and the
uncertainty does not reach the gas at all. **Open: the mean of the
distribution.**

### 1.3 The third body of R12 and R15

`src/modules/lower_atmosphere/mol_rates.f90`.

READ from Cohen & Westberg (1983), J. Phys. Chem. Ref. Data 12, 531, p. 559,
the recommended-rate sheets of this reaction, in the units the sheets print:

| collider | coefficient [cm^6 molecule^-2 s^-1] | range | stated uncertainty |
|---|---|---|---|
| H2 | 2.8e-31 T^-0.6 | 50-5000 K | +/-0.2 in log at 300 K, +/-0.4 at 5000 K |
| H | 8.8e-33, temperature independent | 50-5000 K | +/-0.5 in log throughout |
| Ar | 1.9e-30 T^-1 | 77-5000 K | +/-0.3 in log throughout |

New in `mol_rates`: `k3b_H_H_to_H2_atomic_H`, `k3b_H_H_to_H2_monatomic` and
`h2_association_collider_density(T, n_H2, n_HI, n_He)`, which returns the
H2-EQUIVALENT third-body density, i.e. the collider sum divided by k1(H2). Both
directions of the pair are multiplied by that one density, so the ratio of the
two is the equilibrium constant for every collider separately and the pair
stays an exact detailed balance whatever the mixture is. That is asserted at
four temperatures and three pure-collider mixtures in the test.

**Helium has no coefficient of its own anywhere in hand and is given argon's,
which bounds it from above.** The only published helium third-body calculation,
Paolini, Ohlinger & Forrey (2011), resolves the product states but tabulates no
rate: its total is shown in figures over 0-350 K only (READ, their Figs. 6 and
7), which neither reaches this layer's 200-1600 K nor gives a temperature
dependence. Argon is the nearest tabulated monatomic inert third body of the
same evaluation, and since it is heavier and more polarizable it forms the
longer-lived collision complex, so k1(He) <= k1(Ar). Against k1(H2) the argon
coefficient is 0.43 at 1000 K and 0.86 at 200 K, so the substitution removes the
factor ~2 by which applying k1(H2) to helium overstated the monatomic colliders.
**Open: a helium three-body coefficient with a temperature dependence.**

**R13 is NOT collider-resolved**: its source is a single coefficient and no
evaluation in hand splits it by third body, so it keeps the total density. Said
at the site.

**Every other heavy particle is left out** of the collider sum rather than
counted at the H2 efficiency: no evaluation in hand gives a coefficient for H+,
for the molecular ions or for a metal atom, and in the molecular layer they are
together below 1e-4 of the sum. Above the H2 front they are not, but there the
H2 both reactions act on has gone.

### 1.4 The radiative rate of the cascade

`src/modules/lower_atmosphere/h2_vibrational_relaxation.f90`.

`h2_vibrational_relaxation_init` reduces the H2 line list of
`molecular_infrared_data` (Roueff et al. 2019, 302 levels, 1833
electric-quadrupole and magnetic-dipole lines) to the LARGEST total spontaneous
decay rate over the ladder, grouping lines by the upper term energy each line
carries. **MEASURED: 5.5939e-06 s^-1**, against 8.3e-07 s^-1 for the v = 1 value
the fraction used before, a factor 6.56. The initializer is called from
`h2_thermochemistry_init` in `mol_rates`, the serial prologue every molecular
configuration passes through before the first parallel region opens; the module
variable falls back to Hollenbach & McKee's published v = 1 value if it has not
run, so a caller that skips the prologue gets an older published number and not
a meaningless one.

**The level ladder itself is not read**, only the line list: the one Boltzmann
sum over the ladder stays the single place that touches it
(`h2_partition_sum_uniqueness.py` enforces that, and the first implementation of
this reduction tripped it).

**Cross-checked against the complete quadrupole set.** The tree's line list
leaves fifteen near-dissociation levels with no downward transition, so its
maximum could in principle be the maximum of an incomplete network. It is not.
MEASURED on the complete set of Wolniewicz, Simbotin & Dalgarno (1998), which
connects ALL bound rovibrational levels (4661 lines, read from the copy the
Cloudy c25.00 distribution carries, the published paper being
`references/Wolniewicz_1998_ApJS_115_293.pdf`): the maximum is 5.5945e-06 s^-1,
agreeing to 1.1e-04 relative, and it sits at the SAME level, (v = 1, J = 28), a
high-J rotational level of low v and not a near-dissociation one. Against the
complete set only three levels of the code's ladder have no downward transition
at all, (0,0), (0,1) and (14,4). **The tree's own list is therefore adequate for
this number and no new table is imported.**

### 1.5 The collisional coefficients: three calculations in place of two fits

Same module. Each coefficient is the THERMAL v = 1 -> v' = 0 relaxation: for
each temperature of the source's own grid, the rates out of the initial
rotational levels of v = 1 summed over every final j' of v' = 0 and averaged
over the initial j with the Boltzmann weights `g_I (2j+1) exp(-E(1,j)/kT)` of
the code's own ladder. All three are tabulated in the module and interpolated in
`log k`, linear in T, with the end values held rather than extrapolated.

| collider | source | grid | validation |
|---|---|---|---|
| H | Lique (2015), MNRAS 453, 810, supplementary data `references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat`, 1431 transitions | 100-5000 K, step 100 K | MEASURED 1.781e-13 cm^3 s^-1 at 300 K against his Table 2's 1.8e-13 |
| He | Jozwiak, Thibault, Viel, Wcislo & Lique (2024), A&A 685, A113, VizieR tables `references/Jozwiak_2024J_A+A_685_A113/`, 26 ortho + 27 para levels, 520 + 539 converged transitions | 20-8000 K, 43 points | MEASURED 2.639e-15 at 808 K and 7.015e-15 at 1000 K, reproducing the inventory's independent reduction of the same files to 0.6 per cent |
| H2 | Le Bourlot, Pineau des Forets & Flower (1999), MNRAS 305, 802, `references/LeBourlot_1999_MNRAS_305_802.pdf`, data in the Cloudy c25.00 tables | 100-6000 K, 9 points | MEASURED 3.508e-15 at 1000 K and 2.946e-14 at 1500 K |

**Why the fits had to go, MEASURED.** The Hollenbach & McKee (1979) eq. (6.29)
atomic-hydrogen fit is 1.14 (200 K), 3.47 (300 K), 1.88 (808 K), 1.55 (1000 K)
and 1.03 (1600 K) times the calculation, crossing it near 1750 K. Their H2-H2
fit is 371 (300 K), 107 (808 K), 54 (1000 K) and 21 (1527 K) times it. **This
also settles a discrepancy the module recorded as unchecked**: Burton,
Hollenbach & Tielens (1990) Table 7 reprint the H2-H2 coefficient with 18100 in
the exponent where eq. (6.29) has 12000, and against the calculation that
reading is 6.35, 5.15, 3.37 and 2.26 times too large at the same four
temperatures. Burton et al. are the closer of the two and both are too large;
the question of which fit is right is retired, because neither is used.

**The ordering the three sources give, MEASURED at 1000 K:** atomic hydrogen
7.51e-12, helium 7.01e-15, H2 3.51e-15 cm^3 s^-1. Atomic hydrogen is three
orders above the other two, its exchange channel being what makes it efficient,
and H2 is no better a quencher of H2 than helium is. The fits had the H2
collider three orders too close to the atomic one.

**Helium is now IN the collider sum.** The module header used to say helium was
"deliberately absent ... that understates the de-excitation rate", an
open-ended caveat; it is now a term with its own published coefficient.
`h2_vibrational_heat_fraction` takes `n_He` as an OPTIONAL argument.

**What is NOT used, and why.** The Cloudy distribution also carries the ORNL
H2-H2 set of Wan et al. (2018): MEASURED, all 240 of its rows are purely
rotational, v = 0 to v = 0, so it contains no vibrational transition and cannot
give this coefficient. And `coll_rates_He_ORNL.dat` covers far more of the
ladder (v = 0 to 14, 137 v = 1 -> 0 rows, T = 1 to 1e4 K) than the Jozwiak
tables, but its stated source is "Lee et al. 2007, ApJ, in preparation", which
is unpublished; no number is taken from it. Both are recorded, neither is used.

### 1.6 The exact form the fraction approximates, recorded and not built

Written at `h2_vibrational_heat_fraction` on the Cloudy reading: the exact
statement of the same physics is the net collisional heat of the ladder,
`sum over level pairs and colliders of (n_u C_ul - n_l C_lu) n_M (E_u - E_l)`,
with the upward coefficients from the downward by detailed balance and the
populations from a statistical equilibrium of the whole X ladder. It carries no
n_cr, no single level standing in for the cascade and no fraction at all, and it
is exact in both limits. **It is not buildable here, for data and not for
effort**: the H2-H rate coefficients stop at v = 3 in everything that
distributes numbers, while the molecules this fraction is applied to are born at
v = 10 to 14 (R15) and v = 5 to 6 (R6). The distance from the distinction
mattering in this layer is `1 - f_quench` = 6.8e-07 at the certified base
(MEASURED, section 3).

A second record at the R15 site: depositing the internal share as heat wherever
the ladder is thermalized, with no prompt translational part, is the same
assumption a code that solves the ladder makes without a nascent distribution.
Cloudy spreads a newly formed molecule over the levels in LTE proportions, which
leaves the net collisional heat of the formation exactly zero and the whole bond
energy in the thermal pool, i.e. f_trans = 0 with f_quench = 1 by another name.

### 1.7 One comment outside the recipients: the Lyman-Werner fragment energy

`src/modules/lower_atmosphere/lyman_werner.f90` deposits a constant 0.4 eV of
kinetic energy per Solomon-process dissociation, on Black & Dalgarno's (1977)
"about 0.4 eV", with no statement of how good that is. The check had already
been made: weighting the level-resolved mean kinetic energies of Abgrall,
Roueff & Drira (2000) by the module's own rate weights gives 0.397 eV at 100 K,
0.406 at 1300 K and 0.429 at 2700 K (READ,
`docs/p39_lw_cross_section_sources.md` section 4.1), so the constant is right to
2 per cent over the layer. **The site now says so and points there. The constant
is unchanged and nothing was re-measured.**

### 1.8 The control key

`EXHALE_REACTION_HEAT_RECIPIENTS=0` restores items 1, 2 and 3 TOGETHER: R5 and
R16 deposit their full enthalpy and send nothing to n = 2, R6 deposits its full
enthalpy promptly, and R12 and R15 take the total heavy-particle density as
their third body. It is read once by
`reaction_heat_recipients_corrected()` in `molecular_reaction_heat`, and the two
consumers (the heating assembly and the H(n=2) update) both run serially. `0` is
the only value it reads; anything else, including an unset key, is the corrected
physics.

**It does NOT restore items 4 and 5**, the all-level radiative rate and the
three collider coefficients. Those are corrections of the same item with no key,
and the residual between the control binary and the new binary under the key is
exactly them: 6.1e-07 relative on `heat_mol_chem` at cell 1 (MEASURED,
section 3).

---

## 2. Tests

`src/tests/physics_probe/molecular_energy_recipients.f90`, new, registered in
`src/tests/physics_probe/run.sh`. 62 assertions in five groups:

- **A. the enthalpy is untouched**: q(R5), q(R16) and q(R6) against the values
  the network audit read from the same table; `E(n=2)` is the one table entry
  and is 0.75 IP(H); heat + radiated-or-excited = enthalpy for each of the three
  channels; the kinetic shares 0.749 and 1.554 eV and the R6 internal share
  2.4795 eV.
- **B. R12 and R15 are an exact detailed-balance pair collider by collider**: at
  300, 800, 1500 and 3000 K, for a gas of pure H2, pure H and pure He, the ratio
  of the two directions is `keq_H_H_to_H2` to 1e-13 relative, and each
  collider's rate is its own published coefficient.
- **C. the limit row**: with n(H) = n(He) = 0 and n(H2) the whole density the
  H2-equivalent third body is the total density exactly, so the correction
  replaces a MIXTURE and not a normalization.
- **D. the radiative and collisional data**: the all-level maximum against the
  inventory's 5.594e-06 and against the v = 1 total; the three collider tables
  against the values their own papers publish; the detailed balance of the
  helium table at four temperatures, the excitation and de-excitation rows being
  READ literals in the test and the energy gap the module's own `e_vib_v1_eV`.
- **E. the ledger deposits what the recipients leave it**: `heating_of_composition`'s
  own molecular assembly run on cells that carry ONE channel at a time, so the
  whole output is that channel's: R5 alone, R16 alone, R6 and R7 alone at zero
  collider density (where the internal share is withheld in full), and R15 alone
  in pure atomic hydrogen (where the rate must be `k1(H) n(H)^3`); the H(n=2)
  source is the same two rates with the same densities; and the ledger does not
  depend on `use_excited_H`.

**RED, MEASURED.** Against the entry text -- a scratch copy of the tree with the
v = 1 radiative rate, the two Hollenbach & McKee fits and no helium collider
restored, run with `EXHALE_REACTION_HEAT_RECIPIENTS=0` -- **17 assertions FAIL
and 45 pass**: the two all-level-A rows, the four table rows and the helium
share, the two efficiency-ordering rows, the four ledger rows and the three
n = 2 source rows.

**GREEN, MEASURED.** On the delivered code **all 62 pass**, and the whole
`physics_probe` suite reports PASSED over its 1577 assertions.

**Every other suite whose drivers link the changed modules**, MEASURED with the
private binary, single-threaded:

| suite | result |
|---|---|
| `physics_probe` | PASSED |
| `ionization_imposed_fractions` | PASSED, 31 assertions |
| `coupled_source_step` | PASSED, 30 assertions |
| `spectrum_type` | PASSED, 170 assertions |
| `grid_and_gates` | 220 pass, 1 FAIL: `outer_iteration_ending_is_the_stagnation_one` |
| `a2_m1` | no `run.sh`; it is a driver plus a Python cross-check, not a runnable suite |

**The `grid_and_gates` failure is not this item's.** MEASURED: the CONTROL
binary, built from the tree before any edit of this item, fails the same row with
the same measured value (`pass_budget` against `no_progress`). Reported, not
fixed.

---

## 3. Impact

### 3.1 The certified molecular base state, held as it stands

`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13` reloaded on a scratch
copy with `Load IC? True` and `Restart intent: stationary evaluate`, so the
state is measured and not moved; single-threaded. Control = a private build of
the tree as it stood before this item.

`heat_mol_chem` (column 19 of `output/Heating_breakdown.txt`), MEASURED
[erg cm^-3 s^-1]:

| | cell 1 (808.3 K) | cell 25 (1527.3 K) |
|---|---|---|
| control | 1.11205700e-06 | 6.75251873e-08 |
| corrected | 8.71940708e-07 | 6.12082692e-08 |
| corrected, `EXHALE_REACTION_HEAT_RECIPIENTS=0` | 1.11205632e-06 | 6.75251702e-08 |
| **corrected against control** | **-21.6 %** | **-9.4 %** |
| control key against control | -6.1e-07 relative | -2.5e-07 relative |

The total heating of the same cells falls from 1.127911e-06 to 8.877949e-07
(-21.3 per cent) and from 7.685033e-08 to 7.053341e-08 (-8.2 per cent).

**The control key closes the difference to 6e-07, not to zero, and that residual
is items 4 and 5**, which the key does not cover: the fraction moved from
1 - 5.4e-08 to 1 - 6.8e-07.

### 3.2 Where the movement comes from: the R15 rate

MEASURED at cell 1 from the state's own composition (n(HI) = 1.8680e12,
n(HeI) = 6.1716e12, n(H2) = 5.1474e11 cm^-3, T = 808.3 K):

| | cell 1 | cell 25 |
|---|---|---|
| total heavy-particle density [cm^-3] | 8.5543e+12 | 3.5437e+12 |
| H2-equivalent collider density [cm^-3] | 6.6521e+12 | 3.1171e+12 |
| ratio | 0.7776 | 0.8796 |
| k15 two-body [cm^3 s^-1], old | 4.3132e-20 | 1.2197e-20 |
| k15 two-body [cm^3 s^-1], new | 3.3540e-20 | 1.0729e-20 |
| R15 volumetric rate [cm^-3 s^-1], old | 1.5050e+05 | 7.3188e+03 |
| R15 volumetric rate [cm^-3 s^-1], new | 1.1704e+05 | 6.4378e+03 |
| **change** | **-22.2 %** | **-12.0 %** |

R15 carries 95 per cent of `heat_mol_chem` at cell 1 (READ, predecessor memo),
so 0.95 x 22.2 = 21.1 per cent accounts for essentially the whole 21.6 per cent
fall. The recipient corrections contribute the rest.

### 3.3 The collider sum and the thermalized fraction

MEASURED at the same two cells with the three new coefficients:

| | cell 1 | cell 25 |
|---|---|---|
| share of the collider sum, atomic H | 99.796 % | 99.036 % |
| share, He | 0.198 % | 0.918 % |
| share, H2 | 0.0059 % | 0.0458 % |
| `1 - f_quench`, corrected | 6.79e-07 | 3.79e-07 |
| `1 - f_quench`, entry text (two fits, no He, v = 1 A) | 5.37e-08 | 5.23e-08 |

So the fraction is still 1 to a part in 1e6 and no deposit in this layer is
held back; what the correction buys is that the statement is now true of the
whole cascade and of a complete collider sum. **At 1527 K against an H2 collider
ALONE at n(H2) = 4e9 cm^-3 the fraction would be 0.958**, and the retired fit
would have claimed 0.999 there: that is the configuration in which item 5
matters.

### 3.4 The H(n=2) source

**`use_excited_H` is OFF in this configuration** (the code arms it only when
both a stellar effective temperature and a stellar radius are stated, and this
case loads a numerical SED instead), which the run confirms: `heat_Hpe` and
`heat_Hdx` are identically zero in `Heating_breakdown.txt`. So this state
exercises the "when it is off" branch: the 10.199 eV leaves the ledger as
radiation and no trapping is modeled, the way R23's 153 nm photon does.

Had it been on, MEASURED from the same state: the source would be 1.10e-02 cm^-3
s^-1 at cell 1 and 1.45e-01 at cell 25, carrying 1.80e-13 and 2.37e-12
erg cm^-3 s^-1, i.e. 2.1e-07 and 3.9e-05 of `heat_mol_chem` there. The two
carriers are traces in this layer (n(H2+) = 7.9e-02 and n(HeH+) = 1.6e-05 cm^-3
at cell 1), which is why the recipient correction for R5 and R16 does not move
this base.

### 3.5 The bounded re-solve

Same scratch copy, `Restart intent: stationary`, `EXHALE_OUTER_PASSES=3`,
single-threaded, both binaries.

**The control binary accepts the loaded state at once.** MEASURED: outer pass 1
returns `hydro info = 0` with mass 7.71e-09, momentum 2.21e-12 and the worst
gated species row 5.64e-07 of its 1.0e-05; pass 2 reports
`ACCEPTED -- every active equation of this state is within its own tolerance`
and the flux gate accepts 3.1661e-11. The state is its own root for that binary,
which is what "certified" means and why it is the archive's.

**The corrected binary does not reach a root in the bounded budget, and that is
the finding.** MEASURED: outer pass 1 returns `hydro info = 2` (not converged)
with mass 2.48e-06, momentum 1.76e-10 and **energy 4.27e-01**, in 2178 s; outer
pass 2 the same, mass 1.70e-01, momentum 2.73e-06, energy 4.27e-01, in 2337 s,
with the worst gated species row 2.38e-09 (the elemental He/H partition, well
inside its tolerance). Outer pass 3 was still running at JFNK iteration 48, with
`||R||` oscillating between 1.9e-01 and 8.6e-01 and the worst row the momentum
of cell 1, when the two-hour wall-clock cap of the run stopped it. **No output
state was written, so no base temperature, x2, cell-1 mass row or Mdot can be
quoted from it, and none is.**

**What the two tell us together.** The energy row, and only the energy row, is
what refuses: its residual sits at 4.27e-01 through both completed passes while
mass, momentum and every species row are at or near their tolerances. That is
the direct signature of the correction: `heat_mol_chem` fell 21.6 per cent at
the base, so the energy balance that the archived state satisfied is no longer
satisfied there and the layer has to find a new temperature. **The certified
state is NOT a solution of the corrected equations, and re-solving it is a
campaign and not a bounded pass.** That is the item's largest consequence for
the plan and it is handed over rather than attempted here.

### 3.6 The regression cases

**Matched short runs, not the cases' own step counts.** The three molecular
cases of the matrix are 12000-step relaxation snapshots and, on the machine as
it was loaded, one of them takes hours per binary; six of them, twice over, is
not a measurement this item can complete. What is reported instead is a MATCHED
pair: each case run from its own `output/*_IC.txt` with `EXHALE_MAXSTEPS=500`
under both binaries, single-threaded, on scratch copies, so the two states
compared have taken the same 500 steps from the same initial condition. It is a
lower bound on the movement the full case will show, not the movement itself.

MEASURED, maximum relative change over all 504 rows, control against corrected:

| case | rho | v | p | T | heat | cool | x2 | worst species above 1 cm^-3 |
|---|---|---|---|---|---|---|---|---|
| `mol_base_handoff` | 3.1e-04 | 2.2e-03 | 3.0e-03 | 2.6e-03 | 1.4e-02 | 2.8e-03 | 7.8e-03 | H2 7.8e-03 |
| `mol_sec_ion` | 2.9e-04 | 5.9e-03 | 3.0e-03 | 2.5e-03 | 1.1e-01 | 2.7e-03 | 7.7e-03 | H2 7.7e-03 |
| `mol_ir_bands` | 3.1e-04 | 2.2e-02 | 3.0e-03 | 2.6e-03 | 1.2e-02 | 2.8e-03 | 7.7e-03 | H2 7.8e-03 |

The largest movements sit where the layer is, r = 1.07 to 1.60; the temperature
moves by at most 2.6e-03 and the H2 content by 7.8e-03 after 500 steps. **All
three goldens will move and none of them by much at this depth of relaxation;
the full 12000-step number is not measured here and the advisor's refresh should
take it from the matrix run.** No golden was read or written by this item.

**The atomic case is byte-identical.** `wasp_he23off` run from its own IC with
`EXHALE_MAXSTEPS=300` under both binaries: `Hydro_ioniz.txt`,
`Ion_species.txt`, `Heating_breakdown.txt`, `Cooling_breakdown.txt`,
`Excited_H.txt`, `IC_dump.txt`, `OI_levels.txt`, `Ion_species_adv.txt` and
`OI_levels_adv.txt` compare byte for byte. `Hydro_ioniz_adv.txt` differs in its
HEADER ONLY, and by another item's change: the control build predates item L9's
replacement of the `# coupling:` line by `# derived_from:`. Its NUMBERS are
identical to the bit, maximum relative difference exactly 0.0 over all
504 x 10 entries (MEASURED). **So nothing this item changed executes in an
atomic configuration**, which is what the gates on `thereis_mol` assert.

**One consequence of that header difference, reported because it bears on every
number above**: the control binary is a build of the tree as it stood when this
item opened, and other items landed in it afterwards. The 21.6 per cent fall in
`heat_mol_chem` is nevertheless attributable to this item alone, because the
SAME new binary run with `EXHALE_REACTION_HEAT_RECIPIENTS=0` reproduces the
control binary's value to 6.1e-07 (section 3.1): had another item moved the
molecular heating numerically between the two builds, that agreement could not
hold.

---

## 4. What remains open

1. **The mean of the R6 product distribution.** The peak is settled in both
   primary papers; the mean is not, and the spread is about one vibrational
   level, +/-15 per cent of the internal share. It reaches the gas only where
   `f_quench` < 1.
2. **The l-partition of the n = 2 source.** Statistical, 1 : 3, because the two
   general sources resolve only n. It decides how much of the 10.199 eV leaves
   as Lyman-alpha and how much as the two-photon continuum.
3. **A helium three-body coefficient with a temperature dependence.** Paolini et
   al. (2011) is the only helium calculation and tabulates nothing; argon stands
   in, as an upper bound.
4. **The level-resolved cascade**, and with it the exact net-collisional-heat
   form recorded at the fraction. Blocked on the same gap the inventory
   identified: nothing that distributes numbers reaches past v = 4 for H2-H,
   and the nascent molecule is born at v = 10 to 14. The rate tables of Bossion
   et al. (2018) remain item 1 of the request list.
5. **`gamma_10_H2` rests on a data file whose paper is now in hand but whose
   numbers have not been checked against it line by line.** The published PDF
   arrived during this item; the tabulated reduction should be verified against
   the paper's own tables.
6. **Three call sites do not pass `n(He)` to the heat fraction.** The two
   Lyman-Werner and photoelectric sites in `util_ion_eq.f90` and the
   photoelectron partition in `ionization_equilibrium.f90` are handed no helium
   density and so drop a term worth 0.2 per cent of the collider sum at this
   base. Passing it is a one-line change in files this item does not own.
7. **The certified base state has to be re-solved, and it is a campaign.**
   Section 3.5: the corrected equations leave the archived state with an energy
   residual of 4.27e-01 through two full outer passes while every other row is
   at its tolerance, and three bounded passes did not reach a root within two
   hours. `molecular_scalar_gj1132_kzz1e9/HeH2.13` and every model built on it
   are states of the old ledger.
8. **The chemistry solver still hoists the TOTAL-density third body.**
   `set_mol_coeffs` in `System_HeH_mol.f90` computes k15 and k13 once per cell
   from `(T, n_tot)` because they were invariant across the cell's Newton solve;
   a collider-resolved k15 is a function of the unknowns, so carrying the
   correction into the composition solver means moving that coefficient inside
   the residual and the Jacobian. **The correction therefore reaches the ENERGY
   ledger and not yet the CHEMISTRY**, and the two now use different R15 rates in
   the same run. This is the largest thing this item leaves undone and it is
   outside its file scope. The same holds for the R15/R12 diagnostic in
   `write_output.f90`.

---

## 5. Files changed

| file | what |
|---|---|
| `src/modules/lower_atmosphere/molecular_reaction_heat.f90` | the R5/R16 and R6 recipients, the collider-resolved third body of the R12/R15 pair in the ledger, the control key, the n = 2 source and excitation-energy accessors |
| `src/modules/lower_atmosphere/mol_rates.f90` | the two further Cohen & Westberg coefficients and the H2-equivalent collider density; the vibrational-relaxation initializer called from the serial prologue |
| `src/modules/lower_atmosphere/h2_vibrational_relaxation.f90` | the all-level radiative rate reduced from the line list; the three tabulated collider coefficients; the optional helium term; the v = 5 and v = 6 term energies; the recorded target form |
| `src/modules/radiation/excited_hydrogen.f90` | the chemical n = 2 source in the rate matrix, the populations and the level residual |
| `src/modules/lower_atmosphere/lyman_werner.f90` | comment only: the Abgrall et al. (2000) range behind the adopted 0.4 eV |
| `src/tests/physics_probe/molecular_energy_recipients.f90` | new driver, 62 assertions |
| `src/tests/physics_probe/run.sh` | the new driver registered |
| `docs/input_schema.md` | the control key and what it does and does not restore |
| `docs/PLAN_20260913_lhs_stationary.md` | row L7g |
| `docs/Update_EXHALE_stage2.md` | one bullet in section 11 |

---

## Follow-up: one third-body rate, the helium collider, the Le Bourlot check

Item L7h, 2026-09-17. Three of the items section 4 left open are closed here:
number 8 (the chemistry still hoisted the total-density third body), number 6
(three call sites did not pass the helium collider) and number 5 (the Le
Bourlot reduction had not been checked against the published paper). Every
number below is MEASURED here or READ from the source it names.

### 1. One R15 rate for the composition and for the energy

`set_mol_coeffs` in `src/modules/nonlinear_system_solver/System_HeH_mol.f90`
hoisted `mk15 = rk_R15_3body_H2(T, n_tot)` once per cell while the energy
ledger of section 1.3 formed the collider sum, so the H2 abundance was set by
one rate and its heat charged at another. The collider sum is now formed
inside the rows, at the composition the row is evaluated at, by the same
routine the ledger calls (`mol_rates::h2_association_collider_density`), and
`mk15` carries the M = H2 coefficient k1(H2) alone. Both directions of the
pair take that one density, so R12 and R15 remain an exact detailed balance
collider by collider. The helium collider of the row is the ground singlet
`n_heiSI`, which is the density the ledger passes. `EXHALE_REACTION_HEAT_RECIPIENTS=0`
now restores the total density in the chemistry as well as in the ledger, read
from the one accessor that owns the key.

The turnover scale of the H2 row takes the collider sum evaluated with each
collider at the whole of its element, n(H2) = n_H/2, n(H) = n_H and n(He) =
n_He, which no mixture reaches at once and which is therefore still an upper
bound on the row.

The `write_output.f90` base H2 loss budget, which reported R12 and R15 at the
total density, forms the same collider sum and is gated on the same key.

**What it changes, MEASURED on the certified base state.** The state was
reloaded and HELD (`Restart intent: stationary evaluate`), single-threaded,
against a control binary built from the tree as this item found it. Because
the state is held, every written composition is bitwise identical between the
two binaries: x(He II) and n(H2) at cells 1, 25 and 100 and the H2 column,
2.428387e+30 cm^-2, do not move, and neither do T, rho, `heat_mol_chem`
(8.71940708e-07 and 6.12082692e-08 erg cm^-3 s^-1 at cells 1 and 25) or the
total heating. That is the ledger's number and this item did not touch it. The
same binary with `EXHALE_REACTION_HEAT_RECIPIENTS=0` returns 1.11205632e-06
and 6.75251702e-08, the control key's +27.5 and +10.3 per cent, which is
section 3.1 read the other way round.

**What does move is the H2 carrier balance of that state, and it moves by six
orders.** MEASURED from the certification inventory of the same two runs:

| | control | with the collider-resolved chemistry |
|---|---|---|
| carrier balance H2, gated row measure | 5.645e-07 at cell 244 | **7.817e-01 at cell 243** |
| its tolerance | 1.0e-05 | 1.0e-05 |
| verdict | within | ABOVE |
| refusing entries of the inventory | 2 | 3 |

The reason is that the collider sum and the total density part company ABOVE
the H2 front, in the opposite direction to the base. MEASURED from the held
state's own composition:

| cell | r [R_p] | T [K] | n_tot [cm^-3] | collider sum [cm^-3] | ratio | R15 rate, total density | R15 rate, collider sum | change |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.0002 | 808.3 | 8.5543e+12 | 6.6521e+12 | 0.778 | 1.5050e+05 | 1.1703e+05 | -22.2 % |
| 25 | 1.0048 | 1527.3 | 3.5437e+12 | 3.1171e+12 | 0.880 | 7.3188e+03 | 6.4378e+03 | -12.0 % |
| 243 | 1.3165 | 5635.6 | 1.4431e+09 | 3.5890e+09 | 2.487 | 8.4488e-07 | 2.1012e-06 | +148.7 % |
| 300 | 1.8571 | 5046.7 | 1.5249e+08 | 4.2712e+08 | 2.801 | 1.5860e-09 | 4.4424e-09 | +180.1 % |

(volumetric rates in cm^-3 s^-1). In the wind the third body is atomic
hydrogen, whose coefficient is temperature independent while k1(H2) falls as
T^-0.6, so the collider sum stands 2.5 to 2.8 times the total density at 5000
to 5600 K where at the base it stood 0.78 of it. The two columns are what the
run carried before this item: the left one in the H2 balance row, the right
one in the heat charged for the same reaction. They are now one number, and
the test asserts it to 1e-14 relative.

**Consequence for the archived states, reported and not acted on.** The H2
carrier row of `molecular_scalar_gj1132_kzz1e9/HeH2.13` is no longer within
its tolerance above the front. That state was already not a solution of the
corrected energy equation (section 3.5); it is now not a solution of the
corrected H2 balance either, and the re-solve that section 3.5 handed over has
one more row to satisfy.

**The regression cases**, run from their own initial condition with
`EXHALE_MAXSTEPS=500` under both binaries, single-threaded, on scratch copies.
MEASURED, maximum relative change over all 504 rows:

| case | rho | v | p | T | heat | cool | n(HeII) | n(H2) | worst species above 1 cm^-3 |
|---|---|---|---|---|---|---|---|---|---|
| `mol_base_handoff` | 2.5e-02 | 3.4e+00 | 3.8e-02 | 4.3e-02 | 3.4e-02 | 7.2e-01 | 8.5e-01 | 6.5e+00 | H3+ 9.5e+00 |
| `mol_sec_ion` | 2.4e-02 | 9.4e+00 | 3.0e-02 | 2.6e-02 | 7.7e-02 | 4.0e-01 | 8.3e-01 | 5.8e+00 | H3+ 1.3e+01 |

**These are not small and they sit in one place, the H2 front.** Every maximum
above falls between r = 1.09 and 1.19, and the temperature's median relative
movement over the whole grid is 3.9e-04 and 5.2e-04 with a 95th percentile of
2.2e-02. At r = 1.1097 in `mol_base_handoff`, MEASURED, n(H2) goes from
2.1523e+09 to 1.6226e+10 and n(H3+) from 0.912 to 9.61 cm^-3 while T rises from
838.4 to 858.9 K and rho moves 0.5 per cent; at r = 1.1861 the two states are
within 20 per cent of each other again. The front moves OUTWARD, which is the
direction the rate table above requires: where the third body is atomic
hydrogen at 900 to 5600 K the collider sum stands 2.5 to 2.8 times the total
density, so the association that makes the H2 runs that much faster than the
chemistry was running it. `log10 Mdot` is 7.93 for both binaries in both cases.
All the molecular goldens will move, and by much more than the 500-step
movement the predecessor item measured for the ledger alone. No golden was read
or written by this item.

### 2. The helium collider at the three call sites

`h2_vibrational_heat_fraction` takes the neutral helium density as an optional
third collider (section 1.5) and three callers did not pass it. They do now,
each with the GROUND-SINGLET neutral, which is what the quantum calculation
behind the coefficient is for and what the chemistry and the ledger already
used; a helium ion interacts with the molecule through a different potential
and carries no coefficient here. The sites are the Lyman-Werner fluorescence
channel of `heating_of_composition` and the same channel in
`fuv_band_absorption_ledger` (`src/modules/radiation/util_ion_eq.f90`) and the
photoelectron partition's `f_vib_quench` (`ionization_equilibrium.f90`). The
band ledger had no helium density among its arguments and is given one;
`write_output.f90`, its only caller, passes the same array it passes for every
other density. The ground singlet is now formed once in
`heating_of_composition` and read by both of its molecular channels.

At this base the helium term is 0.20 per cent of the collider sum (READ,
section 3.3), so the fraction does not move here; what the change buys is that
one collider sum decides the thermalized share everywhere it is asked for.

### 3. The Le Bourlot check, and the correction it required

`references/LeBourlot_1999_MNRAS_305_802.pdf`, read with `pdftotext -layout`
and, for Table 1, from the page image.

**The paper publishes no rate coefficient.** Its section 2.2 gives a
three-parameter fit form and states that the numbers are distributed through
the CCP7 library; the only numbers in the paper are Table 1, the critical
densities n_cr = A/q of eleven transitions for each of the three perturbers at
500, 1000 and 2000 K. The comparison therefore goes through that table.

**The file and the paper disagree about how many H2-H2 sets there are, and the
code was using the wrong one.** Cloudy carries two files named for an ortho-H2
and for a para-H2 perturber and fills two separate collider slots from them
(`parse_atom_h2.cpp`), and the module combined them at the equilibrium ortho
fraction. The paper says there is one set. Published text, section 2.1: Flower
& Roueff (1998a) computed the excitation of ortho- and para-H2 "by para-H2 in
its rotational ground state (J = 0)"; test calculations with ground-state
ortho-H2 as the perturber "showed that the rate of v = 1 -> 0 vibrational
relaxation was insensitive to the rotational state of the perturber"; and "we
have, therefore, adopted the rate coefficients calculated by Flower & Roueff
(1998a) for the excitation of ortho- and para-H2 by para-H2 (J = 0), and
applied them also to the case of excitation by ortho-H2". Three further
readings agree that the para file is that set and the other is not:

- the paper states "A total of 51 rovibrational energy levels ... (v = 0,
  J <= 16; v = 1, J <= 13; v = 2, J <= 10; v = 3, J <= 8)". MEASURED, the para
  file holds exactly that, 51 levels to (3,8) in 627 rows; the file named for
  the ortho perturber holds 39 levels reaching only (1,13) in 361 rows, and its
  transition list is a strict subset of the other's.
- MEASURED, the two are not the same data: on the 361 rows they share the
  ortho file stands a median 1.74 above the para file at 1000 K (1.15 to 17.1
  on the purely rotational rows, 0.75 to 197 on the vibrational ones).
- **Table 1, MEASURED against both.** For 1-0 S(1), i.e. (v = 1, j = 3) ->
  (v = 0, j = 1), with A = 3.470e-07 s^-1 READ from the tree's own line list
  (`molecular_infrared_data`, the 2.1219 micron line from the 6951.3 K level):

  | perturber and file | 500 K | 1000 K | 2000 K |
  |---|---|---|---|
  | published Table 1, H2 | 3.7e+11 | 1.8e+10 | 4.2e+08 |
  | para file | 2.33e+11 | 1.65e+10 | 5.80e+08 |
  | ratio, published / para file | 1.59 | 1.09 | 0.72 |
  | ortho file | 1.13e+11 | 7.5e+09 | 2.94e+08 |
  | ratio, published / ortho file | 3.3 | 2.4 | 1.43 |
  | published Table 1, He | 4.2e+10 | 1.5e+09 | 6.8e+07 |
  | helium file of the same paper | 5.59e+10 | 1.56e+09 | 5.50e+07 |
  | ratio, published / helium file | 0.75 | 0.96 | 1.24 |

  (critical densities in cm^-3.) The helium row is what establishes that the
  procedure and not the data is under test: the same arithmetic on the same
  paper's helium file scatters about one within 25 per cent. Against that, the
  para file scatters about one and the ortho file is low by a factor 1.4 to
  3.3 throughout, i.e. its rates are systematically too fast.

**The correction.** `gamma_10_H2` is now the reduction of the para-H2
perturber file alone, for an H2 collider of either form. The reduction itself
is unchanged: thermal v = 1 -> v' = 0, summed over final j' and Boltzmann
averaged over the initial j of v = 1 with the weights g_I (2j+1) exp(-E(1,j)/kT)
of the code's own ladder. Reproducing the delivered table from the two files
and the ladder gives the OLD mixture to 2e-04 relative at every one of the
nine temperatures, which is what says the new column is the same reduction of
one file and not a different procedure:

| T [K] | 100 | 300 | 500 | 1000 | 1500 | 2000 | 3000 | 4500 | 6000 |
|---|---|---|---|---|---|---|---|---|---|
| before, ortho/para mixture | 6.7957e-18 | 2.1946e-17 | 1.1361e-16 | 3.5080e-15 | 2.9462e-14 | 1.1758e-13 | 6.0673e-13 | 2.0890e-12 | 4.0506e-12 |
| after, para file alone | 5.9364e-18 | 1.6258e-17 | 7.9918e-17 | 2.6310e-15 | 2.3015e-14 | 9.5076e-14 | 5.1397e-13 | 1.8230e-12 | 3.5736e-12 |
| ratio | 0.874 | 0.741 | 0.703 | 0.750 | 0.781 | 0.809 | 0.847 | 0.873 | 0.882 |

(cm^3 s^-1.) The coefficient falls by 12 to 30 per cent. What section 1.5 and
the module header state about it moves with it, MEASURED: the Hollenbach &
McKee (1979) eq. (6.29) fit is now 500, 147, 72 and 27 times the calculation at
300, 808, 1000 and 1527 K (was 371, 107, 54 and 21) and the Burton et al.
(1990) reading 8.57, 7.04, 4.50 and 2.89 times it (was 6.35, 5.15, 3.37 and
2.26); the three colliders at 1000 K are 7.51e-12, 7.01e-15 and 2.63e-15
cm^3 s^-1 for H, He and H2; the H2-only critical density is 2.1e+09 cm^-3 at
1000 K, and at 1527 K against an H2 collider alone at n(H2) = 4e9 cm^-3 the
thermalized fraction is 0.947, not the 0.958 section 3.3 quoted.

**Where it does and does not matter.** Not at this base: H2 is 0.006 per cent
of the collider sum at cell 1 and 0.046 per cent at cell 25 (READ, section
3.3), so 1 - f_quench moves by parts in 1e4 of itself and no deposit is held
back. It matters in a cold H2 base with no atomic hydrogen, which is the
configuration section 3.3 named.

### 4. Tests

`src/tests/physics_probe/molecular_third_body_in_the_chemistry.f90`, new,
registered in `src/tests/physics_probe/run.sh`. Nine assertions on a synthetic
cell that is the certified base in order of magnitude: the row's R15 rate
against the ledger's expression to 1e-14; that the total-density coefficient
is NOT that number, so the first assertion has content; that R12 and R15 carry
the same third body, their channel ratio reducing to n(H2)/(K_eq n(H)^2) with
no collider density left in it; the three derivative entries the numerical
Jacobian builds, d(R15)/d n(H2), d(R15)/d n(H) and d(R12)/d n(H), by central
differences of the production rows against the analytic derivative of the
collider sum, with a row asserting that the collider term is 20 per cent of
the second of them and not a rounding of it; and that the H2 row's turnover
bound rises with the cell's helium.

**RED, MEASURED** against an entry text (the delivered closure with
`System_HeH_mol.f90` and the `gamma_10_H2` table restored to what this item
found): 5 of the 9 fail, `one_k15_for_the_row_and_the_ledger` at 4.3132e-20
against 3.3541e-20, `collider_sum_is_not_the_total_density` at exactly 1,
`dR15_dnH2_is_the_H2_collider` and `dR12_dnHI_is_the_atomic_third_body` at
exactly zero, and `dR15_dnHI_has_the_collider_term` at 1.6114e-07 against
1.5602e-07. Two rows of `molecular_energy_recipients.f90` fail with them,
`lebourlot_H2_thermal_v1_v0_at_1000K` at 3.508e-15 against 2.631e-15 and the
1500 K one at 2.946e-14 against 2.302e-14. The four that pass on the entry
text do so for the reason each states: the detailed-balance row because the
frozen coefficient was one-sided in neither direction, the turnover row
because the total density carries helium too, and the two bound rows because
they are properties of the coefficients.

**GREEN, MEASURED**: all 9 pass and the whole `physics_probe` suite reports
PASSED over 1587 assertions.

### 5. What is still open

The three items closed here leave section 4's list at items 1, 2, 3, 4 and 7,
and item 7 grows: the archived base state now fails the H2 carrier row as well
as the energy row, and the two are the same re-solve.
