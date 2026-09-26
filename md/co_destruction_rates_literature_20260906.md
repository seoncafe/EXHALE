# CO destruction and formation rates in the ceiling's regime: what the
# published literature gives, and what it settles

**The ceiling this document is named for no longer exists.** It was deleted on
2026-09-06 (item CEILING-DEL), after the literature below was used to write the
destruction model of B1 T5.1 into the CO carrier row (items B3b-CO, B3b-CO2).
Every phrase below that speaks of "the ceiling's regime", of where "the ceiling
fires" or of "ceiling-active runs" is describing the tree as it stood on
2026-09-06, and B1 decision 7's exclusion was lifted by construction when the
observable it named was removed. MEASURED at the deletion: with the rates of
this document in the row, the ceiling never fired at all on the
`oxygen_chemistry` case, so every output column was byte-identical across its
removal. The literature, the rate values and every timescale measurement below
stand as made.

**Item COlit of PLAN_20260906_rev2. Read-only on the code.** This document
supplies the rate literature that B1 section 5.2 lists as the blocked part of
alternative A (the one-sided instantaneous CO destruction model), and then
uses it to answer the question alternative A cannot be written without: does
the timescale ordering `tau_dest << tau_res << tau_form` hold anywhere the CO
ceiling is active.

Every number below is labeled READ (transcribed from a source file or a
published paper) or MEASURED (evaluated here, by the script named at the
point of use, from data that is itself READ). No number is quoted from a
preprint. Sources that could not be obtained from the publisher are listed in
section 8 as requests, and nothing in this document depends on a number from
one of them.

---

## 1. Verdict, first

> **Superseded in part, and by later evidence in this same document.** The
> verdict below was reached with the neutral channels only, because no source
> for the ion channels was available when it was written. UMIST RATE22
> (section 12) supplies `He+ + CO -> O + C+ + He` as a measured, accuracy-A,
> temperature-independent Langevin rate whose products EXHALE already carries,
> and with it the timescale ordering `tau_dest << tau_res` **does** hold above
> the helium ionization front, by one to four decades (MEASURED, section 12.5).
> **The corrected verdict is section 12.6**: the one-sided destruction model of
> B1 T5.1 is constructible above that front and not below it, so the domain
> argument T5.2 returns a domain with a lower boundary rather than an empty
> one. What survives unchanged is the judgement on the present ceiling: its
> switching surface is the CO equilibrium constant, the destruction domain's
> boundary is `n(He+)`, and these are not the same surface, so B1 decision 7
> (exclusion) stays correct until T5.1 is implemented. Read section 1 as the
> record of what the thermal channels alone support.


**The one-sided destruction model of B1 T5.1 is not constructible in the
regime the ceiling acts in.** The requirement `tau_dest << tau_res` fails
there by two to ten orders of magnitude, and it fails for a reason that is
thermochemical rather than one of missing data: carbon monoxide is bound by
11.16 eV (MEASURED from the code's own Shomate table, section 3), so every
neutral channel that breaks it has a barrier of at least 6.7 eV and is
therefore slow at 3000 K by an amount that no choice of rate compilation can
change. The fast destroyers of CO (photodissociation in the 912 to 1118 A
predissociating bands, charge transfer with He+ and H+, electron impact) all
need partners that are absent at the base and only appear far above the
temperature where the ceiling first fires.

So B1's decision 7 (adopt alternative B, exclusion) is settled on physics and
not on data availability. What the literature does change is the reading of
the ceiling itself: section 7 shows the ceiling is a defensible kinetic
statement only above roughly 5000 K, is wrong by four to ten orders of
magnitude in the 2000 to 3500 K band, and that the band where it is wrong is
exactly the band where it fires first.

**Photodissociation, added in section 10 from the published Visser, van
Dishoeck and Black (2009) model, does not overturn this.** It is by far the
fastest CO destruction channel in the regime, twelve or more decades faster
than the barrier channels at 1456 K, and it is the channel any real CO model
would be built around. But with its published unshielded rate and its
published shielding function evaluated on the columns of an actual EXHALE run
(MEASURED, section 10.7: `log N(CO) = 18.96`, `log N(H2) = 21.80`, giving
`Theta ~ 2.6e-5` and `tau_dest ~ 2.8e8 s`), the ordering `tau_dest << tau_res`
still fails at the base by two to six decades. In the configuration the
molecular regression gate actually runs (`mol_base_handoff`) it fails by twelve,
because that configuration models no radiation longward of 911.8 A and so has
no photodissociating flux in the CO bands at all. And the standard model
cannot simply be evaluated at the temperatures of interest: Visser et al. stop
at `T_ex(CO) = 512 K` and state that no data exist for dissociating
transitions out of the `v'' = 1` level that is populated above about 500 K
(section 10.4).

**No compilation covers the band where the ceiling fires.** Section 11.4
collects the stated validity of all six sources consulted: KIDA 10 to 300 K,
UMIST RATE12 per entry with interstellar provenance, Venot 2012 and 2020
300 to 2500 K and 0.01 to 100 bar, VULCAN 500 to 2500 K, NSRDS-NBS 67's
`CO + M` 7000 to 15000 K, Visser's photodissociation to about 500 K. The
3000 to 4000 K band is a gap between the combustion literature below it and
the shock-tube literature above it. Venot's pressure floor is the sharper
mismatch still: 0.01 bar against EXHALE's `p_base = 1e-6` bar, four decades
below the validated range.

---

## 2. The ceiling as it is coded

READ, `src/modules/lower_atmosphere/diffusive_photochemistry.f90`, subroutine
`limit_to_element_budget` (the CO block sits at lines 660 to 700 of the file
as read on 2026-09-06; other workers are editing the same tree, so the routine
name is the stable reference).

```text
nco_eq = co_equilibrium_density( nC_free(j), nO_free(j), max(TK(j),1) )
if ( nc(ic_CO) > nco_eq )  nc(ic_CO) = nco_eq
```

`co_equilibrium_density` (READ, `src/modules/lower_atmosphere/oxygen_rates.f90`)
solves `n_C n_O / n_CO = K_c(T)` for the CO root and caps it at
`min(n_C_tot, n_O_tot)`. `K_c` comes from `equilibrium_constant_conc`, which
is `exp(-dG/RT) n_std^dn` with `dG` from the Shomate table of the same module
and `n_std = p_std / (k_B T)`, `p_std = 1 bar`.

The ceiling is therefore a pure thermodynamic statement with no kinetics in
it, and the code says so in its own comment: "The network has no CO kinetics
to destroy it with, so the statement made here is the one the thermodynamics
supports" (READ, same routine).

**Where it fires.** MEASURED here (`co2.py`, from the Shomate table above),
for solar C/H = 2.7e-4 and O/H = 4.9e-4, the temperature at which the
equilibrium CO fraction of the carbon falls through one half:

| n_H [cm^-3] | T at x_CO = 0.5 [K] |
|---|---|
| 1e10 | 2970 |
| 1e12 | 3323 |
| 1e13 | 3533 |
| 1e14 | 3772 |

Below those temperatures the ceiling is inactive and the transported value
stands. So the ceiling's active domain begins at the top of the temperature
range the brief states for the hot Uranus and HD 209458 b lower atmospheres
(300 to 3000 K) and continues upward through the wind.

**CO2 is not in the code.** READ, `src/modules/init/species_table.f90`: the
carriers are `isp_H2`, `isp_OH`, `isp_H2O`, `isp_CO` (indices 34, 38, 39, 40);
there is no CO2 index. Its absence is what removes the single fastest CO
destruction channel of the whole combustion literature from consideration
(section 4.3).

---

## 3. Reaction energies from the code's own thermodynamic table

The brief asks for the reaction energy of each channel from the species
enthalpies of `molecular_reaction_heat.f90`. **That file does not carry them.**
READ, `src/modules/lower_atmosphere/molecular_reaction_heat.f90`, subroutine
`species_enthalpies`: the table is H2, H+, He+, H2+, H3+, HeH+ only, zero at
"ground-state H, ground-state He and free electrons at rest". There is no C,
no O, no CO, no OH and no H2O in it. The oxygen and carbon enthalpies of this
code live in one other place, the Shomate table of `oxygen_rates.f90`
(species `ith_H, ith_H2, ith_O, ith_OH, ith_H2O, ith_CO, ith_C`), which is
also what `water_photolysis.f90` uses for its photolysis thresholds ("The
threshold energies are the reaction enthalpies at 298.15 K computed from the
Shomate table of oxygen_rates", READ).

**This is a gap worth recording**: a CO balance row written as B1 T5.1
specifies it would take its `eps_s` from a table that does not contain its
species, and the two tables have different zeros (`molecular_reaction_heat`
uses ground-state neutral atoms at rest; the Shomate coefficients carry the
NIST absolute enthalpy including the enthalpy of formation, with the eighth
NIST coefficient H dropped for exactly that purpose, READ from the
`oxygen_rates.f90` header). Any implementation of T5.1 must first state which
of the two is the reference and convert the other, or the ledger will not
close.

MEASURED here (`shom.py`), from the Shomate coefficients of `oxygen_rates.f90`,
the reaction enthalpy `dH(T)` in eV per reaction (positive means endothermic):

| reaction | 298 K | 1000 K | 2000 K | 3000 K | 4000 K |
|---|---|---|---|---|---|
| CO + H -> C + OH | +6.718 | +6.710 | +6.688 | +6.683 | +6.694 |
| CO + H2 -> C + H2O | +6.067 | +6.049 | +6.052 | +6.077 | +6.112 |
| CO + M -> C + O + M | +11.156 | +11.236 | +11.305 | +11.361 | +11.419 |
| C + OH -> CO + H | -6.718 | -6.710 | -6.688 | -6.683 | -6.694 |
| C + O -> CO | -11.156 | -11.236 | -11.305 | -11.361 | -11.419 |

Consistency of the table with the bond energies it should reproduce, MEASURED
at 298.15 K: D(CO) = 11.156 eV, D(OH) = 4.438 eV, D(H2) = 4.519 eV. The last
is to be compared with `h2_dissociation_energy_eV()` of `mol_rates.f90`, which
is the code's single definition of D0(H2); the 298 K enthalpy difference is
above the 0 K dissociation energy by the thermal terms, as it must be.

---

## 4. The channels, one by one

### 4.1 CO + H -> C + OH, and its reverse

**The forward rate is published and directly tabulated.**

READ, Westley, F. 1980, *Table of recommended rate constants for chemical
reactions occurring in combustion*, NSRDS-NBS 67, National Bureau of
Standards, doi 10.6028/NBS.NSRDS.67 (downloaded to
`references/Westley_1980_NSRDS-NBS_67.pdf`), entry "CO + H -> C + OH",
source code `75 BEN/GOL`, reaction order 2, in the table's units of
cm^3 mol^-1 s^-1:

```text
A = 2.0e13 ,  n = 0.5 ,  E/R = 77755 K       (no temperature range given)
```

which in cm^3 s^-1 is

```text
k(CO + H -> C + OH) = 3.32e-11 T^0.5 exp(-77755/T)   cm^3 s^-1 .
```

The same page carries the second hydrogen channel, "CO + H -> O + CH",
A = 2.5e13, n = 0.5, E/R = 88020 K, and the reverse, "C + OH -> CO + H",
A = 6.3e11, n = 0.5, E/R = 0, that is

```text
k(C + OH -> CO + H) = 1.05e-12 T^0.5   cm^3 s^-1 .
```

That reverse rate is the one both reference networks use. READ, Tsai et al.
2017, ApJS 228, 20 (`references/Tsai_2017_ApJS_228_20.pdf`), Table 5, entry
R115: `OH + C -> CO + H`, `1.05e-12 T^0.500`, source "NSRDS 67". READ, the
network file in this workspace,
`EXHALE_v1.00/VULCAN/thermo/NCHO_photo_network.txt` line 91, entry 111:
`[ OH + C -> CO + H ] 1.05E-12 0.500 0.0` with the temperature-range column
left as `--`.

Both entries trace to the same origin, Benson, Golden et al. 1975, which in
NSRDS-NBS 67's own classification is a thermochemical-kinetics estimate and
not a measurement. The `E/R = 0` and `n = 0.5` of the reverse are the
signature of a capture rate; the forward `E/R = 77755 K` is the reverse
multiplied through the equilibrium constant of the 6.72 eV endothermicity
(6.718 eV / k_B = 77960 K, MEASURED in section 3, against the tabulated
77755 K, a 0.3 per cent difference).

**Independent check by detailed balance.** MEASURED here (`co2.py`): taking
the reverse rate `1.05e-12 T^0.5` and dividing by `K_c` computed from the
code's own Shomate table (the same construction `rate_from_detailed_balance`
uses in `oxygen_rates.f90`) gives a forward rate that agrees with Westley's
tabulated forward entry to 13 to 19 per cent over 2000 to 5000 K:

| T [K] | k_f by detailed balance | k_f from NSRDS-NBS 67 |
|---|---|---|
| 2000 | 1.71e-26 | 1.94e-26 |
| 2500 | 4.47e-23 | 5.16e-23 |
| 3000 | 8.60e-21 | 1.01e-20 |
| 4000 | 6.39e-18 | 7.59e-18 |
| 5000 | 3.49e-16 | 4.14e-16 |

Two routes that share only the endothermicity, one from a 1980 evaluation and
one from the 1998 NIST-JANAF thermochemistry this code carries, agree to
better than 20 per cent. The magnitude of the rate is therefore not in
question, whatever one thinks of the Benson estimate.

**Products and energy.** Products C + OH; the reaction takes 6.68 to 6.72 eV
out of the gas over 298 to 4000 K (section 3). Written into a one-sided model
this is a cooling term of that size per destruction, not a heating term.

**Validity.** NSRDS-NBS 67 gives no temperature range for the `75 BEN/GOL`
entries, which is how the compilation marks an estimate. KIDA's rule for the
reverse-type entry is explicit and is quoted in section 5.2: below 300 K the
capture form is what is recommended, above it the extrapolation is not. For
the argument of section 7 this does not matter: a factor of a few in a
capture rate cannot move a conclusion that turns on four to ten decades.

### 4.2 CO + H2

There is no bimolecular CO + H2 channel in either reference network. What
exists is the association

```text
CO + H2 + M -> H2CO + M
```

READ, Tsai et al. 2017 Table 5 entry R283, source Wang & Frenklach (1997):

```text
k0   = 2.80e-20 T^-3.420 exp(-42445/T)   cm^6 s^-1
kinf = 7.14e-17 T^1.500  exp(-40055/T)   cm^3 s^-1
```

with no temperature range stated in the table. The same coefficients are in
the workspace network file (line 393, entry 707), also with the range column
blank on both limits.

The exponential is 40000 to 42000 K. At 3000 K it is `exp(-14.1)` at best.
H2CO is not an EXHALE species, so the channel would in any case need two more
carriers (H2CO and HCO) before it could be written down.

The direct abstraction `CO + H2 -> C + H2O` is endothermic by 6.05 to 6.11 eV
(MEASURED, section 3) and appears in neither compilation.

### 4.3 CO + OH -> CO2 + H, and what the absence of CO2 implies

This is the fastest CO destruction channel in the whole combustion
literature, and it is the one channel of this section with a real evaluation
behind it.

READ, Baulch, D. L., Cobos, C. J., Cox, R. A., et al. 1992, JPCRD 21, 411
(`references/Baulch_1992JPCRD.pdf`, the summary table page transcribed as
`OH + CO -> H + CO2`):

```text
k = 1.05e-17 T^1.3 exp(250/T)   cm^3 s^-1 ,   300-2000 K ,
uncertainty +/- 0.2 in log10 at 300 K rising at the top of the range.
```

READ, Tsai et al. 2017 Table 5 entry R39: `OH + CO -> H + CO2`,
`1.05e-17 T^1.500 exp(-259.0/T)`, source "NIST 1992BAU/COB411-429", that is
the same Baulch evaluation. The workspace network file carries the same
numbers (line 54, entry 37, range `300-2000`). Note that the two do not agree on the exponent: 1.3 in the
printed Baulch table against 1.5 in Tsai's. The exponential is the same term
in both once the network file's `k = A T^B exp(-C/T)` convention is applied to
its `C = -250.0` (READ, network entry 37). The Baulch reading here is from
`pdftotext` on a scanned page whose optical character recognition is poor
(the printed `1.05` comes out as `l.Oj` and `1.3` as `1'13`), so the 1.3
should be checked against the page image before it is used; the discrepancy is
recorded rather than resolved. Anyone implementing this channel must resolve
it against the Baulch page rather than against either network.

An older recommendation on the same page of NSRDS-NBS 67 (READ) gives the
same reaction over a wider range, 250 to 2500 K, as
`log10 k [cm^3 mol^-1 s^-1] = 10.83 + 3.94e-4 T`, source `76 BAU/DRY`.

**What its absence implies.** CO2 is not a species of this code (section 2),
so this channel cannot be written without adding a carrier. That is not a
bookkeeping inconvenience: it removes the only CO destruction channel in the
300 to 2500 K range that has an evaluated, measured rate coefficient. Every
other channel available in the ceiling's own regime is either an estimate
(section 4.1), an association into a species the code does not carry
(sections 4.2 and 4.4), or valid only far above the regime (section 4.5). So
a CO destruction model built inside the present species set is built entirely
out of estimates and extrapolations, and section 7 shows that even taken at
face value they are far too slow.

A second consequence: `OH + CO -> H + CO2` is also a net **sink of OH** in the
molecular layer, and OH is a transported carrier here (READ,
`diffusive_photochemistry.f90` section 2 of the header: OH "holds up to 14% of
the oxygen where the carriers live"). Its omission is therefore an omission in
the oxygen budget as well as in the carbon one, not only in the CO row.

### 4.4 CO + O + M -> CO2 + M

READ, Tsai et al. 2017 Table 5 entry R265, source "NIST 1986TSA/HAM1087" for
the low-pressure limit and Simonaitis & Heicklen (1972) for the high-pressure
limit:

```text
k0   = 1.70e-33 exp(-1510/T)   cm^6 s^-1
kinf = 2.66e-14 exp(-1459/T)   cm^3 s^-1
```

The workspace network file carries the same coefficients (line 384, entry
689) with the range `300-2500` on `k0`. Products CO2 + M; again CO2 is not a
species of this code. The reverse (`CO2 + M -> CO + O + M`) is a CO source,
not a sink.

### 4.5 CO + M -> C + O + M, thermal dissociation

READ, NSRDS-NBS 67 entry "CO + M -> C + O + M", source `76 BAU/DRY`:

```text
k = 8.8e29 T^-3.5 exp(-128700/T)   cm^3 mol^-1 s^-1 ,
range 7000-15000 K , M = Ar or CO , uncertainty factors f = 0.3, F = 1.8 .
```

that is `1.46e6 T^-3.5 exp(-128700/T)` cm^3 s^-1.

**The stated range starts at 7000 K.** This is the cleanest single statement
in this document: the evaluated literature does not offer a thermal
dissociation rate for CO anywhere inside the ceiling's regime, because there
is nothing to evaluate there. The activation temperature, 128700 K, is the
11.16 eV bond of section 3 (11.156 eV / k_B = 129470 K, MEASURED, 0.6 per
cent from the tabulated value). Extrapolated down to 3000 K the coefficient is
of order 2e-25 cm^3 s^-1, giving a destruction time of 4e11 s at
n = 1e13 cm^-3, and the extrapolation is flagged here only to show how far it
is from mattering.

The formation counterpart is in both networks: READ, Tsai et al. 2017 Table 5
entry R289, `O + C + M -> CO + M`, `k0 = 9.10e-22 T^-3.100 exp(-2114/T)`,
source "NSRDS 67"; the same in the workspace file (line 428, entry 767).

### 4.6 The ion channels: CO + H+ and CO + He+

These are the channels that actually destroy CO in the region where the
ceiling fires, and they are the ones for which no source could be obtained
here. Both compilations that carry them, UMIST RATE12 and KIDA, refused
their data files to this machine (section 8).

What can be stated without a source, from data already in the code:

- **CO + H+ -> CO+ + H is endothermic.** The first ionization energy of CO is
  above that of H, so the charge transfer has a threshold. **Confirmed from
  RATE22 in section 12.3**, which carries only the exothermic direction,
  `H + CO+ -> CO + H+`. The guess offered in the first draft of this
  paragraph, that the networks carry a radiative association
  `H+ + CO -> HCO+ + hv` instead, is **wrong**: RATE22 has no such entry
  either. There is no `H+ + CO` destruction channel at all.
- **CO + He+ -> C+ + O + He is exothermic and is the standard fast destroyer**
  in the astrochemical networks, but He+ requires 24.6 eV photons and is
  absent from the molecular layer entirely.

Neither is quoted with a number here. Both are listed in section 8 as
requests, and section 7's conclusion does not use them: they can only make
destruction faster in the hot region, where section 7 already finds the
ceiling defensible, and they are zero in the cool region, where section 7
finds it is not.

### 4.7 CO + e

Dissociative recombination applies to CO+, not to CO. The neutral channel is
electron-impact dissociation, whose threshold is the 11.16 eV bond of section
3. At 3000 K the mean electron energy is 0.39 eV; the fraction of a Maxwellian
above 11.16 eV is `exp(-43)` to within the usual prefactors. No rate is needed
to conclude that this channel is closed in the ceiling's regime, and none is
quoted.

### 4.8 Photodissociation, and the shielding functions

CO photodissociates entirely through predissociating **lines** between roughly
912 and 1118 A, not through a continuum, which is why it self-shields and why
the shielding cannot be written as `exp(-sigma N)`.

The code already records this. READ,
`src/modules/lower_atmosphere/water_photolysis.f90` section 3:

> CO absorbs below 1118 A (in the LW band) in predissociating lines that
> self-shield. It is not carried: CO is chemically frozen in this network
> (decision D4) so its photodissociation would have no consumer, and adding
> it as an absorber alone would remove photons from the water cycle without
> the compensating oxygen release. A network that lets CO react must add
> both.

The band in question is the code's own LW band, 912 to 1110 A (READ,
`oxygen_rates.f90`, `fuv_band_lo_A`/`fuv_band_hi_A`), which already carries
the H2 Lyman-Werner flux and the H2O and OH continua. So a CO photodissociation
term is not an addition to a free band, it is a fourth absorber sharing an
existing beam, and the "one band, one incident flux, one beam" rule the same
header states would apply to it.

**The standard reference is Visser, van Dishoeck & Black 2009, A&A 503, 323.**
The user supplied the publisher PDF on 2026-09-06
(`references/Visser_2009A&A_503_323.pdf`), and **section 10 below is the
reading of it**: the unshielded rate and its field, the 37-band line list, the
shielding-function tables and their parameters, the model's stated ceiling near
500 K, and how the function would enter the FUV band treatment of this code.

**What the workspace has instead, and why it is not enough.** READ,
`EXHALE_v1.00/VULCAN/thermo/photo_cross/CO/CO_cross.csv`: a three-column
absorption / dissociation / ionization cross section on a wavelength grid
running 6.199 to 166.499 nm in 1604 points. MEASURED from that file: the grid
spacing through the CO band region is 0.1 nm uniformly (250 points between 90
and 115 nm) and the peak dissociation cross section is 2.599e-16 cm^2 at
107.599 nm. A 0.1 nm grid is three to four orders of magnitude coarser than
the CO band lines, so this table is a band-averaged continuum and carries no
shielding information at all. Consistent with that, MEASURED by grep: there is
no self-shielding term anywhere in the VULCAN Python source in this workspace
(no match for "shield" in `VULCAN/*.py`). The accompanying
`CO/note.txt` reads "including ionization/photodissociation <22.8 nm from
PhiDRATES" (READ), and `CO_branch.csv` gives the single branch CO -> C + O
with ratio 1 at all wavelengths (READ).

So the photodissociation channel cannot be implemented from anything in this
workspace. It needs the Visser 2009 shielding functions, which is what makes
that request the one item of section 8 that actually blocks a piece of
alternative A.

---

## 5. The compilations

### 5.1 UMIST RATE12 (McElroy et al. 2013, A&A 550, A36)

The user supplied the publisher PDF on 2026-09-06
(`references/McElroy_2013A&A_550_A36.pdf`); **section 11.1 below is the
reading of it**. In short: the paper describes the database and its format and
prints no reaction entries at all, so the individual CO rates still have to
come from the network file, which is not in this workspace and remains
request 2 of section 8.

### 5.2 KIDA (Wakelam et al. 2012, ApJS 199, 21)

Downloaded, `references/Wakelam_2012_ApJS_199_21.pdf`. Two statements from
the published text bear directly on the question, quoted verbatim:

> Each year, a subset of the reactions in the database (kida.uva) will be
> provided as a network for the simulation of the chemistry of dense
> interstellar clouds with temperatures between 10 K and 300 K.

> By default, the rate coefficients are valid between 10 and 300 K. For
> updated rate coefficients, the temperature range is given by the results of
> experiments or the estimates of experts. Extrapolations outside this range
> of temperature are not recommended unless the rate coefficient is predicted
> to be totally independent of temperature, as occurs for the Langevin model
> of exothermic ion-non-polar-neutral reactions.

The rate form is the Arrhenius-Kooij `k = A (T/300)^B exp(-C/T)` (READ, their
equation 1), and photodissociation is stored as `k = A exp(-C A_v)` (their
equation 5) in the Draine (1978) interstellar field, which is not a form an
exoplanet code can use without recomputing `A` for the stellar spectrum.

**So KIDA answers the brief's question negatively**: it does not give CO
destruction rates valid to 3000 K, and it says in its own text that
extrapolating its entries there is not recommended, with an exception
(Langevin ion-neutral rates) that covers exactly the two ion channels of
section 4.6 and none of the neutral ones.

### 5.3 Venot et al. 2012 (A&A 546, A43) and Venot et al. 2020 (A&A 634, A78)

The user supplied both publisher PDFs on 2026-09-06
(`references/Venot_2012A&A_546_A43.pdf`,
`references/Venot_2020A&A_634_A78.pdf`); **sections 11.2 and 11.3 below are
the reading of them**. The answer to the brief's question is no on both
counts: the stated range is 300 to 2500 K and 0.01 to 100 bar, and the 2020
update does not extend it.

### 5.4 Moses et al. 2011 (ApJ 737, 15)

The user supplied the publisher PDF on 2026-09-06
(`references/Moses_2011_ApJ_737_15.pdf`); **section 13.1 below is the reading
of it**. Nothing was ever quoted from the arXiv copy, so nothing is withdrawn.
The short answer: the paper states no single validity envelope for its
network, its one quoted fitting interval is 500 to 2500 K, and it carries no
CO self-shielding.

### 5.5 Tsai et al. 2017 (VULCAN, ApJS 228, 20)

Downloaded, `references/Tsai_2017_ApJS_228_20.pdf`. The answer to the brief's
question is in the abstract, quoted verbatim:

> It is constructed for gaseous chemistry from 500 to 2500 K, using a reduced
> C-H-O chemical network with about 300 reactions.

and in section 4:

> most of the rate coefficients obtained from the NIST database (validated
> from 500 to 2500 K)

**So VULCAN's own stated validity stops at 2500 K, not 3000 K.** Its reverse
rates are not fitted but computed on the fly:

> Having obtained the reversed rate coefficients, we do not fit them with the
> functional form of the generalized Arrhenius equation, because we find that
> the fitting procedure can fail or produce noticeable errors ... Instead, we
> reverse the rate coefficients, at given values of the temperature and
> pressure, on the fly and before the calculation starts.

This is the same construction `rate_from_detailed_balance` implements in
`oxygen_rates.f90`, so a CO row written for EXHALE would reverse its rates the
same way VULCAN does, and the cross-check of section 4.1 shows the two agree.

### 5.6 Tsai et al. 2021 (ApJ 923, 264)

The user supplied the publisher PDF on 2026-09-06
(`references/Tsai_2021_ApJ_923_264.pdf`); **section 13.2 below is the reading
of it**. The answer to the question this section originally left open is no:
its photochemistry is continuum Beer-Lambert on a constant 0.1 nm wavelength
grid, which cannot represent CO line self-shielding at any column, and its
temperature-dependent cross sections cover eight molecules of which CO is not
one.

### 5.7 What the VULCAN network in this workspace actually carries

READ, `EXHALE_v1.00/VULCAN/thermo/NCHO_photo_network.txt`, which is the file
`vulcan_cfg.py` names (`network = 'thermo/NCHO_photo_network.txt'`, READ).
Every entry in which CO appears as a reactant, with the file's own source and
temperature-range columns:

| id | reaction | A | B | C | source column | T range column |
|---|---|---|---|---|---|---|
| 37 | OH + CO -> H + CO2 | 1.05e-17 | 1.500 | -250.0 | (blank) | 300-2000 |
| 153 | CH3O + CO -> CH3 + CO2 | 2.61e-11 | 0.000 | 5940.0 | (blank) | 300-2500 |
| 499 | CO + HO2 -> CO2 + OH | 2.51e-10 | 0. | 11900.0 | 1986TSA/HAM1087 | 300-2500 |
| 625 | O_1 + CO -> CO2 | 8.00e-11 | 0.000 | 0.0 | 1975TUL1893 | 100-2100 |
| 637 | CH2_1 + CO -> CH2 + CO | 5.00e-11 | 0.000 | 0.0 | Wang&Frenklach(1997) | (blank) |
| 689 | CO + O + M -> CO2 + M | k0 1.70e-33 | 0.000 | 1510.0 | (blank) | 300-2500 |
| 693 | H + CO + M -> HCO + M | k0 5.29e-34 | 0.000 | 370.0 | (blank) | 300-2500 |
| 703 | CO + CH3 + M -> CH3CO + M | k0 3.95e-10 | -7.500 | 5490.0 | (blank) | 300-1700 |
| 707 | CO + H2 + M -> H2CO + M | k0 2.80e-20 | -3.420 | 42445.0 | (blank) | (blank) |
| 799 | CO -> C + O (photo) | (cross section file `CO`, branch 1) | | | | |

Three observations follow directly from that table.

1. **Every product of every bimolecular entry is a species EXHALE does not
   carry**: CO2, CH3, HO2, HCO, CH3CO, H2CO. The only entry whose products
   are inside the code's species set is the photodissociation, entry 799.
2. **There is no `CO + H` entry and no `CO + M -> C + O + M` entry.** The
   `CO + H -> C + OH` reaction of the brief exists in the network only as
   the thermodynamic reverse of entry 111 (`OH + C -> CO + H`), generated on
   the fly by the mechanism section 5.5 describes, because the network file
   lists forward reactions with odd ids and VULCAN builds the even-id
   reverses down to the "# reverse stops" line (READ, the file's own comment
   at line 439).
3. **The network's carbon chemistry is a hydrocarbon chemistry.** CO is
   destroyed in it by being carried into CH3, HCO, H2CO and CH3CO, all of
   which return it. Removing those species, as EXHALE has, does not leave a
   reduced CO destruction network. It leaves none.

---

## 6. Summary table of usable rates

"Usable" means: the products are species EXHALE carries, or the missing
species is named. Rates in cm^3 s^-1 (two-body) or cm^6 s^-1 (three-body).

| channel | rate coefficient | T range | products | dH [eV] | source |
|---|---|---|---|---|---|
| CO + H -> C + OH | 3.32e-11 T^0.5 exp(-77755/T) | none stated (estimate) | C, OH | +6.68 to +6.72 | NSRDS-NBS 67, entry `75 BEN/GOL`; local PDF |
| CO + H -> O + CH | 4.15e-11 T^0.5 exp(-88020/T) | none stated (estimate) | O, CH (CH absent) | not evaluated, CH not in the Shomate table | NSRDS-NBS 67, `75 BEN/GOL` |
| C + OH -> CO + H | 1.05e-12 T^0.5 | none stated (estimate) | CO, H | -6.68 to -6.72 | NSRDS-NBS 67; Tsai 2017 R115; network entry 111 |
| CO + OH -> CO2 + H | 1.05e-17 T^1.3 exp(250/T), exponent to be checked against the page image (section 4.3) | 300-2000 K | CO2 (absent), H | not evaluated, CO2 not in the Shomate table | Baulch et al. 1992 JPCRD 21, 411; local PDF |
| CO + OH -> CO2 + H | log10 k [cm^3/mol/s] = 10.83 + 3.94e-4 T | 250-2500 K | CO2 (absent), H | as above | NSRDS-NBS 67, `76 BAU/DRY` |
| CO + O + M -> CO2 + M | k0 = 1.70e-33 exp(-1510/T) | 300-2500 K | CO2 (absent) | as above | Tsai 2017 R265; network entry 689 |
| CO + H2 + M -> H2CO + M | k0 = 2.80e-20 T^-3.42 exp(-42445/T) | none stated | H2CO (absent) | not evaluated | Wang & Frenklach 1997 via Tsai 2017 R283 |
| CO + M -> C + O + M | 1.46e6 T^-3.5 exp(-128700/T) | **7000-15000 K** | C, O | +11.16 to +11.42 | NSRDS-NBS 67, `76 BAU/DRY`; local PDF |
| C + O + M -> CO + M | k0 = 9.10e-22 T^-3.1 exp(-2114/T) | none stated | CO | -11.16 to -11.42 | NSRDS-NBS 67 via Tsai 2017 R289 |
| CO + hv -> C + O | `k = chi k0 Theta(N_CO, N_H2) exp(-gamma A_V)`, `k0 = 2.6e-10 s^-1` in the Draine field; `Theta` from Tables 5 to 8 | bands 911.75-1117.80 A; `T_ex(CO)` to 512 K only | C, O | +11.16 threshold | Visser et al. 2009, A&A 503, 323; local PDF (section 10) |
| CO + hv -> C + O | `k = alpha exp(-gamma A_V)`, `alpha = 2.40e-10`, `gamma = 3.88`, unshielded | 10-41000 K (default bound) | C, O | as above | UMIST RATE22 entry 8259, `C`/`C`, after Heays et al. 2017 (section 12.3) |
| CO + He+ -> C+ + O + He | 1.60e-09, `beta = gamma = 0` (Langevin) | 10-41000 K (default bound); temperature independent | **C+, O, He, all carried** | exothermic; He+ recombination energy 24.59 eV less the 11.16 eV bond and the C ionization | UMIST RATE22 entry 4068, `M`/`A` (section 12.3) |
| CO + H+ | **no entry exists**; RATE22 carries only the exothermic reverse `H + CO+ -> CO + H+` | | | endothermic | UMIST RATE22 entry 499 (section 12.3) |
| CO + H3+ -> HCO+ + H2 | 1.36e-09 `(T/300)^-0.14 exp(3.4/T)` | 10-400 K only | no (HCO+) | | UMIST RATE22 entry 3497, `C`/`A` |
| CO + e | threshold 11.16 eV, closed below 5000 K | | | | no rate needed |

---

## 7. The timescale assessment

### 7.1 The two inequalities are not equally hard

B1 T5.2 requires `tau_dest << tau_res << tau_form`. The second half of that is
easier than it looks. By detailed balance, at any state,

```text
tau_dest / tau_form  =  P_form / L_dest  =  K_eq / Q ,
```

with `Q` the local reaction quotient. The ceiling fires precisely when the
parcel is CO-supersaturated, `Q > K_eq`, so `tau_dest < tau_form` holds **by
construction wherever the ceiling is active**. That half of the ordering is
never the obstacle and does not need rate data at all.

The obstacle is entirely `tau_dest << tau_res`.

### 7.2 tau_dest

MEASURED here (`co2.py`), `tau_dest = 1 / (k_f n_H)` for the CO + H channel of
section 4.1, using the detailed-balance rate (the NSRDS-NBS 67 tabulated
forward rate gives values 13 to 19 per cent shorter, which changes nothing
below):

| T [K] | n_H = 1e10 | 1e12 | 1e13 | 1e14 |
|---|---|---|---|---|
| 1000 | 6.4e32 | 6.4e30 | 6.4e29 | 6.4e28 |
| 2000 | 5.9e15 | 5.9e13 | 5.9e12 | 5.9e11 |
| 2500 | 2.2e12 | 2.2e10 | 2.2e9 | 2.2e8 |
| 3000 | 1.2e10 | 1.2e8 | 1.2e7 | 1.2e6 |
| 3500 | 2.7e8 | 2.7e6 | 2.7e5 | 2.7e4 |
| 4000 | 1.6e7 | 1.6e5 | 1.6e4 | 1.6e3 |
| 5000 | 2.9e5 | 2.9e3 | 2.9e2 | 2.9e1 |
| 6000 | 1.9e4 | 1.9e2 | 1.9e1 | 1.9 |
| 8000 | 6.3e2 | 6.3 | 6.3e-1 | 6.3e-2 |

(seconds). The other neutral channels are slower still: the thermal
dissociation of section 4.5 gives 4e11 s at 3000 K and 1e13 cm^-3, and the
CO2-forming channels do not exist in this species set.

### 7.3 tau_res

READ, `diffusive_photochemistry.f90` header section 1, measured in an EXHALE
run of HD 189733 b: the base cell has `v = 0` exactly by the lower boundary
condition; the flow time `r/|v|` in the cells above it is 1e8 to 2e8 s; the
diffusive time of the base cell is 8.0e5 s with `K_zz = 0` and 634 s with
`K_zz = 1e9`.

Taking `tau_res = min(L/|v|, L^2/(D + K_zz))` as B1 T5.2 defines it, the
residence time in the molecular layer of that run is **634 s to 8.0e5 s**, set
by diffusion, with the flow time two to five decades longer. Higher in the
wind the flow time falls, so the bracket for the whole region is roughly
1e2 to 1e8 s.

### 7.4 Comparison, and where the ordering fails

Put the two together at the temperatures where the ceiling actually fires
(section 2: onset at 2970 K at n_H = 1e10 rising to 3772 K at n_H = 1e14):

| n_H [cm^-3] | T at ceiling onset [K] | tau_dest there [s] | tau_res [s] | tau_dest / tau_res |
|---|---|---|---|---|
| 1e10 | 2970 | 1.5e10 | 1e2 to 1e8 | 1.5e2 to 1.5e8 |
| 1e12 | 3323 | 8.9e6 | 1e2 to 1e8 | 1e-1 to 9e4 |
| 1e13 | 3533 | 2.2e5 | 6e2 to 8e5 | 0.3 to 4e2 |
| 1e14 | 3772 | 5.2e3 | 6e2 to 8e5 | 6e-3 to 9 |

**The ordering `tau_dest << tau_res` fails at the ceiling's onset at every
density**, marginally at the base (where the high density shortens `tau_dest`
and `K_zz` shortens `tau_res` by comparable amounts) and catastrophically in
the rarefied gas, where the ceiling is nonetheless applied to every cell.

Deeper inside the stated regime it fails by more. At 2000 K and
n_H = 1e13 cm^-3, `tau_dest = 5.9e12 s` against `tau_res` of at most 8e5 s: a
factor of 7e6. At 1000 K it is 24 decades.

The ordering only recovers above roughly 5000 K at n_H = 1e12, or 4000 K at
1e14: `tau_dest` there is 1e2 to 1e3 s and comparable to the shortest
residence times. Above 6000 K it is fast by any measure, which is the regime
the code's own comment appeals to ("the collisional and radiative processes
that take it apart at 2e4 K run far faster than the flow", READ) and where
the ion and photon channels of sections 4.6 and 4.8 take over anyway.

### 7.5 What this settles

> **Point 1 below is superseded by section 12.6.** It was written before the
> UMIST RATE22 ratefile was available, from the neutral channels alone. With
> the He+ channel of section 12.3 the ordering does hold above the helium
> ionization front, and the domain argument T5.2 returns a domain with a lower
> boundary rather than an empty one. Points 2 to 4 stand as written.


1. **Alternative A as B1 T5.1 writes it, a one-sided destruction model with a
   domain argument, cannot be constructed for the ceiling's regime.** Not for
   want of rate data: the rates exist, they were found, and they say the
   destruction is slow. The domain argument T5.2 asks for would return an
   empty domain below about 4000 K.
2. **The present ceiling is not a kinetic statement where it first fires.**
   Between roughly 2500 and 4000 K it imposes instantaneous destruction on a
   channel whose actual time is 1e4 to 1e10 s against a residence time of at
   most 1e5 to 1e8 s. The code's own comment already says the ceiling
   over-suppresses "around the 3000-4000 K turnover"; the numbers above are
   the width and the size of the error it names, and they are larger than
   "the width of the interval it gets wrong cannot be measured without the
   rate the audited set does not contain" suggests. The rate exists
   (section 4.1) and the interval is measurable.
3. **B1 decision 7 (adopt exclusion) is confirmed on physics.** A run in
   which the ceiling fires during accepted physical integration is a run in
   which CO has been removed on a timescale the kinetics does not support, so
   excluding it from validated results is the correct treatment, and it is
   correct independently of whether the rate compilations are ever obtained.
4. **If a CO row is ever written, it must be two-sided, not one-sided.** The
   quenched CO the ceiling destroys is real: a parcel carried out of the
   molecular layer keeps its CO for 1e5 to 1e10 s, which is longer than it
   takes to reach the wind. The physics that eventually removes it is
   photodissociation in the 912 to 1118 A band and charge transfer with He+,
   both of which need the sources of section 8 and the second of which needs
   He+ to be available in the same cell. Neither is a neutral thermal channel,
   so the reaction set of section 6 is not the set such a row would use.

---

## 8. Requests for the user

**All but one item is now supplied.** Six of the seven were provided by the
user on 2026-09-06, and the seventh is not blocking. The list is kept, with
each item marked, so that the record of what was needed and what it answered
stays together.

1. **SUPPLIED.** Visser, R., van Dishoeck, E. F., & Black, J. H. 2009,
   A&A 503, 323, bibcode 2009A&A...503..323V,
   doi 10.1051/0004-6361/200912129, at
   `references/Visser_2009A&A_503_323.pdf`. Read in **section 10**: unshielded
   rate `2.6e-10 s^-1` in the Draine (1978) field, 37 bands and 855 lines per
   isotopologue over 911.75 to 1117.80 A, the `Theta[N(12CO), N(H2)]` tables
   with their four parameter sets, and the statement that the model is not
   extended above `T_ex(CO) = 512 K` for want of data on the `v'' = 1` level.

2. **SUPPLIED, and by the later release.** The RATE12 data file was requested;
   the user pointed instead to `~/RT_Codes/UMIST/` with the **RATE22**
   ratefile (`rate22_final.rates`, 8767 entries) and its paper, Millar, T. J.,
   Walsh, C., Van de Sande, M., & Markwick, A. J. 2024, A&A 682, A109,
   bibcode 2024A&A...682A.109M, doi 10.1051/0004-6361/202346908, at
   `~/RT_Codes/UMIST/UDfA2024.pdf`. Read in **section 12**, which quotes every
   CO reactant entry. This is the item that changed the verdict: entry 4068,
   `He+ + CO -> O + C+ + He`. The RATE12 paper itself was also supplied
   (`references/McElroy_2013A&A_550_A36.pdf`, section 11.1), and it prints no
   reaction entries.

3. **SUPPLIED.** Venot, O., Hebrard, E., Agundez, M., et al. 2012, A&A 546,
   A43, at `references/Venot_2012A&A_546_A43.pdf`. Read in **section 11.2**:
   300 to 2500 K and 0.01 to 100 bar, CO photodissociation as entry J9 with no
   shielding, thermodynamic reversal with three named exceptions.

4. **SUPPLIED.** Venot, O., Bounaceur, R., Dobrijevic, M., et al. 2020,
   A&A 634, A78, at `references/Venot_2020A&A_634_A78.pdf`. Read in
   **section 11.3**: a methanol update of V12 with no new validity range, its
   own experimental base at 800 to 1700 K and 1 to 50 bar.

5. **SUPPLIED.** Moses, J. I., Visscher, C., Fortney, J. J., et al. 2011,
   ApJ 737, 15, at `references/Moses_2011_ApJ_737_15.pdf`. Read in
   **section 13.1**: no stated validity envelope, one quoted fitting interval
   of 500 to 2500 K, a grid running 1000 to 1e-6 bar, and no CO
   self-shielding (the string "shield" occurs once in the paper, about NNH).

6. **SUPPLIED.** Tsai, S.-M., Malik, M., Kitzmann, D., et al. 2021,
   ApJ 923, 264, at `references/Tsai_2021_ApJ_923_264.pdf`. Read in
   **section 13.2**: the answer is no. VULCAN's photochemistry is continuum
   Beer-Lambert on a constant 0.1 nm wavelength grid, which cannot represent
   CO line self-shielding, and its temperature-dependent cross sections cover
   eight molecules of which CO is not one.

7. **CLOSED without the book (2026-09-06; the user could not obtain volume 3).** The two recommendations it covers stay here as second-hand citations through NSRDS-NBS 67 (Westley 1980), which prints the 1976 recommendation and its temperature range; neither enters the revised verdict (`CO + M` is ruled out by its barrier at any transcription; `CO + OH` needs CO2, absent from the species table, and is superseded by Baulch et al. 1992). A file uploaded as `Baulch_1976_Vol3_...pdf` is byte-identical to `Baulch_1981.pdf` (volume 4) and is not volume 3. Original entry: Baulch, D. L., Drysdale, D. D., Duxbury, J.,
   & Grant, S. J. 1976, *Evaluated Kinetic Data for High Temperature
   Reactions, Volume 3* (London: Butterworths). The `76 BAU/DRY` primary
   source behind the `CO + M -> C + O + M` and `CO + OH -> CO2 + H`
   recommendations quoted at second hand from NSRDS-NBS 67 in sections 4.3 and
   4.5. Volumes 1 (1972) and 4 (1981) are already in `references/`; volume 3
   is not. **Confirm**: the stated lower temperature limit of the `CO + M`
   recommendation and whether the evaluation says anything about extrapolating
   below it. **Why it does not block anything**: the NSRDS-NBS 67
   transcription is a published recommendation in its own right, both
   reactions it covers are ruled out for EXHALE on other grounds (`CO + M` by
   its 11.16 eV barrier, `CO + OH` by the absence of CO2 from the species
   table), and neither enters the revised verdict of section 12.6.

**One further item is not a literature request but a measurement**, recorded
here because it is what would close the last open point of section 12.5: the
`tau_dest` against `tau_res` profile of section 12.5 evaluated on a
**converged** run that reaches 3000 K, so that the nesting of the helium
ionization front inside the CO ceiling's firing region is measured rather than
inferred. The `oxygen_chemistry` case used there is a step-pinned relaxation
snapshot whose physical cells never exceed 1456 K.

---

## 9. Files

*(Sections 10 and 11 below were added after the user supplied four of the
requested PDFs; the lists here cover the whole document.)*

**Downloaded into `references/` by this item** (publisher-served, verified as
PDF):

- `Tsai_2017_ApJS_228_20.pdf` (1 965 219 bytes), IOP, ApJS 228, 20.
- `Wakelam_2012_ApJS_199_21.pdf` (475 760 bytes), IOP, ApJS 199, 21.
- `Westley_1980_NSRDS-NBS_67.pdf` (9 136 525 bytes), NIST
  (`nvlpubs.nist.gov`), NSRDS-NBS 67, doi 10.6028/NBS.NSRDS.67. This is the
  source of the CO + H and C + OH rates of section 4.1 and the CO + M rate of
  section 4.5, and it is cited at second hand as "NSRDS 67" by Tsai et al.
  2017 and by the VULCAN network file.

**Supplied by the user on 2026-09-06 and read here** (publisher PDFs, in
`references/`):

- `Visser_2009A&A_503_323.pdf`, A&A 503, 323. Section 10.
- `McElroy_2013A&A_550_A36.pdf`, A&A 550, A36. Section 11.1.
- `Venot_2012A&A_546_A43.pdf`, A&A 546, A43. Section 11.2.
- `Venot_2020A&A_634_A78.pdf`, A&A 634, A78. Section 11.3.
- `Moses_2011_ApJ_737_15.pdf`, ApJ 737, 15. Section 13.1.
- `Tsai_2021_ApJ_923_264.pdf`, ApJ 923, 264. Section 13.2.

**Pointed to by the user on 2026-09-06**, outside `references/`, at
`~/RT_Codes/UMIST/` (source `https://umistdatabase.uk/`, from
`where_UMIST.txt`):

- `UDfA2024.pdf`, Millar et al. 2024, A&A 682, A109, the published RATE22
  paper. Section 12.1 and 12.2.
- `rate22_final.rates` (8767 entries) and `rate22_G_final.rates` (the reduced
  ratefile, 51 entries zeroed). Sections 12.3 and 12.4.
- `rate22_dipole.specs` and `rate22_revised_CtoO_0.44.specs`, not used here.

**Already in `references/` and used here**:

- `Baulch_1992JPCRD.pdf` (JPCRD 21, 411), for CO + OH -> CO2 + H.
- `Moses_2011_ApJ_737_15_arXiv.pdf`, superseded by the publisher PDF above;
  nothing was ever quoted from it (sections 5.4 and 13.1).

**Read in the code**:

- `src/modules/lower_atmosphere/diffusive_photochemistry.f90` (the ceiling,
  the carrier set, the measured chemical and transport times).
- `src/modules/lower_atmosphere/oxygen_rates.f90` (the Shomate table,
  `co_equilibrium_density`, `equilibrium_constant_conc`,
  `rate_from_detailed_balance`, the FUV bands).
- `src/modules/lower_atmosphere/molecular_reaction_heat.f90` (the H/He
  enthalpy table, and the finding of section 3 that it carries no C, O or CO).
- `src/modules/lower_atmosphere/water_photolysis.f90` (the CO absorber note,
  and the continuum-absorber pattern of section 10.6).
- `src/modules/lower_atmosphere/lyman_werner.f90` (the band, the mean photon
  energy, the shielded-rate assembly and the cell-mean integration that
  section 10.6 takes as the pattern for a CO term).
- `src/modules/init/species_table.f90` (the carrier indices, and the absence
  of CO2).
- `VULCAN/vulcan_cfg.py`, `VULCAN/thermo/NCHO_photo_network.txt`,
  `VULCAN/thermo/photo_cross/CO/*`.

**Run outputs read** (existing files, nothing re-run):

- `backup/regression/mol_base_handoff/input.inp` and `base.inp` (the spectrum
  range and `p_base` of section 10.7).
- `backup/regression/mol_lyman_werner/input.inp` (`Stellar LW flux: 343.0`).
- `backup/regression/oxygen_chemistry/output/Oxygen_chemistry.txt`,
  `Lyman_Werner.txt` and `FUV_bands.txt` (the CO and H2 densities, the
  star-ward H2 column, the band fluxes and the dayside dilution of
  section 10.7).

**Evaluation scripts** (scratch, not part of the repository):
`shom.py` and `co2.py` under the session scratchpad. They transcribe the
`c_shom`, `t_shom` and `n_shom` arrays of `oxygen_rates.f90` and reimplement
`gibbs_energy_shomate`, `equilibrium_constant_conc` and
`co_equilibrium_density` in Python; every table labeled MEASURED above comes
from them, together with the CO cross-section and column integrations of
section 10.7.

---

## 10. CO photodissociation as Visser, van Dishoeck and Black (2009) define it

The publisher PDF was supplied by the user on 2026-09-06 and is at
`references/Visser_2009A&A_503_323.pdf` (A&A 503, 323-343). Everything in
this section is READ from it, quoted by section, page or table; the page
numbers are the journal's own running numbers.

### 10.1 The unshielded rate and the radiation field it belongs to

READ, section 3.2 (p. 328):

> We obtain an unshielded CO photodissociation rate of 2.6 x 10^-10 s^-1.
> This rate is 30% higher than that of vDB88, due to the generally larger
> oscillator strengths in our data set.

and, in the same paragraph, the field dependence stated explicitly:

> Clearly, the rate depends on the choice of radiation field. If we adopt
> Habing (1968), Gondhalekar et al. (1980) or Mathis et al. (1983) instead of
> Draine (1978), the photodissociation rate becomes 2.0, 2.0 or
> 2.3 x 10^-10 s^-1, respectively.

The standard field is named in section 3.1 (p. 328): "We adopt Draine (1978)
as our standard unattenuated interstellar radiation field." The four tables of
shielding functions carry the rate again to four digits at the head of each
isotopologue block; for 12CO it is `k0 = 2.592e-10 s^-1` in Table 5,
`2.590e-10` in Table 6, `2.588e-10` in Table 7 and `2.592e-10` in Table 8
(READ). So `k0` itself is essentially independent of the excitation
temperature and the Doppler width, and all of the parameter dependence sits in
the shielding function.

The estimated accuracy is stated in section 3.5 (p. 330): "we estimate the
absolute photodissociation rates to be accurate to about 20%", and the
relative accuracy between isotopologues "about 10%".

**The geometry statement matters for a planet.** READ, section 3.3 (p. 328):

> Equation (2) assumes the radiation is coming from all directions. If this is
> not the case, such as for a cloud irradiated only from one side, k0,i should
> be reduced accordingly.

EXHALE already applies exactly that reduction to its own band: the
`Lyman_Werner.txt` header of the `oxygen_chemistry` regression case reads
"dayside dilution applied: 0.500" (READ), turning a 343.0 erg cm^-2 s^-1 flux
at the planet into a 171.5 erg cm^-2 s^-1 beam. A CO term would take the same
factor from the same place.

### 10.2 Band structure and line list

READ, section 3.1 (p. 328):

> For each of our 37 CO bands, we include all lines originating from the first
> ten rotational levels (J'' = 0-9) of the v'' = 0 level of the electronic
> ground state. That results in 855 lines per isotopologue. In addition, we
> have 48 H lines and 444 H2 lines, for a total of 5622. We use an adaptive
> wavelength grid that resolves all lines without wasting computational time
> on empty regions. For typical model parameters, the wavelength range from
> 911.75 to 1117.80 Angstrom is divided into ~47 000 steps.

**37 bands, not 33.** Table 1 (pp. 325) is numbered 1 to 33 following vDB88,
but its footnote b reads "The numbering follows vDB88. Their bands 1 and 2 are
split into four and two components", so the table has 37 rows: 1A to 1D, 2A,
2B and 3 to 33. The shortest is band 1A at 912.7037 A and the longest is
band 33, the E0 transition, at 1076.0796 A (READ, Table 1). The E0 band is the
strongest contributor at the cloud edge (section 3.3, p. 329).

The H2 and H shielding lines are named in section 2.5 (p. 328): "we include H
Lyman lines up to n = 50 and H2 Lyman and Werner lines (transitions to the
B 1Sigma+u and C 1Piu states) from the v'' = 0, J'' = 0-7 levels of the
electronic ground state", with data from Abgrall et al. (1993a,b) as compiled
for the Meudon PDR code.

**The integration range 911.75 to 1117.80 A is EXHALE's LW band plus 8 A.**
The code's own bands are `fuv_band_lo_A = (912, 1110, 1202, 1231, 1451)` and
`fuv_band_hi_A = (1110, 1201, 1230, 1450, 2304)` (READ, `oxygen_rates.f90`), so
Visser's range covers the whole LW band and the first 7.8 A of B1. The
`lyman_werner.f90` header already names 1110 A as the end of the H2
Lyman-Werner system, following Draine & Bertoldi (1996), and Visser's longest
CO band sits at 1076 A, inside it. A CO photodissociation term therefore
belongs in the existing LW band with the existing beam, and the "one band, one
incident flux, one beam" rule of the `oxygen_rates.f90` header applies to it
without modification.

### 10.3 The shielding functions

READ, section 3.3 (p. 328), equation (2):

```text
k_i = chi * k0_i * Theta_i * exp(-gamma A_V) ,
```

"with chi a scaling factor for the UV intensity and k0,i the unattenuated rate
in a given radiation field. The shielding function Theta_i accounts for
self-shielding and shielding by H, H2 and the other CO isotopologues".

The reduction of the full line-by-line problem to a two-dimensional table is
justified in section 5.1 (p. 333):

> The transition from atomic to molecular hydrogen occurs much closer to the
> edge of the cloud than the C+-C-CO transition, so the column density of
> atomic H is roughly constant at the depths where shielding of CO is
> important. In addition, H shields CO by only a few per cent. Therefore, it
> is a good approximation to compute the shielding functions on a grid of CO
> and H2 column densities, while taking a constant column of H.

**Table structure** (READ, Tables 5 to 8, pp. 334 and online material pp. 1-4).
Each table is a matrix `Theta[N(12CO), N(H2)]` repeated for six isotopologues,
each block headed by that isotopologue's own `k0`. The grid is:

```text
log10 N(12CO) [cm^-2] :  0, 13, 14, 15, 16, 17, 18, 19     (8 columns)
log10 N(H2)   [cm^-2] :  0, 19, 20, 21, 22, 23             (6 rows)
```

The four parameter sets, from the table captions and footnotes:

| table | b(CO) [km/s] | T_ex(CO) [K] | T_ex(H2) [K] | N(12CO)/N(13CO) | other |
|---|---|---|---|---|---|
| 5 | 0.3 | 5 | 51.5 | 69 | the reference set |
| 6 | 0.3 | 50 | 501.5 | 69 | "Additional rotational lines of CO and H2 were included as described in Sect. 4.4" |
| 7 | 3.0 | 5 | 51.5 | 69 | b(H2) = 11.2, b(H) = 15.9 km/s |
| 8 | 0.3 | 5 | 51.5 | 35 | |

Section 5.1 also states the fixed isotope ratios (Wilson 1999:
[12C]/[13C] = 69, [16O]/[18O] = 557, [18O]/[17O] = 3.6), the extra
`5.2e15 cm^-2` of J'' = 4-7 H2 carried throughout to account for UV pumping,
and that a five times finer grid for more parameter sets is available from the
authors' web page.

**Accuracy of the approximation** (READ, section 5.2, p. 333): tested against
the full integration on a 2880-point grid of translucent cloud models, "The
rate from our approximate method is within 10% of the 'real' rate in 98.3% of
all points (Fig. 6). In no cases is the difference between the approximate
rates and the full model more than 40%."

**And the caveat that matters here**, same section:

> For example, the shielding functions from Table 5 (T ex (CO) = 5 K) can
> easily give photodissociation rates off by a factor of two when applied to a
> high-density, high-temperature PDR.

### 10.4 The T_ex and b dependence, and the hard limit at 500 K

The abstract states the size of the effect: "Increasing the excitation
temperature or the Doppler width can reduce the photodissociation rates and
the isotopic selectivity by as much as a factor of three for temperatures
above 100 K" (READ, p. 323).

The direction is not monotonic. READ, section 4.4 (p. 332): "The 12CO rate
increases from 4 to 16 K, as described in Sect. 4.1. At higher temperatures
the increased overlap with H2 lines takes over and the rate goes down."
Section 4.1 (p. 329) gives the mechanism: raising T_ex spreads the absorption
over more rotational lines, which unsaturates the low-J lines (less
self-shielding, higher rate) but also moves the lines into the H2 absorption
(more shielding, lower rate). Section 4.2 (p. 330) does the same for b: "The
integrated intensity in each line remains the same when b(CO) increases, so a
larger width is accompanied by a lower peak intensity. The resulting reduction
in self-shielding then causes a higher 12CO photodissociation rate", with
"Natural broadening is the dominant broadening mechanism up to
b(CO) = 6e-12 A_tot".

The grid of section 4.4 couples the two, with the auxiliary relations of the
Fig. 3 caption (p. 332): `b(H2) = sqrt(14) b(CO)`, `b(H) = sqrt(28) b(CO)` and
`T_ex(H2) = [T_ex(CO)]^1.5`.

**The model has a stated ceiling of about 500 K, and it is a data limit, not a
choice of convenience.** READ, section 4.4 (p. 332):

> T ex (CO) is raised from 4 to 512 K in steps of factors of two. The v'' = 1
> vibrational level of 12CO lies at 2143 cm^-1 above the v'' = 0 level, so it
> starts to be thermally populated at ~500 K. No data are available on
> dissociative transitions out of this level, so we choose not to go to higher
> excitation temperatures.

**This is the decisive statement of the section for B1.** The regime the CO
ceiling acts in is 300 to 3000 K (and the ceiling itself fires only above about
3000 K, section 2). At those temperatures CO is vibrationally excited, and the
standard photodissociation model is explicitly not extended there because the
molecular data for the dissociating transitions out of `v'' = 1` do not exist.
So the shielding functions cannot simply be looked up for EXHALE's gas: Tables
5 and 6 are computed at `T_ex(CO)` of 5 and 50 K against a base at 1140 K, and
the paper's own extension of the grid stops a factor of two to six below the
temperatures of interest.

MEASURED here for scale, the purely thermal Doppler parameter
`b = sqrt(2 k T / m)` of CO and H2:

| T [K] | b(CO) [km/s] | b(H2) [km/s] |
|---|---|---|
| 300 | 0.42 | 1.57 |
| 1140 | 0.82 | 3.07 |
| 1456 | 0.93 | 3.47 |
| 2000 | 1.09 | 4.06 |
| 3000 | 1.34 | 4.98 |

So the gas of interest sits between Table 5 (b = 0.3) and Table 7 (b = 3.0) in
Doppler width and above every table in excitation temperature. Note also that
Visser's `b(H2) = sqrt(14) b(CO)` is the correct thermal ratio
(`sqrt(m_CO/m_H2) = 3.73`), so the tables' auxiliary widths are consistent with
a common kinetic temperature, and Table 7's `b(H2) = 11.2 km/s` corresponds to
a CO Doppler width of 3.0 km/s, i.e. to gas well above 3000 K if the width is
thermal.

### 10.5 Isotopologues

Not needed here, and this is worth stating rather than passing over. The
paper's own aim is isotope-selective photodissociation, and five of the six
blocks of every table (13CO, C17O, C18O, 13C17O, 13C18O) exist for that
purpose. EXHALE carries one carbon species and one oxygen species with no
isotopic structure (`species_table.f90`), so only the 12CO block of each table
is relevant, and the `N(12CO)/N(13CO)` distinction between Tables 5 and 8 is
immaterial to it: the 12CO blocks of the two differ only through the small
back-shielding of 12CO by 13CO. Nothing in EXHALE would ever call the other
five blocks.

### 10.6 How the shielding function would enter EXHALE's FUV band treatment

The pattern already exists in this code, twice, and Visser's Eq. (2) maps onto
the first of the two exactly.

READ, `lyman_werner.f90`, function `lyman_werner_dissociation_rate`:

```text
k = F_LW / e_lw_photon_erg * h2_lw_dissociation_cross_section(N_H2, T, n_H)
    ,  then  k = k * exp(-tau_cont)
```

with `e_lw_photon_erg` the mean photon energy of a flat-`F_lambda`
912-1110 A band (12.2635 eV, READ from the same file) and `tau_cont` the H2O
and OH continuum of the same interval from `water_photolysis.f90`. The header
of that routine says it is "DB96 eq. (40) in structure, a line self-shielding
factor, exp(-tau) for the continuum, no dust term, with the level-resolved
CLOUDY table in place of a closed-form fit".

Visser's Eq. (2) is the same three factors: an unattenuated rate, a line
shielding function of the star-ward columns, and a continuum term. The
translation is

```text
k_CO(j) = ( F_LW_beam / <h nu> ) * sigma_CO_band * Theta(N_CO, N_H2)
          * exp(-tau_cont) ,
```

where `sigma_CO_band` is the band-mean CO dissociation cross section (the
optically thin limit, which is what `k0/photon flux` amounts to), `Theta` is
the Visser table, and `tau_cont` is the same H2O plus OH continuum the H2 rate
already carries, since the two absorbers share the beam. Three practical
points follow from the existing code rather than from the paper:

1. **The cell mean, not a face value.** `lyman_werner.f90` does not evaluate
   its rate at a face; `lyman_werner_dissociation_rate_cell_mean` integrates
   the local rate across the cell with three-point Gauss-Legendre on segments
   refined to 0.05 dex of column and 0.5 in `tau`, "because the dissociation
   cross section falls by four decades across the self-shielding transition,
   and it falls fastest at the H2 front, exactly where one cell can carry a
   large fraction of a decade of H2 column" (READ, the routine header). The
   Visser `Theta` falls by four decades across its own grid, from 1.0 to
   5.24e-4 along the `N(CO)` axis at `N(H2) = 0` and to 3.9e-7 along the
   `N(H2)` axis, so a CO rate needs the same treatment for the same reason.
2. **The column axis needs a clamp.** `lyman_werner.f90` exposes
   `h2_shield_max_column()` and returns the edge value above it. Visser's grid
   ends at `log N(CO) = 19` and `log N(H2) = 23`, and section 5.1 notes that
   beyond the tabulated range "photodissociation at these depths is typically
   already so slow a process that it is no longer the dominant destruction
   pathway for CO", which is the justification for an edge clamp rather than
   an extrapolation.
3. **`Theta` is a function of two columns, not one.** The H2 rate carries one
   column, `N_H2`; a CO rate carries `N_CO` as well, and `N_CO` has to be
   accumulated star-ward through the same integration that already builds
   `N_H2` and the H2O and OH columns in `FUV_bands.txt`. That is an addition
   to the column bookkeeping, not to the beam.

The `water_photolysis.f90` pattern is the other half: a continuum absorber
with one cross section in each band and `exp(-tau)` on the star-ward column, which
is what the H2O and OH terms of `FUV_bands.txt` are. CO is not that: it is a
line absorber, so it takes the `lyman_werner.f90` shape and not the
`water_photolysis.f90` one. The code says so already, in the note quoted in
section 4.8 above.

### 10.7 What CO photodissociation contributes to `tau_dest`

Two numbers are needed, and both can be MEASURED from files already in this
workspace.

**The optically thin rate.** MEASURED here, from the CO photodissociation
cross section of `VULCAN/thermo/photo_cross/CO/CO_cross.csv` (READ; 0.1 nm
grid, peak `2.599e-16 cm^2` at 107.6 nm) integrated against a flat
`F_lambda` over EXHALE's own 912-1110 A band:

```text
k_thin = 8.066e-7 s^-1 per (erg cm^-2 s^-1) of band flux .
```

The flux-weighted mean dissociation cross section over that band is
`1.577e-17 cm^2`. Using the beam flux of the `oxygen_chemistry` regression
case, 171.5 erg cm^-2 s^-1 (READ, `output/Lyman_Werner.txt` header, which is
the 343.0 erg cm^-2 s^-1 of `input.inp` with the 0.500 dayside dilution
applied):

```text
k_thin = 1.383e-4 s^-1  ,   tau = 7.2e3 s .
```

This is a consistency check rather than an independent result: `k_thin` is the
unshielded limit, which is what Visser's `k0` is, and the 0.1 nm table is
adequate for it because the optically thin rate depends only on the integrated
cross section. The table is useless for anything shielded, which is the whole
point of section 4.8.

**The columns of an actual run.** MEASURED here from the `oxygen_chemistry`
regression case (HD 209458 b with `Oxygen chemistry: True`, the only regression
case that carries CO), integrating `n_CO` from `output/Oxygen_chemistry.txt`
star-ward and calibrating the radius scale against the `N_H2` column that
`output/Lyman_Werner.txt` writes:

```text
at the base cell (T = 1439 K, n_CO = 2.30e10 cm^-3) :
    log10 N(CO) = 18.96 ,   log10 N(H2) = 21.80 .
```

Both sit inside the Visser grid, near its `N(CO) = 19` column and between its
`N(H2)` rows 21 and 22. Log-bilinear interpolation in Table 5 (MEASURED here;
the four bracketing entries are `1.150e-3`, `1.941e-4`, `7.329e-5` and
`1.437e-5`) gives

```text
Theta ~ 2.6e-5  ,   k_CO ~ 3.6e-9 s^-1  ,   tau_dest(photo) ~ 2.8e8 s .
```

**Three conclusions, and the third is the one that matters.**

1. **Photodissociation is by far the fastest CO destruction channel in the
   ceiling's regime, and it is not close.** At 1456 K the barrier channel
   `CO + H -> C + OH` of section 4.1 has `tau_dest` beyond 1e20 s; shielded
   photodissociation gives 2.8e8 s. Any CO destruction model for this code is
   a photodissociation model with the neutral channels as a negligible
   correction, until the gas is hot enough (above about 5000 K, section 7.4)
   for the collisional and ionic channels to take over.
2. **It is still not fast enough for the one-sided model at the base.**
   `tau_res` there is 634 s to 8.0e5 s (READ, section 7.3), so 2.8e8 s misses
   the requirement `tau_dest << tau_res` by two to six decades. Adding
   photodissociation does not rescue alternative A in the layer; it moves the
   failure from twelve decades to three.
3. **Where the ceiling fires, the shielding is gone and the rate is the thin
   one, but the model has no flux there to use.** Higher in the wind the CO
   and H2 columns above a cell fall, `Theta -> 1`, and `tau_dest(photo)`
   approaches 7.2e3 s, which is comparable to `tau_res` and would make an
   instantaneous-destruction ceiling arguable. But the flux is a configuration
   input, not a property of the gas, and in the configuration where the
   molecular gate is defined it is zero: READ,
   `backup/regression/mol_base_handoff/input.inp`, the spectrum is
   `Spectrum type: Power-law` over `[E_low,E_mid,E_high] = [13.60, 123.98,
   1.24e3]` eV, that is wavelengths at or below 911.8 A, and the file sets no
   `Stellar LW flux` key. The 912-1118 A interval that carries every CO band
   is outside the modeled spectrum entirely in that case. So in
   `mol_base_handoff` the CO photodissociation rate is identically zero and
   `tau_dest` is the barrier channel's 1e12 s or worse, while the ceiling
   removes CO instantaneously.

The verdict of section 7 is therefore unchanged and is now quantitative on the
one channel that could have overturned it: with photodissociation included at
its published unshielded rate and its published shielding function, the
ordering `tau_dest << tau_res` still fails in the molecular layer by two to six
decades, and in the configuration the molecular regression gate runs it fails
by twelve because that configuration carries no photodissociating flux at all.

**What would be needed to do better than the estimate above**, in order:
(i) `Theta` recomputed at the actual `T_ex(CO)` and `b(CO)` of the gas, which
the paper's own grid does not reach and which its section 4.4 says cannot be
extended past about 500 K without new molecular data for the `v'' = 1` level;
(ii) a CO cross section resolved on the band lines rather than the 0.1 nm
continuum grid of the VULCAN table, which is the same data the shielding
function is built from; (iii) the star-ward `N(CO)` column added to the
existing FUV column integration.

---

## 11. UMIST RATE12 and the Venot networks

The three publisher PDFs were supplied by the user on 2026-09-06:
`references/McElroy_2013A&A_550_A36.pdf`,
`references/Venot_2012A&A_546_A43.pdf` and
`references/Venot_2020A&A_634_A78.pdf`.

### 11.1 UMIST RATE12 (McElroy et al. 2013, A&A 550, A36)

**The paper does not print any reaction entries.** It is a description of the
database and its format: sections 2 and 2.1 define the file layout and the rate
formulae, section 3 lists what changed since RATE06, section 5 shows two
model applications. MEASURED by grep over the full text: there is no line in
the paper containing a CO reaction with its coefficients, and the string "CO"
appears only in discussion of self-shielding, of CSE models and of species
names. So the specific entries the brief asks for, `He+ + CO`, `H+ + CO`,
`CO + H` and the CO photodissociation, **cannot be quoted from this paper**;
they are in the network file. **That file is now available**: the user pointed
to `~/RT_Codes/UMIST/`, which holds the later RATE22 release, and every CO
entry of it is quoted in section 12.3. What is said here about RATE12 stands
as the reading of the RATE12 paper; the numbers come from section 12.

What the paper does give, and it is directly useful, is the **format and the
quality flags each entry carries** (READ, Table 1, p. 2): reaction number,
type code, up to two reactants and four products, `NE` (the number of fitted
temperature ranges), the coefficients `alpha, beta, gamma`, the fitted bounds
`Tl` and `Tu`, a source type `ST` and an accuracy code `ACC`. The two flag
sets are quoted verbatim from that table:

```text
ST  source type :  E  Estimated
                   M  Measured
                   C  Calculated
                   L  A combination of a number of experimental values
                      from the Literature
ACC accuracy    :  A  <25%
                   B  <50%
                   C  within a factor of 2
                   D  within an order of magnitude
                   E  highly uncertain
```

The rate forms (READ, section 2.1, p. 2). Two-body reactions:

```text
k = alpha (T/300)^beta exp(-gamma/T)   cm^3 s^-1        (their eq. 1)
```

and interstellar photoreactions, type `PH`:

```text
k = alpha exp(-gamma A_V)   s^-1                        (their eq. 4)
```

"where alpha represents the rate coefficient in the unshielded interstellar
ultraviolet radiation field, A_V is the dust extinction at visible wavelengths
and gamma is the parameter used to take into account the increased dust
extinction at ultraviolet wavelengths."

**So RATE12's CO photodissociation entry is an unshielded rate times a dust
term, exactly as KIDA's is, and the database does not carry the CO shielding
function.** The paper says so itself, in the introduction (p. 1):

> Other environments for which Rate12 can be used without modification, along
> with careful treatment of molecular hydrogen and CO self shielding, are
> photodissociation regions (PDRs), hydrodynamic shock regions and diffuse
> clouds.

that is, the self-shielding is the user's responsibility and is what section 10
of this document supplies.

**Extrapolation policy** (READ, section 2, p. 2), which is the same question
KIDA answers in section 5.2 above:

> In order to evaluate the rate coefficient outside the given temperature
> range, we recommend that the user chooses the expression that is closest to
> the temperature of interest. While there is no guarantee that this will give
> the correct rate coefficient, we have taken care to ensure that, when
> evaluated at low temperatures (<50 K), no rate coefficient will become
> unphysically large.

The care is taken at the cold end only; nothing is claimed at 3000 K. Combined
with the reaction-type census of their Table 2 (p. 2), which counts 619
neutral-neutral and 2589 ion-neutral entries out of 5399 total, the shape of
the database for this question is clear: RATE12 is where the `He+ + CO` and
`H+ + CO` channels of section 4.6 live, it stores them in a form whose
temperature validity has to be read off each entry's own `Tl`/`Tu`, and it
does not carry a neutral CO destruction channel useful at 3000 K, since
`CO + H -> C + OH` is a combustion reaction and not an interstellar one.

### 11.2 Venot et al. 2012 (A&A 546, A43)

**The validity range is stated in the abstract** (READ, p. 1):

> The network we release is robust for temperatures within 300-2500 K and
> pressures from 10 mbar up to a few hundred bars, for species made of C, H,
> O, and N. It is validated for species up to 2 carbon atoms and for the main
> nitrogen species (NH3, HCN, N2, NOx).

and again in section 2.1.1 (p. 3), for the combustion mechanism the scheme is
built on: "This C0-C2 mechanism has been widely validated in the 300-2500 K,
0.01-100 bar range for several types of reactors, such as shock-tubes,
perfectly stirred reactors, plug-flow reactors, rapid compression machines,
and laminar flames".

**So Venot 2012 does not reach 3000 K.** It stops at the same 2500 K as the
NIST-derived VULCAN network of section 5.5, and for the same reason: both are
combustion mechanisms, and the combustion experiments that validate them stop
there.

**The pressure limit is the more serious mismatch for this code.** The lower
bound, 10 mbar in the abstract and 0.01 bar in section 2.1.1, is `1e-2` bar.
EXHALE's molecular base sits at `p_base = 1.000e-06` bar (READ,
`backup/regression/mol_base_handoff/base.inp`), which is **four orders of
magnitude below the low-pressure end of the stated validity range**. A network
validated against shock tubes and stirred reactors at 0.01 to 100 bar is not
validated at a microbar, where three-body channels are negligible and the
radiative and photochemical ones are not. That limit applies to Venot 2012 and
2020 and, through the shared combustion provenance, to the VULCAN network as
well.

**The scheme's construction**, for the record (READ, section 2.1.1, p. 3): 957
reversible and 6 irreversible reactions, 105 neutral species; rate
coefficients are "those recommended for the individual processes by the main
kinetics databases for combustion (Tsang & Hampson 1986; Manion et al. 2008;
Smith et al. 1999; Baulch et al. 2005)"; the reaction list is released through
KIDA. Reverse rates are computed thermodynamically, with three named
exceptions (READ, section 2.1.3, p. 4): "our nominal network uses
thermodynamical reversal for most of the reactions but not for three important
ones", the unimolecular initiations of methane. This is the same construction
VULCAN uses (section 5.5) and the same one `rate_from_detailed_balance`
implements here.

**CO photodissociation is carried, without shielding.** READ, their
photodissociation table (the `J` list): entry J9 is `CO + h nu -> C + O(3P)`,
with the cross section from Olney et al. (1997) and the branching ratio from
Huebner et al. (1992). MEASURED by grep over the full text: the string
"shield" does not occur anywhere in the paper. So Venot 2012 treats CO
photodissociation as a bare cross section with no self-shielding function,
which for a hot Jupiter lower atmosphere at `N(CO)` of `1e19 cm^-2`
(section 10.7) overestimates the rate by the four decades of `Theta`.

MEASURED by grep, no bimolecular CO destruction reaction is printed in the
paper's own text; the reaction list is in the released file, not in the
article.

### 11.3 Venot et al. 2020 (A&A 634, A78)

This is an update of the methanol sub-mechanism of V12, not a new scheme.
READ, abstract (p. 1): "The updated scheme involves 108 species linked by a
total of 1906 reactions", and section 3 (p. 2): "The full and updated chemical
scheme that we present in this paper, hereafter called V20, contains 108
species, 948 reversible reactions, and 10 irreversible reactions (i.e., 1906
reactions in total). This scheme can be downloaded from the KInetic Database
for Astrochemistry (KIDA)". A reduced scheme of 44 species and 582 reactions
is derived in Appendix D for 3D models.

**No new validity range is stated, and none is claimed.** MEASURED by grep:
the strings "2500", "300 K" and "3000 K" do not occur in a validity statement
anywhere in the paper. The experimental base cited for the update is narrower
than V12's, not wider: the autoignition studies of section 2 (p. 2) are at
"temperatures greater than 1300 K and moderate pressures (5 bar)" and, for the
rapid compression machine work, "temperatures ranging from 800 to 1700 K and
pressures between 1 to 50 bar". So V20 inherits V12's 300-2500 K,
0.01-100 bar envelope and does not extend it.

What V20 changes that touches CO at all is the methanol pathway, and the paper
reports its effect on CO only for the ice giants (READ, abstract, p. 1):
"Concerning Uranus and Neptune, the update of the chemical scheme modifies the
abundance of CO and thus impacts the deep oxygen abundance required to
reproduce the observational data." For hot Jupiters the conclusion (p. 12) is
that "The new updated scheme V20 gives very similar results to the former
scheme for hot Jupiters."

### 11.4 What the three compilations settle

Collecting sections 5.2, 5.5, 11.1, 11.2 and 11.3, the brief's question
"which of them give CO destruction rates valid to 3000 K" has one answer:

| compilation | stated validity | reaches 3000 K? | carries CO shielding? |
|---|---|---|---|
| KIDA (Wakelam 2012) | 10-300 K by default | no | no (Draine field, `A exp(-C A_V)`) |
| UMIST RATE12 (McElroy 2013) | per entry, `Tl`/`Tu`; interstellar provenance | no, and not claimed | no, stated to be the user's responsibility |
| UMIST RATE22 (Millar 2024) | per entry; `Tu = 41000` is an undocumented default on 6562 of 8767 entries | only the Langevin ion entries, defensibly | no, and it names Visser 2009 as where to get it |
| Venot 2012 | 300-2500 K, 0.01-100 bar | no | no (bare cross section, J9) |
| Venot 2020 | inherits V12 | no | not addressed |
| VULCAN (Tsai 2017) | 500-2500 K | no | no (no "shield" in the source) |
| NSRDS-NBS 67 (Westley 1980) | per entry; `CO + M` at 7000-15000 K | only above 7000 K | not applicable |
| Visser 2009 | `T_ex(CO)` to 512 K, hard data limit near 500 K | no | yes, and it is the only source that does |

**No source reaches 3000 K for a channel whose rate depends on temperature.**
Two stop at 2500 K because combustion experiments stop there, two are
interstellar databases whose defaults are 300 K or an undocumented 41000 K,
one covers only above 7000 K, and the photodissociation model stops near
500 K because the molecular data for vibrationally excited CO do not exist.
The temperature band in which the CO ceiling fires, roughly 3000 to 4000 K, is
a gap between the combustion literature below it and the shock-tube literature
above it.

**The exception is the one that matters, and it is why section 12.6 revises
the verdict.** RATE22's `He+ + CO -> O + C+ + He` has `beta = gamma = 0`: it
is temperature independent at the Langevin value, so there is no temperature
range to cover. KIDA's own rule names that case as the one extrapolation it
permits (section 5.2), and UMIST's accuracy code for the entry is `A`, better
than 25%. The channel that removes CO where the ceiling fires therefore does
not need the missing 3000 K data at all: it needs `n(He+)`, which the code
already carries.


---

## 12. UMIST RATE22: the network file, and what it changes

The user pointed to `~/RT_Codes/UMIST/`, which holds the RATE22 release of the
UMIST Database for Astrochemistry: `rate22_final.rates` (8767 lines, the full
ratefile), `rate22_G_final.rates`, `rate22_dipole.specs`,
`rate22_revised_CtoO_0.44.specs`, `UDfA2024.pdf` and `where_UMIST.txt`
(`https://umistdatabase.uk/`). **This supersedes the RATE12 data-file request
of section 8 item 2**, and it changes the assessment; section 12.5 states how.

**The paper is the published version.** READ from the PDF's own header:
"A&A, 682, A109 (2024)", `https://doi.org/10.1051/0004-6361/202346908`,
"(c) The Authors 2024", "Received 15 May 2023 / Accepted 27 October 2023",
title "The UMIST Database for Astrochemistry 2022", authors T. J. Millar,
C. Walsh, M. Van de Sande and A. J. Markwick. ADS gives bibcode
**2024A&A...682A.109M**, doi 10.1051/0004-6361/202346908, properties
`EPRINT_OPENACCESS, OPENACCESS, PUB_OPENACCESS` (READ, ADS API). Note the
mismatch of years that the citation has to carry: the database is RATE22, the
paper is 2024.

### 12.1 File format and rate forms

READ, section 2.1 (p. 2):

> Our basic gas-phase ratefile, RATE22, now contains some 8767 individual rate
> coefficients. These correspond to 737 species involving 17 elements,
> increases of over 40% and 55% in reactions and species, respectively, from
> RATE12. The additional elements are Al, Ar, Ca, and Ti. The basic format is
> that each line of data consists of 18 colon-separated entries: the first two
> are the reaction number and reaction type, defined in Table 1, followed by
> two reactants and up to four products. The ninth entry denotes the number of
> temperature ranges, NTR, over which the rate coefficient is defined, while
> entries 10-12 give the values of alpha, beta and gamma used to calculate the
> rate coefficients. Entries 13-14 give the temperature range over which the
> rate coefficient is defined, entry 15 provides the method by which the rate
> coefficient has been determined (M: measured; C: calculated; E: estimated;
> L: literature).

Entry 16 is the accuracy code, entry 17 the DOI or URL, entry 18 the notes;
"Further details on the file format are given in McElroy et al. (2013)", so
the accuracy codes are the A to E scale already quoted in section 11.1
(A < 25%, B < 50%, C within a factor of 2, D within an order of magnitude,
E highly uncertain).

Rate forms, READ from section 2.2 (pp. 2-3): binary reactions follow the
de Kooij-Arrhenius formula `k = alpha (T/300)^beta exp(-gamma/T)` cm^3 s^-1
(their eq. 1); cosmic-ray ionization `k = alpha` s^-1 (eq. 2); interstellar
photoprocesses `k = alpha exp(-gamma A_V)` s^-1 (eq. 3), "where alpha
represents the rate coefficient in the unshielded interstellar ultraviolet
radiation field"; cosmic-ray induced photoreactions
`k = alpha (T/300)^beta gamma/(1 - omega)` s^-1 (eq. 4).

Provenance census, READ from their Table 1 (p. 2): 8767 reactions and 737
species against RATE12's 6173 and 467; 6990 DOIs, 474 URLs, 1303 with no
reference; 1889 measured, 1508 calculated, 2996 literature and 2374 estimated.
Section 2.2 adds that ten reactions are fitted as sums of formulae, "so that
the number of fully independent reactions in RATE22 is 8757".

### 12.2 What the paper says about self-shielding, and where it sends the reader

READ, section 2.2 (p. 3), quoted in full because it settles the question
section 11.1 could only infer:

> We have not explicitly included the effects of self-shielding in the
> database. Such a process can occur in situations in which dissociating
> photons are absorbed through line rather than continuum processes and acts
> in addition to the extinction caused by dust grains. It can be important for
> molecules including H2, CO, N2, OH, and H2O and depends on the column
> density of the molecule with these species showing almost complete shielding
> once a column density of around 10^15 cm^-2 is reached (Heays et al. 2017).
> Numerical approaches often involve the use of look-up tables. Details on the
> self-shielding of H2, CO, and N2 are discussed by Sternberg et al. (2014),
> Visser et al. (2009), Li et al. (2013) and Heays et al. (2014).

So the database itself directs the user to Visser et al. (2009) for the CO
shielding function, which is section 10 of this document. The two sources fit
together exactly: RATE22 supplies `alpha` (the unshielded rate and the ion
channels), Visser supplies `Theta`.

### 12.3 Every CO reactant entry in `rate22_final.rates`

MEASURED here by parsing the file: **68 entries have CO as one of the two
reactants**, distributed as IN 23, RA 21, NN 13, CE 5, AD 3, CP 1, CR 1,
PH 1. The 21 radiative associations are all with hydrocarbon cations
(`CnHm+ + CO -> ... + PHOTON`, all with `T = 10-300 K`) and none of their
partners or products exists in EXHALE; they are not listed below.

The entries the brief asks for, quoted verbatim from the file (the columns are
`alpha`, `beta`, `gamma`, `Tl`, `Tu`, method, accuracy):

| # | type | reaction | alpha | beta | gamma | T range [K] | ST | ACC | source | covers 300-3000 K? | products in EXHALE? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4068 | IN | `He+ + CO -> O + C+ + He` | 1.60e-09 | 0.00 | 0.0 | 10-41000 | M | A | `06_0948D.notes` | **yes** | **yes** (O, C+, He) |
| 7351 | NN | `H + CO -> OH + C` | 1.10e-10 | 0.50 | 77700.0 | **2590**-41000 | L | C | Mitchell GF, ApJSS, 54, 81 (1984) | no, starts at 2590 | **yes** (OH, C) |
| 8259 | PH | `CO + PHOTON -> O + C` | 2.40e-10 | 0.00 | 3.88 | 10-41000 | C | C | Heays et al., A&A 602, A105 (2017) | yes | **yes** (O, C) |
| 7665 | NN | `OH + CO -> CO2 + H` | 2.81e-13 | 0.00 | 176.0 | 80-**3150** | M | A | Frost MJ et al., JPC, 97, 12254 (1993) | **yes** | no (CO2) |
| 3497 | IN | `H3+ + CO -> HCO+ + H2` | 1.36e-09 | -0.14 | -3.4 | 10-**400** | C | A | `06_0771D.notes` | no | no (HCO+) |
| 3498 | IN | `H3+ + CO -> HOC+ + H2` | 8.49e-10 | 0.07 | 5.2 | 10-**400** | C | A | `06_0771D.notes` | no | no (HOC+) |
| 7140 | NN | `C + CO -> C2 + O` | 2.94e-11 | 0.50 | 58025.0 | 1934-41000 | L | C | Mitchell GF, ApJSS, 54, 81 (1984) | partly | no (C2) |
| 7318 | NN | `CO + O2 -> CO2 + O` | 4.20e-12 | 0.00 | 24000.0 | 300-6000 | L | A | NIST, Tsang & Hampson, JPCRD 15, 1087 | yes | no (CO2, O2) |
| 7319 | NN | `CO + O2H -> CO2 + OH` | 2.51e-10 | 0.00 | 11900.0 | 250-2500 | L | A | NIST, Tsang & Hampson, JPCRD 15, 1087 | no (stops 2500) | no (CO2, O2H) |
| 3195 | IN | `H2+ + CO -> HCO+ + H` | 2.16e-09 | 0.00 | 0.0 | 10-41000 | M | A | | yes | no (HCO+) |
| 456 | CE | `H2+ + CO -> CO+ + H2` | 6.44e-10 | 0.00 | 0.0 | 10-41000 | M | A | | yes | no (CO+) |
| 3322 | IN | `H2O+ + CO -> HCO+ + OH` | 5.00e-10 | 0.00 | 0.0 | 10-41000 | M | A | Jones et al., CPL, 77, 484 (1981) | yes | no (HCO+) |
| 4672 | IN | `OH+ + CO -> HCO+ + O` | 1.05e-09 | 0.00 | 0.0 | 10-41000 | M | A | Jones et al., CPL, 77, 484 (1981) | yes | no (HCO+) |
| 664 | CE | `O+ + CO -> CO+ + O` | 4.90e-12 | 0.50 | 4580.0 | 2000-10000 | C | B | | partly | no (CO+) |
| 822 | CP | `CO + CRP -> CO+ + e-` | 3.90e-17 | 0.00 | 0.0 | 10-41000 | L | C | Prasad and Huntress 1980 | not applicable | no |
| 990 | CR | `CO + CRPHOT -> O + C` | 1.30e-17 | 1.17 | 105.0 | 10-41000 | L | C | | not applicable | yes |

**Four findings from the file, each of which answers a question this document
had to leave open or infer.**

1. **`H+ + CO` does not exist in RATE22.** MEASURED by direct search: the only
   entry pairing H+ with CO is the *reverse*, number 499,
   `CE: H + CO+ -> CO + H+`, `alpha = 7.50e-10`, `beta = 0`, `gamma = 0`,
   `T = 10-41000 K`, `M`, `B`, doi 10.1103/PhysRevLett.52.2084. The database
   carries the exothermic direction only. This is published confirmation of
   what section 4.6 could state only from the ionization energies: the charge
   transfer `H+ + CO -> CO+ + H` is endothermic and is not a CO destruction
   channel. There is also no `H+ + CO -> HCO+ + PHOTON` radiative association
   in the file, so the guess offered in section 4.6 that the networks carry
   one instead is **wrong and is corrected here**: they carry nothing.

2. **`CO + M -> C + O + M` does not exist in RATE22 either.** MEASURED: the
   collisional-dissociation block (`CD`) has 14 entries, unchanged from
   RATE12 (their Table 1), and it covers `H2`, `CH`, `OH`, `H2O`, `O2`,
   `HOC+` and `O2-` as the dissociating partner but **not CO**. The strongest
   bond in the network is simply not broken collisionally anywhere in the
   database. This is consistent with NSRDS-NBS 67 putting the only evaluated
   `CO + M` rate at 7000-15000 K (section 4.5).

3. **`CO + H2` does not exist as a CO destruction channel.** MEASURED: every
   file entry pairing CO with H2 has CO as a *product*.

4. **`He+ + CO -> O + C+ + He` is the entry that matters, and it is the best
   one in the table**: measured (`M`), accuracy `A` (better than 25%),
   temperature independent (`beta = gamma = 0`) at the Langevin value
   `1.6e-9 cm^3 s^-1`, and **every product is a species EXHALE already
   carries**. Its temperature independence is exactly the case KIDA's rule
   (section 5.2) allows to be extrapolated: "Extrapolations outside this range
   of temperature are not recommended unless the rate coefficient is predicted
   to be totally independent of temperature, as occurs for the Langevin model
   of exothermic ion-non-polar-neutral reactions."

**Two caveats on the temperature ranges, both MEASURED from the file.**

- The upper bound `41000` appears on **6562 of the 8767 entries**, and the
  string does not occur anywhere in the paper. It is an undocumented file
  default for entries with no measured upper limit, **not a statement that the
  rate was validated to 41000 K**. Every rate quoted above with `Tu = 41000`
  must be read as "no upper bound was determined", and for the Langevin
  `He+ + CO` entry that is a defensible extrapolation while for others it is
  not.
- The `H + CO -> OH + C` entry has a genuine *lower* bound, `Tl = 2590 K`, and
  it disagrees with NSRDS-NBS 67 (section 4.1). In cm^3 s^-1 the RATE22 form
  is `1.10e-10 (T/300)^0.5 exp(-77700/T) = 6.35e-12 T^0.5 exp(-77700/T)`
  against Westley's `3.32e-11 T^0.5 exp(-77755/T)`: **RATE22 is a factor 5.2
  lower** (MEASURED). The exponentials agree to 0.07 per cent, so the
  difference is entirely in the pre-exponential, and the two trace to
  different chains (Mitchell 1984 against Benson & Golden 1975). The detailed
  balance check of section 4.1 favours Westley's value, since it reproduces it
  to 13-19 per cent from the code's own thermochemistry. Either way the
  channel is far too slow to matter, so the discrepancy is recorded and not
  resolved.

### 12.4 The `_G` ratefile

MEASURED by comparing the two files line by line: `rate22_G_final.rates`
differs from `rate22_final.rates` on **55 lines, 51 of which have `alpha` set
to zero** in the `_G` version (the remaining four differ only in reference
formatting). This is the "reduced ratefile" of the paper's section 2.3
(p. 3), which removes the endothermic reactions that Tinacci et al. (2023)
identified in KIDA: "We have searched for each of these in our database
finding that 53 overlap in terms of reactants and products with our list and
should be removed... Both the full and the reduced ratefiles are made
available to the community" (READ). One CO entry is among them, number 309
`CE: CN+ + CO -> CO+ + CN`, whose `alpha` goes from `6.30e-10` to `0.00e-10`.
None of the entries this document uses is zeroed in `_G`.

### 12.5 The ion channel changes the timescale assessment

Section 7 concluded that `tau_dest << tau_res` fails everywhere, and it
reached that conclusion from the neutral channels alone, because no source for
the ion channels was available. **RATE22 supplies the ion channel, and it
changes the answer.**

MEASURED on the `oxygen_chemistry` regression case, the only case that carries
CO (see section 10.7 for how the columns are built). `tau_photo` uses the
optically thin rate of section 10.7 with the Visser Table 5 shielding function
bilinearly interpolated in `log N(CO)` and `log N(H2)`; `tau(He+)` uses RATE22
entry 4068 with the `HeII` column of `output/Ion_species.txt`; `tau_res` is
`min(tau_adv, tau_diff)` as the run itself writes them in
`output/Oxygen_chemistry.txt`:

| r [R_p] | T [K] | Theta | tau_photo [s] | tau(He+) [s] | tau_dest [s] | tau_res [s] | tau_dest/tau_res |
|---|---|---|---|---|---|---|---|
| 1.0000 | 1439 | 2.6e-05 | 2.8e+08 | 1.3e+10 | 2.7e+08 | 4.0e+05 | 6.7e+02 |
| 1.0078 | 1456 | 7.5e-05 | 9.6e+07 | 2.9e+09 | 9.3e+07 | 1.3e+07 | 7.2e+00 |
| 1.0167 | 1455 | 2.5e-04 | 2.9e+07 | 1.3e+08 | 2.3e+07 | 1.3e+07 | 1.8e+00 |
| 1.0307 | 1452 | 1.3e-03 | 5.8e+06 | 2.7e+07 | 4.8e+06 | 3.2e+06 | 1.5e+00 |
| 1.0530 | 1446 | 4.9e-03 | 1.5e+06 | 3.3e+05 | 2.7e+05 | 1.5e+06 | 1.8e-01 |
| 1.0883 | 1420 | 2.4e-02 | 3.0e+05 | 4.7e+02 | 4.6e+02 | 8.1e+05 | 5.7e-04 |
| 1.1443 | 1279 | 1.8e-01 | 3.9e+04 | 4.7e+01 | 4.7e+01 | 4.9e+05 | 9.6e-05 |
| 1.2334 | 847 | 5.5e-01 | 1.3e+04 | 4.5e+02 | 4.3e+02 | 3.2e+05 | 1.4e-03 |
| 1.3749 | 794 | 6.3e-01 | 1.2e+04 | 1.3e+04 | 6.1e+03 | 2.2e+05 | 2.8e-02 |
| 1.9567 | 794 | 6.6e-01 | 1.1e+04 | 1.3e+04 | 5.9e+03 | 1.2e+05 | 5.0e-02 |
| 3.4247 | 794 | 6.9e-01 | 1.0e+04 | 1.3e+04 | 5.8e+03 | 8.2e+04 | 7.1e-02 |

The supporting composition, MEASURED from the same run: `n(He+)` is
`4.9e-2 cm^-3` in the base cell, rises through the helium ionization front to
a maximum of `1.3e7 cm^-3` near `r = 1.1` and stays near `5e4 cm^-3` in the
outer wind; `n(H+)` goes from `1.6e3` at the base to `7e5 cm^-3` outside.

**What the table says.**

1. **The ordering `tau_dest << tau_res` holds above `r ~ 1.05`, by one to four
   decades**, and it is the He+ channel that carries it: `tau(He+)` falls from
   `1.3e10 s` at the base to `47 s` at `r = 1.14`. Photodissociation alone
   would not do it inside `r ~ 1.2`, because the shielding is still strong
   there.
2. **It fails below `r ~ 1.04`**, by up to 670 at the base, and it fails
   because helium is neutral there and the CO and H2 columns above are thick.
   That is the molecular layer, and it is exactly where CO sits at its
   equilibrium value and the ceiling never fires.
3. **The boundary of the destruction domain is compositional, not thermal.**
   It is the helium ionization front, set by `n(He+)`, and it has nothing to
   do with the CO equilibrium constant that the present ceiling switches on.
   In a converged hot-Jupiter wind the two surfaces are nested the right way
   round (helium ionizes well below the radius where `T` reaches 3000 K), so
   the ceiling fires inside a region where destruction really is fast, and its
   answer is approximately right. But its stated reason, "where the
   equilibrium constant says CO cannot exist, it does not", is not the
   physics that removes the CO. The physics is He+ charge transfer and, higher
   up, unshielded photodissociation.

**Caveat on the measurement, and it is not a small one.** The
`oxygen_chemistry` case is a step-pinned relaxation snapshot, not a converged
solution, and its temperature never exceeds 1456 K over the physical cells
(MEASURED: 793 to 1456 K). So the table above does **not** contain the
ceiling's own firing region, and the statement in point 3 about how the two
surfaces nest is an inference from the structure, not a measurement of a
converged wind. Confirming it needs the same profile on a converged run that
reaches 3000 K, which this item did not run and which is stated here as what
would settle it.

### 12.6 The revised verdict

Sections 1 and 7 said the one-sided destruction model of B1 T5.1 is not
constructible. **With the RATE22 He+ channel that is no longer right, and this
section supersedes that part of the verdict.** The corrected statement:

- **The model is constructible above the helium ionization front.** The
  destruction partners are named and published: `M_1 = He+` with
  `k = 1.6e-9 cm^3 s^-1` (RATE22 entry 4068, measured, accuracy better than
  25%, products O + C+ + He, all carried by EXHALE), and
  `J_CO = k_thin Theta(N_CO, N_H2)` with `k_thin` from the unshielded rate of
  Visser et al. and `Theta` from their Table 5 or 6, products C + O. Both
  inequalities of T5.2 then hold there: `tau_dest << tau_res` by one to four
  decades (MEASURED, section 12.5), and `tau_dest < tau_form` by construction
  wherever the ceiling fires (section 7.1).
- **It is not constructible in the molecular layer**, `r < 1.04` in the run
  measured, where `tau_dest` is `1e8` s or longer against a residence time of
  at most `1e5` s. The domain argument T5.2 therefore returns a domain with a
  lower boundary, and that boundary must be reported per cell rather than
  assumed.
- **What still has to be added to make T5.1 real** is now a short list, not an
  open question: the `N_CO` column in the FUV column integration
  (section 10.6), a `Theta` table at the run's own `T_ex` and `b` rather than
  at 5 or 50 K (section 10.4, and the paper's own limit near 500 K), and the
  He+ density, which the code already has.
- **What does not change**: the present ceiling is still a thermodynamic
  switch standing in for a kinetic process, its switching surface is not the
  destruction domain's boundary, and B1 decision 7 (exclude runs in which it
  fires during accepted physical integration) remains the correct treatment
  until T5.1 is implemented. The reason for the exclusion is now different:
  not that no destruction model exists, but that the one the code uses is not
  the one the physics supports.

---

## 13. Moses et al. (2011) and Tsai et al. (2021)

Both publisher PDFs were supplied by the user on 2026-09-06,
`references/Moses_2011_ApJ_737_15.pdf` (ApJ 737, 15, 37 pp.) and
`references/Tsai_2021_ApJ_923_264.pdf` (ApJ 923, 264, 42 pp.). Section 5.4
above quoted nothing from the arXiv copy and nothing has to be withdrawn; what
follows replaces it.

### 13.1 Moses et al. 2011: no stated validity range, and no self-shielding

**There is no single stated temperature and pressure validity for the
network.** MEASURED by reading section 2.2 and searching the full text: the
paper states the network's provenance rather than its envelope. READ,
section 2.2 (p. 3):

> extensive updates that account for high-temperature kinetics have been
> included based largely on combustion-chemistry literature (e.g., Baulch
> et al. 1992, 1994, 2005; Atkinson et al. 1997, 2006; Smith et al. 2000;
> Tsang 1987, 1991; Dean & Bozzelli 2000)... The model contains ~1600
> reactions, with the rate coefficients of ~800 of the reactions being taken
> from literature values, and the remaining ~800 reactions being the reverse
> of these "forward" reactions, with the rate coefficients of the reverse
> reactions being calculated internally at each pressure-temperature point
> along the grid using the thermodynamic principle of microscopic
> reversibility (e.g., k_for/k_rev = K_eq...).

That is the same construction as VULCAN (section 5.5) and Venot et al.
(section 11.2), and the same one `rate_from_detailed_balance` implements here.
The only explicit temperature interval quoted for a fit in the paper is 500 to
2500 K, for the three methanol channels R864, R862 and R860 (READ, section 2.5,
p. 10: "The predicted rate coefficients for reactions R864-R860 are well
reproduced over the 500-2500 K region by the modified Arrhenius expressions").
So **the answer to the brief's question is the same as for every other
combustion-derived network in this document: it does not reach 3000 K, and it
does not claim to.**

The pressure range is stated, and it is the model's grid rather than the
network's validity (READ, section 2.1, p. 3): profiles come "from the 3D GCM
simulations of Showman et al. (2009) for the region from 170 bar to a few
microbar or the 1D models of Fortney et al. (2006, 2010) for the region from
1000 to 1.0e-6 bar". The three chemical regimes are given in section 3
(p. 10): "Thermochemical equilibrium dominates at pressures greater than a few
bars, transport-induced quenching can dominate for some species in the ~1 to
1e-3 bar region, and photochemistry can dominate at pressures less than
~1e-3 bar". EXHALE's molecular base at 1e-6 bar is the very top of that grid,
and the paper says plainly what it is not doing there (READ, section 2.1,
p. 3): "Although the models extend into the thermosphere, our intent is not to
specifically model thermospheric chemistry. Instead, the thermospheric
temperature profile and escape boundary condition are included in an ad hoc
manner".

**CO destruction.** The paper's own account is that CO photolysis is the
initiating destruction and that the network largely recycles it. READ,
section 3.2 (p. 14):

> Even in the case of CO photolysis, some recycling pathways such as the
> following can operate
>
>     CO + h nu -> C + O
>     C + H2 -> CH + H
>     CH + H2O -> H2CO + H
>     H2CO + H -> HCO + H2
>     HCO + H -> CO + H2
>     Net : H2O -> O + H2 ,
>
> allowing CO to persist. However, these recycling schemes are not 100%
> effective, and some of the CO photodestruction leads to the production of
> other species. Atomic carbon and oxygen released from CO photolysis at high
> altitudes can remain as C and O

and the thermal channel that appears in the text is the CO2 one already
covered in section 4.3, written there as the reversible pair
`OH + CO <-> CO2 + H` (READ, section 3.2, p. 14). Every intermediate of the
recycling scheme (CH, H2CO, HCO) is a species EXHALE does not carry, which is
the same finding as section 5.7 item 3 for the VULCAN network: remove the
hydrocarbons and there is no CO chemistry left, only the photolysis.

**No self-shielding.** MEASURED by searching the full text: the string
"shield" occurs exactly once in the paper, at p. 20, and it refers to NNH
being "shielded from photolysis by a larger surrounding column", not to CO.
Moses et al. therefore compute CO photolysis with a bare cross section, as
Venot et al. (2012) do (section 11.2). This matters for the comparison Tsai
et al. make: they note (READ, Tsai et al. 2021, section 4.1.1) that "CO
photolysis appears to be stronger in M11 and generates more atomic carbon
around microbar" than in VULCAN, which is what an unshielded CO rate does.

The paper's own uncertainty statement on the carbon-oxygen kinetics, for the
record (READ, section 2.5, p. 10): "the kinetics of CO <-> CH4
interconversion in reducing environments remains somewhat uncertain. The rate
coefficients for the above rate-limiting reactions are uncertain by perhaps a
factor of three".

### 13.2 Tsai et al. 2021: the answer to the open question is no

Section 5.6 left open whether the photochemical version of VULCAN implements
CO self-shielding, and called it "the one check that could show a working
implementation of section 10.6 in a code of this kind". **It does not, and its
numerical scheme could not carry one.**

**The radiative transfer is continuum Beer-Lambert plus two-stream diffuse.**
READ, section 2.3 (p. 6), their equations (8) and (9):

```text
J(z, lambda) = J(inf, lambda) exp(-tau(z, lambda)/mu) + J_diff(z, lambda)   (8)

tau = int [ sum_i (sigma_a,i + sigma_s,i) n_i ] dz                          (9)
```

"where sigma_a,i and sigma_s,i are the cross section of absorption and
scattering, respectively. The absorption cross section sigma_a,i can differ
from the photodissociation cross section because absorption is not necessarily
followed by dissociation." The photolysis rate is then their equation (11),
`k = int_lambda q(lambda) sigma_a(lambda) J(z, lambda) dlambda`.

That is a **sum over absorbers of `sigma n`**, with no term of the form
`Theta(N_CO, N_H2)`. Self-shielding is present only to the extent that the
wavelength grid resolves the absorbing lines, which brings in the second
point.

**The wavelength grid is 0.1 nm, and that is the whole answer.** READ,
Appendix B (p. 39):

> ...survey with a constant 0.1 nm resolution, and the UV cross sections from
> the Lieden database have the same resolution of 0.1 nm. Therefore, we
> consider constant resolution of 0.1 nm as [the native resolution]

with a table of the errors from coarser binning (0.2, 0.5, 1 and 10 nm bins,
against the 0.1 nm reference) and a figure captioned "The stellar UV flux of
GJ 436 and cross sections of H2 and H2O, showing the native resolution of
0.1 nm adopted in the model". The cross sections come "from the Leiden
Observatory database (Heays et al. 2017) whenever possible" (READ,
section 2.3, p. 7), which is the same source RATE22's CO photoprocess entry
cites (section 12.3) and the same 0.1 nm table that sits in this workspace as
`VULCAN/thermo/photo_cross/CO/CO_cross.csv` (MEASURED, section 4.8). A 0.1 nm
grid is three to four orders of magnitude coarser than the CO band lines, so
the CO self-shielding transition cannot appear in it at any column.

**The three occurrences of "shield" in the paper are about CO2, not CO.**
MEASURED by searching the full text; the relevant passage is section 5.3
(p. 33): "We confirm the analysis in Venot et al. (2013) that although the CO2
abundance is not directly influenced by the temperature dependence of CO2
photolysis, the shielding effects can impact other species. As CO2 absorbs
more strongly with increasing temperature, the UV photosphere is lifted to
lower pressure." That is continuum shielding by a temperature-dependent cross
section, which is the `water_photolysis.f90` pattern, not the
`lyman_werner.f90` one.

**Temperature-dependent cross sections are carried for eight molecules, and CO
is not among them.** READ, section 2.5 (p. 8): "we have included
temperature-dependent photoabsorption cross sections of H2O (EXOMOL), CO2
(Venot et al. 2018; with 1160 K from EXOMOL), NH3 (EXOMOL), O2 (Frederick &
Mentall 1982; Vattulainen et al. 1997), SH (Gorman et al. 2019), H2S (Gorman
et al. 2019), COS (Gorman et al. 2019), and CS2 (Gorman et al. 2019) in the
current version of VULCAN." The same section quotes the justification for
leaving the rest at room temperature, and it is a statement about wavelength
integrals: "Heays et al. (2017) suggested that as temperatures increased by a
few hundred kelvin, the excitation of vibrational and rotational levels
(limited to v <= 2) in many cases only causes minor broadening of the cross
sections and does not alter its wavelength integration."

**That justification is exactly the one that fails for CO**, and this is worth
stating because it is the same physics as section 10.4. A line-shielded rate
does not depend on the wavelength integral of the cross section; it depends on
where the lines are and how wide they are, which is what "minor broadening"
changes. Visser et al. measure the size of it: raising `T_ex(CO)` and `b(CO)`
changes the rate "by as much as a factor of three for temperatures above
100 K" (section 10.4), while leaving the integrated cross section alone. So
the argument that licenses a fixed cross section for most molecules is
precisely the argument that does not license one for CO.

### 13.3 What this implies for the section 10.6 pattern

Three consequences, and the first is the one that changes how section 10.6
should be read.

1. **No consulted code of this kind implements CO self-shielding.** VULCAN
   2017 has no photochemistry at all; VULCAN 2021 has continuum photochemistry
   on a 0.1 nm grid; Moses et al. 2011 and Venot et al. 2012 use bare cross
   sections; RATE22 stores an unshielded `alpha` and sends the user to Visser
   et al. So section 10.6 does not describe a pattern to be copied from a
   reference implementation. It describes one that would have to be written,
   and **the closest working precedent in reach is EXHALE's own
   `lyman_werner.f90`**, which already does for H2 exactly what CO needs: a
   band flux divided by a mean photon energy, an effective cross section read
   from a two-dimensional table of column and temperature, a cell mean taken
   by refined quadrature, and an edge clamp above the table's maximum column.
2. **The 0.1 nm cross-section table cannot be the input.** Section 10.7 used
   it for the optically thin rate, which is legitimate because the thin rate
   depends only on the integrated cross section. It cannot be used for
   anything shielded, and Tsai et al.'s Appendix B is the published statement
   of why: 0.1 nm is the native resolution of the Leiden data, and the CO
   bands are not resolved in it at any binning.
3. **The temperature argument has to be made explicitly, not inherited.** The
   Heays criterion Tsai et al. quote is about wavelength integrals and is
   sound for continuum absorbers; for CO the relevant statement is Visser
   et al.'s, that `T_ex` and `b` move the rate by up to a factor of three and
   that their own grid stops near 500 K for want of data on the `v'' = 1`
   level. Any CO term in this code operates at 800 to 3000 K, above that
   limit, and has to carry the limitation in the comment beside it.
