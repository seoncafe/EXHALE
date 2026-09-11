# B3b-CO: the one-sided CO destruction model, as an implementable design

Written 2026-09-06 for the advisor's approval before B3b-CO is briefed to a
worker. Read-only item: nothing in `src/` was changed, nothing was built and
nothing was run for it.

Inputs, all read for this document:
`docs/b1_target_system_20260906.md` section 5 (T5.1, T5.2, decision 7 as
revised on 2026-09-06);
`docs/co_destruction_rates_literature_20260906.md` sections 3, 10, 12, 12.5
and 13;
and the source files
`src/modules/lower_atmosphere/diffusive_photochemistry.f90`,
`lyman_werner.f90`, `water_photolysis.f90`, `oxygen_rates.f90`,
`molecular_reaction_heat.f90`,
`src/modules/radiation/util_ion_eq.f90`,
`src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/files_IO/write_output.f90`,
`src/modules/files_IO/write_setup_report.f90`,
`src/modules/time_step/certification.f90`.

Every number below is labeled READ (transcribed from a source file or a
publication as quoted in the literature document) or MEASURED (arithmetic
carried out here on such numbers). No new run was made, so nothing here is
measured on a wind.

Position, in one paragraph. The code today removes CO with a thermodynamic
switch that carries no rate: where the CO to C + O equilibrium constant
forbids CO, the transported CO is set to its equilibrium value
(the thermodynamic-ceiling block of `limit_to_element_budget`,
`diffusive_photochemistry.f90`; that block was deleted on 2026-09-06, item
CEILING-DEL, so every citation of it below is a record of the state this
design was written against). B1 decision 7 excludes any run in
which that switch acts during accepted physical integration. Section 12.5 of
the literature document establishes that the physics which actually removes
CO in that region is `He+` charge transfer, with unshielded photodissociation
taking over higher in the wind, and that both are published with usable rate
data. This design replaces the switch with those two rates, states the domain
in which a destruction-only model is legitimate, and turns the switch's
cumulative record into a record of that domain.

---

## 1. The two channels as equations in the CO carrier's source

### 1.1 Where the CO row lives today, and what it contains

READ. CO is one of the five carriers of the transport operator
(`ic_H2, ic_OH, ic_H2O, ic_CO, ic_Hp`). Its balance row is assembled in
`carrier_residual` (`diffusive_photochemistry.f90`), which is a
backward-Euler discretization of

```text
d n_CO / dt + (1/r^2) d[ r^2 ( n_CO v + Phi_CO ) ] / dr = S_CO ,
```

with `Phi_CO` the molecular-diffusion plus settling flux of
`carrier_face_flux` and `S_CO` the chemical source returned by
`carrier_source` (same file). That source was, when this was written,

```text
src(ic_CO) = 0.0d0        ! chemically frozen (decision D4)
```

so the row was transport only, and the CO the operator produced was then
clamped from above in `limit_to_element_budget`, first by the
thermodynamic ceiling and afterwards by the carbon and oxygen
element budgets. The element clamps are conservation and stay;
the thermodynamic ceiling is what this design replaces.

Two further facts of the present layout that the design has to respect.

- The carrier step holds the ion stages frozen. `n(He+)` is `cbg_nheii(j)`,
  filled once per step by `carrier_state` from the `f_sp` the step was handed
  (the `cbg_*` declarations in the module head, filled in `carrier_state`);
  the ionized stages of carbon and oxygen are
  `cbg_nCion(j)` and `cbg_nOion(j)` in the same block.
- The nuclei a carrier does not hold are shared over that element's stages in
  proportion to what they already held (`carrier_write_back`, the oxygen block
  and the carbon block at the end of the routine). So a nucleus released by CO
  destruction returns to its element's pool inside the carrier step, and the
  ionization sweep that follows re-solves the stage partition of that element
  from its own balance.

Without `Molecular carrier transport` there is no CO row at all: CO is the
algebraic equilibrium value `co_equilibrium_density(nm_el(C), nm_el(O), T)`
assigned in `ioniz_eq` (`ionization_equilibrium.f90`). Section 4 states what this
design does there. In practice the two are one case, because a run with
`Oxygen chemistry: True` turns the carrier transport on by default
(READ, `backup/regression/oxygen_chemistry/EXHALE_resolved.out:18`,
`carrier_transport T`, from an `input.inp` that sets no transport key).

### 1.2 Channel D1: `He+ + CO -> C+ + O + He`

Rate, READ from `docs/co_destruction_rates_literature_20260906.md`
section 12.3, which quotes RATE22 entry 4068 of `rate22_final.rates`
(Millar, Walsh, Van de Sande and Markwick 2024, A&A 682, A109):
`alpha = 1.60e-9`, `beta = 0`, `gamma = 0`, method `M` (measured), accuracy
`A` (better than 25 per cent), `T = 10 to 41000 K`. In the RATE22 rate form
`k = alpha (T/300)^beta exp(-gamma/T)` this is temperature independent:

```text
k_D1 = 1.60e-9 cm^3 s^-1 .
```

The upper bound `41000` is a file default and not a validated ceiling
(literature document section 12.3, MEASURED: it appears on 6562 of the 8767
entries and nowhere in the paper). The extrapolation is defensible for this
entry alone, because the value is the Langevin rate of an exothermic
ion-neutral reaction and KIDA's own rule permits extrapolating exactly that
case (literature document section 5.2). The code comment must say so.

Term in the source:

```text
S_CO  <- S_CO  - k_D1 * n(He+) * n_CO ,      n(He+) = cbg_nheii(j) .
```

`n_CO` is the trial carrier density `nc(ic_CO)`, so the term is a diagonal
contribution to the CO row and the Jacobian sees it; `n(He+)` is frozen
background for the step, like every other ion density the rows use.

Products, and where they go.

- **O** goes to free atomic oxygen. Nothing is added for it: `n_o0` in
  `carrier_source` is the closure
  `max(nO_free - cbg_nOion - nc(ic_OH) - nc(ic_H2O) - nc(ic_CO), 0)`, so the
  oxygen nucleus leaves CO and appears as free atomic O in the same evaluation
  of the same cell. The oxygen element total does not move, which is the
  conservation statement the test asserts.
- **C+** goes to the carbon element pool. Carbon is carried today as the metal
  element `iel_C` with its three stages `nm(:,melem_i0(iel_C)+0..2)`
  (`species_table.f90`, the `iel_C`, `mion_elem`, `mion_ethr` and `melem_i0`
  parameter arrays); inside the carrier step those stages
  are frozen as the sum `cbg_nCion` and `nC_free`, and `carrier_write_back`
  hands the carbon not held by CO back to the stages in proportion to what
  they held. So the carrier step conserves the carbon nucleus and does not
  decide its charge; the ionization sweep that runs next re-solves C I, C II
  and C III from its own balance. The reaction's charge bookkeeping (He+
  neutralized, one carbon ionized) is therefore not resolved inside the
  carrier step, and this must be written at the term rather than left to be
  discovered: the step moves nuclei, the sweep moves charge.
- **He** is background and is not written.

Reaction energy, from the one formation-energy table
(`species_formation_energy`, `molecular_reaction_heat.f90`, whose
reference state is every element as a neutral ground-state free atom at rest
with the free electron at zero):

```text
q_D1 = eps(He+) + eps(CO) - eps(C+) - eps(O) - eps(He) .
```

MEASURED from the table's own entries:
`eps(He+) = e_th_HeI = 24.587389011 eV` (READ, `parameters.f90` through
the ionization-energy table of `species_formation_energy`,
`molecular_reaction_heat.f90`), `eps(CO) = oxygen_formation_energy_eV(ith_CO)
= -11.11569 eV` (READ, `oxygen_rates.f90` header, the 0 K formation energy from
the Shomate `F` coefficients), `eps(C+) = mion_ethr(C I) = 11.26 eV` (READ, the
`mion_ethr` parameter array of `species_table.f90`), `eps(O) = eps(He) = 0` by
the reference state. Hence

```text
q_D1 = +2.2117 eV per event, exothermic.
```

Sign and recipient. Positive `q` is heat released, the convention of
`molecular_reaction_energy_eV`. The recipient is the products'
kinetic energy, so the term is thermal heating and enters `gamma_chem` of
`molecular_chemical_heating`, where every other collisional channel
of the network already deposits. No channel-specific energy is written down:
the value is formed from the species table, which is the rule the module
header of `molecular_reaction_heat.f90` states and which B3b-ENTH extends.

One implementation constraint, found while reading, that the brief must carry.
The reaction table `mreac_react` / `mreac_prod` is a `parameter` array of
`f_sp` keys, and the key of C II is `6 + melem_i0(iel_C) + 1`, which is a
runtime quantity (it depends on which metals the run carries). So this channel
cannot be a compile-time row of that table as written. The design's answer:
add a function in the same module that forms `q_D1` at run time from
`species_formation_energy` calls with the runtime carbon key, so the energy is
still built from the one table and still transcribed nowhere. Making
`mreac_*` runtime-initialized instead is the alternative, and it is open
choice 5.

### 1.3 Channel D2: CO photodissociation on the 912 to 1110 A beam

The rate, in the form Visser, van Dishoeck and Black (2009), A&A 503, 323
write it (their Eq. 2, READ through literature document section 10.3):

```text
k_i = chi k0_i Theta_i exp(-gamma A_V) ,
```

an unattenuated rate, a line shielding function of the star-ward columns, and
a continuum term. EXHALE's translation, in the shape of
`lyman_werner_dissociation_rate` (`lyman_werner.f90`):

```text
k_CO(j) = ( F_LW_beam / e_lw_photon_erg ) * sigma_CO_band
          * Theta( N_CO, N_H2 ) * exp( -tau_cont ) ,
```

with, term by term:

- `F_LW_beam = fuv_band_flux(ib_LW)` (`util_ion_eq.f90`), which is
  `F_LW_star` multiplied by `dayside_dilution()`. This is the fourth absorber
  of a beam that already carries three: H2 in the Lyman and Werner lines, and
  H2O and OH in their continua. The "one band, one incident flux, one beam"
  rule of the `oxygen_rates.f90` band header applies without modification,
  because Visser's integration range 911.75 to 1117.80 A covers the LW band
  and its longest CO band sits at 1076.08 A, inside it (READ, literature
  document section 10.2).
- `e_lw_photon_erg = 1.96483e-11 erg = 12.2635 eV`, the mean photon energy of
  a flat `F_lambda` 912 to 1110 A band (READ, the `e_lw_photon_erg` parameter
  of `lyman_werner.f90`). **Superseded on 2026-09-06 (item LW-NORM-B):** the
  band is 912 to 1201 A and the constant is `1.88021e-11 erg = 11.7354 eV`
  (READ 2026-09-07 from the same declaration).
- `sigma_CO_band = 1.577e-17 cm^2`, the flux-weighted mean CO dissociation
  cross section over that band (MEASURED in literature document section 10.7
  from `VULCAN/thermo/photo_cross/CO/CO_cross.csv`). Equivalently
  `k_thin = 8.066e-7 s^-1` per `erg cm^-2 s^-1` of band flux.
- `Theta(N_CO, N_H2)`, the Visser shielding function, section 1.4 below.
- `tau_cont`, the H2O and OH continuum depth of the same interval,
  `tau_b(:,ib_LW)` (built by `fuv_lw_photon_field`, `util_ion_eq.f90`), which is identically zero for a
  run without the oxygen chemistry. Note that CO photodissociation exists only
  in runs that have the oxygen chemistry, so this term is always live here.

**The conversion of Visser's `k0` to this beam, stated.** Visser's
`k0 = 2.592e-10 s^-1` (Table 5 header, READ) is the unshielded rate in the
Draine (1978) field, and their section 3.2 states plainly that the rate
depends on the choice of field, quoting 2.0, 2.0 and 2.3e-10 for Habing,
Gondhalekar and Mathis. EXHALE has a band flux and not an interstellar field,
so the code does **not** carry `k0` at all: it carries the flux-weighted mean
cross section above, and multiplies it by the band photon flux of the run's
own beam. The two routes must agree, and that is a test rather than an
assumption: `k0 / k_thin = 2.592e-10 / 8.066e-7 = 3.213e-4 erg cm^-2 s^-1`
(MEASURED here), which is the 912 to 1110 A energy flux the Draine field must
have for the two statements to be the same statement. Increment (ii) computes
that integral from the published Draine field and requires agreement inside
Visser's own stated absolute accuracy of about 20 per cent (their section 3.5,
READ). If it disagrees by more than that, the cross-section route is the one
to keep (it is the run's own beam) and the disagreement is recorded at the
term, because it is then a statement about the 0.1 nm cross-section table.

**Why the cross-section route and not `k0` with a color correction.** Visser's
`k0` is an integral of the same molecular data over a specified field; using
it with EXHALE's flux would require the flux's spectral shape inside the band,
which the code does not carry (the band is flat `F_lambda` by construction
everywhere else in this beam). Using the band-mean cross section is the same
approximation the H2, H2O and OH terms of this beam already make, and it keeps
the four absorbers of one beam on one convention.

**The one-sided reduction.** Visser's section 3.3 states that Eq. 2 assumes
radiation from all directions and that "if this is not the case, such as for a
cloud irradiated only from one side, k0,i should be reduced accordingly"
(READ). EXHALE applies that reduction already, in the same place and by the
same factor as every other absorber of the beam: `fuv_band_flux` multiplies by
`dayside_dilution()`, which the `oxygen_chemistry` case reports as 0.500, so
343.0 erg cm^-2 s^-1 at the planet becomes a 171.5 erg cm^-2 s^-1 beam (READ,
`output/Lyman_Werner.txt` header). No further factor is applied to the CO
term, and the comment says why: applying one would give the four absorbers of
one beam two different geometries. This is open choice 4.

**The cell mean, not a face value.** The rate is evaluated as the mean over
the cell, by the quadrature `lyman_werner_dissociation_rate_cell_mean`
(`lyman_werner.f90`, through `column_cell_mean_of_a_cross_section`) uses:
composite three-point Gauss-Legendre on
segments cut geometrically in the column, at most `seg_dex = 0.05` decades of
column and `seg_dtau = 0.5` of continuum depth wide, with a head segment for a
cell whose star-ward column is zero, and at most 256 segments. The reason is
the same reason stated there and it is quantitative: `Theta` falls by four
decades along the `N(CO)` axis of Table 5 and by more than six along the
`N(H2)` axis (READ, literature document section 10.6), so one cell at the CO
front carries a large fraction of a decade and a face value is wrong in one
direction across the whole cell. The CO term therefore gets its own cell-mean
function of the same shape, sharing the quadrature constants with the H2 one.

**Products and energy.** `CO + h nu -> C + O`, the neutral channel of
RATE22 entry 8259 and of Visser. The oxygen goes to `n_o0` and the carbon to
the carbon pool, exactly as in D1, and neither is charged. The photon pays the
bond and the remainder is fragment kinetic energy, which is the B3b-PE
recipient rule and is already how this beam's other absorbers deposit
(the `e_lw_fragment_erg` parameter of `lyman_werner.f90`;
`photolysis_threshold_erg` (`oxygen_rates.f90`) and the `heat_fuv` sum for H2O
and OH in `ioniz_eq` (`ionization_equilibrium.f90`)):

```text
threshold  D0(CO) = eps(C) + eps(O) - eps(CO) = 11.11569 eV   (READ)
deposit    <h nu> - D0(CO) = 12.2635 - 11.1157 = 1.1478 eV per event  (MEASURED)
```

Two facts about that number belong in the code comment. First, the threshold
sits at 1115.4 A (MEASURED), which is longward of the band's red edge at
1110 A, so every photon of the band is above threshold, but only by 0.054 eV
at the edge. Second, CO does not absorb flat across the band: it absorbs in 37
discrete bands between 912.70 and 1076.08 A (READ, Visser Table 1 through
literature document section 10.2), the strongest of which sits at 11.52 eV.
The flat-band mean is therefore an approximation of the deposited energy, not
of the rate; it is the same convention the H2 term uses, and it is open
choice 3.

### 1.4 The shielding function: which Visser table, and its domain

**Table 5** is the reference set: `b(CO) = 0.3 km/s`, `T_ex(CO) = 5 K`,
`T_ex(H2) = 51.5 K`, `N(12CO)/N(13CO) = 69` (READ, literature document
section 10.3). Only its 12CO block is used: EXHALE carries no isotopic
structure, so the five other blocks of every table are never called, and the
isotope-ratio difference between Tables 5 and 8 is immaterial to the 12CO
block (READ, section 10.5). The grid is
`log10 N(12CO) = 0, 13, 14, 15, 16, 17, 18, 19` by
`log10 N(H2) = 0, 19, 20, 21, 22, 23`, so 8 by 6 for the block the code needs.
Interpolation is bilinear in the two logarithms with the edge value returned
outside the grid, which is the `h2_shield_locate` convention
(`h2_shield_locate`, `h2_self_shielding_table.f90`) and which Visser's own section 5.1
justifies at the deep end: beyond the tabulated range "photodissociation at
these depths is typically already so slow a process that it is no longer the
dominant destruction pathway for CO" (READ).

**Why Table 5 and not Table 6, and the honest statement of the domain.**
Neither table is at the excitation temperature of this gas. Visser raise
`T_ex(CO)` from 4 to 512 K and stop there, for a stated data reason: the
`v'' = 1` level of 12CO lies 2143 cm^-1 above `v'' = 0` and starts to be
populated near 500 K, and "no data are available on dissociative transitions
out of this level" (READ, their section 4.4). The molecular layer of this code
runs at 800 to 3000 K. Their section 5.2 says what that costs: the Table 5
functions "can easily give photodissociation rates off by a factor of two when
applied to a high-density, high-temperature PDR" (READ). So:

- Table 5 is adopted as the shipped table, because it is the paper's reference
  set and because the gas sits between Table 5 and Table 7 in Doppler width
  (MEASURED, literature document section 10.4: `b(CO)` is 0.42 km/s at 300 K
  and 1.34 at 3000 K, against 0.3 and 3.0 for the two tables) but above every
  table in excitation temperature, so no table is closer on both axes.
- `T > 512 K` is a **stated out-of-domain state** of the closure, recorded per
  cell as a B6 item of class 4.2 (out-of-domain closure), with the H3+ cooling
  model's own domain record as the precedent (`read_validity_states`,
  `certification.f90`).
  It is informational and does not by itself invalidate a result, on the same
  reading decision 8 gave the H3+ record, but it must be printed, because a
  factor of two in `Theta` is a factor of two in the rate.
- The comment at the table states the limit and its cause, in the wording of
  Visser's own sentence, so that a later reader does not take the table for a
  general-purpose one.

Carrying two tables (5 and 6, that is `T_ex(CO) = 5 K` and 50 K) and reporting
the spread as the model's own uncertainty is open choice 1.

### 1.5 The assembled CO row

```text
S_CO(j) = - [ k_D1 n(He+)_j + k_CO(j) ] n_CO ,
k_D1    = 1.60e-9 cm^3 s^-1 ,
k_CO(j) = < ( F_LW_beam / e_lw_photon_erg ) sigma_CO_band
             Theta(N_CO, N_H2) exp(-tau_cont) >_cell ,
```

with no formation term, which is what "one-sided" means (T5.1), and with the
thermodynamic ceiling of `limit_to_element_budget` removed. Heat deposited per
unit volume per unit time, positive into the gas:

```text
Gamma_CO = [ q_D1 k_D1 n(He+) + ( <h nu> - D0(CO) ) k_CO ] n_CO ,
```

`q_D1` and `D0(CO)` both formed from `species_formation_energy`, the first
through the reaction difference and the second through the same difference for
the photolysis channel, so the two channels of one molecule cannot disagree
about how much energy its bond holds.

### 1.6 The new column in the beam's bookkeeping

`Theta` is a function of two columns, and only one of them exists today.
`fuv_lw_photon_field` (`util_ion_eq.f90`) builds `NH2col`, `NH2Ocol` and
`NOHcol` with `calc_column_dens_one` and carries `tau_lya` beside them. The
CO term needs `NCOcol`, built by the same routine from `n_CO`, on the same
radial points, so that `NCOcol(j)` is the column at the inner face of cell `j`
and `NCOcol(j+1)` the column at its star-ward face, which is the convention
the cell mean depends on. Two consequences:

- `fuv_lw_photon_field` gains one input array (`nCO`) and two outputs
  (`NCOcol`, `k_CO`), and `carrier_photolysis` gains `cph_kco(j)`
  beside `cph_klw(j)`, filled in the same call, so that the CO rate is
  rebuilt from the current carrier columns at the top of each step exactly as
  the H2 rate is. This is what makes the CO layer shield itself inside the
  relaxation rather than only between Picard passes.
- The photon ledger of the band gains a CO term. The LW band's absorbed-photon
  accounting in `write_output.f90` currently splits the interval between the
  H2 lines and the H2O and OH continua; a fourth absorber that removes photons
  from the same beam has to appear in that split or the ledger stops closing.
  Whether the CO line absorption is also subtracted from the beam driving the
  other three (that is, whether CO shields H2 as well as itself) is open
  choice 2.

---

## 2. The domain (T5.2), and what replaces the ceiling's record

### 2.1 The inequality, evaluated on the state

T5.2 requires `tau_dest << tau_res << tau_form`. The first inequality is what
makes a destruction-only model legitimate; it is evaluated per cell from the
same arrays the row uses:

```text
tau_dest(j) = 1 / ( k_D1 n(He+)_j + k_CO(j) ) ,
tau_res(j)  = min( tau_adv(j), tau_diff(j) ) ,
tau_adv(j)  = r(j) R0 / ( |v(j)| v0 ) ,
tau_diff(j) = ( dr_j(j) R0 )^2 / ( D_CO(j) + K_zz(j) ) ,
```

with `D_CO(j) = carrier_diffusion_coefficient(j, ic_CO)` and `K_zz(j) =
kzz_cell(j)`. The two residence times are the definitions
`write_oxygen_chemistry` (`write_output.f90`) already writes for H2, with the CO diffusion
coefficient in place of the H2 one; using the code's own definitions is what
makes the domain record comparable with the timescale columns the run already
prints.

One departure from T5.2's own wording, stated rather than hidden: T5.2 writes
`tau_res = min(L/|v|, L^2/(D + K_zz))` with `L` the local gradient scale,
while the code's `tau_adv` uses `r` and its `tau_diff` uses the cell width
`dr_j`. The measured table of literature document section 12.5 was built with
the code's definitions, so adopting them keeps the design's numbers and the
run's numbers on one convention; the gradient-scale form is open choice 6.

### 2.2 The record, and what "off" means where it fails

The model is **in domain** in a cell when `tau_dest <= f_dom * tau_res`, with
`f_dom` a stated constant. Recommended value `f_dom = 0.1`, one decade, which
is the weakest reading of "much less than" that is still a statement; the
measured profile of literature document section 12.5 puts every cell above
`r = 1.05` at `tau_dest/tau_res` of 1.8e-1 to 9.6e-5 and every cell below
`r = 1.04` at 1.5 to 6.7e2, so a threshold anywhere between 1e-1 and 1 selects
the same boundary on that state. The value is a stated parameter of the model
and not a tuning knob: it is written at the site with the measurement above
beside it.

**What happens where the test fails.** The rates are physics and are evaluated
everywhere; what fails below the helium ionization front is the omission of
formation, not the destruction terms. So the model does not switch off any
term. It records the cell:

- a B6 class 4.2 entry (out-of-domain closure), counted per cell and reported
  with the count of cells, the worst ratio and the radius at which it occurs;
- and the physics of such a cell is stated where it is counted: CO there is
  destroyed on a time far longer than the residence time, so over the run's
  own physical time the omitted formation cannot rebuild what the omitted
  destruction fails to remove, and the transported value stands. That is the
  correct behavior of a transport operator in a quenched layer, and it is the
  opposite of what the ceiling does there.

The second inequality of T5.2, `tau_res << tau_form`, is what the omission of
formation has to be argued from, and it can be recorded with published data
without adding a formation term to the row. RATE22 entry 8597 (READ, from
`~/RT_Codes/UMIST/rate22_final.rates`, parsed here) is the radiative
association `C + O -> CO + PHOTON`, `alpha = 4.69e-19`, `beta = 1.52`,
`gamma = -50.5`, `T = 10 to 14700 K`, calculated, accuracy `B`; entry 7167 is
`C + OH -> CO + H`, `1.00e-10`, but its stated range is 10 to 300 K and it is
out of range here. Evaluating `tau_form = 1/(k_8597 n_O)` for the record, as a
diagnostic that never enters the row, is open choice 7.

### 2.3 The ceiling's removal and the record that replaces it

**Done on 2026-09-06 (item CEILING-DEL).** Everything specified in this
section is in the code, and the measurement that settles the physics is that
in the current tree the ceiling never fired: on `oxygen_chemistry` at 1000
steps the before binary reported `co_ceiling_cells 0`,
`co_ceiling_applications 0` and `co_ceiling_CO_removed 0.0` (MEASURED, from
`EXHALE_resolved.out` and the end-of-run line), so every numerical output
column is byte-identical across the deletion and `log10 Mdot` stays 9.61. With
the published destruction rates in the row and the He+ sink of B3b-CO2 in the
sweep, the carried CO never rises above its chemical equilibrium anywhere on
the grid, and the bound had nothing left to bound. The specification as
written follows.

Removed (done 2026-09-06, item CEILING-DEL): the thermodynamic-ceiling block
of `limit_to_element_budget` (`diffusive_photochemistry.f90`), that is the
`co_equilibrium_density` call, the comparison, the clamp, and the four
cumulative counters `co_ceiling_applications`, `co_ceiling_cells_hit`,
`co_ceiling_CO_removed` and `co_ceiling_seen` with their attempted-ledger
twins and their accessors `carrier_co_ceiling_record`,
`carrier_co_ceiling_record_attempted` and `carrier_co_ceiling_cells`. The
element clamps that follow in the same routine stay: they are conservation,
not thermodynamics.

Kept, and re-pointed: the four consumers of that record.
`read_validity_states` (`certification.f90`) read it into
`rep%n_active_unvalidated_physics`, which is B6 class 4.1;
`write_resolved_config` (`write_setup_report.f90`) printed the three counters
at the end of a run;
`write_oxygen_chemistry` (`write_output.f90`) wrote the `co_ceiling_cells`
header of `Oxygen_chemistry.txt`; `EXHALE_main.f90` and `attempted_step.f90`
carried it through the step ledgers. Each of these takes the domain record in
its place:

- `rep%n_active_unvalidated_physics` loses its CO contribution entirely. The
  CO destruction physics is published and rated, so it is no longer unvalidated
  physics; it becomes a closure with a domain, that is class 4.2, counted with
  the H3+ items. This is the substantive change to the run's validity state
  and it must be reported as such.
- The setup report and the `Oxygen_chemistry.txt` header print, in place of
  the three ceiling counters: the number of cells out of domain, the worst
  `tau_dest/tau_res` and its radius, the number of cells above the 512 K
  `T_ex` limit of the shielding table, and the `f_dom` in force.
- The marks keep the two-family indexing of `ledger_family`
  (initialization or continuation against accepted physical steps) that the
  ceiling record already has, for the same reason: a cell first out of domain
  during relaxation is not a cell the physical history touched.

**B1 decision 7 and the exclusion rule.** T5.3 excludes any configuration in
which the ceiling is active during accepted physical integration. With the
ceiling gone that observable no longer exists, and the exclusion is lifted by
construction, exactly as T1.10 says the oxygen exclusion is lifted when T1.9
exists. What replaces it is not an exclusion but a reported domain: a run in
which the destruction model is out of domain in some cells is still a physical
result, with that record attached to it. Whether the advisor wants a stronger
rule (for example, refusing certification when out-of-domain cells hold more
than a stated share of the carbon) is open choice 8.

---

## 3. Increments and tests

Four increments, each RED before and GREEN after, each with its own test
program under `src/tests/`. The naming follows the existing suites
(`carrier_retry`, `a2_m1`, `spectrum_type`).

### (i) The He+ channel alone

Changes: `carrier_source` gains the D1 term and `nC_free` in its argument
list where the reference density needs it; `molecular_reaction_heat` gains the
runtime `q_D1` function and `molecular_chemical_heating` gains the term. The
photodissociation term is absent in this increment, and the ceiling is still
in place, so the increment is additive and its effect can be measured alone.

Test suite `co_destruction`, part 1:

1. `k_D1` returned at `T = 300, 1500, 3000 K` equals `1.60e-9` to the bit at
   all three, which is what `beta = gamma = 0` means, and is the RED check
   against any accidental Arrhenius factor.
2. `q_D1` from the table equals `+2.2117 eV` to 1e-6 eV, and the same value is
   recovered as `eps(He+) + eps(CO) - eps(C+)` computed independently in the
   test from `species_formation_energy`, so a change in any one entry moves
   both together.
3. Carbon and oxygen budgets: a single cell, no transport, one step of the
   source alone. `n_CO + n_C(all stages)` and `n_CO + n_OH + n_H2O + n_O(all
   stages)` are each constant to round-off (relative change below 1e-14) over
   a step in which `n_CO` falls by half.
4. The energy identity of T1.5 on the same cell: the heat deposited equals
   the drop in `sum_s n_s eps_s` from the CO destroyed, to round-off.

Impact measurement: `oxygen_chemistry` and `mol_carrier` on scratch copies,
before and after. `mol_carrier` runs without the oxygen chemistry and carries
no CO, so it must be bitwise identical; if it is not, the term is being
evaluated where no CO exists.

### (ii) The CO column and the shielded photodissociation

Changes: a new `co_self_shielding_table.f90` carrying the 12CO block of Visser
Table 5 as a `parameter` array with its two axes, in the shape of
`h2_self_shielding_table.f90`; a `co_photodissociation.f90` (or an addition to
`lyman_werner.f90`) with the point rate and the cell mean; `NCOcol` and
`k_CO` in `fuv_lw_photon_field`; `cph_kco` in `carrier_photolysis`; the D2
term in `carrier_source` and its fragment heat.

Test suite `co_destruction`, part 2:

5. **Table transcription.** The 48 entries of the 12CO block are reproduced at
   the grid points to the precision they are printed with. Four of them are
   already independently quoted and serve as the anchors: `1.150e-3`,
   `1.941e-4`, `7.329e-5` and `1.437e-5` at the four corners bracketing
   `log N(CO) = 18.96, log N(H2) = 21.80` (READ, literature document
   section 10.7), together with the axis-end values `Theta = 1.0` at the
   origin of both axes, `5.24e-4` at `log N(CO) = 19` with `N(H2) = 0`, and
   `3.9e-7` at `log N(H2) = 23` (READ, section 10.6).
6. **Interpolation and clamp.** Log-bilinear interpolation at the same point
   returns `Theta = 2.6e-5` to 1 per cent of the value quoted in
   section 10.7; above `log N(CO) = 19` and above `log N(H2) = 23` the edge
   value is returned and not an extrapolation.
7. **The band-flux conversion.** With `Theta = 1` and `tau_cont = 0` the rate
   equals `8.066e-7 s^-1` per `erg cm^-2 s^-1` of band flux to 1e-3 relative,
   and at the `oxygen_chemistry` beam flux of 171.5 erg cm^-2 s^-1 it equals
   `1.383e-4 s^-1` (READ, section 10.7). Separately, the ratio `k0/k_thin =
   3.213e-4 erg cm^-2 s^-1` is compared with the 912 to 1110 A energy flux of
   the published Draine (1978) field, computed in the test, and the two agree
   inside the 20 per cent Visser state for their absolute rate. A larger
   disagreement is reported at the site, not tolerated silently.
8. **The cell mean against a fine reference.** A synthetic cell spanning 0.01
   to 10 decades of CO column, with and without continuum depth, at columns
   covering the whole grid: the composite rule agrees with a 4000-segment
   reference of the same integrand to better than 1e-3 relative, which is the
   accuracy `lyman_werner_dissociation_rate_cell_mean` reports for its own
   integrand (2.5e-4 at `seg_dex = 0.05`, READ). If CO needs a finer
   `seg_dex` than H2, the test is what says so and the constant is set from
   it, not from the H2 value.
9. **The column.** `NCOcol` built by `calc_column_dens_one` on a uniform
   `n_CO` reproduces the analytic column to round-off, and
   `NCOcol(j) - NCOcol(j+1)` equals the CO of cell `j`, which is the face
   convention the cell mean assumes.
10. **The out-of-domain mark.** A cell at `T = 600 K` sets the `T_ex` limit
    flag; a cell at 400 K does not.

### (iii) The domain record and the ceiling's removal

Changes: the ceiling block and its counters deleted; the domain evaluation and
its record added; the four consumers re-pointed; the output headers rewritten.

Test suite `co_destruction`, part 3:

11. A synthetic column with `n(He+)` rising through a front: the cells in
    domain and out of domain are exactly those satisfying and violating
    `tau_dest <= f_dom tau_res` on the arrays the run holds, with no gap and
    no overlap at the boundary cell.
12. The record is family-indexed: activity in initialization mode does not
    appear in the physical family, and a rejected trial leaves both empty
    after restoration (the joint test with A3 and B3a that AT-5 (iii) already
    specifies).
13. `certification.f90` reports the CO items under class 4.2 and reports zero
    CO contribution to class 4.1.
14. AT-5 is rewritten against the new observable and its four parts are
    re-stated in B1 section 5; the old parts (i) and (ii), which name the
    ceiling, no longer have a producer and must not be left asserting one.

### (iv) Impact on the regression cases

Two cases touch this code.

**`oxygen_chemistry`.** READ, `docs/Update_EXHALE_stage2.md` (the A1scale and A3b
entries): the case is refused at step 2 under the physical-terms acceptance,
with a CO row at 5.5e-8 beside a clamped cell, 37 retry attempts to the floor
all refused; under A3b it runs its 1000 steps in initialization mode with 998
exhausted intervals, its carriers frozen at the step-2 values, and the final
certification refusing on four carrier rows at 1.0, described there as "the
real imbalance beside a clamped cell, open physics of the excluded CO model,
not a controller defect".

The design's prediction, stated as a prediction and not as a result, because
this item ran nothing: the clamp that produces that imbalance is the
thermodynamic ceiling, and removing it removes the cause **provided** the
element clamps of the same routine do not fire in the same cells. They are a
different clamp with a different cause (a carrier asking for more nuclei than
its element holds), and the ceiling's own record shows it acting in of order
160 cells at every step of this case (READ, the comment at the
`pct_cell_constrained` declaration of `diffusive_photochemistry.f90`), so the
ceiling is the clamp to beat. The measurement is part
of increment (iii): run the case before and after, and report whether the
refusal at step 2 disappears, whether the case reaches step 1000 in physical
mode, and what its `log10 Mdot` is against the 9.61 recorded for the A3b run.
If the refusal survives the removal, the cause is the element clamp and the
finding is reported rather than worked around.

**Outcome (2026-09-06, item CEILING-DEL).** The prediction was wrong in its
premise and the measurement is kept here for that reason. With the rates of
increments (i) and (ii) in the row, the ceiling **never fired at all** on
`oxygen_chemistry` at 1000 steps: `co_ceiling_cells 0`,
`co_ceiling_applications 0`, `co_ceiling_CO_removed 0.0` (MEASURED, before
binary). The "of order 160 cells at every step" was READ from a comment
written before the CO row carried any rate. Removing the clamp therefore moved
no number: every output column is byte-identical and `log10 Mdot` is 9.61
before and after. Whatever produces the step-2 refusal, it is not this
clamp.

**`mol_carrier`.** No oxygen chemistry, no CO, no FUV bands. Required
bitwise identical after every increment. This is the guard that the new terms
are gated on the same condition the rest of the oxygen path is
(`thereis_oxychem`).

Goldens: not refreshed by this item. The `oxygen_chemistry` golden changes
when the physics changes, and it is refreshed once at the end of the series
that contains B3b-CO, with the movement reported.

---

## 4. What this changes for the user

- **No new input key.** The destruction model replaces the ceiling in every
  run with `Oxygen chemistry: True`; there is nothing to switch on. The
  default is therefore "on wherever CO exists", and that is deliberate: the
  thing it replaces was also unconditional, and a switched-off destruction
  model would leave a transported CO indestructible, which is the failure the
  ceiling was introduced to prevent (READ, the comment that opened the
  ceiling block, deleted with it on 2026-09-06: CO
  advected to 1.6 R_p and 28000 K taking the whole oxygen and carbon
  inventory with it). The stated constant `f_dom` of the domain test is a
  parameter of the model in the source, not an input key.
- **`co_ceiling_cells`, `co_ceiling_applications` and `co_ceiling_CO_removed`
  are retired** (done 2026-09-06, item CEILING-DEL) from the setup report,
  from the `Oxygen_chemistry.txt` header and from the certification's class
  4.1, which is left with no producer at all and reports "not produced". In their place: the out-of-domain
  cell count, the worst `tau_dest/tau_res` with its radius, the count of cells
  above the shielding table's 512 K excitation-temperature limit, and the
  `f_dom` in force. Any script reading the old three fields has to be updated;
  they exist only in EXHALE's own outputs.
- **`FUV_bands.txt` gains a CO column**, `N_CO[cm^-2]`, beside `N_H2O` and
  `N_OH`, and `Lyman_Werner.txt` (or `Oxygen_chemistry.txt`, open choice 9)
  gains `k_CO[1/s]` and `Theta` beside `f_shield` and `k_LW`. Both files carry
  a `# columns` schema header, so `examples/exhale_io.py` adapts without
  change.
- **A run whose spectrum carries no 912 to 1110 A flux loses only the D2
  channel.** The `mol_base_handoff` family sets no `Stellar LW flux` and its
  power-law spectrum stops at 911.8 A (READ, literature document
  section 10.7), so `k_CO` is identically zero there and CO is destroyed by
  the He+ channel alone. Those cases carry no CO in the first place, so this
  is a statement about a configuration a user might build, not about the
  matrix.

---

## 5. Open choices for the advisor

1. **One shielding table or two.** Ship Visser Table 5 alone (adopted in the
   text above), or ship Tables 5 and 6 and report the spread between
   `T_ex(CO) = 5 K` and 50 K as the model's own uncertainty in the out-of-domain
   record.
2. **Does CO shield the other absorbers of the beam.** The CO lines remove
   photons the H2 lines and the H2O and OH continua would otherwise see;
   including that coupling keeps the band's photon ledger exact, and excluding
   it keeps the H2 rate bitwise unchanged in every existing case.
3. **The photon energy of the deposited fragment energy.** The flat-band mean
   12.2635 eV (adopted, and the convention of every other absorber of this
   beam), or an oscillator-strength weighted mean over Visser's 37 bands,
   which needs their Table 1 transcribed.
4. **The dilution convention.** Take the one-sided reduction from
   `dayside_dilution()` alone, as every absorber of this beam does (adopted),
   or apply a separate geometric factor to `Theta` because the tables were
   computed for a slab.
5. **How the D1 reaction energy is formed.** A runtime function beside the
   `mreac_*` tables (adopted, because the C II key is runtime), or converting
   `mreac_react` and `mreac_prod` from `parameter` arrays to runtime-filled
   ones so that every channel lives in one table.
6. **The residence time.** The code's `r/|v|` and `dr^2/(D + K_zz)` (adopted,
   and the definitions the measured table of section 12.5 used), or T5.2's
   gradient-scale `L/|v|` and `L^2/(D + K_zz)`, which needs a gradient scale
   the code does not currently form.
7. **Whether `tau_form` is evaluated for the record** from RATE22 entry 8597
   (`C + O -> CO + PHOTON`), as a diagnostic that never enters the row, so
   that both inequalities of T5.2 are measured rather than one.
8. **What the domain record does to certification.** Informational only, as
   decision 8 made the H3+ record (adopted in the text), or a refusal above a
   stated share of the carbon held in out-of-domain cells.
9. **Where the CO rate is written.** A new `CO_photodissociation.txt`, or
   extra columns in `Lyman_Werner.txt` (which is the beam's file and already
   carries `f_shield` and `k_LW`), or in `Oxygen_chemistry.txt` (which is
   where `n_CO` is).
10. **Whether the non-transport path changes at all.** Without
    `Molecular carrier transport` CO has no balance row and stays at
    `co_equilibrium_density`; the destruction model has nothing to act on
    there. Leave that path as it is and record it as an equilibrium closure in
    B6, or refuse the combination of `Oxygen chemistry: True` with carrier
    transport off, which no shipped case uses.

## 6. Advisor decisions on section 5 (2026-09-06)

1. One shielding table in the code: Visser et al. Table 6 (`T_ex(CO) = 50 K`,
   `T_ex(H2) = 501.5 K`, `b = 0.3 km/s`), the set closest to the warm gas of
   the layer; the test records the Table 5 values at the grid points beside
   it as the spread, and the 512 K `T_ex` limit is a stated out-of-domain
   state above it.
2. CO does not shield the beam's other three absorbers: its lines shield CO
   (and H2's lines shield CO through `N_H2` in `Theta`); no CO term enters
   `tau_cont`. Stated at the code site with Visser's construction as the
   reason.
3. The photon energy of the CO dissociation events is the
   oscillator-strength weighted mean over Visser's 37 bands (one constant,
   its derivation in a comment), not the flat-band mean.
4. Dilution: `dayside_dilution()` alone; Visser's reduction for one-sided
   illumination is that factor and nothing else.
5. `q_D1` is formed at run time from the one formation-energy table
   (`species_formation_energy`), not a compile-time row.
6. `tau_res = min(r/|v|, dr^2/(D_CO + K_zz))` from the code's own arrays, as
   the literature document formed it; T5.2's gradient-scale time is a later
   refinement, noted.
7. `tau_form` is evaluated for the record only, from the RATE22 formation
   entries whose reactants EXHALE carries (8597 `C + O -> CO + photon` and any
   other with C, O, OH, H2O reactants), so the record states both sides of
   the ordering.
8. The domain record is informational (B6 category 2, out-of-domain
   closure): the rates are on everywhere, the record marks the cells where
   the one-sided destruction assumption fails; it does not refuse
   certification by itself.
9. `output/FUV_bands.txt` gains the CO column, `k_CO` and `Theta` per cell.
10. Both paths use the same source terms: with carrier transport the CO row
    carries them; without it the local balance evaluates CO with the same
    rates and is recorded as an equilibrium closure. The ceiling and its four
    counters are deleted in both (done 2026-09-06, item CEILING-DEL).

Sequencing: B3b-CO is briefed after B4-6 hands back `util_ion_eq.f90` (the
FUV beam) and after B3c owns the coupled source step, so the new source terms
enter the energy ledger through the one table from the first day.

---

## 7. Landed (2026-09-06), with the numbers

Increments (i), (ii) and (iii) are in the code, the ceiling's deletion
included (item CEILING-DEL, the same day; section 7.4).
Every number below is MEASURED by the implementing item unless marked READ.

### 7.1 What differs from the design above, and why

- **Section 1.3 is restated on the 912-1201 A beam.** The Lyman-Werner
  interval was renormalized to 912-1201 A on the same day this design was
  written (item LW-NORM-B), so the band mean is taken over that beam and
  `<hv>_band = 11.7354 eV`, not the 12.2635 eV of the old 912-1110 A band.
  CO absorbs only over 912-1118 A (37 bands to 1076.08 A, the dissociation
  continuum ending at the 11.1157 eV bond), so the shipped cross section is
  `sigma_CO_band = 1.0160e-17 cm^2` per photon of the WHOLE beam, against
  `1.5496e-17` over 912-1110 A alone. Both are photon-flux weighted means of
  the same 0.1 nm table over a flat `F_lambda`.
- **The cross-check against Visser's `k0` is not the arithmetic section 1.3
  proposed.** The Draine (1978) field is not integrated here: the code
  already carries that field's 912-1110 A photon flux, 1.232e7 cm^-2 s^-1,
  READ from Draine & Bertoldi (1996) Table 1 through `lyman_werner.f90`.
  With `k0 = 2.590e-10 s^-1` (Table 6 header) that implies a band-mean
  dissociation cross section of `2.102e-17 cm^2`, against our `1.5496e-17`
  on the same band: **ours is 26 per cent lower** (MEASURED), a little
  outside Visser's stated 20 per cent absolute accuracy. The cross-section
  route is kept, as the design says it should be, and the disagreement is
  recorded at the term and asserted in the test.
- **A loaded spectrum gets the same flat-`F_lambda` cross section.** The
  design and the brief allow an SED-resolved integral over Visser's band
  positions for `Spectrum type: Load`. It is not done, for two reasons
  stated at the code: every other absorber of this beam is normalized per
  photon of a flat-`F_lambda` band, so a shape-aware CO would put the four
  absorbers of one beam on two conventions; and the integral would have to
  live in `sed_read.f90`, which was not this item's file.
- **`q_D1` is a row of the reaction table, not a runtime function** (a
  departure from decision 5, whose premise is wrong): `mion_fsp(im_CII)` is
  a `parameter`, not a runtime quantity, so the C II key is available at
  compile time and the channel is `ir_D1` of `oreac_react`/`oreac_prod` like
  every other collisional channel. The energy is still formed at run time
  from `species_formation_energy` and is transcribed nowhere, which is what
  decision 5 is for, and the row additionally gets the nucleus and charge
  balance test the table already has.
- **The photon ledger of the band prints CO on its own row.** Decision 2
  keeps CO out of `tau_cont` and out of the other absorbers' shielding, so
  `beam_loss_ph` carries no CO term; counting CO's absorptions in
  `rated_ph`, which is compared against that beam loss, would compare two
  different beams (MEASURED: it drove the LW budget residual from -1.8e-2 to
  -8.2e-2). CO is therefore printed as `co_absorbed_ph`, `co_frac`,
  `co_absorbed_en`, `co_heat_en`, and `co_frac` is the size of the
  approximation decision 2 makes.

### 7.2 The constants, as they ended up

| quantity | value | provenance |
|---|---|---|
| `k_D1` | 1.60e-9 cm^3 s^-1, T independent | READ, RATE22 4068 |
| `sigma_CO_band` | 1.0160e-17 cm^2 per 912-1201 A photon | MEASURED |
| `<hv>` of a dissociation event | 12.8674 eV | MEASURED, Visser Table 1, f-weighted |
| `D0(CO)` | 11.1157 eV | MEASURED, the one formation-energy table |
| fragment deposit | 1.7517 eV | MEASURED |
| `q_D1` | +2.2117 eV | MEASURED, the same table |
| `f_dom` | 0.1 | stated parameter |
| `T_ex` limit | 512 K | READ, Visser sec. 4.4 |

The design predicted `q_D1 = +2.2117 eV` from the same table and the code
returns 2.211694 eV.

### 7.3 The shielding table's own spread

Table 6 is shipped and Table 5 is carried beside it, not called by the rate.
MEASURED over the 48 grid points of the two 12CO blocks: `Theta(Table 5) /
Theta(Table 6)` runs from 0.549 to 22.07, the maximum at
`log N(CO) = 17, log N(H2) = 22`. That is larger than the factor of two
Visser's section 5.2 quotes, because their sentence is about the rate a whole
cloud model gives and this is the worst single grid point; it is the size of
the model's uncertainty on the one axis on which neither table reaches this
gas, and the run's domain record counts the cells above 512 K so the two can
be read together.

### 7.4 The ceiling, deleted the same day (item CEILING-DEL)

B3b-CO left the thermodynamic ceiling of `limit_to_element_budget` and its
four counters in place, because removing them changes `EXHALE_main.f90`,
`certification.f90` and `write_setup_report.f90`, none of which was that
item's to edit. Item CEILING-DEL did it on 2026-09-06: the
`co_equilibrium_density` call, the comparison, the clamp, the four cumulative
counters with their attempted-ledger twins and their three accessors are gone,
`limit_to_element_budget` applies conservation only, and the domain record
stands in the setup report (seven `co_domain_*` keys), the end-of-run report
and the certification, where the CO items are counted under class 4.2 beside
the H3+ records.

The physics finding of that item: **in the current tree the ceiling never
fired.** MEASURED on `oxygen_chemistry` at 1000 steps with the before binary:
`co_ceiling_cells 0`, `co_ceiling_applications 0`,
`co_ceiling_CO_removed 0.0`. This supersedes the expectation stated above that
the ceiling "still cuts CO back above about 3000 K": with the rates of B3b-CO
in the row and the He+ sink of B3b-CO2 in the sweep, the carried CO never
rises above its chemical equilibrium anywhere on the grid. Every numerical
output column of that case is byte-identical across the deletion, and
`log10 Mdot` is 9.61 before and after.

The design's "of order 160 cells" of section 3 (iv) was measured before the CO
row carried any rate. It is no longer the state of the code, and the comment
that carried it was deleted with the clamp.

### 7.5 The He+ ledger closed (2026-09-06, item B3b-CO2)

Section 7.4's list of what B3b-CO left open began with the P0 the item found
in its own physics: `He+ + CO -> C+ + O + He` deposited 2.2117 eV per event
into heating channel 18 while nothing anywhere consumed the helium ion. Row
(2) of `mol_heh_rows` (`System_HeH_mol.f90`) now carries

```text
- k_D1 n_CO n(He+) ,   k_D1 = rk_D1_Hep_CO() = 1.60e-9 cm^3 s^-1 ,
```

with `n_CO` the cell's own reservoir `ieq_cell%n_co`. Both molecular systems
get it from that one row: `System_HeH_mol_metals` calls `mol_heh_rows` and
adds nothing of its own, and the same call reaches the constrained
chemical-equilibrium residual. There is no analytic Jacobian to extend; both
molecular systems are solved by `hybrd1` on a finite-difference Jacobian, and
the derivative that Jacobian builds is asserted directly by the new test.
The turnover scale of the row (`set_mol_turnover_rates`) gains
`k_D1 n_CO n_He`, the same all-of-the-element bound its other terms carry;
without it the He+ row of the shielded base would be divided by a bound far
below the reaction it is mostly made of.

CO is a background density for that row and not one of its unknowns. That is
not an approximation of convenience: the one-sided model gives CO no
formation channel, so CO has no local equilibrium and cannot be a balance row
at all. The reaction is therefore split across the two operators, and each
half reads the other's frozen value: the carrier row removes the CO at the
sweep's He+, the He+ row removes the ion at the carrier's CO, both from
`rk_D1_Hep_CO`. What the split leaves is a lag of one operator, which is what
`dom_cells_HeP` measures.

The products were already where they belong and stay there. The oxygen atom
is supplied by the free-atomic-oxygen closure `n_o0` of `carrier_source`, and
the carbon nucleus by `nC_free - n_CO` in `carrier_write_back`, which shares
it over C I, C II and C III in proportion to what they held; the ionization
sweep then re-solves that partition from its own balance. The carrier
operator moves nuclei and never the element totals, so the element census
(`element_census.f90`, whose `bsp_nC` and `bsp_nO` already count the nuclei
in CO) closes on carbon and oxygen by construction, and does so in the run.
The charge the reaction moves from helium to carbon is not resolved inside
either operator: the sweep sets the carbon stages from photoionization and
recombination, as it does for every metal.

**What it changes, MEASURED** on `oxygen_chemistry` at 1000 steps, one
scratch run each with the CO physics of B3b-CO alone separating the two
binaries:

| quantity | before | after |
|---|---|---|
| `log10 Mdot` | 9.61 | 9.61 |
| exit status | 0 | 0 |
| largest share of the total heating taken by channel 18 | 0.672 at r = 1.081 | 0.035 at r = 1.117 |
| channel 18 at r = 1.0894 [erg cm^-3 s^-1] | 6.90e-7 | 9.75e-9 |
| `n(He+)` at r = 1.0894, relative | 1 | 0.015 |
| smallest `n(He+)` ratio over 1.0-1.4 Rp | 1 | 0.0073 at r = 1.011 |
| `T` at r = 1.0233 | 1156.4 K | 1088.2 K |
| `T` at r = 1.0894 | 1364.5 K | 1279.6 K |
| `dom_cells_HeP` | 281325 | 281382 |
| `dom_cells_out` | 401374 | 440999 |

The heating peak the item was opened for is gone: two thirds of the local
heating at the helium ionization front was being drawn from a helium ion
population that the reaction should already have destroyed, and with the sink
in place the same channel takes 3.5 per cent at its largest. `dom_cells_HeP`
does not fall, and that is the informative part: the counter compares
`k_D1 n_CO` with the recombination rate, and the comparison still says that
in the molecular layer this reaction is the leading He+ loss. What the
counter now reports is the size of the operator split's lag, not a missing
sink; the wording of its header line in `write_output.f90` still says the
latter and is named for correction in the item's report.

Cases with no CO are untouched: `mol_carrier` at 300 steps and `mol_metals`
at 400 steps are byte-identical across the change in every output file.
