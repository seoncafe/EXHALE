# L6: the C I ionization balance is right, its line cooling is not

Item L6 of `docs/PLAN_20260913_lhs_stationary.md`, opened by the out-of-scope
finding of `docs/lhs1140b_stationary_L4b_20260913.md` section 8: carbon is 94
to 99 percent neutral between 1.0 and 1.5 R_p of the LHS 1140 b rung and its
line emission is 79 to 94 percent of the cooling there. The state audited is
`LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k00/output/`. Every
number is MEASURED on this tree unless it is marked READ.

## 1. Verdict

**The large neutral fraction is right and is not the problem. The C I line
cooling is wrong, by one to three decades.**

*The ionization balance.* The C I photoionization rate this state carries was
recomputed independently in Python from the GJ 1132 SED and the published
Verner et al. (1996) cross section, with the photon grid, the columns and the
cell-mean attenuation of the Fortran reproduced, and it agrees with the rate
the code's own composition implies to 0.5 to 1.3 percent (table A). Feeding
that rate, the Badnell recombination coefficient, the Voronov collisional
ionization and the two carbon charge-exchange rows back into the balance by
hand reproduces the code's x(C I) to between 0.004 and 0.14 percent (table B).
Carbon is 94 to 99 percent neutral because recombination beats photoionization
by a factor 30 to 100 there: `alpha n_e` is 5.1e-5 to 5.2e-6 s^-1 against a C I
photoionization rate of 5.0e-7 to 8.8e-7 s^-1. The 11.26 eV edge is a bin edge
of the photon grid to 1.8e-15 eV, the band that drives the rate is essentially
unattenuated (tau = 1.1e-2 at 12 eV at the deepest cell audited, and it is C I
itself, not hydrogen or helium), and the stellar flux between 11.26 and
13.598 eV is 0.566 erg cm^-2 s^-1 against 19.41 above the H I edge. Nothing in
the chain is broken.

*The cooling.* `cool_CI_ne_func` (`src/modules/radiation/Cool_coeff.f90`)
saturates the ground-term fine structure but leaves the rest of the C I
cooling as a coronal, low-density-limit fit. The channels that carry 99.7
percent of it at these temperatures excite the even-parity metastables
2p2 1D2 and 1S0, whose critical densities are 2.0e4 and 1.3e7 cm^-3 against a
local n_e of 1.2e7 to 3.3e7. The coronal form is therefore evaluated one to
three decades outside its density domain, and the cooling it returns exceeds
the LTE ceiling on what a C I atom can radiate in those lines by a factor 20
at 1.86 R_p rising to 1500 at 1.03 R_p (table C). That bound needs no
collision data: statistical equilibrium is monotone in the collision rate and
saturates at the Boltzmann population, so no collisionally excited line can
radiate more than its LTE emission.

*The fix.* The two metastable terms of the fit now carry the exact two-level
saturation, with the ceiling built from NIST level energies and transition
probabilities written at the code site; the two permitted terms stay coronal,
which is right for them (n_cr ~ 1e16 cm^-3). C I cooling at the audited state
falls by the factors of table D, and `cool/heat` over 1.15 to 1.50 R_p, the
band whose energy row refuses the L4b rung, falls from 0.65 to between 0.04
and 0.13, against 0.042 at the same cell in the certified metals-off case
`atomic_scalar_gj1132_kzz1e9/HeH1.60` (READ from
`docs/lhs1140b_stationary_L4b_20260913.md` section 3.2).

*What this does to L4b.* L4b measured the refusal as an energy row that the
stationary solve cannot take back down, and traced the row to a cooling that
stands at 0.53 to 0.66 of the heating where the composition relaxation moves
X_He. That cooling was C I. Whether the rung certifies with the corrected
cooling was NOT measured here; it is the L4b item's own measurement to redo.

## 2. What the code does, and where

| quantity | routine | source of the data |
|---|---|---|
| C I photoionization cross section | `cross_sec.f90` `metal_photoion_sigma`, column 1 | Verner, Ferland, Korista & Yakovlev 1996 Table 1 below E_max = 291 eV, Verner & Yakovlev 1995 above, plus the 1s/2s subshells above 291 eV |
| threshold | `species_table.f90` `mion_ethr(im_CI)` = 11.26 eV | the same value gates the fit and charges the photoelectron |
| attenuated field | `util_ion_eq.f90` `photoionization_field_at_cell_HHe` | `tau_out` from H I, He I, He II, He 2^3S, H2 and every metal ion's own column; cell mean `exp(-tau_out)(1-exp(-dtau))/dtau` |
| rate | same routine, `Pm_loc(i) = sum(int_f sigma/E dE)` | rectangle rule on the photon grid |
| recombination C II + e | `Cool_coeff.f90` `rec_CII_range` -> `alpha_rec_metal('CI')` | Badnell 2006 radiative fit plus the adf09 dielectronic sum |
| collisional ionization | `Cool_coeff.f90` `ion_coeff_CI_range` | Voronov 1997 Table 1 row C 0 |
| charge exchange | `charge_exchange.f90` rows A15, A16 | Stancil et al. 1998 by way of Huang et al. 2023 Table 4 |
| balance row | `ion_residual_core.f90` `metal_rows` | `n0 P + (n0 beta - alpha n1) n_e = 0`, charge exchange added by `cx_add_to_fvec` |
| line cooling | `Cool_coeff.f90` `cool_CI_ne_func` | ground term in statistical equilibrium; the rest a CHIANTI v11.0.2 coronal fit |

Density convention, READ from the residual: the photoionization term multiplies
the neutral density alone, the collisional ionization and the recombination
terms both multiply `n_e`, and the charge-exchange rows multiply the two
reactant densities. `n_e` in that row is the full electron density, metals
included (`metal_electron_sum`).

The photon grid is built in `sed_read.f90` from the SED rows, split at every
ionization threshold of an active absorber. For this run the grid floor is the
He 2^3S threshold, 4.7678 eV, so the 11.26 eV edge is interior and is inserted
exactly; MEASURED, the nearest bin edge to 11.26 eV is 11.2600000000 eV.

## 3. The photoionization rate (question 1)

Independent calculation: the SED file is read and binned exactly as
`sed_read.f90` bins it (geometric-mean interior edges, the table's own end
points, thresholds inserted), the cross sections are transcribed from the same
published fits, the columns outside each cell are accumulated the way
`advance_starward_columns` accumulates them from the cell widths recovered from
the r column, and the attenuation is the same cell mean. The comparison column
"implied" is what the code's own x(C I) requires, given the recombination,
collisional-ionization and charge-exchange coefficients: it inverts the balance
row rather than reading a diagnostic, so it needs no new output.

TABLE A

| cell | r/R_p | T [K] | n_H | n_e | tau(12 eV) | P(C I) independent [1/s] | P(C I) implied [1/s] | implied/independent |
|---|---|---|---|---|---|---|---|---|
| 100 | 1.0250 | 1869 | 7.190e+10 | 3.234e+07 | 1.07e-02 | 5.0486e-07 | 5.0726e-07 | 1.0047 |
| 150 | 1.0615 | 2355 | 1.568e+10 | 3.338e+07 | 3.52e-03 | 5.1842e-07 | 5.2200e-07 | 1.0069 |
| 200 | 1.1489 | 3869 | 1.586e+09 | 3.074e+07 | 8.58e-04 | 5.8712e-07 | 5.9377e-07 | 1.0113 |
| 233 | 1.2657 | 5273 | 3.669e+08 | 2.340e+07 | 3.72e-04 | 7.0990e-07 | 7.1923e-07 | 1.0131 |
| 240 | 1.3003 | 5478 | 2.763e+08 | 2.101e+07 | 3.12e-04 | 7.3844e-07 | 7.4777e-07 | 1.0126 |
| 250 | 1.3578 | 5700 | 1.859e+08 | 1.756e+07 | 2.42e-04 | 7.7590e-07 | 7.8492e-07 | 1.0116 |
| 260 | 1.4261 | 5834 | 1.260e+08 | 1.431e+07 | 1.86e-04 | 8.0784e-07 | 8.1641e-07 | 1.0106 |
| 267 | 1.4816 | 5873 | 9.626e+07 | 1.224e+07 | 1.54e-04 | 8.2657e-07 | 8.3486e-07 | 1.0100 |
| 300 | 1.8571 | 5423 | 2.661e+07 | 5.311e+06 | 5.65e-05 | 8.8051e-07 | 8.8925e-07 | 1.0099 |

What attenuates the band, MEASURED at 12.0 eV (cross sections in 1e-18 cm^2:
H I 0, He I 0, He 2^3S 1.466, C I 17.62, O I 0, N I 0, because the O I and N I
thresholds are 13.62 and 14.53 eV):

| cell | r/R_p | tau(He 2^3S) | tau(C I) | tau(H I) |
|---|---|---|---|---|
| 100 | 1.0250 | 7.92e-07 | 1.070e-02 | 0 |
| 200 | 1.1489 | 7.87e-07 | 8.570e-04 | 0 |
| 233 | 1.2657 | 7.66e-07 | 3.711e-04 | 0 |
| 300 | 1.8571 | 6.86e-07 | 5.580e-05 | 0 |

So the answer to "what blocks the C I band, since hydrogen does not" is:
nothing does. C I itself carries all of it and reaches 1 percent only at the
deepest cell audited. The unattenuated rate is 9.089e-07 s^-1, of which 56
percent comes from 11.26 to 13.598 eV, and the attenuated rate is 5.0e-07 to
8.8e-07 s^-1 over 1.0 to 1.86 R_p. The rate falls inward not because of
absorption in its own band but because the cells further in see the same weak
field with a slightly larger C I column.

Integrated stellar fluxes at the planet, MEASURED from the grid this run
builds: 45.84 erg cm^-2 s^-1 over the whole 4.77 to 1211 eV span, 0.5656
between 11.26 and 13.598 eV, 19.41 above 13.598 eV.

## 4. Recombination, charge exchange and the balance (question 2)

C II + e -> C I is the Badnell total, radiative plus dielectronic, keyed by the
recombined ion. MEASURED split at the audited cells: the radiative part is
1.39e-12 at 1869 K falling to 6.66e-13 at 5873 K, the dielectronic part
2.00e-13 to 3.07e-13, so the dielectronic channel is 13 to 32 percent of the
total. The fit reproduces the published 1e4 K value 4.7e-13 cm^3 s^-1
(MEASURED from the coefficients in the source: 4.72e-13).

The two carbon charge-exchange rows are radiative charge transfer and are
negligible: `C + H+ -> C+ + H` contributes 1.0e-8 to 5.6e-8 s^-1 against a
photoionization rate 15 to 80 times larger, and the reverse `C+ + H -> C + H+`
carries a 1.7e5 K activation barrier and is dead (1e-44 to 1e-19 s^-1). Their
validity range is the one Huang et al. (2023) Table 4 states for the Stancil
fits; nothing in the audited band is outside it. Collisional ionization of C I
is dead as well (3.9e-32 to 1.8e-11 s^-1).

TABLE B, the hand balance against the code, all terms in s^-1

| cell | T [K] | P(C I) | beta n_e | alpha n_e | k(C+H+) n_H+ | k(C+ + H) n_HI | x(C I) code | x(C I) hand | hand/code |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 1869 | 5.049e-07 | 3.90e-32 | 5.143e-05 | 5.65e-08 | 5.0e-44 | 0.98916 | 0.98920 | 1.0000 |
| 150 | 2355 | 5.184e-07 | 8.93e-26 | 4.540e-05 | 5.37e-08 | 2.5e-36 | 0.98748 | 0.98756 | 1.0001 |
| 200 | 3869 | 5.871e-07 | 3.21e-16 | 3.293e-05 | 4.15e-08 | 1.2e-24 | 0.98107 | 0.98127 | 1.0002 |
| 233 | 5273 | 7.099e-07 | 2.48e-12 | 2.325e-05 | 2.49e-08 | 6.2e-20 | 0.96898 | 0.96936 | 1.0004 |
| 240 | 5478 | 7.384e-07 | 5.77e-12 | 2.072e-05 | 2.12e-08 | 1.7e-19 | 0.96421 | 0.96463 | 1.0004 |
| 250 | 5700 | 7.759e-07 | 1.26e-11 | 1.718e-05 | 1.65e-08 | 4.0e-19 | 0.95543 | 0.95591 | 1.0005 |
| 260 | 5834 | 8.078e-07 | 1.77e-11 | 1.393e-05 | 1.26e-08 | 5.7e-19 | 0.94385 | 0.94440 | 1.0006 |
| 267 | 5873 | 8.266e-07 | 1.76e-11 | 1.191e-05 | 1.03e-08 | 5.3e-19 | 0.93371 | 0.93432 | 1.0007 |
| 300 | 5423 | 8.805e-07 | 1.14e-12 | 5.246e-06 | 3.91e-09 | 1.1e-20 | 0.85452 | 0.85574 | 1.0014 |

Why the neutral fraction is large, in one line: the photons that ionize C I sit
below the H I edge, where this M dwarf puts 0.57 erg cm^-2 s^-1, while the
electrons that recombine it come from hydrogen ionized by the 19.4 erg cm^-2
s^-1 above the edge. The ratio x(C II)/x(C I) = P/(alpha n_e) is then 1e-2 to
1.7e-1 and carbon stays neutral. That is the opposite of an H II region, where
the same 11 to 13.6 eV band is strong and C I is absent; here the band is weak
and the electron density is set by a much harder one.

## 5. The line cooling (question 3)

What `cool_CI_ne_func` returns is `W_FS/n_e` plus a four-term coronal remainder
fit. `W_FS` is the exact statistical equilibrium of the 2p2 3P ground term
(the [C I] 609.1 and 370.4 um lines); MEASURED, it is 8.7e-17 erg cm^-3 s^-1
at cell 233 against a total C I cooling of 4.29e-9, so it is not what carries
the cooling. The remainder is:

```
( 1.68601052e-18 exp(-16400.7/T) + 4.56976084e-18 exp(-22917.8/T)
+ 4.41186120e-17 exp(-60832.4/T) + 2.88720952e-16 exp(-157618.0/T) ) / sqrt(T)
```

MEASURED at 5273 K the four terms are 1.035e-21, 8.13e-22, 5.95e-24 and
2.7e-34 erg cm^3 s^-1, so the two low-excitation terms carry 99.7 percent. The
only C I levels they can be are the even-parity metastables 2p2 1D2
(10192.657 cm^-1 = 14665 K) and 1S0 (21648.030 cm^-1 = 31147 K): the next
level up is 2s2p3 5S* at 33735.121 cm^-1 = 48539 K, and the two high terms
belong to the odd-parity 2s2p3 manifold whose 1657 A multiplet is permitted.
Back-solving the first term for its effective collision strength gives
Upsilon(3P-1D) = 0.34 at 1869 K and 0.65 at 5273 K, which is the right size for
that transition and confirms the identification.

Both metastables decay by forbidden transitions. NIST Atomic Spectra Database
(READ, retrieved 2026-09-13): A(1D2) = 2.2e-4 + 7.3e-5 + 5.9e-8 = 2.93e-4
s^-1 (9850.250, 9824.118, 9808.295 A), A(1S0) = 6.0e-1 + 2.3e-3 + 2.2e-5 =
6.023e-1 s^-1 (8727.131, 4621.569, 4627.344 A). With the back-solved Upsilon
the electron critical density of 1D2 is 2.0e4 cm^-3; 1S0 is at 1.3e7 cm^-3 for
an assumed Upsilon of 0.4. The local n_e is 1.2e7 to 3.3e7 cm^-3.

The bound that settles it needs no collision data at all. In statistical
equilibrium the population of a level rises monotonically with the collision
rate and saturates at its Boltzmann value, so one C I atom cannot radiate more
in these lines than

```
W_LTE(T) = sum_u f_u^Boltz sum_l A_ul h nu_ul ,
```

with f_u over the five levels of the 2p2 configuration. Adding the two
permitted terms at their coronal value, which is correct for them, gives the
ceiling of table C.

TABLE C

| cell | r/R_p | T [K] | n_e | n(C I) | cool(C I) as coded | ceiling | as coded / ceiling | C I share of cooling | cool/heat |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 1.0250 | 1869 | 3.234e+07 | 1.976e+07 | 4.1611e-09 | 2.7690e-12 | 1503 | 93.6% | 0.446 |
| 150 | 1.0615 | 2355 | 3.338e+07 | 4.302e+06 | 5.5141e-09 | 4.0420e-12 | 1364 | 94.1% | 0.562 |
| 200 | 1.1489 | 3869 | 3.074e+07 | 4.321e+05 | 7.8025e-09 | 2.5618e-11 | 305 | 93.6% | 0.649 |
| 233 | 1.2657 | 5273 | 2.340e+07 | 9.875e+04 | 4.2906e-09 | 5.5672e-11 | 77 | 87.0% | 0.652 |
| 240 | 1.3003 | 5478 | 2.101e+07 | 7.399e+04 | 3.2508e-09 | 5.2761e-11 | 62 | 85.0% | 0.631 |
| 250 | 1.3578 | 5700 | 1.756e+07 | 4.935e+04 | 2.0420e-09 | 4.3766e-11 | 47 | 82.1% | 0.591 |
| 260 | 1.4261 | 5834 | 1.431e+07 | 3.305e+04 | 1.1925e-09 | 3.2297e-11 | 37 | 79.0% | 0.539 |
| 267 | 1.4816 | 5873 | 1.224e+07 | 2.497e+04 | 7.8539e-10 | 2.4490e-11 | 32 | 76.7% | 0.497 |
| 300 | 1.8571 | 5423 | 5.311e+06 | 6.316e+03 | 6.8013e-11 | 3.4088e-12 | 20 | 60.3% | 0.266 |

Rates in erg cm^-3 s^-1. The cooling and heating channels are MEASURED from
`output/Cooling_breakdown.txt` and `output/Heating_breakdown.txt` of a
post-processing pass on the same state. The heating there is 88 to 90 percent
He I photoionization, which is what a He/H of 9 makes it.

So the 79 to 94 percent share is not physical: it is one to three decades of
line emission that the atoms cannot produce. Against the same ceiling, C I is
0.06 to 3.0 percent of the cooling, and the metals stop being the dominant
coolant of this band altogether.

The code already knows this failure mode. `Update_EXHALE_stage1` section 48
(READ) fixed it for the ground-term fine structure of C I, C II, N II and O I,
and the Fe II coefficient is a two-dimensional statistical-equilibrium table
for the same reason. What was left coronal is the metastable manifold above
the ground term, and for C I that is where the cooling is.

## 6. The fix, and what it moves (question 4)

`src/modules/radiation/Cool_coeff.f90`, `cool_CI_ne_func`. The remainder fit is
split into the metastable pair and the permitted pair, and the metastable pair
is saturated by the exact two-level solution

```
Lambda_meta,eff = Lambda_meta / (1 + n_e Lambda_meta / W_LTE)
```

which is the coronal rate as n_e -> 0 and W_LTE/n_e as n_e -> infinity, with the
crossover at the critical density because Lambda/W_LTE = q_ul/A identically.
`W_LTE` is the new `lte_emission_CI_metastable`, built from the NIST level
energies and transition probabilities listed above, written at the code site.
H-atom de-excitation of the two metastables is not carried, which can only
leave the result above the true one. The permitted pair is untouched.

TABLE D, the same state before and after (post-processing pass, same input)

| cell | r/R_p | cool(C I) before | after | factor | cool total before | after | cool/heat before | after |
|---|---|---|---|---|---|---|---|---|
| 100 | 1.0250 | 4.1611e-09 | 2.7844e-12 | 1494 | 4.4460e-09 | 2.8767e-10 | 0.446 | 0.029 |
| 150 | 1.0615 | 5.5141e-09 | 4.0428e-12 | 1364 | 5.8611e-09 | 3.5101e-10 | 0.562 | 0.034 |
| 200 | 1.1489 | 7.8025e-09 | 2.5544e-11 | 305 | 8.3319e-09 | 5.5499e-10 | 0.649 | 0.043 |
| 233 | 1.2657 | 4.2906e-09 | 5.5265e-11 | 78 | 4.9336e-09 | 6.9824e-10 | 0.652 | 0.092 |
| 240 | 1.3003 | 3.2508e-09 | 5.2300e-11 | 62 | 3.8230e-09 | 6.2455e-10 | 0.631 | 0.103 |
| 250 | 1.3578 | 2.0420e-09 | 4.3268e-11 | 47 | 2.4863e-09 | 4.8765e-10 | 0.591 | 0.116 |
| 260 | 1.4261 | 1.1925e-09 | 3.1812e-11 | 37 | 1.5088e-09 | 3.4808e-10 | 0.539 | 0.124 |
| 267 | 1.4816 | 7.8539e-10 | 2.4042e-11 | 33 | 1.0239e-09 | 2.6259e-10 | 0.497 | 0.128 |
| 300 | 1.8571 | 6.8013e-11 | 3.2698e-12 | 21 | 1.1271e-10 | 4.7969e-11 | 0.266 | 0.113 |
| 350 | 3.0508 | 2.2985e-13 | 1.2586e-14 | 18 | 2.1709e-12 | 1.9536e-12 | 0.068 | 0.061 |
| 400 | 5.9043 | 1.7384e-17 | 1.4947e-18 | 12 | 3.5529e-14 | 3.5513e-14 | 0.019 | 0.019 |

### Movement on the two regression cases

MEASURED. Both cases were run on scratch copies with two private builds of this
tree, one with the change and one without, `OMP_NUM_THREADS=8` for both. The
control reproduces the golden to 6.4e-9 on `wasp_full` and 3.4e-8 on
`mol_metals` (the residual is the reduction order of eight threads against the
golden's one), so it is a valid control. No golden was refreshed.

| case | quantity | control | with the saturation | shift |
|---|---|---|---|---|
| `wasp_full` | stop | count 17635, du 9.9551e-04 | count 17671, du 9.8531e-04 | same du threshold |
| | log10 Mdot [g/s] | 13.35 | 13.36 | |
| | mass flux rho v r^2, median over r > 1.2 R_p | | | +1.2% (range +1.1 to +1.6%) |
| | cooling at r = 1.0002 R_p | 2.949e-06 | 1.601e-06 | -46% |
| | T at r = 1.0111 R_p [K] | 2913.4 | 2954.0 | +1.4% |
| | rho, largest move (r = 1.2823 R_p) | 9.280e+09 | 9.648e+09 | +4.0% |
| | max abs relative move, `Ion_species.txt` | | | 1.5e-01 |
| `mol_metals` | log10 Mdot [g/s] | 10.58 | 10.58 | |
| | mass flux, median over r > 2 R_p | | | +0.17% |
| | cooling at r = 1.0059 R_p | 1.770e-08 | 4.763e-09 | -73% |
| | T, largest move (r = 4.2103 R_p) [K] | 1357.8 | 1372.8 | +1.1% |
| | max abs relative move, `Ion_species.txt` | | | 2.4e-02 |

Both cases stop on a `du` threshold, whose path dependence the regression notes
put at several percent in the mass-loss rate, so the mass-flux shifts above are
not separable from that spread; the profile shifts are. What is unambiguous is
the cooling: the C I channel falls by half at the WASP-121 b base and by
three quarters at the hot-Uranus base, and both goldens move well past the 1e-3
relative tolerance. They are left as they are.

### Suites

`src/tests/physics_probe` is the suite whose drivers link `Cool_coeff.f90`. Run
with `EXHALE_OBJDIR`/`EXHALE_TEST_OUT` pointing at the private build: 1481
assertions PASS, 0 FAIL. Two drivers abort before their assertions,
`heating_channel_closure` and `co_helium_ion_sink`, both with
`keq_H_H_to_H2 called before h2_thermochemistry_init` raised inside
`mol_rates.f90`. That is a molecular initialization that the driver never
performs; it cannot be reached from a coefficient in the C I cooling, and the
molecular initialization is another item's file. Reported, not touched.

## 7. Noticed, not fixed

- **The same defect is in five more coolants.** The metastable upper levels of
  the C/N/O fits carry the same coronal form. NIST A values (READ, retrieved
  2026-09-13) and the critical densities they imply for an assumed collision
  strength of order unity: N I 2D* A = 2.79e-5 s^-1, n_cr ~ 8e3 cm^-3;
  O II 2D* A = 2.08e-4, n_cr ~ 1e4; N II 1D A = 3.90e-3, n_cr ~ 4e5;
  O I 1D A = 7.48e-3, n_cr ~ 1e6. All are below the n_e of a wind base. N I
  carries 2.5 to 3.0 percent of the cooling of this rung at 1.15 to 1.5 R_p
  and is the next one to matter; the C II 4P metastable, by contrast, decays
  at 50 to 130 s^-1 and stays coronal at these densities. Fixing them needs
  the same treatment and the same split of each fit, which is a separate item.
- **The fits themselves were not re-derived.** `XUVTOP` is not set on this
  machine and the CHIANTI database is not installed, so
  `cooling_data/chianti_cooling.py` and `fit_fs_saturation.py` cannot be run
  and the coronal curves could not be checked against CHIANTI at the audited
  temperatures. What is checked here is the ceiling, which is independent of
  the collision data.
- **The regression case directories use "arm" as a noun** (`armA_LW`,
  `armD_D2`, `armHeH_*`, `arm_heh1_x2matched` under `backup/regression/`).
  Renaming them would break `run_check.sh` and the goldens beside them, so
  they are left alone and reported.
  **DONE 2026-09-16** (PLAN_20260916_rev3 section 10): all sixteen were
  renamed; the old-to-new mapping is `docs/named_case_audit.md` section 6.

## 8. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
D=$EX/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k00

# the cooling and heating channels of the audited state
mkdir -p /tmp/L6/pp/output
cp $D/input.inp $D/base.inp $D/lower_atmosphere_profile.dat /tmp/L6/pp/
cp $D/output/Hydro_ioniz.txt /tmp/L6/pp/output/Hydro_ioniz_IC.txt
cp $D/output/Ion_species.txt /tmp/L6/pp/output/Ion_species_IC.txt
sed -i "s|^Spectrum file: .*|Spectrum file: $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt|" \
    /tmp/L6/pp/input.inp
sed -i -e '/^Restart intent:/d' -e '/^Solver:/d' -e 's/^Do only PP: .*/Do only PP: True/' \
    /tmp/L6/pp/input.inp
cd /tmp/L6/pp && OMP_NUM_THREADS=8 $EX/EXHALE.x
#   -> output/Cooling_breakdown.txt, output/Heating_breakdown.txt

# the private build with the saturation, and the same pass again
cd $EX && make OBJDIR=build_L6 EXE=EXHALE_L6.x
```

The independent photoionization and balance calculation is a Python
transcription of the published cross-section fits and of `sed_read.f90`'s
binning; it needs only `output/Ion_species.txt`, `output/Hydro_ioniz.txt` and
the SED file, and is reproduced by the tables above rather than kept as a
script in the tree.
