# L6b: the metal line cooling at the local electron density

Item L6b of `docs/PLAN_20260913_lhs_stationary.md`, opened by section 7 of
`docs/lhs1140b_stationary_L6_20260913.md`: L6 saturated the C I metastables by
hand and listed five more coolants carrying the same defect. This item audits
every metal and helium cooling term of `src/modules/radiation/Cool_coeff.f90`
against its own critical density and replaces the coronal treatment of the six
C/N/O coolants with the CHIANTI statistical equilibrium at the cell's own
electron density, following the user's instruction (2026-09-13) to adopt the
method of `~/MoCHII/MoCHII_v1.00` rather than bolt a two-level saturation onto
each metastable term. Every number is MEASURED on this tree unless marked READ.

## 1. Verdict

**Six coolants were evaluated above the critical density of the levels they
excite and are now density-resolved: C I, C II, N I, N II, O I and O II.** The
coefficient the cooling assembly takes for them stood at 0.02 to 5.4e4 times
the CHIANTI statistical equilibrium of the same ion over the wind-relevant
(T, n_e) domain; it now stands at 0.969 to 1.009, every ion, every point of a
7 x 6 grid spanning 1000 to 20000 K and 1e4 to 1e9 cm^-3 (section 6).

**What it does to the fixtures.** The LHS 1140 b rung loses 17 to 21 percent
of its cooling over 1.15 to 1.50 R_p, the band whose energy row refuses the
L4b rung, and `cool/heat` there falls from 0.043-0.128 to 0.036-0.105.
`wasp_full` moves a great deal: the cooling below 1.1 R_p falls by up to 72
percent, the temperature above 1.2 R_p rises 7 to 14 percent, and the
mass-loss rate is unchanged at log10 Mdot = 13.36. `mol_metals` barely moves
(T to 8e-4, Mdot unchanged), because at 890 to 2650 K the C/N/O channels are in
the exponential tail and the layer is cooled by H3+ and the ground-term fine
structure. Both goldens move; neither was refreshed (section 10).

**No other cooling term in the file needs it.** Mg I, Mg II, Ca II and Na I are
single permitted resonance lines with critical densities 3.5e13 to 1.1e16
cm^-3; Fe I is built from permitted multiplets and has no CHIANTI model atom;
Fe II already carries its own two-dimensional statistical-equilibrium table;
and the one helium channel that is coronal above its critical density, the
19.82 eV excitation into He 2^3S, is below 2e-7 of the cooling in every fixture
because of its exponential (sections 3 and 7).

**A defect in the fits themselves was found on the way, and the change removes
it from the production path.** The closed-form C/N/O fits were built from the
CHIANTI `.scups` THEORETICAL transition energies, while the level populations
that decide the emission use the observed level energies. The two differ by up
to 14 percent in dE for the O II 4S* - 2D* excitation, which at 2000 K is a
factor 14 in exp(-dE/kT). MEASURED, the table's own coronal column divided by
the fit is 220 at 1000 K and 13.8 at 2000 K for O II, 4.2 and 1.9 for N II,
2.3 and 1.5 for N I, 1.6 and 1.2 for C I, and 0.92 to 1.01 for C II and O I
across 1000 to 20000 K (section 6). The fits are kept as the coronal reference
and for the legacy `cno_cool 0` branch; the assembly no longer uses them for
these six ions.

**The two-level saturation this item first implemented, in the format L6 used
for C I, is gone.** Its measurement is kept in section 8 because it is the
evidence that decided the design: at best it left the coefficient within a
factor 8 of CHIANTI, and it could not do better, because it pairs a ceiling
with an excitation rate whose exponent is a fit parameter and not a level
energy. The C I saturation L6 landed is replaced by the same table as the rest.

## 2. Design: which route, and the form

The instruction is to adopt the method of `~/MoCHII/MoCHII_v1.00`, the same
author's photoionization code: `md/COOLING_LOCAL_NE_PLAN.md`, implemented and
the default there since 2026-07-26 (`cooling_model='local_ne'`). The line
cooling of every ion covered by a level model is taken from the statistical
equilibrium of that model at the cell's own `(T, n_e)`. MoCHII reports that at
fixed state this moves its total cooling from 1.116 times Cloudy c25.00 to
1.011 (READ, that document's section 1).

### Route (a), an offline table, is the one that fits this code

EXHALE already carries exactly this structure for Fe II: `cool_logL_FeII_ne` is
a two-dimensional `(log T, log n_e)` array of `parameter` data built offline by
`cooling_data/fe2_cooling.py`, and `cool_table_value_2d` interpolates it. The
axes exist (`cool_logT`, 41 points over log T 3.0 to 5.0; `cool_logne`, 29
points over log n_e 0 to 14), and the repository rule that atomic data are
fitted offline and only evaluated online is the rule this follows. This is also
MoCHII's own design point D1, for the same reason: a bilinear interpolation is
ten flops, and for a given ion the level populations normalized to the ion
density depend only on `(T_e, n_e)`, so a table in those two variables carries
the full physics and the only error is interpolation error, which section 6
measures.

Route (b), porting `nlevel_mod.f90` for a runtime statistical-equilibrium
solve, buys MoCHII a property EXHALE does not need: there the cooling and the
emitted line list must come from one solver on one level model. EXHALE's
transit forward model reads level populations only for the O I ground term and
the H n=2 system, both of which have their own exact solutions already. Against
that, a solve of a 40 to 200 level atom for each of six ions inside the
ionization sweep, which runs for every cell at every step and inside the Newton
residual, is a cost with no return here. Route (a) is chosen.

### What is tabulated, and what is not

For an ion whose ground term is split the coefficient is

```
Lambda(T, n_e) = W_FS(T, n_e, n_HI, beta)/n_e  +  Lambda_rem(T, n_e)
```

with `W_FS` the code's existing EXACT statistical equilibrium of the split
ground term and

```
Lambda_rem(T, n_e) = ( sum_{u outside the ground term, l} f_u A_ul E_ul )/n_e
```

the tabulated quantity, `f_u` the level populations of the full CHIANTI model.
The ground term is excluded from the table because the code solves it exactly
itself, and solves it BETTER than any such table could: it carries the
neutral-hydrogen de-excitation channel of the fine-structure levels, and the
escape probability of those lines enters inside the solution rather than as a
factor on it. Neither is available from CHIANTI. Including the ground term in
the table would both double count it and throw those two away.

MoCHII keeps its closed-form fits as the carrier and tabulates only the
suppression ratio, so that an ion's truncated level model cannot bias the
total (its design point D2: Fe III is 14 of 34 levels there). That construction
was implemented here first and then rejected on measurement, for two reasons
specific to this code:

- the CHIANTI model used here is the FULL one, and it is the same data the
  fits were made from, so there is nothing for the fit to cover that the table
  does not. MEASURED, the table's coronal column is at or above the fit for
  C I, N I, N II and O II at every temperature, and below it by at most 8
  percent for C II and O I;
- where it is below by those few percent, the reason is the ground-term
  weighting convention of the fit and not a missing transition, and leaving
  that few percent unsuppressed is ruinous once the rest is suppressed by three
  decades: MEASURED for O I at 5000 K and n_e = 1e9, the D2 form returns 11
  times the CHIANTI value, all of it the unsuppressed 4 percent.

So the table carries the cooling outright. The price is that the low-density
limit is no longer the fit; it is the CHIANTI coronal limit of the same
channels, and section 6 reports the two side by side. That difference is the
`.scups` energy defect of section 1, so keeping the old limit would be keeping
an error.

### Scope

Six ions: C I, C II, N I, N II, O I, O II. Mg I, Mg II, Ca II and Na I are
single permitted resonance lines whose critical densities are decades above
anything here (table 1); Fe I has no CHIANTI model atom in v11 and its upper
levels are permitted; Fe II already has its own two-dimensional table. C III,
N III, O III, Si, K and S carry no cooling coefficient at all in
`cool_coeff_of_ion`.

### Provenance

Method: `~/MoCHII/MoCHII_v1.00/md/COOLING_LOCAL_NE_PLAN.md`, design points D1
to D4, and the tools that built MoCHII's version of it,
`tools/fitting/{chianti_cooling.py, fit_nlevel.py, verify_nlevel_pyneb.py}`
and `src/nlevel_mod.f90` (the same author's code). Atomic data: CHIANTI
v11.0.2 at `/nfs/mocafe/kiseon/RT_Codes/CHIANTI/dbase` (VERSION 11.0.2, READ)
through ChiantiPy 0.15.2. Both are recorded in the header of
`cooling_data/metal_cooling_density_resolved.py` and at the table in
`Cool_coeff.f90`.

## 3. The audit

Critical densities are MEASURED as n_cr = A_tot/q_ul with
q_ul = 8.629e-6 Upsilon/(g_u sqrt(T)). Level energies and transition
probabilities: NIST Atomic Spectra Database (READ, retrieved 2026-09-13).

Electron densities of the three fixtures, MEASURED from
`output/Cooling_breakdown.txt`: WASP-121 b (`wasp_full`) 7.7e8 to 4.0e9 cm^-3
over T = 2366 to 11197 K; the hot Uranus (`mol_metals`) 5.1e6 to 7.5e7 over
T = 893 to 2645 K; the LHS 1140 b rung 5.3e6 to 3.3e7 over T = 506 to 5875 K.

TABLE 1, the metal and helium cooling terms

| term | routine | upper level of the channel | A_tot [s^-1] | n_cr [cm^-3] | above n_cr in a fixture? | action |
|---|---|---|---|---|---|---|
| C I above the ground term | `cool_CI_ne_func` | 2p2 1D2, 1S0; 2s2p3 5S*; 2s2p 3s/3p | 2.9e-4 to ~1e8 | 2.4e4 to ~1e16 | yes, all three | **tabulated** |
| C II above the ground term | `cool_CII_ne_func` | 2s2p2 4P; 2s2p2 2D and above | 9.9 to ~1e8 | 6.6e8 to ~1e15 | yes, `wasp_full` | **tabulated** |
| N I above the ground level | `cool_NI_ne_func` (new) | 2p3 2D*, 2P*; 3s 4P and above | 7.6e-6 to ~4e8 | 1.3e4 to ~1e16 | yes, all three | **tabulated** |
| N II above the ground term | `cool_NII_ne_func` | 2p2 1D2, 1S0; 2s2p3 5S*, 3D* | 3.9e-3 to ~3e8 | 1.0e5 to ~1e15 | yes, all three | **tabulated** |
| O I above the ground term | `cool_OI_ne_func` | 2p4 1D2, 1S0; 3s 5S*, 3S* | 7.5e-3 to 5.6e8 | 2.6e6 to ~1e16 | yes, all three | **tabulated** |
| O II above the ground level | `cool_OII_ne_func` (new) | 2p3 2D*, 2P*; 2s2p4 4P and above | 3.1e-5 to ~1e8 | 4.2e4 to ~1e15 | yes, all three | **tabulated** |
| C I, C II, N II, O I ground-term fine structure | same routines | 2p 2P, 2p2 3P, 2p4 3P | 8e-8 to 9e-5 | 1e0 to 1e5 | yes | already exact statistical equilibrium |
| Mg I 2852 A | `cool_MgI_func` | 3s3p 1P*1 | 4.91e8 | 8.8e15 to 1.1e16 | no | coronal correct |
| Mg II h and k | `cool_MgII_func` | 3p 2P* | 2.60e8 | 3.3e14 to 7.0e14 | no | coronal correct |
| Ca II H and K | `cool_CaII_func` | 4p 2P* | 1.47e8 | 1.9e14 to 3.6e14 | no | coronal correct |
| Na I D | `cool_NaI_func` | 3p 2P* | 6.16e7 | 3.5e13 to 8.3e13 | no | coronal correct |
| Fe I | `cool_logL_FeI` table | permitted E1 multiplets | ~1e7 to 1e8 | >= 1e13 | no | coronal correct for the upper levels, section 7 |
| Fe II | `cool_FeII_ne` | multilevel | - | - | - | already a 2-D (T, n_e) table |
| He I 19.82 eV | `lambda_coex_HeI` | He 2^3S | 1.272e-4 | 1.0e5 to 1.6e5 | yes, but negligible | left coronal, section 7 |
| He 2^3S -> 2^3P 10830 A | `coex_rate_HeI23S_10830` | He 2^3P | 1.022e7 | 2.5e13 to 3.7e13 | no | coronal correct |
| H I Ly-alpha, He II Ly-alpha | `lambda_coex_HI`, `lambda_coex_HeII` | 2p | 6.3e8, 1.0e10 | >= 1e17 | no | coronal correct |

Si, K and S carry no line cooling at all in `cool_coeff_of_ion`; that is an
omission, not a density error, and is out of scope here.

## 4. The reproduction (RED, then GREEN)

RED and GREEN are the same quantity: the coefficient the assembly takes,
against the CHIANTI statistical equilibrium of the same ion at the same
(T, n_e). The full table is section 6; the summary over a 7 x 6 grid covering
1000 to 20000 K and 1e4 to 1e9 cm^-3, with the Fortran evaluated at n_HI = 0
for the like-for-like:

| ion | RED, coronal fit / CHIANTI | GREEN, as coded now / CHIANTI |
|---|---|---|
| C I | 1.00 to 5.4e+04 | 0.992 to 1.009 |
| C II | 0.98 to 1.82 | 0.969 to 1.000 |
| N I | 0.99 to 5.2e+04 | 0.979 to 1.007 |
| N II | 0.77 to 1.2e+04 | 0.988 to 1.001 |
| O I | 0.99 to 354 | 0.986 to 1.005 |
| O II | 0.02 to 1.9e+04 | 0.979 to 1.002 |

A second, collision-data-free reproduction of the same failure, MEASURED on the
LHS 1140 b rung state
`LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k00/output/` by a
post-processing pass: the emission of one ion divided by the LTE ceiling of
every level of its ground configuration plus the coronal value of its permitted
terms. Statistical equilibrium is monotone in the collision rate and saturates
at the Boltzmann population, so no value above 1 is physically possible. This
is the delivered tree, BEFORE this item (it already carries L6's C I fix):

| r/Rp | T [K] | n_e | C I | C II | N I | N II | O I | O II |
|---|---|---|---|---|---|---|---|---|
| 1.0250 | 1868.5 | 3.234e+07 | 0.999 | 1.00 | 1.55e+03 | 153 | 1.62 | 561 |
| 1.1489 | 3868.8 | 3.074e+07 | 0.997 | 0.019 | 101 | 341 | 12.0 | 260 |
| 1.2657 | 5273.2 | 2.340e+07 | 0.993 | 0.0122 | 37.1 | 174 | 8.79 | 73.3 |
| 1.3578 | 5699.8 | 1.756e+07 | 0.989 | 0.0088 | 24.2 | 111 | 6.22 | 44.0 |
| 1.4816 | 5872.6 | 1.224e+07 | 0.982 | 0.0060 | 16.0 | 72.6 | 4.21 | 28.3 |
| 1.8571 | 5423.2 | 5.311e+06 | 0.959 | 0.0027 | 7.98 | 37.3 | 1.96 | 15.3 |

## 5. The change

`src/modules/radiation/Cool_coeff.f90`

- The tabulated `log10 Lambda_rem(T, n_e)` for the six ions, on the existing
  `cool_logT` and `cool_logne` axes (41 x 29 each), generated by
  `cooling_data/metal_cooling_density_resolved.py`.
- `remainder_table_value`, bilinear interpolation in
  `(log10 T, log10 n_e)` with edge clamping, and
  `metal_cooling_above_ground_term(i, T, n_e, lam_coronal)`, which selects the
  table for an ion by its canonical `im_*` index and returns `lam_coronal`
  unchanged for an ion with no table. The case list lives there and nowhere
  else.
- `cool_CI_ne_func`, `cool_CII_ne_func`, `cool_NII_ne_func`,
  `cool_OI_ne_func`: the ground term is solved as before; the coronal
  remainder is now the argument `lam_coronal` and the value used is the table's.
- `cool_NI_ne_func`, `cool_OII_ne_func` and their grid wrappers, new. N I and
  O II have a single-level 4S* ground term, so there is no fine structure to
  solve and the whole coefficient is the table.
- Removed: the two-level saturation this item first implemented and L6's C I
  version of it (`metastable_saturated_coefficient`, the
  `lte_emission_*_metastable` functions, and the NIST level and decay
  constants that only served them). The levels they named are named in the
  comment on the table instead.

`src/modules/radiation/util_ion_eq.f90`, four lines: the `cno_chianti`
override block in `eval_cool` now also takes N I and O II from the
density-dependent coefficients. Without this the two ions keep the pure coronal
fit, since they have no fine-structure floor and so were never in that block.
This is the one file outside `Cool_coeff.f90` the physics change touched.

`src/tests/physics_probe/heating_channel_closure.f90`,
`src/tests/physics_probe/co_helium_ion_sink.f90`: `call
h2_thermochemistry_init` before the code that reaches `keq_H_H_to_H2`, with
`use mol_rates`. See section 9.

`cooling_data/metal_cooling_density_resolved.py`, new: the table builder.
Its `statistical_equilibrium_cooling` is the one place the CHIANTI sum is
formed, so the table and every check of it use the same function; passing
`n_ground_term = 0` gives the whole ion, which is what section 6 compares
against.

## 6. Against CHIANTI v11.0.2

The comparison quantity is the one the assembly uses,

```
Lambda(T, n_e) = ( sum_{u>l} f_u(T,n_e) A_ul E_ul ) / n_e   [erg cm^3 s^-1],
```

with `f_u` from ChiantiPy's `populate`, the full multilevel statistical
equilibrium of the ion under electron and proton collisions and every radiative
decay of the `.wgfa` file. `ion.intensity` is an emission measure and is not
used. The Fortran side is a Python transcription of the six coefficients which
reproduces the binary to 5.6e-7 relative or better at thirty cells of the LHS
rung (MEASURED against `output/Cooling_breakdown.txt` divided by n_e n_ion), so
it is the code's coefficient and not a second model of it.

CHIANTI has no neutral-hydrogen collisions, so the grid evaluates the Fortran
at n_HI = 0. The size of that omission is MEASURED separately: the ratio of the
coefficient with the run's own n_HI to the coefficient at n_HI = 0 is 1.000 to
1.005 at every cell of all three fixtures. That is not an accident of small
n_HI (it reaches 2.8e12 cm^-3 at the WASP-121 b base): the fine-structure
levels the H channel acts on are saturated into LTE there, and in that limit
the collider identity drops out of the populations.

TABLE 2, Lambda [erg cm^3 s^-1], n_HI = 0

| ion | T [K] | n_e | CHIANTI | coronal fit | fit/CHIANTI | as coded now | now/CHIANTI |
|---|---|---|---|---|---|---|---|
| C I | 2000 | 1e+05 | 1.6205e-24 | 1.1438e-23 | 7.06 | 1.6190e-24 | 0.999 |
| C I | 2000 | 1e+09 | 2.1142e-28 | 1.1429e-23 | 5.41e+04 | 2.1147e-28 | 1.000 |
| C I | 5000 | 1e+07 | 1.6457e-23 | 1.5607e-21 | 94.8 | 1.6450e-23 | 1.000 |
| C I | 12000 | 1e+09 | 1.0675e-21 | 1.2639e-20 | 11.8 | 1.0625e-21 | 0.995 |
| C II | 5000 | 1e+09 | 1.1102e-24 | 1.8103e-24 | 1.63 | 1.1083e-24 | 0.998 |
| C II | 12000 | 1e+09 | 1.3214e-21 | 1.8110e-21 | 1.37 | 1.3131e-21 | 0.994 |
| N I | 2000 | 1e+09 | 4.0005e-31 | 1.5309e-26 | 3.83e+04 | 3.9952e-31 | 0.999 |
| N I | 5000 | 1e+07 | 4.3490e-24 | 1.2980e-22 | 29.8 | 4.3440e-24 | 0.999 |
| N I | 12000 | 1e+09 | 6.2626e-23 | 5.1917e-21 | 82.9 | 6.2103e-23 | 0.992 |
| N II | 2000 | 1e+09 | 1.8694e-28 | 1.3084e-24 | 6999 | 1.8678e-28 | 0.999 |
| N II | 5000 | 1e+07 | 9.8940e-24 | 9.8174e-22 | 99.2 | 9.8922e-24 | 1.000 |
| N II | 12000 | 1e+09 | 3.1175e-22 | 1.0829e-20 | 34.7 | 3.0968e-22 | 0.993 |
| O I | 2000 | 1e+09 | 1.0397e-27 | 4.6604e-26 | 44.8 | 1.0394e-27 | 1.000 |
| O I | 5000 | 1e+09 | 1.8865e-25 | 6.6761e-23 | 354 | 1.8860e-25 | 1.000 |
| O I | 12000 | 1e+07 | 2.5462e-22 | 1.4866e-21 | 5.84 | 2.5435e-22 | 0.999 |
| O II | 2000 | 1e+05 | 4.8638e-29 | 1.0344e-28 | 2.13 | 4.8518e-29 | 0.998 |
| O II | 5000 | 1e+09 | 1.1113e-26 | 3.7475e-23 | 3372 | 1.1097e-26 | 0.999 |
| O II | 12000 | 1e+07 | 5.1823e-22 | 4.5903e-21 | 8.86 | 5.1378e-22 | 0.991 |

### The interpolation error

MEASURED, the tabulated `Lambda_rem` interpolated off the grid nodes against a
direct CHIANTI solve of the same sum, 80 random `(T, n_e)` points for each ion
with n_e from 1e2 to 1e10:

| ion | 2000 K to 1e5 K, median | worst | 1000 to 2000 K, median | worst | largest Lambda_rem there |
|---|---|---|---|---|---|
| C I | 3.5e-03 | 2.6e-02 | 1.3e-02 | 5.6e-02 | 1.0e-23 |
| C II | 2.8e-03 | 4.7e-02 | 4.3e-01 | 5.8e+09 | 4.9e-33 |
| N I | 5.9e-03 | 4.3e-02 | 2.4e-02 | 6.6e-02 | 1.4e-26 |
| N II | 2.8e-03 | 4.5e-02 | 2.1e-02 | 6.7e-02 | 1.9e-24 |
| O I | 4.6e-03 | 3.3e-02 | 2.4e-02 | 6.0e-02 | 3.1e-26 |
| O II | 7.9e-03 | 4.6e-02 | 3.4e-02 | 7.7e-02 | 8.4e-28 |

The C II entry below 2000 K is the one place the log-log interpolation fails,
and it fails on a number that is zero for every purpose: the C II remainder has
its lowest excitation at 61715 K, so below 2000 K it falls 2.9 decades across
one grid cell in log T, and its largest value anywhere in that range is
4.9e-33 erg cm^3 s^-1 against a C II coefficient of ~1e-25 set entirely by the
[C II] 158 um fine structure, which is solved exactly. The T axis is the one
`cool_logT` already declares; refining it for that corner would be refining a
zero.

### The table's coronal column against the fit it replaces

MEASURED, the ratio of the table's floor column (n_e = 1 cm^-3) to the
closed-form coronal remainder:

| ion | 1000 K | 2000 | 3000 | 5000 | 8000 | 12000 | 20000 |
|---|---|---|---|---|---|---|---|
| C I | 1.555 | 1.218 | 1.122 | 1.042 | 1.016 | 0.988 | 0.965 |
| C II | 0.939 | 0.960 | 0.922 | 1.002 | 1.012 | 0.993 | 1.005 |
| N I | 2.332 | 1.500 | 1.275 | 1.161 | 1.090 | 1.046 | 1.026 |
| N II | 4.153 | 1.941 | 1.511 | 1.277 | 1.143 | 1.062 | 1.017 |
| O I | 0.924 | 0.955 | 0.959 | 0.972 | 0.996 | 0.992 | 0.989 |
| O II | 220.4 | 13.81 | 5.396 | 2.620 | 1.712 | 1.350 | 1.138 |

These are not interpolation or method differences: they are the fit's own
error, and its cause is identified. The fits were built from the `.scups`
theoretical transition energies. MEASURED against the observed level energies
of the same database, the O II 4S* - 2D* excitation carries dE = 30660.6 cm^-1
in the `.scups` file against an observed 26810.8, a ratio 1.1436; N I
4S* - 2D* is 1.0315, N I 4S* - 2P* 1.0525, C I 3P - 1D 1.0333. At 2000 K the
O II mismatch alone is exp(0.1436 x 38574/2000) = 16. The fit exponents follow:
the O II fit's low term sits at 44128.6 K against a 2D* level at 38574 K, a
ratio 1.144. The tabulated populations use the observed energies, so the table
does not carry the defect. The fits are left as they are, as the coronal
reference and for the legacy branch, and the defect is reported.

A separate check, MEASURED with the ground-term-Boltzmann convention the fits
were made with: the fits reproduce the CHIANTI v11.0.2 coronal curve OF THIS
INSTALLATION, computed the same way from the same `.scups` energies, to within
2.4 percent (C I), 4.4 (C II), 0.6 (N I), 0.6 (N II), 2.3 (O I) and 0.2
(O II) over 1e3 to 1e5 K. So the fits are faithful to what they were fitted
to; what they were fitted to used the wrong energy in the exponential.

## 7. Audited and left as it is

- **The 19.82 eV helium channel.** `lambda_coex_HeI` excites He 2^3S
  (A(2^3S -> 1^1S) = 1.272e-4 s^-1, n_cr ~ 1e5 cm^-3) and assumes every
  excitation radiates, which is the coronal limit three to four decades above
  that critical density. MEASURED as `1.1e-19 T^0.082 exp(-2.3e5/T)` times
  n_He I(singlet) n_e against the total cooling, it is below 1.8e-7 of the
  cooling in `wasp_full`, 1.8e-33 in `mol_metals` and 2.9e-11 at the LHS rung,
  because 19.82 eV is 2.3e5 K and no fixture exceeds 11200 K. Correcting it is
  not a table: EXHALE already tracks n(2^3S) as a species, so the right
  treatment is an energy ledger against the tracked population, which belongs
  with the helium metastable and not here. Reported.
- **Fe I.** Its table is built from permitted E1 multiplets whose upper levels
  are coronal at every density here, but its LOWER levels are taken Boltzmann
  over the metastable manifold, which is the LTE assumption. That is right at a
  wind base and overestimates the cooling in thin gas, the opposite failure
  mode, and it is the treatment Huang et al. (2023) section 2.5 adopts (READ).
  CHIANTI v11 has no Fe I model atom, so the route taken here is not available
  for it. Left as it is and reported.
- **The legacy branches.** With `cno_cool 0` the AIOLOS forms are selected and
  none of the density corrections apply, including the ground-term ones already
  in the tree. Those branches are kept for comparison and are documented as
  coronal at the code site.
- **Mg I and Ca II omissions.** The Mg I fit carries the 2852.13 A resonance
  line only, not the 4571.10 A intercombination line (3s3p 3P*1,
  A = 2.54e2 s^-1), and the Ca II fit carries H and K only, not the
  7291.47/7323.89 A forbidden decays of the 3d 2D metastable (A = 1.3 s^-1
  each). Both are omissions that UNDERSTATE the cooling, not density errors;
  reported, not changed.
- **Regression case directories using "arm" as a noun** (`armA_LW`, `armD_D2`,
  `armHeH_*`, `arm_heh1_x2matched` under `backup/regression/`). Renaming them
  breaks `run_check.sh` and the goldens beside them. Reported again, unchanged.
  **DONE 2026-09-16** (PLAN_20260916_rev3 section 10): all sixteen were
  renamed; the old-to-new mapping is `docs/named_case_audit.md` section 6.

## 8. The two-level saturation, and why it was replaced

This item first implemented the construction L6 used for C I, for all six
coolants: split each coronal fit into its metastable and permitted terms and
apply the exact two-level saturation

```
Lambda_meta,eff = Lambda_meta/(1 + n_e Lambda_meta/W_LTE)
```

with `W_LTE` the LTE emission of those levels built from NIST level energies
and transition probabilities. It works, in the sense that it removed two to
four decades of the error: MEASURED over the same grid as section 4, it left
the coefficient at 0.86 to 7.58 of CHIANTI for C I, 0.98 to 1.09 for C II,
0.31 to 8.24 for N I, 0.08 to 3.02 for N II, 0.45 to 2.46 for O I and 0.02 to
3.75 for O II.

It could not do better, and the reason is instructive. The fits are EFFECTIVE
multi-exponential fits; nothing constrained their exponents to be level
energies, and section 6 shows they are not. Saturating a term against a
specific level's ceiling therefore pairs a ceiling with an excitation rate
carrying the wrong exponent. Saturating each term against its own level rather
than the pair against a shared ceiling was also MEASURED and was not uniformly
better (root-mean-square |ln(ratio to CHIANTI)| over the grid fell for C I 0.79
to 0.69, N I 0.87 to 0.82, N II 0.66 to 0.60 and O II 0.86 to 0.79, but rose
for O I 0.29 to 0.31). The table has no such coupling: it is the statistical
equilibrium itself.

## 9. The physics_probe drivers

`heating_channel_closure` and `co_helium_ion_sink` aborted with
`(mol_rates) ERROR: keq_H_H_to_H2 called before h2_thermochemistry_init`.
Both call into the molecular rates directly, `heating_of_composition` and
`set_mol_coeffs`, which bypass the serial prologue of `ioniz_eq` where the
initializer is called (`ionization_equilibrium.f90` line 1411, and
`diffusive_photochemistry.f90` line 1188, both guarded by
`h2_thermochemistry_ready`). Each driver now calls `h2_thermochemistry_init`
itself.

A sweep of `src/tests` for the same gap, MEASURED by
`grep -rn "call set_mol_coeffs\|call mol_heh_rows\|call
molecular_chemical_heating\|call heating_of_composition\|keq_H_H"`, finds no
other driver that reaches the molecular rates outside `ioniz_eq`;
`h2_rovibrational_identity` already calls the initializer.

RED, the delivered drivers built from a copy of this tree: two `ERROR STOP 1`,
one per driver. GREEN, after the change and with the table version of the
cooling: `physics_probe: PASSED`, 1491 assertions, 0 failures, exit status 0.

One note on running the suite: with `PATH=/usr/bin:$PATH`, the build recipe of
the shell notes, the suite picks up `/usr/bin/python3` (3.6.9) and
`metal_photoion_table_transcription.py` fails on a `subprocess.run` keyword
added in 3.7. That is the PATH, not the driver; with the default toolchain
(`/opt/miniconda3/bin/python3`, 3.11.16) the suite passes.

## 10. Movement on the fixtures

MEASURED. The control is a build of the CURRENT tree with only this item's
four files reverted to the state it was handed (other workers moved other
files during the item, so a control built at the start would no longer isolate
this change). Both binaries were run on scratch copies with
`OMP_NUM_THREADS=8`. No golden was refreshed.

### The LHS 1140 b rung, post-processing pass

| r/Rp | T [K] | cool before | cool after | after/before | heat | cool/heat before | after |
|---|---|---|---|---|---|---|---|
| 1.0250 | 1868.5 | 2.8767e-10 | 2.8609e-10 | 0.994 | 9.9785e-09 | 0.0288 | 0.0287 |
| 1.1489 | 3868.8 | 5.5499e-10 | 4.5985e-10 | 0.829 | 1.2844e-08 | 0.0432 | 0.0358 |
| 1.2657 | 5273.2 | 6.9824e-10 | 5.5120e-10 | 0.789 | 7.5669e-09 | 0.0923 | 0.0728 |
| 1.3578 | 5699.8 | 4.8765e-10 | 3.9298e-10 | 0.806 | 4.2085e-09 | 0.1159 | 0.0934 |
| 1.4816 | 5872.6 | 2.6259e-10 | 2.1688e-10 | 0.826 | 2.0587e-09 | 0.1276 | 0.1053 |
| 1.8571 | 5423.2 | 4.7969e-11 | 4.2436e-11 | 0.885 | 4.2438e-10 | 0.1130 | 0.1000 |

Rates in erg cm^-3 s^-1. Over 1.15 to 1.50 R_p, the band whose energy row
refuses the L4b rung, the metal channels fall by these factors: N II 89 to 347,
O II 35 to 257, N I 25 to 123, O I 6.2 to 12.3, C I 1.3 to 1.8, C II 1.00 to
1.06. `cool/heat` there falls from 0.043 to 0.128 down to 0.036 to 0.105.

### `wasp_full`, WASP-121 b with He 2^3S and metals

The case stops on the `du` threshold, whose path dependence the regression
notes put at several percent in the mass-loss rate, so the mass-flux number
below is not separable from that spread; the profile numbers are.

| quantity | before | after |
|---|---|---|
| stop | count 17671, du 9.8531e-04 | count 12251, du 9.7613e-04 |
| log10 Mdot [g/s] | 13.36 | 13.36 |
| mass flux rho v r^2, median over r > 1.2 R_p | | -0.49% (range -5.5% to -0.12%) |
| T at r = 1.0111 R_p [K] | 2954.0 | 2758.8 (-6.6%) |
| T at r = 1.2000 R_p [K] | 8210.4 | 8804.9 (+7.2%) |
| T at r = 1.5006 R_p [K] | 10871.5 | 12411.0 (+14.2%) |
| cool at r = 1.0500 R_p | 2.8144e-06 | 7.7368e-07 (-72.5%) |
| cool at r = 1.2000 R_p | 3.6707e-05 | 4.3995e-05 (+19.9%) |
| rho, largest move (r = 1.1147 R_p) | 1.0374e+11 | 7.4224e+10 (-28.5%) |
| max abs relative move, `Ion_species.txt` | | 1.83 |

The cooling falls by up to 72 percent below 1.1 R_p and the gas there runs
cooler, while above 1.2 R_p the temperature rises 7 to 14 percent and the
cooling there RISES 20 percent: with less cooling at the base the wind arrives
hotter and denser higher up, and the metal channels there respond. Both
goldens move far past the 1e-3 relative tolerance. They are left as they are.

### `mol_metals`, the hot-Uranus gate with solar trace metals

| quantity | before | after |
|---|---|---|
| log10 Mdot [g/s] | 10.58 | 10.58 |
| mass flux, median over r > 2 R_p | | +0.007% |
| largest T move (r = 4.2103 R_p) [K] | 1372.8 | 1373.9 (+0.08%) |
| largest cool move (r = 1.9476 R_p) | 1.0034e-09 | 9.9348e-10 (-0.99%) |
| max abs relative move, `Ion_species.txt`, species above 1e-10 of their peak | | 1.8e-03 |

This gate barely moves, and the reason is physical: at 890 to 2650 K the six
C/N/O coefficients are down in the exponential tail, and what cools that layer
is the H3+ infrared emission and the ground-term fine structure, neither of
which this item touches. The metal channels themselves do move, by the factors
N I 53 to 5.3e3, O II 0.33 to 1.4e3, N II 1.0 to 639, O I 1.0 to 12.4 and C I
0.97 to 1.29; they are simply not what sets the temperature there. That golden
moves by less than the 1e-3 relative tolerance in the hydrodynamic columns and
more than it in `Ion_species.txt`. It is left as it is.

## 11. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00

# the private build
cd $EX && make OBJDIR=build_L6b EXE=EXHALE_L6b.x

# the table
XUVTOP=/nfs/mocafe/kiseon/RT_Codes/CHIANTI/dbase OMP_NUM_THREADS=4 \
    python3 $EX/cooling_data/metal_cooling_density_resolved.py out.txt
#   -> paste out.txt into the table block of Cool_coeff.f90

# the LHS rung post-processing pass (as in the L6 memo section 8)
D=$EX/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k00
mkdir -p /tmp/L6b/pp/output
\cp -f $D/input.inp $D/base.inp $D/lower_atmosphere_profile.dat /tmp/L6b/pp/
\cp -f $D/output/Hydro_ioniz.txt /tmp/L6b/pp/output/Hydro_ioniz_IC.txt
\cp -f $D/output/Ion_species.txt /tmp/L6b/pp/output/Ion_species_IC.txt
sed -i "s|^Spectrum file: .*|Spectrum file: $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt|" \
    /tmp/L6b/pp/input.inp
sed -i -e '/^Restart intent:/d' -e '/^Solver:/d' -e 's/^Do only PP: .*/Do only PP: True/' \
    /tmp/L6b/pp/input.inp
cd /tmp/L6b/pp && OMP_NUM_THREADS=8 $EX/EXHALE_L6b.x

# the test suite
EXHALE_OBJDIR=$EX/build_L6b EXHALE_TEST_OUT=/tmp/L6b/probe \
    bash $EX/src/tests/physics_probe/run.sh
```

ChiantiPy 0.15.2 reads a `chiantirc` file and raises `TypeError: 'bool' object
is not iterable` when there is none; a `$HOME/.config/chiantirc` with a
`[chianti]` section carrying `abundfile`, `ioneqfile`, `wavelength`, `flux` and
`gui` settles it. Limit the BLAS thread count when running the builder: left
alone, `populate` takes every core on the machine.
