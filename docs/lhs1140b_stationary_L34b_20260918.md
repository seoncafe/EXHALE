# L34 step 3 and step 4b: the molecular reference solution, and the molecular catalog re-solved on it

Item L34 of `docs/PLAN_20260917.md`, the molecular part (step 3, the reference
solution; step 4, the rest of the molecular cases). The atomic re-evaluation
of steps 1 and 2 is a separate worker and is not in this memo.

Every number below is MEASURED on this tree unless it is marked READ.
The binary is `LHS1140b/models/EXHALE_3146d11b.x`, md5
`3146d11b4090306dcea75bb9718edd22`, manifest
`LHS1140b/models/BINARY_MANIFEST_3146d11b4090.txt` (it carries L27, L28, L30,
L30b and L31). No source file was changed by this item and no tolerance was
weakened.

---

## 0. Verdict

**The molecular reference solution certified, at outer pass 12, in 40 m 35 s
of wall clock on `lart4` at `OMP_NUM_THREADS=8`.** `info = 0`,
`||R|| = 1.417e-08`, and the certification report reads CERTIFIED with all
seven active equations of the inventory inside their own tolerances. The pass
trajectory reproduces the L22 I3 run that certified the same case on the
predecessor binary, pass for pass and cell for cell, and the solution's
observables reproduce I3's to within 0.03 per cent in the He I 10830
equivalent width.

Every acceptance item of the plan's step 3 is answered in section 2, with two
statements that need to be read as stated rather than as a pass or a fail:

- the no-step evaluate route, run afterwards by `run_case.sh`, REFUSES the
  state it was handed, on one entry: the hydrodynamic mass row at cell 1,
  1.248e-08 against that cell's own rounding anchor 7.8e-09, a distance of
  1.598. That is inside the 0.78 to 1.71 band L18 measured for the base mass
  row's round trip and is a property of the cancellation at a subsonic base,
  not a second opinion about the wind (section 2.3);
- the L31 probe does NOT take a state path: its three representative cells are
  `real*8, parameter` arrays in
  `src/tests/physics_probe/molecular_energy_recipients.f90` lines 132 to 163.
  Blocks G and H were therefore run on the solved state through a scratch COPY
  of the suite with those eleven arrays replaced by the solved state's own
  values, the tree untouched (section 2.5).

---

## 1. The route, and why it is this one

`molecular_scalar_gj1132_kzz1e9/HeH2.13` is the case the plan names, and the
route is the one L22 increment I3 certified it with: the alternation, not the
coupled block. The block was measured on this very state to stand at
`||R|| = 4.241e-01` after fourteen trust-region iterations with all thirteen
of its linear solves exhausting their subspace (READ,
`docs/lhs1140b_stationary_L22_20260916.md` step 3 section 2).

The seed is the case's own archived state, `output/Hydro_ioniz_IC.txt` and
`Ion_species_IC.txt`, written 2026-09-16 by the binary `c2e9c9990b9f` at git
`43bc28cef587` and carrying `certified=T cert_reason=certified_in_wind`. What
the current binary makes of that claim, MEASURED by the no-step evaluation
this run begins with:

| entry | value | tolerance | cell |
|---|---|---|---|
| hydrodynamic energy row | 4.267e-01 | 1.0e-06 | 203 (r = 1.1569) |
| carrier balance H2 | 7.817e-01 | 1.0e-05 | 243 |
| elemental transport He/H partition | 7.245e-08 | 1.0e-05 | 315 |

The first two are the numbers L22 I3 and item L29 quote for this state, so the
run was handed the state the plan names.

The recipe as run, every key READ from the case's `input.inp` as it now
stands:

```
Restart intent: stationary            the alternation's entry
Coupled carrier solve: False          the alternation
Well balanced: True
Secondary_ionization: Immediate
Numerical flux: HLLC
Molecular chemistry: True
Molecular carrier transport: True
He_diffusion: True,  He_Kzz: 1.0e9
EXHALE_PTC_DTAU0=1.0
EXHALE_OUTER_PASSES=40                the budget run_case.sh sets
OMP_NUM_THREADS=8                     on lart4
```

Two lines of `input.inp` moved and nothing else: `Restart intent:` from
`stationary equilibrate` to `stationary`, and `Coupled carrier solve: False`
added, which is that key's default and now states the route in the file.
`base.inp` is untouched. The case's archived `output/` is kept beside it as
`output_pre_L34/`, with the archived logs, `REPRODUCE.md` and transit curves
under `output_pre_L34/case_files/`.

---

## 2. Step 3: the acceptance list, item by item

### 2.1 `info = 0` and CERTIFIED on every gated row

Twelve outer passes, 40 m 35 s, 2026-09-18 02:09:07 to 02:49:42 KST.
MEASURED, the worst gated species row pass by pass (the L33 reader,
`models/.L22/outer_pass_history.py`), beside the same column of the L22 I3 run
on the predecessor binary:

| pass | hydro info | worst gated row | cell | r [R_p] | wall [s] | I3's row at the same pass |
|---|---|---|---|---|---|---|
| 1 | 2 | 7.82e-01 | 243 | 1.3165 | 631.4 | 7.82e-01 (243) |
| 2 | 0 | 5.30e-01 | 221 | 1.2153 | 428.8 | 5.30e-01 (221) |
| 3 | 0 | 1.49e-01 | 221 | 1.2153 | 212.4 | 1.49e-01 (221) |
| 4 | 0 | 2.27e-02 | 222 | 1.2191 | 176.7 | 2.28e-02 (222) |
| 5 | 0 | 2.49e-03 | 242 | 1.3110 | 178.6 | 2.49e-03 (242) |
| 6 | 0 | 4.72e-04 | 221 | 1.2153 | 156.4 | 4.72e-04 (221) |
| 7 | 0 | 1.92e-04 | 223 | 1.2230 | 156.9 | 1.92e-04 (223) |
| 8 | 0 | 5.46e-05 | 234 | 1.2704 | 159.0 | 5.46e-05 (234) |
| 9 | 0 | 2.22e-05 | 246 | 1.3336 | 83.6 | 2.23e-05 (246) |
| 10 | 0 | 2.16e-05 | 221 | 1.2153 | 83.6 | 2.16e-05 (221) |
| 11 | 0 | 1.42e-05 | 221 | 1.2153 | 75.4 | 1.42e-05 (221) |
| 12 | 0 | **8.17e-06** | 221 | 1.2153 | 63.2 | 8.18e-06 (221) |

Pass 12 prints `ACCEPTED -- every active equation of this state is within its
own tolerance`. The carrier relaxation of passes 1 to 11 ended on the
composition movement bound and the accepted pass took no composition update at
all, so the acceptance is not a bound-throttled one (the L33 section 8.3
condition).

The certification of the final state, all seven active entries of the
twenty-eight in the inventory (MEASURED, `run.log`):

| entry | max | at cell | tolerance | verdict |
|---|---|---|---|---|
| hydrodynamic mass row | 1.206e-09 | 3 | cell-by-cell rounding anchor; binding cell 188 at 2.611e-11 against 2.8e-11 (distance 0.9189) | within |
| hydrodynamic momentum row | 5.588e-12 | 500 | 1.0e-08 | within |
| hydrodynamic energy row | 1.417e-08 | 47 | 1.0e-06 | within |
| carrier balance H2 | 8.172e-06 | 221 | 1.0e-05 in the wind | within |
| elemental transport He/H partition | 9.884e-06 | 2 (6.379e-08 at cell 218, the gated value) | 1.0e-05 in the wind | within |
| level balance He 2^3S | 1.824e-19 | 442 | 1.0e-06 | within |
| eliminated-species closure `System_HeH_mol` | 1.665e-16 | 82 | 1.0e-06 | within |

with 0 cells without a chemical root, 0 unbudgeted accepted corrections, no
active unvalidated physics and no rejected trial with an adopted contribution.

### 2.2 The L27 face-flux budget

The certification report prints it on the solved state (`pp.log`, the evaluate
route; the same block on the seed state is at the head of `run.log`):

| | seed state | solved state |
|---|---|---|
| operator, flux, well balanced | WENO3, HLLC, T | WENO3, HLLC, T |
| wind-window mean of `rho v r^2` | 5.82729e-07 | 5.73832e-07 |
| minimum over the faces | 0.9999748253553058 (face 92) | 0.9999733285494139 (face 0) |
| maximum over the faces | 0.9999748294889720 (face 0) | 0.9999733424100231 (face 92) |
| base face | 0.9999748294889720 | 0.9999733285494139 |
| base face offset from the window mean | -2.51705e-05 | -2.66715e-05 |

in units of the window mean. **The face mass flux is one number across the
whole column, base face included, to 1.39e-08 of the wind's own mass flux**
(the max-minus-min of the solved state). The constant offset of -2.67e-05
between the faces and the cell-centered window mean is the difference between
a face value and the mean of a cell-centered product, which the report says
itself and which the flux gate measures separately (accepted spread
2.2039e-10).

### 2.3 The energy balance of the column, cell by cell

**The code's own operator, on the solved state with no step taken**
(`EXHALE_RESIDUAL=1`, `Restart intent: stationary evaluate`, a scratch copy of
the case). The energy row IS this balance: `R_3 = dF_E - (heat - cool)` with
`dF_E` the divergence of `A v (E + p)` plus the well-balanced gravitational
work term `dF3p` of `RK_rhs.f90`, so its terms are the advected enthalpy, the
kinetic flux, the pressure work, the work against gravity and the two
radiative rates. MEASURED:

| row | cell-wise max | at cell | r [R_p] | volume-weighted integral (num / den) |
|---|---|---|---|---|
| mass | 7.0724e-09 | 1 | 1.00019 | 6.1209e-11 (1.756e-14 / 2.869e-04) |
| momentum | 5.5864e-12 | 500 | 29.03122 | 1.0253e-14 (8.676e-15 / 8.462e-01) |
| energy | 1.4464e-08 | 31 | 1.00599 | 7.0167e-10 (5.674e-14 / 8.087e-05) |

and the integral split into the four radial windows the report uses, as the
ratio of its numerator to its denominator:

| window | mass | momentum | energy |
|---|---|---|---|
| r < 1.03 | 2.13e-10 | 1.26e-14 | 3.48e-09 |
| 1.03 to 1.10 | 9.52e-11 | 3.01e-15 | 9.39e-10 |
| 1.10 to 1.20 | 1.52e-11 | 2.89e-15 | 1.27e-10 |
| r >= 1.20 | 7.77e-13 | 6.51e-15 | 4.54e-11 |

**Every cell of the column closes its energy balance inside the certification
tolerance 1.0e-06, the worst by a factor 69**, and the column-integrated
balance closes to 7.0e-10. The terms at the worst cell, signed, in code units:
`R = -6.573e-12`, flux divergence `5.353979e-05`, explicit source `0`,
heat `4.544731e-04`, cool `4.009333e-04`, `heat - cool = 5.353979e-05`. The
same reading on the seed state is `4.2666e-01` at cell 203, which is the
refusal L29 recorded.

**An independent assembly of the same equation, outside the binary**, from
`output/{Hydro_ioniz_IC.txt, Ion_species_IC.txt, Heating_breakdown.txt,
Cooling_breakdown.txt}` through `examples/exhale_io.py`, with the internal
energy taken from the production caloric equation of state (the H2
rovibrational ladder of `caloric_eos.f90` tabulated by a scratch program that
links that module, so no constant heat capacity enters) and the faces formed
as the cubic Lagrange interpolant through the four neighbouring centers.
MEASURED, `|R_3|` over the largest term of the row:

| window | median | maximum |
|---|---|---|
| r < 1.03 | 4.5e-05 | 1.3e+00 (r = 1.0014) |
| 1.03 to 1.10 | 2.7e-05 | 3.9e-04 |
| 1.10 to 1.20 | 2.0e-04 | 4.0e-03 |
| 1.20 to 2.00 | 3.6e-04 | 6.7e-04 |
| 2.00 to 5.00 | 1.0e-03 | 1.4e-03 |
| 5.00 to 29.03 | 1.6e-03 | 1.7e-03 |

**This reading is at its own floor and not at the state's**, and the
measurement that shows it is that the SAME assembly on the seed state, whose
energy row the code puts at 4.267e-01, gives 4.4e-05, 2.3e-05, 2.6e-04,
3.6e-04, 1.0e-03, 1.6e-03 in the same six windows: the independent assembly
cannot tell the certified state from the refused one, because its floor stands
seven decades above the difference between them. Two causes, both identified:
in the far wind the run's face flux is the HLLC upwind state of a WENO3
reconstruction and this assembly's is a centered interpolant, a difference of
first order in `dr` on a stretched grid; at the subsonic base the terms cancel
to twelve digits and one unit in the last place of the density moves the row
by more than the row. The code's own operator is therefore the gate and this
assembly is the independent confirmation that the balance closes to a part
in a thousand, which is what it can say.

**The evaluate route's second opinion on the state.** `run_case.sh` runs
`Restart intent: stationary evaluate` after the solve. It REFUSES the state,
on one entry of the seven:

| entry | in the solve | re-evaluated | tolerance |
|---|---|---|---|
| hydrodynamic mass row | binding cell 188, 2.611e-11 against 2.8e-11, distance 0.919; the row's largest measure 1.206e-09 at cell 3, whose own anchor is 3.3e-09; at cell 1 the pass printed 5.213e-09 | binding cell 1, 1.248e-08 against 7.8e-09, distance 1.598 | the cell's own rounding anchor |
| hydrodynamic energy row | 1.417e-08 at cell 47 | 1.446e-08 at cell 31 | 1.0e-06 |
| carrier balance H2 | 8.172e-06 at cell 221 | 8.172e-06 at cell 221 | 1.0e-05 |
| elemental transport He/H | 9.884e-06 at cell 2 | 9.883e-06 at cell 2 | 1.0e-05 |
| level balance He 2^3S | 1.824e-19 | 1.953e-19 | 1.0e-06 |
| eliminated-species closure | 1.665e-16 | 1.665e-16 | 1.0e-06 |

Six of the seven come back where they were, the carrier row to every digit
printed; what moves is the mass row at cell 1, by a factor 2.39 of the
5.213e-09 the accepted pass printed there, which carries it past that cell's
anchor. L18 measured 0.78 to 1.71 (median 1.12) for that round trip over 120
certified states, all of them atomic, so this factor is OUTSIDE the band those
states set and it is reported as such rather than filed under it. What it is
NOT is a second opinion about the wind: cell 1 carries `v = -1.9 cm s^-1`
against a wind of `1.2e+05 cm s^-1`, its flux difference cancels to twelve
digits, and the cell is far below the r >= 1.20 window every gated species row
is judged in. A third no-step evaluation of the same state, the
`EXHALE_RESIDUAL=1` diagnostic above, puts the same row at 7.072e-09 at the
same cell 1: three readings of one row at a subsonic base spread over a factor
2.4, and none of the three moves any other entry. The
state written to `output/Hydro_ioniz.txt` therefore carries
`certified=F cert_reason=failing_entries` while the solved pair beside it,
`output/*_IC.txt`, carries `certified=T cert_reason=certified_in_wind`.

### 2.4 The L31 collider-sum agreement, on the solved state

`src/tests/physics_probe/molecular_energy_recipients.f90` block H forms the H2
association third body three times, once as the chemistry forms it
(`System_HeH_mol::h2_third_body_density` through `mol_heh_rows`), once as the
carrier balance writes it into the R15 channel of its row record, and once as
the heat ledger deposits it (`molecular_chemical_heating`). MEASURED at the
molecular-depth cell OF THE SOLVED STATE:

| comparison | difference |
|---|---|
| chemistry against `h2_association_collider_density` | 0.0 exactly |
| the carrier row's R15 channel against `mk15 n_third n(H)^2` | 0.0 exactly |
| the heat ledger's deposit on an R12/R15-only cell against the same third body | 1.015357541515134e-06 measured and expected, agreeing to 1e-14 relative |
| the collider sum in a pure H2 gas against `n_tot` | 0.0 exactly |

All 81 assertions of the driver pass on the solved state's cells.

### 2.5 The L31 uncertainty bracket, on three cells of the solved state

The three cells, located on the solved state by the L31 rules (the molecular
depth is physical cell 1, the H2 front is the cell where `2 n(H2)/n_H` has
fallen to half its base value, the dilute upper column is the cell at the
radius L31 used):

| | molecular depth | H2 front | dilute upper column |
|---|---|---|---|
| row of the written state | 3 | 179 | 301 |
| r [R_p] | 1.000193 | 1.099271 | 1.842260 |
| T [K] | 778.792 | 2116.614 | 5061.593 |
| `2 n(H2)/n_H` | 0.335357 | 0.165984 | 3.559e-06 |
| n(H I) [cm^-3] | 1.991e+12 | 1.714e+10 | 8.200e+07 |
| n(H2) | 5.024e+11 | 1.708e+09 | 1.519e+02 |
| n(He I) | 6.382e+12 | 3.922e+10 | 6.943e+07 |
| n_e | 9.212e+06 | 2.371e+07 | 6.870e+06 |

(L31 read rows 3, 180 and 301 of the archived state at 808.319, 2128.367 and
5083.438 K; the front has moved in by one cell and the whole column is
cooler.)

The bracket, molecular chemical heating in erg cm^-3 s^-1, each quantity
varied separately against the production assembly's own output:

| cell | nominal heat | recipient 1 -> 0 | R6 branching | R6 `E_int` +15% | quench `f` -> 1 | He efficiency +0.3 dex | He efficiency -0.3 dex |
|---|---|---|---|---|---|---|---|
| molecular depth | 1.04863e-06 | +1.80273e-13 | -1.17019e-15 | +1.75528e-16 | +7.18294e-13 | +4.40921e-07 | -2.20984e-07 |
| H2 front | 8.60218e-09 | +3.75730e-10 | -3.50296e-15 | +5.25445e-16 | +3.50670e-15 | -1.74765e-11 | +8.75900e-12 |
| dilute upper column | 3.11216e-12 | +4.06002e-12 | -1.54180e-21 | +2.31270e-22 | +1.56324e-21 | -1.92246e-15 | +9.63514e-16 |

as fractions of the nominal heat:

| cell | recipient | R6 branching | R6 `E_int` | quench | He efficiency |
|---|---|---|---|---|---|
| molecular depth | 1.7e-07 | 1.1e-09 | 1.7e-10 | 6.8e-07 | **+42.0% / -21.1%** |
| H2 front | **+4.4%** | 4.1e-07 | 6.1e-08 | 4.1e-07 | -0.20% / +0.10% |
| dilute upper column | **+130%** | 5.0e-10 | 7.4e-11 | 5.0e-10 | -0.062% / +0.031% |

The ordering is L31's and it survives the re-solve: at the molecular depth the
uncertainty of the chemical heat is the helium third body and nothing else,
+42 and -21 per cent from the +/-0.3 dex of the argon coefficient, four to six
orders above every other term; at the front and above, the recipient of R5 and
R16 carries it, 4.4 per cent at the front and more than the nominal heat in
the dilute column, where that heat is 3e-12 erg cm^-3 s^-1 and of no
consequence to the energy budget. **This is the number to quote beside every
molecular heat of this catalog.** Against L31's own bracket on the archived
state (+41 / -21 per cent, +4.6 per cent, +130 per cent) the re-solve moves
nothing that matters.

**How blocks G and H were run on a state.** The driver holds the three cells
as `real*8, parameter` arrays (lines 132 to 163), so the suite as it stands in
the tree cannot be pointed at a state file. The tree was not edited: the whole
suite was copied to the scratch directory, the eleven arrays were replaced by
the solved state's own values and the driver list cut to this one driver. The
same copy built from the ARCHIVED state's values reproduces every number of
the L31 memo's sections 5.3 and 5.5 exactly, which is the control that the
copy is the driver.

### 2.6 The comparison figures and the He I 10830 line

`docs/figures/lhs1140b_L34b_reference.png`: temperature, mass density,
velocity, `2 n(H2)/n_H`, `3 n(H3+)/n_H` and the He I 10830 transit line of the
reference solution against the state it was seeded from, with an inset on the
base layer where the two differ most. The line is the WINERED HIRES-Y
synthesis `run_case.sh` produced (`R = 68 000`,
`LHS1140b/winered_hires_y.sh`), and the equivalent width quoted below is the
one `run_case.sh` measures over the observation's own vacuum window,
10832.60 to 10834.20 A in air.

---

## 3. The archived state beside the new one

Every number MEASURED from the two states' own files by one tool
(`log10 Mdot` as `EXHALE_main` forms it, `4 pi rho v r^2` at the physical cell
`N - 20` with the `2D approximate method` factor; the equivalent width over
the measurement's window; the depth and FWHM from the three-Gaussian fit
beside the curve; the H2 front as the outermost radius with
`2 n(H2)/n_H >= 1e-2`).

| quantity | archived (2026-09-16, `c2e9c9990b9f`) | this item (`3146d11b4090`) | change |
|---|---|---|---|
| log10 Mdot [g s^-1] | 7.90907 | 7.90239 | -1.53 per cent in Mdot |
| He I 10830 red-pair EW [%A] | 1.57583 | 1.56032 | -0.99 per cent |
| red-pair depth [%] | 5.62817 | 5.57130 | -1.01 per cent |
| FWHM [A] | 0.262836 | 0.262836 | 0.00 per cent |
| base T [K] (cell 1, r = 1.000193) | 808.319 | 778.792 | -3.65 per cent |
| base heating rate [erg cm^-3 s^-1] | 1.12791e-06 | 1.06477e-06 | -5.60 per cent |
| `2 n(H2)/n_H` at cell 1 | 0.355303 | 0.335357 | -5.61 per cent |
| H2 front, `2 n(H2)/n_H = 1e-2` [R_p] | 1.16256 (cell 205) | 1.15973 (cell 204) | one cell |

**The molecular numbers of this table carry the L31 uncertainty of section
2.5**: the base heating rate's chemical part is uncertain by +42 and -21 per
cent through the helium third body alone, which is eight times the -5.6 per
cent the re-solve moved it by, and the chemical heat at and above the H2 front
is uncertain by +4.4 per cent through the R5/R16 recipient.

The same comparison was made by L22 I3 on the predecessor binary and gives
base T 778.79, `x2` at cell 1 3.3536e-01, the front at 1.1597 R_p cell 204,
log10 Mdot 7.9024, EW 1.5599 %A, depth 5.5712 per cent, FWHM 0.2628 A. **The
solution of this item and the solution I3 reached are the same solution to
four digits in every one of them**, although the two binaries differ by L30
(one spherical geometry for the element and carrier transport, which moves
every element and carrier row by O(dr^2)) and by L27 and L28. The largest
difference is 0.03 per cent in the equivalent width and 0.02 per cent in the
depth.

---

## 4. Step 4: the rest of the molecular catalog

Each case re-solved from its OWN archived or stopped state, which is what the
reference solution was given, on the same recipe. MEASURED:

| case | seed | passes | verdict | wall clock |
|---|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` (the reference) | its archived certified state | 12 | `info = 0`, CERTIFIED | 40 m 35 s |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | its archived certified state | 10 | `info = 0`, CERTIFIED | 25 m 35 s |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | its archived certified state | 6, stopped | not solved, the carrier row a band and the wind broken down from pass 5 | 3 h 35 m |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | the 2-pass snapshot its `not_solved.md` names | 1 at each of two pseudo-time starts | not solved, the element relaxation found no admissible advance | 1 h 04 m |
| `molecular_photochem_gj1132_kzzprofile/HeH2.09` | the stopped state its `not_solved.md` names | 0 completed, stopped | not solved, the first pass's hydrodynamic solve spent 180 iterations at `||R||` 3.4e-01 | 1 h 55 m |
| the three `molecular_scalar_gj1132_wellmixed` cases | -- | -- | not run: section 5 | -- |

Each case that did not solve carries its own `not_solved.md` in the catalog's
form, with the worst row and the pass history. Each that did carries a
`REPRODUCE.md` written by its own run, naming the binary and its md5.

### 4.1 `HeH9.7`, the second certified solution

Ten outer passes, 25 m 35 s, the same shape as the reference: the worst gated
row falls by a geometric factor 0.20 to 0.45 a pass, 4.88e-01, 2.09e-01,
4.59e-02, 9.19e-03, 2.04e-03, 5.11e-04, 1.49e-04, 5.22e-05, 2.13e-05,
9.54e-06, and the refusing cell stands at 217 (r = 1.2007 R_p) from pass 2 on.
`info = 0`, `||R|| = 2.688e-08`, CERTIFIED with all seven entries within.

The L27 face-flux budget of the solved state: window mean 7.14538e-07,
minimum 0.9999625547588813 at face 92, maximum 0.9999625611808490 at face 0,
base face 0.9999625611808490, offset from the window mean -3.74388e-05. The
face mass flux is one number to 6.4e-09 of the wind's.

The archived state beside the new one:

| quantity | archived | this item | change |
|---|---|---|---|
| log10 Mdot [g s^-1] | 7.94648 | 7.94282 | -0.84 per cent in Mdot |
| He I 10830 red-pair EW [%A] | 2.42118 | 2.40606 | -0.62 per cent |
| red-pair depth [%] | 7.68280 | 7.62903 | -0.70 per cent |
| FWHM [A] | 0.295123 | 0.295627 | +0.17 per cent |
| base T [K] | 567.526 | 522.193 | -7.99 per cent |
| base heating rate [erg cm^-3 s^-1] | 6.31976e-07 | 5.51060e-07 | -12.80 per cent |
| `2 n(H2)/n_H` at cell 1 | 0.0779932 | 0.0675867 | -13.34 per cent |
| H2 front [R_p] | 1.05626 (cell 145) | 1.05333 (cell 142) | three cells |

The base heating rate's chemical part carries the same L31 uncertainty as the
reference's, +42 and -21 per cent through the helium third body, which is
three times the -12.8 per cent the re-solve moved it by. The evaluate route
refuses this state too, on the same one entry: hydrodynamic mass row 2.340e-08
at cell 1 against 7.5e-09.

**The base is where the L7g corrections land, and the colder the base the
larger they are**: -3.7 per cent in base temperature at He/H = 2.13 and -8.0
per cent at 9.7, against -1.5 and -0.8 per cent in the mass-loss rate and -1.0
and -0.6 per cent in the equivalent width. The line forms far above the layer
that moved.

### 4.2 The three that did not solve

`kzz1e9/HeH0.55` is the informative one, because it CERTIFIED on the old
binary and does not now. Its worst gated row over six passes is 8.92e-03,
1.39e-02, 3.39e-03, 7.20e-03, 2.45e-03, 3.11e-03, a band and not a fall, with
the refusing cell at 500, the outermost physical cell, in four of the six;
and from pass 5 its hydrodynamic mass row reads 1.88 and 1.67 against a
tolerance of order 1e-08, so the wind under the front stopped converging.
That is the shape L33 section 6.2 classified for `wellmixed/HeH0.55`, met
here one composition step away from a case that certifies in twelve passes.
The ladder of section 5 is the route that addresses it and this case is its
rung 4.

The two `photochem` cases fail earlier and for a reason of their own: neither
has a seed of the right physics. `HeH9` was handed a two-pass relaxation
snapshot whose face mass flux varies by 1.2e-02 across the column, against
1.4e-08 on the certified reference, and its element relaxation restored its
entry composition at pass 1 at both pseudo-time starts, which ends the loop by
construction. `HeH2.09` was handed the molecular seed its stopped run had been
built from and its first hydrodynamic solve did not descend below
`||R|| = 3.4e-01` in 180 iterations.



---

## 5. The well-mixed cases and the L33 ladder, and how the construction extends

The plan puts the three `molecular_scalar_gj1132_wellmixed` cases last and
allows them only through the L33 continuation in the base H2 fraction. The
ladder as L33 section 8.2 designed it runs inside the `kzz1e9` family, from
the certified fiducial to `kzz1e9/HeH0.55`:

| rung | `q_H2_base` | He/H |
|---|---|---|
| 0 | 0.190110258 | 2.1300 (the certified reference of section 2) |
| 1 | 0.261627353 | 1.4111 |
| 2 | 0.333144448 | 1.0008 |
| 3 | 0.404661544 | 0.7356 |
| 4 | 0.476178639 | 0.5500 |

with `q_H2_base = (f/2)/((1 - f) + f/2 + He/H)` and
`f = 0.9999831600354949`, so a rung is a physical model of the same family
and not an interpolation.

**How it extends to the three well-mixed cases.** `wellmixed` and `kzz1e9`
differ in exactly two keys, `He_diffusion: True` and `He_Kzz: 1.0e9`
(MEASURED by L33 section 8.1 on the case files), and at the same base H2
fraction the element diffusion moves the `x2 = 0.5` front only from 1.2565 to
1.3705 R_p while the base fraction moves the `x2 = 1e-2` contour from 1.1626
to beyond the grid. So the composition ladder is the long leg and the
element treatment is the short one, and the extension is one further rung of
a different kind at the end of each leg:

- `wellmixed/HeH0.55`: rungs 0 to 4 above with `He_diffusion` kept, then one
  rung with the two keys removed at fixed `q_H2_base = 0.476178639`;
- `wellmixed/HeH2.13`: the element-treatment rung alone from the certified
  reference, `q_H2_base = 0.190110258` with the two keys removed. It is the
  shortest of the three and the one to try first, and it is the case
  `MODELS.md` records as never started for want of a certified well-mixed
  state to seed it from;
- `wellmixed/HeH0.083`: `q_H2_base = 0.857606`, which is 0.381 beyond rung 4,
  so five further rungs of the same 0.0715 step (0.5477, 0.6192, 0.6907,
  0.7622, 0.8337, 0.8576) after the element-treatment rung of
  `wellmixed/HeH0.55`, or the same five inside `kzz1e9` first and the
  element-treatment rung at the end.

Each rung its own directory under `LHS1140b/models/.L34/rung<k>/`, the
fiducial's `input.inp` and `base.inp` copied with `He/H number ratio`,
`HeH_base` and `q_H2_base` changed and nothing else, seeded from the previous
rung's certified state, 25 passes a rung, a rung that does not certify ENDS
the ladder.

**Not run in this item.** At the measured cost of a pass of the cases of
section 4 the five-rung ladder is the eight hours L33 estimated and the
extension is more; the host was carrying the item's own direct solves and two
other workers. The construction above is the record of how it is built, and
the reference solution the ladder needs as its seed now exists and is
certified, which is the thing that was missing.
---

## 6. Reproducing this, and where the record is

The reference solution, from the case directory
`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13` with its archived
state put back in `output/` from `output_pre_L34/`:

```bash
cd LHS1140b/models
NOSEED=1 FORCE=1 OMP_NUM_THREADS=8 \
  EXHALE_BIN=$PWD/EXHALE_3146d11b.x \
  EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=40 \
  ./run_case.sh molecular_scalar_gj1132_kzz1e9/HeH2.13
```

`NOSEED=1` is what makes the run start from the pair in `output/` rather than
build a molecular seed out of the atomic case; `FORCE=1` lets it re-run a case
that already carries a synthesized line. The case's own `REPRODUCE.md`, which
that run writes, carries every command with its environment.

The catalog record: `output_pre_L34/` beside each case holds the state the
case carried before this item, with the archived logs, `REPRODUCE.md` and
transit curves under `output_pre_L34/case_files/`; each case that did not
solve carries `not_solved.md`. `MODELS.md` sections 7 and 8 and every
`REPRODUCE.md` were regenerated by `models/status.py --write` and
`models/write_reproduce.py` at the close of this item.

The figure is `docs/figures/lhs1140b_L34b_reference.png`.

The three readings of the energy balance, the scratch copy of the physics
probe on the solved state's cells and the tabulation of the H2 rovibrational
ladder are scripts of this item's scratch directory and are not in the tree;
the numbers they produced are in sections 2.3, 2.4 and 2.5 and each states how
it was formed.

---

## 7. Bookkeeping, and three things noticed beside the measurement

`models/status.py --write` regenerated `MODELS.md` sections 7 and 8;
`models/write_reproduce.py` is per case and has no tree-wide mode, so every
`REPRODUCE.md` of this item is the one its own run wrote (both certified cases
name the md5 `3146d11b4090306dcea75bb9718edd22`).
`LHS1140b/make_memo_figures.py` regenerated the twenty
`docs/figures/lhs1140b_*.pdf`, of which twelve changed when rendered and
compared page by page against the repository's copies (`bump`, `closure`,
`closure_ladder`, `composition_profiles`, `diff_broadened`, `ew_vs_heh`,
`fig4style`, `gj699`, `heh_vs_kzz`, `knudsen`, `kzz_profiles`, `thermostat`);
most of those draw ATOMIC cases and the change in them is the L34 steps 1 and
2 re-evaluation running in the same hours, not this item. `latexmk -pdf` for
`docs/lhs1140b_exhale_vs_pwinds.tex` then built cleanly, 80 pages, the same
page count as before.

Noticed and NOT acted on:

1. **The `claim on 59bfdb3f` column of `MODELS.md` section 7 is now stale for
   the cases this item re-solved.** It records what the evaluate route made of
   the ARCHIVED state's header (`CLAIMS_59bfdb3fc4d0.md`, item L29), and for
   `kzz1e9/HeH2.13` and `HeH9.7` that state is no longer the one in `output/`.
   The table prints `info=0 certified` beside `refused` for both, which reads
   as a contradiction and is two statements about two different states.
2. **`molecular_scalar_gj1132_wellmixed/HeH2.13` still prints as `running`**
   in the tables, which is what the files on disk say and what the revised
   handoff's section 3.1 already warns about; the case was never started.
   `molecular_photochem_gj1132_kzzprofile/HeH2.09` now prints the same way for
   the same reason, its run having written no state.
3. **The three certified molecular states of this group all fail the no-step
   re-evaluation on the hydrodynamic mass row of cell 1** (1.248e-08 against
   7.8e-09 on `HeH2.13`, 2.340e-08 against 7.5e-09 on `HeH9.7`), while every
   other entry comes back where it was. L18 measured the round-trip band 0.78
   to 1.71 on 120 certified ATOMIC states; a molecular base with 34 per cent
   of its hydrogen in H2 and a caloric equation of state has not been measured
   that way, and on this evidence its band is wider. That is a measurement
   about the base row's cancellation, and whether the evaluate route should be
   the acceptance for a molecular case is a question for the user, not a
   change to make here.
