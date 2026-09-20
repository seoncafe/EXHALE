# D7d: the helium state the metal + He and the advected He <-> H reactions react from

Item D7d of `docs/PLAN_20260918_rev2.md`: the two sites the D7c report listed
as noticed outside its own scope.

## Verdict

Both are real and both are fixed.

`System_implicit_adv_HeH_TR`, the advection-corrected ionization system that
`post_process_adv` integrates upwind across each cell, handed the He <-> H
pair of Huang et al. (2023), ApJ 951, 123, Table 4 group B the SUMMED He I
fraction while its row 2 is the ground-singlet balance. It now passes the
singlet `xheiS`, so the advection correction reacts the pair from the same
reservoir the equilibrium systems do.

The metal + He reactions of Table 4 group C were handed the summed He I by
all three systems that assemble them. Those rates are ground-state rates, by
the same argument and the same arithmetic as group B: all three now pass the
ground singlet, so one reservoir feeds the whole reaction set.

The group C correction changes nothing unless `cx_full 1` is set in
`metals.inp`, because without it those rows are not in the active set and the
helium density handed to the generic assembly is never read. The advection
correction changes the `*_adv.txt` profiles wherever the He 2^3S metastable
is tracked, which is every triplet run, and those profiles are what the
transit tool reads.

## 1. Which helium state the rates are for

READ from the journal PDF, `references/Huang_2023_ApJ_951_123.pdf`, Table 4.
Every group C reactant is written as the bare element symbol, `Si + He+`,
`Si+ + He`, `C + He+`, `C+ + He`, `O + He+`, `O+ + He`, exactly as the group
B rows write `He + H+` and `He+ + H`.

The barriers of the three reverse rows settle which state that is. A
reaction X+ + He -> X + He+ is endothermic by the difference of the two
ionization potentials, so the exponential of the fit is that difference in
temperature units:

| Table 4 row | barrier in the fit | IP(He) - IP(X) | in kelvin |
|---|---|---|---|
| C2, Si+ + He | exp(-19.1/T4) | 24.587 - 8.152 = 16.436 eV | 1.907e5 K |
| C4, C+ + He  | exp(-15.5/T4) | 24.587 - 11.260 = 13.327 eV | 1.547e5 K |
| C6, O+ + He  | exp(-12.7/T4) | 24.587 - 13.618 = 10.969 eV | 1.273e5 K |
| B1, He + H+  | exp(-12.75/T4) | 24.587 - 13.598 = 10.989 eV | 1.275e5 K |

(barriers READ from Table 4 through `docs/charge_exchange_table4.md`;
ionization potentials READ from the standard values; the division by
8.6173e-5 eV/K is arithmetic done here, MEASURED.) Every one of them matches
the ground-state value to the two digits the fit carries. He(2^3S) lies
19.82 eV above the singlet, so its ionization potential is 4.77 eV and each
of those collisions is exothermic for it by 3.4 to 8.8 eV, with no barrier
and a rate of a wholly different size. That reaction is carried nowhere in
this code, and charging the metastable at a ground-state rate is not an
approximation of it.

## 2. What each site passed, and what it passes now

READ from the source at the entry text of this item.

| site | helium density at the entry text | now |
|---|---|---|
| `System_implicit_adv_HeH_TR.f90:155`, to `he_h_cx_fvec_adv` | `xhei = xheiS + xheiTR`, the summed He I fraction | `xheiS`, the ground singlet unknown (`:161`) |
| `System_HeH_TR_metals.f90:131`, to `cx_add_to_fvec` | `n_hei = (1-x2-x3) n_He`, the summed He I | `n_heiSI` (`:135`) |
| `System_HeH_mol_metals.f90:261`, to `cx_add_to_fvec` | `n_hei`, the free neutral helium | `n_heiSI` (`:266`) |
| `constrained_chemical_equilibrium.f90:1203`, to `cx_add_to_fvec` | `n_hei = sden(is_HeI_SI) + sden(is_HeITR)` | `sden(is_HeI_SI)` (`:1229`) |

The three group C sites were corrected in one change, not one at a time:
correcting one of them would have recreated exactly the divergence between
systems that item D7c removed.

`System_implicit_adv_HeH`, the advection system without a metastable, already
passed the singlet; all of its neutral helium is the singlet and its own
comment said so. There is no analytic Jacobian to match at any of the four
sites: the advection systems are solved by `hybrd1` with a numerical one, and
`cx_add_to_jac` is reached only from `System_HeH_metals`, whose helium is the
singlet because that system has no metastable.

The comment at `constrained_chemical_equilibrium.f90:1148` called the summed
reservoir a deliberate choice, "the He I reservoir the metal charge exchange
reacts with (the metastable is a level inside the neutral stage)". That
sentence is gone; the reactant is now stated where the assembly is called,
with the reason.

## 3. Tests

Two new rows of the D7a suite `src/tests/charge_exchange_rows/`. Both isolate
the reaction they are about by evaluating the same residual twice at one
state with the reaction switched on and off, so that the difference is that
reaction's contribution to each row and nothing else.

`advection_pair_reservoir.f90` runs `System_implicit_adv_HeH_TR` at the state
D7c's `he_h_pair_reservoir` uses, with c1 = dr/v = 1e4 s, and converts its
rows to physical sources by S(H I) = d(1) n_H/c1 and S(He I) = d(2) n_He/c1.
It asserts that those sources are the closed form built from the two Table 4
rates and the ground singlet, that adding a metastable at fixed singlet,
He+ and He2+ densities moves neither of them, and that they are the sources
`System_HeH_TR` makes out of the same physical state through its own row
basis.

`metal_helium_reservoir.f90` runs the three metal systems and the network's
balance rows, differencing `cx_full` on against off with oxygen as the only
element carrying atoms, which leaves C5 and C6 as the whole difference. The
state then isolates C6, the one group C row whose reactant is the helium: all
the oxygen is ionized so C5 has no neutral oxygen to start from, the hydrogen
is ionized down to 1 atom cm^-3 of H0, every photoionization, recombination
and collisional coefficient is zero and the He <-> H pair is off, so the rows
the sources are read from carry charge exchange alone. C6 is suppressed by
exp(-12.7/T4), so in an ordinary cell it is a part in 1e5 of the group C
source and the reservoir it reacts from cannot be resolved at all.

RED, MEASURED against a private build of the entry text: all six assertions
of `advection_pair_reservoir` fail, the advection-corrected He I source
sitting 1.22e-2 cm^-3 s^-1 below the equilibrium one, which is
k1 n(He 2^3S) n(H+) at that state; and six of the thirteen assertions of
`metal_helium_reservoir` fail at 1.2903e-1 relative in all three systems,
which is the metastable's share 1.111e7/8.611e7 of the neutral helium there.

GREEN, MEASURED with the delivered text: the whole suite is 151 assertions,
0 failures. The advection-corrected and the equilibrium system agree on the
pair's sources to 7e-17 of the rows; the three metal systems agree on the
group C sources bitwise, the printed difference being exactly zero for two of
them and 2.4e-16 relative for the network.

The network's rows had to be made reachable from a test: they are assembled
by a private routine that needs the module's own layout, field and row
scales. `constrained_chemical_equilibrium` therefore gained one public
diagnostic entry, `constrained_network_balance_rows_of_cell`, which sets the
field to unity and the row scales to one and calls the production assembly,
together with the species slots the test indexes. No production path calls
it, which is the arrangement `constrained_network_layout_of_cell` already
had.

## 4. Impact

Control and measured are two private builds taken back to back, with only the
four reactants between them: an md5 list of every `src/**/*.f90` before the
control build and after the measured one shows no production source moving
between them (MEASURED; two test drivers of other suites did move, and
neither is linked into either binary). Both run single-threaded
(`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`) on scratch copies.

### `wasp_full` (He 2^3S + metals), 400 steps, `cx_full` off

MEASURED: the data lines of all ten output files are byte-identical, the
`*_adv.txt` profiles included.

The group C part is byte-identical by construction: without `cx_full` those
rows are not in the active set, so the helium density handed to the generic
assembly is never read. The advection part is byte-identical for a reason
that belongs to this snapshot and not to the correction: the advection
correction is applied to 51 cells, all of them between 1.0000 and 1.0099 Rp
at 2194 to 2476 K, where the metastable is at most 9.2e-12 of the neutral
helium (MEASURED from the run's own `Ion_species_adv.txt`) and where
exp(-12.75/T4) is 1e-23 and below. The 452 cells that carry a metastable
worth anything, up to 0.2355 of the neutral helium at 1.5437 Rp, are
`adv_retained`: this is a 400-step relaxation snapshot, its mass row does not
hold there, and those cells keep the equilibrium composition.

### `mol_metals` (molecular chemistry + metals), 400 steps, `cx_full` off

MEASURED: the data lines of all nine output files are byte-identical. Five
cells are corrected, between 1.1749 and 2.7140 Rp at 472 to 896 K, one of
them carrying a metastable of 7.1e-2 of its neutral helium; at those
temperatures exp(-12.75/T4) is below 1e-60 and the pair contributes nothing
to the row at all.

(The golden of this case runs 12000 steps. Both builds were run at 400 here,
the same cap for each, which is a like-for-like comparison and not the
golden's state.)

### A converged, certified state: LHS 1140 b, `atomic_scalar_gj1132_kzz1e9/HeH2.13`

The case above corrects the composition of 51 cells out of 504, so it says
little about the advection change. This one is a Newton-converged stationary
model of this session's planet, re-run through `Restart intent: stationary
evaluate` from its own state so that no step and no solve is taken and only
the post-process is redone. It corrects 477 cells of 504, from 1.0050 to
30 Rp at 547 to 5730 K.

MEASURED, with log10(Mdot) = 7.87 g/s in both:

| file | verdict |
|---|---|
| `Hydro_ioniz.txt`, `Ion_species.txt`, and both `_IC` files | byte-identical |
| `Cooling_breakdown.txt`, `Heating_breakdown.txt` | byte-identical |
| `Ion_species_adv.txt` | moves; largest relative movement 1.413e-13 (`HII` at 1.0524 Rp), then `HeII` 8.98e-14 at 1.1249, `HeITR` 4.37e-14 at 1.1271, `HeIII` 2.98e-14, `HeI` 1.28e-15, `HI` 5.37e-16 |
| `Hydro_ioniz_adv.txt` | moves; `cool` 2.43e-13 at 1.0524 Rp, `heat` 8.80e-15, `p` 1.46e-15, `T` 1.19e-15 |

The equilibrium state is untouched, which is the expected division: this item
changes only the advection-corrected post-process and the group C rows.

The size is set by two measured factors. The metastable is at most 1.109e-5
of the neutral helium in this model, so that is the fraction of the He + H+
rate the correction removes; and that rate is itself barrier-suppressed at
these temperatures, exp(-12.75/T4) being 2e-10 at the hottest corrected cell,
so a part in 1e5 of it is a part in 1e13 of the composition.

### `wasp_full` with `cx_full 1`, 400 steps: the group C change on its own

`cx_full 1` appended to the case's `metals.inp` (the setup echo confirms
`charge-exchange full Table 4 mode = T`, MEASURED). This configuration also
runs the advection change, but that change is byte-identical on this case, as
the paragraph above measured, so everything below is group C.

`wasp_full` carries C, N, O, Mg, Ca, Na and Fe but no Si, so four of the six
group C rows are active: C + He+, C+ + He, O + He+ and O+ + He.

MEASURED, log10(Mdot) = 7.20 g/s in both, over the physical cells:

| column | max relative | where | the density there |
|---|---:|---:|---|
| `FeI` | 9.579e-5 | 1.2783 Rp | 1.8e-12 cm^-3 |
| `CaI` | 4.113e-5 | 1.3416 Rp | 2.6e-13 cm^-3 |
| `OI`  | 3.848e-5 | 1.4803 Rp | 6.3e-11 cm^-3 |
| `NI`  | 2.454e-5 | 1.4545 Rp | 1.4e-11 cm^-3 |
| `CI`  | 2.297e-5 | 1.2783 Rp | 6.6e-11 cm^-3 |
| `MgI` | 6.129e-6 | 1.2783 Rp | 1.8e-11 cm^-3 |
| `HeI` | 4.751e-9 | 1.2328 Rp | |
| `HI`, `HeII`, `HeITR` | 8.5e-11, 5.5e-11, 5.6e-11 | 1.217 to 1.218 Rp | |
| `v`, `heat`, `cool` | 8.0e-11, 3.3e-11, 6.7e-12 | 1.22 to 1.35 Rp | |
| `T`, `p`, `rho` | 5.0e-13, 3.8e-13, 3.4e-13 | 1.39 to 1.46 Rp | |

The neutral metals are where it shows, because they are what the reverse
group C rows X+ + He -> X + He+ make: their rate loses exactly the
metastable's share of the neutral helium, and this run's metastable is
7.217e-6 of it at the He 2^3S peak, 1.0588 Rp (MEASURED). The columns that
move at 1e-5 all sit at 1e-11 to 1e-13 cm^-3; the bulk state moves at 1e-13.
This is a 400-step relaxation snapshot, so the numbers bound this
configuration's sensitivity and are not the error of a quoted mass-loss rate.

## 5. Files changed

| file | change |
|---|---|
| `src/modules/nonlinear_system_solver/System_implicit_adv_HeH_TR.f90` | passes `xheiS` to `he_h_cx_fvec_adv`; `xhei` had no other reader, so its declaration and its assignment are gone; the comment states the reactant, the barrier that identifies it, and that row 4 gives the metastable no charge-exchange term |
| `src/modules/nonlinear_system_solver/System_HeH_TR_metals.f90` | passes `n_heiSI` to `cx_add_to_fvec`; `n_hei` had no other reader here either and is gone; the comment names the group C reactant |
| `src/modules/nonlinear_system_solver/System_HeH_mol_metals.f90` | the same one-argument change and comment; `n_hei` stays, as the intermediate the singlet is formed from |
| `src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90` | passes `sden(is_HeI_SI)`; the summed-reservoir local and the comment that called it deliberate are gone; a dead assignment of the same name in `hydrogen_helium_row_terms` removed; new public diagnostic `constrained_network_balance_rows_of_cell` and the species slots the row test indexes |
| `src/modules/radiation/charge_exchange.f90` | header only: the group C helium reactant with its three barriers, and the same statement for `he_h_cx_fvec_adv` |
| `src/tests/charge_exchange_rows/advection_pair_reservoir.f90` | new driver |
| `src/tests/charge_exchange_rows/metal_helium_reservoir.f90` | new driver |
| `src/tests/charge_exchange_rows/run.sh`, `README.md` | the two new test names, in the default set, and their sections |
| `docs/lhs1140b_stationary_D7d_20260918.md` | this memo |

No `Makefile` edit: `make test` discovers `src/tests/<name>/run.sh` in the
recipe, and both drivers live in a suite that already exists.

## 6. Noticed outside scope

- `hydrogen_helium_row_terms` in `constrained_chemical_equilibrium.f90`, the
  diagnostic that decomposes the H+ and He+ rows term by term, assigned a
  local `n_hei = sden(is_HeI_SI) + sden(is_HeITR)` that nothing in the
  routine read. It was the summed reservoir of the entry text, left behind
  when the decomposition was written, and it is removed: an unread
  assignment naming a reservoir the routine does not use is a statement
  about the physics that is no longer true anywhere else in the file.
- `docs/audit_20260905/plan_20260918_rev1_source_probe.f90` still does not
  compile against the current `cx_add_to_fvec` signature, as the D7a and D7c
  reports said. Unchanged here, for the reason they gave: it is a dated
  record of that review.

## 7. Housekeeping

The private builds `build_D7d`, `build_D7d_ctl` and the binaries
`EXHALE_D7d.x`, `EXHALE_D7d_ctl.x` were deleted at the end. Nothing under
`backup/regression/` was written, no golden was refreshed, and the tree's own
`EXHALE.x` and `build/` were never touched.
