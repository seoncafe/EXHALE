# charge_exchange_rows

Tests of item D7 of `docs/PLAN_20260918_rev2.md`: **the equation basis of the
helium row, who owns the cell's rate coefficients, and the same reaction set
seen by a solver whose unknowns are the stage densities themselves.**

The generic charge-exchange assembly `cx_add_to_fvec` writes its increments
in the ionization-positive form, +R on the donor's lower-stage boundary row
and -R on the acceptor's. Every row of every metal solver is written that way
except one: the summed He I balance of the He 2^3S systems
(`ion_residual_core::heh_tr_rows`, row 2) is written He I-gain positive, so a
reaction that destroys He II and makes He I has to enter it with the opposite
sign. The assembly therefore takes `he_row_sign`, which multiplies the
increments of the He I <-> He II boundary row alone. The reactions that reach
that row are the metal + He / He+ pairs of Huang et al. (2023), ApJ 951, 123,
Table 4 (group C, Si/C/O with He and He+), active only when `cx_full 1` is set
in `metals.inp`.

`cx_kc` and `cx_metal_base` are thread-local and hold one cell at a time.
`cx_set_cell(T)` fills them and records the temperature; `cx_add_to_fvec`,
`cx_add_to_jac`, `charge_exchange_stage_sources` and `cx_add_to_turnover` take
the cell temperature and stop the run when it is not the one the rates were
filled at, so no caller can assemble a cell on the rates of whatever cell its
thread handled last. The turnover's cell is optional in the interface while
two of its call sites still do not state one, and a call that states none is
refused when the thread holds no cell at all.

A transported ionization stage carries the stage DENSITY as its unknown, so
its charge-exchange source is the stoichiometric number-density source of that
stage and not a boundary flow. `charge_exchange_stage_sources` returns it, and
both it and `cx_add_to_fvec` read one rate vector from one reaction table, so
the two bases describe the same events; `stage_sources` asserts that equality
reaction by reaction and at a launch-region composition.

Run everything with

    src/tests/charge_exchange_rows/run.sh [reaction_sources] [stage_sources] \
                                          [pair_reservoir] [advection_pair] \
                                          [metal_helium] [stale_cell] \
                                          [stale_cell_stage] \
                                          [stale_cell_turnover] \
                                          [turnover_rates] [transport_refusal]

`EXHALE_OBJDIR` and `EXHALE_EXE` point the suite at another build (a private
`OBJDIR`/`EXE` pair); with either of them set the `make -q` staleness check is
skipped. Output is one `PASS|FAIL <name> measured=... reference=... tol=...`
line per assertion; the suite exits nonzero if any of them fails.

## reaction_sources (`reaction_source_orientation.f90`)

Links the production objects and assembles ONE reaction at a time, with the
two reactants at 1 cm^-3 and every other species absent, at T = 10000 K. The
rows are then converted to physical number-density sources by the conversion
of the basis they were assembled in, and the assertions are on those sources,
never on the rows.

| basis | solver | metal row base | `he_row_sign` | row 2 is |
|---|---|---|---|---|
| `standard` | `ion_system_HeH_metals` | 4 | +1 | the He I -> He II boundary flow |
| `triplet` | `ion_system_HeH_TR_metals` | 5 | -1 | the summed He I source |
| `molecular` | `ion_system_HeH_mol_metals` | 9 | +1 | the He II source |

The reactions, written from Table 4 and from the group-E source rather than
read from the module's descriptor tables, are the six group-C pairs
(`C1_Si_Hep`, `C2_Sip_He`, `C3_C_Hep`, `C4_Cp_He`, `C5_O_Hep`, `C6_Op_He`),
one group-A reaction (`A9_Si_Hp`) and the group-E electron capture
(`E1_O2p_H`, O2+ + H0 -> O+ + H+, Barragan et al. 2006).

| assertion (per reaction, per basis) | reference | tolerance |
|---|---|---|
| `<rxn>_rate` | a strictly positive rate | 0 |
| `<rxn>_<basis>_stoichiometry` | -R, +R on the donor's two stages and -R, +R on the acceptor's, and nothing on any other stage of any element, divided by R | 1e-13 |
| `<rxn>_<basis>_nucleus_balance` | the stages of each element sum to zero, divided by R | 1e-13 |
| `<rxn>_<basis>_charge_source` | sum over species of stage times source is zero, divided by R | 1e-13 |

and, once each,

| assertion | reference | tolerance |
|---|---|---|
| `group_C_absent_without_cx_full` | the metal + He group assembles nothing when `cx_full` is off | 0 |
| `group_E_active_without_cx_full` | the O2+ + H0 capture is assembled anyway, gated by its own scale factor | 0 |
| `group_E_scale_is_its_own_gate` | `cx_o2p_h_scale` = 1, the published rate | 0 |
| `group_B_He0_Hp_not_in_generic_set`, `group_B_Hep_H0_not_in_generic_set` | the He <-> H pair is excluded from the generic set in both directions, so it is counted once | 0 |
| `group_B_He0_Hp_ionizes_He`, `group_B_He0_Hp_neutralizes_H`, `group_B_He0_Hp_triplet_basis` | `he_h_cx_fvec` puts +k1 n(He I) n(H II) on the He row and -k1 n(He I) n(H II) on the H row, and flips the He row in the He I-gain basis | 1e-14 relative |
| `zero_abundance_<basis>` | every row is zero when every density is | 0 |
| `jacobian_vs_fd_<basis>` | `cx_add_to_jac` against a central difference of the rows `cx_add_to_fvec` assembles in the same basis, at n_H = 1e8, n_He = 1e7, x(H II) = 0.30, x(He II) = 0.20, x(He III) = 0.05 and every metal at 0.25/0.10; the rate is bilinear in densities linear in the unknowns, so the difference is exact to rounding | 1e-9 of the largest derivative |
| `cell_order_independence` | eight cells from 2000 to 16000 K, assembled in increasing and then in decreasing order, give the same rows both times | 0 |

It also prints, as a DIAGNOSTIC, what the triplet basis produced before the
orientation argument existed: the generic assembly asked for `he_row_sign` =
+1 there, on Si I + He II -> Si II + He I, gives a He II source of
+7.561e-11 cm^-3 s^-1 where the reaction requires -7.561e-11, and a total
ionic-charge source of +1.5122e-10 cm^-3 s^-1 where an exchange requires 0.

## pair_reservoir (`he_h_pair_reservoir.f90`)

Which helium state the He <-> H pair of Table 4 (group B) reacts from, asked
of the four systems that carry helium. The forward rate, from Glover &
Jappsen (2007), carries the barrier exp(-12.75/T4), and 12.75e4 K = 10.99 eV
is the ionization-potential difference 24.587 - 13.598 eV of GROUND-STATE
helium against hydrogen, so the reactant is He(1^1S) and not the sum over
He I. He(2^3S) lies 19.82 eV above the singlet and its ionization potential
is 4.77 eV, so the same collision is exothermic for it: a different reaction
with a different rate, carried nowhere in this code, and in particular not in
the metastable balance `ion_residual_core::tr_triplet_row`.

Each system's residual is evaluated twice at one state, with
`he_h_charge_exchange` true and false; the difference is the pair's
contribution to each row and nothing else. The state is n_H = 1e9 cm^-3,
n_He = 1e8 cm^-3, x(H II) = 0.30, x(He II) = 0.20, x(He III) = 0.05, at
T = 10000 K, with the canonical element list present and every element empty.
The comparison state adds a metastable of 10 percent of the helium nuclei
while holding n(H0), n(H+), n(He 1^1S), n(He+) and n(He2+) at the same
absolute values, which the fraction unknowns allow only by enlarging the
helium nucleus density by exactly the metastable density.

| assertion (per system) | reference | tolerance |
|---|---|---|
| `pair_He_source_singlet_<system>`, `pair_H_source_singlet_<system>` | the sources are k2 n(He+) n(H0) - k1 n(He 1^1S) n(H+) and its negative, with the ground singlet and not the summed He I | 1e-13 of the assembled rows |
| `pair_He_source_agrees_<system>_with_HeH_TR`, `pair_H_source_agrees_...` | the four systems make the same pair sources out of one state | 1e-13 of the assembled rows |
| `pair_He_source_no_metastable_channel_<system>`, `pair_H_source_...` | adding a metastable at fixed singlet moves neither source | 1e-13 of the assembled rows |
| `pair_jacobian_column<k>_<system>`, k = 1..4, triplet systems | the pair's derivative with respect to each unknown of the triplet layout against a central difference of the rows; column 4 is nonzero exactly because the reservoir is the singlet | 1e-5 relative |

## advection_pair (`advection_pair_reservoir.f90`)

The same question of `System_implicit_adv_HeH_TR`, the advection-corrected
system `post_process_adv` solves upwind across each cell, whose row 2 is the
ground-singlet balance and whose row 1 tracks the neutral hydrogen. The pair
is isolated the same way, by differencing the residual with
`he_h_charge_exchange` true and false, and the rows are converted to physical
sources by S(H I) = d(1) n_H/c1 and S(He I) = d(2) n_He/c1, with c1 = dr/v.

The state is the one `pair_reservoir` uses, at c1 = 1e4 s. The comparison
state adds a metastable of 10 percent of the helium nuclei while holding
n(H0), n(H+), n(He 1^1S), n(He+) and n(He2+) at the same absolute values.

| assertion | reference | tolerance |
|---|---|---|
| `adv_pair_He_source_singlet`, `adv_pair_H_source_singlet` | the sources are k2 n(He+) n(H0) - k1 n(He 1^1S) n(H+) and its negative, with the ground singlet and not the summed He I | 1e-11 of the assembled rows |
| `adv_pair_He_source_no_metastable_channel`, `adv_pair_H_source_...` | adding a metastable at fixed singlet moves neither source | 1e-11 of the assembled rows |
| `adv_pair_He_source_agrees_with_HeH_TR`, `adv_pair_H_source_...` | the advection-corrected and the equilibrium triplet system make the same pair sources out of one physical state, each through its own row basis and normalization | 1e-11 of the assembled rows |

## metal_helium (`metal_helium_reservoir.f90`)

Which helium state the metal + He reactions of Table 4 group C react from,
asked of the three systems that assemble them: `System_HeH_TR_metals`,
`System_HeH_mol_metals` and the balance rows of
`constrained_chemical_equilibrium`. Those rates are ground-state rates: the
three reverse rows Si+ + He, C+ + He and O+ + He carry exp(-19.1/T4),
exp(-15.5/T4) and exp(-12.7/T4), and those are the ionization-potential
differences of ground-state helium, 24.587 - 8.152 = 16.436 eV = 19.07e4 K,
24.587 - 11.260 = 13.327 eV = 15.47e4 K and
24.587 - 13.618 = 10.969 eV = 12.73e4 K.

Group C is active only with `cx_full`, so each system's residual is evaluated
with `cx_full` true and false and the difference taken. Oxygen is the only
element carrying atoms, so every group D row has a reactant at zero, and
groups A and E are active in both evaluations and cancel.

The state isolates C6, the one group C row whose reactant is the helium: all
the oxygen is ionized (x(O II) = 0.95, x(O III) = 0.05, n_O = 1e5 cm^-3), so
C5 has no neutral oxygen to start from; the hydrogen is ionized down to 1 atom
cm^-3 of H0; every photoionization, recombination and collisional coefficient
of the cell is zero and the He <-> H pair is off, so the rows carry charge
exchange alone. The helium is n_He = 1e8 cm^-3 with x(He II) = 0.20 and
x(He III) = 0.05 at T = 10000 K, and the comparison state adds a metastable of
10 percent of the helium nuclei at fixed singlet, He+ and He2+ densities.

| assertion (per system) | reference | tolerance |
|---|---|---|
| `groupC_He_source_no_metastable_channel_<system>`, `groupC_O_source_...` | adding the metastable at fixed singlet moves neither source | 1e-9 of the group C source |
| `groupC_He_source_agrees_<system>_with_HeH_TR_metals`, `groupC_O_source_...` | the three systems make the same group C sources out of one state | 1e-9 of the group C source |
| `groupC_charge_balance_<system>` | the helium and oxygen sources sum to zero: C6 moves one electron between the two elements and nothing else | 1e-9 of the group C source |

## stale_cell (`stale_cell_rates_refused.f90`)

Loads one cell, assembles it, loads a second cell, assembles it, then
assembles the FIRST cell again without reloading its rates. The program is
expected to stop with a nonzero status and the `(charge_exchange) ERROR`
message; a clean exit is a failing verdict.

| assertion | reference | tolerance |
|---|---|---|
| `stale_cell_rates_refused` | nonzero exit and the named refusal | 0 |

## stage_sources (`transported_stage_sources.f90`)

The same reaction set in the basis a transported stage is written in. Each
reaction is assembled once through `cx_add_to_fvec` in each of the three
bases above and once through `charge_exchange_stage_sources`; the rows are
converted to physical sources by the basis's own conversion and the two are
compared stage by stage. The reaction list adds one group-D row (`D1_C_Sip`,
C + Si+ -> C+ + Si), which touches neither hydrogen nor helium, to the eight
of `reaction_sources`.

| assertion | reference | tolerance |
|---|---|---|
| `<rxn>_stage_rate` | a strictly positive rate read off the donor's own stage | 0 |
| `<rxn>_stoichiometry` | -R, +R on the donor's two stages and -R, +R on the acceptor's, divided by R | 1e-13 |
| `<rxn>_nucleus_balance` | the stages of each element sum to zero, divided by R | 1e-13 |
| `<rxn>_charge_source` | sum over species of stage times source is zero, divided by R | 1e-13 |
| `<rxn>_gross_split` | gross production minus gross loss is the net source of every stage, divided by R | 1e-13 |
| `<rxn>_gross_nonnegative` | neither gross magnitude is negative | 0 |
| `<rxn>_<basis>_stage_sources_match` | the stage sources equal the sources the basis's rows represent, divided by R | 1e-13 |
| `production_state_<basis>_stage_sources_match` | the same with the whole table active at a launch-region composition, relative to the largest stage source | 1e-13 |
| `group_C_absent_without_cx_full`, `group_D_absent_without_cx_full` | exactly zero | 0 |
| `group_E_active_without_cx_full`, `group_E_off_at_zero_scale` | the group-E capture follows `cx_O2p_H` and not `cx_full` | 0 |
| `group_B_He0_Hp_not_in_stage_sources`, `group_B_Hep_H0_not_in_stage_sources` | exactly zero: the He <-> H pair is applied by `he_h_cx_fvec` alone and is counted once | 0 |
| `cell_order_independence` | eight cells from 2000 to 16000 K, evaluated in increasing and then in decreasing order, bitwise identical | 0 |
| `thread_count_independence` | the same eight cells on one thread and on four, bitwise identical | 0 |

## stale_cell_stage (`stale_cell_stage_sources_refused.f90`)

`stale_cell`'s assertion for the stage-source interface: two cells are loaded
and evaluated, then the first is evaluated again without reloading its rates.

| assertion | reference | tolerance |
|---|---|---|
| `stale_cell_stage_sources_refused` | nonzero exit and `(charge_exchange) ERROR` in the output | 0 |

## stale_cell_turnover (`stale_cell_turnover_refused.f90`)

`stale_cell`'s assertion for the acceptance normalization: the turnover bound
of a row is `kc(T)` times two element counts, so it belongs to the cell whose
temperature filled the rates. Two cells are bounded, then the first is bounded
again without reloading its rates.

| assertion | reference | tolerance |
|---|---|---|
| `stale_cell_turnover_refused` | nonzero exit and `(charge_exchange) ERROR` in the output | 0 |

## turnover_rates (`turnover_cell_rates.f90`)

That the turnover bound is the loaded cell's and carries no history, on one
thread and on four. The two cells are 6000 K and 10000 K, and they are
alternated.

| assertion | reference | tolerance |
|---|---|---|
| `turnover_moves_with_the_cell_temperature` | strictly positive | 0 |
| `turnover_of_a_reloaded_cell_is_its_own` | 0 | 0 |
| `turnover_on_several_threads_is_each_cells` | 0 | 0 |

The first row is the non-vacuity guard of the other two: a bound that did not
move with the temperature would make them hold for any implementation.

## transport_refusal (`ionization_transport_with_metals.sh`)

Runs `examples/14_diffusion` plus `Ionization transport: True` twice, capped
at three steps: once with the example's `metals.inp` and once without it. The
transported stage sources now carry the metal charge exchange of Table 4
through `charge_exchange_stage_sources`, the same reaction set and the same
rate coefficients the local sweep's rows take through `cx_add_to_fvec`, so
both mixtures run; while those terms were missing the metal-bearing one was
refused by name, and this fixture is what that refusal was asserted by.

| assertion | reference | tolerance |
|---|---|---|
| `ionization_transport_with_metals_runs` | exit 0 and no `(input_read) ERROR` | 0 |
| `ionization_transport_metal_free_runs` | the same key without metals runs and exits 0 | 0 |
