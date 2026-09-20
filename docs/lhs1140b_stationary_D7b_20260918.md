# D7b: the metal charge exchange of a transported ionization stage

Item D7b of `docs/PLAN_20260918_rev2.md`. The transported ionization stages
took their source from the H/He rows of `ion_residual_core` alone, while the
local sweep of the same cell adds the metal charge exchange of Huang et al.
(2023), ApJ 951, 123, Table 4 to those same rows. With metals in the mixture
the stage the flow carried and the composition the sweep returned answered two
different H+ and He+ balances of one cell, and the combination was refused
(item D7a). The stage sources now carry those reactions and the refusal is
lifted.

## Verdict

The stage sources take the metal terms through one reaction-source interface,
`charge_exchange::charge_exchange_stage_sources`, which returns the
stoichiometric number-density source of each hydrogen, helium and metal stage.
It and `cx_add_to_fvec` read one rate vector built from one reaction table, so
there is one spelling of every rate and of the stoichiometry, and the two
solvers differ only in the equation basis they map it into. For the transport
operator that map is the identity: its unknowns ARE the stage densities.

MEASURED on a metal-bearing cell at 8000 K with the whole table active, the
term that was missing is 1.81e-3 of the H II source and 2.38e-6 of the He II
source; after the change the two agree to 2.2e-16 and to 0 exactly. On the
solved 500-cell column of the key-on `wasp_full` run the transported stage
source and the sweep's row-implied source of the metal charge exchange are
equal at EVERY cell, difference exactly 0.

## 1. Which reactions reach a transported stage

READ from the header and the descriptor tables of
`src/modules/radiation/charge_exchange.f90`, and from
`references/Huang_2023_ApJ_951_123.pdf` Table 4 (page 25).

| group | reactions | reaches | active when |
|---|---|---|---|
| A | metal + H and metal + H+, 23 rows | the H I <-> H II boundary | metals present (the default set) |
| B | He + H+ and He+ + H | H and He | applied by `he_h_cx_fvec`, excluded from the generic set, counted once |
| C | metal + He and metal + He+, 6 rows (Si, C, O) | the He I <-> He II boundary | `cx_full 1` in `metals.inp` |
| D | metal + metal, 32 rows | no H or He stage | `cx_full 1` |
| E | O2+ + H0 -> O+ + H+ (Barragan et al. 2006), not in Table 4 | the H I <-> H II boundary | `cx_O2p_H` scale, positive by default |

None of group C touches He2+: `cx_acc_stg` is 1 and `cx_don_stg` 0 for every
one of its rows, so `cx_fvidx(cx_He, 0)` is the only helium row they reach
(READ). Group A is on in every metal-bearing run, which is why the refusal was
not conditional on `cx_full`.

## 2. The interface, and why it is shaped this way

A local ionization system carries BOUNDARY FLOWS as its rows: row 1 is the net
upward flow across the H I <-> H II boundary, rows 2 and 3 the two helium
boundaries (with the one orientation difference item D7a fixed), and the metals
follow from `cx_metal_base`. A transported stage carries the STAGE DENSITY as
its unknown, so what it needs is the stoichiometric source of that stage. The
two are related exactly: for an element whose boundary flows are u_0 and u_1,

```
   S_0 = -u_0 ,   S_1 = u_0 - u_1 ,   S_2 = u_1 .
```

`charge_exchange_stage_sources` returns S directly, with the optional gross
production and loss of each hydrogen and helium stage beside it, and
`cx_add_to_fvec` assembles u as it always did. Both call `cx_reaction_rates`,
which is the one place a rate `R = k n(donor) n(acceptor)` is formed, and both
address the reaction table `cx_don_el/cx_don_stg/cx_acc_el/cx_acc_stg`, which
is the one place the stoichiometry is written.

**Why `cx_add_to_fvec` is not written as a caller of the stage-source
routine.** The plan allowed either shape. Making the residual assembly build
stage sources first and convert them into its rows afterwards changes the
ORDER of the floating-point additions into `fvec`: the entry text adds each
reaction's rate to the row it is already accumulating, while a conversion adds
a completed difference. The two differ by rounding in the rows of every
metal-bearing key-off run, so byte identity on `wasp_full` and `mol_metals`
would have been given up for a code shape, not for physics. The shared object
is therefore the rate vector and the reaction table, which is what "one
evaluator" protects; that the two bases then describe the same events is not
assumed but asserted reaction by reaction in
`src/tests/charge_exchange_rows/transported_stage_sources.f90`, and measured on
a whole solved column (sections 6 and 7).

The routine refuses an evaluation whose cell is not the one this thread's rate
coefficients were filled at, through the same `cx_require_cell` guard item D7a
added, and `carrier_source` calls `cx_set_cell` for the cell it is assembling
before it. `cx_kc` is thread-local, `carrier_source` runs inside an OpenMP
loop, and a source that inherited a worker's previous cell would depend on the
schedule; the suite measures eight cells from 2000 to 16000 K in increasing and
in decreasing order and on one thread against four, bitwise identical in all
four arrangements.

## 3. Where the metal stage densities come from

A metal stage is not a transported unknown: this operator moves the hydrogen
and helium stages and the molecular carriers, and the ionization sweep that
follows re-solves the metal partition from its own balance. The metal stage
densities are therefore part of the frozen background of one transport step,
exactly as the molecular ions and the electrons are, and they are filled in
`carrier_state` from the same `f_sp` the rest of the background comes from:

```
   cbg_nm0(:,e) = f_sp(:, mion_fsp(melem_i0(e)))     * nd
   cbg_nm1(:,e) = f_sp(:, mion_fsp(melem_i0(e)+1))   * nd
   cbg_nm2(:,e) = f_sp(:, mion_fsp(melem_i0(e)+2))   * nd   (melem_top >= 2)
```

with `nd = rho n0`. They join `save_carrier_module_state` and
`restore_carrier_module_state`, so measuring a carrier balance still leaves the
module in the state it was in.

The helium reactant of group C is the ground singlet He(1^1S), the reservoir
all three metal systems pass since item D7d, so `carrier_source` passes
`n_heiSI` and not the summed He I.

## 4. The He <-> H pair of the atomic branch, corrected with it

`carrier_source`'s atomic branch passed the SUMMED He I to `he_h_cx_fvec`
while `System_HeH_TR`, whose rows it is evaluating, passes the ground singlet
(item D7c). That is the same defect D7c fixed in the local systems, left in the
transport operator because that file was not D7c's. It is fixed here: the two
now react the pair from one reservoir, which the acceptance of this item
requires anyway, since the stage sources and the sweep's rows must be one
expression. The Table 4 rate of He + H+ carries the barrier exp(-12.75/T4), and
12.75e4 K = 10.99 eV is the ionization-potential difference 24.587 - 13.598 eV
of GROUND-state helium against hydrogen; He(2^3S) lies 19.82 eV above the
singlet and the same collision is exothermic for it, a different reaction that
this code carries nowhere.

No regression case reaches that branch: the four cases that set `Ionization
transport: True` (`carrier_model_a_newton`, `hp_trace_seed`, `hp_zero_seed`,
`hp_front`) are all `Molecular chemistry: True` and take the molecular branch,
which already passed the singlet (READ from their `input.inp`).

## 5. The channel ledger

The 26-channel layout of item D3a gains six channels, the gross production and
the gross loss of H II, He II and He III under the metal charge exchange
(`ich_metal_*`, 27 to 32). They are filled by whoever applies the reactions,
exactly as the four channels of the He <-> H pair are, and
`charge_exchange_stage_sources` hands out those magnitudes directly rather than
having them rebuilt. The He III pair carries the zero its own gross rates
return, because no reaction of the present set captures an electron onto He III
or promotes He II; it exists so that the ledger of that stage stays complete if
the set grows.

With the channels present, `stage_source_from_channels(chan, irow)` reproduces
the assembled source of a metal-bearing cell in both gas branches to 0.0
exactly (MEASURED, tolerance 8 eps = 1.78e-15).

The atomic kernel also gains the call-site argument D3a's item 3 named.
`mol_heh_rows` can tell the transport operator's evaluation from the local
sweep's because the operator is its only caller that asks for the
production/loss split; `heh_tr_rows` is asked for nothing that would
distinguish them, so it wrote every record `A-`. It now takes an optional
`transport_operator`, which `carrier_source` sets, and the record carries `AC`
for the operator's evaluation and `AS` for the sweep's, matching the `MC`/`MS`
of the molecular branch that the reader of the record already keys on. The
rows themselves do not read it.

## 6. What was measured

All numbers MEASURED on 2026-09-18 with a private build of the delivered tree
(`EXHALE_D7b.x`, md5 `1abf66c7024f0c2a066814e9749f2f21`) against a private
build of this item's entry text (`EXHALE_D7b_ctl.x`, md5
`7bcd4f59ded4c6e069f6ce897ad652e3`), both from the same `make`, both
single-threaded for the byte-identity runs.

### The term that was missing

One cell at 8000 K, n(H I) = 6e8, n(H II) = 4e8, n(He 1^1S) = 7e7,
n(He II) = 2e7, n(He III) = 1e6 cm^-3, every metal present at solar order and
spread 0.55 / 0.35 / 0.10 over its three stages, `cx_full` on.

| quantity [cm^-3 s^-1] | entry text | delivered |
|---|---:|---:|
| S(H II), transported | 1.707322700000000e+07 | 1.710423815876733e+07 |
| S(H II), sweep row | 1.710423815876733e+07 | 1.710423815876733e+07 |
| gap, relative | **1.813e-03** | 2.178e-16 |
| S(He II), transported | -1.197031965000000e+07 | -1.197034813802725e+07 |
| S(He II), sweep row | -1.197034813802725e+07 | -1.197034813802725e+07 |
| gap, relative | **2.380e-06** | 0 exactly |

### The solved column

`solved_state_sources.x`, which reads a run's `Ion_species.txt` and
`Hydro_ioniz.txt` and evaluates both bases at every physical cell: the maximum
absolute difference over the 500 cells of the key-off `wasp_full` state is
0.000e+00 for H II (scale 1.380e+05 cm^-3 s^-1) and 0.000e+00 for He II (scale
5.082e+02), and the same on the key-on solved state (section 7).

### Key-off byte identity

Control and measured are two private builds of the two source texts, the
runs single-threaded on scratch copies, the data lines compared with the
provenance headers stripped (those carry a timestamp).

| case | steps | files | verdict |
|---|---:|---:|---|
| `wasp_full` (He 2^3S + metals, the triplet system with `cx_add_to_fvec`) | 400 | 4 | every data line byte-identical; log10(Mdot) 7.20 g/s both ways |
| `mol_metals` (molecular chemistry + metals) | 12000 | 9 | every data line byte-identical; log10(Mdot) 10.58 g/s both ways |
| `mol_carrier` (molecular carrier transport, ionization key off) | 12000 | 7 | every data line byte-identical; log10(Mdot) 10.57 g/s both ways |

### Bounds-checked

A build with `-O1 -fcheck=bounds,do,mem` (`EXHALE_D7b_fc.x`): the two suites
run clean under it (250 and 10 assertions, 0 failures), and a metal-bearing
key-on run (`examples/14_diffusion` plus `Ionization transport: True`, 60
steps) finishes with no runtime trap.

### The suites that reach the changed modules

MEASURED with the delivered build, 0 failures in every one:
`charge_exchange_rows` 250, `stage_row_balance` 10, `carrier_retry` 156,
`steady_species_rows` 197, `ionization_imposed_fractions` 48,
`ionization_stage_flux` 45, `element_operator` 40,
`carrier_returned_state_acceptance` 36, `coupled_source_step` 30,
`carrier_reference_scales` 16, `carrier_outer_boundary` 16,
`carrier_helium_inventory` 15, `carrier_boundary_jacobian` 12,
`carrier_constraint_attribution` 10, `species_masses` 9,
`constrained_network_layout` 8, and `grid_and_gates` 284, which carries the
atomic ionization-transport fixture.

## 7. The key-on solve

`wasp_full` with `Ionization transport: True`, seeded by its own key-off
400-step state (`Load IC? True`, `Restart option change: iontrans`, the loader
noting that the H+ column of the seed is a local equilibrium) and run through
the partitioned stationary route: `Restart intent: stationary`,
`Solver: Newton 100.0`, `Resid tol: 1.0e-8`, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=90`, 8 threads on lart4. The run starts: the refusal is
lifted and the metal-bearing transported configuration is accepted.

**It does not certify.** The alternation is refused at outer pass 6 because the
joint distance has not fallen in three consecutive passes, and the state is
written `certified=F`. Pass history (MEASURED; the pass allowance of the
`transported_ionization` class is 90, D9):

| pass | hydro info | mass | momentum | energy | worst gated species row | s |
|---:|---:|---:|---:|---:|---|---:|
| 1 | 2 | 1.16 | 1.97 | 1.96 | 1.00 at cell 398, carrier balance He+ | 65.0 |
| 2 | 2 | 1.57 | 1.98 | 1.97 | 1.00 at cell 401, carrier balance He+ | 82 |
| 3 to 6 | 2 | 1.57 | 1.98 | 1.97 | the same | 80 to 24 |

The seven entries that refuse the written state: the hydrodynamic mass row
(9.998e-1 against 3.0e-12 at cell 400), momentum (1.980 against 1.0e-8 at 364),
energy (1.975 against 1.0e-6 at 110), the three carrier balances H+, He+ and
He++ (9.998e-1, 1.000, 9.990e-1 against 1.0e-5, all at cell 401, a wind cell)
and the H(n=2) level balance (1.000 against 1.0e-6 at cell 395).

**The refusal is the case, not the metal charge exchange.** The same
configuration with `metals.inp` removed, seeded by its own metal-free 400-step
state, refuses at outer pass 5 with the SAME seven entries and the same
saturated measures (mass 1.197 at cell 396, momentum 1.954 at 488, energy 1.974
at 138, the three carrier balances at 1.000 at cell 396, H(n=2) 1.000 at cell
393) (MEASURED). That run reaches none of the code this item added
(`thereis_metals` is false), so the alternation's failure to approach a joint
fixed point on this case is independent of it. `wasp_full` is a 400-step
relaxation snapshot of WASP-121b whose hydrodynamic rows stand at order unity
from the start; it was never a stationary solution, and restarting one with the
stationary intent does not make it one.

**What the run does establish**, which is what this item is accountable for:
on the state it wrote, the transported stage source and the source the local
sweep's rows imply for the metal charge exchange are equal at EVERY one of the
500 physical cells, difference exactly 0.000e+00 against scales of 1.353e+05
(H II) and 5.174e+02 (He II) cm^-3 s^-1 (MEASURED). The accepted flux gate of
the written state is 1.6414e+01.
