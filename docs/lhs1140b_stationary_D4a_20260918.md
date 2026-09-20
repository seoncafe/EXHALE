# D4a: what the thermochemical closure hands back, measured

Item D4a of `docs/PLAN_20260918_rev2.md` (the measurement; D4b is the
decision and is not taken here). 2026-09-18 (KST). Every number below is
MEASURED with the driver named in section 1 unless it is marked READ, in
which case the source line it is read from is given. No production source
was changed, no tolerance was moved, and no golden was touched.

The question D4 leaves open is which of two things the failing
`carrier_retry` row means: production's closure exits before its promise is
met, or the test assumes a bound the residual contract does not give. The
measurements answer it; section 8 states the answer with the numbers.

---

## 1. Provenance

| item | value |
|---|---|
| source snapshot | `Makefile` + `src/` of the working tree copied 2026-09-18 14:07 UTC to the session scratch directory `D4a/tree/` |
| repository HEAD at the snapshot | `3c73905ca8a7fe2af92a5c2b014c225b2feedff3` |
| working tree against HEAD | 89 files changed, 12171 insertions, 1130 deletions (`git diff --stat`) |
| compiler | GNU Fortran (conda-forge gcc 16.2.0-4) 16.2.0, `/opt/miniconda3/bin/gfortran` |
| production link, `make` in the snapshot | `-O3 -fopenmp`, `-L/opt/miniconda3/lib -lopenblas -ldl`; `EXHALE.x` md5 `44afa4848d143bd6570f7be86ee72e88` |
| source manifest (md5 of the md5 list of every `.f90`, `.f`, `.sh`, `.py` under `src/` plus `Makefile`) | `f3ca58742e0027ed07ed2d6713d29051`, taken after the build and re-taken at the end of the work; **identical**, so nothing the measurement used moved under it |
| measurement driver | `src/tests/carrier_retry/closure_returned_state.f90`, md5 `c9439efbe8a220f8d926ea7a2ffecaeb`; runner `src/tests/carrier_retry/run_closure_returned_state.sh`, md5 `c55ce007facadf97b86836bb8109025e` |
| driver build | the same compiler, `-O0 -g -fcheck=all -fbacktrace -fopenmp`, linking the snapshot's production sources through `src/tests/physics_probe/source_closure.py` plus MINPACK and OpenBLAS; executable md5 `b149cbafb91ce58396782791dba3cf98` |
| threads | `OMP_NUM_THREADS=1` |
| saved return | `D4a_returned_state.txt` in the object directory of `run_closure_returned_state.sh` (`build/tests/closure_returned_state/` by default) |

The existing suite `src/tests/carrier_retry/run.sh` was also built and run
unmodified on the same snapshot, as the anchor of section 2. The driver is
a NEW program; `carrier_retry.f90` was not edited.

Every source line number quoted below is the snapshot's. Those in
`ionization_equilibrium.f90`, `caloric_eos.f90`, `ion_residual_core.f90` and
the two tolerance lines of `diffusive_photochemistry.f90` were checked
against the working tree at the end of the work and still hold; the rest of
`diffusive_photochemistry.f90` is under edit by another item and its
numbering will move.

Tolerances, READ from source:

| name | value | where |
|---|---|---|
| `ieq_res_tol` | `1.0d-6` | `ionization_equilibrium.f90:509` |
| `chem_cycle_tol` | `1.0d-5` | `diffusive_photochemistry.f90:829` |
| `chem_cycles_max` | `12` | `diffusive_photochemistry.f90:828` |
| the fixture's presence rule | `f_sp > 1.0d-20` | `carrier_retry.f90`, the section-11 row |

---

## 2. The state, and how closely this driver reproduces the suite's

The column is the twelve-cell molecular hydrogen and helium column of
`carrier_retry.f90`: 80 per cent of the hydrogen nuclei in H2, helium at
`HeH = 0.0793`, no metal abundance, `Ionization transport` on, the
molecular chemistry on, the carrier transport on, the wind `v = 20 r`, and
one relaxation pass at trust `1e-2` entered from the same checkpoint (one
transport interval taken from the entry composition).

The unmodified suite reports, on this snapshot (MEASURED):

```
one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance
    measured = 1.100521341448987E-06   tol = 1.000000000000000E-06   FAIL
the largest further movement is f_sp column 37 at cell 1,
    from and to 1.328597E-12 1.328596E-12
```

The driver reports `1.110026074407983E-06`, the same species (column 37,
HeH+) at the same cell (1). The 0.9 per cent difference is accounted for
and is not noise: the suite's section 11 re-seeds the mass row, the frozen
background, the composition and the carrier checkpoint but NOT the
pressure, so `conserved_of` builds the conserved state from whatever `p_col`
the program's earlier passes left. Repeating the sequence with the pressure
each repetition returns gives (MEASURED):

| repetition | `p(1)` entering | one further sweep moves | species, cell |
|---|---|---|---|
| 1 | 7.270331132884705e-1 | 1.104145377496377e-6 | HeH+, 1 |
| 2 | 7.329694837946903e-1 | 1.094853731742269e-6 | HeH+, 1 |
| 3 | 7.389498536586548e-1 | 1.293075021802678e-6 | HeH+, 1 |
| 4 | 7.449759681774529e-1 | 1.083573868280137e-6 | HeH+, 1 |

So the row's value carries a +-10 per cent dependence on the run history
the sequence is entered with, and the suite's 1.1005e-6 and this driver's
1.1100e-6 both sit inside that spread. **What does not vary: the row is
above `ieq_res_tol` in every repetition, and the binding species and cell
are HeH+ at cell 1 in every one.**

---

## 3. Table 1: the production return, and the module state beside it

The composition, pressure, temperature, electron and total particle
density, the imposed fractions with their flags, the radiation field and
the rate state of every physical cell are written to
`D4a_returned_state.txt` (128 lines, `ES23.15E3`). The binding cell:

| quantity | value |
|---|---|
| `r` | 1.020000000000000 |
| `p` [code] | 7.211569230996558e-1 |
| `T` returned | 1.497358208954703e+3 K |
| `n_e` | 1.145524008649463e+5 cm^-3 |
| `n_tot` | 4.816080522277298e+9 cm^-3 |
| `x_h2_fixed`, `x_h2_fix` | true, 7.865001382104033e-1 |
| `x_hp_fixed`, `x_hp_fix` | true, 6.952990865123550e-6 |
| `P_HI`, `P_HeI`, `P_H2` [1/s] | 3.859070138347077e-8, 6.859510884395538e-7, 1.019690636596275e-7 |
| `k_LW` [1/s] | 0 |
| `rchiiB` [cm^3/s] | 1.113552969549173e-12 |

The pass ends on the movement bound with the last step inside it, after 4
kept steps; the closure of its last trial reports `converged` with a last
temperature increment of 2.368608466342650e-7.

**The module state the next real operation reads**, READ from the source
and listed because measurement 4 is about whether a refused trial puts it
back:

| item | module | read by |
|---|---|---|
| `bg_cell(1-Ng:N+Ng)` (rates, `T_K`, `n_tot`, `x_*_fixed`/`x_*_fix`, `nh`, `nhe`) | `ionization_equilibrium` | the carrier transport operator (`carrier_source`) and the seeds of the next sweep |
| `ieq_rate_cell`, `ieq_ne_cell`, `ieq_TK_cell`, `ieq_ntot_cell`, `ieq_met_coef`, `ieq_neq_stored`, `ieq_mbase_stored`, `ieq_iox_stored`, `ieq_rates_ready` | `ionization_equilibrium` (lines 307 to 318, written at 2138 to 2150) | `ionization_closure_residual_cell` and `ionization_closure_residual_profile`, i.e. the certification's chemical residual |
| `ieq_nonroot_streak` | `ionization_equilibrium:551` | the non-root persistence stop |
| `nk_per_mass`, `x_h2`, `molecular_cell` (private) | `caloric_eos:159-161` | `pressure_from_energy_density` and `energy_density_from_pressure`, so `Apply_BC`, the reconstruction and the residual assembly; refreshed only by `get_species_densities` |
| `dp_bc`, `n_part_cell1` | `global_parameters` (written inside the sweep, `ionization_equilibrium.f90:3167`) | the base ghost of `Apply_BC` |
| the carrier checkpoint (composition, carrier fractions, history, attempt statistics) | `diffusive_photochemistry` | the next transport interval |

---

## 4. Table 2: the chemical residual of the returned composition

Two evaluations, both with `ionization_closure_residual_profile`, which
solves nothing and writes nothing:

- **res (stale rates)**: the residual of the returned composition against
  the rate state left in `ieq_rate_cell` when the pass returned.
- **res (at T returned)**: the residual of the same composition against a
  rate state rebuilt at the RETURNED temperature from the returned
  composition (one sweep at `T_ret` on a copy, which fills `ieq_rate_cell`
  and `ieq_ne_cell` from `(T_ret, f_ret)` before it solves anything), the
  composition itself untouched.

| cell | res (stale rates) | res (at T returned) | `ieq_TK_cell` [K] | `bg_cell%T_K` [K] | T returned [K] |
|---|---|---|---|---|---|
| 1 | 5.332326629e-4 | 1.610675168e-16 | 1497.252852 | 1497.358574 | 1497.358209 |
| 2 | 3.075043380e-5 | 3.948886135e-17 | 1500.196155 | 1500.190124 | 1500.190127 |
| 3 | 2.106969082e-5 | 5.932800158e-17 | 1500.127846 | 1500.123712 | 1500.123715 |
| 4 | 1.905834360e-5 | 9.130940273e-17 | 1500.125222 | 1500.121486 | 1500.121488 |
| 5 | 1.986009216e-5 | 1.438365547e-16 | 1500.134273 | 1500.130380 | 1500.130383 |
| 6 | 2.063189524e-5 | 2.320921613e-16 | 1500.144300 | 1500.140256 | 1500.140259 |
| 7 | 2.136700667e-5 | 4.135943490e-16 | 1500.155744 | 1500.151556 | 1500.151561 |
| 8 | 2.207624332e-5 | 1.134245472e-15 | 1500.169333 | 1500.165008 | 1500.165013 |
| 9 | 2.276074582e-5 | 1.296143055e-15 | 1500.186331 | 1500.181875 | 1500.181881 |
| 10 | 2.342173942e-5 | 2.849548096e-15 | 1500.209370 | 1500.204788 | 1500.204795 |
| 11 | 2.392681726e-5 | 7.718291970e-15 | 1500.244118 | 1500.239445 | 1500.239454 |
| 12 | 4.076411851e-5 | 2.592940096e-14 | 1500.391369 | 1500.383404 | 1500.383434 |

Binding cell 1; largest residual of any cell at the returned temperature
2.592940096e-14, which is **eight decades below `ieq_res_tol`**.

Two readings follow.

1. **The returned composition IS a root of the network at the temperature
   it is handed back with**, to 2.6e-14 over the column and 1.6e-16 at the
   binding cell. `bg_cell%T_K`, which the refusal branch does restore and
   which is the temperature the accepted trial's sweep ran at, stands
   2.4e-7 relative from the returned temperature at the binding cell (the
   closure's last increment, 2.37e-7).
2. **`ieq_rate_cell` does not describe the returned composition.** It stands
   7.0e-5 relative away in temperature at cell 1. The reason is measurement
   4: the pass ends by REFUSING a trial on the movement bound, the refusal
   branch of `relax_photochemical_composition` restores `f_sp`, `bg_cell`,
   `p`, `T`, `heat`, `cool` and `eta`, but not the rate state the closure
   evaluator reads. The 5.3e-4 in the first column is therefore the residual
   of one state measured against another state's rates, and it is
   the quantity the certification would read if it were asked for this
   state's chemical residual right now.

### Table 2b: row by row at the binding cell, at the returned temperature

| row | species | `x` | normalized residual |
|---|---|---|---|
| 1 | H II | 6.95299087e-6 | 0 |
| 2 | He II | 2.55521375e-7 | 0 |
| 3 | He III | 6.93675220e-11 | 0 |
| 4 | H2 | 7.86500138e-1 | 0 |
| 5 | H2+ | 9.82731724e-9 | 9.54659707e-18 |
| 6 | H3+ | 2.80180424e-5 | 1.61067517e-16 |
| 7 | HeH+ | 1.67948624e-12 | 5.98727848e-19 |
| 8 to 27 | the metal stages | 0 | 0 |

Rows 1 to 4 are exactly zero because they are identity rows: with
`Ionization transport` on, `impose_transported_ionization_fractions`
(`ion_residual_core.f90:260-274`) replaces rows 1, 2 and 3 by
`x - x_fix`, and the transported H2 carrier replaces row 4 the same way.
Only three genuine balances remain at this cell, H2+, H3+ and HeH+, and
**the HeH+ row of the state the test fails on carries a residual of
5.99e-19, thirteen decades below `ieq_res_tol`.**

---

## 5. Table 3: complete fixed-conserved-energy closure increments

One increment is one cycle of the production routine (a chemistry sweep AND
the equation-of-state temperature update at the unchanged conserved
variables), taken by `equilibrate_chemistry_at_fixed_conserved_state`
itself with `chem_cycles_cap_for_test = 1`, on an isolated copy.

| increment | thermal `dT/T` | largest chemical movement | where | H2+ | H3+ | HeH+ | H I |
|---|---|---|---|---|---|---|---|
| 1 | 1.2574e-11 | 1.1100e-6 | HeH+, cell 1 | 1.5842e-8 | 1.0743e-7 | 1.1100e-6 | 1.4101e-11 |
| 2 | 6.3637e-15 | 1.5724e-10 | HeH+, cell 12 | 1.0146e-11 | 5.2804e-12 | 5.7890e-12 | 8.5487e-16 |
| 3 | 4.4405e-16 | 5.5730e-13 | H3+, cell 1 | 2.5944e-13 | 5.5730e-13 | 1.8167e-13 | 3.4195e-16 |

(The last four columns are relative movements at the binding cell.)

The thermal part of the first increment is 1.26e-11, six decades below
`chem_cycle_tol`; the chemical part is 1.11e-6. The state settles by about
four decades a step, which is a contraction and not arithmetic noise, and
confirms L36d's third-sweep reading (1.0085e-10 there, 1.5724e-10 here).

The chain, from the first increment at the binding cell (MEASURED):
H2+ 1.584e-8 -> H3+ 1.074e-7 is a factor 6.78, H3+ 1.074e-7 -> HeH+
1.110e-6 a factor 10.33. This is L36d's "factor 7 to 10 per link",
measured here on the complete closure increment rather than on a chemistry
sweep alone.

---

## 6. Table 4: what a refused trial puts back

A refusal was forced through the closure's own cycle cap
(`chem_cycles_cap_for_test = 1`): the closure then spends its budget on
every trial, returns `ok = .false.` with the reason `cycle budget spent`,
and the pass refuses each trial on the chemistry and ends on
`carrier_relax_chemistry_refused` with 0 kept steps (both asserted, both
PASS). **A nonfinite or inadmissible trial state cannot be forced from
outside**: the public test hooks of `diffusive_photochemistry` are
`carrier_headroom_set_for_test`, `carrier_co_domain_perturb_for_test`,
`carrier_perturb_checkpointed_state_for_test`,
`carrier_reject_leading_attempts_for_test`,
`carrier_history_reset_for_test`, `carrier_initial_substep_for_test` and
`chem_cycles_cap_for_test`, and none of them manufactures a composition
that is not a number or a thermal state that is not a gas.

Difference between the state at the entry of the pass and the state after
the pass refused every trial (0 = bit for bit):

| item | difference |
|---|---|
| composition `f_sp` | 0 |
| frozen cell state `bg_cell` | 0 |
| heating, cooling, eta | 0 |
| pressure `p` | 2.4869e-14 |
| temperature `T` | 6.8834e-14 |
| rate state `ieq_rate_cell` | 4.7238e+7 |
| `ieq_ne_cell` | 2.0310e+7 cm^-3 |
| `ieq_TK_cell` | 1.4870e+1 K |
| `ieq_ntot_cell` | 4.7238e+7 cm^-3 |
| `ieq_met_coef` | 5.6610e-7 |
| caloric state, as `p` from the unchanged `rho e` | 1.0711e-3 (p is of order 0.72, so 1.5e-3 relative) |
| base ghost electron term `dp_bc` | 1.0525e-6 |
| base ghost particle count `n_part_cell1` | 1.8047e-4 |
| carrier checkpoint matches | false |

Readings:

- The composition and the frozen cell state ARE restored bit for bit, which
  is what the refusal branch of `relax_photochemical_composition` writes
  back explicitly.
- `p` and `T` differ at 1e-14 because the pass's first act is to recompute
  them from the conserved state through the caloric equation of state,
  while the entry values were seeded from `comp_p_from_T`. This is the
  round trip, not a failure to restore.
- **The rate state the closure evaluator and the certification read, the
  caloric state, and the base ghost quantities are NOT restored.** The
  caloric state is the consequential one: `pressure_from_energy_density`
  maps the unchanged thermal energy to a pressure 1.5e-3 relative away from
  the one the entry state had, and that map is an inverse of the caloric
  equation of state, not an arithmetic identity, so nothing downstream
  cancels it. The module's own header states the rule it breaks: "a path
  that installs a different `f_sp` must refresh these arrays BEFORE it maps
  an energy density to a pressure" (`caloric_eos.f90`, the block above line
  159). The refusal branch restores `f_sp` and does not refresh them.
- The carrier checkpoint does not match. It covers the attempt statistics
  and the history as well as the state, so this is not by itself a
  statement about the physical state; the composition row above is.

---

## 7. Table 5: the conditioning

### 5a. What a row residual of `ieq_res_tol` admits in abundance

The returned composition is a root at the returned temperature (section 4),
so the normalized residual of a row grows from zero as its own unknown is
displaced. `|R_k|` was read at five relative displacements of `x_k` at the
binding cell; `|R_k|/eps` is the sensitivity of the row to its own unknown,
and `ieq_res_tol` divided by it is the RELATIVE abundance error the
residual tolerance admits. The displacement is one-sided, because the code
reports the absolute value of the row residual and a central difference
would be taken through the kink at the root.

| row | `x` | `|R|` at 1e-8 | 1e-7 | 1e-6 | 1e-4 | 1e-2 | `|R|/eps` | relative abundance admitted at `ieq_res_tol` |
|---|---|---|---|---|---|---|---|---|
| H II | 6.953e-6 | 6.953e-14 | 6.953e-13 | 6.953e-12 | 6.953e-10 | 6.953e-8 | 6.953e-6 | 1.4382e-1 |
| He II | 2.555e-7 | 2.555e-15 | 2.555e-14 | 2.555e-13 | 2.555e-11 | 2.555e-9 | 2.555e-7 | 3.9136e+0 |
| He III | 6.937e-11 | 6.937e-19 | 6.937e-18 | 6.937e-17 | 6.937e-15 | 6.937e-13 | 6.937e-11 | 1.4416e+4 |
| H2 | 7.865e-1 | 7.865e-9 | 7.865e-8 | 7.865e-7 | 7.865e-5 | 7.865e-3 | 7.865e-1 | 1.2715e-6 |
| H2+ | 9.827e-9 | 2.221e-17 | 1.362e-16 | 1.276e-15 | 1.266e-13 | 1.266e-11 | 1.266e-9 | 7.8981e+2 |
| H3+ | 2.802e-5 | 1.784e-16 | 3.346e-16 | 1.896e-15 | 1.737e-13 | 1.742e-11 | 1.742e-9 | 5.7418e+2 |
| HeH+ | 1.679e-12 | 6.042e-19 | 6.533e-19 | 1.144e-18 | 5.516e-17 | 5.456e-15 | 5.456e-13 | **1.8328e+6** |

The four identity rows scale exactly as `x_k` by construction, which is the
check that the measurement is reading what it is meant to.

**This is the decisive number of D4a.** `ieq_res_tol = 1e-6` bounds the
relative abundance of the transported H2 at 1.3e-6 and of H II at 0.14, but
it bounds the relative abundance of HeH+ only at 1.8e+6. A relative HeH+
error of 1 per cent produces a row residual of 5.5e-15, nine decades
below the tolerance. The residual tolerance carries no information at all
about the relative abundance of HeH+ at this cell, because the row's
normalized residual is proportional to an abundance thirteen decades below
the gas.

### 5b. A controlled perturbation of the transported parent

`x(H2)` of the binding cell is a quantity the sweep holds fixed there, so a
relative change of it propagates rather than being erased. One sweep after
the change, relative response divided by the perturbation:

| eps | d ln H2+ | d ln H3+ | d ln HeH+ | d ln n_e |
|---|---|---|---|---|
| 1e-8 | 6.7571e-1 | 1.0703e+0 | 8.3615e-1 | 6.1279e-1 |
| 1e-7 | 6.7569e-1 | 1.0703e+0 | 8.3614e-1 | 6.1279e-1 |
| 1e-6 | 6.7568e-1 | 1.0703e+0 | 8.3614e-1 | 6.1279e-1 |

The responses are flat in the perturbation size to four digits, so this is a
derivative and not a numerical artifact, and it is of order one: the chain
does NOT amplify a change of the transported parent. The factors 6.78 and
10.33 of section 5 are therefore not a sensitivity of HeH+ to H2; they are
how the residual of the inner solve, whose stopping tolerance is a relative
STEP criterion at `sqrt(dpmpar(1)) = 1.49e-8` (`ionization_equilibrium.f90`,
the `ieq_inner_tol` block), distributes itself over the three coupled
eliminated ions when they are re-solved together with the electron closure.

---

## 8. Table 6: the stopping rule against the residual at the exit

The closure ends on `dT/T < chem_cycle_tol` alone
(`diffusive_photochemistry.f90`, the exit test inside
`equilibrate_chemistry_at_fixed_conserved_state`). Restarted from the state
it handed back, on an isolated copy:

| quantity | value |
|---|---|
| verdict | `ok = .true.`, reason `converged`, **1 cycle** |
| temperature increment of that cycle | 1.2574e-11, against `chem_cycle_tol` 1e-5 |
| composition after that one cycle, against the returned one | 1.1100e-6, HeH+ at cell 1 |
| HeH+ row residual of the returned composition at the returned temperature | 5.9873e-19, against `ieq_res_tol` 1e-6 |
| largest row residual of any cell at the returned temperature | 2.5929e-14 |

### 6b. The claim in the routine's own header

The header of `equilibrate_chemistry_at_fixed_conserved_state` states: "The
tolerance is the sweep's own reaction-residual tolerance (ieq_res_tol,
1e-6): a temperature change below it cannot move a composition the sweep
has converged to that tolerance." The returned composition was swept again
at displaced temperatures:

| `dT/T` | largest relative movement | species, cell | HeH+ at the binding cell |
|---|---|---|---|
| 0 | 1.1100e-6 | HeH+, 1 | 1.1100e-6 |
| 1e-7 | 7.3240e-7 | HeH+, 12 | 6.5912e-7 |
| 1e-6 | 4.7634e-6 | HeH+, 12 | 3.3990e-6 |
| 1e-5 | 4.5074e-5 | HeH+, 12 | 4.3981e-5 |

**RED for the claim as it is written about abundance.** A temperature change
of exactly `ieq_res_tol` moves the composition by 4.76e-6, about five times
that tolerance; a change at `chem_cycle_tol`, which is the tolerance the
code actually stops on, moves it by 4.51e-5, 45 times `ieq_res_tol`.

**GREEN for the same claim read about the residual.** The largest movement,
4.5e-5 relative in HeH+, corresponds through table 5a to a row residual of
about 2.5e-17, eleven decades below `ieq_res_tol`. The claim is true of the
quantity the tolerance measures and false of the quantity the test measures.

Two further defects of that header, both READ:

1. It names `ieq_res_tol, 1e-6` as the tolerance the cycle uses. The code
   tests `chem_cycle_tol`, whose value is `1.0d-5`, ten times larger. The
   comment states a tolerance the routine does not use.
2. It says "the state handed back is one whose composition, temperature,
   pressure and rate coefficients belong together". The composition handed
   back is the one the last sweep produced at the temperature of the cycle
   BEFORE the exit, and the temperature handed back is the one computed
   after it; the two are `dT` apart, and `dT` is bounded only by
   `chem_cycle_tol = 1e-5`. On this fixture the realized gap is 2.4e-7
   (table 2), well inside, but the contract admits ten times `ieq_res_tol`.

---

## 9. The proof that the diagnostics put the state back

Every diagnostic that advances the state runs on a copy of the composition,
with the module list of section 3 saved and restored. The caloric state is
private to `caloric_eos` and is restored the way the production code
restores it, by refreshing it from the composition being reinstated
(`get_species_densities`, which calls `caloric_state_from_composition`).

The restore is proven by re-evaluating the residual profile of table 2
after the diagnostics and requiring the same bits. Three assertions, all
PASS on this build:

```
PASS the_as_solved_residual_profile_is_reinstated_bit_for_bit          0.0  tol 0
PASS the_as_solved_residual_profile_survives_three_isolated_increments 0.0  tol 0
PASS and_the_returned_composition_is_untouched                         0.0  tol 0
```

The first is taken after the two isolated sweeps of measurement 2, the
second after the three complete closure increments of measurement 3. The
third compares the saved return with the array the production pass wrote.
Because the residual profile is a function of the composition AND of
`ieq_rate_cell`, `ieq_ne_cell`, `ieq_met_coef` and `ieq_mbase_stored`, a
bitwise identical profile is a statement about all of them together.

The two remaining assertions of the driver state the forced refusal
(`a_forced_closure_refusal_ends_the_pass_on_the_chemistry`,
`and_keeps_no_step`) and the two items the production refusal branch does
restore (`a_refused_trial_returns_the_entry_composition`,
`and_the_frozen_cell_state`). The rest of table 4 is printed and NOT
asserted, because it is a measurement of the production code and not a gate
this driver may set.

---

## 10. Which of the two outcomes the measurements support

The measurements support **the second outcome: the test assumes a bound the
residual contract does not give.** The production closure did not exit
before its own promise was met on this fixture.

The numbers, in the order they settle it:

1. The composition production returned is a root of its own network at the
   temperature it was returned with: largest normalized residual over the
   column 2.59e-14, at the binding cell 1.61e-16, and the HeH+ row that
   binds the failing test 5.99e-19, all against `ieq_res_tol = 1e-6`
   (section 4).
2. Restarted from that state the closure reports `converged` in one cycle
   with a temperature increment of 1.26e-11, six decades inside
   `chem_cycle_tol` (section 8). There is no early exit to repair here.
3. `ieq_res_tol` admits a relative HeH+ abundance error of 1.83e+6 at this
   cell (section 7). The assertion
   `one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
   compares a relative abundance movement of 1.11e-6 with a RESIDUAL
   tolerance that bounds that same movement only at 1.8e+6, which is twelve
   decades looser. The comparison is dimensionally wrong, not marginal, and
   no choice of the numerical value of `ieq_res_tol` repairs it.
4. What the 1.11e-6 does measure is the inner solve's own repeatability
   distributed over the three coupled eliminated ions: the inner tolerance
   is a relative step criterion at 1.49e-8, the response of the chain to its
   transported parent is of order one (section 7, table 5b), and the
   movement contracts by about four decades a step (section 5). It is the
   distance of the returned state from the sweep's fixed point at frozen
   temperature, and that distance is not bounded by any residual tolerance.
5. The row's value is not even a function of the state alone: it spreads
   from 1.084e-6 to 1.293e-6 with the pressure the sequence is entered with
   (section 2), while the binding species and cell never move.

Three findings sit beside that verdict and are NOT settled by it. They are
reported because they are real and because D4b's scope touches them:

- **The rate state the certification reads is left by a refused trial.**
  The refusal branch of `relax_photochemical_composition` restores `f_sp`,
  `bg_cell`, `p`, `T`, `heat`, `cool` and `eta`, and does not restore
  `ieq_rate_cell`, `ieq_ne_cell`, `ieq_TK_cell`, `ieq_ntot_cell` or
  `ieq_met_coef`. Anything that asks for the chemical residual of the
  returned state after such a pass is answered against another state's
  rates: 5.33e-4 instead of 1.61e-16 at cell 1, a gap of twelve decades
  (sections 4 and 6).
- **The caloric state is not restored by a refused trial either**, so the
  energy-to-pressure map stands 1.5e-3 relative away from the one the entry
  state had, against the rule the module's own header states. `dp_bc` and
  `n_part_cell1`, the base ghost quantities, are not restored either.
- **The routine's header names `ieq_res_tol` where the code tests
  `chem_cycle_tol`**, a factor ten, and asserts an implication from
  temperature to composition that measurement 6b refutes for abundance
  (section 8).

No tolerance change is recommended here, and none is needed for the verdict:
`ieq_res_tol`, `chem_cycle_tol` and the 1e-20 presence rule were not moved,
and the numbers above are decades away from any of them.
