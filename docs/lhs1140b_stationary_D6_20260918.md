# D6: the helium held outside the atomic stages is reserved from the stage simplex

Item D6 of `docs/PLAN_20260918_rev1.md`, with the refinements of
`docs/PLAN_20260918_review.md` section 8 (the inventory) and
`docs/PLAN_20260918_review1.md` sections 7.5 (the metastable level) and 7.3
(the proton source ledger). Written 2026-09-18 (KST).

## 1. The defects, as they stood in the source

1. **The helium inventory.** `carrier_state` takes the helium nucleus count
   from `element_nucleus_counts`, which counts the nucleus inside HeH+ by
   stoichiometry, and the two ionized stages are carried as fractions of
   that count. `carrier_ionization_stage_projection` bounded their sum by
   one alone, and `limit_to_element_budget` reserved nothing of the helium
   for the molecule. `carrier_source` and `carrier_write_back` then formed
   the neutral stage as `max(N_He - n_HeII - n_HeIII - n_HeH+, 0)`, so a
   partition with a tenth of the helium in HeH+ and the two stages at 0.60
   and 0.35 passed the bound at 0.95 and was written back as a cell holding
   1.05 of the helium it has.

2. **The metastable level.** `carrier_source` reads the He 2^3S population
   from the frozen background and formed the ground singlet as
   `max(n_HeI - n_He(2^3S), 0)`. A trial leaving less neutral helium than
   the frozen level holds has a negative singlet, and the clip returned
   zero and let the trial stand: the rows were evaluated at a partition no
   cell can be in. He 2^3S is a level inside He I, `0 <= n(2^3S) <= n(HeI)`,
   not a fourth stage and not a nucleus of its own.

3. **The proton source ledger.** `mol_heh_rows` returns the proton row
   already split into `p_Hp` and `l_Hp`; `carrier_source` then moves that
   row through `he_h_cx_fvec` (the He <-> H charge-exchange pair) and takes
   the MOVED row as `src(H+)` while the record `sprod(H+)`, `sloss(H+)`
   kept the split from before the move. `sprod - sloss` was therefore not
   the source the residual was assembled from.

## 2. What the admissible set is now

The two ionized stages of a cell may hold the helium nuclei of that cell
less the helium held outside the atomic stages:

    n(He II) + n(He III) <= n_He,avail
    n_He,avail = n_He,nuc - n_He(in molecules) - n_He(2^3S)

`carrier_helium_available_to_stages` is that density, and
`carrier_helium_neutral_of_partition` is the closure that follows from it:
the ground singlet is `n_He,avail - n(He II) - n(He III)` and the neutral
stage is that singlet plus the frozen metastable. Both differences are
nonnegative on the admissible set, so neither is clipped anywhere.

The molecular part is read from the species table (`bsp_nHe` over the
species that are not He I, He II or He III and are not excited levels), so
a further helium-bearing molecule enters it without a coefficient written
at a call site. Today it is HeH+ alone.

The bound is applied in one place and reaches every user of the partition:

| where | how |
|---|---|
| the projection | `stage_simplex_projection` takes the bound cell by cell (`xsum_max`); its diagnostic is now `stage_simplex_sum_over_limit`, the amount by which the sum exceeded the nuclei available to it |
| the element limiter | calls that projection, as before |
| every Newton trial | `solve_carriers` projects the damped trial after the nonnegativity clip and the outer-ghost fill, before its rows are assembled |
| the difference step of the Jacobian | a helium stage at the bound is differenced INWARD, so the derivative is of the chemistry and not of the closure outside the set |
| the source evaluation | `carrier_source` closes the neutral helium through the routine above |
| the write-back | the remainder subtracts the same frozen molecular count; it is nonnegative on the admissible set and is not floored |

The normalization was NOT changed: the stage unknown is still a fraction of
the element's nucleus count, which is the count the stage face flux divides
by (`ionization_stage_transport`), so the flux and source definitions stand
as they are. The reference density of a helium stage
(`carrier_element_reference_density`) also stays the nucleus count: it sets
the floor of the difference step and the absolute floor of the row, a SCALE
rather than a bound, and a scale that can fall to zero where an element is
wholly bound in a molecule would leave the row's diagonal derivative at the
1e-300 guard. The two agree to the molecular and metastable fractions of
the helium, which are trace; the comment at that table now says so instead
of claiming the stages can reach every nucleus.

`carrier_ionization_stage_projection` refuses to project the helium stages
when the frozen helium count of the step has not been built: the bound is a
statement about the composition, and the only bound available without it is
the one this reservation exists to remove.

The proton ledger takes the charge-exchange pair's net on the side of its
own sign, measured as the change the pair makes to the row rather than
rebuilt from the two rate expressions, which live in one place. The helium
entries were already the positive and negative parts of a net source, and
the comments now say so where a reader meets them.

## 3. What was measured

### 3.1 The suites

New suite `src/tests/carrier_helium_inventory/` measures the species the
production routines return. Against a reconstruction of the entry text
(verified: a build of it gives md5 `ebc6449b22c8ba8c5a78fbf2a9946785`,
which is the md5 of the control build made from the untouched tree before
any edit) and against the corrected text:

| assertion | entry text | corrected |
|---|---|---|
| stage sum within the helium available to it | +7.000e-2 | -1.110e-16 |
| ground singlet of the exhausted cell (of its helium nuclei) | -3.500e-2 | +8.327e-17 |
| ground singlet of the counterexample cell | -7.000e-2 | +5.551e-17 |
| helium nuclei of the written state, relative | 5.000e-2 | 1.110e-16 |
| hydrogen nuclei of the written state, relative | 2.220e-16 | 2.220e-16 |
| mass of the written state, relative | 1.618e-2 | 0 |
| source ledger, molecular gas, worst carrier | H+, 6.489e-1 | 0 |
| source ledger, atomic gas, worst carrier | H2, 0 | 0 |
| lower ghosts not written | 0 | 0 |
| electron density of the returned species | 0 | 0 |

The 5.000e-2 is the review's counterexample: the written state held 1.05 of
the helium the cell had. Charge is not conserved by the write-back and is
not asserted to be: the operator transports the ionization partition, which
is what a charge is. What is asserted is that the electron budget the next
sweep reads is the one the written species imply.

`src/tests/ionization_stage_flux/` gained the same statement at the level
of the projection alone (case 4b): a background with a tenth of the helium
in a molecule and a fiftieth in the metastable projects a stage sum of 0.95
onto 0.88 and keeps the ratio of the two stages. 40 of 40 assertions pass.

### 3.2 The suites that link the changed objects

`carrier_boundary_jacobian`, `carrier_constraint_attribution`,
`carrier_outer_boundary`, `carrier_reference_scales`,
`carrier_returned_state_acceptance`, `certification`, `element_operator`,
`species_masses`, `steady_species_rows`, `attempted_step`,
`ionization_stage_flux` and `carrier_helium_inventory` pass.

`carrier_retry` fails on
`one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
at 1.100521341448987e-06 against 1e-06, and `coupled_block_jacobian` fails
on `the_assembled_action_reproduces_the_central_difference_carrier`
(8.9540842e-01) and `the_species_rows_scale_is_the_certifications`
(3 of 13). Every one of those numbers is identical on the entry text, so
they are the standing items D4 and the L22 coupled block, untouched here.

`ionization_imposed_fractions` aborted on both texts: its driver assembles
a metals row without loading the metal charge-exchange rate coefficients of
its cell, and the assembly refuses a set that belongs to another
temperature. A `cx_set_cell(ieq_cell%T_K)` in its setup, where the sweep
does the same, makes it run: 49 assertions, all passing.

### 3.3 Impact

Key-off cases (no `Ionization transport`, so no helium stage is carried):
`wasp_full` (400 steps), `mol_diffusion` (12000 steps) and `mol_carrier`
(12000 steps) give byte-identical data lines in `Hydro_ioniz.txt`,
`Ion_species.txt` and both `*_adv.txt` files between a control build of the
entry text and the corrected one. The provenance comment carries the
wall-clock time of the run and differs; nothing else does.

Key-on: on `hp_front` (100 steps) the largest relative movement of any
column is 1.2e-14 (the He 2^3S column; the hydrodynamic columns move
9.1e-16 or less), which is the round-off of a differently associated
closure. On `carrier_model_a_newton` (`EXHALE_OUTER_PASSES=3`,
`EXHALE_JFNK_MAXIT=5`) the columns move by

| file | column | max relative |
|---|---|---|
| Hydro_ioniz | rho | 1.181e-6 |
| | v | 4.100e-5 |
| | p | 7.390e-8 |
| | T | 1.161e-6 |
| | heat | 1.165e-6 |
| | cool | 8.153e-6 |
| Ion_species | H I | 3.205e-5 |
| | H II | 2.306e-6 |
| | He I | 1.181e-6 |
| | He II | 4.821e-6 |
| | He III | 1.182e-6 |
| | H2 | 1.192e-6 |
| | H2+ | 5.304e-6 |
| | H3+ | 6.165e-6 |
| | HeH+ | 1.014e-5 |

**The constraint did not bind on either case, and the movement is not the
constraint.** On the state `carrier_model_a_newton` returns, the helium
held outside the atomic stages is at most 5.61e-6 of the element and the
stage sum reaches 0.445, so the sum stands 0.55 below its bound; on
`hp_front` the reservation reaches 1.58e-3 and the sum 0.675, 0.32 below
the bound. Three further measurements pin the mechanism:

- A build of the corrected text with the trial projection taken back out
  gives data lines byte-identical to the corrected build on
  `carrier_model_a_newton` (3 passes, maxit 5), and so does a build with
  the inward difference step also taken out. Neither fires on this case.
- The same solve at `EXHALE_OUTER_PASSES=1` is byte-identical between the
  two builds at `EXHALE_JFNK_MAXIT=1` and at 5; the difference enters at
  the second outer pass.
- The four gated carrier rows of the three-pass solve end at the same
  measures and the same cells in both builds: H2 1.208e-3 at 248, H+
  2.422e-2 at 267, He+ 9.908e-1 at 248, He++ 8.130e-2 at 248.

What remains is the different floating-point association of the neutral
helium closure -- `(N_avail) - n_HeII - n_HeIII` in place of
`max(N_nuc - n_HeII - n_HeIII - n_HeH+, 0)` followed by
`max(n_HeI - n_TR, 0)` -- a round-off-level change of the rows, which the
100-step marching case shows at 1.2e-14 and which the stationary iteration
amplifies to 1e-5 over three passes.

That the reservation is inactive on these two states is a measurement of
the states, not a reason to leave the projection as it was: conservation of
helium nuclei is an invariant the operator claims, and a projection that
admits a state violating it is wrong whether or not a run has reached one.

## 4. Noticed outside the item, not fixed

- `steady_newton.f90` line 14313 builds the carrier densities of a
  diagnostic row as `nc = fc(jcell,:)*nrho_c(jcell)`, the density of the
  mass row, for EVERY carrier. An ionization stage's unknown is a fraction
  per element NUCLEUS (`carrier_nucleus_reference`), so the H+, He+ and
  He++ entries of that diagnostic are wrong by the ratio of the element
  nucleus density to `rho n0` whenever a stage row is reported. The file
  belongs to another item.
- `element_inventory_report` skips helium in its element loop
  (`if (ie .eq. ien_He) cycle`), so the one report that could have shown a
  helium breach directly shows only the He/H ratio departure.
- The hydrogen and oxygen closures of `carrier_source` still clip a
  negative remainder (`n_hi`, `n_o0`), and the difference step of a
  molecular carrier is not turned inward at its budget. The same treatment
  as the helium closure would need the hydrogen and oxygen budgets in the
  trial projection, which is a larger change than this item.
