# Element inventories by calculation context (N4a deliverable)

The three-context assignment of every species that carries a hydrogen,
helium, oxygen or carbon nucleus, and the constraints each context states.
Advisor draft of 2026-09-09 for PLAN_20260909_rev1 item N4a, verified row by
row against the source on 2026-09-09 by the N4a worker; every correction the
verification made is marked **corrected** with what the source says. The map
this document specifies is `src/modules/lower_atmosphere/element_inventory.f90`
and its rows are `src/tests/steady_species_rows/`.

Sources READ (2026-09-09 working tree): `diffusive_photochemistry.f90`
(`carrier_state`, `hydrogen_available_to_carriers`,
`carrier_element_headroom`, `carrier_source`, `limit_to_element_budget`,
`carrier_write_back`, `carrier_steady_residual`); `steady_newton.f90`
(`set_transported_species_rows`, `freeze_species_unknown_box`,
`eval_residual`, `write_species_rows_into_composition`,
`arm_carrier_log_unknown`); `System_HeH_mol.f90` and
`System_HeH_mol_metals.f90` (the unknown list of the local solve);
`ionization_equilibrium.f90` (the element totals and the CO reservoir);
`species_table.f90` (`bsp_nH`, `bsp_nHe`, `bsp_nO`, `bsp_nC`, `bsp_charge`,
the metal stage tables); `element_census.f90` (the full-state inventory and
its gate); `certification.f90` (`cert_tol_carrier`, `cert_tol_element`).

## 1. What the code calls "free" and "available"

| Symbol | Definition in the source | What it is |
|---|---|---|
| `nH_free(j)` | sum over base species with `bsp_nH > 0` (excited levels skipped) of `bsp_nH * n`, minus `n(HeH+)`, floored at 0 (`element_inventory.f90: carrier_element_totals`, called by `carrier_state`) | the TOTAL hydrogen nuclei of the cell outside HeH+ (H I, H II, H2, H2+, H3+, OH, H2O all included). Not "free" in the sense of unbound. |
| `nO_free(j)` | sum over base species with `bsp_nO > 0` of `bsp_nO * n`, plus **all three oxygen stages including the neutral one** (`do k = 0, melem_top(iel_O)`), floored at 0 | the TOTAL oxygen of the cell (molecular and atomic stages). **corrected**: the earlier wording said "every oxygen ion stage", which reads as excluding O I; the loop starts at `k = 0`. The distinction matters because `cbg_nOion` uses `do k = 1, melem_top` and therefore excludes O I. |
| `nC_free(j)` | same for carbon, all stages from the neutral up | the TOTAL carbon. |
| `nH_avail(j)` | `max(nH_free - n(H II) - 2 n(H2+) - 3 n(H3+), 0)` without proton transport; `max(nH_free - 2 n(H2+) - 3 n(H3+), 0)` with it (`element_inventory.f90: carrier_hydrogen_budget`, called by `hydrogen_available_to_carriers` with the step's frozen background) | the hydrogen the carriers may hold: the total less what the step FREEZES and rescales (H II when it is not a carrier; H2+, H3+ always). H I is not subtracted: the write-back gives H I the remainder. **corrected**: the floor at zero was missing, and the `cbg_*` densities are the frozen background `carrier_state` wrote from the `f_sp` it was called with, which is the last sweep's state on the marching path and the outer iterate's in the stationary solve. |
| `carrier_element_headroom(ic)` | H2: `nH_avail/2`; H+: `nH_avail`; OH, H2O: `nO_free` each; CO: `min(nO_free, nC_free)` | one box side per carrier from ONE element's budget, divided by the nuclei of that element the carrier holds; no joint bound. |

## 2. Species by context

Legend: T = transported unknown (a Newton unknown of the stationary system,
or a carrier of the transport operator); E = eliminated at the current state
(solved by the local equilibrium sweep, re-solved at every evaluation);
C = closure of an element budget (its density is the remainder, not a row);
F = fixed data for the step (frozen at the last sweep and rescaled by the
write-back, or computed before the sweep); R = a reservoir condition at the
base cell.

| Species | local source solve (`System_HeH_mol*`, one cell) | carrier operator (`carrier_steady_residual`, `carrier_source`, write-back) | global stationary solve (`eval_residual`, species rows) | certification |
|---|---|---|---|---|
| H I | **C**, not a row: the unknowns are `x(1) = n(H II)/n_H`, `x(4..6) = 2 n(H2)/n_H, 2 n(H2+)/n_H, 3 n(H3+)/n_H`, `x(7) = n(HeH+)/n_H`, and "Neutral atomic H and neutral He close the element budgets" (`System_HeH_mol.f90` header) | F, rescaled to the hydrogen remainder | C through the sweep; its density is the write-back remainder | closure check of H nuclei |
| H II | E (row 1, the complete proton balance) | F, rescaled (no proton transport) / T (`ic_Hp`, key `Ionization transport: True`, which `input_read` refuses without carrier transport) | E / T likewise | same |
| H2 | E (row 4) unless `ieq_cell%x_h2_fixed`, which replaces row 4 by `x(4) - x_h2_fix` | T (`ic_H2`, always in the carrier set) | T when `carrier_transport` and `thereis_mol`; else E (refused with a coupled solve since B5i) | carrier row, `cert_tol_carrier = 1e-8` |
| H2+, H3+ | E (rows 5, 6) | F, rescaled with H I | E | none of their own |
| HeH+ | E (row 7) | F, its H nucleus removed from `nH_free` (rescaling it would move helium) | E | none of its own |
| He I | **C**, not a row: neutral He closes the helium budget, which includes the He nucleus of HeH+ | F (helium not a carrier) | C through the sweep | element row |
| He II, He III | E (rows 2, 3) | F | E; helium ELEMENT total T when `He_diffusion` (`srow_element_he`, a mass fraction, base R) | element row, `cert_tol_element = 1e-8` |
| He 2^3S | E (row 8 when `thereis_HeITR`); an excited LEVEL of He I, so every budget sum skips it (`bsp_is_excited_level`) | F | E | none of its own |
| OH, H2O | E (rows `iox`, `iox+1`, `iox = 8` or 9 with the triplet) | T (`ic_OH`, `ic_H2O`) under the oxygen chemistry | T when carried; else E | carrier rows |
| CO | **F**, not a row and not eliminated by the network: `ioniz_eq` takes CO from the CO to C + O chemical equilibrium of the local (n, T) (`co_equilibrium_density`), or, with the carriers transported, from the transported reservoir clamped to the C and O element totals, and then removes it from both element totals before the cell sweep (`ionization_equilibrium.f90` 1039-1062). **corrected**: the destruction expression `-(k_D1 n(He+) + k_ph) n(CO)` is the CARRIER operator's CO row (`carrier_source`), not the local solve's. | T (`ic_CO`), destruction only: `He+ + CO` (UMIST RATE22 4068) and the shielded 912-1201 A photodissociation | T when carried | carrier row |
| O I .. O III, C I .. C III | E (the metal block), solved against the FREE family `nm_tot = nm_el - n(CO)`, with `nm0(iel_O)` further reduced by `n(OH) + n(H2O)`, so O I means FREE ATOMIC oxygen | F, rescaled to the oxygen / carbon remainder after the carriers | E; the ELEMENT totals O and C are T when `He_metal_diffusion` (`srow_element_trace`, base R) | element rows |
| N, Mg, Ca, Na, Fe stages | E | not touched by the carriers | E; element totals T under `He_metal_diffusion` | element rows |
| Si, K, S (not in `metals.inp`) | absent | absent | **NOT registered**: `set_transported_species_rows` skips an element whose reservoir abundance is zero (`if (melem_ab(im) .le. 0.0d0) cycle`, `steady_newton.f90`). **corrected**: the earlier row recorded the state before that skip, when an absent element was registered with an identity row on its own lower bound (1500 unknowns of the atomic fixture). The skip is in the 2026-09-09 working tree and is the N1 owner's change; the map treats an element with zero reservoir abundance as absent either way. | none |
| electrons | E (charge closure `calc_ne`) | recomputed from the written species | E | charge closure |

## 3. Constraints the code states, and the ones it does not

**Stated as box sides** (`carrier_element_headroom`, made a face of the
coupled solve's unknown box by `freeze_species_unknown_box`):
`n(H2) <= nH_avail/2`, `n(H+) <= nH_avail`, `n(OH) <= nO_free`,
`n(H2O) <= nO_free`, `n(CO) <= min(nO_free, nC_free)`; and `carrier_log_floor`
below in `ln n`.

**Stated as shared sums, on the marching path only**
(`limit_to_element_budget`, run inside `solve_carriers` before the
write-back). **corrected**: the earlier draft said the shared constraints
were "enforced only after the fact, by the write-back". They are enforced
before it, by that limiter, which
- clamps `n(CO)` to `nC_free`;
- scales OH and H2O together into `max(nO_free - n(CO), 0)`;
- scales the whole H-bearing carrier set (H2, OH, H2O and, where it is
  carried, H+) by one factor when
  `2 n(H2) + n(OH) + 2 n(H2O) [+ n(H+)] > nH_avail`.

What has NO shared statement is the coupled stationary solve: its unknown box
carries the five separate sides above, and the screen of `eval_residual`
compares the NUMBER of cells outside a side with the iterate's own number.
Neither bounds a shared sum, and neither bounds a magnitude.

**Redistribution, and where the element total moves**
(`carrier_write_back`): the hydrogen remainder
`rest = nH_free - 2 n(H2) - n(OH) - 2 n(H2O) [- n(H+)]` is shared over H I
(and H II, H2+, H3+ by proportional rescaling); with `rest <= 0` and
`held > 0` those stages are ZEROED; the oxygen remainder
`nO_free - n(OH) - n(H2O) - n(CO)` and the carbon remainder
`nC_free - n(CO)` are shared over their stages through `max(rest, 0)`.

**corrected**, and the sign matters: a negative remainder does not destroy
nuclei, it CREATES them. The carriers keep what they hold and the closure
species cannot go negative to absorb the excess, so the element total of the
written state exceeds the total it was formed from by `|rest|`. The
write-back's invariant is that the element totals it is handed come back
unchanged, and `max(rest, 0)` is the one place that invariant breaks.
`element_census_verify` around the call is what sees it, at its gate of
`1e-9` relative, and only when `EXHALE_ELEMENT_ASSERT` is 1 or 2: with the
census off, which is the default, nothing reports it at all.

**Not stated anywhere** (the physical feasible set of the stationary solve):
- hydrogen: `2 n(H2) + n(OH) + 2 n(H2O) + n(H+) + n(H I) + 2 n(H2+) + 3 n(H3+) = nH_free`,
  so the carriers alone obey `2 n(H2) + n(OH) + 2 n(H2O) [+ n(H+)] <= nH_avail`;
- oxygen: `n(OH) + n(H2O) + n(CO) <= nO_free`;
- carbon: `n(CO) <= nC_free` (stated as a side, but `min(nO_free, nC_free)`
  for CO ignores that OH and H2O also draw on `nO_free`);
- OH and H2O carry hydrogen that their own side (`nO_free`) does not mention;
- the oxygen budget of the SIDES is the element total, while the chemistry
  rows close free atomic oxygen against `nO_free - n(O ions above neutral)`
  (`carrier_source`: `n_o0 = max(nO_free - cbg_nOion - n(OH) - n(H2O) - n(CO), 0)`).
  The box therefore admits states whose own source rows have to clamp free
  atomic oxygen at zero. Hydrogen has no such gap: the side and the closure
  both go through `hydrogen_available_to_carriers`.

**A second constraint the sides do not touch: the ratio of two elements.**
The box bounds the partition of one element among its species. It says
nothing about the amounts of two elements standing in the ratio the run was
given, and nothing restores that ratio when a solve writes a species density
directly: `write_species_rows_into_composition` writes the carrier column
AFTER the element projection, so the Newton's `n(H2)` is not rescaled onto an
element total, and with helium diffusion off no operator re-imposes one.
MEASURED under the map (2026-09-09, rows
`inventory_reference_*_keeps_its_He_over_H` and
`inventory_coupled_state_He_over_H_is_seen`): the marching goldens hold the
He/H nucleus ratio to their base cell's within 3.4e-13 to 3.8e-11, and two
saved coupled-carrier JFNK states depart from it by 4.0e-2 and 6.7e-1 while
every carrier budget of the same states holds exactly.  Seen live inside one
such solve under `EXHALE_INVENTORY_REPORT=1`, the departure grows from
1.4e-11 at the first iterate through 4.9e-5 at the third, the first whose
`n(H2)` came from a Newton step, to 1.1e-1 at the thirteenth, while the
hydrogen occupancy stays at 0.985 and the worst carrier breach stays at
zero.

## 4. The map (N4a), and what it does not do

1. One stoichiometric source: `bsp_nH`, `bsp_nHe`, `bsp_nO`, `bsp_nC` and
   the metal stage tables carry the nuclei counts, and
   `element_inventory.f90` reads them. `carrier_element_totals`,
   `carrier_hydrogen_budget` and `element_box_side` are that module's, and
   `carrier_state`, `hydrogen_available_to_carriers` and
   `carrier_element_headroom` call them, so no call site of the operator
   carries a stoichiometric coefficient of its own.
2. Three contexts of one constraint: (a) the full inventory of a state,
   grouped as carried / locked / reserved / closure, with
   `n_total = n_carried + n_locked + n_reserved + n_closure` for each
   element; (b) the conditional budget of a calculation that holds some
   densities fixed, evaluated against a budget FROZEN at another state
   (`element_inventory_of_candidate`), which is what makes the constraint
   more than an identity; (c) the remainder a consistent redistribution
   leaves, `element_free - n_carried`.
3. Violation magnitude, not count:
   `max(0, 2 n(H2) + n(OH) + 2 n(H2O) [+ n(H+)] - nH_avail)/nH_avail` and the
   oxygen and carbon analogues, cell by cell, with their maxima and the
   element that carries the worst. The occupancy `n_carried/carrier_budget`
   is reported beside them, so a state's distance from the constraint is a
   number and not a verdict. `eval_residual`'s count comparison stays where
   it is and is a diagnostic beside them.
4. Absent elements are not part of the feasible set and are not unknowns;
   the map treats an element with zero reservoir abundance as absent.
5. The write-back's `max(rest, 0)` clamps are where the element total moves.
   Under the map a negative remainder is an infeasible state reported by
   magnitude; the clamp itself is unchanged by this item.
6. What the map does NOT do in this item: it reports and does not constrain.
   Trial construction, the chemistry rows and certification still enforce
   the five separate sides, so nothing here establishes that the three
   enforce one constraint. That integration is N4b's, and the Stage B gate
   may not claim it from this item.

> **Note 2026-09-10 (decision 22 a, item N30):** the element and carrier row
> tolerances named here as 1e-8 are 1e-5 in the wind (r >= 1.20) and reported
> below it since N30.
