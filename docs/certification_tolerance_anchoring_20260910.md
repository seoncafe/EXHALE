# Anchoring the species-row certification tolerances

PLAN_20260909_rev1 items N8b (the anchoring) and N8c (the same five anchors
re-measured after N29); `docs/ISSUES_20260909_review.md` section 5.4 and
`docs/ISSUES_20260909.md` 3.1 (d) and (e). 2026-09-10.

`cert_tol_element` and `cert_tol_carrier` both stand at `1e-8` and both are
UNCONFIRMED: no configuration with a species row has converged under them, so
neither value has ever been measured against a state. This memo produces the
measurements the review asks for and a PROPOSAL. **No tolerance was changed by
this item.** The proposal is section 7 and is the user's decision.

Every table carries a second reading, "after N29". N8b measured the anchors
on candidate states whose species did not carry their own density: N29 found
the element projection creating mass (6.4e-3 on the atomic candidate) and
fixed it. Item N8c produced the two candidates again with the corrected
operator and re-measured all five anchors on them. Section 9 gathers what
moved and what did not; the reading of section 7 is unchanged by it.

Every number below is labelled MEASURED (N8b or N8c ran it) or READ (from
source or from another item's report). The two candidate states are the ones
the reloads of N26 handed back:

- **atomic element candidate**: HD 209458 b, `He_diffusion`, `He_metal_diffusion`,
  eight element rows. Before N29 it was handed back at `||R||` 3.719e-4 with
  the helium row of cell 246 at 6.2e-4 of its scale (READ, N26/N27); after
  N29 the same reload is handed back at `||R||` **2.271e-4** in 275
  iterations with the helium row of cell 245 at **5.932e-4** (MEASURED, N8c).
  Re-entered here through `Restart intent: stationary evaluate`, which
  measures the state as loaded and takes no step.
- **carrier candidate**: the hot-Uranus molecular column with the H2 carrier
  transported, handed back at `||R||` 0.25 with the H2 row of cell 205 at
  0.48 (READ, N26). Re-entered the same way; the H2 row reproduces at
  **4.7996747301264420e-1** at cell 205 (MEASURED, and the same seventeen
  digits after N29: the carrier reload carries no metals and no element row, so
  the projection N29 corrected is never reached and the handed-back state is
  byte-identical to N26's).

Before N29 both reloads carried a caveat, measured in section 6: the writer
and the loader did not return the same state. N29 closed it.

## 1. The regimes, and why one number cannot serve both

An element or carrier balance is a cancellation between transport terms whose
ratio is set by the flow. Above the sonic point the advective divergence
carries the row; below the homopause the eddy and settling terms carry it and
each is orders of magnitude larger than their difference, so the same
equation is resolvable to very different precision in the two places. The
operator says so itself (READ, `src/modules/functions/binary_element_diffusion.f90`
line 1869-1873):

> the reachable floor is set by the conditioning of the tridiagonal solve and
> varies by orders of magnitude across the cases here -- 1e-12 in a wind,
> ~1e-5 in the K_zz = 2e12 homopause column of test T7, where the time term
> is negligible against the eddy term.

The certification now reports the two separately. `cert_regime_layer_r = 1.10`
and `cert_regime_wind_r = 1.20` (`certification.f90`) bound two REPORTING
windows; no row is accepted or refused by them, and the whole-column maximum
remains the only number the verdict reads.

MEASURED on the two candidates, `EXHALE_CERT_ANCHOR=1`:

| row | state | whole column | layer, r < 1.10 | wind, r >= 1.20 |
|---|---|---|---|---|
| element He/H (atomic) | N8b | 5.974e-4, cell 245 | 5.602e-4, cell 197 | 2.983e-4, cell 292 |
| element He/H (atomic) | after N29 | 5.932e-4, cell 245 | 5.793e-4, cell 196 | 3.048e-4, cell 291 |
| element Fe (atomic) | N8b | 5.837e-4, cell 106 | 5.837e-4, cell 106 | 3.086e-6, cell 292 |
| element Fe (atomic) | after N29 | 5.922e-4, cell 114 | 5.922e-4, cell 114 | 5.419e-5, cell 479 |
| carrier H2 | N8b | 4.800e-1, cell 205 | 4.800e-1, cell 205 | 7.186e-2, cell 291 |
| carrier H2 | after N29 | 4.800e-1, cell 205 | 4.800e-1, cell 205 | 7.186e-2, cell 291 |

The wind reading is 2 to 200 times smaller than the layer reading on the same
equation and the same state (2 to 11 times after N29, whose better resolved
iron row raises the iron wind reading by a factor 18). Averaging them into one threshold would either
refuse every wind state a code of this discretization can produce or admit a
layer that is not resolved at all.

## 2. The five anchors, MEASURED

### (1) Reproducibility

| what | element rows (atomic candidate) | after N29 | carrier row | after N29 |
|---|---|---|---|---|
| two evaluations of one state | **0 exactly**, all 17 digits of every row | **0 exactly** | **0 exactly** | **0 exactly** |
| the closure's seed family (READ, N5 contract 2) | 3.07e-3 `cert_tol_element` at a relative seed 1e-6, linear, i.e. **3.07e-5** in the row measure per unit relative seed | unchanged (READ) | 4.27e-4 `cert_tol_carrier` at 1e-6, i.e. **4.27e-6** per unit relative seed | unchanged (READ) |
| the seed the state's own sweep carries | 1.736e-10 | **1.645e-14** | 2.617e-14 | 2.617e-14 |
| so the closure floor of this state | **5.3e-15** | **5.0e-19** | **1.1e-19** | **1.1e-19** |
| a one-ulp change of the stored state | see below | see below | see below | see below |

The atomic state's own sweep seed falls by four decades because the loader no
longer has to repair a composition whose mass sum missed the file's density:
the reload equilibrates in one pass at 1.6e-14 instead of three passes at
1.7e-10 (MEASURED).

The state-perturbation ladder (relative perturbations of the transported
composition, drawn independently per cell; three sizes, linear in all three):

| where | amplification, d(row measure)/d(relative composition) | after N29 | floor at eps = 2.2e-16 | after N29 |
|---|---|---|---|---|
| element He/H, base region r ~ 1.001 | **1.1e3** | **1.15e3** | 2.4e-13 | 2.5e-13 |
| element He/H, r ~ 1.010 | **2.0e3** | **2.03e3** | 4.4e-13 | 4.5e-13 |
| element He/H, binding cell r = 1.152 | 0.30 | 0.50 | 6.6e-17 | 1.1e-16 |
| element He/H, layer maximum r = 1.084 | not read separately | 0.70 | | 1.6e-16 |
| element He/H, wind r = 1.269 (1.266 after N29) | 0.23 | 0.39 | 5.1e-17 | 8.7e-17 |
| carrier H2, layer r = 1.098 | 0.26 | 0.26 | 5.7e-17 | 5.7e-17 |
| carrier H2, wind r = 1.287 | 0.53 | 0.53 | 1.2e-16 | 1.2e-16 |

**Neither reproducibility nor double-precision representation forbids 1e-8
anywhere except in the base region**, where the row amplifies a composition
change by three decades and the floor is 4e-13. A run that only knows its
composition to 1e-12 there, which is what the closure delivers, cannot read
its base element row below 1e-9.

### (2) Derivative convergence

The banded Jacobian action against the directional difference of the same
map, on a direction supported on one cell, seven probe lengths
(`EXHALE_SPECIES_JAC_TEST=1`, `EXHALE_JAC_TEST_CELL=j`; the cell selector is
new in this item). MEASURED on the atomic candidate at its first iterate:

| cell | radius | least relative error | after N29 | at probe | shape |
|---|---|---|---|---|---|
| 1 | 1.0002 | **4.52e-7** | **7.63e-8** | 1.49e-8 | a clean V |
| 245 | 1.1533 (the binding cell) | **6.41e-5** | **6.62e-5** | 1.49e-10 | a PLATEAU at 6.6e-5 to 6.9e-5 over four decades of probe length |
| 292 | 1.2692 (wind) | **2.10e-5** | **2.12e-5** | 1.49e-6 | a shallow V |
| whole system | | 2.26e-7 all rows, 1.03e-7 species rows | 5.37e-8 all rows, 5.57e-8 species rows | 1.49e-8 | |

Both readings are taken at the first JFNK iterate of a run entered with
`Restart intent: relaxation`, which is two CFL steps past the loaded state
and not the candidate itself; the two are compared to each other and not to
the candidate's own rows.

The plateau matters. At the base and over the whole system the error falls
like the probe and then rises like its reciprocal, which is the signature of
a correct pair, and the least error is the arithmetic optimum. At the binding
cell it does neither: it sits at 6.7e-5 whatever the probe is, so what is left
is not truncation and not cancellation but the part of the true derivative
that the band does not hold. **A Newton step there is built on a model that is
6.4e-5 wrong**, so it can remove at most 1 - 6.4e-5 of the row it is aimed at
and the smallest row it can resolve in one step from this state is
6.4e-5 x 6.0e-4 = **3.8e-8**; in the wind, 2.1e-5 x 3.0e-4 = **6.3e-9**. Those
two numbers sit within a factor 4 of `cert_tol_element` itself. This is the
first of the five anchors that is anywhere near 1e-8. After N29 the same
product reads 6.62e-5 x 5.93e-4 = **3.9e-8** at the binding cell and
2.12e-5 x 3.05e-4 = **6.5e-9** in the wind.

### (3) A manufactured steady column, three nested grids

`src/tests/certification` now builds the measurement so that it can be
repeated without a planet run. The column: helium and ten trace elements
transported, WENO3, no gravity, no ambipolar field, thermal diffusion
`he_alphaT = 0.15` against a temperature gradient, `r^2 F_rho` constant (a
steady mass row), uniform grid over [1, 2] with 801, 401 and 201 cells whose
cell centers are nested, so a fine solution restricts onto a coarse grid with
NO interpolation. The composition is written as mass fractions that sum to
one, so the column carries its own density exactly.

Two states are read:

- **the flat column with the drift off**, which is analytically a steady
  state of the operator (uniform composition, constant `r^2 F_rho`): the row
  reads **2.72e-14**, and that is the assembly's own round-off floor. After
  N29: 2.72e-14, the same number.
- **the relaxed column with the drift on**, whose steady state has a real
  gradient (helium contrast 3.73e-3 over the column, 3.02e-3 after N29,
  because the projection no longer adds mass to the helium it moves):

| cells | dr | row measure | after N29 |
|---|---|---|---|
| 801 | 1.25e-3 | 5.41e-13 (converged; the relaxation's own floor) | 7.90e-13 |
| 401 | 2.50e-3 | **2.466e-6** | **2.469e-6** |
| 201 | 5.00e-3 | **7.419e-6** | **7.429e-6** |

observed order **1.589** before and after. Extrapolated to the production
spacing at the binding radius (dr = 1.8447e-3 R_p, MEASURED on the atomic
fixture's grid, which N29 does not move) this column would read 1.5e-6
(1.52e-6 after N29). The column's own mass closure after the relaxation
falls from 3.37e-3 to 2.27e-15, which is section 3.

That is a LOWER bound on the production discretization error, because this
column's gradient is far gentler than a real ionization front's. The
production number is section (4).

### (4) Grid refinement of the candidate itself

**[corrected 2026-09-11]** This anchor and the remapping control below are
CONSISTENCY DIAGNOSTICS of one candidate on two grids, not a measurement of
a floor. Two reasons, both from review section 3.2 of
`docs/solver_approach_analysis_20260910_review.md`. First, the state
evaluated on the N = 1000 grid is a remap of an unconverged N = 500
candidate, not a solution of the N = 1000 discrete equations, so the
Richardson-style split below ("about 78 percent is discretization") is an
estimate of how the evaluated row responds to refinement and not a division
of the residual into an unavoidable and an avoidable part. A nearly exact
round trip through the interpolant bounds the round trip, not the derivative
error on the intermediate grid. Second, and whatever the split: a truncation
error of this size does not bound the algebraic residual of the discrete
equations from below (review section 3.1), and MEASURED on this same
500-cell grid the existing H2 transport routine took the wind H2 row from
7.380744e-2 to **2.449021e-13** on a held background
(`docs/solver_partition_experiment_20260911.md` section 4). The numbers
below stand as measured; what is withdrawn is reading them as a floor.

The atomic candidate mapped onto the N = 1000 grid of the same configuration
and evaluated with no step (`Restart intent: stationary evaluate`), at one
fixed physical radius r = 1.1515:

| grid | dr at the probe | He/H row at the probe | after N29 |
|---|---|---|---|
| N = 500 (the state's own) | 1.84e-3 | **5.9678e-4** | **5.9273e-4** |
| N = 1000 | 9.2e-4 | **2.4801e-4** | **2.4677e-4** |

ratio 2.41 (2.402 after N29), i.e. the row of THIS state reads about 2.4
times smaller on the finer grid. The remap control (N = 500 -> 1000 -> 500 with
the same interpolant) returns 5.96770e-4 against the original 5.96779e-4, a
relative difference of 1.6e-5 (5.92725e-4 against 5.92734e-4, 1.5e-5, after
N29), so **at that radius the remap contributes at most 0.02 percent and the
factor 2.4 is a grid effect**. Richardson at order 2 puts the h -> 0 limit at
1.32e-4 (1.31e-4 after N29): the estimate puts **about 78 percent of the
5.9e-4 the state is refused for in the direction of the discretization on
its own grid**, and about 1.3e-4 in the h -> 0 limit of the same estimate.
Both readings agree to two figures. [corrected 2026-09-11: read as a
diagnostic, per the note opening this anchor. The wording was "is
discretization ... and about 1.3e-4 is a genuine imbalance", which claims a
split this experiment does not establish: the fine-grid evaluation is a
remap of an unconverged candidate and the assumed order is imported from a
smooth manufactured column, review section 3.2.]

Two limits of this test, stated rather than smoothed over:

- the interpolant has to be C2. A linear remap of log(composition) injects a
  residual of order 1/h at the coarse nodes and gives 1.158e-3 on the fine
  grid, i.e. it reads LARGER instead of smaller (MEASURED, first attempt).
- below r ~ 1.03 the remap itself dominates whatever the interpolant: the
  base rows of the remapped state read 0.99 against 5.6e-4 for the state on
  its own grid (0.99 on the N = 1000 remap and 0.20 on the round trip after
  N29). **The refinement statement holds for r >~ 1.03 only.** The
  reason is section 2 (1): the base region amplifies a composition change by
  1e3, so a remap accurate to 1e-6 leaves a base row of 1e-3.

### (5) Conserved fluxes

The face fluxes the element operator itself carries
(`output/element_flux_profile.txt`, `EXHALE_DIFFUSION_CHECK=1`), reduced to a
median and a relative radial spread the way `reduce_element_flux_windows`
does (READ). MEASURED on the atomic candidate:

| window | F_He spread | after N29 | F_H spread | after N29 | face mass flux spread | after N29 |
|---|---|---|---|---|---|---|
| layer, r < 1.10 | **39.3** | **9.79** | 37.3 | 9.09 | 37.8 | 9.25 |
| wind, r >= 1.20 | **2.64e-2** | **2.71e-2** | 1.55e-2 | 1.36e-2 | 5.69e-3 | 3.99e-3 |
| escape window, r >= 2.0 | **2.17e-3** | **2.01e-3** | 1.67e-3 | 2.14e-3 | 1.31e-3 | 1.55e-3 |

beside the certification's own cell-centered mass-flux spread of the state,
**5.497e-7** (READ, the N26 run; **9.849e-7** after N29). The element flux is
conserved to 2e-3 where the same faces conserve mass to 1.5e-3, and in the
diffusion-dominated layer neither is conserved at all: the standing base wave
puts the spread at 40 times the median before N29 and still at 9 times after
it. **The 1e-6 the candidate is quoted with is a cell-centered mass flux and says
nothing about the element fluxes**, which are three to four decades worse in
the same state. N29 buys a factor 4 in the layer and nothing in the wind,
which is what a mass-closure fix should do: the leak was largest where the
composition changed most.

## 3. A conservation error the inventory has no row for (CLOSED by N29)

Found while measuring (5) and reported here because it bounds everything
above. **N29 fixed it**; the paragraphs below are the measurement as N8b
made it, with the post-N29 reading beside each number. The candidate state's
species did not carry its own density:

    max_j | sum_i m_i n_i(j) - rho(j) | / rho(j)  =  6.68e-3   (MEASURED)

with the code's own mass policy (`calc_rho`), 1.9e-9 at the base and growing
monotonically outward. The three-unknown control with no element transport
(`backup/regression/wasp_full_newton`) reads **1.93e-9** on the same measure,
so this was specific to the element-transport path.

It is reproduced in the suite on a synthetic column that starts exactly
closed: **6.2e-16 as built, 3.37e-3 after the element transport relaxation**,
with a helium contrast of 3.7e-3, i.e. the closure error was the same size as
the composition change the operator made. `project_elements` did not give
back to hydrogen the mass it takes from helium.

**After N29** (MEASURED, N8c, the same three readers):

| reader | N8b | after N29 |
|---|---|---|
| the atomic candidate's own state | 6.68e-3 (6.4e-3 on N26's file) | **1.95e-9** |
| the carrier candidate's own state | not read | **1.93e-9** |
| the suite's synthetic column after the relaxation | 3.37e-3 | **2.27e-15** |
| the three-unknown control, `golden/wasp_full_newton` | 1.93e-9 | 1.93e-9 (unchanged) |

1.9e-9 is the precision of the written file itself, which the three-unknown
control has always read, so the element-transport path is now closed as well
as any state the code writes.

The two consequences the closure had are gone with it: (a) a reload
reconstructs rho from the species (`load_IC`, `calc_rho`) and now loads the
density the file's own `rho` column carries, which is section 6; (b) no row
of the stationary inventory measures this quantity, so a state could be
certified while its composition and its density described two different
atmospheres. (b) is still true and is still not this memo's to fix: what
measures it now is the element operator's own test suite, not the run-level
census.

## 4. What binds, per row class and regime

Each cell reads "N8b -> after N29"; one number means the two agree.

| anchor | element, wind | element, layer | element, base r < 1.03 | carrier, wind | carrier, layer |
|---|---|---|---|---|---|
| repeated evaluation | 0 | 0 | 0 | 0 | 0 |
| representation (state x eps) | 5.1e-17 -> 8.7e-17 | 6.6e-17 -> 1.6e-16 | 4.4e-13 -> 4.5e-13 | 1.2e-16 | 5.7e-17 |
| closure seed | 5.3e-15 -> 5.0e-19 | 5.3e-15 -> 5.0e-19 | 5.3e-15 -> 5.0e-19 | 1.1e-19 | 1.1e-19 |
| derivative | **6.3e-9 -> 6.5e-9** | **3.8e-8 -> 3.9e-8** | 4.5e-7 -> 7.6e-8 x the row | not measured | not measured |
| discretization | **~4.7e-4 -> ~4.6e-4** (from (4)) | ~4.7e-4 -> ~4.6e-4 | not separable from the remap | not measured | not measured |
| flux conservation | **2.2e-3 -> 2.0e-3** | **39 -> 9.8** | 40 -> 9.8 | not measured | not measured |
| the tolerance today | 1e-8 | 1e-8 | 1e-8 | 1e-8 | 1e-8 |
| the state stands at | 3.0e-4 -> 3.0e-4 | 5.6e-4 -> 5.8e-4 | ~1e-6 -> ~1e-6 | 7.2e-2 | 4.8e-1 |

Read down each column: the arithmetic floors are five to eleven decades below
1e-8 and forbid nothing; the derivative sits ON 1e-8 to within a factor 4;
the discretization and the flux conservation are four to six decades ABOVE
it. **`cert_tol_element = 1e-8` is not refused by any floor of the arithmetic
or of the closure.** [corrected 2026-09-11] The two sentences that followed
read "What refuses it is that the discrete equation itself only describes the
continuum equation to about 5e-4 on this grid at this radius, and that the
state's element fluxes are only conserved to 2e-3. A state certified at 1e-8
on this grid would be a state whose residual is four decades below the error
of the discretization it is a residual of." The first half of that is
withdrawn: the residual of the discrete equations and the accuracy of those
equations against the continuum are different quantities, and an algebraic
error decades below the discretization error is the normal and desirable
situation, not a contradiction (review section 3.1; MEASURED, the
frozen-background H2 wind row of 2.449021e-13 on the 500-cell grid,
`docs/solver_partition_experiment_20260911.md` section 4). What the table
does establish about 1e-8 is what the last two rows say: the derivative the
solve is steered by resolves only 3.8e-8 at the binding cell, and the states
this implementation produces stand at 3.0e-4 to 4.8e-1, four to eight
decades away. The flux-conservation row remains a statement about the
STATE handed back (its element fluxes are conserved to 2e-3 in the wind and
not at all in the layer), which is a property of that state and of the
layer's standing base wave, not of the tolerance. **N29
moves no column of this table by as much as a decade except the closure seed
(four decades down, and it was already the least binding of the floors) and
the layer flux conservation (a factor 4 down, and still nine decades above
the tolerance).**

## 5. What the two candidates ARE, in the certification's language

With the tolerances as they stand (nothing changed), on the post-N29 states
(MEASURED, N8c):

- **atomic candidate**: 11 of 13 active equations refuse it. Three
  hydrodynamic rows (mass 1.494e-4 at cell 106, momentum 2.227e-5 at cell
  500, energy 2.271e-4 at cell 108) and all eight element rows: He/H
  5.932e-4 (cell 245), Fe 5.922e-4 (cell 114), Ca 1.169e-4 (cell 131), Mg
  2.767e-5 (cell 132), O 2.218e-5 (cell 500), C 1.806e-5 (cell 499), Na
  1.739e-5 (cell 487), N 1.044e-5 (cell 133). The level balance (1.643e-19)
  and the closure (3.051e-13) are within tolerance. Zero cells without a
  chemical root, zero residual samples refused of 29948 admitted.
  The state's own mass closure is 1.95e-9, the precision of the file.
- **carrier candidate**: 4 entries refuse it, of which the carrier H2 row at
  4.800e-1 (cell 205, r = 1.098) is the species one; the other three are
  hydrodynamic (mass 4.007e-3 at cell 211, momentum 8.241e-2 at cell 499,
  energy 2.519e-1 at cell 51). Mass closure 1.93e-9. Byte-identical to the
  state N26 handed back.

The reading is unchanged from N8b in kind and in size. What N29 changed is
the atomic hydrodynamic rows: the momentum row of the written state falls
from 4.860e-1 as reloaded (1.812e-5 as written) to 2.227e-5 written and
2.227e-5 reloaded, because the reload now reconstructs the density the file
carries (section 6).

With the proposal of section 7 (element and carrier rows at 1e-5 in the wind,
layer rows reported and not gating), the atomic candidate would still be
refused by its element rows (5.932e-4 against 1e-5, and 3.048e-4 in the wind
window) and by its hydrodynamic rows; the carrier candidate would still be
refused by 7.2e-2 in the wind. **No proposal here would certify either
state, before or after N29.** That is the point of stating it: the anchoring
is not a way to admit the states that exist.

## 6. The reload was not the state (CLOSED by N29)

MEASURED, and before N29 it limited every number above that was taken
through a reload. Evaluating the handed-back atomic state through
`Load IC? True` gave:

| quantity | handed back (READ, N26) | as reloaded (MEASURED) |
|---|---|---|
| `\|\|R\|\|` | 3.719e-4 | **4.860e-1** |
| worst hydrodynamic row | momentum 1.812e-5, cell 499 | momentum 4.860e-1, cell 500 |
| element He/H row | 6.061e-4, cell 245 | 5.974e-4, cell 245 |

The element rows reproduced to 1.5 percent; the momentum row did not, by
four decades. The cause was section 3: the loader rebuilds rho from the
species and the species did not carry the file's rho, so the reloaded
density was up to 0.67 percent below the written one and the outermost
momentum row, a cancellation at the 1e-10 level in code units, was destroyed
by it.

**After N29** the same round trip on the post-N29 candidate (MEASURED, N8c):

| quantity | as the writer reported it | as reloaded |
|---|---|---|
| `\|\|R\|\|` | 2.271e-4 | **2.2714e-4** |
| hydrodynamic mass row | 1.4940801779e-4, cell 106 | 1.4940797199e-4, cell 106 |
| hydrodynamic momentum row | 2.2271068257e-5, cell 500 | 2.2271089825e-5, cell 500 |
| hydrodynamic energy row | 2.2713430987e-4, cell 108 | 2.2713617755e-4, cell 108 |
| element He/H row | 5.9318678537e-4, cell 245 | 5.9318678540e-4, cell 245 |

Every row returns on its own cell, the element row to eleven digits and the
hydrodynamic rows to six. **The anchoring above is therefore about the
candidate itself.** The carrier reload was never affected: its writer reports
`||R||` 2.519e-1 and its reload 2.5576e-1, on the energy row of a base cell
in both.

## 7. PROPOSAL, for the user's decision

This is item 22 for `docs/To_be_determined_by_user_20260906.md`, written here
because this item does not own that file. Nothing below is implemented.

> **Item 22. The species-row certification tolerances.** `cert_tol_element`
> and `cert_tol_carrier` stand at 1e-8, chosen from `newton_floor` of the
> composition Newton and never measured against a state. N8b measured the
> five anchors of `ISSUES_20260909_review` 5.4 (memo
> `docs/certification_tolerance_anchoring_20260910.md`) and N8c re-measured
> all five on the states the corrected element projection of N29 hands back.
> The floors of the arithmetic (0, 9e-17, 5e-19) forbid nothing; the
> derivative sits at 6.5e-9 to 3.9e-8; the discretization of the operator on
> the production grid is ~4.6e-4 at the binding radius and the element flux
> is conserved to 2e-3 in the wind and to 10 below r = 1.10.
>
> **(a) Adopt the anchored values.** `cert_tol_element` and
> `cert_tol_carrier` become **1e-5** for cells at r >= `cert_regime_wind_r`,
> a decade above the discretization error a smooth manufactured column shows
> at the production spacing (1.52e-6 after N29) and one and a half decades
> below what the candidate states carry in that window (3.048e-4 atomic,
> 7.19e-2 carrier), so nothing that exists is admitted by it. Cells below
> `cert_regime_layer_r` are REPORTED and do not gate until the base wave is
> removed, because their element flux conservation is 9.8 after N29 (39
> before it) and no tolerance is meaningful against that. Cost: a state can be certified in the wind while
> its layer is not certified, and the report has to say so in one line.
> [corrected 2026-09-11: the clause "a decade above the discretization error
> a smooth manufactured column shows at the production spacing" is a choice
> of where to put a gate given the ACCURACY the operator has at that
> spacing, and it is kept as that. It is not an argument that a smaller
> residual is unattainable there, which review section 3.1 refutes and the
> MEASURED 2.449021e-13 of
> `docs/solver_partition_experiment_20260911.md` section 4 contradicts
> directly. Decision 22 (a) was taken on this option and stands.]
> **(b) Keep 1e-8 and state that the fixtures cannot certify.** The value
> stays where it is, the report says that no state this implementation has
> produced reaches it at the radii that bind, and the species-row solves are
> documented as uncertified by every method measured so far.
> This is honest and it makes the certification a measurement rather than a
> gate for these rows.
> [corrected 2026-09-11: option (b) read "no state of this discretization can
> reach it at the radii that bind, and the species-row solves are documented as
> permanently uncertified until the operator's order or the grid changes".
> That inference is WITHDRAWN. Spatial truncation error is not a lower bound
> on the algebraic residual of the discrete equations (review section 3.1),
> so neither the operator's order nor the grid has to change for these rows
> to fall further, and MEASURED on the same 500-cell grid the existing H2
> transport routine took the wind H2 row from 7.380744e-2 to 2.449021e-13 on
> a held background (`docs/solver_partition_experiment_20260911.md` section
> 4). What is measured is that no state produced by the methods tried so far
> reaches 1e-8, or 1e-5, at those radii. The choice between (a) and (b) is
> unaffected: decision 22 (a) was taken and stands.]
> **(c) A two-regime tolerance, both values gating.** 1e-5 in the wind and a
> separate, larger number in the layer measured after the base wave is fixed.
> This is (a) with the layer promoted from reported to gating, and it cannot
> be written until that measurement exists.
>
> The advisor's reading: (a) or (b). (a) makes the gate say something true
> about the wind; (b) says nothing and claims nothing. (c) needs work that has
> not been done. What must NOT happen is a single loosened number covering
> both regimes, which is what review 5.4 warns against.
>
> Whatever is chosen, three things were named as prerequisites and are not
> tolerance questions. Two of them are now discharged: the mass closure of
> section 3 (6.7e-3 -> 1.9e-9) and the reload of section 6 (four decades ->
> eleven digits), both by N29. The third stands: the element flux spread in
> the layer, 39 before N29 and **9.8** after it, which is the standing base
> wave and not the tolerance.

**The numbers after N29 do not move this proposal.** Every anchor is within
a factor 2 of the value N8b measured except the closure seed, which falls
four decades and was already the least binding floor, and the layer flux
spread, which falls by a factor 4 and remains nine decades above any
candidate tolerance. Option (a) is supported exactly as it was written: 1e-5
gating in the wind, the layer reported. The layer's element flux spread, the
number the advisor asked for, is **9.79 for helium, 9.09 for hydrogen and
9.25 for the face mass flux** over r < 1.10 (MEASURED, N8c).

## 8. How to repeat any of it

- the regime split and the full-precision row measures: `EXHALE_CERT_ANCHOR=1`
  on any run; `EXHALE_CERT_ANCHOR_R=<r>` adds every row's reading at the cell
  nearest that radius, which is how one physical radius is followed across
  two grids. Both are off by default and change no other line of the report.
- the derivative curve at one cell: `EXHALE_SPECIES_JAC_TEST=1` with
  `EXHALE_JAC_TEST_CELL=<j>`.
- the five numbers on a synthetic column, with no planet run:
  `src/tests/certification/run.sh`, the rows of
  `the_five_anchoring_numbers_of_a_synthetic_column`.
- the two candidate states themselves: the atomic one from
  `backup/regression/atomic_elem_newton` with `IC/` copied into `output/`,
  `Load IC? True` and `OMP_NUM_THREADS=8 EXHALE_JFNK_MAXIT=300`; the carrier
  one from the N26 carrier fixture at a cap of 150. Both take about an hour
  on a loaded machine.
- the mass closure of any written state, from the files alone:
  `sum_i m_i n_i` against the `rho` column, helium at 3.9715259 m_H, a metal
  stage at its `melem_A`, a molecule at the nuclei it carries and the
  He 2^3S column inside He I.

## 9. After N29: the two candidates produced again

Item N8c, 2026-09-10. Everything in this section is MEASURED with a private
build of the tree that carries N29 (`make OBJDIR=build_n8c EXE=EXHALE_n8c.x`,
conda-forge gfortran 16.2.0, md5 `a9627a14`). No tolerance and no code value
was moved by N8c; the certification suite is 42 rows, 0 FAIL.

### The two states

| | atomic element candidate | carrier candidate |
|---|---|---|
| recipe | `backup/regression/atomic_elem_newton/IC`, `Load IC? True`, `EXHALE_JFNK_MAXIT=300`, 8 threads | the N26 carrier fixture, cap 150, 8 threads |
| outer iterations | 275 (N26: 252) | 63 (N26: 63) |
| trust-region steps accepted / rejected | 254 / 22 (N26: 235 / 18) | 53 / 11 (N26: 53 / 11) |
| model refused | 0 | 0 |
| stop | stagnation detector, `info = 2` | stagnation detector, `info = 2` |
| `\|\|R\|\|` handed back | 2.271e-4 (N26: 3.719e-4) | 2.519e-1 (N26: 2.519e-1) |
| binding row | element He/H of cell 245, r = 1.153 (N26: element He of cell 246) | carrier H2 of cell 205, r = 1.098 (N26: the same) |
| entries that refuse | 11 of 13 (N26: 11) | 4 of 6 (N26: 4) |
| residual samples admitted / refused | 29948 / 0 on all six names | 4406 / 0 on all six names |
| trial or probe states refused by an admissibility screen | 0 | 0 |
| flux spread at the gate window | 9.849e-7 (N26: 5.497e-7) | 4.629e-2 (N26: 4.629e-2) |
| flux spread r >= 1.03 / r >= 1.10 | 3.554e-4 / 3.297e-6 | 9.339e-2 / 6.662e-2 |
| mass closure of the state | **1.95e-9** (N26: 6.68e-3) | **1.93e-9** (N26: 1.93e-9) |
| `log10 Mdot` | 10.10 | 10.60 |

The carrier state is **byte-identical** to the one N26 handed back, over
every numeric line of `Hydro_ioniz.txt` and `Ion_species.txt`. That is what
the code says it should be: the case carries no `metals.inp` and no element
row, so the projection N29 corrected is never reached. Every carrier number
in this memo therefore stands unchanged, and it was re-measured rather than
assumed.

The atomic state is a different state, as the ulp-level sensitivity of this
reload (N26c) guarantees it must be after a change of this size. It is a better
one on every reading that is not a tolerance: `||R||` 1.6 times lower, the
mass closure six decades lower, the reload returning the writer's rows.

### What moved and what did not

- **closure of the state**: 6.68e-3 -> 1.95e-9, the precision of the written
  file. This is the finding of section 3 and N29 closed it.
- **the reload**: `||R||` 3.719e-4 written against 4.860e-1 reloaded, on
  different rows, becomes 2.271e-4 against 2.2714e-4 on the same rows.
  Section 6 is closed.
- **the closure seed the state's own sweep carries**: 1.736e-10 -> 1.645e-14,
  so the closure floor of the element rows falls from 5.3e-15 to 5.0e-19.
- **the layer element flux spread**: 39.3 -> 9.79. A factor 4, in the window
  where the composition changes most, which is where the leak was.
- **everything else is within a factor 2**: the row measures themselves
  (He/H 5.974e-4 -> 5.932e-4), the perturbation amplifications, the
  derivative plateau (6.41e-5 -> 6.62e-5), the grid-refinement ratio (2.41 ->
  2.402) and its Richardson limit (1.32e-4 -> 1.31e-4), the manufactured
  column at three grids and its observed order 1.589, and the wind and escape
  flux spreads.

### The reading

Unchanged. The arithmetic floors are now five to fifteen decades below 1e-8
instead of five to eleven, so they forbid even less; the derivative still
sits on 1e-8 within a factor 4; the discretization at the binding radius is
still ~4.6e-4 and the element fluxes are still conserved only to 2e-3 in the
wind and to 10 in the layer. Neither 1e-8 nor 1e-5 is reached by either
candidate, and option (a) of section 7 still admits neither.
[corrected 2026-09-11: the last sentence read "`cert_tol_element = 1e-8` is
still refused by the discretization and the flux conservation and by nothing
else". The discretization does not refuse a tolerance on the algebraic
residual (review section 3.1; MEASURED, 2.449021e-13 for the
frozen-background H2 wind row on this grid,
`docs/solver_partition_experiment_20260911.md` section 4). What refuses
these two candidates is the states themselves: they stand at 3.0e-4 and
4.8e-1, and the derivative the solve is steered by resolves 3.8e-8 at the
binding cell.]

## Anchor (6), 2026-09-11: the mass row's rounding floor

`PLAN_20260911_partitioned_solver.md` item P16, on the user's decision of
2026-09-11: re-anchor the hydrodynamic mass tolerance ABOVE the measured
rounding floor of the state, not by a new hand constant. What P7 and P10
left open was the last refusing entry of the HD 209458 b element reload,
the mass row at 8.893E-10 against a fixed `cert_tol_mass = 3e-12`, where
N33 had already measured that a base-layer mass row carries no signal at
that level.

Every number below is MEASURED by this item unless it says READ.

### What the row's rounding floor is

The continuity row of cell `j` is a difference of two face fluxes over the
cell volume,

```
R_1(j) = ( A_+ F_+  -  A_- F_- ) / dV_j ,     s_1(j) = max|A F| / dV_j ,
```

and the certification judges `|R_1|/s_1`, the fractional change of the mass
flux across the cell. Two candidate estimates of the smallest such fraction
the arithmetic can resolve were measured against the row's actual
non-smoothness.

- **the flux-difference estimate**, the textbook rounding of a difference of
  two correctly rounded addends, `eps (|A_+ F_+| + |A_- F_-|) / dV_j` over
  `s_1`. Because `s_1` IS the larger addend, this is algebraically between
  one and two ulps at every cell of every state, and MEASURED it is
  4.4409E-16 to within a part in a thousand throughout all three fixtures.
  It is the expression the item was briefed with, and it cannot serve: it is
  a constant, so `c_round` times it never rises above 3e-12 for any margin
  the measurement supports, and the step the assembly actually takes stands
  up to 2.4E+05 times above it.
- **the signal estimate**, `eps rho(|v| + c_s) A / dV_j` over `s_1`, which
  is `eps` times `max_faces[rho(|v|+c_s) A] / max_faces|A F|` and so about
  `eps/Mach` where the flux is a small residue of the momentum the cell
  carries at its own signal speed. This is the one adopted
  (`mass_row_rounding_floor`, `steady_residual.f90`). Its justification is
  N33's measurement: the interface flux is assembled from the JUMPS of the
  two reconstructed states, those jumps stand 3.2e5 to 5.9e6 below the
  states in the base layer (READ, N33), and the first quantity of the
  assembly that steps under a one-ulp perturbation is the interface flux, by
  4.803E-16 of an O(1) state where its own value is smaller by that factor
  (READ, N33). So the flux carries the last bit of `rho c_s` and not of
  `rho v`.

### The measurement of `c_round`

`EXHALE_MASS_FLOOR_SCAN=1` (default off, `steady_residual.f90`) adds one ulp
to the density of every physical cell, runs the whole flux assembly again,
and reports the STEP each mass row takes against both estimates. The state
as given is assembled once more afterwards and reproduces its own mass rows
to 0.00E+00 on all three fixtures, so the step is the perturbation's and not
the measurement's.

Each fixture is entered with `Restart intent: stationary evaluate`, which
measures the state as loaded and takes no step.

| fixture | base Mach | signal estimate at cell 1 | flux estimate | largest step / signal estimate | largest step / flux estimate |
| --- | --- | --- | --- | --- | --- |
| `wasp_full_newton` | 3.937E-03 | 5.7049E-14 | 4.4409E-16 | **1.5931** (cell 207) | 1.109E+02 (cell 52) |
| `atomic_elem_newton` (HD 209458 b) | 5.062E-06 | 3.2031E-11 | 4.3701E-16 | **1.8857** (cell 82) | 2.4365E+05 (cell 152) |
| `carrier_elem_newton` (hot Uranus) | 2.400E-04 | 9.3076E-13 | 4.4275E-16 | **2.6806** (cell 17) | 5.0915E+03 (cell 17) |
| the HD 209458 b state this item certifies (below) | 2.767E-06 | 1.0041E-09 | 4.4409E-16 | **1.7733** (cell 212) | 1.4201E+06 (cell 28) |

The signal estimate tracks the measured step within a factor 1.6 to 2.7 on
states whose base Mach numbers span three decades and whose floors span five
(5.7E-14 to 1.0E-09); the flux estimate is a constant and is short by two to
six decades.

`c_round = 10` (`cert_mass_round_margin`, `certification.f90`): a round
number 3.7 times above the largest ratio measured here, and inside N33's
independent reading of the same floor against its own bound, which was 5.5
to 12 (READ). Nothing above that is measured and nothing above that is
taken.

### The tolerance

```
tol_mass(j) = max( 3e-12, min( 1.0, c_round * floor(j) ) )
```

`cert_tol_mass = 3e-12` stays, as the FLOOR of the tolerance: where the
cell's arithmetic permits 3e-12 the row is held to it. The ceiling of 1.0 is
the row measure at which the mass flux changes by the whole of itself across
the cell, which no arithmetic argument admits; it is not an operative branch
on any state measured (it needs a base Mach below 1e-15). The verdict on the
mass row is taken CELL BY CELL, because the cell whose measure stands
furthest outside its own tolerance is not in general the cell of the largest
measure; `row_max` and its cell are untouched, and the report states which
cell took the verdict, against which number, and whether that number was the
fixed value or the rounding anchor.

What is given up where the anchor binds: the mass row of that cell is judged
at `eps/Mach` of its own flux instead of at 3e-12, which on the HD 209458 b
base layer is about 3e-10. What the row's scale was introduced for is
undamaged, the flux errors it was given to catch being fractions of a
percent to tens of percent (the 30 percent at 1.03 R_p of
`docs/p54_base_layer_mass_flux.md`), seven decades above any floor measured
here, and the mass-flux spread gate reads the same quantity independently.

### Impact, MEASURED against a private build of the entry text of the two files

Three runs, 8 threads, 500 cells, reloaded from each fixture's `IC/`, with
`Restart intent: stationary`. The control is a build of the entry text of
`certification.f90` and `steady_residual.f90` in the same tree; the machine
was shared, so wall times are reported and are not a claim.

**1. `wasp_full_newton` reloaded.** Unmoved. `info = 0`, CERTIFIED, mass row
1.212E-12 at cell 500 with the verdict taken there against the FIXED 3.0E-12
at distance 0.404, momentum 7.055E-10, energy 1.920E-09 identical to the
control's digit for digit, and every written output byte-identical apart from
the provenance timestamp of `Hydro_ioniz.txt` (one second apart).
`log10 Mdot` 13.30 on both, from the written profile under the run's
`Rate/2 + Mdot/2` convention. The anchor cannot reach this case: its mass
row binds in the WIND, where the estimated floor is below 3e-13, and
`cert_tol_mass` is the floor of the tolerance.

The cold-start `wasp_full_newton` of `make check` is unaffected for the same
reason and for one more: READ from `backup/regression/wasp_full_newton/run.log`,
that run ends `-> converged: JFNK steady solution (||R|| < 1.000E-08 AND flux
spread < 2.000E-05)` with a final mass row of 2.624E-13 at cell 500 against
3.0E-12, every hydrodynamic row within its tolerance and the state CERTIFIED,
and the log carries no `the acceptance gate is met and the certified ... is
above its own ...` line at any iteration. So the mass row never bound that
stop; a tolerance that only ever RISES cannot move a stop the mass row did
not hold.

**2. The HD 209458 b element reload, partitioned** (`Load IC? True`,
`Coupled carrier solve: False`, `Restart intent: stationary`, `metals.inp`,
`EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80 EXHALE_DIFF_OMEGA=0.5`).
**The state CERTIFIES, at pass 12**, which is the entry the plan's "what
remains" list had open:

```
(EXHALE_main) outer pass 12: ACCEPTED -- every active equation of this
                             state is within its own tolerance.
```

The certification of that state: mass row 7.703E-10, largest at cell 2, the
verdict taken at cell 34 (r = 1.00665) where 6.542E-10 stands against
**4.6E-09, the rounding anchor of that cell**, at distance 0.141; momentum
1.562E-14 against 1E-08; energy 9.981E-09 against 1E-06; every elemental
wind row within 1E-05, so the state is certified IN THE WIND with the layer
rows reported (worst 8.966E-05, Fe at cell 145, r = 1.043). Twelve passes,
1047 s as the sum of the pass timings, on a machine also running a test
suite.

The mass row of each pass, and which gate took its verdict (MEASURED, the
`verdict at cell` line of each pass's certification):

| pass | mass row (cell) | verdict cell, r | measure there | tolerance | gate | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| as loaded | 1.509E+00 (191) | 285, 1.2477 | 1.936E-01 | 3.1E-12 | anchor | ABOVE |
| 1 | 4.106E-09 (1) | 117, 1.0294 | 1.904E-10 | 8.1E-10 | anchor | within |
| 2 | 1.322E-09 (39) | 117, 1.0294 | 1.904E-10 | 8.1E-10 | anchor | within |
| 3 | 5.484E-09 (12) | 160, 1.0530 | 2.270E-09 | 3.2E-11 | anchor | ABOVE |
| 4 | 8.407E-10 (34) | 221, 1.1145 | 4.764E-13 | 3.0E-12 | fixed | within |
| 5 | 7.019E-10 (20) | 97, 1.0219 | 1.847E-10 | 7.8E-10 | anchor | within |
| 6 | 6.734E-10 (19) | 66, 1.0132 | 3.049E-10 | 2.3E-09 | anchor | within |
| 7 | 9.305E-10 (35) | 128, 1.0344 | 6.888E-11 | 2.4E-10 | anchor | within |
| 8 | 9.241E-10 (5) | 79, 1.0165 | 3.692E-10 | 1.5E-09 | anchor | within |
| 9 | 7.807E-10 (3) | 92, 1.0203 | 1.625E-10 | 1.0E-09 | anchor | within |
| 10 | 7.050E-10 (1) | 97, 1.0219 | 1.497E-10 | 8.6E-10 | anchor | within |
| 11 | 8.147E-10 (15) | 97, 1.0219 | 1.655E-10 | 8.5E-10 | anchor | within |
| 12 | 7.703E-10 (2) | 34, 1.0067 | 6.542E-10 | 4.6E-09 | anchor | within |

Three things this says. The row's largest measure sits in the base layer at
1e-9 to 5e-9 throughout and is admitted there, because the base cells' own
tolerances are larger still; the VERDICT then moves out to r = 1.006 to
1.115, where the anchor is 3e-12 to 5e-9 depending on the cell, so the gate
has not been removed, only moved to where the arithmetic still resolves the
row; and the row is refused at one of the twelve passes (pass 3, distance
69.9 at cell 160), so the anchored tolerance still refuses states.

**Does 8.893E-10 at cell 3 fall inside the tolerance there?** Yes, by a
factor of 10.9. The certified state was re-entered with `Restart intent:
stationary evaluate` and the scan armed: at cell 3 its estimated floor is
9.6526E-10, so `tol_mass(3)` is 9.65E-09, and the row that state itself
carries at cell 3 is 8.9806E-10, which is the same quantity P10 reported at
8.893E-10. The reloaded certified state certifies again, with its mass row
largest at cell 1 (1.650E-09) whose own tolerance is 1.0E-08 and the verdict
taken at cell 18 (1.143E-09 against 6.7E-09, distance 0.170).

The tolerance MOVES with the state, which is the point of anchoring it
there: the base Mach number of a converged column is not the base Mach
number of the reloaded one. At cell 3 it is 7.867E-06 on the state as loaded
and 1.190E-06 on the certified one, and the floor there accordingly rises
from 3.0E-11 to 9.7E-10 (MEASURED).

The control, a build of the entry text of the same two files run on the same
recipe, 12 passes and 1040.83 s: **NOT CERTIFIED, one entry refusing, the
hydrodynamic mass row at 7.159E-10 against 3.0E-12 at cell 13**. Everything
else about the two runs agrees: the worst gated species row of every one of
the twelve passes is identical to three digits and so is its cell, and the
eight elemental rows of the last pass agree to every digit printed (He/H
2.793E-05 at cell 238, C 1.000E-05, O 1.238E-05, N 1.496E-05, Mg 4.377E-05,
Ca 1.942E-05, Na 1.941E-05, Fe 8.966E-05). Passes 1 to 5 are identical in
every reported number in the two runs; they part at pass 6, the first pass
whose hydrodynamic solve stops on the loop-top test because the mass row is
inside its own tolerance (`info = 0` where the control's is `info = 2`).

The wall times of the two runs (1047 s against 1041 s) are not a comparison:
the P16 run shared the machine with a test suite and the control did not.
The two other passes that end `info = 0` in the new run are 6 and 12.

Note that this control is the CURRENT tree's entry text and not P10's binary,
so its ladder is not P10's ladder number for number (P10 read a mass row of
8.893E-10 at cell 3 after twelve passes; this control reads 7.159E-10 at cell
13). Items P11 and P12 and two advisor changes landed in between. What the
two controls agree on is the finding: after twelve passes the ONE entry
refusing the state is the hydrodynamic mass row of a base-layer cell.

**3. The hot-Uranus carrier reload, partitioned**, `EXHALE_OUTER_PASSES=3
EXHALE_JFNK_MAXIT=40`. The mass row's verdict per pass, and what it costs:

| pass | hydro info | mass row (cell) | verdict cell, r | tolerance, gate | verdict | seconds |
| --- | --- | --- | --- | --- | --- | --- |
| control 1 | 1 | 3.375E-12 (80) | one number for the column | ABOVE 3.0E-12 | ABOVE | 28.42 |
| P16 1 | **0** | 6.944E-12 (1) | 500, 4.6215 | 3.0E-12, fixed | **within** | **8.75** |
| control 2 | 1 | 1.370E-11 (1) | one number for the column | ABOVE 3.0E-12 | ABOVE | 20.50 |
| P16 2 | 2 | 1.256E-11 (1) | 180, 1.0707 | 3.0E-12, fixed | ABOVE | 3.13 |
| control 3 | 1 | 3.645E-12 (72) | one number for the column | ABOVE 3.0E-12 | ABOVE | 20.47 |
| P16 3 | 2 | 5.519E-12 (19) | 180, 1.0707 | 3.0E-12, fixed | ABOVE | 3.59 |

The base cells of this fixture are admitted by their own anchor at every
pass (cell 1 at 6.944E-12 and 1.256E-11 against a floor of 9.3E-13 times
ten), which is what N33's reading of this row predicted, and the verdict
moves out to **cell 180, r = 1.0707**, where the estimated floor is below
3e-12 and the row stands at 5.1E-12 to 9.5E-12 against the FIXED tolerance.
So the mass row still refuses this fixture and now refuses it at a cell
where the arithmetic resolves the row: the refusal has moved from rounding
to signal.

The H2 carrier row is unchanged pass by pass, as it must be, the carrier
relaxation sitting outside the hydrodynamic solve: 7.38E-02 at cell 290,
7.22E-02 at 289, 7.06E-02 at 288 in both runs. It still refuses the state
(2 entries refuse the last pass in both runs). The three passes cost 15.47 s
against the control's 69.39 s, because each hydrodynamic solve now stops when
its rows are within their tolerances instead of spending its budget on a
mass row below the rounding of the base layer.

### What the code carries now

- `steady_residual.f90`: `mass_row_rounding_floor(j, u)`, the estimate, public
  and the only expression of it; `mass_row_rounding_floor_scan`, the
  measurement, armed by `EXHALE_MASS_FLOOR_SCAN` and off by default.
- `certification.f90`: `cert_mass_round_margin = 10` and
  `cert_tol_mass_ceiling = 1`; `cert_tol_mass_at(j, u)` and the pure
  `mass_row_cell_verdict(q, floor, tol, dist, within, anchored)`, which is
  the rule and the one thing the tests drive; `mass_row_verdict`, the loop
  over the column that takes the verdict and records where it bound;
  `cert_mass_gate_name`, the two words the report prints, read by the report
  and by the test from one place. `cert_tol_mass` keeps its value and is now
  read only through the accessor.
- The report gains, under the mass row: the cell that took the verdict with
  its radius, the measure there, the tolerance there, the distance, and which
  of the two gates applied; and, when the largest measure sits elsewhere, the
  tolerance that applied at that cell. The refusing-entry list names the cell
  that refused rather than the cell of the largest measure, which on this row
  may be inside its own tolerance.
- Tests, `src/tests/certification/`: the block
  `the_mass_row_tolerance_is_anchored_on_its_rounding_floor`, 11 rows on the
  floors MEASURED here (5.7E-14, 9.3E-13, 3.0E-11), which assert that a floor
  below the fixed value leaves the fixed value, that a floor above it gives
  `c_round` times the floor, that a row at half the anchor is within where
  the superseded fixed rule (transcribed in the test) refuses it, that a row
  above the anchor still refuses, that the ceiling binds, that the report
  names which gate applied, and that the margin stands above every ratio
  measured. 57/0 at entry, **70/0** after, and 0 PASS / 1 FAIL on the control
  objects, where the compile names `cert_mass_round_margin`,
  `cert_tol_mass_ceiling` and `cert_mass_gate_name` as absent.

Reported by this item and CLOSED the same day by items P17 and P18 (READ,
`docs/Update_EXHALE_stage2.md` section 8): at the time of writing the best-iterate
ledger still formed its hydrodynamic distance from the three row maxima over
the fixed tolerances (`rc(1)/cert_tol_mass`), so it ranked the mass row
against 3e-12 while the stop test and the acceptance read the cell's own
tolerance. P17 made the ledger read
`hydrodynamic_distance_from_certification_by_cell`, the maximum over cells
and rows of `|R_kj|/scale_kj/tol_k(j)` with `tol_1(j)` this anchor, and P18
removed the fixed-tolerance form, so one expression of the distance exists.

### How to repeat any of it

- the rounding floor of the mass row on any state, with both estimates and
  the measured step beside them: `EXHALE_MASS_FLOOR_SCAN=1`, which fires
  once, on the first state a run assembles a stationary residual for. Off by
  default; it restores the row terms it found and decides nothing.
- the floor of a state already written: re-enter its directory with
  `Restart intent: stationary evaluate`, which measures the state as loaded
  and takes no step, with the scan armed.
- the tolerance rule itself, with no planet run:
  `src/tests/certification/run.sh`, the rows of
  `the_mass_row_tolerance_is_anchored_on_its_rounding_floor`.
- the three impact runs: the recipes in the dated sections of
  `backup/regression/atomic_elem_newton/README.md` and
  `backup/regression/carrier_elem_newton/README.md`, and
  `wasp_full_newton` with `Load IC? True` and `Restart intent: stationary`
  added.
