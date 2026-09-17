# LHS 1140 b stationary series, item L22 step 3 (iii): the design of the coupled block

Written 2026-09-16, closed 2026-09-17 (KST). Item L22 of `docs/PLAN_20260916_rev3.md` section 3
names three algorithms in its step 3 and keeps them apart: (i) cell-local
pseudo-time preconditioning of the unchanged carrier residual, (ii)
conservative physical subcycling, and (iii) a coupled block. The user approved
this document on 2026-09-16 ("fix the ghost-rule inconsistency first, then
design the coupled block"). It is a design and nothing else: **no source file
is edited here, nothing is built and nothing is run.** Every number below is
READ, from the memos named beside it or from the source at the line cited.
Nothing was re-measured for this document.

The source identities this design was read against (READ, `md5sum` at the time
of writing, 2026-09-16): `src/EXHALE_main.f90` `8ab093e892d98077`,
`src/modules/time_step/steady_newton.f90` `33cbd5065c7eea75`,
`src/modules/lower_atmosphere/diffusive_photochemistry.f90` `4f9687fc1c64d045`,
`src/modules/functions/binary_element_diffusion.f90` `8ce3a0244d0715719`. Other
workers are editing `diffusive_photochemistry.f90` in parallel (L22 step 2c,
the outer ghost rule), so its line numbers move; every statement about it is
therefore made by ROUTINE NAME, and a line number beside a routine name is the
line that routine begins at in the identity above.

**One outer ghost rule is assumed.** Step 2b section 1 of
`docs/lhs1140b_stationary_L22_20260916.md` measured that the certification's
outer ghost is DATA (the ionization sweep solves the ghost cells,
`ionization_equilibrium.f90` loops `1-Ng, N+Ng`) while every line-search trial
inside `solve_carriers` refills `ftry(N+1:N+Ng,ic) = ftry(N,ic)`, a COPY, and
that the copy reduces the limited slope of cell N to exactly zero
(`the_copied_outer_ghost_leaves_no_slope_on_the_outflow_face` = 0.0, MEASURED
there) while a ghost that is not a copy leaves a slope of 3.62e-02 of the cell
value. L22 step 2c is fixing that. **This design assumes the outcome: ONE
outer ghost rule shared by the carrier operator, by the certification that
judges it and by the relaxation that advances it.** Where the coupled block
below evaluates a carrier row it evaluates the same rule; if step 2c ends with
two rules the coupling check of section 7 increment I1 will fail on the last
two cells and that failure is the correct signal, not a tolerance to widen.

---

## 1. The problem statement, from the measurements

**The alternation fails because the composition update the carriers need moves
the base thermodynamics out of the range in which the wind held at that
composition is a useful approximation, and the wind's response in that layer is
not small.** The measurements that say so, in order.

**1.1 Every pass ends on the bound, and the cell that sets the bound is not the
cell that refuses.** READ, `docs/lhs1140b_stationary_L7e_20260915.md` section
23.3: on both refusing well-mixed states every carrier relaxation of every pass
ended on the composition movement bound and not one ended on its own residual;
the bound was cut 5.0e-03, 2.5e-03, 1.25e-03, then the floor 1.0e-03, and the
composition movement of a pass fell in exact proportion (1.67e-02, 8.41e-03,
4.21e-03, 2.10e-03, 1.68e-03). The cell that ATTAINS the bound is the H2 front
(cell 231 at 1.2565 R_p with row measure 6.6e-03 on `wellmixed/HeH0.083`; cell
225 at 1.2309 R_p with 5.1e-03 on `wellmixed/HeH0.55`); the cell that REFUSES
the certification is four to five times worse and decades away in abundance
(cell 500 at 29.0 R_p, 2.46e-02; cell 306 at 1.95 R_p with x2 = 2.7e-06,
2.48e-02). The gate those are refused against is
`cert_tol_carrier_wind` = 1.0e-05 at `r >= cert_regime_wind_r` = 1.20 R_p
(READ, `certification.f90` lines 248 and 333).

**1.2 The refusal is not the rows and not the Jacobian.** READ, L22 step 2b.
The two cell-500 rows are a one-sided diffusive influx balanced against an
advective divergence, with the chemistry 2.9 and 11.9 percent of the residual,
and every intervention that STATES the outer diffusive face makes the row worse
by 1.12x to 20.6x on `wellmixed/HeH0.083` and 2.7x to 182x on
`kzz1e9/HeH0.083`; cell 306 of `HeH0.55` is an interior cell whose row is the
residue of a near cancellation of two live diffusive faces against the
He+ + H2 destruction, and the boundary does not move it by one digit. The
deferred reconstruction terms of the carrier Jacobian are measured harmless:
the assembled action reproduces a central difference of the full carrier
residual to a mean relative error of 1.542e-03 over the direction's support,
the pass assembled with the reconstruction terms is the control pass to every
printed digit (17 kept steps, 22 refused trials, kept displacement 2.904e-02,
the same ending and the same carrier steady residual), and every transport
substep of the control already ends at a row residual between 1e-13 and 1e-17
after 2 to 5 Newton iterations with every interval covered in one attempt, so
the Jacobian is far better than the iteration needs. The one entry the carrier
relaxation's block-tridiagonal matrix cannot hold is
`d res(314)/d f_c(312)`, two cells off the diagonal, where the assembled action
is exactly zero and the operator's is 1.679e-03; over the whole column the
largest advective entry outside that band is 9.435e-02 of the band of its own
row. Section 4 returns to that entry, because the coupled block's band does
hold it.

**1.3 What the bound is for, and what happens when it is lifted.** The bound
exists because the wind is held fixed while the carriers relax, so a pass may
move the composition only as far as the fixed wind remains a useful background;
since item L7e section 17 it is written on the relative change of the cell's
PARTICLE COUNT `n_tot + n_e` against the composition the pass was entered with
(`carrier_particle_count_change`), which is the quantity the wind feels through
the pressure at the conserved thermal energy. L22 step 2 waived it over a
selection set on the frozen `wellmixed/HeH0.55` state. MEASURED there (step 2b
section 6, `EXHALE_L22B_DISPLACEMENT=1`):

| what the pass handed back | control | the set below 5 R_p waived |
|---|---|---|
| kept transport steps | 17 | 70 |
| carrier displacement of the pass | 2.904e-02 | 1.599e-01 |
| max abs(dp)/p | 8.862e-03 at cell 230, r = 1.2521 | 5.016e-02 at cell 10, r = 1.0019 |
| max abs(dT)/T | 4.971e-03 at cell 275, r = 1.5539 | 1.870e-02 at cell 1, r = 1.0002 |
| max abs(dmbar)/mbar | 9.804e-03 at cell 217, r = 1.2007 | 6.448e-02 at cell 3, r = 1.0006 |
| max abs(d(n_tot+n_e))/(n_tot+n_e) | 1.000e-02 at cell 221, r = 1.2153 | 6.893e-02 at cell 3, r = 1.0006 |

The control attains the bound EXACTLY (1.000e-02 against `trust` = 1.0e-02);
the waived pass stands at 6.9 times it, and the cell that carries the
displacement moves from the front (1.20 to 1.55 R_p) to the BASE (1.0006 to
1.0019 R_p). With the front free of the bound the carrier relaxation reached
its OWN fixed point for the first time on this route: 70 kept steps, `grow`
from 1.0 to 1.4e+12, not one trial refused, carrier steady residual 7.32e-10
volume-weighted with a worst cell of 1.05e-09 at cell 309, against 7.01e-03 and
2.62e-01 in the control. The base H2 partition fell from 0.934 to 0.784 over
those 70 steps and x(H2) at 1.2007, 1.3578 and 1.6044 R_p was multiplied by
1.103, 1.855 and 8.833 (READ, the step 2 log of that run).

**1.4 And the wind could not follow it.** The pass-2 hydrodynamic solve of that
run stood at `||R||` = 1.118 after 205 JFNK iterations with the line-search
damping at 9.5e-07 and the pseudo-time step at 7.3e-05; from iteration 100 the
residual norm oscillated between 1.100 and 1.162 with the line-search factor
alternating 1.0 and 9.54e-07, all 500 MASS rows outside their tolerance with
the worst at cell 123, r = 1.038, the momentum row with 483 to 492 cells
outside and the energy row with all 500. That is a limit cycle and not a slow
descent (READ, step 2b section 6), and it sits in the same inner layer, r near
1.0 to 1.04, that the composition moved by 5 to 7 percent.

**1.5 The statement of the problem.** Put the four measurements together. At
the bound the carrier block does not close (1.1, 1.2: the refusing rows stay at
2.4e-02 to 2.5e-02 against 1.0e-05 and the bound is cut to its floor while the
worst row does not fall); off the bound the carrier block does close at a fixed
wind (1.3: 1.05e-09), and the composition it closes at is one whose base p, T
and mean particle mass stand 5.0e-02, 1.9e-02 and 6.4e-02 away from the state
the hydrodynamic rows were last solved at, seven times the bound; and the wind
that has to carry that base does not converge (1.4). **The alternation has no
step between those two: any composition good enough for the carrier rows is
already outside the range in which holding the wind is admissible, and any
composition inside that range leaves the carrier rows where they are.** That
is a property of the splitting, not of the bound's value, and neither a
different bound nor a different pseudo-time distribution inside the carrier
half addresses it: L22 step 2 measured that relocating the bound's controlling
cell only relocates it to the edge of the mask (in five of six comparisons the
cell that took over is the first cell of the selection set or its immediate
neighbor) and that the refusing row RISES in every run and rises MORE where the
bound was waived.

What the measurements do NOT establish, and this design does not claim, is that
a coupled block is the only remedy. Step 2b section 6 states it: a shorter
alternation with the same bound on the inner layer would move the same way.
What the coupled block is, is the one algorithm of the three in step 3 that
addresses the splitting itself rather than the length of a half step.

---

## 2. The unknowns of the coupled block

### 2.1 The decision

**The block's unknowns are the three conserved hydrodynamic variables and every
registered transported species balance, in every cell of the column, solved as
one nonlinear system; the block IS the existing `Coupled carrier solve: True`
route of `solve_steady_jfnk`, not a new solver.**

This is the first finding of the source reading and it changes what the item
has to write. READ, `steady_newton.f90` `set_transported_species_rows` (line
2301): with the key on, the registry adds one row and one unknown for the
helium/hydrogen partition where `he_diffusion .and. thereis_He`, one for each
transported trace element that has a reservoir, and one for each solved carrier
where `thereis_mol .and. carrier_transport`; then
`nvar_jac = 3 + nspec_row` (line 2374). `eval_residual` (line 3815) assembles
the hydrodynamic rows through `assemble_residual`, the carrier rows through
`carrier_steady_residual` (the same routine, in the same module, that the
certification and the alternation's progress measure read; called at line 4321)
converted to code time by `tscale_code = R0/(v0*n0)`, and the elemental rows
through `element_transport_residual` (line 4340) converted by
`elem_he_to_code()` and `elem_tr_to_code()`, all into ONE vector `Fvec`. The
species unknowns are packed from and written back into the composition by
`pack_species_rows` (line 3453) and `write_species_rows_into_composition`.

For the three refusing cases the registry is short. READ,
`LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.083/input.inp`:
`Molecular chemistry: True`, `Molecular carrier transport: True`,
`Well balanced: True`, `Reconstruction scheme: PLM`, and NO `He_diffusion`, no
metals and no ionization transport. So `nspec_row` = 1 (the H2 carrier),
`nvar_jac` = 4, and with N = 500 the block has **2000 unknowns, four per
cell: rho, rho v, E and n(H2)**.

### 2.2 Why the whole column and not a band

The alternative the plan asks to be discussed is a separate base-layer solve
with the wind above held: take the base out of the alternation, since 1.3 and
1.4 put the displacement and the failure both at r < 1.04. **Rejected**, for
three reasons, each of them a measurement or a code fact and not a preference.

1. **The wind's response is not local** (review row R5, accepted in rev3
   section 0.1). The pressure of a cell couples to its neighbors through the
   momentum row, and the radiation field of a cell couples to every cell below
   it through the column density, which is why the coupled route's Jacobian
   ACTION is a difference of the FULL residual and not of the banded model:
   `jv_product` (line 10514) states it, "captures the non-local radiation
   coupling the banded Jacobian omits". A band of cells with the wind above
   held would need a boundary condition at its outer edge that the code does
   not have and that no measurement supports; the one boundary the code does
   have at a moving interface, the base, took items L21 and L26 to state.
2. **The failure is not confined to the band a base solve would carry.** In
   1.4 all 500 mass rows are outside tolerance and the energy row is outside in
   all 500; the worst cell is at r = 1.038 but the column is not close
   anywhere. A band solve would return a converged band inside a column that
   still refuses.
3. **The measured displacement is not confined to the base either.** In the
   CONTROL column of the table in 1.3 the displacement maxima sit at cells 217
   to 275, r = 1.20 to 1.55, which is the front and not the base. The base
   carries the displacement only once the bound is off. A base-only block would
   be designed on one of the two regimes.

The band idea is not worthless; it is the natural SECOND experiment if the
whole-column block converges too slowly, and section 8 keeps it as a named
fallback measured against a cost, not as a design.

### 2.3 Why the element fractions stay out of these three cases, and stay in generally

The three refusing cases do not set `He_diffusion`, so no elemental balance is
transported in them and the registry adds no element row. **Nothing is added
for them.** Where a run does diffuse elements (the atomic catalog, the
`lower_profile` and `mol_diffusion` regression cases) the registry already adds
those rows to the SAME block, and the design keeps that: one block for every
balance the run alternates, because the reason for the block, the splitting
between a held wind and a moving composition, is identical for an element and
for a carrier. The registry's own ordering rule is kept as it stands (the
elements are registered first because their write-back is a projection of the
whole species vector, `set_transported_species_rows` says so, and a carrier
slot written before it would be rescaled).

---

## 3. The residual of the block

### 3.1 The rows

| row | how it is formed | unit conversion into the block | scale the certification judges it on |
|---|---|---|---|
| mass, momentum, energy of cell j | `assemble_residual` (`steady_residual.f90` line 217) on the conserved state with the ghosts filled by `Apply_BC` and the composition of that state | none (code time) | `residual_row_scale(k, j, u)`; tolerances `cert_tol_mass_at(j,u)` with floor `cert_tol_mass` = 3.0e-12, `cert_tol_momentum` = 1.0e-8, `cert_tol_energy` = 1.0e-6 (`certification.f90` lines 193-195, 800) |
| carrier c of cell j | `carrier_steady_residual`, which calls `carrier_state`, `carrier_geometry`, `carrier_advective_state`, `carrier_diffusivities`, `carrier_face_coefficients`, `carrier_photolysis` and `carrier_residual` with `dt_big = 1.0d30`, "which kills the time term exactly" | `tscale_code = R0/(v0*n0)` | the sum of the row's own physical terms with the two faces counted separately, `abs(F_dif,in) + abs(F_dif,out) + abs(F_adv,in) + abs(F_adv,out) + abs(net source)`; tolerance `cert_tol_carrier_at(r)` = 1.0e-5 at r >= 1.20 R_p |
| element a of cell j | `element_transport_residual` (`binary_element_diffusion.f90` line 2773) on the face mass flux of THIS state, `face_mass_flux_of_state(rho, Frho_elem)` | `elem_he_to_code()`, `elem_tr_to_code()` | `escale_he`, `escale_tr`; tolerance `cert_tol_element_at(r)` = 1.0e-5 at r >= 1.20 R_p |
| the base cell of an element row | the Dirichlet reservoir statement `Y(base) - srow_base_value(i)` where `species_row_base_equation(i) .eq. base_row_reservoir_condition` | none | the same |
| an element with no reservoir | `Y(is)`, the row that says the element stays absent, so that the column is not empty and the band factorization not singular | none | not gated |

The scales are the certification's own, which is the requirement of section 3
step 1 (c) of the plan: each half normalized on its own physical scale, the
scale independent of the iteration length. Nothing in this design invents a
scale.

### 3.2 The coupling terms, and why this is a block

**How a hydrodynamic row depends on the composition.**

1. **Through the pressure and the temperature at the conserved state.** The
   contract is `pressure_and_temperature_at_fixed_conserved_state`
   (`diffusive_photochemistry.f90` line 6162): the pressure of a cell is
   `pressure_from_energy_density(j, rho, e_th)`, the inverse of the caloric
   equation of state, and the temperature then follows from `comp_T_from_p` at
   the particle count of the composition. Review row R25 is the condition that
   makes this a real coupling and not a formality: `p = (gamma-1) e_th` is
   independent of the particle count only at FIXED gamma, and the molecular
   caloric equation of state moves the heat capacity with the composition, so
   at a fixed conserved state BOTH p and T move when the composition moves.
   The measurement of 1.3 is that dependence evaluated: 5.0e-02 in p and
   1.9e-02 in T at the base for a composition whose base H2 partition fell from
   0.934 to 0.784.
2. **Through the mean particle mass and the particle count.** These set the
   temperature at a given pressure and the scale height of the layer; 1.3
   measures 6.4e-02 and 6.9e-02.
3. **Through the heating and the cooling.** `eval_residual` runs
   `excited_H_update` and `ioniz_eq` at every evaluation and hands the heat and
   cool it used into `assemble_residual`, so the energy row's source is the one
   of the composition the evaluation itself produced. That is also why the
   banded model must be the FULL-residual one, `build_banded_jac_full` (line
   4615): its own comment says the frozen version "omits the dominant local
   d(heat-cool)/dE coupling, so Newton has no descent direction even near the
   solution".
4. **Through the opacity.** The optical depth of a cell is an integral over the
   cells below it, so a composition change at the base moves the radiation
   field, hence heat and cool, everywhere above it. This is the non-local term
   `jv_product` exists for.

**How a carrier row depends on the hydrodynamic state.** `carrier_steady_residual`
re-derives everything it needs from `(rho, v, f_sp)` at each evaluation:

1. **The advective face flux**, through `carrier_advective_state`, which forms
   the face mass flux `Frho` from `rho` and `v`, and through the reconstruction
   of the carrier mass fraction `Y = m_c f_c / msum` in `species_face_fraction`.
2. **The temperature in every rate**, through `carrier_state`, which returns
   `TK` and `mbar`; the H2 chemistry of `carrier_source` (line 4119) is the
   Koskinen et al. (2022) Table-1 network of `mol_rates` used by
   `System_HeH_mol.f90`, whose three-body and charge-transfer channels are
   strongly temperature dependent, and whose largest loss channels the L7e row
   dump names (He+ + H2 at 97 percent of the loss at cell 306, H+ charge
   transfer at the front).
3. **The density in the three-body terms and in every two-body rate**, through
   `nd = rho*n0` and the frozen background `cbg_*` of the same state.
4. **The diffusion and settling coefficients**, through `carrier_diffusivities`
   and `carrier_face_coefficients(ntot, TK, mbar, gphys, ...)`: the settling
   drift is `G = (m_c - mbar) g/(kT)`, which depends on the mean particle mass
   of the mixture and on T.
5. **The photolysis**, through `carrier_photolysis(rho, TK, f_sp)`, hence on
   the column density of the gas below.

**So the block's residual is one map of one vector.** In the alternation the
first list is evaluated at a composition the hydrodynamic solve holds and the
second at a wind the carrier relaxation holds; in the block both are evaluated
at the same iterate, and the root is a state at which the wind and the
composition are steady states of each other. That is the sentence the outer
loop already carries for the coupled route (`EXHALE_main.f90` line 6770: "With
the carrier row among the unknowns the outer loop is not an alternation at all:
one solve returns a wind and a carrier partition that are steady states of each
other, and the loop below runs once").

### 3.3 The statement rev3 asks to be answered

Section 3 step 3 (iii) of the plan records that "a linearized pressure
derivative at fixed hydro is not the wind response". Agreed, and the block does
not use one. The block does not linearize the wind's response to a composition
change at all: it carries the hydrodynamic rows as equations and lets the
nonlinear solve find the state at which they and the carrier rows vanish
together. The only linearization anywhere is the Jacobian ACTION used to build
a Krylov direction, which is a directional difference of the full residual and
is discarded at the end of each outer iteration; the acceptance is on the
nonlinear residual (section 5).

---

## 4. The linearization

### 4.1 What the action must carry, and what may be lagged

**The criterion.** A term may be lagged (left out of the Jacobian's action)
where its measured size in the action is below the linear tolerance the Krylov
solve is asked for. The model for measuring that is L22 step 2b section 5: form
the action of the assembled operator on a direction concentrated on named
cells, form a central difference of the full residual along the same direction
with the ghosts filled by the rule the trial uses, and report the relative error
row by row and the size of the entries the band cannot hold. That measurement
is what excluded the carrier reconstruction terms (mean relative error
1.542e-03 with them lagged, and the pass identical to every digit with them in)
and what found the one entry the band drops (cell 314, assembled action exactly
zero against the operator's 1.679e-03).

**In the block nothing of the coupling is lagged in the action.**
`jacobian_action_of_direction` (line 10864) calls `jv_product` (line 10514),
which forms `J v = (F(Y + eps v) - F0)/eps` with `F` the FULL residual of
`eval_residual`, that is, with the equilibrium sweep, the heating and cooling,
the radiation field, the carrier rows and the element rows all recomputed at
the perturbed state. So every cross term of section 3.2, including the
non-local radiation one, is in the action to first-order finite-difference
accuracy. This is the second finding of the source reading, and it is the main
technical argument for (iii) over a hand-assembled block: **the coupling
derivatives do not have to be written down at all.** What has to be measured is
that the difference quotient is a usable one, which is increment I1.

**What IS lagged is the preconditioner**, and deliberately:
`build_banded_jac_full` probes the full residual by graph coloring with stride
`ncolor_jac = kl_jac + ku_jac + 1`, so same-color columns whose row supports
overlap through the non-local column density contaminate each other; the
routine's own comment calls that "preconditioner error". A preconditioner may
be wrong; it changes the number of Krylov products and not the root.

### 4.2 Can the matrix-free JFNK absorb the carrier unknowns

It already does. The question the plan puts is really the other one: the
`Coupled carrier solve: True` route existed before the alternation was chosen,
and the alternation was preferred; **what is different now?**

**What the coupled route did, and why it was set aside.** READ:

- The refusal at startup for a molecular configuration with the carriers
  ELIMINATED (`input_read.f90` line 2220) is the finding B5h: with n(H2) not an
  unknown, the fast chemistry of the shielded layer cycles
  H2 to H2+ to H3+ to H2 without changing the number of H2 nuclei, so the local
  balance rows fix only the partition and leave the content where the seed put
  it; 0.77 of any seed perturbation of the layer's H2 content survives every
  pass of the sweep, at perturbations of 1e-6 and of 1e-2 alike, and the base
  cell's energy row amplifies that by about 4e3, so two evaluations of one state
  differ by 3.2e5 times `Resid tol` per unit relative seed change. With
  `Molecular carrier transport: True` the same measurement reads 9.7e-10, below
  `Resid tol` = 1e-8. **That refusal is a guard the design keeps, and the three
  refusing cases satisfy it already** (they set the key).
- N7 and the items of `docs/PLAN_20260909_rev1.md` Stage C found the route
  taking no step at all on the atomic element fixture: an element with no
  reservoir put 1500 unknowns exactly on the lower face of their own box, the
  fraction-to-the-boundary rule then returned zero for any direction with a
  component pointing below zero, the first Krylov cycle managed 0 products of
  40, the initial trust radius came out zero and eleven outer iterations
  reported a zero dogleg step at `||R||` = 1.888 (READ, the comment at
  `set_transported_species_rows`).
- The probe step of the matrix-free difference, `sqrt(eps)(1 + ||Y||)/||v||`, is
  set by the whole vector while the component of `v` at one species unknown is
  set by the banded preconditioner: on the coupled `mol_carrier` reload it put
  the carrier of cell 217 at -2.4e-16 with 26 species unknowns within round-off
  of their bounds, the sample was refused at every halving, and 12 of 20 outer
  iterations ended with no Krylov direction at all (READ, the comment at
  `largest_step_inside_the_species_box`, line 6919).
- On `mol_diffusion`, one blocked component cost the whole direction and 18
  Krylov cycles were truncated for that reason (same comment).
- N26: the upper ghost of a transported column was frozen at the pre-solve
  composition, leaving the outermost element row unconstrained from outward;
  making it follow the iterate moved the atomic element reload from `||R||`
  1.36e-2 to 3.7e-4 (READ, `docs/code_status_20260910.md` line 107).

**What is different now.** Every item in that list has been answered in the
source, and three further things landed after the coupled route was set aside:

1. **The three refusals above are fixed.** Absent elements are skipped by the
   registry; `largest_step_inside_the_species_box` computes the
   fraction-to-the-boundary length instead of searching for it by halving, and
   counts and skips a component with no room instead of refusing the whole
   direction; the upper ghost follows the iterate (N26).
2. **The trust region is the step control wherever the system carries a species
   row** (`steady_newton.f90` line 16085). MEASURED there on the 1 microbar hot
   Uranus reloaded with `He_diffusion` and `Coupled carrier solve: True`
   (nvar = 4, the same width as the three refusing cases), from one marching
   snapshot with everything else equal: the line search took `||R||` 1.986 to
   1.723 and aborted at outer iteration 18 on "no descent direction exists for
   the banded model", with accepted steps running 1.72, 6.44, 1.82, 25.2; the
   trust region took `||R||` 1.986 to 4.345e-09 against a target of 1e-08 in 54
   accepted steps of 66 outer iterations, with the elemental transport row of
   the returned state at 4.633e-11. **A coupled solve at exactly this width has
   already converged on this tree.**
3. **The L4h pseudo-time guard** (READ, `docs/PLAN_20260913_lhs_stationary.md`
   row L4h, `docs/lhs1140b_stationary_L4h_20260915.md`). The diagnosis was that
   the pin at the floor is the step-size rule and not the direction: at the
   floor the shift `I/dtau` stands above every entry of the Jacobian, the Krylov
   cycle reaches its tolerance in one product and the step is `dY = -dtau F`, an
   explicit Euler step at the CFL interval. The rules now default on: the cut of
   an accepted step is clipped at `max(lam, 0.1)`; it stops at 1e-3 of
   `dtau_earned`, the largest pseudo-time at which a WHOLE step was taken; eight
   consecutive damped steps return the solve to its best iterate with `dtau/10`,
   and two such returns with no full step between them stop it there with
   `info = 2`. MEASURED on the element reload with the earned anchor: at
   iteration 27 the guarded run stands at `||R||` 4.035e-02 and merit 9.71e-04
   against 1.543e+00 and 1.11e-02 unguarded, a factor 38 and 11. **The limit
   cycle of 1.4, line-search factor alternating 1.0 and 9.54e-07 at
   dtau 7.3e-05, is exactly the shape rule (4) was written for**, and that run
   predates nothing: it is a run of the alternation's hydrodynamic half, which
   carries the same guard. So the guard alone does not explain the failure of
   1.4 away, and this design does not claim it does; what it claims is that the
   pseudo-time control the block will use is the guarded one and not the one
   the old coupled experiments met.
4. **The well-balanced momentum scale** (P4 of
   `docs/PLAN_20260911_partitioned_solver.md`). Under `Well balanced: True`,
   which all three refusing cases set, `store_row_terms` used to scale the
   momentum row by `max(abs(dF(2)), abs(S(2)), abs(Smom))` with `S(2) = 0` by
   construction, so every scaled momentum residual was exactly 1. The reference
   is now the pressure force the well-balanced option cancels analytically.
   A block whose merit is a 2-norm over row kinds cannot be steered by a row
   family that reads 1 everywhere.
5. **The R4 carrier Jacobian corrections** (L7e section 23.1): the outward
   face's derivative at the last cell with `Frho(N) < 0`, the inward face's
   derivative at the first cell when the base composition is not imposed, and a
   spurious division by `carrier_mass_amu(ic)` that made the coded coefficient a
   factor `m_c` too small (2.0088 measured for H2). Those are corrections to
   the carrier relaxation's OWN matrix; in the block they matter only through
   the preconditioner, since the action is matrix-free.
6. **The thread question is closed for this path** (L15 stage 2): the whole
   105-iteration trajectory of the element reload is bitwise identical at 1, 8
   and 16 threads, and the molecular path reproduces over two 8-thread runs and
   one 1-thread run. So a coupled block may be measured at 8 threads without
   the measurement being a measurement of the schedule; the controls of this
   design are nonetheless single-threaded (section 7).

### 4.3 The preconditioner

**The banded preconditioner as it stands, with the carrier rows in it.**
`set_transported_species_rows` sets `kl_jac = ku_jac = 3*nvar_jac - 1` and
`ncolor_jac = kl_jac + ku_jac + 1`, from the statement that WENO3 reaches two
cells and within a cell every variable couples, so a column touches rows within
`2*nvar + (nvar - 1)`. For the three refusing cases `nvar_jac = 4`, so
`kl_jac = ku_jac = 11` and `ncolor_jac = 23`.

**And that band holds the entry the carrier relaxation's matrix could not.**
The dropped entry of L22 step 2b is `d res_carrier(j+2)/d f_c(j)`, at a flat
distance of `2*nvar_jac = 8` in the block's ordering, inside a band of 11. The
block therefore needs no special treatment for the reconstruction stencil: it
is inside the band by construction, and `build_banded_jac_full` probes it by
finite difference of the operator itself, so nothing has to be differentiated
by hand. The measurement of I1 states it rather than assuming it.

**Block structure.** The ordering is `Y(nvar*(j-1)+k)`, so the matrix is block
banded with a 4 by 4 dense block on the diagonal (the hydro-hydro, the
hydro-carrier, the carrier-hydro and the carrier-carrier couplings of one cell)
and two further cells of blocks each side. The factorization is `dgbtrf` on the
general-band storage `ab(kl+ku+1 + i - j, j)` with `ldab = 2*kl+ku+1`; the
route also carries a normal-equations variant at `kl_jac_normal = 2*kl_jac`.
No new structure is proposed. The one decision is which banded builder the
block uses: **`build_banded_jac_full`, not `build_banded_jac`**, because the
frozen version omits `d(heat-cool)/dE`, which its own comment says leaves
Newton without a descent direction even near the solution; that cost is
`ncolor_jac + 1` = 24 full residual evaluations per preconditioner build at
nvar = 4 (section 8).

### 4.4 The pseudo-time

**One `dtau` for all rows, the PTC ramp with the L4h guard; NOT the carriers'
own `dt_code`.** Reasons:

1. A single `dtau` with the row scaling that already exists is a shift
   `I/dtau` of the SCALED system, so it does not change the root; the
   root-preserving cell-local variant `M_j dU_j/dtau_j + R_j(U) = 0` is the
   subject of step 3 (i) and is being tested there on the carrier operator
   alone. Writing a second, different pseudo-time distribution inside the block
   before (i) reports would make the two experiments inseparable, which is the
   error review row R45 warns against.
2. `dt_code(j) = min(t_diff, t_adv)/t_scale` (the diffusion and advection
   scales of the carrier operator, formed in `relax_photochemical_composition`)
   is the explicit stability interval of the carrier transport at a FIXED wind.
   In the block the wind is not fixed and that interval is not the stability
   interval of anything the block integrates.
3. The measured pathology the block must avoid is the pseudo-time collapsing to
   the explicit interval (L4h: "an explicit march at the CFL interval by then").
   Handing the block an interval formed cell by cell from an explicit stability
   condition walks into it by construction.

If the block's pseudo-time turns out to be held by the carrier rows alone, the
evidence for that is the row-kind ledger the route already prints
(`write_row_kind_shares`, which reports which row kind holds `||R||` and which
holds the merit), and the remedy at that point is a row scaling and not a second
clock. That is an outcome of the measurement, not a design decision taken now.

---

## 5. Globalization: what replaces the movement bound

**Inside the block the composition movement bound is not evaluated at all.** It
is a property of the alternation, where it bounds how far a composition may move
at a held wind; in the block the wind is not held and the quantity it bounds,
the particle count at fixed conserved state, is a function of unknowns the block
is solving. The bound stays, unchanged, in the alternation.

What takes its place is already in the coupled route, and this design adopts it
as it stands:

1. **Positivity and the species box.** `species_box_lo` and `species_box_hi` per
   species unknown; `largest_step_inside_the_species_box` (line 6919) returns
   the largest multiple of a direction that keeps every species unknown inside
   its box, the fraction-to-the-boundary rule; `zero_the_blocked_components_of`
   removes a component with no room at any step length so that one blocked
   component does not cost the whole direction; `fix_active_species_bounds` and
   `freeze_species_unknown_box` declare the active set, and `jv_product` applies
   the SAME declaration to the operator it samples (`hold_the_active_bounds_of`),
   so the difference quotient is a directional derivative of the residual
   composed with the projection onto the box, which is the map the trial itself
   evaluates. A probe may use at most `box_safety` = 0.9 of the room to a face,
   because a step to the face itself puts an unknown where the reaction network
   has no interior to linearize about.
2. **Simplex feasibility, as element constraint rows.** The route forms, per
   cell, an elemental budget from `carrier_headroom(j, ic)` times
   `nuclei_per_particle(isp, ien)` for hydrogen, oxygen and carbon, the
   constraint value `c(ie)` and its derivatives with respect to every unknown of
   the cell through `element_constraint_derivatives`, and it declines to form
   the rows at all when a budget is not finite. That is the statement "the
   carriers of a cell may not claim more nuclei of an element than the cell
   has", which is the simplex feasibility the reviews ask for, carried as
   constraint rows with gradients rather than as a scalar excursion budget.
3. **The trust region on the joint step.** `use_tr` is the default wherever the
   system carries a species row, with the dogleg's Newton leg built by
   `pgmres_with_restarts` at a stated relative tolerance `tr_gm_rtol` (a loose
   leg makes the model predict an increase at full step length, section 159.3),
   and the model of a PROJECTED step recomputed with its own product. The
   measurement of 4.2 item 2 is the case for it over the line search.
4. **The admissibility screen of the nonlinear trial.** `eval_residual` returns
   `admissible = .false.` when the equilibrium sweep met a cell whose reaction
   residual is not finite, when any cell's composition was accepted without
   being a certified root of the network (acceptance class 4 or 6), or when a
   row of the assembled residual is not finite; `trial_state_is_admissible` and
   `probe_direction_is_usable` (`certification.f90`) are the screens the step
   control and the probe read. A trial that is inadmissible is rejected, not
   damped.
5. **The judged merit of L4e**, if it applies: the route already distinguishes
   the convergence measure `||R||` (a maximum over cells of each hydrodynamic
   row against its own largest term) from the merit (a 2-norm of every row on
   the Newton's own row scales), and prints which kind holds each. The design
   does not change either; it requires the ledger to be ON for every run of
   section 7, because a block whose merit is held by one row kind while `||R||`
   is held by another is the first thing to know about it.

**The acceptance is unchanged.** A state is accepted only when the
certification of the FULL set of active equations passes on the REFRESHED
state, from a residual assembled at that state: the three hydrodynamic rows,
the transported carrier balances, the elemental transport balances, the level
populations, the eliminated-species closure and the conservation records, each
against its own tolerance (`certification.f90`; the contract is stated at
`steady_wind_with_element_diffusion`, `EXHALE_main.f90` line 6495). In
particular the carrier rows are judged by `carrier_steady_residual` at
`cert_tol_carrier_wind` = 1.0e-5 over `r >= 1.20 R_p` exactly as now, the
element and charge identities are the ones the certification already records,
and **no tolerance and no gate is changed by this item.** A hydrodynamic
`info = 0` with an alternated species row refusing is not an accepted state,
and neither is a coupled `info = 0` with a refusing row: the flag is a statement
about the rows in the Newton registry and the certification is the acceptance.

Since the block's pseudo-time is a shift of the scaled system and its
intermediate iterates are not physical-time updates, **no inventory constraint
is imposed on intermediate iterates** (review row R33, rev3 acceptance for (i)).
The inventory equation `E_a` of the plan's acceptance belongs to (ii),
conservative physical subcycling, and is not part of this design.

---

## 6. The contract with the rest of the solver

### 6.1 When the block is entered

**Recommendation: the alternation runs first and the block is entered when the
composition's progress is held by the bound, not always.** The reason is the
alternation's one measured virtue: it is a globalization. The state the block
would start from at pass 1 of a campaign case is a mapped seed, and the coupled
experiments of 4.2 are the record of what a coupled solve does from a state far
from its root. The alternation's passes at the bound move the front and leave
the wind stationary at each step, which is a continuation; what they cannot do
is finish (section 1.5).

The entry condition, stated on quantities the outer loop already forms:

- the carrier relaxation of the pass returned `carrier_relax_movement_bound`
  (the ending `carrier_ending_block_closed` is its complement), AND
- the pass did not count as progress under the joint rule of L14 section 3.1 (a
  pass counts as progress when the joint distance from certification fell, or
  when the state's distance is carried by an entry other than a hydrodynamic
  row and the composition's own distance to the fixed point of its operator
  fell), AND
- that has now happened on `outer_no_fall_max` = 3 consecutive passes, which is
  today the condition on which the loop REFUSES the state ("the joint distance
  of the state has not fallen in 3 consecutive passes; the alternation is not
  approaching a joint fixed point", `EXHALE_main.f90` line 7612).

So the block is entered exactly where the alternation gives up, and the change
to the outer loop is that this ending becomes a handover rather than a refusal.
That keeps every run in which the alternation certifies bit for bit what it is
today, which is the standing rule for a new algorithm.

The key: `Coupled carrier solve:` already exists and already selects the block
for the WHOLE solve. The handover needs a third state, so the key becomes
three-valued, `False` (the alternation only, today's default and the default
after this item), `True` (the block only, today's meaning, unchanged) and
`On stall` (the alternation, then the block from the state it stalled at). The
new value is opt-in; `docs/input_schema.md` and the setup report echo it.

### 6.2 What the outer loop does with the block's result

The block returns a state and an `info`. The outer loop then does what it does
with any refreshed state: it refreshes the primitives (`U_to_W`), the species
densities and the temperature of the composition beside them, and certifies. A
block that returns `info = 0` and certifies ends the solve with the state
accepted; a block that returns `info = 0` and does NOT certify is a refusal
named as such, exactly as the alternation's is; a block that returns nonzero
hands back its best iterate under the L4h rules and the loop reports which
entries refuse. There is no path on which the block's own flag accepts a state.

### 6.3 How the certification reads it

Unchanged, and it must be unchanged, because a block that changed the
acceptance would be a knob and not an algorithm (the L22 row's own condition:
"how any of the three is prevented from becoming a knob that certifies a state
by loosening what holds it"). One point deserves a statement because it is a
real difference: after the block the composition is the one the block's own last
residual evaluation produced, and the certification evaluates its rows on the
state it is handed. Today, after a carrier relaxation, the same is true. So
nothing about the certification's input changes; what changes is which operator
produced it, and the certification names the operator it certifies against
(`docs/restart_contract_design_20260909.md` section 7).

### 6.4 The restart intent and the written state

Unchanged (`restart_contract_design_20260909.md` section 7 and 7.1). The work
state holds the conserved density, momentum and total energy, and p and T follow
the composition at that conserved energy through the caloric equation of state,
which is the contract `pressure_and_temperature_at_fixed_conserved_state`
states; the block's unknowns are exactly those conserved variables plus the
species, so a state it hands back is written and reloaded by the same rules as
one the alternation hands back. `assert_written_state_is_the_accepted_one`
applies unchanged. The one thing the block adds to the written metadata is which
route produced the state, which belongs in `cert_reason` under the L23 metadata
contract and is that item's to carry, not this one's.

---

## 7. The tests, in order

**T1, the coupling terms of the action against a finite difference of the full
residual, on a frozen state.** The model is L22 step 2b section 5. A new suite
`src/tests/coupled_block_jacobian/` reloads a frozen molecular state with the
block registered (`set_transported_species_rows(.true.)`), forms
`jacobian_action_of_direction` on directions that are (a) purely hydrodynamic,
(b) purely carrier and (c) mixed, concentrated on named cells, and compares each
against a central difference of `eval_residual` along the same direction with
the ghosts filled by the rule the trial uses. Rows:

- `the_action_of_a_carrier_direction_moves_the_energy_row` (the cross term of
  section 3.2 item 1, nonzero and of a stated size);
- `the_action_of_a_hydrodynamic_direction_moves_the_carrier_row` (the cross term
  of 3.2 items 1 to 5);
- `the_assembled_action_reproduces_the_central_difference` to a stated relative
  error over the direction's support, the number the criterion of 4.1 is read
  against;
- `the_band_holds_the_two_cell_reconstruction_entry` (the claim of 4.3: the
  entry at a flat distance of `2*nvar_jac` is inside `kl_jac`);
- `the_outer_ghost_rule_is_the_same_in_both_evaluations` (the assumption of the
  preamble; it fails if L22 step 2c leaves two rules).

RED before the increment is not available for the new module entities the suite
names, exactly as L22 step 2b recorded for its own suite; the rows that can be
stated against the entry text are stated against it and the ones that cannot are
marked as such in the suite's README.

**T2, a manufactured column the block must reproduce.** The certified
`molecular_scalar_gj1132_kzz1e9/HeH2.13` state, reloaded, solved by the block,
and required to return a state certifying on the same entries to the same
tolerances, with the profiles agreeing with the alternation's certified solution
to the certification tolerance. Review row R37 is the condition under which this
is read: the old certified state is not the oracle of a changed operator, so the
test is that the block's state CERTIFIES, and the profile comparison is reported
beside it and not used as the gate. A byte comparison is not expected and is not
asked for.

**T3, the three refusing cases**, in the order the measurements put them:
`molecular_scalar_gj1132_wellmixed/HeH0.083` (stalled, the refusing cell at the
domain edge), `molecular_scalar_gj1132_kzz1e9/HeH0.083` (the single-cell spike
at cell 500, and the case whose pass-2 hydrodynamic solve is the weaker half)
and `molecular_scalar_gj1132_wellmixed/HeH0.55` (still travelling; its front
moved about one cell per pass). Each from the state the campaign left, each
single-threaded, each with the row-kind ledger on. The outcomes to report per
case: the certification entry by entry, the worst gated carrier row and its
cell, `||R||` and the merit with the row kind holding each, the trust-region
accepted and rejected counts, the pseudo-time trajectory against the L4h rules,
the Krylov cycles that exhausted their subspace (the L24 failure mode), the
admissibility refusals and the blocked species components, and the wall time.

`molecular_scalar_gj1132_wellmixed/HeH2.13` is NOT in this list: it has no
seed, which is a different statement (rev3 section 3, review row R8).

**T4, the goldens that move.** None is expected to move, and the reasoning is a
code-path one that must be measured and not assumed:

- `mol_carrier`, `mol_base_handoff`, `mol_diffusion`, `mol_metals`,
  `mol_lyman_werner`, `mol_ir_bands`, `mol_sec_ion` are MARCHING cases; L22
  step 2 measured that `mol_carrier` never prints an `outer pass` line, so the
  outer loop and the block are not entered at all.
- `lower_profile` and the `*_newton` fixtures reach the stationary route. With
  the new key value absent they take the alternation, unchanged.
- The fixtures that would change are the ones a run of this item POINTS at the
  block: `backup/regression/carrier_elem_newton`, `carrier_model_a_newton`,
  `atomic_elem_newton`. Those are fixtures, not goldens of the physics matrix,
  and their movement is reported as a measurement.

The rule stands: goldens are refreshed once, at the end of the change series,
with the movement reported, and never to make a new algorithm pass.

**Suites to run** (the rule of `docs/worker_rules.md`: every suite whose driver
links a module the increment changed, not only the one the brief names):
`krylov_and_dogleg`, `attempted_step`, `steady_completion_flag`,
`steady_species_rows`, `certification`, `grid_and_gates` (which carries
`coupled_carrier_h2`, the startup refusal rows of the key), `carrier_retry`,
`carrier_returned_state_acceptance`, `carrier_reference_scales`,
`carrier_constraint_attribution`, `carrier_boundary_jacobian`,
`carrier_outer_boundary`, `element_operator`, `molecular_seed`,
`steady_selfconsistent_residual`.

---

## 8. Cost, increments, files and risks

### 8.1 Cost

At `nvar_jac` = 4 and N = 500 (the three refusing cases): 2000 unknowns,
`kl_jac` = `ku_jac` = 11, `ncolor_jac` = 23.

| item | count | what one costs |
|---|---|---|
| unknowns | 2000 (4 per cell) | |
| full residual evaluation | 1 | one `Apply_BC`, one `excited_H_update`, one `ioniz_eq` sweep over the column, one `carrier_steady_residual` (which re-forms the diffusivities, the face coefficients and the photolysis), one `assemble_residual` |
| banded preconditioner build | `ncolor_jac` + 1 = 24 full residual evaluations | the dominant cost of an outer iteration |
| band factorization `dgbtrf` | 1 | about `neq * kl * ku` = 2.4e5 multiply-adds, negligible beside one residual evaluation |
| GMRES product | 1 full residual evaluation each | at `EXHALE_GM_M` = 40 a full cycle is 40 |
| outer iteration | 24 + (products actually taken) + trial evaluations | |

So an outer iteration costs between about 30 and 70 full residual evaluations
at nvar = 4, against 8 + products at nvar = 3 for a hydrodynamic-only solve
(`ncolor` = 17 there, so 18 + products). The block is therefore roughly 1.3
times a three-unknown iteration in preconditioner cost and identical per
product; the width of the vector is 4/3 of it. Against the alternation the
comparison is not a ratio of iterations but of outcomes: the campaign runs of
the three cases spent 24 and 40 outer passes each and did not certify.

The reference point for how many outer iterations to expect is the one coupled
solve measured at this width: 54 accepted steps in 66 outer iterations from
`||R||` 1.986 to 4.345e-09 (4.2 item 2). A budget of 200 outer iterations per
solve is the starting cap, with the L4h best-iterate return as the stopping
rule.

### 8.2 Increments, each with its own test

| # | what lands | files touched | its test | entry condition |
|---|---|---|---|---|
| I1 | the coupling measurement only: a diagnostic key that writes the block's Jacobian action against a central difference of the full residual on a frozen state, and the suite that reads it | `src/tests/coupled_block_jacobian/` (new), `src/modules/time_step/steady_newton.f90` (the diagnostic key and its writer, no change to any solve path) | T1; byte identity of a marching case with the key unset | L22 step 2c landed, so the outer ghost rule is one |
| I2 | the three-valued `Coupled carrier solve:` key and the handover from the alternation, with its log and its refusal wording | `src/modules/files_IO/input_read.f90`, `src/modules/files_IO/write_setup_report.f90`, `src/EXHALE_main.f90`, `docs/input_schema.md` | `grid_and_gates/coupled_carrier_h2` extended with the new value; a run with the key absent reproduces the alternation bit for bit | I1 |
| I3 | the block on a state that already has an answer | none (a measurement on I2's binary) | T2 on `kzz1e9/HeH2.13` | I2 |
| I4 | the block on the three refusing cases | none | T3 | I3 certifies |
| I5 | the record, the fixture movement, the golden decision | `docs/lhs1140b_stationary_L22_20260916.md` (a step 3 section), `docs/PLAN_20260913_lhs_stationary.md` (row L22), `docs/audit_20260905/L22_step3_20260916/` | the suites of section 7 | I4 |

Nothing in I1 to I5 changes a default; the block is reached only by the new key
value until a measurement says otherwise, and that decision is the user's.

### 8.3 Risks, each with the measurement that meets it BEFORE it is met

| risk | where it was met before | the measurement that meets it first |
|---|---|---|
| species unknowns on their bounds, no Krylov direction | N7 (1500 unknowns on a bound, 0 products of 40, zero radius, eleven zero dogleg steps at `||R||` 1.888) | I1 reports the blocked-component count and the fraction-to-the-boundary length on the frozen state, before any solve is attempted |
| the probe step set by the whole vector leaves a species unknown outside its box | the `mol_carrier` reload (carrier of cell 217 at -2.4e-16, 12 of 20 outer iterations with no direction) | the same I1 report, at the probe length the route would actually use |
| the linear solve exhausting its subspace on one species row | L24 (62 of 91 Krylov cycles at 0.91-0.999 of the right-hand side on the helium element row of cell 335) | T3 reports Krylov cycles that reached the cap and the row they were held by, per solve |
| the pseudo-time collapsing to the explicit interval | L4h; and the limit cycle of 1.4 | T3 reports the pseudo-time trajectory against `dtau_earned` and counts the firings of L4h rules (1), (2) and (4) |
| the carrier row's local chemistry having no root at a trial state, so the trial is inadmissible in a whole layer | the acceptance classes of `eval_residual`; the `n_trial_no_chem_root` counter of `solve_refusal_statistics` | T3 reports the refusal statistics, which the route already accumulates |
| the base cell's mass row, which is the lower-atmosphere handoff and the wind network disagreeing across one face | L7e section 21, open as L7f; and 1.4's worst cell at r = 1.038 | T3 reports the mass row cell by cell in the inner layer; if the block converges everywhere except there, the item is handed to L7f and says so |
| the outer ghost rule differing between the block's evaluation and the certification's | L22 step 2b section 1 | the last row of T1 |
| a golden or a fixture moving through a path the key does not gate | the `mol_carrier` comparison of L22 step 2b section 7, which had no control because the tree moved under it | I2 builds a control binary from the same live tree with its own files at their entry text, as `docs/worker_rules.md` requires, and the comparison is made against that and not against the tree's `EXHALE.x` |

One risk has no measurement that precedes it, and it is named as such: the
block may converge to a state whose H2 front sits outside the domain. P21 of
`docs/PLAN_20260911_partitioned_solver.md` measured that at a FIXED wind, with
the movement bound disabled, the carrier operator reaches a genuine fixed point
whose column is fully molecular, x2 = 0.98 at 1.1 R_p and 0.33 at 3.0 R_p, with
the front outside the domain; the front's position is then set by the
hydrodynamic response and not by the carrier operator. The block is exactly the
algorithm that lets the hydrodynamic response set it, so the answer it gives is
the answer to that question, and it is a physical result to be compared against
Koskinen et al. (2022) Model A and not a solver failure. It is stated here so
that a front leaving the domain is not read as a defect when it appears.

---

## 9. The decisions for the user

Everything not listed here is decided in this document.

| # | decision | recommendation |
|---|---|---|
| D1 | The coupled block is the EXISTING `Coupled carrier solve: True` route of `solve_steady_jfnk`, extended with an entry condition, rather than a new solver written for L22 | **Adopt.** The route already carries the hydrodynamic, carrier and element rows in one residual, a matrix-free action of the full residual that needs no coupling derivative written by hand, a banded preconditioner whose band holds the entry the carrier relaxation's matrix drops, a trust region measured to converge at this width, and the positivity and simplex feasibility the reviews require. Writing a second block would duplicate all of it. |
| D2 | The block's unknowns are the whole column, every registered balance; the base-band solve with the wind above held is rejected | **Adopt the whole column.** Three measured reasons in section 2.2. Keep the band as a named fallback if the whole-column block is too slow, to be opened only on a measured cost. |
| D3 | Entry: the alternation first, the block on stall, rather than the block always | **Alternation first.** The alternation is the continuation that gets the state near enough for a coupled solve; the block is what finishes. Every run in which the alternation certifies is then unchanged. The alternative, block always, is available today by setting the key and needs no work from this item. |
| D4 | The pseudo-time inside the block is one `dtau` with the L4h guard, not the carriers' `dt_code` | **Adopt one `dtau`.** An interval formed cell by cell from an explicit stability condition is the pathology L4h diagnosed, and a second pseudo-time distribution would make this item and step 3 (i) inseparable. |
| D5 | The composition movement bound is not evaluated inside the block; it stays in the alternation unchanged | **Adopt.** Its replacement inside the block is the species box with the fraction-to-the-boundary rule, the element constraint rows, the trust region and the admissibility screen, all of which exist. |
| D6 | The acceptance, the tolerances and the certification gate are unchanged by this item | **Adopt.** Any other answer makes the block a knob, which is what the L22 row forbids. |
| D7 | Whether a failure confined to the base cell's mass row is handed to L7f rather than chased here | **Hand it to L7f.** That row is the lower-atmosphere handoff disagreeing with the wind network across one face, it is open since L7e section 21, and it is not a property of the splitting this item addresses. |
| D8 | Whether the three refusing cases are run at 8 threads or 1 | **1 thread for every control and comparison, 8 permitted for a first exploratory solve.** L15 stage 2 closed the reproducibility question for this path, but the series' own rule is that controls are single-threaded. |
| D9 | Whether to open increments I1 and I2 now, with I3 to I5 gated on I1's measurement | **Open I1 and I2.** They change no default and no number; I1 is the measurement that decides whether the block is worth I3. |
