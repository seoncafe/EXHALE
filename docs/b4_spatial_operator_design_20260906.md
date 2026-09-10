# B4: the conservative spatial operator, boundary count and hydrostatic residual

**Status: design for the advisor's approval. Nothing here is implemented.**
Step B4 of `docs/PLAN_20260906_rev2.md`, Phase 4 of
`docs/development_plan_20260905_rev3.md` section 4.5 (D1, D2, D3). It takes
its target from `docs/b1_target_system_20260906.md` sections 1, 2, 4 and 9.2,
its inventory of what exists from `docs/d0_governing_system_20260906.md`
sections 2.2, 2.6, 4.4, 4.6 and `docs/b1a_active_equation_inventory_20260906.md`
sections 2 and 3.

**Provenance.** Every statement about the present code is READ from the tree at
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` on 2026-09-06, with
file and routine. Nothing was MEASURED for this document: no build, no run, no
test. Numbers quoted from other documents are marked as READ from that
document, with the document named. Line numbers drift while parallel steps of
the same plan edit the tree; the routine name is the citation.

**Owners.** Every increment carries `Owner: (unassigned)`; the advisor assigns
them with the implementation briefs.

---

## 1. The operators as they stand today

### 1.1 The table

`fluxes` means a divergence of quantities defined at cell faces `r_edg(j)`;
`advective` means a difference of cell values multiplied by a cell-centered
velocity, which conserves nothing by construction. "Same faces as mass" asks
whether the operator's face flux is built from the hydrodynamic face mass flux
`face_flux(1,j)` that `RK_rhs` stores.

| row / class | discretization | reconstruction | time level | boundary treatment | same faces as mass |
|---|---|---|---|---|---|
| mass `rho` | face fluxes, `dF(1,j) = (A_p F_p(1) - A_m F_m(1))/dV`, `RK_rhs.f90` cell loop | PLM (MC/minmod, `theta = 2`, `PLM_rec.f90`) or ESWENO3 (`Reconstruction.f90`) of the primitive vector `W = (rho, v, p)`, one `gamma_eff` per side | SSP-RK3, three stages, `Apply_BC` after each | base: the characteristic face state at `r_edg(0)` is put into `WL_out(:,0)` by `Rec_BC` and the interior reconstruction is left to the Riemann solver; outer: zero-gradient copy, linear extrapolation under WENO3 | definitional |
| momentum `rho v` | face fluxes plus a cell source: `S(2) = -0.5(rho_L + rho_R)(Gphi_i(j) - Gphi_i(j-1))/dr` and, under PLM, `+(A_p - A_m) p_C/dV` (`source`, `Source.f90`); under WENO3 `dF(2)` gains `(p_R - p_L)/dr` (`RK_rhs.f90`) | as mass | as mass | as mass | n/a |
| energy `E = e_kin + u_th` | face fluxes `v(E+p)` plus the gravity work `A_p F_p(1)(Gphi_i(j) - Gphi_c(j)) - A_m F_m(1)(Gphi_i(j-1) - Gphi_c(j))` | as mass | as mass; sources at a later split stage | as mass | the gravity work term uses `F(1)`; the enthalpy flux is bulk, carries no species |
| equilibrium species (every stage of H, He, the metals, the molecular ions) | **none**: not transported at all. `f_sp` is re-solved from scratch in every cell by `ioniz_eq` each step, and `rho` was then overwritten by `calc_rho` of that composition (D0 C2) until item B3c on 2026-09-06 made the density `intent(in)` and `calc_rho` a check | none | after the hydro, before the energy step (`EXHALE_main` marching order, D0 section 6.1) | the two base ghosts re-solve their ionization locally (D0 C18) | no flux exists |
| the `_adv` post-process species (`System_implicit_adv_H/HeH/HeH_TR`) | cell-centered upwind marching of an ODE in the mixing ratio: `x_j` from `x_{j-1}` with `c1 = dr/v(j-1)` (`post_process_adv.f90` around the `adv_cell%c1 = As` assignments), point values, backward Euler in space | none | post-process only, off the marching trajectory | the innermost cell keeps the input fractions; cells outside the validity gate keep the equilibrium value | no: it is a non-conservative fraction ODE and never sees a mass flux |
| transported carriers (H2, OH, H2O, CO, and H+ under `Ionization transport`) | **split**: the diffusive/drift part is a face flux `Jf(j)` on `r_edg`, the advective part is a cell-centered upwind difference `cadv (f_j - f_{j-1})` with `cadv = n_tot(j) v(j)/(r_j - r_{j-1})` (`carrier_residual`, `carrier_geometry` sets `ntv = ntot*v*v0`) | donor cell, plus a deferred van Leer harmonic correction `adv_corr` frozen inside the Newton (`carrier_slope`) | backward Euler, one implicit step per marching step, split after the element diffusion | base: direction from the base face mass flux, never `v(1)`; H2 Dirichlet inflow applied **as an advective term only**; face 0 carries no diffusive flux; the proton gets no Dirichlet value | **no** when this was written: the advective term was not a flux at all, and `ntv` is a cell product `n_tot v`, unrelated to `face_flux(1,j)`. **The cell-velocity advection was removed on 2026-09-07 (item B4-1c):** in the marching loop the carriers are advected as `F_rho Y_c/m_c` inside the Runge-Kutta stages and the carrier rows carry the diffusive half alone. **Item B5b removed the cell-centered form entirely:** the stationary rows and the fixed-wind relaxation take the same face-flux divergence, and `ntv`, `adv_corr` and `carrier_slope` no longer exist |
| element diffusion, He/H binary | **split**: diffusive face flux `J` on `r_edg` (`element_face_flux`), advection in **advective form** `rho v dX/dr` one-sided on the **cell** velocity (`composition_residual`, header "solved here in the ADVECTIVE form") | one-sided upwind on the cell velocity, deliberately not face-averaged | backward Euler, one step per marching step, split before the carriers | inner Dirichlet reservoir `X_base` written into cells `1-Ng .. 1`; outer zero gradient | **no** when this was written, and the code said why: the routine receives `rho, v, T, dt` and not the hydro face mass fluxes (D0 C28). **Changed on 2026-09-06 (item B4-1):** the marching route takes the advection on the hydro's face mass fluxes inside the Runge-Kutta stages and `solve_mass_fraction` is called with `advect` false; the advective form survives in the fixed-wind relaxation |
| element diffusion, trace metals | same split, on the mixing ratio `f_X = n_X/n_H` (`solve_trace_element_in_hydrogen`), against a **fixed** hydrogen background | one-sided upwind on the cell velocity | with the binary step | inner Dirichlet `fbase`, outer zero gradient | **no**, and additionally **no counter-flux**: `sum_s rho_s w_s /= 0` by the amount the metals carry (D0 C29) |
| chemical (formation/excitation) energy `u_form` | **none**: not a row, not a flux, not a source of the transport step | none | n/a | n/a | absent |
| species enthalpy in the energy flux | **none**: the face energy flux is `v(E+p)` with `E` thermal plus kinetic only | none | n/a | n/a | absent |
| low-Mach dissipation | added to the numerical flux before the divergence (`RK_rhs` adds `Ddis`), so the marching right-hand side and the steady residual carry the same term; the mass component is identically zero | fourth difference of `v` | inside the RK stage | none of its own | it is a face flux, and it is conservative |
| Shapiro filter | cell-value 1-2-1 low pass `u_j += 0.25 eps (u_{j-1} - 2u_j + u_{j+1})` on the conservative variables (`Apply_BC.f90` `shapiro_filter`), applied after the energy step | none | after the whole split step | uses the ghosts `j = 0` and `j = N+1` | **no**: it is not in flux form and carries no cell volume, so on a stretched grid it does not conserve `sum_j dV_j u_j` even in the interior. Off by default |
| radiation column | rectangle rule on cell centers, `N(j) = N(j+1) + n(j) dr_j(j) R0` (`calc_column_dens`, `calc_column_dens_one`, `calc_column_dens_metals`) | none | before the sweep | outermost cell seeded with its own `n dr` | it is an inner-face column used as a cell mean; rev 3 section 10.4 assigns it to Phase 4 (D1) |

Grid note, READ: `define_grid.f90` (`define_grid`, both cell-width
assignments) now reads
`dr_j(2-Ng:N+Ng) = r_edg(2-Ng:N+Ng) - r_edg(1-Ng:N+Ng-1)`, so the one-cell
shift that `docs/audit_20260905/source_audit/group_A_report.md` reported as a
P0 (the width of cell `j+1` stored in `dr_j(j)`, 0.7 to 1.5 percent over the
stretched region, READ from that report) is corrected in the present tree.
B4 inherits the corrected grid and must not re-derive that defect.

### 1.2 Consistency verdict

**No. The species mass fluxes are not consistent with the mass flux today, and
for most of the species there is no species mass flux at all.** Stated by
class:

1. **Equilibrium species**, which is every stage of every element in every
   configuration, are not transported. Their spatial operator is the identity:
   `ioniz_eq` re-solves the local equilibrium in each cell, and then
   `calc_rho` of the accepted composition **overwrites** `rho`, which
   `EXHALE_main` pushes back into the conserved mass row. So the composition
   does not ride on the mass flux; the mass row is instead modified by the
   composition (D0 C2, B1 T2.1). `sum_s m_s F_s = F_rho` cannot even be
   evaluated: only one side of it exists.
2. **Carriers.** Only the diffusive and drift part is a face flux. The
   advective part is `n_tot v df/dr` on cell values, so the sum over carriers
   of the advected mass is not `F_rho` at any face, and a face-by-face
   comparison is not defined because the operator forms no advective face
   flux.
3. **Element diffusion**, both the H/He binary and the trace metals, is
   advective by construction and says so at the code site: helium mass is
   conserved only to the order of the hydro truncation in a transient, exactly
   at steady state. The diffusive half is conservative on the same face
   positions `r_edg(0:N)` as the hydro, but it is not built from the hydro's
   face states.
4. **Chemical energy** is transported by nothing, and the energy face flux
   carries no species. The bulk enthalpy flux `v(E+p)` is the only material
   energy flux in the code.

The one place where a species operator used the hydro face was the
**direction** of the carrier base inflow, taken from the base face mass flux
rather than from `v(1)`. That was a sign, not a flux. Since item B5b every
carrier row that carries advection reads the face mass flux itself, and the
direction at the base is the sign of that flux.

Two further facts a brief must carry, because they are the reason the present
inconsistency is not visible in the accepted states:

- the element census (`element_census.f90`) checks the element **ratio**
  `n_El/n_H` always and the **absolute** nucleus density only when the caller
  declares that its operator held `rho` fixed, and the module is off unless
  `EXHALE_ELEMENT_ASSERT` is set. So a common factor moving every density
  passes;
- `project_elements` meets both element totals exactly after the diffusion
  step by depositing each shortfall into the neutral ground stage, and
  `limit_to_element_budget` rescales the carriers back onto the element
  simplex after the transport step. Both restore the element constraint after
  a non-conservative update, so what a run reports is the constraint, not the
  conservation.

---

## 2. The target operator

### 2.1 One set of faces

**T-B4.1.** There is exactly one set of faces, `r_edg(1-Ng:N+Ng)`, and exactly
one mass flux per face, `F_rho(j)`, the one the Riemann solver returns and
`RK_rhs` stores in `face_flux(1,j)`. Every material transport in the code is a
divergence of a flux defined at those faces, with the same spherical volume
`dV = (r_edg(j)^3 - r_edg(j-1)^3)/3` and the same areas `r_edg^2`.

**T-B4.2 (the consistency identity).** For every face `j` and every accepted
state,

```text
sum_s m_s F_s(j) = F_rho(j)      exactly, by construction,
F_s(j) = F_rho(j) Y_s^face(j) / m_s  +  A(j) rho^face(j) w_s^face(j) / m_s ,
sum_s Y_s^face(j) = 1 ,    sum_s rho_s^face(j) w_s^face(j) = 0 .
```

`Y_s` is the species mass fraction with the electron mass carried on its ion
(B1 decision 2), so `sum_s Y_s = 1` is the same statement as the species mass
sum of B1 T2.2 and the free electrons contribute no mass flux. The identity is
by construction and not by tolerance: the advective part is one scalar
`F_rho(j)` multiplied by a vector that sums to one, and the diffusive part is a
vector that sums to zero.

**T-B4.3 (how the face composition is built).** `Y_s^face(j)` is the
reconstruction of `Y_s` with the **same** limiter as the hydrodynamic
primitive variables, evaluated on the side the Riemann solver's contact wave
selects, and then normalized to sum to one. Three properties are required and
are the increment-1 tests: a uniform `Y_s` is preserved exactly for any
velocity field; `Y_s` stays in `[0, 1]`; and the normalization moves no
species by more than the reconstruction's own truncation error, which is
measured and reported rather than assumed.

Deciding the upwind side from the contact speed rather than from a cell
velocity is what removes the two defects the present operators each worked
around in their own way: the cell whose two face velocities straddle zero
(which froze the composition of the first free cell of a breathing base, the
reason the element operator refuses face-averaged upwinding) and the cell that
is outflowing at both faces (which evacuated the metals, the reason the trace
arm moved to the mixing ratio). Neither arises when the species ride on the
same face mass flux the density rides on, because a cell can then only lose
what the mass row says it loses.

### 2.2 Thermal and chemical energy on the same faces

**T-B4.4.** The evolved third row stays `U(3) = e_kin + u_th` (B1 T1.1): no
flux, no EOS inversion and no boundary condition is redefined. The material
energy transport is therefore two terms on the same faces:

```text
thermal:   F_E(j) = F_rho(j) [ (E + p)/rho ]^face(j)      (the present flux)
chemical:  F_form(j) = sum_s eps_s F_s(j)
```

with `eps_s` the formation and excitation energy of B1 T1.2, from the single
declared reference (every element a neutral ground-state free atom at rest).
`F_form` is **not** added to the conserved row, and (advisor correction of
2026-09-06, replacing the first draft of this paragraph) it is **not a source
of the thermal energy equation either**. With thermal plus kinetic energy as
the evolved variable the two balances are

```text
d_t u_th   + div(F_E)    = Q_th      (photoheating h nu - I, cooling, reaction heats)
d_t u_form + div(F_form) = S_form    (I supplied by photoionization, I emitted by
                                      recombination, collisional exchanges with u_th)
```

and the ledger identity of B1 T1.5 is their sum: `d_t (u_th + u_form) +
div(F_E + F_form) = Q_external`. Advecting ionized gas out of a cell moves
`u_form`, not `u_th`; a `-div(F_form)` term in the thermal equation would
heat a cell by the ionization energy leaving it, which no process does.
`F_form` therefore belongs to the reservoir equation, which the code does not
carry as an evolved quantity: B4-2 builds it as a DIAGNOSTIC balance (the
reservoir's transport and its sources, evaluated on the state) so that the
identity above can be tested to round-off in a closed column and reported at
every step in `phys` mode (B3a decision 7). The thermal energy flux itself is
unchanged by `F_form`; what B4-2 changes in the thermal flux is T-B4.5 (the
face composition in the enthalpy).

**T-B4.5 (species enthalpy in the thermal flux).** The thermal face flux uses
the mixture enthalpy of the **face** composition, so that a face across which
`gamma_eff` or the mean molecular weight changes carries the enthalpy of what
actually crosses it and not of the cell it is charged to. This is the
composition dependence the caloric EOS already carries per cell
(`caloric_eos`, `adiabatic_index_from_state`); the face composition of T-B4.3
is what makes it available at a face.

### 2.3 The diffusive fluxes and the background counter-flux

**T-B4.6.** The diffusive velocities `{w_s}` come from the (B9) family of
Koskinen et al. (2022) written for **every** diffusing species, neutrals and
ions together, closed by the derived zero-current condition (J0) of B1 T4.1
and the ambipolar field (F) of B1 T4.2, and solved by one matrix inversion per
cell under the mass condition (B10). Three consequences for the operator:

1. **One inversion, not two arms.** The H/He binary and the trace metals are
   rows of the same matrix (B1 decision 5), so `sum_s rho_s w_s = 0` holds by
   construction at every face and the counter-flux is derived rather than
   assumed. The present sequential solve against a frozen hydrogen background
   survives only as an interim, and only with its residual
   `|sum_s rho_s w_s| / sum_s rho_s |w_s|` measured and reported.
2. **The proton is in the set with its charge** (B1 decision 6). Holding the
   most abundant ion at `w_s = 0` is not consistent with a field derived from
   a closure that sums over all charged species. Every species for which the
   implementation still sets `w_s = 0` carries the Peclet diagnostic
   `Pe_s = L |w| / (D_s + K_zz)` of B1 T4.6 in every cell, reported through B6.
3. **Faces, not cells.** `w_s` is evaluated at faces, from face values of the
   gradients (`d ln(n_e T)/dr` for the field, `d ln x_s/dr` for the Fickian
   term, `d ln psi/dr` for the carrier-count term), on the same `r_edg`
   positions. The present central-difference field on cell centers
   (`binary_element_diffusion`, `eEf`) becomes a face quantity, which is what
   makes the (J0) and (B10) sums face identities rather than cell identities.

### 2.4 One operator or two

**Recommendation: one operator with one face-flux assembly, two solvers.**

The advective part of every species is the same expression, `F_rho Y_s^face`,
and it belongs in `RK_rhs` beside the mass flux: it is explicit, it is at the
hydrodynamic time level, and it costs one multiplication per species per face.
The diffusive part is stiff (the base cell's molecular diffusion time is READ
from the `diffusive_photochemistry` header as 8.0e5 s with `K_zz = 0` against
a step of order 1 s) and must stay implicit.

So the target is:

- **one face-flux module** that owns `Y_s^face`, `F_s^adv = F_rho Y_s^face`,
  `w_s^face` from the (B9) inversion, `F_s^diff = A rho^face w_s^face / m_s`,
  and `F_form = sum_s eps_s F_s`. It is called once per RK stage for the
  advective part and once per implicit step for the diffusive part, and it is
  the only place any of these quantities is defined;
- **the explicit advective divergence inside the RK stages**, so that the
  species and the mass are advected by the same operator at the same time
  level and a uniform composition is preserved to the bit;
- **one implicit diffusive solve** for the whole species vector, replacing the
  present two (the element binary plus its trace arm, and the carrier block
  tridiagonal). The chemistry source stays inside that solve for the carriers,
  which is what the carrier Newton exists for; species with no kinetic model
  of their own contribute a row with a zero source.

The carrier operator does not disappear: it keeps its chemistry source, its
headroom constraint and its base handoff. What it loses is its own advection
term and its own notion of a velocity. The element diffusion operator loses
its advective form entirely, and with it the reason the code gives for that
form ("this routine receives only the new rho, v, T and dt"), because the face
mass flux becomes an argument.

### 2.5 Boundary equations by cell class, with the constraint count

The count is set by the characteristics for the hydrodynamic rows and by the
sign of the face mass flux for the material rows. `Ng` ghosts are written per
side; the ghosts are volume averages of the hydrostatic isentrope through the
face state and not copies of it, which the target keeps.

| class | rows | imposed at the face | count | checked, never imposed |
|---|---|---|---|---|
| base face `r_edg(0)`, hydrodynamic | 3 | LODI `C^-` compatibility plus the `(p, s)` reservoir, split by the eigenvalue signs: `0 < v < c` gives 2 reservoir and 1 interior; `-c < v < 0` gives 1 and 2; `v < -c` gives 0 and 3 | 3 | that the inflow is subsonic. `v > c` needs 3 reservoir conditions and 0 interior: the closure is then outside its stated domain, the state carries a category (2) validity record and certification is refused (B1 section 9.2, replacing the present diagnosis-only C33) |
| base face, species | one per independent species | **inflow** (`F_rho(0) > 0`): the full composition vector `Y_s^base` of the handoff, imposed on the advective flux, with `sum_s Y_s^base = 1` and the element and charge constraints checked on it; **outflow**: none, the interior composition is carried out | `n_indep` on inflow, 0 on outflow | that the imposed vector satisfies the element totals and charge neutrality |
| base face, diffusive | one per species | zero diffusive flux at face 0 when a handoff states the partition (the partition boundary is where the two models meet), or the reservoir Dirichlet where the handoff states a gradient | equal to the advective count, and never additional to it | `sum_s rho_s w_s = 0` and `sum_s z_s n_s w_s = 0` at face 0 |
| base ghosts, composition | | **no local re-solve.** The ghosts carry the imposed inflow composition, not a locally equilibrated one. Where the code today lets the two lower ghosts re-solve their ionization because the proton has no Dirichlet value (D0 C18), the target states the ghost closure explicitly and gives A2 its residual | 0 additional | |
| physical cells `1..N` | 3 hydrodynamic plus the active rows of B1 T1.7 | the equations of the physics | | the element totals, charge neutrality, the species mass sum against the supplied density. An anchored row (`set_base_fix` replacing the first `nfix` rows by `Y - Yfix` in the stationary solve) is an imposed condition and is declared as one in the certification record |
| outer face `r_edg(N)` | 3 plus the species | free outflow: no incoming characteristic, no condition imposed, on any row | 0 | that the state is outflowing. A state whose outer boundary is not outflowing is outside the closure's domain and carries a validity record; the WENO3 extrapolation dropping back to zero gradient stays a category (4) count |

**The count identity B4 must print** is: for each face, the number of imposed
conditions equals the number of characteristics entering the domain there for
the hydrodynamic rows, plus the number of species rows whose face flux is
inflowing. A configuration whose imposed set exceeds that count is
over-specified and is refused, which is D3 of rev 3 section 4.5.

### 2.6 The hydrostatic residual as an operator identity

**T-B4.7.** For a state that is hydrostatic in the discrete sense, the
operator returns zero velocity change to round-off, or the operator's residual
is stated as a number with its order in `dr`.

The discrete statement, on the cell `j` with the present assembly, is

```text
R_j = (A_p F_p(2) - A_m F_m(2))/dV  -  S(2)_j  -  [geometric pressure term]_j
```

with `F(2)` the momentum flux the Riemann solver returns from the reconstructed
face states, `S(2)_j = -0.5(rho_L + rho_R)(Gphi_i(j) - Gphi_i(j-1))/dr` and the
geometric term `(A_p - A_m) p_C/dV` under PLM or `(p_R - p_L)/dr` under WENO3.
Three things make `R_j` nonzero today and each is a separate object:

1. **the gravity source and the pressure divergence are divided by different
   lengths**: the source by `dr_j(j)`, the flux difference by `dV`. On a
   stretched grid `A_p - A_m` over `dV` is not `1/dr_j`, so the pair does not
   cancel at the same order the reconstruction claims;
2. **the face pressure comes from a reconstruction that is not
   well-balanced**: PLM of `p` on a stratified column reproduces a linear
   profile, not the exponential one hydrostatic balance requires, so the face
   pressure carries a truncation error of the order of `(dr/H)^2` that the
   gravity source does not carry;
3. **the gravity source uses the average of the two face densities** while the
   flux carries the reconstructed face pressures; the two are consistent only
   to the reconstruction's order.

**What B4 must produce**, in this order: the residual `R_j` written as a
diagnostic on the existing `hydrostatic_column` case at `Grid cells:` 250,
500, 1000 and 2000 on the same input; the observed order of
`max_j |R_j| / max_j |rho_j g_j|` and of the spurious velocity
`max_j |v_j| / c_s`; and, only if the observed order is below the
reconstruction's design order, a well-balanced form of the pressure-gravity
pair. **A spatial imbalance is not repaired by a time-step device**, which is
the standing instruction of rev 3 section 4.5 (D2): the base cycle on
`mol_sec_ion` is measured only after the residual test passes.

The present state of that case is READ from
`backup/regression/hydrostatic_column/README.md`: at 300 steps
`max |v|` over the physical cells is 2.8e5 cm s^-1 against a sound speed of
3.5e5 cm s^-1 at 1140 K, that is 0.8 of `c_s`, and the column drains rather
than settles (at 3000 steps `du` grows from 4.20 to 2.13e3 and the base
velocity reverses). **The case does not today have a hydrostatic solution to
measure a residual against**, and the README records as an open question
whether one exists under the present boundary treatment. B4's first act on
this case is therefore not the refinement ladder but the question the ladder
presupposes: does the discrete operator plus the base and outer conditions
admit a discrete hydrostatic state at all. Two candidate answers, both of
which the increment must distinguish by measurement rather than by argument:
the outer free-outflow condition drains the column because it imposes nothing
and the column has nowhere to rest, or the pressure-gravity pair itself drives
the flow. The discriminator is the residual `R_j` evaluated on the **analytic**
hydrostatic column (the isentrope the base reservoir is carried along,
evaluated as cell volume averages), which needs no converged run: if `R_j` is
already of the order of `rho g` there, the operator is the source and the
boundary question is secondary.

---

## 3. Increments

Each increment is bounded, carries its own test, and states which regression
cases it touches. The transport-active cases of the default matrix are READ
from their `input.inp`: `mol_carrier` (`Molecular carrier transport: True`),
`mol_diffusion` (`He_diffusion: True`, `He_Kzz: 1.0e9`) and `lower_profile`
(`He_diffusion: True`, `He_metal_diffusion: True`, the only case with the
trace-metal arm and therefore the only one carrying C29). No case in the
matrix sets an oxygen key, so the conditional carrier default
(`carrier_transport = thereis_oxychem`) never fires there. `wasp_full`,
`wasp_he23off`, the five `mol_*` snapshot cases and `hydrostatic_column` run
no transport operator at all.

### B4-1. Species face fluxes consistent with the mass flux

**Landed in two increments: the element rows 2026-09-06 (B4-1), the carrier
rows 2026-09-07 (B4-1c). Both are in the tree.** The first was briefed without
`diffusive_photochemistry.f90`, which another item held at the time, so the
element half landed alone and the carriers followed. The transported set is
now: helium against the hydrogen component (normalized), each trace metal when
`He_metal_diffusion` is on, and every carrier the photochemical transport
operator solves (H2 always, OH/H2O/CO under `Oxygen chemistry`, H+ under
`Ionization transport`) when `Molecular carrier transport` is on. A
configuration with neither element diffusion nor carrier transport has an
empty transported set and is byte-identical.

Build `Y_s^face` and `F_s = F_rho Y_s^face / m_s` in the face-flux module;
replace the carrier advective term and the element-diffusion advective term by
the divergence of `F_s`; delete the cell-velocity upwind selections and the
mixing-ratio formulation that existed to work around them.

Tests (all RED before, GREEN after):
- **passive scalar**: a uniform `Y_s` on a compressing and an expanding flow
  is preserved to round-off at every cell and every step, for both PLM and
  WENO3, at any velocity field including one whose two face velocities
  straddle zero;
- **face identity**: `|sum_s m_s F_s(j) - F_rho(j)| / |F_rho(j)| <= ` round-off
  at every face and every stage, asserted in the operator;
- **element budget in a transient**: each element's nucleus total over the
  domain changes only by its boundary fluxes, to round-off, over a transient
  that today conserves it only to hydro truncation (this is B1 AT-4d);
- **order**: a smooth composition step advected across the domain converges at
  the reconstruction's design order.

What landed, and where. A new module `src/modules/flux/species_face_flux.f90`
(`species_face_fraction`, `species_face_flux`, `species_advective_update`)
holds the face composition, the face flux and the Runge-Kutta stage update,
and asserts the face identity at every face of every stage. The scalar
reconstruction is `Reconstruct_scalar` in `Reconstruction.f90` and
`PLM_rec_scalar` in `PLM_rec.f90`, sharing the limiter `minmod_mc_slope` and
the ESWENO3 grid coefficients `weno3_geometry_coefficients` with the
primitive reconstruction, so the two cannot drift apart.
`binary_element_diffusion.f90` owns the transported set through
`species_advection_active`, `species_advection_begin_step`,
`species_advection_stage` and `species_advection_project`; the marching path
calls the first three inside the Runge-Kutta stages of `EXHALE_main.f90` and
the fourth at the element-transport split point, immediately before the
diffusive step. The element operator now advects only when it is GIVEN an
advecting mass flux, which the fixed-wind relaxation
`relax_element_composition` supplies and the marching path does not; the
stationary balance `element_transport_residual` keeps its own advective term,
because the two forms agree where `r^2 rho v` is constant and that is what a
stationary state means. `positivity_limited_fluxes` now writes a replaced
first-order flux back into `face_flux`, so the species rows ride on the flux
the mass row of the same stage actually used.

The trial operation. No new operation code was needed and none was added:
the stage update is part of the Runge-Kutta stages and is refused with them at
the `as_op_hydro` injection point, and the projection into `f_sp` sits in the
element-transport split slot and is refused at `as_op_diffusion`. The
composition is rebuilt from `f_sp` at the top of every attempt, so a restored
checkpoint restores it and no checkpoint field was added either. MEASURED on
`mol_diffusion` with `EXHALE_REJECT_AFTER_OP=4` at step 20: the run reports
`step 20 REFUSED at the element diffusion step` and goes on.

Golden movement, MEASURED 2026-09-06 on scratch copies with
`OMP_NUM_THREADS=1`, against a binary built from the same tree with the stage
calls removed and the element operator's cell-velocity advection restored.
The `mol_*` and `hp_*` cases were run at 200 steps and `mol_diffusion` and
`lower_profile` at 1000, not at the 12000 their `maxsteps` pins: one pair of
12000-step runs takes about 40 minutes at this grid and the matrix is 16
cases, so the step counts were cut and are stated rather than the runs being
skipped. `wasp_*` was capped at 300 as the brief directed.

BYTE-IDENTICAL, 14 cases, every one of them without element transport:
`wasp_full`, `wasp_he23off`, `wasp_full_newton`, `mol_base_handoff`,
`mol_metals`, `mol_lyman_werner`, `mol_ir_bands`, `mol_sec_ion`,
`mol_carrier`, `oxygen_chemistry`, `hp_zero_seed`, `hp_trace_seed`,
`hp_front`, `hydrostatic_column`. `mol_carrier` is in this list because the
carriers are B4-1c.

MOVED, the two cases that run element transport, as the largest relative
movement of any column of any row of each output file:

| case | file | max relative movement | where |
|---|---|---|---|
| `mol_diffusion` (1000 steps) | `Cooling_breakdown.txt` | 3.59e-3 | |
| | `Ion_species.txt` | 7.10e-5 | H3+ at `r = 1.104` |
| | `Heating_breakdown.txt` | 5.79e-5 | |
| | `Hydro_ioniz.txt` | 4.54e-5 | `v` at `r = 1.014` |
| `lower_profile` (1000 steps) | `Cooling_breakdown.txt` | 1.07e-2 | |
| | `element_flux_profile.txt` | 2.30e-3 | |
| | `Ion_species.txt` | 1.23e-3 | O III at `r = 1.077` |
| | `OI_levels.txt` | 1.21e-3 | |
| | `Heating_breakdown.txt` | 9.09e-4 | |
| | `Hydro_ioniz.txt` | 1.43e-4 | `cool` at `r = 1.076` |

The `_adv` files move with their parents and are not listed separately.

The transported quantity itself, the helium nucleus ratio `n_He/n_H`, moves by
1.9e-10 at the base, 3.5e-6 at the molecular front and 4.5e-8 at the top of
`mol_diffusion`, and by 2.7e-8, 1.3e-5 and 5.1e-10 at the same three places of
`lower_profile` (its largest movement, 6.0e-5, is at `r = 1.30`). The
base-layer mass flux `rho v r^2` moves by at most 9.2e-7 relative over
`1.00 <= r <= 1.16` in `mol_diffusion` and 2.5e-7 in `lower_profile`.

So the increment moves the two transport cases much less than the estimate it
replaces. The reason is READ from the runs and not from the estimate: at 1000
steps neither case has separated its elements far from the reservoir ratio, so
the composition the two forms advect is nearly uniform, and a uniform
composition is what both forms carry exactly. The estimate remains the right
statement of where the difference LIVES (the base layer, where the face flux
is not flat), and a converged 12000-step pair is what would measure it at the
separation those cases reach; that measurement is not in this increment and is
named here as open.

#### B4-1c, the carrier rows (landed 2026-09-07)

What landed. The carriers of `photochemical_transport_step` are declared to
the transported set by `carrier_set_init` once the input keys are parsed
(`advected_carrier_register` in `binary_element_diffusion.f90`), and are
advanced by a second call to `species_advective_update` inside the same three
Runge-Kutta stages, with `n_norm = 0`: a carrier is one species INSIDE an
element, the element is already carried by the closing member of the
normalized set, and normalizing the carrier beside it would count that mass
twice. The face identity is therefore a statement about the element rows, and
what keeps a carrier consistent with its element is the transport operator's
own write-back, which restores the element totals cell by cell.

The carrier rows are advanced from cell 1, not from cell 2 as the element rows
are: cell 1 is the element operator's Dirichlet reservoir, but the partition
of an element among its carriers there is a state the base face carries into.
The inflow composition is the inner ghosts', which is the handoff value where
one is stated (`q_H2_base`, pinned on the ghosts by the ionization sweep) and
the base cell's own partition where none is, so the one rule reproduces both
branches of decision D6 (B4-4a).

`diffusive_photochemistry.f90` lost the cell-velocity advective term, the
deferred correction `carrier_advection_correction` and the base Dirichlet
branch from its MARCHING rows, and with them the two named defects those
forms worked around. Two paths keep an advective term, each for a stated
reason: `carrier_steady_residual`, because the balance a stationary state
satisfies is the whole equation and not the remainder of an operator split;
and the fixed-wind relaxation `relax_photochemical_composition`, which takes
many transport steps with no hydrodynamic stage between them, so the advection
the marching loop takes with the mass row is not being taken anywhere else
there. The switch is `carrier_rows_advect`, and it is the same division the
element operator makes.

**Item B5b removed the last difference between those two paths and the
stages.** Where a row carries the advective term it is now the divergence of
the same face species fluxes: `carrier_mass_fractions` for the face variable,
`species_face_fraction` for the face value, `species_face_flux` for the flux
and `species_flux_divergence` for the divergence, which
`species_advective_update` also calls. The cell-velocity form, the deferred
correction `adv_corr`, `carrier_slope` and the base Dirichlet record
(`carrier_base_state`, `fc_base`, `base_dirichlet`, `base_inflow`) are gone
from the module: the direction at the base is the sign of the face mass flux
and the value is the ghost's. The block-tridiagonal Jacobian of the
fixed-wind relaxation carries the donor-cell linearization of the same term
(`carrier_advective_face_coefficients`); the limiter and the second-order
part of the reconstruction stay out of it, as the deferred correction did,
so what changed is where the Newton converges from and not what it converges
to. The stationary ELEMENT rows of `element_transport_residual` were moved
to the same face fluxes. `element_diffusion_step` and
`relax_element_composition` followed in item B5c: they take the face mass
flux itself (`Frho_in`) in place of the smoothed steady mass flux
`mdot/(4 pi r^2)` and form the term through `element_advective_divergence`,
which is the same three routines, so the elemental Picard alternation and the
elemental stationary row are one operator and the relaxation's fixed point is
the row's zero. `docs/steady_solver_design.md` section 15.

Tests. Two rows added to `src/tests/species_face_flux/`, MEASURED on a private
build: `uniform_carrier_partition_preserved` 0.83 ulp against a 8 ulp bound on
a sign-alternating face mass flux, and
`h2_nucleus_total_changes_only_by_boundary_flux` 1.36e-15 against 1e-12 over
the range that begins at the base cell. RED before in the strongest sense: the
driver does not compile against the pre-increment module. The suite is 15
driver rows plus 2 whole-binary rows, all GREEN. `carrier_retry` (101),
`carrier_returned_state_acceptance` (37), `carrier_constraint_attribution`
(11), `carrier_reference_scales` (15) and `acceptance_classes` (22) are GREEN
with no reference re-derived: their assertions are about the retry, the
restoration and the classification of a returned state, none of which was
built on the advective term.

Movement, MEASURED 2026-09-07 on scratch copies, `OMP_NUM_THREADS=1`, one run
at a time, against a binary built from the SAME tree snapshot with
`binary_element_diffusion.f90` and `diffusive_photochemistry.f90` at their
pre-increment text. Step counts: `mol_*` and `lower_profile` 200,
`mol_carrier` also 1000, `wasp_*` and `hydrostatic_column` 300, `hp_*` 100
(their own `maxsteps`), not the 12000 the `mol_*` cases pin.

BYTE-IDENTICAL, 12 cases, every one without `Molecular carrier transport`:
`wasp_full`, `wasp_he23off`, `wasp_full_newton`, `mol_base_handoff`,
`mol_metals`, `mol_lyman_werner`, `mol_ir_bands`, `mol_sec_ion`,
`mol_diffusion`, `lower_profile`, `oxygen_chemistry`, `hydrostatic_column`.

`oxygen_chemistry` is in that list for a reason that is not an absence of
effect and is worth stating: `carrier_transport` defaults to on with the
oxygen chemistry, but at 200 steps in initialization mode NOT ONE carrier
interval is covered, before or after (201 of 201 uncovered in both), so the
operator's write-back never runs and the carriers never reach `f_sp`. The
rows themselves did change, and the run log says so: the worst CO row at
`j = 310` of step 1 reads `|res| 6.03e-5` against physical terms `1.433e2`
before and `5.68e-5` against `1.431e2` after.

MOVED, the four cases that transport carriers, as the largest relative
movement of any column of any row of each file:

| case | file | max relative movement | where |
|---|---|---|---|
| `mol_carrier` (1000) | `Cooling_breakdown.txt` | 6.21e-1 | N I at `r = 1.100` |
| | `Ion_species.txt` | 4.97e-1 | H3+ at `r = 1.100` |
| | `Heating_breakdown.txt` | 3.67e-1 | `r = 1.100` |
| | `Hydro_ioniz.txt` | 2.65e-2 | `v` at `r = 1.186` |
| `mol_carrier` (200) | `Cooling_breakdown.txt` | 4.48e-1 | O II at `r = 1.000` |
| | `Ion_species.txt` | 2.29e-1 | H3+ at `r = 1.108` |
| | `Hydro_ioniz.txt` | 1.07e-2 | `cool` at `r = 1.000` |
| `hp_zero_seed` / `hp_trace_seed` (100) | `Cooling_breakdown.txt` | 1.68e-2 | N I at `r = 4.708` |
| | `Ion_species.txt` | 8.47e-3 | H3+ at `r = 4.708` |
| | `Hydro_ioniz.txt` | 1.35e-3 | `cool` at `r = 1.102` |
| `hp_front` (100) | `Cooling_breakdown.txt` | 1.13e-1 | O II at `r = 4.622` |
| | `Ion_species.txt` | 1.13e-2 | H3+ at `r = 4.708` |
| | `Hydro_ioniz.txt` | 1.80e-3 | `heat` at `r = 4.708` |

The transported quantity, `n(H2)`: on `mol_carrier` at 1000 steps it moves by
2.0e-3 at the base, 3.7e-1 at `r = 1.100` and 3.2e-4 at the top; on the `hp_*`
cases by 1e-15 at the base, 3.7e-4 at `r = 1.100` and 8e-3 to 1e-2 at the top.
The base-layer mass flux `rho v r^2` moves by at most 9.7e-4 relative over
`1.00 <= r <= 1.16` on `mol_carrier` (1000) and 3e-5 to 1e-4 on the `hp_*`
cases.

**Where the movement is, and it is the place the design named.** The H2 front
itself does not move: the `x_H2 = 0.5` crossing of `mol_carrier` sits at
`r = 1.0698` at 1000 steps before and after, and `x_H2(1.05)` is 0.7308
against 0.7310. What moves is the TAIL beyond the front, which the
cell-velocity form was carrying by its own numerical diffusivity:
`x_H2(1.10)` falls from 3.19e-3 to 2.02e-3, a 37 per cent reduction, and that
one number is what every large entry in the table above is a consequence of,
the molecular column beyond the front being what sets H3+, the metal cooling
of those cells and the local heating. Every carrier interval was covered in
both runs and the element limiter moved no cell in either.

No golden was refreshed.

### B4-2. Thermal and chemical energy on the same fluxes

Add the face composition to the thermal enthalpy flux (T-B4.5) and the
chemical flux `F_form = sum_s eps_s F_s` as a source of the thermal energy
equation (T-B4.4). The species energy table `eps_s` is B1 T1.2 and is defined
at exactly one site, shared with B3b's local ledger, so the transport and the
source cannot disagree about a species' formation energy.

Tests:
- **closed column material energy**: over a domain with no external exchange,
  `sum_j dV_j (u_th + u_form)` changes only by the boundary fluxes
  `F_E + F_form` at the two end faces, to round-off (this is the B4 row of the
  rev 2 test matrix, "material energy close");
- **advected front**: a composition discontinuity advected at constant
  velocity through a uniform-pressure gas produces no temperature change other
  than the one the mixture enthalpy implies;
- **the ledger meets B3b**: on one accepted step, the sum of the transport
  contribution and the source contribution to `u_form` equals the change in
  `sum_s n_s eps_s`, to round-off.

Expected golden movement: **every case with a spatially varying composition,
which is all of them**, including the atomic `wasp_*` pair. This is the
increment that changes numbers everywhere, and it does so because a term that
is physically present is absent today, not because a discretization is
refined. Size, ESTIMATED: the ratio of the chemical to the thermal enthalpy
carried by one hydrogen nucleus is `eps(H II)/((5/2) k T)`, that is
`13.6 eV / 2.15 eV = 6.3` at `T = 1e4 K` (READ: `eps(H II) = I(H)`, B1 T1.2).
The **divergence** is what enters, so the effect concentrates where the
ionization fraction changes along the flow, and it is a sink of `u_th` where
ionization increases outward. **The increment's first act is a diagnostic
measurement of `div(F_form)` against `heat - cool` on the existing accepted
states, reported before the term is switched on**, so that the advisor sees
the size before any golden moves. Whether the term defaults on is a physics
judgment on that measurement, not a golden-stability judgment.

### B4-3. The charged closure and the background counter-flux

One (B9) matrix inversion per cell face over all diffusing species, closed by
(J0) and (B10); the proton enters with its charge; the ambipolar field becomes
a face quantity; the trace metals leave their frozen-background arm.

Tests (B1 AT-4a to AT-4c):
- **zero net diffusive mass** on a static isothermal column with a trace
  metal: `|sum_s rho_s w_s| / sum_s rho_s |w_s|` at or below the stated order
  at every cell and step. RED today for a metals-on run by C29;
- **the current condition** on a column with several ion stages:
  `|sum_s z_s n_s w_s| / sum_s |z_s| n_s |w_s|` at or below the stated order;
- **the three ambipolar limits**: `dmeff` reproduces 2.9715, 2.4715 and 1.6477
  at `m_He/m_H = 3.9715` (this one is a check on the existing routine and is
  expected GREEN before and after);
- **the omitted order printed**: electron inertia, the thermal force and the
  resistive term evaluated as diagnostics, not assumed small (B1 T4.2).

Expected golden movement: `mol_diffusion` and `lower_profile` move; every
ionized-region profile of a diffusion-active run is expected to change because
the proton stops being held at `w_s = 0` (B1 decision 6 states this
explicitly and states that it is not a reason to defer). `mol_carrier` moves
only if it carries a charged carrier, which it does not
(`Ionization transport` is absent from its input). The eight cases with no
diffusion are expected byte-identical.

### B4-4. The boundary constraint count and the ghost composition

Print the count identity of section 2.5 for the active configuration; refuse
an over-specified set; replace the local re-solve of the two base ghosts'
ionization by the imposed inflow composition, with the ghost closure stated
and its residual available to A2; make the supersonic base inflow and the
non-outflowing outer boundary validity records rather than diagnostics.

Tests:
- **count**: a configuration with a transported proton and an inflowing base
  reports the expected number of imposed conditions, and an artificially
  over-specified set is refused with the face and the row named;
- **ghost composition**: the imposed vector satisfies the element totals and
  charge neutrality, and the two base ghosts no longer differ from the imposed
  partition (RED today by C18);
- **domain records**: a run driven supersonic at the base carries the category
  (2) record and is refused certification, instead of only printing a warning.

Expected golden movement: cases whose base face is inflowing with a handoff,
that is the `mol_*` family and `lower_profile`. Size, ESTIMATED: the ghost
composition enters through the reconstruction into the base face state, so the
change is of the order of the difference between the locally equilibrated and
the imposed ionization in the two ghost cells, which is largest in the
molecular base where the local equilibrium is the thing the handoff exists to
replace. `wasp_full`, `wasp_he23off` and `hydrostatic_column` have no
composition handoff and are expected byte-identical.

### B4-5. The hydrostatic residual

The diagnostic `R_j` of section 2.6, the analytic-column evaluation, the
refinement ladder on `hydrostatic_column`, and a well-balanced pressure and
gravity pair if the ladder says one is needed.

Tests:
- **analytic column**: `R_j` on the cell-volume-averaged hydrostatic isentrope,
  with its order in `dr` over the four grids;
- **the ladder**: `max_j |R_j| / max_j |rho_j g_j|` and
  `max_j |v_j| / c_s(T_eq)` both fall at the reconstruction's design order as
  `N` doubles, with no odd-even pattern (the pass condition the case README
  already states, and does not assert);
- **no regression in the wind**: a well-balanced form must reproduce the
  uniform-flow and shock-tube behavior of the present flux unchanged.

Expected golden movement: if a well-balanced correction lands, **every case
moves**, since the pressure-gravity pair is in every run. If the ladder shows
the present pair already at design order and the drain is the outer boundary,
nothing moves and the increment ends in a measurement and a statement.
`hydrostatic_column` is not in `golden/` today (its reference sits in
`baseline_post170_20260905/`), so it can be re-measured freely.

### B4-6. The cell-mean photon field

Rev 3 section 10.4 assigns to D1 the cell-mean photon field, that is the
inner-face column of `calc_column_dens` and the `k_lw` cell mean, as one
discretization for every absorber of one beam. It is listed here as the sixth
increment because it shares the face and volume conventions of the others, and
it is independent of B4-1 to B4-5: it can be briefed in parallel.

Expected golden movement: every case, since the column enters every optical
depth. Size, ESTIMATED: the difference between an inner-face column and a
cell-mean column is half a cell of absorber, so of the order of
`n_s dr_j / 2` against `N_col`, largest at the base where the column is built
and the cells are thin relative to nothing else, and negligible in the
optically thin wind.

### Order and dependencies

`B4-1 -> B4-2 -> B4-3` is a chain: the energy flux needs the species fluxes,
and the counter-flux needs a place to put a diffusive face flux. `B4-4` needs
`B4-1` (the inflow composition is imposed on a face flux that does not exist
yet). `B4-5` and `B4-6` are independent of all of them and of each other.

---

## 4. Test matrix (the B4 rows of PLAN rev 2 section 5)

| row | fixture | quantity | acceptance |
|---|---|---|---|
| B4-a diffusion-only column, mass | static isothermal column, no wind, diffusion on | `sum_j dV_j rho_j` and each element's nucleus total | change only by the boundary fluxes, to round-off |
| B4-b diffusion-only column, elements | the same column through a transient | each element's nucleus total | to round-off, not to hydro truncation (B1 AT-4d) |
| B4-c current condition | the same column with several ion stages | `|sum_s z_s n_s w_s| / sum_s |z_s| n_s |w_s|` | at or below the stated order at every cell and step (B1 AT-4b) |
| B4-d zero net diffusive mass | the same column with a trace metal | `|sum_s rho_s w_s| / sum_s rho_s |w_s|` | at or below the stated order (B1 AT-4a); RED today |
| B4-e reacting transport balance | a molecular base with an active chemistry source | per species, the sum of the time term, the flux divergence and the source | closes to round-off per accepted step, with the boundary fluxes named |
| B4-f closed column material energy | a domain with no external exchange, composition varying in space | `sum_j dV_j (u_th + u_form)` | changes only by `F_E + F_form` at the two end faces, to round-off |
| B4-g hydrostatic residual | `hydrostatic_column`, `Grid cells:` 250/500/1000/2000 | `max_j |R_j| / max_j |rho_j g_j|` | falls at the reconstruction's design order; no odd-even pattern |
| B4-h spurious velocity | the same ladder | `max_j |v_j| / c_s(T_eq)` | falls with the residual. Today, READ, 0.8 at 300 steps on the production grid |
| B4-i boundary count | one configuration per cell class, plus one deliberately over-specified | number of imposed conditions per face against the entering characteristics plus the inflowing species rows | equal; the over-specified set refused with the face and row named |
| B4-j face identity | every case that runs a transported species | `|sum_s m_s F_s - F_rho| / |F_rho|` per face per stage | round-off, asserted in the operator, not in the test |

Tolerances, units and floors are set in each implementation brief. "Round-off"
means round-off and is not used for a relative agreement of 1e-3, per the
vocabulary of rev 2 section 3.2.

---

## 5. Open choices for the advisor

1. **Do the equilibrium species keep no transport at all, or do they ride on
   the mass flux with the same reconstruction?** Today they are re-solved
   locally every step and the mass row is rewritten from them; the minimal B4
   makes the composition ride on the face mass flux and leaves the local solve
   to redistribute within an element, which is the only form in which
   `sum_s m_s F_s = F_rho` is a statement about all species rather than about
   the transported subset.
2. **Is the chemical energy flux `F_form` on by default when it lands, or
   measured first and switched on by a separate decision?** The recommendation
   above is to measure `div(F_form)` against `heat - cool` on the existing
   accepted states and report before switching, because the term is physically
   required and its size is not known.
3. **Does the advective species divergence go inside the RK stages, or stay a
   split operator at the end of the step?** Inside the stages gives one time
   level and a bitwise-preserved uniform composition; split keeps the RK stage
   cost and the existing operator boundaries unchanged.
4. **Where does the low-Mach dissipation stand in a conservative operator?**
   It is already a face flux with an identically zero mass component, so it
   needs no change to satisfy T-B4.2; the question is whether the species
   fluxes should carry a corresponding term at all, given that the mass one
   does not.
5. **Does the Shapiro filter stay?** It is not in flux form, carries no cell
   volume, and therefore breaks every conservation statement this design makes
   whenever it is armed. Either it becomes a conservative face-flux
   fourth-difference, or arming it becomes a validity record that excludes the
   state from conservation certification.
6. **Does the `_adv` post-process keep its own cell-centered fraction ODE, or
   read the transported composition?** Once the species ride on the face mass
   flux the marching solution already carries the advected composition the
   `_adv` files were built to approximate, and B1 section 8 governs the
   product; keeping both means two spatial operators for one quantity.
7. **Is a well-balanced pressure and gravity pair in B4's scope, or a separate
   item?** It touches every case and is the one increment that can move an
   atomic golden for a purely discrete reason.
8. **Interim or nothing for the trace metals?** B1 decision 5 allows the
   sequential frozen-background solve as an interim with its residual measured
   and its validity domain stated; the alternative is to refuse metal
   diffusion until the joint inversion lands, which takes `lower_profile` out
   of the matrix in the meantime.

## 6. Advisor decisions on section 5 (2026-09-06)

1. The composition rides on the face mass flux with the same reconstruction
   (the minimal B4 of section 2); the local equilibrium solve redistributes
   within each element and no longer rewrites `rho`. This is the form in
   which `sum_s m_s F_s = F_rho` is a statement about every species.
2. `F_form` is measured first (done, section 7): increment B4-2 builds the
   reservoir balance as a diagnostic that the ledger identity T1.5 tests,
   NOT as a thermal source (advisor correction above: with a thermal energy
   variable the formation flux is the reservoir's transport and enters no
   thermal equation). The thermal flux change of B4-2 is T-B4.5 alone (face
   composition in the enthalpy), measured and switched on with the B4 series.
3. The advective species divergence goes inside the RK stages: one time
   level, a uniform composition preserved bitwise.
4. The species fluxes carry no low-Mach dissipation term, matching the mass
   flux's identically zero component.
5. The Shapiro filter stays as it is (off by default); arming it becomes a
   validity record (B6 category 4) that excludes the state from conservation
   certification. A conservative rewrite is not scheduled.
6. Once B4-1 lands the `_adv` post-process reads the transported
   composition; its own cell-centered fraction ODE is retired under B1
   section 8 (a corrected product is either the same assembly or refused).
7. The well-balanced pressure and gravity pair is a separate increment
   (B4-5b), measured first on the analytic column of B4-5; it moves every
   case and is reported on its own.
8. Trace metals: the sequential frozen-background solve stays as the interim
   of B1 decision 5 with its residual measured and its domain stated;
   `lower_profile` stays in the matrix.

Sequencing note: B4-1 touches `EXHALE_main.f90`, the flux and RK modules, the
`System_implicit_adv_*` systems, the carrier operator and the element
diffusion operator; it is briefed after A0-impl (main), A1scale (carriers)
and A2 step 4 (diffusion residual) have handed those files back, and after
B3a's controller owns the adoption boundary, so that the new species step is
a trial operation from its first day.

---

## 7. B4-2 pre-measurement (2026-09-06)

Increment B4-2 states that its first act is a diagnostic measurement of
`div(F_form)` against `heat - cool` on the existing accepted states, reported
before the term is switched on, and the advisor's decision 2 of section 6
repeats it. This section is that measurement. Everything in it is MEASURED by
`src/utils/formation_energy_flux_diagnostic.py` on saved output directories;
no production source was changed, no binary was built and no golden was
refreshed. The `eps_s` values are READ, each with its site.

### 7.1 The `eps_s` table used

Reference of B1 T1.2: every element a neutral, ground-state, free atom at
rest; the free electron carries zero and the ionization energy sits on the ion.

| species | `eps_s` [eV] | READ from |
|---|---|---|
| H I, He I, neutral metals, e- | 0 | the reference |
| H II | +13.598435 | `e_th_HI`, `parameters.f90` (NIST ASD) |
| He II | +24.587389 | `e_th_HeI` |
| He III | +79.005154 | `e_th_HeI + e_th_HeII` |
| He 2^3S | +19.819614 | `parameters.f90`, 159855.9743 cm^-1 above the ground singlet (NIST ASD). The `HeI` output column is total He I with the triplet inside it, so this is excitation on top of a nucleus already counted |
| H(n=2) | +10.198826 | `e_th_HI (1 - 1/4)`; available only where `Excited_H.txt` exists |
| H2 | -4.478075 | `-D0(H2)`, `D0_H2_cm = 36118.11 cm^-1` (Huber & Herzberg 1979) through `mol_rates::h2_dissociation_energy_eV` |
| H2+ | +10.947852 | `IP(H2) - D0(H2)`, `IP(H2) = 15.425927 eV` (`e_th_H2`, NIST Chemistry WebBook) |
| H3+ | +4.771490 | `h(H2) + h(H+) - D0(H3+)`, `D0(H3+) = 35076 cm^-1` (Mizus et al. 2019, Mol. Phys. 117, 1663) |
| HeH+ | +11.753281 | `h(H+) - D0(HeH+)`, `D0 = 16448.84 - 1566.6764 cm^-1` (Coxon & Hajigeorgiou 1999, J. Mol. Spectrosc. 193, 306, Tables 9 and 8) |
| metal stage `k` | cumulative thresholds | `mion_ethr`, `species_table.f90` |

The four molecular entries are exactly `species_enthalpies` of
`molecular_reaction_heat.f90`, the code's single site for those energies,
written against this same reference. No case measured here runs the oxygen
network, so OH, H2O and CO have no entry, which is B1 T1.9's increment.

### 7.2 The operator and its control

The output files carry cell centers, so the diagnostic places the face at the
midpoint of two centers and takes the face flux as the average of the two
cell-centered values, then forms `(A_+F_+ - A_-F_-)/dV` with `A = r^2` and
`dV = (r_+^3 - r_-^3)/3`, the two ghost rows supplying the end cells'
neighbours. This is not the operator B4-1 builds, so what is measured is the
SIZE of the term, not its discrete value in an implementation.

The control is the same stencil on the mass flux, whose exact divergence
vanishes in a steady state. Reported as the relative non-flatness
`|A_+F_+ - A_-F_-| / (|A_+F_+| + |A_-F_-|)`, median over the energetic cells:

| case | floor (mass flux) | signal (`F_form`) | ratio |
|---|---|---|---|
| `wasp_full_newton` | 1.6e-7 | 7.2e-3 | 4.6e4 |
| `wasp_full` | 6.4e-4 | 5.3e-3 | 8.3 |
| `mol_base_handoff` | 1.1e-3 | 9.5e-3 | 8.8 |
| `mol_carrier` | 1.4e-3 | 1.1e-2 | 8.4 |
| `hd209` | 5.7e-4 | 1.1e-2 | 19 |
| `wasp52` | 3.7e-6 | 5.7e-3 | 1.5e3 |
| `lower_profile` | 1.1e-2 | 1.1e-2 | 1.0 |

Every case separates except `lower_profile`, whose saved mass flux is as
non-flat as its formation flux; its numbers are reported but discriminate
nothing and the case wants re-measuring on a Newton-finished state.
Statistics are taken over the cells with `heat >= 1e-4 heat_max`, since the
far field of `hd209` and `lower_profile` otherwise divides by heating rates of
order 1e-16 erg cm^-3 s^-1.

### 7.3 What was measured

The reservoir against the thermal energy density, at its peak: `u_form` is
5.7 times `(3/2) p` on `wasp_full_newton`, 33 on `hd209`, 30 on `wasp52`, 40
on `mol_base_handoff`, 17 on `lower_profile`. The section 3 estimate
`eps(H II)/((5/2) k T) = 6.3` at 1e4 K was the right size.

`|div F_form|` against the local rates, over the energetic cells:

| case | median /`heat` | median /`\|heat-cool\|` | > 0.1 of `\|heat-cool\|` | > 1 of `\|heat-cool\|` | max /`heat` |
|---|---|---|---|---|---|
| `wasp_full_newton` | 0.141 | 0.382 | 66% | 33% | 0.517 |
| `wasp_full` | 0.121 | 0.361 | 85% | 23% | 0.527 |
| `mol_base_handoff` | 8.40 | 10.0 | 92% | 76% | 329 |
| `mol_carrier` | 8.50 | 12.7 | 91% | 78% | 331 |
| `lower_profile` | 0.103 | 0.219 | 58% | 36% | 1168 |
| `hd209` | 0.145 | 0.263 | 55% | 44% | 109 |
| `wasp52` | 1.33 | 1.67 | 74% | 60% | 9.21 |

On `wasp_full_newton` the signed ratio `div F_form / heat` is **positive in
every one of the 500 physical cells** (minimum +1.06e-2, maximum +0.517).
Positive is a sink of `u_th`: the ionization fraction rises outward, the
outflow carries more formation energy out of a cell than into it, and the
heating term pays for that transport today with no equation saying so. The
same sign holds through the `wasp52` wind. The molecular snapshots change
sign several times, which is the H2 fraction rising and falling across the
dissociation front.

The ionization front is where it is largest, as section 3 expected. On
`wasp_full_newton` the steepest `x_HII` gradient is at `r = 1.2345 R_p`
(`x_HII = 0.495`) and there `|div F_form|/|heat - cool| = 4.41` and
`/heat = 0.515`, the maximum of both. By radial band, median `|div|/heat`:

| band | `wasp_full_newton` | `hd209` | `wasp52` | `mol_base_handoff` |
|---|---|---|---|---|
| 1.0-1.1 | 0.052 | 0.008 | 0.039 | 26.0 |
| 1.1-1.3 | 0.256 | 0.136 | 2.98 | 0.111 |
| 1.3-2.0 | 0.132 | 2.05 | 1.44 | 2.00 |
| > 2.0 | -- | 4.71 | 2.67 | 8.12 |

Small in the dense, nearly neutral base of an atomic run; largest at the
ionization front; of order unity or above in the outer wind, where the heating
has fallen off faster than the advected chemical energy.

Integrated over the physical cells (the same `dV` in all three sums, the
common `4 pi` omitted):

| case | `int div(F_form) dV` | `/ int heat dV` | `/ int cool dV` | `/ int (heat-cool) dV` |
|---|---|---|---|---|
| `wasp_full_newton` | 3.577e+25 | 0.284 | 0.329 | 2.05 |
| `wasp_full` | 3.736e+25 | 0.282 | 0.325 | 2.15 |
| `mol_base_handoff` | 5.597e+22 | 4.46 | 112 | 4.65 |
| `mol_carrier` | 4.386e+22 | 3.49 | 88.5 | 3.63 |
| `lower_profile` | 4.691e+19 | 9.6e-4 | 2.8e-3 | 1.5e-3 |
| `hd209` | 4.179e+22 | 1.01 | 0.814 | -4.15 |
| `wasp52` | 1.395e+24 | 2.31 | 14.5 | 2.75 |

`lower_profile`'s heat integral is dominated by its deep base, where the term
is 3 per cent of the heating; that and the control failure above are why its
column ratio is small.

The figure `heat`, `cool` and `|div F_form|` against `r` for `wasp_full` and
`mol_base_handoff` is written by the same script (`--plot-dir`). On
`wasp_full` the formation-flux divergence tracks the heating curve through the
ionization peak a factor 2 below it, then flattens and crosses above the
cooling curve beyond 1.4 R_p.

### 7.4 The present thermal flux, READ

`Num_Fluxes.f90`, `Phys_flux`, builds
`E = 0.5 rho v^2 + energy_density_from_pressure(jcell, rho, p)` and
`PF(3) = v (E + p)`. `energy_density_from_pressure` (`caloric_eos.f90`)
returns `p/(gamma_ad-1)` for an atomic cell and `n_k u(p/n_k)` through the H2
rovibrational ladder for a molecular one, so **the enthalpy of the internal
degrees of freedom of H2 is already in the present flux**. What is not in it
is a face composition: the caloric map is addressed by the owning cell
`jcell`, so a face across which `gamma_eff` changes carries the internal
energy of the cell it is charged to. That is the T-B4.5 defect, and it is a
different statement from the H2 degrees of freedom being absent, which they
are not.

The difference between `v(E+p)` and a monatomic `(5/2) p v`, as a divergence
against `|heat - cool|`:

| case | median | max | kinetic (max) | H2 rovibrational (max) | max `u_rv/(3/2 p)` |
|---|---|---|---|---|---|
| `wasp_full_newton` | 8.1e-3 | 0.735 | 0.735 | 0 | 0 |
| `hd209` | 6.4e-6 | 5.62 | 5.62 | 0 | 0 |
| `wasp52` | 4.2e-2 | 67.2 | 67.2 | 0 | 0 |
| `mol_base_handoff` | 0.747 | 202 | 31.7 | 207 | 0.577 |
| `mol_carrier` | 0.398 | 1456 | 34.1 | 1490 | 0.576 |

In the atomic runs the whole difference is the kinetic flux `(1/2) rho v^3`,
which the conserved row carries and a bare enthalpy flux would not. In the
molecular snapshots the rovibrational share reaches 58 per cent of `(3/2) p`,
so a monatomic enthalpy flux there would be wrong by a factor 1.6 in what the
flux carries. The `u_rv(T)` used in the comparison is read out of the same
302-level ladder the caloric EOS sums (`molecular_infrared_data.f90`,
`h2_lev_T` and `h2_lev_g`, Roueff et al. 2019, A&A 630, A58), parsed from the
Fortran source so that no second ladder exists.

### 7.5 Verdict for the increment

`div(F_form)` is the size of the formation-energy reservoir transport, a term
of the total-energy ledger that the code does not carry today (see 7.6: it is
NOT a thermal source). On a converged atomic state it is outward transport of
ionization energy in every cell, 14 per cent of the local heating in the median, half the
heating at the ionization front, and twice the net `heat - cool` over the
column; in the outer wind of `wasp52` and in the molecular base it exceeds the
local heating rate. This is not a refinement of a discretization: it is a
transport of energy the gas performs and the equation does not record, and the
heating term is silently paying for it. B4-2 switches it on, which is the
advisor's decision 2 of section 6, and the movement it produces is the
consequence of adding a term that belongs there.

### 7.6 Advisor reading of the measurement (2026-09-06)

The numbers of 7.1 to 7.5 are the size of the formation-energy RESERVOIR
transport, not of a missing thermal source. On `wasp_full_newton` 28 percent
of the absorbed heating leaves the column as advected ionization energy
(`div(F_form)` integrated), which is the familiar statement that part of the
photon energy an escaping atmosphere absorbs is carried away as ionization
potential rather than thermalized; it is a term of the total-energy budget
and of the ledger identity T1.5, and it is what B4-2's diagnostic balance
has to reproduce. The phrase "a physically required term the energy equation
omits today" in 7.x applies to the LEDGER (the code carries no reservoir
balance), not to the thermal equation, whose photoheating already carries
`h nu - I` per event. Where a molecular snapshot shows `div(F_form)` several
times the heating, the reservoir is the H2 bond energy advected through a
relaxing front, a transient of the snapshot.

---

## 8. B4-5 hydrostatic residual of the pressure and gravity pair, measured

**Status: measurement only. No production source was changed for this
section.** It is the first act section 2.6 asks for: the residual `R_j` on the
**analytic** hydrostatic column, which needs no converged run, taken before any
refinement ladder on the marched case and before any well-balanced correction
(B4-5b).

**Test.** `src/tests/grid_and_gates/hydrostatic_residual.f90`, registered as
`hydrostatic_residual` in that suite's `run.sh`. It links the production
objects and drives the production `define_grid`, `set_gravity_grid`,
`Apply_BC`, `Reconstruct`, `Num_flux` and `RK_rhs` + `Source`; nothing in it
re-implements a flux, a source or a grid.

**What is measured.** `R_mom(j) = dF(2,j) - S(2,j)`, the momentum the operator
adds in one right-hand-side evaluation to a state that should stay at rest,
normalized by the local weight `rho_j Dphi(r_j)` with the production `Dphi`.
The state is the isothermal hydrostatic column of the mechanical
`hydrostatic_column` case (0.0457 M_J, 0.49 R_J, `T_eq` 1140 K, He/H 0.0793,
`Domain mode: Spherical`, `r_max` 3, base grid 2.0e-4 x 50), in code units

```text
T = 1,   p = nhat rho,   rho(r) = exp[-(phi(r) - phi(1))/nhat],   v = 0
```

with `nhat = 0.820798` the particle count per unit mass of the neutral mixture
and `phi` the code's own potential, evaluated as cell volume averages by an
8-point Gauss-Legendre rule. MEASURED `b0 = 17.5734`, `H/R_p = 0.0467`,
`v0 = 3.0667e5 cm/s`, `c_s = 3.5869e5 cm/s`.

All numbers below are MEASURED 2026-09-06, the driver at `-O0` against the
production objects at `-O3`, `OMP_NUM_THREADS=1`. The two toolchains of the
tree give the same printed digits: conda-forge gcc 16.2 (the compiler `make`
and `run.sh` pick) and the system gfortran 13.1 agree on every number in this
section. ROE and HLLC agree to 1 percent or better in this state, except
`R(r=2)` under WENO3 on the two coarsest grids (6 and 2 percent, where that
residual is at its smallest), so one Riemann solver is quoted in the tables
and the other is in the run's own output.

### 8.1 The prescribed ladder: `Grid cells:` 250 / 500 / 1000 / 2000

`R/rho g` (ROE; the interior maximum is over cells 3..N-2, the physical cells
whose stencil reads no ghost):

| scheme | N | max\|R\| | rms | interior max | at r | R(1) | R(r=2) | R(N-1) | R(N) | max dv/c_s |
|---|---|---|---|---|---|---|---|---|---|---|
| PLM | 250 | 7.593e-1 | 4.99e-2 | 3.331e-3 | 1.0097 | 9.67e-5 | 1.624e-3 | -2.121e-1 | 7.593e-1 | 1.322e-4 |
| PLM | 500 | 7.528e-1 | 3.53e-2 | 1.194e-3 | 1.0100 | 1.000e-4 | 2.209e-4 | -2.358e-1 | 7.528e-1 | 1.296e-4 |
| PLM | 1000 | 7.510e-1 | 2.50e-2 | 4.439e-4 | 1.0101 | 1.013e-4 | 3.385e-5 | -2.445e-1 | 7.510e-1 | 1.288e-4 |
| PLM | 2000 | 7.504e-1 | 1.77e-2 | 1.566e-4 | 1.0102 | 1.018e-4 | 5.423e-6 | -2.480e-1 | 7.504e-1 | 1.285e-4 |
| WENO3 | 250 | 6.447e-2 | 4.29e-3 | 3.732e-3 | 1.0097 | 1.537e-4 | 1.803e-3 | 1.617e-2 | -6.447e-2 | 1.122e-5 |
| WENO3 | 500 | 2.322e-2 | 1.08e-3 | 1.333e-3 | 1.0100 | 1.590e-4 | 2.430e-4 | 6.120e-3 | -2.322e-2 | 3.999e-6 |
| WENO3 | 1000 | 8.791e-3 | 2.90e-4 | 4.952e-4 | 1.0101 | 1.610e-4 | 3.677e-5 | 2.364e-3 | -8.791e-3 | 1.508e-6 |
| WENO3 | 2000 | 3.282e-3 | 7.64e-5 | 1.749e-4 | 1.0102 | 1.618e-4 | 5.720e-6 | 8.804e-4 | -3.282e-3 | 5.620e-7 |

Observed orders (least squares over the four grids):

| quantity | length it is quoted in | PLM | WENO3 |
|---|---|---|---|
| `R` at `r = 2` (deep interior) | the local `dr_j(r=2)` | **1.961** | **1.977** |
| interior maximum | `1/N` | 1.466 | 1.468 |
| `max\|R\|` over the domain | `1/N` | **0.006** | 1.429 |
| `max dv/c_s` | `1/N` | 0.013 | 1.437 |
| `R(1)`, base cell | `dr_j(1)` | see 8.2 | see 8.2 |

The two lengths are separate because this ladder refines only one region:
`Grid cells:` grows while `Base grid [dr,cells]: 2.0e-4 50` is held, so
MEASURED `dr_j(1)` moves only from 1.894e-4 to 1.994e-4 R_p (5 percent) while
`dr_j(r=2)` falls from 2.873e-2 to 1.565e-3 (18x). **The ladder the case
README states does not refine the base**, which is why the base cell is
measured separately below.

### 8.2 A base ladder, `Base grid [dr,cells]:` 2e-4 50 / 1e-4 100 / 5e-5 200 / 2.5e-5 400 at N = 1000

`R(1)/rho g`: PLM 1.013e-4, 5.027e-5, 2.498e-5, 1.237e-5 (order **1.004** in
`dr_j(1)`); WENO3 1.610e-4, 8.007e-5, 3.983e-5, 1.974e-5 (order **1.002**).
The base cell is first order in its own width, and the WENO3 value at
`drc = 2e-4` (1.61e-4) sits next to the 1.749e-4 that
`docs/phaseC_characteristic_base_bc.md` section 7.2 measured for the
characteristic condition on the hot-Uranus geometry by a different route.
The characteristic base condition is consistent and is not where the imbalance
is.

### 8.3 The verdicts

**(i) The pressure and gravity pair is at PLM's design order and below
WENO3's.** In the deep interior, away from both boundaries and away from the
grid junction, `R` falls as `dr^1.96` under PLM and as `dr^1.98` under WENO3.
Second order is the design order of the PLM pair; it is one order below the
third that section 2.6 sets for WENO3. The pair is therefore **not** the
source of an O(1) hydrostatic imbalance anywhere in the interior, and a
well-balanced correction (B4-5b) would buy accuracy in the WENO3 branch, not
consistency.

**(ii) [Measured before B4-4a, 2026-09-06; closed the same day: the outer ghost is now the isothermal hydrostatic continuation of cell N and the outer-cell residual falls at order 1.6 to 1.8, see `Update_EXHALE.md` section 7 item B4-4a.] The O(1) residual is at the OUTER free-outflow boundary, and under PLM
it does not converge.** `max|R|` sits at cell N at every grid and for every
scheme. Under PLM it is 0.750 of the local weight and grid independent (order
0.006), with cell N-1 at 0.248; under WENO3 it starts at 6.4e-2 and falls at
order 1.43. The mechanism is READ from the code and is consistent with that
split: `Apply_BC_W`'s free-outflow ghost is a zero-gradient copy of cell N
under PLM, so the forward difference entering the MC limiter of `PLM_rec` is
exactly zero, a zero enters the limiter's argument list and the limited slope
of cell N is exactly zero, leaving that cell with no pressure gradient of its
own; under WENO3 the ghost is the linear extrapolation, which keeps one.

**(iii) The largest INTERIOR residual is at the grid junction, not at the
scheme's smoothest place.** It sits at `r = 1.0100` at every N, which is
`1 + N_low*drc = 1 + 50*2.0e-4`, the cell where the uniform base region meets
the stretched region and the width jumps. It falls at order 1.47 in `1/N`
because the first stretched cell approaches `drc` as the stretched region
refines, not because the scheme's order is 1.47.

**(iv) The spurious velocity of one CFL step.** `max_j |dv_j|/c_s` with
`dv = R_mom dt_CFL/rho` is 1.32e-4 to 1.28e-4 under PLM over the whole ladder
(47.4 to 46.1 cm/s), and 1.12e-5 to 5.62e-7 under WENO3 (4.03 to 0.20 cm/s).
Its location is the outermost cell in every case. Under PLM it does not fall
with resolution, which is the same statement as (ii). For scale, the case
README records `max |v| = 2.8e5 cm/s` (0.8 `c_s`) on the marched state at 300
steps: a 46 cm/s kick per step, applied at the outer cells and never
converging away, is a candidate source of that drain, and the base cell at
1e-4 of its own weight is not.

**(v) The discriminator of section 2.6 is answered.** The two candidates it
names were "the outer free-outflow condition drains the column" and "the
pressure-gravity pair itself drives the flow". On the analytic column, MEASURED:
the pair is second order everywhere in the interior and the base cell is first
order in its own width, while the outer cell carries 0.75 of its own weight
at every resolution under the production PLM stage. **The outer condition is
where the O(1) residual is.** B4-5b (a well-balanced pair) is not what removes
it; the outer boundary closure of B4-4 / D3 is.

Two limits of this measurement, stated so it is not read for more than it
says. It is the truncation error of the operator on the exact solution, one
right-hand-side evaluation, not the marched trajectory: the drain also
involves the time integration and the two conditions in feedback. And the
column carries no source term, by construction (the case switches the
irradiation down 18 orders and the cooling has no key), so nothing here speaks
about the thermal balance the same case's README records.

### 8.4 An observation outside this item's scope

`define_grid` gives the Mixed grid's stretch factor a Newton iteration whose
convergence flag is declared `real*8 :: tol = 1.0`. An initializer in a
subroutine declaration implies SAVE, so on the SECOND and every later call the
loop is skipped and the stretch parameter keeps its initial guess 1.01,
whatever `N` and `r_max` ask for. Production calls `define_grid` exactly once
(`init.f90`), so no run is affected; anything that builds two grids in one
process is. MEASURED, before the ladder was moved to one grid per process: a
second grid at N = 2000 came back with `dr_j(1) = 7.3e-11 R_p` instead of
2.0e-4. Not fixed here: `define_grid.f90` is not this item's file.
