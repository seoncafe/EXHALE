# LHS 1140 b, item L12b: the ionization-stage flux, its discrete identity, and its boundaries

Item L12 of `docs/PLAN_20260913_lhs_stationary.md`, section 8 of
`docs/PLAN_20260916_rev3.md` (rows R17, R18, R41 of its section 0).  This
memo corrects section 1.3 of `docs/lhs1140b_stationary_L12a_design_20260913.md`
and prepares stages B to F.  **No file under `src/modules/` or `src/EXHALE_main.f90`
is changed by it.**  The one code deliverable is the acceptance suite
`src/tests/ionization_stage_flux/`, which tests the formula derived here and
the face routines the production operator itself calls.

Every number is MEASURED (computed or run here) unless marked READ.

---

## 0. Verdict, before the detail

1. **The corrected stage flux is**

   ```
      F_k(f) = x_k(f) N_el(f)  -  n_el(f) K(f) [ x_k(j+1) - x_k(j) ] / dr(f)
   ```

   for every stage k of an element, neutral stages included, with `x_k` the
   fraction of that element's nuclei in stage k, `N_el` the element NUCLEUS
   face flux the element operator forms (advection plus gradient plus eddy
   plus drift), and `n_el(f)`, `K(f)`, `dr(f)` the face element nucleus
   density, the face eddy coefficient and the face spacing.  The design's
   `Phi_x = -n_tot K d(x n_el/n_tot)/dr` is the whole mixing-ratio eddy flux
   and counts `x_k J_el,eddy` a second time, because `x_k N_el` already holds
   it.  Section 2 derives this in the operator's own variable.

2. **The identity `sum over all stages k of F_k(f) = N_el(f)` holds to
   rounding at every face if and only if three conditions hold**: one element
   face flux `N_el(f)` multiplies every stage; the face fractions close the
   simplex exactly, which they do when the closing (neutral) stage is not
   reconstructed but taken as one minus the reconstructed others; and the
   closing stage's eddy term is taken as MINUS the sum of the others'.  The
   third condition is not cosmetic: forming it from the closing stage's own
   cell values is exact in exact arithmetic and loses the gradient to the
   rounding of `1 - sum` wherever the carried stages are trace, MEASURED at
   2.9e-3 of the term (section 4).  Section 3 states the conditions, section
   4 the measurements.

3. **The test is implemented**: `src/tests/ionization_stage_flux/`, 15 rows,
   15 passed and 0 failed (MEASURED).  It tests the formula, the face
   fractions and the telescoping, NOT stage B, which does not exist.  Its
   four RED rows measure 5.0e-3 (independent reconstruction, PLM), 8.6e-4
   (the same, WENO3), 1.99 (a second copy of the element face flux) and
   2.0e-2 (a stage-dependent face density).

4. **Stage B needs one accessor the element operator does not expose.**
   `drift_and_gradient_face_coefficients`, `element_face_flux` and
   `settling_coefficient` are private to `binary_element_diffusion`, so no
   caller can read `Agrd`, `Bdrf`, `updrf` or `J(f)`, and the element NUCLEUS
   face flux is formed nowhere: the operator forms two divergences and never
   a single face flux.  Section 7 names what has to be added.

5. **The neglect of drift between charge stages is not small on this wind,
   and this memo does not certify it.**  MEASURED on the certified 45 R_p
   state with the code's own non-resonant ion-neutral coefficients and the
   code's own ambipolar field, the H II against H I drift reaches 0.72 of the
   bulk velocity above 1.2 R_p and the He II against He I drift 0.53.
   Resonant charge exchange, the dominant momentum-transfer channel for both
   pairs, is absent from those coefficients and would roughly halve both
   (READ: resonant H+ + H at 1e4 K about 2e-15 cm^2 against the about
   1e-15 cm^2 rigid core, `docs/collisional_validity.md` section 1).  The
   approximation is safe where the local root is returned anyway, that is
   inside about 1.2 R_p where the Damkohler number is 1e4 and above, and
   degrades outward in step with the continuum assumption itself, which the
   run's own tool already reports as unvalidated.  Section 6.

---

## 1. The operator's own convention, read from the source

`src/modules/functions/binary_element_diffusion.f90` (all READ):

- **The transported variable is a MASS fraction**, `X = rho_He/rho`, the
  helium mass over the mixture mass, with the element counted over every
  species through the `bsp_nH`/`bsp_nHe` weights (module header, STATE).
  Component 1 is hydrogen together with the trace metals and the heavy nuclei
  bound in the molecular carriers, and `rho_1 + rho_He = rho` exactly.
- **The transport equation** is
  `d(rho X)/dt + (1/r^2) d/dr [ r^2 ( rho X v + J ) ] = 0` with
  `J = -rho (D_12 + K_zz) dX/dr - rho D_12 G X (1-X)` (module header,
  TRANSPORT and DIFFUSIVE MASS FLUX).
- **The face coefficients** (line 2009 and its neighbors,
  `drift_and_gradient_face_coefficients`) are

  ```
     dr_f    = max(rp(j+1) - rp(j), 1)
     rhof    = 0.5 (rho(j) + rho(j+1))
     Df      = 0.5 (D(j) + D(j+1))
     Gf      = 0.5 (G(j) + G(j+1))
     Kf      = 0.5 (kzz_cell(j) + kzz_cell(j+1))
     Agrd(j) = rhof (Df + Kf)/dr_f        gradient + eddy, >= 0
     Bdrf(j) = rhof Df Gf                 drift, either sign
  ```

  every face quantity an arithmetic average of its two cells, and
  `updrf(j)` the Peclet switch: 0 (central) where
  `|Bdrf| dr_f <= 2 rhof (Df + Kf)`, -1 where the drift is inward and +1
  where it is outward.  **Faces 0 and N are left at zero**, which is the
  zero-diffusive-flux boundary at both ends.
- **The diffusive face flux** (`element_face_flux`) is
  `J(f) = -Agrd (X_r - X_l) - Bdrf [X(1-X)](f)`, the two factors of the drift
  product taken from opposite sides of the face.
- **The advective half is not in this operator's face coefficients.**  It is
  the divergence of the SAME face mass fluxes the density rides on, formed by
  `species_face_fraction`, `species_face_flux` and `species_flux_divergence`
  of `src/modules/flux/species_face_flux.f90`
  (`element_advective_divergence`).  The face mass fraction is the
  reconstruction (PLM or WENO3) evaluated on the side the face MASS FLUX
  selects, `F_rho >= 0` takes the left state, scaled back onto [0,1] against
  its own donor cell average.
- **There is no single element face flux anywhere in the operator.**  The
  cell row adds two divergences with two different discrete geometries: the
  diffusive half through `Kj (sR J(j) - sL J(j-1))` with
  `Kj = 1/(rp(j)^2 (rep(j) - rep(j-1)))`, the advective half through
  `species_flux_divergence`, whose volume is `(r_+^3 - r_-^3)/3`
  (`composition_residual`).  The two agree to `O(dr^2)` and are not the same
  number.
- **The eddy coefficient acts on the MASS fraction**, not on the mole
  fraction: "eddy mixing transports the mixture as a whole and has no
  preferred species, so it acts on `dX/dr` alone" (module header, the
  paragraph on `-dln(psi)/dr`).  This is a stated choice of the operator and
  the stage term below is written to match it.

The present proton carrier (`src/modules/lower_atmosphere/diffusive_photochemistry.f90`,
all READ) carries `f_sp(:,isp_HII)`, a fraction per unit mass, on the BULK
face mass flux; its molecular diffusion is set to zero
(`carrier_diffusivities`, the `ic .eq. ic_Hp` branch) while `kzz_cell` still
enters through `Agrd(j,ic) = ntf (Df + Kf)/dr_f`, and `carrier_face_flux`
turns the mass fraction into the particle mixing ratio `X = wc f` before
differencing it.  So the proton's present eddy flux is exactly
`-n_tot(f) K(f) d(n_HII/n_tot)/dr`, the `Phi_x` of the L12a design.  **It is
correct today** because the bulk mass flux carries no eddy term of its own.
It stops being correct the moment the advective half becomes the element
flux, which does.

---

## 2. The correction, derived in the operator's variable

### 2.1 In number densities

Write `n_k` for the density of stage k of an element, `n_el = sum_k n_k` the
element nucleus density, `x_k = n_k/n_el`, `y = n_el/n_tot`.  The eddy flux of
a species is the mixing-ratio flux `-n_tot K d(n_k/n_tot)/dr`, so with
`n_k/n_tot = x_k y`,

```
   -n_tot K d(x_k y)/dr = -n_tot K y dx_k/dr  -  x_k n_tot K dy/dr
                        = -n_el K dx_k/dr     +  x_k Phi_el,eddy ,
```

where `Phi_el,eddy = -n_tot K dy/dr` is the element's own eddy flux.  The
element flux `N_el` the stage rides on already contains `Phi_el,eddy`, so the
stage flux is

```
   F_k = x_k N_el  -  n_el K dx_k/dr ,
```

and the additional stage term is `-n_el K dx_k/dr` alone.  Summing over all
stages, `sum_k x_k = 1` and `sum_k dx_k/dr = 0`, so `sum_k F_k = N_el`
identically: the ionization operator moves charge and cannot move a nucleus,
which is the property section 1.2 of the L12a design asks for.

### 2.2 In the operator's mass-fraction variable

The operator's gradient coefficient acts on `X`, so the same statement has to
be made with mass fractions.  All stages of one element share a nucleus mass
to within the electron mass, which is the mass the species table assigns
(`bsp_mass`), so the mass fraction of stage k is `X_k = x_k X_el` exactly in
the code's own weights.  Then

```
   -rho K dX_k/dr = x_k ( -rho K dX_el/dr )  -  rho X_el K dx_k/dr
                  = x_k J_el,eddy            -  rho_el K dx_k/dr ,
```

and dividing by the nucleus mass `m_el` returns section 2.1 with
`n_el = rho_el/m_el`.  **The mass-fraction and nucleus-fraction conventions
agree on the stage term**, and they agree for the same reason: the stages of
one element have one mass.  So the stage term is

```
   E_k(f) = -( rho_el(f)/m_el ) K(f) [ x_k(j+1) - x_k(j) ] / dr(f)
          = -n_el(f) K(f) [ x_k(j+1) - x_k(j) ] / dr(f)
```

with `n_el(f) = 0.5 (n_el(j) + n_el(j+1))`, `K(f) = 0.5 (kzz(j) + kzz(j+1))`
and `dr(f) = max(rp(j+1) - rp(j), 1)`, the same three face rules
`drift_and_gradient_face_coefficients` uses.

### 2.3 The three carried fractions of stage B

```
   hydrogen   x(H I)   x(H II)                      N_H  = element H nucleus flux
   helium     x(He I)  x(He II)  x(He III)          N_He = element He nucleus flux

   F(H II)   = x(H II)(f)   N_H(f)  + E(H II)(f)
   F(H I)    = x(H I)(f)    N_H(f)  + E(H I)(f)
   F(He II)  = x(He II)(f)  N_He(f) + E(He II)(f)
   F(He III) = x(He III)(f) N_He(f) + E(He III)(f)
   F(He I)   = x(He I)(f)   N_He(f) + E(He I)(f)
```

The three carried rows are `x(H II)`, `x(He II)` and `x(He III)`; `x(H I)` and
`x(He I)` are the closing stages and are not rows.  The He 2^3S level is a
sublevel inside He I and is not a stage of this partition (`species_table.f90`
`bsp_is_excited_level`, READ); it stays with the sweep, for the Damkohler
reason of L12a section 3.

The element nucleus flux is the element MASS flux divided by the element's
mass per nucleus:

```
   N_He(f) = [ F_rho(f) Y_He(f) + J(f) ] / m_He
   N_H(f)  = [ F_rho(f) Y_1(f)  - J(f) ] / m_1(f)
```

with `Y_1 = 1 - Y_He` the closing member of the mass-fraction pair and
`J_1 = -J_He` the binary closure of the module header.  `m_He` is a constant;
`m_1(f)` is the mass of component 1 per hydrogen nucleus at the face, which is
constant only in an atomic, reservoir-composition gas and is a face quantity
elsewhere.  **`m_1(f)` has to be ONE face value shared by every hydrogen
stage**, for the same reason `n_el(f)` does (section 3).

---

## 3. The exact discrete identity, and the three conditions

**Statement.**  At every face f, with the stage fluxes of section 2.3,

```
   sum over ALL stages k of F_k(f)  =  N_el(f)          to rounding.
```

**Condition C1, one element flux.**  The same number `N_el(f)` multiplies
every stage.  Two copies of it, for instance one rebuilt per stage from its
own donor rule, do not return one element flux.  MEASURED break: 1.99 of the
term (section 4).

**Condition C2, the face fractions close the simplex exactly.**  The face
interpolation must be a linear map whose coefficients sum to one and are the
same for every stage: one stencil, one donor side taken from the sign of
`N_el` and not of each stage's own gradient, one limiter factor.  The
production reconstruction is NOT such a map, because
`species_face_fraction` applies a bound-preserving scaling toward the donor
cell average whenever a reconstructed face value leaves [0,1], and that
scaling is chosen species by species.  **What restores C2 is the rule the
code already follows for the mass fractions**: the closing member is not
reconstructed, it is one minus the reconstructed others
(`species_advective_update`, "The closing member of the normalized set is not
reconstructed ... which is what makes the sum exact", READ).  With that rule
`sum_k x_k(f) = 1` to one subtraction's rounding, MEASURED 2.2e-16.
Reconstructing every stage on its own instead leaves the sum away from 1 by
the reconstruction's own truncation, MEASURED 5.0e-3 (PLM) and 8.6e-4
(WENO3).

**Condition C3, the stage term telescopes.**  `n_el(f)`, `K(f)` and `dr(f)`
are one face value shared by every stage, and **the closing stage's term is
formed as minus the sum of the others**, `E_close = -sum_k E_k`.  Both halves
are needed and for different reasons.  A stage-dependent face density (a
donor cell chosen from each stage's own gradient) breaks the sum outright,
MEASURED 2.0e-2.  Forming the closing stage's term from its own cell values
by the same formula is exact in exact arithmetic, and in floating point it
differences numbers near 1 while the physical difference can be far below
`eps`: in the neutral-background limit the residual is MEASURED at 2.9e-3 of
the term.  This is the same failure the face fractions have and it takes the
same remedy.

**What the identity does NOT require.**  It is linear in the diffusive half
of the element flux, so it holds for any value of `J(f)`, including a
reversed and a vanishing one.  That is what lets the test below state the
whole mathematical content without reading the operator's private face
coefficients, and it is also why C1 is the only condition the diffusive half
imposes.

**At the divergence level.**  The cell row must apply the face fractions to
each of the element row's two divergences separately, with the geometric
factors each already carries (section 1, last bullet): `x_k(f)` on the
advective face flux inside `species_flux_divergence`, `x_k(f)` on `J(f)`
inside `Kj (sR ... - sL ...)`, and the stage term `E_k` through whichever of
the two the implementation puts it in.  Summing over stages then returns the
element row's own two divergences and element conservation is exact by
construction and not by cancellation.  Forming one face element flux by
adding the two halves and then taking one divergence of it would be a third
discretization of the element balance, and the operator's own header records
what two spellings of one term cost.

**Charge consistency.**  After transport,

```
   n_e = sum_k Z_k n_k = x(H II) n_H + [ x(He II) + 2 x(He III) ] n_He
         + the frozen metal and molecular-ion charges
```

is a sum over charges of the transported partition, not a separate variable,
which is what makes the replacement of one term by its trial value exact.
The proton case of this already exists (`carrier_source`,
`n_e = cbg_ne - cbg_nhii + n_hii`, READ) and the generalization is the same
sum over the three charges.  MEASURED consistency of the rebuilt sum: 2.6e-16.

**The simplex, and what an upwind stage flux does at a simplex face.**  Each
`x_k` belongs to [0,1] and they sum to one, so the admissible set is the
simplex, not a box, and the element headroom of L12a section 2.2 cannot be
reused: an ionization fraction is bounded by the sum rule and not by a free
nucleus density.  The closing-member rule buys the SUM exactly and buys
nothing else: where the reconstructed members are separately scaled back into
[0,1], the closing member can still leave it.  MEASURED largest negative
closing face value on the simplex-face column: -8.6e-4.  An upwind stage flux
at such a face donates from one cell, so a donor at `x_k = 0` sends nothing
of stage k and an acceptor at `x_k = 1` receives nothing more, exactly as the
element operator's drift product does at the ends of its own axis; what the
upwind rule does NOT give is nonnegativity of the closing stage, because that
stage is a difference and not a donated quantity.  **Stage B therefore owes a
simplex projection on the cell state** (scale the carried fractions by the
largest factor that keeps their sum at or below one), measured and reported
rather than asserted, in the shape `he_fraction_over_one` and
`species_fraction_excursion` already have.

---

## 4. The test

`src/tests/ionization_stage_flux/ionization_stage_flux_tests.f90` and its
`run.sh`.  It links the production objects and sets the grid itself, so no
`input.inp`, no ionization solve and no hydrodynamics are involved.  The face
fractions come from `species_face_fraction` and the element advective face
flux from `species_face_flux`, the two routines the element operator's own
advective half calls; the diffusive half is supplied by the driver, for the
reason of section 3, and is run positive, reversed and absent.

**Columns.**  A helium element that separates from hydrogen by a factor 4
over the column, as the LHS 1140 b wind does, crossed with four compositions
(all stage fractions varying; all constant; the neutral-background limit at
1e-12; helium driven onto the simplex faces, `x(He I)` to 0 and `x(He III)` to
1), four face mass fluxes (outflowing, reversed, alternating in sign face by
face, identically zero), three diffusive halves and both reconstructions.
`K_zz = 1e9 cm^2/s`, nonzero everywhere, filled through the production
`eddy_diffusion_on_grid`.

**Result (MEASURED, `OMP_NUM_THREADS=1`, private build `build_L12`): 15 rows,
15 passed, 0 failed.**

| row | measured | reference | tolerance |
|---|---|---|---|
| `stage_fractions_close_the_simplex_at_every_face_PLM` | 2.2204e-16 | 0 | 8 eps |
| `stage_fractions_close_the_simplex_at_every_face_WENO3` | 2.2204e-16 | 0 | 8 eps |
| `sum_of_stage_fluxes_is_the_element_flux_PLM` | 3.7323e-16 | 0 | 100 eps |
| `sum_of_stage_fluxes_is_the_element_flux_WENO3` | 2.8128e-16 | 0 | 100 eps |
| `stage_eddy_term_sums_to_zero` | 0.0000e+00 | 0 | 8 eps |
| `charge_consistency_of_the_transported_partition` | 2.5941e-16 | 0 | 4 eps |
| `the_closing_member_excursion_is_finite` | 8.6300e-04 | reported | not gating |

**The RED rows, each a deliberately mismatched rule that must fail the
identity (MEASURED):**

| row | measured |
|---|---|
| `independently_reconstructed_stages_break_the_simplex_PLM` | 5.0251e-03 |
| `independently_reconstructed_stages_break_the_simplex_WENO3` | 8.6300e-04 |
| `independently_reconstructed_stages_break_the_identity_PLM` | 5.0251e-03 |
| `independently_reconstructed_stages_break_the_identity_WENO3` | 8.6300e-04 |
| `a_second_element_face_flux_breaks_the_identity_PLM` | 1.9899e+00 |
| `a_second_element_face_flux_breaks_the_identity_WENO3` | 1.9899e+00 |
| `differencing_the_closing_member_loses_the_trace_gradient` | 2.9061e-03 |
| `a_stage_dependent_face_density_breaks_the_eddy_sum` | 2.0325e-02 |

The RED rows are not decoration: **the first form of the closing stage's eddy
term written here failed the identity at 8.1e-3** and the 2.9e-3 row is what
that failure was localized to.  The rule `E_close = -sum_k E_k` of condition
C3 is a result of this test and not an assumption that went into it.

**What the suite does not cover.**  That stage B multiplies the operator's
own `J(f)` rather than a second copy of it, and that it applies the face
fractions to the element row's two divergences with their own geometric
factors.  Both need the accessor of section 7, and both become rows of this
suite when it exists.

---

## 5. Boundaries, and the matching to the inner equilibrium

**The base face, f = 0.**  The element operator sets the diffusive flux to
zero at faces 0 and N by its face coefficients (the loop runs `j = 1, N-1`),
and cell 1 is the Dirichlet reservoir that carries no element equation.  The
stage rows inherit that: cell 1 and the two lower ghosts stay with the sweep,
exactly as the proton does today (`ionization_equilibrium.f90`, the
`j .ge. 1` guard, READ), so no handoff has to state an ionization fraction.
The face flux at f = 0 is then purely advective, `N_el(0) = F_rho(0) Y_el(0)/m_el`,
and the stage fractions it carries are those of the cell the gas comes from:
with `F_rho(0) > 0` (inflow from the reservoir) the reservoir composition,
which for a scalar handoff is the base cell's own equilibrium partition and
for a `base.inp` or profile handoff is the handed-over one; with
`F_rho(0) < 0` (the breathing base, which the code admits and `docs/` records
as normal) the domain's own, taken upwind from cell 1.  The rule is the sign
of the face flux and nothing else, which is the rule `species_face_fraction`
already applies.

**The top face, f = N.**  The diffusive flux is zero there by the same face
coefficients, and the outer ghosts carry the outermost physical cell's
composition, so the advective reconstruction sees no compositional step across
the top (`element_transport_residual`, THE BOUNDARIES ARE THE CALLER'S, READ).
For the stage rows this is the outgoing-characteristic statement: with
`N_el(N) >= 0` the boundary is outflow and needs no data, the column's own
composition leaving; with `N_el(N) < 0` the zero-gradient ghost brings the
same composition back in, so no stage is created or destroyed at the top.  The
stage term `E_k` vanishes at that face with the ghost equal to cell N.

**Zero element flux.**  `N_el(f) = 0` does not make the stage flux zero: the
stage term `-n_el K dx_k/dr` survives, and the identity holds trivially
(`0 = 0`).  Only with `K(f) = 0` as well does the face carry nothing.  The
test runs both.

**Reversed element flux.**  The donor side flips for every stage together,
because it is taken from the sign of `N_el` and not from any stage's own
gradient.  A stage-by-stage donor rule is one of the RED rows.

**The matching to the inner equilibrium approximation.**  Below the binding
radius the composition is the local root and above it the transported one.
In a finite-volume method the face flux is single valued and shared by the two
cells that meet at it, so **the flux is continuous across the gate by
construction**; what the matching supplies is the inflowing stage FRACTION,
taken from the equilibrium side.  The physical condition for that to be a
matching and not a jump is that the equilibrium side really carries the local
root, that is `Da >> 1` there.  The L12a Damkohler table starts at 1.99 R_p
and gives `Da_H(r/v) = 1.81`, `Da_He = 2.11` there (READ), which is O(1) and
not a gate; inside that radius the table has no entry and the gate radius is
therefore an open measurement stage C owes, taken on the same construction the
table used.  Until it exists the honest default is the one the certification
already uses, `cert_regime_wind_r = 1.20` (READ), with the cells below it
reported and not gating.

---

## 6. The collisional regime of the neglected stage drift

Setting the drift between charge stages of one element to zero says that all
stages of an element move at one velocity.  What decides it is the drift
velocity an ion would acquire against its own neutral under the ambipolar
field, against the velocity at which the flow carries them both.

Within one element the masses are equal to within the electron mass, so the
settling coefficient of section 1 reduces to the field term alone,

```
   G_stage = -(Z_ion - Z_neutral) eE/(kT) = dln(n_e T)/dr ,
   w       = D_in G_stage ,
```

using the operator's own ambipolar field `eE = -kT dln(n_e T)/dr` and its own
`ion_neutral_pair_diffusion` for `D_in`.  MEASURED on the certified 45 R_p
state `LHS1140b/models/.L8/r45/output/` (the state L12a measured, `He_Kzz` =
1e9, `He_diffusion: True`):

| r [R_p] | T [K] | v [cm/s] | x(H II) | D_in(H+,H) [cm^2/s] | \|w_H\| [cm/s] | \|w_H/v\| | \|w_He/v\| |
|---|---|---|---|---|---|---|---|
| 1.05 | 1578 | 5.72e+00 | 4.69e-04 | 2.70e+08 | 1.79e+00 | 3.14e-01 | 2.01e-01 |
| 1.20 | 4466 | 2.10e+02 | 8.93e-03 | 2.43e+10 | 4.18e-01 | 1.99e-03 | 1.19e-03 |
| 2.00 | 4111 | 3.13e+03 | 3.69e-02 | 7.61e+11 | 1.28e+03 | 4.09e-01 | 2.46e-01 |
| 4.02 | 2038 | 1.59e+04 | 1.14e-01 | 8.81e+12 | 6.65e+03 | 4.17e-01 | 2.63e-01 |
| 9.98 | 831 | 6.46e+04 | 3.12e-01 | 1.11e+14 | 3.30e+04 | 5.11e-01 | 3.44e-01 |
| 20.17 | 447 | 1.27e+05 | 5.22e-01 | 5.40e+14 | 7.11e+04 | 5.61e-01 | 3.95e-01 |
| 44.22 | 252 | 2.03e+05 | 7.92e-01 | 2.58e+15 | 9.69e+04 | 4.76e-01 | 3.49e-01 |

Largest over `r >= 1.2 R_p`: `|w_H/v| = 0.72`, `|w_He/v| = 0.53` (MEASURED).
Over 1.0 to 1.2 R_p the ratio reaches 2.50 and 1.77, at velocities of a few
cm/s where the composition is the local root anyway.

Two corrections, both toward a smaller drift, neither applied:

1. Resonant charge exchange is not a channel of `ion_neutral_pair_diffusion`
   (the module header states it is excluded because both partners carry the
   same element).  Frictions add, so adding it lowers `D_in` and lowers `w`.
   READ, `docs/collisional_validity.md` section 1: resonant H+ + H at 1e4 K is
   about 2e-15 cm^2 against the about 1e-15 cm^2 rigid core, so the friction
   roughly triples and `w` falls to roughly a third.  That leaves
   `|w_H/v|` of order 0.15 to 0.25 in the wind, which is still not a
   negligible term.
2. The trace-ion picture behind `w = D_in G` fails where the gas is mostly
   ionized, and `x(H II)` reaches 0.79 at 44 R_p (MEASURED): there the ion is
   the bulk and the neutral is the trace that drifts, with the same relative
   velocity but the opposite attribution.

**Where the approximation is safe, measured.**  The gas stops being
collisional for the neutrals long before it does for the ions.  MEASURED with
`src/utils/collisional_validity.py` on the same state:

| r [R_p] | Kn_bulk | Kn(H I) | Kn(H II) | Kn(He I) | Kn(He II) |
|---|---|---|---|---|---|
| 1.20 | 4.74e-04 | 6.30e-04 | 2.65e-06 | 4.38e-04 | 1.34e-06 |
| 4.02 | 1.89e-02 | 2.78e-02 | 1.24e-06 | 1.94e-02 | 6.26e-07 |
| 9.98 | 8.06e-02 | 1.73e-01 | 6.37e-07 | 1.36e-01 | 3.21e-07 |
| 20.17 | 1.21e-01 | 4.39e-01 | 4.07e-07 | 3.82e-01 | 2.05e-07 |
| 40.42 | 1.06e-01 | 8.85e-01 | 2.89e-07 | 8.35e-01 | 1.46e-07 |

`Kn(H I)` crosses 0.1 at 7.25 R_p and `Kn(He I)` at 8.42 R_p (MEASURED); the
ions stay Coulomb-collisional to 1e-7 everywhere.  The tool's own verdict on
this state is HYDRODYNAMIC RESULT UNVALIDATED, max Kn 0.88 in the
heating and acceleration region (MEASURED, `Kn_bulk = 0.1` at 12.91 R_p,
critical point 40.06 R_p).

**Statement to be written at the site in stage B.**  All stages of an element
are given one velocity.  This is exact where the ion-neutral momentum transfer
is fast against the flow, which on this wind means inside about 7 R_p for
hydrogen and about 8 R_p for helium by the neutral Knudsen number, and it is
immaterial inside about 1.2 R_p where the Damkohler number returns the local
root whatever the transport does.  Between those radii and above them the
neglected drift is a systematic of order 0.2 to 0.5 of the advective stage
flux, measured here, and it is not smaller than the ionization transport the
option exists for.  Above about 13 R_p the continuum stage equation is itself
unvalidated by the run's own Knudsen measure, and the drift is not the leading
approximation there.  **Stage B carries this as a stated approximation with
its measured size, and an ambipolar stage drift is a separate item, not a
refinement folded into this one.**

---

## 7. What stage B needs that the operator does not expose

All four are in `src/modules/functions/binary_element_diffusion.f90` and all
four are private today (READ, the `public ::` list of lines 313 to 429):

1. **`drift_and_gradient_face_coefficients`**, or an accessor returning
   `Agrd(0:N)`, `Bdrf(0:N)`, `updrf(0:N)` for a given state.  Without it no
   caller can form the diffusive half of the element face flux, and a second
   copy of it inside the ionization operator is the "two spellings of one
   term" the module header records the cost of.
2. **`element_face_flux`**, or the assembled `J(f)` array.  The drift product
   takes its two factors from opposite sides of the face, so `J(f)` cannot be
   rebuilt from `Agrd` and `Bdrf` without also rebuilding that rule.
3. **The element NUCLEUS face flux itself**, which exists nowhere: the
   operator forms two divergences and never a single face flux, and the two
   carry different geometric factors (section 1).  The clean addition is one
   routine that returns, for a given state and face mass flux, the four
   arrays the stage rows need: the advective face element mass flux
   `F_rho(f) Y_el(f)` for each element, the diffusive face flux `J(f)`, the
   face element nucleus density `n_el(f)` and the face mass per nucleus
   `m_1(f)` of component 1.  It is the same object
   `element_advective_divergence` and `composition_residual` already build
   internally, exposed rather than duplicated.
4. **`settling_coefficient`'s `Gco`**, if the stage rows are ever to carry a
   drift of their own (section 6 says they should not yet).  Only `dmeff` is
   exposed today, through `relative_settling_mass`.

Until 1 to 3 exist, stage B cannot be written without a second copy of the
element flux, and the identity of section 3 would then be a statement about
two operators instead of one.  **This is the blocking item for stage B and it
is a change to `binary_element_diffusion.f90`, which this item does not
own.**

---

## 8. Stage F: the line from the certified composition, and what is retired

**What the line reads today.**  `EXHALE_transit.py` reads
`output/Hydro_ioniz_adv.txt` and `output/Ion_species_adv.txt` (lines 79 and
80, READ) and takes the He 2^3S density from column 7 of the ion file
(`nheiTR`, line 323, READ), with the temperature and velocity from the hydro
file.  So every He I 10830 equivalent width of the atomic catalog is quoted
on the post-processed composition and not on the certified state, which is
the factor 13.3 of L10 (15.5116 against 1.1702 percent A, READ).

**What stage F changes.**  Once `x(H II)`, `x(He II)` and `x(He III)` are
carried and certified, the composition of the SOLUTION is the transported
one, and the line is synthesized from it:

- the transit tool's two input files become `output/Hydro_ioniz.txt` and
  `output/Ion_species.txt`, selectable rather than hard-coded, with the file
  pair named in the products' headers so a curve always says which state it
  stands on;
- the He 2^3S density stays column 7 of the ion file, that is, the sweep's
  own local solution of row 4 on the transported and certified
  `(n(H I), n(He I), n_e, T)`.  The level is not carried, by the Damkohler
  argument of L12a section 3, and the 5 percent promotion rule of that section
  is evaluated at stage D and reported;
- the line profile itself is unchanged: impact-parameter Voigt integration,
  disk average, instrument and rotation convolution.  Only the composition it
  integrates moves.

**What is retired, and what is not.**  The `_adv` advection correction is
retired AS THE SOURCE OF THE CARRIED FRACTIONS, that is, for `x(H I)`,
`x(H II)`, `x(He I)`, `x(He II)` and `x(He III)`: those are then solved and
certified in the run, and a second, differently discretized answer for them in
the file the line reads would be the state-mixing L10 objects to.  What is NOT
retired is `post_process_adv` itself:

- it stays the independent discretization the acceptance test of L12a section
  4.2 is measured against (first-order upwind marching at the bulk velocity,
  no eddy term, no element drift, against the simultaneous solve), and its
  measured difference from the solution is written to `pp.log` and to the file
  header so the identity is stated by the run;
- it keeps every quantity it alone produces, the advection-corrected
  temperature and internal-energy balance among them, for the metal stages
  while those are not carried, and for any run with the key off, which is the
  default.

---

## 9. Noticed outside this item's scope, reported and not changed

- **`docs/input_schema.md` K15f is stale on two counts**, both already
  recorded in L12a section 7 and both still present (READ, line 147): it says
  `Ionization transport: True` with `Solver: Newton` also needs
  `Coupled carrier solve`, and that "Each of the three is a fatal `error stop`
  in `input_read`".  That refusal was removed on 2026-09-12.  The same stale
  sentence is repeated in `Update_EXHALE_stage1.md` section 169.1.  Not
  edited: that file belongs to the input-schema item and stage B rewrites K15f
  anyway.
- **The two discrete geometries of the element row** (section 1, last
  bullet) are an `O(dr^2)` inconsistency inside one operator: the diffusive
  divergence uses `1/(rp^2 (rep_+ - rep_-))` and the advective one
  `(rep_+^3 - rep_-^3)/3`.  It is not a defect of the identity derived here,
  which is a face statement, but it means "the element face flux" is not a
  single well-defined number in the present code, and section 7 item 3 is
  written around that.  Reported, not changed: it moves every element row of
  every state.
