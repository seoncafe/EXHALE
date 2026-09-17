# LHS 1140 b, item L12a: design of ionization transport on the atomic partitioned stationary route

Item L12a of `docs/PLAN_20260913_lhs_stationary.md`, opened by item L10
(`docs/lhs1140b_stationary_L10_20260913.md`), which measured that the
certified LHS 1140 b wind is solved with a photoionization-equilibrium
composition its own post-process contradicts by a factor 4.5 in the neutral
hydrogen density and 6.0 in the He 2^3S column, and that the He I 10830
equivalent width differs by a factor 13 between the two compositions. This
is a DESIGN memo. No file under `src/` is changed by it.

Every number is MEASURED (computed here from the files named) unless marked
READ. The state measured throughout is the certified 45 R_p solution,
`LHS1140b/models/.L8/r45/output/`, with `R_0 = 1.1273716464e9 cm`.

---

## 0. Verdict, before the detail

1. **Three fractions must be carried: x(H II) per H nucleus, x(He II) and
   x(He III) per He nucleus.** Their Damkohler numbers on this wind fall
   below 0.3 above 4 R_p and reach 1.1e-2 (H) and 6e-2 (He) at the outer
   boundary, so the local root is not the composition of that gas.
2. **The He 2^3S level must NOT be carried.** Its destruction rate is 15.7
   to 2.2e4 times the flow rate everywhere in the domain, on the gradient
   scale of its own equilibrium fraction; it is a local function of the
   transported (n(H I), n(He I), n_e, T) to within 1/Da, at worst 6 percent
   at 44 R_p and 5e-5 inside 4 R_p.
3. **The variable has to be the fraction per element nucleus, not the mass
   fraction per unit mass that the present proton carrier uses.** With
   element diffusion on, which is the LHS 1140 b configuration, the helium
   element separates from hydrogen by a factor 4.1 between the base and
   2 R_p (n(He)/n(H) 1.602 -> 0.499 -> 0.387), so an ion stage advected as a
   mass fraction on the bulk mass flux and an element advected on its own
   element flux do not compose into a consistent partition. Per element
   nucleus, the ionization operator moves charge only and the element
   operator moves nuclei only, and element conservation is exact by
   construction rather than by cancellation.
4. **The refusal in `input_read.f90` is a statement about the OPERATOR, not
   about the physics**: `Ionization transport: True` requires
   `Molecular chemistry` and `Molecular carrier transport` only because the
   proton was made the fifth carrier of the molecular transport and
   chemistry solve (`ic_Hp`, `diffusive_photochemistry.f90`), which does not
   run in an atomic gas. Nothing in the proton's own equation needs a
   molecule.
5. **The acceptance test is that the post-process becomes an identity**, and
   its tolerance cannot be the 1e-5 of the certification: the post-process
   solves a first-order upwind recursion at the bulk velocity with no eddy
   term and no element drift, so the two equations differ by terms this memo
   bounds and the implementation must measure. The proposed gate is
   1e-2 on the composition in the wind, against the factors 4.5 and 6.0 it
   stands at now.

---

## 1. What is carried, and the equation each carried quantity satisfies

### 1.1 The measurement that decides the set

For a stage fraction x of an element whose nucleus flux is conserved, the
stationary balance is

```
   n_el ( v dx/dr )  =  P  -  L x n_el ,
```

with no geometric divergence term: the divergence of the nucleus flux is
what cancels it. The local root x_eq = P/(n_el L) is the composition of the
gas only where the relaxation rate L is fast against the rate at which the
flow carries the parcel through a change of x_eq, so the number that decides
is

```
   Da = L * l / |v| ,      l = |x_eq / (dx_eq/dr)| ,
```

and `l = r` is the coarse version of it. Both are given below. `L` is the
full loss rate of the stage as the code writes it: for hydrogen
`P(H I) + alpha_B n_e`, for helium `P(He I) + alpha(He II) n_e`, with the
photoionization rate read back from the run's own equilibrium
(`P = n_ion n_e alpha / n_neutral`, the construction L10 used, so the rate is
the run's and not a second estimate of it) and alpha a Case B power law used
only to set a time scale; for the He 2^3S level

```
   L(2^3S) = A31 + Q31 n(H I) + (q31a + q31b) n_e + k_ci(2^3S) n_e
```

with `A31 = 1.272e-4 s^-1` and the four coefficients as
`util_ion_eq.f90::HeITR_coeffs` and `Cool_coeff.f90` write them (READ:
`ioniz_HeI23S_H`, `coex_HeI_23S_21S`, `coex_HeI_23S_21P`, `ci_HeI23S`).
Photoionization of the metastable is not in this sum, so the level's
Damkohler numbers below are LOWER bounds.

| r [R_p] | t_flow = r/v [s] | Da_H (r/v) | Da_H (l/v) | Da_He (r/v) | Da_He (l/v) | Da_2^3S (r/v) | Da_2^3S (l/v) | n(He)/n(H) |
|---|---|---|---|---|---|---|---|---|
| 1.99 | 7.18e+05 | 1.81e+00 | 3.10e+01 | 2.11e+00 | 8.33e-01 | 1.10e+05 | 2.16e+04 | 0.4990 |
| 4.02 | 2.84e+05 | 2.14e-01 | 1.08e+00 | 2.93e-01 | 2.10e-01 | 1.74e+03 | 4.07e+02 | 0.3914 |
| 5.97 | 2.11e+05 | 8.95e-02 | 3.14e-01 | 1.40e-01 | 1.25e-01 | 2.45e+02 | 7.64e+01 | 0.3879 |
| 8.06 | 1.84e+05 | 5.05e-02 | 1.66e-01 | 9.25e-02 | 1.04e-01 | 8.00e+01 | 4.21e+01 | 0.3873 |
| 9.98 | 1.74e+05 | 3.51e-02 | 1.06e-01 | 7.50e-02 | 1.08e-01 | 4.61e+01 | 5.12e+01 | 0.3872 |
| 15.00 | 1.71e+05 | 2.05e-02 | 3.78e-02 | 6.06e-02 | 1.85e-01 | 2.70e+01 | 4.87e+01 | 0.3871 |
| 20.17 | 1.79e+05 | 1.61e-02 | 2.09e-02 | 5.85e-02 | 5.13e-01 | 2.46e+01 | 2.13e+01 | 0.3871 |
| 24.93 | 1.90e+05 | 1.46e-02 | 1.56e-02 | 5.97e-02 | 2.24e+00 | 2.50e+01 | 1.74e+01 | 0.3870 |
| 28.74 | 2.00e+05 | 1.41e-02 | 1.35e-02 | 6.16e-02 | 5.02e-01 | 2.59e+01 | 1.62e+01 | 0.3870 |
| 44.22 | 2.45e+05 | 1.43e-02 | 1.13e-02 | 7.30e-02 | 1.75e-01 | 3.13e+01 | 1.57e+01 | 0.3870 |

The share of the 2^3S loss rate carried by the 1.272e-4 s^-1 radiative decay
alone is 0.001 at 2 R_p, 0.110 at 6, 0.481 at 10, 0.926 at 20 and 0.996 at
44 R_p; below 10 R_p the Penning channel Q31 n(H I) dominates. So the level's
Damkohler number never approaches unity from either end, and it is the decay
constant that holds it up in the outer wind where everything else has thinned
out.

**Set carried: H II, He II, He III. Not carried: the He 2^3S level, every
metal stage (the metals are absent from this configuration and their
Damkohler numbers are not measured here), and the molecular ions.**

The helium column of the table also states why He III cannot be left out:
x(He III) rises from 1.2e-4 at 2 R_p to 0.301 at 44 R_p on the solution, so
helium is a three-stage problem in the outer wind, not a two-stage one.

### 1.2 The variable, and how element conservation is kept

`System_HeH_TR.f90` already solves the local problem in exactly the
variables proposed here: `x(1) = n(H II)/n_H`, `x(2) = n(He II)/n_He`,
`x(3) = n(He III)/n_He`, `x(4) = n(He 2^3S)/n_He` (READ). The transported set
is then rows 1 to 3 of that system, and the sweep keeps row 4 and the
metals.

The present proton carrier does not use that variable. It carries
`f_sp(:,isp_HII)`, a fraction per unit mass, through the same face mass flux
as H2 (`carrier_state`, `species_face_flux`; stage 1 section 169 states the
choice, READ). With a uniform elemental composition the two are proportional
by a constant, because the H nuclei per unit mass are then constant, and the
choice is free. **With `He_diffusion: True` they are not**: n(He)/n(H) runs
1.602 at the base, 1.483 at 1.05, 0.988 at 1.20, 0.499 at 2.0 and 0.387
above 4 R_p (MEASURED, last column of the table), so the hydrogen mass
fraction varies by a factor 2.3 over the inner two radii. Advecting
Y(H II) on the bulk mass flux while the element operator advects Y(H) on the
mass flux PLUS its own diffusive flux gives a ratio x = Y(H II)/Y(H) that
satisfies neither equation. The defect is in the region 1 to 2 R_p, which is
where the He 2^3S population peaks.

**Judgement: carry the fraction per element nucleus, and give its advective
half the ELEMENT nucleus flux, not the bulk mass flux.** The element flux is
the object `binary_element_diffusion` already forms and
`relax_element_composition` already balances; the ionization rows read it the
way the element relaxation reads the face mass flux from its caller
(`face_mass_flux_of_state`). Then

- the element operator moves nuclei and does not touch the partition;
- the ionization operator moves the partition and cannot move a nucleus;
- `n(He II) = x(He II) n(He)` is unambiguous, and the He element total of a
  cell is untouched by the ionization step, exactly as the carrier operator's
  header claims for its own species.

### 1.3 The transport equation and where each term comes from

**CORRECTED 2026-09-16 by `docs/lhs1140b_stationary_L12b_derivation_20260916.md`:
the `Phi_x` written below double counts the eddy flux, because the element
flux `x` rides on already carries it; the additional stage term is
`-n_el K dx/dr` alone.  Read that memo for the corrected flux, its discrete
identity, the boundary and matching conditions and the collisional regime.**

For each carried fraction, in a cell j,

```
   d(n_el x)/dt + (1/r^2) d/dr [ r^2 ( n_el x u_el  +  Phi_x ) ]  =  S(x; local state)
   Phi_x = - n_tot K_zz d(x n_el/n_tot)/dr
```

- `u_el` is the element's own velocity, that is, the element nucleus face
  flux the element operator forms (bulk advection plus the element diffusive
  and settling flux). Where `He_diffusion` is off it is the bulk mass flux
  and nothing changes.
- **No stage-resolved diffusive flux.** An ion drifting relative to its own
  neutral is an ambipolar problem and needs the (Z e E)/(kT) term the
  molecular carrier header states is absent for its four neutral species;
  the element operator builds a `settling_coefficient` with that field
  already, but it applies it to an ELEMENT and not to a stage. Setting the
  stage-resolved drift to zero is the same statement the present proton
  carrier makes by setting its molecular diffusion to zero (READ,
  `docs/input_schema.md` K15f), written once instead of twice.
- **The eddy term is kept** because `K_zz` is bulk parcel exchange and does
  not distinguish charge. Where it matters at all, below the homopause,
  Da is 1e4 and larger and the row returns the local root regardless, so
  keeping it is a consistency statement rather than a measurable one.
- `S` is taken WHOLE from the rows the local solve already uses, at the trial
  fractions, with everything else frozen: rows 1 to 3 of `heh_tr_rows`
  (`ion_residual_core.f90`) in an atomic gas and rows 1 to 3 of
  `mol_heh_rows` (`System_HeH_mol.f90`) in a molecular one. That is the rule
  the carrier operator already follows for H2 ("P_i - L_i is not written
  here", module header, READ), and it is what keeps photoionization,
  secondary ionization, radiative recombination, collisional ionization, the
  He/H charge exchange, the Penning source of protons and, in a molecular
  gas, the R10 and R13 channels from existing in two transcriptions.
- The electron density inside the rows is the frozen background with the
  moving stages replaced by their trial values, `n_e = n_e,bg - sum_k Z_k
  n_k,bg + sum_k Z_k n_k,trial`. The proton case of this already exists
  (`carrier_source`, `cbg_ne - cbg_nhii + n_hii`, READ) and the generalization
  is the same sum over the three charges.
- The He 2^3S population that appears in those rows (the Penning term of the
  H row, the q13/q31 terms of the He rows) is frozen background inside the
  transport step and refreshed by the sweep between steps, which is what the
  Damkohler numbers of section 1.1 license.

---

## 2. Where this sits in the partitioned route

### 2.1 The route, and why nothing has to be reordered

The LHS 1140 b configuration takes `steady_wind_with_element_diffusion`
(`EXHALE_main.f90`) with `species_alternated = he_diffusion` true and
`carrier_in_newton` false. One pass is, in order (READ from the routine):

```
 solve_steady_jfnk / solve_steady_ptc   (wind at fixed composition)
 U_to_W, get_species_densities, comp_T_from_p
 ioniz_eq                               (the sweep)
 assemble_residual, certification_evaluate   (the joint test on THIS state)
 progress control
 relax_element_composition  + ioniz_eq       (element rows, fixed wind)
 relax_photochemical_composition             (carrier rows, fixed wind)
```

The ionization rows enter as a third relaxation **after** the element
relaxation, because the flux they ride on is the element flux that
relaxation just set, and before the pass line. In a molecular run they are
rows of the same solve as the carriers (section 2.2), so the order is
unchanged there too. `species_alternated` gains `.or. ionization_transport`,
which is what makes an atomic run with the key on take the alternation
instead of the single solve it takes today.

### 2.2 One operator or two

The proton row lives today inside `diffusive_photochemistry`, whose Newton
block is dense within a cell precisely so that H2 and H+ move together: on
the hot Uranus the dominant proton sink is charge exchange with H2 (R10, 70
times radiative recombination in the outer wind, READ from stage 1 section
168), so splitting them into two alternated operators would replace an
implicit coupling by a Picard one at the H2 front.

**Judgement: one operator, with the transported set resolved from the
configuration.** The set becomes

```
   {H2, OH, H2O, CO}         where the run carries them   (mass fraction, own diffusion)
   {H II, He II, He III}     where the key is on          (nucleus fraction, element flux)
```

and the `thereis_mol` and `carrier_transport` gates on the ionization part
are lifted, so the set may consist of the ionization rows alone. That keeps
one Newton block per cell, one residual, one row-scale table, one headroom
rule and one certification shape; two operators would give the same physics
two spellings, which is the situation the code has already paid for once
(`ATES/ATES_extended`, the decoupled metal solver).

What has to be generalized inside it, concretely:

- `carrier_set_init`: indices `ic_Hp`, and new `ic_Hep`, `ic_Hepp`; the set
  no longer implies H2.
- `carrier_state`: the ionization members are filled as
  `n(stage)/n(element)` rather than from `f_sp` directly.
- `carrier_element_reference_density`: He II and He III take `nHe_free`, as
  H II takes `nH_free`. The module already records what happens when a
  carrier is charged to the wrong element (the proton on `nO_free`: the
  diagonal derivative comes back exactly zero,
  `src/tests/carrier_reference_scales/`, READ).
- `carrier_face_flux` / `carrier_face_coefficients`: for an ionization row
  the molecular `Agr` is zero, `Bst` is zero, and the advective half comes
  from the element flux.
- `hydrogen_available_to_carriers` and the element headroom: an ionization
  fraction is bounded by 1 and by the sum rule x(He I) + x(He II) +
  x(He III) = 1, not by a free-nucleus density. The headroom of an
  ionization row is therefore the simplex face, and it needs stating
  separately rather than reusing the nucleus budget.
- The base: no handoff states an ionization fraction, so the two lower
  ghosts stay with the sweep, exactly as the proton does now
  (`ionization_equilibrium.f90`, the `j .ge. 1` guard, READ).

### 2.3 How the sweep receives the transported composition

`ioniz_eq` already has the mechanism, for one row. `ieq_cell%x_hp_fixed`
replaces row 1 of the molecular system by `x(1) - x_hp_fix`
(`System_HeH_mol.f90:585`, `System_HeH_mol_metals.f90:290`), and
`ionization_equilibrium.f90` sets it inside a block gated on
`thereis_mol .and. carrier_transport .and. (bg_ready .or. do_load_IC)`.

The change is: gate the ionization part of that block on
`ionization_transport .and. (bg_ready .or. do_load_IC)` alone, and set
three fixed fractions instead of one. The substitution itself must be
written **once**, in `ion_residual_core.f90`, and called by every
`System_*` module that carries these rows, replacing the two ad-hoc lines
that exist today. Seven systems can reach these rows
(`System_H`, `System_HeH`, `System_HeH_TR`, `System_HeH_metals`,
`System_HeH_TR_metals`, `System_HeH_mol`, `System_HeH_mol_metals`), and
seven transcriptions of one substitution is the drift the `g2s`/`G2s`
incident of 2026-08-12 is the standing warning about.

The sweep then solves, at the imposed H and He partition: the He 2^3S level
(row 4), the metal stages, the molecular partition where there is one, the
electron density, the temperature and the heating and cooling. **The heating
and cooling of the energy equation are therefore evaluated at the transported
composition**, which is the point of the whole item: photoheating is
proportional to the neutral density, and L10 measured the equilibrium
composition carrying up to 5.9 times too little of it at 44 R_p.

### 2.4 The certification rows

One entry per carried ionization fraction, on the shape
`carrier_row_entry` already has: `regime_gated = .true.`, measure
`|res(j)| / terms(j)` with `terms` the sum of the row's own terms, verdict
taken over `r >= cert_regime_wind_r = 1.20`, the rest of the column
reported.

**The 1e-5 of `cert_tol_carrier_at` may not simply be reused.** Its five
anchors (READ, `certification.f90` lines 326 to 360 and
`docs/certification_tolerance_anchoring_20260910.md`) are all measurements of
a TRANSPORT row, whose balance is a cancellation between an advective and a
diffusive divergence: the element operator's own two floors, the 1.52e-6
discretization error of a smooth manufactured column at the production
spacing, and the 3.05e-4 and 7.19e-2 the two candidate states carried. An
ionization row is a cancellation between an advective divergence and a
reaction rate, and none of those five numbers is about it.

So the implementation owes one anchoring measurement before the tolerance is
fixed: a smooth manufactured ionization column (a chosen `x_eq(r)` and the
`P`, `L` that produce it), refined by a factor two at the production
spacing, giving the discretization error of the row itself. The gate is then
a decade above that, as the carrier gate is a decade above its own anchor.
Until that number exists the entry is REPORTED and not gating, which is what
the report already does for the cells below 1.20 R_p, and the memo records
the expectation that it will land near 1e-5 rather than asserting it. No
tolerance is to be chosen so that a state passes.

---

## 3. The He 2^3S level

**Judgement: the level is solved locally, by the sweep, on the transported
and certified (n(H I), n(He I), n_e, T). It is not carried.**

The measurement is section 1.1: `Da(2^3S)` on the gradient scale of its own
equilibrium fraction is 2.2e4 at 2 R_p, 4.2e1 at 8, 4.9e1 at 15 and 1.57e1
at 44 R_p, and the smallest value anywhere in the physical domain is 15.7.
A linear level balance `dn/dt = P - L n` relaxes on 1/L, so the stationary
solution departs from the local root by O(1/Da): at worst 6 percent at the
outer boundary, 2 percent at 20 R_p, 2.4e-2 at 8 R_p and 5e-5 inside
4 R_p. Repeating the calculation on the `_adv` composition instead of the
solution's (more neutral hydrogen, so a faster Penning channel; fewer
electrons) gives the same verdict with the smallest value 2.7e1 on `r/v`.

What the level is NOT insensitive to is the composition that feeds it. On
this state `n(He 2^3S)` differs between the two closures by a factor 33 at
10 R_p, 79 at 28.7 and 93 at 44.2 (MEASURED, `_adv` over solution 3.05e-2,
1.27e-2, 1.07e-2), and its radial column over 1 to 30 R_p by 6.045
(4.073e11 against 6.738e10 cm^-2, MEASURED). That factor follows the
electron density, which the recombination source of the metastable is
proportional to: `n_e(_adv)/n_e(solution)` is 0.195, 0.127 and 0.117 at the
same radii (MEASURED). **The line is a statement about the ionization
closure, and the level equation on top of it is not where the error is.**

The cost of the alternative is the reason to state the judgement rather
than take the safe option. Carrying the level would make the post-process an
exact identity in all four unknowns instead of three, because
`System_implicit_adv_HeH_TR` advects all four. It would also add a row whose
loss rate is 2.2e4 times the flow rate at 2 R_p, that is, the stiffest row
in the system by four decades, whose unknown is 1e-7 to 3e-4 of its element,
to a Newton block whose column and row scales had to be separated once
already for exactly this reason (`carrier_column_scale` against
`carrier_row_term_scale`, READ). It buys a correction of at most 6 percent
of the level, in one place, against a factor 13 in the equivalent width that
the three ionization rows buy.

**The decision rule is fixed here, before the measurement**: if the
implemented state shows the post-process moving the He 2^3S column by more
than 5 percent when the three ionization fractions already agree to their
own gate, the level joins the carried set. Otherwise it does not, and the
measured difference is reported as the systematic it is.

---

## 4. Acceptance

### 4.1 The reproduction (RED), on the certified 45 R_p state

| quantity | solution | `_adv` | ratio |
|---|---|---|---|
| n(H I) at 10.0 R_p | | | 1.382 |
| n(H I) at 28.7 R_p | | | 2.720 |
| n(H I) at 44.2 R_p | | | 4.512 |
| n_e at 10.0 / 28.7 / 44.2 R_p | | | 0.195 / 0.127 / 0.117 |
| n(He 2^3S) radial column, 1 to 30 R_p [cm^-2] | 4.073e11 | 6.738e10 | 6.045 |
| photoheating at 44.2 R_p | | | 5.89 (READ, L10 section 2) |
| He I 10830 equivalent width [percent A] | 15.5116 | 1.1702 | 13.3 (READ, L10 section 5) |

(All ratios MEASURED here from the two file pairs except the two marked
READ. The item brief quotes "n(H I) 5.9x"; 5.9 is the photoheating ratio,
and the neutral hydrogen ratio is 4.51, both at 44.2 R_p.)

The measured equivalent width for comparison is 1.108 +/- 0.030 percent A
(Cherubim et al. 2026, READ from `LHS1140b/MODELS.md` through L10 section 5).

### 4.2 The test (GREEN)

Run the post-process on the certified state produced with the ionization
transport on, and compare its `_adv` composition with the solution's own, on
the cells where `adv_correction_valid` is true:

```
   max_j  | x_adv(j) - x_solve(j) | / max( x_solve(j), 1e-8 )     for x(H I), x(He I), x(He III)
```

**Gate: 1e-2 over the cells with r >= 1.20 R_p.** It is not 1e-5, and the
reason is that the two equations are not the same equation. What separates
them, each to be measured separately during the implementation rather than
assumed:

1. **Discretization.** The post-process is a first-order upwind marching
   recursion cell by cell (`As = dr/(v(j-1) v0)`, `x_old` from cell j-1,
   `post_process_adv.f90`, READ); the transported rows are a simultaneous
   solve of the same balance on the same grid with the flux reconstruction
   the hydrodynamics uses. First order against second order on a smooth
   profile at this spacing is the leading term, and it is bounded by
   refining the grid.
2. **The eddy flux** the post-process has no term for. Above 1.2 R_p the
   ionization rows are advection-dominated by the Damkohler argument of
   section 1.1 in reverse (the eddy time is `dr^2/K_zz`); the measurement is
   to re-run the comparison with `He_Kzz` set to zero and read the change.
3. **The element drift** the post-process has no term for. It is bounded by
   the measured helium separation: n(He)/n(H) changes by 1.1 percent between
   4 and 44 R_p and by 3.4e-4 between 10 and 44 R_p (MEASURED), so above
   10 R_p this term is below 1e-3 and it is only inside 4 R_p that it can
   matter, where Da is of order unity and the two closures nearly agree
   anyway.
4. **The He 2^3S level**, not carried, compared separately against the
   5 percent rule of section 3.

Beside the composition gate, two physical readings that the item exists for:

- the He 2^3S radial column computed from the solution and from `_adv` agree
  to the same 1e-2, against 6.045 now;
- `EXHALE_transit.py` run on the two pairs gives the same equivalent width to
  within the composition gate, against a factor 13.3 now.

### 4.3 What this changes that is not a test

Adopting the transported closure changes every composition-dependent number
of the LHS 1140 b campaign, and the critical point of this wind sits at
40.06 R_p (READ, L8), that is, INSIDE the region where the closure is wrong,
so the mass-loss rate itself is at stake and not only the line. The
certified 30, 45 and 60 R_p solutions become states of a superseded closure,
and the He/H that the campaign crosses to meet the measured equivalent width
moves with them. This is to be reported as a consequence, not guarded
against.

---

## 5. Implementation plan

### 5.1 Files, with what each change is

| file | change |
|---|---|
| `src/modules/nonlinear_system_solver/ion_residual_core.f90` | ONE routine that imposes the transported fractions on rows 1 to 3 of an ionization residual, replacing the two ad-hoc `x_hp_fixed` lines; the only place the substitution is written |
| `src/modules/nonlinear_system_solver/ion_cell_state.f90` | the imposed He II and He III fractions beside the existing `x_hp_fix`, and one flag for the set |
| `System_HeH.f90`, `System_HeH_TR.f90`, `System_HeH_metals.f90`, `System_HeH_TR_metals.f90`, `System_HeH_mol.f90`, `System_HeH_mol_metals.f90`, `System_H.f90` | one call each to the routine above, in place of nothing (five of them) or of the existing line (two) |
| `src/modules/radiation/ionization_equilibrium.f90` | the imposed-fraction block gated on `ionization_transport` instead of on `thereis_mol .and. carrier_transport`; three fractions set instead of one |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | the transported set generalized: `ic_Hep`, `ic_Hepp`; ionization members held as fractions per element nucleus; their reference density, headroom (the simplex face), face coefficients (no molecular diffusion, no settling) and advective flux (the element flux); the `thereis_mol` and `carrier_transport` early returns lifted for a set that is ionization only; `carrier_source` extended to rows 2 and 3 with the trial electron sum; the module header rewritten, since it is no longer a molecular statement |
| `src/modules/functions/binary_element_diffusion.f90` | `advected_carrier_register` gains the element a carried species belongs to, so an ionization row is advected on its element's nucleus flux and a molecular carrier on the bulk mass flux, as now |
| `src/EXHALE_main.f90` | `species_alternated` includes the ionization rows; the third relaxation placed after the element relaxation; the pass diagnostics extended |
| `src/modules/time_step/certification.f90` | one gated entry per ionization row, with its own tolerance accessor; the anchoring of section 2.4 |
| `src/modules/time_step/steady_newton.f90` | the row registry for `Coupled carrier solve: True`, so the coupled route can carry the same rows; last, and only if the alternation is certified first |
| `src/modules/files_IO/input_read.f90` | the two refusals removed; the resolved set reported; the remaining refusals stated for configurations in which the option means something else |
| `src/modules/files_IO/load_IC.f90` | the restart pairing, which already carries `iontrans`, extended to the set actually carried |
| `src/modules/files_IO/write_setup_report.f90` | the carried set echoed at startup |
| `src/modules/init/parameters.f90` | the key comment: what it now means, which stages it carries, and the validity range |
| `src/modules/post_process/post_process_adv.f90` | the measured difference between the two closures written to `pp.log` and to the file header, so the identity is stated by the run and not only by a test |
| `docs/input_schema.md` | K15f rewritten (it is stale already, section 6) |
| `docs/EXHALE_user_manual.tex`, `README.md`, `README_HOWTO.md` | the key's requirements |

### 5.2 Stages, each with its own measurement

| stage | what it does | what is measured |
|---|---|---|
| A | the single imposed-fraction routine, no behaviour change | `hp_front` byte-identical; `carrier_retry`, `carrier_returned_state_acceptance`, `carrier_reference_scales`, `carrier_constraint_attribution`, `steady_species_rows`, `certification` suites |
| B | the proton row's variable changed to the fraction per H nucleus, advected on the H element flux | `hp_front`, `hp_zero_seed`, `hp_trace_seed` byte-identical up to the round-off of one global rescale, because those cases carry no element diffusion and a uniform elemental composition makes the two variables proportional (MEASURED premise: neither input carries `He_diffusion`) |
| C | the molecular gates lifted; an atomic run carries H II alone | the 45 R_p state re-solved; x(H I) against `_adv` by the measure of section 4.2 |
| D | He II and He III rows added | the same measure on all three; the He 2^3S column and the equivalent width; the 5 percent rule of section 3 |
| E | the certification entries and the tolerance anchoring of section 2.4 | the manufactured-column refinement; the gate set from it |
| F | the new fixture | below |

### 5.3 Regression reach

- **Every golden of a key-off case is untouched**, and the key stays off by
  default. The changes of stages C to F live inside branches the key opens;
  stage B changes a variable inside those branches. The one change outside a
  branch is the imposed-fraction routine of stage A, which is a refactor
  carrying a byte-identity claim.
- **Three cases carry the key and ARE goldens**: `hp_front`,
  `hp_zero_seed` and `hp_trace_seed` are all in `DEFAULT_CASES` and all have
  a `backup/regression/golden/` directory (MEASURED by reading the script and
  listing the directory). They are goldens OF this option, so they are what
  each stage is measured on:
  - stage A must leave them byte-identical;
  - stage B changes the row's variable by a constant factor per run, because
    those cases carry a uniform elemental composition and no element
    diffusion, so the equation is the same equation and only its arithmetic
    moves. The expectation is agreement at round-off, which the standing rule
    of 2026-09-05 counts as identical below 0.1 percent and which the
    harness admits at `REGRESSION_REL_TOL = 1e-3`. The measured relative
    movement is reported; no refresh.
  - stage C does not reach them: they carry molecules, so lifting the
    molecular gate changes nothing there.
  - **stage D moves them deliberately.** `hp_front` carries helium at
    He/H = 0.0793 with `Include He23S? True`, so giving the key its full
    meaning changes what those runs transport. That is a golden refresh at
    the end of the series, taken once, reported with the measured movement
    and with the physical reason: the key names the ionization state and
    helium is part of it.
- `carrier_model_a_newton` also carries the key and is a diagnostic entry
  point, not a golden (READ, its README); it is re-pinned with the series.
- **A new fixture**, `backup/regression/lhs1140b_ionization_transport/`: the
  certified 45 R_p state as `IC/`, the L8 input with
  `Ionization transport: True`, and a `run.sh` of the shape
  `carrier_model_a_newton/run.sh` has. It is the case in which the
  ionization rows, the element rows and the He 2^3S level all run together
  in an atomic gas, which no existing case does.

### 5.4 Expected cost

Stage A one day; B one day including the byte-identity measurement; C two to
three days, since it is where the operator stops being molecular and the
first LHS 1140 b solve has to certify; D two days; E two days, most of it the
manufactured column; F one day. The largest single risk is C: the
alternation has to certify a state with three participants (wind, elements,
ionization) where it certifies two today, and the progress control
(`outer_no_fall_max`, the halving of omega and of the movement bound) was
tuned on two.

---

## 6. The alternative (b): local equilibrium in the solve, the line from the post-process

The user's decision of 2026-09-13 keeps (b) open. This section states what it
is and what separates it from (a), so that it can be picked up without
re-deriving the comparison.

**What (b) is.** The stationary solve keeps the local photoionization
equilibrium it has today. The post-process goes on solving the advected
ionization system it already solves, the transit tool goes on reading
`_adv`, and the He 2^3S level is solved with advection on the corrected
composition, which `System_implicit_adv_HeH_TR` already does. Nothing in the
solver changes.

**What it costs to build: nothing.** It is the present code. The only work
would be to state in the outputs which state each file is, and L10 has
already made that change for `collisional_validity.py`.

**What it cannot give.**

1. The wind is still solved at the wrong composition. Photoheating at the
   corrected composition is 1.5 times the run's own at 10 R_p, 2.7 at 20,
   3.7 at 29 and 5.9 at 44 R_p (READ, L10 section 2), and the energy
   equation of the solve never sees it. Temperature, density, velocity and
   the mass-loss rate are all the wrong composition's.
2. For THIS planet the critical point is at 40.06 R_p (READ, L8), inside
   the region where the closure fails, so the mass-loss rate is set in that
   region and not below it. On a hot Jupiter, where the sonic point sits at
   2 to 4 R_p and Da there is large, this objection is weak and (b) is a
   reasonable engineering answer; on LHS 1140 b it is the whole question.
3. The state is not self-consistent in the sense L10 states: `_adv` carries
   the solution's density and velocity beside a temperature and a
   composition for which the momentum equation was never re-solved, so
   anything formed by mixing the two files belongs to neither.
4. Nothing certifies the corrected composition. The certification measures
   the equations the solve holds; the post-process holds different ones, and
   its only verdict is a per-cell status flag.

**What (b) is still good for after (a) lands.** It becomes the acceptance
test of section 4: on a state solved with ionization transport, the
post-process solves the same balance by an independent discretization, and
the difference is the measure. It is worth keeping working for exactly that
reason, and this design does not retire it.

---

## 7. Noticed outside the scope, reported and not changed

- **`docs/input_schema.md` K15f is stale against the code, on one count,
  and miscounts the refusals.** It states that `Ionization transport: True`
  with `Solver: Newton` also needs `Coupled carrier solve`, and that "Each
  of the three is a fatal `error stop` in `input_read`". That refusal was
  removed on 2026-09-12: the comment at the site says so and gives the
  reason, the partitioned alternation having replaced the plain JFNK finish
  the refusal was written for. Two requirements remain in `input_read.f90`,
  `Molecular chemistry` and `Molecular carrier transport`, and both are
  about to change with this design. The same stale sentence is repeated in
  stage 1 section 169.1 ("The key is refused with `Solver: Newton`,
  `EXHALE_PTC`, `Coupled carrier solve` ..."). What K15f says about
  `EXHALE_PTC=1` is CORRECT: that route is refused, in `EXHALE_main.f90`
  rather than in `input_read.f90`.
- **The comment at that PTC refusal is stale** (`EXHALE_main.f90`, the
  `EXHALE_PTC` block): it justifies itself as "refused for the reason
  input_read refuses `Solver: Newton` with the ionization state carried",
  and `input_read` no longer refuses that. The refusal itself stands on its
  own second sentence, which is the physical one: that route solves
  `F(Y) = 0` with the composition on its own local root. Reported and not
  edited, this being a design item and that file being another worker's.
- **`src/tests/` has no suite for the imposed-fraction substitution.** The
  `x_hp_fixed` path is exercised only through `carrier_retry` and the
  fixtures. Stage A should bring one, since it is about to be the single
  definition seven systems share.
- **A quarantined fixture directory is named `_quarantined/armD_D1`**
  (`backup/regression/`, named in the `run_check.sh` header). The noun
  "arm" is not a name this tree uses any more. Renaming a quarantined
  fixture is not this item's work; it is recorded here so the next pass over
  `backup/regression/` can take it.
  **DONE 2026-09-16**: renamed to
  `_quarantined/molecular_base_no_chemistry`, for the combination its input
  carries and the startup check refuses.
