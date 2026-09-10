# Design: binary H/He diffusion for a helium-rich wind (Phase D)

Status: **implemented through milestone M4.** Written 2026-08-25
for Phase D of `lhs1140b_lower_atmosphere_plan_new.md`; revised the same day
against the external review `binary_diffusion_design_review.md` (section 9
records each finding and what was done with it). The operator is
`src/modules/functions/binary_element_diffusion.f90`; section 6 tracks which
milestones and tests are closed, and section 8 which decisions are still
open.

Supersedes the transport formulation of `design_hehe_diffusion.md`
(sections 2-5), whose kernel is the *trace-helium-in-hydrogen* limit of what
is written here. That document remains the record of the Phase-1/2 build
and validation (its sections 7b-7e).

## 1. Why the present kernel cannot be extended

[2026-08-27: the kernel this section reads is
`src/modules/functions/species_diffusion.f90`, deleted at milestone M3
(section 6), so what follows is a record of code that no longer exists.]

`species_diffusion.f90` (read in full for this memo; every statement below
is from the code, not from its comments) solves one conservative implicit
tridiagonal step for the helium *number density* `n_He` (`solve_1elem`,
lines 280-360), with the hydrogen density `n_H` held fixed during the solve
as a lagged background (lines 93, 149), a diffusion coefficient written for
helium *in* hydrogen, flux coefficients that divide by that background
(`nHf/nHl(j)`, line 294 ff.), a Dirichlet base at the global `HeH`, an
outflow-only top face with zero diffusive flux baked into the coefficients
(lines 326-328), and -- after the solve -- a cap `f_He <= HeH` over the
whole domain (line 166). The comments' "sub-cycled explicit" description
and the old design's `Y` variable do not describe this code. Three
consequences follow, all recorded in `design_hehe_diffusion.md` without
being recognized as one problem:

1. **The variable diverges in the limit the target needs.** `f -> infinity`
   as hydrogen becomes the minor element. LHS 1140 b is retrieved at
   H:He ~ 1e-3 (`f ~ 1e3`), and the EXHALE scan runs to `f = 1e4`.
2. **The self-consistent coupling diverged (P2a, "not adopted").** Making
   `n_H` respond to the helium it loses drove `n_H -> 0` near the base and
   the `1/n_H` flux coefficients blew up; the fix adopted was a cap
   `f <= HeH`. That is the trace assumption failing, not a numerical
   accident: a flux written as `n_H D d(n_He/n_H)/dr` has no reason to
   vanish when hydrogen runs out, whereas the true binary flux does. The
   cap is not dormant: in the stored `examples/14_diffusion` output the
   profile sits *exactly* on it between 1.2 and 2 R_p.
3. **No closure on the mass flux.** The hydrogen flux is *implied*
   (`F_H = -F_He`, section 4 of the old design) but never constructed, so
   nothing in the discretization guarantees that the diffusive mass flux
   through a face sums to zero, which is what keeps `rho` -- owned by the
   hydro -- consistent with the composition it is told.

Extending the kernel ("swap the background for the total gas") was ruled
invalid in the plan review; the reason is (3): a background-species
formulation has no second flux to balance against.

## 2. Formulation

### 2.1 Variables

Two *components* in a single-fluid gas of density `rho` and mass-averaged
velocity `v` (the hydro's velocity *is* the mass-averaged one, since the
hydro evolves total momentum):

- component 1: hydrogen **together with the trace metals slaved to it**
  (the present treatment: metal nuclei move with hydrogen nuclei at fixed
  metal/H). Its mass per hydrogen nucleus is `m_1 = m_H + sum_Z A_Z m_Z`,
  the hydrogen part of `comp_mass_per_H()` (`composition.f90:131`), which
  is `m_H` when `eos_metals` is off or metals are absent;
- component 2: helium, all stages.

With `eos_metals` on, `rho` contains the metal mass, and defining `X`
against the H+He subdensity would leave a third mass flux unaccounted
for. Assigning the metals to component 1 keeps the two-component algebra
exact: `rho_1 + rho_He = rho`, and the diffusive mass fluxes close,
`J_1 = -J_He`, with the metal share of `J_1` being `(rho_met/rho_H) J_H`.
Define

```
X   = rho_He / rho                  helium mass fraction, 0 <= X <= 1
x   = n_He / (n_1 + n_He)           helium fraction of the collision partners
```

where `n_1` counts the hydrogen *carriers* (atoms and ions in the atomic
region; see section 5 for the molecular region) -- the metal nuclei are
left out of `x` and of `D_12` as a trace approximation stated at the
coefficient (their number fraction is ~1e-4 at solar abundance). The
relation `x(X)` uses `m_1` and `m_He`.

The *element* number densities used for the mass budget and the write-back
count nuclei over every species: `n_H = n_HI + n_HII + 2 n_H2 + 2 n_H2+ +
3 n_H3+ + n_HeH+`, `n_He = n_HeI + n_HeII + n_HeIII + n_HeH+` (the `bsp_nH`,
`bsp_nHe` weights of `species_table.f90`). He 2^3S carries no term of its
own: it is a metastable level of He I and the `n_HeI` the code carries is
the total He I population, triplet included, so its nucleus is already in
`n_HeI`. Every budget sum skips the flagged column
(`bsp_is_excited_level` in `species_table.f90`;
`element_nucleus_counts` in `binary_element_diffusion.f90`).

`X` is the transported variable: it is bounded, smooth at both ends of the
composition axis, and `rho X` is the helium mass density.

### 2.2 Transport equation

```
d(rho X)/dt + (1/r^2) d/dr [ r^2 ( rho X v + J ) ] = 0                (1)
```

`J` is the diffusive *mass* flux of helium relative to the mass-averaged
velocity. The hydrogen element obeys the same equation with `1 - X` and
`-J`; adding the two gives the hydro's continuity equation exactly, so no
mass is created or lost by the composition update -- this is the
`sum_i rho_i w_i = 0` condition of the plan, built into the variable choice
rather than enforced afterwards.

### 2.3 The binary diffusive flux

For a binary mixture the diffusion velocities relative to the mass-averaged
velocity satisfy `rho_1 w_1 + rho_He w_He = 0`, so both are set by their
difference:

```
J = rho_He (w_He - v) = (rho_1 rho_He / rho) (w_He - w_1)             (2)
```

and the difference is the binary case of the Chapman & Cowling (1970)
diffusion equation, their eq. (18.2,6) -- the same equation Koskinen et al.
(2013, section 2.1) solve in multicomponent form with the ambipolar field
`eE = -(1/n_e) dp_e/dr` added as a force term -- written for two species
with no assumption about which is minor:

```
w_He - w_1 = - D_12 [ (1/(x (1-x))) dx/dr
                      + ( (m_He - m_1) g - (Zbar_He - Zbar_1) e E ) / (k T)
                      + alpha_T dlnT/dr ]                              (3)
```

with `r` increasing outward, `g > 0` the inward gravity and `E` the
ambipolar field (outward positive). **Sign check:** for neutral gas at
uniform `x`, `w_He - w_1 = -D_12 (m_He - m_1) g/kT < 0`: helium drifts
inward relative to hydrogen, i.e. settles. Every discrete form below must
reproduce this sign (the first draft of section 3 did not; review
finding 1).

- `D_12` is the **binary** coefficient, symmetric in the two species. The
  Banks & Kockarts hard-sphere form already in the code,
  `D_12 = 1.52e18 (1/m_H + 1/m_He)^(1/2) T^(1/2) / n`, *is* symmetric --
  the same number serves He-in-H and H-in-He -- with `n = n_H + n_He` the
  total nucleus density. (The ionized-gas correction remains the
  documented caveat it is now; see 8, D3.)
- **The ambipolar field is computed, not assumed.** The present P2b
  form `Delta m_eff = 3 - 0.5 (Zbar_He - Zbar_H)` is the *hydrogen-plasma*
  value `eE = m_H g/2`, and it is wrong in the limit this phase targets:
  in a fully ionized helium plasma (He++, `n_e = 2 n_He`, `p = 3 n_He kT`,
  `rho = 4 m_H n_He`) the massless-electron balance `dp_e/dr = -n_e eE`
  gives `eE = 4/3 m_H g`, so a proton (charge 1, weight `m_H g`) is
  pushed *outward* and the relative settling mass is `(4 - 2*4/3) -
  (1 - 4/3) = 5/3`, not `2.5`. Equation (3) therefore carries the field
  explicitly, evaluated from the solved electron pressure of the cell,

  ```
  e E = - (1/n_e) d(n_e k T)/dr                                        (3a)
  ```

  (Koskinen et al. 2013, section 2.1; `n_e` and `T` are per cell in
  `get_species_densities`), with the mean charges `Zbar_He`, `Zbar_1`
  from the local ionization state as now. Limits, which become test T9:
  H+ plasma `-> eE = m_H g/2`, `Delta m_eff = 2.5` (the P2b value, so the
  validated He-trace behavior is unchanged); He++ plasma `-> 4/3 m_H g`,
  `5/3`; neutral gas `-> 0`, `3`. Validity: equal electron and ion
  temperatures, zero current -- the assumptions the ionization solver
  already makes. `He_ambipolar` keeps its meaning (on/off).
- `alpha_T` thermal diffusion, default 0 as now.

Eddy diffusion `K_zz` (new in this phase) mixes toward uniform *mass*
fraction and carries no settling term; it adds to (2) as

```
J_total = J - rho K_zz dX/dr                                           (4)
```

which is the aeronomy `(D + K)` form (Banks & Kockarts 1973): the gradient
term carries `D + K`, the settling term `D` alone, because eddy mixing
transports the mixture as a whole and has no preferred species. Its use is
the homopause test of section 6 and the Phase-E matching region; the default is `K_zz = 0`, which reproduces (2).

### 2.4 Both dilute limits

Substituting `rho_1 rho_He / rho` and `1/(x(1-x))`:

- **Helium trace** (`x -> 0`, `rho_1 -> rho`): `J -> rho_He (w_He - w_1)`
  and `(1/(x(1-x))) dx/dr -> (1/x) dx/dr = (1/f) df/dr`. Equation (3)
  becomes the present kernel's drift exactly (`design_hehe_diffusion.md`,
  section 2, "Diffusive drift"). This is the limit the HD 209458 b
  validation exercised, and it must reproduce those numbers (section 6,
  test T2).
- **Hydrogen trace** (`x -> 1`, `rho_He -> rho`): `J -> -rho_1 (w_He - w_1)`,
  i.e. the *hydrogen* diffuses through a helium background with the same
  `D_12`, and the settling term now lifts hydrogen (lighter element rises).
  Nothing diverges: the settling prefactor `rho_1 rho_He / rho -> rho_1
  -> 0` turns the *settling* flux off as hydrogen is exhausted, which is
  the behavior the P2a Picard iteration could not obtain from the trace
  form. The *gradient* flux does not vanish at the ends (section 3) --
  helium must be able to diffuse into a helium-free cell -- so a single
  cell at `X = 0` is not an invariant state; only the uniform states
  `X = 0` and `X = 1` over the whole domain are.

### 2.5 What is deliberately not in this phase

- **Multicomponent momentum** (separate velocity fields, drag on metals):
  out of scope by the plan. Crossover mass and O/C/N drag stay analytic
  diagnostics on the converged wind.
- **Metals in the mass budget of the diffusion**: metals keep the present
  treatment (element-by-element trace diffusion against `n_H` when
  `He_metal_diffusion` is on, otherwise slaved). In a helium-dominated
  gas the "background" for a trace metal should be the mixture; this is
  a follow-on once the two-element core is validated (8, D4).

### 2.6 Stage-resolved friction in the ionized wind

Sections 2.3-2.5 fix the *form* of the transport equation; what is left is
the number `D_12` in it. Milestones M2/M3 used one coefficient everywhere,
the neutral Banks & Kockarts hard-sphere form. That is the wrong friction
above the ionization front, and it is wrong by orders of magnitude in the
direction that matters: a helium *ion* moving through a proton gas is held
by Coulomb collisions, whose momentum-transfer cross section at `T ~ 1e4 K`
is `~ pi (Z_s Z_t e^2/kT)^2 ln(Lambda) ~ 1e-12 cm^2` against the
`~1e-15 cm^2` of the hard sphere. Settling that the neutral coefficient
allows in the ionized wind is therefore suppressed by a factor of order
`1e2-1e3`, and the measured consequence of *not* suppressing it was
`+0.18 dex` in the HD 209458 b mass-loss rate once the `f_He <= HeH` cap was
removed (section 9.3). Koskinen et al. (2013, section 2.1) solve the same
Chapman & Cowling equation with collision terms that "account for
neutral-neutral, resonant and non-resonant ion-neutral, and Coulomb
collisions", and report (their section 3.2.2) that Coulomb collisions "are
much more efficient in preventing diffusive separation than collisions with
neutral H."

This is a correction to the physics, not an option: there is no input key
for it and no way to select the neutral-everywhere behavior.

**Pair coefficients.** Each element is resolved into its ionization stages
and the friction is built pair by pair. Every pair coefficient is the
Chapman-Enskog first approximation

```
D_st = 3 k T / (16 n mu_st Omega_st^(1,1))                              (7)
```

with `n` the total (nucleus) density, `mu_st` the reduced mass and
`Omega^(1,1)` the standard collision integral, so the three cases below sit
in one framework and are directly comparable.

- **Neutral-neutral** -- unchanged, the hard-sphere form already in the
  code (Banks & Kockarts 1973):

  ```
  D_st = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2) / n                      (8)
  ```

  Putting the rigid-sphere `Omega^(1,1) = pi d^2 (kT/2 pi mu)^(1/2)` into
  (7) reproduces (8) with `d = 2.7 Angstrom`, which is the check that (8)
  and the two coefficients below are the same approximation.

- **Ion-neutral, non-resonant** -- two channels of one interaction, added.
  The long-range channel is the polarization (Langevin) interaction
  `V = -alpha_n e^2/(2 r^4)`. For that potential `g Q^(1)(g)` is independent
  of the relative speed, `g Q^(1) = 2.21 pi e (alpha_n/mu_st)^(1/2)`, which
  is the constant behind the non-resonant ion-neutral collision frequency of
  Schunk & Nagy (*Ionospheres*, eq. 4.88),
  `nu_in = 2.21 pi (n_n m_n/(m_i+m_n)) (gamma_n e^2/mu_in)^(1/2)`. Inserting
  it into (7) -- the Chapman-Enskog numerical factors cancel exactly for a
  Maxwell-molecule potential -- gives

  ```
  D_pol = k T / (2.21 pi e n (alpha_n mu_st)^(1/2))                     (9a)
  ```

  with `alpha_n` the static dipole polarizability of the neutral partner in
  `cm^3` and `e` in esu. That the imported `2.21 pi` really is a
  momentum-transfer cross section, and not a collision frequency in some
  other convention, is checked against the orbiting (capture) cross section
  of the same potential: `b_L^4 = 8 C_4/(mu g^2)` with `C_4 = alpha e^2/2`
  gives `g sigma_L = 2 pi e (alpha/mu)^(1/2)`, the classical Langevin rate,
  and the constant above is `2.21/2 = 1.105` times it -- capture plus the
  ~10% the glancing collisions add.

  (9a) keeps no repulsive core, so its friction falls as `T^-1` where a
  rigid core would hold it at `T^-1/2`, and taken alone it would let an
  ion-neutral pair become *more* mobile than a neutral-neutral one of the
  same masses above ~1.5e3 K. The two are momentum-transfer cross sections
  of the same encounter -- the long-range attraction and the short-range
  repulsion of one potential -- so to first order their `Q^(1)` add, hence
  their collision integrals add, hence their **frictions** add. Since
  `D = 3kT/(16 n mu Omega^(1,1))` is linear in `1/Omega`, adding frictions
  is adding inverse coefficients:

  ```
  1/D_st = 1/D_pol + 1/D_hs                                             (9b)
  ```

  the same "frictions add" rule the stage mixture (11) uses. The physical
  content is that opening a second channel cannot make a pair more mobile:
  `D_st` never exceeds either channel and reproduces each limit exactly --
  `-> D_pol` where the polarization friction dominates, `-> D_hs` where the
  core does. For an H/He pair the two channels cross at **1504 K**
  (`D_pol/D_hs = 0.0258 T^(1/2)`), so the base is in the polarization limit
  and the `1e4 K` transition layer is core-dominated: there `D_pol` is 2.6
  to 4.7 times `D_hs` and (9b) sits at 0.82 of `D_hs`, i.e. 0.18-0.28 of the
  polarization value alone. Tests T12c and T12d; measured in 9.3.

- **Ion-ion** -- Coulomb, with the momentum-transfer cross section
  `Q^(1) = 4 pi b_90^2 ln(Lambda)`, `b_90 = Z_s Z_t e^2/(mu_st g^2)`.
  Carrying it through the `Omega^(1,1)` integral gives

  ```
  D_st = 3 (k T)^(5/2) / [ 4 (2 pi mu_st)^(1/2) n (Z_s Z_t e^2)^2 ln(Lambda) ]
                                                                       (10)
  ln(Lambda) = ln( 3 k T lambda_D / (Z_s Z_t e^2) ),
  lambda_D   = ( k T / (4 pi n_e e^2) )^(1/2)
  ```

  the `T^(5/2)/(n Z^2 Z^2 ln Lambda mu^(1/2))` scaling of the standard
  Chapman-Enskog Coulomb result (Paquette et al. 1986; Schunk & Nagy give
  the equivalent collision frequency, their eq. 4.142). Neither textbook is
  in `references/`, so the numerical constants of (9) and (10) were
  **re-derived here from the Chapman-Enskog collision integral** rather than
  copied; the derivation is recorded in the code comment beside each
  coefficient, and (8) reproducing itself through (7) is the check that the
  chain is consistent.

- **Resonant charge exchange `H+ + H`** is deliberately *absent*. It is a
  collision between two carriers of the same element -- both belong to
  component 1 -- and the transport equation (3) is driven by the friction
  *between* the two elements only. Internal friction within a component
  does not enter a binary diffusion coefficient. (The same argument removes
  `He+ + He` resonant exchange.) So every ion-neutral pair that does appear
  in the He-H friction is non-resonant, and (9) is the right form for all of
  them.

**Mixture average.** The element still moves as one body (a multi-fluid
treatment, one velocity per stage, remains out of scope by 2.5), so the
element-element friction is the stage-fraction-weighted sum of the pair
frictions -- the Blanc's-law structure already adopted for the molecular
carriers in section 5, now applied across ionization stages as well:

```
1/D_eff(He,1) = sum_{s in He} sum_{t in 1} y_s y_t / D_st              (11)
```

with `y_s` the fraction of the element's *carriers* in stage `s`
(`sum_s y_s = 1` within each element). Frictions add, so it is the inverse
coefficients that are averaged. Limits: an all-neutral gas leaves exactly
one term, `D_eff = D(HeI,HI)` of (8); a fully ionized one leaves exactly
`D_eff = D(He++,H+)` of (10); a half-ionized state is the four-term average.
Those three are test T12. Because every pair coefficient carries the same
`1/n`, so does `D_eff`.

The metal loop has the same structure with the metal's own stages against
the hydrogen carriers, so a metal ion in the ionized wind is held to the
protons by (10) exactly as helium is. `HeH+` carries both elements and is
left out of both carrier lists (as it already is out of the mean charges),
a trace approximation stated at the coefficient.

**Size of the three coefficients.** Measured by T12 on a prescribed column
at `T = 1e4 K`, `n = 2.1e8 cm^-3` of nuclei: `D(HeI,HI) = 2.09e12`,
`D(HeI,H+) = 1.72e12`, `D(He+,HI) = 1.51e12`, `D(He++,H+) = 9.2e9 cm^2/s`
(the polarization channel alone would give 9.73e12 and 5.39e12 for the two
ion-neutral pairs). So an ion-ion pair is `~2.3e2` times slower than the
neutral one while an ion-neutral pair is slightly slower than neutral-
neutral rather than several times faster, and the net effect of (11) is
that settling is nearly frozen wherever both elements are ionized and
essentially unchanged in the neutral base. The measured `D_eff/D_neutral`
across a real front is in section 9.3.

## 3. Discretization

Finite volume on the existing grid (`r`, `r_edg`), faces carry `r^2`
areas, exactly as the current kernel. The composition step is operator
split from the hydro and runs once per relaxation step, after the hydro
has produced the new `rho`, `v`, `T` (`EXHALE_main.f90:601`).

**Transient advection: the advective form.** The routine does not receive
the hydro's Runge-Kutta face mass fluxes or the old density, so a second,
independent conservative advection of `rho X` would not in general keep
a uniform `X` uniform in a compressing or expanding flow (review finding
2). Rather than re-plumb the hydro, this phase transports `X` in its
*advective* form,

```
dX/dt + v dX/dr = -(1/(rho r^2)) d/dr [ r^2 J ]                        (5)
```

discretized with a first-order upwind difference on `X` taken with the
**cell** velocity `v_j`, inside the same tridiagonal as the diffusive part.
(The first implementation used the face velocities `v_f = 0.5 (v_j +
v_{j+1})`; section 9 records why that is wrong and what it did.)
Consequences, stated so
they are tested rather than assumed: (i) a uniform `X` is preserved
*exactly* whatever the hydro does (T6, by construction); (ii) at a steady
state `rho v r^2` is constant, so (5) and the conservative form (1) are
the same equation and the elemental face flux is constant (T8); (iii) in
a transient, helium mass is conserved only to the order of the hydro's
own truncation -- acceptable because the composition is driven to a
*steady* state with the wind (D2), and stated as the limitation of this
phase. Evolving `rho X` inside the Runge-Kutta stages as a fourth
conserved variable is the exact alternative; it is recorded as the
follow-on if a transient application ever needs it.

**Diffusion, settling, eddy** as a conservative implicit (backward-Euler,
tridiagonal) solve in `X`, with face coefficients

```
A_f = (m_1 m_He n / mbar)_f D_12,f (dx/dX)_f              [gradient]
B_f = (rho_1 rho_He / rho)_f D_12,f
      ( ((m_He - m_1) g - (Zbar_He - Zbar_1) eE) / kT
        + alpha_T dlnT/dr )_f                             [settling]
E_f = rho_f K_zz                                           [eddy]
```

so that

```
J_f = -(A_f + E_f) dX/dr|_f  -  B_f                                     (6)
```

-- with the **minus** on `B_f`, which is what equations (2)-(3) give and
what the sign check in 2.3 requires (the first draft had `+ B_f`; review
finding 1). Writing `B_f = beta_f [X(1-X)]_f` with
`beta_f = rho_f D_12,f G_f`, the settling part keeps the Peclet-based
central/upwind hybrid of P2 (central where `|beta_f| dr <= 2 (A_f + E_f)`,
donor-cell upwind otherwise).

**The two factors of `X(1-X)` come from opposite sides of the face**
(2026-08-28; `docs/Update_EXHALE_stage1.md` section 86). The settling flux is a
counter-flow -- the helium flux is matched by an equal and opposite
hydrogen flux, since the two components close -- so a donor-cell rule has
to take each element's mass fraction from the cell *that element* leaves,
and the two elements leave from opposite sides:

```
[X(1-X)]_f = X(j+1) (1 - X(j))     beta_f >= 0   (helium drifts inward)
[X(1-X)]_f = X(j)   (1 - X(j+1))   beta_f <  0   (helium drifts outward)
```

with both factors at the new time level. In the central branch the
product is taken at the face average, `X_f (1 - X_f)`. The first draft of
this section lagged the `(1-X)` factor and solved a step linear in `X`
with one Picard sweep; that is what lost the upper bound (below). The
step is now nonlinear in `X^{new}` and is solved by Newton, whose
tridiagonal Jacobian carries the two face slopes `dJ_f/dX(j) >= 0` and
`dJ_f/dX(j+1) <= 0` -- the M-matrix condition -- for any iterate in
`[0,1]`.

**Bounds `0 <= X <= 1`, from the shutoff of the drift flux at both ends
of the composition axis** (rewritten 2026-08-28; `docs/Update_EXHALE_stage1.md`
section 86 -- the two earlier versions of this paragraph are described at
the end of it). The continuum drift flux vanishes at `X = 0` and at
`X = 1`, because a cell with no helium has none to send and a cell with
no hydrogen has nothing to send it in exchange; the donor/acceptor rule
above is that statement in the discrete flux. A donor at `X = 0` sends
nothing and an **acceptor** at `X = 1` receives nothing. In the central
branch the Peclet condition does the same work: with `X(j) = 1`,
`X_f (1-X_f) = (1 + X(j+1))(1 - X(j+1))/4 <= (1 - X(j+1))/2`, so
`|beta_f| <= 2 (A_f + E_f)/dr` bounds the drift flux by the gradient flux
`(A_f + E_f)(1 - X(j+1))/dr` that carries helium *out* of that cell at
the same face.

The bound then follows cell by cell. Let `X` solve the implicit step and
suppose it first touches 1 in cell `m`, every other cell still in
`[0,1]`. The time term `rho (X_m - X_m^old)/dt >= 0`; the gradient flux
leaves `m` at both faces because `X_m` is the maximum; the upwind
advection contributes `rho|v| (X_m - X_donor)/dr >= 0`; and the drift
flux is an outflow or exactly zero at both faces. Every term of the row
has the same sign and their sum is the row residual, which is zero, so
none of them can be strictly positive and `X_m > 1` is impossible.
`X >= 0` is the same statement, because the scheme is exactly symmetric
under `X -> 1 - X`, `beta -> -beta`: donor and acceptor exchange roles
and the central product is even about `X_f = 1/2`. That symmetry is the
binary symmetry of (2) itself, and it is what the lagged form broke --
written for `Y = 1 - X` the lagged flux carried the constant `W_f` and
not a multiple of `Y`, which is why one end held in practice and the
other did not.

The **row sum is not `rho/dt`**, and the argument above does not need it
to be: it is made on the residual of the nonlinear step, not from a
comparison principle for a linear system. The Newton Jacobian is
nevertheless an M-matrix on `[0,1]`, with diagonal
`rho/dt + (outflow face terms) > 0` and off-diagonals `<= 0` for the
upwind advection, for the gradient term (`A_f + E_f >= 0`) and for both
drift branches; that is what makes the correction well posed.

What is measured, before the clip (`he_fraction_over_one` /
`he_fraction_under_zero`, exposed so the tests read the solve and not its
clip): the excursion outside `[0,1]` is **negative in every case tested**,
i.e. the solve stays strictly inside the range and does not reach the
clip -- `-9.2e-4` on the strong-drift column of T14, `-8.6e-3` on the
pure-helium band of T13, `-1.05e-1` on the LHS 1140 b wind. On the lagged
form the same T14 column overshot by `+0.514` (`X = 1.514`), and an
8e3-step LHS 1140 b solve exceeded 1 in 140 of 7962 steps by up to 0.520
(section 84). The clip at 0 and 1 is kept and is now an assertion at both
ends.

   The gradient coefficient looks like `0/0` at the ends of the
   composition axis but is not: with `rho_He = n_He m_He`, `rho_H = n_H
   m_H`, `x = n_He/n`, `1 - x = n_H/n`,

   ```
   rho_1 rho_He / (rho x (1-x)) = m_1 m_He n / mbar,     mbar = rho / n
   ```

   finite and smooth for all `X` in `[0,1]`. The implementation uses this
   closed form, so no division by a vanishing element density occurs
   anywhere (the P2a divergence divided by `n_H`; here nothing does).
   What *does* vanish at the ends is the settling flux, through its
   `rho_1 rho_He / rho` prefactor: no helium, no helium to settle. The
   gradient flux does **not** vanish there -- a helium-free cell next to
   a helium-bearing one receives helium, as it must -- so only the
   uniform states `X = 0` and `X = 1` over the whole domain are
   invariant. The bounds on `X` do come from the vanishing of the drift
   flux, but at the FACE and not over the domain: it is the shutoff of
   the donor and of the acceptor separately (above) that keeps a cell
   from being pushed past either end. Both factors of the product are
   evaluated at the new time level; the Newton iteration is what carries
   them there.
**Write-back to the species vector, nonnegative by construction.** After
the step the cell has new element totals `n_He^new = rho X/m_He` and
`n_H^new = rho (1 - X)/m_1` (metal nuclei follow at fixed metal/H). The
species vector `f_sp` is projected onto these totals so that (a) no
species goes negative and (b) both element totals are met exactly:
helium-only species are scaled by `r_He = n_He^new/n_He^old`,
hydrogen-only species (including `H2`, `H2+`, `H3+` in the molecular
region) by `r_H`, and a species carrying both elements (`HeH+`) by
`min(r_H, r_He)`; the nuclei that this under-counts for the element with
the larger factor are deposited into that element's neutral ground
species (`HI` or `HeI`). Every operation is a multiplication by a factor
`>= 0` or an addition, so no negative intermediate can arise (review
finding 7). This projection is the *initial guess* handed to `ioniz_eq`,
which re-solves the partition from the element totals on the next call
(`EXHALE_main.f90:609`); the split within an element is what `ioniz_eq`
owns, the element totals are what the diffusion step owns.

**The metals are a ratio, not a density** (2026-08-28; the defect this
states is `docs/Update_EXHALE_stage1.md` section 84). The metals are part of
component 1 at a fixed metal/H, so what the projection has to preserve
for them is `n_X/n_H`, not `n_X`. Multiplying them by `r_H` does that --
in a cell that *had* hydrogen. Where the cell had none, `r_H` is `0/0`;
the code took it as `0`, and that is where the metals were lost. The
sequence is: `X` reaches `1`, so `n_H^new = 0`, so `r_H = 0` and every
hydrogen species and every metal ion in the cell is multiplied by zero.
Zeroing them is *correct at that instant* -- no hydrogen, no
hydrogen-slaved metals, and component 1 carries no mass. What is not
correct is the step after. Hydrogen returns to the cell through the
shortfall deposit above; the metals had no such deposit, so they stayed
at zero while every later call multiplied that zero by something, and a
band of cells at exactly zero `C`, `N` and `O` sat between cells at the
reservoir ratio -- a state with no sink for it in these equations.

The metals therefore return **with** the hydrogen. When `n_H^old` is
zero the cell holds no metal/H of its own, and the ratio is taken from
`melem_ab`: it is the metal/H this operator's own mass budget already
assumes, since `m_1 = ` `mass_per_H_nucleus_without_He()` `= m_H + sum_X
melem_ab_X A_X`, and it is the same reservoir ratio `set_IC` and
`load_IC` build an absent element from (section 80). All of it goes into
the neutral ground stage, for the same reason the hydrogen and helium
shortfalls go into `HI` and `HeI`: a cell that held no hydrogen held no
ionization split either, and `ioniz_eq` re-solves the split from the
element total on the next call. Where `n_H^old > 0` nothing changes, so
a state that never empties a cell of hydrogen is untouched (`make check`
byte-identical, 7/7). The repair is to the *creation* of such a band;
it does not reconstruct one that a restart file already carries, which
is why the operator now also counts them (below).

**The metal census.** No convergence test in this pipeline could see the
band: the steady residual is a residual of the hydro and energy
equations, and their solution with the metals removed is a perfectly
good solution of the equations as posed; the elemental-flux closure
measures a window that need not contain the cells concerned. The
operator therefore checks its own elemental bookkeeping on every step
(`metal_hydrogen_ratio_departure`): the departure of `n_X/n_H` from
`melem_ab` is reported under `EXHALE_DIFFUSION_CHECK=1`, and a
`(cell, element)` pair holding no metal nuclei at all while hydrogen is
present -- which these equations cannot produce -- raises a warning
whether or not the diagnostic is on. T13 is the acceptance test.

The cap `f <= HeH` and the base pile-up limiter of Phase 1 are
**removed**; if a pile-up reappears it is physics (helium settling against
a slow base) to be answered by the base boundary condition, not clamped.

## 4. Boundary conditions

- **Base:** Dirichlet on `X` at the reservoir composition, `X = X(HeH)`
  (or, with a `base.inp` handoff, `X(HeH_base)`): the lower atmosphere is
  the well-mixed reservoir. Unchanged in kind from Phase 1; what changes is
  that it is imposed on a bounded variable.
- **Top:** zero diffusive flux, and for advection the *zero-gradient*
  outer ghost `X(N+k) = X(N)`: an outflowing top face carries the
  domain's own composition out; an inflowing face (`v(N) < 0`, transient
  only) brings in gas of the same composition as the top cell, so no
  element enters or leaves through the top by diffusion and the
  advective elemental flux through the top is `rho X v r^2` with the
  local `X`. This is one boundary condition, stated once (the first draft
  said "outflow-only" in one place and "outer-ghost composition" in
  another; review finding 8). The mass budget of T4 accounts for it.

## 5. Removing the molecular-chemistry / `He_diffusion` exclusion

`input_read.f90` refuses the pair (line ~909). The survey (section 7) will
state the concrete reason; on the physics the resolution is not in doubt:

**Option A (recommended architecture) -- one element-transport
formulation through the molecular region.** Equation (1) transports
*elements*; chemistry moves hydrogen between H, H2, H2+, H3+, HeH+ and
helium between its stages without changing either element's mass, so the
transport equation is the same above and below the molecular front. What
changes with the species set is the *closure* -- the collision partners
of helium are then a mixture of carriers, and a single mean carrier mass
put into the He-H formula is not a derived closure (review finding 6).
The closure adopted here, with its validity stated:

- **Friction: Blanc's law.** Helium against a background made of carriers
  `c` in {H, H+, H2, H2+, H3+} with carrier fractions `y_c`
  (normalized within the hydrogen component) has the effective binary
  coefficient `1/D_{He,1} = sum_c y_c / D_{He,c}`, each `D_{He,c}` the
  Banks & Kockarts hard-sphere form with the carrier mass. Blanc's law is
  exact for a trace species in a mixture and the standard approximation
  otherwise; it reduces to `D_{He,H}` in the atomic region and to
  `D_{He,H2}` in the fully molecular one, and interpolates between (the
  two limits are the analytic checks of T7).
- **Thermodynamic and force terms per carrier.** In (3) the fraction `x`
  counts helium against *carriers* (`n_1 = sum_c n_c`, so an H2 counts
  once), and the settling difference uses the mean carrier mass `<m_c>`
  and mean carrier charge `<Z_c>` of the hydrogen component in place of
  `m_1`, `Zbar_1`. This is what makes gravity act on the particles that
  actually collide (per nucleus the weight of H2 is the same as that of
  H; per collision partner it is twice).
- **The mole-fraction driver the chemistry carries** (added at M4; the
  paragraph above is not complete without it). The transported variable
  is the mass fraction `X`, while the force in (3) is the gradient of the
  *mole* fraction `x`. Writing `psi = n_1/n_H` for the collision partners
  each hydrogen nucleus is spread over (1 atomic, 1/2 fully H2, so
  `<m_c> = m_1/psi`),

  ```
  logit(x) = logit(X) + ln( m_1/(m_He psi) ),
  d logit(x)/dr = d logit(X)/dr - dln(psi)/dr                          (12)
  ```

  The gradient coefficient itself is unaffected -- `A = (m_1 m_He n/mbar)
  D_12 (dx/dX)` collapses to `rho D_12` with the carriers exactly as it
  does with the nuclei, for any `psi`, because the two carrier factors
  cancel -- so the whole content of (12) is one extra term in the
  settling coefficient,

  ```
  G -> G - dln(psi)/dr .                                               (13)
  ```

  It is a real driver, not a change of variable: where hydrogen turns
  molecular going down, each nucleus is spread over fewer collision
  partners, helium's mole fraction rises, and helium diffuses down that
  gradient -- **outward across the molecular front** -- even at a uniform
  mass fraction. With no other force at all the steady state of (13) is
  the statement the whole closure rests on, that diffusion levels the
  mole fraction and not the mass fraction; that is test T7d, and across a
  front where `psi` runs 1/2 -> 1 it moves the He/H *nucleus* ratio by
  the same factor 2. The eddy term does not carry it (eddy mixing has no
  preferred species and acts on `dX/dr` alone), and the same substitution
  applies to the trace-metal loop, whose transported variable `n_X/n_H`
  differs from its mole fraction `n_X/n_1` by the same `psi`.
  `dln(psi)/dr` vanishes identically wherever the hydrogen is atomic.
- `HeH+` is counted as a helium carrier for the friction (its abundance
  is negligible wherever the closure matters) and handled in the
  write-back as in section 3.
- **The carrier density in the pair coefficients.** `n` in (7)-(10) is
  the Chapman-Enskog density of *colliding particles*, so it is the
  carrier density `n_1 + n_He`, not the nucleus density. The two are the
  same number in the atomic region; below a molecular front the carrier
  density is smaller and every pair coefficient correspondingly larger.

**Staging.** The core of this phase (milestones M1-M3 in section 6) is
built and validated in the *atomic* region, where the closure is the
plain binary one and no approximation beyond `D_12` enters. The molecular
closure is milestone M4: implemented behind the same key, gated by T7,
and the `error stop` is lifted only when T7 passes.

*Status (2026-08-26): done. T7a-T7d pass, and the `error stop` of
`input_read.f90` is gone -- the pair is an ordinary configuration.*

**Option B -- a transition radius.** Keep the atomic kernel above the
molecular front and freeze the composition below it at the reservoir
value. Simpler to code, but it puts the homopause *at* the molecular
front by construction, which is exactly the quantity Phase E needs to
compute rather than assume, and it makes the answer depend on a
user-chosen radius. Not recommended.

Under A the `error stop` is replaced by nothing: the pair becomes an
ordinary configuration, exercised by test T7.

## 6. Tests (acceptance)

Each is a short run or a unit check with a stated pass criterion; all
default-off paths must leave `make check` byte-identical (T0).

| id | test | pass criterion |
|---|---|---|
| T0 | regression with `He_diffusion` off | the five off cases byte-identical; since 2026-08-26 the matrix also carries `mol_diffusion` with the flag on (7.5), so `make check` is 6/6 [2026-08-27: 7/7 -- `lower_profile` was added as the seventh default case] |
| T1a | **diffusive equilibrium**: `v = 0`, isothermal column, no eddy, zero diffusive flux at *both* ends; start uniform | `X(r)` relaxes to the barometric separation `dx/dr = -x(1-x) (m_He - m_1) g/kT` (neutral) to grid order; total He mass constant to round-off (the closed-column version of T4) |
| T1b | same with the Dirichlet reservoir base of section 4 | the integrated base flux `int 4 pi r_b^2 J_b dt` accounts for the change of total He mass to round-off |
| T2a | **convergence to the trace equation**: `examples/14_diffusion` configuration (HD 209458 b) at He/H = 0.0833, 1e-2, 1e-3, 1e-4, run with the new operator and with the present kernel | the difference between the two (He/H)/HeH profiles decreases with He/H at the expected first order in `x`; at 1e-4 it is below 0.1% |
| T2b | **report** at the tutorial abundance He/H = 0.0833 (not a gate: the present profile is capped over 1.2-2 R_p and the abundance is not asymptotically small) | (He/H)/HeH profile and He 10830 depth tabulated against a **fresh, Newton-finished run of the present code made before any Phase-D change** -- the stored `examples/14` output is unconverged (JFNK info=2, du 1.3e-2) and in an older column schema, and the tables in `design_hehe_diffusion.md` 7d/7e are from other runs; neither is a usable reference |
| T3 | **H-trace limit**: same planet with He/H = 1000 | finite everywhere, `X` in `[0,1]`, hydrogen rises (X decreasing outward relative to the reservoir), no cap active |
| T4 | **elemental conservation, closed column**: zero flux at both ends, `v = 0`, arbitrary initial `X(r)` | `int rho X r^2 dr` constant to 1e-12 relative over 1e4 steps |
| T5 | **zero net diffusive mass flux at every face**: instrument `J_He + J_H` | identically zero (by construction; the test guards the discretization) |
| T6 | **uniform mixture preserved**: uniform `X`, `g = 0`, any `v(r)` | `X` stays uniform to round-off (advection of a uniform field is continuity) |
| T7 | **molecular closure and homopause**, milestone M4. T7a/T7b: the Blanc carrier coefficient in a forced-atomic and a forced-molecular column. T7c: a prescribed molecular column at a series of constant `K_zz`. T7d: the same column with `g = 0`, where `-dln(psi)/dr` of (13) is the only driver left | T7a/T7b reproduce `D_{He,H}` and `D_{He,H2}` to round-off. T7c: the He/H profile is flat where `K_zz >> D_{He,1}` and separates above, and the homopause read off the profile (where the measured `-d logit(X)/dr` is half the diffusion-limited `G`, i.e. where `D/(D+K) = 1/2`) sits at the radius where `D_{He,1} = K_zz`, within a cell. T7d: the **mole** fraction levels out, not the mass fraction |
| T8 | **moving-wind face flux**: converged LHS 1140 b wind, diffusion on | elemental face flux `4 pi r^2 (rho X v + J)` constant with radius above the base to the same tolerance as the mass flux (`du`) |
| T9 | **ambipolar limits** from (3a) on prescribed states: neutral; fully ionized H+ plasma; fully ionized He++ plasma; a partially ionized front | relative settling mass `3`, `2.5`, `5/3` to round-off in the three limits; finite and continuous through the front |
| T10 | grid and time-step convergence on T1 | second-order in `dr` for the gradient term, first-order where the settling hybrid upwinds |
| T12 | **stage-resolved friction** (2.6) on a prescribed isothermal column at three ionization states, plus the two temperature limits of the ion-neutral pair | `D_eff` reproduces, to round-off, the Banks & Kockarts hard sphere when all neutral, the He++/H+ Coulomb coefficient when fully ionized, and the four-term Blanc average when half ionized -- each compared against a closed form written out independently in the test, not against the module's own pair routines. The combined ion-neutral form (9b) reduces to the polarization channel at 1 K and to the hard sphere at 1e9 K and never exceeds either |
| T13 | **the metals come back with the hydrogen** (2026-08-28): a closed metal-bearing column seeded with a band of pure helium, so `n_H = 0` and `X = 1` exactly in those cells; diffusion refills the band with hydrogen | no cell that holds hydrogen holds zero metals, and `n_X/n_H` equals `melem_ab` to 1e-12 over the whole column. Fails on the pre-fix projection (21 of 21 band cells with hydrogen and no metals, departure 1.0). Also asserts the excursion outside `[0,1]` before the clip, measured where the solve starts ON the bound |
| T14 | **the bounds, read before the clip** (2026-08-28): a closed column at 20 Jeans parameters, steeply stratified, driven 400 steps at `1e6` cell diffusion times, so the drift crosses many cells per step and the base is driven to pure helium; the test reads `he_fraction_over_one` / `he_fraction_under_zero`, which the operator sets from the solve itself | the run must REACH the boundary (max `X >= 0.99`, else the bound asserted is vacuous) and the excursion outside `[0,1]` must stay at or below `1e-12` at both ends (measured `-9.2e-4` and `-8.2e-11`). Fails on the lagged linearization at `+0.514`, i.e. `X = 1.514` |

T2a is the gate that protects what already works (the old kernel is the
comparison, not the truth standard); T3 is the one that shows the new
formulation does what the old one could not.

**Milestones.** M1 baseline capture (D7) and the two prerequisite defects
of 7.3-7.4. M2 the atomic binary operator with (3a), tests T0-T6, T8-T10.
M3 the LHS 1140 b He-rich runs with diffusion (the Phase F precursor,
Newton-finished through the repaired direct-steady route). M4 the
molecular closure and T7; only then the exclusion is lifted.

*Status (2026-08-26): M4 complete -- the molecular closure of section 5 is in
`binary_element_diffusion.f90` (carrier lists, mean carrier mass and charge,
carrier density in the pair coefficients, and the mole-fraction driver (13)),
T7a-T7d pass in the driver, and the `error stop` that refused molecular
chemistry together with `He_diffusion` is removed from `input_read.f90`. The
`He_diffusion` regression case of 7.5 is **added and golden-snapshotted** the
same day (`backup/regression/mol_diffusion`, in the `make check` default set;
`docs/Update_EXHALE_stage1.md` section 73). What follows is the M1-M3 record.*

*Status (2026-08-25): M1, M2 and M3 complete. The three defects of
7.3-7.4 are fixed (`docs/Update_EXHALE_stage1.md` section 67); the atomic binary
operator is `src/modules/functions/binary_element_diffusion.f90`, with T0 and
T1a-T6, T9, T10 passing in `make diffusion_tests && ./diffusion_tests.x`
(section 68) and the T2a sequence measured there (3.05e-6 at He/H = 1e-4,
first order in the helium fraction). M3 (section 69): the baseline capture of
D7 is complete in `backup/phase_d_baseline/`, T2b is tabulated against it,
**T8 passes** (elemental face flux as flat as the mass flux in the code's own
escape window, 5.1e-3 against 7.2e-3), the three LHS 1140 b helium-rich cases
run to `info = 0` with diffusion on, and the trace kernel
`species_diffusion.f90` is deleted -- so T2a stands as a recorded measurement
and is no longer runnable. Decision D3 is
closed as of the same date: the friction is stage-resolved (section 2.6),
T12 passes, and the measured effect is section 9.3.*

## 7. Code integration (from the read-only survey of 2026-08-25)

### 7.1 What already works per cell, and stays

- Composition lives only in `f_sp(1-Ng:N+Ng, n_species)` (mass-normalized,
  `sum_s m_s f_s = 1`, masses in `species_table.f90:90-92`). Densities,
  electrons, particle counts (`composition.f90:35-88`), `T <-> p`
  (`:92-106`), the ionization solver (`ionization_equilibrium.f90:236-241`
  builds the nucleus totals `nh`, `nhe` of each cell *including the molecular
  species*, and hands them to `ieq_cell`, `:409, :494-495`) and the
  advection-corrected post-processing (`post_process_adv.f90:509-515`,
  `heh_loc` per cell when `he_diffusion`) all work cell by cell already. None
  of this changes.
- The base scalars `rho_bc`, `ntot_bc`, `dp_bc` are built from the global
  `HeH` (`composition.f90:131,148,190`, used in `Apply_BC.f90:41,76-78`).
  With the Dirichlet base of section 4 the base composition *is* `HeH`,
  so these stay consistent; they must not be touched.

### 7.2 The module: rewrite, not patch

`species_diffusion.f90` is replaced by a module named for what it computes
(`binary_element_diffusion`, routine `element_diffusion_step`), keeping:
the Thomas solve, the Peclet central/upwind hybrid, the ambipolar mean
charges (lines 109-118), `alpha_T`, and the trace-metal loop over elements
(lines 203-268) unchanged behind `He_metal_diffusion`. Replacing:

- the state: `X` from element mass fractions built with the nucleus
  counts `bsp_nH`, `bsp_nHe` of `species_table.f90:98-` over *all* species
  (atomic and molecular), not `HI+HII` / `HeI..HeTR` (lines 83-87) -- this
  is the first of the three molecular conflicts;
- the flux coefficients: section 3, the closed form `m_1 m_He n/mbar` for
  the gradient and `(rho_1 rho_He/rho)` for the settling, instead of
  `nHf/nHl`; `n` in `D_12` is `n_H + n_He` nuclei (metals
  stay excluded from `n` as now -- a documented trace approximation);
- the top face: keep outflow-only advection with zero diffusive flux, but
  write it as a boundary condition on the face, not inside the
  coefficient table, so that an inflow face (`v(N) < 0`) carries the
  outer-ghost composition instead of a silent zero flux;
- the write-back (lines 177-195): the nonnegative projection of section
  3 (scale by element, `HeH+` by the smaller factor, deposit the shortfall
  into `HI`/`HeI`) -- this is the second and third conflict (`cmet`
  mistaking molecular mass for metals, line 87; one factor cannot serve a
  two-element molecule, the reason `load_IC.f90:204-211` also refuses
  HeH+). The projection is only a first guess: `ioniz_eq` re-partitions
  the species from the element totals on the very next call
  (`EXHALE_main.f90:601 -> :609`);
- the cap (line 166) and the base pile-up limiter: removed;
- `he_kzz`: keep the key, change the default from `1e9` to `0`
  (`parameters.f90:152`) -- with the present default every diffusion run
  silently carries a constant eddy term, and the T1/T3 tests need the
  pure molecular case. User-visible change: an input that relied on the
  default must now say `He_Kzz: 1.0e9` (`examples/14_diffusion` already
  does).

### 7.3 Coupling to the steady solver

The JFNK/PTC residual has no diffusion (`steady_newton.f90:195-246`);
coupling is the outer Picard loop `it_diff` in `EXHALE_main.f90:941-975`
(<= 5 passes x 500 relaxation steps, exit when the He/H field moves by
< 1e-3). Decision D2 keeps this.

**The relaxation step size is a composition time scale, not the hydro CFL
step** (revised in the post-review pass of section 9; the M2/M3 code passed
`dt_loc` down and is what the paragraph above described). The implicit
operator is unconditionally stable, so nothing ties its step to a
sound-crossing time, while the quantities it relaxes evolve on
`dr^2/(D_12 + K_zz)` and `dr/max(|v|, w_s)` with `w_s = D_12 |G|` the
settling drift speed. At the base of HD 209458 b the first of those is
~10^6 s against a ~1 s CFL step. `relax_element_composition` therefore
starts at the smaller of the two scales, per cell, and grows the step
geometrically (x1.5, capped at 10^12 times the start), so the late steps are
direct steady solves; the coefficients are Picard-updated at every step.
The convergence measure is **absolute**,

```
max_j |X_j^new - X_j^old| / X_base
```

with `X_base` the reservoir helium mass fraction, for both the inner stop
(1e-12 per step) and the outer pass stop (1e-3 over a pass, <= 20 passes).

**The outer loop is damped.** Nothing in a Picard iteration of two solves
keeps them from chasing each other, so the composition update handed back is
under-relaxed, `X <- X_old + omega (X_relaxed - X_old)`, with `omega` from
0.5, halved (floor 0.125) on any pass whose drift failed to fall, and the
blend put back through the same nonnegative projection the operator ends on.
The drift the convergence test reads is the *undamped* distance, so a small
`omega` cannot buy a false convergence. `EXHALE_DIFF_OMEGA` pins the starting
value; `EXHALE_DIFFUSION_CHECK=1` writes one composition and mass-flux
profile per pass to `output/diffusion_pass_profiles.txt`. Section 9.2 records
what those profiles then showed about the `Kzz = 0` case.
The relative measure `max |dX|/X` it replaces is meaningless in a cell the
transport has emptied, which is exactly where the old loop's drift lived.

**The relaxation advects on the face mass fluxes of the state (B5c,
2026-09-07).** The non-conservative form (5) is the conservative equation (1)
only where `rho` and `v` satisfy continuity: the two differ by
`X d(r^2 rho v)/dr / (rho r^2)`. In the converged states that is satisfied in
the wind and violated outright at the base -- measured on HD 209458 b and
LHS 1140 b, the spread of `r^2 rho v` below 1.02 R_p is 10^2 to 10^4 times its
own median, because the base carries a standing sound wave whose sign
alternates from cell to cell. Relaxed to convergence on that field, (5)
converges to the composition of a flow that neither conserves mass nor exists,
which is why `relax_element_composition` used to impose a smoothed
`rho v = mdot_steady/r^2` with `mdot_steady` the median of `r^2 rho v` over
`[j_min:N]`.

That repair belonged to the non-conservative form and went with it. The
relaxation now advects with the CONSERVATIVE term `div(F_rho X)` on the face
mass fluxes `F_rho` of the state itself -- the same faces, the same three
routines (`species_face_fraction`, `species_face_flux`,
`species_flux_divergence`) and the same expression the Runge-Kutta stages
integrate and the stationary elemental row of `element_transport_residual`
balances. `div(F_rho X) = X div(F_rho) + F_rho grad X`: the cell-velocity form
was the second term alone, so it had the mass row's own error subtracted out
and a base that does not conserve mass had to be hidden from it; the
divergence form carries that error, and a cell can lose only the fraction of
its mass the mass row loses. The fixed point of the relaxation is then the zero
of the row that judges it, which is what makes a converged elemental Picard
alternation a state the row reads as stationary. The marching path passes no
flux at all: there the advection is the same divergence, taken inside the
stages. `docs/steady_solver_design.md` section 15.

Two defects on that path are fixed in
this phase because Phase F cannot run without them:

- the direct-steady route (`EXHALE_main.f90:399/401`, the `EXHALE_PTC`
  path that every LHS 1140 b case uses through `finish_case.sh`) has **no
  diffusion loop at all** -- a diffused wind cannot be Newton-finished
  there. The `it_diff` loop moves into one routine called from both
  routes;
- the drift metric hard-codes `f_sp(:,3..5)/f_sp(:,1..2)` (`:951-958`),
  omitting `HeTR` and every molecular species; it becomes the same
  element count as 7.2.

### 7.4 Restart

`load_IC.f90:180-231` rescales every cell to the global `HeH` when the
file's composition differs by > 1e-6, so a diffused state cannot be
restarted -- and `finish_case.sh` restarts from a snapshot. With
`He_diffusion` on, `load_IC` must keep the element split of each cell (the
base cell still pinned to `HeH`); the HeH+ refusal there is lifted by the
same element-count logic. Without this, T8 and Phase F are impossible.

*Done. The element split is kept from M2. The HeH+ part was closed at M4
(2026-08-26): with `He_diffusion` on, the only cells rescaled at all are the
base and its inner ghosts, and those are now projected onto their two element
totals exactly as the operator's write-back does -- `HeH+` by the smaller of
the two factors, the shortfall deposited into the neutral ground species.
Without `He_diffusion` the whole column has to be rescaled and the refusal
stands. The same reading found that the helium nucleus count `load_IC` builds
for the composition-match test omitted `HeTR`; it is now in it.*

### 7.5 Output and tests

- `Ion_species.txt` already lets He/H(r) be reconstructed; add nothing to
  the schema. `exhale_io.py` gains a derived `heh_profile` (nuclei, via
  the same `bsp_nH/nHe` weights) so the tests read one definition.
- `make check` contains no diffusion case (grep over
  `backup/regression/*/input.inp`): today the only guard is the early
  `return` when the flag is off. Phase D adds one regression case with
  `He_diffusion: True`, golden-snapshotted at the end of the series so the
  operator is pinned from then on.

  **Done (2026-08-26): `backup/regression/mol_diffusion`**, the sixth default
  case of `run_check.sh`. It is `mol_base_handoff` -- the Tier-2 hot-Uranus
  gate, molecular chemistry on, its pinned `base.inp` copied unchanged --
  plus `He_diffusion: True` and `He_Kzz: 1.0e9`, a 12000-step relaxation
  snapshot (`maxsteps`) like its siblings. One case rather than the atomic
  HD 209458 b run first planned, because after M4 this single run crosses the
  whole operator: the molecular-carrier closure of section 5 below the front,
  the stage-resolved friction pairs of section 2.6 through it, the projection
  back onto `f_sp`, and the Coulomb-suppressed friction of the ionized region
  above the front. Measured against `mol_base_handoff` at the same step
  count, the two differ in `Hydro_ioniz.txt` and `Ion_species.txt` while
  `du` moves only 2.6292 -> 2.6290 and `log10 Mdot` stays 10.58: the case
  pins the operator's arithmetic, not a large physical effect, because on
  this planet the wind sweeps the composition along faster than diffusion
  separates it (`docs/Update_EXHALE_stage1.md` section 73).

### 7.6 Traps the survey flagged

- Name shadowing: the module's own history (`t0` vs `T0`, 5800x error).
  No new local may reuse `T0 R0 v0 n0 p0 N Ng count info` or any other
  `global_parameters` name, case-insensitively.
- The metal loop's `fXbase` and `rX` cap logic was kept as is at M2; it is
  removed in the post-review pass of section 9 (it is the same limiter the
  memo rejects for helium in section 1, and the same face-velocity advection
  that emptied the base cell of helium emptied it of every metal).
- `x_h2 = 2 q (1+HeH)/(1+q)` at `ionization_equilibrium.f90:1305` uses the
  global `HeH` for the molecular base fraction: consistent with the
  Dirichlet base, left alone.

## 8. Decisions requested before coding

- **D1. Adopt Option A as the architecture**, with the Blanc's-law /
  mean-carrier closure of section 5 as milestone M4, gated by T7; the
  exclusion stays in force until T7 passes. Recommended. **RESOLVED
  (2026-08-26):** adopted and built; T7 passes and the exclusion is gone.
  One term had to be added to the closure as section 5 stated it -- the
  mole-fraction driver (13), which is what makes the carrier bookkeeping
  a force and not just a bookkeeping.
- **D2. Keep the operator-split + outer co-convergence** (`it_diff` loop
  with JFNK) rather than adding diffusion to the JFNK residual, and
  transport `X` in the advective form (5) so that the split is
  well-defined without the hydro's face fluxes; transient helium-mass
  conservation is then to hydro truncation order, exact at the steady
  state. Recommended for this phase; the fourth-conserved-variable
  alternative is the recorded follow-on.
- **D3. Ambipolar field from the electron pressure gradient, (3a)**,
  replacing the hydrogen-plasma constant; it reduces to the validated
  P2b value in the He-trace limit and gives `5/3` in the He++ limit.
  Recommended. **RESOLVED, and the friction with it (2026-08-25).** The
  field was adopted at M2. The neutral-everywhere `D_12` that D3 left
  standing as a caveat is now replaced by the stage-resolved friction of
  section 2.6: hard sphere for neutral-neutral, polarization for
  non-resonant ion-neutral, Coulomb for ion-ion, combined by the
  stage-fraction harmonic sum, with the ion-neutral pair carrying the
  polarization and rigid-core channels added as frictions (9b). There is no
  key -- it is a correction to
  the physics, not an option -- and the same treatment runs in the trace
  metal loop, which also stops assuming the hydrogen-plasma `eE = m_H g/2`
  and now reads the computed field. Measured in section 9.3.
- **D4. Metals**: slaved to hydrogen as component 1 (section 2.1), so the
  mass-fraction algebra closes with `eos_metals` on or off; the
  `He_metal_diffusion` element loop stays as is (trace, validated). T4
  runs with both `eos_metals` settings. Recommended.
- **D5. Input keys**: `He_diffusion`, `He_ambipolar`, `He_alphaT`,
  `He_metal_diffusion`, `He_Kzz` keep their names and meaning. The base
  cap has no key and is removed. **One user-visible change:** the
  `He_Kzz` default goes from `1e9` to `0` (7.2); an input relying on the
  old default must state it. A radius-dependent `K_zz` profile is not
  added in this phase (T7 uses the constant). [implemented 2026-08-27 in
  Phase E: a `Lower atmosphere profile:` file supplies `K_zz(p)`, which
  `src/modules/files_IO/lower_atmosphere_profile.f90` interpolates onto the
  grid as `kzz_cell` (`parameters.f90`); `He_Kzz` remains the constant used
  when no profile is given.]
- **D6. The direct-steady route and restart** (7.3, 7.4) are in scope:
  both are prerequisites for running the LHS 1140 b cases with diffusion,
  and both are defects of the present code independent of the
  formulation. Recommended: yes.
- **D7. Baseline capture first.** Before any code change, re-run
  `examples/14_diffusion` with the present binary to a Newton finish and
  keep that output as the T2 reference (section 6). Recommended: yes; it
  is the only way T2 can mean anything.

## 9. Review record (2026-08-25)

`binary_diffusion_design_review.md` (external, read-only against the memo
and the code) raised eight findings. Disposition, each checked against the
memo text and, where it applies, the code:

| # | finding | verdict | done |
|---|---|---|---|
| 1 | settling sign reversed in the face flux of section 3 | **correct** -- (2)-(3) give `-B`, the draft wrote `+B` | (6) rewritten with `-B_f`; sign check added to 2.3 |
| 2 | conservative `rho X` advection is not consistent with the hydro's continuity step through the present interface, so T6 cannot be guaranteed | **correct** on the interface (the routine gets only new `rho, v, T, dt`) | advective form (5) adopted; T6 exact by construction, transient conservation stated as hydro-order, steady state exact (T8); alternative recorded |
| 3 | `X = rho_He/rho` with `1 - X - X_met` for hydrogen breaks the two-component identities when metals are in `rho` | **correct** (`eos_metals` default on) | components redefined: hydrogen + slaved metals vs helium, `m_1` from `comp_mass_per_H`; T4 in both `eos_metals` settings |
| 4 | endpoint states are not invariant: the gradient coefficient is finite at `X = 0, 1`; positivity needs an M-matrix argument for the whole operator | **correct** | claim removed (2.4); M-matrix argument for advection + gradient + hybrid settling, and for `1 - X`, written out in section 3 |
| 5 | the P2b ambipolar constant is the hydrogen-plasma value and cannot give `2.5` in the He-rich limit | **correct**; derived here: He++ plasma gives `eE = 4/3 m_H g`, relative settling mass `5/3` | field computed from `-(1/n_e) dp_e/dr` (3a); T9 rewritten with three analytic limits |
| 6 | a mean carrier mass in the He-H formula is not a derived molecular closure | **correct** | Blanc's-law friction + carrier-based thermodynamic/force terms specified; molecular region staged as M4 behind T7; exclusion lifted only then |
| 7 | HeH+ write-back could drive `HI` negative | **correct** | nonnegative projection (scale by element, two-element species by the smaller factor, shortfall deposited into the neutral ground species) |
| 8 | T1 conflicts Dirichlet with conservation; T2 gates on a capped, non-asymptotic reference; T7 assumes a profile the scope lacks; top BC stated two ways | **correct** on all four | T1a/T1b split; T2a convergence sequence + T2b report; T7 as a constant-`K_zz` series; top BC stated once in section 4 |

Two review remarks were not adopted: restricting the first implementation
to a steady-state-only formulation (the advective form keeps the transient
path usable and is exact at the steady state, which is what D2 needs), and
deferring the ambipolar field (it is cheap to compute from the solved
state and its omission would make T3 wrong by construction). The review's
survey assessment agreed with section 7 on every point.

### 9.1 Post-review findings (2026-08-25, after M3)

Three defects found on the converged HD 209458 b runs of milestone M3
(`backup/phase_d_baseline/new_kzz0`), where the helium ratio of the first
free cell above the base sat 550x below the reservoir while the barometric
He-H separation scale there is `H = kT/(dm g) = 0.0045 R_p`, 24 times the
cell width -- so a 550x drop across one cell is not barometric physics.

**(a) The advection was built from face-averaged velocities.** With
`v_f = 0.5 (v_j + v_{j+1})` and an upwind selection per face, a cell whose
two face velocities straddle zero (`v_f(j-1) < 0 < v_f(j)`) donates at both
faces and receives at neither: its advective coefficients vanish
identically, even though its own `v_j` is large. That is the standing
configuration of the breathing base, where the base sound wave alternates
the sign of `v` from cell to cell -- measured at the first free cell of
`new_kzz0`: `v = (-796, +63, +28, -8, -9, +5) cm/s`, giving face velocities
`-366` and `+45`. The cell was then left with only the molecular diffusion
time `dr^2/D_12 = 7e5 s` against a `~2 s` step, so whatever composition it
held was frozen there. The advection is now the one-sided upwind difference
on the **cell** velocity, which keeps the cell coupled to its donor -- here
the Dirichlet base itself, since `v_2 > 0`. The transport coefficient of
that coupling, `rho v/dr`, is 200 times the diffusive one at this cell.
Regression: T11.

**(b) The trace-metal solver transported the density conservatively with the
same face velocities.** In the conservative form the configuration above is
not merely a decoupling but a drain: the cell has an outflow term at both
faces and an inflow term at neither. It has no counterpart in the hydro's
own `rho`, whose face fluxes come from the Riemann solver and carry no such
divergence. Measured on `new_kzz0`: every metal emptied from the same cell
by ~10^3 (`C I` 2.3e7 against 2.6e10 in its neighbours), and with the metals
the metal-line cooling of that cell. The solver now transports the **mixing
ratio** `f_X = n_X/n_H` in the same advective form as helium, so a uniform
`f_X` is preserved under any velocity field and no spurious divergence can
create or destroy the element. The `f_X <= f_X(base)` cap in the write-back
went with it: settling piles an element up as readily as it depletes one,
and the cap is the limiter section 1 rejects for helium.

**(c) The relaxation was driven at the hydro CFL step, and on a velocity
field with no steady mass flux.** Both are recorded in section 7.3, where the
coupling to the steady solver is specified: the step size is now a
composition time scale grown geometrically, the convergence measure is
absolute, and the advecting flow of the relaxation is the face mass flux of
the state (it was the wind's smoothed steady mass flux until item B5c; the
measurements quoted below were made with that form). Measured effect of (c) alone on the HD 209458 b `Kzz = 0` wind, solving
the operator's steady state offline on the M3 converged state: with the cell
velocities the base takes a five-fold cliff at the fourth cell and a plateau
at 0.10 of the reservoir ratio (where the barometric separation scale is 24
cells); with the steady flux the first free cells are flat to 1% and the
profile leaves the base on the barometric scale (0.97 at 1.05 R_p, 0.90 at
1.10 R_p, 0.85 at 1.20 R_p).

### 9.2 What the damped loop showed (2026-08-25)

The under-relaxation of 7.3 was added to close the `Kzz = 0` HD 209458 b case,
whose outer loop had limit-cycled through all twenty passes. It does not close
it, and the profiles written each pass say why: **the cycle is in the wind, not in the
composition.** Over passes 3-20 the median `r^2 rho v` in the escape window
alternates between 2.1-2.2e-7 and 5.9-6.8e-7 (a factor 2.8 in `Mdot`,
`log10 Mdot` 9.867 against 10.376) while each branch is flat to 1e-3 and every
JFNK solve returns `info = 0` with `||R||` between 9.0e-5 and 9.6e-4, inside
the 1e-3 target. Across the same two states the element ratio differs by at
most 0.0117 above 1.05 R_p, and each parity is stationary to 0.003 from one
occurrence to the next. Damping the composition therefore cannot help, and the
run demonstrates it: `omega` fell to its 0.125 floor by pass 11 and the cycle
held at full amplitude. Started at 0.5 instead, the pass-3 wind solve stagnates
(`info = 2`, line search collapsing at r = 1.009 on the energy row) and the
loop exits on its solver-failure guard.

Two discrete steady states 2.8x apart in mass flux that the residual test
accepts equally is a property of the steady solver at this configuration, not
of the transport operator; the `Kzz = 1e9` case with the same operator
converges to 1.18e-4 in ten passes. It is left open here.

### 9.3 The friction the ionized wind actually has (2026-08-25, decision D3)

The friction is now stage-resolved (section 2.6). What follows is measured,
not estimated.

**The coefficient.** `EXHALE_DIFFUSION_CHECK=1` writes `D_eff`, the neutral
hard-sphere coefficient of the same cell and their ratio into
`output/element_flux_profile.txt`, with the stage pair carrying the largest share of
the friction. On the HD 209458 b `Kzz = 1e9` wind and the LHS 1140 b
`He/H = 0.55` wind:

| r [R_p] | HD 209458 b `D_eff/D_neutral` | dominant pair | LHS 1140 b `D_eff/D_neutral` | dominant pair |
|---|---|---|---|---|
| 1.02 | 1.000 | HeI-HI | 1.000 | HeI-HI |
| 1.05 | 1.000 | HeI-HI | 0.988 | HeI-HI |
| 1.10 | 0.918 | HeI-HI | 0.744 | HeI-HI |
| 1.20 | 0.161 | HeII-HII | 0.163 | HeII-HII |
| 1.50 | 9.56e-3 | HeII-HII | 4.13e-2 | HeII-HII |
| 2.00 | 2.13e-3 | HeII-HII | 1.33e-2 | HeII-HII |
| 3.00 | 3.08e-4 | HeII-HII | 2.87e-3 | HeII-HII |
| 4.00 | 7.59e-5 | HeIII-HII | 1.02e-3 | HeII-HII |

(HD 209458 b from the final `new_kzz1e9_d3b` state, LHS 1140 b from
`heh0p55_diff_d3`, which the ion-neutral change does not reach -- see the
end of this section.)

The transition is the ionization front, and above it the suppression grows
with radius because the Coulomb coefficient scales as `T^(5/2)/n` while the
hard sphere scales as `T^(1/2)/n`: the falling temperature of the outer wind
weakens the ionized friction far more slowly than it weakens the neutral
one. At 4 R_p helium is held to the protons `1.3e4` times more strongly than
the neutral coefficient claimed.

**Limits (T12).** `D_eff` reproduces the Banks & Kockarts hard sphere in an
all-neutral column to `1.2e-16` relative, the He++/H+ Coulomb coefficient in
a fully ionized one to round-off, and the four-term Blanc average of a
half-ionized one to round-off, each against a closed form written
independently in `src/tests/diffusion_tests.f90`. The combined ion-neutral
coefficient (9b) sits within `2.5e-2` of the polarization channel at 1 K and
within `1.2e-3` of the hard sphere at 1e9 K, and never exceeds either.
`make diffusion_tests && ./diffusion_tests.x`: 21 passed, 0 failed.
`make check`: 5/5 byte-identical (the operator is entered only with
`He_diffusion` on, which no regression case sets).

**HD 209458 b** (`backup/phase_d_baseline/new_kzz1e9_d3` against
`new_kzz1e9_ctrl`). The control is *not* the stored `new_kzz1e9_fix`
directory: that run was made with a binary older than the damped outer loop
of 9.2, so it differs from the new run in more than the friction.
`new_kzz1e9_ctrl` is the same configuration re-run from the same converged
state with a binary built from the commit this change sits on, everything
except the friction identical. Both Newton-finished on the direct-steady
route, `info = 0`, `||R||` 8.1e-4, 7.0e-4 and 7.3e-4; all three left the
outer loop on the drift criterion, in 4, 8 and 7 passes.

`new_kzz1e9_d3` is the intermediate state in which the ion-neutral pair
carried the polarization channel alone; `new_kzz1e9_d3b` is the adopted
form with (9b). Both are kept because the difference between them is the
whole content of the ion-neutral decision.

| quantity | pre-D3 control | D3, pol. only | **D3, (9b)** |
|---|---|---|---|
| `log10 Mdot` [g/s] | 10.0106 | 10.0656 | **10.0561** |
| He 10830 red depth [%] | 8.2044 | 7.8094 | **7.8542** |
| He 10830 blue depth [%] | 1.2671 | 1.1980 | **1.2065** |
| `(He/H)/HeH` at 1.10 R_p | 0.9282 | 0.9238 | **0.9219** |
| at 1.20 R_p | 0.8991 | 0.8526 | **0.8502** |
| at 2.00 R_p | 0.8989 | 0.8425 | **0.8403** |
| at 4.00 R_p | 0.8327 | 0.8419 | **0.8391** |

and the metal mixing ratios, each normalized to its own base cell:

| element | ctrl 1.10 | 4.00 | | pol. only 1.10 | 4.00 | | **(9b)** 1.10 | 4.00 |
|---|---|---|---|---|---|---|---|---|
| C | 0.675 | 0.438 | | 0.588 | 0.524 | | **0.628** | **0.563** |
| O | 0.504 | 0.259 | | 0.547 | 0.311 | | **0.536** | **0.301** |
| Mg | 0.124 | 0.032 | | 0.264 | 0.240 | | **0.322** | **0.296** |
| Na | 0.165 | 0.047 | | 0.199 | 0.110 | | **0.199** | **0.110** |
| Fe | 0.000 | 0.000 | | 0.034 | 0.028 | | **0.131** | **0.107** |
| Ca | 0.000 | 0.000 | | 0.0145 | 0.0085 | | **0.0321** | **0.0193** |

Two things are visible and they are the two halves of section 2.6.

1. **Above the front the mixing ratios are flat.** With the Coulomb
   coefficient the diffusive term is negligible against advection, so every
   element is carried at whatever ratio it had when it crossed the front,
   which is what a steady wind must do. The pre-D3 profiles kept settling all
   the way out -- iron and calcium to *zero*, magnesium to 2.5% of the base
   ratio at 4 R_p -- because the neutral coefficient never turned the
   settling off. Iron and calcium are the elements this matters most for:
   they went from completely removed to depleted-but-present, which is the
   difference between having and not having Ca II H&K and Fe II in the
   transmission spectrum.
2. **The plateau each element freezes into is set by the ion-neutral pair,
   and only for the elements that are already ionized where hydrogen is
   not.** Between 1.0 and 1.2 R_p hydrogen and helium are still neutral --
   the measured `x(H+)` is 0.000, 0.004 and 0.072 at 1.05, 1.10 and
   1.20 R_p -- so the He-H friction there is neutral-neutral and the
   ion-neutral form barely enters it: `D_eff` at 1.10 R_p changes by 1%
   between the two columns and the helium profile does not recover
   (0.8425 -> 0.8403 at 2 R_p; the control is 0.8989). The metals are a
   different case, because the low first ionization potentials put them in
   ionized stages while hydrogen around them is still neutral: at 1.10 R_p
   the measured ionized fractions are Fe 0.999, Mg 0.943, Ca 0.269,
   C 0.329, Na 0.053, O 0.003. Going from the polarization channel alone to
   (9b) raises the frozen plateau in exactly that order -- Fe x3.8,
   Ca x2.3, Mg x1.23, C x1.07, Na x1.00, O x0.97 -- which is the check
   that the change acts where the physics says it should. Adopting (9b)
   therefore recovers between a quarter and a third of the way back to the
   control for the low-potential metals, and leaves helium, sodium and
   oxygen where they were.

The mass-loss rate moves **+0.046 dex**. It moves *up*, not back toward the
pre-Phase-D 9.79: the total radiative cooling of the new state is *higher*
(`3.59e-7` against `3.27e-7` in the volume integral, the metal channels
`3.09e-7` against `2.79e-7`, since more metal is left aloft), and the wind is
nonetheless hotter -- by up to 3.3% at 1.11 R_p -- and faster. The rate is a
property of the whole coupled state, and no single-channel account of the
0.046 dex is offered here.

**LHS 1140 b** was not re-run for (9b). Its comparison below is
`heh0p55_diff_d3` (polarization only) against the control, and it already
shows the wind untouched to 9e-16 in temperature: the separation happens in
the neutral first cell, where neither the Coulomb nor the ion-neutral
channel is reached. (9b) changes only ion-neutral pairs, so it cannot move
a case that the larger Coulomb change did not move. That last step is
argued, not measured -- the only claim in this section that is.

(`LHS1140b/exhale/heh0p55_diff_d3` against
`heh0p55_diff_ctrl`, built the same way, same restart-and-finish,
`info = 0`, `||R||` 9.96e-4 both): `log10 Mdot` 7.6065 both, He 10830 red
depth 9.969e-5 against 9.988e-5 % (+0.2%), the He/H profile different by at
most 1.4e-4 relative and the temperature by 8.7e-16 -- the wind itself is
untouched. The composition relaxation exited on its first pass with a drift of
1.9e-4, i.e. the old converged state is a fixed point of the new operator as
well. The reason is in the table above: this wind separates its helium
*inside the first cell above the base*, where the gas is neutral and
`D_eff = D_neutral` by construction -- `(He/H)/HeH` is already 0.72 at
1.0013 R_p and 0.17 at 1.0074 R_p -- and above that there is no helium left
for the Coulomb suppression to hold. The change is a correction to the
ionized wind, and this case does its separation in the neutral base.

### 9.4 The eddy coefficient LHS 1140 b needs (2026-08-25)

The `Kzz = 0` LHS 1140 b case of 9.3 is the reason the open item "those cases
need a `He_Kzz` argued from the planet, not the default 0" was raised. It is
now answered with a literature survey and a seven-point scan, both written up
outside this memo: **`LHS1140b/kzz_decision.md`** (the proposal, with the
recommendation and the alternatives) and **`LHS1140b/kzz_literature.md`** (the
survey, with the verification status of every quoted number). The scan runs
are `LHS1140b/exhale/heh0p55_diff_kzz{1e6,...,1e11}` against
`heh0p55_diff_ctrl`, all `info = 0` with composition drift under 1.3e-3, and
the numbers come out of `LHS1140b/exhale/kzz_scan_table.py`. In summary:
`D_eff` at this wind base is 4.6e5 to 1.3e6 cm^2/s but rises to 2.5e10 by
1.2 R_p, so `K_zz` above ~1e6 lifts the homopause off the base while
`K_zz >= 1e10` is needed before the element ratio is back at the reservoir
value at 1.05 R_p, and no value makes the outer wind well mixed --
`(He/H)/HeH` freezes at 0.43 above 5 R_p even at 1e11. The He 10830 red-pair
equivalent width therefore rises from 2.5e-5 to 0.605 %A between `K_zz` = 0
and 1e10 and then saturates, against 1.109 %A with diffusion off at the same
reservoir composition, while `log10 Mdot` moves only 7.61 to 7.84 over the
whole range. The recommendation on the table is `He_Kzz: 1.0e9`, following
Taylor et al. (2025) at the same 1e-6 bar boundary and matching
`examples/14_diffusion`; the value is a user decision. Two solver facts
belong here: restarting a large-`K_zz` case from the `Kzz = 0` state stalls
the line search, so the scan was built one decade at a time; and on this wind
JFNK settles at `||R||` of 2.7e-3 to 4.7e-3 rather than 1e-3, so the upper
rows state a looser `Resid tol` and carry their achieved `||R||`.

**Adopted 2026-08-25 (user decision): `He_Kzz = 1.0e9` for the LHS 1140 b
runs.** With the operator active at that value the He 10830 equivalent width
crosses the measurement at `He/H = 2.09` rather than the 0.55 of the
diffusion-off calibration (scan in `LHS1140b/kzz_decision.md` section 6).
