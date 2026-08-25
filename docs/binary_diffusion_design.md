# Design: binary H/He diffusion for a helium-rich wind (Phase D)

Status: **design memo, not implemented.** Written 2026-08-25 for Phase D of
`lhs1140b_lower_atmosphere_plan_new.md`; revised the same day against the
external review `binary_diffusion_design_review.md` (section 9 records
each finding and what was done with it). Implementation waits for the
user's go-ahead; the open decisions are listed in section 8.

Supersedes the transport formulation of `design_hehe_diffusion.md`
(sections 2-5), whose kernel is the *trace-helium-in-hydrogen* limit of what
is written here. That document remains the record of the Phase-1/2 build
and validation (its sections 7b-7e).

## 1. Why the present kernel cannot be extended

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
3 n_H3+ + n_HeH+`, `n_He = n_HeI + n_HeII + n_HeIII + n_He(2^3S) + n_HeH+`
(the `bsp_nH`, `bsp_nHe` weights of `species_table.f90`).

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
finding 1). The settling part keeps the Peclet-based central/upwind
hybrid of P2 (central where `|B_f| dr < 2 (A_f + E_f)`, upwind toward the
direction of the settling velocity otherwise).

**Bounds `0 <= X <= 1`, argued for the whole operator.** Coefficients are
frozen within a step at the Picard level, so the step is linear in
`X^{new}`. For each interior row the backward-Euler matrix has diagonal
`1/dt + (outflow face terms) > 0`, off-diagonals that are `<= 0` for the
upwind advection, for the gradient term (`A_f + E_f >= 0`), and for the
settling term whenever the hybrid switches to upwind, while the central
settling branch is used only when its contribution cannot flip the sign
of an off-diagonal (that is the Peclet condition above). The matrix is
therefore an M-matrix: its inverse is nonnegative, and with `X^{old} >=
0` and boundary values `>= 0` the update gives `X^{new} >= 0`. The same
operator applied to `1 - X` (the hydrogen equation is (5) with `-J`, and
`-J` has the identical structure with the settling direction reversed --
the binary symmetry of (2)) gives `1 - X^{new} >= 0`. Hence the bounds
hold per step without clipping; a clip at round-off (`1e-15`) is kept
only as an assertion, not as a limiter.

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
   invariant, and the bounds on `X` come from the M-matrix argument
   below, not from vanishing fluxes. The prefactors are evaluated at the
   face at the new time level through one Picard sweep.
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
- `HeH+` is counted as a helium carrier for the friction (its abundance
  is negligible wherever the closure matters) and handled in the
  write-back as in section 3.

**Staging.** The core of this phase (milestones M1-M3 in section 6) is
built and validated in the *atomic* region, where the closure is the
plain binary one and no approximation beyond `D_12` enters. The molecular
closure is milestone M4: implemented behind the same key, gated by T7,
and the `error stop` is lifted only when T7 passes. Until then the pair
remains refused -- the exclusion is removed by a validated closure, not
by deleting the check.

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
| T0 | regression with `He_diffusion` off | 5/5 byte-identical |
| T1a | **diffusive equilibrium**: `v = 0`, isothermal column, no eddy, zero diffusive flux at *both* ends; start uniform | `X(r)` relaxes to the barometric separation `dx/dr = -x(1-x) (m_He - m_1) g/kT` (neutral) to grid order; total He mass constant to round-off (the closed-column version of T4) |
| T1b | same with the Dirichlet reservoir base of section 4 | the integrated base flux `int 4 pi r_b^2 J_b dt` accounts for the change of total He mass to round-off |
| T2a | **convergence to the trace equation**: `examples/14_diffusion` configuration (HD 209458 b) at He/H = 0.0833, 1e-2, 1e-3, 1e-4, run with the new operator and with the present kernel | the difference between the two (He/H)/HeH profiles decreases with He/H at the expected first order in `x`; at 1e-4 it is below 0.1% |
| T2b | **report** at the tutorial abundance He/H = 0.0833 (not a gate: the present profile is capped over 1.2-2 R_p and the abundance is not asymptotically small) | (He/H)/HeH profile and He 10830 depth tabulated against a **fresh, Newton-finished run of the present code made before any Phase-D change** -- the stored `examples/14` output is unconverged (JFNK info=2, du 1.3e-2) and in an older column schema, and the tables in `design_hehe_diffusion.md` 7d/7e are from other runs; neither is a usable reference |
| T3 | **H-trace limit**: same planet with He/H = 1000 | finite everywhere, `X` in `[0,1]`, hydrogen rises (X decreasing outward relative to the reservoir), no cap active |
| T4 | **elemental conservation, closed column**: zero flux at both ends, `v = 0`, arbitrary initial `X(r)` | `int rho X r^2 dr` constant to 1e-12 relative over 1e4 steps |
| T5 | **zero net diffusive mass flux at every face**: instrument `J_He + J_H` | identically zero (by construction; the test guards the discretization) |
| T6 | **uniform mixture preserved**: uniform `X`, `g = 0`, any `v(r)` | `X` stays uniform to round-off (advection of a uniform field is continuity) |
| T7 | **homopause (molecular + eddy)**, milestone M4: molecular chemistry on + `He_diffusion` on, a series of constant `K_zz` values | runs (no `error stop`); the He/H profile is flat where `K_zz >> D_{He,1}` and separates above; the homopause radius moves with `K_zz` as `D_{He,1}(r_h) = K_zz`; the Blanc's-law coefficient reproduces `D_{He,H}` in a forced-atomic run and `D_{He,H2}` in a forced-molecular one |
| T8 | **moving-wind face flux**: converged LHS 1140 b wind, diffusion on | elemental face flux `4 pi r^2 (rho X v + J)` constant with radius above the base to the same tolerance as the mass flux (`du`) |
| T9 | **ambipolar limits** from (3a) on prescribed states: neutral; fully ionized H+ plasma; fully ionized He++ plasma; a partially ionized front | relative settling mass `3`, `2.5`, `5/3` to round-off in the three limits; finite and continuous through the front |
| T10 | grid and time-step convergence on T1 | second-order in `dr` for the gradient term, first-order where the settling hybrid upwinds |

T2a is the gate that protects what already works (the old kernel is the
comparison, not the truth standard); T3 is the one that shows the new
formulation does what the old one could not.

**Milestones.** M1 baseline capture (D7) and the two prerequisite defects
of 7.3-7.4. M2 the atomic binary operator with (3a), tests T0-T6, T8-T10.
M3 the LHS 1140 b He-rich runs with diffusion (the Phase F precursor,
Newton-finished through the repaired direct-steady route). M4 the
molecular closure and T7; only then the exclusion is lifted.

*Status (2026-08-25): M1, M2 and M3 complete; M4 open. The three defects of
7.3-7.4 are fixed (`docs/Update_EXHALE.md` section 67); the atomic binary
operator is `src/modules/functions/binary_element_diffusion.f90`, with T0 and
T1a-T6, T9, T10 passing in `make diffusion_tests && ./diffusion_tests.x`
(section 68) and the T2a sequence measured there (3.05e-6 at He/H = 1e-4,
first order in the helium fraction). M3 (section 69): the baseline capture of
D7 is complete in `backup/phase_d_baseline/`, T2b is tabulated against it,
**T8 passes** (elemental face flux as flat as the mass flux in the code's own
escape window, 5.1e-3 against 7.2e-3), the three LHS 1140 b helium-rich cases
run to `info = 0` with diffusion on, and the trace kernel
`species_diffusion.f90` is deleted -- so T2a stands as a recorded measurement
and is no longer runnable. T7 and the molecular closure are M4, and the
`He_diffusion` regression case of 7.5 is still to be added.*

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

**The relaxation advects with the steady mass flux.** The advective form (5)
is the conservative equation (1) only where `rho` and `v` satisfy continuity:
the two differ by `X d(r^2 rho v)/dr / (rho r^2)`. In the converged states
that is satisfied in the wind and violated outright at the base -- measured
on HD 209458 b and LHS 1140 b, the spread of `r^2 rho v` below 1.02 R_p is
10^2 to 10^4 times its own median, because the base carries a standing sound
wave whose sign alternates from cell to cell. Relaxed to convergence on that
field, (5) converges to the composition of a flow that neither conserves mass
nor exists. `relax_element_composition` therefore imposes

```
rho v = mdot_steady / r^2,   mdot_steady = median of r^2 rho v over [j_min:N]
```

-- the wind's own mass flux, on the window the solver itself uses to declare
the wind steady -- and the base then relaxes to the barometric profile.
The marching path keeps the cell values `rho_j v_j`: there the wind is
genuinely transient and the cell velocity is the consistent choice.

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

### 7.5 Output and tests

- `Ion_species.txt` already lets He/H(r) be reconstructed; add nothing to
  the schema. `exhale_io.py` gains a derived `heh_profile` (nuclei, via
  the same `bsp_nH/nHe` weights) so the tests read one definition.
- `make check` contains no diffusion case (grep over
  `backup/regression/*/input.inp`): today the only guard is the early
  `return` when the flag is off. Phase D adds one regression case with
  `He_diffusion: True` (the HD 209458 b tutorial size, Newton-finished),
  golden-snapshotted at the end of the series so the He-trace limit is
  pinned from then on.

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
  exclusion stays in force until T7 passes. Recommended.
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
  Recommended. The neutral Banks & Kockarts `D_12` is kept for the
  friction with its caveat recorded at the coefficient (a Coulomb ion-ion
  coefficient changes *where* separation turns on, not the formulation);
  revisit after T3.
- **D4. Metals**: slaved to hydrogen as component 1 (section 2.1), so the
  mass-fraction algebra closes with `eos_metals` on or off; the
  `He_metal_diffusion` element loop stays as is (trace, validated). T4
  runs with both `eos_metals` settings. Recommended.
- **D5. Input keys**: `He_diffusion`, `He_ambipolar`, `He_alphaT`,
  `He_metal_diffusion`, `He_Kzz` keep their names and meaning. The base
  cap has no key and is removed. **One user-visible change:** the
  `He_Kzz` default goes from `1e9` to `0` (7.2); an input relying on the
  old default must state it. A radius-dependent `K_zz` profile is not
  added in this phase (T7 uses the constant).
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
absolute, and the advecting flow of the relaxation is the wind's steady mass
flux. Measured effect of (c) alone on the HD 209458 b `Kzz = 0` wind, solving
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
