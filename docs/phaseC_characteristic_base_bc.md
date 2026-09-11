# Phase C: a characteristic lower boundary at the face

Design note and baseline measurements. 2026-09-03.

**Nothing here is applied to the tree.** The measurements were made with a
scratch build; the design in sections 4 and 5 is a proposal, and section 6
lists the decisions that have to be taken before any of it is written, because
they change the mathematical boundary-value problem the code solves.

> **Sources checked against the published originals, 2026-09-07.** The three
> papers this note rests on were read in their journal versions
> (Thompson 1987, *J. Comput. Phys.* **68**, 1; Thompson 1990, *J. Comput.
> Phys.* **89**, 439; Poinsot and Lele 1992, *J. Comput. Phys.* **101**, 104),
> and section 4.1 now carries their equation numbers instead of the secondary
> sources it used to quote. Three outcomes, in the order that matters for the
> decisions of section 6:
>
> 1. **The characteristic count of section 4.3 is confirmed unchanged**, branch
>    by branch, against Thompson (1990) sections 3.2.4 to 3.2.7 and Poinsot and
>    Lele (1992) Tables III and IV. So is the relation `(C-)`, and so is the
>    two-condition inflow closure of section 4.2. No constant was found wrong.
> 2. **The claim that pair (D) is singular at `v_b -> 0` is withdrawn.** The
>    published constant mass flux condition, Thompson (1990) Eqs. (69) and
>    (70), is regular at `u_1 = 0` at a lower boundary. The singularity belongs
>    to the particular closure this note wrote for it, not to the condition.
>    Section 4.2 carries the correction, and it reopens (D) as a candidate for
>    decision D1.
> 3. **The impossibility argument in D6 is too strong.** Two published forms
>    are transparent to sound and are algebraic functions of the instantaneous
>    face state, so a steady residual can express both: Thompson (1990)
>    Eq. (64), the nonreflecting condition in a gravitational field, and
>    Poinsot and Lele (1992) Eq. (40), the partially reflecting relaxation.
>    Neither changes the number of conditions. Section 6 D6 carries the
>    correction; the measurement that a pinned static quantity reflects
>    (section 7.3) is untouched by it.
>
> A fourth outcome is not a correction but a gap this note never stated: the
> source terms. Section 4.1 now says where gravity and the spherical divergence
> enter the characteristic relations, and why they are absent from `(C-)` as it
> is written.

## 1. Scope, build, and what is measured

The task is item 4 of `docs/open_defects_20260903.md` and Phase C of
`docs/open_defects_20260903_review.md` section 6: replace the component-wise
ghost-cell closure at the lower boundary with a condition imposed at the face
`r_edg(0)`, with the characteristic count stated explicitly.

**Build.** The sources are a copy of the `p_ab` working tree
(`.../scratchpad/p_ab/src`), i.e. the tree of 2026-09-03 plus the (AB) change
in which `Apply_BC` writes the ghosts only and `eval_residual` fills them after
its ionization sweep (`.../scratchpad/p_ab/docs/p_ab_ghost_order.md`). That
build is the right baseline for this work: (AB) moves the *steady root's* cell-1
velocity from `-15` to `-252 cm/s`, i.e. onto the marching attractor's own
`-251 cm/s`, and cuts the 3000-step departure by 44 percent at `r = 1.03` and
73 percent at `1.005`. What (AB) leaves behind is what a characteristic
boundary has to remove, so the baselines below are measured on it.

For the specific probe used here the (AB) change is inert: the probe runs an
atomic gas, where the interior round trip `Apply_BC` used to perform is a
last-bit operation, and the measured numbers are the same on the pre-(AB) tree
to every digit printed (checked on the test A table of section 3.1).

**Probe.** `src/tests/base_bc_probe.f90` (scratch only). It links the
production `define_grid`, `set_gravity_grid`, `Apply_BC`, `PLM_rec`,
`Reconstruction`, `Num_Fluxes`, `Source` and `RK_rhs` verbatim and drives them
with the SSP-RK3 update of `EXHALE_main`, with the radiation and chemistry
source terms absent. What it measures is therefore exactly the
Euler + gravity + boundary discretization, with nothing else able to hide or
imitate a boundary error. Every option is an environment variable
(`BASEBC_*`), so no production path is touched. One gate was added to the
scratch `Apply_BC` (`probe_base_reference`, default false) so the probe can
substitute a reference boundary; see section 3.1.

**Background state.** Isothermal hydrostatic, `b0 = 39.87` (the hot-Uranus gate:
`0.0457 M_J`, `0.49 R_J`, `T_0 = 1140 K`, `mu = 2.27 m_H`), `Mixed` grid,
`Base grid [dr,cells]: 2.0e-4 50`, `Grid cells: 500`, HLLC, WENO3 unless
stated. The realized geometry is `r(0) = 1.000000000`,
`r_edg(0) = 1.000095887`, `r(1) = 1.000191773`, `dr_j(1) = 1.9177e-4 R_p`,
`H/dr = 130.8`. Cell averages of the analytic profile are used as the initial
state (8-point Gauss-Legendre, volume weighted), so the initial data is the
exact discrete representation of the atmosphere and the reconstruction
reproduces the analytic face value to its own order.

Velocities are quoted in code units `v0 = sqrt(k T_0/mu) = 2.03602e5 cm/s`.

## 2. The present boundary, exactly

### 2.1 What is written, and where

`BC_component_constrho(W_in, index)` (`src/modules/states/Apply_BC.f90`) writes
three primitive components at one index:

```fortran
W_in(1,index) = rho_bc
if (base_v_massflux .and. base_flux_const .gt. 0.0d0) then
   W_in(2,index) = base_flux_const/(rho_bc*r(index)**2)
else if (valve_eps .gt. 0.0d0) then
   W_in(2,index) = 0.5d0*(W_in(2,1) + sqrt(W_in(2,1)**2 + valve_eps**2))
else
   W_in(2,index) = max(W_in(2,1),0.0)
endif
if (hydrostatic_base) then
   W_in(3,index) = W_in(3,1) + (W_in(3,2) - W_in(3,1))/(r(2) - r(1))*(r(index) - r(1))
else if (base_ghost_T_continuous) then
   W_in(3,index) = (ntot_bc + dp_bc)*W_in(3,1)/n_part_cell1
else
   W_in(3,index) = ntot_bc + dp_bc
endif
```

It is called from two places, with two different meanings and the same code:

* `Apply_BC_W` calls it with `index = 0, -1`, and there `W` holds **cell
  averages**. The ghost written is a cell average at the ghost centre `r(0)`.
* `Rec_BC` calls it with `index = 0, -1` on `WL_out`, and with `index = -1` on
  `WR_out`. There `WL_out(:,j)` is the **face state** at `r_edg(j)`. So
  `WL_out(:,0)` is the left state of the base face, and it is filled with the
  same three expressions.

`WR_out(:,0)`, the right state of the base face, is *not* overwritten: it is
the interior reconstruction of cell 1 down to the face. `RK_rhs` then calls
`Num_flux(WL(:,0), WR(:,0), ...)` and the HLLC solver mediates the two.

So the boundary is already imposed at the face, through a Riemann problem with
a fully specified left state. What is wrong is not the location of the
*condition*; it is that the left state is written component by component from
quantities that belong somewhere else.

### 2.2 Characteristic count

In this code inflow at the base is `v > 0` (the radial coordinate increases
outward, so gas entering the domain moves outward). At `r_edg(0)` the three
eigenvalues are `v - c`, `v`, `v + c`, and a wave enters the domain when its
eigenvalue is positive.

| base state | `v-c` | `v` | `v+c` | entering | reservoir conditions | interior relations |
|---|---|---|---|---|---|---|
| subsonic inflow `0 < v < c` | out | in | in | 2 | **2** | 1 |
| flow reversal `-c < v < 0` | out | out | in | 1 | **1** | 2 |
| supersonic inflow `v > c` | in | in | in | 3 | **3** | 0 |
| supersonic outflow `v < -c` | out | out | out | 0 | **0** | 3 |

The default branch prescribes `rho` and `p` from the reservoir and takes `v`
from the interior through the valve. **For subsonic inflow the count is
correct**, and the review's section 3.2 is right to say that two prescribed
thermodynamic quantities are not by themselves an over-specification. The same
count is already stated in `docs/supersonic_molecular_base.md` section 2.2,
and `check_base_inflow_is_subsonic` reports (does not enforce) the supersonic
case.

What the code does not have is the other two rows. There is no reversal
branch: `max(v_1, 0)` clips the velocity but does not change how many
conditions are imposed, so at `v < 0` the boundary still states two
thermodynamic quantities where only one is admissible.

### 2.3 The `massflux` branch

`Base velocity: massflux` writes `rho`, `v` and `p`, all three. For subsonic
inflow that is one condition too many by the table above, and it is the
combination the standard reference calls out as inadmissible: Carlson (2011,
NASA/TM-2011-217181, Table 1 footnote) allows `rho` and `T`, `rho` and `p`, or
`rho` and `U`, and states "but not `U` and `p`".

Two qualifications. First, CETIMB does something that looks the same and is
not: Koskinen et al. (2013a, Icarus 226, 1678, section 2.1.1) specify `T_0` and
`p_0`, derive `rho_0` from the ideal gas law, and then "the steady state
continuity equation `rho_0 v_0 r^2 = F_c`, where `F_c` is the flux constant,
was used to calculate the velocity `v_0` at the lower boundary during each time
step", with `F_c` "solved self-consistently by the model". Their `v_0` is a
consistency relation with the model's own output, not an independent datum, and
their scheme is an operator-split van Leer advection with a Crank-Nicholson
diffusion step, not a Riemann-solver finite-volume method. Whether that
constitutes an over-specification in the characteristic sense is not settled by
their papers, and this note does not claim it is.

Second, EXHALE's `base_flux_const` is a running average over `[j_min:N]`, so
the third condition it imposes is not even the local one; it is the wind's flux
constant divided by `rho_bc r^2`.

### 2.4 The location mismatch, and its size

`Rec_BC` writes the face state with expressions evaluated at ghost-cell
quantities:

* the `massflux` branch divides by `r(index)**2`: `r(0)`, the ghost **centre**,
  while the state is consumed at `r_edg(0)`;
* the `hydrostatic_base` branch extrapolates to `r(index)` using `r(1)`, `r(2)`
  as the abscissae of `W_in(3,1)`, `W_in(3,2)`, which in `Rec_BC` are the face
  pressures at `r_edg(1)` and `r_edg(2)`;
* the default branch writes constants, which is the same statement in its most
  visible form: the ghost cell average and the face value are set to the same
  two numbers.

The size is one half-cell of hydrostatic stratification. On the gate grid,

```
rho(r_edg(0))/rho(r(0)) - 1 = -3.815e-3
```

i.e. the face state is 0.38 percent too heavy and too high in pressure. That
ratio is `dr/(2H)`, so it is **first order in `dr`**, but the momentum
residual it produces is that error divided by `dr`, so the residual itself is
**zeroth order**: refining the grid does not remove it. Section 3.1 measures
exactly that.

### 2.5 Two paths read two different interior velocities

The valve reads `W_in(2,1)`. In `Apply_BC_W` that is the cell-1 **average**
velocity. In `Rec_BC` it is `WL_out(2,1)`, the velocity reconstructed at
`r_edg(1)`: the *top* face of cell 1, a full cell above the boundary. The
`hydrostatic_base` and `base_ghost_T_continuous` branches read `W_in(3,1)` and
`W_in(3,2)` with the same ambiguity, and `base_ghost_T_continuous` divides a
face pressure by `n_part_cell1`, a cell-averaged particle count.

So the ghost cell average that `eval_dt`, the ionization sweep and the column
integrals see, and the face state that the Riemann solver sees, are built from
different interior data. On the converged hot Uranus of
`docs/p44_base_sawtooth.md` section 2 the two velocities are `v(1) = -250.2`
and a reconstructed face value between `v(1)` and `v(2) = +14.7 cm/s`; under
the softplus valve (`valve_eps = 1e-4`) the ghost velocity is
`eps^2/(4|v|)`, so the two paths differ by roughly a factor two in the ghost
velocity. Both are small; the point is that they are not the same boundary.

## 3. Baseline measurements with the present boundary

### 3.1 (a) A stationary hydrostatic atmosphere

**Method.** Initialize the exact isothermal hydrostatic atmosphere at rest and
(A1) evaluate the discrete steady residual `R = dF - S` cell by cell, and (A2)
march it to a fixed physical time and read the drift. `R_mom` is normalized by
the local weight `rho g`, so 1 means the boundary leaves the whole weight of
the cell unbalanced.

Five closures are compared. Four are the code's own. The fifth,
`reference`, is not implementable (it hands the boundary the exact analytic
solution, both as ghost **cell averages** and as the **face** state) and is
included as an upper bound on what any face-consistent construction can buy.
`pinnedvalues` is the legacy pair `(rho_bc, p_bc)` held at the face *and* the
ghost with a **fixed zero velocity**, i.e. the legacy boundary with the
interior-velocity copy removed; it separates the valve from the pinned numbers.

A1, cells 1-12, legacy closure, gate grid:

```
#  j    r-1 [R_p]      R_mass          R_mom        R_mom/(rho*g)
   1    1.917733E-04   -7.712567E+00   -9.941287E+00   -2.513517E-01
   2    3.835465E-04   -7.774505E-04   -8.124199E-04   -2.070642E-05
   3    5.753198E-04    1.608493E-06    1.919556E-04    4.931843E-06
   ...
  12    2.301279E-03    1.482092E-06    1.773914E-04    4.898214E-06
```

Cell 1 is out by 25 percent of its own weight; cell 2 by `2.1e-5`; cells 3-12
sit at the interior truncation error, `4.9e-6`. **The boundary error is
confined to two cells and is four to five orders above the scheme's own.**

Full table, `t_end = 1.0 R_p/v0` (the drift is saturated by `t = 0.25`; see
below), all quantities in code units:

| closure | `drc` | `R_mom(1)/rho g` | `R_mass(1)` | `v(1)` | base-face `rho v r^2` | max `|rho/rho_hs - 1|` |
|---|---|---|---|---|---|---|
| legacy isothermal | 2.0e-4 | 2.5135e-1 | -7.7126 | -1.2137e-3 | 1.62e-6 | 4.389e-3 |
| | 1.0e-4 | 2.5067e-1 | -7.7167 | -6.0516e-4 | 1.90e-7 | 2.189e-3 |
| | 5.0e-5 | 2.5033e-1 | -7.7187 | -3.0106e-4 | 1.55e-7 | 1.097e-3 |
| `pinnedvalues` (no valve) | 2.0e-4 | 2.5135e-1 | -7.7126 | -1.2137e-3 | 1.59e-6 | 4.392e-3 |
| | 1.0e-4 | 2.5067e-1 | -7.7167 | -6.0479e-4 | 4.04e-7 | 2.190e-3 |
| | 5.0e-5 | 2.5033e-1 | -7.7187 | -3.0126e-4 | 3.72e-8 | 1.098e-3 |
| `Hydrostatic base: True` | 2.0e-4 | 4.6668e-3 | 6.087e-2 | -3.2015e-5 | -2.23e-5 | 8.706e-4 |
| | 1.0e-4 | 2.3259e-3 | 3.056e-2 | -7.9269e-6 | -5.54e-6 | 2.351e-4 |
| | 5.0e-5 | 1.1582e-3 | 1.527e-2 | -1.9632e-6 | -1.38e-6 | 8.155e-5 |
| `Base ghost temperature: continuous` | 2.0e-4 | 7.5037e-1 | 7.6333 | -6.4995e-3 | -2.55e-3 | 1.095e-1 |
| | 1.0e-4 | 7.5010e-1 | 7.6796 | -3.3358e-3 | -1.44e-3 | 5.824e-2 |
| | 5.0e-5 | 7.4997e-1 | 7.7027 | -1.6879e-3 | -7.64e-4 | 3.002e-2 |
| **`reference`** | 2.0e-4 | **6.2469e-5** | **-2.063e-3** | -1.80e-7 | 5.6e-9 | 3.074e-5 |
| | 1.0e-4 | 3.1437e-5 | -1.004e-3 | -4.30e-8 | 8.3e-9 | 2.663e-5 |
| | 5.0e-5 | 1.5731e-5 | -4.941e-4 | -1.13e-9 | 1.1e-8 | 2.553e-5 |

Five statements come straight off this table.

**(i) The legacy boundary is inconsistent, not merely inaccurate.**
`R_mom(1)/rho g` reads 0.2514, 0.2507, 0.2503 across a 4x refinement: it
converges to `1/4` and does not go to zero. `R_mass(1)` is likewise constant at
`-7.72`. This is the arithmetic of section 2.4: a face error of `dr/(2H)`
divided by `dr` is a fixed fraction of `rho g`, at any resolution.

**(ii) The velocity artifact is first order in `dr`, and it is the production
one.** `v(1)` reads `-1.2137e-3, -6.0516e-4, -3.0106e-4, -1.4999e-4` over a
`8x` refinement (ratios 2.006, 2.010, 2.007), which in physical units is
`-247.1, -123.2, -61.3, -30.5 cm/s`. `docs/p44_base_sawtooth.md` section 7
measures `-248.3, -122.2, -60.0 cm/s` on the full production code at
`dr = 1.93e-4, 9.60e-5, 4.78e-5 R_p`. **The probe reproduces the converged
hot-Uranus base sawtooth to within 1-2 percent with no chemistry, no radiation, no
molecules and no wind**: an isothermal monatomic atmosphere at rest, the same
grid, the same boundary. Whatever else is true of that artifact, it is a
property of the Euler + gravity + boundary discretization alone.

The drift saturates: `v(1)` at `t_end = 0.05, 0.1, 0.25, 0.5, 1.0, 2.0` reads
`+5.49e-4, -8.65e-4, -1.350e-3, -1.269e-3, -1.214e-3, -1.216e-3`. It is a
fixed error the layer settles onto, not a growing instability.

**(iii) The valve is not the cause of the hydrostatic defect.**
`pinnedvalues` reproduces the legacy `R_mom(1)`, `R_mass(1)` and `v(1)` to six
digits. Removing the interior-velocity copy entirely changes nothing here. The
defect is the pinned numbers evaluated at the wrong location.

**(iv) A face-consistent boundary removes it.** `reference` cuts
`R_mom(1)/rho g` by 4024x and `R_mass(1)` by 3738x on the gate grid, drops
`v(1)` from `-247 cm/s` to `-0.04 cm/s`, and cuts the base-layer density drift
by 143x. Its own residual is first order in `dr` and its `v(1)` is at the
interior truncation floor. This is the size of the prize, measured, not
argued.

**(v) `Hydrostatic base: True` is well balanced and `continuous` is worse than
the default.** The pressure extrapolation restores the discrete balance almost
completely (`4.67e-3` of `rho g`, first order in `dr`; `v(1)` second order).
The continuous-`T` ghost is three times worse than the default and, like it,
zeroth order.

That last result must be read with its limit. **This test cannot see a thermal
runaway**, because it has no radiative source term. `Hydrostatic base: True`
runs the hot-Uranus base away to 10200 K in the full code and
`Base ghost temperature: continuous` collapses the H2 front onto the boundary
(`TO_BE_DONE.md` item (W)). Both let a thermodynamic variable float with no
statement of what sets it, and the source terms then set it. The lesson for
Phase C is that a well-balanced boundary and a thermodynamically closed
boundary are two requirements, and the present key set satisfies at most one at
a time.

PLM instead of WENO3 changes none of this: legacy `2.5154e-1 / -7.7206`,
`reference` `3.507e-4 / -1.091e-2`. The `reference` floor is higher on PLM
because the reconstruction is second order; the legacy value is unmoved because
it is a boundary error, not a reconstruction error.

### 3.2 (b) Reflection of an acoustic pulse leaving through the base

**Method.** A Gaussian pressure pulse of relative amplitude `A`, centred 30
cells above the base and 6 cells wide, is launched **downward** (a pure
`v - c` wave: `p' = A p`, `rho' = p'/c^2`, `v' = -p'/(rho a)`). The boundary's
own drift (section 3.1) is larger than a linear pulse, so the same march is run
twice, with and without the pulse, and every amplitude is read on the
**difference** of the two states. The two acoustic amplitudes are

```
f_dn = p'/(rho a) - v'      (toward the base)
f_up = p'/(rho a) + v'      (away from it)
```

normalized by `sqrt(rho a) r`, which removes the geometric and stratification
growth of a linear wave. `|R|` is the maximum normalized `|f_up|` after the
pulse has reached the boundary, divided by the incident `|f_dn|`.

| closure | incident | reflected | **`|R|`** | peak `dp/p` | peak `v` |
|---|---|---|---|---|---|
| legacy isothermal (hard valve) | 1.5796e-4 | 1.5056e-4 | **0.953** | -8.76e-5 | -6.81e-5 |
| `Hydrostatic base: True` | 1.5796e-4 | 1.5904e-4 | **1.007** | +9.31e-5 | +7.15e-5 |
| `Base ghost temperature: continuous` | 1.5796e-4 | 1.3508e-4 | **0.855** | +7.95e-5 | +5.99e-5 |
| `pinnedvalues` (legacy pair, `v` fixed) | 1.5796e-4 | 4.902e-6 | **0.031** | +1.56e-5 | -1.64e-5 |
| `Base velocity: massflux` (`F_c` fixed) | 1.5796e-4 | 4.902e-6 | **0.031** | +1.56e-5 | -1.64e-5 |
| **`reference`** | 1.5796e-4 | 4.348e-7 | **0.0028** | +2.93e-6 | -2.66e-6 |

`|R|` is amplitude independent (0.95401, 0.95317, 0.95303, 0.95304 at
`A = 1e-3 … 1e-6`), so this is the linear reflection coefficient, and it is
grid independent (0.95317 at both `dr = 2e-4` and `1e-4`).

Three statements.

**(i) The legacy base reflects essentially everything, with the sign
inverted.** `|R| = 0.95` and the returning pressure perturbation is negative
where the incident one was positive: the boundary is a **pressure node**, the
acoustic behaviour of a free end. An acoustic disturbance launched at the base
does not leave; it comes back.

**(ii) The reflection follows the velocity, not the pinned numbers.**
Holding the same legacy pair `(rho_bc, p_bc)` at the face with a *fixed* zero
velocity drops `|R|` from 0.953 to 0.031: a boundary whose left state is held
entirely at the unperturbed background is 97 percent transmitting, because the
incoming characteristic then carries no perturbation.

> **Corrected 2026-09-03, by the implementation of section 7.** This
> measurement is right and the inference first drawn from it was wrong. What
> makes the frozen state transmitting is that it prescribes all THREE
> variables, which for a subsonic inflow is one condition too many; the low
> reflection is a symptom of over-specification, not of correctness. A
> boundary that prescribes the admissible two and takes the velocity from the
> outgoing characteristic reflects just as completely as the valve does --
> measured `|R| = 0.956` for the characteristic condition of section 7 against
> `0.953` here. Section 7.3 gives the reason and section 6's D6 carries what
> is left open.

**(iii) The two defects are independent and have opposite causes.** Section
3.1 (iii): the valve contributes nothing to the hydrostatic imbalance. Section
3.2 (ii): the pinned numbers contribute nothing to the acoustic reflection.
This is the same separation `docs/p44_base_sawtooth.md` section 9 reached by a
different route (raising `T_0` closes the thermal gap and leaves the flux error
untouched), and `TO_BE_DONE.md` item (W) records: "the thermal mismatch and the
flux mismatch are two separate defects".

### 3.3 Why the earlier `J^-` experiment failed

`docs/heitr_metals_and_convergence_notes.md` section 2.5 records three reverted
base experiments: CFL `0.6 -> 0.2` (`du` floor `~2.8 -> ~1.9`), zero-gradient
base velocity (`du ~ 2.3`), and "velocity from outgoing Riemann invariant
`J^- = v - 2c/(gamma-1)`", which "did not help HD189733b **and** broke
WASP-121b's convergence (drove a large spurious base inflow when cell 1 is
hot)".

That failure is diagnostic rather than discouraging, and the measurements above
say why. The characteristic count of that experiment was right (two reservoir
conditions plus one interior relation). What was wrong is the *pair*: with
`rho` and `p` both pinned at `T_0`, the face sound speed is the reservoir's,
and when the interior is much hotter the compatibility relation converts a
thermodynamic mismatch the boundary is over-stating into a velocity. **A
compatibility relation is only as good as the reservoir pair it is closed
with**, which is item (W) restated. Section 6 is where that has to be decided.

No characteristic boundary has ever been implemented in this code: the `J^-`
experiment set one component and kept the rest, and
`docs/p55_base_mode.md` section 8 records "No characteristic base boundary was
tried". Neither validation (a) nor (b) had been run before; the closest
existing measurements are the imposed `2 dr` entropy perturbation of
`docs/hd189_base_checkerboard.md` section 5.3 and the base momentum-balance
audit of `docs/newton_scaling_and_base_wall.md` section 1.

## 4. Candidate designs at the face

### 4.1 The face solve

Let `W_i = (rho_i, v_i, p_i)` be the interior reconstruction of cell 1 down to
`r_edg(0)`. That state already exists as `WR_out(:,0)` and is already left
alone by `Rec_BC`. Let `W_b` be the boundary face state to be solved for, and
let `rho_ref`, `c_ref` be reference values taken from the **interior** state
(section 4.5 explains why the interior and not the ghost).

The outgoing acoustic relation is the linearized `v - c` compatibility
condition,

```
    p_b - rho_ref c_ref v_b  =  p_i - rho_ref c_ref v_i .          (C-)
```

This is the wave amplitude of the `u - c` characteristic, and it is verified
against the published originals as follows.

**The wave amplitudes and their normalization.** Poinsot and Lele (1992)
Eqs. (15) to (17) give the characteristic velocities of the `x_1` direction,
`lambda_1 = u_1 - c`, `lambda_2 = lambda_3 = lambda_4 = u_1`,
`lambda_5 = u_1 + c`, with `c^2 = gamma p / rho` (their Eq. (18)), and their
Eqs. (19) to (23) give the amplitude variations, of which the two acoustic ones
are

```
    L_1 = lambda_1 ( dp/dx_1 - rho c du_1/dx_1 )                    P&L (19)
    L_5 = lambda_5 ( dp/dx_1 + rho c du_1/dx_1 )                    P&L (23)
```

**The normalization includes the eigenvalue factor `lambda_i`**, and this is
the same in both originals: Thompson (1990) Eqs. (44) to (47) are identical
expressions in the same primitive ordering `U = (rho, p, u_1, u_2, u_3)` with
left eigenvectors `l_1 = (0, 1, -rho c, 0, 0)` and `l_5 = (0, 1, rho c, 0, 0)`
(his Eq. (46)), and Thompson (1987) Eq. (19) defines `L_i = lambda_i l_i dU/dx`
for outgoing waves. There is no factor of two, no `1/2` and no `1/(rho c)`
inside the `L_i` themselves; those appear only when the `L_i` are assembled
into the primitive equations.

**The relation `(C-)` is the linearized amplitude of that wave.** Poinsot and
Lele (1992), page 109, in the paragraph following their Eq. (23): for the
upstream-propagating wave associated to `lambda_1 = u_1 - c`, with `p'` and
`u'` the pressure and velocity perturbations, "the wave amplitude
`A_1 = p' - rho c u'` is conserved along the characteristic line
`x + lambda_1 t = const`". Written between two states at the same radius, one
of them the interior extrapolation and one the boundary state, that
conservation is exactly `(C-)`.

**The LODI relations.** Poinsot and Lele (1992) section 2.3 derives them from
their Eqs. (9) to (13) and "neglecting transverse and viscous terms". In
primitive variables the system is their Eqs. (24) to (28); the three that a
one-dimensional problem keeps are

```
    d(rho)/dt + (1/c^2) [ L_2 + (1/2)(L_5 + L_1) ] = 0               P&L (24)
    d(p)/dt   + (1/2)(L_5 + L_1)                   = 0               P&L (25)
    d(u_1)/dt + (1/(2 rho c)) (L_5 - L_1)          = 0               P&L (26)
```

and their Eqs. (29) to (32) recast the same system in `T`, `m_1 = rho u_1`,
`s` and `h`. Their Eqs. (33) to (36) invert it for the normal gradients; the
first three of those are Thompson (1990) Eq. (48), term for term.

**The source terms, which the LODI form drops and this problem does not have
the right to drop.** Thompson writes the conservation law with an
inhomogeneous vector from the start: Thompson (1987) Eqs. (1) and (2) carry a
term "which often arises from divergence terms in nonrectangular geometries",
and Thompson (1990) Eq. (1) carries the same vector `D`, whose footnote 1 on
page 441 names both contents that this code has,

```
    (1/r^2) d(r^2 rho u_r)/dr  ->  d(rho u_r)/dr + (2/r) rho u_r
```

for the spherical divergence, and "source terms, such as heating, cooling, or
gravitational forces". The projection of that vector onto the left eigenvector
stays inside the characteristic equation: Thompson (1987) Eq. (13) is
`l_i dU/dt + lambda_i l_i dU/dx + l_i C = 0`, his nonreflecting condition
Eq. (17) is `(l_i dU/dt + l_i C) = 0` and not `l_i dU/dt = 0`, and Thompson
(1990) Eqs. (50) to (54) carry the gravitational acceleration explicitly, the
`u_1` row being

```
    d(u_1)/dt + (1/(2 rho c)) (L_5 - L_1) + (transverse) - g_1 = 0   T90 (52)
```

against Poinsot and Lele's (26), which has no `g_1` because their derivation is
Cartesian with no body force. **In spherical one-dimensional flow there are no
transverse terms at all, so "drop transverse and viscous" removes nothing, and
the source projection is the entire difference between the two papers' momentum
relation.** The gravitational part of it is of the same order as the quantity
`(C-)` constrains: over a half cell the neglected `rho c g_r dr` is the
hydrostatic pressure change `dp` itself.

`(C-)` escapes this only because **both of its states are evaluated at the same
radius** `r_edg(0)`, so the source projection is common to the two sides and
cancels. The hydrostatic and geometric terms are carried instead by the
continuation that produces `W_i` from the cell 1 average (section 4.4, and the
implemented `continue_hydrostatic_isentrope` of section 7.1). Any future form
of this boundary that relates the cell 1 CENTRE to the face, or that time
advances the face state in the Thompson manner, must put `l_1 C` back in
explicitly; the steady form of Thompson (1990) Eq. (63) is
`L_1 = -rho c g_1` in a gravitational field, not `L_1 = 0`.

`(C-)` is also preferred over the closed-form invariant
`J^- = v - 2c/(gamma-1)` because the latter assumes a constant `gamma`, which
this code does not have (section 4.5).

For **subsonic inflow** `0 < v_b < c_b`, (C-) plus two reservoir conditions
closes the face state. That count is the published one: Poinsot and Lele
(1992), section 3.1, "For a subsonic three-dimensional flow, four
characteristic waves are entering the domain (Fig. 1), `L_2`, `L_3`, `L_4`, and
`L_5`, while one of them (`L_1`) is leaving the domain at the speed
`lambda_1 = u_1 - c`. Therefore, the density `rho` (or the pressure `p`) has to
be determined by the flow itself." Removing the two transverse waves leaves two
entering and one leaving, which is two conditions and one interior relation.
Thompson (1990) Eq. (68), left column, is the same statement as a recipe for
the lower boundary `x_1 = a_1` with `0 < u_1 < c`: `L_1` computed from its
definition, `L_2` and `L_5` specified. Sections 4.2-4.4 are the three choices
that remain.

### 4.2 Reservoir pairs

Admissibility is Carlson (2011) Table 1 and its footnote: for a subsonic inflow
face one may specify `(p_t, T_t, flow direction)`, or `rho` and `U`, or `rho`
and `p`, "but not `U` and `p`".

The primaries were checked for this list on 2026-09-07 and do not contain it:
what Poinsot and Lele (1992) Table III gives are four inflow sets they used,
SI 1 `(u_1, u_2, u_3, T)`, SI 2 `(u_1, u_2, u_3, rho)`, SI 3
`(u_1 - 2c/(gamma-1), u_2, u_3, s)` and SI 4 the non-reflecting
`L_2 = L_3 = L_4 = L_5 = 0`, with the well-posedness column pointing at
Strikwerda and at Oliger and Sundstrom rather than proving it there. The
exclusion of `(U, p)` is therefore still on Carlson's authority and not on the
originals'. Two things the originals do say and this table did not: velocity is
an admissible member of an inflow pair (SI 1 and SI 2 both use it), and in one
dimension SI 3 and SI 4 coincide, since "they both express the conservation of
entropy and the non-reflection of acoustic waves at the inlet section".

| pair | what it means physically | solve | notes |
|---|---|---|---|
| **(A) `(p, T)`** | the handoff level states both, as a lower/middle-atmosphere model does. Koskinen et al. (2013a) sec. 2.1.1: "We specified `T0` and `p0` at the lower boundary, and used them to calculate `rho0` from the ideal gas law"; Koskinen et al. (2022) app. B: "The lower boundary conditions are the temperature, pressure (`p1 = 1 μbar`), and species mixing ratios that are based on our lower atmosphere models." | `p_b = p_res`; `rho_b` from `(p_b, T_res)` and the base composition; `v_b` from (C-) | the closest to what EXHALE already has (`base.inp`'s `T_base`/`p_base`, the `Lower atmosphere profile:` reader). Admissible: `(p, T)` is `(p, rho)`. **The measured objection stands**: on the gate the interior accepts `p` to 0.1 percent and rejects `T_0` by 5-14 percent (`docs/p44_base_sawtooth.md` sec. 8.3, item (W)). This pair does not fix that; it makes it explicit, and moves the question to where `T_res` comes from. |
| **(B) `(p, s)`** | pressure plus specific entropy. Identical to (A) for a fixed composition and caloric EOS, but stated so the boundary is closed by a quantity the caloric EOS defines rather than by a temperature the composition has to convert | `p_b = p_res`; `rho_b` from `(p_b, s_res)`; `v_b` from (C-) | better behaved than (A) when the base composition changes (`x_H2` handoff), because `s` carries the rovibrational heat capacity with it instead of leaving the conversion to `n_part_cell1`. Costs an `s(p, rho, x_H2)` inversion; the caloric module already has every piece. |
| **(C) `(T, rho)`** | the legacy pair in another parameterization | `rho_b = rho_res`; `p_b` from `(rho_b, T_res)`; `v_b` from (C-) | admissible, and it is what the code effectively does now except that `v` is copied rather than solved. It keeps item (W)'s over-statement intact. Listed for completeness. |
| **(D) `(p, rho v)`** | pressure plus the mass flux, i.e. CETIMB's velocity condition rewritten as a reservoir datum | `p_b = p_res`; `v_b` from (C-); `rho_b = F_c/(r_b^2 v_b)` | solvable and not the forbidden `(U, p)` combination, since `F_c` constrains `rho v` and not `v`. ~~It is singular at `v_b -> 0`, which is exactly the state the hot-Uranus base sits in.~~ **Corrected 2026-09-07 against the original.** The singularity is in the closure written in the solve column above, not in the condition. Thompson (1990) section 3.2.7.2 imposes a mass flux at a subsonic inflow face by requiring `d(rho u_1)/dt = 0`, which with his Eqs. (50) and (52) and no transverse terms gives `2 u_1 L_2 + (u_1 - c) L_1 + (u_1 + c) L_5 = 2 rho c^2 g_1` (his Eq. (69)), and at a LOWER boundary `x_1 = a_1` with `0 < u_1 < c` this is solved for the entering acoustic amplitude as `L_5 = [ 2 rho c^2 g_1 - 2 u_1 L_2 - (u_1 - c) L_1 ] / (u_1 + c)` (his Eq. (70), left column). The denominator is `u_1 + c`, so the prescription is regular through `u_1 = 0` and through sign changes of `u_1`; it is his Eq. (71), which solves the same relation for the entropy amplitude `L_2` instead and carries `1/(2 u_1)`, that is singular at stagnation. Note also that the published form carries `rho c^2 g_1`: the mass flux condition in a gravitational field is not the zero-gravity one. What remains true is the second half of the original objection, that (D) leaves the base entropy as an output rather than a datum. |

**Recommendation, for decision.** (B), falling back to (A) where no entropy is
available. It is the pair the lower-atmosphere handoff can actually state
(`p_base` fixes the level; the Koskinen-2022-style handoff states `T` and the
composition at that level), it is admissible by the characteristic count, and
it puts the caloric EOS inside the boundary condition rather than beside it.
This is a decision for the user, not a conclusion: see section 6.

### 4.3 Flow reversal and the supersonic branches

The rule has to be explicit, evaluated on the **face** state, and it changes
the number of conditions, not the value of one of them.

```
M_b = v_b / c_b  evaluated at the face
```

| branch | conditions | interior relations |
|---|---|---|
| `M_b >= 1` (supersonic inflow) | all three: `rho_res`, `v_res`, `p_res` | none |
| `0 < M_b < 1` (subsonic inflow) | two, per section 4.2 | (C-) |
| `-1 < M_b <= 0` (reversal) | **one**: `p_b = p_res` | (C-) **and** the entropy relation `s_b = s_i` |
| `M_b <= -1` (supersonic outflow) | none | all three: extrapolate `W_b = W_i` |

**The four rows are the published count**, checked branch by branch on
2026-09-07. Thompson (1990) section 3.2.4 (supersonic inflow): "so that all
`lambda_i < 0`. Consequently we must specify all values of `L_i` from boundary
conditions", the row above with three conditions. Section 3.2.5 (supersonic
outflow): "we can specify no boundary conditions at all, and the evolution of
the flow at the boundary is determined completely by interior data", the
extrapolation row. Section 3.2.6 (subsonic outflow, his `0 < u_1 < c` at the
upper face, which is the reversal branch at ours): one amplitude specified,
the rest computed, i.e. one condition and two interior relations. Section
3.2.7 with his Eq. (68): at a lower face with `0 < u_1 < c`, `L_1` computed
from its definition and `L_2` and `L_5` specified, i.e. two conditions and one
interior relation. Poinsot and Lele (1992) Tables III and IV give the same
counts for the subsonic rows, four Euler conditions at a three-dimensional
subsonic inflow (two once the transverse velocities are removed) and one at a
subsonic outflow. That the count itself changes along a boundary is stated in
Thompson (1990) section 3.2.3: "The number and type of boundary conditions
required may vary from time to time or place to place along the boundaries,
even for one particular problem."

The reversal row is the one the code does not have. Physically it is gas
falling back into the lower atmosphere: the reservoir can still state the
pressure it presents, but its temperature and density are no longer its to
state, because the gas arriving carries the interior's entropy. That is the
correct replacement for `max(v_1, 0)` and for the softplus valve.

**The branch is not differentiable at `M_b = 0`**, and the steady Newton solver
is why `valve_eps` exists. Two ways out, and this is a decision item:

* blend the reversal and inflow closures over a narrow `|M_b| < eps_M` window,
  so the residual stays differentiable at the price of a small window in which
  neither closure is exactly imposed;
* keep the branch exact and give the solver a boundary-aware line search.

The first is the direct analogue of `valve_eps` and is the cheaper start; note
that the softplus valve's own floor `eps^2/(4|v|)` is currently what the ghost
velocity reads on the converged hot Uranus (item (V)), so any blending
parameter has to be reported in the setup report and read in diagnostics.

Since the count depends on `M_b` and `M_b` depends on the solved state, the
branch must be selected on a **predictor** (the interior state `W_i`) and the
selection frozen for the solve, then checked. `check_base_inflow_is_subsonic`
already measures the Mach number and reports crossings; it should become the
selector rather than a warning.

### 4.4 Ghost cell averages consistent with the face state

The face solve gives `W_b` at `r_edg(0)`. Two things then need it:

1. the Riemann flux, which takes `WL(:,0) = W_b` directly;
2. **the ghost cell averages**, which are read by `eval_dt`, the ionization
   sweep, `calc_column_dens`, the WENO3 reconstruction of `WR(:,0)` (through
   `dW(:,0)` and `dW(:,-1)`), and the residual's conservative state.

Writing `W_b` into both is the present defect with a better `W_b`. The ghost
averages have to be the cell averages of a continuation of `W_b` **below** the
face, to the order of the reconstruction. The construction that section 3.1
measured is:

```
rho_g(r) = rho_b * exp( -(phi(r) - phi(r_edg(0))) / c_T^2 )
p_g(r)   = p_b * (rho_g(r)/rho_b)^k
v_g(r)   = rho_b v_b r_edg(0)^2 / (rho_g(r) r^2)
```

with `c_T^2 = p_b/rho_b` and `k = 1` for the isothermal continuation (pair A
with `T` held), or the polytropic exponent of the chosen entropy for pair (B);
`phi` is the code's own `grav_func::phi`, so the continuation uses the same
potential as the source term. The ghost average is then the volume-weighted
integral of each over `[r_edg(j-1), r_edg(j)]`; the probe uses 8-point
Gauss-Legendre, which is exact to far beyond the reconstruction's order and
costs nothing at two cells per step.

The mass-flux form for `v_g` is deliberate: it makes the ghost carry a constant
`rho v r^2`, so the reconstruction of `WR(:,0)` does not see a spurious
velocity gradient at the boundary. With `v_b = 0` it reduces to `v_g = 0`.

`reference` in section 3.1 is this construction with the exact analytic profile
substituted for the continuation, so its numbers are the **floor**, not the
prediction. How close a real continuation gets depends on how well the
isothermal/polytropic continuation matches the true sub-face structure over one
cell; on the gate grid that is `dr/H = 7.6e-3` of a scale height, so the
continuation error should be far below the `3.8e-3` half-cell error it replaces.
That has not been measured and is the first thing an implementation should
measure.

**A cheaper alternative exists and should be measured against it:** solve the
face state, then set the ghost averages by requiring the *reconstruction
operator itself* to return `W_b` at the face, i.e. invert the two-cell WENO3
stencil for `W_0` given `W_b`, `W_1`, `W_2`. That is exact by construction and
carries no continuation assumption, but it makes the ghost average a function
of the interior state (and of the WENO weights), which the steady solver then
has to differentiate. The continuation form above is the more conservative
start.

### 4.5 Composition and the caloric equation of state

Three couplings, each of which the present component-wise closure is silent
about.

**(i) `c_ref` in (C-) belongs to the interior.** The outgoing wave is an
interior wave: it is the state of cell 1, not the reservoir, that determines
what leaves. So `c_ref = sqrt(gamma_eff(1) p_i / rho_i)` with
`gamma_eff` from `adiabatic_index_from_state(1, rho_i, p_i)`. Getting this
wrong is precisely the mechanism by which the earlier `J^-` experiment turned a
hot cell 1 into a spurious inflow (section 3.3). `Num_flux` already keeps
`gam_L` and `gam_R` separate across the same face (`Num_Fluxes.f90` lines
64-75), so a composition jump at the boundary is already handled there; the
boundary condition has to use the same distinction.

**(ii) The face composition is the base composition.** `rho_b` and `p_b` are
related through `nk_per_mass` and the H2 fraction of the base, which come from
`comp_ntot_bc` / the `x_H2` handoff (`q_H2_base`, or the profile). The face
state must be built with the **base** composition (`x_h2(0)`, `nk_per_mass(0)`)
and the interior relation with cell 1's. `caloric_state_from_composition` fills
`x_h2`, `nk_per_mass` and `molecular_cell` over `1-Ng:N+Ng`, so the ghost
entries exist; what does not exist is a statement of which one each half of the
boundary condition uses.

**(iii) `gamma_eff` is not constant, so no closed-form invariant exists.** With
`e = p/(gamma_eff - 1)` and `gamma_eff = gamma_eff(x_H2, T)`, the Riemann
invariant `v - 2c/(gamma-1)` is not integrable in closed form and is wrong
wherever the H2 fraction or the temperature varies across the stencil, which
is the whole molecular base. This is the decisive reason to use the
differential relation (C-) instead. The general statement is Thompson (1987)
Eqs. (14) and (15): a wave amplitude `V_i` with `dV_i = l_i dU + l_i C dt`
exists only if the coefficients satisfy Pfaff's condition for the integrability
of differential forms, "a condition not met for the fluid equations", while
"the characteristic form of (13) holds true independent of (14)". A constant
`gamma` and no source term is one of the special cases in which the integration
does succeed, and it is the case this code is not in. Where `caloric_mixture_active` is false the
two coincide to the order of the linearization, so an atomic run is unaffected
by the choice.

**(iv) The composition/energy ordering.** The (AB) change established that the
residual's ghosts must be filled *after* the ionization sweep, from the same
composition the fluxes use. A face solve inherits that requirement exactly: the
face state's `rho_b <-> p_b <-> T_b` conversion needs a composition, and it must
be the converged one of the residual evaluation, not the previous trial's. The
face solve therefore belongs at the same point in `eval_residual` where the
ghost fill now sits.

## 5. Validation set

Tests (a) and (b) are implemented and their present-boundary baselines are in
section 3; (c) and (d) were run for section 7 and (e) was not. The status
column records the outcome.

| test | what it decides | acceptance | status |
|---|---|---|---|
| **(a) stationary hydrostatic atmosphere** | is the boundary consistent? | `R_mom(1)/rho g` must fall with `dr`, and be within a factor of a few of the interior | **PASSED** (sec. 7.2). Legacy `2.51e-1`, zeroth order; Phase C `1.75e-4`, first order, 1437x lower |
| **(b) outgoing acoustic pulse** | does the boundary trap acoustic energy? | `\|R\| < 0.1` at the linear amplitude, grid and amplitude independent | **FAILED, and shown to be unreachable by any static reservoir pair** (sec. 7.3). Legacy `0.953`, Phase C `0.956`. Carried as D6 |
| **(c) manufactured smooth inflow** | is the boundary the design order? | drive an exact steady solution and require cell 1 within a small factor of the interior truncation error | **PASSED** (sec. 7.2). Mass row 1.8x the interior, momentum row 0.68x |
| **(d) hot-Uranus reload and march** | does it move the root toward the marching attractor? | peak-to-peak face mass flux at `r = 1.005, 1.03, 1.10, 1.20` against the (AB) baseline `0.4443, 0.2517, 0.0403, 0.0084` | **RUN, and NEGATIVE** (sec. 7.4): `0.9284, 0.3074, 0.0473, 0.0096`. The departure grows. `G_dt` versus `assemble_residual` (Phase B's B2) not measured |
| **(e) three base grids** | does the answer converge? | `Base grid [dr,cells]` `2e-4 50`, `1e-4 100`, `5e-5 200` with `Grid cells` adjusted; `log10 Mdot`, the H2 front and the closure ratio of item (AA) must converge, and `v(1)` must fall faster than first order | **not run.** It was gated on (a)-(c) passing, which they now do, so it is the next measurement |

Two remarks on (d) and (e). First, they must be run **after** (a)-(c) pass, not
alongside: a boundary that fails (a) will move (d) for reasons no one can
attribute. Second, (e) is the test on which `Hydrostatic base: True` and
`Base ghost temperature: continuous` are known to fail in ways (a) cannot see
(section 3.1), so it is the gate that a well-balanced-but-thermally-open
closure will not pass.

Probe (a)/(b) reruns are one command each:

```
cd .../scratchpad/phaseC
BASEBC_MODE=<isothermal|hydrostatic|continuousT|massflux|pinnedvalues|reference> \
BASEBC_TEST=<A|B|AB> BASEBC_TEND=1.0 BASEBC_DRC=2.0e-4 BASEBC_NLOW=50 \
OMP_NUM_THREADS=4 ./base_bc_probe.x
```

## 6. Decisions

D1-D5 were put to the user and answered on 2026-09-03; the answers are
recorded inline. D6 was raised BY the implementation and is open.

This is a public numerical-model change: it alters the boundary-value problem,
and it will move every converged state. The following are the user's, not the
implementer's.

**D1. The reservoir pair. ANSWERED: (B) `(p, s)`.** (A) `(p, T)`,
(B) `(p, s)`, (C) `(T, rho)`, or (D) `(p, rho v)`, section 4.2. The answer was
given when (D) was believed singular at stagnation; section 4.2 now records
that it is not, so if the answer rested on that, it is worth revisiting.
Note also that in the characteristic reading (B) is Thompson (1990)'s
`L_2 = 0` (constant inflow entropy, his Eq. (68)) together with his constant
pressure condition `L_5 = -L_1` (his Eq. (67)), which is the fully reflecting
choice among the three he lists; that is the same fact as section 7.3's
measurement. The
recommendation was (B) with (A) as the fallback, and (B) is what section 7
implements. The consequence of any of them is that `T_0` stops being an
independent boundary datum in the sense it is now: with the pair fixed, the
face temperature is whatever the pair and the composition give, and cell 1 will
no longer sit against a pinned `T_0` contact. Item (W)'s "a ghost that pins `p`
and lets `rho` and `T` float together" is *not* one of these: that is one
condition, not two, and it is admissible only in the reversal branch.

**D2. The reversal rule, and its smoothing. ANSWERED: explicit branches on the
characteristic count, blended.** Section 4.3. Implemented as a cubic
smoothstep of width `base_face_mach_blend = 1e-6` in the face Mach number,
a named module constant rather than a user key (section 7.1), reported in the
setup report. Its influence is measured in section 7.4.

**D3. What happens to the existing keys. ANSWERED: retire the four, keep
`Base BC`, and REJECT a retired key at startup rather than ignoring it.**
One kind of thing under one rule (the project's own standard) argues for
retiring rather than accumulating:

| key | proposal |
|---|---|
| `Valve eps:` | **retire**. Its job is taken over by the reversal branch of D2. |
| `Hydrostatic base:` | **retire**. The face condition is well balanced by construction (section 3.1 (iv)); this key achieves well-balancing by leaving `T` unstated, which is the failure mode of item (W). |
| `Base ghost temperature: isothermal\|continuous` | **retire**. `isothermal` is pair (A)/(C); `continuous` is the worst closure in section 3.1 and has no reservoir meaning. |
| `Base velocity: massflux\|valve` | **retired**. With pair (B) it is over-specification. It was off by default and `docs/base_breathing_progress.md` records it as giving the cleanest base at 5x the convergence cost. |
| `Base BC: density\|pressure` | **keep**. It sets the base *level*, not the closure, and the 1 microbar level is a physics choice (Koskinen et al. 2022 sec. 3.1) independent of this work. |
| `Shapiro filter:` | **keep as is, off**. It touches the marching state only and is measured 42x worse than nothing on the molecular hot Uranus (`docs/p55_base_mode.md` sec. 5). Not part of this change. |

**D4. What is expected to move, qualitatively.** No number below is predicted;
these are the directions the measurements above make likely, and they are
stated so a surprise is recognizable.

* The base sawtooth should go. Section 3.1 (ii) shows it is entirely this
  boundary's, and (iv) shows a face-consistent boundary removes it in the same
  test. `v(1)` should stop being `-248 cm/s` at `dr = 2e-4` and stop being
  first order in `dr`.
* The cell-1 residual rows should stop being `O(1)`. Item (V) records
  `R_1/s_1 = -0.22`, `R_2/s_2 = 0.14`, `R_3/s_3 = -0.99` on the molecular hot
  Uranus, unmoved by three row scales, and section 143's measure "cannot drive
  down" cell 1. That is the same cell and the same order as the `0.25` of
  section 3.1.
* ~~The base-layer acoustic mode should damp instead of ringing.~~
  **Withdrawn 2026-09-03 (section 7.3).** A reservoir that holds a static
  quantity is an acoustic node whatever else it does, so the characteristic
  condition reflects as completely as the closure it replaces and this
  prediction had no basis. Item (AA)'s layer swings and
  `docs/p54_base_layer_mass_flux.md` section 5's `2050-2300` step period are
  not addressed by Phase C.
* **The mass-loss rate should not move much.** The wind above `1.4 R_p` is
  insensitive to every base variant measured so far (item (V): flux spread
  `2.4e-3` in all five runs, `log10 Mdot` `10.28-10.34`), and section 3.1 (ii)
  shows the artifact is confined to two cells. If `Mdot` moves by more than the
  cycle noise, that is a result to investigate, not to accept.
* Goldens will not be byte-identical, for every case. This is a physics/model
  change, so a refresh is expected and belongs at the end of the series.

**D5. Scope. ANSWERED: all branches in the first implementation.** Whether the
reversal branch and the supersonic branches are part of the first
implementation or a second one. The inflow branch alone fixes
sections 3.1 and 3.2 for every case measured here; the reversal branch is what
makes `max(v_1, 0)` unnecessary and what a breathing base needs.

**D6. The acoustic reflection (raised by the implementation, section 7.3).**
The characteristic condition reflects sound as completely as the closure it
replaces, and section 7.3 shows this follows from stating two static reservoir
quantities and cannot be cured by choosing a different pair. That measurement
stands.

> **Corrected 2026-09-07 against the originals.** The sentence that used to
> follow, "the only form that is transparent to sound without corrupting the
> steady state is a relaxation condition, which is not a function of the
> instantaneous state and so has no fixed point the steady residual can
> express", is wrong on its second half. Two published forms are transparent
> and are algebraic in the instantaneous face state, so a steady residual can
> carry either as one more algebraic row:
>
> * **Thompson (1990) Eq. (64), the nonreflecting condition in a gravitational
>   field.** His Eq. (63), the characteristic equation with the transverse
>   terms dropped, is
>   `(dp/dt - rho c du_1/dt) + L_1 + rho c g_1 = 0`, so the amplitude of the
>   wave stays constant not when the entering amplitude vanishes but when it
>   equals the source projection: `L_5 = rho c g_1` at a lower boundary
>   `x_1 = a_1`, `L_1 = -rho c g_1` at an upper one. Every symbol on the right
>   is a face quantity. Thompson's own caveat is the one this note reached by
>   measurement from the other side: for a problem in pressure equilibrium
>   "the nonreflecting boundary conditions destroy the pressure balance"
>   (his section 3.2.6.2), which is why he offers the force-free condition
>   (his Eq. (65)) and the constant pressure condition (his Eq. (67),
>   `L_5 = -L_1`) beside it. The constant pressure condition is the fully
>   reflecting one, and it is what pair (B) amounts to.
> * **Poinsot and Lele (1992) Eq. (40), the partially reflecting relaxation.**
>   `L_1 = K (p - p_inf)` with `K = sigma (1 - M^2) c / L`, where `M` is the
>   maximum Mach number in the flow, `L` a characteristic size of the domain
>   and `sigma` a constant; "when `sigma = 0`, Eq. (40) sets the amplitude of
>   reflected waves to 0". This is a function of the instantaneous face
>   pressure alone, so it has a fixed point, and it exists precisely to stop
>   the mean drift that a purely nonreflecting condition leaves: their section
>   2.1 warns that conditions which "impose only values for derivatives and no
>   constraint for the mean values" "may lead to a drift of the mean
>   quantities", and Eq. (40) is their answer to it. `sigma` interpolates
>   between transparent and pinned, which is the dial this note assumed did
>   not exist.
>
> Neither form changes the number of conditions: each replaces one of the two
> admissible conditions of section 4.2, not adds to them. What is NOT settled
> by the literature is whether either keeps the base level this code needs,
> since `p_base` is a physics datum here and not a far-field constant; that is
> what the decision below is now about.

Three ways forward, none taken
here: leave it (the steady state is what this code computes, and a converged
state carries no wave); replace one of the two static conditions with
Poinsot and Lele's Eq. (40) at a `sigma` that has to be measured, which keeps
the count and gives up an exact `p_base`; or reformulate the steady residual to
carry a boundary state of its own. This decides whether item (AA)'s base-layer
ringing is in scope at all.

## 7. Implementation and measurements

Section 152 of `docs/Update_EXHALE_stage1.md`. The patch is
`docs/phaseC_base_bc.diff` (8 files, 28 hunks: one new module, plus the
retirements); the probe is `docs/phaseC_probe.diff`, test code only. Nothing
is applied to the tree.

### 7.1 What was built

`src/modules/states/base_boundary.f90`, a new module, holds the whole
boundary. It has one entry point, `base_boundary_states`, which takes the
first interior CELL AVERAGE and returns three things at once: the primitive
state at the face `r_edg(0)`, the ghost cell averages, and the point state at
the face below the ghost. `Apply_BC_W` calls it and writes the ghosts;
`Rec_BC` reads the face state it stored. That is the "one kind, one rule" the
brief asked for: the two paths can no longer disagree, because there is only
one answer.

Inside it:

* `characteristic_base_face_state` continues cell 1 to the face along **cell
  1's own hydrostatic isentrope** (not a linear extrapolation, and not the
  reconstruction stencil, which reaches into the ghost and would make the
  boundary read its own output), applies the linearized `v-c` compatibility
  relation `p_b - rho_i c_i v_b = p_i - rho_i c_i v_i` with `c_i` at cell 1's
  own `gamma_eff`, and closes it with the two reservoir conditions.
* `continue_hydrostatic_isentrope` is the shared continuation: `h + phi`
  uniform, solved by Newton on `T` with `dh/dT = c_v + 1` from the caloric EOS,
  then `d ln rho = c_v d ln T` with a midpoint `c_v`. Its validity range
  (increments below 0.1 in `ln T`, error below `1e-7`) is stated at the
  routine.
* `base_ghost_averages` builds the ghost cell averages as the 8-point
  Gauss-Legendre volume averages of that continuation below the face, with the
  velocity carrying a constant `rho v r^2`.
* The reservoir is `(p, s)` at the base LEVEL `r = 1`, carried as the isentrope
  through `(p = ntot_bc + dp_bc, T = T0)` at the base composition, and
  `base_boundary` continues it to the face itself. So `Base BC: pressure`
  still names a level that does not move when the grid does, which is the
  half-cell error of section 2.4 removed at the level as well as at the face.

Branches, on the interior predictor Mach number at the face:

| `M_i` | reservoir | interior |
|---|---|---|
| `>= 1` | would owe three; a `(p, s)` reservoir has two | flagged, not invented |
| `0 < M_i < 1` | `p` and `s` | compatibility |
| `-1 < M_i <= 0` | `p` | compatibility and the interior entropy |
| `<= -1` | none | all three |

Only the ENTROPY source changes between the inflow and reversal rows, and only
the whole state changes at the supersonic-outflow row, so both handovers are
one cubic smoothstep each rather than a branch. The window is
`base_face_mach_blend = 1e-6`, a named module constant with its justification
at the declaration: the physical base Mach number of the gate is `6e-6`
(`rho_b v_b r^2 = F_wind` at the base density gives `v_b ~ 1.6 cm/s` against
`c ~ 2.6e5 cm/s`), so `1e-6` leaves every intended operating point wholly
inside the inflow branch.

Supersonic inflow is reported, not repaired: a `(p, s)` reservoir states two
conditions and the model has no lower-atmosphere velocity to state a third, so
`check_base_inflow_is_subsonic` says the state is outside the stated problem
rather than a velocity being invented for it.

Retired, and now rejected at startup with a message naming the replacement
(`retired_base_key` in `input_read.f90`, `error stop`): `Valve eps`,
`Hydrostatic base`, `Base ghost temperature`, `Base velocity`. With them go
`valve_eps`, `hydrostatic_base`, `base_ghost_T_continuous`, `base_v_massflux`
and `base_flux_const`, the JFNK hand-off's `valve_eps = 1e-4` arming, the
mass-flux running average in the marching loop and in `init`, and the two
restart-header fields that carried them. `n_part_cell1` is KEPT and reused: it
is the cell-1 particle count the face solve needs for `T(1)` and `nhat(1)`.
`Base BC: density|pressure` is kept unchanged.

### 7.2 (a), (b) and (c) with the characteristic condition

Same probe, same background, same grids as section 3.

**(a) Stationary hydrostatic atmosphere**, `t_end = 1.0 R_p/v0`, WENO3:

| `drc` | `R_mom(1)/rho g` | vs legacy | `R_mass(1)` | `v(1)` [cm/s] |
|---|---|---|---|---|
| 2.0e-4 | 1.749e-4 | **1437x** | 3.084e-3 | -0.278 |
| 1.0e-4 | 8.638e-5 | 2902x | 1.561e-3 | -0.040 |
| 5.0e-5 | 4.282e-5 | 5846x | 7.832e-4 | -0.0069 |
| 2.5e-5 | 2.127e-5 | 11762x | 3.913e-4 | -0.00009 |

The residual is now **first order in `dr`** (ratios 2.02, 2.02, 2.01) where the
old closure was zeroth order and stuck at `0.2503`. `v(1)` falls faster than
first order and is at the interior truncation floor: the base sawtooth is gone,
not merely reduced, and the `-247, -123, -61, -30 cm/s` ladder of section 3.1
(ii) becomes `-0.28, -0.04, -0.007, -0.0001 cm/s`.

PLM is a different story and is reported as it is: `R_mom(1)/rho g` is smaller
still (2.79e-5 at `drc = 2e-4`), but the base-layer density drift is `6.5e-3`
and grows with refinement, as it did with the old closure (`9.7e-3`). PLM's
base layer is limited by the second-order reconstruction, not by the boundary,
and Phase C does not change that. PLM is stage 1 of a `PLM+WENO3` run.

**(b) Acoustic reflection**: `|R| = 0.956`, amplitude independent
(`1e-3 … 1e-5` all give `0.956`), falling slowly with `dr` (0.956 at `2e-4`).
**The target was not met**, and section 7.3 says why it could not be.

**(c) Manufactured smooth inflow** (new): an EXACT steady spherical isentropic
solution -- `rho v r^2 = F` and `v^2/2 + h + phi = B`, which solves the Euler
equations with gravity, so nothing has to be forced -- driven through the
boundary at base Mach `1.0e-3`, with the reservoir set to that solution's own
state at `r = 1`. The steady residual of the exact solution then measures the
boundary against the scheme it bounds:

| | cell 1 | interior max, cells 3..50 | ratio |
|---|---|---|---|
| `|R_mass| dr/(rho c)` | 6.397e-7 | 3.481e-7 | **1.8** |
| `|R_mom|/(rho g)` | 1.970e-4 | 2.877e-4 | **0.68** |

The boundary is at the accuracy of the interior scheme: it is no longer a
distinguished cell.

### 7.3 The reflection, and why the design note's prediction was wrong

Section 3.2 measured that freezing the ghost velocity dropped `|R|` from 0.953
to 0.031 and concluded that the velocity copy was the reflector. The
implementation shows the inference was wrong. A frozen state prescribes all
three variables, which for subsonic inflow is one too many; its transparency is
over-specification, not correctness. A boundary that prescribes the admissible
two and takes the velocity from the outgoing characteristic reflects
completely: `|R| = 0.956` against the old `0.953`.

The reason is general and worth stating once. Two characteristics enter a
subsonic inflow face; the reservoir states two things; if both are STATIC
quantities then the incoming acoustic amplitude is pinned, and a wave arriving
from inside is reflected. That is true of `(p, s)`, of `(rho, s)` and of
`(rho, p)` alike -- the pair does not matter, the staticness does.

The alternative was implemented and measured, as
`base_incoming_invariant_weight` in `base_boundary.f90`: state the incoming
acoustic invariant at rest, `p_b + rho c v_b = p_res`, instead of the pressure.
Both forms are admissible, both are identical at rest, and:

| | weight 0 (`p_b = p_res`) | weight 1 (invariant) |
|---|---|---|
| (b) `|R|` | 0.956 | **0.0024** |
| (a) `R_mom(1)/rho g` | 3.53e-4 | 1.75e-4 |
| **(c) `R_mom(1)/rho g`** | **1.970e-4** | **1.323e-1** |

Test (c) decides it. `p_b + rho c v_b = p_res` puts the static pressure at the
face below the reservoir's by `rho c v_b`, a fraction `M` of the pressure,
where a real reservoir's Bernoulli drop is `rho v^2/2`, a fraction `M^2`. It
over-states the drop by `1/M`, turning the wind's own mass flux into a base
pressure error that the momentum residual then multiplies by `H/dr`: at base
Mach `1e-3` and `H/dr = 131` that is the measured `0.13`, 460 times the
interior truncation error. **The invariant form buys acoustic transparency by
corrupting the steady state, so it is rejected.** The constant stays at 0 with
this measurement recorded at its declaration.

Getting both would need a RELAXATION form -- the incoming amplitude driven
toward the reservoir at a finite rate (Rudy & Strikwerda 1980, in the
formulation Poinsot & Lele 1992 use) -- which is not a function of the
instantaneous state and therefore has no fixed point the steady residual can
express. That is a real open item and it is D6 below, not a setting.

### 7.4 Production

All runs are the Phase C build against the (AB) build of
`.../scratchpad/p_ab`, same cases, same inputs.

**The hot-Uranus gate rung, reloaded** (`Load IC? True`, `Solver: Newton
100.0`, `Secondary_ionization: Immediate`), the case (AB) measured, so it is
the one comparison with a baseline:

| | (AB) | Phase C |
|---|---|---|
| `info` | 0 | 0 |
| JFNK iterations | 19 | **14** |
| `\|\|R\|\|` accepted | 8.961e-06 | 7.848e-06 |
| flux gate (`r >= r_flux`) | 2.501e-04 | 2.501e-04 |
| `r >= 1.03` | 3.432e-03 | 3.402e-03 |
| `r >= 1.10` | 2.375e-03 | 2.372e-03 |
| `log10 Mdot` | 10.34 | 10.34 |
| **`v(cell 1)`** | **-251.88 cm/s** | **+14.40 cm/s** |
| `T(cell 1)` | 1116.68 K | 1116.72 K |
| H2 front (`2 n_H2 = n_H`) | 1.0449 `R_p` | 1.0449 `R_p` |
| base Mach, max | - | 3.364e-03, subsonic |

**The base sawtooth is gone in the production code**, which is what the probe
predicted: `v(1)` stops being a `-250 cm/s` artifact and becomes `+14.4 cm/s`,
outward, of the order the wind's own mass flux through the base density gives.
Everything the wind is judged by is unmoved: the same flux gate to four
figures, the same `Mdot`, the same H2 front, the same base temperature to
0.04 K. The H3+ column moves 0.9 percent. The solve is 26 percent shorter.

**The same rung from a COLD start** (`Load IC? False`, `du_th 0.5 1.0e-3`)
lands on the same root: `info = 0` at `||R|| = 1.405e-06`, the same three gate
numbers to four figures (3.402e-03 / 2.372e-03 / 2.501e-04), the same
`log10 Mdot = 10.34`, and a state that agrees with the reloaded root to
**3.9e-5 in `rho`, 1.9e-4 in `v` and 1.3e-5 in `T`**, with `v(1) = +14.393`
against `+14.395 cm/s` and `T` at `r - R_p = 8e-3` equal to 818.01 K in both.

**That is a direct repair of the reproducibility defect of
`docs/p44_base_sawtooth.md` section 8.4**, which measured three `info = 0`
states of this configuration putting the H2 front anywhere in 1.053-1.089
`R_p` and `T` at the same radius anywhere in 583-949 K, on this same default
base grid. Two independent paths to the root now agree to five figures.

The cold start is not cheap: 820 JFNK iterations across three attempts
(`info = 1`, then `info = 2`, then `info = 0`) against the reload's 14. Whether
that is harder than the same cold start under (AB) is **not established**:
(AB) did not run one, so there is no baseline.

**WASP-121b, `wasp_full_newton`** (atomic, He 2^3S + metals, Newton):

| | (AB) | Phase C |
|---|---|---|
| `info` | 0 | 0 |
| JFNK iterations | 96 | 108 |
| `\|\|R\|\|` accepted | 9.026e-06 | 9.910e-06 |
| flux gate | 5.242e-05 | 5.238e-05 |
| `r >= 1.03` | 6.751e-05 | 6.733e-05 |
| `log10 Mdot` | 13.21 | 13.21 |
| `v(cell 1)` | 1372.5 cm/s | 1705.9 cm/s |
| max rel. state change | - | `rho` 5.9e-3, `p` 8.8e-3, `T` 3.2e-3 |
| base Mach, max | - | 0.331, subsonic |

A strongly irradiated atomic hot Jupiter converges to the same gates and the
same mass-loss rate, with the state moving by well under a percent and the
base velocity moving by 24 percent: the boundary doing its job on the one cell
that is its own.

**The reload-and-march departure** (root reloaded, `Solver:` removed,
`du_th 1.0e9 1.0e-30`, 3000 steps, `EXHALE_P54_TS=10`; peak-to-peak of the
Riemann face mass flux in units of the root's own `F_0`):

| root, binary | 1.005 | 1.03 | 1.10 | 1.20 | `du` at 3000 |
|---|---|---|---|---|---|
| (AB) root, (AB) binary (control, rerun here) | 0.4443 | 0.2517 | 0.0403 | 0.0084 | 3.2250e-04 |
| **Phase C root, Phase C binary** | **0.9284** | **0.3074** | **0.0473** | **0.0096** | 3.2904e-04 |

The control reproduces `p_ab_ghost_order.md` section 5 exactly (0.4443, 0.2517,
0.0403), so the harness and the normalization are the same ones.

**This is a negative result and it is the honest headline of the production
set: Phase C does not repair the marching departure of item (AD), and it makes
it worse, 2.09x at `r = 1.005` and 1.22x at 1.03.** It is consistent with
section 7.3: the boundary still reflects `|R| = 0.956` of what the base layer
launches at it, so the layer still rings, and the steady residual and the
production update still do not share a fixed point. One plausible mechanism for
the increase is that the retired one-way valve, by clipping `v` at zero, was
also a nonlinear amplitude limiter on that oscillation, and removing it lets
the layer swing more freely; **that is a hypothesis, not a measurement.** `du`
at 3000 steps is unchanged at 2 percent.

**The identity of that increase was measured afterwards** with the section 145
update-map tool (`EXHALE_UPDATE_MAP=1`), each build on its own accepted root,
and **it is not the boundary**. The residual of the accepted root falls by
three orders of magnitude -- mass `9.916e-3 -> 6.067e-6`, momentum
`7.473e-3 -> 5.453e-6`, energy `1.486e-2 -> 8.462e-6` -- while the update-map
defect `|G_dt + R|` at `dt/100` falls only from `2.460e-5` to `1.351e-5`
(mass), `2.643e-5` to `9.840e-6` (momentum), and **does not move at all in the
energy row, `1.892e-4` against `1.891e-4`**. The energy row also stops falling
with `dt` in both builds, and its breakdown columns put it in the
operator-split chemistry/energy update rather than in `hydro`. Defect over that
row's own residual at `dt/100` therefore goes from `2.5e-3, 3.5e-3, 1.3e-2` to
`2.2, 1.8, 22`.

So the root became 1600x more accurate and the splitting error did not move.
The (AB) root's own residual was the same size as that splitting error, so
marching it moved it little; the Phase C root is a far better solution of a
discrete steady problem the production update does not quite solve, so marching
moves it further. **The departure grew because the residual shrank.** The
remaining obstruction is Phase B's B2, not this boundary. `Apply_BC` moves the
interior by a relative `0.000E+00` in both builds, so the (AB) B1 invariant
holds here too. Full table: `docs/phaseC_prep_measurements.md` section 2.

**Regression matrix**, ten cases, Phase C against (AB), single-threaded, same
step caps:

| case | max rel diff rho / v / p / T | worst r | stop step, du (AB -> C) | log10 Mdot | v(cell 1) [cm/s] |
|---|---|---|---|---|---|
| wasp_full | 4.4e-01 / 2.6e+00 / 2.9e-01 / 1.7e-01 | 1.0278 | 11289,9.9989e-04 -> 18026,9.9091e-04 | 13.26 -> 13.24 | -726.20 -> 2383.28 |
| wasp_he23off | 4.4e-01 / 2.7e+00 / 2.9e-01 / 1.7e-01 | 1.0284 | 11276,9.8873e-04 -> 18015,9.8369e-04 | 13.26 -> 13.24 | -755.85 -> 2380.44 |
| mol_base_handoff | 2.3e-01 / 4.3e-01 / 3.1e-01 / 1.6e-01 | 1.1355 | 12000,1.1472e+00 -> 12000,1.2914e+00 | 10.43 -> 10.52 | -139.22 -> 48.99 |
| mol_metals | 2.4e-01 / 5.1e-01 / 3.2e-01 / 1.6e-01 | 1.1321 | 12000,1.2754e+00 -> 12000,1.3237e+00 | 10.43 -> 10.52 | -139.13 -> 53.63 |
| mol_lyman_werner | 1.7e-01 / 2.2e-01 / 2.4e-01 / 1.3e-01 | 2.6316 | 12000,5.5038e-01 -> 12000,8.6940e-01 | 10.47 -> 10.50 | -114.22 -> 42.50 |
| mol_diffusion | 2.3e-01 / 4.3e-01 / 3.1e-01 / 1.6e-01 | 1.1355 | 12000,1.1470e+00 -> 12000,1.2915e+00 | 10.43 -> 10.52 | -139.21 -> 49.02 |
| mol_ir_bands | 2.4e-01 / 5.1e-01 / 3.2e-01 / 1.6e-01 | 1.1321 | 12000,1.2753e+00 -> 12000,1.3236e+00 | 10.43 -> 10.52 | -139.21 -> 53.53 |
| mol_sec_ion | 2.0e-01 / 1.0e+00 / 3.2e-01 / 1.5e-01 | 1.0508 | 12000,4.0963e-01 -> 12000,1.2205e+00 | 10.38 -> 10.47 | -601.15 -> 7.11 |
| mol_carrier | 2.3e-01 / 1.9e-01 / 3.1e-01 / 1.6e-01 | 4.6810 | 12000,5.9036e-01 -> 12000,1.2315e+00 | 10.42 -> 10.50 | -139.52 -> 30.05 |
| lower_profile | 9.6e+00 / 5.9e+01 / 3.6e+00 / 7.5e-01 | 1.1879 | 12000,2.8228e+01 -> 12000,8.1893e+02 | 9.56 -> 8.44 | -581.62 -> 542.74 |

Every case moves, which is expected of a model change, and **every case flips
the sign of `v(cell 1)`**: `-726, -756, -139, -139, -114, -139, -139, -601,
-140, -582 cm/s` become `+2383, +2380, +49, +54, +43, +49, +54, +7, +30,
+543`. That single column is the boundary defect of section 3.1 (ii) being
removed on ten configurations at once.

Three qualifications, because most of these rows do not mean what they look
like:

* **The eight molecular rows are fixed-step snapshots, not answers.** Their
  `du` is 0.4-1.3 in both builds: they stop at step 12000 in the middle of a
  relaxation, by design. A `log10 Mdot` of `10.43 -> 10.52` at the same step
  count is a different point in a transient, not a different mass-loss rate.
  The one molecular case that is actually solved -- the gate rung above --
  gives `10.34 -> 10.34`.
* **The two WASP rows are `du`-stops**, and a `du = 1e-3` stop is known to
  carry a several-percent path-dependent spread in `Mdot`
  (`docs/auto_ic_design.md`; the standing rule is to Newton-finish before
  quoting one). `13.26 -> 13.24` is inside that. The Newton-finished
  `wasp_full_newton` gives `13.21 -> 13.21`. They also take 60 percent more
  steps to reach the same threshold (11289 -> 18026, 11276 -> 18015), which
  is a real cost and is reported as one.
* **`lower_profile` is not converged in either build** (`du = 28.2` before,
  `819` after) and its row compares two transients. Its `9.56 -> 8.44` should
  not be read as a mass-loss rate at all.


**`-fcheck=bounds,do,mem`**, the bounded HD 209458 b case of
`backup/regression/run_fcheck.sh` (1500 steps, cold start, He 2^3S + metals):
**clean**, no runtime check fired, `log10 Mdot = 9.61`.

## 8. What was not established

* ~~The size of the continuation error in section 4.4.~~ **Measured, section
  7.2.** The implemented continuation gives `R_mom(1)/rho g = 1.75e-4` against
  the analytic `reference`'s `6.25e-5` on the same grid, i.e. it keeps
  1437x of the 4024x the reference bounded, and unlike the reference it is
  first order rather than stuck.
* **Everything about the radiative source term, in the probe.** The probe has
  no heating, cooling, ionization or molecular chemistry, so it cannot see a
  thermal runaway and it cannot say what the face condition does to the base
  energy budget of item (AA). The production measurements of section 7.4 are
  where that is decided, and test (e) after them.
* **The grid-convergence test (e).** Not run.
* **Whether pair (D) is well posed near stagnation.** Section 4.2, corrected:
  the published constant mass flux prescription is regular through `u_1 = 0`,
  so the question is no longer the singularity but whether a base whose entropy
  is an output rather than a datum behaves on the hot Uranus. Not measured.
* ~~The Poinsot & Lele normalization.~~ **Settled, section 4.1** (2026-09-07).
  The three journal versions were read. The amplitudes carry the eigenvalue
  factor `lambda_i` in both papers, `(C-)` is their `u - c` amplitude, the
  branch counts of section 4.3 are theirs, and no constant was found wrong.
  Two statements of this note were corrected instead: pair (D)'s singularity
  (section 4.2) and D6's impossibility argument (section 6). What section 4.1
  now records as open is not the citation but the physics: the source
  projection `l_1 C`, which the spherical geometry and the gravitational field
  both feed, cancels out of `(C-)` only because both of its states sit at the
  same radius, and any reformulation that moves off that has to carry it.
* **Whether CETIMB's three boundary data over-specify its scheme.** Section
  2.3: their papers state what is prescribed but not the discrete location of
  the boundary (grid level, not ghost cell) nor the flux stencil that consumes
  it, so the question cannot be settled from the publications.
* **Nothing was measured on the production binary.** Every number in section 3
  is from the isolated Euler probe. The agreement with
  `docs/p44_base_sawtooth.md` section 7 in section 3.1 (ii) is a comparison
  against that document's measurements, not a rerun of them.

## 9. Files

Published sources, read in the journal versions in `../references/`:

```
Thompson_1987JCP_68_1.pdf          K. W. Thompson, J. Comput. Phys. 68, 1 (1987)
                                   "Time Dependent Boundary Conditions for
                                   Hyperbolic Systems". Eqs. (13), (17)-(19).
Thompson_1990JCP_89_439.pdf        K. W. Thompson, J. Comput. Phys. 89, 439 (1990)
                                   "Time-Dependent Boundary Conditions for
                                   Hyperbolic Systems, II". Eqs. (44)-(52),
                                   (62)-(71), footnote 1 p. 441.
Poinsot_1992JCP_101_104.pdf        T. J. Poinsot and S. K. Lele, J. Comput.
                                   Phys. 101, 104 (1992) "Boundary Conditions
                                   for Direct Simulations of Compressible
                                   Viscous Flows". Eqs. (15)-(40),
                                   Tables II-IV.
```

```
docs/phaseC_characteristic_base_bc.md   this note
src/tests/base_bc_probe.f90             the (a)/(b) probe
src/modules/states/Apply_BC.f90         + probe_base_reference gate (scratch only)
logs/final_A.txt                        the section 3.1 table as produced
logs/final_B.txt                        the section 3.2 table as produced
```

Build:

```
PATH=/usr/bin:$PATH make -j8
PATH=/usr/bin:$PATH gfortran -O2 -fopenmp -Jbuild -Ibuild src/tests/base_bc_probe.f90 \
  build/parameters.o build/species_table.o build/molecular_infrared_data.o \
  build/caloric_eos.o build/UW_conversions.o build/grav_field.o \
  build/define_grid.o build/set_gravity_grid.o build/Apply_BC.o build/PLM_rec.o \
  build/Reconstruction.o build/speed_estimate_ROE.o build/Num_Fluxes.o \
  build/low_mach_dissipation.o build/Source.o build/RK_rhs.o -o base_bc_probe.x
```
