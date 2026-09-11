# The carrier marching step and the carrier steady residual (P158)

Measured 2026-09-04 on a scratch copy of the tree at `EXHALE.x`
md5 `063b2673eb0d4f1ca342cc07bf942b4c` (regression 11/11, `run_fcheck.sh`
clean).  Nothing in this note was applied to the tree.

The question was whether the marching loop's carrier stage and
`carrier_steady_residual` are the same discrete operator, because
`docs/p23_transport_on_state.md` section 10 measured that
`|G_dt + R|` for the carrier row does not fall when the step is refined,
while every hydrodynamic row falls first order.

## 1. What the two paths actually share

`photochemical_transport_step` and `carrier_steady_residual` already call one
assembly, in the same order and with the same arguments:

| stage | marching step | steady residual |
|---|---|---|
| carrier fractions, background, element headroom | `carrier_state` | `carrier_state` |
| base face composition | the inner ghosts of `f_sp` | the inner ghosts of `f_sp` |
| grid and gravity | `carrier_geometry(dt_code, ...)` | `carrier_geometry(1, ...)`, then `dt = 1e30` |
| the face mass flux and the mixture mass the advective term rides on | not read: the marching rows carry no advective term | `carrier_advective_state` |
| molecular diffusivities | `carrier_diffusivities` | `carrier_diffusivities` |
| face gradient / drift / donor switch | `carrier_face_coefficients` | `carrier_face_coefficients` |
| photolysis and Lyman-Werner columns | `carrier_photolysis` | `carrier_photolysis` |
| material advection | `species_advective_update` inside the Runge-Kutta stages | `carrier_advective_divergence`, the same face fractions, face fluxes and flux divergence (item B5b) |
| signal-crossing rate | `carrier_signal_rate` | `carrier_signal_rate` |
| rows | `solve_carriers` -> `carrier_residual` | `carrier_residual` |

The base inflow condition is built in exactly one place: the inner ghosts of
the species vector, which the ionization sweep pins to the handoff partition
where a handoff states one.  Both paths reconstruct the face composition from
those ghosts wherever the face mass flux flows inward, so the base condition
is identical between them by construction.  **The base Dirichlet was not the
inconsistency**, and the face fluxes were not duplicated.

Amended twice.  Increment B4-1c moved the marching advection into the
Runge-Kutta stages as the divergence of the hydrodynamic face mass flux and
removed the cell-velocity term from the marching rows, which left the
stationary balance carrying a different discretization of the same term.
Item B5b removed that difference: the stationary rows and the fixed-wind
relaxation now form their advective term from the SAME face fractions, face
fluxes and flux divergence (`species_face_flux.f90`), so the two paths differ
in no term.  The record below was measured before B5b and describes the
cell-velocity form.

Three differences exist and only the last of them changes an answer:

1. the step passes its own `dt_code`, the residual substitutes `dt = 1e30`
   (this is the definition of the steady limit, not a difference);
2. the step runs a Newton and then `limit_to_element_budget`;
3. the step ends in `carrier_write_back`, which converts the solved carrier
   fractions back into `f_sp`.

## 2. The operator IS consistent, measured

A probe inside `photochemical_transport_step` evaluates `carrier_residual`
on the entry state with `dt = 1e30` -- the same coefficients the step is about
to use -- and compares it with the rate of change the step produced, at three
cut points: after the Newton, after the element clamp, after the write-back.
On the hot-Uranus carrier gate marched 12 000 steps, four time steps spanning
four decades:

| dt/dt_CFL | after Newton | after clamp | after write-back |
|---|---|---|---|
| 1     | 9.518e-5 | 9.518e-5 | 9.518e-5 |
| 0.1   | 7.145e-6 | 7.145e-6 | 7.145e-6 |
| 0.01  | 8.160e-7 | 8.160e-7 | 8.160e-7 |
| 0.001 | 8.262e-8 | 8.262e-8 | 8.262e-8 |

(largest over cells of `|G + R|` divided by the row's own terms; worst cell
`j = 445`, `r = 2.911`.)  First order, over four decades, and the clamp and the
write-back add nothing.

## 3. Where the dt-independent term comes from

Section 10 measured `G` on the carrier DENSITY, `n(H2) = f_sp * rho * n0`,
across the whole production step.  Splitting that map by stage, at the cell
carrying the worst relative imbalance (`j = 497`, `r = 4.59`, all four steps):

| stage | rate [cm^-3 s^-1] |
|---|---|
| hydrodynamics | +7.154e-1 |
| carrier transport | -5.470e-2 |
| ionization sweep | 0 |
| steady residual R of the carrier row | +5.474e-2 |

The transport stage reproduces `-R` to three digits.  The whole dt-independent
imbalance is the hydrodynamic stage: it changes `rho` at frozen `f_sp`, so the
carrier density it hands on moves by `n_c d(ln rho)/dt`, and no carrier row
contains that term.  It does not fall with the step because it is a
hydrodynamic rate, not a discretization error.

## 4. The physics underneath: the wrong advected variable

That is not only a measurement artifact.  Write the two stages out.  The
hydrodynamic stage gives

    dn_c|hydro = -dt n_c [ div v + v dln(rho)/dr ]

and the carrier stage, transporting the PARTICLE mixing ratio
`X = n_c/n_tot` with `n_tot` frozen, gives

    dn_c|transport = -dt [ v dn_c/dr - n_c v dln(n_tot)/dr ] + dt S.

Their sum is

    dn_c = -dt div(n_c v) + dt S  -  dt n_c v dln(mbar)/dr ,

so the composed map is the conservation law plus a spurious advective term
`n_c v dln(mbar)/dr`.  It is nonzero wherever the mean mass per particle
varies -- that is, across the H2 dissociation front and the ionization front,
which is the whole region the carriers exist in.

The reason is that `n_tot` is not a conserved density: ionization and
dissociation create particles, so `d(n_tot)/dt + div(n_tot v)` is not zero and
a mixing ratio formed against it is not advected.  `rho` IS conserved -- it is
the hydrodynamic mass row -- so the scalar that rides with the flow is the
fraction per unit mass, `y = n_c/(rho n0)`, which is the variable `f_sp`
already stores and the variable the hydrodynamic stage leaves alone.

Molecular diffusion is a different matter and is genuinely per particle: its
flux is `-n_tot (D + K_zz) dX/dr` with `X` the particle mixing ratio.  So the
correct row keeps the diffusive face flux in `X` and everything else in `y`,
with `X = wfac * y`, `wfac = rho n0 / n_tot = mbar / m_amu`.

## 5. The change

- `carrier_state` returns the carrier fraction per unit mass (it is `f_sp`,
  with no conversion), together with `nrho = rho n0` and `wfac`.
- The time term, the advective coefficient `ntv`, the source conversion, the
  Newton diagonal, the element clamp, the column scale and the write-back are
  all weighted by `nrho` in place of `ntot`.
- `carrier_face_flux` takes the two face values of `wfac` and forms the
  diffusive and drift fluxes on `X = wfac y`, so the diffusion is unchanged
  physics.
- The lower ghosts are no longer written by the carrier operator.  They are the
  inflow reservoir: the ionization sweep pins their composition from the
  handoff, and `carrier_base_state` reads it as the base Dirichlet value.
  `solve_carriers` had been mirroring the base cell into them on every Newton
  trial and `carrier_write_back` had been writing that mirror into `f_sp`, so
  the base Dirichlet degenerated into a zero-gradient condition.  The marching
  loop hid it -- its ionization sweep re-pins the ghost every step -- but
  `relax_photochemical_composition` takes many transport steps between two
  sweeps, so from its second pass on it was relaxing against the base cell's
  own H2 rather than against the handoff.

## 6. The verdict measurement

The production map -- the whole marching body, hydrodynamics included --
measured in the row's own unknown, against the carrier steady residual of the
state it started from.  Largest over cells of `|G + R|` divided by the row's
own terms, and the same quantity at three named cells:

Before (carrier as `n_c/n_tot`):

| dt/dt_CFL | base cell j=1 | r = 1.0984 | r = 1.16 | max over cells |
|---|---|---|---|---|
| 1     | 1.892e-5 | 3.843e-3 | 1.125e-5 | 3.843e-3 (j=205) |
| 0.1   | 3.715e-4 | 3.846e-4 | 1.895e-6 | 1.027e-3 (j=446) |
| 0.01  | 5.871e-5 | 3.846e-5 | 1.897e-7 | 1.022e-4 (j=446) |
| 0.001 | 3.058e-4 | 3.846e-6 | 1.930e-8 | 3.058e-4 (j=1)   |

After (carrier as `n_c/(rho n0)`, lower ghosts left to the boundary):

| dt/dt_CFL | base cell j=1 | r = 1.0984 | r = 1.16 | max over cells |
|---|---|---|---|---|
| 1     | 2.000e-5 | 3.836e-3 | 1.100e-5 | 3.836e-3 (j=205) |
| 0.1   | 1.790e-5 | 3.838e-4 | 1.667e-6 | 3.838e-4 (j=205) |
| 0.01  | 2.885e-6 | 3.838e-5 | 1.667e-7 | 3.838e-5 (j=205) |
| 0.001 | 3.112e-7 | 3.838e-6 | 1.665e-8 | 3.838e-6 (j=205) |

Measured on the carrier DENSITY instead, before the change, the same run gives
0.1912, 0.1913, 0.1912, 0.1912 -- flat to four figures over four decades, which
is what section 10 reported.

Two things separate in that table.  The interior cells were already first order
before the change: the two paths were one operator, and refining the step
refines the imbalance.  The BASE CELL was not -- 1.9e-5, 3.7e-4, 5.9e-5,
3.1e-4, with no trend -- and it is first order after it.  The base value the old
row imposed was the ghost's carrier DENSITY divided by the ghost's PARTICLE
density, and the ghost's particle count is a quantity the sweep derives for a
cell that is not part of the domain; the new row imposes the ghost's mass
fraction, which the handoff states directly.

## 7. The coupled solve was not reached

Both operators were restarted from the same 12 000-step hot-Uranus state with
`Solver: Newton` and the JFNK hand-off at `du < 1e-2`, and both marched to a
plateau instead of reaching it.  Sampled every 2000 steps to step 27 914:

| step | du, new operator | du, old operator |
|---|---|---|
| 13 942 | 1.105e-1 | 1.056e-1 |
| 17 934 | 3.305e-2 | 3.274e-2 |
| 21 926 | 3.667e-2 | 3.683e-2 |
| 27 914 | 3.660e-2 | 3.665e-2 |

`dtu` fell to 4.3e-6 in both, so the state has stopped moving while the mass
flux keeps a 3.7% radial spread.  **The variable change does not move that
plateau**, and the coupled Newton was therefore never offered a state in either
run.  The three gates on the new operator's state and which row is blocked are
still unmeasured.  The plateau is the section 157 wall and it is a separate
obstruction from the one this section closes.

## 8. What was not measured

- Whether the front position and the mass-loss rate of a converged hot-Uranus
  model move under the new variable.  The regression case `mol_carrier` is a
  12 000-step snapshot, not a converged state.
- The oxygen carriers (OH, H2O, CO).  The gate runs H2 alone; the change is
  written for all four but only the H2 row has been exercised.
- Any planet other than the hot Uranus.
- `run_fcheck.sh`.  The change alters argument lists and array uses in the
  carrier module, so it should be run before the change is taken further.
  `EXHALE_ELEMENT_ASSERT=1` over 400 steps of the carrier gate reported no
  element-budget violation, which covers the write-back's ghost change but not
  the bounds checking.
