# LHS 1140 b, item L5c: the base layer's discrete energy flux, and why the hydrodynamic solve stops on it

Item L5c of `docs/PLAN_20260913_lhs_stationary.md`: why the hydrodynamic
solve returns `info = 2` at outer pass 1, what the base cells' discrete
energy flux divergence is made of, and what the two existing keys
`Well balanced` (K46) and `Low-Mach damping` (K31b) do to it.

Every number is MEASURED on the binary of 2026-09-13 03:20 (md5
`97e10317a710b9ccc63addbedde3586a`, the copy L5a used, kept at
`<scratch>/L5a_bin/EXHALE.x`) unless it is marked READ. **No file in `src/`
was changed**; this is a diagnosis. The seeds are L5a's, unchanged.

---

## 1. Verdict

**The base layer's O(1) energy residual is a property of the
DISCRETIZATION, not of the state, and one existing key removes it. With
`Well balanced: True` the same S3 seed, the same input, the same binary
reaches a CERTIFIED stationary solution of the LHS 1140 b atomic wind at
outer pass 8 in 10 minutes** (mass row 4.53e-10, momentum 2.29e-13, energy
5.08e-9, the gated element row 8.76e-6 of 1e-5, mass-flux spread 1.5e-11,
`info = 0`), where every seed of L5a stalls at pass 1 with `info = 2` and a
cellwise row of order one.

The mechanism, measured face by face:

- **The discrete energy flux of a base face is not the enthalpy flux the
  flow carries; it is the enthalpy flux carried at the HLLC contact speed,
  and that speed is set by the reconstruction's pressure jump, not by the
  velocity.** `S* = S*_vel + S*_pres` with
  `S*_pres = (p_R - p_L)/(phi_L - phi_R) ~ -(p_R - p_L)/((rho_L+rho_R) c)`.
  In the base layer `S*_vel/c` is 3.4e-7 to 1.1e-6, the flow's own Mach
  number, and `S*_pres/c` is 5e-6 to 9e-4: one to three decades above it
  (table T2).
- **Split that way, the part of `dF_3` the velocity carries is 1.0 to 1.9
  times the net heating in every base cell -- exactly what a stationary
  balance asks for -- and the part the pressure jump carries is 10^2 to
  10^4 times it, with the alternating sign L5a reported** (table T1). The
  kinetic-energy part is 1e-7 to 1e-16 of the face flux and plays no role
  at all.
- **The pressure jump is the WENO3 reconstruction's own third-order
  truncation error on a steeply stratified column, not an odd-even mode and
  not the base ghost.** Three measurements: it survives eight 1-2-1 passes
  over `ln rho, ln v, ln p` and is reproduced within a factor 4 by a ninth-order
  polynomial fit of the same column (table T3); on a refinement of that
  analytic column it converges at order **2.97 to 3.00** (table T4); and it
  is the same size at faces 3 to 50, whose stencils never touch a ghost
  cell. The seed's own 2 dr residue (a 5-point quadratic residual of
  `ln rho` of -1.0e-2 at cell 1, 1e-5 to 1e-4 from cell 5 out) dominates only in cells 1
  to 3 and modulates the sign further out.
- **The row therefore reads its own scale.** `residual_row_scale` for the
  energy row is `max(|dF_3|, |S_3|, heat, cool)`; with `|dF_3|` four decades
  above the heating the row is `|dF_3 - (heat-cool)|/|dF_3| ~ 1` whatever
  the state does, which is L5a's "the row and its own scale are the same
  term" measured at its source.
- **The conditioning is 1/Mach.** A relative pressure perturbation of 1e-8
  at a base cell moves that cell's `dF_3` by 0.1 of the net heating
  (table T5): the energy row's sensitivity to `p` is 9.9e6 in the base
  against 5.8e3 at 1.86 R_p and 2.3e2 at 5.9 R_p, while its sensitivity to
  `v` is 0.06 to 2 everywhere. Holding the row below the solver's 5e-5
  target with the plain flux would need the pressure field specified to
  5e-12 relative.
- **The well-balanced option removes the term algebraically.** It forms the
  Riemann pressure jump as `dp_eq + dev_R - dev_L`, all three of the size of
  the departure from the cell's own hydrostatic equilibrium. On a column
  that IS the discrete equilibrium of that operator the jump is 1e-15 of the
  pressure and `S*_pres/c` is 1e-15, against -8.3e-5 for the plain flux on
  the same column (table T6). Nothing about the physical Mach number
  changes: `S*_vel/c` stays 4.1e-7.

**`Low-Mach damping` is inert on this state**, at both documented
coefficients. It never touches the mass flux by construction, and on the S3
seed it moves the energy row by at most 2.3e-7 relative and the momentum row
by 2.1e-4 (section 5): the velocity field of a stationary seed is smooth, so
the fourth difference the stress is built from is negligible. The key was
written for the HD 209458 b stagnant layer's 2 dr VELOCITY mode; what
obstructs LHS 1140 b is a systematic PRESSURE truncation error, and a
velocity hyper-viscosity cannot reach it.

**Two of L5a's findings are corrected by this.**

1. L5a's "the pass is refused by the element relaxation, not by the
   hydrodynamics" is right about which test fires and wrong about the cause.
   With `Well balanced: True` the very same element operator, on the very
   same seed, **reaches its fixed point in 32 steps at every pass** (drift
   2.46e-1, 1.23e-1, 6.21e-2, ...), exactly as the HD 209458 b control does.
   The `element_step_out_of_bounds` of L5a is downstream of a hydrodynamic
   state the solve could not converge, not an obstruction of its own. This
   is item **L5b**'s answer, arrived at from the other side.
2. L5a's table T3 reads the base face mass flux off the ghost, which carries
   a constant `rho v r^2` by construction (+3.67 F_wind for S3). That is the
   boundary condition's own face state, not the flux the scheme puts through
   that face. MEASURED from the code's own `R_mass(1)` and the face-1 flux,
   the HLLC mass flux through the base face of the S3 seed is
   **-1746 F_wind**, and through face 1 it is +2290 F_wind, for the same
   reason every other base face is wrong: the pressure jump, not the
   velocity, sets `S*`.

**What this does not settle.** The certified state is the solution of a
DIFFERENT discretization from the one every earlier LHS 1140 b number was
produced with, and it is one solve from one seed. Its `Mdot` is 1.95e8 g/s
(log10 8.29) against the archived 6.38e7 (7.81), a factor 3.1; its first
three cells still carry a velocity artifact of 0.3 cm/s against a stationary
base velocity of 0.150 cm/s (P44's cell-1 defect, 1.6x the signal rather
than 1.8e4x); and its base temperature is 405 K at cell 1 against the 226 K
reservoir. Section 7 lists what L4 has to check before that `Mdot` is
quoted.

---

## 2. What the pass-1 solve does, and why `info = 2`

READ from the L5a run logs (`<scratch>/L5a_{S3,S1c,S1}/run.log`); the
runs are L5a's, not re-run.

### T0. The hydrodynamic JFNK of outer pass 1

| seed | iterations | stop | `\|\|R\|\|` handed back | `\|\|Fs\|\|2` start -> end | lam = 1 at | gm = 1 at |
|---|---|---|---|---|---|---|
| S3 | 604 | **no descent direction exists for the banded model at this state** (`\|\|grad merit\|\|` 1.392) | 1.078 | 3.05e+1 -> 1.771e-2 | 592 of 604 | 598 of 604 |
| S1c | 95 | STAGNATED, no improvement of the best iterate in 40 iterations | 1.842 | 2.92e+1 -> 5.75e-2 | -- | 89 of 95 |
| S1 | 67 | STAGNATED, same test | 1.989 | 1.64e+2 -> 4.76e-2 | -- | 64 of 67 |

**The outer iteration cap is 3000, from the caller** (`(JFNK) outer
iteration cap 3000, from the caller`, READ from every one of the three
logs), and no run comes within a factor 4 of it. `EXHALE_JFNK_MAXIT`
REPLACES the caller's cap (`steady_newton.f90` line 416 and 6506, READ, and
marked TEST ONLY there), so the 400 the brief names would LOWER it below the
604 S3 already takes, and naming a larger number changes nothing: neither of
the two stops in the table is the cap.

**More iterations would not have helped, measured.** Over S3's iterations
400 to 604 -- 204 iterations, a third of the run -- the merit `||Fs||2` moved
from 1.780e-2 to 1.771e-2, 0.5 percent, and `||R||` from 1.203 to 1.078;
every one of those iterations was a damped Gauss-Newton escape at
`mu = 1e-4` to `1e+2` with a full step `lam = 1`. The solve then reported no
descent direction at all along either the Newton or the damped Gauss-Newton
leg and aborted. The valley floor is at a cellwise maximum near 1.1, not at
the tolerance.

### T0b. Where the worst row sits along S3's pass 1

| iterations | worst row | `\|\|R\|\|` | `\|\|Fs\|\|2` |
|---|---|---|---|
| start | mass, cell 150 (r = 1.0615) | 1.829 | 3.05e+1 |
| 1 to 3 | momentum/energy, cells 2 to 3 (r = 1.000 to 1.001) | 1.69 to 2.00 | 8.6e-1 to 5.1e+0 |
| 4 to 40 | energy, cells 3 to 5 | 1.26 to 1.99 | 4.4e-2 to 2.7e+0 |
| 50 to 400 | **mass, cell 3** (r = 1.001) | 1.20 to 1.89 | 1.86e-2 to 1.78e-2 |
| 448 to 604 | momentum, cell 500 (r = 29.03) alternating with mass, cell 3 | 1.08 to 1.18 | 1.771e-2 |

The solve spends its whole life in cells 2 to 5, which is where section 3
measures the discrete energy flux to be four decades above the heating. The
preconditioned Krylov leg reaches its `rtol = 1e-1` in one matrix-vector
product at 598 of the 604 iterations; the obstruction is not the linear
solve.

### T0c. The same pass with `Well balanced: True`

MEASURED, `<scratch>/L5c_rt_S3_wb/run.log`, the same seed, the same input
file plus the one key:

| iteration | 1 | 22 | 26 | 29 | 30 | 31 | 32 | 33 | 34 |
|---|---|---|---|---|---|---|---|---|---|
| `\|\|R\|\|` | 1.947 | 1.125 | 5.358e-1 | 9.185e-2 | 7.406e-3 | 3.766e-5 | 5.061e-7 | 2.658e-8 | 4.622e-9 |
| `\|\|Fs\|\|2` | 1.47 | 2.58e-1 | 3.66e-3 | 2.58e-4 | 1.40e-5 | 2.92e-7 | 3.48e-9 | 2.88e-10 | 5.08e-11 |

`done info=0 ||R|| = 4.621E-09, flux spread = 1.406E-12`, 34 iterations,
85 s. From iteration 29 the convergence is quadratic, which is what a Newton
solve on a residual whose terms are the same size does.

---

## 3. The discrete energy flux divergence of the base cells

### 3.1 What is being decomposed

`RK_rhs` assembles the energy row's flux divergence as

```
dF(3,j) = [ A+ F+(3) - A- F-(3)
          + A+ F+(1) (Gphi_i(j)   - Gphi_c(j))
          - A- F-(1) (Gphi_i(j-1) - Gphi_c(j)) ] / dV
```

(`src/modules/time_step/RK_rhs.f90`, READ), so the potential rides on the
MASS flux and the face flux `F(3)` is the total-energy flux `v (E + p)`.
Under HLLC in a subsonic cell the face flux is exactly

```
F(3) = S* (E* + p*)                     (MEASURED to 2.4e-9 relative)
S*   = (p_R - p_L + phi_L v_L - phi_R v_R)/(phi_L - phi_R)
     = S*_vel + S*_pres ,
S*_vel = (phi_L v_L - phi_R v_R)/(phi_L - phi_R)     the velocity average
S*_pres = (p_R - p_L)/(phi_L - phi_R)                the pressure jump
```

with `phi_K = rho_K (S_K - v_K)`, so at low Mach
`S*_pres ~ -(p_R - p_L)/((rho_L + rho_R) c)` and, as a Mach number,
`S*_pres/c ~ -(p_R - p_L)/(2 gamma p)`. Splitting `F(3)` and the mass flux
`F(1) = rho* S*` the same way and taking the divergence of each part gives
the table below.

### T1. The energy row of the S3 seed, term by term

Units: `(heat - cool)` of the same cell, which is what a stationary energy
balance has to equal. `S_3 = 0` and conduction is off, so
`R_3 = dF_3 - (heat - cool)` exactly.

| cell | r [R_p] | `dF_3`/(heat-cool) | velocity part | pressure-jump part | kinetic part of `F_3` |
|---|---|---|---|---|---|
| 1 | 1.00019 | +4.155e+04 | -1.582e+01 | +4.157e+04 | 2.9e-07 |
| 2 | 1.00039 | -4.612e+04 | +3.382e+00 | -4.612e+04 | 8.4e-08 |
| 3 | 1.00058 | +1.607e+04 | +1.029e+00 | +1.607e+04 | 1.5e-11 |
| 5 | 1.00097 | +2.744e+03 | +1.122e+00 | +2.743e+03 | 1.7e-10 |
| 10 | 1.00193 | -8.301e+02 | +1.057e+00 | -8.312e+02 | 6.5e-10 |
| 20 | 1.00387 | -5.157e+01 | +1.181e+00 | -5.275e+01 | 7.1e-12 |
| 41 | 1.00793 | +1.698e+00 | +1.456e+00 | +2.418e-01 | 2.3e-13 |
| 68 | 1.01372 | -1.697e+00 | +1.667e+00 | -3.363e+00 | 2.2e-12 |
| 100 | 1.02496 | +1.857e+00 | +1.727e+00 | +1.305e-01 | 1.5e-11 |
| 150 | 1.06151 | -1.857e+03 | +2.530e+00 | -1.859e+03 | 3.4e-08 |
| 200 | 1.14888 | +1.977e+00 | +1.936e+00 | +4.067e-02 | 8.3e-09 |
| 300 | 1.85709 | +5.311e-01 | +5.325e-01 | -1.415e-03 | 1.5e-06 |

Read it as: the column the velocity carries is 1.0 to 2.5 everywhere, which
is the stationary balance (L3's Bernoulli mismatch is that departure from
1); the column the pressure jump carries is 10^2 to 10^4 in the base layer
and falls below the velocity part at cell ~41, r = 1.008 R_p. Cell 150 is
the seam of S3's own construction (L5a section 5) and cell 68 is where the
two parts are of the same size and opposite sign, which is why L5a's T5 read
`-1.0` there.

**The kinetic energy plays no part.** `S* rho* S*^2/2` is 1e-7 to 1e-16 of
the face energy flux everywhere, which is `M^2` as it must be. The
decomposition has two terms, not three.

### T2. The face, in Mach numbers

| face | r_edg [R_p] | `S*_vel/c` (the flow) | `S*_pres/c` (the jump) | ratio | `rho v r^2/F_wind` of the face |
|---|---|---|---|---|---|
| 1 | 1.00029 | 3.369e-07 | +9.462e-04 | 2808 | +2.290e+03 |
| 2 | 1.00048 | 4.465e-07 | -5.106e-04 | -1143 | -1.142e+03 |
| 3 | 1.00068 | 4.796e-07 | -7.196e-06 | -15.0 | -1.400e+01 |
| 5 | 1.00106 | 5.480e-07 | -2.351e-05 | -42.9 | -4.188e+01 |
| 10 | 1.00203 | 7.196e-07 | -4.532e-05 | -63.0 | -6.195e+01 |
| 20 | 1.00396 | 1.101e-06 | -5.731e-06 | -5.20 | -4.203e+00 |
| 30 | 1.00590 | 1.541e-06 | -2.961e-06 | -1.92 | -9.213e-01 |
| 50 | 1.00976 | 2.574e-06 | -2.625e-06 | -1.02 | -1.995e-02 |
| 68 | 1.01385 | 3.864e-06 | -1.308e-06 | -0.34 | +6.614e-01 |
| 100 | 1.02519 | 8.383e-06 | -1.707e-06 | -0.20 | +7.963e-01 |
| 200 | 1.15020 | 1.623e-04 | -3.982e-06 | -0.02 | +9.759e-01 |
| 300 | 1.86464 | 2.121e-03 | -1.699e-05 | -0.008 | +9.936e-01 |

The crossover is at face ~50, **r = 1.0098 R_p**: below it the scheme's
contact wave does not move at the flow's velocity, and the mass flux it
transports has nothing to do with `F_wind`. Above r = 1.15 the face carries
the wind to within 2 percent, which is why every gate anchored at 1.2 R_p
passed on states whose base layer was this far out.

### 3.2 Which of the three candidates it is

**(a) The reconstruction: yes, but as its design-order truncation error and
not as an odd-even mode.**

### T3. The face pressure jump `(p_R - p_L)/p`, against two smoothed versions of the same column

| face | r_edg | as loaded | after eight 1-2-1 passes on `ln rho, ln v, ln p` | on a ninth-order fit of `ln p, ln rho` |
|---|---|---|---|---|
| 1 | 1.00029 | -3.011e-03 | -9.762e-05 | +6.945e-05 |
| 2 | 1.00048 | +1.636e-03 | +1.883e-06 | +6.410e-05 |
| 3 | 1.00068 | +2.313e-05 | +9.094e-05 | +5.922e-05 |
| 5 | 1.00106 | +7.584e-05 | +1.352e-04 | +5.066e-05 |
| 10 | 1.00203 | +1.473e-04 | +4.281e-05 | +3.476e-05 |
| 20 | 1.00396 | +1.874e-05 | +1.612e-05 | +1.732e-05 |
| 30 | 1.00590 | +9.717e-06 | +8.764e-06 | +9.247e-06 |
| 68 | 1.01385 | +4.309e-06 | +4.157e-06 | +3.980e-06 |
| 100 | 1.02519 | +5.616e-06 | +5.601e-06 | +5.495e-06 |

Removing every 2 dr component of the state leaves the jump where it was from
face 3 outward, and an exactly smooth analytic column reproduces it within a
factor 4 at every face beyond the third. Only faces 1 and 2 are dominated by
the seed's own 2 dr content (5-point quadratic residual of `ln rho`:
-1.00e-2 at cell 1, -1.29e-3 at cell 2, +1.03e-3 at cell 3, -2.1e-5 at cell 5),
and that content is the base boundary's, which is P44.

### T4. The same jump under refinement of the analytic column

The analytic column of T3's last column, sampled on uniform grids of spacing `h`, the
jump read at the face nearest a fixed radius. The production base spacing is
`dr_base = 2.0e-4 R_p` (the default, this run states no `Base grid` key).

| r | h = 2.0e-4 | 1.0e-4 | 5.0e-5 | 2.5e-5 | order |
|---|---|---|---|---|---|
| 1.0010 `(p_R-p_L)/p` | 5.745e-05 | 7.314e-06 | 9.210e-07 | 1.155e-07 | 2.97, 2.99, 3.00 |
| 1.0010 `\|S*_pres/S*_vel\|` | **33.5** | 4.22 | 0.53 | 0.07 | |
| 1.0050 `(p_R-p_L)/p` | 1.358e-05 | 1.711e-06 | 2.146e-07 | 2.687e-08 | 2.99, 2.99, 3.00 |
| 1.0050 `\|S*_pres/S*_vel\|` | **3.10** | 0.39 | 0.05 | 0.01 | |
| 1.0140 `\|S*_pres/S*_vel\|` | 0.13 | 0.02 | 0.00 | 0.00 | |

**Order 2.97 to 3.00: it is the WENO3 reconstruction's own truncation error,
and it converges away.** A 4x refinement of the base spacing
(`Base grid [dr,cells]: 5.0e-5 400`) would put the spurious contact speed
below the physical one at every radius above 1.001 R_p. This is the same
statement P44 section 7 made from the other end for the velocity artifact
("first order in `dr` and converges away"), now for the quantity that
actually obstructs the solve.

**(b) HLLC at Mach 1e-6: it is the converter, not the source.** The flux
does what a contact-resolving flux is built to do -- it resolves the middle
wave exactly and moves the entropy jump at `S*` -- and `S*` is a solution of
the face Riemann problem, so it must respond to the face pressure jump. Any
flux with that property (ROE has it too, `input_schema.md` K46 core line 17,
READ) converts the jump into a contact speed the same way. What makes it
fatal here is only that the flow's own Mach number is 3e-7.

**(c) The characteristic base condition: no, beyond the first two faces.**
`Rec_BC` replaces the left state of face 0 only (`WL_out(:,0) = base_face_W`,
`Apply_BC.f90` line 343, READ), so face 0 reads the boundary, face 1 reads
the ghost through the WENO3 stencil of cell 1, and faces 2 and beyond touch
no ghost at all. The pattern of T2 runs unchanged out to face 50. What the
boundary does own is the 2 dr content of cells 1 to 3 in T3, which is the
P44 defect and is a second-order contribution to the same quantity.

### 3.3 Why the row cannot be driven down, quantified

### T5. Sensitivity of `dF_3(j)` to a relative perturbation at cell j

In units of `(heat - cool)(j)`, one-sided difference at `eps = 1e-8`.

| cell | r [R_p] | `d dF_3 / (eps dp/p)` | `/ (eps drho/rho)` | `/ (eps dv/v)` |
|---|---|---|---|---|
| 1 | 1.00019 | +9.526e+06 | -4.608e+03 | +6.465e-02 |
| 5 | 1.00097 | +9.917e+06 | -4.575e+02 | +5.400e-01 |
| 20 | 1.00387 | +9.929e+06 | -2.399e+01 | +5.743e-02 |
| 68 | 1.01372 | +6.922e+06 | -5.283e+00 | +9.727e-01 |
| 150 | 1.06151 | +8.067e+05 | +7.759e+02 | +8.618e-01 |
| 300 | 1.85709 | +5.783e+03 | +5.136e-01 | +1.033e+00 |
| 400 | 5.90434 | +2.295e+02 | +4.530e-01 | +2.048e+00 |

The pressure sensitivity is of order `1/M` (measured 9.5e6 where
`M = 3.4e-7`, whose reciprocal is 2.9e6) and the velocity sensitivity is 1,
which is the ratio of the acoustic to the advective part of an upwind energy
flux. A solver asked to hold the row
below the run's own `Resid tol: 5.0e-5` with this flux has to place the base
pressure field to 5e-12 relative, in a residual assembled from differences of
O(1) face pressures. That is the "no descent direction" of T0 stated as a
condition number.

---

## 4. What the well-balanced option changes

Under `Well balanced: True` the pressure component is reconstructed in the
coordinate of the departure from the cell's own constant-density hydrostatic
equilibrium and the Riemann jump is
`dp_wb = dp_eq + dev_R - dev_L`, where `dp_eq` is the mismatch of the two
neighboring equilibria at the shared face
(`Reconstruction.f90`, `Num_Fluxes.f90`, READ).

### T6. The pressure jump under the two operators

"eq column" is the same T and the same base cell, with `rho` and `p`
re-integrated so that `dp_eq = 0` cell by cell, i.e. the discrete
hydrostatic equilibrium of the well-balanced operator itself.

| face | S3: `dp_eq/p` | S3: `dp_wb/p` | S3: `S*_pres/c` (WB) | eq column: `dp_eq/p` | eq column: `S*_pres/c` (WB) | eq column: `S*_pres/c` (plain) | `S*_vel/c` |
|---|---|---|---|---|---|---|---|
| 1 | -1.018e-02 | -9.409e-03 | +2.952e-03 | +1.0e-14 | -1.9e-06 | -8.3e-05 | 4.100e-07 |
| 3 | -5.705e-04 | -1.240e-04 | +3.860e-05 | +1.9e-15 | +7.8e-16 | -3.0e-05 | 4.738e-07 |
| 5 | +1.641e-04 | +1.175e-05 | -3.643e-06 | -8.7e-16 | +2.5e-16 | -2.9e-05 | 5.417e-07 |
| 10 | +1.337e-04 | +1.021e-04 | -3.142e-05 | +8.1e-15 | -9.9e-16 | -1.2e-05 | 7.101e-07 |
| 20 | -1.570e-04 | +2.750e-06 | -8.410e-07 | -2.4e-15 | +1.1e-15 | -4.9e-06 | 1.085e-06 |
| 30 | -1.290e-04 | +1.032e-06 | -3.146e-07 | -4.4e-15 | +8.3e-16 | -2.6e-06 | 1.516e-06 |
| 50 | -9.982e-05 | +2.282e-07 | -6.923e-08 | +5.1e-15 | -5.5e-16 | -2.6e-06 | 2.526e-06 |
| 68 | -1.071e-04 | +1.764e-07 | -5.355e-08 | -1.2e-15 | +1.6e-16 | -1.2e-06 | 3.785e-06 |
| 100 | -1.285e-04 | +5.514e-08 | -1.676e-08 | -2.1e-15 | +1.6e-16 | -1.7e-06 | 8.179e-06 |

Two readings.

1. **On a state the operator recognizes as its own equilibrium the spurious
   contact speed is zero to round-off**, while the physical one is
   untouched. The plain flux on the identical column still reads 202 times
   the physical Mach number at face 1, 64 at face 3, 16 at face 10, 4.5 at
   face 20, and reaches 1.0 only at face 50.
2. **The S3 seed is not that equilibrium**: `dp_eq/p` is 1e-2 at face 1 and
   1e-4 through the layer, which as a Mach number is 3e-5, still 60 times
   the flow's. That is why the seed's own residual is WORSE under the
   well-balanced option (section 5) and the solve nonetheless converges from
   it: the option does not make the entry state better, it makes the target
   exist and the Newton model point at it.

---

## 5. The two keys, measured

`EXHALE_RESIDUAL=1` on the S3 seed, four runs differing by one line of
`input.inp` (`<scratch>/L5c_res_S3_{base,wb,lm2e-2,lm5e-3}`). The `base` row
reproduces L5a's S3 row of T1/T2 to every printed digit, which is the
control.

### T7. The residual of the S3 seed

| option | mass max @ cell | momentum max @ cell | energy max @ cell | mass r<1.03 num/den | energy r<1.03 num/den |
|---|---|---|---|---|---|
| none (control) | 1.8287 @ 150 | 0.2915 @ 499 | 1.5894 @ 68 | 7.061e-3 / 5.030e-3 | 7.361e-3 / 7.363e-3 |
| `Well balanced: True` | 1.8940 @ 6 | 0.4333 @ 499 | 1.9177 @ 41 | 1.126e-2 / 1.150e-2 | 1.169e-2 / 1.169e-2 |
| `Low-Mach damping: 2.0e-2 1.0e-3` | 1.8287 @ 150 | 0.2915 @ 499 | 1.5894 @ 68 | 7.061e-3 / 5.030e-3 | 7.361e-3 / 7.363e-3 |
| `Low-Mach damping: 5.0e-3 1.0e-3` | 1.8287 @ 150 | 0.2915 @ 499 | 1.5894 @ 68 | 7.061e-3 / 5.030e-3 | 7.361e-3 / 7.363e-3 |

The two coefficients are the ones the documents recommend:
`docs/p44_base_sawtooth.md` section 6 names `2.0e-2` as the best of the four
prescriptions for the interior part of the hot-Uranus pattern, and
`docs/hd209_metal_stagnation.md` section 9.2 scans `5.0e-3`, `1.0e-2`,
`2.0e-2`. `M_th` is the default 1e-3, and the base Mach of 3e-7 leaves the
gate `[max(0, 1 - M^2/M_th^2)]^2` fully open, so the term is armed
everywhere in the layer.

**Cell by cell, `Low-Mach damping` changes the mass row by exactly zero (it
does not touch the mass flux, by construction), the momentum row by at most
2.1e-4 relative (cell 148) and the energy row by at most 2.3e-7 relative
(cell 204), at `eps4 = 2.0e-2`; at `5.0e-3` the two are 5.1e-5 and 5.4e-8.**
It is not that the gate is shut; it is that a fourth difference of a smooth
stationary velocity field is nothing.

### T8. Outer pass 1 of the partitioned route, same seed, 16 threads, `EXHALE_PTC_DTAU0=1.0`

| option | hydro | mass | momentum | energy | gated species row | s | outcome |
|---|---|---|---|---|---|---|---|
| none (L5a control, re-READ) | info=2 | 1.08E+00 | 3.24E-02 | 6.26E-01 | 1.69e-3 of 1e-5 | 1325 | REFUSED, `element_step_out_of_bounds` at step 25 |
| `Well balanced: True` | **info=0** | **2.55E-10** | **1.77E-14** | **4.62E-09** | 1.15e-3 of 1e-5 | 85 | element fixed point in 32 steps; pass 2 |
| `Low-Mach damping: 2.0e-2` | info=2 | 1.32E+00 | 1.35E-02 | 1.24E+00 | 1.73e-3 of 1e-5 | 1232 | REFUSED, `element_step_out_of_bounds` at step 21 |

The `Low-Mach damping` run reproduces the control's outcome in every
respect: 543 JFNK iterations, the same `no descent direction exists for the
banded model at this state` abort, `||R||` 1.324, and the same refusal by the
element relaxation at pass 1. The key changes nothing about this obstruction.

### T9. The well-balanced route to the end

| outer pass | hydro | mass | momentum | energy | gated element row | element relaxation | s |
|---|---|---|---|---|---|---|---|
| 1 | info=0 | 2.55e-10 | 1.77e-14 | 4.62e-09 | 1.15e-03 | fixed point, 32 steps, drift 2.46e-01 | 85.3 |
| 2 | info=0 | 5.26e-10 | 5.16e-13 | 7.30e-09 | 5.15e-04 | fixed point, 32 steps, drift 1.23e-01 | 74.4 |
| 3 | info=0 | 4.18e-10 | 1.31e-13 | 6.21e-09 | 2.37e-04 | fixed point, 31 steps, drift 6.21e-02 | 58.4 |
| 4 | info=0 | 4.19e-10 | 1.81e-13 | 1.40e-08 | 1.11e-04 | fixed point | 54.6 |
| 5 | info=0 | 2.39e-10 | 1.03e-13 | 5.55e-09 | 5.43e-05 | fixed point | 53.4 |
| 6 | info=0 | 2.70e-10 | 1.49e-12 | 6.87e-09 | 2.88e-05 | fixed point | 50.9 |
| 7 | info=0 | 3.93e-10 | 4.05e-14 | 6.61e-09 | 1.58e-05 | fixed point | 55.2 |
| 8 | info=0 | 4.53e-10 | 2.29e-13 | 5.07e-09 | **8.76e-06** of 1e-5 | fixed point | 49.1 |

### T9b. The same route from the MAPPED seed, which is the entry of the plan

The mapped seed is the archived state on the current grid, the one that
carries the P44 sawtooth (`v(1) = -538 cm/s`) and the 2.92 percent base
pressure offset L5a found -- none of L1's, L3's or L5a's repairs applied.
With `Well balanced: True` it runs the SAME ladder and certifies at the same
pass: gated element row 1.15e-3, 5.15e-4, 2.37e-4, 1.11e-4, 5.43e-5,
2.88e-5, 1.58e-5, **8.76e-6** at passes 1 to 8, `ACCEPTED`, `CERTIFIED`,
flux gate 1.4087e-11, `info = 0`, hydro rows at pass 8
3.52e-10 / 2.13e-13 / 5.34e-09.

**S1c does the same**: the anchored velocity seed certifies at outer pass 8
on the identical element-row ladder, flux gate 1.4953e-11.

**All three certified states are the same solution.** MEASURED over all 504
rows of `Hydro_ioniz.txt`:

| pair | max rel `rho` | max rel `v` | max rel `p` | max rel `T` |
|---|---|---|---|---|
| S3 vs mapped | 3.4e-10 | 2.5e-08 | 5.0e-11 | 3.5e-10 |
| S3 vs S1c | 1.7e-10 | 9.6e-09 | 2.2e-11 | 1.8e-10 |

and `Mdot` agrees to seven digits between S3 and the mapped seed
(1.94830e8 g/s both). Seeds that differ by four decades in the base velocity
land on the same fixed point, so this is the solver's root and not a seed
artifact -- and the seed engineering of L1, L3 and L5a turns out not to have
been what the solve needed.

`outer pass 8: ACCEPTED -- every active equation of this state is within its
own tolerance`, then
`CERTIFIED: every active equation was evaluated and is within its
tolerance`, `certified IN THE WIND, r >= 1.20`; flux gate 1.4922e-11,
`flux spread by window: r>=1.03 5.97e-11, r>=1.10 3.14e-11`. The element row
halves at every pass, which is the HD 209458 b control's behaviour exactly.

### T10. What the key does to a case that already converged

The HD 209458 b element-diffusion reload of the plan's section 1, the
control of the whole route, run with `Well balanced: True` and nothing else
changed (`<scratch>/L5c_hd209_wb` against `<scratch>/ctrl_hd209_elem`).

| | control | with `Well balanced: True` |
|---|---|---|
| certified at | outer pass 12 | **outer pass 12** |
| hydro mass / momentum / energy | 1.15e-09 / 1.66e-14 / 1.28e-08 | 1.16e-09 / 3.78e-14 / 9.05e-09 |
| worst gated species row | 7.73e-06 at cell 289 | **7.73e-06 at cell 289** |
| `Mdot` [g/s] | 2.65688e10 (log10 10.4244) | 2.65724e10 (log10 10.4244) |

`Mdot` moves by **1.3e-4** and the two solutions agree to 2.3e-4 in `T`,
1.1e-4 in `v` and 2.8e-4 in `rho` beyond 2 R_p; the departures grow inward
(1.8e-3 in `v` beyond 1.2 R_p, 1.8e-2 beyond 1.1) and are confined to the
base: 12 percent in `rho` and 14 percent in `T` at cell 1, where the
option's whole point lies. The cell-1 velocity artifact goes from
-0.968 cm/s to +0.150 cm/s and the first six cells from
(-0.97, -1.11, -0.41, -0.13, +0.12, +0.23) to
(+0.15, -0.15, +0.13, +0.06, +0.09, +0.08) cm/s.

**The key therefore does what it is advertised to do on a case that already
worked: it rewrites the base cells and leaves the wind alone.** It is not
byte-identical and cannot be -- K46 says so -- so adopting it as a default
would be a deliberate golden refresh, not a free change.

---

## 6. Judgment

**The base layer's discrete energy residual is a property of the scheme.**
The state can be made as stationary as one likes -- L5a's S3 has the exact
continuity velocity, a Bernoulli-consistent temperature and an anchored base
face -- and the energy row still enters at O(1), because the quantity the
row differences is not the enthalpy flux of that state. At cell 2 the part
the flow's own velocity carries is 7.3e-5 of `dF_3` and all the rest is an
enthalpy flux carried at a contact speed the reconstruction's third-order
truncation error puts there; `|pressure part|/|dF_3|` is above 0.975 in
every one of the 20 cells inside r = 1.004 R_p. No seed can remove a term
that is a functional of the discretization and the grid.

**The minimal change is a key that already exists and needs no source
edit: `Well balanced: True`.** It replaces the O(1) pressure jump of the
Riemann problem by the departure from each cell's own hydrostatic
equilibrium, which is the quantity that is actually small in this layer; the
spurious contact speed then vanishes with it (T6) and the Newton solve
converges quadratically (T0c) to a certified state (T9). The reproduction
condition is one line:

```
Well balanced: True
```

added to `<scratch>/t_map_stat/input.inp`, with `OMP_NUM_THREADS=16
EXHALE_PTC_DTAU0=1.0`. **Which seed is in `output/` does not matter**: the
mapped, S1c and S3 seeds all certify at outer pass 8 on the same element-row
ladder and land on the same state to 1e-8 (T9b), so the LHS 1140 b recipe
L4 needs is the archived state mapped onto the grid plus this one line.

**The alternative, if the plain flux has to be kept, is base refinement.**
T4 gives the exchange rate: the jump is third order, and
`Base grid [dr,cells]: 5.0e-5 400` (4x the default spacing, the same 0.01
R_p of base) takes the spurious contact speed below the physical one above
r = 1.001 R_p. That costs 350 extra cells and does nothing for cell 1.

**`Low-Mach damping` is the wrong instrument for this layer** and its being
inert here is not a failure of the key: it damps a 2 dr VELOCITY mode, and
what obstructs this wind is a smooth PRESSURE truncation error. It should
not be carried into the LHS 1140 b recipe.

---

## 7. What L4 has to check before the certified state is used

Reported, not resolved; all MEASURED on
`<scratch>/L5c_rt_S3_wb/output/Hydro_ioniz.txt`.

1. **`Mdot` moves by a factor 3.1.** The certified state gives
   1.948e8 g/s (log10 8.29) against the archived 6.379e7 (log10 7.81) and
   the S3 seed's 8.888e7 (7.95). L5a section 5 put the mass flux the base
   layer's energy budget demands at 1.10 to 1.95 times the archived,
   depending on how far out the integral runs; 3.1 is above that band, and
   the heating moves with the column, so the two are not directly
   comparable. It is not a seed artifact: the mapped seed gives the same
   `Mdot` to seven digits (T9b). The profile comparison the user's rule asks
   for -- T, v, rho and the composition overlaid against the archived state
   and against Koskinen-style expectations -- has not been made here.
2. **A three-cell base velocity artifact survives**, much reduced. The
   certified state carries `v` = -0.092, -0.008, +0.200, +0.186, +0.223,
   +0.234 cm/s in cells 1 to 6, against the stationary 0.150 cm/s that its
   own `F_wind` and cell-1 density give, so the cell-centered
   `rho v r^2/F_wind` reads -0.614, -0.046, +1.085, +0.923, +1.020, +0.992
   and is 1.000 to four digits from cell 8 outward. On the archived state
   the same ratio was -1.8e4 at cell 1 and took to cell 233 to reach the
   band. The residual rows are within tolerance throughout (mass 4.5e-10 at
   cell 3), so the faces carry a constant flux and only the cell-centered
   product does not: P44 section 3's reading, now at 1.6 times the signal
   rather than 1.8e4 times it. The base ghost carries -2.79 cm/s,
   -20.6 F_wind.
3. **The base thermal jump is unchanged in character**: cell 1 sits at 405 K
   against the 226 K reservoir at `r = 1`. That is P44 section 8's
   independent defect and nothing here touches it.
4. **Grid convergence of the well-balanced solution has not been measured.**
   The discrete equilibrium the option preserves is a second-order
   approximation of the hydrostatic column, and the column that satisfies it
   exactly differs from a continuum-hydrostatic one by 1.2e-2 at cell 5 and
   1.2e-1 at cell 300 on this grid (section 4's "eq column"). Whether the
   certified `Mdot` is grid converged is the first L4 question.
5. **Turning the key on changes every other run**, though by far less than
   this planet's base needs. `input_schema.md` K46 says the option moves
   every result, and the goldens were taken with it off. Measured on the
   HD 209458 b element reload (T10): the case still certifies at outer pass
   12 with the same worst species row, `Mdot` moves by 1.3e-4, the wind
   beyond 2 R_p by 2 to 3e-4, and the base cells by 12 to 14 percent.
   Adopting the key beyond the LHS 1140 b recipe is therefore a deliberate
   golden refresh and a user decision, not a side effect.

---

## 8. Method and reproduction

`EX` is the repository, `S` the session scratch directory
`.../scratchpad/L`.

```
# the residual of the S3 seed under each key
for tag in base wb lm2e-2 lm5e-3; do
  mkdir -p $S/L5c_res_S3_$tag/output
  \cp -f $S/t_atomic/input.inp $S/L5c_res_S3_$tag/input.inp
  \cp -f $S/seeds/S3/*.txt     $S/L5c_res_S3_$tag/output/
done
printf 'Well balanced: True\n'               >> $S/L5c_res_S3_wb/input.inp
printf 'Low-Mach damping: 2.0e-2 1.0e-3\n'   >> $S/L5c_res_S3_lm2e-2/input.inp
printf 'Low-Mach damping: 5.0e-3 1.0e-3\n'   >> $S/L5c_res_S3_lm5e-3/input.inp
cd $S/L5c_res_S3_<tag> && OMP_NUM_THREADS=4 EXHALE_RESIDUAL=1 $S/L5a_bin/EXHALE.x > run.log

# the state with its ghosts, for the face decomposition
mkdir -p $S/L5c_pp_S3_base/output
sed 's/^Do only PP: False/Do only PP: True/' $S/t_atomic/input.inp \
    > $S/L5c_pp_S3_base/input.inp
printf 'CFL: 1.0e-12\n' >> $S/L5c_pp_S3_base/input.inp
\cp -f $S/seeds/S3/*.txt $S/L5c_pp_S3_base/output/
cd $S/L5c_pp_S3_base && OMP_NUM_THREADS=4 $S/L5a_bin/EXHALE.x > run.log

# the partitioned route
mkdir -p $S/L5c_rt_S3_wb/output
\cp -f $S/t_map_stat/input.inp $S/L5c_rt_S3_wb/input.inp
printf 'Well balanced: True\n' >> $S/L5c_rt_S3_wb/input.inp
\cp -f $S/seeds/S3/*.txt $S/L5c_rt_S3_wb/output/
cd $S/L5c_rt_S3_wb && OMP_NUM_THREADS=16 EXHALE_PTC_DTAU0=1.0 $S/L5a_bin/EXHALE.x > run.log
```

**The decomposition is a reimplementation of the code's own operator,
validated against the code before it was read.** `$S/l5c_scheme.py` rebuilds
the WENO3 reconstruction (`Reconstruction.f90` case `WENO3`, including
`weno3_geometry_coefficients` and the `Rec_BC` upper extrapolation), the
HLLC flux (`Num_Fluxes.f90` case `HLLC`) and the flux difference of
`RK_rhs`, on the state the run writes -- the two ghost cells at each end
included, so no boundary condition has to be rebuilt. The grid is rebuilt
from the written cell centres by the code's own identities
(`r_edg(j) = (r(j)+r(j+1))/2`, `dr_j(j) = r_edg(j)-r_edg(j-1)`,
`define_grid.f90`, READ) and the potential from `phi = -b0/r` with
`b0 = Gc Mp mu/(kb T0 R0)` (`grav_field.f90`, `input_read.f90` line 2050,
READ). The gas is atomic here, so `caloric_mixture_active` is false and the
EOS is `gamma_ad = 5/3` exactly (`caloric_eos.f90`, READ).

Two checks, both against the binary's own output and both passed before any
number above was read off:

- **the mass row**, which has no source, so `R_mass(j) = dF(1,j)` exactly:
  the reimplementation reproduces `output/residual_profile.txt` to a median
  relative difference of **2.7e-6** and a maximum of 3.3e-6 over cells 2 to
  500, the floor being the 9 significant digits the output file carries;
- **the energy row**, rebuilt as `R_3 + (heat - cool)` from the same file
  and the `Do only PP` pass: median **2.6e-6**, maximum 7.2e-4 over the same
  cells.

Cell 1 is excluded from both: its lower face takes `base_face_W` from the
characteristic boundary and not the reconstruction, which is exactly the
17.3 percent difference the check shows there and is the sign that the check
is measuring what it claims to. The code's own base face flux is recovered
from `R_mass(1)` and the reproduced face-1 flux instead
(section 1, -1746 F_wind).

The same script carries the well-balanced operator
(`hydrostatic_equilibrium_of_each_cell`, the departure stencil,
`well_balanced_face_departures` and the `dp_wb` form of `S*`), validated the
same way against `<scratch>/L5c_res_S3_wb`: median 2.7e-6, maximum 6.5e-5.

Code units, as in L1, L3 and L5a: `n_0 = 3.204863e13 cm^-3` (the run's own
`derived n0` line), `T_0 = 226.0 K`, `v_0 = sqrt(kb T0/mu) = 1.365459e5
cm/s`, `p_0 = n_0 mu v_0^2 = 1.000 dyn cm^-2`, volumetric rate unit
`q_0 = n_0 mu v_0^3/R_0 = 1.211192e-4 erg cm^-3 s^-1`, `R_0 = 1.127370e9
cm`, `b_0 = 106.2132`. `F_wind = 7.551987e-7` in code units, the median of
`rho v r^2` over cells 1 to 150 of the S3 seed.

Run directories: `$S/L5c_res_S3_{base,wb,lm2e-2,lm5e-3}` (residual),
`$S/L5c_pp_S3_{base,wb}` (the state with ghosts),
`$S/L5c_rt_S3_{wb,lm2e-2}`, `$S/L5c_rt_{S1c,map}_wb`, `$S/L5c_hd209_wb`
(the route; the plain-flux control of T8 is L5a's `$S/L5a_S3`). Scripts: `$S/l5c_scheme.py`, `$S/l5c_decomp.py`, and the
consolidated table dump `$S/l5c_numbers.txt`.

---

## 9. Observations outside the item, reported and not fixed

1. **`docs/lhs1140b_stationary_L5a_20260913.md` table T3 is labelled "the
   base face, which is where the O(1) residual of the base cells lives" and
   its second column is the ghost's `rho v r^2`, not the flux the scheme
   transports through that face.** On the S3 seed the two differ in sign and
   by a factor 476 in magnitude (+3.67 F_wind against -1746 F_wind). The table's numbers are right
   for what they measure -- the boundary condition's own face state -- and
   the heading over-reads them. L5a's conclusion that the base pressure
   anchor is what moves the integrated norms is unaffected.
2. **The `Low-Mach damping` entry in `docs/input_schema.md` (K31b) and in
   `docs/hd209_metal_stagnation.md` section 8 states the gate and the
   stability bound but not the one property that decides whether the key can
   help: the stress is a fourth difference of the VELOCITY, so it is blind
   to a pressure-reconstruction error however large.** Adding that sentence
   at K31b would have saved this item a run. Not edited, because the key's
   own file header (`low_mach_dissipation.f90` section 2) does say the mass
   and energy fluxes are untouched and the schema is a summary of it.
