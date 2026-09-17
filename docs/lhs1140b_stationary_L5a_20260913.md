# LHS 1140 b, item L5a: a seed with the stationary velocity field

Item L5a of `docs/PLAN_20260913_lhs_stationary.md`, the experiment L1 and L3
left: rebuild the seed from the archived thermodynamics with the velocity a
stationary wind carries, `v = F_wind/(rho r^2)`, and give it to the
partitioned stationary route.

Every number here is MEASURED on the binary of 2026-09-13 03:20 (git 43bc28c,
tree dirty, md5 97e10317a710b9ccc63addbedde3586a, the copy used for every run
of this item is `<scratch>/L5a_bin/EXHALE.x`) unless it is marked READ. No
source file was changed. The seed builder is
`LHS1140b/models/make_stationary_seed.py`.

## 1. Verdict

**The route does not accept any of the seeds: every one ends at outer pass 1
with `hydro info = 2` and the pass REFUSED, exactly as the mapped seed does.**
What the experiment settles is why, and it is not what item L5 supposed.

- **The stationary velocity field is not enough, and it is not where most of
  the base-layer residual was.** Setting `v = F_wind/(rho r^2)` below
  1.266 R_p takes the cell-centered mass flux of the base layer from 1.8e4
  times the wind's value to 0.998 of it in every cell (table T4), and takes
  the base layer's ratio of its own energy flux divergence to its net
  heating,
  `|dF_3|/|heat - cool|`, from a median of 15.6 over cells 1 to 150 to 1.00,
  with the sign of the continuum value in 104 of those 150 cells against 28
  before. It lowers the integrated base-layer numerators by 1.25 (mass) and
  1.27 (energy) and leaves every cellwise row maximum at O(1).
- **The larger defect the velocity field was hiding is the base PRESSURE.**
  The archived column stands 2.92 percent above the 1 microbar reservoir at
  the base face (MEASURED: the column scaling that nulls the base face
  velocity is 0.971649, whose reciprocal is 1.02918), and the characteristic
  base condition turns that into an inflow the interior velocity has no say
  in: `v_b = v_i + (p_b - p_i)/(rho_i c_i)`
  (`src/modules/states/base_boundary.f90` line 352, READ). A relative
  pressure offset d therefore enters as a face Mach number of about d/gamma:
  2.92 percent of pressure predicts Mach 1.75e-2, and the run's own boundary
  report prints `max Mach 1.745E-02` for the seed whose interior velocity is
  already the stationary one, against the stationary Mach 2.7e-7 of that
  face. The ghost the boundary writes carries -2382 cm/s on the mapped seed
  and -1858 cm/s on that one. Scaling the whole column by 0.971649 so that
  its base face meets the reservoir takes the boundary's own Mach number to
  1.168e-6, the ghost to +0.12 cm/s and the base face mass flux from
  8.1e4 times the wind's to 4.7 times it, and takes the
  integrated base-layer numerators down by 18 (mass) and 15 (energy). The
  offset is the same one L1 saw in the momentum row: the archived column is
  hydrostatic under a gravitational parameter 2.91 percent from this run's
  (MEASURED `|dp/dr|/(rho dPhi/dr)` = 1.02908 median over cells 4 to 300),
  which is the 2.26 percent change of R_0 plus the rest L1 did not attribute.
- **The pass is refused by the element relaxation, not by the
  hydrodynamics.** It is ended by the
  He/H element relaxation, which on this configuration takes 24 to 26 steps

  and then finds no admissible one, `element_step_out_of_bounds`
  (`binary_element_diffusion.f90` line 387, READ), restores its entry
  composition and refuses the pass. This happened for all six seeds,
  including the one whose base-layer residual is 18 times smaller, and it is
  what the HD 209458 b control does NOT do (there the relaxation reaches the
  fixed point of the element operator in 30 steps at every pass). The
  relaxation is given the state the hydrodynamic solve hands back, so this is
  six different states refused and not a proof that a converged one would be;
  what it does show is that no seed of this family reaches pass 2, and that
  the He/H operator on this configuration is a candidate obstruction in its
  own right.
- **No seed of this family enters with an energy row below O(1), and the
  reason is the discretization, not the state.** (The solve does bring it
  down afterwards: S3's handed-back row is 0.626.) At the worst cell of
  each seed the code's discrete `dF_3` is -1.0 times the net heating while
  the continuum `(G/r^2) d(h+Phi)/dr` of the same state is +0.49 to +0.53 of
  it (table T5). The row's own scale is the heating, so the row reads about
  2. The enthalpy-flux divergence the stationary balance is made of is not
  what the discrete operator returns there. The reading offered for that is
  L1's reading of the mass row carried over, an interpretation and not a
  measurement of this item: at Mach 1e-6 to 1e-5 the reconstructed interface
  energy flux is dominated by its acoustic part, which stands 1/Mach above
  the advective one.

**The state that goes furthest is S3, the Bernoulli seed** (section 6): its
handed-back momentum numerator is 26 times below the mapped seed's, its
energy numerator 1.6 times below, its flux-spread gate 7.3e-2 against
1.8e-1, and its cellwise energy row 0.626, the only one of the six that is
neither O(1) nor near 2. It is refused at the same place.

**What this hands back to L5.** Two obstructions are named and neither is the
velocity field: the element relaxation's refusal, which gates the route at
pass 1 on this planet, and the discrete energy flux divergence in a layer at
Mach 1e-6. The seed L5a asked for is built, is 18 times closer in the
base-layer norms than the mapped one, and is kept at `<scratch>/seeds/S1c`;
`--base-anchor` in the seed builder is the part of it L4 would need whatever
is done next.

## 2. The seeds

| seed | what it is |
|---|---|
| mapped | the archived state mapped onto this grid, the entry of L1 and L3 |
| S1 | mapped, with `v = F_wind/(rho r^2)` below cell 233 (10-cell blend) |
| S2b | S1, with the odd-even component of ln rho and ln T removed over cells 1 to 30 (tapered to zero by cell 45, four passes of a half-odd-even subtraction) |
| S1c | S1, with the whole column scaled by 0.971649 so its base face meets the 1 microbar reservoir |
| S2c | S1c with the same odd-even filter |
| S3 | the Bernoulli seed: T over cells 1 to 150 from the stationary energy equation, rho over cells 1 to 232 from hydrostatic equilibrium under the current gravity, base face anchored, `v = F/(rho r^2)` |

`F_wind = 2.375219e12 m_H cm s^-1 R_p^2 = 5.05210e6 g s^-1 sr^-1`
(`Mdot = 6.3487e7 g/s`), the median of `rho v r^2` over the 192 cells at
r >= 2 R_p. The stationary base velocity that follows is 0.02971 cm/s,
which reproduces the 0.0297 cm/s L3 derived from the cell-400 face flux.

The matching radius is **cell 233, r = 1.26567 R_p**: the innermost cell from
which `|rho v r^2/F_wind - 1| <= 0.10` holds in every cell outward. The first
crossing of that band is at cell 15, r = 1.0029, and it is accidental: the
archived base layer oscillates across the wind's value by decades while
approaching it, so the innermost-from-which definition is the one used.

### 2.1 What the restart actually reads, and a correction to L1

`load_IC` builds the primitive state as `W = (rho, v, p)` with **rho from
`calc_rho` of the species densities of `Ion_species_IC.txt`**
(`src/modules/files_IO/load_IC.f90` near line 799, READ), `v` and `p` from
`Hydro_ioniz_IC.txt`, and the temperature not read at all: T follows as p
over the particle count the composition carries.

L1 section 7 observation 2 states that the density is rebuilt from p, T and
the composition and that "a hand-made seed has to change the density through
p and T". That is not what the code does, and the consequence is measurable:
a seed whose ln p and ln T columns were filtered over cells 1 to 30 loads
with a density IDENTICAL to the unfiltered one in every cell (odd-even of
ln rho 4.750e-2 in both) and with the odd-even of ln T slightly WORSE than
the unfiltered state's (4.818e-2 against 4.725e-2), because T is p over a
particle count that kept its own odd-even. L1's fourth seed therefore
changed the pressure alone; its measured effect is the pressure smoothing's.
`make_stationary_seed.py` writes both files: the species of a cell are
scaled by one factor, which carries the density and the particle count
together, and the pressure is then set from that particle count and the
temperature the seed is to hold.

## 3. The residual of each seed (`EXHALE_RESIDUAL=1`)

### T1. Cellwise maxima

| seed | mass max @ cell (r) | momentum max @ cell | energy max @ cell (r) |
|---|---|---|---|
| mapped | 1.9108 @ 4 (1.00077) | 0.5995 @ 1 | 1.9791 @ 137 (1.04877) |
| S1 | 1.5798 @ 43 (1.00831) | 0.6971 @ 1 | 1.9426 @ 71 (1.01453) |
| S2b | 1.8790 @ 10 (1.00193) | 0.4962 @ 1 | 1.9426 @ 71 (1.01453) |
| S1c | 1.7281 @ 2 (1.00039) | 0.0637 @ 500 | 1.9732 @ 62 (1.01223) |
| S2c | 1.8790 @ 10 (1.00193) | 0.1648 @ 1 | 1.9732 @ 62 (1.01223) |
| S3 | 1.8287 @ 150 (1.06151) | 0.2915 @ 499 | 1.5894 @ 68 (1.01372) |

No cellwise maximum falls: the row and its own scale are the same term in
this layer and they fall together, which is what L1 measured for the three
velocity seeds and what still holds here. The momentum row is the one
exception and it is the one whose scale is NOT the row: it is the local
weight, and S1c takes its maximum from 0.60 at cell 1 to 0.064 at the outer
boundary, leaving the 2.83 percent the whole column carries (S1c's momentum
window below 1.03 R_p is 2.680e-2 over 9.465e-1, which is 2.83 percent).

### T2. Window numerator / denominator

| seed | mass r<1.03 | mass 1.03-1.10 | momentum r<1.03 | energy r<1.03 | energy 1.03-1.10 | energy r>=1.20 |
|---|---|---|---|---|---|---|
| mapped | 4.882e-2 / 4.788e-2 | 7.599e-6 / 1.930e-4 | 5.884e-2 / 9.799e-1 | 4.764e-2 / 4.764e-2 | 1.778e-5 / 1.652e-5 | 3.156e-6 / 3.891e-5 |
| S1 | 3.908e-2 / 3.479e-2 | 1.792e-7 / 3.100e-5 | 6.493e-2 / 9.788e-1 | 3.752e-2 / 3.752e-2 | 5.632e-7 / 1.016e-5 | 4.058e-6 / 4.017e-5 |
| S2b | 1.771e-2 / 1.792e-2 | 1.792e-7 / 3.100e-5 | 5.509e-2 / 9.731e-1 | 1.656e-2 / 1.656e-2 | 5.632e-7 / 1.016e-5 | 4.058e-6 / 4.017e-5 |
| S1c | **2.730e-3** / 2.342e-3 | 1.742e-7 / 3.012e-5 | **2.680e-2** / 9.465e-1 | **3.093e-3** / 3.094e-3 | 6.487e-7 / 1.027e-5 | 4.011e-6 / 3.913e-5 |
| S2c | 1.981e-2 / 1.742e-2 | 1.742e-7 / 3.012e-5 | 3.999e-2 / 9.516e-1 | 1.848e-2 / 1.849e-2 | 6.487e-7 / 1.027e-5 | 4.011e-6 / 3.913e-5 |
| S3 | 7.061e-3 / 5.030e-3 | 4.027e-5 / 6.964e-5 | 1.682e-2 / 9.044e-1 | 7.361e-3 / 7.363e-3 | 3.820e-4 / 3.864e-4 | 5.022e-5 / 9.922e-5 |

Three readings.

1. **Above 1.03 R_p the velocity field alone settles the wind**: the mass
   numerator of the 1.03-1.10 window falls by 42 (7.599e-6 to 1.792e-7) and
   the energy numerator of the same window by 32 (1.778e-5 to 5.632e-7), and
   both stay there for every later seed. That window is where L1 measured the
   archived state's 4.5 percent per cell flux decay.
2. **Below 1.03 R_p the base pressure anchor is what moves the norms**, not
   the velocity and not the odd-even filter: 4.882e-2 to 3.908e-2 (velocity)
   to 2.730e-3 (anchor) in the mass row, and 4.764e-2 to 3.752e-2 to
   3.093e-3 in the energy row.
3. **The odd-even filter is not an improvement on an anchored column.** S2c
   is 7 times worse than S1c in both rows. Two reasons, both measured: the
   filter moves rho at cell 1 by up to 4.3 percent, which throws the base
   anchor off (the ghost velocity goes from +0.12 cm/s on S1c to +777 cm/s on
   S2c), and the odd-even component of the archived layer is not noise laid
   over a discrete balance, it is part of the discrete balance the archived
   solver found, so removing it from rho and T independently breaks that
   balance.

### T3. The base face, which is where the O(1) residual of the base cells lives

The base ghost the boundary condition writes carries a constant `rho v r^2`
by construction (`base_ghost_averages`, READ), so the ghost's `rho v r^2` IS
the base face mass flux. MEASURED on a step-free post-processing pass
(`Do only PP: True`, `CFL: 1.0e-12`), which loads the seed, applies the
boundary condition and writes the state it then holds:

| seed | base ghost v [cm/s] | base face `rho v r^2` / F_wind | base face Mach (the run's own line) | v(cell 1) [cm/s] |
|---|---|---|---|---|
| mapped | -2382 | -8.11e+04 | 2.146e-02 | -538.07 |
| S1 | -1858 | -6.90e+04 | 1.745e-02 | +0.029715 |
| S2b | -1079 | -4.03e+04 | 1.016e-02 | +0.029297 |
| S1c | **+0.1245** | **+4.74** | **1.168e-06** | +0.029715 |
| S2c | +777 | +2.98e+04 | 7.308e-03 | +0.029297 |
| S3 | +0.1372 | +3.66 | 1.293e-06 | +0.042815 |

The fourth column is the code's own
`base boundary: inflow stayed subsonic, max Mach ...`, printed by
`check_base_inflow_is_subsonic`; it is the boundary's verdict on itself and
needs no arithmetic outside the code. The stationary value of that face is
Mach 2.7e-7, so the anchored seeds stand a factor 4.4 above it and the
others four decades above it.

The base face columns are the whole point of the item: the interior velocity
of cells 1 to 232 can be set exactly stationary and the base face still
carries 6.9e4 times the wind's mass flux, because the face velocity is the
base pressure mismatch divided by `rho c`. The anchor removes it.

### T4. The cell-centered mass flux of the loaded state

`rho v r^2` in units of the same state's far-wind median, from
`output/residual_profile.txt` (the code's own loaded rho and v):

| cell | r [R_p] | mapped | S1 | S1c | S3 |
|---|---|---|---|---|---|
| 1 | 1.00019 | -1.807e+04 | 0.9979 | 0.9979 | 0.9979 |
| 61 | 1.01199 | 112.9 | 0.9979 | 0.9979 | 0.9979 |
| 100 | 1.02496 | 25.21 | 0.9979 | 0.9979 | 0.9979 |
| 139 | 1.05054 | 4.576 | 0.9981 | 0.9981 | 0.9980 |
| 200 | 1.14888 | 1.282 | 0.9984 | 0.9984 | 0.9986 |
| 232 | 1.26106 | 1.099 | 1.0893 | 1.0893 | 0.9986 |
| 309 | 2.00290 | 1.009 | 1.0092 | 1.0092 | 0.9996 |

The 0.9979 rather than 1.0000 is the seed builder's own residue and is
reported rather than removed: it sets v from the density column of the file
while the run rebuilds the density from the species columns, and the ratio
of the two is not the same at the base as in the wind. The run says why in
its own load line: `He_diffusion: the diffused He/H profile of the restart
file is kept (top cell He/H = 3.2814E-01); base cells set to the input
He/H = 1.6261E+00`, so the base cells are re-composed and the wind is not.
MEASURED, the difference is 0.21 percent. The bump at cell 232
is the blend into the archived wind, which is inside the band the matching
radius was chosen with; S3, whose velocity is continuity's everywhere, does
not have it.

### T5. The energy row: the discrete flux divergence against the continuum one

At the worst energy cell of each seed, `dF_3` as the code assembles it
(from its own worst-cell block) against `(G/r^2) d(h + Phi)/dr` differenced
from the same loaded state:

| seed | worst cell | r [R_p] | continuum dF_3 / (heat-cool) | code dF_3 / (heat-cool) |
|---|---|---|---|---|
| mapped | 137 | 1.04877 | +4.83 | -1.017 |
| S1 | 71 | 1.01453 | +0.527 | -0.967 |
| S1c | 62 | 1.01223 | +0.488 | -1.028 |

Over cells 1 to 150 as a whole, the code's `dF_3` (rebuilt exactly as
`R_3 + (heat - cool)`, which holds because `S(3) = 0` and conduction is off):

| seed | cells of 150 where the code's sign is the continuum's | median \|code dF_3\|/(heat-cool) |
|---|---|---|
| mapped | 28 | 15.6 |
| S1 | 104 | 1.047 |
| S1c | 104 | 1.004 |

The seeds do exactly what L5 asked of them on this measure: the base layer's
energy flux divergence goes from 15.6 times the net heating to 1.0 times it,
and from the wrong sign in 122 of 150 cells to the wrong sign in 46. It is
not enough: the continuum value at the worst cell is +0.5 of the heating
(this is L3's Bernoulli mismatch, table T2 there) and the discrete one is
-1.0, so what the discretization returns for the enthalpy-flux divergence in
this layer is not the quantity the stationary balance is made of. The reading
offered for that is in section 1 and is an interpretation.

## 4. The partitioned route on each seed

`Load IC? True`, `Restart intent: stationary`, `Solver: Newton`,
`EXHALE_PTC_DTAU0=1.0`, `OMP_NUM_THREADS=16`. The input file is the failing
control's, `<scratch>/t_map_stat/input.inp`, unchanged except for the seed,
so the comparison is like for like. That file carries
`Reconstruction scheme: PLM+WENO3` with `du_th [PLM,WENO3]: 0.5 1.0e-3`
while the item's brief named PLM; the key makes no difference here, because
`stationary_state_of_the_loaded_restart` sets `rec_method = 'WENO3'` for
both the measurement and the solve whatever the input says
(`src/EXHALE_main.f90` line 5174, READ), and a stationary restart takes no
marching step (the runs report `steps: 0 accepted of 0 attempted`).

### T6. Outer pass 1, which is the only pass any of them reaches

| seed | hydro | mass | momentum | energy | worst gated species row | s | outcome |
|---|---|---|---|---|---|---|---|
| mapped | info=2 | 1.47E+00 | 1.05E-02 | 1.33E+00 | 1.45E-03 of 1.0E-05 at cell 217 | 321.8 | REFUSED, solve info = 1 |
| S1 | info=2 | 1.11E+00 | 3.21E-02 | 1.99E+00 | 1.66E-03 of 1.0E-05 at cell 217 | 130.3 | REFUSED, solve info = 1 |
| S2b | info=2 | 1.62E+00 | 6.03E-02 | 1.17E+00 | 1.42E-03 of 1.0E-05 at cell 217 | 112.9 | REFUSED, solve info = 1 |
| S1c | info=2 | 1.08E+00 | 2.64E-02 | 1.84E+00 | 1.79E-03 of 1.0E-05 at cell 217 | 220.7 | REFUSED, solve info = 1 |
| S2c | info=2 | 1.41E+00 | 4.39E-02 | 1.80E+00 | 1.68E-03 of 1.0E-05 at cell 217 | 167.8 | REFUSED, solve info = 1 |
| S3 | info=2 | 1.08E+00 | 3.24E-02 | **6.26E-01** | 1.69E-03 of 1.0E-05 at cell 217 | 1324.9 | REFUSED, solve info = 1 |

For comparison, the HD 209458 b element reload on the same route is CERTIFIED
at outer pass 12 with the hydrodynamic rows at 1.2e-9, 1.7e-14, 1.3e-8
(control of PLAN section 1, re-READ from `<scratch>/ctrl_hd209_elem/run.log`).

### T7. The hydrodynamic JFNK of pass 1

| seed | iterations | best judged distance | \|\|R\|\| handed back | flux spread |
|---|---|---|---|---|
| mapped | 77 | 5.972E+08 | 1.466E+00 | 1.807E-01 |
| S1 | 67 | 8.944E+08 | 1.989E+00 | 3.077E-01 |
| S2b | 50 | 6.636E+08 | 1.619E+00 | 2.830E-01 |
| S1c | 95 | 1.072E+09 | 1.842E+00 | 3.636E-01 |
| S2c | 55 | 1.011E+09 | 1.797E+00 | 3.916E-01 |
| S3 | 604 | 5.061E+08 | 1.078E+00 | 7.315E-02 |

Every one ends on STAGNATED, no improvement of the best iterate in 40
iterations, monotone search included. The last iterations of S1 are the
shape all of them have: the damped Gauss-Newton escape fires at
`mu = 1e-5` and moves `||Fs||2` by 0.13 percent per iteration
(4.840e-2 at it 57 to 4.752e-2 at it 68) while `||R||`, the cellwise
maximum, sits between 1.94 and 1.97 and its worst cell walks between 114 and
115. S2c and the mapped seed do the same at `mu = 1e-4`. The solve is not
diverging and it is not descending: it is moving along a valley whose floor
is at a cellwise maximum of about 1.8.

What the solve DOES achieve is in the integrated norms of the state it hands
back, which are three to four decades below the seed's: on S1 the mass
numerator below 1.03 R_p goes from 3.908e-2 to 4.993e-6 and the energy one
from 3.752e-2 to 2.815e-5. The cellwise maxima do not follow, because the
scales fall with them.

### T8. Why the pass is refused, which is not the hydrodynamics

| seed | element relaxation |
|---|---|
| mapped | 26 steps, then no admissible one, outcome 3 at step length 3.612E-02 |
| S1 | 24 steps, then no admissible one, outcome 3 at step length 1.605E-02 |
| S2b | 25 steps, then no admissible one, outcome 3 at step length 2.408E-02 |
| S1c | 24 steps, then no admissible one, outcome 3 at step length 1.605E-02 |
| S2c | 24 steps, then no admissible one, outcome 3 at step length 1.605E-02 |
| S3 | 25 steps, then no admissible one, outcome 3 at step length 2.408E-02 |
| HD 209458 b control | the fixed point of the element operator in 30 steps, at every pass (drift 9.85E-01, 5.67E-01, 2.71E-01, 1.41E-01 at passes 1 to 4) |

Outcome 3 is `element_step_out_of_bounds` (`binary_element_diffusion.f90`
line 387, READ): after 24 to 26 accepted steps the next step puts an element
mass fraction outside its bounds at every admissible step length, the
relaxation restores its entry composition, and `EXHALE_main` then refuses the
pass because "no further pass could differ from this one". That refusal is
what stops the route at pass 1 on LHS 1140 b, it is the same for the seed
whose base-layer residual is 18 times smaller, and it is not something a seed
can move. It is the first item L5 should take next.

## 5. The Bernoulli seed S3, and what its construction measures

S3 is the variation the item's step 5 names: T over cells 1 to 150 from the
stationary energy equation, p by hydrostatic re-integration. Building it
measured two things worth keeping whatever the run does.

- **At the archived mass flux the integration runs away.** Holding
  `F = F_wind` and integrating `dB/dr = (heat - cool) r^2/F` outward from
  cell 1 gives T = 3636 K at cell 150 against the archived 2241 K and a
  column 8.19 times more massive at 1.26 R_p. Integrating inward from cell
  150 instead drives h negative before the base. The archived thermal
  structure and the archived mass flux cannot both carry the current
  heating.
- **The mass flux the base layer's energy budget demands.** Holding both ends
  of B at the archived state's own values fixes F:
  `F = int (heat-cool) r^2 dr / [B(end) - B(1)]`. MEASURED, over cells 1 to
  N with N running outward:

  | cells 1..N | r(N) [R_p] | F_required [g s^-1 sr^-1] | F_required/F_wind |
  |---|---|---|---|
  | 50 | 1.0097 | 9.858e+06 | 1.951 |
  | 100 | 1.0250 | 9.452e+06 | 1.871 |
  | 137 | 1.0488 | 7.702e+06 | 1.524 |
  | 150 | 1.0615 | 7.073e+06 | 1.400 |
  | 200 | 1.1489 | 5.963e+06 | 1.180 |
  | 232 | 1.2611 | 5.734e+06 | 1.135 |
  | 300 | 1.8571 | 5.569e+06 | 1.102 |

  This is L3's table T1 column `G_ener/F_wind` as one number per window. The
  reading: with the current heating, the stationary mass flux the heated
  column demands is about 1.10 times the archived wind's when the integral
  runs out to 1.86 R_p, and 1.9 times it when it stops at 1.025 R_p, so the
  archived temperature profile does not carry the enthalpy gradient the
  current heating deposits in the innermost cells. Whether the true
  stationary Mdot of this planet is 1.1 times the archived one is not settled
  by this arithmetic, because the heating itself moves with the column; it
  is a statement about the archived state, not a prediction.

S3 was built at F = 7.073e6 g/s/sr (the cells 1 to 150 budget), which gives
T = 795 K at cell 50 against the archived 569 K, T = 1399 K at cell 100
against 944 K, and a column 2.90 times more massive at 1.26 R_p. Its base
anchor was calibrated the same way as S1c's, to a base ghost velocity of
+0.137 cm/s. Its residual is better than the mapped seed's and worse than
S1c's (T1, T2), and its worst mass cell is cell 150, the seam where the
energy-equation integration stops: that seam is a defect of the
construction, not of the state.

## 6. The S3 run, which is the one that goes furthest

S3 is refused like the others, and it is nonetheless the best state of the
six by every integrated measure the solve produces. It ran 604 JFNK
iterations in 1325 s, five to ten times the others, because its
`||Fs||2` kept falling by fractions of a percent and the stagnation counter
kept resetting; `||Fs||2` was flat at 1.79e-2 from iteration 80 onward while
`||R||` walked between 1.2 and 2.0 and then came down to 1.078.

### T9. The state each solve hands back, integrated numerator over denominator

| seed | mass | momentum | energy | flux spread gate |
|---|---|---|---|---|
| mapped | 7.448e-06 / 5.571e-04 | 4.132e-07 / 1.004e+00 | 2.345e-05 / 1.432e-04 | 1.807e-01 |
| S1 | 5.704e-06 / 4.651e-04 | 1.795e-06 / 1.004e+00 | 3.810e-05 / 1.325e-04 | 3.077e-01 |
| S1c | 9.121e-06 / 4.736e-04 | 7.157e-07 / 1.004e+00 | 3.885e-05 / 1.285e-04 | 3.636e-01 |
| S3 | **5.943e-06 / 6.699e-04** | **1.603e-08 / 9.987e-01** | **1.446e-05 / 1.560e-04** | **7.315e-02** |

The momentum numerator of S3's handed-back state is 45 times below S1c's and
26 times below the mapped seed's, which is the hydrostatic re-integration
under the current gravity showing up where it should; its energy numerator
is 2.7 times below S1c's and 1.6 times below the mapped seed's; and its flux
spread, 7.3e-2, is the only one of the six below the mapped seed's own
1.8e-1 (the other four are 2.8e-1 to 3.9e-1). Its cellwise energy row, 0.626 at cell 3, is the only one of the six
that is not O(1) or near 2, and it sits at the level L3's Bernoulli mismatch
predicts for a state whose thermal structure now carries the required
enthalpy gradient.

What it does not do is change the ending: 25 element-relaxation steps, then
no admissible one, entry composition restored, pass refused.

## 7. Method, and the commands

`EX` is the repository, `S` the session scratch directory
`.../scratchpad/L`.

```
# the seeds
python3 $EX/LHS1140b/models/make_stationary_seed.py $S/t_atomic/output $S/seeds/S1
python3 $EX/LHS1140b/models/make_stationary_seed.py $S/t_atomic/output $S/seeds/S2b --odd-even
python3 $EX/LHS1140b/models/make_stationary_seed.py $S/t_atomic/output $S/seeds/S1c \
        --base-anchor 0.971649
python3 $EX/LHS1140b/models/make_stationary_seed.py $S/t_atomic/output $S/seeds/S2c \
        --base-anchor 0.971649 --odd-even
python3 $S/mkseed_bernoulli.py $S/t_atomic/output $S/L5a_pp_S1/output $S/seeds/S3 0.971649

# the residual of a seed
mkdir -p $S/L5a_res_<tag>/output
\cp -f $S/t_atomic/input.inp $S/L5a_res_<tag>/input.inp
\cp -f $S/seeds/<tag>/*.txt $S/L5a_res_<tag>/output/
cd $S/L5a_res_<tag>
OMP_NUM_THREADS=8 EXHALE_RESIDUAL=1 EXHALE_MASS_FLOOR_SCAN=1 $S/L5a_bin/EXHALE.x > run.log 2>&1

# the boundary condition's own base face state, with no step taken
sed 's/^Do only PP: False/Do only PP: True/' $S/t_atomic/input.inp > $S/L5a_pp_<tag>/input.inp
printf 'CFL: 1.0e-12\n' >> $S/L5a_pp_<tag>/input.inp
cd $S/L5a_pp_<tag> && OMP_NUM_THREADS=8 $S/L5a_bin/EXHALE.x > run.log 2>&1

# the partitioned stationary route
mkdir -p $S/L5a_<tag>/output
\cp -f $S/t_map_stat/input.inp $S/L5a_<tag>/input.inp
\cp -f $S/seeds/<tag>/*.txt $S/L5a_<tag>/output/
cd $S/L5a_<tag> && . $EX/LHS1140b/winered_hires_y.sh
OMP_NUM_THREADS=16 EXHALE_PTC_DTAU0=1.0 $S/L5a_bin/EXHALE.x > run.log 2>&1
```

`Do only PP: True` takes one time step before it writes (L1 section 7
observation 1); `CFL: 1.0e-12` makes that step a no-op, and the resulting
state reproduces the seed's own v at cell 1 to every printed digit, which is
the check that the pass measures the loaded state.

The base pressure anchor was calibrated against the code's own boundary
condition rather than by rebuilding `continue_hydrostatic_isentrope` outside
it: two post-processing passes at two column scalings give the base ghost
velocity, which is linear in the scaling
(-1857.8 cm/s at 1.000 and +111.5 cm/s at 0.970, a slope of -6.56e4 cm/s per
unit scaling), and one secant step lands the ghost at +0.12 cm/s.

Code units used in the arithmetic outside the code, calibrated as in L1 and
L3: `n_0 = 3.204863e13 cm^-3`, `v_0 = 1.365459e5 cm/s`, `p_0 = 1 microbar`,
volumetric rate unit `p_0 v_0/R_0 = 1.211188e-4 erg cm^-3 s^-1`,
`R_0 = 1.1273716e9 cm`, `M_p = 0.0176220 M_J`. Constants READ from
`src/modules/init/parameters.f90`. The gravity is spherical for this run
(`Domain mode: Spherical`), checked by `|dp/dr|/(rho GM/r^2)` = 1.023 flat
over cells 16 to 300 of the archived state.

The S3 seed builder is `$S/mkseed_bernoulli.py` and is kept in the scratch
directory rather than the repository: its construction has a seam at the cell
where the energy-equation integration stops, and its mass flux is set by the
base layer's energy budget rather than measured from the wind, so it is a
measurement of this item and not a tool.

The run directories are `$S/L5a_res_{map,S1,S2b,S1c,S2c,S3}` (residual),
`$S/L5a_pp_{map,S1,S2b,S1c,S2c,S3}` (the boundary state and the heat and
cool), `$S/L5a_{S1,S2,S2b,S1c,S2c,S3}` (the route), `$S/seeds/*` (the seeds)
and `$S/L5a_bin` (the copy of the binary). The analysis scripts are
`$S/collect.py`, `$S/collect_route.py` and `$S/mkseed_bernoulli.py`.
S1c and S2c were originally made in two steps, `$S/scale_column.py` followed
by the seed builder on `$S/scaled_a`; `--base-anchor` folds that into the one
script and reproduces both files to 4.4e-16 (MEASURED).

## 8. Two observations outside the item, reported and not fixed

1. **`docs/lhs1140b_stationary_L1_20260913.md` section 7 observation 2 is
   wrong about where the restart's density comes from**, and section 4's
   instruction that "the density seed has to be made through p and T"
   follows from it. The density comes from the species columns
   (section 2.1 above). L1's fourth seed changed the pressure only. The
   memo was left as it stands rather than edited, because its own
   measurements are unaffected: what it measured is what that seed did.
2. **The odd-even filter and a base pressure anchor cannot be applied
   independently.** The filter moves rho at cell 1 by up to 4.3 percent,
   which is 150 times the anchor's whole correction, so any filter that
   touches cell 1 has to re-anchor afterwards. `make_stationary_seed.py`
   does not: `--odd-even` on an anchored column undoes the anchor (S2c,
   base ghost +777 cm/s). A filter that holds cell 1 fixed, or an anchor
   applied after the filter, would be the fix; neither was needed for the
   verdict of this item, and `--base-anchor` carries the warning in its help
   text.
