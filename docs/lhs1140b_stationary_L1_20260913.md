# LHS 1140 b, item L1: the stationary residual of the mapped seed, taken apart

Item L1 of `docs/PLAN_20260913_lhs_stationary.md`. Every number here is
MEASURED on the current binary (`EXHALE.x` of 2026-09-13 03:20, git
e8eb6e6ebbbd, tree dirty) unless it is marked READ. No source file was
changed; the whole measurement is the `EXHALE_RESIDUAL=1` diagnostic, a
post-processing pass at a step of zero length, and arithmetic on the files
those two write.

## 1. Verdict

**The O(1) residual of the mapped seed is not the base-layer velocity
sawtooth and not the change in the heating. It is the archived state itself:
that state was never stationary on the present row scale, anywhere between
the base and about 1.1 R_p, and its momentum row is out by a constant 2.8
percent of the local weight over the whole column.** Row by row:

- **Energy row, which sets `||R||` = 1.979.** The sawtooth contributes
  **zero**: the seed whose base velocity is replaced by the constant mass
  flux of cells 31 to 40 gives the same cellwise maximum, 1.979069 at cell
  137, to every digit the code prints, and the same window contributions
  above 1.03 R_p. The heating change contributes **2.5 percent**: putting
  the archived heating back into the same hydrodynamic state lowers the
  worst cell from 1.979 to 1.930 and the base window from 1.798 to 1.287.
  What is left is the state. In every one of cells 1 to 137, that is
  everywhere below r = 1.0488, the energy row is dominated by its own flux
  divergence and not by its radiative terms: |dF_3|/|heat - cool| is 7.1e5
  at cell 1, 3.5e2 at cell 4, 4.9e1 at cell 16 and 4.4 at cell 119. Nothing
  in the row balances that divergence, so the row reads 1.00 of itself
  through the whole layer. At cells 130 to 144 the divergence has fallen to
  the size of the heating and carries the opposite sign, so the two add and
  the row reads twice the local heating rate; that is the 1.979 maximum at
  cell 137. Only from cell 145 outward does the flux divergence take the
  sign of the heating, and there the state comes within 7 to 19 percent of
  balance (dF_3/(heat - cool) = 0.81 at 1.15 R_p, 0.93 at 2.4 R_p) and the
  row falls to 0.07 to 0.30.
- **Mass row, 1.911 at cell 4.** The sawtooth contributes about **30
  percent** of the maximum: the flux-matched seed brings it to 1.338 and the
  3-point average to 1.739. The rest is not in the cell-centered velocity at
  all. The reconstructed face mass fluxes of the base cells stand one to
  five decades above the state's far-wind mass flux and alternate in sign,
  and they still do when every cell-centered rho v r^2 of cells 1 to 30 has
  been set to one constant. The cause is a percent-level odd-even component
  of the density and temperature profile (5.2e-2 in ln rho at cell 1 and
  4.8e-2 at cell 3) in a layer whose Mach number is 2.5e-6: a relative
  density jump d at a face enters the interface flux through its acoustic
  part as d rho c_s, which is d/Mach times rho v. Reducing that component by
  5.2 times lowers the base numerator of the mass row by 2.8 times (4.882e-2
  to 1.747e-2) and of the energy row by 2.9 times (4.764e-2 to 1.622e-2),
  which is the direction and roughly the size the argument predicts, but it
  leaves both cellwise maxima at O(1), because the row and its own scale are
  the same term and fall together. The measurement identifies the face mass
  flux as the quantity that carries the residual; how much of that face flux
  is the density jump and how much is the reconstructed velocity is bounded
  by these two seeds and not resolved by them.
  Above 1.03 R_p the mass row is a flat 4.5 percent of its own scale, which
  is not small when it is followed outward: the state's mass flux is 115
  times its far-wind value at r = 1.012 and reaches that value only near
  1.3 R_p (section 4.1). The archived run's flux gate began at 1.2 R_p and
  reported a spread of 0.077 for this state.
- **Momentum row.** Away from cell 1 it is not a base-layer quantity at all:
  R_2/(rho dPhi/dr) is -2.79 percent (median) at every cell from r = 1.004
  to 10 R_p, and the independently differenced |dp/dr|/|rho dPhi/dr| of the
  same state is 1.0285. A gravity 2.79 percent stronger would null the row.
  The constant changes listed in `LHS1140b/MODELS.md` section 5 are of that
  size and sign: R_0 grew by 2.26 percent when the Jupiter radius became
  7.1492e9 cm (READ; the archived setup report prints 1.73 R_earth for this
  planet and the current one prints 1.77), and b_0 = G M_p m_H/(k_B T_0 R_0)
  is inversely proportional to R_0, so a column hydrostatic under the old
  b_0 is too steep by that factor under the new one. The remaining few
  tenths of a percent are not attributed: I did not recover the constants the
  archived binary carried.

**What this means for L5.** No damping of the base layer and no change of the
base condition will bring this seed below O(1), because two of the three
rows are out for reasons that have nothing to do with the base: the momentum
row by a global rescaling of the gravitational parameter, and the energy row
because the archived state's energy flux divergence has the wrong sign
inside 1.05 R_p. The seed is a starting point for a relaxation, not a state
near a root.

## 2. What the residual is measured against

`residual_row_scale` (`src/modules/time_step/steady_residual.f90`, READ)
divides each row by the largest term that row itself contains:

| row | scale s_k(j) | what R_k/s_k then means |
|---|---|---|
| mass | max(\|F r^2\|_{j+1/2}, \|F r^2\|_{j-1/2})/dV_j | the fractional change of the mass flux across the cell |
| momentum | max(\|ram\|, \|dp/dr\|, \|rho dPhi/dr\|, \|S_visc\|) | the imbalance as a fraction of the largest term of the momentum equation |
| energy | max(\|dF_3\|, \|S_3\|, heat, cool) | the imbalance as a fraction of the largest term of the energy equation |

So "the residual is 1.9" means: the mass flux changes by 191 percent of the
larger of the cell's two face fluxes, which at cell 4 is 18 times the wind's
mass flux, so the flux reverses sign across that cell. "The energy residual
is 1.98" at cell 137 means the row is out by 1.98 times the local
photoionization heating rate, the largest term it holds. "The energy
residual is 1.0000" in cells 1 to 17 means the row is out by 100 percent of
its flux divergence, which there is 50 to 7e5 times the heating: nothing in
the row balances it.

The mass row has no source (`Source.f90` sets S(1) = 0) and the energy row
has none either (S(3) = 0), so R_1 is the face-flux difference alone and
R_3 = dF_3 - (heat - cool) exactly. That identity is what the term tables
below are built from.

## 3. Where the residual sits

### T1. The four seeds and the marched snapshot

`R_k/s_k`, cellwise maximum with its cell, and the window numerator over
window denominator of the integrated norm. All from the code's own
`write_residual_breakdown`.

| state | mass | momentum | energy | mass r<1.03 | momentum r<1.03 | energy r<1.03 | energy 1.03-1.10 |
|---|---|---|---|---|---|---|---|
| mapped seed | 1.9108 @ 4 (1.00077) | 0.5995 @ 1 | 1.9791 @ 137 (1.04877) | 4.882e-2 / 4.788e-2 | 5.884e-2 / 9.799e-1 | 4.764e-2 / 4.764e-2 | 1.778e-5 / 1.652e-5 |
| v 3-point averaged, cells 1-30 | 1.7392 @ 2 | 0.7085 @ 1 | 1.9791 @ 137 | 3.779e-2 / 3.686e-2 | 6.593e-2 / 9.799e-1 | 3.744e-2 / 3.744e-2 | 1.778e-5 / 1.652e-5 |
| v from the constant mass flux of cells 31-40 | 1.3384 @ 4 | 0.6919 @ 1 | 1.9791 @ 137 | 3.937e-2 / 3.766e-2 | 6.459e-2 / 9.786e-1 | 3.756e-2 / 3.756e-2 | 1.778e-5 / 1.652e-5 |
| odd-even p and T removed, cells 1-30, v flux-matched | 1.9768 @ 6 | 0.4859 @ 1 | 1.9791 @ 137 | 1.747e-2 / 2.018e-2 | 5.348e-2 / 9.709e-1 | 1.622e-2 / 1.622e-2 | 1.778e-5 / 1.652e-5 |
| CFL 0.05 march, step near 7227 | 1.7074 @ 165 (1.08028) | 0.7611 @ 185 (1.11432) | 1.9573 @ 327 (2.37300) | 1.499e-3 / 3.115e-2 | 3.269e-3 / 9.460e-1 | 4.469e-3 / 4.464e-3 | 4.247e-3 / 4.248e-3 |

The three velocity seeds differ only in the velocity of cells 1 to 30. The
energy row's cellwise maximum, its cell, and every energy window above
1.03 R_p are identical across all four seeds, which is the cleanest single
statement that the energy row is not a base-layer quantity.

### T2. The mapped seed, cells 1 to 40

| cell | r [R_p] | v [cm/s] | R_1/s_1 | R_2/s_2 | R_3/s_3 | R_1 | R_2 | R_3 | \|F r^2\|/F_far at the outer face |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.00019 | -538.1 | 1.082 | 0.5988 | 1.0000 | 2.329e+02 | 1.575e+02 | 2.264e+02 | 6.28e+03 |
| 2 | 1.00039 | 32.2 | 0.973 | 0.0795 | 1.0000 | -1.725e+01 | -1.938e+01 | -1.652e+01 | 1.67e+02 |
| 3 | 1.00058 | 29.0 | 0.891 | 0.0354 | 1.0007 | -4.197e-01 | -7.782e+00 | -4.483e-01 | 1.82e+01 |
| 4 | 1.00077 | -14.8 | 1.909 | 0.0258 | 1.0029 | -9.780e-02 | -5.401e+00 | -1.073e-01 | 1.65e+01 |
| 5 | 1.00097 | -6.5 | 1.434 | 0.0267 | 0.9982 | 1.540e-01 | -5.216e+00 | 1.729e-01 | 3.81e+01 |
| 6 | 1.00116 | 2.8 | 0.761 | 0.0270 | 1.0033 | -8.172e-02 | -4.737e+00 | -8.931e-02 | 9.09e+00 |
| 7 | 1.00135 | 1.5 | 0.616 | 0.0281 | 0.9953 | 4.115e-02 | -4.377e+00 | 5.847e-02 | 2.37e+01 |
| 8 | 1.00155 | 21.9 | 0.318 | 0.0302 | 0.9950 | 3.111e-02 | -4.220e+00 | 5.103e-02 | 3.47e+01 |
| 9 | 1.00174 | 3.0 | 0.944 | 0.0286 | 1.0020 | -9.233e-02 | -3.717e+00 | -1.302e-01 | 1.95e+00 |
| 10 | 1.00193 | 1.0 | 1.181 | 0.0265 | 1.0050 | -3.581e-02 | -3.287e+00 | -5.165e-02 | 1.08e+01 |
| 11 | 1.00213 | -11.7 | 1.124 | 0.0262 | 0.9947 | 3.406e-02 | -3.059e+00 | 4.829e-02 | 1.33e+00 |
| 12 | 1.00232 | 13.8 | 0.958 | 0.0288 | 0.9982 | 8.571e-02 | -3.119e+00 | 1.356e-01 | 3.18e+01 |
| 13 | 1.00251 | 7.8 | 1.865 | 0.0298 | 1.0009 | -1.668e-01 | -3.036e+00 | -2.615e-01 | 2.75e+01 |
| 14 | 1.00271 | -6.3 | 0.492 | 0.0266 | 0.9958 | 3.808e-02 | -2.628e+00 | 5.825e-02 | 1.40e+01 |
| 15 | 1.00290 | 0.1 | 1.650 | 0.0271 | 0.9984 | 9.968e-02 | -2.576e+00 | 1.563e-01 | 2.15e+01 |
| 16 | 1.00309 | 2.9 | 0.087 | 0.0281 | 0.9795 | 5.748e-03 | -2.553e+00 | 1.180e-02 | 2.35e+01 |
| 17 | 1.00329 | 2.2 | 0.458 | 0.0282 | 1.0053 | -3.025e-02 | -2.445e+00 | -4.692e-02 | 1.28e+01 |
| 18 | 1.00348 | 1.4 | 0.016 | 0.0280 | 0.5477 | -5.824e-04 | -2.313e+00 | 2.970e-04 | 1.26e+01 |
| 19 | 1.00367 | 2.7 | 0.050 | 0.0280 | 0.9471 | 1.850e-03 | -2.216e+00 | 4.368e-03 | 1.32e+01 |
| 20 | 1.00387 | 1.8 | 0.206 | 0.0280 | 1.0211 | -7.635e-03 | -2.122e+00 | -1.177e-02 | 1.05e+01 |
| 21 | 1.00406 | 1.5 | 0.036 | 0.0280 | 1.3724 | -1.050e-03 | -2.029e+00 | -8.959e-04 | 1.01e+01 |
| 22 | 1.00425 | 1.8 | 0.011 | 0.0280 | 0.8571 | 3.251e-04 | -1.944e+00 | 1.454e-03 | 1.02e+01 |
| 23 | 1.00445 | 1.7 | 0.059 | 0.0280 | 0.9436 | 1.804e-03 | -1.869e+00 | 4.047e-03 | 1.09e+01 |
| 24 | 1.00464 | 2.2 | 0.036 | 0.0281 | 1.3011 | -1.107e-03 | -1.796e+00 | -1.041e-03 | 1.05e+01 |
| 25 | 1.00483 | 1.6 | 0.118 | 0.0280 | 1.0469 | -3.461e-03 | -1.721e+00 | -5.380e-03 | 9.25e+00 |
| 26 | 1.00503 | 1.4 | 0.063 | 0.0280 | 1.1233 | -1.620e-03 | -1.652e+00 | -2.191e-03 | 8.68e+00 |
| 27 | 1.00522 | 1.5 | 0.096 | 0.0280 | 0.9576 | 2.577e-03 | -1.591e+00 | 5.419e-03 | 9.60e+00 |
| 28 | 1.00541 | 2.2 | 0.024 | 0.0282 | 0.8959 | 6.669e-04 | -1.537e+00 | 2.057e-03 | 9.84e+00 |
| 29 | 1.00561 | 1.9 | 0.273 | 0.0280 | 1.0184 | -7.504e-03 | -1.473e+00 | -1.317e-02 | 7.15e+00 |
| 30 | 1.00580 | 0.8 | 0.074 | 0.0279 | 0.9360 | 1.587e-03 | -1.412e+00 | 3.485e-03 | 7.72e+00 |
| 31 | 1.00599 | 2.4 | 0.171 | 0.0281 | 0.9746 | 4.450e-03 | -1.368e+00 | 9.093e-03 | 9.31e+00 |
| 32 | 1.00619 | 1.8 | 0.255 | 0.0281 | 1.0202 | -6.617e-03 | -1.316e+00 | -1.194e-02 | 6.94e+00 |
| 33 | 1.00638 | 1.2 | 0.060 | 0.0280 | 0.9243 | 1.240e-03 | -1.263e+00 | 2.881e-03 | 7.39e+00 |
| 34 | 1.00657 | 2.1 | 0.081 | 0.0281 | 0.9464 | 1.827e-03 | -1.221e+00 | 4.148e-03 | 8.04e+00 |
| 35 | 1.00677 | 1.9 | 0.197 | 0.0280 | 1.0299 | -4.405e-03 | -1.176e+00 | -8.067e-03 | 6.46e+00 |
| 36 | 1.00696 | 1.2 | 0.043 | 0.0280 | 0.8960 | 8.010e-04 | -1.134e+00 | 2.016e-03 | 6.75e+00 |
| 37 | 1.00715 | 2.1 | 0.039 | 0.0281 | 0.8974 | 7.692e-04 | -1.098e+00 | 2.039e-03 | 7.02e+00 |
| 38 | 1.00735 | 1.6 | 0.161 | 0.0281 | 1.0411 | -3.153e-03 | -1.058e+00 | -5.894e-03 | 5.89e+00 |
| 39 | 1.00754 | 1.4 | 0.010 | 0.0280 | 0.7511 | 1.581e-04 | -1.019e+00 | 7.008e-04 | 5.95e+00 |
| 40 | 1.00773 | 1.7 | 0.074 | 0.0281 | 0.9314 | 1.321e-03 | -9.867e-01 | 3.139e-03 | 6.42e+00 |

### T3. The five largest cells of each row above r = 1.03, mapped seed

| row | cell | r [R_p] | R/s | R | s | v [cm/s] |
|---|---|---|---|---|---|---|
| mass | 121 | 1.0366 | 0.0446 | -3.3868e-04 | 7.5917e-03 | 54.1 |
| mass | 122 | 1.0372 | 0.0446 | -3.1745e-04 | 7.1185e-03 | 54.4 |
| mass | 120 | 1.0359 | 0.0446 | -3.6093e-04 | 8.0959e-03 | 53.9 |
| mass | 118 | 1.0346 | 0.0446 | -4.1037e-04 | 9.2062e-03 | 53.4 |
| mass | 119 | 1.0353 | 0.0446 | -3.8482e-04 | 8.6333e-03 | 53.6 |
| momentum | 500 | 29.0312 | 0.0551 | -6.0638e-12 | 1.0996e-10 | 168578.2 |
| momentum | 499 | 28.5489 | 0.0485 | -5.5097e-12 | 1.1366e-10 | 165521.4 |
| momentum | 498 | 28.0729 | 0.0443 | -5.3519e-12 | 1.2070e-10 | 162768.3 |
| momentum | 497 | 27.6050 | 0.0442 | -5.6582e-12 | 1.2793e-10 | 160137.2 |
| momentum | 496 | 27.1453 | 0.0426 | -5.7872e-12 | 1.3589e-10 | 157661.6 |
| energy | 137 | 1.0488 | 1.9791 | -2.7682e-04 | 1.3987e-04 | 60.7 |
| energy | 136 | 1.0479 | 1.8836 | -2.9335e-04 | 1.5574e-04 | 60.1 |
| energy | 138 | 1.0496 | 1.8715 | -2.6109e-04 | 1.3951e-04 | 61.3 |
| energy | 135 | 1.0471 | 1.7965 | -3.1111e-04 | 1.7318e-04 | 59.5 |
| energy | 139 | 1.0505 | 1.7698 | -2.4621e-04 | 1.3912e-04 | 61.9 |

### T4. The energy row of the mapped seed, term by term

| cell | r [R_p] | heat | cool | heat - cool | dF_3 | R_3 | dF_3/(heat-cool) | R_3/s_3 | dominant term |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.00019 | 3.193e-04 | 1.463e-06 | 3.1784e-04 | 2.2644e+02 | 2.2644e+02 | 712433.330 | 1.0000 | flux divergence |
| 2 | 1.00039 | 3.065e-04 | 1.499e-06 | 3.0497e-04 | -1.6524e+01 | -1.6525e+01 | -54182.889 | 1.0000 | flux divergence |
| 4 | 1.00077 | 3.097e-04 | 1.564e-06 | 3.0809e-04 | -1.0699e-01 | -1.0730e-01 | -347.262 | 1.0029 | flux divergence |
| 8 | 1.00155 | 2.573e-04 | 1.682e-06 | 2.5561e-04 | 5.1284e-02 | 5.1028e-02 | 200.629 | 0.9950 | flux divergence |
| 16 | 1.00309 | 2.487e-04 | 1.849e-06 | 2.4685e-04 | 1.2048e-02 | 1.1802e-02 | 48.809 | 0.9795 | flux divergence |
| 20 | 1.00387 | 2.455e-04 | 1.920e-06 | 2.4356e-04 | -1.1523e-02 | -1.1766e-02 | -47.310 | 1.0211 | flux divergence |
| 30 | 1.00580 | 2.402e-04 | 2.073e-06 | 2.3813e-04 | 3.7236e-03 | 3.4854e-03 | 15.637 | 0.9360 | flux divergence |
| 40 | 1.00773 | 2.336e-04 | 2.201e-06 | 2.3135e-04 | 3.3703e-03 | 3.1389e-03 | 14.568 | 0.9314 | flux divergence |
| 60 | 1.01176 | 2.168e-04 | 2.400e-06 | 2.1441e-04 | -1.0127e-02 | -1.0341e-02 | -47.232 | 1.0212 | flux divergence |
| 80 | 1.01722 | 1.919e-04 | 2.526e-06 | 1.8942e-04 | -4.6361e-03 | -4.8255e-03 | -24.475 | 1.0409 | flux divergence |
| 100 | 1.02496 | 1.637e-04 | 2.552e-06 | 1.6114e-04 | -1.8138e-03 | -1.9749e-03 | -11.256 | 1.0888 | flux divergence |
| 119 | 1.03528 | 1.458e-04 | 2.568e-06 | 1.4327e-04 | -6.3047e-04 | -7.7375e-04 | -4.400 | 1.2272 | flux divergence |
| 130 | 1.04301 | 1.419e-04 | 2.591e-06 | 1.3935e-04 | -2.7689e-04 | -4.1624e-04 | -1.987 | 1.5033 | flux divergence |
| 137 | 1.04877 | 1.399e-04 | 2.609e-06 | 1.3726e-04 | -1.3955e-04 | -2.7682e-04 | -1.017 | 1.9791 | flux divergence |
| 145 | 1.05626 | 1.361e-04 | 2.630e-06 | 1.3350e-04 | -4.0389e-05 | -1.7388e-04 | -0.303 | 1.2774 | heat - cool |
| 160 | 1.07347 | 1.232e-04 | 2.647e-06 | 1.2059e-04 | 4.4501e-05 | -7.6091e-05 | 0.369 | 0.6174 | heat - cool |
| 180 | 1.10467 | 9.772e-05 | 2.572e-06 | 9.5148e-05 | 6.6063e-05 | -2.9086e-05 | 0.694 | 0.2976 | heat - cool |
| 200 | 1.14888 | 7.037e-05 | 2.347e-06 | 6.8020e-05 | 5.5038e-05 | -1.2982e-05 | 0.809 | 0.1845 | heat - cool |
| 250 | 1.35776 | 2.168e-05 | 1.134e-06 | 2.0551e-05 | 1.8377e-05 | -2.1744e-06 | 0.894 | 0.1003 | heat - cool |
| 300 | 1.85709 | 3.956e-06 | 1.711e-07 | 3.7850e-06 | 3.5131e-06 | -2.7191e-07 | 0.928 | 0.0687 | heat - cool |
| 327 | 2.37300 | 1.256e-06 | 4.038e-08 | 1.2154e-06 | 1.1314e-06 | -8.3959e-08 | 0.931 | 0.0669 | heat - cool |
| 400 | 5.90434 | 2.779e-08 | 2.505e-10 | 2.7542e-08 | 2.5534e-08 | -2.0079e-09 | 0.927 | 0.0722 | heat - cool |
| 500 | 29.03122 | 6.002e-11 | 1.479e-13 | 5.9876e-11 | 4.2431e-11 | -1.7445e-11 | 0.709 | 0.2906 | heat - cool |


Reading T2: the velocity alternates in sign through cell 15 and the mass row
is O(1) there; from cell 16 outward the velocity is smooth and positive and
the mass row falls to a few percent, while the energy row stays at 1.00
because its own flux divergence still has nothing to balance it. The last
column is the reconstructed face mass flux in units of F_far, the same
state's mass flux at the face of cell 400 (r = 5.90 R_p, 5.4636e-7 in code
units), and it is the quantity that makes both rows O(1): 7.6e4 F_far at the
base face, 6.3e3 at the first interior face, and still 6 F_far at cell 40,
where the state's own cell-centered flux is 0.8 F_base = 7 F_far.

## 4. The base layer: the sawtooth and the physics change separated

The archived base layer carries an odd-even mode in the density and the
temperature, not only in the velocity. MEASURED on the seed, with the
odd-even component of ln rho defined as ln rho_j - (ln rho_{j-1} + ln
rho_{j+1})/2:

| cell | r [R_p] | n [cm^-3] | T [K] | dln n/dln r | odd-even of ln rho | \|F r^2\|/F_base, mapped | \|F r^2\|/F_base, flux-matched |
|---|---|---|---|---|---|---|---|
| 1 | 1.00019 | 7.941e+13 | 252.25 | -607.2 | 5.17e-02 | 6.9e+02 | 7.6e+02 |
| 2 | 1.00039 | 7.062e+13 | 269.01 | -593.4 | -2.67e-03 | 1.8e+01 | 7.8e+01 |
| 3 | 1.00058 | 6.313e+13 | 286.63 | -333.8 | -4.75e-02 | 2.0e+00 | 1.2e+01 |
| 4 | 1.00077 | 6.207e+13 | 277.85 | -228.4 | 2.71e-02 | 1.8e+00 | 3.6e+01 |
| 6 | 1.00116 | 5.167e+13 | 303.71 | -607.9 | 5.27e-03 | 9.9e-01 | 1.2e+00 |
| 8 | 1.00155 | 4.029e+13 | 358.72 | -470.2 | -3.53e-02 | 3.8e+00 | 2.1e+01 |
| 12 | 1.00232 | 3.151e+13 | 396.39 | -405.1 | -2.24e-02 | 3.5e+00 | 1.4e+01 |
| 20 | 1.00387 | 2.234e+13 | 429.78 | -219.8 | -3.64e-04 | 1.1e+00 | 5.8e-01 |
| 30 | 1.00580 | 1.501e+13 | 475.55 | -199.0 | 1.52e-03 | 8.4e-01 | 1.1e+00 |

F_base is the state's own mass flux at cells 31 to 40, 5.000e-6 in code
units, which is 9.2 times the flux the same state carries at 5.9 R_p; the
column's flux is not constant and section 4.1 says where it goes. The
logarithmic density gradient swings between -150 and -650 from one cell to
the next, and the odd-even component of ln rho reaches 5.2e-2 at cell 1 and
-4.8e-2 at cell 3. The sound speed there is 1.1e5 cm/s and the flux-matched
velocity is 0.28 cm/s, a Mach number of 2.5e-6, so a density jitter of a few
percent produces a face mass flux two to three decades above rho v. The last
two columns show that this survives the flux match: with every cell-centered
rho v r^2 of cells 1 to 30 set to one constant, the face fluxes are 12 to
760 times it and alternate in sign. That is why the velocity seeds change so
little.

### 4.1 The seed's mass flux is not constant below 1.15 R_p

The mass row above 1.03 R_p is a flat 4.5 percent of its own scale in every
cell (T3), which is small cell by cell and enormous once it is followed
outward. MEASURED, the state's own rho v r^2 in units of its value at the
face of cell 400 (r = 5.90 R_p):

| r [R_p] | 1.004 | 1.008 | 1.012 | 1.017 | 1.025 | 1.036 | 1.051 | 1.073 | 1.149 | 1.36 | 5.90 | 29.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rho v r^2 / its far-wind value | 17 | 7 | 115 | 56 | 25 | 10 | 4 | 2 | 1.3 | 1.05 | 1.00 | 1.07 |

The archived state's wind is a wind only above about 1.3 R_p. Below that it
carries one to two decades more mass flux, the excess decaying at 4.2
percent per cell, which is exactly the 0.045 the mass row reads there. The
velocity steps from 2 cm/s at cell 50 to 44 cm/s at cell 60 (r = 1.008 to
1.012) and then sits at 50 cm/s to cell 100. The flux gate of the archived
run began at 1.2 R_p and reports a spread of 0.077 for this state, so none
of this was visible to it.

T5 shows the three velocity seeds and the seed in which the odd-even
component of ln p and ln T has been removed over cells 1 to 30 (four passes
of a half-odd-even subtraction, tapered to zero by cell 45; the largest
change to p is 1.3 percent and to T 4.1 percent, and the odd-even amplitude
of ln rho over cells 1 to 30 falls from 4.75e-2 to 9.2e-3). Note that
`load_IC` reconstructs the density from p, T and the composition, so the
density seed has to be made through p and T.

### T5. The four seeds, cells 1 to 16 of the base layer

| cell | v mapped | v 3-point | v flux-matched | R_1/s_1 mapped | 3-point | flux-matched | p,T smoothed | R_3/s_3 mapped | 3-point | flux-matched | p,T smoothed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | -538.1 | -168.6 | 0.28 | 1.082 | 1.062 | 1.127 | 1.001 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 2 | 32.2 | -159.0 | 0.31 | 0.973 | 1.739 | 1.102 | 0.899 | 1.0000 | 1.0000 | 1.0000 | 0.9997 |
| 3 | 29.0 | 15.4 | 0.35 | 0.891 | 1.054 | 0.843 | 0.241 | 1.0007 | 1.0000 | 0.9998 | 1.0011 |
| 4 | -14.8 | 2.5 | 0.35 | 1.909 | 0.417 | 1.338 | 0.566 | 1.0029 | 0.9990 | 0.9998 | 1.0006 |
| 5 | -6.5 | -6.2 | 0.38 | 1.434 | 1.224 | 0.786 | 0.733 | 0.9982 | 1.0004 | 1.0004 | 1.0011 |
| 6 | 2.8 | -0.7 | 0.42 | 0.761 | 1.738 | 1.152 | 1.976 | 1.0033 | 0.9993 | 1.0011 | 1.0012 |
| 7 | 1.5 | 8.7 | 0.48 | 0.616 | 1.518 | 0.946 | 0.576 | 0.9953 | 1.0007 | 1.0004 | 1.0014 |
| 8 | 21.9 | 8.8 | 0.54 | 0.318 | 0.051 | 0.046 | 0.014 | 0.9950 | 1.0166 | 1.0842 | 1.0218 |
| 9 | 3.0 | 8.6 | 0.57 | 0.944 | 1.506 | 0.928 | 0.471 | 1.0020 | 0.9989 | 0.9996 | 0.9983 |
| 10 | 1.0 | -2.6 | 0.60 | 1.181 | 0.710 | 1.129 | 0.682 | 1.0050 | 0.9987 | 0.9995 | 0.9978 |
| 11 | -11.7 | 1.0 | 0.63 | 1.124 | 0.579 | 0.948 | 0.040 | 0.9947 | 1.0017 | 1.0006 | 1.0403 |
| 12 | 13.8 | 3.3 | 0.69 | 0.958 | 1.354 | 1.041 | 0.378 | 0.9982 | 1.0005 | 1.0004 | 1.0057 |
| 13 | 7.8 | 5.1 | 0.73 | 1.865 | 1.190 | 0.896 | 0.276 | 1.0009 | 0.9994 | 0.9995 | 0.9914 |
| 14 | -6.3 | 0.5 | 0.74 | 0.492 | 0.467 | 1.280 | 0.954 | 0.9958 | 0.9958 | 0.9991 | 0.9964 |
| 15 | 0.1 | -1.1 | 0.77 | 1.650 | 1.031 | 0.863 | 1.098 | 0.9984 | 1.0019 | 1.0014 | 0.9936 |
| 16 | 2.9 | 1.7 | 0.81 | 0.087 | 1.055 | 0.804 | 0.105 | 0.9795 | 0.9968 | 1.0104 | 1.0902 |

## 5. Does the heating change drive the base energy row

No. Two measurements say so.

**(a) The row is not made of heating there.** T4 gives dF_3 and (heat - cool)
cell by cell. Their ratio is 7.1e5 at cell 1, 3.5e2 at cell 4, 4.9e1 at cell
16, 1.6e1 at cell 30 and 4.4 at cell 119. The heating could be multiplied or
divided by 3.4 in cells 1 to 100 and the energy row would not move by one
percent. The ratio reaches unity at cell 137, r = 1.0488, which is the worst
cell of the whole table, and there the heating has not changed at all: the
current heating is 0.95 times the archived one at that radius.

**(b) Restoring the archived heating leaves the row where it is.** The
hydrodynamic state fixes dF_3, so the row the archived heating would give on
this same state is R_3 = dF_3 - (heat_arch - cool_arch), with the scale
max(\|dF_3\|, heat_arch, cool_arch):

| window | R_3/s_3 max with the current heating | with the archived heating |
|---|---|---|
| r < 1.03 | 1.798 | 1.287 |
| 1.03 to 1.10 | 1.979 | 1.930 |
| 1.10 to 1.20 | 0.316 | 0.370 |
| r >= 1.20 | 0.291 | 0.279 |

The heating change accounts for 2.5 percent of the worst cell and for 28
percent of the base window's maximum; nothing else in the table moves.

**Which channel carries the factor.** MEASURED, current code on the mapped
state against the archived `Heating_breakdown.txt` interpolated to the
current cell centers. The total ratio is 3.62 at cell 1 falling to 3.31 by
cell 5 (the 3.36 quoted in the plan is cell 4), 1.94 at r = 1.0101, 1.00 at
1.0353, 0.86 at 1.0496 and 0.92 to 0.98 in the wind. The factor is carried
by the three photoionization channels together, and the H I and He I logarithms
of the ratio agree with each other to 0.005 at every radius:

| r [R_p] | ln(current/archived): H I | He I | He II | total | He 2^3S | He recombination photons |
|---|---|---|---|---|---|---|
| 1.0002 | 1.308 | 1.311 | 1.206 | 1.287 | -0.097 | -0.120 |
| 1.0010 | 1.220 | 1.223 | 1.104 | 1.197 | -0.115 | -0.161 |
| 1.0050 | 0.916 | 0.921 | 0.740 | 0.884 | -0.189 | -0.563 |
| 1.0101 | 0.713 | 0.719 | 0.495 | 0.664 | -0.240 | -1.055 |
| 1.0200 | 0.454 | 0.463 | -0.020 | 0.376 | -0.350 | -1.580 |
| 1.0353 | 0.142 | 0.071 | -0.463 | -0.004 | -0.468 | -2.301 |
| 1.0496 | -0.092 | -0.097 | -0.545 | -0.157 | -0.533 | -2.858 |
| 1.2007 | -0.006 | 0.028 | -0.319 | -0.084 | -0.509 | -4.465 |
| 5.0 | 0.064 | 0.049 | 0.057 | -0.021 | -0.158 | -8.177 |

He I is 97 percent of the total at the base (3.738e-8 of 3.867e-8 erg
cm^-3 s^-1 at cell 1), as it must be at He/H = 1.6261 in a nearly neutral
layer. Read as an effective optical depth, the base cell now sees a column
smaller by dtau = 1.2 to 1.3 in all three bands at once, and the difference
falls monotonically to zero by r = 1.035 and reverses sign above it. A
common dtau across bands whose thresholds are 13.6, 24.6 and 54.4 eV is the
signature of a change in the overlying column, or in the quadrature of the
attenuated field, and not of the rate coefficients, which would move each
channel by its own factor. The code items in that area are the photon grid
whose bin edges are now the ionization thresholds of H I, He I and He II
(item 2c-QUAD, `src/tests/grid_and_gates/photon_grid_threshold_edges.f90`
states the rectangle-rule argument and the 9.5 percent overcount it removed
at zero optical depth, READ), the cell width `dr_j` and the sub-Lyman band
of `J_inc`. I did not run a controlled build with one of them reverted, so
this report names them and does not apportion the factor among them. The
9.5 percent of the quadrature item is the right size for the far wind, where
the total ratio is 0.98 to 1.02 and the optical depth is small; it is not
the right size for the base, where the ratio is 3.6.

The He 2^3S channel is 0.91 at the base and 0.55 to 0.86 outward, and the
helium recombination-photon channel falls from 0.89 at the base to 0.10 at
r = 1.035 and 3e-4 at 5 R_p, both as the plan states.

## 6. The CFL 0.05 snapshot

The snapshot copied at 10:32 on 2026-09-13 is near step 7227 of
`scratchpad L/t_cfl_0.05` (the run was left untouched and was still running
at step 11121 when this was written). Its base layer has relaxed a long way
and its wind has not:

- the momentum row below 1.03 R_p fell from 0.599 to 0.021 cellwise and from
  5.884e-2 to 3.269e-3 in the window numerator, a factor 18;
- the mass row below 1.03 fell from 1.911 to 1.205 cellwise and from
  4.882e-2 to 1.499e-3 in the numerator, a factor 33;
- the energy row is unchanged in character: 1.00 in every window from the
  base to 1.20 R_p, with the cellwise maximum now at r = 2.373;
- the base face mass flux fell from 4.16e-2 to 1.00e-4 in code units, a
  factor 416, and is now 1.8e2 to 7.9e2 times that state's own far-wind flux
  (5.583e-7) against 7.6e4 for the seed;
- the odd-even amplitude of ln rho over cells 1 to 30 did NOT fall, 5.17e-2
  in the seed against 5.83e-2 in the snapshot: what the march removed is the
  velocity and the face flux the jitter was driving, not the jitter;
- the march has driven a large inward flow through 1.01 to 1.15 R_p
  (v = -5.0e3 cm/s at cell 120 and -7.7e3 at cell 140, against +54 and +63
  in the seed; the cell-centered mass flux there is -474 and -956 times the
  far-wind flux), and that is where the mass and momentum rows now sit
  (1.703 at cell 165, 0.761 at cell 185).

So the low-CFL march is not converging on the seed's wind. Its mass-flux
spread rises monotonically: 0.0773 at step 2, 0.159 at 602, 0.699 at 3002,
1.313 at 6002, 1.543 at 7202, 1.991 at 9602 and 2.272 at 11121. The CFL 0.1
run does the same and faster: 0.0773 at step 2, 1.312 at 3001, 2.435 at
6001, 3.625 at 9001 and 4.818 at 11349. The plan's section 1 entry saying
these marches "settle" and that du is "falling" at step 499 does not survive
the full history: du at step 499 is about 0.13 on its way up from 0.077, not
on its way down.

### T7. The CFL 0.05 snapshot against the mapped seed

| cell | r [R_p] | v seed | v snapshot | R_1/s_1 seed | snapshot | R_2/s_2 seed | snapshot | R_3/s_3 seed | snapshot |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.00019 | -538.1 | -25.2 | 1.082 | 0.585 | 0.5988 | 0.0022 | 1.0000 | 0.9995 |
| 2 | 1.00039 | 32.2 | 46.4 | 0.973 | 0.448 | 0.0795 | 0.0003 | 1.0000 | 0.9998 |
| 4 | 1.00077 | -14.8 | 9.5 | 1.909 | 0.407 | 0.0258 | 0.0008 | 1.0029 | 1.0005 |
| 8 | 1.00155 | 21.9 | 37.6 | 0.318 | 0.075 | 0.0302 | 0.0009 | 0.9950 | 1.0303 |
| 16 | 1.00309 | 2.9 | 16.5 | 0.087 | 0.032 | 0.0281 | 0.0018 | 0.9795 | 1.0619 |
| 30 | 1.00580 | 0.8 | -6.7 | 0.074 | 0.318 | 0.0279 | 0.0028 | 0.9360 | 1.0026 |
| 60 | 1.01176 | 49.9 | -207.4 | 0.034 | 0.029 | 0.0279 | 0.0061 | 1.0212 | 1.0016 |
| 100 | 1.02496 | 51.0 | -1986.3 | 0.042 | 0.008 | 0.0279 | 0.0181 | 1.0888 | 1.0012 |
| 120 | 1.03592 | 53.9 | -4997.0 | 0.045 | 0.011 | 0.0279 | 0.0196 | 1.2428 | 1.0020 |
| 140 | 1.05146 | 62.6 | -7747.0 | 0.038 | 0.047 | 0.0278 | 0.0274 | 1.6752 | 0.9956 |
| 165 | 1.08028 | 92.6 | 50.3 | 0.021 | 1.703 | 0.0276 | 0.0823 | 0.5012 | 0.9972 |
| 185 | 1.11432 | 144.9 | 3583.3 | 0.011 | 0.854 | 0.0275 | 0.7612 | 0.2588 | 0.9996 |
| 200 | 1.14888 | 209.8 | 20978.5 | 0.007 | 0.036 | 0.0274 | 0.0406 | 0.1845 | 1.0110 |
| 250 | 1.35776 | 732.2 | 14357.5 | 0.002 | 0.035 | 0.0269 | 0.0373 | 0.1003 | 1.0365 |
| 300 | 1.85709 | 2494.7 | 9123.5 | 0.000 | 0.031 | 0.0262 | 0.0297 | 0.0687 | 1.2086 |
| 327 | 2.37300 | 4861.3 | 8788.7 | 0.000 | 0.020 | 0.0259 | 0.0275 | 0.0669 | 1.9573 |
| 400 | 5.90434 | 29255.0 | 29891.2 | 0.000 | 0.001 | 0.0264 | 0.0268 | 0.0722 | 0.1294 |
| 500 | 29.03122 | 168578.2 | 168653.1 | 0.003 | 0.002 | 0.0551 | 0.0227 | 0.2906 | 0.1682 |

## 7. Two observations outside the item, reported and not fixed

1. **`Do only PP: True` takes one time step before it writes.** The exit test
   `if (do_only_pp) exit` sits at the end of the marching loop body
   (`src/EXHALE_main.f90` line 3928), so the run reports "steps: 1 accepted
   of 1 attempted" and every file it writes, `Heating_breakdown.txt`
   included, describes the state after that step and not the state that was
   loaded. On a state that is not stationary this is not small: on this seed
   the default-CFL step moved the base velocity from -538 to -1174 cm/s and
   the velocity at r = 1.05 from 61 to 97 cm/s, and it changed the
   certification's own row measures (mass 1.318 at cell 2 against 1.131 at
   cell 1 for the loaded state). `EXHALE_MAXSTEPS=0` does not prevent it.
   The measurement here used `CFL: 1.0e-12`, which makes the step a no-op:
   the resulting state matches the `EXHALE_RESIDUAL` state to 5e-9 in n and
   T, and its heat and cool columns reproduce the code-unit values the
   residual diagnostic prints at cell 137 to seven digits.
2. **`load_IC` discards the density column of `Hydro_ioniz_IC.txt`.** Both
   read statements (`src/modules/files_IO/load_IC.f90` lines 348 and 351)
   read that column into the scratch variable `tmp`; the loaded state is
   built from v, p, T and the composition of `Ion_species_IC.txt`. This is
   consistent and not a defect, but a hand-made seed has to change the
   density through p and T, which is what the fourth seed above does.

## 8. Method and how to reproduce

The diagnostic is `EXHALE_RESIDUAL=1` (`src/EXHALE_main.f90` near line
1093): it equilibrates the loaded composition, applies the boundary
condition, evaluates the ionization sweep and the steady residual in WENO3,
writes `output/residual_profile.txt` (r, n, v, T, R_mass, R_mom, R_energy in
code units) and the breakdown, and stops. It does not write the row scales,
so they were rebuilt outside the code and checked against the code's own
numbers:

- **mass.** R_1 = (A_+ F_+ - A_- F_-)/dV with S(1) = 0, so G_j = A_j F_j
  obeys G_j = G_{j-1} + R_1(j) dV_j. The recursion was anchored at cell 400,
  where the row is 3e-4 and G is the wind's rho v r^2, and run in both
  directions. The cell faces are r_edg(j) = (r_j + r_{j+1})/2 and
  dV_j = (r_edg(j)^3 - r_edg(j-1)^3)/3, both READ from
  `src/modules/init/define_grid.f90`, evaluated on the ghost-inclusive radius
  column of the output file.
- **energy.** S(3) = 0, so dF_3 = R_3 + (heat - cool) exactly, and
  s_3 = max(\|dF_3\|, heat, cool). heat and cool come from the step-free
  post-processing pass, converted to code units by R_0/(p_0 v_0).
- **momentum.** s_2 was estimated as max(\|ram\|, \|dp/dr\|, rho dPhi/dr)
  with rho dPhi/dr = rho b_0/r^2 analytic and dp/dr differenced from the
  output. This is the only estimated scale in the report.

Code units were calibrated from the two states, not assumed:
n_0 = 3.204863e13 cm^-3 (the code prints it), and the conversion of a
volumetric rate, R_0/(p_0 v_0) = 8.256357e3, measured from the heat and cool
the breakdown prints at cell 137 against the same cell of the
post-processing file, which agrees with R_0/(p_0 v_0) = 8.256393e3 computed
from R_0 = 0.157692 x 7.1492e9 cm, v_0 = sqrt(k_B T_0/m_H) = 1.365449e5
cm/s and p_0 = n_0 k_B T_0 = 1 microbar to 4e-6.

The reconstruction reproduces the code's own cellwise maxima:

| state | mass, code / rebuilt | momentum | energy |
|---|---|---|---|
| mapped seed | 1.9108 / 1.9092 @ 4 | 0.59954 / 0.5988 @ 1 | 1.979069 / 1.9791 @ 137 |
| v 3-point | 1.7392 / 1.7392 @ 2 | 0.70848 / 0.7076 @ 1 | 1.979069 / 1.9791 @ 137 |
| flux-matched | 1.3384 / 1.3384 @ 4 | 0.69192 / 0.6940 @ 1 | 1.979069 / 1.9791 @ 137 |
| CFL 0.05 | 1.7074 / 1.7029 @ 165 | 0.76115 / 0.7612 @ 185 | 1.957302 / 1.9573 @ 327 |

### Commands

With `EX` the repository and `S` the session scratch directory
`.../scratchpad/L`:

```
# the residual of a seed
mkdir -p $S/L1_mapped/output
cp $S/t_atomic/input.inp $S/L1_mapped/
cp $S/t_atomic/output/{Hydro_ioniz_IC.txt,Ion_species_IC.txt} $S/L1_mapped/output/
cd $S/L1_mapped && . $EX/LHS1140b/winered_hires_y.sh
EXHALE_RESIDUAL=1 OMP_NUM_THREADS=8 $EX/EXHALE.x > run.log 2>&1

# the heat and cool of that same state, with no step taken
sed 's/^Do only PP: False/Do only PP: True/' $S/t_atomic/input.inp > $S/L1_pp_seed/input.inp
printf 'CFL: 1.0e-12\n' >> $S/L1_pp_seed/input.inp
cd $S/L1_pp_seed && . $EX/LHS1140b/winered_hires_y.sh
OMP_NUM_THREADS=8 $EX/EXHALE.x > run.log 2>&1
```

The seed builders and the analysis are `$S/mkseed.py` (`smooth`,
`fluxmatch`), `$S/mkseed2.py` (the p and T odd-even removal), `$S/rows.py`
(the row scales) and `$S/tables*.py`. The four test directories are
`$S/L1_mapped`, `$S/L1_smooth`, `$S/L1_fluxmatch`, `$S/L1_pTsmooth`, plus
`$S/L1_cfl005` and the two post-processing directories `$S/L1_pp_seed` and
`$S/L1_cfl005_pp`.
