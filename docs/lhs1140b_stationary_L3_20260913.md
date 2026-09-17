# LHS 1140 b, item L3: can the base layer be stationary under the present verdict

Item L3 of `docs/PLAN_20260913_lhs_stationary.md`, the question L1 left open:
the base cells of this planet are the ones whose residual is O(1), and
`Update_EXHALE_stage1` section 155 warns that the residual there is not a
function of the state alone. If those cells cannot carry a balance that the
certification admits, no solver can be asked to find one.

Every number is MEASURED on the binary of 2026-09-13 03:20 (git 43bc28c, tree
dirty, md5 97e10317a710b9ccc63addbedde3586a) unless marked READ. No source
file was changed; the whole measurement is three environment diagnostics the
code already carries (`EXHALE_RESIDUAL`, `EXHALE_MASS_FLOOR_SCAN`,
`EXHALE_UPDATE_MAP`) and arithmetic on the files they write, with every
rebuilt quantity checked against a number the code itself prints (section 7).

## 1. Verdict

**Double precision is not what stops this planet. A state that satisfies all
three hydrodynamic tolerances in the LHS 1140 b base layer exists: at the
thermodynamics of the mapped seed and a mass flux equal to the far wind's,
the continuity row's rounding floor is 8.43e-10 at cell 1 against a tolerance
of 8.43e-9, the energy row's is 4.56e-9 against 1e-6, and the momentum row's
is 1.90e-14 against 1e-8. No cell of the column is unresolved.** The three
findings behind that, in the order they were measured:

- **The stationary state of this layer is a wind at Mach 2.7e-7.** With
  rho v r^2 equal to the seed's far-wind value, the base velocity is
  0.0297 cm/s, 1.8e4 times smaller than the -538 cm/s the seed carries there.
  This is why L1 found the energy row dominated by its own flux divergence at
  7e5 times the net heating: the divergence is the seed's velocity, not a
  property of the layer.
- **The layer's thermodynamics is already close to a stationary one.**
  Continuity and the energy equation are simultaneously satisfiable only if
  `F_wind d(h + Phi)/dr = (heat - cool) r^2`. The seed's profile carries 42 to
  100 percent of the required Bernoulli gradient, and the mass flux the energy
  equation alone demands rises to 1.85 times the far-wind value near 1.02 R_p
  and falls back to 1.05 to 1.08 beyond 1.15 R_p. A factor 1.85, not a factor
  1e4.
- **The continuity row is admissible only through the anchored rule.**
  `cert_tol_mass_at` returns `max(3e-12, 10 x floor)`, and the floor of a
  stationary state here is exactly `eps/Mach`. Against the fixed 3e-12 alone
  the floor of cell 1 stands 281 times too high, and the anchored gate is the
  one that binds out to cell 232 (r = 1.261). The ceiling of that rule,
  `cert_tol_mass_ceiling = 1`, would be reached only below Mach 2.2e-15.

**What this hands to L5.** The obstruction is dynamical, not arithmetic, and
the marching side of it is now a measured number rather than a suspicion: the
base cell's thermal time E/heat is 3.78e7 s, which is 4.78e7 steps at the run's
own step size, and the transit time at the required velocity from the base to
1.15 R_p is 1.44e8 s, or 1.82e8 steps. The marches of L2 ran 1e4 to 3e4 steps.
On WASP-121 b the same two numbers are 2.9e4 and 5.7e4 steps. The update map
of this seed is clean and linear in dt with no floor (section 5), so the
marching operator's fixed point IS the steady residual's root; it is simply
three to four decades of wall time away.

## 2. What the stationary state of this layer requires

Conduction is OFF for this run: `cond_on` defaults to `.false.`
(`src/modules/init/parameters.f90` line 1100, READ) and `input.inp` carries no
`Conduction` key, so `transport_active()` is false and the energy row is
`R_3 = dF_3 - (heat - cool)` with no viscous or conductive source. The code's
`dF_3` is not the plain enthalpy-flux divergence: `RK_rhs` adds
`[A_+ F_1^+ (phi_i(j) - phi_c(j)) - A_- F_1^- (phi_i(j-1) - phi_c(j))]/dV`
(line 272, READ), which is the gravitational work written so that it vanishes
on a discrete equilibrium. In the continuum that term is `rho v dPhi/dr`, so
the stationary energy equation the code imposes is

    (1/r^2) d/dr[ r^2 (E + p) v ] + rho v dPhi/dr = heat - cool .

With `G(r) = rho v r^2` and `h = (E + p)/rho` this is a linear first-order
equation for G,

    d(G h)/dr + G dPhi/dr = (heat - cool) r^2 ,

integrated outward from cell 1 with `G(1) = F_wind`. Continuity separately
requires `G = F_wind` everywhere, so the two hold together only if

    F_wind d(h + Phi)/dr = (heat - cool) r^2 .

`F_wind = 5.0118e6 g s^-1 sr^-1` (Mdot = 6.298e7 g/s), the seed's own mass
flux at the face of cell 400, r = 5.90 R_p.

### T1. The velocity a stationary state requires, LHS 1140 b

`v_cont` is what continuity alone gives, `F_wind/(rho r^2)`; `v_ener` is what
the energy equation alone gives from the integration above; `G_ener/F_wind` is
how far the two disagree.

| cell | r [R_p] | T [K] | v_cont [cm/s] | Mach_cont | v_ener [cm/s] | Mach_ener | G_ener/F_wind |
|---|---|---|---|---|---|---|---|
| 1 | 1.00019 | 252.3 | 2.966e-02 | 2.684e-07 | 2.966e-02 | 2.684e-07 | 1.000 |
| 2 | 1.00039 | 269.0 | 3.334e-02 | 2.922e-07 | 3.417e-02 | 2.994e-07 | 1.025 |
| 3 | 1.00058 | 286.6 | 3.728e-02 | 3.165e-07 | 3.873e-02 | 3.288e-07 | 1.039 |
| 4 | 1.00077 | 277.8 | 3.790e-02 | 3.268e-07 | 4.362e-02 | 3.761e-07 | 1.151 |
| 5 | 1.00097 | 284.3 | 4.068e-02 | 3.468e-07 | 4.894e-02 | 4.172e-07 | 1.203 |
| 8 | 1.00155 | 358.7 | 5.830e-02 | 4.424e-07 | 6.531e-02 | 4.956e-07 | 1.120 |
| 10 | 1.00193 | 366.5 | 6.411e-02 | 4.813e-07 | 7.644e-02 | 5.739e-07 | 1.192 |
| 15 | 1.00290 | 402.0 | 8.341e-02 | 5.979e-07 | 1.079e-01 | 7.732e-07 | 1.293 |
| 20 | 1.00387 | 429.8 | 1.047e-01 | 7.254e-07 | 1.458e-01 | 1.010e-06 | 1.393 |
| 30 | 1.00580 | 475.6 | 1.552e-01 | 1.023e-06 | 2.432e-01 | 1.602e-06 | 1.567 |
| 40 | 1.00773 | 521.3 | 2.217e-01 | 1.395e-06 | 3.739e-01 | 2.353e-06 | 1.686 |
| 60 | 1.01176 | 619.8 | 4.257e-01 | 2.455e-06 | 7.725e-01 | 4.455e-06 | 1.814 |
| 80 | 1.01722 | 755.8 | 8.843e-01 | 4.613e-06 | 1.633e+00 | 8.518e-06 | 1.846 |
| 100 | 1.02496 | 944.4 | 2.016e+00 | 9.382e-06 | 3.599e+00 | 1.675e-05 | 1.785 |
| 119 | 1.03528 | 1232.1 | 4.877e+00 | 1.975e-05 | 7.878e+00 | 3.189e-05 | 1.615 |
| 137 | 1.04877 | 1730.5 | 1.218e+01 | 4.114e-05 | 1.686e+01 | 5.697e-05 | 1.385 |
| 150 | 1.06151 | 2240.7 | 2.309e+01 | 6.785e-05 | 2.897e+01 | 8.513e-05 | 1.255 |
| 165 | 1.08028 | 2921.5 | 4.507e+01 | 1.145e-04 | 5.242e+01 | 1.331e-04 | 1.163 |
| 180 | 1.10467 | 3627.8 | 8.144e+01 | 1.832e-04 | 9.072e+01 | 2.041e-04 | 1.114 |
| 200 | 1.14888 | 4493.9 | 1.630e+02 | 3.239e-04 | 1.763e+02 | 3.502e-04 | 1.081 |
| 250 | 1.35776 | 5469.9 | 6.910e+02 | 1.177e-03 | 7.298e+02 | 1.243e-03 | 1.056 |
| 300 | 1.85709 | 4287.4 | 2.455e+03 | 4.326e-03 | 2.597e+03 | 4.577e-03 | 1.058 |

Two readings. First, the required velocity is nowhere near what the seed
carries: 0.0297 cm/s at cell 1 against -538 cm/s, and 12.2 cm/s at cell 137
against 60.7. Second, the two requirements never separate by more than a
factor 1.85. The seed's temperature and density profile is a factor of order
unity away from one that can carry this wind, and the state's distance from a
root is in its velocity field.

### T2. How far the seed's thermal structure is from the Bernoulli condition

`B = h + Phi`. The required gradient is `(heat - cool) r^2/F_wind`.

| cell | r [R_p] | heat - cool [erg cm^-3 s^-1] | required dB/dr [erg g^-1 cm^-1] | the seed's dB/dr | ratio actual/required |
|---|---|---|---|---|---|
| 1 | 1.00019 | 3.772e-08 | +9.570e+03 | +7.324e+03 | 0.765 |
| 2 | 1.00039 | 3.618e-08 | +9.185e+03 | +7.624e+03 | 0.830 |
| 5 | 1.00097 | 3.633e-08 | +9.233e+03 | +8.224e+03 | 0.891 |
| 10 | 1.00193 | 3.112e-08 | +7.924e+03 | +3.353e+03 | 0.423 |
| 20 | 1.00387 | 2.947e-08 | +7.533e+03 | +3.303e+03 | 0.439 |
| 50 | 1.00967 | 2.698e-08 | +6.976e+03 | +3.428e+03 | 0.491 |
| 100 | 1.02496 | 1.944e-08 | +5.181e+03 | +3.351e+03 | 0.647 |
| 137 | 1.04877 | 1.660e-08 | +4.635e+03 | +4.502e+03 | 0.971 |
| 150 | 1.06151 | 1.570e-08 | +4.491e+03 | +4.483e+03 | 0.998 |
| 200 | 1.14888 | 8.161e-09 | +2.738e+03 | +2.663e+03 | 0.973 |

## 3. What double precision permits in this layer

The definitions are the code's own. The continuity row's floor is
`mass_row_rounding_floor` (`steady_residual.f90` line 521, READ),

    f_1(j) = eps max_faces[ rho (|v| + c_s) A ] / max_faces|A F_1| ,

which for a state with `rho v r^2 = F_wind` is `eps/Mach`, up to the maximum
over the cell's two faces: MEASURED, `f_1/(eps/Mach)` is 1.019 at cell 1 of
LHS 1140 b and 1.013 at cell 1 of WASP-121 b. The
energy and momentum rows have no such function in the code, so the same
definition was carried over to them, with the row's own conserved variable in
place of rho and the row's own scale in the denominator:

    f_3(j) = eps max_faces[ E (|v| + c_s) A ] / dV / max(|dF_3|, heat, cool)
    f_2(j) = eps (p_+ A_+ + p_- A_-) / dV / max(|dp/dr|, rho dPhi/dr)

and for a stationary state `|dF_3| = |heat - cool|`, so the energy scale is the
heating rate. `f_3` is an estimate, not a code measurement: the one entry point
that would measure it directly is closed (section 6). Its amplification is
bounded by the continuity row's, which IS measured: `EXHALE_MASS_FLOOR_SCAN`
reports step/signal at most 1.87 on this seed, 1.85 on the certified
WASP-121 b root and 2.70 on the hot Uranus seed, and item P16 READ the same
quantity at 1.59 to 2.68 over three states.

### T3. LHS 1140 b, a stationary state at the seed's thermodynamics

`tol_mass` is `cert_tol_mass_at = max(3e-12, 10 f_1)`; `floor/tol` is 0.100
wherever the anchored gate binds, by construction, and lower where the fixed
3e-12 binds.

| cell | r [R_p] | Mach | f_1 | tol_mass | f_1/tol | f_3 | f_3/1e-6 | f_2 | f_2/1e-8 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.00019 | 2.684e-07 | 8.427e-10 | 8.427e-09 | 0.100 | 4.374e-09 | 4.37e-03 | 8.44e-15 | 8.4e-07 |
| 2 | 1.00039 | 2.922e-07 | 8.275e-10 | 8.275e-09 | 0.100 | 4.430e-09 | 4.43e-03 | 8.72e-15 | 8.7e-07 |
| 3 | 1.00058 | 3.165e-07 | 7.601e-10 | 7.601e-09 | 0.100 | 4.520e-09 | 4.52e-03 | 9.19e-15 | 9.2e-07 |
| 4 | 1.00077 | 3.268e-07 | 7.018e-10 | 7.018e-09 | 0.100 | 4.221e-09 | 4.22e-03 | 9.22e-15 | 9.2e-07 |
| 5 | 1.00097 | 3.468e-07 | 6.796e-10 | 6.795e-09 | 0.100 | 3.973e-09 | 3.97e-03 | 9.42e-15 | 9.4e-07 |
| 6 | 1.00116 | 3.752e-07 | 6.404e-10 | 6.404e-09 | 0.100 | 4.015e-09 | 4.02e-03 | 9.99e-15 | 1.0e-06 |
| 7 | 1.00135 | 4.074e-07 | 5.919e-10 | 5.919e-09 | 0.100 | 4.221e-09 | 4.22e-03 | 1.08e-14 | 1.1e-06 |
| 8 | 1.00155 | 4.424e-07 | 5.451e-10 | 5.451e-09 | 0.100 | 4.541e-09 | 4.54e-03 | 1.16e-14 | 1.2e-06 |
| 9 | 1.00174 | 4.632e-07 | 5.020e-10 | 5.020e-09 | 0.100 | 4.559e-09 | 4.56e-03 | 1.20e-14 | 1.2e-06 |
| 10 | 1.00193 | 4.813e-07 | 4.795e-10 | 4.795e-09 | 0.100 | 4.389e-09 | 4.39e-03 | 1.21e-14 | 1.2e-06 |
| 11 | 1.00213 | 5.023e-07 | 4.615e-10 | 4.615e-09 | 0.100 | 4.238e-09 | 4.24e-03 | 1.24e-14 | 1.2e-06 |
| 12 | 1.00232 | 5.374e-07 | 4.421e-10 | 4.421e-09 | 0.100 | 4.337e-09 | 4.34e-03 | 1.29e-14 | 1.3e-06 |
| 13 | 1.00251 | 5.617e-07 | 4.133e-10 | 4.133e-09 | 0.100 | 4.373e-09 | 4.37e-03 | 1.32e-14 | 1.3e-06 |
| 14 | 1.00271 | 5.763e-07 | 3.954e-10 | 3.954e-09 | 0.100 | 4.166e-09 | 4.17e-03 | 1.32e-14 | 1.3e-06 |
| 15 | 1.00290 | 5.979e-07 | 3.853e-10 | 3.853e-09 | 0.100 | 3.983e-09 | 3.98e-03 | 1.33e-14 | 1.3e-06 |
| 16 | 1.00309 | 6.227e-07 | 3.714e-10 | 3.714e-09 | 0.100 | 3.884e-09 | 3.88e-03 | 1.35e-14 | 1.3e-06 |
| 17 | 1.00329 | 6.474e-07 | 3.567e-10 | 3.567e-09 | 0.100 | 3.798e-09 | 3.80e-03 | 1.36e-14 | 1.4e-06 |
| 18 | 1.00348 | 6.727e-07 | 3.430e-10 | 3.430e-09 | 0.100 | 3.711e-09 | 3.71e-03 | 1.38e-14 | 1.4e-06 |
| 19 | 1.00367 | 6.993e-07 | 3.301e-10 | 3.301e-09 | 0.100 | 3.633e-09 | 3.63e-03 | 1.40e-14 | 1.4e-06 |
| 20 | 1.00387 | 7.254e-07 | 3.176e-10 | 3.176e-09 | 0.100 | 3.552e-09 | 3.55e-03 | 1.42e-14 | 1.4e-06 |
| 52 | 1.01006 | 1.958e-06 | 1.165e-10 | 1.165e-09 | 0.100 | 1.859e-09 | 1.86e-03 | 1.89e-14 | 1.9e-06 |
| 110 | 1.02996 | 1.377e-05 | 1.678e-11 | 1.678e-10 | 0.100 | 2.575e-10 | 2.57e-04 | 1.34e-14 | 1.3e-06 |
| 138 | 1.04965 | 4.281e-05 | 5.402e-12 | 5.402e-11 | 0.100 | 9.008e-11 | 9.01e-05 | 1.45e-14 | 1.5e-06 |
| 177 | 1.09927 | 1.674e-04 | 1.370e-12 | 1.370e-11 | 0.100 | 3.070e-11 | 3.07e-05 | 1.70e-14 | 1.7e-06 |
| 217 | 1.20069 | 5.086e-04 | 4.496e-13 | 4.496e-12 | 0.100 | 1.353e-11 | 1.35e-05 | 1.58e-14 | 1.6e-06 |

Column extremes: `f_1` is largest at cell 1, 8.427e-10; `f_3` at cell 9,
4.559e-09; `f_2` at cell 51, 1.899e-14. `10 f_1` reaches the ceiling 1 in no
cell, so no cell is unresolved. The anchored gate binds out to cell 232
(r = 1.261) and the fixed 3e-12 takes over beyond it.

### T4. WASP-121 b, the certified stationary root, the same quantities

`backup/regression/wasp_full_newton/output`, which carries `certified=T` in
its own header. Re-measured here: mass 2.7010e-13, momentum 1.5347e-10,
energy 6.5452e-09, against tolerances 3e-12, 1e-8, 1e-6. `(h-c)/s_3` is the
net heating over the energy row's own scale, for the state as it stands.

| cell | r [R_p] | Mach | f_1 | tol_mass | f_1/tol | f_3 | f_3/1e-6 | (h-c)/s_3 |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.00020 | 3.942e-03 | 5.705e-14 | 3.000e-12 | 0.019 | 3.928e-12 | 3.93e-06 | +0.797 |
| 2 | 1.00039 | 3.981e-03 | 5.656e-14 | 3.000e-12 | 0.019 | 3.910e-12 | 3.91e-06 | +0.795 |
| 3 | 1.00059 | 4.021e-03 | 5.601e-14 | 3.000e-12 | 0.019 | 3.899e-12 | 3.90e-06 | +0.794 |
| 4 | 1.00079 | 4.061e-03 | 5.545e-14 | 3.000e-12 | 0.018 | 3.889e-12 | 3.89e-06 | +0.792 |
| 5 | 1.00099 | 4.102e-03 | 5.491e-14 | 3.000e-12 | 0.018 | 3.879e-12 | 3.88e-06 | +0.790 |
| 6 | 1.00118 | 4.142e-03 | 5.437e-14 | 3.000e-12 | 0.018 | 3.869e-12 | 3.87e-06 | +0.788 |
| 7 | 1.00138 | 4.183e-03 | 5.384e-14 | 3.000e-12 | 0.018 | 3.860e-12 | 3.86e-06 | +0.786 |
| 8 | 1.00158 | 4.224e-03 | 5.332e-14 | 3.000e-12 | 0.018 | 3.850e-12 | 3.85e-06 | +0.784 |
| 9 | 1.00178 | 4.265e-03 | 5.280e-14 | 3.000e-12 | 0.018 | 3.841e-12 | 3.84e-06 | +0.782 |
| 10 | 1.00197 | 4.306e-03 | 5.230e-14 | 3.000e-12 | 0.017 | 3.832e-12 | 3.83e-06 | +0.780 |
| 11 | 1.00217 | 4.347e-03 | 5.180e-14 | 3.000e-12 | 0.017 | 3.822e-12 | 3.82e-06 | +0.778 |
| 12 | 1.00237 | 4.389e-03 | 5.131e-14 | 3.000e-12 | 0.017 | 3.813e-12 | 3.81e-06 | +0.776 |
| 13 | 1.00256 | 4.430e-03 | 5.083e-14 | 3.000e-12 | 0.017 | 3.804e-12 | 3.80e-06 | +0.775 |
| 14 | 1.00276 | 4.472e-03 | 5.035e-14 | 3.000e-12 | 0.017 | 3.795e-12 | 3.79e-06 | +0.773 |
| 15 | 1.00296 | 4.514e-03 | 4.989e-14 | 3.000e-12 | 0.017 | 3.785e-12 | 3.79e-06 | +0.771 |
| 16 | 1.00316 | 4.556e-03 | 4.942e-14 | 3.000e-12 | 0.016 | 3.776e-12 | 3.78e-06 | +0.769 |
| 17 | 1.00335 | 4.599e-03 | 4.897e-14 | 3.000e-12 | 0.016 | 3.767e-12 | 3.77e-06 | +0.768 |
| 18 | 1.00355 | 4.641e-03 | 4.852e-14 | 3.000e-12 | 0.016 | 3.758e-12 | 3.76e-06 | +0.766 |
| 19 | 1.00375 | 4.684e-03 | 4.808e-14 | 3.000e-12 | 0.016 | 3.749e-12 | 3.75e-06 | +0.765 |
| 20 | 1.00395 | 4.726e-03 | 4.764e-14 | 3.000e-12 | 0.016 | 3.740e-12 | 3.74e-06 | +0.763 |
| 51 | 1.01006 | 6.147e-03 | 3.665e-14 | 3.000e-12 | 0.012 | 3.443e-12 | 3.44e-06 | +0.728 |
| 128 | 1.02995 | 1.248e-02 | 1.822e-14 | 3.000e-12 | 0.006 | 1.517e-12 | 1.52e-06 | +0.736 |
| 179 | 1.04997 | 2.287e-02 | 1.007e-14 | 3.000e-12 | 0.003 | 7.124e-13 | 7.12e-07 | +0.792 |
| 261 | 1.10038 | 8.417e-02 | 2.913e-15 | 3.000e-12 | 0.001 | 6.912e-14 | 6.91e-08 | +0.660 |
| 352 | 1.19995 | 2.618e-01 | 1.084e-15 | 3.000e-12 | 0.000 | 5.922e-15 | 5.92e-09 | +0.099 |

The certified root's own continuity row, 2.70e-13, sits 4.7 times above its
base-layer floor and 11 times below its tolerance. Its energy row, 6.55e-09,
sits at cell 380 in the wind (r = 1.245) where the estimated floor is
2.7e-15, so that row is at the solve's target and nowhere near the
arithmetic.

### T5. Where the net heating stands against the rounding of the energy flux

At cell 1 of the LHS column, with `rho v r^2 = F_wind`:

| quantity | value |
|---|---|
| `\|A F_1\|` at both faces | 5.0118e+06 g s^-1 sr^-1 |
| `max[rho(\|v\|+c_s)A]/\|A F_1\|` | 3.795e+06 (this is 1/Mach) |
| `\|A F_3\| = h_tot F_wind` | 9.181e+16 erg s^-1 sr^-1 |
| `(heat - cool) dV` | 1.067e+16 erg s^-1 sr^-1 |
| `(heat-cool)dV / (2 eps \|A F_3\|)` | 2.62e+14 |
| `1/f_3` (the same against the signal rounding) | 2.29e+08 |

There is no cell of the LHS column in which the net heating falls below the
rounding of the energy flux, on either estimate. The energy row can be held
to 1e-6 everywhere from the base outward, with a margin of 219 at its worst
cell and at least 81 once the largest amplification ever measured for the
continuity row (2.70) is applied to it.

## 4. Judgment

**(a) Does a state satisfying the present tolerances exist in this layer?
Yes.** All three rows have room, and the ranking of the three margins is

| row | tolerance | worst floor in the column | margin |
|---|---|---|---|
| continuity | `max(3e-12, 10 f_1)` = 8.43e-09 at cell 1 | 8.43e-10 | 10 by construction |
| energy | 1e-6 fixed | 4.56e-09 at cell 9 | 219 |
| momentum | 1e-8 fixed | 1.90e-14 at cell 51 | 5.3e+05 |

The continuity row is the binding one and it is admissible ONLY through the
anchored rule: the fixed 3e-12 alone stands 281 times below the floor of
cell 1, and the anchored gate is what binds every cell out to r = 1.261. This
is the same rule and the same regime the HD 209458 b element reload already
exercises (base Mach 5.1e-6, floor 3.2e-11; item P16, READ); LHS 1140 b is
one further decade down the same line, at Mach 2.7e-7 and floor 8.4e-10, and
the rule covers it without any change.

**(b) Is any row above its arithmetic floor? No,** so the fallback the item
contemplated is not needed. One caveat is worth recording rather than
hidden: the energy row's tolerance is a fixed number and its floor here is a
carried-over estimate, not a code measurement. Its margin, 219 before
amplification and at least 81 after, is the thinnest of the three, and it is
the one number in this report that a direct measurement would improve. What
would measure it is an ulp scan of the energy row of the kind
`mass_row_rounding_floor_scan` already performs for continuity, or the
quadruple-precision control, which cannot be reached from the diagnostic
entry point as the code stands (section 6).

**(c) Why the marching and the solvers do not reach it.** Four candidates,
three of them now measured.

1. **The marching operator is not the problem.** `EXHALE_UPDATE_MAP` on the
   mapped seed gives, as the largest departure over the column of the
   production update map from the steady residual, `|G_dt + R|`:

   | dt/dt_CFL | dt [code] | mass | momentum | energy |
   |---|---|---|---|---|
   | 1 | 9.5675e-05 | 3.223e+01 | 3.997e+01 | 3.211e+01 |
   | 0.1 | 9.5656e-06 | 3.310e+00 | 6.077e+00 | 3.277e+00 |
   | 0.01 | 9.5673e-07 | 3.368e-01 | 6.367e-01 | 3.338e-01 |
   | 0.001 | 9.5675e-08 | 3.375e-02 | 6.397e-02 | 3.345e-02 |

   against `|R|` = 2.329e+02, 1.575e+02, 2.264e+02. The departure is linear in
   dt over three decades with no floor, the chemistry column is 2.4e-03 and
   the transport and filter columns are exactly zero, and `Apply_BC` moves the
   interior by exactly 0. So the section-161 kick is present and is exactly
   what section 161 says it is, 13.8 percent of the mass row at the full step
   and 0.14 percent at CFL/100, and the fixed point of the march is the root
   of the steady residual to that order. Lowering CFL removes the kick and
   does not bring the state closer, which is the L2 measurement.
2. **The march cannot reach the state on any step budget yet run.** The base
   cell's own relaxation times against the run's step size, `dt = 0.790 s`
   (0.401 of a base-cell sound crossing):

   | quantity | LHS 1140 b | WASP-121 b |
   |---|---|---|
   | thermal time E/heat of cell 1 | 3.777e+07 s = 4.78e+07 steps | 1.066e+05 s = 2.92e+04 steps |
   | transit at the stationary velocity, base to 1.15 R_p | 1.439e+08 s = 1.82e+08 steps | 2.077e+05 s = 5.70e+04 steps (to 1.567 R_p) |

   The marches of L2 ran 1e4 to 3e4 steps. They are three to four decades
   short of the layer's own clock, and no reduction of CFL closes that.
3. **The seed is not near the root.** The required base velocity is 1.8e4
   times below the seed's, and L1 MEASURED the seed's base face mass flux at
   7.6e4 times its own far-wind value. A Newton or Gauss-Newton direction
   computed at that state is not a local model of the root.
4. **The acceptance band narrows as the solve converges.** The continuity
   row's scale is the face mass flux, so as the state approaches the root that
   scale falls by four decades toward `F_wind/dV` and the rounding floor rises
   with it, from 1.1e-14 on the seed (MEASURED by the floor scan at cell 1) to
   8.4e-10 at the root. The tolerance rises by the same factor, so the solver
   must finish inside one decade of the residual's own non-smoothness, against
   53 on WASP-121 b (3e-12 over 5.7e-14). This is the state-dependence section
   155 warns about, stated as a number; it is a candidate for L5 and is not
   demonstrated here to be what refuses the solves.

## 5. Why LHS 1140 b is harder than the hot Uranus

### T6. The base cell of three configurations

Cell 1 of each, all at a 1 microbar base. LHS 1140 b and the hot Uranus are
the mapped seeds, WASP-121 b the certified root.

| quantity | LHS 1140 b | hot Uranus (Model A) | WASP-121 b |
|---|---|---|---|
| T_base [K] | 252.3 | 1127.1 | 2358.5 |
| c_s [cm/s] | 1.105e+05 | 2.632e+05 | 5.104e+05 |
| scale height H [cm] | 4.173e+06 | 8.801e+07 | 2.599e+08 |
| first cell width dr [cm] | 2.179e+05 | 6.823e+05 | 3.114e+06 |
| H/dr | 19.2 | 129.0 | 83.5 |
| F_wind [g s^-1 sr^-1] | 5.012e+06 | 1.122e+09 | 3.202e+12 |
| Mdot [g/s] | 6.298e+07 | 1.410e+10 | 4.024e+13 |
| stationary Mach at cell 1 | 2.684e-07 | 1.455e-05 | 3.942e-03 |
| heat(cell 1)/heat_max | 1.000 | 0.375 | 0.106 |
| radius of heat_max [R_p] | 1.0002 (cell 1) | 1.0259 | 1.2466 |
| cool/heat at cell 1 | 4.58e-03 | 3.25 | 0.203 |
| E/heat at cell 1 [s] | 3.777e+07 | 1.755e+07 | 1.066e+05 |
| continuity floor at cell 1, stationary | 8.43e-10 | 1.53e-11 | 5.70e-14 |
| the same, of the state as loaded | 1.11e-14 (seed) | 1.58e-12 (seed) | 5.70e-14 (root) |

The LHS 1140 b base layer is the hardest of the three for two reasons that
compound. First, its heating peaks in cell 1: the 1 microbar level of this
planet is already the top of the absorbing column, so there is no shielded
base and the full XUV deposition lands in the cells whose velocity is
smallest. Second, nothing local absorbs it: the cooling rate there is
0.46 percent of the heating, against 325 percent on the hot Uranus, whose
H3+ and H2 bands put its base in radiative balance and let it be nearly
static without being wrong. On LHS 1140 b every erg deposited at the base has
to leave as enthalpy flux carried by a wind whose mass flux is 224 times
smaller than the hot Uranus's and 6.4e5 times smaller than WASP-121 b's. The
velocity that carries it is Mach 2.7e-7, and that single number produces both
of the difficulties this item was asked about: the continuity row's rounding
floor, which is `eps/Mach` and therefore 1.5e4 times the WASP-121 b value,
and the relaxation time, which is the layer's thermal time and therefore
5e7 steps of the marching scheme.

## 6. Two observations outside the item, reported and not fixed

1. **`EXHALE_RESID_QUAD` cannot fire from the `EXHALE_RESIDUAL` entry point.**
   `assemble_residual` reads the environment only when
   `ieq_sweep_state_kind .ne. ieq_state_marching` (`steady_residual.f90` line
   239, READ), and the diagnostic block of `EXHALE_main.f90` (near line 1093)
   never moves that flag off its initial `ieq_state_marching`
   (`ionization_equilibrium.f90` line 387, READ). MEASURED:
   `EXHALE_RESID_QUAD=1` and `=2` on the certified WASP-121 b root return rows
   identical to the production ones in every printed digit (2.7010e-13,
   1.5347e-10, 6.5452e-09) and print no quadruple-assembly line. The control
   experiment the module header describes is therefore unreachable from the
   one entry point that evaluates a single state and stops, which is also the
   entry point that would measure the energy row's rounding floor asked for in
   section 4(b).
2. **The residual profile's second column is a mass density, not a particle
   density.** `output/residual_profile.txt` is headed `n[cm-3]` and is written
   as `W(1,jj)*n0` (`EXHALE_main.f90`, the block near line 1141, READ), which
   is rho in units of the hydrogen atom mass. MEASURED on this state: the
   column reads 7.94094e+13 at cell 1 while `p/(k_B T)` of the same cell is
   2.79612e+13, a ratio of 2.840, the mean mass per particle in m_H at
   He/H = 1.6261. A reader taking the column as a particle density is off by
   that factor.

## 7. Method, and how the numbers were checked

Everything below the diagnostics is arithmetic on their output files,
performed in cgs, and every rebuilt quantity was required to reproduce a
number the code itself prints.

- **State.** The mapped seed is `scratchpad L/t_atomic/output/{Hydro_ioniz_IC,
  Ion_species_IC}.txt`, the same seed L1 used. Its heat and cool in cgs come
  from `L1_pp_seed/output/Hydro_ioniz.txt`, the step-free post-processing pass
  of L1 section 8 (`Do only PP: True` with `CFL: 1.0e-12`).
- **Rows.** `EXHALE_RESIDUAL=1` writes `output/residual_profile.txt` in code
  units. The face mass flux was rebuilt by the continuity recursion
  `G_j = G_{j-1} + R_1(j) dV_j` anchored at the face of cell 400 on the
  state's own `rho v r^2`, with `r_edg(j) = (r_j + r_{j+1})/2` and
  `dV_j = (r_edg(j)^3 - r_edg(j-1)^3)/3`, and the energy row from
  `dF_3 = R_3 + (heat - cool)`, which is exact because `S(3) = 0`
  (`Source.f90` line 51, READ) and conduction is off.
- **Code units.** Calibrated from the breakdown line of each run rather than
  assumed. LHS: cell 137 prints `rho = 5.489268E-03`, `rho v = 2.438791E-06`,
  giving rho_0 = 5.363429e-11 g cm^-3 (n_0 = 3.20485e13 m_H cm^-3, against the
  3.2049e13 the run prints) and v_0 = 1.365459e5 cm/s (against
  sqrt(k_B T_0/m_H) = 1.365459e5), hence p_0 = rho_0 v_0^2 = 0.9999997
  microbar and R_0/(p_0 v_0) = 8.256359e3 (L1 measured 8.256357e3).
  WASP-121 b: cell 380 prints `rho = 3.471104E-03`, `rho v = 3.637313E-03`.
- **Checks passed.** The rebuilt continuity row reproduces the code's
  `|R_1|/s_1` at cells 1, 2, 3, 4, 5, 10, 16, 20, 30 to five digits
  (1.08241/1.08240, 0.97342/0.97342, 0.89126/0.89126, 1.91073/1.91080,
  1.43406/1.43410, 1.17937/1.17930, 0.08697/0.08697, 0.20610/0.20611,
  0.07374/0.07375); the rebuilt energy row reproduces L1's T4 at cells 1, 4,
  16, 137, 145, 200, 327 (1.0000, 1.0029, 0.9795, 1.9791, 1.2774, 0.1845,
  0.0669); the rebuilt continuity floor reproduces the code's own
  `EXHALE_MASS_FLOOR_SCAN` column to five digits at every one of cells 1 to 30
  on both planets; the rebuilt WASP-121 b rows reproduce 2.7009e-13 against
  the code's 2.7010e-13 and 6.5452e-09 against 6.5452e-09.
- **Constants.** R_J = 7.1492e9 cm, M_J = 1.8982e30 g, G = 6.67430e-8,
  k_B = 1.380649e-16, m_H = 1.67353284e-24 g, all READ from
  `src/modules/init/parameters.f90`; gamma = 5/3 (atomic, the caloric mixture
  is inactive for LHS 1140 b and WASP-121 b).

### Commands

With `EX` the repository and `S` the session scratch directory
`.../scratchpad/L`:

```
# the residual and the continuity row's rounding floor of the mapped seed
mkdir -p $S/L3_lhs_res/output
\cp -f $S/t_atomic/input.inp $S/L3_lhs_res/
\cp -f $S/t_atomic/output/{Hydro_ioniz_IC.txt,Ion_species_IC.txt} $S/L3_lhs_res/output/
cd $S/L3_lhs_res
OMP_NUM_THREADS=8 EXHALE_RESIDUAL=1 EXHALE_MASS_FLOOR_SCAN=1 $EX/EXHALE.x > run.log 2>&1

# the same on the certified WASP-121 b root (its final state loaded as the IC)
mkdir -p $S/L3_wasp_res/output
sed 's/^Load IC? False/Load IC? True/' $EX/backup/regression/wasp_full_newton/input.inp \
    > $S/L3_wasp_res/input.inp
\cp -f $EX/backup/regression/wasp_full_newton/metals.inp $S/L3_wasp_res/
\cp -f $EX/backup/regression/wasp_full_newton/output/Hydro_ioniz.txt \
       $S/L3_wasp_res/output/Hydro_ioniz_IC.txt
\cp -f $EX/backup/regression/wasp_full_newton/output/Ion_species.txt \
       $S/L3_wasp_res/output/Ion_species_IC.txt
cd $S/L3_wasp_res
OMP_NUM_THREADS=8 EXHALE_RESIDUAL=1 EXHALE_MASS_FLOOR_SCAN=1 $EX/EXHALE.x > run.log 2>&1

# the same on the hot Uranus Model A seed
mkdir -p $S/L3_hotu_res/output
\cp -f $EX/backup/regression/carrier_model_a_newton/{input.inp,base.inp} $S/L3_hotu_res/
\cp -f $EX/backup/regression/carrier_model_a_newton/IC/*.txt $S/L3_hotu_res/output/
cd $S/L3_hotu_res
OMP_NUM_THREADS=8 EXHALE_RESIDUAL=1 EXHALE_MASS_FLOOR_SCAN=1 $EX/EXHALE.x > run.log 2>&1

# the update map of the mapped seed at four step sizes
mkdir -p $S/L3_lhs_umap/output
\cp -f $S/t_atomic/input.inp $S/L3_lhs_umap/
\cp -f $S/t_atomic/output/{Hydro_ioniz_IC.txt,Ion_species_IC.txt} $S/L3_lhs_umap/output/
cd $S/L3_lhs_umap
OMP_NUM_THREADS=8 EXHALE_UPDATE_MAP=1,0.1,0.01,0.001 $EX/EXHALE.x > run.log 2>&1
```

The analysis is `$S/l3lib.py` (the file reader and the grid), `$S/l3.py` (the
state, the continuity recursion, the row scales and the floors), `$S/tab.py`
(the two states with their code-unit calibration) and `$S/final.py` (the
energy integration of section 2). The run directories are `$S/L3_lhs_res`,
`$S/L3_wasp_res`, `$S/L3_hotu_res`, `$S/L3_lhs_umap`, plus `$S/L3_wasp_q1`
and `$S/L3_wasp_q2` for the observation of section 6, and `$S/L3_bin` holds
the copy of the binary the measurements were made with.
