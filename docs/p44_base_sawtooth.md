# The hot-Uranus base velocity sawtooth: what it is, and what four candidate prescriptions do to it

Investigation date: 2026-09-02. Every number below was measured in this
investigation unless it is explicitly marked as taken from another document.
No file in `src/`, `build/`, `backup/` or any pre-existing document was
modified; the binary used is a copy of the tree's `EXHALE.x` taken at
2026-09-02 21:48:04 KST (`gfortran` 13.1.0, 1457856 bytes, md5
`204dffe6fa4d91168dbb070cd38cfa8b`), run from a scratch directory
(`.../scratchpad/p44/`, paths in section 8). This is a measurement note, not a
recommendation: the choice of a prescription is left open.

The working tree was being edited by another worker while these runs were made,
so the binary predates the sources now sitting in `src/`. Every source statement
quoted below was re-checked against the tree as of 2026-09-02 22:20 KST and is
unchanged there: `Apply_BC.f90` (the base ghost closures, the valve and
`shapiro_filter`) has not been touched since 2026-09-01, the `Shapiro filter` /
`Base velocity` / `Base ghost temperature` branches of `input_read.f90`, the
Shapiro call and the `valve_eps = 1e-4` assignment in `EXHALE_main.f90`, the
mass row of `assemble_residual`, and the setup-report gap of section 11 all read
the same before and after.

---

## 1. Summary

**The alternation is very likely a discrete odd-even velocity mode of the
collocated scheme in a very low Mach number layer, not physical infall and not
the HD 189733 b entropy checkerboard.** Three measurements point the same way:

* The pressure is smooth. Over cells 1-12 the alternating amplitude of
  `ln(p/rho^gamma)` is 3.81e-3 against 3.28e-5 for `ln p`, a factor 116. But,
  unlike HD 189733 b, `rho`, `p` and `T` are also *monotone* through the same
  window: the alternation is essentially confined to `v` (section 2).
* The base layer is very well resolved. `H(T_eq)/dr = 288` cells in the code's
  own report (127 cells with the run's actual mean molecular weight
  `mu = 2.27 m_H`), against 2.5 for converged HD 189733 b. The
  `docs/hd189_base_checkerboard.md` trend puts that firmly in the regime where
  the entropy mode decays; the alternating density amplitude here is 2.3e-3
  over cells 1-12 and 3.5e-5 over cells 3-12, i.e. it is a two-cell feature
  rather than the extended checkerboard of HD 189733 b.
* **The discretization is not transporting the flux the cell-centred product
  suggests.** Integrating the steady mass residual `R_mass` (dumped by
  `EXHALE_RESIDUAL=1`) down from the smooth wind gives a face mass flux at the
  base of `+0.98 F_wind` in the baseline, while the cell-centred `rho v r^2` of
  cell 1 is `-196 F_wind` (section 3, and panel (c) of the figure). The base
  face is supplying the wind at essentially the right rate, outward, in every
  converged configuration measured.

**Base refinement settles it: the artifact is first order in `dr` and converges
away, while the wind does not move** (section 7, added 2026-09-03). Across a 4x
base refinement of chain-converged cold starts, `v(1)` reads -248.3, -122.2 and
-60.0 cm/s (ratios 2.03 and 2.04 for successive halvings of `dr`), the cell-1
mass-flux ratio follows, and on the finest grid cells 2-12 are all positive and
smooth. Over the same ladder `log10 Mdot` reads 10.28, 10.29, 10.29 and the H2
front 1.0527, 1.0634, 1.0636 `R_p`. So it is neither an undamped null mode
(which would keep its amplitude) nor a physical oscillation: it is a
first-order-in-`dr` boundary discretization error occupying a few cells, and it
does not contaminate the wind. **What does not converge is the base thermal
jump**: the first cell's excess over the pinned `T0` ghost doubles at every
refinement (29.1, 61.0, 132.0 K) and the base temperature gradient steepens with
it (section 7.3). Refining the grid removes the sawtooth without removing the
boundary condition problem underneath it.

**What holds the first cell above `T0` is a boundary flux mismatch, not a heating
rate** (section 8, added 2026-09-03). At cell 1 the total photoheating is
1.34-1.47e-8 erg cm^-3 s^-1 on all three grids -- volumetric and
`dr`-independent -- while the cooling is 2.4-3.4e-6, H3+ infrared to five
figures: **radiation is a net sink there**, 166 to 256 times the heating. The
gravitational work at the same cell falls as `dr` (it is proportional to the
artifact velocity) while the cell-1 steady residual **grows as `1/dr`** in both
the mass and the energy row, which is what a fixed face-flux error looks like in
a shrinking cell; multiplied back by the cell volume it is a mass-flux error of
order `F_wind` and an energy-flux error of order the wind's own enthalpy flux.
Extrapolated to the base face the interior agrees with the pinned ghost on `p` to
0.1% and disagrees on `T` and `rho`, so the isothermal density-pinned ghost fixes
two thermodynamic variables where the interior supplies one. **`T0 = 1140 K`
is Koskinen et al. (2022) Model A's lower-boundary temperature at 1 microbar**,
applied here nine times deeper; moving the boundary deeper still widens the
disagreement (+29.1 -> +35.3 K), the sign that reading predicts. Section 7.3's
"the base thermal structure does not converge" is corrected in section 8.3:
`T(r)` at fixed physical radius does converge outside the innermost `~1e-3 R_p`.

**Base depth turns out to be the strongest lever, and the two defects are
independent** (section 9, added 2026-09-03). The gate planet's base radius is
Koskinen's 1 microbar radius, while the hand-chosen `n0` puts the base pressure
nine times deeper; moving it to 1 microbar (`Base BC: pressure 1.0`) cuts the
base-face mass-flux error 9.7x (-202.8 -> -21.0 `F_wind`), the thermal mismatch
3.3x (+29.1 -> +8.9 K), removes the alternation beyond cell 2, makes the
molecular layer reproducible to five digits **on the default grid**, converges
from cold in 27 minutes where the 9 microbar base did not converge in 1 h 41 m,
and leaves `log10 Mdot` at 10.27 against 10.28. Across a 28x range in base
pressure `v(1)` is constant to 8% (-228, -248, -248 cm/s) while the flux error
tracks `rho_bc` (ratios 9.7 and 2.9 against 9.0 and 3.16), so the artifact is a
fixed velocity error whose flux impact is set by how heavy the base is.
Meanwhile a run with `T0` raised to the 1200 K the interior asked for closes the
thermal gap (+29.1 -> -8.7 K) and leaves the flux error untouched
(-202.8 -> -203.8): **the thermal mismatch and the flux mismatch are two
independent defects**. Section 8.5's "`T0` is applied at the wrong level" is
withdrawn -- Koskinen's lower atmosphere is isothermal at `T_eff`, so 1140 K is
their value at any depth below their boundary (section 9.1).

**The clearly actionable measurement is section 8.4**: on the default base grid
the converged molecular layer is **not reproducible** -- three `info = 0` states
of the same configuration put the H2 front anywhere in 1.053-1.089 `R_p` and `T`
at `r - R_p = 8e-3` anywhere in 583-949 K -- while one base refinement
(`Base grid [dr,cells]: 1.0e-4 100` with `Grid cells: 550`) makes it reproducible
to 0.0008 `R_p` and 1.3% and puts it in agreement with the 4x grid. `Mdot` never
notices (10.28-10.37 throughout).

**From a cold start the mode is absent from the initial condition and from the
whole breathing transient**, and appears between steps ~24000 and ~30000 as `du`
falls through 0.5 to 0.13; it then holds `v(1) = -248 +- 2 cm/s` for 150000
further steps, the same value the JFNK-finished restarts give (section 7.1).

**The pattern separates into two parts that respond to completely different
things.**

* *Cells 2-12* carry a genuine 2`dr` mode that dissipation removes. A Shapiro
  filter removes it while the marching loop runs (alternating density amplitude
  2.28e-3 -> 2.84e-4, sign flips 3 -> 1), and `Low-Mach damping` removes it in
  the converged state (5.06e-4 -> 1.09e-4, flips 5 -> 3).
* *Cell 1* does not respond to any of that. Its velocity is -249, -249, -248,
  -246 and -247 cm/s in the baseline, with the Shapiro filter, with
  `Base velocity: massflux`, with `Low-Mach damping` and with the two combined.
  It is set by the lower boundary closure itself, and it survives every
  dissipative remedy tested.

**The Shapiro filter is a diagnostic here, and the diagnosis is negative.** It
acts on the marching state only, so the JFNK finish that follows re-establishes
the pattern: the filtered and unfiltered converged states differ by 0.06 cm/s in
`v(1)` and by 0.05 dex in `log10 Mdot`, and the filtered one is slightly *worse*
on every base metric (section 5).

**On the premise that motivated this note.** The reading that "the base face is
an outflow, so no inflow composition condition can fire on this planet" is not
what the measurement shows, on two counts. First, in EXHALE's own convention
inflow at the base is `v > 0` (outward, into the domain), and the ghost velocity
is `+0.94 cm/s` in the converged state -- but that number is an artifact: it is
exactly the softplus floor `eps^2/(4|v_1|)` of the smooth valve that
`EXHALE_main.f90` switches on (`valve_eps = 1e-4`, i.e. 30.7 cm/s) before the
JFNK finish, and it is `0.000` in the marching state where the hard valve
`max(v_1, 0)` applies. Second, and more to the point, the *face* mass flux the
scheme transports there is positive and of order the wind flux. What reads as an
outflow is the cell-centred `v(1)`, which is the artifact. A composition
condition gated on `v(1)` or on the ghost velocity would be reading a
discretization mode; one gated on the face mass flux would read the physical
state.

---

## 2. The pattern, as measured

Converged reference state, `.../scratchpad/ladder/runs/heh0p0793c_r1/output/`
(hot Uranus, `He/H = 0.0793`, `Molecular chemistry: True`, no metals, no
Lyman-Werner field, `Grid type: Mixed` with the default `Base grid` of 50
uniform cells of `dr = 2.0e-4 R_p`; the realized `dr = 1.9303e-4 R_p`).
`F_wind` is the median of `rho v r^2` over `r > 1.4 R_p`, where the profile is
flat to 4e-3.

| `j` | `r` [R_p] | `v` [cm/s] | `rho` [mH/cm^3] | `p` [cgs] | `T` [K] | `rho v r^2` | `/F_wind` |
|---|---|---|---|---|---|---|---|
| ghost -1 | 0.999807 | 0.9361 | 1.2204e+14 | 8.9939 | 1140.0 | 1.142e+14 | 0.542 |
| ghost 0 | 1.000000 | 0.9361 | 1.2204e+14 | 8.9939 | 1140.0 | 1.142e+14 | 0.543 |
| 1 | 1.000193 | -250.2 | 1.3398e+14 | 8.9612 | 1100.3 | -3.354e+16 | -159 |
| 2 | 1.000386 | 14.73 | 1.3951e+14 | 8.8835 | 1047.5 | 2.056e+15 | 9.77 |
| 3 | 1.000579 | -4.312 | 1.3933e+14 | 8.8090 | 1040.1 | -6.015e+14 | -2.86 |
| 4 | 1.000772 | 2.029 | 1.3829e+14 | 8.7336 | 1038.9 | 2.810e+14 | 1.33 |
| 5 | 1.000965 | -0.3724 | 1.3680e+14 | 8.6595 | 1041.3 | -5.104e+13 | -0.242 |
| 6 | 1.001158 | 0.5745 | 1.3509e+14 | 8.5860 | 1045.5 | 7.779e+13 | 0.369 |
| 7 | 1.001351 | 0.3200 | 1.3326e+14 | 8.5136 | 1051.0 | 4.276e+13 | 0.203 |
| 8 | 1.001544 | 0.5137 | 1.3147e+14 | 8.4422 | 1056.3 | 6.774e+13 | 0.322 |
| 9 | 1.001737 | 0.4816 | 1.2987e+14 | 8.3718 | 1060.4 | 6.276e+13 | 0.298 |
| 10 | 1.001930 | 0.4411 | 1.2854e+14 | 8.3021 | 1062.4 | 5.691e+13 | 0.270 |
| 11 | 1.002123 | 0.3254 | 1.2746e+14 | 8.2332 | 1062.5 | 4.165e+13 | 0.198 |
| 12 | 1.002316 | 0.2082 | 1.2655e+14 | 8.1648 | 1061.2 | 2.647e+13 | 0.126 |

`F_wind = 2.106e+14` (in the units of the table). Taking
`Mdot = 0.5 * 4 pi R_p^2 m_H F_wind` -- the halving is the run's
`Rate/2 + Mdot/2` setting -- gives `log10 Mdot = 10.43` against the 10.42 the
run reports, so this is the flux the run's own mass-loss rate is built on.

> The task brief quotes `9.7e13` for the wind `rho v r^2` of this run. That
> value could not be reproduced: the flux is flat at 2.106e14 from `r = 1.33`
> outward, and it is 2.106e14 that reproduces the reported `log10 Mdot`. The
> brief's cell-1 to cell-5 numbers agree with the file to 1-3%, so the two
> readings are of the same state.

**Alternating amplitudes.** Defined here as the projection of `ln x` onto
`(-1)^j` after a cubic detrend, over a 10- or 12-cell window (this is close to,
but not identical with, the second-difference projection used in
`docs/hd189_base_checkerboard.md` section 11.3; the two are compared only
qualitatively).

| field | cells 1-12 | cells 3-12 | HD 189733 b converged, cells 1-12 (that document) |
|---|---|---|---|
| `ln(p/rho^gamma)` | 3.81e-3 | -- | 0.0875 |
| `ln rho` | 2.27e-3 | 3.48e-5 | 0.0537 |
| `ln T` | 2.29e-3 | -- | 0.0443 |
| `ln p` | 3.28e-5 | -- | 0.00197 |

**Reading.** The pressure is smooth, as in HD 189733 b, but everything else is
23 times smaller and 65 times more localized (the 1-12 to 3-12 window drop is
65x here against ~6x there). `rho`, `p` and `T` are monotone through cells 2-12;
only `v` alternates, and it does so strictly over cells 1-5 (successive amplitude ratios 17.0, 3.4, 2.1, 5.5) and
then stops -- cells 6-12 are all positive and smoothly
decaying. So this is **a boundary-localized, geometrically damped odd-even
oscillation of the velocity**, not the sustained entropy checkerboard of
HD 189733 b, in which `rho` and `T` alternated in anti-phase at 5% amplitude
out to cell ~20.

**Grid margin.** `EXHALE_setup.out` reports `H(T_eq)/dr = 288.3` cells
(`1/(b0 dr_1)`, i.e. with `mu = m_H`). With the run's own mean molecular weight,
`mu = rho_1/(p_1/k T_1) = 2.27 m_H`, the physical value is 127 cells. Either
number is far above the 8-17 transition that
`docs/hd189_base_checkerboard.md` section 4.4 identifies, and far above the
10-cell threshold at which `write_setup_report.f90` warns. Underresolution of
the base scale height -- the controlling parameter on HD 189733 b -- is not the
mechanism here.

---

## 3. What the scheme actually transports through the base

`EXHALE_RESIDUAL=1` evaluates the finite-volume steady residual of the loaded
state and writes `output/residual_profile.txt`, whose `R_mass` column is
`dF(1,j)` in code units, i.e. the discrete divergence of the numerical mass
flux. Integrating it downward from a cell in the smooth wind (`j = 400`,
`r = 2.1 R_p`, where the cell-centred product and the face flux agree to 1e-4)
recovers the face flux the scheme is using:

`Phi_(j-1/2) = Phi_(j+1/2) - R_mass(j) r_j^2 dr_j` (code units).

Baseline (`A_base`), with `F_wind` the same flux in code units
(`5.18e-6`):

| `j` | `v` [cm/s] | cell-centred `rho v r^2 / F_wind` | face flux `Phi / F_wind` |
|---|---|---|---|
| 1 | -248.9 | **-196** | **+0.98** |
| 2 | +21.9 | +18.1 | +0.70 |
| 3 | -4.91 | -4.24 | +0.95 |
| 4 | +3.06 | +2.73 | +0.95 |
| 5 | -0.103 | -0.095 | +0.95 |
| 6 | +1.06 | +1.00 | +0.95 |
| 8 | +0.773 | +0.76 | +0.95 |
| 16 | +0.586 | +0.63 | +0.95 |
| 201 | +1000 | +1.01 | +1.00 |
| 401 | +3.47e5 | +1.00 | +1.00 |

The same reconstruction at the base face for every state measured:

| run | cell-centred `rho v r^2 / F_wind` at cell 1 | reconstructed face flux `/ F_wind` |
|---|---|---|
| `heh0p0793c_r1` (reference) | -159 | +0.36 |
| `A_base` | -196 | +0.98 |
| `B_shapiro` | -181 | +0.34 |
| `C_massflux` | -202 | +0.94 |
| `E_lowmach` | -198 | +0.44 |
| `D_ghostT` (unconverged) | -26.9 | +14.6 |

The reconstruction is approximate -- it assumes the exact divergence form and it
accumulates the residual of the whole quasi-hydrostatic layer, which is why the
converged runs spread over 0.34-0.98 rather than sitting on 1.00. What is not
approximate is the sign and the order of magnitude: the base faces of every
converged run carry an **outward** mass flux **of order the wind flux**, while
the cell-centred product at cell 1 is 200 times the wind flux **inward**. The
cell-centred `rho v` at the base is therefore not a flux the discretization
transports, which is the signature of a collocated odd-even velocity mode rather
than of physical infall.

*Related observation, not exercised by this configuration.*
`diffusive_photochemistry.f90` states in its header that its advection is "a
one-sided upwind difference on the cell velocity, not on face-averaged
velocities". At the base of this planet the cell velocity is the artifact above,
so the upwind direction that operator would choose at cell 1 is set by the mode.
That operator needs `Oxygen chemistry: True`, which is off here, so nothing in
these runs depends on it; it is recorded because it is the same quantity a base
composition condition would naturally be gated on.

---

## 4. The four configurations, measured

All runs restart from the same converged reference state
(`ladder/runs/heh0p0793c_r1/output/`) with `Load IC? True`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-3` and `Solver: Newton 100.0`, i.e. the chain
settings that produced it, so the comparison is like for like. The converged runs each march 2002
steps (the secondary-ionization flip at step 2 plus the `N_stall = 2000` hold)
and then hand off to the JFNK finish. `OMP_NUM_THREADS=8`,
`EXHALE_MAXSTEPS=20000`.

`A(ln rho)` is over cells 1-12; "flips" counts sign changes of `v` over cells
1-12; `r_front` is where `2 n(H2) = n(H I)`; the spread is
`(max-min)/median` of `rho v r^2` over `r > 1.4 R_p`.

| run | key added | JFNK | `\|\|R\|\|` | `log10 Mdot` | `v(1)` [cm/s] | `v_ghost` | cell 1 `/F_wind` | `A(ln rho)` | flips | `r_front` | max `n(H3+)` | spread |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| reference `heh0p0793c_r1` | -- | info=0 | 9.89e-4 | 10.42 | -250.2 | +0.94 | -159 | 2.27e-3 | 5 | 1.1012 | 3.09e5 | 4.1e-3 |
| **A** `A_base` (control) | -- | info=0 | 8.67e-4 | **10.29** | -248.9 | +0.94 | -196 | 5.06e-4 | 5 | 1.0614 | 3.54e5 | 2.9e-3 |
| A' `A_base2` (2nd cycle) | -- | info=0 | 1.70e-4 | 10.29 | -247.1 | +0.95 | -197 | 4.91e-4 | 5 | 1.0607 | 3.56e5 | 2.4e-3 |
| **B** `B_shapiro` | `Shapiro filter: 0.02 4` | info=0 | 9.61e-4 | 10.34 | -248.7 | +0.94 | -181 | 7.21e-4 | 7 | 1.0839 | 3.23e5 | 2.4e-3 |
| **C** `C_massflux` | `Base velocity: massflux` | info=0 | **2.66e-4** | 10.28 | -247.5 | +1.71 | -202 | 1.90e-4 | 5 | 1.0536 | 3.71e5 | 2.4e-3 |
| **D** `D_ghostT` | `Base ghost temperature: continuous` | **info=2** | 2.28e-3 | (10.35) | -39.2 | +5.28 | -26.9 | 7.45e-2 | 3 | 1.0003 | 4.51e5 | 2.4e-3 |
| E `E_lowmach` | `Low-Mach damping: 2.0e-2` | info=0 | 7.62e-4 | 10.29 | -246.2 | +0.95 | -198 | **1.09e-4** | 3 | 1.0699 | 3.39e5 | 2.4e-3 |
| E' `E_lowmach2` (2nd cycle) | as E | info=0 | 4.02e-4 | 10.30 | -245.6 | +0.95 | -191 | 3.73e-4 | 3 | 1.0675 | 3.43e5 | 2.4e-3 |
| F `F_hydrobase` | `Hydrostatic base: True` | **info=2** | 2.41e-3 | (10.47) | +64.9 | +68.4 | +32.8 | 2.45e-3 | 0 | 1.0001 | 9.74e5 | 2.5e-3 |
| G `G_mf_lm` | C and E together | info=0 | 9.87e-4 | 10.29 | -246.9 | +1.71 | -199 | 1.04e-3 | 7 | 1.0692 | 3.40e5 | 2.4e-3 |

E, F and G are additions, not part of the four configurations asked for: E and G
because `Low-Mach damping` is the key the code documents for exactly this
symptom (`README_HOWTO.md`: "a shell has stopped flowing ... the local Mach
number falls to 1e-5 and a 2 dr mode grows there that the contact-resolving HLLC
flux no longer damps"), F because it is the other closure that writes the base
ghost pressure and is mutually exclusive with D.

**The control matters.** The first restart cycle moves `log10 Mdot` from 10.42
to 10.29 and the H2 front from 1.1012 to 1.0614 with no key added at all; the
second cycle then moves them by 0.00 dex and 0.0007 `R_p`. So the reference
state was not yet a fixed point of the restart-plus-JFNK workflow, and **A, not
the reference, is the control**. Cycle-to-cycle drift on either baseline is
about 0.01 dex in `Mdot` and 0.002-0.003 `R_p` in the front, which is the noise
floor against which the runs must be read.

### (a) Amplitude of `rho v r^2` in cells 1-8, against the wind value

| run | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| A baseline | -196 | +18.1 | -4.24 | +2.73 | -0.10 | +1.00 | +0.60 | +0.76 |
| B Shapiro | -181 | +15.3 | -3.84 | +0.46 | -0.63 | +0.39 | -0.02 | +0.39 |
| C massflux | -202 | +18.0 | -3.46 | +2.63 | -0.03 | +0.97 | +0.62 | +0.79 |
| D ghost `T` | -26.9 | +46.2 | -16.7 | +13.0 | +2.51 | +8.62 | +5.29 | +6.33 |
| E Low-Mach | -198 | +13.9 | -1.16 | +1.10 | +0.18 | +0.18 | +0.39 | +0.27 |

### (b) `log10 Mdot`

Quoted only for the runs that reached `info = 0`: A 10.29, B 10.34, C 10.28,
E 10.29, G 10.29, against the control's own second-cycle 10.29. B's +0.05 dex is
five times the cycle noise; C, E and G move `Mdot` by at most 0.01 dex. D and F
never converged, so their 10.35 and 10.47 are not quotable.

### (c) Does the base face become an inflow?

In no run does `v(1)` change sign except D (-39 cm/s, still negative) and F
(+65 cm/s, in a state that has already destroyed the layer -- see below). The
ghost velocity is positive in every run, but in A, B and E it is only the
softplus floor of the smooth valve (`0.5(v_1 + sqrt(v_1^2 + eps^2))` with
`eps = 1e-4` code units `= 30.7 cm/s`, which for `v_1 = -248.88` evaluates to
`0.940` against the `0.9412` in the file). C is the one run where the ghost
velocity is a real number: `F_c/(rho_bc r^2) = +1.71 cm/s`. And as section 3
shows, the *face* flux is already outward at about the wind rate in all of them.

### (d) H2 front and H3+

The front sits at 1.053-1.084 `R_p` across A, B, C, E and G, against 1.0614 for
the control and a 0.002 `R_p` cycle noise; B is the outlier at 1.0839 (+0.023),
in the same direction as its +0.05 dex in `Mdot`. Peak `n(H3+)` moves with the
front, 3.2-3.7e5 cm^-3, and `n(H3+)` in cell 1 is 2.91e5 cm^-3 in every one of
them to 0.5%. D collapses the front onto the base (1.0003) and F destroys it
altogether (1.0001, peak `n(H3+)` 9.7e5 at the ghost).

### (e) Residual and wind-region spread

`||R||` is the volume-weighted JFNK norm at exit. C is the best (2.66e-4, 3.3
times below the control), E next (7.62e-4), B the worst of the converged runs
(9.61e-4). The wind-region spread of `rho v r^2` is 2.4e-3 in every converged
run, i.e. the wind above `r = 1.4 R_p` is insensitive to all of this.

### D and F in detail

**D (`Base ghost temperature: continuous`)** does not work on this planet. The
JFNK finish fails three times (`info = 2`, then 1, then 2) and the run exits on
`du`, and the state it leaves is a different solution, not a cleaner version of
the same one: the base ghost heats to 1775 K against `T0 = 1140 K`, the
alternating density amplitude rises by 147x to 7.45e-2, the cell-2 flux reaches
+46 `F_wind`, and the molecular front collapses to 1.0003 `R_p`. This is the
opposite of what the key does on HD 189733 b
(`docs/hd189_base_checkerboard.md` section 13.2: amplitude down 2.9x, `||R||`
down 1.9x, `Mdot` +0.01 dex). The difference is plausibly that HD 189733 b had a
5x ghost-to-cell-1 temperature jump for the key to remove, whereas here the jump
is 3.5% (1140 vs 1100 K) and there is nothing for it to fix -- what it does
instead is release the base temperature from its `T0` pin in a layer whose
molecular chemistry depends steeply on it. This reading is inferred from the two
measurements, not demonstrated.

**F (`Hydrostatic base: True`)** fails harder and in the way section 4.2 of that
same document recorded for HD 189733 b (its E1: alternation gone, base draining).
Here the base runs away to 10200 K in the ghost and 5800 K in cell 1, the
pressure rises 9x to 80 cgs, the whole first 12 cells settle into a uniform
+66 cm/s outflow at 33 `F_wind`, and the molecular layer is gone. The
alternation is indeed removed (0 sign flips) but of a state that is no longer
the problem's solution.

---

## 5. The Shapiro filter, isolated

`shapiro_filter` is called from the marching loop only
(`EXHALE_main.f90`, every `shapiro_every` steps); the JFNK residual does not see
it, which `write_setup_report.f90` states in a comment. To separate the two,
both the baseline and the filtered run were also stopped at
`EXHALE_MAXSTEPS = 1900`, i.e. before the hand-off, so the marching state itself
could be read.

| | `v` cells 1-8 [cm/s] | `A(ln rho)` 1-12 | flips | `r_front` | `log10 Mdot` (du-stop, not quotable) |
|---|---|---|---|---|---|
| A2 marching only, no filter | -247, +16.1, -3.00, +3.38, +0.96, +1.92, +1.66, +1.85 | 2.28e-3 | 3 | 1.1012 | 10.41 |
| B2 marching only, Shapiro | **-219**, +34.8, +9.53, +16.4, +13.9, +15.1, +14.8, +15.0 | **2.84e-4** | **1** | 1.1013 | 10.41 |

So the filter does work, on the part of the pattern that is a 2`dr` mode: cells
2-5 stop alternating and cells 6-12 settle on a nearly uniform +14.5 to
+15.1 cm/s, and the alternating density amplitude falls 8x. It does not touch the front (1.1012 vs
1.1013) or the wind. **But it does not remove cell 1** (-219 against -247, 11%),
and after the unfiltered JFNK finish the whole pattern is back: B against A
differs by 0.06 cm/s in `v(1)`, and B is worse than A on `A(ln rho)` (7.21e-4
against 5.06e-4) and on sign flips (7 against 5).

Note also that, filtered, cells 2-12 carry a *uniform* cell-centred
`rho v r^2` of about 9 `F_wind`. Smoothing the mode does not make the
cell-centred product equal the transported flux; it only makes the error
smooth. That is consistent with section 3 and with the reading that the two are
different quantities in this layer.

---

## 6. Candidate judgment (tentative)

Against the three readings the brief poses:

**(i) A grid null mode that only a filter removes, with no physical content --
this appears to be right for cells 2-12, and wrong for cell 1.** Cells 2-12
respond to every dissipative operator tried (Shapiro during marching, Low-Mach
damping in the flux), the pressure is smooth across them, and the face flux the
scheme transports through them is flat at about `F_wind`. Cell 1 responds to
none of them.

**(ii) A reflection created by the boundary condition -- this appears to be the
right reading for cell 1, but no boundary variant tested removes it without
destroying the solution.** `v(1)` is -249, -249, -248, -246, -247 cm/s across
the baseline and four different interior/boundary treatments; the two keys that
actually rewrite the base ghost (D, F) both move it, and both fail to converge
and both collapse the molecular layer. So the boundary is where cell 1 comes
from, and none of the existing closures is a usable prescription for it here.

**(iii) A physical base oscillation that survives every prescription -- this
does not appear to be supported.** A physical infall of `-250 cm/s` at
`rho = 1.34e14` would be a mass flux 200 times the wind, and the steady mass
residual says the faces are not carrying it; the layer is also 127 scale-height
cells thick, so there is nothing for a physical 2-cell oscillation to be.

> **See also section 9.6** for the final judgment table, and **section 8.7**,
> which settles the intermediate question -- whether the
> first cell is raised above `T0` by a real heating rate or by a boundary flux
> mismatch -- in favour of the flux mismatch, and corrects section 7.3.

> **Update after section 7 (2026-09-03).** The refinement ladder sharpens (i)
> and (ii) into a single statement that neither quite made: the alternation is a
> **boundary discretization error that is first order in `dr`**. It is not an
> undamped null mode -- a null mode keeps its amplitude under refinement, and
> this one halves each time `dr` halves. It is not a fixed-width physical
> boundary layer either -- its cell count falls (6 -> 4) rather than doubling.
> And it is not physical: `Mdot` and the H2 front are unchanged across the same
> 4x refinement. Reading (iii) is now excluded on grid evidence as well as on
> the flux evidence of section 3. What survives, and what none of the four
> prescriptions addresses, is the base *thermal* mismatch of section 7.3, which
> refinement makes worse rather than better.

**Which prescription cleans the base without moving `Mdot` or the front?** On
these measurements, `Low-Mach damping: 2.0e-2` is the best of the four for the
interior part of the pattern: it is the only key that reduces the converged
alternating amplitude (4.6x on the first cycle, though only 1.3x on the
second -- E' 3.73e-4 against A' 4.91e-4, so the amplitude gain is not stable
across cycles while the sign-flip count is), it cuts the sign flips 5 -> 3, it leaves `Mdot`
at the control value to 0.00-0.01 dex and the front within 0.009 `R_p`, and it
enters the numerical flux, so the JFNK residual and the marching loop see the
same equation. `Base velocity: massflux` is second: 2.7x on the amplitude, the
best `||R||` of all (2.66e-4), `Mdot` -0.01 dex, front -0.008 `R_p`; it is also
the only key that gives the ghost a velocity with a meaning. Combining the two
(G) is not additive and was worse than either alone on `A(ln rho)`, which is not
explained here. The Shapiro filter changes nothing that survives the JFNK
finish, and `Base ghost temperature: continuous` and `Hydrostatic base` both
destroy this solution.

**None of them turns the base face into an inflow in the cell-centred sense, and
on the evidence of section 3 that may be the wrong target.** The face already
carries an outward flux of order `F_wind`. If a base composition condition needs
to know the direction of the flow across the base face, the face mass flux -- or
the ghost velocity under `Base velocity: massflux`, which is derived from it --
looks like the quantity to gate on, and the cell-centred `v(1)` does not.

---

## 7. Cold start and grid refinement (2026-09-03)

Section 10 listed two things the restart-based measurements could not settle: what
a run started from scratch does, and how the pattern responds to the base grid
(which needs a cold start, because `load_IC` does not interpolate onto a
different grid). Both were run. Predictions were written down before the refined
grids left their cold-start transient and are reproduced verbatim in section 9.

### 7.1 Cold start: the mode is not in the initial condition and not in the transient

Cold runs use the reference `input.inp` with `Load IC? False`,
`du_th [PLM,WENO3]: 0.5 1.0e-3`, `Solver: Newton 2.0e-2`, `OMP_NUM_THREADS=8`
and a 200000-step budget, with a snapshot archived at every thousandth step.
Two quantities track the pattern: the odd-even purity of the velocity,
`oe = |sum_j (-1)^j v_j| / sum_j |v_j|` over cells 1-6 (1 = a pure 2`dr`
alternation, 0 = a smooth profile), and `v(1)` itself.

**(a) It is not there at the start.** The initial condition
(`EXHALE_DUMP_IC=1`) has all velocities positive and monotone
(+29.6, +59.2, +88.8, +118.4, +148.0, +177.6 cm/s over cells 1-6), zero sign
changes and an alternating density amplitude of 1.5e-13 -- numerically zero.

**(b) It appears late, when the base breathing dies.** Steps 1-1000 are the
known cold-start base breathing: `v(1)` swings over +8.8e3, +9.0e3, -1.3e3,
-6.0e3, -1.5e4 cm/s at steps 5, 10, 100, 200 and 500, with `oe` below 0.15 and
no strict alternation anywhere. That continues to step ~23000. Then:

| step | `v(1)` [cm/s] | `oe` | `du` |
|---|---|---|---|
| 0 (IC) | +29.6 | 0.00 | -- |
| 12091 | +96.1 | 0.10 | 4.75 |
| 17141 | -144.0 | 0.89 | 0.71 |
| 21180 | -477.2 | 0.32 | 0.087 |
| 24118 | -263.1 | 0.87 | 0.52 |
| 27102 | -251.5 | 0.96 | 0.29 |
| 30224 | -244.9 | **1.00** | 0.13 |
| 36146 | -242.3 | 1.00 | 0.026 |
| 37857 | *secondary ionization activated* | | 0.0200 |
| 39081 | +1015.9 | 0.03 | 0.014 |
| 60107 | -248.7 | 0.98 | -- |
| 198041 | -248.3 | 0.99 | 0.0143 |

So the mode **appears between step ~24000 and ~30000**, as `du` falls through
0.5 to 0.13 and the base settles into its quasi-hydrostatic state; the
secondary-ionization flip at step 37857 knocks it out for ~10000 steps and it
re-forms at the same value. From step ~50000 to the end of the budget it is
stationary: `v(1) = -248 +- 2 cm/s` over 150000 steps.

**(c) The cold-start value is the converged value.** `v(1) = -248.3` at 198638
marching steps, against `-248.9` in the JFNK-finished restart of section 4 and
`-250.2` in the reference chain state: the same number to 1%. The amplitude is
therefore not something the JFNK finish creates, nor something the marching
loop leaves behind -- both routes land on it.

**Low-Mach damping makes no difference on a cold start.** `CS_E`
(`Low-Mach damping: 2.0e-2`, otherwise identical) tracks `CS_A` step for step:
at the 200000-step budget the two give `v(1) = -213.5` against `-216.0`,
`F_wind` 2.3977e14 against 2.3987e14, and an H2 front at 1.11116 against
1.11122 `R_p` -- differences of 1% or less. It is only once each is
JFNK-finished that the key separates them (section 7.3), which is consistent
with its stated purpose: the term is gated on `M < 1e-3`, and a marching state
whose base is still adjusting is not yet there everywhere.

Neither 200000-step cold run reached `info = 0` (three JFNK attempts each,
`info = 1`, best `||R|| = 2.07e-3`), so their `log10 Mdot = 10.47` is not
quotable. Chained restarts of both did converge, and are used below.

### 7.2 Base refinement: prediction, then measurement

`Base grid [dr,cells]` sets the uniform base region at fixed 0.01 `R_p` extent;
`Grid cells` was raised in step so the stretched region keeps its 450 cells and
the outer grid is unchanged:

| label | key | realized `dr` [R_p] | total cells |
|---|---|---|---|
| 1x | (default) | 1.930e-4 | 500 |
| 2x | `Base grid [dr,cells]: 1.0e-4 100` + `Grid cells: 550` | 9.601e-5 | 550 |
| 4x | `Base grid [dr,cells]: 5.0e-5 200` + `Grid cells: 650` | 4.776e-5 | 650 |

**Predictions, recorded before the refined runs left their transient**
(full text in section 9): P1 -- if it is a 2`dr` grid null mode it occupies a
fixed *number of cells*, so halving `dr` halves its *physical* width and does
not by itself reduce `v(1)`. P2 -- if it is a boundary reflection with a
physical length scale it occupies a fixed *physical width*, so halving `dr`
*doubles* the number of alternating cells and reduces the per-cell amplitude.
P3 -- on cell 1: a fixed error in the mass *flux* keeps `rho v r^2` near
-200 `F_wind` at every `dr`; a fixed momentum deposited in the first cell makes
`v(1)` scale with `dr`. P4 -- `Mdot` and the front should already be near
grid-converged.

**Measurement.** All three states below are chain-converged (`info = 0`), each
the second restart cycle of its own cold start, so they sit at the same point of
the workflow. `nalt` is the number of cells over which `v` alternates strictly
from cell 1; `D2(ln p)` is the `(-1)^j` projection of the discrete second
difference of `ln p`, which is grid-independent for a true null mode and falls
with `dr` for a converging discretization.

| grid | `dr` [R_p] | `v(1)` [cm/s] | ratio | cell 1 `/F_wind` | `nalt` | width [R_p] | `D2(ln p)` | `\|\|R\|\|` | `log10 Mdot` | `r_front` |
|---|---|---|---|---|---|---|---|---|---|---|
| 1x | 1.930e-4 | **-248.3** | -- | -202.8 | 6 | 1.16e-3 | 7.14e-5 | 2.48e-4 | **10.28** | 1.0527 |
| 2x | 9.601e-5 | **-122.2** | 2.03 | -94.5 | 4 | 3.84e-4 | 1.55e-5 | 3.57e-4 | **10.29** | 1.0634 |
| 4x | 4.776e-5 | **-60.0** | 2.04 | -43.8 | 4 | 1.91e-4 | 3.23e-6 | 5.66e-4 | **10.29** | 1.0636 |

**`v(1)` is first order in `dr`, to 2%.** Halving `dr` halves it twice over
(ratios 2.03 and 2.04); the cell-1 mass-flux ratio follows (2.15 and 2.16), and
the 2`dr` content of the pressure falls faster still (4.6x then 4.8x, i.e. as
`dr^2.2`). On the 4x grid the alternation is essentially gone from the interior:
cells 2-12 of `CS_A4x` read +2.36, +0.12, +1.38, +0.97, +1.13, +1.06, +1.08,
+1.06, +1.05, +1.04, +1.04 cm/s, all positive and smooth, with only cell 1
left negative.

So P2 is refuted (the cell count falls from 6 to 4, it does not double, and the
physical width falls by 6x rather than staying fixed) and P1 is only half right:
the disturbance does stay a few cells wide, but its amplitude does **not** stay
put -- it converges away at first order. P3's second branch is the one that
matches.

**Mdot and the H2 front are grid-converged (P4 confirmed).** `log10 Mdot` reads
10.28, 10.29, 10.29 across a 4x refinement -- inside the 0.01 dex restart-cycle
noise of section 4. The front moves 0.011 `R_p` from 1x to 2x and then 0.0002
`R_p` from 2x to 4x, so it is converged from the 2x grid on; peak `n(H3+)` is
3.73e5, 3.51e5, 3.50e5 cm^-3 on the same ladder. The base artifact does not
contaminate the wind.

### 7.3 What does *not* converge: the base thermal jump

The same three states, at the lower boundary:

| grid | `dr` [R_p] | `T_ghost` [K] | `T(1)` [K] | `T(1) - T_ghost` | `rho(1)/rho_ghost` | `T(2)` | `T(3)` |
|---|---|---|---|---|---|---|---|
| 1x | 1.930e-4 | 1140.0 | 1169.1 | **29.1** | 1.0333 | 1113 | 1060 |
| 2x | 9.601e-5 | 1140.0 | 1201.0 | **61.0** | 1.0076 | 1159 | 1123 |
| 4x | 4.776e-5 | 1140.0 | 1272.0 | **132.0** | 0.9521 | 1218 | 1173 |

The excess of the first cell over the pinned ghost **doubles every time `dr`
halves** (ratios 2.10 and 2.16), and the difference quotient `(T_2 - T_1)/dr`
steepens with it: -2.9e5, -4.4e5 and -1.13e6 K per `R_p`. The density ratio
across the boundary moves the same way, from 1.033 through 1.008 to 0.952.

> **Corrected 2026-09-03 (section 8.2).** The first reading written here -- that
> the base thermal structure "does not converge" -- was too strong, and the way
> it was measured is what made it look that way. `T(1)` is sampled at
> `r - R_p = dr`, so a finer grid reads a *different point* of the profile, and
> the profile rises steeply toward the face; comparing the three grids at the
> same physical radius instead shows `T(r)` converging for `r - R_p > 1e-3 R_p`
> (2x and 4x agree to 0.5% there). What survives, and is measured in section 8,
> is narrower and better posed: the innermost `~1e-3 R_p` does not converge,
> the first cell carries a steady residual that grows as `1/dr`, and the
> converged interior extrapolates to a base-face temperature well above the
> pinned `T0`. Read section 8 in place of the paragraph below.

Note what this means for the base condition: the velocity artifact of sections
2-3 vanishes as `dr -> 0`, but the thermal mismatch it sits next to does not.
Refining the grid removes the sawtooth without removing the boundary condition
problem underneath it.

### 7.4 Low-Mach damping in the JFNK-converged cold-start state

`CS_E_r1` is the chained, `info = 0` finish of the Low-Mach cold start, to be
read against `CS_A_r2`/`A_base` on the same grid: `v(1) = -245.2` against
-248.3, `nalt = 4` against 6, `D2(ln rho) = 2.01e-4` against 6.6e-4, `||R|| =
4.46e-4` against 2.48e-4, `log10 Mdot = 10.29` and `r_front = 1.0644` --
i.e. the same wind as the 2x and 4x grids give, with a base narrower by two
cells. That reproduces the section 4 result from an independent (cold) starting
point: the key shortens the interior tail of the pattern and leaves the wind
alone, and it does not touch cell 1.

---

## 8. The base thermal mismatch: term-by-term budget, a deeper base, and where 1140 K comes from (2026-09-03)

Section 7.3 reported that the first cell's temperature excess over the pinned
ghost doubles at every refinement and read that as a base structure that does not
converge. This section measures the energy budget that produces it, tests a
deeper lower boundary, and traces the boundary temperature to its source. The
first result is a correction to section 7.3 itself.

### 8.1 No photo channel heats the first cell -- radiation there is a net sink

`Heating_breakdown.txt` and `Cooling_breakdown.txt` at cell 1 of the three
chain-converged grid states (each the second restart cycle of its own cold
start; `Viscosity` and `Conduction` are off by default and were not given, so
there is no conduction term):

| quantity at cell 1 [erg cm^-3 s^-1] | 1x | 2x | 4x |
|---|---|---|---|
| `r - R_p` [R_p] | 1.930e-4 | 9.601e-5 | 4.776e-5 |
| `T` [K] | 1169.1 | 1201.0 | 1272.0 |
| **photoheating, total** | **1.467e-8** | **1.422e-8** | **1.342e-8** |
| `heat_HeI` | 8.879e-9 (60.5%) | 8.607e-9 (60.5%) | 8.124e-9 (60.5%) |
| `heat_H2` | 5.781e-9 (39.4%) | 5.604e-9 (39.4%) | 5.289e-9 (39.4%) |
| `heat_HI` | 5.708e-12 (0.04%) | 5.704e-12 | 5.811e-12 |
| `heat_HeII`, `heat_He23S`, `heat_He_recomb`, both Penning channels | <= 1.3e-18 | <= 1.3e-18 | <= 1.3e-18 |
| `heat_metals`, `heat_H2_LW`, `heat_FUV_photolysis`, `heat_Hpe`, `heat_Hdx` | 0 | 0 | 0 |
| **radiative cooling, total** | **2.437e-6** | **2.718e-6** | **3.440e-6** |
| `H3p_IR` | 100.0000% | 100.0000% | 100.0000% |
| `reco`, `brem` | 1.4e-16, 4.5e-17 | 1.4e-16, 4.6e-17 | 1.5e-16, 5.0e-17 |
| **net radiative (heat - cool)** | **-2.423e-6** | **-2.704e-6** | **-3.426e-6** |
| cooling / heating | 166 | 191 | 256 |

Three readings. **(i)** Every photo channel is dr-independent to 9% -- these are
volumetric rates and they behave like it. **(ii)** They are irrelevant in size:
the total photoheating is 0.4-0.6% of the cooling, and the cooling is H3+
infrared to five figures. This run carries no metals (`abundance_* = 0` in
`EXHALE_resolved.out`), no Lyman-Werner field and no FUV photolysis, so those
channels are identically zero; secondary ionization is active (staged, the
default) but SvS85 partitions the photoelectron *inside* the photo channels
rather than as a column of its own, so it is contained in the `heat_HeI`,
`heat_H2` and `heat_HI` numbers above. **(iii) Radiation is a net sink at cell 1
on every grid**, and the sink deepens with refinement -- which follows the
temperature rather than causing it, since H3+ cooling is steeply temperature
dependent. **Nothing in the radiative budget raises the first cell above `T0`.**

### 8.2 The hydrodynamic terms: `dr`-proportional work, a `1/dr` residual

With no volumetric source available, the energy that holds cell 1 above `T0` has
to arrive through the base face. Two measurements separate the two ways that can
happen. `EXHALE_RESIDUAL=1` was run on each converged state to get the steady
residual row by row (`R = dF - S - (heat - cool)`, code units converted with
`q0 = n0 mu v0^3 / R_p = 1.409e-3 erg cm^-3 s^-1`); the gravitational work
`S_E = -rho v g` is computed from the profile with `g = GM_p/r^2` (the tidal term
is 1% of it at the base and is neglected).

| cell 1 [erg cm^-3 s^-1] | 1x | 2x | 4x | scaling |
|---|---|---|---|---|
| net radiative | -2.42e-6 | -2.70e-6 | -3.43e-6 | volumetric, `dr`-independent |
| gravitational work `-rho v g` | +2.58e-5 | +1.24e-5 | +5.76e-6 | **`~ dr`** (2.08, 2.15) |
| steady residual `R_energy` | +2.84e-5 | +5.22e-5 | +2.31e-4 | **`~ 1/dr` or steeper** (1.84, 4.43) |

The gravitational work is proportional to `v(1)`, so it inherits the first-order
scaling of section 7.2 and vanishes with the grid. The residual does the
opposite. Both restart cycles of each grid agree on that (`R_energy/u_3` reads
4.54e-6, 1.09e-6, 2.11e-6 s^-1 on 1x; 1.06e-6, 3.88e-6 on 2x; **1.89e-5,
1.72e-5** on 4x; the mass row behaves the same way, `R_mass` = 3.7e-2, 9.0e-3,
1.7e-2 / 6.3e-3, 3.1e-2 / **1.39e-1, 1.33e-1**), so the growth is not cycle
scatter.

**A residual that grows as `1/dr` in a cell of width `dr` is a fixed flux
mismatch at the face**, not a volumetric rate: multiplying back by the cell
volume gives a face-flux error that is roughly grid-independent. For mass,
`R_mass * dr` = 3.3e-6, 2.9e-6, 6.3e-6 in code units against a wind flux of
5.2e-6, i.e. **the base face carries a mass-flux error of order the wind flux
itself on every grid** -- the same conclusion section 3 reached by integrating
the residual from the wind down. For energy, `R_energy * dr * R_p` = 19, 17,
38 erg cm^-2 s^-1 (9.7, 4.7, 42 on the other cycle) against roughly
20 erg cm^-2 s^-1 of enthalpy flux that the wind itself carries through the base
(`Mdot` = 2.7e10 g s^-1 spread over `4 pi R_p^2`, at `5kT/2 mu` with
`T ~ 1200 K`). **So the answer to the discriminator is the face-flux one: what
sets the first cell is a boundary flux mismatch of the same order as the wind's
own flux, not a heating rate.**

### 8.3 Correction to section 7.3: `T(r)` does converge; the innermost 1e-3 `R_p` does not

Section 7.3 compared `T` at *cell 1*, which sits at `r - R_p = dr` and therefore
moves with the grid. Comparing the same three states at the same **physical**
radius instead:

| `r - R_p` [R_p] | 1e-4 | 2e-4 | 5e-4 | 1e-3 | 2e-3 | 4e-3 | 8e-3 | 1.5e-2 | 6e-2 |
|---|---|---|---|---|---|---|---|---|---|
| `T` 1x [K] | -- | 1167.1 | 1082.0 | 976.9 | 846.7 | 709.6 | 582.8 | 483.4 | 279.9 |
| `T` 2x [K] | 1199.2 | 1155.6 | 1057.3 | 948.4 | 813.3 | 725.7 | 688.4 | 634.6 | 346.8 |
| `T` 4x [K] | 1213.6 | 1138.2 | 1042.2 | 943.2 | 815.4 | 727.6 | 689.8 | 636.4 | 354.2 |

2x and 4x agree to 0.3% for `r - R_p >= 2e-3` and to 0.5% at 1e-3; they still
differ by 1.5% at 2e-4 and by 1.2% at 1e-4, and the 1x column is a different
solution altogether beyond `r - R_p ~ 4e-3` (see section 8.4). So `T(r)` **is**
grid-converging outside the innermost `~1e-3 R_p`, and section 7.3's "the base
thermal structure does not converge" was an artifact of sampling a steep profile
at a moving point. What does not converge is narrower: the innermost `~1e-3 R_p`
(5 cells at 1x, 20 at 4x), where the face-flux error of section 8.2 is absorbed.

The mismatch the ghost imposes is best stated at the face itself. Extrapolating
the first three interior cells quadratically to the base face
(`r - R_p = dr/2`), for both restart cycles of each grid:

| grid | `T_face` [K] | `rho_face` [1e14 cm^-3] | `p_face` [cgs] |
|---|---|---|---|
| pinned ghost | **1140.0** | **1.2204** | **8.9939** |
| 1x | 1199, 1209 | 1.235, 1.219 | 9.000 |
| 2x | 1225, 1268 | 1.208, 1.164 | 8.997 |
| 4x | 1302, 1506 | 1.135, 0.952 | 8.996 |

**The pressure at the base face matches the pinned ghost to better than 0.1% on
every grid; the temperature and the density do not, and the gap widens with
refinement** (and scatters between cycles on the finest grid, 1302 against 1506,
so the extrapolated face temperature is itself not converged). The interior
solution agrees with the boundary on *one* thermodynamic variable and disagrees
on the split of it: it wants a hotter, thinner base at the same pressure. That is
the shape of an over-specified lower boundary -- `rho` pinned to `rho_bc` and `T`
pinned to `T0` fixes two variables where the interior only supplies one.

### 8.4 What base refinement actually buys: a reproducible molecular layer

The restart-cycle scatter of section 4 turns out to be a function of the base
grid. Every state below is JFNK-converged (`info = 0`) in the same configuration;
`A_base` is the state of section 4, the others are cycles of the cold starts:

| grid | state | `log10 Mdot` | `r_front` | `T` at `r-R_p` = 2e-4 | 1e-3 | 4e-3 | 8e-3 |
|---|---|---|---|---|---|---|---|
| 1x | `A_base` | 10.29 | 1.0614 | 1178.4 | 970.9 | 715.8 | 657.1 |
| 1x | cold r1 | 10.37 | 1.0893 | 1142.9 | 993.0 | 979.3 | 949.2 |
| 1x | cold r2 | 10.28 | 1.0527 | 1167.1 | 976.9 | 709.6 | 582.8 |
| 2x | cold r1 | 10.28 | 1.0642 | 1147.0 | 943.9 | 733.2 | 697.1 |
| 2x | cold r2 | 10.29 | 1.0634 | 1155.6 | 948.4 | 725.7 | 688.4 |
| 4x | cold r1 | 10.29 | 1.0638 | 1131.7 | 942.9 | 729.1 | 691.5 |
| 4x | cold r2 | 10.29 | 1.0636 | 1138.2 | 943.2 | 727.6 | 689.8 |

**On the default base grid the converged molecular layer is not reproducible.**
Three `info = 0` states of the same configuration put the H2 front anywhere in
1.053-1.089 `R_p` (spread 0.037) and `T` at `r - R_p = 8e-3` anywhere in
583-949 K (a factor 1.63). One base refinement removes it: the 2x pair agrees to
0.0008 `R_p` in the front and 1.3% in `T`, the 4x pair to 0.0002 `R_p` and 0.3%,
and **the two grids agree with each other** (front 1.0638 against 1.0636, `T` at
8e-3 within 0.4%). `Mdot` is insensitive to all of it (10.28-10.37, i.e. +-0.05
dex, the du-stop path spread the project already documents).

This is the practically useful result of the refinement study: `Base grid
[dr,cells]: 1.0e-4 100` (with `Grid cells: 550`) is what makes the molecular
layer's thermal structure and H2 front a property of the physics rather than of
the restart cycle. The default grid gets `Mdot` right and the layer wrong.
### 8.5 Where 1140 K comes from, and at what level it is valid

`T0` is read from the `Equilibrium temperature` key
(`input_read.f90:143`), so EXHALE uses one number both as the planet's
equilibrium temperature and as the temperature the base ghost is pinned to. For
this planet that number is 1140 K, and `base.inp` repeats it (`T_base 1140.0`)
with the stated intent that "the handoff does not move the planet".

**It is not a `T_eq` guess, and it is not a VULCAN result.** The Tier-2 gate
planet is defined in `docs/lower_atmosphere_coupling.md` as "hot-Uranus-like:
0.0457 M_J, R_p = 0.49 R_J, T_eq = 1140 K, HD209 orbit/spectrum", and the same
document's Tier-1 validation gate names the source: "reproduce Koskinen 2022
Model A (hot Uranus, 0.05 au): r_0 = 1.34 R_p, T_0 = 1140 K, q_H2 ~ 0.84 at
1 microbar". That traces to the published paper, which states it directly
(Koskinen et al. 2022, ApJ 929, 52, section 3.2, "Numerical Modeling"):

> "The lower boundary of the model is at r0 = 1.34Rp (p0 = 10−6 bar, T0 = 1,140 K)
> where the boundary conditions are consistent with our models of the lower and
> middle atmosphere"

So 1140 K is a **lower-atmosphere model temperature at the 1 microbar level**, and
`q_H2 = 0.84` is its companion (`base.inp` uses 0.75, a representative
photochemical value, deliberately below the equilibrium fit -- that is what the
regression case exists to test).

> **Corrected 2026-09-03 (section 9.1).** The "level" half of what follows does
> not survive. Koskinen et al. state in the same paper (section 2.2) that their
> lower and middle atmosphere is **isothermal at the effective temperature**, so
> 1140 K is their value at 1 microbar, at 9 microbar and anywhere below their
> boundary -- there is no profile to be higher up on. Our own Tier-1 column says
> the same. What survives is the *radius* half: the gate planet's base radius is
> Koskinen's 1 microbar radius while its hand-chosen `n0` puts the base pressure
> nine times deeper. Read section 9.1 in place of item 1 below.

**Two mismatches follow, and they both point the same way.**

1. **Level.** Koskinen's 1140 K is quoted at `p0 = 1e-6 bar`. This run's base sits
   at `p = 8.994 dyn cm^-2 = 8.99e-6 bar`, i.e. **nine times deeper**, and
   `base.inp` declares `p_base 1.0e-5 bar` for exactly that reason. Koskinen's
   own profile has the temperature *decreasing* with radius just above the lower
   boundary ("the temperature first slightly decreases with the radius, and then
   increases rapidly to about 4400 K around r ~ 3 Rp"), which is the same sense
   the EXHALE solution shows below `r - R_p ~ 1e-2`. Going nine times deeper
   along a profile that falls outward means the level in question should be
   **hotter** than 1140 K, not equal to it.
2. **Radius.** Koskinen's lower boundary is at `r0 = 1.34 R_p`; EXHALE's is at
   `1.00 R_p`, because `Planet radius: 0.49 R_J` is the base radius here rather
   than the 1-bar radius. The two boundaries are not the same surface.

**Judgment, tentative.** The interior solution asking for 1200-1300 K at the base
face (section 8.3) is not obviously in conflict with the published physics -- it
is what one would expect from applying a 1-microbar temperature at a
nine-microbar level on a profile that decreases outward. What the code does is
pin the published 1-microbar value at a deeper level, and the interior then
disagrees with it by 5-14%. That makes the "physically correct base temperature"
question answerable in principle: it is the temperature of the Koskinen lower and
middle atmosphere model **at the level EXHALE actually places its boundary**, not
the value published for 1 microbar. Reading that number off Koskinen's Figure 7
was not attempted here; it is the obvious next step and it does not need a code
change.
### 8.6 A deeper lower boundary makes the mismatch worse, not better

If the interior wants a hotter base than `T0`, and `T0` is a value published for a
shallower level, then moving the boundary *deeper* should widen the gap rather
than close it. That is what happens. `Log10 lower boundary number density` was
raised 14.00 -> 14.50 (base pressure 8.99e-6 -> 2.844e-5 bar, a factor 3.16)
with everything else unchanged -- same default 1x grid, same `base.inp`, so the
composition handoff is untouched (`EXHALE_resolved.out` differs only in the last
digits of `ntot_bc_per_H`). Cold start to the 200000-step budget (which did not
reach `info = 0`, `du = 3.5e-2`), then two chained restarts, both `info = 0` and
agreeing with each other to four digits.

| quantity | shallow (`n0 = 1e14`) | deep (`n0 = 10^14.5`) |
|---|---|---|
| base pressure at the ghost | 8.99e-6 bar | 2.844e-5 bar |
| `T_ghost` (pinned) | 1140.0 K | 1140.0 K |
| `T(1)` | 1169.1 K | 1175.3 K |
| **`T(1) - T_ghost`** | **+29.1 K** | **+35.3 K** |
| `rho(1)/rho_ghost` | 1.0333 | 1.0289 |
| `v(1)` | -248.3 cm/s | -247.9 cm/s |
| cell 1 `rho v r^2 / F_wind` | -202.8 | **-588.7** |
| `v` at cells 5, 10 | -0.05, +0.76 cm/s | +7.9, -18.0 cm/s |
| odd-even purity, cells 1-6 | 1.000 | 0.825 |
| radiative net at cell 1 | -2.42e-6 (cool/heat 166) | -1.96e-6 (cool/heat 222) |
| `\|\|R\|\|` | 2.48e-4 | 3.59e-4 |
| **`log10 Mdot`** | **10.28** | **10.31** |
| **`r_front`** | **1.0527** (cycles 1.053-1.089) | **1.08302** (both cycles) |

**(a) The temperature mismatch does not shrink; it grows**, +29.1 -> +35.3 K.
Deepening the boundary by a factor 3.16 in pressure does not bring the pinned
1140 K into agreement with what the interior wants. (The explanation offered
here at the time -- that a deeper level should be hotter along Koskinen's
profile -- is withdrawn in section 9.1: their lower atmosphere has no profile.
The measured trend is simply monotone in base pressure, +35.3, +29.1, +8.9 K at
28, 9 and 1 microbar; see section 9.3.) The radiative picture there is
unchanged: still a net sink, still H3+ to five figures, cooling/heating 222.

**(b) The sawtooth is unchanged at cell 1 and worse above it.** `v(1)` is the
same to 0.2% (-247.9 against -248.3), but the base is 3.16 times denser, so the
cell-1 mass-flux ratio is 2.9 times worse (-589 against -203), and the
alternation no longer dies at cell 5: it is still +-18 cm/s at cell 10, against
below 1 cm/s in the shallow run. This is the dense-base behaviour
`docs/base_breathing_progress.md` already records from the other direction (a
*lighter*, 1-microbar-anchored base gave the cleanest base of its 2x2 sweep).

**(c) The wind moves a little, and becomes reproducible.** `log10 Mdot` goes
10.28 -> 10.31 and the H2 front 1.053 -> 1.083 `R_p`. Both are modest, but the
notable part is that the deeper base gives the **same front to five digits in two
independent restart cycles** (1.08302, 1.08302) on the default grid, where the
shallow base scattered over 1.053-1.089. So the deeper boundary and the finer
base grid buy the same thing -- a molecular layer that is a property of the
physics rather than of the restart cycle -- by different routes, and the deeper
one costs a wider sawtooth.
### 8.7 Candidate judgment, updated

The question section 7 left open was whether the first cell is raised above `T0`
by a real heating rate or by a boundary flux mismatch. **On these measurements it
is the flux mismatch.**

* **Not a volumetric heating rate.** Every photo channel at cell 1 is
  `dr`-independent to 9%, and their sum is 0.4-0.6% of the local cooling.
  Radiation is a **net sink** of 2.0-3.4e-6 erg cm^-3 s^-1 at cell 1 on every
  grid and at both base depths. There is no heating term to blame.
* **A face-flux mismatch, of order the wind's own flux.** The cell-1 steady
  residual grows as `1/dr` in both the mass and the energy row while the
  gravitational work at the same cell falls as `dr`; multiplied back by the cell
  volume the residual is a roughly grid-independent face-flux error, about the
  wind mass flux for mass and about the wind's enthalpy flux for energy. The
  base face is delivering an error of the same size as the thing it is supposed
  to deliver.
* **The over-specification is visible directly.** Extrapolated to the base face,
  the interior agrees with the pinned ghost on `p` to better than 0.1% on every
  grid and disagrees on `T` (by 5-14%, growing with refinement) and on `rho` (by
  the compensating amount). The ghost fixes two thermodynamic variables; the
  interior supplies one.
* **~~`T0` itself is applied at the wrong level.~~** *Withdrawn 2026-09-03 --
  see section 9.1.* 1140 K is Koskinen et al. (2022) Model A's lower-boundary
  temperature at `p0 = 1e-6 bar`, but their lower and middle atmosphere is
  isothermal at that value, so it is equally their number at 9e-6 bar. What is
  inconsistent is the pressure the gate puts at Koskinen's own base radius, not
  the temperature; and the thermal mismatch, while monotone in base pressure
  (+35.3, +29.1, +8.9 K at 28, 9, 1 microbar), turns out to be independent of the
  flux error (section 9.4).
* **Corrected from section 7.3**: `T(r)` at fixed physical radius *does*
  converge outside the innermost `~1e-3 R_p`; the apparent divergence was the
  first cell centre moving inward along a steep profile. The innermost `1e-3
  R_p` does not converge, and that is where the face-flux error is absorbed.
* **The one clearly actionable measurement** is section 8.4: on the default base
  grid the converged molecular layer is not reproducible across restart cycles
  (front 1.053-1.089 `R_p`, `T` at `r - R_p` = 8e-3 spanning 583-949 K), and one
  base refinement (`Base grid [dr,cells]: 1.0e-4 100` with `Grid cells: 550`)
  makes it so, in agreement with the 4x grid. `Mdot` never notices.

No prescription is proposed here.

---

## 9. What the published base temperature actually is, and what a different base condition does (2026-09-03)

Section 8.5 argued that `T0 = 1140 K` is Koskinen et al. (2022) Model A's value
for the 1 microbar level, applied here at 9 microbar, and read the interior's
1200-1300 K as what a deeper level should have. **The first part of this section
refutes the second half of that argument.**

### 9.1 The published temperature at 9 microbar is 1140 K -- their lower atmosphere is isothermal

Koskinen et al. (2022) state the assumption directly (section 2.2, "Lower and
Middle Atmosphere"):

> "We assume that the temperature is constant throughout the lower and middle
> atmosphere and equal to the effective temperature of the planet. In order to
> calculate the effective temperature, we use a generic Bond albedo AB = 0.3 and
> assume uniform redistribution of energy around the planet."

So their lower and middle atmosphere has **no depth dependence at all**: the
1140 K quoted for `p0 = 1e-6 bar` is the same number at 9e-6 bar, at 1e-3 bar,
anywhere below their escape-model boundary. (1140 K is reproduced by that recipe:
a Sun-like star at 0.05 au with `A_B = 0.3` and full redistribution gives
`T_eff = 5772 (R_sun/2a)^(1/2) (1-0.3)^(1/4) = 1138 K`.) Our Tier-1 column
implements the same section and says so in its own header
(`lower_column.f90`: "isothermal at `T_eq`", "Approximations (documented):
isothermal column at `T_eq`"), so it returns 1140 K at 9 microbar too.

**Digitization of Figure 7** (their *upper* atmosphere solution, whose domain
begins at the 1 microbar boundary) was done anyway, because it says what their
converged interior does just above the same boundary. Method: `pdftoppm -r 600
-f 9 -l 9`, page 9, the single panel of Figure 7 (`T` solid on the left axis, `V`
dashed on the right). Plot box columns 1325 and 3747 carry `r/R_p` = 1 and 10;
the left-axis label rows are 464, 703, 942, 1181, 1419, 1658 for 6000 down to
1000 K, i.e. 238.8 px per 1000 K, so **the box bottom (row 1776) is 506 K, not
zero** -- reading it as zero is the trap, and it produced one wrong pass of this
measurement before the label rows were located. Validation against the paper's
own text: the digitization gives 4384 K at `r = 3.0 R_p` against "about 4400 K
around r ~ 3 Rp", 5538 K at the upper boundary against "about 5540 K", and
1149 K at the first plotted point against the stated `T0 = 1140 K` -- better
than 1% on all three.

| source | what it is | `T` at 1e-6 bar | `T` at 9e-6 bar (EXHALE's base level) |
|---|---|---|---|
| Koskinen+2022 sec. 2.2, lower and middle atmosphere | isothermal at `T_eff` (`A_B` = 0.3, uniform redistribution) | 1140 K | **1140 K** (isothermal by construction) |
| Koskinen+2022 Fig. 7, upper atmosphere (digitized here) | their converged escape solution; domain starts at 1 microbar | 1149 K at the first plotted point (`r` = 1.331); **minimum 860 K at `r` = 1.372**, then rising | not defined (below their domain) |
| EXHALE Tier-1 column, `lower_column.f90` | the same section 2.2 equations | 1140 K | **1140 K** |
| `docs/lower_atmosphere_figs/data_g2` | EXHALE 12000-step relaxation snapshot; `base.inp` records that no VULCAN run exists for this planet | -- | `T`(cell 1) = 1152.8 K against a 1155.1 K ghost, i.e. still at `T0` and not yet relaxed |
| **EXHALE converged interior (section 8.3)** | face extrapolation of the JFNK-converged states | -- | **1199, 1209 (1x); 1225, 1268 (2x); 1302, 1506 (4x)** |

**Reading.** No published source supports 1200-1300 K at this level, and section
8.5's "the value is applied at the wrong level" argument does not survive: within
Koskinen's framework 1140 K is the right number at *any* level below their
boundary, because the model has no gradient there. Worse for the earlier
argument, their own escape solution moves the **opposite way** just above the
boundary -- down to 860 K at `r = 1.372 R_p` (a 25% dip, which is what the paper
calls "the temperature first slightly decreases with the radius") -- while EXHALE
goes *up*, to 1169-1272 K in its first cell. The two codes disagree in sign about
what happens immediately above a 1140 K lower boundary.

**The disagreement is not confined to the first cell.** Putting both codes on
`(r - r_base)/r_base` (Koskinen's radii are in units of their 1-bar radius,
2.5559e9 cm, and their boundary is at 1.34 of it, so their axis is scaled by
2.5559/3.4256) gives:

| `(r - r_base)/r_base` | 1e-3 | 3e-3 | 1e-2 | 3e-2 | 1e-1 |
|---|---|---|---|---|---|
| Koskinen+2022 Fig. 7 | 1137 | 1113 | 1044 | **862** | 1554 |
| EXHALE, 1 microbar base | 1030 | 849 | 609 | **394** | 959 |
| EXHALE, 9 microbar base (default) | 977 | 766 | 545 | **400** | 746 |

Both codes dip and then rise, but EXHALE's converged molecular layer is
100-450 K colder than Koskinen's over the same range, and its minimum is ~400 K
against their 860 K. The project already records the likely reason on the other
side of the ledger -- H3+ carries more than 99% of the radiative cooling from the
base to the front here (`docs/lower_atmosphere_coupling.md`, and section 8.1 of
this note measures 100.0000% at cell 1) -- so the comparison above is a
quantitative handle on that, and it is not a base-condition question. It is
recorded here and not pursued.

One genuine inconsistency does survive, and it is about pressure rather than
temperature. Koskinen's lower boundary is `r0 = 1.34 R_1bar` = 3.425e9 cm with
`R_1bar` = 2.5559e9 cm (`lower_column.f90` header), and the gate planet's
`Planet radius: 0.49 R_J` = 3.4256e9 cm is **the same surface to 0.2%**. So
EXHALE's base radius is Koskinen's 1 microbar radius -- but the hand-chosen
`Log10 lower boundary number density: 14.00` puts the pressure there at
8.99e-6 bar, nine times their value. The self-consistent choice is
`Base BC: pressure 1.0`, which derives `n0` = 1.1119e13 cm^-3; that run is
section 9.3.
### 9.2 `Base BC: pressure` is a base-depth key, not a change of ghost closure

The over-specification of section 8.3 -- `rho` pinned to `rho_bc` *and* `T`
pinned to `T0`, where the interior only agrees on `p` -- was expected to be
testable with `Base BC: pressure`. It is not. In the current source the key does
exactly one thing (`input_read.f90`, the only two places `base_bc_mode` is read
outside the setup report):

```
if (base_bc_mode .eq. 1) then
   n0 = base_p_ubar/(kb_erg*T0*ntot_bc)
```

The lower ghost is still written by `BC_component_constrho`, which pins
`W(1) = rho_bc` and `W(3) = ntot_bc + dp_bc` exactly as before. So the key
**derives `n0` from a target base pressure at `T0`** and changes nothing else:
it moves the boundary in depth, it does not change what is imposed there.

Measured confirmation: `Base BC: pressure 8.994` (the pressure the density
anchor already produces) derives `n0 = 1.0000e14`, and a 1000-step cold run
under it reproduces the density-anchored run of the same length to the round-off
in `n0` -- ghost `rho` 1.220429e14 against 1.220421e14, and at most 0.1% in
`rho`, `p` and `T` anywhere in the transient state.

The keys that *do* change the closure are `Base ghost temperature: continuous`
(`T` floats, `rho` pinned) and `Hydrostatic base: True` (`p` extrapolated, `rho`
pinned); both were run in section 4 from a converged restart and both failed.
Section 9.4 runs the first of them from a cold start, which is the fair test.
### 9.3 Base depth is the dominant lever, and 1 microbar is the self-consistent choice

Section 9.1 ends with an inconsistency in the gate setup: the base *radius*
(`Planet radius: 0.49 R_J` = 3.4256e9 cm) is Koskinen's 1 microbar radius
(`1.34 R_1bar` = 3.425e9 cm) to 0.2%, while the hand-chosen
`Log10 lower boundary number density: 14.00` puts the *pressure* there at
8.99e-6 bar. `Base BC: pressure 1.0` is the self-consistent choice; it derives
`n0 = 1.1119e13 cm^-3`. Run from a cold start on the default 1x grid and chained
twice, together with the 9 microbar baseline of section 8 and the 28 microbar run
of section 8.6 (all `info = 0`, all second restart cycles unless noted):

| base | `p` at the ghost | `rho_bc` [cm^-3] | `T(1)` | **`T(1) - T_ghost`** | `rho(1)/rho_ghost` | `v(1)` [cm/s] | **cell 1 `/F_wind`** | `log10 Mdot` | `r_front` |
|---|---|---|---|---|---|---|---|---|---|
| 1 microbar | 1.000e-6 bar | 1.357e13 | 1148.9 K | **+8.9 K** | 1.0151 | -227.8 | **-21.0** | 10.27 | 1.03551 |
| 9 microbar (default) | 8.994e-6 bar | 1.220e14 | 1169.1 K | **+29.1 K** | 1.0333 | -248.3 | **-202.8** | 10.28 | 1.0527 |
| 28 microbar | 2.844e-5 bar | 3.859e14 | 1175.3 K | **+35.3 K** | 1.0289 | -247.9 | **-588.7** | 10.31 | 1.08302 |

Three things fall out.

**(a) The velocity artifact is a fixed velocity, and its flux impact scales with
the base density.** `v(1)` is -228, -248, -248 cm/s across a 28x range in base
pressure -- constant to 8% -- while the cell-1 mass-flux ratio goes -21.0,
-202.8, -588.7. The ratios of the flux error (9.7 and 2.9) match the ratios of
`rho_bc` (9.0 and 3.16). So the base-face error is `v_err * rho_bc` with `v_err`
set by the discretization, and **the only way to shrink it without refining the
grid is to make the base lighter.**

**(b) The temperature mismatch shrinks 3.3x at 1 microbar**, +29.1 -> +8.9 K, and
`T(1) = 1148.9 K` there lands on the 1149 K that the Figure 7 digitization reads
at Koskinen's own first plotted point. That agreement may be partly coincidence
-- the two codes disagree in sign about what happens above the boundary (section
9.1) -- but the residual mismatch at the correct level is 0.8% of `T0`, against
2.6% at 9 microbar.

**(c) It also fixes the section 8.4 reproducibility problem, on the default
grid.** Two independent chained restarts of the 1 microbar run give
`r_front = 1.03551` and `1.03551`, `log10 Mdot` 10.27 and 10.27, `v(1)` -228.12
and -227.80, peak `n(H3+)` 2.1636e5 both -- five digits, where the 9 microbar
base on the same grid scattered over 1.053-1.089 `R_p`. And it converges from
cold in **26 m 42 s / 66712 steps**, against the 9 microbar run failing to reach
`info = 0` in 200000 steps and 1 h 41 m.

The cost is that the molecular layer thins with the base: peak `n(H2)` falls
6.4e13 -> 6.15e12 and peak `n(H3+)` 3.73e5 -> 2.16e5 cm^-3, and the H2 front moves
in from 1.053 to 1.0355 `R_p`. `Mdot` is almost untouched by any of it (10.27,
10.28, 10.31 over a 28x pressure range, i.e. 0.04 dex).
### 9.4 Changing `T0`, and changing the closure: the temperature gap and the flux error are independent

Two further cold starts on the default 1x grid at the default 9 microbar base,
each chained to `info = 0` where it could reach it.

**`Equilibrium temperature: 1200.0`** (with `T_base 1200.0` in `base.inp` so the
handoff does not put 1140 K back). 1200 K is what the interior's own face
extrapolation asked for at 1x/2x (section 8.3), so this asks whether pinning the
ghost at the temperature the interior wants removes the mismatch and the flux
error. **Caveat**: `T0` is also the code's temperature scale and the equilibrium
temperature, so this changes more than the ghost -- the IC, `v0`, `q_H2` ceiling
and the whole thermal normalization move with it.

**`Base ghost temperature: continuous`** from a cold start -- the closure that
lets `T_ghost` float to `T(1)` and pins only `rho`. Section 4 ran it from a
converged restart, where it failed; a cold start is the fair test.

| state | `T_ghost` | `T(1)` | **`T(1) - T_ghost`** | `rho(1)/rho_ghost` | `v(1)` | **cell 1 `/F_wind`** | `log10 Mdot` | `r_front` |
|---|---|---|---|---|---|---|---|---|
| baseline, 9 microbar | 1140.0 | 1169.1 | **+29.1** | 1.0333 | -248.3 | **-202.8** | 10.28 | 1.0527 |
| `T0 = 1200`, cycle 1 | 1200.0 | 1020.7 | -179.3 | 1.2456 | -241.1 | -234.9 | 10.28 | 1.0538 |
| **`T0 = 1200`, cycle 2** | 1200.0 | 1191.3 | **-8.7** | 1.0675 | -241.2 | **-203.8** | 10.28 | 1.0533 |
| ghost `T` continuous, cycle 1 (`info = 0`) | 1501.3 | 1501.3 | 0 by construction | 1.0555 | +112.5 | +88.2 | 10.31 | **1.0009** |
| ghost `T` continuous, cycle 2 (`info = 2`) | 1839.9 | 1839.9 | 0 by construction | 0.9985 | +4113 | +3131 | -- | 1.0384 |
| 1 microbar base (section 9.3) | 1140.0 | 1148.9 | **+8.9** | 1.0151 | -227.8 | **-21.0** | 10.27 | 1.0355 |

**The `T0 = 1200` result is the decisive one.** Raising the pin to the
temperature the interior asked for **does** close the thermal gap -- +29.1 K
becomes -8.7 K, a 3.3x reduction and a sign flip -- and **does not touch the
base-face mass-flux error at all**: -202.8 becomes -203.8, i.e. 0.5%, with
`v(1)` moving 3% and `Mdot`, the front and peak `n(H3+)` unchanged to three
digits. So the temperature mismatch of section 8.3 and the flux mismatch of
section 8.2 are **two independent defects**, and the first is not the cause of
the second. (Note also that the interior does not simply track the pin: with
`T0 = 1200` it asks for 1191 K, i.e. it now wants slightly *less* than the
boundary offers, so the +29 K of the baseline is not a fixed offset either.)

**The continuous ghost is not a prescription for this planet.** From a cold start
it converges once, but to a base at 1501 K with the H2 front collapsed onto the
boundary (1.0009 `R_p`), a velocity field alternating +112, -98, +48, -106,
+701 cm/s over cells 1-5 and a 4116 K spike at cell 5; the next restart cycle
does not converge at all (`info = 2`, `||R||` = 1.55e-2, ghost at 1840 K,
`v(1)` = +4113 cm/s). Removing the temperature pin removes the only thing
holding the base thermal state, and the layer runs away -- the same failure the
restart test of section 4 found, reached from the opposite direction.

Worth recording: the continuous-ghost cycle-1 state is the **only** configuration
in this whole study in which `v(1)` is positive, i.e. the base face is a
cell-centred inflow. It costs the molecular layer.
### 9.5 Cost and consequences of making `1.0e-4 100` the molecular default

Section 8.4 is the case for it: on the default base grid the converged molecular
layer is not reproducible across restart cycles (H2 front 1.053-1.089 `R_p`, `T`
at `r - R_p` = 8e-3 spanning 583-949 K), and one refinement makes it reproducible
to 0.0008 `R_p` and 1.3% and puts it in agreement with the 4x grid. What follows
is the cost, measured.

**Per step: +12%.** Three alternating pairs of 2000-step cold runs, same binary,
`OMP_NUM_THREADS=8`, run back to back so they see the same machine load:
1x 24.95, 24.22, 26.38 ms/step; 2x 25.92, 27.98, 28.21 ms/step. Medians 24.95
against 27.98, i.e. **+12%**, which is the cell-count ratio (500 -> 550, +10%)
plus noise. (Two earlier pairs at 3000 steps gave +23% and +3%; the machine was
loading and unloading, which is why the alternating protocol was used.)

**Per unit physical time: exactly 2x the steps.** The base cells set the CFL
minimum, so halving `dr_base` halves `dt`. Measured on the cold starts by the two
milestones the log prints, which are physical states rather than step counts:

| milestone | 1x | 2x | 4x | ratio |
|---|---|---|---|---|
| PLM -> WENO3 at `du` = 0.5 | step 17652 | 35471 | 71217 | 2.009, 2.008 |
| secondary ionization at `du` = 2e-2 | step 40557 | 81731 | 156005 | 2.015, 1.909 |

So a *fixed-step* run covers half the physical time and a *fixed-time* run costs
**2.01 x 1.12 = 2.25** times the wall clock.

**For a run driven to convergence the finer grid was cheaper.** The 1x cold start
did not reach `info = 0` in its 200000-step budget (1 h 41 m); the 2x cold start
converged at step 83731 (47 m) and the 4x at 158424 (1 h 37 m). On this planet
the refinement paid for itself.

**Which goldens move: all six molecular cases, and only those.** The default
matrix is `wasp_full wasp_he23off mol_base_handoff mol_metals mol_lyman_werner
mol_diffusion mol_ir_bands mol_sec_ion lower_profile`; exactly six carry
`Molecular chemistry: True` (the six `mol_*`), and `lower_profile`, `wasp_full`
and `wasp_he23off` do not. **None of the six sets `Base grid` or `Grid cells`**,
so all six would inherit a changed default and all six goldens move.

Two things make this more than a re-snapshot:

* **All six are 12000-step relaxation snapshots** (`<case>/maxsteps` = 12000),
  not converged solutions. With `dt` halved, 12000 steps cover **half the
  physical time**, so the golden would move for two reasons at once. Preserving
  the snapshot convention means going to 24000 steps, i.e. **2.24x the wall time
  for those six cases**; keeping 12000 steps means the cases now pin a different
  physical state and the convention has to be restated.
* **The row count changes**, 500 -> 550 interior cells (554 -> 604 output rows
  with ghosts), so the goldens cannot be diffed against the old ones at all --
  they have to be re-snapshotted with `run_check.sh golden` after a `check`.

**One design consequence to weigh.** Gating a *grid* default on a *physics* key
means the same planet is run on different grids with molecules on and off, which
is exactly the comparison Gate 0 of `docs/lower_atmosphere_coupling.md` makes
("mol-off byte-equivalent: molecular chemistry off reproduces the atomic
systems"). That gate would no longer be a like-for-like comparison. Putting the
key explicitly in the molecular `examples/` and regression inputs, rather than
making it a conditional default, keeps the grid a property of the input file.
Which way to go is not decided here.
### 9.6 Candidate judgment, updated

**Corrections to section 8 first.** Section 8.5 said `T0 = 1140 K` is applied
"at the wrong level". That does not survive section 9.1: Koskinen's lower and
middle atmosphere is isothermal at `T_eff`, so 1140 K is their value at 1
microbar, at 9 microbar and everywhere below their boundary; our Tier-1 column
says the same, and no source supports the 1200-1300 K the interior asks for. The
real inconsistency in the gate setup is not the temperature but the **pressure**:
the base radius is Koskinen's 1 microbar radius while the hand-chosen `n0` puts
the base pressure nine times deeper.

**Where the defects now stand, and what moves each.**

| defect | what it is | what moves it | what does not |
|---|---|---|---|
| cell-1 velocity artifact, `v(1) ~ -248 cm/s` | first-order-in-`dr` boundary discretization error (section 7.2) | base refinement (`~ dr`) | base depth (-228 to -248 over 28x in pressure), `T0`, the ghost closure |
| base-face **mass-flux** error, -21 to -589 `F_wind` | that velocity times `rho_bc` | **base depth** (`~ rho_bc`, measured 9.7x and 2.9x against 9.0x and 3.16x) and base refinement | `T0 = 1200` (-202.8 -> -203.8), Shapiro, `Base velocity: massflux`, `Low-Mach damping` |
| base **thermal** mismatch, `T(1) - T_ghost` | the isothermal pin disagreeing with the interior | `T0` itself (+29.1 -> -8.7 K at `T0` = 1200), base depth (+35.3, +29.1, +8.9 K at 28, 9, 1 microbar) | base refinement (it grows, though section 8.3 shows that is mostly sampling) |
| molecular layer not reproducible across restart cycles | front 1.053-1.089 `R_p`, `T` at 8e-3 spanning 583-949 K on the default grid at 9 microbar | **either** 2x base refinement **or** a 1 microbar base -- both give five-digit reproducibility | nothing else tested |

**The two independent defects.** The `T0 = 1200` run separates them cleanly: it
closes the thermal gap by 3.3x and leaves the mass-flux error at 0.5% of its
former value -- unchanged. Any prescription aimed at one of them should not be
expected to fix the other.

**Base depth is the strongest single lever, and 1 microbar is also the
self-consistent choice.** Moving the base to 1 microbar -- the level the gate
planet's own base radius corresponds to, and the level the published `T0` refers
to -- cuts the base-face mass-flux error by 9.7x, the thermal mismatch by 3.3x,
removes the velocity alternation beyond cell 2, makes the molecular layer
reproducible to five digits **on the default grid**, converges from cold in
27 minutes instead of not converging in 1 h 41 m, and leaves `log10 Mdot` at
10.27 against 10.28. Its cost is a thinner molecular layer (peak `n(H2)` 6.4e13
-> 6.15e12 cm^-3, front 1.053 -> 1.0355 `R_p`) -- which is a physical
consequence of a lighter base, not an error.

**What does not work**: `Base ghost temperature: continuous` (section 9.4 --
converges once to a 1501 K base with the front collapsed, then stops converging),
`Hydrostatic base: True` (section 4), and the Shapiro filter (section 5).
`Base BC: pressure` is not a closure change at all (section 9.2), so the
over-specification of section 8.3 has **no key in the current code that tests it
directly**; a ghost that pins `p` and lets `rho` and `T` float together is the
missing option, and it does not exist.

**On the grid default**: section 9.5 gives the measured cost (+12% per step,
2.01x the steps per unit physical time, all six molecular goldens move and all
six are 12000-step snapshots whose physical duration would halve). Note that a
1 microbar base buys the same reproducibility at no step cost -- it converges
faster -- so the two candidate defaults are alternatives, not complements, and
the choice between them is a physics choice about where the boundary belongs.

No prescription is chosen here.

---

## 10. The base level now has one source, and the gate moved to 1 microbar (2026-09-03)

Section 9 ended with the inconsistency that survives: the gate planet's base
radius is Koskinen's 1 microbar radius while its hand-chosen `n0` put the base
pressure nine times deeper, and `Base BC: pressure` -- the key that looked like
it would settle it -- only moves the level, it does not change the closure. The
user approved putting the base at 1 microbar and letting `base.inp` state it
(2026-09-03). This section is what was changed and what it moved.

### 10.1 The input change: `p_base` fixes the level, and the pair is checked

`input_read.f90` now treats a `base.inp` that carries `p_base` as what it is --
a lower-atmosphere handoff written AT that pressure, whose composition, `q_H2`
and temperature all refer to that level -- so the level is the handoff's to
state:

```
n0 = p_base/(k_B T0 ntot_bc)
```

the same relation `Base BC: pressure` already used. Three rules follow, and all
three are refusals rather than warnings, because a run that marches on a base
level nobody chose is not a model of anything:

* `Log10 lower boundary number density` became **optional**. If it is given
  alongside `p_base` the two must agree to **1%**; if they do not, startup stops
  and names both numbers, the `n0` each implies, `ntot_bc`, and the two ways to
  fix it (delete one or the other).
* `p_base` and `Base BC: pressure` together must be the same level (round-off
  tolerance), or startup stops.
* If no source states the level at all, startup stops and lists the three keys
  that can.

A lower-atmosphere **profile** is inside the same rule. It states the level
itself: `p_base_bar` is the profile's matching pressure `p_match_bar`, the level
`apply_lower_atmosphere_profile` read `T0`, `R0`, `q_H2` and every elemental
ratio at, and `n0` follows from it by the same relation. A density key given
beside a profile therefore has to agree with it to 1% or startup stops, exactly
as beside `p_base` (`input_read.f90`, `base_level_from_handoff`). The
`lower_profile` case had a density key that disagreed; the key was removed from
`examples/17_lower_profile/input.inp` and from the case, so the run now starts
at the profile's own level and its solution changed (base density 1.0e14 ->
5.0e12 cm^-3, printed log10 Mdot 8.44 -> 7.92 for that change alone; measured
in `docs/Update_EXHALE_stage2.md` section 6, item 2c-BASE, not re-measured here).

`write_setup_report.f90` carries one line naming both numbers and the source,
which is `base.inp p_base`, `the profile matching level p_match`,
`"Base BC: pressure"` or `the density key`:

```
 - Base level: n0 =  1.1690E+13 cm^-3 -> p =  1.0000E-06 bar (level from base.inp p_base)
```

### 10.2 Every `p_base` in the repository disagreed with its density key

The consistency check was run over all 649 `base.inp` files in the tree. 597
carry no `p_base`; those of them that hand a **profile** over (the LHS 1140 b
examples, `lower_profile`) state the level through the profile's `p_match_bar`
and are inside the same rule, and the rest state it through the density key.
**Of the 52 that do carry `p_base`, all 52 disagreed** -- the key had
never actually located the base level anywhere:

| group | files | `p_base` | `p` from the density key | ratio |
|---|---|---|---|---|
| `mol_*` gate cases (metals off) | 4 | 1.0e-5 bar | 8.994e-6 | 0.899 |
| `mol_metals`, `mol_ir_bands` (metals on) | 2 | 1.0e-5 | 9.007e-6 | 0.901 |
| `armHeH_0p3` | 1 | 1.0e-5 | 1.014e-5 | 1.014 |
| `armA_*`, `armD_*`, `arm_heh1_x2matched` (He/H = 1) | 12 | 1.0e-5 | 1.210e-5 | 1.210 |
| `armHeH_3` / `_10` / `_30` | 3 | 1.0e-5 | 1.392e-5 / 1.508e-5 / 1.551e-5 | 1.39-1.55 |
| `vulcan_work/**` HD 209458 b runs | 30 | 1.0e-6 or 1.0e-4 | 1.4e-5 to 3.2e-5 | **0.14-32.2** |

The 22 `backup/regression` cases were retargeted (section 10.3). The 30
`vulcan_work` directories were **not touched** -- they are records of past
VULCAN comparison runs, not part of the matrix -- but every one of them will now
be refused if it is re-run, and their stated `p_base` is wrong by factors of 7
to 32 against the level they actually used. That is a finding about those runs,
not about this change.

### 10.3 The gate moved to the level its own radius already implied

The six molecular regression cases and the sixteen diagnostic runs now read

```
q_H2_base 0.84            (gate, He/H = 0.0793)
p_base    1.000e-06
```

with `Log10 lower boundary number density` deleted from `input.inp`. 0.84 is
Koskinen et al. (2022, ApJ 929, 52, section 3.2) at the same level: "The volume
mixing ratio of H2 at the lower boundary is about q0 = 0.84 while H is a minor
species with q0 = 0.026." It passes the ceiling `0.5/(0.5 + He/H)` = 0.863110651
at 97.3% of it.

The He/H ladder keeps its own rule -- **fix the fraction of H nuclei bound into
H2 and vary only He/H** -- so its `q_H2_base` values were recomputed from the new
gate value through `x2 = 2 q (1+He/H)/(1+q)` and its inverse
`q = x2/(2(1+He/H) - x2)`, with `x2 = 0.985447826` now (it was 0.925114286 when
the gate said 0.75):

| He/H | `q_H2_base` | ceiling `0.5/(0.5+He/H)` | fraction of ceiling | `x2` check |
|---|---|---|---|---|
| 0.0793 (gate) | **0.840000000** | 0.863110651 | 97.3% | 0.985447826 |
| 0.3 | 0.610353658 | 0.625000000 | 97.7% | 0.985447826 |
| 1.0 | 0.326896922 | 0.333333333 | 98.1% | 0.985447826 |
| 3.0 | 0.140486207 | 0.142857143 | 98.3% | 0.985447826 |
| 10.0 | 0.046893592 | 0.047619048 | 98.5% | 0.985447826 |
| 30.0 | 0.016151029 | 0.016393443 | 98.5% | 0.985447826 |

All six pass the ceiling. The `x2` column is the code's own
`base_h2_nuclei_fraction()` formula evaluated on each pair, so the ladder is one
experiment in He/H and not two.

**The 12000-step snapshot convention survives the move, measured.** At the same
12000 steps the retargeted gate reaches a *smaller* `du` than the old one
(`mol_base_handoff` 2.56 against 3.75; `mol_sec_ion` 2.23 against 4.68) and runs
faster (5 m 35 s against 6 m 59 s; 10 m 12 s against 11 m 39 s, single-threaded),
so the lighter base does not cost the convention physical time -- it gains some.

### 10.4 What moved, and what did not

`make check` in the tree, after the change and with the goldens deliberately
**not** refreshed:

| case | result |
|---|---|
| `wasp_full`, `wasp_he23off` | **PASS**, data identical |
| `lower_profile` | **PASS**, data identical |
| the six `mol_*` | **FAIL** -- moved, as intended |

Byte-identity was also checked directly, by building the same source snapshot
with and without the change and running the three cases with both binaries:
`Hydro_ioniz.txt` and `Ion_species.txt` identical in all six comparisons.

The six molecular cases move as follows (row count unchanged, 504, so the
goldens can be re-snapshotted rather than restructured; no NaN in any of them):

| case | `rho_bc` old -> new | `r_front` old -> new | peak `n(H3+)` old -> new | `N(H2)` new/old | `log10 Mdot` old -> new |
|---|---|---|---|---|---|
| `mol_base_handoff` | 1.220e14 -> 1.427e13 | 1.1546 -> 1.0741 | 6.68e4 -> 4.40e4 | 0.094 | 10.56 -> 10.42 |
| `mol_metals` | 1.234e14 -> 1.441e13 | 1.1562 -> 1.0748 | 3.89e1 -> 4.66e1 | 0.095 | 10.57 -> 10.41 |
| `mol_lyman_werner` | 1.220e14 -> 1.427e13 | 1.1243 -> 1.0350 | 6.67e4 -> 4.42e4 | 0.067 | 10.55 -> 10.51 |
| `mol_diffusion` | 1.220e14 -> 1.427e13 | 1.1546 -> 1.0741 | 6.68e4 -> 4.40e4 | 0.094 | 10.56 -> 10.42 |
| `mol_ir_bands` | 1.234e14 -> 1.441e13 | 1.1562 -> 1.0748 | 3.89e1 -> 4.66e1 | 0.095 | 10.57 -> 10.41 |
| `mol_sec_ion` | 1.220e14 -> 1.427e13 | 1.1401 -> 1.0622 | 2.89e5 -> 1.68e5 | 0.085 | 10.51 -> 10.39 |

The base is 8.6 times lighter, the H2 column falls to 7-10% of what it was, the
H2 front moves in by 0.08-0.09 `R_p`, peak `n(H3+)` falls by a third where H3+ is
the coolant, and `Mdot` falls by 0.04-0.16 dex. All of that is the direct
consequence of putting the base nine times higher; none of it is a change of
physics. **The goldens were not refreshed** -- that is left to be done once in
one pass.

---

## 11. Two smaller findings, not acted on

* **`EXHALE_setup.out` does not echo `Shapiro filter` or `Base velocity`.**
  `write_setup_report.f90` reports `Low-Mach damping`, `Hydrostatic base` and
  `Base ghost temperature` in the human-readable report, and carries all of them
  plus `shapiro_eps`, `shapiro_every` and `base_v_massflux` in
  `write_parse_dump` (`parse_dump.txt`), but the two keys are absent from the
  setup report itself. `input_read` does echo them to stdout, so they are not
  invisible; but the report is what the project documentation points at for
  "did my key take effect?".
* **The base ghost velocity reported in the output files is not a physical
  quantity in a JFNK-finished run.** It is `eps^2/(4|v_1|)` from the softplus
  valve that `EXHALE_main.f90` switches on at the hand-off (`valve_eps = 1e-4`),
  which is why it reads `+0.94 cm/s` on a base whose cell-1 velocity is
  `-249 cm/s`, and `0.000` in a marching-only run of the same state. Anything
  that reads the base flow direction from `v(0)` will read that number.

Neither was changed: `src/` and `docs/` outside this file were left alone for a
concurrent worker.

---

## 12. Reproduction

Scratch root `S = .../scratchpad/p44/`, binary `$S/EXHALE.x` (copy described at
the top). All runs restart from
`.../scratchpad/ladder/runs/heh0p0793c_r1/output/` copied to
`output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt`, with the reference
`input.inp` and `base.inp` and one key appended.

```
$S/runs/A_base                  (control, no key)
$S/runs/A_base2                 (control, second restart cycle)
$S/runs/B_shapiro               Shapiro filter: 0.02 4
$S/runs/C_massflux              Base velocity: massflux
$S/runs/D_ghostT                Base ghost temperature: continuous
$S/runs/E_lowmach               Low-Mach damping: 2.0e-2
$S/runs/E_lowmach2              (E, second restart cycle)
$S/runs/F_hydrobase             Hydrostatic base: True
$S/runs/G_mf_lm                 C and E together
$S/runs/A2_march1900            control, EXHALE_MAXSTEPS=1900 (marching only)
$S/runs/B2_shapiro_march1900    Shapiro, EXHALE_MAXSTEPS=1900 (marching only)
$S/resid_*/                     EXHALE_RESIDUAL=1 dumps of the states above
```

Each run: `OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=20000 $S/EXHALE.x`, 70-75 s except
D (7 min) and F (10 min), which never converged.

Section 7 (cold starts and base refinement), same binary, `Load IC? False`,
`du_th [PLM,WENO3]: 0.5 1.0e-3`, `Solver: Newton 2.0e-2`,
`OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=200000`, each with the snapshot archiver
running:

```
$S/cs/CS_A      cold start, default grid                      1 h 41 m, no info=0
$S/cs/CS_E      cold start + Low-Mach damping: 2.0e-2         1 h 41 m, no info=0
$S/cs/CS_A2x    cold start, Base grid 1.0e-4 100 / cells 550    47 m, info=0
$S/cs/CS_A4x    cold start, Base grid 5.0e-5 200 / cells 650  1 h 37 m, info=0
$S/cs/CS_*_r1, $S/cs/CS_*_r2   chained restarts of each (Load IC, Newton 100.0)
$S/birth/s0, s1 ... s1000      bounded cold runs, EXHALE_DUMP_IC=1 and MAXSTEPS
```

The step-resolved cold-start series of section 7.1 is read from the snapshot
archive of the pre-existing `ladder/runs/heh0p0793` (the same configuration,
198638 steps) and reproduced by `$S/cs/CS_A`; it is exported to
`$S/birth_series.txt`. Analysis scripts: `$S/anal.py` (base table and
alternating amplitudes), `$S/front.py` (H2 front and H3+), `$S/summary.py` (the
section 4 table), `$S/oe.py` + `$S/scan.py` (the odd-even purity series),
`$S/grid.py` and `$S/grid2.py` (the refinement table), `$S/make_fig.py` (the
figure). The predictions of section 7.2 are stored verbatim in
`$S/predictions.txt`, written at 2026-09-02 22:33 KST. The Figure 7
digitization of section 9.1 is calibrated on the axis label rows located in the
600 dpi render (left axis 464/703/942/1181/1419/1658 px for 6000-1000 K, x box
edges 1325/3747 px for `r/R_p` = 1/10) and validated against three numbers the
paper states in words (4384 vs "about 4400 K" at `r ~ 3`, 5538 vs "about
5540 K" at the upper boundary, 1149 vs the stated `T0` = 1140 K at the first
plotted point).

Section 8 (energy budget, deeper base, `T0` provenance), same binary:

```
$S/cs/CS_Adeep       cold start with Log10 lower boundary number density: 14.50
                     (base 2.844e-5 bar), default grid, 1 h 18 m, no info=0
$S/cs/CS_Adeep_r1, _r2   chained restarts of it, both info=0, identical to 4 digits
$S/resid_g1x, g2x, g4x   EXHALE_RESIDUAL=1 dumps of the three chain-converged grids
$S/resid_g1xa, g2xa, g4xa, g1xb   the same for the other restart cycle of each
```

Section 9 (published `T0`, base depth, closure and `T0` variants), same binary:

```
$S/cs/CS_P1ubar        cold start, Base BC: pressure 1.0 (n0 = 1.1119e13)
                       26 m 42 s, info=0 straight from cold
$S/cs/CS_P1ubar_r1,_r2 chained restarts, both info=0, identical to 5 digits
$S/cs/CS_Psame         Base BC: pressure 8.994, 1000-step equivalence check
$S/cs/CS_T1200         cold start, Equilibrium temperature 1200.0 + T_base 1200.0
                       1 h 28 m, no info=0
$S/cs/CS_T1200_r1,_r2  chained, both info=0
$S/cs/CS_GhostT        cold start, Base ghost temperature: continuous
                       1 h 30 m, budget exhausted at du = 3.1e-2
$S/cs/CS_GhostT_r1     chained, info=0 (base at 1501 K, front at 1.0009)
$S/cs/CS_GhostT_r2     chained again, info=2 -- does not stay converged
$S/time_t1x, t2x       2000-step timing pairs for section 9.5
$S/fig/p9-09.png       Koskinen page 9 at 600 dpi; $S/k22_fig7_T.txt the
                       digitized solid curve (2253 points, r/Rp 1.331-9.699)
```

Analysis for section 8: `$S/budget.py` (channel tables), `$S/hydro_terms.py`
(the cell-centred hydro terms, kept for the record -- they are dominated by the
artifact velocity and are not what section 8.2 uses). The Koskinen quotation in
section 8.5 is from the publisher PDF at
`references/Koskinen_2022_ApJ_929_52.pdf`, read with `pdftotext -layout`.

Figure: `docs/lower_atmosphere_figs/fig_p44_base_sawtooth.png`, twelve panels.
(a) `rho v r^2 / F_wind` over cells 1-12 for the four configurations asked for;
(b) the same with `Low-Mach damping` and the two marching-only states;
(c) the baseline's cell-centred product against the face flux reconstructed from
the steady mass residual; (d) the cold start, odd-even purity of `v` and `du`
against marching step, with the window in which the mode appears and the
secondary-ionization flip marked; (e) the three base grids on a common physical
radial axis; (f) `|v(1)|` against `dr` on log-log, with `dr` and `dr^2`
reference slopes; (g) the cell-1 energy terms against `dr`, showing that no
source term scales as `1/dr` except the steady residual itself; (h) `T(r)` on a
common physical axis for three 1x restart cycles and the 2x and 4x grids, with
the pinned `T0`; (i) the shallow and the deeper base compared over cells 1-12; (j) the three
base depths over cells 1-12; (k) the thermal gap and the base-face flux error
side by side for the baseline, the `T0 = 1200` run and the 1 microbar base,
showing that they move independently; (l) the digitized Koskinen Figure 7 and
the two EXHALE bases on a common `(r - r_base)/r_base` axis.

---

## 13. Scope

Checked: the eleven runs listed above, all from the same converged reference
state on one planet, and the steady residual dumps of six of them.

Added 2026-09-03: four cold starts to a 200000-step budget (baseline,
`Low-Mach damping`, and 2x and 4x base refinement), their chained restarts to
`info = 0`, eleven bounded cold runs for the first thousand steps, and the
snapshot series of the pre-existing 198638-step cold run. Later the same day: the
channel-by-channel heating and cooling budget of cells 1-5 on the three grids,
steady-residual dumps of seven converged states, a fifth cold start with the base
0.5 dex deeper and its two chained restarts, and the Koskinen et al. (2022)
lower-boundary statement read from the publisher PDF.

Not checked: a cold start under `Shapiro filter`, `Base velocity: massflux`,
`Base ghost temperature: continuous` or `Hydrostatic base` (only the baseline and
`Low-Mach damping` were run from scratch); any other planet; `Base BC: pressure`;
the transit observables; whether cell 1 responds to a Shapiro filter applied for
many more than 2000 marching steps; whether the 0.05 dex `Mdot` offset of the
Shapiro run persists over further restart cycles; an 8x base grid, which would
test whether `v(1)` stays first order below `dr = 5e-5 R_p`; and the term-by-term
budget of cell 1 that section 7.3 says the non-converging base temperature
needs.

Added later on 2026-09-03: the Koskinen Figure 7 digitization, four further
cold starts (1 microbar base, matched-pressure equivalence check, `T0` = 1200 K,
continuous ghost) and their chained restarts, and six alternating 2000-step
timing runs.

Not checked in section 9: whether the 1 microbar base behaves the same way on a
refined grid, or on any other planet; whether its thinner molecular layer changes
the transit observables; a ghost that pins `p` and lets `rho` and `T` float
together, which is the closure the over-specification of section 8.3 calls for
and which does not exist in the code; base depths between 1 and 9 microbar; and
whether the 100-450 K temperature deficit against Koskinen through the molecular
layer (section 9.1) is the H3+ cooling, which would need its own study. The
`T0` = 1200 K run changes the code's whole temperature normalization, not only
the ghost, so it bounds rather than isolates the effect of the pin.

Not checked in section 8: the split of the energy flux divergence itself into
its enthalpy and kinetic parts (the cell-centred estimate is dominated by the
artifact velocity, so only the code's own residual was used, and that gives the
net); the temperature Koskinen's lower and middle atmosphere model actually has
at 9 microbar, which is the number that would settle what the base should be
pinned to (their Figure 7 was not digitized); whether the deeper base behaves the
same way on a refined grid (only the default grid was run deep); any base depth
between 8.99e-6 and 2.844e-5 bar; and `Base BC: pressure`, which anchors the
base by pressure instead of density and is the natural third test of the
over-specification of section 8.3.
