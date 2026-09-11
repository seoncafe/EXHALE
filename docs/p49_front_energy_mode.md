# P49. What blocks the JFNK in the molecular configuration: it is not the H2 front's
# chemistry, it is the positivity guard of the reconstruction

**Judgment first, and it contradicts both hypotheses the item was opened to
test. The energy row at the front is smooth in the temperature, its
finite-difference directional derivative converges over six decades of step
size, and a two-sided scan of that row is exactly antisymmetric down to
`h = 1e-9` -- so a descent direction along it EXISTS. What has no descent
direction is the state as a whole, because `F(Y)` is genuinely DISCONTINUOUS
at the iterate: the WENO3 reconstruction at the face of cell 273
(`r = 1.2324 R_p`) returns a left density of `2.50e-15` against cell averages
of `1e-4` -- eleven orders below, and still positive by a hair -- so the
iterate sits exactly on the switching surface of `positivity_limited_faces`
(`src/modules/states/Reconstruction.f90`). Raising the density of cell 272 by
one part in `1e9` tips that face state through zero, the guard replaces the
face pair by the first-order cell averages, the face density jumps from
`2.50e-15` to `4.17e-6`, and the merit jumps from `202.5` to `945.1` -- by the
SAME factor 4.67 at every step size from `1e-2` to `1e-9`, which is what a
step discontinuity looks like and what a kink does not. Every direction the
solver can build crosses that surface, so every one of them ascends: the
Newton/PTC step ascends linearly in `lam` over six decades and all ten damped
Gauss-Newton steps ascend by factors 4.7 to 27.8. The section 126 escape is
not defeated by an indefinite Jacobian here; it is defeated by a residual that
is not continuous. Separately, and smaller: the residual's cooling is
evaluated from the composition handed IN, not the one the sweep returns, so
the temperature derivative the Jacobian carries at the front cell has an
effective exponent of 0.26 where the self-consistent one is about 2.0.**

Everything below is measured on this tree, reproduced from the archived state
of the P43 measurement, and reproduces its printed digits exactly.

## 1. The state, and that it is the same one

Source: the P43 working copy `.../scratchpad/p42/src` (second-order carrier
advection and the automatic Lyman-Werner output; not yet in the repository
tree), copied to `.../scratchpad/frontmode/src` and built with
`PATH=/usr/bin:$PATH make`. Case: the hot Uranus of `p42/E1d` -- 500 cells,
`Molecular chemistry: True`, `Molecular carrier transport: True`,
`Stellar LW flux: 343.0`, `Secondary_ionization: Immediate`,
`Solver: Newton 4.0e-2`, `Load IC? True` from the 150,000-step marched state.
`OMP_NUM_THREADS=8`.

The reproduction lands on the archived failure to every printed digit:

```
 (JFNK) it  62  ||R||=  2.782E-02  ||Fs||2=  1.85E+02  lam= 1.00E+00  gm=  3  worst j= 255 r=  1.186 k=3
 (JFNK) it  65  ||R||=  2.775E-02  ||Fs||2=  2.00E+02  lam= 7.63E-06  gm=  4  worst j= 255 r=  1.186 k=3
 (JFNK) it  66  no descent along the Newton or the damped Gauss-Newton direction; ||grad merit||= 4.583E+11
 (JFNK) no descent direction exists for the banded model at this state -- aborting
 (JFNK) below r_esc [1:391] |R|: mass= 4.617E-03 mom/grav= 3.445E-04 energy= 2.375E-02  max cell j=259 r=  1.1952 k=3
 (JFNK) done info=2 ||R||=  2.775E-02  non-monotone accepts=46
```

The convergence measure of this source is the section 127/131 one, not the
signal-speed scale of section 133 (that change is in the repository tree and
not in this copy), so `||R||` is on the old scale throughout.

**How the blocking state is produced, which is not by marching.** The run
reaches the hand-off at step 2 and the FIRST steady solve converges,
`info = 0` at `||R|| = 5.028e-4`. What follows is the carrier-transport outer
pass:

```
 (EXHALE_main) steady-wind diffusion outer pass 1: 63 relaxation steps, omega = 0.500, composition drift =  5.53E-01
    carrier steady residual max =  8.65E-03 at cell  273 r=  1.232 carrier H2
```

and the SECOND solve is handed a state with `||R|| = 8.166e-2`, 160 times
worse, whose worst row is already the energy row at `r = 1.198`. Comparing the
two states the two solves were handed, over `1.19 - 1.28 R_p`:

| at `r = 1.195` | handed to solve 1 | handed to solve 2 |
|---|---:|---:|
| `rho` [g cm^-3] | 6.99e-13 | 3.53e-14 |
| `v` [cm s^-1] | 227 | 6663 |
| `p` [code] | 3.94e-2 | 2.63e-3 |
| `f(H2)` | 4.87e-2 | 3.18e-1 |
| `f(H2)` at `r = 1.2324` | 1.02e-5 | 2.39e-1 |

Both are SMOOTH profiles. The cell-to-cell oscillation the abort state carries
(section 2) is therefore built by the 66 iterations of the second solve, not
handed to it: the solve drives the density at `1.19 - 1.23 R_p` down by a
further factor 45 and the temperature from about 900 K to 65,000 K while
lowering `||R||` from 8.2e-2 to 2.8e-2. **The convergence measure rewards the
excavation.**

## 2. The energy row, term by term

At the abort state, cells around the blocking one. Code units; multiply
`heat`, `cool` and the row entries by `q0 = 1.4090e-3 erg cm^-3 s^-1` for cgs
(the factor is measured, section 3). `dFE` is the flux-divergence-minus-source
part, `R3 = dFE - (heat - cool)`. `D3 = |E_j|` is the merit's own energy scale.
Times are in seconds, `t_therm = E/|heat-cool|`, `t_flow = dr_j/|v|`.
`class` is the acceptance class of that cell's chemical solve (4 = non-root
under the relaxation amnesty), `res` its normalized reaction residual.

| j | r | T [K] | rho [g/cm3] | `R3` | `dFE` | `heat` | `cool` | `|R3|/D3` | `t_therm` | `t_flow` | class | `res` | `f(H2)` | `f(H3+)` | `x(HII)` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 254 | 1.1833 | 664 | 1.18e-13 | 4.476e-4 | 9.093e-4 | 4.651e-4 | 3.358e-6 | 1.75 | 6.21e3 | 7.65e2 | 1 | 9.1e-18 | 0.327 | 1.5e-7 | 6.8e-4 |
| 255 | 1.1856 | 2191 | 2.96e-14 | 1.402e-2 | 1.360e-2 | 1.512e-4 | 5.673e-4 | 65.9 | 5.71e3 | 2.43e2 | 1 | 3.5e-22 | 0.325 | 1.8e-6 | 7.1e-4 |
| 256 | 1.1880 | 6258 | 8.78e-15 | 9.843e-3 | -7.449e-3 | 4.690e-5 | 1.734e-2 | 54.8 | 1.16e2 | 1.54e2 | 1 | 8.8e-19 | 0.324 | 9.5e-5 | 4.1e-6 |
| 257 | 1.1903 | 13673 | 3.85e-15 | 8.077e-3 | -2.715e-3 | 2.173e-5 | 1.081e-2 | 47.2 | 1.77e2 | 2.71e2 | 1 | 1.8e-17 | 0.322 | 1.9e-4 | 1.4e-6 |
| 258 | 1.1927 | 28235 | 1.81e-15 | 5.655e-3 | -6.210e-4 | 1.088e-5 | 6.287e-3 | 34.1 | 2.95e2 | 3.03e2 | 1 | 5.2e-21 | 0.321 | 3.6e-4 | 1.5e-6 |
| **259** | **1.1952** | **65173** | **7.73e-16** | **1.105e-2** | **2.736e-4** | **2.381e-5** | **1.080e-2** | **67.6** | **1.70e2** | **2.69e2** | **4** | **2.4e-3** | **0.313** | **4.5e-3** | **1.3e-4** |
| 260 | 1.1976 | 57245 | 8.70e-16 | -4.148e-3 | -4.144e-3 | 4.727e-6 | 5.034e-8 | 25.7 | 3.86e5 | 4.06e2 | 4 | 5.6e-6 | 0.315 | 0 | 0 |
| 261 | 1.2001 | 53033 | 9.43e-16 | 1.351e-3 | 1.356e-3 | 5.146e-6 | 3.771e-8 | 8.37 | 3.53e5 | 4.93e2 | 4 | 5.2e-6 | 0.315 | 0 | 0 |
| 262 | 1.2026 | 62047 | 8.01e-16 | 1.392e-3 | 1.395e-3 | 4.345e-6 | 7.917e-7 | 8.70 | 5.03e5 | 3.70e2 | 4 | 8.0e-6 | 0.313 | 0 | 0 |
| 263 | 1.2052 | 66462 | 7.48e-16 | 1.190e-3 | 1.187e-3 | 4.075e-6 | 6.386e-6 | 7.45 | 7.72e5 | 3.33e2 | 4 | 1.3e-5 | 0.312 | 0 | 0 |
| 264 | 1.2078 | 69022 | 7.19e-16 | 9.725e-4 | 9.659e-4 | 3.965e-6 | 1.058e-5 | 6.10 | 2.69e5 | 3.03e2 | 4 | 1.7e-5 | 0.310 | 0 | 0 |

**Which term carries the row.** At 254 the row is heating (`heat` is 99% of
the balance); at 255 it is the flux divergence; from 256 to 259 it is the
COOLING -- `cool` is 176%, 134%, 111% and 97.7% of `|R3|` there, i.e. a
radiative loss that nothing in the row balances; from 260 outward `cool`
collapses by five orders and the row is the flux divergence again.

**Which channel that cooling is.** From `Cooling_breakdown.txt` written at
this state (cgs):

| r | T [K] | `cool` total | largest channel | second |
|---:|---:|---:|---|---|
| 1.1833 | 664 | 4.731e-9 | H3+ IR 4.275e-9 (90%) | reco 3.58e-10 |
| 1.1856 | 2191 | 7.994e-7 | H3+ IR 7.993e-7 (100%) | reco 3.13e-11 |
| 1.1880 | 6258 | 2.443e-5 | H3+ IR 2.443e-5 (100%) | Lya 2.26e-14 |
| 1.1927 | 28235 | 8.858e-6 | H3+ IR 8.851e-6 (99.9%) | Lya 3.64e-9 |
| **1.1952** | **65173** | **1.699e-5** | **H3+ IR 1.585e-5 (93%)** | coex HeI 8.14e-7 |
| 1.1976 | 57245 | 7.090e-11 | coex HeI 3.05e-11 | Lya 2.21e-11 |
| 1.2052 | 66462 | 8.994e-9 | coex HeI 4.70e-9 | coio 2.66e-9 |

So the row that blocks the solve is 93% **H3+ infrared cooling at 65,000 K**,
computed from an `f(H3+) = 4.5e-3` that the chemical solve did not certify as
a root (class 4, reaction residual 2.4e-3). H3+ does not survive at 65,000 K;
the coolant is an artefact of a non-root acceptance, and it is the largest
single term in the residual the solver cannot reduce.

**Timescales.** At 256-259 `t_therm/t_flow` is 0.75, 0.65, 0.97, 0.63 -- the
cells sit exactly at the radiative/advective transition, which is where a
steady solution is hardest to pin. From 260 outward the ratio is 715 to 2319,
because the class-4 acceptance there sets `n(H3+) = x(HII) = 0` exactly and
removes every coolant: those cells are effectively adiabatic.

## 3. How the residual's temperature derivative is built, and what it misses

Scanning the energy unknown of cell 259 alone, in 61 steps of `1e-3` relative
(`eps` is the relative change of `E_259`; the composition seed is the iterate's
and is held fixed, exactly as a Jacobian column or a line-search trial does):

| quantity | at `eps = -3.0e-2` | at `eps = +1.9e-2` | `d ln X / d ln T` |
|---|---:|---:|---:|
| `T` [K] | 63253.10 | 66385.39 | -- |
| `cool` | 1.0720e-2 | 1.0857e-2 | **0.263** |
| `heat` | 2.3424e-5 | 2.4046e-5 | 0.542 |
| returned `f(H3+)` | 4.2641e-3 | 4.6356e-3 | **1.728** |
| `R3` | 5.354e-3 | 1.470e-2 | -- |

The two boldface exponents cannot both belong to the same rate, because
`h3p_cooling_rate` is exactly proportional to `n(H3+)`
(`util_ion_eq.f90:1319-1332`). They do not: within one `ioniz_eq` call the
molecular densities the cooling uses are built at the top of the routine from
the composition handed IN (`ionization_equilibrium.f90:570-573`,
`nmol_eq(:,3) = f_sp_io(:,isp_H3p)*n_in_dim`), `eval_cool` is called at line
779 with them, and the per-cell chemical sweep that changes the composition
runs afterwards and writes `f_sp_io` back at lines 1918-1921. There is one
`eval_cool` call per `ioniz_eq` and no inner self-consistency loop
(`grep "call eval_cool"`, `grep "call ioniz_eq"`).

So the residual's cooling responds to temperature only through the emissivity
at fixed densities -- exponent 0.263 -- while the self-consistent response
adds the density's own 1.728, about 2.0 in total. **The finite-difference
Jacobian, banded or matrix-free, sees roughly one seventh of the cooling's
true temperature sensitivity at this cell.**

The size of the lag at the iterate itself is small everywhere but there.
Recomputing `eval_cool` from the RETURNED composition (what
`write_cool_breakdown_eq` does) and dividing by the `cool` the residual used
gives the same constant `q0 = 1.4090e-3` at cells 254, 255, 256, 258, 260,
261 and 263 to five digits, and `1.5727e-3` at cell 259: **the lag is under
0.1% at every cell except the blocking one, where it is 11.6%.**

## 4. Why no direction descends

### 4.1 The two hypotheses, measured

**(B) "the finite-difference Jacobian is meaningless / the chemistry is
exponentially stiff in T": REFUTED along the blocking coordinate.** The
directional derivative `J*v` with `v` the relative perturbation of the energy
unknown of cell 259 converges cleanly:

| `h` | `(Jv)` at that row | `||Jv||` |
|---:|---:|---:|
| 1e-4 | 1.910992e-1 | 2.625400e-1 |
| 1e-5 | 1.910957e-1 | 2.625360e-1 |
| 1e-6 | 1.910954e-1 | 2.625357e-1 |
| 1e-7 | 1.910953e-1 | 2.625473e-1 |
| 1e-8 | 1.910953e-1 | 2.637013e-1 |
| 1e-9 | 1.910953e-1 | 3.609337e-1 |
| 1e-10 | 1.910954e-1 | 2.490732e+0 |

Six digits stable over four decades, round-off taking over below `h = 1e-8` in
the norm only. And the effective exponent the residual carries is 0.26 (section
3), not an exponential.

**(A) "the chemistry root selection is non-smooth in T": PRESENT, but it is
not what blocks the solve.** The T scan does show the accepted state jumping
between two branches at isolated points -- at `eps = 2.0e-2`, `2.1e-2` and
`2.9e-2` the solve returns `x(HII) = 0` exactly, `f(H3+) = 1.21e-4` instead of
`4.60e-3`, reaction residual `1.37e-5` instead of `3.2e-3`, and `T` 450 K below
the trend -- and returns to the first branch in between, which is seed- and
path-dependent selection rather than a threshold. **But the energy row does
not notice**: `cool` reads 1.0859e-2 on the jumped points against 1.0857e-2 on
the trend, because the cooling uses the seed's molecular densities (section 3)
and the seed does not jump. And a two-sided scan of the merit on the energy
rows of cells 255, 259 and 260 is smooth and antisymmetric:

| cell | `h` | `df2/f2` at `+h` | at `-h` | cells with class >= 3, `+h` / `-h` |
|---:|---:|---:|---:|---:|
| 259 | 1e-2 | +2.448183e-2 | -1.317689e-2 | 16 / 16 |
| 259 | 1e-4 | +1.898731e-4 | -1.887398e-4 | 16 / 16 |
| 259 | 1e-6 | +1.893122e-6 | -1.893009e-6 | 16 / 16 |
| 259 | 1e-9 | +1.893066e-9 | -1.893066e-9 | 16 / 16 |
| 255 | 1e-9 | +2.453093e-10 | -2.453090e-10 | 16 / 16 |
| 260 | 1e-9 | -1.946960e-9 | +1.946960e-9 | 16 / 16 |

Identical to nine digits with the WENO weights recomputed (`weno_mode = 0`)
and frozen at `Y` (`weno_mode = 2`), so the nonlinear reconstruction weights
are not involved either -- the same conclusion section 126 reached on the
He/H = 0.3 state. **Lowering `E` at cell 259 lowers the merit. A descent
direction along the front's energy rows exists.**

### 4.2 What actually blocks it: `F` is discontinuous

Measuring the merit gradient by forward differences of `f2` itself over cells
244-274, in the scaled coordinates, one component is the whole norm:

```
 (probe) ||grad f2||(window, scaled) =  7.42579E+08
 (probe) largest component: cell  272 row 1  g =  7.42579E+08
```

against 1e0 to 1e3 for every other component, including the front energy rows
(cell 259 row 3 reads 3.83e2). Stepping along that direction:

| step `t` | `df2/f2` along `+g` | along `-g` |
|---:|---:|---:|
| 1e-2 | +3.665786 | +8.281928e-1 |
| 1e-3 | +3.666659 | +1.651184e-1 |
| 1e-4 | +3.666745 | +1.897242e-2 |
| 1e-6 | +3.666755 | +2.195950e-4 |
| 1e-9 | +3.666755 | +8.412381e-7 |

The `+` column is **constant to seven digits at every step size**. That is not
a large gradient; it is a step discontinuity, and the iterate is on its
negative side. (The `-` column decays to zero, so the other side is
continuous.)

Dumping both sides of that switch at `h = 1e-9` on the density of cell 272:
every cell's temperature, class, heating, cooling and composition is
**identical to every printed digit** on the two sides, and so is every
reaction residual except cell 272's own, which moves in its last digits
(2.72e-19 / 1.72e-19 / 9.55e-19) as the solver reseeds. What
changes is the reconstruction:

| | `f2` | faces dropped to first order | `WL(rho)` at face 273 | `WR(rho)` at face 273 |
|---|---:|---:|---:|---:|
| at `Y` | 202.5166 | 2 | **2.50212e-15** | 9.05825e-7 |
| at `Y + h` | **945.0952** | **3** | **4.16758e-6** | 5.06403e-5 |
| at `Y - h` | 202.5167 | 2 | 2.28712e-14 | 9.05825e-7 |

and, with it, the mass-row residual of cell 273: `|R1|/D` goes from 9.501 to
918.3, a factor 97, with cell 274 following (20.90 to 52.64).

`positivity_limited_faces` (`src/modules/states/Reconstruction.f90:139-206`)
replaces a face pair by the piecewise-constant cell averages whenever the
reconstructed `rho` or `p` is not strictly positive. At this iterate the
WENO3 left density at face 273 is `2.50e-15` where the neighbouring faces
carry `1e-4` -- eleven orders below, positive by a hair. An arbitrarily small
increase of the neighbouring density tips it negative, the guard fires, and
the face density jumps by nine orders of magnitude. **The residual is
discontinuous at `Y`, and `Y` sits on the discontinuity.**

That accounts for every direction the solve tried, since all of them put their
largest scaled component at exactly that place (`max|dZ_newton| = 7.42e-3` at
cell 273, row 1):

| direction | `df2/f2` |
|---|---:|
| Newton/PTC, `lam = 1` | +1.332256e-1 |
| Newton/PTC, `lam = 1e-2` | +1.562406e-3 |
| Newton/PTC, `lam = 1e-6` | +6.992441e-7 |
| damped Gauss-Newton, `mu = 1e2` | +3.665784 |
| `mu = 1e1` | +3.661376 |
| `mu = 1e0` | +3.716514 |
| `mu = 1e-1` | +4.257652 |
| `mu = 1e-2` | +6.614426 |
| `mu = 1e-3, 1e-4, 1e-5, 1e-7` | not admissible (NaN) |
| `mu = 1e-6` | +2.682886e1 |

The Newton column falls linearly with `lam` over six decades and never changes
sign -- the section 126 signature of an uphill direction, reproduced here.
The Gauss-Newton column lands on `+3.6658` at `mu = 1e2`, the same plateau the
`+g` direction reaches: it is the post-jump value, not a gradient effect. So
the section 126 escape is intact and is simply being asked to descend a
function that is not continuous. `||grad merit|| = 4.583e11` in the banded
model is not a statement that the point is non-stationary; it is that model
reading a jump as a slope.

## 5. Verdict on the two hypotheses

* **(A) non-smooth chemistry root selection: present, measured, not the
  blocker.** 16 cells per sweep are accepted as class 3 or 4 (projection onto
  the element budget, or non-root under the amnesty), and the accepted branch
  at cell 259 does jump discontinuously with `T` at isolated points. It does
  not reach the energy row, because the cooling in the residual is built from
  the seed composition, which does not jump.
* **(B) genuine stiffness / a meaningless finite-difference derivative:
  refuted as stated.** `J*v` converges to six digits over four decades of
  step size, and the residual's temperature exponent at the blocking cell is
  0.26. The stiffness that IS there -- the composition's own exponent 1.73 --
  is absent from the residual by construction, so it can neither help nor
  hurt the Newton model.
* **The blocker is a third thing: a step discontinuity of `F` at the
  positivity guard of the reconstruction**, at a face whose WENO3 density has
  been driven to `1e-11` of the cell average by the same solve.

## 6. What was implemented, and what it did

Both changes live in the scratch working copy
`.../scratchpad/frontmode/src` (a copy of the P43 working copy
`.../scratchpad/p42/src`). Nothing in the repository tree was touched; the
tree application waits on sections 136 and 134 landing first.

### 6.1 (i) The positivity repair is continuous

`src/modules/states/Reconstruction.f90`. `positivity_limited_faces` no longer
swaps a face pair between two schemes. Each face state is scaled toward ITS
OWN cell average by the largest factor that keeps rho and p above their floor,

    W_face  <-  W_avg + theta ( W_rec - W_avg )
    theta    =  min over rho and p of  (q_avg - eps)/(q_avg - q_rec)

after the linear scaling limiter of Zhang and Shu (2010, J. Comput. Phys.,
229, 8918; doi:10.1016/j.jcp.2010.08.016; the publisher PDF is
`references/Zhang_2010JCP_229_8918.pdf`), their Section 2.2. Two new functions
carry the arithmetic, `positivity_scaling` (the two variables' minimum) and
`positive_variable_scaling` (one variable's share).

**It is not their limiter, and the differences are recorded here** rather than
left to the citation. Checked against the published version and measured in
`src/tests/physics_probe/positivity_limiter_scaling.f90`:

* They scale the **conserved** vector `w = (rho, m, E)` in two steps: density
  about the cell average, their Eqs. (2.4) and (2.5) with
  `theta_1 = min{(rho_avg - eps)/(rho_avg - rho_min), 1}`, then the whole
  vector, their Eqs. (2.10) and (2.11), with `theta_2` the root `t_eps` of the
  **quadratic** of their Eq. (2.12), `p[(1-t) w_avg + t w_rec] = eps`. This
  code scales the **primitive** vector `(rho, v, p)` in one step, so its
  `theta` is the root of a linear equation in `p`. The two numbers differ.
  Along their path `rho`, `m` and `E` are linear in `t`, so
  `p(t) = (gamma-1)(E - m^2/(2 rho))` is concave in `t` and never falls below
  the straight line between `p_avg` and `p_rec` that this code interpolates:
  at equal `eps` this `theta` is the smaller, that is the more restrictive,
  never the less safe. MEASURED at a Mach 60 state with `eps = 1e-3 p_avg`:
  `theta = 0.85629` here against `t_eps = 0.94910` from their Eq. (2.12).
* Their `theta_2` is one number per **cell**, the minimum over that cell's
  quadrature points, which is what keeps the cell average of the limited
  polynomial equal to `w_avg` (the conservativity property their Section 2.2
  requires). Here `theta` is one number per **face state**, so a cell whose two
  ends are scaled by different factors no longer reconstructs to its own
  average. MEASURED departure on a cell with both ends limited: `8.3e-2` of the
  cell average, against round-off when the cell's minimum `theta` is used for
  both ends.
* Their floor is absolute, `eps = 1e-13` in their computations and
  `eps = min_j {1e-13, rho_avg, p(w_avg)}` in their implementation flowchart.
  Here it is relative, one unit in the last place of the cell average.

What carries over is their Lemma 2.5: the limited value is a convex combination
of an admissible cell average with the reconstruction, so it stays admissible,
and no quadrature enters that argument. Their **Theorem 2.1 does not** carry
over. It concludes that the next *cell average* is admissible, and its
hypotheses are the `N`-point Legendre Gauss-Lobatto set of their Eq. (1.7) with
`2N - 3 >= k` (a finite-volume reconstruction to two face values supplies only
`N = 2`, which covers the `k = 1` of PLM and not the `k = 2` of WENO3; their
Remark 2.6 names finite-volume WENO as the open implementation case) and the
CFL condition of their Eq. (2.1), which nothing here imposes.

Pressure positivity of the scaled state is not at issue in this form: `p` is
itself one of the scaled variables, so the limited face carries the pressure
the scaling put there and its internal energy `p/(gamma-1)` has the same sign.
The nonlinear map that forces their quadratic has no counterpart here.

**The floor.** The switch this replaces tested `q > 0`, so its floor was zero,
and the instruction was to keep it. **Zero is not usable and the code does not
use it**: with `eps = 0` the scaling puts the face value at EXACTLY zero and
hands `Num_Fluxes` a `sqrt(gamma p/0)`. What is used instead is the smallest
floor that is not a chosen number -- one unit in the last place of the cell
average, `eps = epsilon(1.0d0)*q_avg`. No dimensional constant and no tuning
enter. The residual jump that survives at the crossing is `2.2e-16` of the
cell average against the factor `1.7e9` the hard switch made, so the
discontinuity is reduced to the level of the round-off the `theta = 1` branch
already carries.

**The floor is not delivered, and that is an open defect.** `eps` is one unit
in the last place of `q_avg`, so `q_avg + theta (q_rec - q_avg)` is a
cancellation of the same size as the floor and the result is decided by the
rounding of the product rather than by `theta`. MEASURED over 200000
reconstructed values that cross zero, using the production
`positive_variable_scaling` and the production update: the delivered value is
**exactly zero for 1013 of them** and never negative. Zero is the value this
floor exists to avoid, since a zero face density divides in `v = m/rho` and in
the sound speed of `Num_Fluxes.f90`. The same fact seen from the paper's side:
at Mach 60 a floor of one ulp of `p` is nineteen decades below `E`, so the
admissible set of their Eq. (2.6), which is stated on the conserved vector,
cannot be tested at that floor at all. Mapping the limited face to `(rho, m, E)`
and back returns `p = 0` where the limited pressure was `1.16e-10`.

The repair is to clamp the scaled variable to the floor it was solved for,
`q <- max(q_avg + theta (q_rec - q_avg), eps)`, which is the same number in
exact arithmetic and delivers `eps` in all 200000 of the sweep. It is not
applied here: it changes the bits of every run in which the limiter fires, so
it is a proposed change waiting on a golden refresh. Runs in which
`n_faces_positivity_limited` ends at zero are unaffected either way.

**Byte-identity.** A face with `theta = 1` is left untouched rather than
rewritten as `W_avg + 1*(W - W_avg)`, which is not bitwise `W`. A run in which
the guard never fires is therefore unchanged to the bit, and that is measured
(section 6.3).

The counter now counts face STATES rather than face pairs, because the limiter
scales the offending side only; `EXHALE_main.f90`'s end-of-run line says so.

### 6.2 (ii') A trial with an uncertified composition is refused

`src/modules/time_step/steady_newton.f90`, in `eval_residual`'s section 121
admissibility. A state whose equilibrium sweep accepted any cell as class 4
(the relaxation amnesty) or class 5 (the constrained continuation) is not a
state the residual describes, so it is refused where `admissible` is asked for.
Refusals are counted in the solver, separately from the `ioniz_eq` probe
ledger, and reported at the end of the solve:

```
 (JFNK) trial states refused for an uncertified composition: 164; worst cell count 13, worst reaction residual  1.455E-02
```

**Two refinements, both forced by measurement rather than chosen.**

1. **The test is a comparison, not a threshold.** The absolute form -- refuse
   any trial with an uncertified cell -- was written first and measured. On
   this state the ITERATE itself carries 11 such cells, so all 17 Jacobian
   colors and every Krylov direction were refused, the solve had no Newton
   model at all (`||grad merit|| = 0.000E+00`) and aborted at iteration 1 with
   `||R||` unchanged at 8.166e-2. A rule that rejects the neighbourhood of the
   point it stands on cannot be used to leave that point. The test is
   therefore `n_uncertified(trial) <= n_uncertified(iterate)`, the reference
   being refreshed from the iterate at the top of every outer iteration. At a
   fully certified iterate the two forms are the same rule.
2. **It applies to states that could be ADOPTED, not to derivative samples.**
   With the comparison in place the solve ran 18 iterations and then aborted
   again with 13 of 17 Jacobian colors zeroed and the GMRES cycle truncated:
   the refusal was now falling on the finite-difference probes. A Jacobian
   color and a `J*v` sample are points where the residual is read, never
   states anyone keeps, and a model built from an imperfect derivative is
   still a model. `eval_residual` therefore takes `may_be_adopted` (default
   `.true.`), and `build_banded_jac_full` and `jv_product` pass `.false.`

The log for both intermediate forms is kept: `E1d_fix_absolute_gate.log` and
`E1d_fix_relative_gate.log`.

### 6.3 What the two changes did

**The excavation is gone, and the blocking mode with it.** The second steady
solve of the reproduction case, same state, same source, three builds:

| | p42 as it was | + (i) + (ii') |
|---|---|---|
| outer iterations | 66 | 47 |
| `\|\|R\|\|` at exit | 2.775e-2 | **1.872e-2** |
| below `r_esc`, mass row | 4.617e-3 | **9.127e-3** |
| below `r_esc`, energy row | 2.375e-2 | **1.489e-2** |
| worst cell at exit | j=259, `r = 1.1952`, energy | j=266, `r = 1.2130`, energy |
| `T` at `r = 1.195` in the worst state | **65,173 K** | no such cell |
| `rho` at `r = 1.195` | 7.7e-16 g cm^-3 (45x below the state handed in) | not evacuated |
| trials refused for chemistry | -- | 164 (worst 13 cells, worst residual 1.5e-2) |
| `info` | 2 | 2 |

The 65,000 K, thousand-fold-evacuated layer the old solve built no longer
appears: the solve now works on the front itself, at 1.21-1.24 R_p, and its
below-escape energy residual falls by a factor 1.6 rather than being held up
by an H3+ population that cannot exist.

**Given the run's own retries, the hand-off then converges.** The reproduction
case pins `Stall [tol,N]: 1.0e-14 2000000`, so the retry of section 126(d)
would wait two million marching steps and never fire; re-run with that hold at
20,000 (`E1d_fix2`, the only input change, and it is a run parameter, not
physics), the three attempts are:

| attempt | at step | `\|\|R\|\|` in -> out | `info` | how it ended |
|---:|---:|---|---:|---|
| 1 | 2 | 8.166e-2 -> 1.872e-2 | 2 | no descent, 47 iterations |
| 2 | 44,120 | 5.393e-2 -> 2.766e-3 | 1 | iteration limit, no abort |
| 3 | 64,120 | 4.678e-3 -> **4.711e-4** | **0** | converged in 3 iterations |

**"no descent direction exists" is printed once in the whole run**, on the
state the carrier pass produces at marching step 2, and never again. At
attempt 3 the accepted solution's below-escape rows are mass 1.978e-4,
mom/grav 1.268e-5, energy 2.964e-4, and its worst cell is the base, j=1 --
not the front. On this source that is the residual gate met; the flux gate does
not exist here (section 8).

**The temperature hole is gone.** In the state the run finally writes the
front region is smooth and monotone -- `T` 1477 K at 1.152 R_p rising to
2455 K at 1.298, `rho` falling smoothly from 1.03e11 to 5.42e9 mH cm^-3, `v`
from 832 to 1.15e4 cm s^-1 -- against the 129-69,000 K cell-to-cell swing and
the 4.6e8 mH cm^-3 evacuation the blocked solve carried. The maximum
temperature anywhere in the wind is 6469 K at 1.955 R_p, an ordinary
photoionized peak.

**The atomic cases move, and by how much is measured.** The guard fires ~1200
times in the WASP-121b cases, so they cannot be byte-identical; the coordinator
asked for the count first, and it is:

| case | guard firings, p42 base | with (i) | trials refused by (ii') |
|---|---:|---:|---:|
| `wasp_full` | 1195 | 1300 | 0 |
| `wasp_he23off` | 1208 | 1302 | 0 |
| `wasp_full_newton` | 1195 | 1300 | 0 |
| `wasp_he23off_newton` | 1208 | 1302 | 0 |
| `lower_profile` | 0 | 0 | 0 |

(The two counts are not directly comparable: the old one counted face pairs,
the new one face states.) `refused = 0` in all of them, so **the whole move is
(i); (ii') changes nothing outside the molecular configuration.**

How far they moved, against the same source without the change:

| case | `log10 Mdot` base -> fix | `info`, `\|\|R\|\|` base -> fix | max relative change, `r > 1.2 R_p` |
|---|---|---|---|
| `wasp_full` | 13.26 -> **13.26** | (stops on `du`) | rho 4.6e-5, v 4.8e-5, T 5.9e-6 |
| `wasp_he23off` | 13.26 -> **13.26** | (stops on `du`) | rho 2.2e-5, v 2.2e-5, T 2.5e-6 |
| `wasp_full_newton` | 13.21 -> **13.21** | 0, 9.886e-4 -> 0, **8.039e-4** | rho 2.9e-3, v 3.0e-3, T 8.9e-4 |
| `wasp_he23off_newton` | 13.22 -> **13.22** | 0, 8.525e-4 -> 0, **6.985e-4** | rho 2.3e-3, v 2.3e-3, T 6.6e-4 |
| `lower_profile` | -- | -- | **byte-identical** |

The rate is unchanged to the printed digit in all four, both Newton solves
still return `info = 0` and reach a SMALLER residual than before, and the one
case in which the guard never fires is byte-identical -- which is the
byte-identity property the change was required to have. The largest single
relative change anywhere in `wasp_full` is 2.4e-1 in the velocity of the base
cell at `r = 1.00235`, where `v` is -1.9 cm s^-1: that is the collocated
odd-even mode of `docs/p44_base_sawtooth.md`, not a change in the wind.

### 6.4 The lagged cooling is not a defect (prescription (iii-a), judged)

Within one `ioniz_eq` the cooling is evaluated from the composition handed in
and the composition is solved afterwards (section 3). **At the fixed point the
two are the same object**: the steady solve converges when `f_sp` is
reproduced by the sweep it seeds, and there the cooling is evaluated on the
composition it returns. The lag is therefore a property of the ITERATION, not
of the solution, and re-evaluating `eval_cool` after the sweep would change no
converged answer -- only the path to it. It is not corrected. What it does cost
is the Newton model: the residual's temperature exponent at the front is 0.26
where the self-consistent one is about 2.0 (section 3), so the Jacobian
understates the cooling's response by about a factor 7. That is prescription
(iii-b) -- add the analytic `d ln cool/d ln T` of the dominant channel to the
banded Jacobian's energy diagonal -- and it is NOT done here: it is held until
(i) and (ii') are shown to be insufficient.

### 6.5 What was tried and rejected

**A log-temperature energy row is not indicated.** The two-sided scan of
section 4.1 shows the energy rows smooth and well scaled in `E` itself, with
derivatives of order 2 per unit scaled step; the obstruction was in a mass row
two cells away. Not implemented.

## 7. Reproduction

* **Both changes are now in the tree** (`docs/Update_EXHALE_stage1.md` section 138),
  where they were re-measured against the tree's own section 133 two gates and
  the refreshed block-F goldens. The numbers in sections 6.3 and 6.5 below were
  taken on the P43 working copy, which is where the blocked state exists; the
  tree numbers are in section 138.
* Investigation source: `.../scratchpad/frontmode/src`, built with
  `PATH=/usr/bin:$PATH make`; binary kept as `EXHALE_fix.x`. The unmodified
  P43 copy is built beside it as `EXHALE_base.x` from
  `.../scratchpad/frontmode/basebuild`.
* Diagnostic build used for sections 2 to 4, and NOT part of the deliverable:
  `.../scratchpad/frontmode/probe_build` (`front_energy_mode_report` in
  `steady_newton.f90`, per-cell `ieq_class_cell`/`ieq_res_cell` in
  `ionization_equilibrium.f90`), gated by `EXHALE_FRONT_PROBE`.
* Measurement runs: `E1d` (the anatomy), `E1dB`-`E1dF` (the direction,
  gradient, kink, switch and hand-off dumps), `E1d_fix` (both changes),
  `E1d_fix2` (the same with the retry hold shortened to 20,000 steps),
  `reg/` and `regbase/` (the case matrix on the two builds).
* Channel breakdowns of the abort state: `.../scratchpad/frontmode/snap1/`.

## 8. What was NOT measured

* **The three-layer table is NOT produced, and nothing is quoted for it.**
  The `info = 0` solve of attempt 3 is a converged wind, but the outer loop
  then took a SECOND carrier pass (400 relaxation steps, composition drift
  0.328, worst at `r = 2.119` in the wind rather than at the front), which
  threw the state to `||R|| = 3.619e-1`; the fourth solve aborted, and the
  state the run writes is that one. `Hydro_ioniz.txt` and the run's
  `log10 Mdot = 10.11` therefore do NOT describe the converged solution and
  are not quoted as one. Capturing it needs one more run that stops at the
  accepted hand-off -- either an output write at the `info = 0` return, or the
  carrier outer passes capped at one. **The carrier/wind alternation not
  converging is a separate open question**, and it is the one
  `docs/supersonic_molecular_base.md` section 13.7 left open ("the transported
  front is still propagating"); it is not what P49 was about and nothing here
  addresses it.
* **The two-gate acceptance cannot be evaluated on this source.** The flux
  gate of section 133 (`flux_spread_of_state`, `steady_gates_met`) is in the
  repository tree and NOT in the P43 working copy this work is built on, which
  still carries the section 127/131 residual measure alone.
* **Six of the nine regression cases cannot run on this source at all**, before
  and after the change: the P43 copy's `input_read.f90` still requires
  `Log10 lower boundary` (`req(...)`, line 135) where the tree has made it
  optional, and `mol_base_handoff`, `mol_metals`, `mol_lyman_werner`,
  `mol_diffusion`, `mol_ir_bands` and `mol_sec_ion` all carry that key in
  `base.inp` instead. They error identically on both builds. The matrix
  actually exercised is `wasp_full`, `wasp_he23off`, `lower_profile` and the
  two `*_newton` cases, compared against the SAME source without the change
  rather than against the tree goldens -- the tree goldens do not describe this
  source.
* `run_fcheck.sh` (the `-fcheck=bounds,do,mem` build) was not run.
* Whether the He/H = 1 reference moves. It was not re-run.
