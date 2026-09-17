# L25: the low-Mach base of the seven low-XUV atomic cases, measured

Item L25 of `docs/PLAN_20260913_lhs_stationary.md`, specified in
`docs/PLAN_20260916_rev3.md` section 5. Step 1 is written here; step 2 (the
Rieper velocity-jump factor on the Roe branch) is a separate worker's and is
appended below this section when it is done.

Every number is MEASURED on this tree on 2026-09-16 at one thread unless it is
marked READ. **No production source was changed by step 1.** What was added is
one test, `src/tests/low_mach_stress_energy/`, which is standalone and links
no production object.

## Step 1: the existing damping option

### 1.1 Verdict

**The nonnegative-dissipation claim of the `low_mach_dissipation.f90` header
is FALSE as the term is written. It is true only where the face coefficient
`eps4 g rho_f lambda_f r_edg^2` is constant.** The header's integration by
parts is correct about the boundary terms, which vanish exactly because the
stress is set to zero at the faces `j = 0` and `j = N`; what it drops is the
variable coefficient. The discrete kinetic-energy rate is the quadratic form
`E' = d^T W L d` with `d` the face differences of `v`, `W` the diagonal of the
face coefficients and `L` the second difference; the symmetric part of `W L`
is not negative semidefinite unless `W` is a multiple of the identity.
MEASURED (section 1.3):

- constant coefficient, uniform grid, N = 64 and N = 500: `lambda_min` of the
  symmetric part is at the rounding of the eigensolver, so the claim holds
  there;
- a coefficient falling four decades smoothly over 500 cells: `lambda_min` is
  48 times the eigensolver's rounding, so the form is already indefinite;
- **the term's OWN GATE closing**, i.e. the coefficient stepping to zero at one
  face, which is what `g = [max(0, 1 - M_f^2/M_th^2)]^2` does wherever the Mach
  number crosses `M_th`: `lambda_min/lambda_max = -6.9e-03`, and the negative
  mode sits on the cells astride the gate edge;
- **on the actual LHS 1140 b grids and states**: `lambda_min = -4.147e+02`
  (0.02 XUV) and `-6.958e+02` (0.03 XUV) against an eigensolver rounding of
  8.0e-02, and the negative mode lives at r = 4.0 to 5.1 R_p and 3.5 to
  4.3 R_p, which is exactly where the gate closes (the gate is open on faces
  1-388, r <= 4.98, and 1-376, r <= 4.23). The same `lambda_min` comes out of
  the block that no ghost value reaches, so the sign failure is interior and
  no boundary closure produces or removes it.

Because the energy pair conserves total energy exactly, a positive `E'` is a
term that takes heat OUT of the gas and puts it into the velocity field. The
mode that does this is at the gate edge, in the wind at 3.5-5 R_p, and not at
the base the item is about.

**The correctly signed form rev3 names is derived in section 1.2 and measured
in section 1.3: it is nonpositive on every grid tested, including the two
actual ones, and it is algebraically IDENTICAL to the present stress wherever
the coefficient is constant.** It is reported and NOT inserted into the
production code.

**On the frozen states the damping shrinks the alternating base velocity by
about a third and changes nothing else that is measurable** (section 1.4).
With `eps4 = 0.02`, `M_th = 1e-3`, 200 Newton iterations of one outer pass: the
2 dr amplitude of `v` over cells 2-8 ends at 0.664 (0.02 XUV) and 0.685
(0.03 XUV) of the undamped run, but at three times the value the loaded state
carried, so the solve GROWS the mode in both branches and the damping only
slows it; `Mdot` is the same to six significant digits (the gate closes at
r ~ 4-5 R_p and `Mdot` is read at 20.8 R_p); the base rows end above their
tolerances by five to six decades in both branches; neither state certifies
with the damping on or off; and the `ratio = -26.8` trust-region model failure
of L17 does not occur in either branch (1.007 off, 1.006 on at the first
radius), so its disappearance is not the damping's doing.

**A finding outside L25 blocks the acceptance question of the item as the
brief states it** (section 1.5): on today's tree neither `.L14` state
reproduces the residual it was certified or measured at on 2026-09-14/15. The
physical cells are loaded bit for bit; the two base GHOST cells are not, and
the base energy row at cell 2 is now its own flux divergence (`R = -9.078e-02`
against `-1.719e-09` in L17). The tree binary `EXHALE.x` gives the same
numbers as the private build, so this is the tree and not this build. The
`0.03` control therefore does not re-certify with the damping OFF either, and
"does the control re-certify with the damping on" cannot be answered today.

### 1.2 The discrete energy form, derived on the grid the code uses

**What the code's divergence and face averaging actually give.**
`contact_mode_dissipation_flux` fills, for the interior faces `j = 1..N-1`
only,

```
D_j = eps4 g_j rho_f,j lambda_f,j ( v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1} )
D_E,j = v_f,j D_j ,        D_0 = D_N = 0 ,
```

with `rho_f`, `lambda_f`, `v_f` the arithmetic face averages and `g_j` the
gate. `RK_rhs` adds `D` to the numerical flux of the face loop and forms the
cell update with the spherical weights

```
d(rho v)_j/dt |_D = - ( A_j D_j - A_{j-1} D_{j-1} ) / dV_j ,
A_j = r_edg(j)^2 ,    dV_j = ( r_edg(j)^3 - r_edg(j-1)^3 ) / 3 .
```

Under `Well balanced: True` the same `A_j D_j - A_{j-1} D_{j-1}` appears, in
both the PLM and the WENO3 branch, so the weights below are the weights of the
LHS 1140 b runs as well. The mass flux is not touched, so at frozen density
the discrete kinetic energy `sum_j dV_j (rho v)_j^2/(2 rho_j)` changes at

```
E' = sum_j dV_j v_j d(rho v)_j/dt |_D
   = - sum_{j=1}^{N} v_j ( A_j D_j - A_{j-1} D_{j-1} )
   = sum_{j=1}^{N-1} A_j D_j ( v_{j+1} - v_j ) .                     (E)
```

**(E) is exact and carries NO boundary term**, because `D_0 = D_N = 0`: on this
point the header is right. With `d_j = v_{j+1} - v_j` the undivided third
difference is the plain second difference of `d`,

```
v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1} = d_{j+1} - 2 d_j + d_{j-1} = (L d)_j ,
```

so with `w_j = A_j eps4 g_j rho_f,j lambda_f,j >= 0`

```
E' = sum_{j=1}^{N-1} w_j d_j (L d)_j = d^T W L d ,   W = diag(w) .      (F)
```

`L` is negative semidefinite, so for `W = w I` the form is
`-w sum_j (d_{j+1} - d_j)^2 <= 0` and the header's conclusion follows. For a
varying `W` the form is governed by the symmetric part `(W L + L W)/2`, which
is NOT a product of `W` with a definite matrix and is not sign definite: one
summation by parts turns (F) into `-sum_j q_j (B^T (W d))_j` with
`q_j = 2 v_j - v_{j+1} - v_{j-1}`, whereas a dissipative form needs
`-sum_j w_j q_j^2`. The two differ by the terms in which the difference
operator falls on `W`, i.e. by the coefficient's own variation, and the gate
makes that variation a step.

Neither the face averaging of `rho_f`, `lambda_f` nor the spherical areas can
repair this: they enter only through `w >= 0`, which is precisely the thing
that varies.

**The energy pair.** With `D_E,j = v_f,j D_j` the total-energy flux difference
telescopes over the physical cells and vanishes, so total energy is conserved
exactly and the internal-energy change is exactly `-E'`. Wherever `E' > 0` the
pair is therefore an anti-dissipative term that converts heat into kinetic
energy, not the reverse; the header's last sentence of its section 2 is the
claim that this cannot happen.

**The correctly signed replacement (rev3 `dv/dt = -M^-1 B^T K B v`).** Let

```
M = diag( rho_j dV_j )                  the volume (mass) matrix, M > 0
(B v)_j = q_j = 2 v_j - v_{j+1} - v_{j-1}   the cell second difference
K = diag( k_j ) >= 0 ,   k_1 = k_N = 0
```

and take the conservative momentum flux, at the interior faces only,

```
A_j D_j = k_j q_j - k_{j+1} q_{j+1} ,   j = 1..N-1 ,   D_0 = D_N = 0 ,  (G)
```

with the same energy pair `D_E,j = v_f,j D_j`. Substituting (G) into (E) and
summing by parts once gives, with no boundary term left because `k_1 = k_N = 0`,

```
E' = - sum_{j=2}^{N-1} k_j q_j^2 <= 0                                   (H)
```

for every velocity field and every nonuniform grid. In matrix form the update
is `dv/dt = -M^-1 B^T K B v` with `B^T = B`, so
`d(v^T M v/2)/dt = -(Bv)^T K (Bv) <= 0` at fixed `M`. No ghost value enters
(G) at all: the first interior face reads `q_2`, which needs `v_1, v_2, v_3`.

Two properties make (G) the natural replacement rather than a different
scheme:

1. **For constant `k` it IS the present stress.** Expanding (G) with
   `k_j = k_{j+1} = k` gives `k ( v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1} )`, the
   third difference of the header's equation (1). The two forms differ only
   through the coefficient's variation, which is exactly what breaks the sign.
   The order of accuracy, the 2 dr decay rate and the explicit-stability bound
   of the header's sections 4 and 6 carry over unchanged.
2. **It is a flux**, so momentum and, with the same pair, total energy stay
   conserved to round-off, and the stencil is still `j-2 .. j+2`.

The cell-centered coefficient that matches the present one is
`k_j = eps4 g(M_j^2) rho_j lambda_j A^c_j` with `A^c_j` the cell's own area
(`r_j^2`, or the mean of its two face areas) and `k_1 = k_N = 0`; the gate is
then read at the cell, which is where `M_j^2` is defined in the first place.

### 1.3 The measured quadratic forms and eigenvalues

`src/tests/low_mach_stress_energy/run.sh [Hydro_ioniz.txt]`. The driver
assembles the matrix `A` of `E' = -v^T A v` over the N physical cells,
symmetrizes it and calls LAPACK `dsyev`; nonnegative dissipation is
`lambda_min >= 0`. The reference on every row is `-100 eps lambda_max`, the
accuracy of the symmetric eigenvalue problem itself. `jst` is the present
stress, `bkb` the form (G). The ghost values the present stress reaches at the
first and last interior face are closed by zero gradient; the row
`jst_state_interior` repeats the measurement on the block `3..N-2`, which no
ghost value reaches.

| configuration | operator | lambda_min | lambda_max | lambda_min/lambda_max | verdict |
|---|---|---|---|---|---|
| uniform grid, constant coefficient, N = 64 | present | +1.150e-15 | 1.598e+01 | +7.2e-17 | dissipative |
| uniform grid, constant coefficient, N = 500 | present | -3.896e-16 | 1.600e+01 | -2.4e-17 | dissipative |
| uniform grid, coefficient falling 4 decades smoothly, N = 500 | present | **-1.521e-11** | 1.431e+01 | -1.1e-12 | **indefinite** (48 x the eigensolver rounding 3.2e-13) |
| uniform grid, the gate closed over the outer half, N = 200 | present | **-1.096e-01** | 1.599e+01 | **-6.9e-03** | **indefinite** |
| LHS 1140 b, 0.02 XUV state, N = 500 | present | **-4.147e+02** | 3.610e+12 | -1.1e-10 | **indefinite** (5.2e+03 x the rounding 8.0e-02) |
| LHS 1140 b, 0.02 XUV, interior block only | present | **-4.147e+02** | 3.441e+12 | -1.2e-10 | **indefinite**, no ghost involved |
| LHS 1140 b, 0.03 XUV state, N = 500 | present | **-6.958e+02** | 3.613e+12 | -1.9e-10 | **indefinite** (8.7e+03 x the rounding) |
| uniform grid, constant coefficient, N = 64 / 500 | (G) | -6.5e-16 / -1.3e-15 | 1.598e+01 / 1.600e+01 | -4e-17 / -8e-17 | dissipative |
| uniform grid, smooth 4-decade coefficient | (G) | -2.552e-16 | 1.409e+01 | -1.8e-17 | dissipative |
| uniform grid, gate closed over the outer half | (G) | -7.251e-16 | 1.599e+01 | -4.5e-17 | dissipative |
| LHS 1140 b, 0.02 XUV state | (G) | -2.106e-05 | 3.563e+12 | -5.9e-18 | dissipative |
| LHS 1140 b, 0.03 XUV state | (G) | -1.124e-04 | 3.566e+12 | -3.2e-17 | dissipative |

Where the negative mode of the present stress lives, and the same number
evaluated twice (once through the assembled matrix, once face by face from the
stress itself, so that the assembly is checked):

| configuration | E' on the negative mode, matrix | E' face by face | mode support | where the gate closes |
|---|---|---|---|---|
| gate closed over the outer half, N = 200 | +1.096e-01 | +1.096e-01 | cells 96-102 | face 100 |
| LHS 1140 b, 0.02 XUV | +4.147e+02 | +4.147e+02 | cells 372-390, r = 4.010-5.120 | face 388, r = 4.978 |
| LHS 1140 b, 0.03 XUV | +6.958e+02 | +6.958e+02 | cells 362-378, r = 3.528-4.342 | face 376, r = 4.227 |

On fields that are not the negative mode the stress does remove energy, as its
size suggests: five pseudo-random fields on the 0.02 state give `E'` between
-1.85e+13 and -2.56e+13, and the alternating base mode
`v_j = (-1)^j` in cells 1-6 gives -1.923e+13 for the present stress and
-1.814e+13 for (G). **The alternating base mode is damped by both forms**; the
indefiniteness is not a statement about that mode.

The gate on the two states (MEASURED from the states' own Mach numbers with
`M_th = 1e-3` and gamma = 5/3): open on faces 1-388 of 499 (r <= 4.978) at
0.02 XUV and 1-376 (r <= 4.227) at 0.03 XUV, one contiguous block in both, and
the face coefficient `w` runs over eight decades inside it
(2.20e+03 to 2.67e+11 at 0.02 XUV).

The grid the driver uses is the state's own cell centers with the face radii
reconstructed as their midpoints, since `r_edg` is not written to
`Hydro_ioniz.txt`; the sound speed is `sqrt(5/3 p/rho)` (both runs are
atomic). Neither approximation touches the sign question, which is a statement
about the variation of a nonnegative coefficient.

### 1.4 The damping on the frozen states

Fixture: read-only copies of `LHS1140b/models/.L14/x002_HeH2.13/output/` (the
stalled 0.02 transient, `certified=F cert_reason=no_stationary_claim`) and
`.L14/x003_HeH2.13/output/` (the 0.03 control, `certified=T` as its own run
wrote it) with their `input.inp`, placed as `output/*_IC.txt` in a scratch
directory, `Load IC? True`, `Restart intent: stationary`, `Well balanced:
True`, `Numerical flux: HLLC`, `Reconstruction scheme: PLM`,
`EXHALE_PTC_DTAU0=1.0` as the case's `REPRODUCE.md` records,
`EXHALE_OUTER_PASSES=1`, `EXHALE_JFNK_MAXIT=200`, `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`. The damping is switched on by the single added line
`Low-Mach damping: 0.02 1.0e-3` (the classical JST range, inside the explicit
bound `eps4 < 1/(16 CFL)` the header states). Binary `EXHALE_L25a.x`,
md5 `f23da8d84fafc2f0587b610add56e714`, built as
`make OBJDIR=build_L25a EXE=EXHALE_L25a.x` from the tree of 2026-09-16 and
deleted afterwards; `EXHALE.x` and `build/` were not touched.

**What the damping does, in one line: it takes a third off the alternating
velocity of cells 2-8, leaves `Mdot` unchanged to six digits, and leaves the
base rows, the certification verdict and the solver `info` where they were.**

#### The base rows at the end of the pass

Both states run the full 200 Newton iterations of the single outer pass and
both return `info = 1`. The worst cell of every hydrodynamic row is inside
cells 1-6 in every run, so the rows below ARE the base rows.

| | 0.02 transient, off | 0.02 transient, on | 0.03 control, off | 0.03 control, on |
|---|---|---|---|---|
| mass row, max (cell) | 2.924e-01 (1) | 4.916e-01 (1) | 5.055e-01 (1) | 5.061e-01 (1) |
| mass tolerance at that cell | 8.3e-07 | 4.4e-07 | 3.0e-07 | 3.0e-07 |
| momentum row, max (cell) | 1.896e-08 (2) ABOVE | 1.744e-08 (1) ABOVE | 3.488e-09 (1) within | 4.220e-09 (1) within |
| energy row, max (cell) | 7.458e-01 (1) | 4.033e-01 (2) | 2.023e-01 (3) | 2.025e-01 (3) |
| certification | NOT CERTIFIED | NOT CERTIFIED | NOT CERTIFIED | NOT CERTIFIED |
| solver `info` | 1 | 1 | 1 | 1 |

#### The alternating velocity of cells 1-6

`v` in cm/s, and the amplitude of the 2 dr component measured as the RMS of
the second difference of `v` divided by four (on a pure `(-1)^j a` mode that
is `a`).

| cell | 0.02: loaded | off | on | 0.03: loaded | off | on |
|---|---|---|---|---|---|---|
| 1 | -1.28048e-01 | -2.67210e+00 | -2.55482e+00 | -1.27875e-01 | -2.77708e+00 | -2.69295e+00 |
| 2 | -1.52989e-01 | +3.14554e-01 | +2.15860e-01 | -1.52670e-01 | +3.29632e-01 | +2.38182e-01 |
| 3 | +2.11780e-02 | -2.51379e-01 | -1.78430e-01 | +2.18352e-02 | -2.61326e-01 | -1.88961e-01 |
| 4 | -1.35775e-02 | +7.09577e-02 | +3.31774e-02 | -1.29416e-02 | +7.43774e-02 | +3.60439e-02 |
| 5 | +5.60928e-03 | -3.07235e-02 | -1.34321e-02 | +6.30735e-03 | -3.16149e-02 | -1.39667e-02 |
| 6 | -4.56185e-04 | +1.17967e-02 | +4.29429e-03 | +2.63183e-04 | +1.27054e-02 | +4.96099e-03 |
| 2 dr amplitude, cells 1-8 | 3.009e-02 | 3.766e-01 | 3.300e-01 | 3.013e-02 | 3.920e-01 | 3.504e-01 |
| 2 dr amplitude, cells 2-8 | 2.431e-02 | 1.114e-01 | **7.403e-02** | 2.435e-02 | 1.162e-01 | **7.955e-02** |
| on / off, cells 2-8 | | | **0.664** | | | **0.685** |
| on / off, cells 1-8 | | | 0.876 | | | 0.894 |

**The damping shrinks the mode and does not remove it, and the solve grows it
in both branches.** Over cells 2-8 the amplitude leaves the loaded state at
2.43e-02, and 200 iterations take it to 1.11e-01 with the damping off and
7.40e-02 with it on: the damped run ends at 0.66 of the undamped one and at
three times the value it started from. Over cells 1-8 the reduction is only
11-12 per cent, because cell 1 is set by the base boundary and not by the
stress.

#### Movement of the base state and of the mass-loss rate

Relative to the state as loaded; `Mdot` at cell N-20, the cell the
documentation evaluates it at.

| | 0.02, off | 0.02, on | 0.03, off | 0.03, on |
|---|---|---|---|---|
| `d rho / rho` at cell 1 | +3.67e-02 | +5.13e-02 | +5.30e-02 | +5.37e-02 |
| `d rho / rho`, worst of cells 2-6 | -2.00e-03 | -2.89e-03 | -2.85e-03 | -2.91e-03 |
| `d T` at cell 1 [K] | -17.53 | -24.17 | -24.91 | -25.20 |
| `d T`, worst of cells 2-6 [K] | +0.46 | +0.70 | +0.65 | +0.68 |
| `Mdot` [g s^-1] | 1.294097e+06 | 1.294095e+06 | 1.957191e+06 | 1.957191e+06 |
| `Mdot` as loaded [g s^-1] | 1.294101e+06 | | 1.957199e+06 | |
| `d Mdot / Mdot` | -3.1e-06 | -4.6e-06 | -4.1e-06 | -4.1e-06 |

**`Mdot` is the same number with and without the damping, to six significant
digits**, which is what the gate is designed to give: `M_th = 1e-3` closes the
term at r ~ 4-5 R_p, and cell 480 is at r = 20.8 R_p. The base cell moves
1.5 per cent more in density and 6.6 K more in temperature with the damping on
at 0.02 XUV, and by less than the difference between the two states at
0.03 XUV.

#### The Newton trajectory

| | 0.02, off | 0.02, on | 0.03, off | 0.03, on |
|---|---|---|---|---|
| `\|\|R\|\|` at iteration 1 | 1.7180e+00 | 1.6070e+00 | 1.7070e+00 | 1.7890e+00 |
| `\|\|R\|\|` minimum (iteration) | 4.8520e-01 (171) | 4.9160e-01 (85) | 5.0550e-01 (199) | 5.0610e-01 (194) |
| `\|\|R\|\|` at iteration 200 | 7.4580e-01 | **4.9180e-01** | 5.0650e-01 | 5.1010e-01 |

The minima are the same to two per cent in all four; on the 0.02 transient the
damped run ENDS at its minimum while the undamped one wanders back up by 54
per cent, and on the 0.03 control the two are indistinguishable. Neither
branch descends: both leave a state whose rows are five to six decades above
the certification, which on this tree is where the state arrives (section 1.5).

#### The trust-region model

`EXHALE_TRUST_REGION=1`, `EXHALE_PTC_DTAU0=1.0e8`, 20 iterations, the 0.02
transient, which is the measurement L17 made.

| | off | on |
|---|---|---|
| first radius | 3.873e-05 | 3.873e-05 |
| `pred` / `actual` / `ratio` at iteration 1 | 8.533e-03 / 8.590e-03 / **1.007** | 9.449e-03 / 9.504e-03 / **1.006** |
| radius cuts at iteration 1 | 0 | 0 |
| iterations accepted of 20 | 14 | 16 |
| first refused iteration | 14 | 12 |
| `pred` where the model first fails | 4.634e-11 | 7.581e-11 |
| `\|\|R\|\|` at iteration 20 | 8.861e-01 | 5.472e-01 |

**The `ratio = -26.8` model failure L17 measured at the first radius does not
occur here, with the damping off OR on: the model is right to 0.7 per cent at
the first radius in both.** The refusals that do occur start at `pred` ~ 1e-11,
which is the merit's own rounding, and L17 recorded the same floor
(`below ||s|| ~ 1e-09 the measured change stops shrinking`). The disappearance
of the -26.8 is therefore not attributable to the damping: it is present in the
control branch as well, and what changed between L17 and today is the base
boundary (section 1.5).

#### The same comparison at the pre-L21 blend window

`EXHALE_BASE_MACH_BLEND=1.0e-6`, otherwise identical, 200 iterations:

| | 0.02, off | 0.02, on | 0.03, off | 0.03, on |
|---|---|---|---|---|
| mass row, max (cell) | 3.527e-01 (1) | 3.533e-01 (1) | 3.520e-01 (1) | 3.534e-01 (1) |
| energy row, max (cell) | 1.520e-01 (3) | 1.518e-01 (3) | 1.537e-01 (3) | 1.531e-01 (3) |
| `\|\|R\|\|` at iteration 200 | 3.569e-01 | 3.575e-01 | 3.535e-01 | 3.594e-01 |

The rows are lower than at the default window and the off/on difference is
under half a per cent, so the reading of the previous tables does not depend
on the window: with the damping on or off, and at either window, neither state
certifies and neither solve descends.

#### The third state of the brief

**NOT DONE.** `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` was set up from the
seed the brief names,
`.stopped/atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7_dtau0-1_20260916055251/output/*_IC.txt`,
and `load_IC` refuses that file on this tree:

```
(load_IC) ERROR: output/Hydro_ioniz_IC.txt: the key "certified" is stated twice.
  # mapped-from-coupling: sec_ion=T sec_ion_step=0 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
```

The refusal is not this state's: it belongs to every state
`src/utils/map_state_to_grid.py` writes, i.e. every mapped seed of
`run_case.sh` pass 0. See section 1.5.


### 1.5 What could not be measured, and why

READ, then MEASURED. `.L14/x003_HeH2.13` was written `certified=T` with its
hydrodynamic mass row at 3.764e-08 (cell 7), momentum 1.960e-14 (cell 1) and
energy 8.255e-07 (cell 7) by `EXHALE_L14.x`
(md5 `a11038c245050d4f11852af829e13fc7`) on 2026-09-14, and L17 re-measured the
same pair on 2026-09-15 at mass 1.619e-04 / momentum 8.045e-12 / energy
5.216e-04 for the transient. On this tree, on 2026-09-16, the same two states
as loaded read

| row | 0.02 transient | 0.03 control | L17 / the state's own run |
|---|---|---|---|
| hydrodynamic mass | 9.989e-01 at cell 2 | 9.984e-01 at cell 2 | 1.619e-04 at cell 1 / 3.764e-08 at cell 7 |
| hydrodynamic momentum | 1.863e-03 at cell 1 | 1.955e-03 at cell 1 | 8.045e-12 at cell 492 / 1.960e-14 at cell 1 |
| hydrodynamic energy | 1.000e+00 at cell 2 | 1.000e+00 at cell 2 | 5.216e-04 at cell 1 / 8.255e-07 at cell 7 |

and neither certifies. Three measurements locate this, and none of them is
L25:

1. **It is not this build.** The tree binary `EXHALE.x`
   (md5 `c2e9c9990b9f14f1be8cd77abca68945`) evaluating the 0.03 control gives
   the same three numbers to every printed digit.
2. **The state is loaded faithfully.** The state the run writes back agrees
   with the state it read to the last bit in every physical cell; the only rows
   that move are the two base GHOST cells, whose density moves by +118 and
   +117 per cent and whose velocity by 7.1e-02 and 2.5e-02 cm/s. The base
   boundary writes a different ghost state today.
3. **The base row is now its own flux divergence.** At the worst cell of the
   0.02 transient (energy, cell 2) the signed terms are `R = -9.0777e-02`,
   flux divergence `-9.0773e-02`, heat - cool `+3.250e-06`; in L17 the same
   cell had `R = -1.7186e-09` against a flux divergence of `+3.264e-06`. The
   base flux divergence is 2.8e+04 times what it was, which is what puts every
   row at its own scale.

**A second refusal, found while setting up the third state of the brief, and
reported because it is wider than this item.** `load_IC` on this tree refuses
every state file written by `src/utils/map_state_to_grid.py`, which is the seed
mapper `LHS1140b/models/run_case.sh` uses in its pass 0:

```
(load_IC) ERROR: output/Hydro_ioniz_IC.txt: the key "certified" is stated twice.
  # mapped-from-coupling: sec_ion=T sec_ion_step=0 recon=PLM certified=F cert_reason=no_stationary_claim mode=init
```

The mapper appends a `# mapped-from-coupling:` provenance line
(`map_state_to_grid.py` line 468) carrying the source file's `certified=` and
`cert_reason=`; the new duplicate-key refusal of `load_IC.f90`
(`parse_certification_claim`, uncommitted) is reached through the line selector
`index(line,'coupling:') > 0` at lines 397 and 481, which matches
`# mapped-from-coupling:` as well as `# coupling:`, and the claim is
accumulated across both lines. 25 of the seeds under
`LHS1140b/models/.stopped/*/output/` carry the line. The selector wants the
line's own label; `load_IC.f90` belongs to another item and was not touched.

`base_boundary.f90` carries uncommitted work in this tree (item L21 moved the
entropy-branch handover from the cell-centered product `rho_1 v_1 r_1^2` to the
Riemann face flux and narrowed the blend window from a face Mach number of
1e-6 to 1e-8, and item L26 measured the boundary's remaining discontinuity on
2026-09-16). Restoring the window alone with `EXHALE_BASE_MACH_BLEND=1.0e-6`
does NOT restore the rows (mass 9.976e-01 and 9.964e-01 at cell 2 on the two
states), so the width is not the whole of it.

**Consequence for this item.** The acceptance clause of rev3 section 5 for an
UNCHANGED operator branch, "the certified 0.03 state re-certifies to
tolerance", cannot be evaluated on this tree: the control does not re-certify
with the damping off. What section 1.4 measures is therefore the contrast
between the two branches under today's boundary, which is the question step 1
was opened for, and not the absolute standing of either state.

### 1.6 What step 1 says about the candidates, neutrally

- **The damping as written is not a dissipative operator on this grid.** That
  is a property of the term and not of the case: it holds on the certified
  state too, and its negative mode is created by the term's own gate. Any use
  of the present form puts an anti-dissipative term in the wind at 3.5-5 R_p,
  inside `r_esc = 2 R_p`'s exterior but well below the outer boundary, whose
  size on this state is small against what the term removes elsewhere but is
  not bounded by anything in the term.
- **The correctly signed form exists, is a flux, and costs nothing where the
  coefficient is constant.** Adopting it would be a change of the operator, so
  by R37 the control would be re-solved and judged on its own residual rather
  than against the old state. It is written down here and not inserted.
- **The damping acts on the base mode, and the base mode is not what holds the
  case.** Both forms damp the alternating mode as a quadratic form (section
  1.3), and on the frozen states the damped solve ends with a third less of it
  over cells 2-8 (section 1.4). It still ends with three times what the loaded
  state carried, the base rows still stand five to six decades above their
  tolerances, `Mdot` is unchanged to six digits and neither state certifies.
  Step 1 therefore does not support the reading that the damping repairs the
  low-XUV base, and it does not support the reading that the base mode is
  merely cosmetic either: what it shows is that removing a third of the mode
  leaves the rows where they were.
- **Where the kink is remains open after step 1.** L26 measured the base
  BOUNDARY to be discontinuous at a zero window flux by 37 per cent of the face
  density and to have a slope kink of 1.9e+04 when the window's extremal cell
  changes; L17 attributed a kink to the HLLC star-flux selection at the base
  cells. Step 1 adds a third object with a kink, the damping's own gate, but
  its gate sits at r ~ 4-5 R_p and not at the base, so it is not a candidate
  for the base kink. Nothing here chooses between the boundary and the flux;
  step 2's three-way flux comparison is the measurement that separates them.

## Step 2: the velocity-jump factor on the Roe branch

Item L25 step 2 of `docs/PLAN_20260916_rev3.md` section 5, with rows R9, R10,
R11, R35, R37 and R47 of its section 0. Written by the step 2 worker; step 1
above is another worker's and is untouched. Every number below is MEASURED on
this tree on 2026-09-16 at one thread unless it is marked READ.

### 2.1 Verdict

**The two contrasts separate cleanly, and they point in opposite directions.**

- **Roe minus HLLC (the flux family).** On the stalled target
  `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` the unmodified Roe flux takes the
  three hydrodynamic rows from `||R||` 1.000 to **2.744e-08 in 23 Newton
  iterations of one outer pass**, where HLLC stands at **1.779 after 60**. On
  the certified 0.10 control the unmodified Roe flux re-solves to `||R||`
  8.566e-08 in 22 iterations and leaves only the gated He/H partition row
  above its tolerance. The flux-family change is where the movement is.
- **Corrected Roe minus Roe (the Mach factor).** The factor leaves the mass and
  the energy rows of the frozen states unchanged to four digits and multiplies
  the judged momentum residual of base cells 3-6 by three to five decades; it
  turns a trust-region ratio of 1.00 into ratios of -74, -22 and -47 with four
  of five trials refused; and on both the control and the target it **prevents**
  the solve that the unmodified Roe flux completes, ending with the base
  velocity in an almost pure odd-even mode (normalized Nyquist projection over
  cells 1-6: 0.99 on the control, 0.90 on the target, against 0.09 and 0.14 in
  the states they started from) and, on the target, with a mass flux at cell
  400 that has changed sign.

The mechanism is the one the unit test isolates and the paper predicts: the
factor is the only damping the odd-even velocity mode has at these faces, and
where the two neighboring velocities cancel in the Roe average the factor is
exactly zero, so that mode is left entirely undamped. Rieper's asymptotic
premise (his section 3.2, the local and the global Mach number of comparable
order) is not met by this column: the global Mach number is 1e-2 at the outer
boundary and 1e-8 at the base, so `min(Ma_face, 1)` is applied eight decades
apart on one grid.

**Nothing here chooses a production design.** The option is delivered off, and
the Roe-minus-HLLC contrast is a separate finding that belongs to L25's own
decision, not to this worker's.

### 2.2 What was implemented

`src/modules/flux/Num_Fluxes.f90`, ROE branch, immediately after the jumps are
formed and before the expansion coefficients:

```
if (low_mach_velocity_jump)                                  &
   dvel = min(abs(v_avg)/a_avg, 1.0d0)*dvel
```

`v_avg` and `a_avg` are the Roe-averaged velocity and sound speed the branch
has already formed, so this is Rieper (2011, J. Comput. Phys. 230, 5263)
eq. 3.15 with the local Mach number of his eq. 3.16, `Ma = (|U|+|V|)/a`, in the
one-dimensional case where the tangential component `V` is absent. The scaled
jump enters `a1` and `a3` and nothing else: `a2` carries `drho` and `dp` only.
The eigenvalues, the eigenvectors, the central flux, the well-balanced pressure
departure `dp_wb`, the entropy correction and the admissibility test
(`speed_estimate_ROE`, the HLLE branch) are untouched.

The module header states the physics, the validity (his section 3.2 premise;
his section 5, an accuracy correction and not a stiffness one, HLL and Rusanov
not suited, an HLLC adaptation needing its own derivation) and the default.

The input key is `Low Mach velocity jump: True|False`, default False, parsed in
`input_read.f90` beside `Well balanced` and echoed by
`write_setup_report.f90` (those two files were edited for this key alone, as the
brief allows). The flag is the module variable
`Numerical_Fluxes::low_mach_velocity_jump`; the setup report warns that the key
has no effect when `Numerical flux:` is not `ROE`. One defect of this worker's
own first draft is recorded because it would recur: the value is `get_word(line, 5)`,
not word 3, because the key is four words, and with word 3 the parser read the
word `velocity` and the option silently stayed off. The setup report echo is
what caught it.

### 2.3 The unit test, RED and GREEN

`src/tests/physics_probe/roe_low_mach_velocity_jump.f90`, added to the driver
list of `src/tests/physics_probe/run.sh`. It calls the production `Num_flux` on
the ROE branch and compares against the analytic three-wave Roe flux formed in
the test from its own Roe averages, so the reference is independent of the
branch and not a snapshot of it. Twenty-four assertions:

| group | what it asserts |
|---|---|
| `roe_key_off_row_*` | with the key off the returned flux is the uncorrected Roe flux of the analytic reference, 1e-12 relative |
| `roe_key_off_repeatable_row_*` | a corrected call between two uncorrected ones moves the uncorrected flux by exactly zero |
| `roe_supersonic_*` | at `Ma_Roe` = 4.039 the corrected flux is bitwise the uncorrected one (tolerance 0) |
| `roe_low_mach_face_mach`, `roe_low_mach_scaled_row_*` | at `Ma_Roe` = 9.99999652e-04 the flux equals the analytic reference with the velocity jump scaled by that Mach number, 1e-12 relative |
| `roe_pressure_jump_part_unchanged_row_*` | the flux at factor zero implied by the two calls (the flux is affine in the factor) equals the analytic flux with the velocity jump removed and the pressure-jump and entropy parts intact, 1e-9 relative |
| `roe_alternating_*` | for `v_L = -v_R` at equal density and pressure the Roe-averaged velocity is exactly zero, the corrected flux is the central flux to 1e-15 in every row, and the damping removed is the momentum row's `-rho_avg a_avg dvel / 2` |

RED, MEASURED: with the factor replaced by `1.0d0` in a scratch copy of the
tree (the option present, the scaling absent), 7 of the 24 assertions fail,
`roe_low_mach_scaled_row_{m,p,e}`,
`roe_pressure_jump_part_unchanged_row_{m,p,e}` and
`roe_alternating_momentum_damping_removed` and
`roe_alternating_central_flux_row_p`; the six that describe the key off and
the supersonic face pass, which is what separates the two halves of the
option. GREEN, MEASURED: with the delivered source the whole
`physics_probe` suite passes, 1515 PASS and 0 FAIL, this driver included; the
`grid_and_gates` driver `free_outflow_boundary`, which also links `Num_flux`,
passes as well. (`riemann_wave_speeds` and `positivity_limiter_scaling` are in
the same suite and are covered by that run.)

**What assertion (4) measures, which is the probe rev3 asked for.** At a face
whose two velocities cancel in the Roe average the factor is exactly 0, and the
velocity-jump dissipation is removed in full. The MASS row of that pair is zero
with and without the factor: with `dp = 0` the two acoustic coefficients are
equal and opposite and their mass-row contributions cancel, and the central
mass flux is `rho(-u + u)/2 = 0`. What the factor removes is entirely in the
MOMENTUM row, `-rho_avg a_avg dvel / 2`, measured as -1.290994450027e-04
against the analytic -1.290994450027e-04 for `rho = p = 1`, `u = 1e-4`,
`gamma = 5/3`. That term is the only damping the 2 dr velocity mode has at such
a face, so with the factor on the mode is undamped there. Sections 2.7-2.9
below are that statement at the scale of a column.

### 2.4 Byte identity with the key off

Two binaries built from ONE snapshot of the tree taken at 19:17 KST, differing
only by this worker's three hunks: `EXHALE_ctl.x` md5 `40d1614319d8` (the
hunks reverted) and `EXHALE_new.x` md5 `3d4c08e7b86b`. Both regression cases
were run on scratch copies at `OMP_NUM_THREADS=1` with the case's own
`maxsteps`, MEASURED:

| case | `Numerical flux:` in its `input.inp` | final marching count | data rows |
|---|---|---|---|
| `wasp_he23off` | HLLC | 12236, `du` 9.7535E-04, identical | `Hydro_ioniz.txt`, `Ion_species.txt`, `Hydro_ioniz_adv.txt`, `Ion_species_adv.txt` all byte-identical |
| `mol_base_handoff` | HLLC | 12000, `du` 7.9006E-01, identical | the same four files byte-identical |

The only difference in any output file is the `run=` timestamp of the
provenance header line, which is what `run_check.sh` already excludes when it
reports "data identical".

**Both cases select HLLC, so neither exercises the branch the option lives on.**
The Roe-branch identity was measured separately and is stronger: the frozen
0.02 state was solved for 200 Newton iterations under `Numerical flux: ROE`
with the key off by both binaries, and the two runs agree in **every printed
`(JFNK)` line of the whole trajectory (zero differing lines over 108
iterations) and in every data row of `Hydro_ioniz.txt` and `Ion_species.txt`**.
With the key off the delivered branch is the entry text to the bit.

### 2.5 The comparison, and which states carry which measurement

Three flux settings, everything else fixed (boundary, reconstruction `PLM`,
`Well balanced: True`, seed, source terms, solver keys, `EXHALE_PTC_DTAU0` as
the case had, `EXHALE_OUTER_PASSES=1`, one thread): unmodified HLLC, unmodified
ROE, ROE with `Low Mach velocity jump: True`. The only edits to a case's
`input.inp` are the `Numerical flux:` line, the added key, and an absolute path
for `Spectrum file:` (the same file).

| label | state | role | `dtau0` |
|---|---|---|---|
| A | `.L14/x002_HeH2.13/output` (the stalled 0.02 transient) | FROZEN state only | 1.0e8 |
| B | `.L14/x003_HeH2.13/output` | FROZEN state only | 1.0e8 |
| C | `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13`, its `output/*_IC.txt` (`certified=T`) | the control, RE-SOLVED | 1.0 |
| D | `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7`, seeded from `.stopped/...dtau0-1_20260916055251/output/*_IC.txt` | the stalled target | 1.0 |

A and B are used only as frozen states: both were certified before L21 changed
the base boundary and neither re-certifies on this tree, so no re-solve of them
is an acceptance statement (the step 1 worker measured that; READ). C is the
control certified on the tree binary `c2e9c9990b9f`, and by R37 it is re-solved
under each operator and judged on its own residual. D is the case the item is
about; it was run with `EXHALE_JFNK_MAXIT=60` rather than 200 and on a binary
built from a later snapshot of the tree (md5 `bd96d010ad88`, which carries the
`load_IC` repair of the mapped-seed refusal), so its three settings are
comparable with each other and not digit-for-digit with A, B and C.

### 2.6 Spatial accuracy: the frozen-state residual of base cells 1-6

The residual `F` at outer iteration 1, before any Newton step, dumped by
`EXHALE_L15_TRACE` with `EXHALE_L15_DUMP_AT=1` and `EXHALE_JUDGED_ROWS=1`, and
divided by the row scale `Drow` the certification judges that row at. MEASURED.

A, the 0.02 frozen state:

| row | setting | cell 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| mass | HLLC | +6.94e+07 | -3.02e+09 | +5.51e+01 | +1.54e+02 | +7.12e+01 | -7.97e+01 |
| mass | ROE | -3.11e+09 | -3.02e+09 | -3.58e+04 | +6.27e+01 | -4.08e+01 | +1.95e+02 |
| mass | ROE + factor | -3.11e+09 | -3.02e+09 | -3.58e+04 | +6.89e+01 | -4.02e+01 | +1.96e+02 |
| momentum | HLLC | +1.86e+05 | -4.10e+04 | +3.80e-04 | -2.23e-04 | -5.00e-05 | -8.83e-05 |
| momentum | ROE | +1.28e+05 | -4.10e+04 | -2.49e-02 | -6.54e-02 | +4.61e-02 | -1.39e-02 |
| momentum | ROE + factor | +8.14e+04 | -3.12e+04 | -2.23e+03 | +1.15e+03 | -3.84e+02 | +1.58e+02 |
| energy | HLLC | +1.00e+06 | -1.00e+06 | +1.18e+02 | +1.89e+02 | +5.50e+01 | -2.11e+02 |
| energy | ROE | -1.00e+06 | -1.00e+06 | -2.25e+05 | -1.67e+01 | +2.17e+03 | -1.12e+03 |
| energy | ROE + factor | -1.00e+06 | -1.00e+06 | -2.24e+05 | +7.58e+01 | +2.18e+03 | -1.12e+03 |

C, the certified control (the HLLC column is the state's own operator, so it is
the state's distance from its own root):

| row | setting | cell 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| mass | HLLC | -5.69e-01 | +1.34e-01 | +1.11e-01 | -9.37e-02 | -4.97e-02 | +1.16e-01 |
| mass | ROE | -6.18e+08 | -6.50e+05 | -3.31e+04 | +1.02e+05 | -1.19e+04 | +1.66e+04 |
| mass | ROE + factor | -6.17e+08 | -6.51e+05 | -3.29e+04 | +1.02e+05 | -1.19e+04 | +1.66e+04 |
| momentum | HLLC | +1.56e-06 | +1.28e-06 | +2.71e-07 | -7.80e-07 | +4.25e-07 | +4.76e-07 |
| momentum | ROE | -3.39e+03 | +2.03e+01 | +1.76e+00 | -3.26e-01 | +3.12e-01 | -7.03e-02 |
| momentum | ROE + factor | -8.01e+04 | +1.75e+04 | -3.28e+03 | +2.08e+03 | -6.75e+02 | +2.98e+02 |
| energy | HLLC | -3.89e-01 | +1.30e-01 | +3.32e-02 | -4.49e-02 | -3.89e-02 | +7.76e-02 |
| energy | ROE | -1.00e+06 | +8.04e+05 | +4.33e+04 | -7.16e+04 | +1.20e+04 | -5.43e+03 |
| energy | ROE + factor | -1.00e+06 | +8.04e+05 | +4.36e+04 | -7.16e+04 | +1.20e+04 | -5.43e+03 |

D, the target's seed:

| row | setting | cell 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| mass | HLLC | -3.32e+08 | -4.03e+09 | +1.61e-01 | -1.96e-01 | +2.50e-01 | -6.56e-03 |
| mass | ROE | -4.79e+09 | -4.03e+09 | -3.59e+04 | +3.15e+02 | -7.08e+02 | +3.36e+02 |
| mass | ROE + factor | -4.79e+09 | -4.03e+09 | -3.59e+04 | +3.25e+02 | -7.09e+02 | +3.36e+02 |
| momentum | HLLC | +2.68e+05 | -5.30e+04 | -1.55e-07 | -1.73e-06 | +1.14e-06 | -4.03e-07 |
| momentum | ROE | +1.79e+05 | -5.30e+04 | -5.12e-02 | -6.10e-02 | +4.71e-02 | -1.45e-02 |
| momentum | ROE + factor | +1.29e+05 | -4.29e+04 | -2.30e+03 | +1.19e+03 | -3.96e+02 | +1.63e+02 |
| energy | HLLC | +1.00e+06 | -1.00e+06 | +4.64e-01 | -2.75e-01 | +1.57e-01 | +2.09e-02 |
| energy | ROE | -1.00e+06 | -1.00e+06 | -4.96e+04 | +8.11e+02 | +4.15e+02 | -2.17e+02 |
| energy | ROE + factor | -1.00e+06 | -1.00e+06 | -4.95e+04 | +8.40e+02 | +4.14e+02 | -2.15e+02 |

The two contrasts, read off the three tables:

- **Roe minus HLLC** moves every row by decades on all three states, which is
  what a different operator at a fixed state does and not by itself a defect of
  either. On the control, where the state IS the HLLC root, it is the whole
  story: 1e-1 to 6e8 in the mass row, 1e-6 to 3e3 in the momentum row.
- **Corrected Roe minus Roe** leaves the mass and energy rows unchanged to
  three or four digits on all three states, and acts only on the MOMENTUM row:
  at cells 3-6 it multiplies it by 3e+04 to 1e+05 on A, by 2e+03 to 6e+03 on C
  and by 2e+04 to 4e+04 on D; at cells 1-2 it lowers it by a factor 1.5 to 2 on
  A and D and raises it by a factor 24 to 860 on C.

The same contrast in the two summary quantities the solve prints at that
iterate, MEASURED, where the volume-weighted norm and the cell count disagree
in sign:

| state | setting | `mom/grav` below `r_esc` | cells outside the momentum tolerance, of 500 | cells outside the mass tolerance |
|---|---|---|---|---|
| A | HLLC | 1.863e-03 | 2 | 496 |
| A | ROE | 1.278e-03 | 7 | 497 |
| A | ROE + factor | 8.139e-04 | 63 | 495 |
| B | HLLC | 1.955e-03 | 2 | 2 |
| B | ROE | 1.339e-03 | 9 | 21 |
| B | ROE + factor | 8.719e-04 | 73 | 178 |
| C | HLLC | 3.400e-14 | 0 (the gate is met at the loop top) | 0 |
| C | ROE | 3.392e-05 | 60 | 24 |
| C | ROE + factor | 8.005e-04 | 97 | 210 |
| D | HLLC | 2.678e-03 | 2 | 2 |
| D | ROE | 1.787e-03 | 64 | 21 |
| D | ROE + factor | 1.288e-03 | 97 | 223 |

On A, B and D the factor LOWERS the volume-weighted fractional violation of
hydrostatic balance below the escape radius by a third while RAISING the number
of cells whose momentum row is outside its own tolerance by an order of
magnitude; on C, where the comparison starts from a converged HLLC state, both
rise. A volume-weighted norm and a worst-cell count are different functionals
and this is where they part.

**The alternating velocity of cells 1-6** is a property of the state and not of
the flux, so at iteration 0 it is the same number for all three settings: the
normalized projection of `v(1:6)` on its odd-even mode is 0.1335 for A, 0.1334
for B, 0.0897 for C and 0.1364 for D. What each operator does to it is in
section 2.9.

### 2.7 Jacobian-action regularity

`EXHALE_TRUST_REGION=1 EXHALE_TR_TRACE=1 EXHALE_JFNK_MAXIT=5`, one outer pass,
one thread, the same frozen states. The trust region's first radius is measured
on the Krylov step, and `ratio` is the measured change of the squared merit
over the change the banded model predicts; 1 is a model that describes the
state. MEASURED:

| state | setting | first radius | ratio, iterations 1-5 | trials refused |
|---|---|---|---|---|
| A | HLLC | 3.873e-05 | 1.007, 1.000, 1.000, 1.000, 1.000 | 0 |
| A | ROE | 6.434e-01 | 1.011, 1.043, 1.034, 1.108, 1.106 | 0 |
| A | ROE + factor | 1.263e+02 | 6.157, -2262.8, -362.6, -84.6, -16.0 | 4 |
| B | HLLC | 2.108e-01 | 1.028, 1.002, 1.003, 1.004, 0.162 | 0 |
| B | ROE | 4.767e-01 | 1.019, 1.052, 1.179, -1.995, 1.225 | 1 |
| B | ROE + factor | 3.873e+04 | overflow, overflow, overflow, overflow, 38593 | 4 |
| C | HLLC | 3.873e-05 | 0.486, 0.699, 0.652, 0.620, 0.563 (squared merit 1e-22: rounding) | 4 |
| C | ROE | 3.873e-05 | 1.000, 1.000, 1.000, 1.001, 1.002 | 0 |
| C | ROE + factor | 1.682e-04 | -7.330, -0.864, 1.000, 1.000, -2.602 | 3 |
| D | HLLC | 3.873e-05 | 1.009, 0.947, 0.998, 0.997, 0.999 | 0 |
| D | ROE | 1.094e-03 | 1.001, 1.001, 1.002, 1.001, 1.002 | 0 |
| D | ROE + factor | 1.676e-03 | 2201.3, -73.97, -21.59, -21.82, -46.84 | 4 |

The GMRES forcing was 1.00e-01 at every one of these iterations and the cycle
reached it, in 1 to 11 products of a subspace of 40; no setting exhausted the
subspace here, so the linear tolerance is not what separates them.

**The contrast is entirely in the second one.** Roe minus HLLC leaves the ratio
at 1 on every state (B iteration 4 is the single exception, one refused trial).
Corrected Roe minus Roe turns it into two- to four-digit numbers of the wrong
sign with most trials refused, on all four states. The first trust radius,
which is the length of the Krylov step, grows by a factor 200 (A), 8e4 (B) and
1.5 (D) when the factor is switched on: the linearized system's Newton step
becomes long, which is what the removal of the damping of a mode does to the
operator's smallest singular direction.

For reference, L17 reported `ratio` -26.804 with `EXHALE_TRUST_REGION=1` on
`x0.10/HeH9.7` under HLLC (READ, `docs/lhs1140b_stationary_L17_20260915.md`
section 5). On this tree and from this seed, HLLC on that case gives 1.009,
0.947, 0.998, 0.997, 0.999 in the same five iterations: the state and the base
boundary have both moved since (L21), so the two are not the same measurement
and the older number is not reproduced or refuted here. What carries the
signature now is the corrected Roe.

### 2.8 Nonlinear convergence

One outer pass, `EXHALE_JFNK_MAXIT` 200 on A and C and 60 on D, one thread.
MEASURED:

| state | setting | `||R||` at entry | at iteration 60 | at the end | iterations | how it ended | worst row and cell at the end |
|---|---|---|---|---|---|---|---|
| A | HLLC | 1.000 | 1.823 | 1.119 | 173 | no descent along either direction | energy of cell 54, 1.119 |
| A | ROE | 1.000 | 1.773 | 1.403 | 108 | no descent | energy of cell 2, 9.806e-01 |
| A | ROE + factor | 1.000 | not reached | 1.915 at iteration 25 | 25 | the leg was stopped at a bounded horizon (see below) | momentum of cell 1 |
| C | HLLC | 3.886e-07 | - | 3.886e-07 | 0 | the acceptance gate met at the loop top | CERTIFIED |
| C | ROE | 1.004 | - | 8.566e-08 | 22 | the gate met, the energy row still 4.114e-06 of 1e-06 at cell 1, then converged | only the gated He/H partition, 3.454e-04 at cell 253 |
| C | ROE + factor | 1.004 | 1.795 | 1.848 | 142 | damped Gauss-Newton escape, no progress | mass 1.601 at cell 2, energy 1.571 at cell 6 |
| D | HLLC | 1.000 | 1.779 | 1.779 | 60 (the cap) | the cap | mass 5.616e-01 at cell 1, energy 1.095 at cell 2 |
| D | ROE | 1.000 | - | 2.744e-08 | 23 | converged | only the gated He/H partition, 2.507e-04 at cell 217 |
| D | ROE + factor | 1.000 | 1.440 | 1.440 | 60 (the cap) | the cap | mass 4.081e-03 at cell 322, momentum 1.080e-02 at cell 392, energy 8.973e-01 at cell 3 |

The A leg with the factor is a BOUNDED HORIZON and not a completed run: its
first attempt was discarded because two processes of this worker wrote into one
directory (the first was launched before the key-parsing defect of section 2.2
was found and its binary outlived the shell that was stopped), and the clean
rerun was itself stopped at iteration 25. What the three A legs say at the one
iteration they all reached, MEASURED at iteration 20:

| setting | `||R||` | merit `||Fs||2` |
|---|---|---|
| HLLC | 1.974 | 5.69e-04 |
| ROE | 1.523 | 8.93e-04 |
| ROE + factor | 1.965 | 4.15e-02 |

The merit is two decades worse with the factor at the same iteration, which is
the same direction as C and D; the A state is not a certified control on this
tree and no acceptance statement is made from it.

- **Roe minus HLLC**: on C and on D the flux family decides the outcome. HLLC
  leaves D at 1.779 after 60 iterations, which is the stall L25 is about; Roe
  reaches 2.744e-08 in 23. On the control Roe re-solves in 22 from a state 1.004
  away in its own residual. On A neither converges.
- **Corrected Roe minus Roe**: on C the factor turns a 22-iteration convergence
  into 142 iterations with no progress; on D it turns a 23-iteration convergence
  into the iteration cap at 1.440. The factor removes the convergence the flux
  family gained.

### 2.9 The control and the target, re-solved (R37)

The control is judged on its own residual and on the movement it undergoes, not
against the old state. MEASURED, cells 1-6 and the mass flux at the outermost
cell:

| quantity | control under HLLC | control under ROE | control under ROE + factor |
|---|---|---|---|
| `d rho` cells 1-6 | 0 (the state is its root) | +7.5, +14.8, +17.1, +16.1, +12.4, +9.0 per cent | -1.7, +1.1, +0.4, +1.2, -0.6, +0.3 per cent |
| `d T` cells 1-6 | 0 | -7.2, -13.7, -16.1, -16.1, -14.1, -11.9 per cent | +1.8, -1.1, -0.4, -1.2, +0.6, -0.3 per cent |
| `v` cell 1 [cm/s] | -6.572e-01 | -3.523e-01 | +2.675e+00 |
| odd-even projection of `v(1:6)` | 0.0897 | 0.2022 | 0.9929 |
| `Mdot` at cell 500 [g/s] | 6.56161e+06 | 6.52184e+06 (-0.61 per cent) | 6.56283e+06 (+0.02 per cent) |
| `Mdot` at cell 400 [g/s] | 6.57705e+06 | 6.53717e+06 | 3.74006e+06 (-43 per cent) |
| certification | CERTIFIED | only the gated He/H partition row refuses | four rows refuse |

The target, from its `.stopped` seed:

| quantity | seed | HLLC, 60 iterations | ROE, converged | ROE + factor, 60 iterations |
|---|---|---|---|---|
| `d rho` cells 1-6 | - | +14.7, -0.3, -0.8, -0.7, -0.7, -0.7 per cent | +154.8, +142.7, +121.2, +97.0, +71.8, +52.3 per cent | +34.0, +3.1, -2.2, -1.8, -1.8, -1.8 per cent |
| `T` cell 1 [K] | 578.5 | 503.1 | 221.8 (-61.7 per cent) | 429.6 |
| `v` cell 1 [cm/s] | -1.293e-01 | -3.629e+00 | -4.852e-01 | -2.520e-04 |
| odd-even projection of `v(1:6)` | 0.1364 | 0.5143 | 0.1997 | 0.9049 |
| `Mdot` at cell 500 [g/s] | 6.68095e+06 | 6.68062e+06 | 6.41484e+06 (-3.98 per cent) | 6.66126e+06 |
| `Mdot` at cell 400 [g/s] | 6.69155e+06 | 6.69109e+06 | 6.42492e+06 | -2.70404e+07 (sign reversed) |

Two things are worth stating plainly. First, R9 is confirmed: a base-local
change moves `Mdot`, by 0.61 per cent on the control and 3.98 per cent on the
target under the unmodified Roe flux, and the Roe solution is a base one to two
hundred Kelvin colder and up to 155 per cent denser than the HLLC one. Whether
that is the better solution is a physics question this measurement does not
answer, and grid convergence of it was not run. Second, the corrected Roe
states are not solutions at all: the outer wind mass flux is 43 per cent low on
the control and has changed sign on the target, and the base velocity is an
almost pure 2 dr mode.

### 2.10 Where the base kink is, given L26

L26 measured the base BOUNDARY to be discontinuous at a zero window flux by 37
per cent of the face density and to carry a slope kink of 1.9e+04 when the
window's extremal cell changes (READ, `docs/lhs1140b_stationary_L26_20260916.md`).
This step's contribution to that question is narrow and is stated as such:

- The boundary was IDENTICAL in all twelve runs. Every difference in this
  section is therefore a flux difference, and none of it can be attributed to
  the boundary.
- Under the unmodified Roe flux the control and the target both solve, the
  target in 23 iterations, with that same boundary. **So the boundary's kink, on
  these two states, is not what prevents the hydrodynamic rows from converging.**
  That is the sharpest statement here and it belongs to L26 as much as to L25.
- It does not follow that the boundary is smooth, or that the HLLC star-flux
  selection is the cause of the HLLC stall. The two operators differ in more
  than the contact selection, and a state reached by a different flux may simply
  sit where the boundary's kinks are not active. Separating those needs the
  frozen-state directional-derivative probe of L26 under both fluxes, which was
  not run here.

### 2.11 A neutral reading

- **Accuracy changed, and by both contrasts.** Roe minus HLLC changes every row
  of the frozen residual by decades, which is the operator change. Corrected Roe
  minus Roe changes only the momentum row, lowering its volume-weighted norm
  below the escape radius by a third on three of four states and raising its
  worst cells by three to five decades on all four.
- **Regularity changed, and by the second contrast alone.** The banded model
  keeps `ratio` ~1 under HLLC and under unmodified Roe and loses it under the
  factor on every state tried.
- **Convergence changed, and by both, in opposite directions.** The flux family
  is what lets the control and the target solve; the factor is what stops them.
- **The factor's validity condition is not met here.** Rieper's section 3.2 asks
  for a local Mach number of the same order as the global one; this column spans
  Mach 1e-8 to 1e-2, and his section 5 says the correction is for accuracy and
  removes no stiffness. The measurements are consistent with the factor being
  applied outside the regime it was derived in, and the mechanism the unit test
  isolates (factor exactly 0 where the neighboring velocities cancel, so the
  odd-even mode keeps no damping) is visible in the states that come out.
- **No production design is chosen here.** What this step delivers is the
  option, off by default, with its test and its measurement. The Roe-minus-HLLC
  result is a finding for L25's decision and is deliberately not acted on:
  adopting a flux family would change every certified state in the catalog, and
  by R37 each would have to be re-solved and judged on its own residual,
  conservation and grid convergence, none of which was done here.

### 2.12 What was not measured

- Grid convergence of any of these states, under any flux.
- The standing L25 programme's other required cases (a stationary contact, a
  resolved low-speed flow with a fixed background, a compressible hot-Uranus or
  WASP control). The unit test covers the stationary-contact and cancelling-
  velocity limits at the level of a single face and nothing beyond that.
- A gate on the factor. The plan asks for a dimensionless one, the face Mach
  number against a power of `delta p / p_ref` with the exponent derived for the
  gravitational equilibrium and the caloric EOS. None is implemented: the
  delivered option applies the published factor at every face of the ROE branch
  and nothing else.
- Anything with more than one outer pass, so the He/H partition row, which is
  alternated outside the Newton solve, is nowhere converged; on C and D under
  the unmodified Roe flux it is the only row that refuses certification.
- `EXHALE_TR_TRACE=1` was set on the convergence runs as well, and
  `EXHALE_TRUST_REGION=1` only on the five-iteration probes, so the step control
  of the convergence table is the tree's default and not the trust region.

## Step 3: the Roe flux on the low-XUV cases

Item L25 step 3 of `docs/PLAN_20260916_rev3.md` section 5, with rows R9, R35
and R37 of its section 0. Written by the step 3 worker; steps 1 and 2 above are
other workers' and are untouched. **No source was changed and nothing was
built.** Every run used the tree binary `EXHALE.x`, md5
`c2e9c9990b9f14f1be8cd77abca68945`, verified before the series. Every number
below is MEASURED on the host `lart3` on 2026-09-16/17 at 8 threads unless it
is marked READ. The run directories are `LHS1140b/models/.L25/`; no catalog
directory and no golden was touched.

### 3.1 Verdict

- **Three of the seven unsolved low-XUV cases CERTIFY under the unmodified Roe
  flux, with the catalog recipe otherwise unchanged**, and a fourth reaches a
  stationary wind whose three hydrodynamic rows are all inside their tolerances
  and is refused only by the gated elemental He/H partition row at 10.85 R_p.
  Three were still running at the end of the window and are reported with their
  state. The headline case of the item,
  `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7`, is one of the three: it certifies
  at `||R||` 4.143e-07 on outer pass 14, where HLLC from the same seed spent
  its forty passes without a root.
- **The flux-family systematic where both solve is small in the wind and a few
  per cent at the base.** On the two certified controls that finished, `Mdot`
  moves by -0.31 and -0.62 per cent, the He 10830 equivalent width by -0.31 and
  -0.44 per cent, and the first-cell temperature and mass density by -5.7 / +5.9
  and -7.2 / +7.5 per cent. The mass-flux spread over the wind window is
  unchanged to three digits.
- **Resolution answers the question the two fluxes disagreed on, and it answers
  it in the Roe solution's favour on convergence and against a large base
  difference.** On the doubled grid the UNMODIFIED HLLC flux certifies the same
  `x0.10/HeH9.7` case that it cannot solve at the catalog resolution, and the
  state it reaches is the Roe 500-cell state: the two agree to 1.0 per cent in
  `Mdot`, 0.4 per cent in the equivalent width and 0.1 to 3 per cent in T and
  rho everywhere above 1.005 R_p, and they part only in the bottom two or three
  cells, where the two grids also have different first centers. So the 500-cell
  HLLC stall is a resolution failure of that flux on that grid and not a
  different wind.
- **Nothing here chooses a production design.** What a decision to adopt ROE
  would require is stated in section 3.6.

### 3.2 What was run, and with what

The catalog recipe of `LHS1140b/models/run_case.sh` with `Numerical flux: ROE`
as the ONLY change to the case's `input.inp` (the `Spectrum file:` path was
made absolute, to the same file). `OMP_NUM_THREADS=8`,
`EXHALE_OUTER_PASSES=40`, `SEED_ATTEMPTS=1`, `FORCE=1`. The seven low-XUV cases
were taken at `EXHALE_PTC_DTAU0=1.0e8`, which is what `models/run_lowxuv_dtau0.sh`
prescribes for them and what the catalog attempt used; the certified controls
were taken at `EXHALE_PTC_DTAU0=1.0`, which is what their own `REPRODUCE.md`
records (READ).

Seeds, stated as the brief asks. For `x0.10_kzz1e9/HeH9.7` and the three
photochem cases at 0.15, 0.20 and 0.25 the seed is the same one the catalog
attempt used, `models_20260914_preL21/<the same case>/output`, which is what
`pick_seed.py` returns today as tier0 (MEASURED by running it). The three cases
at 0.01 of the fiducial spectrum have no tier0 state at all, so they take
`pick_seed.py`'s current first choice, which is a certified case at another XUV
normalization: `x0.10_kzz1e9/HeH2.13` for `x0.01/HeH2.13`,
`x0.15_kzz1e9/HeH9.7` for `x0.01/HeH9.7` and
`photochem x0.10_kzzprofile/HeH9.7` for `photochem x0.01/HeH9.7`. The catalog
attempt of 2026-09-14 recorded no seed path in the files it left, so that
choice is stated as today's and not claimed to be the identical file.

The certified controls were re-solved from their own certified states
(`models/<case>/output`), which is what R37 asks for.

### 3.3 Item 1: the seven low-XUV cases under ROE

MEASURED. `Mdot` is 4 pi rho v r^2 of the state's outermost physical cell,
computed from the state; the equivalent width is the red-pair EW over the
measurement's vacuum window 10832.60-10834.20 A, read the way
`LHS1140b/make_memo_figures.py` reads it. The odd-even entry is the normalized
Nyquist projection of `v` over physical cells 1-6.

| case | outcome | pass | `||R||` | refusing row | `Mdot` [g/s] | EW [%A] |
|---|---|---|---|---|---|---|
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | **CERTIFIED** | 14 | 4.143e-07 | none | 6.4102e+06 | 0.2353 |
| `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7` | **CERTIFIED** | 6 | 8.045e-08 | none | 1.3560e+07 | 0.6711 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | **CERTIFIED** | 39 | 5.749e-07 | none | 6.1609e+05 | 1.370e-06 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | not certified, pass budget spent | 41 | 1.855e-07 | elemental He/H partition, 2.234e-04 of 1.0e-05 at cell 440, r = 10.85 | 5.5041e+05 | 3.358e-05 |
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | still running at the end of the window | 6 | 1.29 | - | - | - |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | still running | 9 | 1.67 | - | - | - |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | still running | 10 | 1.94 | - | - | - |

The fourth case is worth stating plainly: its three hydrodynamic rows are all
WITHIN their tolerances at the end (mass 9.314e-08 of 1.3e-06 at cell 1,
momentum 1.039e-14 of 1.0e-08 at cell 79, energy 1.855e-07 of 1.0e-06 at cell
52), and what refuses it is the composition in the outer wind. That is the
situation `run_case.sh` already describes as "no seed can be nearer to the wind
than the wind the solve already found"; under HLLC this case never reached it.

The base of the certified states, MEASURED, cells 1-6:

| case | T(1..6) [K] | rho(1) [cm^-3] | odd-even of v(1:6) | Mach at cell 1 |
|---|---|---|---|---|
| `x0.10_kzz1e9/HeH9.7` | 221.78 225.98 241.58 265.69 300.17 335.60 | 1.1177e+14 | 0.2367 | -6.89e-06 |
| `photochem x0.20/HeH9.7` | 181.72 187.57 207.08 235.35 271.49 308.35 | 1.3451e+14 | 0.2172 | -1.01e-05 |
| `x0.01_kzz1e9/HeH2.13` | 222.48 222.75 228.78 241.38 263.50 287.75 | 9.2437e+13 | 0.2554 | -4.35e-06 |
| `photochem x0.01/HeH9.7` (uncertified) | 181.23 181.06 187.80 203.92 236.47 275.74 | 1.3485e+14 | 0.2520 | -9.29e-06 |

For the target case the mapped seed carried T(1) = 578.54 K and
rho(1) = 4.3859e+13 (MEASURED from the seed), so the Roe solve moved the first
cell to 221.78 K and 1.1177e+14: a base 62 per cent colder and 155 per cent
denser than the state it started from. That is the same direction and the same
size as the movement step 2 measured on this case at one outer pass (READ,
section 2.9), so the two measurements agree.

The alternating base velocity is NOT removed by the flux change: the odd-even
projection of the certified states is 0.22 to 0.26, against 0.15 in the seeds
and 0.09 to 0.11 in the certified HLLC catalog states. The Roe solutions carry
MORE of that mode than the HLLC states of the catalog do, and certify anyway.

### 3.4 Item 2: the flux family where both solve

The control is re-solved under ROE from its own certified HLLC state and
judged on its own residual (R37). Two of the four finished inside the window.
MEASURED:

| case | flux | outcome | `Mdot` [g/s] | EW [%A] | T(1) [K] | rho(1) [cm^-3] | `du` | odd-even |
|---|---|---|---|---|---|---|---|---|
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | HLLC (catalog) | certified | 7.37309e+07 | 1.3749 | 238.00 | 8.6581e+13 | 1.940e-04 | 0.0858 |
| the same | ROE | **certified**, pass 3, 2.187e-09 | 7.35009e+07 | 1.3707 | 224.44 | 9.1654e+13 | 1.945e-04 | 0.1835 |
| movement | | | **-0.31 %** | **-0.31 %** | **-5.70 %** | **+5.86 %** | +0.26 % | |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | HLLC (catalog) | certified | 6.55794e+06 | 9.1473e-05 | 239.92 | 8.5911e+13 | 2.940e-03 | 0.1077 |
| the same | ROE | **certified**, pass 15, 4.702e-08 | 6.51761e+06 | 9.1069e-05 | 222.63 | 9.2374e+13 | 2.936e-03 | 0.2387 |
| movement | | | **-0.62 %** | **-0.44 %** | **-7.21 %** | **+7.52 %** | -0.12 % | |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | ROE | still running, pass 4, 1.9e-03 | - | - | - | - | - | - |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7/k04` | ROE | still running, pass 8, 4.0e-03 | - | - | - | - | - | - |

**That is the systematic of the flux choice at the XUV levels where both
solve**: a few tenths of a per cent in the mass-loss rate and in the He 10830
equivalent width, about 6 to 8 per cent in the first cell's temperature and
density with the two of opposite sign at nearly equal magnitude (the base
pressure is held by the boundary, so a colder first cell is a denser one), and
nothing measurable in the mass-flux spread. The HLLC catalog reference values
for the two cases that did not finish are `Mdot` 2.0882e+07, EW 0.1894 (x0.30)
and `Mdot` 7.4646e+07, EW 2.2206 (photochem fiducial k04), MEASURED from the
catalog states for whoever completes those two legs.

### 3.5 Item 3: grid convergence on `x0.10_kzz1e9/HeH9.7`

The catalog grid is `Grid cells: 500` (the default; the shared grid file states
`N 500 r_min 1.0001933 r_max 29.0312`). The doubled grid was made by adding
`Grid cells: 1000` to the same input; the code's Mixed construction then gives
`N 1000 r_min 1.0001971 r_max 29.5816`, so **the doubled grid is not a
refinement of the 500-cell grid: both ends move**, by 3.8e-06 R_p at the base
and 1.9 per cent at the outer edge. The seed was mapped onto it with
`src/utils/map_state_to_grid.py --extrapolate-beyond 28.0` and the run was made
with `NOSEED=1`, because `run_case.sh` maps only onto the 500-cell grid.

MEASURED:

| flux | cells | outcome | `Mdot` [g/s] | EW [%A] | `du` | T(1) [K] | rho(1) [cm^-3] | odd-even of v(1:6) |
|---|---|---|---|---|---|---|---|---|
| ROE | 500 | **certified**, pass 14, 4.143e-07 | 6.4102e+06 | 0.2353 | 2.00e-03 | 221.78 | 1.1177e+14 | 0.2367 |
| HLLC | 500 | never solved (the catalog stall; step 2 measured 1.779 after 60 Newton iterations, READ) | - | - | - | - | - | - |
| HLLC | 1000 | **certified**, pass 11, 3.552e-07 | 6.4761e+06 | 0.2344 | 9.21e-04 | 247.70 | 1.0032e+14 | 0.2201 |
| ROE | 1000 | still running at the end of the window, 30 Newton iterations, `||R||` 1.99, state still the seed | - | - | - | - | - | 0.1326 |

**The doubled grid is what decides this.** HLLC, which has no solution of this
case at 500 cells, certifies it at 1000, and the state it reaches is the state
the Roe flux reaches at 500: `Mdot` differs by +1.03 per cent, the equivalent
width by -0.40 per cent. Profile against profile, at matched radius, MEASURED
(ROE 500 against HLLC 1000):

| r [R_p] | T ratio | rho ratio |
|---|---|---|
| 1.0005 | -21.1 % | +23.4 % |
| 1.0010 | -18.1 % | +15.3 % |
| 1.0020 | -7.9 % | -0.5 % |
| 1.0050 | -1.1 % | -9.0 % |
| 1.0100 | +1.0 % | -10.8 % |
| 1.0500 | +1.6 % | -8.6 % |
| 1.2000 | +1.2 % | -3.0 % |
| 2.0000 | -0.1 % | -1.0 % |
| 5.0000 | -0.3 % | -0.5 % |
| 10.000 | +0.1 % | -0.1 % |

So the disagreement the two fluxes have at the base is confined to the bottom
two or three cells, where the two runs also do not share a grid, and it is gone
by 1.2 R_p in temperature and by 2 R_p in density. The wind, the mass-loss rate
and the line are the same object. The base odd-even mode is 0.24 at 500 cells
and 0.22 at 1000 and does not fall with resolution in either flux, which is
what a 2 dr mode of the discretization does.

**What this does not establish.** Only one case was refined, and only by a
factor two; the ROE leg at 1000 cells had not converged when the window closed,
so the pair (ROE 500, ROE 1000) is missing and the Roe flux's own convergence
rate is not measured here. The two grids differ at both ends, so the base
comparison at r - 1 < 0.002 R_p is not a pure resolution difference.

### 3.6 Item 4: conservation and the critical point of every certified Roe state

MEASURED. `du` is the mass-flux spread over the wind window the run prints;
the two closures are the ones the certification prints. The L8 check is the
Mach number at the outermost physical cell: a value below 1 says the critical
point is outside the domain, which is where L8 put it for this planet (40.06
R_p against a 30 R_p domain, READ,
`docs/lhs1140b_stationary_L8_20260913.md`).

| state | `du` (window) | chemistry mass closure | composition closure | Mach at r_max | critical point |
|---|---|---|---|---|---|
| `x0.10/HeH9.7` ROE 500 | 2.00e-03 (309..500) | 3.882e-16 | 4.208e-16 | 0.0815 | outside the domain |
| `photochem x0.20/HeH9.7` ROE | 6.94e-04 (309..500) | 4.100e-16 | 8.038e-16 | 0.1773 | outside |
| `x0.01/HeH2.13` ROE | 4.09e-03 (309..500) | 3.598e-16 | 5.374e-16 | 0.00599 | outside |
| `photochem x0.01/HeH9.7` ROE (uncertified) | 4.03e-03 (309..500) | 2.788e-16 | 7.662e-16 | 0.00649 | outside |
| `fiducial/HeH2.13` ROE | 1.945e-04 (309..500) | 4.738e-16 | 4.880e-16 | 0.6541 | outside |
| `x0.10/HeH2.13` ROE | 2.936e-03 (309..500) | 3.262e-16 | 4.355e-16 | 0.0682 | outside |
| `x0.10/HeH9.7` HLLC 1000 | 9.21e-04 (546..1000) | 3.824e-16 | 4.795e-16 | 0.0844 | outside |

Every closure is at the rounding of the element budget, as it is for the
catalog states. No state reaches its sound speed inside the domain, so every
one of them is the subsonic part of a transonic escape whose critical point
sits outside a 30 R_p domain, which is the L8 situation and carries the L8
conditions on what may be quoted from such a run.

### 3.7 A neutral reading

- **What the Roe results say about the seven.** The flux family is what
  separates a solved case from a stalled one at these XUV levels: three of the
  seven certify under ROE with nothing else changed, a fourth reaches a
  stationary wind and is held only by the composition in the outer wind, and
  three had not finished when the window closed. That does not make the Roe
  flux right; it makes the HLLC stall a property of that flux on that grid, and
  section 3.5 shows it is a resolution property, because HLLC solves the same
  case on twice the cells.
- **What the systematic is where both solve.** A few tenths of a per cent in
  `Mdot` and in the He 10830 equivalent width, and 6 to 8 per cent in the
  first cell. R9 is confirmed again: a base-local change moves `Mdot`, and it
  moves it by far less than the base itself moves.
- **Whether resolution decides which base is right.** For this one case, yes in
  the wind and no at the base. The refined HLLC solution and the coarse Roe
  solution are the same wind to 1 per cent in `Mdot` and 0.4 per cent in the
  line, so neither base is a different physical state; the 10 to 20 per cent
  difference that remains is in the bottom two cells of grids whose first
  centers are not the same, and deciding it needs a refinement of one flux at
  fixed domain, which was not run.
- **What adopting ROE for the catalog would require.** Every certified case of
  `models/` re-solved under it and judged on its own residual, conservation and
  a grid check, not compared against the HLLC state it replaces (R37); the
  transit products and the memo figures rebuilt from the new states; the
  regression goldens refreshed deliberately at the end of the series with the
  refresh reported, and the byte-identity of the HLLC-selecting cases measured
  first, since `Numerical flux:` is an input key and the matrix cases select
  HLLC; and the numbers of this memo re-measured on whatever binary the
  adoption is made with. None of that was done here, and no production design
  is chosen by this step.

### 3.8 What was not measured, and one defect found on the way

- The ROE leg of the grid study at 1000 cells, and any refinement beyond a
  factor two; the three cases of the seven and the two controls that were still
  running when the window closed (section 3.9).
- A case whose two grids share both domain ends, which is what a clean
  resolution statement about the base needs; the code's Mixed grid moves
  `r_min` and `r_max` with the cell count, and there is no key to hold them.
- **A defect outside this item, found and NOT fixed** (it is not in a file this
  item owns): with the tree binary `c2e9c9990b9f`, the post-processing pass that
  `run_case.sh` prescribes since item L9, `Restart intent: stationary evaluate`,
  **writes no `Hydro_ioniz_adv.txt` / `Ion_species_adv.txt`**, so the transit
  synthesis that follows it exits 1 with `FileNotFoundError: ./output/Hydro_ioniz_adv.txt`.
  MEASURED on every case of this series that reached the products. The current
  `src/EXHALE_main.f90` does call `post_process_adv` on that route (READ, lines
  5617-5652), and the message the binary prints on that route is not in the
  current source, so the tree binary predates that repair; the catalog's own
  `pp.log` files state `Restart intent: relaxation` (READ), meaning the catalog
  products were made with the older one-step recipe at CFL 1e-12. This series
  therefore made its products with that older recipe, which is the one the
  catalog numbers it is compared against were made with, through
  `models/.L25/post.sh`. Anyone re-running `run_case.sh` end to end on this
  binary will hit the same failure after a successful solve.

### 3.9 The runs still going when this was written (02:00 KST, 2026-09-17)

On `lart3`, all under `LHS1140b/models/.L25/`:

| pid | directory | state |
|---|---|---|
| 6775 | `roe_p015_HeH9.7` | pass 6, `||R||` 1.29 |
| 8999 | `roe_p025_HeH9.7` | pass 9, `||R||` 1.67 |
| 8250 | `roe_s001_HeH9.7` | pass 10, `||R||` 1.94 |
| 10606 | `roe_s010_HeH9.7_d1` | pass 8, `||R||` 3.6e-03 (the target at `dtau0` 1.0 from the stopped attempt's state, the configuration step 2 measured) |
| 15057 | `roe_s030_HeH2.13` | pass 4, `||R||` 1.9e-03 |
| 18290 | `roe_pfid_HeH9.7` | pass 8, `||R||` 5.8e-03 |
| 17371 | `g1000_roe` | 30 Newton iterations, `||R||` 1.99 |
| 24423, 28424 | `g500_hllc` | restarted, the HLLC 500-cell leg of the grid study |

Each will write its own `run.log`, and the solved ones need
`models/.L25/post.sh <dir>` for their products, for the reason in section 3.8.

## Step 4: continuation in XUV for the helium-rich cases

Item L25 step 4, the user's instruction of 2026-09-17: for the low-XUV atomic
cases that step 3 left unsolved, start from a SOLVED case at a higher XUV with
the same composition and lower the XUV in small steps, each rung seeded from
the previous certified rung. Written by the step 4 worker; steps 1 to 3 above
are other workers' and are untouched. **No source was changed and nothing was
built.** Every run used the tree binary `EXHALE.x`, md5
`c2e9c9990b9f14f1be8cd77abca68945`, verified on both hosts before the series.
Every number below is MEASURED on the host `lart4` on 2026-09-17 at 8 threads
unless it is marked READ. The run directories are
`LHS1140b/models/.L25/ladder/<ladder>/x<tag>/`; no catalog directory and no
golden was touched.

### 4.1 What a rung is

A rung is the catalog `input.inp` of the ladder's TARGET case with ONLY
`Spectrum file:` changed to that rung's scaled spectrum. So the flux (`ROE`,
the user's decision of 2026-09-17), the base boundary, `Resid tol`, the
`He/H number ratio` and the lower-atmosphere profile are the target's at every
rung and the XUV is the only thing that moves. VERIFIED by diff: each rung's
`input.inp` against its target's differs in that one line and nothing else.
The lower-atmosphere profile of the photochem group is the SAME file at every
XUV level (md5 `3f96da22717e79667a5c93e932ded7f3` at 0.01, 0.10, 0.15, 0.20,
0.25 and 0.30, MEASURED), so carrying the target's copy onto a rung changes no
boundary.

The recipe is the catalog's low-XUV one: `EXHALE_PTC_DTAU0=1.0e8` (item L14,
`models/run_lowxuv_dtau0.sh`), `EXHALE_OUTER_PASSES=40`, `SEED_ATTEMPTS=1`,
`FORCE=1`, `OMP_NUM_THREADS=8`, through `models/run_case.sh`. The seed of the
first rung is the certified state named in the table below; the seed of every
later rung is the `output/` of the rung above it, mapped onto this case's grid
by `src/utils/map_state_to_grid.py` the way `run_case.sh` maps any seed. A rung
that does not certify STOPS its ladder, and the driver
(`models/.L25/ladder/run_ladder.sh`) records the host, the wall time, the
`info` and the certification verdict of every rung in
`models/.L25/ladder/<ladder>.ladder.log`.

### 4.2 The spectra

Four scalings did not exist and were written into `LHS1140b/sed/`: 0.18, 0.16,
0.28 and 0.26 of the fiducial GJ 1132 spectrum at the b orbit. They were made
by the rule the existing grid was made by (L14 section 3.3): the whole flux
column of `lhs1140_sed_gj1132_at_b.txt` times the scale, the wavelength column
untouched, written `%.6e`, the base file's header kept with one line added
stating the scale. CHECK, MEASURED: regenerating the EXISTING `0.20`, `0.15`,
`0.10`, `0.03` and `0.015` files by that same routine reproduces all 29992 data
rows of each of them with no mismatch, so the four new files are on the grid the
old ones are on. The new files' own flux columns are the base column times the
scale to 5.0e-07 relative, which is the `%.6e` rounding of the format.

### 4.3 The definitions the tables use

MEASURED from `output/Hydro_ioniz.txt`, physical cells only (rows 3 to 502 of
the file):

- `Mdot` is `4 pi rho v r^2` of the outermost physical cell, with `rho` the
  mass density (the file's `rho[mH/cm3]` column times the hydrogen mass) and
  the radial scale `R0` the state's own header states. NB for the photochem
  group that is the profile's matching-level radius, not the planet radius; the
  step 3 tables of this memo used the planet radius there, which is why the
  photochem `Mdot` of the same state is 5.8 per cent lower in section 3.3 than
  it is here. This section's value is the one the binary's own
  `Log10 of steady-state Mdot` agrees with (MEASURED: photochem `x0.20`,
  1.4347e+07 here against `10^7.16` = 1.45e+07 from `pp.log`).
- T(1..3) and rho(1..3) are the first three physical cells.
- the odd-even entry is `|sum_i (-1)^i v_i| / sum_i |v_i|` over physical cells
  1 to 6, which is not the normalization section 3.3 used, so the two columns
  are not to be compared number for number.
- "ACCEPTED at pass N" is the outer pass the run log declares
  `ACCEPTED -- every active equation of this state is within its own tolerance`
  on.

Every rung certifies IN THE WIND (`cert_reason=certified_in_wind`, r >= 1.20
R_p), which is what every certified Roe state of step 3 does as well.

### 4.4 The products these runs do NOT carry

The transit synthesis of `run_case.sh` fails on every rung, for the reason
already stated in section 3.8: on this binary the `Restart intent: stationary
evaluate` route writes no `*_adv.txt`, so `EXHALE_transit.py` exits 1 with
`FileNotFoundError: ./output/Hydro_ioniz_adv.txt`. The SOLVE is unaffected
(`info = 0`, certified, the state written), and the ladder is gated on the
solve's verdict and not on the products. So these rungs carry a certified state
and no He 10830 equivalent width. Making the line of a rung needs
`models/.L25/post.sh <rung>`, which takes a step at CFL 1e-12 and therefore
writes a different state back; it was NOT run on any rung IN PLACE, so that the
state a later rung was seeded from, and the state named as a campaign seed
below, is the solved one. The line of the four TARGET rungs was made afterwards
on COPIES of them and is in section 4.10; the intermediate rungs have none and
need none.

### 4.5 Verdict

**All four targets certify by continuation in XUV, and every one of the
eighteen rungs certified: not one ladder was stopped.** The four cases step 3
left open are open no longer, and none of them needed a source change, a
different flux, a different pseudo-time start or a larger pass budget. What
they needed was a nearer seed.

| target | ladder | rungs | outcome | the rung the campaign seeds from |
|---|---|---|---|---|
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | L1 | 0.18, 0.16, 0.15 | **CERTIFIED**, pass 7 | `models/.L25/ladder/L1/x0p15/output` |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | L2 | 0.28, 0.26, 0.25 | **CERTIFIED**, pass 5 | `models/.L25/ladder/L2/x0p25/output` |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | L3 | 0.07, 0.05, 0.03, 0.02, 0.015, 0.01 | **CERTIFIED**, pass 23 | `models/.L25/ladder/L3/x0p01/output` |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | L4 | 0.07, 0.05, 0.03, 0.02, 0.015, 0.01 | **CERTIFIED**, pass 24 | `models/.L25/ladder/L4/x0p01/output` |

The `not_solved.md` of each of the four catalog directories was removed and
each carries a one-line `seed_from_ladder.txt` naming the state above. Nothing
else was written into a catalog directory.

### 4.6 The rung tables


#### ladder L1
| rung | outcome | ACCEPTED at pass | `Mdot` [g/s] | T(1..3) [K] | rho(1..3) [cm^-3] | odd-even of v(1:6) |
|---|---|---|---|---|---|---|
| x0.18 | CERTIFIED | 7 | 1.2819e+07 | 181.63 187.14 205.93 | 1.3457e+14 1.1931e+14 9.9572e+13 | 0.2801 |
| x0.16 | CERTIFIED | 8 | 1.1299e+07 | 181.55 186.68 204.69 | 1.3463e+14 1.1959e+14 1.0013e+14 | 0.2830 |
| x0.15 | CERTIFIED | 7 | 1.0543e+07 | 181.50 186.44 204.02 | 1.3466e+14 1.1974e+14 1.0043e+14 | 0.2845 |

#### ladder L2
| rung | outcome | ACCEPTED at pass | `Mdot` [g/s] | T(1..3) [K] | rho(1..3) [cm^-3] | odd-even of v(1:6) |
|---|---|---|---|---|---|---|
| x0.28 | CERTIFIED | 6 | 2.0519e+07 | 181.99 188.98 210.84 | 1.3431e+14 1.1822e+14 9.7444e+13 | 0.2683 |
| x0.26 | CERTIFIED | 6 | 1.8968e+07 | 181.93 188.66 210.00 | 1.3436e+14 1.1841e+14 9.7801e+13 | 0.2705 |
| x0.25 | CERTIFIED | 5 | 1.8194e+07 | 181.89 188.50 209.56 | 1.3438e+14 1.1851e+14 9.7990e+13 | 0.2716 |

#### ladder L3
| rung | outcome | ACCEPTED at pass | `Mdot` [g/s] | T(1..3) [K] | rho(1..3) [cm^-3] | odd-even of v(1:6) |
|---|---|---|---|---|---|---|
| x0.07 | CERTIFIED | 17 | 4.3861e+06 | 221.56 224.87 238.67 | 1.1187e+14 1.0208e+14 8.9312e+13 | 0.3082 |
| x0.05 | CERTIFIED | 26 | 3.2042e+06 | 221.37 223.81 235.80 | 1.1197e+14 1.0254e+14 9.0320e+13 | 0.3130 |
| x0.03 | CERTIFIED | 24 | 1.8810e+06 | 221.31 222.82 232.96 | 1.1200e+14 1.0297e+14 9.1349e+13 | 0.3184 |
| x0.02 | CERTIFIED | 23 | 1.2315e+06 | 221.45 222.36 231.27 | 1.1193e+14 1.0319e+14 9.1983e+13 | 0.3218 |
| x0.015 | CERTIFIED | 23 | 9.1318e+05 | 221.55 222.07 230.17 | 1.1188e+14 1.0332e+14 9.2403e+13 | 0.3239 |
| x0.01 | CERTIFIED | 23 | 6.0016e+05 | 221.65 221.66 228.66 | 1.1183e+14 1.0351e+14 9.2979e+13 | 0.3266 |

#### ladder L4
| rung | outcome | ACCEPTED at pass | `Mdot` [g/s] | T(1..3) [K] | rho(1..3) [cm^-3] | odd-even of v(1:6) |
|---|---|---|---|---|---|---|
| x0.07 | CERTIFIED | 15 | 4.6022e+06 | 181.04 183.96 197.34 | 1.3499e+14 1.2125e+14 1.0354e+14 | 0.2993 |
| x0.05 | CERTIFIED | 27 | 3.2593e+06 | 180.87 183.05 194.85 | 1.3511e+14 1.2181e+14 1.0475e+14 | 0.3045 |
| x0.03 | CERTIFIED | 25 | 1.8980e+06 | 180.90 182.19 192.11 | 1.3509e+14 1.2237e+14 1.0614e+14 | 0.3109 |
| x0.02 | CERTIFIED | 24 | 1.2311e+06 | 181.05 181.75 190.40 | 1.3498e+14 1.2266e+14 1.0704e+14 | 0.3150 |
| x0.015 | CERTIFIED | 23 | 9.0372e+05 | 181.15 181.47 189.30 | 1.3491e+14 1.2285e+14 1.0763e+14 | 0.3175 |
| x0.01 | CERTIFIED | 24 | 5.8235e+05 | 181.23 181.06 187.80 | 1.3485e+14 1.2312e+14 1.0843e+14 | 0.3207 |

### 4.7 What the continuation changed, and what it did not

**It did not change the wind.** The clearest statement of that is the pair at
0.01 of the photochem column. Step 3's direct solve of that case
(`models/.L25/roe_p001_HeH9.7`) spent its forty passes and was refused by the
gated elemental He/H partition row at 2.234e-04 against 1.0e-05 at cell 440,
with its three hydrodynamic rows all inside their tolerances. The final rung of
ladder L4 is the SAME state: cell by cell over all 500 physical cells, MEASURED,
the two agree to 1.9e-07 in mass density, 1.8e-07 in velocity and 4.4e-07 in
temperature, on identical cell centers. What moved is the composition's
distance to its own fixed point: the same gated row reads 7.20e-06 in the rung
and is accepted at outer pass 24.

So the reading of these four cases is not that the low-XUV wind was hard to
find. It is that the ELEMENTAL PARTITION relaxes toward its fixed point by a
fixed factor per outer pass, and a seed one large XUV step away enters so far
out that forty passes are not enough decades. Shortening the step is what the
budget needed. Ladder L1 crossed 0.20 to 0.15 in three steps and each rung cost
5 to 8 passes; ladders L3 and L4 crossed 0.10 to 0.01 in six steps and each
rung cost 15 to 27. Both are inside the forty the catalog allows, while the single step
from 0.10 to 0.01 was not.

**The gated row rises before it falls on every rung.** MEASURED on
`L3/x0p05`, every pass: 1.25e-02 at the first, a maximum of 4.36e-02 at the
sixth, then a fall that steepens until it settles at almost exactly a factor
0.5 a pass from the eighteenth on (1.60e-03, 8.18e-04, 4.14e-04, 2.09e-04,
1.05e-04, 5.26e-05, 2.63e-05, 1.32e-05, 6.60e-06 accepted at the
twenty-sixth). That is the behaviour item L14 recorded, and it is why a rung
must not be judged on its early passes.

**The base keeps its alternating velocity, and slightly more of it at lower
XUV.** The odd-even projection of the first six cells rises monotonically down
every ladder: 0.2683 to 0.2716 over L2, 0.2801 to 0.2845 over L1, 0.3082 to
0.3266 over L3 and 0.2993 to 0.3207 over L4. So the mode the L17 diagnosis
identified at the low-Mach base does not go away with the XUV and does not stop
these cases from certifying. It is a 2 dr mode of the discretization, as
section 3.5 found it to be under grid refinement.

**The profiles down each ladder are smooth and monotone**, which is what a
continuation in a single parameter should give and is a check on the states
themselves. `Mdot` falls monotonically with the XUV on all four ladders, and
over the scalar ladder it falls by a factor 7.31 while the XUV falls by a factor
7.0 (from 4.3861e+06 at 0.07 to 6.0016e+05 at 0.01), so the mass-loss rate is
very nearly linear in the XUV at this end, MEASURED, which is the
energy-limited behaviour and not a signature of the solver. Over that same
factor 7 in XUV the FIRST cell's temperature moves by 0.09 K on the scalar
ladder and 0.19 K on the photochem one, while the THIRD cell's falls by 10.0
and 9.5 K, so what the XUV moves at the base is the gradient and not the
boundary value the pressure condition fixes.

### 4.8 What this does NOT establish

- Nothing here is a statement about the flux family. Every rung was run with
  `Numerical flux: ROE`, which is the user's decision of 2026-09-17 for these
  cases, and no HLLC continuation was run, so it is NOT shown that a
  continuation would rescue these cases under HLLC as well. Whether the
  continuation or the flux is what closes `x0.01_kzz1e9/HeH9.7`, which item L17
  diagnosed as needing a low-Mach flux, is not separated by this series.
- The rungs carry no He 10830 equivalent width, for the reason of section 4.4,
  so the line of these four cases is still to be made, in the campaign re-run
  or by `models/.L25/post.sh`.
- The certification of every rung is IN THE WIND (r >= 1.20 R_p); the element
  and carrier rows below that radius are reported by the run and do not gate,
  exactly as in the certified states of step 3 and of the catalog.
- No intermediate rung is a catalog case and none is claimed as a physical
  result: the XUV scalings 0.18, 0.16, 0.28 and 0.26 exist only to shorten the
  step.

### 4.9 Where the runs ran, and what is left running

Every rung of every ladder ran on `lart4` (72 cores, the same NFS tree, the
same binary, md5 verified there before the series), at 8 threads, at most four
jobs at once, by the user's instruction of 2026-09-17 to move off `lart3`.
Nothing was started on `lart3`. **Nothing of this series is still running:** at
06:12 KST on 2026-09-17 no `EXHALE.x` and no ladder driver of this series
remained on either host, MEASURED by `pgrep` with each pid's `/proc/<pid>/cwd`
checked. The record of every rung is in
`models/.L25/ladder/L{1,2,3,4}.ladder.log`, one line a rung with the host, the
wall time, the `info` and the verdict.

Wall times, MEASURED: L1 and L2 took 231 to 354 seconds a rung, L3 and L4 took
540 to 1285. The whole series took 1 hour 38 minutes from the start of the
first rung (04:33:46) to the end of the last (06:12:17).

### 4.10 The He 10830 line of the four targets

Asked for after the ladders closed, and made the way step 3 made its lines so
that the numbers sit beside the catalog's: each of the four TARGET rungs was
COPIED to a sibling `<rung>_pp/` and `models/.L25/post.sh` was run on the COPY,
then `EXHALE_transit.py` through the WINERED HIRES-Y kernel
(`EXHALE_TRANSIT_RES_HETR=68000`, `MPLBACKEND=Agg`), on `lart4` at 8 threads.
**The recipe is the OLDER one**, `Do only PP: True` with `CFL: 1.0e-12`, which
is what the catalog's own `pp.log` files record (READ) and what writes the
`*_adv.txt` profiles the transit tool consumes; it takes one small step and
writes a different state back, which is why it was run on a copy and not on the
rung. VERIFIED after the fact: the four solved rung directories carry no
`*_adv.txt` and no `tpm_*` and their `output/Hydro_ioniz.txt` still has the
mtime of the ladder run, so the states named in section 4.5 and in each
`seed_from_ladder.txt` are unchanged.

`EW` is the red-pair equivalent width over the measurement's vacuum window
10832.60 to 10834.20 A, `red_depth` and `fwhm_A` are the three-Gaussian fit of
`tpm_He10830_metrics.txt`, all read exactly as `models/run_case.sh` reads them.
MEASURED:

| XUV | case | state | `EW` [%A] | `red_depth` [%] | `fwhm_A` [A] | red/blue |
|---|---|---|---|---|---|---|
| 0.30 | photochem HeH9.7 | catalog, HLLC | 1.0097 | 3.5415 | 0.2664 | 7.075 |
| **0.25** | **photochem HeH9.7** | **L2 target rung** | **0.84287** | **3.0020** | **0.2623** | **7.088** |
| 0.20 | photochem HeH9.7 | step 3, `roe_p020` | 0.67112 | 2.4326 | 0.2578 | 7.101 |
| **0.15** | **photochem HeH9.7** | **L1 target rung** | **0.48218** | **1.7829** | **0.2527** | **7.117** |
| 0.10 | photochem HeH9.7 | catalog, HLLC | 0.26737 | 1.0007 | 0.2487 | 7.152 |
| **0.01** | **photochem HeH9.7** | **L4 target rung** | **3.3576e-05** | **0.0001** | **0.2987** | **7.490** |
| 0.10 | scalar HeH9.7 | step 3, `roe_s010` | 0.23532 | 0.8787 | 0.2487 | 7.166 |
| **0.01** | **scalar HeH9.7** | **L3 target rung** | **3.2456e-05** | **0.0001** | **0.3461** | **5.477** |

**The line falls monotonically with the XUV and is gone by 0.01.** Down the
photochem column the equivalent width falls by a factor 3.8 from 0.30 to 0.10
while the XUV falls by 3, so it is slightly steeper than linear there, and then
by a further factor 7965 from 0.10 to 0.01, four orders of magnitude for one
order in XUV. At 0.01 of the fiducial spectrum the He I 10830 line is not a
detection of anything: the red depth is 1e-4 per cent. That is the same
statement the two `not_solved.md` of the 0.01 cases already made from the
ladder of item L17 (READ), now made on solved states of the cases themselves.

**The three new lines join their neighbours smoothly**, which is a check on the
rung states independent of the residual: `EW`, `red_depth` and `fwhm_A` are each
monotone across the six photochem entries with the ladder rungs interleaved
between catalog and step 3 states, and the red-to-blue ratio drifts slowly from
7.075 to 7.152 over the range where the line is measurable.

**The 0.01 pair is a consistency check that came out exact.** The L4 target rung
and step 3's uncertified direct solve of the same case give `EW` 3.3576e-05,
`red_depth` 0.0001, `fwhm_A` 0.2987 and red/blue 7.490, agreeing in every digit
printed. That is what section 4.7 predicted from the two states agreeing to
4.4e-07, and it confirms that the certification the continuation bought changed
the composition's distance to its fixed point and not the observable.

Caveats. The two catalog entries at 0.30 and 0.10 are HLLC states while every
other row is a Roe state, so the flux-family systematic of section 3.4 (the
equivalent width lower by 0.3 to 0.4 per cent under ROE) sits between them and
is smaller than any step in this table. The `fwhm_A` of the two 0.01 rows is a
fit to a line of depth 1e-4 per cent and is not meaningful. The products live
in `models/.L25/ladder/L{1,2,3,4}/x<tag>_pp/`, not in the rung directories and
not in any catalog directory.
