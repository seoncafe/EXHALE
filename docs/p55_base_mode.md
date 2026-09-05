# P55 -- Why the flat steady root is not a fixed point of the marching loop

**Short answer: the marching loop does not have that state as a fixed point,
and it never did.** The departure is not a growing mode released by round-off:
it starts at step 1 at a finite size and grows *linearly* in step number for the
first twenty steps, at a rate the accepted state's own residual predicts. Two
independent things produce it, both measured here:

1. **The state the JFNK solve accepts is a zero of the residual only in the
   rows and cells the acceptance norm can see.** Judged by its own largest term
   -- the P54b construction applied to rows 2 and 3 instead of row 1 -- the
   accepted hot-Uranus state's energy row is out by `1.2e-3` at `r = 1.03` and
   `6.0e-3` at `r = 1.05`, and its momentum row by `2.5e-5` at `r = 1.03`;
   cell 1 is out by `1.4e-1` (momentum) and `9.9e-1` (energy). The gate reads
   `3.0e-9`, `-3.6e-7` and `7.8e-4 / -1.5e-4` for the same cells. WASP-121b is
   the same picture with different numbers.
2. **The marching contains one operator-split stage the steady residual has no
   counterpart for** -- the pressure re-derived from `T` at the post-sweep
   particle count (`comp_p_from_T` after `ioniz_eq`). On a `Load IC` restart
   this is a **`dt`-INDEPENDENT jump**: the same `dE` for CFL `0.6`, `0.06`,
   `0.006` and `0.0006`, reaching `6.2e-3` of `E` in one step at `r = 1.10` on
   WASP-121b, twelve times that step's whole hydro term. It decays by three
   decades over the first ten steps and is `O(dt)` thereafter.

Neither is a discrete instability of the scheme, and the variant scan says so
directly: **halving and quartering the time step changes the excursion by 3
percent at equal model time**, dropping the base cells to first-order
reconstruction changes nothing on WASP-121b and makes the hot Uranus 12 times
worse, the periodic Shapiro filter makes the hot Uranus **42 times worse**, and
freezing the first 20 cells outright does not help WASP-121b at all. The one
thing that does help is starting from a state the marching itself produced:
3 times smaller on WASP-121b, 41 times on the hot Uranus.

Figures: `p55_base_mode.png` -- (a) the linear-in-step departure, (b) the
3000-step series, (c) the variant scan. `p55_operator_terms.png` -- (a) the
accepted state judged by its own row terms against what the gate reads, (b) why
a `1e-5` momentum residual moves the mass flux, (c) the restart kick.

---

## 1. What was measured, on what, with what

**States.** The two `EXHALE_CONT_SCALE=gate` solutions of P54b section 10.4,
copied read-only from `.../scratchpad/p54flux/n1G` and `.../w1G`:

| | hot Uranus G | WASP-121b G |
|---|---|---|
| JFNK | `info = 0`, 9 iterations | `info = 0`, 17 iterations |
| `\|\|R\|\|` the solve reported | 8.636e-06 | 3.138e-06 |
| face-flux spread `r >= 1.01 / 1.03 / 1.20` | 1.043e-04 / 1.043e-04 / 3.778e-06 | 3.133e-05 / 2.325e-05 / 3.990e-06 |
| `log10 Mdot` | 10.34 | 13.21 |
| base BC | `rho = rho_bc`, hard one-way valve `max(v_1,0)`, isothermal ghost | same |
| Shapiro / low-Mach damping / transport / diffusion / carrier | all off | all off |

The two spread columns reproduce P54b section 10.4 digit for digit on reload,
so the hydro state the runs below start from is the state the solver accepted.
The residual does not reproduce as tightly (`1.385e-05` volume-weighted
re-measured against `8.636e-06` accepted on the hot Uranus, `8.457e-06` against
`3.138e-06` on WASP-121b): the composition in the restart file is written at
finite precision and is not the equilibrium composition of the state, which is
P54's open item and which section 4 below turns into a number.

**Binary.** `.../scratchpad/p54flux/src` + `Makefile` + `build/` copied on
2026-09-03, i.e. tree `ea9e7492` plus `p54_probes.diff` plus `p54b_scale.diff`.
Three files changed, all env-gated, all off unless their variable is set
(`p55_probes.diff`):

| gate | what it does |
|---|---|
| `EXHALE_P55_DU=<n>` | `output/p55_dU.txt`: for the first `n` steps, per cell, the increment each operator-split stage of the step adds to each conserved variable (hydro RK3, composition pressure projection, semi-implicit energy, everything after it) |
| `EXHALE_P55_TS=<k>` | `output/p55_ts.txt`: the P54 face-flux series at `r = 1.005, 1.03, 1.10, 1.20` with the base-cell state beside it, every `k` steps and every step for the first 20 |
| `EXHALE_P55_MOM=1` | `output/p55_row_terms.txt` (inside the `EXHALE_RESIDUAL` branch): `dF_k`, `S_k`, `heat`, `cool`, `R_k` and `s_k` per row per cell, so a residual can be read as a fraction of its own row's largest term |
| `EXHALE_P55_NOCHEM=1` | skips `comp_p_from_T` after the ionization sweep, so the marching keeps the pressure the hydro stage produced |
| `EXHALE_P55_FLATREC=<m>` | piecewise-constant face states in the first `m` cells |
| `EXHALE_P55_FIXBASE=<m>` | holds the first `m` physical cells at the state the run started from |

**Default off, checked.** With every `EXHALE_P55_*` variable unset, a 25-step
hot-Uranus run of this build and of the P54b build it was copied from produce
byte-identical `output/Hydro_ioniz.txt` and `output/Ion_species.txt`
(`md5 025132cf...`, `d7644001...`; `runs/nullchk_p54`, `runs/nullchk_p55`), and
the face-flux spreads of section 1 reproduce P54b section 10.4 digit for digit.

Runs: `OMP_NUM_THREADS=4`, `Reconstruction scheme: WENO3` with
`du_th [PLM,WENO3]: 1.0e-30` so that the marching and the residual use the same
reconstruction, `Solver:` deleted so nothing can hand back to the JFNK, and
`EXHALE_MAXSTEPS` as the only stop.

---

## 2. Are the steady residual and the marching RHS the same discrete operator?

**Read in the code first.** `assemble_residual` (`steady_residual.f90`) calls
`Reconstruct` and `RK_rhs` -- the same two routines the marching loop calls --
so the spatial operator is identical by construction, not by measurement. The
differences are all in what surrounds them:

| candidate difference | status in these two configurations |
|---|---|
| ghost fill order / 1-sweep lag of `n_part_cell1` | **not exercised.** `n_part_cell1` is read only by the continuous-temperature base ghost (`base_ghost_T_continuous`); both runs use the isothermal ghost. `eval_residual` refreshes it before `Apply_BC` anyway |
| `base_flux_const`, the 100-step exponential moving average behind `Base velocity: massflux` | **not exercised.** Both runs use the hard one-way valve |
| position of the ionization sweep | `ioniz_eq` takes `T` as `intent(in)`; it does not move the temperature. The sweep itself is therefore not a difference |
| **the pressure re-derived at the post-sweep particle count** | **a real difference.** The marching does `comp_T_from_p` (pre-sweep composition) -> `ioniz_eq` -> `comp_p_from_T` (post-sweep composition) -> `W_to_U`, so `u(3,:)` is rewritten. `eval_residual` never re-derives `p`. Section 4 |
| semi-implicit energy step | consistent: it solves `e(T_new) - e(T_old) = dt (heat - cool(T_new))` and writes `u(3) = rho v^2/2 + e(p_new)`, i.e. it adds `+dt (heat - cool)`, which cancels the hydro stage's `-dt (dF_3 - S_3)` exactly when `R_3 = 0` |
| the one PLM step at the head of a `PLM+WENO3` run | real but **measured to be irrelevant** (variant G, section 5) |
| Shapiro, low-Mach damping, viscosity, conduction, element diffusion, carrier transport | all off; their calls are guarded |

**Measured, cell by cell, on one step from the root.** The increments of the
four stages of one marching step, in code units of `1/t_s`, largest over the
column, at step 1 from each accepted state:

| | hydro (RK3) | chem projection | semi-implicit energy | after it |
|---|---|---|---|---|
| **hot Uranus** mass | 9.164e-01 | 3.04e-12 | 0 | 0 |
| momentum | 2.550e+00 | 3.71e-16 | 0 | 0 |
| energy | 1.359e+00 | **5.705e-03** | 5.092e-03 | 0 |
| **WASP-121b** mass | 8.338e-01 | 7.66e-12 | 0 | 0 |
| momentum | 1.561e+00 | 5.13e-14 | 0 | 0 |
| energy | 1.821e+00 | **3.868e+00** | 1.058e+00 | 0 |

The chem projection touches only the energy row -- as it must, since it rewrites
`p` at fixed `rho` and `v` -- and on WASP-121b it is the largest single term of
the step, twice the hydro RK3 and four times the heating.

---

## 3. The `dt` scan: which terms are rates and which are jumps

The same first step taken at CFL `0.6`, `0.06`, `0.006`, `0.0006`
(`dt = 1.014e-4 ... 1.014e-7` on WASP-121b). A term of the *equation* has
`dU/dt` independent of `dt`; a *projection* has `dU` independent of `dt`, so its
`dU/dt` grows as `1/dt`.

WASP-121b, `dU/dt` in code units, at four radii:

| r | row | CFL 0.6 | CFL 0.06 | CFL 0.006 | CFL 0.0006 | `-R` |
|---|---|---|---|---|---|---|
| 1.00509 | mass | 1.977e-06 | 1.306e-06 | 1.300e-06 | 1.294e-06 | 1.300e-06 |
| 1.03004 | mass | 1.354e-06 | 1.429e-06 | 1.430e-06 | 1.431e-06 | 1.430e-06 |
| 1.10007 | mass | 8.671e-08 | 1.440e-07 | 1.271e-07 | 1.254e-07 | 1.250e-07 |
| 1.10007 | **energy** | -3.294e+00 | -3.795e+01 | -3.846e+02 | -3.851e+03 | -1.571e-04 |
| 1.10007 | of which chem | -3.851e+00 | -3.851e+01 | -3.851e+02 | -3.851e+03 | -- |

**The mass row converges to the steady residual to three digits.** The energy
row does not converge at all: it is `1/dt` to four decades, and the whole of it
is the chem projection, whose `dU` is the same number `-3.905e-4` at every time
step. That is a **finite jump taken once per step, of a size that has nothing to
do with the equation being solved**, and it is `-6.2e-3` of the cell's energy.

**It is a restart artifact, not a permanent term.** The same measurement at the
39th step, when the composition the sweep is handed is the sweep's own output
rather than the file's:

| r | CFL 0.6 | CFL 0.06 | CFL 0.006 |
|---|---|---|---|
| 1.005 | 6.21e-10 | 5.62e-11 | 1.89e-14 |
| 1.03 | 2.07e-08 | 2.09e-09 | 1.40e-13 |
| 1.10 | 5.59e-06 | 6.65e-07 | -2.98e-11 |

(`dE_chem/E` per step; each decade of CFL divides it by about ten). So it scales
with `dt` once the composition and the state agree, i.e. it is then a legitimate
rate. **The kick exists only because a `Load IC` restart hands the hydro step a
composition that is not the equilibrium composition of the state it belongs
to.**

---

## 4. Where the departure actually comes from: the accepted state's own residual

`EXHALE_P55_MOM=1` writes both terms of each row, so the residual can be divided
by the row's own largest term instead of by `rho (|v| + c_s)/dr`. This is the
P54b construction, applied to rows 2 and 3.

**Hot Uranus G:**

| r | `R_2/s_2` (the gate) | `R_2/max(\|dF_2\|,\|S_2\|)` | `R_3/s_3` (the gate) | `R_3/`own |
|---|---|---|---|---|
| 1.00019 (cell 1) | 7.80e-04 | **1.43e-01** | -1.52e-04 | **-9.94e-01** |
| 1.00058 (cell 3) | -6.3e-09 | -1.1e-06 | -1.2e-08 | -1.16e-02 |
| 1.0100 | -3.6e-10 | -4.8e-08 | -3.5e-10 | -8.2e-04 |
| 1.0302 | -3.63e-07 | **-2.48e-05** | 2.98e-09 | **1.23e-03** |
| 1.0501 | 1.96e-07 | 1.05e-05 | 4.00e-08 | **6.00e-03** |
| 1.1499 | 2.67e-09 | 1.27e-07 | -1.86e-06 | -2.37e-03 |
| 1.2001 | 5.95e-10 | 3.62e-08 | -3.12e-06 | -2.12e-03 |

**WASP-121b G:**

| r | `R_2/s_2` | `R_2/`own | `R_3/s_3` | `R_3/`own |
|---|---|---|---|---|
| 1.00020 (cell 1) | 1.02e-04 | **1.87e-02** | -2.19e-05 | **-7.75e-01** |
| 1.0100 | 2.12e-11 | 3.92e-09 | -4.18e-10 | -9.91e-06 |
| 1.1001 | -1.56e-10 | -3.13e-08 | 9.69e-07 | 1.68e-04 |
| 1.1995 | -8.43e-10 | -3.29e-07 | 1.10e-04 | **1.71e-03** |

The gate's tolerance is `1e-5`. **In its own units every one of these states
passes; in the units of the physics each row is balancing, the energy row of
the layer is out by one part in a thousand and cell 1 is out by all of it.**

**Why one part in `1e5` of the momentum row is enough to destroy the flux
gate.** The momentum row in a quasi-hydrostatic layer is the near-cancellation
of two terms that are enormous compared with the wind's mass flux:

| r | `max(\|dF_2\|,\|S_2\|) / F_0`, hot Uranus | WASP-121b |
|---|---|---|
| 1.005 | 4.41e+05 | 7.69e+03 |
| 1.03 | 1.12e+05 | 2.24e+03 |
| 1.10 | 4.46e+03 | 1.08e+02 |
| 1.20 | 1.92e+02 | 1.71e+01 |

so a residual `epsilon` of that row displaces the mass flux by
`epsilon * (|dF_2|/F_0)` of the wind per crossing time. On the hot Uranus at
`r = 1.03`, `epsilon = 2.5e-5` and `|dF_2|/F_0 = 1.1e5`:

```
predicted   |R_2| r^2 dt / F_0 = 4.314e-04 of F_0 per step
measured    slope of F(1.03) over steps 1-20 = 4.714e-04 of F_0 per step
                                                         ratio 1.09
```

**Nine percent agreement, with no free parameter.** The accepted state was flat
to `1.043e-04` over `r >= 1.01`; **one marching step moves the flux at
`r = 1.03` by 4.5 times that whole spread**, and it does so at a constant rate,
because it is simply the time integration of a residual the acceptance norm
divided away.

At the other radii the local residual does not account for the drift (ratios 6
to 500 on the hot Uranus, `5e4` to `2e5` everywhere on WASP-121b), so the rest
of the layer is being moved by its neighbours and by cells 1-2, not by its own
imbalance.

**The departure is linear, not exponential.** Slopes over steps 1-20, in `F_0`
per step:

| | r = 1.005 | 1.03 | 1.10 | 1.20 |
|---|---|---|---|---|
| hot Uranus | 4.85e-05 | 4.71e-04 | 1.86e-06 | -8.02e-07 |
| WASP-121b | 2.80e-04 | 7.58e-05 | 1.81e-05 | 2.29e-05 |

Figure (a) shows the fit: on the hot Uranus the first twenty steps at `r = 1.03`
lie on a straight line to the width of the marker. A linear rise from step 1 is
the signature of a non-zero `du/dt`, not of an instability amplifying round-off.

---

## 5. The variant scan: what damps it (nothing) and what does not

Each variant is 3000 steps from the same accepted state; the number is the
peak-to-peak excursion of the face mass flux at `r = 1.03` in units of `F_0`,
compared at the **same model time** (the end of the shortest run,
`0.109 t_s` / `0.075 t_s`), so a smaller `dt` is not rewarded for covering less
history.

| variant | hot Uranus | vs baseline | WASP-121b | vs baseline |
|---|---|---|---|---|
| **A** baseline, WENO3 throughout | **0.534** | 1.00 | **2.148** | 1.00 |
| **G** the standard `PLM+WENO3` two-stage start (one PLM step) | 0.561 | 1.05 | 2.148 | 1.00 |
| **B** Shapiro filter `0.5`, every 4 steps | 22.55 | **42.3** | 1.976 | 0.92 |
| **F** composition pressure projection removed | 0.600 | 1.12 | 2.209 | 1.03 |
| **C** CFL 0.3 | 0.588 | 1.10 | 2.115 | 0.98 |
| **D** CFL 0.15 | 0.593 | 1.11 | 2.078 | 0.97 |
| **E2** piecewise-constant reconstruction, `j <= 2` | 2.99 | 5.6 | 2.152 | 1.00 |
| **E5** piecewise-constant reconstruction, `j <= 5` | 6.69 | 12.5 | 2.151 | 1.00 |
| **H2** cells `j <= 2` frozen | 0.817 | 1.53 | 2.231 | 1.04 |
| **H5** cells `j <= 5` frozen | **0.043** | **0.081** | 2.197 | 1.02 |
| **H20** cells `j <= 20` frozen | 0.043 | 0.081 | 2.083 | 0.97 |
| **R** restart from the marching's own 3000-step state | **0.013** | **0.025** | **0.675** | **0.31** |

Read across the rows:

* **The time step is irrelevant.** Quartering `dt` changes the excursion by 3
  percent at equal model time. Whatever this is, it is not a stability limit of
  the explicit integrator, and lowering CFL will not buy anything but wall
  clock. (SSP-RK3 is the only integrator in the code; there is no key for the
  order, so only `dt` was varied.)
* **The reconstruction at the base is not carrying it.** Making cells 1-2 or
  1-5 piecewise constant leaves WASP-121b unchanged to three digits and makes
  the hot Uranus 6 to 12 times worse -- first order there is more diffusive but
  also less accurate, and the layer notices the accuracy more than the
  diffusion.
* **The Shapiro filter is not a remedy here.** It removes 8 percent on
  WASP-121b and multiplies the hot Uranus by 42. CETIMB applies it to a base
  that is genuinely ringing; applied to a state that is drifting because its
  residual is not zero, it is a low-pass filter on the wrong signal, and on the
  molecular base it drives its own excursion (`v(1)` reaches `-2.79e-3` against
  the baseline's `-8.19e-4`).
* **The composition pressure projection is not the driver** of the later
  excursion, even though it is the largest term of step 1: removing it moves
  the number by 3 to 12 percent.
* **Freezing the base helps only where the base is the source.** On the hot
  Uranus, freezing `j <= 5` divides the excursion by 12 (and `j <= 2` is not
  enough -- it makes it worse, since a frozen cell 1-2 next to a free cell 3
  is its own discontinuity). On WASP-121b, freezing 20 cells does nothing:
  its layer's departure is not launched at the base.
* **Only a marched start helps on both.** A state the marching itself produced
  is 3 to 41 times quieter under exactly the same operator. That is the
  statement in its cleanest form: the marching has a fixed point, the JFNK has
  a different one, and what is being released is the distance between them.

---

## 6. Verdict

**(a) A discrete fixed-point mismatch. Confirmed, and it is the answer** -- but
the mismatch is mostly in the STATE, not in the OPERATOR.

* The *spatial* operator is identical by construction: `assemble_residual`
  calls the same `Reconstruct` and `RK_rhs` the marching calls (section 2).
* The operator does carry one extra stage, the composition pressure projection,
  and on a restart it is a `dt`-independent jump of up to `6.2e-3` of the cell
  energy (sections 2-3). It is real, it is the largest term of step 1 on
  WASP-121b, and it should be fixed -- but removing it changes the outcome by
  3 to 12 percent (variant F), so it is not what makes the state leave.
* What makes the state leave is that **it is not a steady state**. Its energy
  row is out by `1.2e-3` to `6.0e-3` of that row's own largest term through the
  layer, its momentum row by `2.5e-5` at `r = 1.03`, and its first cell by
  `1.4e-1` and `9.9e-1`. The acceptance measure divides rows 2 and 3 by
  `rho(|v|+c_s)/dr` and volume-weights the cells, and neither survives that:
  cell 1 carries `5.5e-06` (hot Uranus) and `1.9e-04` (WASP-121b) of the
  column's volume, and in the layer `|v| + c_s` is set by `c_s`, not by `v`.
  **P54b removed the Mach factor from the continuity row and left it in the
  other two.**
* The consequence is quantitative and needs no new physics: the flux drift at
  `r = 1.03` on the hot Uranus is `4.71e-4 F_0` per step and the state's own
  momentum residual predicts `4.31e-4` (section 4).

**(b) A genuine discrete instability of the scheme or the base BC. Refuted for
the first twenty steps, and not needed afterwards.** The rise is linear in step
number, not exponential; it is `dt`-independent at equal model time to 3
percent; it survives every reconstruction and filter variant tried; and it is
3 to 41 times smaller from a state the same operator produced. What *is* a
discrete artifact is cell 1's two-cell sawtooth (P44), and it shows up here as
that cell's O(1) residual -- but freezing it does not fix WASP-121b.

**(c) A physical instability. Not supported as the cause.** The later behaviour
-- the excursion of order `F_0` and its slow oscillation (autocorrelation peaks
at 2060 to 2560 steps at `r = 1.03`-`1.20` on WASP-121b, `0.20`-`0.25 t_s`,
consistent with P54's `2100` steps) -- is the layer's damped acoustic *response*
to being handed a finite imbalance, not a mode growing out of nothing. A
physical instability would not care which of two nearby states it started from;
this one is 41 times weaker from the marched state.

---

## 7. What to do, physics first

**(1) Give rows 2 and 3 the measure P54b gave row 1.** Scale each row by the
largest term it differences,

```
   s_1(j) = max( |F_{j-1/2}| r^2 , |F_{j+1/2}| r^2 ) / dV_j      (P54b, done)
   s_2(j) = max( |dF_2(j)| , |S_2(j)| )
   s_3(j) = max( |dF_3(j)| , |S_3(j)| , heat(j) , cool(j) )
```

all of which `RK_rhs` and the ionization sweep already produce, so the cost is
the `max` itself. On the two states here that measure reads `6.0e-3` and
`1.7e-3` where the present one reads `4.0e-8` and `1.1e-4`, i.e. **the gate
would have refused both states**, which is the point. P54b's own experience says
to put it in the ACCEPTANCE test only and leave the Newton system on the
existing state scales (`EXHALE_CONT_SCALE=gate`, not `face`): rescaling the
solver by a quantity that varies by a factor 9 across the first two cells is
what destroyed variant B there, and `s_2`/`s_3` vary by more than that.
Untested here: whether the JFNK can actually reach a state that passes it.

**(2) Stop the restart from injecting the composition kick.** Two ways, both
physics-preserving: re-equilibrate the composition once (`ioniz_eq` at the
loaded `T`) *before* the first hydro step of a `Load IC` run, so the projection
is taken at `t = 0` where it belongs instead of inside step 1; or write the
species file with enough digits that the reload is the state. The first is
better -- it also makes the reloaded residual reproduce the accepted one, which
P54 records as an open item.

**(3) Make the composition energy step conservative, or put it in the
residual.** `comp_T_from_p` -> sweep -> `comp_p_from_T` holds `T` and lets the
particle count move, so the internal energy changes by an amount no heating or
cooling term accounts for. That is a defensible operator split only if the
reaction energetics in `heat`/`cool` are the *complement* of it; whether they
are was not checked here. Either way the steady residual should contain the
same stage, otherwise the two operators cannot have the same fixed point except
where the sweep is exactly idempotent.

**(4) A characteristic (non-reflecting) base boundary is the next experiment,
not the first.** The base is where the hot Uranus's problem lives (freezing
`j <= 5` divides its excursion by 12) and it is where every row of both states
is worst. But it is *not* where WASP-121b's lives (freezing 20 cells does
nothing), so a base BC alone cannot be the whole answer, and it should be tried
after (1), on states the new measure accepts.

**(5) Filters last, and not this one.** The Shapiro filter as configured is
measured to be 42 times worse than doing nothing on the molecular hot Uranus.
If a filter is wanted later it should be justified on a state whose residual is
actually zero, where what is left to damp is a wave and not a drift.

---

## 8. What was not established

* **Only two states, one configuration of the base.** Both use
  `rho = rho_bc` + the hard one-way valve + the isothermal ghost. The
  `Base velocity: massflux` exponential moving average and the
  continuous-temperature ghost (the two places where the marching carries state
  the residual does not have) were therefore **never exercised**, and the
  section-2 finding that the ghost fill is not a difference is a statement about
  these runs only.
* **Cross-run comparison of `dU/dt` against `R`.** The residual is measured in a
  separate run whose ionization sweep has already corrected the reloaded
  composition, while the marching's step 1 uses the file's. The caloric EOS is
  composition-dependent, so the two runs' `p` differ slightly and the mass and
  momentum rows disagree by factors up to about 2 at cells where both are
  `~1e-7`. That is why section 3 uses the `dt` scan (a single run, internally
  consistent) for the operator question and section 4 uses the predicted vs
  measured *slope* rather than a cell-by-cell subtraction.
* **No characteristic base boundary was tried**, and no grid refinement. The
  variants tested are reconstruction order, filter, time step, the chem
  projection, and freezing cells.
* **Whether the marching's attractor is the physical steady state** is still
  open, exactly as P54b section 10.5 left it. Variant R shows it is much closer
  to being a fixed point of the marching; it says nothing about which of the two
  is right.
* **Whether the JFNK can converge on the row-2/row-3 measure of section 7** was
  not tested. P54b's variant B is a warning that a new scale can make the solve
  unreachable.
* The autocorrelation periods quoted are weak (`r = 0.14` to `0.22` on
  WASP-121b, `0.01` to `0.09` on the hot Uranus, whose baseline barely
  oscillates), so "the period is consistent with P54's 2100 steps" is one
  measured number against another, not an identification of a mode.

---

## 9. Files

```
docs/p55_base_mode.md         this note
docs/p55_base_mode.png        the departure, the series, the variant scan
docs/p55_operator_terms.png   the accepted state by its own terms; the restart kick
docs/p55_probes.diff          the six env gates, against the P54b working copy
docs/p55_reload.png           the reload figure (section 10)
src/, Makefile, EXHALE.x      the probe build (copy of .../p54flux + the diff)
states/hotU_G, states/WASP_G  read-only copies of .../p54flux/n1G and /w1G
   (the probe build and the run directories stayed in the scratch tree; only
    this note, its two figures and the reload figure are carried here)
runs/hotU_res, hotU_resD, wasp_res   EXHALE_RESIDUAL (+P54_FLUX, +P55_MOM) on each state
runs/<planet>_du40            EXHALE_P55_DU=41, 40 steps
runs/<planet>_cfl<c>          the one-step dt scan
runs/<planet>_c40_<c>         the 39-step dt scan
runs/<planet>_{A,G,B,F,C,D,E2,E5,H2,H5,H20,R}   the variant scan, 3000 steps each
runs/nullchk_p54, runs/nullchk_p55   the default-off byte-identity check
an_du.py an_du2.py an_cmp.py an_cfl.py an_ts.py an_ts2.py fig1.py fig2.py
mkrun.sh mkvar.sh runvar.sh patch_p55.py
```

---

## 10. Addendum -- the restart kick removed, and what it turned out to be

Written after sections 1-9, with the fix of `p55_reload.patch` in the build.
**It amends the WASP-121b half of the verdict.** Section 5 says no knob damps
the departure; one does, and section 5 did not try it.

### 10.1 The fix

`equilibrate_loaded_composition` puts a state read from a file on the fixed
point of its own ionization sweep before the first hydro step, holding `rho`,
`v` and `p` and iterating

```
   f_sp  <-  ioniz_eq( rho, T ),      T = p / ( n_tot(f_sp) + n_e(f_sp) )
```

with the same sweep the loop calls. On by default for a state read from a file;
`EXHALE_RELOAD_EQ=0` restores the previous behaviour bit for bit. Full account
in `docs/Update_EXHALE.md` section 144; figure `p55_reload.png`.

**Verification (i) -- the projection becomes a rate.** `dE_chem/E` in the first
step, at four time steps:

| | CFL 0.6 | 0.06 | 0.006 | 0.0006 |
|---|---|---|---|---|
| WASP-121b `r=1.10`, before | -6.2295e-03 | -6.2297e-03 | -6.2297e-03 | -6.2297e-03 |
| WASP-121b `r=1.10`, after | **-3.2041e-06** | -3.2067e-07 | -3.2070e-08 | -3.2072e-09 |
| hot Uranus `r=1.03`, before | 5.4549e-07 | 5.5078e-07 | 5.5131e-07 | 5.5136e-07 |
| hot Uranus `r=1.03`, after | **-5.8808e-09** | -5.8712e-10 | -5.8546e-11 | -5.6966e-12 |

Before the fix the number does not move with `dt` at all -- four decades, four
identical values, which is what makes it a projection. After it, each decade of
`dt` divides it by ten to four digits, which is what makes it a rate. The
first-step magnitude falls by **1944** (WASP-121b at `r=1.10`) and by **93**
(hot Uranus at `r=1.03`).

**Verification (ii) -- it does not change how the state leaves.** The 20-step
slopes and the 3000-step excursion reproduce the "chem projection removed"
variant (F) of section 5 to two digits: hot Uranus p2p at `r=1.03` is 0.530
against the baseline's 0.534, WASP-121b 2.219 against 2.148. **The kick was
real, it is removed, and it was not what drove the layer** -- which is what
section 6 concluded from variant F and is unchanged.

### 10.2 What the departure on WASP-121b actually was

Measuring the fix exposed it. WASP-121b's input carries no
`Secondary_ionization:` key, so the SvS85 coupling is STAGED. The run that
produced the state armed it at step 2242 and the JFNK finished with it on. The
run that restarts the state re-stages it, and the flip never fires again inside
3000 steps -- so **every marching experiment of sections 5 and 6 on WASP-121b
was run at a photoionization rate the state was not converged under.**

The restart file now says what coupling it was written under -- a `# coupling:`
line written by `write_coupling_state_header` and read by `load_IC` -- and the
restart starts there. Peak-to-peak `\Delta(Fr^2)/F_0` over the full 3000 steps:

| | `r=1.005` | `1.03` | `1.10` | `1.20` | `v(1)` end | `du` end |
|---|---|---|---|---|---|---|
| baseline (section 5) | 7.9689 | 4.2961 | 1.1976 | 0.6349 | 5.766e-03 | 1.178e-02 |
| composition equilibrated | 8.0041 | 4.3197 | 1.1970 | 0.6341 | 5.880e-03 | 1.130e-02 |
| **+ the coupling header** | **0.0382** | **0.0174** | **0.0051** | **0.0040** | **3.150e-03** | **5.187e-05** |

A factor 250 at `r = 1.03`, and `du` ends at `5.19e-05` against the `5.03e-05`
the accepted state carries. **On WASP-121b the JFNK root IS a fixed point of the
marching, once the marching runs the physics the root was found under.** The
hot Uranus is unchanged to five digits, as it must be: its input sets
`Secondary_ionization: Immediate`, so it never had the discontinuity.

### 10.3 The verdict, amended

* **Hot Uranus: sections 4 and 6 stand unchanged.** No configuration
  discontinuity exists there, the departure is linear from step 1, and the
  state's own momentum-row residual at `r = 1.03` predicts the measured flux
  drift to nine percent. The accepted state is not a steady state, and the
  acceptance measure cannot see it. Prescription (1) of section 7 -- give rows 2
  and 3 the measure P54b gave row 1 -- is unaffected.
* **WASP-121b: the departure was overwhelmingly a physics-configuration
  discontinuity at restart, not a discrete mode.** Sections 5 and 6 attributed
  it to the operator; it belonged to the run's own staging. What remains after
  the arming (p2p 0.017 at `r = 1.03`) is 130 times smaller than what section 5
  measured and is consistent with the state being a fixed point.
* The section 6 sentence "the marching has a fixed point, the JFNK has a
  different one" is therefore **true on the hot Uranus and not established on
  WASP-121b**. So is the reading of P54b section 10.5, whose WASP-121b marching
  experiments were re-staged in exactly the same way.
* Section 5's "no knob damps it" should be read as "no knob among those tried",
  and the knob that was missing is not a numerical one.

### 10.4 Regression scope

All ten cases of the default matrix are cold starts (`Load IC? False`) and none
of them reaches the JFNK hand-off, so neither half of the fix can touch them --
checked by running all ten for 600 steps against the previous binary and
comparing the way the harness does (`grep -v '^ *#'`): all PASS. The one added
`#` line is the only change to those files' bytes. Every named Load-IC case
moves; the census, and the one bounded Load-IC Newton case run to its cap, are
in `docs/Update_EXHALE.md` section 144.

The addendum was measured on a build WITHOUT the G23 row-scale change, which was
not published when it was made; only the `||R||` comparison reads a row scale and
would move under it. The patch touches `EXHALE_main.f90`, `utilities.f90`,
`parameters.f90`, `write_output.f90`, `load_IC.f90` and `write_setup_report.f90`
-- none of G23's files.

### 10.5 What the addendum did not establish

* The coupling restore is **driven by the file and by nothing else**. An IC that
  carries no `# coupling:` line -- `jfnk_hd189`, written 2026-08-09 -- starts
  staged exactly as before, bit for bit, which is what it needs: arming that one
  blindly makes its composition mismatch worse (4.12e-3 staged -> 6.26e-3 armed)
  and produces a `du = 1.0e2` transient by step 5500. There is no environment
  variable for it any more.
* Whether the hot Uranus has its own configuration discontinuity of another kind
  was not looked for; its restart mismatch is `1.8e-05`, which leaves little
  room for one.
* The WASP-121b run with the arming was not carried past 3000 steps, so
  "fixed point" here means "did not leave in 3000 steps at `du = 5.2e-05`",
  not a proof.
* The residual of the reloaded WASP-121b state improves by 14 to 17 in the hydro
  rows but its energy row still reads `1.0e-05` against the `3.1e-06` the solve
  accepted; that remainder was not chased.

---

## 11. The norm the gate reads, and where the layer's 3.4e-3 comes from

Written after section 10, on the tree as it stands (section 143's row scales
plus section 144's restart handling). It resolves an apparent contradiction:
section 144.5 reports a root with `||R|| = 9.958e-06` whose mass flux is
`3.433e-03` out over `r >= 1.03`. If every row is scaled by its own largest
term, how can a converged root have a non-flat layer?

**Answer: it does not. The conserved mass flux of that root is exactly
constant. The `3.4e-03` is the difference between the cell-centred product the
flux gate measures and the Riemann face flux the scheme conserves -- two
different functionals, differing here by five orders.** The norm does hide
something, but not the layer: it hides the base cells on the hot Uranus and
fifty WIND cells on WASP-121b.

Variant and probes: `p55_cellmax.diff` (three files, all env-gated, all off by
default) as measured. **Section 145 has since adopted the cellwise-max norm and
the face-flux gate as the tree defaults, and the section 133 denominator naming
is corrected there; read this section as the measurement that led to them, not
as a description of the current code.**

### 11.1 What the gate's `||R||` is, in the code

`steady_gates_met(rnorm, ...)` tests `rnorm < resid_tol`. `rnorm` is built by
`resid_relnorm` (`steady_newton.f90`) as

```
   rc_wind  = relnorm_over_cells(R, u, j_min, N,     resid_vol)
   rc_layer = relnorm_over_cells(R, u, 1,   j_min-1, resid_vol)
   rnorm    = maxval( max(rc_wind, rc_layer) )
```

and `relnorm_over_cells` returns, per row `k`,

| `resid_vol` | what it returns |
|---|---|
| `.true.` (**the default**) | `sum_j |R_kj| V_j / sum_j s_kj V_j`, `V_j = r_j^2 dr_j` |
| `.false.` | `max_j |R_kj| / max_j s_kj` -- a **ratio of maxima** |

**Neither is a maximum over cells of a cell's own scaled residual.** The first
is a ratio of volume-weighted sums; the second takes the two maxima
independently, and they are not in the same cell -- `s ~ 1/dr` is largest in the
smallest cells at the base, while `|R|` need not be. So a cell whose own
`|R|/s` is large can be paid for by other cells under either.

`residual_norms` / `residual_norms_vol` (the `EXHALE_RESIDUAL` diagnostic) are
the same two forms.

### 11.2 The residual the gate actually read, cell by cell

Measured by dumping the solver's own `F` at the acceptance point
(`EXHALE_P55_CELL=1`), divided by each cell's own `residual_row_scale`. This is
the quantity `rnorm` is a norm of -- not a reload, and not a re-assembly (see
11.5 for why that distinction is not pedantic).

Hot Uranus, the section 144.5 root:

| r | `R1/s1` | `R2/s2` | `R3/s3` |
|---|---|---|---|
| 1.00019 (cell 1) | 9.07e-06 | 9.26e-07 | 3.16e-05 |
| 1.00039 (cell 2) | 8.83e-06 | 4.50e-08 | 9.88e-06 |
| 1.00058 (cell 3) | **3.10e-05** | 2.91e-07 | **2.17e-03** |
| 1.00097 | 7.69e-10 | 1.21e-08 | 1.57e-06 |
| 1.00502 | 3.00e-11 | 1.01e-08 | 1.72e-06 |
| 1.03016 | 3.73e-09 | 2.43e-08 | 1.77e-06 |
| 1.09969 | 3.37e-11 | 2.06e-09 | 5.13e-06 |
| 1.20011 | 1.04e-11 | 4.35e-11 | 1.21e-06 |
| 2.00271 | 1.66e-12 | 1.27e-11 | 3.60e-07 |

**Cells above the `1e-5` tolerance, by region:**

| state | 1.00-1.30 | of which 1.03-1.30 | whole column |
|---|---|---|---|
| hot Uranus | 1 mass, 0 mom, 6 energy | **0, 0, 0** | 1, 0, 6 (all in cells 1-3) |
| WASP-121b | 0, 0, 20 | 0, 0, 18 | 0, 0, **70** |

**The layer is not hidden by the norm.** Between `r = 1.03` and `1.30` not one
cell of either row exceeds `1e-5` on the hot Uranus; the layer's own scaled
residual there is `1e-11` to `5e-6`. The dilution is real but it is elsewhere:

* on the **hot Uranus** the seven offending cells are **cells 1-3**, the base
  (item (V)), where the mass row is `3.10e-05` and the energy row `2.17e-03`
  against a gate of `1e-5` -- 218 times the accepted `||R||` of `9.958e-06`.
  The run's own log even names the cell (`max cell j=3 r=1.0006 k=3`); it prints
  the diluted number beside it.
* on **WASP-121b** the seventy offending cells are mostly in the **WIND**: 50 in
  `r = 1.30-2.00`, 18 in `1.10-1.30`, 2 at the base, at `1.0e-05` to
  `6.4e-05`. There the volume weighting is diluting a wind, not a layer.

### 11.3 Where the 3.4e-03 comes from: the two flux functionals

The mass row's scale is `max(|F_in| r^2, |F_out| r^2)/dV`, so `R1/s1` is
**exactly** the fractional change of the conserved face flux across the cell.
Summing it over a window bounds the total relative variation of that flux. Both
bound and direct measurement, on the same accepted states:

| case | window | face-flux spread | `sum R1/s1` (the bound) | cell-centred spread (**the gate**) |
|---|---|---|---|---|
| hot Uranus | `r >= 1.03` | 1.21e-08 | 9.33e-08 | **3.433e-03** |
| | `r >= 1.10` | 0 (one unique value) | 4.38e-09 | 2.375e-03 |
| | `r >= 1.20` | 0 | 3.17e-10 | 2.501e-04 |
| WASP-121b | `r >= 1.03` | 0 | 1.29e-09 | 6.751e-05 |
| | `r >= 1.20` | 0 | 3.16e-10 | 5.242e-05 |

Above `r = 1.10` the face mass flux takes **one unique double-precision value
over three hundred cells** (`4.94565219e-05` on the hot Uranus). It varies only
in cells 1-3, where it falls `4.80e-04 -> 3.49e-04 -> 4.95e-05` -- the base
sawtooth of `docs/p44_base_sawtooth.md`.

**So the wind carries one mass flux to the last bit, and the flux gate reports
`2.5e-04` at its own window.** What the gate measures is the cell-centred
`rho v r^2`, which is not the conserved quantity; the number it returns on this
state is the reconstruction difference between the two, not a failure of
conservation. P54's recommendation 7(a) -- measure the flux at the faces -- was
never carried in, and this is a far sharper argument for it than the one that
was made there.

### 11.4 The `1e-4` of P54 against the `3.4e-3` of section 144.5: neither norm nor physics

Both numbers are right and they are **different functionals**. Measured with one
script on all three states:

| state | cell-centred `r>=1.03` | `r>=1.10` | `r>=1.20` |
|---|---|---|---|
| hot Uranus, pre-P53 G root (`p54flux`) | 3.523e-03 | 2.365e-03 | 2.489e-04 |
| hot Uranus, section 143 root, default norm | 3.433e-03 | 2.375e-03 | 2.501e-04 |
| hot Uranus, section 143 root, **cellmax** norm | 3.433e-03 | 2.375e-03 | 2.501e-04 |
| WASP-121b, pre-P53 G root | 7.747e-05 | 5.627e-05 | 5.627e-05 |
| WASP-121b, section 143 root, default norm | 6.751e-05 | 5.242e-05 | 5.242e-05 |
| WASP-121b, section 143 root, **cellmax** norm | 6.751e-05 | 5.243e-05 | 5.243e-05 |

The pre-P53 root's **cell-centred** spread at `r >= 1.03` is `3.523e-03`, within
3 percent of the section 143 root's `3.433e-03`. P54 section 10.4's `1.043e-04`
is that same state's **face-flux** spread. **Nothing changed between P53 and now,
and nothing changes with the norm**: the difference is which of the two flux
measures is quoted. Section 144.5's sentence "the new root's own flux is still
3.433e-03 out ... fourteen times the gate's number" compares a cell-centred
spread with a face spread and should be read that way.

On the functional that IS the conserved one, the two roots are not equal, and the
newer one is better by four orders:

| face-flux spread | `r>=1.03` | `r>=1.10` | `r>=1.20` |
|---|---|---|---|
| hot Uranus, pre-P53 G root | 1.043e-04 | 6.077e-06 | 3.778e-06 |
| hot Uranus, section 143 root | **1.213e-08** | **0** | **0** |

That is a real improvement and it is invisible to the gate, which reads the
other functional and reports `2.501e-04` on the state whose conserved flux is a
single double-precision number.

### 11.5 An acceptance norm that is the maximum over cells: it converges, it is cheap, and it changes nothing else

`EXHALE_RESID_NORM=cellmax` (default off) makes `relnorm_over_cells` return
`max_j |R_kj|/s_kj` -- "every cell is steady on the scale its own physics sets",
which is what `residual_norms`' own comment says the norm is for. It changes the
ACCEPTANCE measure only: `cell_state_scales` and the line-search merit are built
from `state_scales_of_cell` and are untouched.

| | hot Uranus, default | hot Uranus, **cellmax** | WASP-121b, default | WASP-121b, **cellmax** |
|---|---|---|---|---|
| JFNK `info` | 0 | **0** | 0 | **0** |
| JFNK iterations | 7 | 9 | 8 | 12 |
| `||R||` reported (its own measure) | 9.958e-06 | 2.002e-06 | 7.991e-06 | 9.764e-06 |
| cell-max, mass row | 3.100e-05 | **2.021e-07** | 3.232e-07 | **4.374e-08** |
| cell-max, momentum row | 9.255e-07 | 2.077e-07 | 2.540e-06 | 1.987e-07 |
| cell-max, energy row | 2.174e-03 | **2.002e-06** | 6.397e-05 | **9.764e-06** |
| cells above `1e-5` | **7** | **0** | **70** | **0** |
| flux spread `r>=1.03` / `r>=1.20` | 3.4327e-03 / 2.501e-04 | 3.4327e-03 / 2.501e-04 | 6.7506e-05 / 5.242e-05 | 6.7507e-05 / 5.243e-05 |
| `log10 Mdot` | 10.34 | 10.34 | 13.21 | 13.21 |
| 3000-step departure, p2p at `r=1.03` | 0.5333 | 0.5334 | 0.0103 | 0.0104 |
| `du` at step 3000 | 3.363e-04 | 3.377e-04 | 5.073e-05 | 5.087e-05 |

**It converges on both, in 2 to 4 extra JFNK iterations, and it does what it
says**: no cell of any row is left above the tolerance, including the base cells
the volume weighting was hiding, and including WASP-121b's fifty wind cells.

**And it changes nothing else.** The mass-loss rate is identical to the printed
digit, the flux spreads agree to five digits, and the 3000-step departure from
the new root is the same to 0.2 percent. So the norm is not what makes a root
leave the marching loop either -- the third candidate for item (AD), after the
restart (section 144.2, 144.3) and the row scale (section 144.5), is excluded.

### 11.6 A separate finding: the accepted state's residual is not reproducible

Re-assembling the residual on the SAME accepted `u` in memory, after the
post-solve refresh has run one more ionization sweep, does not give the residual
the gate read:

| r | mass, solver `F` | mass, re-assembled | energy, solver `F` | energy, re-assembled |
|---|---|---|---|---|
| 1.00019 | 9.07e-06 | **2.72e-01** | 3.16e-05 | **9.95e-01** |
| 1.00039 | 8.83e-06 | **8.58e-01** | 9.88e-06 | **9.98e-01** |
| 1.00058 | 3.10e-05 | 4.50e-05 | 2.17e-03 | 3.16e-03 |
| 1.03016 | 3.73e-09 | 1.02e-08 | 1.77e-06 | 8.14e-06 |
| 2.00271 | 1.66e-12 | 1.66e-12 | 3.60e-07 | 5.73e-06 |

Volume-weighted, row 3 goes from `~1e-05` to `4.554e-01`. Cells 3 and outward
move by factors 2 to 7 -- the ordinary non-idempotency of a sweep at
`xtol = sqrt(eps)`. **Cells 1 and 2 move by 30000**, and in the MASS row, which
does not read `heat` or `cool` at all: what moved is `rho`, which `ioniz_eq`
rewrites through `calc_rho`. The base cells are so close to the reconstruction's
switching surfaces that re-solving the chemistry on the same conserved state
moves their residual from `1e-5` to order unity.

This is item (V) in a new form, and it is why every number in 11.2 and 11.3 is
taken from the solver's own `F` and not from a re-assembly.

### 11.7 What section 11 did not establish

* **Only the two states of section 144.5**, both under the tree's current row
  scales. The cellmax variant was not run on the regression matrix.
* **Whether the flux gate should move to the faces.** The measurement says the
  functional it uses is not the conserved one; it does not say what threshold a
  face-flux gate should carry, and every threshold in the code is calibrated
  against the present one.
* **Why cells 1-3 sit where they do.** Item (V). The cellmax norm forces the
  solve to bring them below `1e-5` and it succeeds on both cases, which is new,
  but nothing here says whether the state it reaches there is physical.
* The `EXHALE_RESIDUAL` diagnostic's printed header read `max|R|/max|u|` and
  `sum|R|V/sum|u|V`, whereas since section 133 the denominator has been
  `residual_row_scale`, not `|u|`. Corrected in section 145's rewrite of the
  breakdown report.
