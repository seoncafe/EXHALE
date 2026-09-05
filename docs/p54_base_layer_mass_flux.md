# P54 -- Is the base layer's mass flux flat on a state both steady gates accept?

**Short answer: no, and it is not a diagnostic artifact.** The quantity the
finite-volume scheme actually conserves -- the Riemann face flux
`F_{j+1/2} r_{j+1/2}^2` -- is as non-flat as the cell-centred `rho v r^2` the
open item (AA) quotes, to within one percent of the spread at every window
tested above `r = 1.002 R_p`. The accepted state carries a smooth 30 percent
mass-flux deficit centred on `r = 1.027` that recovers only by `r = 1.15`, the
continuity row has no source that could produce it, and marching the accepted
state forward moves it immediately: a coherent damped oscillation of period
about 2100 steps (`3.5e3` s) and an amplitude of order the wind flux itself,
which after 20000 steps has left the layer on a *different* non-flat profile
that both gates would now reject. The residual gate does not see any of this
because the continuity row's scale exceeds the flux error by exactly the
reciprocal of the local Mach number, which in the layer is `4e-5` to `5e-4`.

Figure: `p54_base_layer_mass_flux.png` -- (a) the two flux measures on the
accepted state and the profile 20000 steps later, (b) the Mach-number identity
that makes the residual gate blind, (c) the time series.

Reading (i) of item (AA) -- "the finite-volume state is steady in its own
discrete sense while the cell-centre products are not" -- is therefore **not
supported**. Reading (ii) -- "the gate's window genuinely does not cover a
region in which the solution is not steady" -- is what the measurement shows.

---

## 1. What was measured, on what, with what

**State.** The 1 microbar hot-Uranus `He/H = 0.0793` rung of the P23 ladder,
the run that meets both steady gates:

```
.../scratchpad/ladder2/runs/heh0p0793cc_r1/
  JFNK done info=0  ||R|| = 2.944e-06  flux spread = 2.460e-03
  log10 Mdot = 10.30 g/s
  molecular chemistry on, He 2^3S on, Secondary_ionization: Immediate,
  base.inp p_base = 1e-6 bar, q_H2_base = 0.84, 500 cells (base 2e-4 x 50)
```

Read only. Every measurement below is made by re-loading that run's
`Hydro_ioniz.txt` / `Ion_species.txt` as an initial condition.

**The reload is faithful for the hydro state.** The re-measured mass row of the
volume-weighted residual over the wind is `2.9440e-06` against the solver's
accepted `2.944e-06`, and the re-measured cell-centred flux spread over
`r >= 1.2` is `2.4600e-03` against the accepted `2.4600e-03`. (The energy row
does not reproduce as tightly: the solver accepted `2.944e-06` as its maximum
over rows, while on the reload the energy row reads `2.68e-05`. The most likely
cause is the finite precision of the species file the ionization solve is
restarted from; it was not chased, and it does not enter any number below,
because the mass row is what P54 is about.)

**Binary.** `src/` + `Makefile` copied from the repository working tree on
**2026-09-03 13:25:27 KST**, when the copy was byte-identical to the tree
(`md5` over all `src/*.f90|*.py|*.txt`: `a6b3b040a6f2111345f0b6037e739adc` on
both). It compiled unmodified (`PATH=/usr/bin:$PATH make`, gfortran -O3), so
the P53 completion tree was not needed. The two probe files were later
re-derived from the tree originals and the patch reproduces the working copy
byte for byte, so the diff in `p54_probes.diff` is exact against the tree as it
stood.

**Probes** (`p54_probes.diff`, both env-gated, nothing else touched):

| gate | what it writes |
|---|---|
| `EXHALE_P54_FLUX=1` (inside the existing `EXHALE_RESIDUAL=1` branch) | `output/p54_mass_flux.txt`: per cell, the cell-centred `rho v r^2`, the two Riemann face fluxes `F r^2` bounding it, the continuity residual `R_1` unscaled, its scale `s_1 = residual_row_scale(1,j,u)`, the ratio, `dr` and `c_s` |
| `EXHALE_P54_TS=<interval>` (marching loop) | `output/p54_flux_ts.txt`: one line every `<interval>` steps with the face flux at `r = 1.005, 1.03, 1.10, 1.20` and the two cell-centred window spreads |

`face_flux` is the array `RK_rhs` already fills; the probe only reads it, after
`assemble_residual` has run on the state, so the flux written is the flux the
update used.

Units: the code's own. `F0` below is the mean face flux over `r >= 1.2`,
`4.4642e-05` in code units. Physical scales for the same run:
`R0 = 3.5031e9 cm`, `v0 = 3.0667e5 cm/s`, `t_s = R0/v0 = 1.1423e4 s`,
`dt = CFL min_j dr_j/(|v_j|+c_j) = 1.68 s` on this state.

---

## 2. H1 -- measurement artifact? Refuted above `r = 1.002`, confirmed only in cell 1

**Windowed spreads, `(max - min)/|mean|`, on the accepted state.**

| `r >=` | cells | spread, cell centre | spread, Riemann face | min `F/F0` | max `F/F0` |
|---|---|---|---|---|---|
| 1.000 | 500 | 2.722e+01 | 1.316e+00 | 0.6963 | 1.9187 |
| 1.002 | 490 | 3.247e-01 | 3.282e-01 | 0.6963 | 1.0004 |
| 1.005 | 475 | 3.248e-01 | 3.284e-01 | 0.6963 | 1.0004 |
| 1.010 | 449 | 3.244e-01 | 3.280e-01 | 0.6963 | 1.0004 |
| 1.020 | 409 | 3.193e-01 | 3.229e-01 | 0.6963 | 1.0004 |
| 1.030 | 383 | **3.105e-01** | **3.137e-01** | 0.6997 | 1.0004 |
| 1.050 | 347 | 2.152e-01 | 2.170e-01 | 0.7874 | 1.0004 |
| 1.100 | 294 | **2.920e-02** | **3.041e-02** | 0.9701 | 1.0004 |
| 1.150 | 262 | 6.554e-03 | 6.710e-03 | 0.9937 | 1.0004 |
| 1.200 | 240 | **2.460e-03** | **2.513e-03** | 0.9979 | 1.0004 |
| 1.500 | 165 | 7.812e-04 | 7.867e-04 | 0.9996 | 1.0004 |

The three bold cell-centre numbers are item (AA)'s `3.10e-1 / 2.91e-2` and the
gate's `2.46e-3`, reproduced. **The face column is the same number to within
one to four percent of itself at every window from 1.002 outward.** Whatever
makes the flux non-flat is in the conserved flux, not in the reconstruction of
it at cell centres.

**Where the two do differ is cell 1 and only cell 1** -- the odd-even mode of
`docs/p44_base_sawtooth.md`:

| j | r | v (code) | `F_cell/F0` | `F_inner face/F0` | `F_outer face/F0` |
|---|---|---|---|---|---|
| 1 | 1.000193 | -7.641e-04 | **-21.1411** | 1.7476 | 1.9187 |
| 2 | 1.000386 | +1.065e-04 | 2.9914 | 1.9187 | 0.8331 |
| 3 | 1.000579 | +1.595e-05 | 0.4549 | 0.8331 | 0.9693 |
| 4 | 1.000772 | +4.026e-05 | 1.1650 | 0.9693 | 0.9341 |
| 5 | 1.000965 | +2.919e-05 | 0.8560 | 0.9341 | 0.9493 |
| 6 | 1.001158 | +3.252e-05 | 0.9653 | 0.9493 | 0.9440 |
| 10 | 1.001930 | +3.038e-05 | 0.9367 | 0.9461 | 0.9462 |

Cell 1's centred product is `-21 F0` -- negative, and 21 times the wind flux in
magnitude -- while both of its faces carry `+1.75` and `+1.92 F0`. This is the
P44 finding in its cleanest form: the cell-1 velocity is a two-cell artifact and
the face flux there is inflow of order the wind flux, not `-21` times it. It is
also the whole reason the cell-centred spread over all cells (27.2, with a sign
change, so `max/min` is not even defined) is meaningless while the face spread
over all cells (1.316, single-signed) is at least a number.

**Verdict H1: false above `r = 1.002 R_p`, true only for cell 1.** The face
flux removes the sawtooth and nothing else. It does not rescue the gate.

---

## 3. H2 -- does the residual scale hide it? Yes, by exactly the Mach number

The mass row has no source (`Source.f90`: `S(1) = 0.0`), so the unscaled
continuity residual is the flux divergence itself,

```
R_1(j) = [ F_{j+1/2} r_{j+1/2}^2 - F_{j-1/2} r_{j-1/2}^2 ] / dV_j ,
```

and the gate reads `R_1/s_1` with `s_1 = rho_j (|v_j| + c_j)/dr_j`
(`residual_row_scale`). Writing `dV_j ~ r_j^2 dr_j` and `F ~ rho v r^2`,

```
   R_1/s_1  =  [ D(F r^2)/F0 ] * ( F0 / (rho r^2 (|v| + c_s)) )
            ~  [ D(F r^2)/F0 ] * |v|/(|v| + c_s)   ~   (dF/F) * Mach .
```

**Measured, on the accepted state**, the ratio of `|R_1|/s_1` to
`(dF/F0) * Mach` is `1.109, 1.004, 1.015, 1.007, 1.008, 1.010, 1.006, 0.979,
0.908, 0.656` at `r = 1.001, 1.004, 1.010, 1.023, 1.032, 1.048, 1.093, 1.177,
1.329, 2.119`. The identity holds to one percent everywhere the flow is
subsonic; it degrades only where `|v|` approaches `c_s` and the approximation
`|v|/(|v|+c_s) ~ Mach` stops being one.

| r | Mach | `s_1` | `\|R_1\|` | `\|R_1\|/s_1` | `dF/F0` across the cell | `dF/F0` the gate would still pass at `10^{-5}` |
|---|---|---|---|---|---|---|
| 1.00019 | 9.72e-04 | 5.03e+03 | 3.96e-02 | 7.86e-06 | 1.71e-01 | 1.03e-02 |
| 1.00058 | 2.07e-05 | 5.07e+03 | 3.15e-02 | 6.21e-06 | 1.36e-01 | 4.82e-01 |
| 1.00097 | 3.88e-05 | 5.10e+03 | 3.50e-03 | 6.87e-07 | 1.52e-02 | 2.58e-01 |
| 1.00405 | 4.38e-05 | 4.94e+03 | 1.33e-04 | 2.69e-08 | 5.79e-04 | 2.28e-01 |
| 1.00985 | 5.25e-05 | 3.59e+03 | 8.46e-04 | 2.35e-07 | 3.80e-03 | 1.91e-01 |
| 1.02336 | 9.72e-05 | 8.46e+02 | 1.05e-04 | 1.24e-07 | 8.89e-04 | 1.03e-01 |
| 1.03151 | 1.49e-04 | 4.27e+02 | 9.01e-05 | 2.11e-07 | 9.85e-04 | 6.70e-02 |
| 1.04811 | 3.70e-04 | 1.28e+02 | 2.32e-04 | 1.82e-06 | 3.77e-03 | 2.70e-02 |
| 1.09346 | 4.51e-03 | 6.54e+00 | 5.86e-05 | 8.95e-06 | 1.90e-03 | 2.22e-03 |
| 1.17655 | 3.42e-02 | 4.34e-01 | 2.59e-06 | 5.98e-06 | 1.78e-04 | 2.93e-04 |
| 1.32880 | 1.15e-01 | 5.97e-02 | 1.30e-07 | 2.19e-06 | 2.09e-05 | 8.70e-05 |
| 2.11892 | 5.42e-01 | 2.05e-03 | 3.29e-09 | 1.60e-06 | 4.50e-06 | 1.85e-05 |
| 4.72534 | 1.57e+00 | 7.37e-05 | 3.76e-10 | 5.10e-06 | 8.34e-06 | 6.36e-06 |

The last column is the statement in the form that matters: **at `r = 1.01` the
residual gate at `10^{-5}` still passes a mass flux that changes by 19 percent
from one cell to the next**, and at `r = 1.03` by 6.7 percent. At `r = 1.3` and
above the same tolerance constrains the flux to `1e-4` per cell. The gate is not
weak in the wind; it is weak in the layer, in proportion to `1/Mach`.

Two further points about how the norm is formed, both checked and neither the
dominant effect here:

* `relnorm_over_cells` with `volume_weighted = .false.` forms
  `max_j |R| / max_j s`, a ratio of maxima rather than a maximum of ratios. On
  this state the two agree closely for the mass row in the layer
  (`4.90e-05` vs `4.97e-05`), so the ratio-of-maxima construction is not what
  hides the layer here; the Mach factor is.
* The below-escape "layer" region is `[1 : j_min-1]` with `j_min = 392`, i.e.
  everything below `r_esc = 2.0 R_p`. The whole non-flat region is inside it,
  and the volume-weighted mass row over that region reads `1.08e-06`. It is not
  that the layer is averaged in with the wind; the layer's own number is small.

**Verdict H2: confirmed, and quantified.** The scaled continuity residual is the
fractional flux error per cell times the local Mach number. A layer at
`Mach ~ 5e-5` is invisible to it by a factor `2e4`.

---

## 4. H4 -- is anything adding mass? No

Enumerated in the source, for the configuration this state runs:

| candidate | status on this state | writes to `u(1,:)`? |
|---|---|---|
| hydro source term | `Source.f90`, `S(1) = 0.0` unconditionally | no |
| base ghost / valve | `BC_component_constrho` pins `rho = rho_bc` and sets `v_ghost` by the smooth valve; a boundary condition, not a source | ghosts only |
| Shapiro filter | `shapiro_eps <= 0`, reported "Shapiro filter: off" | no-op |
| low-Mach contact damping | `low_mach_damping_active()` false, reported off | no-op |
| positivity flux limiter | `positivity_limited_fluxes`; conservative substitution of a first-order flux, and the run reports no limited faces | conservative |
| ghost positivity fallback | zero-gradient fallback of the outer ghost; none reported | ghosts only |
| binary element diffusion | `he_diffusion = F`; call is guarded | `rho` is `intent(in)` |
| molecular carrier transport | `carrier_transport = F`; the operator moves `f_sp` only | `rho` is `intent(in)` |
| viscosity / conduction | `transport_active()` false; rows 2 and 3 only in any case | no |
| semi-implicit energy | `energy_semi_implicit.f90` writes `u(3,:)` only | no |
| ionization / composition | changes `f_sp`, `T`, `p`; `W(1,:)` passes through untouched | no |
| Newton frozen base | `set_base_fix` is called only from the `EXHALE_PTC` branch; `nfix_base = 0` here, all 500 cells are unknowns | n/a |

**Verdict H4: there is no mass source and no non-flux modification of `u(1,:)`
in this configuration.** The mass equation is a pure conservative flux
divergence, so a non-flat `F r^2` is `d(rho)/dt != 0`, full stop. The sign says
where: `R_1 < 0` (density rising) over `1.004 <= r <= 1.027` and `R_1 > 0`
(density falling) over `1.028 <= r <= 1.15`, i.e. mass is piling up at the foot
of the deficit and draining above it.

---

## 5. H3 -- is it time dependent? Yes, and it is the dominant term

The accepted state was reloaded and marched forward with the JFNK hand-off
removed (`Solver:` line deleted, `du_th [PLM,WENO3]: 1.0e9 1.0e-30`, so nothing
can stop the run), 20000 steps, `EXHALE_P54_TS=50`. Everything else -- planet,
`base.inp`, secondary ionization, molecular chemistry -- as in the accepted run.
`dt = 1.68 s` on the accepted state (`CFL = 0.6`; the minimum is taken at cell 1
throughout, whose `dr` and `c_s` barely move), so 20000 steps is `9.3` hours of
model time.

**The state moves at once, and it moves a lot.** Face flux in units of the
accepted state's `F0`, over the whole 20000 steps:

| radius | min | max | mean | peak to peak |
|---|---|---|---|---|
| 1.005 | 0.221 | 2.251 | 1.342 | 2.030 |
| 1.03 | 0.702 | 1.667 | 1.213 | 0.964 |
| 1.10 | 0.976 | 1.230 | 1.133 | 0.254 |
| 1.20 | 0.998 | 1.060 | 1.035 | 0.062 |

**It is one coherent, damped oscillation, not noise.** The autocorrelation of
each of the four signals peaks at essentially the same lag:

| radius | autocorrelation peak lag | in seconds | `r` at the peak |
|---|---|---|---|
| 1.005 | 2050 steps | 3446 | 0.38 |
| 1.03 | 2100 steps | 3530 | 0.37 |
| 1.10 | 2150 steps | 3614 | 0.20 |
| 1.20 | 2300 steps | 3866 | 0.20 |

A period of `3.4-3.9e3` s. For comparison, the sound round trip between the
base and `1.10 R_p` measured on this same state, `2 * sum(dr/c_s)`, is
`3288` s, and to `1.20 R_p` it is `5234` s. **The period is consistent with a
standing acoustic mode trapped between the base boundary and the sonic-transit
region**, which is the mode `docs/p44_base_sawtooth.md` and the CETIMB Shapiro
filter are about. That is a consistency statement from one measured period
against one measured crossing time, not an identification of the mode.

**Where it goes.** The amplitude decays through the run (correlation 0.37 at
the first repeat) and the layer settles by about step 12000 onto a **different
and still non-flat profile**, with the sign of the deviation reversed:

| | `F/F_outer` at 1.005 | 1.03 | 1.10 | 1.20 | 2.0 |
|---|---|---|---|---|---|
| the accepted state | 0.943 | **0.700** | 0.969 | 0.998 | 1.000 |
| after 20000 steps | 1.317 | **1.189** | 1.107 | 1.018 | 0.998 |

and the run is not converged there either: `du = 6.22e-03`,
`||R|| = 2.03e-05` (L-inf) / `3.08e-05` (volume weighted), flux spread over the
gate window `2.69e-02` cell centred and `2.66e-02` at the faces. **Both gates
would reject the state 20000 steps after they accepted its ancestor**, by a
factor 11 on the flux gate and 2 to 3 on the residual gate. The wind level also
drifted, `F_outer` `4.4642e-05 -> 4.5493e-05`, `+1.9` percent.

**Verdict H3: confirmed.** The accepted state is one phase of a damped
oscillation of the layer, and the acceptance happened at a phase in which the
gate window looked flat. The oscillation is a real time dependence of the
solution, not a measurement of it.

---

## 6. Verdict

1. **The non-flatness is physical to the discretization, not a diagnostic
   definition.** Cell-centred and Riemann-face spreads agree to within one to
   four percent at every window above `r = 1.002 R_p` (section 2). Reading (i)
   of item (AA) is not supported.
2. **The exception is cell 1**, where the centred product is `-21 F0` and the
   faces carry `+1.75` and `+1.92 F0`. Everything about the "spread over all
   cells is 27 and changes sign" statement is that one cell; the face flux is
   single signed over the whole domain.
3. **The continuity equation has no source and nothing else writes `u(1,:)`**
   in this configuration (section 4), so a non-flat `F r^2` is exactly
   `d(rho)/dt != 0`.
4. **The residual gate cannot see it, by construction.** The scaled continuity
   residual equals the fractional flux error per cell times the local Mach
   number, verified to one percent (section 3). At `Mach = 5e-5` a tolerance of
   `1e-5` admits a 19 percent flux change per cell.
5. **The state is genuinely non-stationary** (section 5), and 20000 further
   steps leave it further from both gates than it started.

So the accepted state is steady in the wind and unsteady in the layer, and both
gates -- the residual one by its scale, the flux one by its window -- are
measurements of the wind only.

---

## 7. Recommendation

**(a) Measure the flux at the faces, not at the cell centres.** It costs
nothing (`face_flux` is already assembled by `RK_rhs` for every residual
evaluation) and it is the quantity the scheme conserves, so it is the one a
steady state makes constant. Concretely it buys two things: the gate stops
reporting the cell-1 sawtooth as a flux (`-21 F0`), and the spread becomes a
well-defined single-signed number over the whole column instead of one that
changes sign. It does **not** buy a smaller number in the layer: the face
spread is 1 to 4 percent *larger* than the cell-centred one at every window.
This is the one place where switching the measure changes the answer, and it
changes it by removing an artifact, not by removing the problem.

**(b) Do not move `r_flux` down to 1.0 and expect the gate to hold.** Whatever
window is chosen, the layer of this state is not flat, and a whole-column flux
gate at the present `flux_spread_th = 5e-3` would reject the entire hot-Uranus
ladder including every rung P23 converged. The honest options are:

* keep `r_flux = 1.2` as the *acceptance* gate and **report the face-flux
  spread over `r >= 1.03` and `r >= 1.10` beside it**, which is (AA)'s own
  "cheapest next step" and costs one extra pass. On the accepted state that
  would have printed `3.14e-01` and `3.04e-02` next to `2.51e-03`, and the
  disagreement is the item;
* or gate the whole column on a threshold appropriate to a layer that is
  hydrostatic to `1e-4` -- but there is at present no measurement saying what
  such a threshold should be, and inventing one would certify the same states
  by a different route.

**(c) The residual gate needs a companion in the layer, not a smaller
tolerance.** Tightening `resid_th` does not help: at `Mach = 5e-5` even
`1e-9` would admit a 2 percent flux change per cell, and the wind rows would
stop converging long before the layer rows started to bite. What is missing is
a statement that does not carry the Mach factor. The flux flatness itself is
exactly that statement, which is the reason the flux gate exists; the finding
here is that its *window* excludes the region where the residual gate is weak,
so between them the two gates leave `1.0 <= r <= 1.15` uncovered by anything.

**(d) The oscillation is the thing to fix, not the measure.** Section 5 shows a
coherent mode of period `~3.5e3` s and an amplitude of order the wind flux
living in `1.0 <= r <= 1.15`. Until it is damped or resolved away, any single
state drawn from it will be accepted or rejected by its phase. Two existing
handles that address the mode rather than the measurement of it: the periodic
Shapiro filter (off by default; CETIMB uses it for exactly this) and base-cell
refinement, which (AA) records as repairing the `He/H = 0.3` rung. Neither was
tested here.

---

## 8. What was not established

* **The cause of the mode.** The period is consistent with a base-to-`1.1 R_p`
  acoustic round trip; it was not shown to be that mode, and no boundary-
  condition variant was tried. `docs/p44_base_sawtooth.md` reports that no BC
  variant moved `Mdot`, which is a different question from whether one damps the
  oscillation.
* **Grid dependence.** Everything here is on the 500-cell default grid
  (`base 2e-4 x 50`). (AA)'s refinement result -- `1e-4 x 100`, 550 cells --
  repairs the `He/H = 0.3` rung and spoils `He/H = 30`; that experiment was not
  repeated on this rung, and it is the single most informative missing
  measurement. It needs a cold start (88000-120000 steps on the two rungs that
  were run that way), which did not fit in this session.
* **Whether the state 20000 steps on ever converges.** It was still drifting
  (`du = 6.2e-3`) when the step cap stopped it. Whether the layer settles onto a
  flat flux given `1e5` steps, or onto a persistent non-flat one, is open.
* **The energy row's disagreement on reload** (`2.68e-05` re-measured against
  the `2.944e-06` the solver accepted, while the mass row reproduces to five
  digits). The likely cause is the precision of the species file the ionization
  solve restarts from. It was not chased; it does not enter any mass-flux
  number here, but it means a reloaded state is not residual-identical to the
  state the solver accepted, which is worth knowing before any gate is
  re-measured this way.
* **Other planets and rungs.** One state, one planet. The identity of section 3
  is algebraic and will hold anywhere; the size of the layer's flux deficit is a
  property of this state.

---

## 9. Files

```
docs/p54_base_layer_mass_flux.md    this note
docs/p54_base_layer_mass_flux.png   the figure: cell-centred vs face flux, the
                                    Mach identity, and the marching time series
docs/p54b_scale.png                 the same for the row-scale variants
docs/p54g23_row_scale_scan.md       rows 2 and 3 given the same measure
docs/p54g23_census.md               the 34 named cases under both measures
docs/p54g23_nonconvergent.md        which row and cell blocks each failure
```

Section 10 below is the P54b experiment the change came out of; sections 1-9
are the diagnosis it rests on.

The change these measurements led to is `Update_EXHALE.md` section 143. The
probe builds, the run directories (`h1/`, `h3long/`, `n1*/`, `w1*/`, `ts*/`,
`ws*/`, `rr_*/`) and the analysis scripts stay in the campaign scratch
(`.../scratchpad/p54flux/` and `.../scratchpad/p54g/`); nothing in them is
needed to read this note.

The source of the state, read only:
`.../scratchpad/ladder2/runs/heh0p0793cc_r1/`.

## 10. P54b -- scaling the continuity row by the face mass flux

**Read this section with its outcome in hand.** It is the experiment that led to
`Update_EXHALE.md` section 143, and section 143 went further than section 10.6
recommends: all three conservation rows are now scaled by their own largest
term, not the continuity row alone. What section 10 settles, and what section
143 rests on, is *which way* to do it -- the measure and not the solver -- and
that finding is 10.4. The three variants below are experiment labels, switched
in the campaign builds by an environment variable that is **not** in the code:
`A` is the measure section 143 replaced, `B` rescales measure and solver both,
`G` rescales the measure only. The figure is `docs/p54b_scale.png`.

Section 3 says the dimensionless continuity residual is the fractional flux
error per cell times the local Mach number, so the residual gate cannot see a
quasi-hydrostatic layer. P54b removes the Mach factor: it replaces the
continuity row's scale by the row's own largest term, the larger of the two
Riemann face mass fluxes the row differences,

```
   s_1(j) = max( |F_{j-1/2}| r_{j-1/2}^2 , |F_{j+1/2}| r_{j+1/2}^2 ) / dV_j ,
   dV_j   = ( r_{j+1/2}^3 - r_{j-1/2}^3 ) / 3      (the dV of RK_rhs)
```

with a floor `1e-6 * rho c_s/dr` as a zero guard -- unused on every state
measured, since the face mass flux is single signed throughout (section 2).
`R_1/s_1` is then the fractional flux error per cell, with no Mach factor.

The mass unknown's state scale `D_1` is derived from `s_1` rather than left at
`rho`, so that the relation the code states between the two, `s_k = D_k
(|v|+c_s)/dr`, keeps holding and the JFNK merit `||D^-1 F||_2` and the
acceptance test stay the same weighting of the same residual. Rows 2 and 3 are
untouched, and every measurement below confirms they are: their numbers are
bit-identical between the two builds on the same state.

Switched in the campaign build by `EXHALE_CONT_SCALE=face`, default off, the
unset build reproducing the section 3 numbers digit for digit
(`4.902541E-05 / 2.944050E-06` for the mass row on the accepted state). Neither
the variable nor this variant is in the code: what was carried in is the
measure-only form of 10.4, generalized to three rows, as section 143.

**Staleness guard.** The face fluxes are the ones `RK_rhs` assembled, cached by
`assemble_residual`; a scale request for a cell whose density is not bit-for-bit
the cached one belongs to a different state and falls back on the density
scale, counted and reported. On WASP-121b the count is **0**; on the hot Uranus
it is **3628** out of order 3e5 requests (about one percent), all from
line-search trial states, so about one percent of the hot-Uranus numbers below
are the density scale rather than the face scale. That is not what decides the
hot-Uranus outcome -- the blocking row is the same row and the same cells at
every one of the 179 iterations -- but it is a contamination and it is stated.

### 10.1 The measure itself, on states already accepted

Volume-weighted relative residual per row, re-measured on the written state
(`EXHALE_RESIDUAL=1`), with and without the new scale. Only row 1 changes.

| state | row | density scale | face scale | ratio |
|---|---|---|---|---|
| hot Uranus, ladder2 accepted | mass | 2.944e-06 | 5.584e-03 | 1900 |
| hot Uranus, A (density-scale Newton) | mass | 4.857e-06 | 1.957e-02 | 4030 |
| WASP-121b, A (density-scale Newton) | mass | 6.533e-06 | 3.770e-04 | 58 |
| any of them | momentum, energy | unchanged | unchanged | 1.000 |

The factor is the reciprocal Mach number of the cells that carry the residual,
and it is 70 times larger on the molecular hot Uranus than on the atomic hot
Jupiter because the hot Uranus's layer is that much more subsonic.

### 10.2 Does the solve converge on the new measure?

Both runs start from the same state and use the same procedure -- `Load IC`
from that state, `Solver: Newton 100.0`, `du_th [PLM,WENO3]: 1.0e9 1.0e-3`,
`Secondary_ionization: Immediate` -- and differ only in which variant the build carries.
The atomic case is `backup/regression/wasp_full_newton`'s own `input.inp` and
`metals.inp` (md5 `14930836...`, `21c1d02d...`), run as that case defines it:
a **cold start** (`Load IC? False`), so the two builds are compared on the whole
run and not only on a polish.

| | hot Uranus A: density | hot Uranus B: face | WASP-121b A: density | WASP-121b B: face |
|---|---|---|---|---|
| JFNK `info` | **0** | **2** | **0** | **0** |
| JFNK iterations | 5 | 179 (aborted) | 14 | 14 |
| residual evaluations | 115 | **8238** | 338 | **340** |
| `\|\|R\|\|` the solve reported | 2.031e-06 | 7.798e-04 | 7.398e-06 | 1.835e-06 |
| flux spread at the gate (`r >= 1.2`) | 2.106e-03 | **9.520e-03** | 1.357e-03 | **1.278e-04** |
| how it ended | both gates met | `NEITHER gate met`; no descent direction exists for the banded model | both gates met | both gates met |
| `log10 Mdot` | 10.34 | 10.33 | 13.21 | 13.22 |

**On the atomic hot Jupiter the change is free and it works.** Same iteration
count, 340 residual evaluations against 338 (+0.6 percent), the same mass-loss
rate to 0.01 dex, and the mass flux an order of magnitude flatter at the gate.

**On the molecular hot Uranus it does not converge.** The solve grinds from
2.23e-2 down to a plateau at 7.8e-4, spends most of its 179 iterations in
damped Gauss-Newton escapes, and aborts on `no descent direction exists for the
banded model at this state` with `||grad merit|| = 4.9e8`. It costs 72 times
the residual evaluations of the run that succeeds, and the state it hands back
is worse at the existing gate than the state it started from (9.52e-3 against
2.46e-3, i.e. it now fails the flux gate it entered passing).

**Which row and which cells block it.** The worst scaled residual is the MASS
row at every one of the 179 iterations, and it sits in the uniform base region:

| worst cell | iterations it was worst | radius |
|---|---|---|
| j = 16..22 | 129 of 179 | 1.003 - 1.004 |
| j = 1..5 | 24 | 1.000 - 1.001 |
| j = 13..15 | 14 | 1.003 |
| j = 224..242 | 5 | 1.13 - 1.16 |
| others | 7 | |

The reason is visible in the face fluxes of the state the solve is handed. On
the hot Uranus the base carries 7 to 9 times the wind flux in its first two
cells and then a smooth 7 percent rise over the next forty, so hundreds of
cells each hold a flux error of 1e-3 to 2.5e-3 of `F_0`:

| j | r | `F_in/F_0` | `F_out/F_0` | face-scaled `R_1/s_1` |
|---|---|---|---|---|
| 1 | 1.00019 | 9.067 | 7.122 | 2.14e-01 |
| 2 | 1.00039 | 7.122 | 1.016 | 8.57e-01 |
| 3 | 1.00058 | 1.016 | 1.017 | 1.21e-03 |
| 11 | 1.00212 | 1.020 | 1.021 | 9.20e-04 |
| 17 | 1.00328 | 1.028 | 1.030 | 1.95e-03 |
| 26 | 1.00502 | 1.048 | 1.051 | 2.48e-03 |
| 41 | 1.00791 | 1.086 | 1.088 | 1.95e-03 |

On WASP-121b, by contrast, the same state has its whole continuity imbalance in
**two** cells, and once the solve is asked to flatten the flux it does:

| j | r | A: `F_out/F_0` | A: `R_1/s_1` | B: `F_out/F_0` | B: `R_1/s_1` |
|---|---|---|---|---|---|
| 1 | 1.00020 | 1.1377 | 6.86e-03 | 1.0138 | 4.38e-03 |
| 2 | 1.00039 | 1.1086 | 2.56e-02 | 0.9999 | 1.37e-02 |
| 3 | 1.00059 | 1.1274 | 1.67e-02 | 1.0001 | 1.03e-04 |
| 5 | 1.00098 | 1.1257 | 4.64e-04 | 1.0000 | 2.28e-07 |
| 11 | 1.00216 | 1.1247 | 4.31e-05 | 1.0000 | 1.31e-07 |
| 21 | 1.00411 | 1.1247 | 1.49e-04 | 1.0000 | 1.67e-07 |
| 41 | 1.00803 | 1.1151 | 5.66e-04 | 1.0001 | 4.01e-08 |

From the third cell outward the WASP-121b B state carries **one** mass flux to
eight significant figures. The two cells that remain are the base ghost and its
neighbour, and they pass only because the norm is volume weighted and their
volume is negligible -- which is worth saying plainly: the new measure does not
solve the base cell either, it just stops the base cell from being the whole
column's excuse.

### 10.3 The layer, and the wind it feeds

Measured on the written states. Spreads are of the cell-centred `rho v r^2`
(the face numbers are within a few percent of them everywhere above 1.002 and
are given for the two windows that matter).

| state | v(1) [cm/s] | T(1) [K] | H2 front | spread `r>=1.03` | `r>=1.10` | `r>=1.20` | `F/F_0` at 1.03 | at 1.10 |
|---|---|---|---|---|---|---|---|---|
| hot U, ladder2 accepted (the input) | -234.3 | 1097.9 | 1.0419 | 3.105e-01 | 2.920e-02 | 2.460e-03 | 0.703 | 0.970 |
| hot U, A: density | -15.3 | 1111.5 | 1.0425 | 6.065e-02 | 1.589e-02 | 2.106e-03 | 1.061 | 1.016 |
| hot U, B: face | -259.3 | 1111.8 | 1.0419 | 8.191e-02 | 1.685e-02 | **9.523e-03** | 1.065 | 1.010 |
| WASP, A: density | +1555.6 | 2360.5 | -- | 5.469e-02 | 6.952e-03 | 1.357e-03 | 1.055 | 1.007 |
| WASP, B: face | +1373.3 | 2356.8 | -- | **1.538e-04** | **1.278e-04** | **1.278e-04** | 1.000 | 1.000 |

Face-flux spreads on the same states: WASP A `5.44e-02` at `r >= 1.03`, WASP B
`9.17e-05` -- a factor **593**; over the whole column including the base cells,
`1.34e-01` against `1.38e-02`, a factor 10. The hot Uranus moves by less than a
factor 1.4 either way and in the wrong direction at the gate.

Neither the H2 front (1.0425 against 1.0419) nor the base temperature (1111.5
against 1111.8 K) nor the rate (10.34 against 10.33; 13.21 against 13.22)
distinguishes the two builds. What moves is the base velocity, `-15` cm/s on
the state the density scale converges to against `-259` cm/s on the one the
face scale stalls at, and on WASP-121b `+1556` against `+1373` cm/s.

### 10.4 Separating the measure from the solver: it was the solver

Variant B rescales two things at once -- what the solve is JUDGED by
(`residual_row_scale`) and what the Newton system and the line-search merit are
SCALED by (`state_scales_of_cell`, through `cell_state_scales`). The hot-Uranus
failure could be either: the state cannot be found, or this solver cannot find
it. A third build separates them: `G` puts the face flux
in the measure alone and leaves the Newton system and the merit on the density
scale, so the iterates are the ones variant A produces and only the acceptance
test changes.

**It converges, on both cases, and it gives the flattest mass flux of the
three.**

| | A: density | B: face, measure and solver | G: face, measure only |
|---|---|---|---|
| **hot Uranus** | | | |
| JFNK `info` | 0 | **2** | **0** |
| iterations | 5 | 179 | 9 |
| residual evaluations | 115 | 8238 | **206** |
| `\|\|R\|\|` reported (its own measure) | 2.031e-06 | 7.798e-04 | 8.636e-06 |
| flux spread at the gate `r >= 1.2` | 2.106e-03 | 9.520e-03 | **2.489e-04** |
| cache misses | -- | 3628 | **0** |
| `log10 Mdot` | 10.34 | 10.33 | 10.34 |
| **WASP-121b** | | | |
| JFNK `info` | 0 | 0 | 0 |
| iterations | 14 | 14 | 17 |
| residual evaluations | 338 | 340 | **406** |
| `\|\|R\|\|` reported | 7.398e-06 | 1.835e-06 | 3.138e-06 |
| flux spread at the gate | 1.357e-03 | 1.278e-04 | **5.627e-05** |
| cache misses | -- | 0 | **0** |
| `log10 Mdot` | 13.21 | 13.22 | 13.21 |

So the hot Uranus's layer **can** be made flat: the state exists, the existing
solver reaches it in nine iterations at 1.8 times the baseline cost, and the
thing that stopped variant B was rescaling the Newton system by a quantity that
varies by a factor 9 across the first two cells, not the measure.

**What the flat state looks like.** Face-flux spread `(max-min)/|mean|` of the
conserved face flux on the written states:

| window | hot U input | hot U A | hot U **G** | WASP A | WASP **G** |
|---|---|---|---|---|---|
| `r >= 1.01` | -- | -- | **1.043e-04** | -- | -- |
| `r >= 1.03` | 3.137e-01 | 5.648e-02 | **1.043e-04** | 5.438e-02 | **2.326e-05** |
| `r >= 1.10` | 3.041e-02 | 1.354e-02 | 5.651e-05 | 6.885e-03 | 1.725e-05 |
| `r >= 1.20` | 2.513e-03 | 1.880e-03 | **3.778e-06** | 1.312e-03 | **3.990e-06** |
| whole column | 1.316e+00 | 5.906e+00 | 6.056e+00 | 1.339e-01 | 1.322e-02 |

From `r = 1.01 R_p` outward the hot Uranus now carries **one** mass flux to one
part in `1e4`, where the state the ladder accepted was 30 percent out at 1.03.
The whole-column number is unchanged because the first two cells are unchanged:
the base ghost and its neighbour still carry 6 to 9 times the wind flux, and no
variant here touches them.

Nothing physical moves with it. `log10 Mdot` is 10.34 for both A and G and
13.21 for both A and G; the H2 front moves from 1.0425 to 1.0449 `R_p`, the
base temperature from 1111.5 to 1116.9 K, and the base velocity is -15.25
against -15.09 cm/s. **The 30 percent flux deficit was not carrying any of the
answer -- it was a residual the old measure could not see, and removing it costs
nine iterations and changes no reported quantity.**

### 10.5 Is the flat state a steady state of the marching?

Each solved state was reloaded and marched with the JFNK hand-off removed
(`Solver:` deleted, `du_th [PLM,WENO3]: 1.0e9 1.0e-30`, `EXHALE_P54_TS=50`),
the same procedure section 5 used.

**On WASP-121b the answer is that the marching does not preserve any of them,
and the scale makes no difference to that.** Face flux in units of the solved
state's own `F_0`:

| | A: density | B: face |
|---|---|---|
| peak-to-peak at `r = 1.005`, first 200 samples | 7.90 | 7.88 |
| at `r = 1.03` | 4.26 | 4.31 |
| at `r = 1.10` | 1.17 | 1.19 |
| largest difference between the two at any sample, `r = 1.03` | 0.082 | 0.082 |
| state after 20000 steps: `v(1)` | 2153.2 cm/s | 2154.0 cm/s |
| spread `r >= 1.03` | 1.254e-01 | 1.247e-01 |
| spread `r >= 1.20` | 2.972e-02 | 2.963e-02 |
| `log10 Mdot` | 13.35 | 13.35 |

Both runs leave their solved state within a hundred steps, swing to seven times
the wind flux at `r = 1.005` around step 1250, and after 20000 steps land on the
**same** state to three digits -- a state 22 times less flat at the gate window
than either solve produced, and 0.14 dex higher in the rate. The two series
track each other to about one percent of excursions of order four, so the
residual scale plays no part in it: **the marching operator has an attractor of
its own, and it is not the JFNK solution.** That is a pre-existing property of
the hand-off (it is the base mode `docs/p44_base_sawtooth.md` and the HD 189733 b
note describe), measured here, not something P54b introduces -- but it does mean
that "converged" has two different meanings in this code and they do not agree
on this case.

**On the hot Uranus the oscillation is halved and not removed.** Compared over
the first 8000 steps, which all three runs cover and which is about four periods
of the section 5 mode:

| marched from | p2p at `r=1.005` | at 1.03 | rms at 1.03 | at 1.10 | rms at 1.10 | at 1.20 |
|---|---|---|---|---|---|---|
| A: density scale | 2.386 | 1.126 | 0.137 | 0.136 | 0.019 | 0.060 |
| B: face, both | 0.535 | 0.433 | 0.068 | 0.086 | 0.012 | 0.015 |
| **G: face, measure only** | 1.601 | 0.508 | 0.081 | 0.098 | 0.014 | 0.018 |

The flat state is quieter -- half the excursion at 1.03 and a quarter at 1.20 --
but it still leaves itself, and after 20000 steps the A and B runs land within
0.03 percent of each other in `v(1)` (-251.36 against -251.29 cm/s) on a state
whose spread at `r >= 1.03` is 6.2e-2 either way. **The same attractor, again.**

So the answer to "is the flat state a true steady state" is: it is a steady
state of the residual the solver solves, and it is not a fixed point of the
marching loop. The two differ by everything that is operator split out of the
residual, and that difference is larger than the improvement P54b makes. Which
of the two is the physical steady state is not settled by anything measured
here.

### 10.6 Recommendation, and what was adopted

**Adopted 2026-09-03 as `Update_EXHALE.md` section 143, and carried further than
this recommendation asks.** The recommendation below is the measure-only form
for the continuity row; what went in is that form for all three conservation
rows (`mass_flux_row_scale`, `momentum_row_scale`, `energy_row_scale` in
`steady_residual.f90`), after `docs/p55_base_mode.md` section 4 measured the
momentum and energy rows to be mis-scaled the same way. The one thing this
section settles and section 143 keeps verbatim is the second half of the
sentence: the measure, not the solver. The extension to the other two rows is
measured in `docs/p54g23_row_scale_scan.md` and its case census in
`docs/p54g23_census.md`.

**Adopt the face mass flux as the continuity row's SCALE in the convergence
measure, and do NOT rescale the Newton system with it.** On the two cases tested this is the only variant that converges on both,
it is the flattest of the three on both, it costs 1.8 and 1.2 times the
baseline's residual evaluations, and it moves no reported physical quantity
(`log10 Mdot` 10.34 and 13.21, unchanged from the baseline to the digit
printed). It removes the section 3 blindness at its root: the measure no longer
carries a Mach factor, so a layer at Mach 5e-5 is held to the same fractional
flux error as the wind.

What that buys concretely on the case item (AA) was opened about: the 1 microbar
hot-Uranus rung goes from a mass flux 30 percent out at `r = 1.03` to one flat
to `1.0e-4` from `r = 1.01` outward, in nine JFNK iterations, with the H2 front
at 1.0449 instead of 1.0425 `R_p` and the same mass-loss rate.

Two things this does NOT do, and neither should be claimed for it:

* it does not touch the first two cells, whose face flux is still 6 to 9 times
  the wind flux; they pass because the norm is volume weighted;
* it does not make the solved state a fixed point of the marching loop
  (section 10.5). The oscillation is halved, not removed.

If it is carried in, the flux gate `flux_spread_th` and `r_flux` should be
revisited in the same change: with this measure the solved states are flat to
`1e-4` over `r >= 1.03`, so the existing `5e-3` at `r >= 1.2` stops being the
binding constraint and the window can start much lower. **Every golden will
move**, because the solve now stops at a different state -- a physically better
one by the criterion the gate encodes, but a different one -- so this is a
deliberate golden refresh, not a byte-identical change.

*What actually happened when it went in* (section 143.4): of the 34 named cases,
**seven goldens moved and twenty-seven were byte-identical**, because a run that
stops on `du` or on `maxsteps` never reads a row scale. The default ten-case
matrix was 10/10 identical and needed no refresh at all. The flux gate's own
defaults were left alone; the window question is section 143.5 and still open.

### 10.7 What P54b did not establish

* **Only two cases here.** One molecular hot Uranus and one atomic hot Jupiter;
  the regression matrix could not be run while this was measured, because
  `backup/regression/` was in use by another worker. *Closed afterwards*: all 34
  named cases were run under the adopted measure and the one it replaced, and
  the table is `docs/p54g23_census.md`.
* **Why variant B destroys the solve.** The evidence is that it does (179
  iterations, no descent direction, `||grad merit|| = 4.9e8`) and that removing
  the rescaling of `state_scales_of_cell` fixes it. The mechanism -- presumably
  that `D_1` then varies by a factor 9 across the first two cells and
  conditions the scaled Jacobian badly -- was not measured.
* **The staleness fallback** fired 3628 times in variant B on the hot Uranus
  and zero times everywhere else, including in the variant that is recommended.
  The B numbers are therefore about one percent contaminated. G and the WASP
  runs are clean.
* **The marching attractor.** Sections 10.5 shows the JFNK solution and the
  marching fixed point are different objects on both planets, by more than the
  difference P54b makes. That is the larger open question and it is untouched
  here.
* `tsG` was stopped at 10600 steps rather than 20000 (the machine was carrying
  another worker's regression matrix all session); the comparison in 10.5 is
  over the 8000 steps all three cover.
