# The complete update-map diagnostic, and what it decides

Review of 2026-09-03, B2: "Expose a diagnostic that applies one production step
without stop logic or file output and returns `G_dt`. Compare it with the steady
residual at the same state for at least three decreasing time steps. Break down
differences after the Euler, chemistry, energy, carrier, diffusion, boundary,
and filter stages. This test decides item 3."

**It decides it: `G_dt` approaches the steady residual at first order in `dt`,
on both planets. The production map and the residual share a fixed point, and
the difference at a CFL step is the ordinary truncation error of an explicit
step, not an operator mismatch** -- unless the Shapiro filter is on, in which
case the two provably cannot share one, and the diagnostic shows it diverging.

## What it is

`EXHALE_UPDATE_MAP=<f1,f2,...>` (env, the convention every other diagnostic in
this code follows; `=1` means the default triple 1, 0.1, 0.01). The run takes
one production step per requested factor, at `dt = f x dt_CFL`, restoring the
state between them, and writes

```
   G_dt(q) = [ Phi_dt(q) - q ] / dt          against       R(q)
```

per row and per cell, with the increment of every operator-split stage beside
it, to `output/update_map.txt` and as a summary to stdout.

**`Phi_dt` is the marching loop's own body, not a copy of it.** The diagnostic
sets the stop thresholds out of the way, caps the loop, restores the state at
the top of each step and snapshots `u` at the stage boundaries; the step that
runs between the snapshots is the step the run takes. A second implementation of
the step would not have been able to answer the question it is asked.

Two things it fixes by construction, both of which produced wrong answers while
this note was being written and are recorded so they are not repeated:

* the state must be restored BEFORE `u_old = u`, because the positivity retry
  loop restarts every attempt from `u_old`; a restore after that line is undone.
* the run must be held in ONE reconstruction. A `PLM+WENO3` run takes its first
  step in PLM and switches at the end of it, so the three time steps would
  otherwise be compared across two different spatial operators.

## The result

Hot Uranus, the section 144.5 root, `max` over cells in code units of `1/t_s`:

| `dt / dt_CFL` | `|R|` mass | `|G_dt|` mass | `|G_dt + R|` mass | mom | energy |
|---|---|---|---|---|---|
| 1 | 1.568 | 0.9164 | **0.6521** | 0.5332 | 0.9604 |
| 0.1 | 1.568 | 1.479 | **0.0893** | 0.0527 | 0.1400 |
| 0.01 | 1.568 | 1.559 | **0.00923** | 0.00517 | 0.0234 |

WASP-121b, the same:

| `dt / dt_CFL` | `|G_dt + R|` mass | mom | energy |
|---|---|---|---|
| 1 | 9.63e-02 | 1.33e-01 | 1.93e-01 |
| 0.1 | 1.09e-02 | 1.69e-02 | 2.45e-02 |
| 0.01 | 1.10e-03 | 1.73e-03 | 2.34e-01 |

The mass and momentum rows fall by a clean factor ten per decade of `dt` on both
planets -- first order, which is what a consistent explicit step gives. **The
42 percent defect section 144.1 measured at a CFL step is that truncation
error**, and it is not evidence of a different operator.

The energy row converges on the hot Uranus (0.960, 0.140, 0.0234) and does not
converge monotonically on WASP-121b, where the composition pressure projection
contributes a RATE of 1.124 against a residual of 1.133 -- the same size as the
whole row. That is the one stage the steady residual has no counterpart for
(section 144.3), and on WASP-121b it is not small. It is the open end of this
measurement.

## The stage breakdown

Hot Uranus, per stage, as a rate (`max` over cells, `1/t_s`):

| stage | `dt/dt_CFL = 1` | 0.1 | 0.01 | what that means |
|---|---|---|---|---|
| hydro (SSP-RK3 of `dF - S`) | 1.359 | 2.183 | 2.301 | converging to `|R|` |
| chem (sweep + pressure re-derived at the post-sweep particle count) | 5.80e-03 | 9.38e-03 | 9.86e-03 | a rate, converging |
| energy (semi-implicit heating-cooling) | 5.092e-03 | 5.091e-03 | 5.091e-03 | a rate, exactly |
| transport (viscosity, conduction) | 0 | 0 | 0 | off |
| filter (Shapiro) | 0 | 0 | 0 | off |
| of which `Apply_BC` | 6.8e-12 | 9.1e-11 | 7.6e-10 | `dU` of order `1e-15`: the conservative-primitive round trip is round-off |

The `Apply_BC` row answers review item 10 for this state directly: the
round trip the boundary fill performs moves the INTERIOR by a relative
**4.2e-16** at most. It is not a source of the motion; the review's B1
recommendation to preserve the interior bytes remains right on its own terms
(it removes an unnecessary operation), but nothing measurable rests on it here.

## With the Shapiro filter on, the two cannot share a fixed point, and here is
## the number

Hot Uranus, `Shapiro filter: 0.5 1` (every step), everything else identical:

| `dt / dt_CFL` | filter stage, as a rate | `|G_dt + R|` mass |
|---|---|---|
| 1 | 1.187e+01 | 1.185e+01 |
| 0.1 | 1.189e+02 | 1.189e+02 |
| 0.01 | 1.188e+03 | 1.188e+03 |

The filter's contribution to `dU` is **the same at every time step** -- it is
applied per STEP, not per unit time -- so as a rate it grows as `1/dt` without
bound, and `G_dt` diverges from `-R` instead of converging to it. The residual
contains no such term and cannot. **A filtered marching map and an unfiltered
residual do not share a fixed point, and the divergence is exactly the filter's
own increment.** Either the filter's steady operator goes into the residual, or
the filter is prohibited in steady-root validation; the diagnostic now measures
which of those a given run needs.

## What this does not settle

* One state per planet, both with `Secondary_ionization` already matched to the
  file (section 144.2) and with transport and diffusion off, so the `carrier`
  and `diffusion` stages are measured as zero rather than as small.
* The energy row on WASP-121b, where the composition pressure projection is a
  rate as large as the residual. Whether the residual should carry that term or
  the step should not is the question section 144.3 left open, and this
  measurement sharpens it without answering it.
* The diagnostic restores `u` and `f_sp` between steps but not the module-level
  ledgers the sweep accumulates; those are counters, and no measured quantity
  reads them.

## Where the snapshots may sit

The snapshots were first taken between the STATEMENTS of one RK3 stage, around
the `Apply_BC` call inside it. Every one of them was inside
`if (upmap_n .gt. 0)` and did nothing with the diagnostic off, and the ten
regression cases all failed at a relative `3e-14` to `9e-11`. `a - b*c` is a
fused multiply-add candidate and `-ffp-contract=fast` is gfortran's default, so
an inserted statement changes the scheduling and with it whether the FMA is
taken; an FMA rounds once where a multiply and a subtract round twice, and one
ulp per step becomes `1e-11` over eleven thousand. The boundary contribution is
now taken once, out of line, in `update_map_begin_step`, and every remaining
snapshot sits between operator-split stages and never between the statements of
one. Anything added to this diagnostic must keep that rule. Section 145.8.
