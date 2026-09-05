# The coupled carrier solve after the characteristic base boundary (2026-09-04)

**Judgement first. The coupled solve still does not converge, and the
obstruction has moved.** With the carrier as a fourth Newton unknown the JFNK
returns `info = 2` on all three attempts, all three stagnating with no
improvement of the best iterate in 40 iterations. The floor is no longer
`||R|| = 1.00` in base cell `j = 1`: in **111 of the 120 JFNK iterations the
worst row is the CARRIER row `k = 4`, in cells `j = 221-240`, `r = 1.121-1.154`**
-- the H2 front. Only 5 iterations are worst in the energy row and 4 in the mass
row. At the hand-off state the base momentum row reads `mom/grav = 4.852E-04`
against `energy = 1.440E+00` in the same print, so section 152 did remove the
momentum imbalance it was built to remove; what it did not remove is the
molecular layer.

**And a control arm says the step control is not the whole story.** The same
hand-off state solved by the Picard alternation instead (no carrier unknown at
all) also returns `info = 2` on all three attempts, stagnating at
`||R|| = 3.705E-01` with its worst cell in the **energy** row at `r = 1.0117`
to `r = 1.0494` -- inside the same molecular layer -- and leaving the carrier
row at `7.450E-03` against its target `1.0E-03`. Neither route reaches a steady
state of the molecular layer, and the two fail in different rows of the same
place.

Because the solve does not converge, **the P23 re-verdict below is taken on a
MARCHED state, not on an accepted steady state**, and is labelled as such
throughout.

Binary: the tree's `EXHALE.x`, md5 `063b2673eb0d4f1ca342cc07bf942b4c`, copied
into scratch and run there. `OMP_NUM_THREADS=8`. Nothing in the repository was
written.

---

## 1. What was run

All three runs are the `mol_carrier` regression case (hot Uranus, 1 microbar
base level from `base.inp`, `q_H2_base = 0.84`), copied into scratch, with
`Molecular carrier transport: True` and:

| key | marching | coupled finish | Picard finish |
|---|---|---|---|
| `Secondary_ionization` | `Immediate` | `Immediate` | `Immediate` |
| `Stall [tol,N]` | `0.0 2000` | `0.0 2000` | `0.0 2000` |
| `Solver` | `Newton 0.0` (hand-off disabled) | `Newton 1.0e2` | `Newton 1.0e2` |
| `Load IC?` | `False` | `True` | `True` |
| `Coupled carrier solve` | absent | `True` | absent |
| `EXHALE_MAXSTEPS` | 200000 | none | 7000 |

The two finishes are fed the SAME state: the 200,000-step marched
`Hydro_ioniz.txt` / `Ion_species.txt`, copied in as `*_IC.txt`. They differ in
one input line.

## 2. The marching trajectory

200,000 steps in 1 h 37 m. Final `du = 5.8899E-03`, `dtu = 6.7832E-07`,
`log10 Mdot = 10.43`. Every ionization-equilibrium root converged
(`hybrd1 info: 1 100800002`), 2 constrained-continuation roots, base inflow
stayed subsonic at `max Mach 1.633E-01`.

| step | flux spread `r>=1.20` (gate) | `r>=1.10` | `r>=1.03` | `\|\|R\|\|`(ref) |
|---:|---|---|---|---|
| 41,500 | 3.194E-01 | 5.966E-01 | 8.834E-01 | 1.030 |
| 62,000 | 2.948E-01 | 5.758E-01 | 8.634E-01 | 0.9745 |
| 88,000 | 2.513E-01 | 5.363E-01 | 8.248E-01 | 0.8793 |
| 127,500 | 1.803E-01 | 4.849E-01 | 7.749E-01 | 0.7982 |
| 167,500 | 1.274E-01 | 4.388E-01 | 7.279E-01 | 0.7518 |
| 187,500 | 1.097E-01 | 4.183E-01 | 7.066E-01 | 0.7344 |
| 200,000 | 9.696E-02 | 4.027E-01 | 6.930E-01 | -- |

**This differs from the earlier campaign in one way that matters.** In
`p23_transport_on_state.md` section 8.1 the gate window fell to a minimum of
`1.502E-02` near 200,000 steps and then rose by a factor 4.7. Here it is
**monotone decreasing through 200,000 steps and has not reached a minimum**.
The two runs are not the same experiment -- this one has
`Secondary_ionization: Immediate` and the section 152 boundary -- so the
comparison says only that the earlier turnaround does not reproduce here, not
why. A consequence for the work: **the hand-off state used below is not
demonstrably the flattest state this march can reach**, and marching further
was not tested.

## 3. The coupled solve: three attempts, three stagnations

`Coupled carrier solve: True`, JFNK hand-off fired at steps 2, 2002 and 4002
(2000 marching steps between attempts, the `N_stall` hold).

| | attempt 1 | attempt 2 | attempt 3 |
|---|---|---|---|
| entry `\|\|R\|\|` | 1.440 | 1.120 | 1.071 |
| entry worst cell | `j=1  r=1.0002  k=3` | `j=93  r=1.0206  k=3` | `j=95  r=1.0213  k=3` |
| entry `mass / mom-grav / energy` below `r_esc` | 1.508E-01 / **4.852E-04** / 1.440 | 2.354E-02 / 7.616E-03 / 1.120 | 2.415E-02 / 6.978E-03 / 1.071 |
| returned `\|\|R\|\|` (own composition) | 1.282 | 1.064 | 1.065 |
| flux spread, gate window | 1.104E-01 | 9.319E-02 | 9.514E-02 |
| carrier row, volume-weighted | 5.45E-03 | 6.43E-03 | 6.37E-03 |
| carrier row, worst cell | 2.00E-01 @ `j=231` | 2.23E-01 @ `j=229` | 2.47E-01 @ `j=229` |
| trials refused, element budget | 102 | 78 | 93 |
| worst reaction residual of a refused trial | 7.403E+06 | -- | 1.548E+00 |
| residual evaluations | 2081 | -- | 1794 |
| stop | STAGNATED, 40 it | STAGNATED, 40 it | STAGNATED, 40 it |
| `info` | 2 | 2 | 2 |

**All three gates fail, and by wide margins:**

```
gate NOT met: residual ||R|| = 1.065E+00 >= 1.000E-05
gate NOT met: flux spread   9.514E-02 >= 2.000E-05
gate NOT met: carrier row   6.367E-03 >= 1.000E-03
```

The residual is the section 155 self-consistent one -- the print says
`residual at the state's own composition, 10 sweeps` -- so the number is
measured at the composition the written state carries, which is now the
default.

### 3.1 Where it is blocked

Worst-row histogram over all 120 JFNK iterations of the three attempts:

| row | count | cells |
|---|---:|---|
| `k = 4` carrier | **111** | `j = 221-240`, `r = 1.121-1.154` |
| `k = 3` energy | 5 | -- |
| `k = 1` mass | 4 | -- |
| `k = 2` momentum | 0 | -- |

The most frequent single cells are `j = 230` (`r = 1.135`, 15 iterations),
`j = 233` (`r = 1.141`, 12), `j = 229` (`r = 1.134`, 11), `j = 231`
(`r = 1.137`, 10), `j = 224` (`r = 1.126`, 10). **That band is the H2 front**:
in the same state `f(H2)` crosses `1e-2` at `r = 1.1333` and `1e-4` at
`r = 1.2897`, and `n(H3+)` last exceeds `1 cm^-3` at `r = 1.1389`.

The handed-back states put their worst cell in a hydrodynamic row instead, and
in two different places:

| attempt | worst row and cell | signed terms |
|---|---|---|
| 1 | energy, `j=175`, `r=1.06649` | `R = -3.595E-04`, flux div `-9.844E-05`, heat `2.805E-04`, cool `1.940E-05`, scale `2.805E-04` |
| 3 | energy, `j=2`, `r=1.00039` | `R = 7.018E-02`, flux div `-2.386E-02`, heat `8.704E-04`, cool `5.128E-03`, scale `6.592E-02` |

Attempt 1's worst cell, `r = 1.0665`, is the `f(H2) = 0.5` front of that same
state (`r = 1.0656`) to within four cells.

### 3.2 Against the previous number

The obstruction that sections 153 and 154 measured was `||R|| = 1.00` in base
cell `j = 1`, energy row, with the base momentum row contributing
`7.24E-03`/`1.21E-01` to the below-escape ledger. In the entry ledger of this
run's attempt 1, at the same print and the same normalization:

```
below r_esc [1:391] |R|: mass= 1.508E-01  mom/grav= 4.852E-04  energy= 1.440E+00
```

**The momentum row at the base is now three digits below the energy row.** That
is section 152 doing what it was built to do. The energy row at cell 1 is not
fixed by it -- it is the entry worst cell of attempt 1 -- but it is no longer
where the solve spends its iterations. In attempts 2 and 3 the entry worst cell
has already moved off the boundary, to `j = 93` / `j = 95` at `r ~ 1.021`.

The window decomposition of the returned state (attempt 3) puts the remaining
weight in two places, not one:

```
energy   r<1.03    :  5.244E-05 / 1.042E-04   (0.50 of its own scale)
energy   1.10-1.20 :  4.514E-05 / 1.180E-04   (0.38)
mass     1.10-1.20 :  3.961E-05 / 1.927E-03
momentum r<1.03    :  6.092E-05 / 3.648E-01   (1.7E-04)
```

## 4. Control arm: the Picard alternation on the same state

Same hand-off state, same input except the `Coupled carrier solve` line
removed, so the carrier is relaxed at fixed wind instead of solved.

| | attempt 1 | attempt 2 | attempt 3 |
|---|---|---|---|
| entry `\|\|R\|\|` | 1.440 | 3.831E-01 | -- |
| returned `\|\|R\|\|` | 5.928E-01 | 3.831E-01 | 3.705E-01 |
| flux spread, gate window | 9.222E-03 | 2.609E-02 | 2.521E-02 |
| worst cell of returned state | energy `j=61  r=1.0119` | energy `j=60  r=1.0117` | energy `j=153  r=1.04944` |
| `info` | 2 | 2 | 2 |

Carrier row of the final state: `7.450E-03` against `Carrier resid tol`
`1.000E-03`. Marched state after the three attempts: `log10 Mdot = 10.39`.

**The comparison is the useful part.** On the same state the Picard route
reaches `||R|| = 0.37` and a gate window of `2.5E-02`, against the coupled
route's `1.065` and `9.5E-02` -- roughly 3x and 4x better. Making the carrier a
Newton unknown makes both hydrodynamic gates worse, and it does so by moving
the iteration into the front, where 111 of 120 steps are then spent. But the
Picard route does not converge either, and its own worst cell sits in the
energy row of the same molecular layer (`r = 1.012` to `1.049`), with the
carrier row 7.5x its target. **So the molecular layer is not a steady state of
this system under either route, and that is not a step-control fact.**

## 5. The front still creeps, at the same rate

Snapshots taken from the marching run's `output/` (rewritten every 1000 steps,
so each row lags its label by up to 1000 steps). Front definition as in
`p23_transport_on_state.md` section 9: `f(H2) = 2 n(H2)/n_H`, log-linear
interpolation of the first downward crossing.

| label step | `f = 0.5` | `f = 1e-2` | `f = 1e-4` | max r with `n(H3+) > 1` | gate window |
|---:|---|---|---|---|---|
| 40,000 | 1.0586 | 1.1054 | 1.1494 | 1.1063 | 3.149E-01 |
| 80,000 | 1.0602 | 1.1108 | 1.1735 | 1.1133 | 2.606E-01 |
| 120,000 | 1.0618 | 1.1170 | 1.2886 | 1.1193 | 1.899E-01 |
| 160,000 | 1.0637 | 1.1243 | 1.2910 | 1.1288 | 1.315E-01 |
| 200,000 | 1.0656 | 1.1333 | 1.2897 | 1.1389 | 9.696E-02 |

`f = 0.5` advances 0.0070 over 160,000 steps, `0.0018` per 40,000 steps.
`p23_transport_on_state.md` section 8.1 gives `1.0485 -> 1.0562` over the same
160,000 steps, `0.0019` per 40,000. **The rate is unchanged to two digits by
the characteristic boundary**, which is what section 9 predicted when it showed
the front does not respond to freezing the base cells. The `f = 1e-4` column
jumps between 80,000 and 120,000 because the outer tail of `f(H2)` flattens and
the crossing moves to a different shoulder; the two inner thresholds are the
ones to read.

## 6. P23, re-taken on the transport-on marched state

**Caveat stated first: this is the 200,000-step marched state, `du = 5.9E-03`,
gate window `9.7E-02`. It is not an accepted steady state, and no run in this
campaign produced one.** The numbers are quoted so the transport-on layer can
be compared against the block-I accepted state of `p23_thermal_budget.md`
section 13.1 and against Koskinen et al. (2022) Model A, not as a converged
result.

`r/r_base` is the model's own radius in `R_p`, the base being at `r = 1`.

| `r/r_base` | 1.00 | 1.05 | 1.10 | 1.15 | 1.20 | 1.30 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|---|
| **T [K], this state** | 1114 | 1091 | 867 | 757 | 901 | 1560 | 2389 | 2697 | 2253 |
| T, block I (sec. 13.1) | 1117 | 630 | 593 | 986 | 1451 | 2117 | 2684 | 2686 | 2164 |
| T, Koskinen Model A | 1082 | 1180 | 1566 | 1717 | 1901 | 2267 | 2871 | 4120 | 4742 |
| **deficit K22 - this** | -32 | +89 | +699 | +960 | +1000 | +707 | +482 | +1423 | +2489 |
| **f(H2), this state** | 0.949 | 0.622 | 0.215 | 9.11e-4 | 2.75e-4 | 8.62e-5 | 3.48e-5 | 1.13e-5 | 4.44e-6 |
| f(H2), block I | 0.958 | 0.393 | 3.12e-3 | 2.62e-4 | 7.84e-5 | 3.95e-5 | 1.98e-5 | 5.08e-6 | 1.33e-6 |
| f(H2), Koskinen A | 0.984 | 0.982 | 0.976 | 0.963 | 0.936 | 0.843 | 0.688 | 0.508 | 0.403 |
| **q(H2), this state** | 0.785 | 0.405 | 0.111 | 4.22e-4 | 1.27e-4 | 3.99e-5 | 1.61e-5 | 5.21e-6 | 2.06e-6 |
| **n(H3+) [cm^-3], this state** | 1.33e5 | 1.72e3 | 1.70e2 | 1.25e-1 | 2.86e-2 | 7.81e-3 | 1.14e-3 | 3.51e-5 | 2.16e-6 |
| n(H3+), block I | 1.40e5 | 7.87e2 | 1.21 | 3.14e-2 | 8.08e-3 | 2.41e-3 | 3.26e-4 | 9.95e-6 | 5.20e-7 |
| n(H3+), Koskinen A | 1.90e4 | 3.14e3 | 4.99e3 | 7.02e3 | 9.10e3 | 1.16e4 | 1.06e4 | 7.16e3 | 1.88e3 |
| **n_e [cm^-3], this state** | 8.06e5 | 2.75e7 | 4.96e7 | 6.89e7 | 7.19e7 | 6.54e7 | 6.41e7 | 4.63e7 | 2.01e7 |

H2 front of this state, three levels: `f = 0.5` at `1.0656`, `f = 1e-2` at
`1.1333`, `f = 1e-4` at `1.2897`. `log10 Mdot = 10.43` (marched; the Picard arm
after 7000 further steps gives 10.39, and neither is a Newton-grade number).

**The P23 verdict does not change.** Transport-on marching keeps H2 alive
noticeably further out than the block-I accepted state -- `f(H2) = 0.215` at
`1.10 r_base` against block I's `3.1e-3`, a factor 70 -- and it lifts `n(H3+)`
at that radius by a factor 140. Both are real improvements in the right
direction. They are nowhere near enough: Koskinen Model A still has
`f(H2) = 0.976` there, and by `1.15 r_base` this state has fallen to `9e-4`
while Model A is still at `0.963`. The temperature deficit is unchanged in
character -- 700 to 1000 K through the molecular layer, growing to 2500 K at
`3 r_base` -- and the `n(H3+)` profile still has the wrong SHAPE: ours peaks in
the first cells and falls monotonically, Model A rises outward to `1.16e4` at
`1.30` and holds `1.9e3` at `3.00`. The attribution in section 7 of
`p23_thermal_budget.md` -- thermal starvation for want of the absorber -- is
untouched by this campaign.

## 7. Recommendation on sections 153 and 154

**Re-measure 153 (the scaled trust region) on this hand-off state. Do not adopt
154 (the acceptance-measure merit).** Three reasons for the first and two
against the second, and one reason to rank both below an operator-level
question.

**For re-measuring 153:**

1. *Its stated reason for being reserved is now false.* Section 153-154 says
   the two controls were kept default-off because "with both in place the two
   controls bottom out at the same place and the same number, `||R|| = 1.00` in
   BASE CELL `j = 1` ... so the obstruction had left the solver before either
   section was written, and section 152 is where it went." Section 152 has now
   landed, and the obstruction did not follow it: the base momentum row fell by
   three digits, the iteration moved to the carrier row at the front, and the
   floor is now `1.065` in a different place. Whether a step control helps is
   therefore an open question again, and the sections themselves say they are
   "kept default-off pending a re-measurement on top of the characteristic
   boundary". This is that re-measurement.
2. *The failure mode now visible is the one 153 was measured to fix.* This
   campaign's line search refuses **78 to 102 trial states per attempt for the
   element budget / an uncertified composition** (worst reaction residual
   `7.403E+06` in attempt 1, `1.548` in attempt 3, worst cell count 38-47), and
   collapses to `lam = 9.54E-07` repeatedly before each damped Gauss-Newton
   escape. Section 153's one robust measured advantage, in both its own runs
   and in section 154's head-to-head, is **zero element-budget refusals against
   the line search's 18 and 19**, because it enforces admissibility inside the
   step instead of discovering it at the trial point. That is a direct match to
   what is failing here.
3. *The patch still applies.* `p51d2_trust_region.patch` applies to the current
   tree with 13/13 hunks and no fuzz (offsets 65-75 lines). It touches only
   `steady_newton.f90` and only the four-unknown branch, gated by
   `EXHALE_TRUST_REGION=1`, so the three-unknown path and every atomic golden
   are untouched by construction.

**Against 154:**

1. Its own head-to-head at identical merit found the **line search better**:
   `||R||` `4.551E-01` vs `1.004` (state A) and `3.295E-01` vs `1.004`
   (state B), flux spread 3x better, carrier row 2x better. The trust region
   won only the element-budget column.
2. Changing the line search's merit to the acceptance measure made the
   `mol_carrier` case **worse** -- `||R||` `3.430E-01 -> 1.279`, spread
   `1.528E-02 -> 2.729E-01`, carrier `1.53E-03 -> 6.65E-03` -- because the
   acceptance measure weights the cell nothing can balance and the minimizer
   walks toward it. That diagnosis was made when the unbalanceable cell was
   `j = 1`; the cell has moved, so the effect might have changed sign, but
   nothing here says it has, and `p51d2b` is a superset that would drag the
   merit change in with the trust region. Apply `p51d2_trust_region.patch`
   alone.

**Why both rank below an operator-level question.** The Picard control arm
stagnates on the same state with no carrier unknown at all, in the energy row
of the same molecular layer, and leaves the carrier row 7.5x its target. A step
control cannot repair a state that is not near a root. Section 10.2 of
`p23_transport_on_state.md` measured exactly that: the carrier marching step and
the carrier residual do not agree in the `dt -> 0` limit -- `|G+R|` is identical
to four digits across `dt` spanning three decades, worst at the base cells
(1.9e-2 of the row's own terms) and 3.4% at the front -- and named the one
unisolated candidate, the **base inflow Dirichlet**. Every coupled hand-off in
this campaign starts from a marched state, so that mismatch is upstream of the
solver. **Isolate the base inflow Dirichlet first; re-measure 153 second.**

> **Corrected 2026-09-04, `Update_EXHALE.md` section 158.** The recommendation
> above was right about the order -- the operator before the step control --
> and wrong about what the operator defect was. The two paths do **not**
> construct the base inflow condition differently: they call the same assembly
> sequence in the same order, and `carrier_base_state` builds that condition in
> one place for both. Evaluated in the variable the carrier residual is written
> in, the marching step converges on the residual cleanly at first order
> (`9.5e-5`, `7.2e-6`, `8.2e-7`, `8.3e-8` over four decades of `dt`). The flat
> `|G+R|` of section 10.2 was the diagnostic's own choice of variable: it formed
> `G` in the carrier DENSITY, so the whole `dt`-independent part was the HYDRO
> stage's term arriving as a constant. The real defect it led to is the
> ADVECTED VARIABLE -- the carrier was transported as a particle mixing ratio
> against `n_tot`, which ionization and dissociation do not conserve, leaving a
> spurious `n_c v d ln(mbar)/dr` exactly across the two fronts -- and section
> 158 moves the row onto the fraction per unit mass that the mass row does
> conserve.

**A cheap third item.** The gate window was still falling monotonically at
200,000 steps here (`9.7E-02`), where the earlier campaign turned around near
that step count. Marching to 400,000 and re-offering the hand-off costs about
1.6 hours of wall clock and would say whether the coupled floor of `1.06` is a
property of the solve or of how far the march had got. It was not tested.

## 8. What was not established

* Whether marching past 200,000 steps changes any of this. Not tested.
* Whether the coupled floor `||R|| ~ 1.06` is the same number on a different
  planet or a different base level. One case only.
* The third Picard attempt of the first launch was lost when that process ended
  without writing a completion line; the arm was re-run from the same restart
  files and the numbers in section 4 are from the re-run, whose attempts 1 and
  2 reproduce the first launch to all printed digits (`5.928E-01` /
  `3.831E-01`, spread `9.222E-03` / `2.609E-02`). The cause of the first
  launch's early exit was not diagnosed.
* Sections 153 and 154 were NOT applied or run in this campaign. The
  recommendation above rests on their recorded measurements plus a dry-run
  check that the patch still applies, not on a new measurement of them.
* No transmission-spectrum or post-processed quantity was computed; the P23
  table is read directly from `Hydro_ioniz.txt` and `Ion_species.txt`.
