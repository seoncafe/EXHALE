# L4h -- the pseudo-time ramp, both ends

Plan item L4h of `docs/PLAN_20260913_lhs_stationary.md`, approved 2026-09-14
15:05. The item is the sentence "the ramp's shrink on damped steps can pin the
pseudo-time at its floor", and the opposite end recorded beside it: a
pseudo-time that climbs to 1e+14 while the residual stands at 1.9.

Code: `src/modules/time_step/steady_newton.f90` only, in
`solve_steady_jfnk`. Two option keys:

| key | default | what it does |
|---|---|---|
| `EXHALE_PTC_RAMP_GUARD` | **on** (`=0` turns it off) | statements (1), (2) and (4) of section 4: the cut of an accepted step clipped at ten, its floor a fixed span below the pseudo-time the solve has earned at a whole step, and a cap on a run of damped steps. With it off the binary is the arithmetic of L4g/L7d, bitwise (section 5.1). |
| `EXHALE_PTC_RAMP_GROWTH_GATE` | **off** (`=1` turns it on) | statement (3): the doubling of a full step held where the stagnation window bought less than `merit_fall_min`. Off because it was MEASURED to hold the near states the doubling exists for (section 5.5). |

`docs/input_schema.md` documents the `input.inp` keys and carries no table of
solver environment variables; these two are documented here and at their
declarations in `steady_newton.f90`.

**THE DELIVERED RULES ARE IN THE TREE BINARY `EXHALE.x`, md5
`c2e9c9990b9f14f1be8cd77abca68945`**, built 2026-09-16 from the finished
source of the 2026-09-13..16 series; its manifest
`LHS1140b/models/BINARY_MANIFEST_c2e9c9990b9f.txt` records
`src/modules/time_step/steady_newton.f90` at md5
`b907b67d63f1725f2cf5271a30d6310b`, which is the file this memo describes.
**Everything named `EXHALE_L4h*.x` below is a MEASUREMENT build of this item
and will be deleted at the end of the series**; the numbers in section 5 were
taken with them and are quoted against the build each was run on.

| measurement build | md5 | what it was |
|---|---|---|
| `EXHALE_L4h_ctl.x` | `9a500b4089fdefdb59b49bd0c03540ad` | the tree before this item's edit, used once, in section 2, to show that the current source reproduces the recorded L14 run |
| `EXHALE_L4h.x`, final | `aaf762175e5481547de12285e814683a` | the delivered rules; its own control is itself with `EXHALE_PTC_RAMP_GUARD=0` |
| `53bcbe13648bfc19cdf8337245d23c84` | | the per-step residual ratio refuted in section 5.5 |
| `f4a50f708b37f41907a6b77c0c26bc81`, `8ab3b72d4c03bc4223242e2a69b7f119` | | the earned-pseudo-time floor refuted in section 5.7.1 |
| `06fde0d649c1947042c2dc303497eafb`, `9a944cdf573541a47a50ec9599703c07` | | the CFL-discriminator rules, on which sections 5.2, 5.3 and 5.7.3 were measured |

A control is always the SAME binary with `EXHALE_PTC_RAMP_GUARD=0`, never a
binary built from an earlier state of the tree: other items land in `src/`
while this one runs. The tree's own `build/` was not touched; this item built
into `build_L4h/` and `build_L4h_ctl/`.

---

## 1. The verdict first

**The pin is a defect of the step-size rule, not of the Newton direction.**
At the iterations the residual stops moving, the Krylov cycle reaches the
linear tolerance it is asked for -- in ONE product -- and the step it returns
solves the system it was given. What that system has become is the problem.
MEASURED at the explicit-stable floor `dtau` = 7.93e-05: the shift `1/dtau` =
1.3e+04 stands above every entry of the assembled band, one matrix-vector
product reaches the requested tolerance, and sixty iterations leave the merit
and the state at the printed digits. Read as an explicit march at the CFL
interval, which is what it looks like, such a march needs
`O(t_flow/dt_CFL) ~ 1e+05` steps on this column against a cap of 200.

**That reading is a description of a measurement, not a theorem**, and
section 4 states where the code leans on it and where it does not: the system
is `(M/dtau + J) dY = -F`, an explicit-Euler limit `dY -> -dtau F` needs
`dtau M^-1 J` small in an operator sense, and a hydrodynamic CFL estimate does
not establish that for a `J` that also carries the source, radiation and
chemistry couplings. What is certain is the measurement: at that floor the
iteration buys nothing, and the residual is pinned by the step-length rule and
not by the direction.

**How the pseudo-time got there is the asymmetry of the ramp.** The rule of
L7d multiplies `dtau` by the line-search fraction `lam` on a damped step and
by two on a full one. `lam` is a halving sequence with no lower bound (its
20th backtrack is 9.5e-07), so ONE short accepted step can undo twenty
doublings, and there is no floor on the accepted-step cut other than the
explicit-stable interval itself, which is twelve decades below where these
solves start. The fall is not the rejection cut: every step in the measured
collapse was ACCEPTED and lowered the merit.

**At the other end the same rule reads only the step length.** It doubles on
every `lam >= 0.5` step whatever the residual did, so on a state whose
residual is flat it removes its own shift in about forty iterations and hands
the solve to the unshifted Newton on a state it is not near. The obvious
repair is the rule of the pseudo-transient literature, a factor read from the
residual's own progress, and **that repair was implemented and MEASURED not
to work here**: the near states the doubling exists for move their residual
as little per iteration as the flat ones do, so every progress test tried
holds their pseudo-time down as well. The flat end is therefore left open,
with the gate delivered off and the measurement recorded.

**The collapse end is repaired** by bounding the cut of a step that
succeeded, by stopping that cut a fixed span below the pseudo-time the solve
has earned at a whole step, and by capping a run of damped steps. **The flat end is
not**: a gate on the growth was implemented and MEASURED not to separate the
two regimes the ramp has to serve, and it is delivered off. Section 4 states
the four rules and section 5 measures them; sections 3.3 and 5.5 are the flat
end.

---

## 2. What was reproduced, and with which binary

`LHS1140b/models/.L14/x002_HeH2.13` is the cleanest of the three recorded
cases: the 0.02-XUV rung of the He/H = 2.13 column of LHS 1140 b, seeded from
the CERTIFIED 0.03 state beside it, `EXHALE_PTC_DTAU0=1.0e8`,
`EXHALE_OUTER_PASSES=40`, 8 threads.

The control build of this item was run on a fresh copy of that case with the
same seed and the same environment
(`LHS1140b/models/.L4h/ctl_x002_HeH2.13`). It reproduces the L14 run **to
every printed digit** over the iterations it was allowed to reach -- `||R||`,
the merit, `lam`, `dtau`, the Krylov product count and the worst row, at
iterations 1 through 41, over a pseudo-time range 2.0e+08 to 1.05e+14 that
includes full steps, damped steps, the ceiling and the cut. The tables of
section 3 are therefore quoted from the L14 log of the same case, which is
the same arithmetic on the current source.

The two other recorded cases are `.L14/x002_HeH9.7` (energy row of cell 12,
`lam` 5.0e-01 at `dtau` 2.56e+10) and `.L14/p002_HeH9.7` (energy row of cell
26, `lam` 6.25e-02 walking down from `dtau` 3.91e+05); copies with their own
seeds are in `LHS1140b/models/.L4h/`.

The opposite end is `LHS1140b/models/.stopped/atomic_scalar_gj1132x0.01_
kzz1e9_HeH{9.7,2.13}_202609140938`. **Those two logs were written by a binary
that predates the L7d rule**: their `dtau` doubles on EVERY accepted step,
including `lam` = 1.6e-02, and falls by exactly four only on the rejected
ones (`lam` = 9.54e-07, the search's twentieth halving). That is visible in
the logs and is stated here because it changes what the pair is evidence
for: it is evidence about growth keyed on acceptance alone, which the
current rule has inherited in the milder form "growth keyed on step length
alone".

---

## 3. The diagnosis

### 3.1 The ratchet, in the first pass of `x002_HeH2.13`

Accepted steps only; no iteration in this window was rejected.

| it | `||R||` | merit | `lam` | `dtau` | worst row |
|---|---|---|---|---|---|
| 131 | 1.279e-02 | 1.98e-07 | 1.00e+00 | 7.81e+19 | energy 79 |
| 134 | 1.029e-02 | 1.68e-07 | 1.00e+00 | 7.81e+19 | mass 32 |
| 135 | 1.174e-02 | 1.85e-07 | **3.91e-03** | 3.05e+17 | energy 46 |
| 138 | 1.815e-02 | 1.92e-07 | 3.12e-02 | 3.81e+16 | energy 5 |
| 141 | 1.284e-02 | 1.93e-07 | 3.12e-02 | 5.96e+14 | energy 9 |
| 145 | 9.139e-03 | 1.69e-07 | 6.25e-02 | 2.33e+12 | energy 44 |
| 149 | 1.253e-02 | 1.36e-07 | 3.12e-02 | 5.82e+11 | energy 10 |
| 156 | 1.124e-02 | 1.20e-07 | 6.25e-02 | 7.28e+10 | energy 4 |
| 159 | 1.122e-02 | 1.20e-07 | 2.50e-01 | 2.84e+08 | energy 4 |
| 167 | 6.817e-03 | 1.04e-07 | 1.56e-02 | 7.11e+07 | energy 36 |
| 171 | 5.950e-03 | 1.03e-07 | 3.12e-02 | 1.73e+04 | mass 1 |
| 172 | 6.965e-03 | 1.04e-07 | **3.91e-03** | **6.78e+01** | energy 64 |

Eighteen decades of pseudo-time in 41 iterations, against a merit that falls
by 47 percent and a `||R||` that falls by a factor 1.8. The doubling can
return one factor 2 per full step; a single `lam` = 3.9e-03 step takes eight.

### 3.2 The pin, in outer pass 39 of the same run

The pass enters at `||R||` = 5.082e-04 and hands back the state it was given.

| it | `||R||` | merit | `lam` | `dtau` | Krylov products | linear tolerance |
|---|---|---|---|---|---|---|
| 47 | 2.158e-04 | 1.25e-08 | 1.00e+00 | 1.91e+02 | 3 | reached |
| 51 | 2.886e-04 | 1.50e-08 | 3.12e-02 | 2.98e+00 | 3 | reached |
| 53 | 1.763e-04 | 1.44e-08 | 7.81e-03 | 4.66e-02 | 3 | reached |
| 57 | 2.019e-04 | 1.42e-08 | 1.00e+00 | 7.45e-01 | 2 | reached |
| 58 | 2.025e-04 | 1.42e-08 | 2.44e-04 | 1.82e-04 | 3 | reached |
| 59 | 2.028e-04 | 1.42e-08 | 3.91e-03 | **7.93e-05** | **1** | reached |
| 60 | 2.024e-04 | 1.42e-08 | 1.95e-03 | 7.93e-05 | 1 | reached |
| 61 | 2.024e-04 | 1.42e-08 | 6.10e-05 | 7.93e-05 | 1 | reached |
| 62 | 2.015e-04 | 1.42e-08 | 3.91e-03 | 7.93e-05 | 1 | reached |
| 63 | 1.480e-03 | 1.23e-08 | (damped Gauss-Newton escape) | 1.59e-04 | 1 | reached |
| 64-68 | 1.46e-03 | 1.23e-08 | 7.8e-03 ... 2.4e-04 | 7.93e-05 | 1 | reached |
| 69 | -- | -- | no descent along either direction; `||grad merit||` 1.657e-06 | | | |

7.93e-05 is `dtau_floor`, the explicit-stable interval `CFL min_j dr_j /
(|v_j| + c_j)` of the state the solve was entered at. The pass is 69
iterations long, 60 of them spent walking the pseudo-time down and 10 sitting
on the floor, and the merit falls 1.62e-08 -> 1.23e-08.

**The direction is not what fails.** Every one of those iterations reports
`the requested linear tolerance was reached`, and at the floor it is reached
in ONE matrix-vector product, which is what a preconditioned system equal to
the identity looks like. The judged rows are unmoved throughout (hydrodynamic
distance 5.67e+02, 495 of 500 cells outside the mass tolerance), and the
`lam` the search ends on is the outcome of backtracking against a merit
change that the step, of size `dtau |F|` with `dtau` = 7.9e-05, cannot make:
the ramp then multiplies the pseudo-time by that number and the loop closes
on itself.

The same shape, at the same place, ends the other two rungs:
`x002_HeH9.7` at the energy row of cell 12 with `lam` 5.00e-01 and `dtau`
2.56e+10, and `p002_HeH9.7` at the energy row of cell 26 with `lam` 6.25e-02
and `dtau` 3.91e+05 falling.

### 3.3 The other end

On the 0.01-XUV pair the pseudo-time runs from 1.00e+00 to its `1e14*dtau0`
ceiling in about fifty iterations while `||R||` stands at 1.9 and the merit
falls a few percent an iteration (6.07e-05 -> 2.77e-05 over twenty-five
iterations, He/H = 2.13). Once there the shift is 1e-14 and the iteration is
the unshifted Newton on a state it is not near; it neither converges nor
falls back: the cut is taken only on an iteration that admitted no step at
all, and the search keeps finding short steps that pass the non-monotone
test, so the growth outruns the cut.

The residual ratio these iterations carry is ~1. **A switched-evolution-
relaxation ramp would not have moved `dtau` at all on them**, which is the
whole content of statement (3) below.

---

## 4. The rules

All of them are in `solve_steady_jfnk`, guarded by `ptc_ramp_guard_on`
(`EXHALE_PTC_RAMP_GUARD=0` restores the entry text), and they act only where
the L7d ramp acts, i.e. with `pseudo_time_doubles_on_an_accepted_step` on,
which is the default.

**THE DISCRIMINATOR THE RULES REST ON, and it took three attempts to find.**
A damped step means two different things, and the two have to be told apart
because the right answer is opposite in each:

- **`dtau >> dt_CFL`.** The shift `I/dtau` is one term among the Jacobian's,
  the step is a Newton step the model shaped, and a damped step that bought
  nothing says the pseudo-time is too large. **Cutting is right here** -- it
  is what pseudo-transient continuation IS. The LHS 1140 b 0.02 rung's useful
  march from 1e+08 down to ~1e+02, six decades of `dt_CFL`, lives entirely in
  this regime.
- **`dtau <~ dt_CFL`.** **Cutting further was MEASURED to be inert** -- on
  the element reload, sixty iterations below that scale leave the merit and
  the state at the printed digits while `dtau` falls twelve decades, and at
  `dtau` = 7.93e-05 the Krylov cycle reaches its tolerance in ONE product with
  the shift 1.3e+04 above every entry of the assembled band. **This half of
  the discriminator is an empirical control, not a proved limit** (section
  4.3), and `ptc_explicit_march_factor` is a knob for that reason.

Two other discriminators were implemented and MEASURED to fail, and both
failures are in section 5.7: **how far the pseudo-time has travelled** (a
floor a fixed span below the largest pseudo-time a whole step was taken at)
serves the element reload and costs the LHS rung a factor 3.2 on its first
pass; **whether this step progressed** (hold the pseudo-time where a damped
step bought less than a fraction of the merit) is worse still, because
non-progressing damped steps happen in BOTH regimes -- it freezes `dtau` at
4.10e+03 through nine of them on the element reload and ends that pass worse
than the control.

**(1) The cut of an ACCEPTED step is bounded at ten.**
`dtau_ramp = max(lam, ptc_shrink_clip)`, `ptc_shrink_clip = 0.1`. The reason
the pseudo-time falls at all on a damped step is L7d's and stands: the linear
model was good over the fraction `lam` and no further. What does not stand is
an unbounded fall on a step that SUCCEEDED -- `lam` is a halving sequence
whose twentieth backtrack is 9.5e-07, so one short accepted step can undo
twenty doublings. The factor of ten is the one the implicit CFD practice
reports, stated of the timestep devices generally and not of the SER ratio
alone: "All such devices are "clipped" into a range about the current timestep
in practice. Typically, the timestep is not allowed to more than double in a
favorably converging situation, or to be reduced by more than a factor of ten
in an unfavorable one, unless feasibility is at stake, in which case the
timestep may be drastically cut" (Gropp, Keyes, McInnes and Tidriri,
NASA CR-1998-208435, section 2.1, verbatim).

**That sentence is the only published counterpart this rule has, and it is a
report of practice, not a hypothesis of any theorem** (section 4.2).

**(2) And nothing else bounds it.** The cut proceeds normally, floored at the
explicit-stable interval `dtau_floor` exactly as it always was. The march the
continuation is for is not forbidden at any depth.

**(3) A gate on the growth, WHICH IS OFF** (`ptc_growth_gate_on`,
`EXHALE_PTC_RAMP_GROWTH_GATE=1`). With it on, a full step doubles the
pseudo-time only if the stagnation window bought at least `merit_fall_min`.
It is off because no residual-progress test separates the two regimes of the
ramp (section 5.5): the near states the doubling exists for move their merit
as little per iteration as the flat ones do.

**(4) A run of iterations below that measured scale is capped.**
`n_damped_run` counts iterations that took no full step and no progressing
damped step, AND only while `dtau <= ptc_explicit_march_factor*dtau_floor`
(`ptc_explicit_march_factor` = 10): above that a cut still changes the step
and there is nothing for the count to be about. At `ptc_damped_run_max` = 8
the solve returns to its best iterate with `dtau = max(0.1*dtau_best,
dtau_floor)` -- a tenth of the pseudo-time that iterate was REACHED at, so
the step after the return is not itself an explicit one -- and after
`ptc_damped_restarts_max` = 2 such returns with no full step between them it
stops and hands back its best iterate with `info = 2` and a line naming the
reason.

A stronger form of (4), counting EVERY iteration below that threshold
so that a full explicit step does not reset the count, was implemented
(`march_counts_every_iteration`, `EXHALE_PTC_MARCH_COUNT_ALL=1`) and MEASURED
to be worse: it fires on the element reload where the delivered form does not,
and ends that pass at `||R||` 8.730e-01 against the delivered 6.121e-02
(section 5.7.3). It is off.

On `info`: no new value was introduced. Every caller of `solve_steady_jfnk`
branches on `info == 0` alone (`src/EXHALE_main.f90`), so a new code would be
read exactly as 2 is; the reason is named in the log line instead.

The trust-region route is untouched in substance: it sets `lam = 1` on every
iteration, so (1) never fires there, and the damped-step count advances only
on its rejected iterations and only below that threshold.

**The default is on. The fixed point of the iteration is `F(Y) = 0` at every
`dtau`** -- the pseudo-time enters the step and not the residual -- so what a
ramp rule can change is which states are reached and how fast, and that is
what section 5 measures.

### 4.3 The one place the rules lean on a reading, and how far

Statement (4) acts only below `ptc_explicit_march_factor * dt_CFL`, and that
threshold is the single place where the design leans on an interpretation
rather than on a number.

**What was measured** (element reload, one thread, section 5.2): at
`dtau` = 7.93e-05 the shift `1/dtau` = 1.3e+04 stands above every entry of the
assembled band; the Krylov cycle reaches its requested tolerance in ONE
matrix-vector product; and sixty consecutive iterations leave the merit and
the state unchanged at the printed digits while the pseudo-time falls twelve
decades. Cutting there buys nothing.

**What was NOT established.** The shifted system is `(M/dtau + J) dY = -F`.
An explicit-Euler limit `dY -> -dtau F` needs `dtau M^-1 J` to be small in an
OPERATOR sense. A hydrodynamic CFL estimate -- `min_j dr_j/(|v_j| + c_j)`,
which is what `dt_CFL` is here -- does not establish that for a `J` that also
carries the source, radiation and chemistry couplings of this system, and a
finite multiple of a stability limit is not an asymptotically small parameter
in any case. "The iteration is the explicit march below `K dt_CFL`" is
therefore a NAME for what was measured on one fixture, not a derivation.

**What follows for the code.** The threshold is an empirical control and is
exposed as `EXHALE_PTC_EXPLICIT_MARCH_FACTOR` so that a state where the
measurement comes out otherwise can be met without a rebuild. `K` = 3, 10 and
100 were run on the element reload and give the same pass to the digit,
because statement (4) never fires on that fixture at all -- so the sensitivity
study is vacuous there and the threshold is, so far, untested by a case that
exercises it.

### 4.1 What was checked at the source, and what was not

| reference | checked | how |
|---|---|---|
| Kelley and Keyes 1998, SIAM J. Numer. Anal. 35, 508-523, "Convergence Analysis of Pseudo-Transient Continuation", doi 10.1137/S0036142996304796 | **OBTAINED and read** (sections 1.1, 1.2, 1.3) | `references/Kelley_1998SIAMJNA_35_508.pdf`, JSTOR scan; formula (1.7) of Assumption 1.1 is not legible in it and is not quoted |
| Coffey, Kelley and Keyes 2003, SIAM J. Sci. Comput. 25, 553, "Pseudotransient Continuation and Differential-Algebraic Equations", doi 10.1137/S106482750241044X | **OBTAINED and read** (sections 1, 2, 3) | `references/Coffey_2003SIAMJSC_25_553.pdf`, publisher PDF with a text layer |
| Mulder and van Leer 1985, J. Comput. Phys. 59, 232-246, "Experiments with Implicit Upwind Methods for the Euler Equations" | **OBTAINED and read** (section 2, equations (2), (3), (4) and the limit-cycle discussion) | `references/Mulder_1985JCP_59_232.pdf`, scanned; the scan renders Greek letters inconsistently (alpha as "c1"/"CI"/"a", epsilon as "E"/"c") and they are restored in the quotations below |
| Gropp, Keyes, McInnes and Tidriri, "Globalized Newton-Krylov-Schwarz Algorithms and Software for Parallel Implicit CFD", NASA CR-1998-208435 / ICASE Report 98-24, 40 pp, section 2.1 | **OBTAINED and read** | NTRS id `19980233244`; `references/Gropp_1998_NASA_CR_208435.pdf` |

All four were read. What follows is what they do and do not say about the
rules of section 4.

### 4.2 What the papers support, and what is ours alone

**The name.** Kelley and Keyes section 1.2 writes "the "switched evolution
relaxation" (SER) method, so named in [21]"; Gropp et al. section 2.1 writes
"successive evolution-relaxation (SER) [56]". Both references are Mulder and
van Leer 1985: the two names are the same rule from the same source, and an
earlier revision of this memo wrongly reported the second as a correction of
the first. **Neither name appears in Mulder and van Leer's own text** as
scanned; the paper states the rule and the later literature named it.

**What the primary source actually prescribes is a TARGET ON THE CHANGE OF
THE STATE, not a residual ratio.** Mulder and van Leer section 2 monitors
convergence with a scaled max-norm residual (their equation (3),
`RES^n = max_{k,i} |g_ki| / (|w_ki| + h_ki)`, `h_ki` "some positive constant
that prevents division by zero") and then, verbatim: "The time-step is derived
from this quantity according to" their equation (4),

    dt^n = epsilon / RES^n,

whose purpose they state as "Equation (4) guarantees that in the explicit case
(alpha = 0) the relative change of the state quantities per time-step will
nowhere exceed epsilon". The successive-ratio form of SER is a reading of (4)
between two steps, not what (4) says; in the taxonomy of Gropp et al. section
2.1 this is their third device, the one that "sets target maximum magnitudes
for change in each component of the state vector and adjusts the timestep so
as to bring the last measured change to the target". No cap on `dt` is stated
in that passage.

**AND THE ORIGINAL DOES NOT COUPLE THE TIMESTEP TO A RELAXATION FACTOR.** The
paper does contain an under-relaxation, and it is a different device in a
different place. Their scheme is `L^n Delta_t W = [I/dt^n + alpha M^n]
Delta_t W = G^n` (equation (2)), where "For alpha = 1 and M = dG/dW (the
Jacobian of G with respect to W), we have the backward Euler or implicit Euler
scheme"; against a two-vector limit cycle in one test problem they write "The
most obvious remedy is to make alpha > 1, implying under-relaxation. It turned
out that for alpha >= 1.2 the scheme converged: the convergence slowed down
with increasing alpha." **The under-relaxation is applied to the OPERATOR, by
over-weighting the implicit term, and equation (4) for the timestep is
untouched by it.** The timestep is never reduced on account of the iteration
being damped.

So the search through the primary source closes the question it was opened
for: **our `lam`-keyed cut has no published counterpart in any of the three
papers** -- not in Kelley and Keyes' or Coffey et al.'s `phi`, which only
caps; not in Gropp et al., whose factor of ten is a report of practice; and
not in Mulder and van Leer, whose under-relaxation acts on `L` and leaves the
timestep rule alone. Statements (1), (2) and (4) rest on section 5 and on
nothing in the literature.

**The form is cumulative in both published statements.** Kelley and Keyes
(1.5) is `delta_n = delta_0 ||F(x_0)||/||F(x_n)||`, described verbatim as
"In its simplest, unprotected form, SER increases the timestep in inverse
proportion to the residual reduction"; Coffey, Kelley and Keyes (1.3) writes
both forms and shows them equal,
`delta_n = delta_{n-1} ||F(u_{n-1})||/||F(u_n)|| = delta_0 ||F(u_0)||/||F(u_n)||`.
The per-step ratio this item implemented and refuted (section 5.7.2) is the
derived relation, not the primary definition; **the cumulative form was never
implemented here** and the sections below make no claim about it.

**The clip on the GROWTH is supported, explicitly.** The phi of Coffey,
Kelley and Keyes (1.5) is `phi(xi) = xi` for `xi <= xi_t` and `delta_max`
above it -- a cap and nothing else -- and their section 3 reports of their own
run, verbatim: "In this implementation, the maximum increase of delta_n from
one time step to the next is limited to a factor of 2". Kelley and Keyes
Assumption 1.1 admits the same family ("Other choices could either limit the
growth of delta or allow delta to become infinite after finitely many steps.
Our formal assumption on phi accounts for all of these possibilities").

**The CUT is not.** Neither paper's update ever reduces the pseudo-time below
what SER gives: phi caps from above. The factor of ten this item's statement
(1) uses has its counterpart in the Gropp report's sentence about practice and
in no theorem.

**And the situation the rules are about is absent from both papers.** Their
algorithm takes the FULL step -- Kelley and Keyes Algorithm 1.1 step 2(b) is
"Set x = x + s" -- with no line search anywhere in it, and Gropp section 2.1
says why: the continuation "does not require reduction in ... at each step, as
do typical linesearch or trust region globalization strategies ...; it can
climb hills". **The quantity `lam` has no counterpart in either paper, and
neither does a cut keyed on it.** Statements (1), (2) and (4) are therefore
OURS: the citations name where the continuation and the growth clip come from
and are not evidence for them. What justifies them is section 5 and nothing
else.

**The explicit-stable interval as a floor has no counterpart either.** In
Kelley and Keyes section 1.1 the CFL enters as the SCALING MATRIX `V`, which
"is typically diagonal and approximately equilibrates the local CFL number
(based on local cell diameter and local wave speed) throughout the domain" --
not as a bound on `delta`.

One thing the papers do bear on that this item met independently: Coffey,
Kelley and Keyes section 3 reports that "the success of the method is
relatively sensitive to the initial time step", which is the `dtau0`
sensitivity L4g and L7d record and which this item's reproduction cases run
at `EXHALE_PTC_DTAU0=1.0e8` for.

---

## 5. Controls and verification

### 5.0 Every control here is SINGLE-THREADED, and that is not a detail

**This solve is not reproducible at eight threads.** Two runs of the same
binary with the same input and `OMP_NUM_THREADS=8` agree for about twenty
iterations and then part: measured on the element reload, the control build
and the measured binary WITH THE GUARD OFF -- the same arithmetic -- printed
identical lines through iteration 21 and different ones from iteration 22,
and ended pass 1 at mass rows 3.08e-02 and 1.20e-02. The reduction order of
the threaded sweeps moves the last bits, and a Newton iteration on these
states amplifies them.

At `OMP_NUM_THREADS=1` the same case is bitwise reproducible: two runs of the
control build print identical lines for the 28 iterations compared. Every
control/measured comparison in 5.1, 5.2 and 5.5 is therefore taken at one
thread. **An eight-thread comparison of this fixture measures the thread
schedule, not the change**, and a first attempt at 5.2 that was taken that
way was discarded.

**NOTE ADDED 2026-09-17 (item L15).** The divergence recorded in this section
is NOT reproducible on the present tree. Traced quantity by quantity, the whole
105-iteration trajectory of this fixture, its 152 Krylov products and 923
recorded quantities, is bitwise identical at 1, 8 and 16 threads (one trace md5
`bd12e7c8ad917beda538d789080e10e3`), and so are the two molecular
configurations and the closing run on the catalog binary. The same recipe that
ended pass 1 at mass rows 3.08e-02 and 1.20e-02 here now ends it at 6.11e-10
and 4.57e-10. The 2026-09-15 observation is left as taken and is attributed to
nothing; `docs/lhs1140b_stationary_L15_20260916.md` stages 2, 5 and 6. The rule
this section states, that a control and a measured run are compared at one
thread, is kept: non-reproduction on these fixtures is not a proof that no race
exists.

(This is the regression harness's own rule -- `run_check.sh` runs every case
at `OMP_NUM_THREADS=1` "so results are deterministic" -- applied to the
fixtures that are not in the harness.)

### 5.1 With the option off

`EXHALE_PTC_RAMP_GUARD=0`, the delivered binary against `EXHALE_L4h_ctl.x`
(the tree as it stood before this item's edit), on the element reload of 5.2,
same seed, same environment, one thread: **every printed JFNK iteration
identical to every digit** -- `||R||`, the merit, `lam`, `dtau`, the Krylov
product count, the worst cell and the worst row -- over the 28 iterations
compared, which include the whole of the collapse this item is about
(iterations 12 to 23, where the pseudo-time goes from 2.05e+03 to the
explicit-stable interval). The guard adds no arithmetic to the path when it
is off.

### 5.2 The element reload, `backup/regression/atomic_elem_newton`

HD 209458 b, atomic with solar trace metals and binary element diffusion,
reloaded from `IC/` with `Load IC? True`, `Coupled carrier solve: False`,
`Restart intent: stationary`; `EXHALE_OUTER_PASSES=1 EXHALE_JFNK_MAXIT=80
EXHALE_DIFF_OMEGA=0.5`, **one thread** (section 5.0), on scratch copies. The
control is the SAME BINARY with `EXHALE_PTC_RAMP_GUARD=0`: a binary built
from an earlier state of the tree is not a control here, because other items
land in `src/` while this one runs.

| end of pass 1, 80 iterations | guard off | guard on |
|---|---|---|
| `\|\|R\|\|` | 8.336e-01 | **6.121e-02** |
| mass row | 3.08e-02 | **1.65e-02** |
| momentum row | 1.90e-05 | **1.44e-05** |
| energy row | 8.34e-01 | **6.12e-02** |
| worst gated species row | 2.83e-03 at cell 301 (elemental transport O) | the same, at the same cell |

**Thirteen times better on `||R||` and on the energy row**, the element
relaxation the same row at the same cell.

**And the whole of it is statement (1).** Statement (4) fired ZERO times, in
all of `K` = 3, 10 and 100 -- which is why those three runs are bit-identical
to each other, and why the `K` sensitivity is vacuous on this fixture. The
pseudo-time still reaches the explicit-stable interval 9.68e-05 in the
guarded run (iteration 24) as it does in the control (iteration 23); the
x 0.1 clip slows the descent by one iteration, it does not stop it. Once at
the floor the steps alternate `lam` = 1 and damped, so every `lam` = 1 resets
the count and it never reaches eight.

What the clip buys is the PATH: from iteration 28 the guarded run's `||R||`
falls -- 1.485, 1.236, 1.137 -- where the control oscillates between 1.5 and
1.9 for the rest of the pass. That is a weaker kind of evidence than a
mechanism, and it is stated as what it is.

### 5.3 The three 0.02 rungs: the ramp is repaired and the case is still not solved

The item's own reproduction, with the delivered binary at
`EXHALE_PTC_DTAU0=1.0e8`, 40 outer passes, 8 threads, each rung seeded from
the certified 0.03 state beside it.

**The collapse is gone.** On `x002_HeH2.13` the pseudo-time makes the march
the continuation is for -- 1e+08 down to ~5e+01 over the first pass -- and
stays there taking FULL steps: at iteration 175, `lam` = 1.00, `dtau` 1.05e+02
and doubling, the Krylov cycle reaching its tolerance in 3 products. Compare
the control at its own iteration 172: `lam` 3.91e-03, `dtau` 6.78e+01 having
fallen eighteen decades, and from pass 39 onward pinned at the explicit floor
7.93e-05 with the state handed back unchanged.

**And the first pass is the control's, to within the run-to-run spread:**

| pass 1 | control | delivered |
|---|---|---|
| mass row | 2.80e-04 | 2.41e-04 |
| momentum row | 5.30e-12 | 1.14e-11 |
| energy row | 6.45e-04 | 6.95e-04 |
| elemental row | 2.95e-02 at cell 374 | the same, at the same cell |

Eight per cent on the energy row, against an eight-thread run-to-run spread
that section 5.0 measures at a factor 2.6 on this class of fixture. **No
claim is made that the guard beats the control on this pass**; the claim is
that it no longer costs it, which the earlier designs did (section 5.7).

**THE 0.02 WALL IS NOT THE RAMP.** With the collapse removed the hydrodynamic
solve still stalls, decades above the anchored tolerance, and it stalls in the
same place on every rung:

| rung | passes | hydrodynamic rows over those passes | worst row at the stall |
|---|---|---|---|
| `x002_HeH2.13` | 1-7 | mass 2.41e-04, 2.73e-04, 2.11e-04, 1.87e-04, ..., 1.93e-04, 1.93e-04; energy 6.95e-04, 8.62e-04, 6.77e-04, 4.47e-04, ..., 5.91e-04, 5.99e-04 | mass of **cell 1** |
| `x002_HeH9.7` | 1-2 | mass 4.30e-05, 4.56e-05; energy 4.64e-05, 3.86e-05 | |
| `p002_HeH9.7` | 1-8 | mass 5.08e-02 and momentum 2.86e-09 IDENTICAL to every digit over all eight, energy 2.24e-01 -> 2.25e-01 | mass of **cell 31** |

The two scalar rungs are not flat pass by pass -- `x002_HeH2.13`'s fourth
pass falls on both rows, 6.77e-04 -> 4.47e-04 and 2.11e-04 -> 1.87e-04 -- but
over seven passes they go nowhere: mass 2.41e-04 -> 1.93e-04 and energy
6.95e-04 -> 5.99e-04, a factor 1.2 while the elemental row falls by sixteen,
2.95e-02 -> 1.81e-03. `p002_HeH9.7` does not move at all, two of its three
rows identical to every printed digit over eight passes.

The elemental row falls throughout (2.95e-02 -> 1.44e-02 on the 2.13 rung
over three passes, 2.94e-02 -> 1.78e-02 on the profile rung over eight), so
the composition is travelling and it is the WIND half that refuses -- the same
division L14 recorded, now with the ramp excluded as its cause.

The JFNK state at the stall says the solve is healthy in every respect except
the one that matters. On `x002_HeH2.13` at iteration 175:

```
||R|| 2.249E-04   ||Fs||2 1.69E-08   lam 1.00E+00   dtau 1.05E+02   gm 3
forcing term 1.00E-01: the cycle reached 1.314E-02 in 3 product(s) of 40;
                       the requested linear tolerance was reached
judged rows: distance 7.430E+02   hydrodynamic 7.430E+02
cells outside the tolerance of their row, of 500:
                       mass 491 (worst cell 1), momentum 0, energy 137 (2)
```

Full steps, a doubling pseudo-time, a linear solve that reaches its tolerance
in three products -- and a judged distance of **743**, held by the mass row of
cell 1 with 491 of 500 cells outside the mass tolerance and 137 outside the
energy tolerance while the momentum row is inside everywhere. The profile
rung is the same shape four decades worse: at its iteration 172, `lam` 1.00,
`dtau` 5.50e+02, 4 products, judged distance **1.501e+05**, worst row the mass
of cell 31, 500 of 500 mass cells and 494 energy cells outside.

The certification anchors these rows at mass 3.00e-12 as a FLOOR -- ten times
the cell's own rounding of the flux difference where that stands above it --
momentum 1.00e-08 and energy 1.00e-06.

**So L4h's verdict on this case is:** the ramp guard removes the collapse and
restores the march, and the three 0.02 rungs still do not certify, because
the hydrodynamic solve stalls at `||R|| ~ 2e-04` (and ~1e-01 on the profile
rung) with the mass row of cell 1 and the energy rows near the base binding.
That is a different defect. It was raised as plan item **L17** and closed by
the user the same day -- diagnosed, not solved: at 0.02 XUV the base is
hydrostatic to Mach 1e-8 with alternating-sign velocities, and the HLLC
contact selection puts every base cell on a kink of the residual (the line is
gone by 0.05). A smooth low-Mach flux would be a discretization change and
was not opened. Not a ramp matter, and not L4h's to fix. All three rungs were
stopped and their logs kept in `LHS1140b/models/.L4h/`.

**And the wall is a failure of the routes tried, not a proof that no route
exists.** What has been shown is that the ramp is not the obstruction and that
the base-adjacent rows refuse under the seeds, the pseudo-time starts and the
flux this campaign has used. Nothing here shows that the 0.02 and 0.01 columns
have no solution the code can reach, nor that a different route -- a smoother
low-Mach flux at the base, a different base boundary condition, a continuation
in a variable other than XUV, a tolerance anchored on something other than the
flux-difference rounding -- would fail. The three 0.01 cases are recorded as
NOT SOLVED BY THE ROUTES TRIED, which is a weaker and more accurate statement
than unsolvable.

### 5.4 The suites

Run from the item's own object directory `build_L4h/`, i.e. against the
delivered objects, with `EXHALE_SPECIES_EXE` pointed at the delivered binary.

| suite | PASS | FAIL |
|---|---|---|
| `krylov_and_dogleg` | 334 | 0 |
| `attempted_step` | 70 | 0 |
| `certification` | 84 | 0 |
| `steady_species_rows` | 199 | 0, the element and carrier reload rows included |

`steady_species_rows` and `steady_selfconsistent_residual` are the two suites
whose drivers link `steady_newton.o`; the other three were run because the
item names them. `steady_selfconsistent_residual` needs a run log
(`EXHALE_STEADY_RUNLOGS`) and was not run.

### 5.5 The growth gate, and why it is off

The fixture is the LHS 1140 b C/N/O column, `scratchpad/L4e/red_A` -- the
state L4g and L7d measured the ramp on -- reloaded with `Restart intent:
stationary`, `EXHALE_OUTER_PASSES=1 EXHALE_JFNK_MAXIT=40`, ONE thread, from
copies in the scratch directory (the item's record was not written to).

`dtau` and `||R||` over the first fourteen iterations:

| it | doubling (the delivered rule) | with the growth gate on |
|---|---|---|
| 1 | 2.00e+00, `\|\|R\|\|` 4.403e-01 | 2.00e+00, 4.403e-01 |
| 2 | 4.00e+00, 3.563e-01 | **2.00e+00**, 3.563e-01 |
| 3 | 8.00e+00, 3.294e-01 | 2.00e+00, 3.332e-01 |
| 4 | 1.60e+01, 3.105e-01 | 2.00e+00, 3.255e-01 |
| 5 | 3.20e+01, 2.799e-01 | 2.00e+00, 3.216e-01 |
| 6 | 6.40e+01, 2.775e-01 | 5.00e-01, 3.216e-01, `lam` 9.5e-07 |
| 7-13 | 1.28e+02 ... 8.19e+03, 2.9e-01 -> 3.1e-01 | 1.25e-01 down to 6.93e-05, `lam` 9.5e-07 at every one |
| 14 | 1.64e+04, **2.790e-02** | 1.00e-03, **1.959e+00** |

The gate holds the pseudo-time at 2 while the doubling carries it to 1.6e+04.
By iteration 6 the line search on the gated run finds no admissible step at
all (`lam` = 9.54e-07 is its twentieth halving), the cut takes `dtau` to the
explicit floor, and at iteration 14 the state is at `||R||` 1.959, **worse
than the 6.979e-01 it was given**, while the ungated run is at 2.790e-02 and
converging. This is the same shape L7d section 3.1 recorded for the
threshold at `lam` = 1.

The reason is in the numbers beside it: this state's merit hovers at 1.28e-01
to 1.42e-01 against an entry value of 1.44e-01 -- inside the tenth the gate
asks for -- while `||R||` falls from 6.979e-01 to 2.799e-01. **The residual
these states reduce is not the functional the gate reads, and the functional
the gate reads does not move on them.** A per-step ratio has the same
property and was measured to leave `dtau` at 1.13 on the same state.

So: no test of residual progress available in this solve separates the near
states the doubling exists for from the flat ones it damages. The gate is
delivered off, and the flat end stays open (section 6).

With the gate off, the delivered rules reproduce the doubling on this state
**iteration for iteration**: `dtau` 2, 4, 8, 16, 32, 64, 128, 256, 512,
1.02e+03, 2.05e+03, 4.10e+03, 8.19e+03, 1.64e+04 and the same `||R||` to
every digit over the fourteen iterations compared. Every accepted step of
this state is `lam` = 0.5 or 1, so statements (1) and (2) never fire on it --
which is the same reason L7d's threshold at 0.5 left it bitwise unchanged.

### 5.6 `wasp_full_newton`: the golden is already stale on the CONTROL

`wasp_full_newton` is the only default regression case that reaches
`solve_steady_jfnk`, and it is the one this item could move that carries a
golden. Run through a copy of `backup/regression` with `REGRESSION_EXE`, one
thread, the harness's own recipe.

**The control binary -- the tree as it stood BEFORE this item's edit, with no
guard in it at all -- already fails that golden on all four files:**

```
FAIL Hydro_ioniz.txt      max rel 6.532e-01 at row 1 col 7
FAIL Ion_species.txt      max rel 6.187e-01 at row 482 col 23
FAIL Hydro_ioniz_adv.txt  max rel 1.000e+00 at row 291 col 10 (0 against 1.95e-14)
FAIL Ion_species_adv.txt  max rel 6.187e-01 at row 482 col 23
==> REGRESSION FAIL
```

Sixty-five per cent on the continuity row of cell 1 and sixty-two on a
species column: **this case's golden was left behind by items that landed in
`src/` before this one**, and no verdict about this item can be read off
`golden/`. (L4g section 3.5 records the same for its own pair, "Both cases
fail against `golden/` with the control binary as well".) The
`Hydro_ioniz_adv.txt` row is the known `adv_mass_row` artefact: exact zero
against 1.95e-14, which the relative measure reports as a full difference.

The verdict that CAN be read is control against measured, and it is small.
Both binaries end the marching phase at the same step with the same numbers
(`final: count=4709 du= 9.9192E-03 dtu= 1.7604E-04`), and the outputs differ
by:

| file | max relative difference, control against measured |
|---|---|
| `Hydro_ioniz.txt` | **2.409e-09** (row 393, col 7) |
| `Ion_species.txt` | **5.880e-09** (row 395, col 23) |
| `Ion_species_adv.txt` | **5.880e-09** (the same entry) |
| `Hydro_ioniz_adv.txt` | 1.000e+00 -- the `adv_mass_row` artefact again, 9.96e-15 against an exact zero |

Six decades inside the harness's 1e-3, and the same order L7d measured for
its own change on this case (4.76e-09 and 7.98e-09). **This item does not
move `wasp_full_newton`.**

That comparison was taken with the superseded build `53bcbe13`, whose
statements (1) and (2) differ from the delivered ones; the run of the
delivered binary `88721ef4` against the same control is in
`scratchpad/regw_measF` and was still marching when the memo was closed.
`carrier_model_a_newton`, the four-unknown route L4g measured at 1.8e-10, was
not started.

---

### 5.7 The two designs and one refinement that were refuted

All on the same two fixtures, one thread, the same binary against itself with
the guard off.

**5.7.1 A floor a fixed span below the pseudo-time the solve had EARNED**
(the largest at which it took a whole step, span 1e-3; `dtau0` alone was
measured not to serve, being 1.0 on the stationary restart). It does exactly
what it says on the element reload: the guarded run's `dtau` is held at 2.05,
a thousandth of the 2.05e+03 of its last whole step, where the unguarded one
falls to 9.68e-05, and the pass ends at `||R||` 4.032e-02 against 8.336e-01.
**And it costs the LHS rung a factor 3.2 on its first pass** (energy row
2.08e-03 against the control's 6.45e-04), because that rung's march down to
~1e+02 is legitimate and the floor forbids it. Serves one fixture, breaks the
other.

**5.7.2 A hold conditioned on the step's own progress** -- cut where the
accepted damped step lowered the merit by at least a fraction `f`, hold where
it did not. Worse than either. MEASURED at `f` = 1e-2 on the element reload,
the pass ends at energy row **1.21e+00**, worse than the control's 8.34e-01,
and the trace says why: iterations 12 to 20 are nine damped steps in a row
that buy nothing, and the pseudo-time is FROZEN at 4.10e+03 through all of
them while `||R||` climbs 1.67 -> 1.95.

The premise fails because **non-progressing damped steps occur in BOTH
regimes**. At large `dtau` one is the statement that `dtau` is too large and
the cut is the right answer; at `dtau ~ dt_CFL` the same step means the
opposite. Progress-on-this-step no more separates them than
distance-travelled does.

The literature says the same thing about the premise, in a sentence that was
not read until after the measurement (NASA CR-1998-208435 section 2.1,
verbatim apart from the norm symbols the scan does not carry): "We emphasize
that pseudo-transient continuation does not require reduction in ... at each
step, as do typical linesearch or trust region globalization strategies ...;
it can climb hills." A rule that conditions the pseudo-time on per-step
residual reduction is asking the continuation for something it is explicitly
not required to give.

**5.7.3 Counting EVERY iteration below the measured threshold**
(`march_counts_every_iteration`), so that a `lam` = 1 there -- a full
EXPLICIT step, not a full Newton step -- does not reset the count. The
reasoning is sound and the mechanism works: on the element reload it fires
where the delivered form never does, three returns taken at useful
pseudo-times (5.12e+01 and 4.10e+01, not at the floor) and then the named
stop. **But the pass ends at `||R||` 8.730e-01 against the delivered
6.121e-02**, fourteen times worse and slightly worse than the control,
because counting at eight interrupts exactly the grinding at the floor that
earns the delivered form its 6.12e-02. The C/N/O near state is untouched by
it (identical to the control over 18 iterations), so the damage is confined
to the regime it acts in. Off.

With the flag off the delivered binary is identical to the build before the
variant existed over the full 99-iteration pass, so the variant costs nothing
where it is not asked for.

---

## 6. What is NOT settled here

- **The direct banded route `solve_steady_ptc` was not changed.** Its ramp is
  already the residual ratio with the cut clipped at ten
  (`dtau*max(lam,0.1)*(f2/f2_try)`) and its rejection cut floors at `dtau0`,
  so both statements of the guard are already true there. It is not the route
  any LHS 1140 b case takes.

- **The return of rule (4) does not release `dtau_earned`.** The rejection
  cut does, because an iteration that admits no step says the scale it
  records is no longer the state's; a run of eight damped steps is weaker
  evidence than that, and the pseudo-time is meant to come back once a whole
  step is taken again. So after a return the floor of the accepted-step cut
  is still a thousandth of the earned value and `dtau` climbs back to it on
  the first full step. That is deliberate and it is what the shift coming
  back looks like; if the return is ever measured to need a lasting cut, the
  anchor is the line to change.

- **One guard-off run of an intermediate build did not reproduce the
  entry-text build and was not explained.** `scratchpad/reg/elem_goff2`,
  build `8ab3b72d` with `EXHALE_PTC_RAMP_GUARD=0`, one thread, parted from
  `elem_ctl1a` at iteration 22, while the delivered binary under the same
  setting reproduces it for 28 (section 5.1) and the earlier `elem_goff1`
  reproduced it for the 11 it was compared over. The two directories were
  copied from `atomic_elem_newton_ctl` at different moments, and that
  directory was being written by an eight-thread run at the time, so the
  copies may not have carried the same auxiliary files. **It was not
  chased**; the statement of 5.1 rests on the delivered binary's own run.

- **One exit skips the damped-step count**: the iteration that finds the
  banded factorization singular cuts `dtau` by four and `cycle`s past the
  bottom of the loop. That is a factorization failure and a separate
  fall-back, and it is left as it is.

- **The best-iterate return does not re-evaluate the residual.** `Y`, `u` and
  `rnorm` are the best iterate's after the return; `F` is still the previous
  iterate's until the top of the next iteration rebuilds it, so the worst-row
  line printed in the iteration the return happens in describes the iterate
  before it. This is the behaviour the stagnation restart beside it has had
  since section 126 and it is not new here.

- **`backup/regression/atomic_elem_newton` does not reproduce the pass-1
  numbers `docs/lhs1140b_stationary_L4g_20260914.md` section 7.3 records for
  its measured build** (2.16e-09 / 4.24e-13 / 1.55e-08). The control build of
  THIS item, which carries the L4g rule, gives 3.08e-02 / 1.90e-05 /
  8.34e-01. Either the source moved between L4g and 2026-09-15 or the fixture
  configuration reconstructed here (`Load IC? True`, `Coupled carrier solve:
  False`, `Restart intent: stationary`, `EXHALE_OUTER_PASSES=12
  EXHALE_JFNK_MAXIT=80 EXHALE_DIFF_OMEGA=0.5`) differs from the one L4g ran.
  **This was not chased.** The control/measured comparison of section 5.2 is
  unaffected: both builds ran that same configuration.

- **The regression case directories `backup/regression/arm*`** (`armA_LW`,
  `armD_D2`, `armHeH_*`, `arm_heh1_x2matched`, ...) carry the noun "arm" that
  the naming rule forbids. Renaming them would move the reference data the
  goldens are keyed to, which is the kind of change that is asked about
  first; it is recorded here and not acted on.
  **DONE 2026-09-16** (PLAN_20260916_rev3 section 10, user's decision): all
  sixteen were renamed. None of them is in `DEFAULT_CASES` or has a golden
  twin, so no reference data moved; the mapping is
  `docs/named_case_audit.md` section 6.

---

## 7. What is kept, and what is not solved

All runs of this item are stopped. What is kept:

| directory | what it holds |
|---|---|
| `LHS1140b/models/.L4h/{x002_HeH2.13,x002_HeH9.7,p002_HeH9.7}` | the three 0.02 rungs with the delivered binary, 7, 2 and 8 outer passes, the evidence of section 5.3 and of plan item L17 |
| `LHS1140b/models/.L4h/ctl_x002_HeH2.13` | the control of section 2 |
| `scratchpad/reg/elem_*` and `scratchpad/cno_*` | the element reload and the C/N/O near state in every configuration sections 5.1, 5.2, 5.5 and 5.7 quote |
| `scratchpad/regw_{ctl,meas}` | `wasp_full_newton` control and measured, section 5.6 |

**The three 0.01 catalog cases are NOT SOLVED**:
`atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`,
`atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` and
`atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`. The XUV ladder of L14
reaches 0.03 on all three and refuses at 0.02. **L4h removes the pseudo-time
collapse from that refusal and does not remove the refusal**: the cause is
the one plan item L17 records, the base being hydrostatic to Mach 1e-8 at
0.02 XUV with alternating-sign velocities and the HLLC contact selection
putting every base cell on a kink of the residual. L17 is closed as
diagnosed, not solved; a smooth low-Mach flux would be a discretization
change and was not opened.

They are **not solved by the routes tried** -- these seeds, these pseudo-time
starts, this flux -- and that is all the evidence supports. No result here
shows the columns are unreachable.
