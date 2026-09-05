# The scaled trust region re-measured on the section 158 tree (2026-09-04)

**Judgement first. The trust region is the better-behaved control by every
measure section 153 claimed for it, and it is not the obstruction.** Both
controls stall in the same place, on the same row, at the same radius: the
carrier row at the H2 dissociation front. Section 158 fixed a real defect in
how that row is transported and did **not** move where the coupled solve fails.

## What was run

One march, two arms, one binary. The arms differ by the environment variable
`EXHALE_TRUST_REGION` and by nothing else -- same executable, same `input.inp`,
same restart files -- so any difference between them is the step control and
not the state, the physics or the build.

* **march** -- `mol_carrier` physics with `Coupled carrier solve: True`,
  `Secondary_ionization: Immediate`, `Stall [tol,N]: 0.0 2000`, no `Solver:`
  line so that it marches only. 200,000 steps, 8 threads, 1 h 36 m 41 s,
  ending at `du = 6.51e-3` and `log10 Mdot = 10.43`.

  Section 157.1 marched the same case on the pre-158 tree to `du = 5.89e-3`
  and `log10 Mdot = 10.43` in 1 h 37 m. **The hand-off state reproduces**: the
  carrier operator of section 158 changed where a 12,000-step snapshot lands
  (`mol_carrier`'s golden moved `10.50 -> 10.52`) without changing where a
  200,000-step march arrives.

* **ls** -- restart of that state, `Load IC? True`, `Solver: Newton 1.0` so the
  JFNK takes it at once. Trust region off: the backtracking line search.
* **tr** -- the same restart with `EXHALE_TRUST_REGION=1`.

Both arms hand off at step 2 and, on failure, march 2,000 steps and retry,
three attempts in all -- the existing recovery, unchanged.

Logs: `.../scratchpad/p159/{march,ls,tr}/`. The acceptance criteria were fixed
in `.../scratchpad/p159/PLAN.md` **before** the arms were run.

## The ledgers

Per JFNK attempt, in order:

| | ls (line search) | tr (trust region) |
|---|---|---|
| `info` | 2, 2, 2 | 2, **1**, 2 |
| best `\|\|R\|\|` | 1.345, 1.110, 1.095 | 1.462, **0.916**, 1.117 |
| flux spread | 1.07e-1, 9.14e-2, 9.49e-2 | 9.53e-2, 1.06e-1, 9.09e-2 |
| carrier row (volume-weighted) | 3.97e-3, 4.83e-3, 4.67e-3 | 5.51e-3, 5.67e-3, 5.76e-3 |
| **trials refused, element budget** | **78, 167, 148** | **0, 0, 0** |
| **trials refused, uncertified composition** | **78, 167, 148** | **0, 0, 0** |
| steps rejected | (not a line-search concept) | **0, 0, 0** |
| model refused by the ray test | -- | **0, 0, 0** |
| steps accepted | -- | 40, 500, 42 |
| residual evaluations | 2010, 2085, 2082 | 1307, 16258, 1371 |
| iterations with `lam <= 1e-6` | **25** | **0** |

Worst-row histogram over every JFNK iteration of the run:

| | k = 3 (energy) | k = 4 (carrier) |
|---|---|---|
| ls | 5 | 115 |
| tr | 0 | 580 |

## What the numbers say

**1. Section 153's claim is confirmed, and by a wider margin than the draft
claimed.** The line search refuses 78 to 167 trial states per attempt on the
element budget, and the same count again for an uncertified composition, with
the step length collapsing to `lam = 9.54e-7` on 25 iterations -- the exact
signature section 157 recorded. The trust region refuses **nothing**: no trial
outside the element budget, no step rejected, no model refused by the ray test,
in 582 accepted steps. Admissibility enforced inside the step does what it was
built to do.

**2. Its model is better than the draft measured.** Section 153 reported the
ratio of actual to predicted reduction at median 0.93 with 42 of 47 iterations
in [0.5, 2]. On this tree it is `ratio = 1.000` with `pred` and `actual` equal
to four digits on every iteration printed.

**3. Neither control converges, and both fail in the same place.** The worst
row is the carrier row (`k = 4`) in 115 of 120 line-search iterations and in
**580 of 580** trust-region iterations, at `j = 224-234`, `r = 1.126-1.143`.
That is the H2 dissociation front. It is where section 157 found the
obstruction on the pre-158 tree, and it has not moved.

**4. So section 158 fixed a defect that was not the obstruction.** That is
worth stating plainly rather than leaving implied. The carrier row was being
transported in a variable the mass row does not conserve; that was physically
wrong and is fixed; and the coupled solve fails in the same cells it failed in
before. The two facts are independent and both are results.

**5. The trust region reaches the lowest residual seen on either arm** --
`||R|| = 0.916`, in the attempt that ran to the 500-iteration cap accepting
every one of its 500 steps -- at eight times the cost (16,258 residual
evaluations against about 2,000). It is progress, and it is slow.

Per P17, no magnitude from either arm is quoted as physics: neither reached
`info = 0`, and the mass-loss rates the two arms end on (`10.41`, `10.42`) are
recorded here only as the two arms' end states, not as a rate.

## Why the trust radius contracts: the leg belongs to a different model

The first version of this memo recorded the contraction as unexplained; the
second blamed the Krylov tolerance of the dogleg's Newton leg. **Both were
wrong, and the second was refuted by its own measurement.** The full account,
with the trace lines, the four-rung `EXHALE_TR_GMRTOL` ladder that changed
nothing (mean GMRES iterations 1 -> 5, `pred <= 0` fixed at 22 of 40, radius
trajectory and finishing residual identical), the shift identity

```
[TR-trace] ||r0||= 8.684E-01  idtau*||(D/Drow)sN||= 1.101E+00  ||r0+AsN||= 1.100E+00
```

and the unshifted-leg fix, is section 159.3 of `docs/Update_EXHALE.md`. In one
paragraph: the leg solved the pseudo-transient-shifted system while the
predicted reduction was measured in the unshifted model, so for an accurate leg
the model residual at full length is exactly the shift term -- which sat just
above `||r0||`, made the model predict a rise on half the iterations, and
walked the radius down by the safeguard-shrink-then-boundary-double compound.
With the leg unshifted (`A sN = -r0`, preconditioner keeps its shift),
`pred <= 0` drops from 22 of 40 to **0 of 40**, the radius grows
`3.96e-2 -> 1.98`, and the one rejection that occurs is the textbook kind. The
unshifted leg is now the trust region's default; `EXHALE_TR_SHIFTED_LEG=1`
restores the old construction for measurement.

**And the solve still stagnates** (`info = 2`, best `||R|| = 1.462` on the
first attempt): the merit falls at every accepted step while the acceptance
norm's best iterate does not improve. A third, differently-repaired step
control stalling at the same carrier row confirms again that the obstruction
is the molecular layer, not the step control.

## Status

Section 153 is in the tree, environment-gated (`EXHALE_TRUST_REGION`), default
off, and adds no input key. With it off the regression is 11 of 11
byte-identical, so the physics of every existing run is untouched.

Its patch predated section 147, which split `eval_residual`'s composition seed
from its result; the five call sites it adds were translated onto the current
interface, following the idiom the neighbouring line search already uses
(`seed = f_sp`, result into a workspace). The accepted evaluation is the one
exception: its result goes to `f_sp_acc`, so that the composition of the
accepted state still reaches the caller, which is what the old inout workspace
did.
