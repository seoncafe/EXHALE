# A stagnant metal-cooled layer below the escape radius, and what it did to the JFNK solve

HD 209458 b with molecular chemistry AND solar trace metals
(`examples/16_molecular_metals`) does not reach a Newton-grade steady state
with the prescription that works for the metals-off twin
(`examples/15_molecular`): at `Resid tol: 1.0e-5` the JFNK finish aborts with
`info = 2` after twelve consecutive line searches that find no descent step,
from a cold start and from a warm one alike. Sections 1 and 3 record what was
measured about why, section 6 where the case stood after the first three
2026-08-13 changes, section 7 what is left, and section 8 the fourth change —
a conservative artificial dissipation of the stagnant layer that the marching
loop and the steady residual both see, under which the warm restart reaches
`info = 0` at `Resid tol: 1.0e-5` for the first time. Section 9 (2026-08-15)
closes three of section 7's items by measurement: the cold path does not
respond to the damping coefficient and is abandoned, the warm prescription is
pinned to a seed converged with the current binary, and widening the acceptance
window to include the layer is rejected — with what the layer's residual
actually is, now printed by every run, and what its weak determination costs in
the Balmer lines.

This memo records what was measured, what was fixed on the strength of it, and
what is still open. Diagnostic runs live under the scratch tree
`mm16/` referred to below (`base/` cold, `warmstall/` a warm restart from the
stalled state, `warm15/` a warm restart from the converged metals-off solution,
`resid_nm/` and `resid_wm/` single residual dumps).

## 1. The layer

Turning solar C/N/O on inverts the local energy balance in a thin shell just
above the base. At `r = 1.02 R_p`, from the two runs' own `Hydro_ioniz.txt`:

| run | heat [erg cm^-3 s^-1] | cool | cool/heat | T [K] | v [cm/s] |
|---|---|---|---|---|---|
| `examples/15` (metals off, converged) | 1.019e-7 | 3.095e-9 | 0.030 | 2206 | +159 |
| `mm16/base` (metals on, cold, stalled) | 2.851e-8 | 2.769e-7 | 9.71 | 1872 | -47.9 |

Cooling exceeds heating by a factor 10 where without metals heating exceeded
cooling by a factor 33, and the flow there is infalling rather than outflowing.

The layer that results is nearly hydrostatic, and it carries a `2 dr` odd-even
mode. Over `r = 1.015-1.030`, with the alternating component of a profile
`f` measured as `max_j |f_j - (f_{j-1} + f_{j+1})/2|`:

| state | alternating `v` / mean `\|v\|` | alternating `T` / mean `T` | mean `\|v\|` [cm/s] |
|---|---|---|---|
| `examples/15` converged (metals off) | 0.001 | 0.0001 | 2.07e2 |
| `mm16/warm15` (16 from the 15 solution) | 0.068 | 0.0015 | 1.57e1 |
| `mm16/warmstall` | 4.70 | 0.057 | 6.16e0 |
| `mm16/base` (16 cold) | 9.98 | 0.154 | 2.36e1 |

The mean speed in that shell drops by one to two orders relative to the
metals-off solution, and the alternating component of the velocity grows to
several times the mean. The inviscid HLLC flux does not damp a contact-mode
oscillation, so nothing in the scheme removes it once the layer stops flowing.
[2026-08-15: "nothing in the scheme removes it" is no longer true. The gated
fourth-difference dissipation of section 8 (`Low-Mach damping`,
`src/modules/flux/low_mach_dissipation.f90`) damps exactly this contact mode;
it is off by default.]

## 2. The solve accepted on one functional and stepped on another

`resid_relnorm` — the quantity `Resid tol` is compared against, and the one
that decides `info` — is a volume-weighted ratio over `[j_min:N]`, i.e. over
`r >= r_esc` only (2.00 R_p here, out of a domain reaching 4.15 R_p). The
JFNK line search minimized a different functional: the scaled 2-norm
`\|\|F/D\|\|_2` over the WHOLE domain `1..N`, which in this configuration is
dominated by the stagnant shell of section 1 and by the first two cells.

On the warm restart from the stalled state (`mm16/warmstall/run.log`) the two
move in opposite directions:

```
 (JFNK) start ||R||=  1.886E-01  ||Fs||2=  5.70E+00
 (JFNK) it   1  ||R||=  4.025E-03  ||Fs||2=  1.57E-01
 (JFNK) it   2  ||R||=  3.958E-05  ||Fs||2=  5.67E-01
 (JFNK) it   3  ||R||=  5.270E-04  ||Fs||2=  2.43E-01
 ...
 (JFNK) it 117  ||R||=  9.548E-05  ||Fs||2=  1.58E-02
 (JFNK) line search found no descent step in 12 consecutive iterations -- aborting
```

Iteration 2 produced a state with `||R|| = 3.958e-5`. Over the following 115
iterations the line search cut the whole-domain merit from 5.67e-1 to 1.58e-2,
a factor 36, while `||R||` rose to 9.5e-5. The solve reached a good wind
solution and then walked away from it, by design.

For contrast, on the metals-off `examples/15` run the two functionals move
together and the minimum of both is the last iterate (`it 272`,
`||R|| = 9.542e-6`, `info = 0`).

## 3. The equilibrium solve inherited a cell's previous composition

When the coupled molecular + metals equilibrium produced no admissible root in a
cell, `ioniz_eq` kept that cell's incoming fractions and warned. Counting those
warnings inside the JFNK phase of each log, normalized by the number of outer
iterations:

| run | JFNK iterations | "molecular equilibrium failed" messages | per iteration |
|---|---|---|---|
| `examples/15` (metals off) | 272 | 214 | 0.8 |
| `mm16/warm15` | 125 | 13035 | 104 |
| `mm16/warmstall` | 117 | 43804 | 374 |

With metals on, essentially every residual evaluation carried at least two cells
whose composition came from an earlier evaluation rather than from a solve at
the state being evaluated. That composition is not the equilibrium of the cell
at that state, and the error is not small: it is whatever the last accepted
iterate happened to leave there.

**What that does and does not break.** It does not corrupt the Newton model
within an outer iteration: `solve_steady_jfnk` evaluates the frozen-weight
residual, every colored Jacobian probe (`build_banded_jac_full` restores
`fwork = f_sp_base` per probe), every `J*v` product and every line-search trial
from the SAME base `f_sp`, so the Jacobian columns and the merits being compared
all sit at one composition. (An earlier version of this section said they mixed
two states; that was wrong.) What it does break is the identity of the map
between outer iterations: `f_sp` is overwritten by the accepted trial's
composition, so the residual the next iteration differentiates and searches
differs from the previous one, by a finite amount, at exactly those cells. The
`||R||` jump from 3.958e-5 at iteration 2 to 5.270e-4 at iteration 3 in
section 2 is that effect.

**What the cells were actually failing on.** Measured with a counting build on
`mm16/cb16` (HD 209458 b, molecular + solar C/N/O, warm restart, 40 marching
steps, single-threaded): 224 cell solves ended below the solver tolerance,
223 of them keeping a physical but unconverged root and 1 inheriting its
previous composition. The roots that were rejected were not on a different
branch: at step 1, cell `j = 117` (`r = 1.0294`, `T = 1652 K`) the first attempt
converged (`info = 1`) to a root whose only violations were
`x(He II) = -3.8e-6` and `x(H2) = -4.5e-4`, and the molecular-basin retry
returned a strongly molecular root, `x(H2) = 0.998`, violating the simplex by
`x(He II) = -1.2e-8` — round-off around a fully neutral helium. Both were
rejected and the cell kept its previous composition. The failure is a root
resolved onto a face of the simplex, not an unphysical branch.

Section 5 item 3 records what was done about it.

## 4. Aligning the merit's region was tried and rejected

The obvious response to section 2 is to make the line search minimize the same
functional the solve is accepted on: the volume-weighted RMS of `F/D` over
`[j_min:N]` instead of `||F/D||_2` over `1..N`. That was implemented and run
against the configurations that converge today. It breaks them:

| case | whole-domain `\|\|F/D\|\|_2` | `[j_min:N]` volume-weighted RMS |
|---|---|---|
| `examples/15` (HD 209458 b, molecular, cold) | `info = 0`, `\|\|R\|\| = 9.542e-6`, 272 iterations | `info = 2`, best 1.171e-4, 30 iterations |
| `photo_deep_secion_cont` (HD 209458 b, `base.inp` handoff, warm) | `info = 0`, `\|\|R\|\| = 7.284e-4`, 4 iterations | `info = 2`, best 3.249e-3, 15 iterations |

This reproduces, on two further cases, the measurement recorded in
`docs/newton_scaling_and_base_wall.md` section 4 on 2026-08-11, which rejected
the same change for the same reason. The reading given there survives the
present case: the Newton step is computed from the full residual over ALL rows,
so a merit that ignores most of those rows rejects the steps that step takes.
The region mismatch between the merit and the convergence measure is therefore
kept, and its one real cost — that the iterate with the smallest `||R||` need
not be the last one — is handled in the bookkeeping instead.

## 5. What was changed (2026-08-13)

Items 1 and 2 are in `docs/Update_EXHALE.md` section 58, item 3 in section 59.

1. **A solve that reached the tolerance is reported as converged.** The solver
   already kept the best-`||R||` iterate, already returned it on failure, and
   already tested the tolerance at the top of every outer iteration; it now
   also judges `info` on the iterate it returns, which closes the two
   remaining holes (the restore was gated on `info /= 0`, and a solve whose
   last allowed iteration met the tolerance fell out with `info = 1`). This is
   a narrow fix: measured on the stalled warm restart at `Resid tol: 5.0e-5`,
   the binaries before and after hand off at the same marching step, take the
   same 2 outer iterations, print the same iterate lines, and both end
   `info = 0` at `||R|| = 3.958e-5`. It is not what makes section 2's state
   usable — the top-of-loop test already caught it — and at
   `Resid tol: 1.0e-5` the configuration genuinely does not converge.

2. **`load_IC` no longer accepts a metal column that is present but zero.** It
   decided whether the restart file carried an element by looking only at
   whether the column existed. The schema-2 writer emits the metal columns
   unconditionally, so a metals-OFF output has them present with the value
   zero, and a metals-ON restart from such a file ran with **zero metal density
   everywhere** while the base boundary condition counted metals in the mass
   and particle budget. Such a restart reported `||R|| = 9.54e-6` on the
   `examples/15` solution — it was reproducing the metals-off solution under a
   metals-on label. `load_IC` now also tests that the loaded stages are not
   identically zero, rebuilds such an element from the abundance exactly as a
   cold start does, and warns on stdout and in the setup report.

3. **The molecular equilibrium no longer inherits a cell's previous
   composition.** Two things replace it, both in
   `src/modules/radiation/ionization_equilibrium.f90`. The molecular-basin
   retry seed is now a state: the H2 dissociation equilibrium of the local
   `(p, T)` together with every element in its own ionization balance at the
   incoming `n_e` (`dissociation_ionization_balance_at_fixed_ne`), where it
   used to combine a chemical-equilibrium H2 fraction with the *previous
   state's* ionization fractions, a pair that can put more H nuclei in H2 and
   H II together than the cell has. And when none of the three attempts lands
   inside the simplex, the cell now keeps the root that lies closest to it
   (`element_budget_violation`) and clamps that root onto the element budget
   (`clamp_fractions_to_element_budget`) instead of restoring its incoming
   fractions. Both outcomes are computed from the state being evaluated;
   neither carries anything over from an earlier evaluation.

   *Measured*, on the same `mm16/cb16` run as section 3: cell solves ending
   below the solver tolerance fall from 224 to 77, cells inheriting a previous
   composition from 1 to 0. On `examples/15` the solve gets better, not just
   different: `info = 0` at `||R|| = 6.105e-6` in **89** outer iterations
   against `info = 0` at 9.542e-6 in 272, hand-off at marching step 25839
   against 25837, `log10 Mdot = 10.24` against 10.23, H2 = H I front at
   1.0098 R_p against 1.0084.

   *One thing was tried and reverted.* The simplex test rejects a root for a
   component below `-1e-10`, and the rejected roots of section 3 sit at
   `-1e-8`, i.e. inside the solver's own `xtol = sqrt(machine epsilon)`.
   Widening the test to that band is defensible on paper and it is worse in
   practice: it lets a cell stop on a root sitting just outside a face instead
   of retrying from another starting point, and the retried root is the better
   one. With the band widened, `examples/15` goes from `info = 0` at 6.105e-6
   in 89 outer iterations to `info = 2` at 1.112e-4 in 133. The band stays at
   `1e-10`, and the reason is now recorded at the site.

   *The two changes are needed together*, and `examples/15` is a marginal solve
   that responds to either one alone by failing. Same case, same cold start,
   4 threads:

   | simplex band | attempt-2 seed | no admissible root | `examples/15` |
   |---|---|---|---|
   | 1e-10 | previous-state mix | inherit (before) | `info = 0`, 9.542e-6, 272 it |
   | 1e-10 | balance | clamp (**adopted**) | `info = 0`, 6.105e-6, 89 it |
   | 1e-10 | previous-state mix | clamp | `info = 2`, 3.550e-4, 114 it |
   | 1e-10 | balance | inherit | `info = 2`, 1.359e-5, 166 it |
   | xtol | balance | inherit | `info = 2`, 1.359e-5, 166 it |
   | xtol | balance | clamp | `info = 2`, 1.112e-4, 133 it |

   The band is the only one of the three that never helps. Section 6 runs the
   same matrix on `examples/16`, where the ordering is different again.

## 6. Where `examples/16` now stands

With the item-1 and item-2 changes, warm-started from the converged
`examples/15` state — the restart the `load_IC` fix repairs — with
`Resid tol: 5.0e-5`, `Solver: Newton 5.0e-2`, `du_th [PLM,WENO3]: 0.5 1.0e-3`,
the case converged: hand-off at marching step 50668, `info = 0`,
`||R|| = 3.337e-5` in 30 outer iterations, `log10 Mdot = 9.65`. That was the
first Newton-grade steady state for this configuration, and it is reproducible
(re-run 2026-08-13 with that binary: the same 50668, 3.337e-5, 30, 9.65).

**Item 3 costs it.** With the equilibrium no longer inheriting a previous
composition, the same warm restart bottoms at `||R|| = 7.931e-5` in 67 outer
iterations and reports `info = 2`; `Resid tol: 2.0e-5` and `1.0e-5` produce the
identical trajectory, so nothing is gained by tightening. The cold path bottoms
at 1.081e-3 in 500 outer iterations against 5.038e-4 in 281 before. Neither
path converges now.

This is not an argument for the inheritance: it is what a marginal solve does
under an O(1) change to a handful of cells' composition. The two HD 209458 b
cases respond in opposite directions, which is what "marginal" means here:

| variant | `examples/15` (cold) | `examples/16` warm, `Resid tol: 5.0e-5` |
|---|---|---|
| before item 3 | `info = 0`, 9.542e-6, 272 it | `info = 0`, 3.337e-5, 30 it |
| new retry seed only, inheritance kept | `info = 2`, 1.359e-5, 166 it | `info = 0`, 3.059e-5, 30 it |
| clamp only, old retry seed | `info = 2`, 3.550e-4, 114 it | `info = 0`, 4.853e-5, 93 it |
| **adopted** (new seed + clamp) | `info = 0`, 6.105e-6, 89 it | `info = 2`, 7.931e-5, 67 it |

Each half of the change on its own converges `examples/16` and fails
`examples/15`; together they converge `examples/15`, better than before, and
fail `examples/16`. `log10 Mdot` is 9.64-9.65 in every row that converges.

What did not change is that the residual for `examples/16` is dominated by the
stagnant shell of section 1, at `j = 93`, `r = 1.021`, in every one of these
runs. That shell, not the equilibrium bookkeeping, is what the case is waiting
on.

## 7. Still open

- **Damping the stagnant layer — addressed, see section 8.** Nothing in the
  inviscid scheme removes the `2 dr` contact mode once the metal cooling stops
  the flow at `r ~ 1.02`. The Shapiro filter and explicit viscosity were both
  measured against a different symptom (the base momentum row) and refuted
  there (`TO_BE_DONE.md` item (A)). A conservative artificial dissipation that
  the marching loop and the steady residual both see is now available as
  `Low-Mach damping` (default off), and section 8 records what it does to this
  configuration.
- **`examples/16` itself — the warm path settled, the cold path closed, see
  sections 9.2 and 9.3.** Section 3's inheritance is gone (section 5 item 3),
  and the case is no better for it. With `Low-Mach damping` at
  `eps4 = 5e-3`-`1e-2` the warm restart now reaches `info = 0` at
  `Resid tol: 1.0e-5` (section 8.2), which it never did before; the cold path
  and the larger coefficients still do not. Measured 2026-08-15: no coefficient
  in `off`-`5e-3`-`1e-2`-`2e-2` rescues the cold path, which is therefore
  closed (section 9.2), and the warm prescription of section 8.2 reproduces
  only from a seed converged with the *current* binary (section 9.3).
- **Whether the acceptance window should reach below `r_esc` — measured and
  rejected, see section 9.4.** The residual of the layer is now printed by
  every run, and at convergence it is 0.2 to 10 in the relative norm the
  window uses, on cases that converge today including the metals-off
  `examples/15`. Widening the window on this norm would fail them all rather
  than determine the layer.
- **The layer is only weakly determined, and it reaches the Balmer lines —
  documented, not fixed, see section 9.5.** Two solves that both return
  `info = 0` at `Resid tol: 1.0e-5` differ inside the layer by tens of percent
  in `rho` and `T`, and that propagates to 4.8% relative in the H-alpha peak
  and 8.4% in H-beta. He I 10830 (0.2%) and the wind are insensitive. Nothing
  measured so far pins the layer down; the honest handling is to quote the pair
  and carry the caveat.
- **Cells that end below the solver tolerance — still open.** The equilibrium
  still returns,
  in a few cells per evaluation, a root the solve did not resolve to `xtol`
  (77 cell solves over the 40-step `mm16/cb16` run, down from 224). Those roots
  are computed at the state being evaluated, so they are not the path
  dependence of section 3, but they do depend on which starting point the
  ladder reached them from. `hybrd1` returns `info = 4` there — no progress —
  on unknowns spanning `1e-14` to `1`, with `mode = 2` and `diag = 1`, i.e. no
  variable scaling at all. Scaling the molecular unknowns has not been tried.

## 8. Damping the layer, and what it did (2026-08-13)

`Low-Mach damping: <eps4> [<M_th>]` (default off) adds a gated fourth-difference
stress to the numerical momentum flux, and its work term to the energy flux,
inside `RK_rhs` — the one routine `assemble_residual`, and therefore the JFNK
residual, also calls. Statement, derivation and stability bound:
`src/modules/flux/low_mach_dissipation.f90`; write-up:
`docs/Update_EXHALE.md` section 60.

Where the term is applied is the point of it. A Shapiro filter smooths the
marching state; the Newton residual never sees it, so the mode it suppresses
while marching is still an undamped mode of the system Newton solves, and the
two disagree about what a steady state is. A term inside the flux keeps one
equation: the fixed point of the marching loop is the zero of the residual.

### 8.1 How big the term is, and where it is exactly zero

Measured on the converged profile of each run. The first three columns are what
the code itself prints at the end of any run that uses the key; recomputing them
in cgs from `Hydro_ioniz.txt` independently of the Fortran agrees to the printed
digits.

| run | `eps4` | peak `\|D_p\|/\|rho v^2+p\|` | at | gate open out to | first `M > 0.1` |
|---|---|---|---|---|---|
| `examples/15` (metals off) | 2e-2 | 1.83e-4 | 1.0005 | 1.043 R_p | 1.76 R_p |
| `examples/16` warm | 1e-2 | 8.36e-5 | 1.0005 | 1.108 R_p | 1.99 R_p |
| `examples/16` warm | 5e-3 | 3.94e-5 | 1.0005 | 1.105 R_p | 1.99 R_p |
| hot Uranus (molecular + solar metals) | 2e-2 | 3.70e-5 | 1.0005 | 1.083 R_p | — |

The peak sits at the first interior face, where the base boundary condition
forces the steepest velocity gradient in the domain, and it is 0.02% or less of
the physical momentum flux there. Through the stagnant shell itself
(`r = 1.005-1.04`) the ratio is `1e-7` or below. Above `M_th = 1e-3` the gate is
identically zero, and the first cell with `M > 0.1` is at 1.76-2.00 R_p, so the
wind, the `r_esc = 2 R_p` window over which `||R||` is measured, and the `N-20`
cell at which `Mdot` is evaluated all sit where the term does not exist.

### 8.2 `examples/16`, warm restart from the converged `examples/15` state

All rows: `Resid tol: 1.0e-5`, `Solver: Newton 5.0e-2`,
`du_th [PLM,WENO3]: 0.5 1.0e-3`, `Max steps: 150000`, 4 threads, same binary,
same restart files.

| `eps4` | hand-off step | `info` | best `\|\|R\|\|` | outer it | `log10 Mdot` | H2 = H I front |
|---|---|---|---|---|---|---|
| off | 48060 | 1 | 5.216e-5 | 500 (cap) | — | — |
| 5.0e-3 | 48062 | **0** | **9.934e-6** | 162 | 9.65 | 1.0112 |
| 1.0e-2 | 48007 | **0** | **9.987e-6** | 202 | 9.67 | 1.0117 |
| 2.0e-2 | 65836 | 2 | 2.529e-3 | 237 | — | — |
| 4.0e-2 | 66015 | 2 | 2.870e-3 | 66 | — | — |

This is the first `info = 0` at `Resid tol: 1.0e-5` this configuration has ever
produced. The two converged rows stop where they do because the solver tests the
tolerance at the top of each outer iteration and exits on the first iterate
below it — 9.93e-6 and 9.99e-6 are "the first crossing", not a knife edge. The
key-off control, with the same restart and the same binary, spent its full
500-iteration budget and plateaued at 5.2e-5. Nothing diverges in the rows that
do not converge: each falls back to time-marching, as it does with the key off,
and runs out the 150000-step cap at a `du` plateau of 1.0e-3 to 2.6e-3,
`log10 Mdot` 9.58-9.66. Those are marching states, not Newton-grade ones, which
is why the table leaves their `Mdot` blank.

**The coefficient has a working range and this case sits near its top.** At
`eps4 >= 2e-2` the solve is worse, not better, and the reason is visible before
the Newton phase: the stress changes the base solution enough to slow the `du`
descent (`du` at step 40000 is 0.139 with the key off, 0.360 at `eps4 = 2e-2`),
so the hand-off comes 18000 steps later, from a different state. The two
coefficients that work are the bottom of the classical JST range; the top of it
is already too much for this configuration. **Whoever uses the key on a new
configuration has to run the pair.**

**The two converged solutions are not the same solution.** Halving `eps4`
from 1e-2 to 5e-3 moves `log10 Mdot` by 0.02 (5%) and the H2 = H I front by
0.0005 R_p, and in the wind (`r > 1.2 R_p`) `rho` by <= 7.4%, `v` by <= 3.7% and
`T` by <= 0.9%. Inside the stagnant layer they differ a great deal — `rho` by
137% at `r = 1.013`, `T` by 38% at `r = 1.015` — and they differ in what the
layer looks like: over `r = 1.015-1.030` the alternating component of `v`
relative to the mean `|v|` is 0.024 at `eps4 = 1e-2` and 4.29 at 5e-3. Both
states satisfy `||R|| < 1e-5` on the wind window. The layer itself is therefore
still weakly determined, and the residual measure does not see it.

### 8.3 The cold path, and the two configurations that must not move

The cold `examples/16` start (`Resid tol: 5.0e-5`) improves and still fails:
`info = 1` at `||R|| = 9.307e-4` with `eps4 = 2e-2` against `info = 1` at
1.081e-3 with the key off, from the identical hand-off step 112920 — the term is
inert through the cold relaxation, when the flow is fast enough everywhere to
close the gate, and only acts at the end.

The two molecular configurations that do converge today are not disturbed
(`Resid tol: 1.0e-5`, cold, `eps4 = 2e-2` where on):

| case | key | hand-off | `info` | `\|\|R\|\|` | outer it | `log10 Mdot` | front |
|---|---|---|---|---|---|---|---|
| `examples/15` | off | 25839 | 0 | 6.105e-6 | 89 | 10.24 | 1.0098 |
| `examples/15` | on | 25839 | 0 | 3.642e-6 | 100 | 10.24 | 1.0095 |
| hot Uranus | off | 34537 | 0 | 5.826e-6 | 23 | 10.30 | 1.0774 |
| hot Uranus | on | 34537 | 0 | 6.252e-6 | 23 | 10.30 | 1.0774 |

Both hand off at the same marching step and reach the same `Mdot`. On the hot
Uranus the whole profile moves by at most 7.3e-4 relative in `rho` and 7.1e-4 in
`T`. On `examples/15` the wind above `r = 1.2 R_p` moves by at most 1.2% in
`rho`, 0.5% in `v` and 0.16% in `T`, while the first few cells above the base
move much more (`rho` by up to 16%, `T` by 6.9% below `r = 1.05`) — those cells
carry an odd-even velocity oscillation of 86 times the local mean with the key
off, which the term reduces to 77, so they are the cells the term is for and
they are not pinned down to begin with.

## 9. The cold path closed, the acceptance window question closed, and the seed pinned (2026-08-15)

Section 7 left three measurable questions: whether any damping coefficient
rescues the cold start, whether the acceptance window should be widened to
include the layer, and what the weak determination of the layer costs. All
three were measured on 2026-08-15 and are recorded here; the section also
records a diagnostic added to make the third visible in every run, and a
Lyman-Werner A/B on the metals-on configuration.

Every run in this section: the `examples/16` configuration with
`Solver: Newton 5.0e-2`, `Resid tol: 1.0e-5`, `du_th [PLM,WENO3]: 0.5 1.0e-3`,
`Max steps: 150000`, 4 threads, the same binary (the one carrying the
diagnostic of section 9.1), in scratch run directories outside the repository.

### 9.1 The residual of the layer is now printed

`resid_relnorm_below_escape` in `src/modules/time_step/steady_newton.f90`
applies the *same* volume-weighted relative norm the convergence measure uses
to the complementary region `[1:j_min-1]`, i.e. below the escape radius, and
also locates the cell with the largest cell-wise relative residual
`|F_k,j| / |u_k,j|` in that region, reporting `j`, `r(j)` and which of the
three rows (mass, momentum, energy) it is. `write_resid_below_escape` prints
one line at four points: the start and the end of the PTC solve and the start
and the end of the JFNK solve. The JFNK end line is emitted *before* the
best-iterate restore, because that is the last point at which `F` and `u` are a
consistent pair, so it refers to the last iterate visited — which is also the
returned state whenever no restore happens.

The value is printed and never tested; nothing in the solve depends on it. The
change is +95 lines in that one file and `make check` is **5/5 byte-identical**
(`wasp_full`, `wasp_he23off`, `mol_base_handoff`, `mol_metals`,
`mol_lyman_werner`) — the output files are unchanged and only stdout grows.

### 9.2 The cold path does not respond to the coefficient, and is closed

Cold start (`Load IC? False`), scanning the `Low-Mach damping` coefficient:

| `eps4` | JFNK outcome |
|---|---|
| off (section 8.3) | `info = 1`, best 1.081e-3 |
| 5.0e-3 | `info = 1`, best 9.953e-4 |
| 1.0e-2 | `info = 1`, best 1.025e-3; falls back to marching and runs out the 150000-step cap at a `du` plateau of ~7e-3 |
| 2.0e-2 (section 8.3) | `info = 1`, best 9.307e-4 |

All four bottom at the same ~1e-3 floor, a factor 100 above the tolerance, and
the two coefficients that work on the warm restart (5e-3, 1e-2) do nothing here
that `off` does not. The cold hand-off state appears to be outside the Newton
basin of the solution the warm restart finds, and the coefficient is not the
lever. **The cold path for this configuration is closed**: the warm restart is
the prescription.

### 9.3 The warm prescription needs a seed converged with the current binary

Section 8.2 says "warm restart from the converged `examples/15` state" without
saying *which* converged state, and that turns out to matter. Two seeds were
run A/B, same binary, same everything else:

- **Seed A** — the copy `examples/15_molecular/output` held until 2026-08-15
  (written 2026-08-13 10:15, i.e. by the binary from *before* the section 5
  item-3 change; H2 = H I front at 1.0089 R_p; now kept as
  `examples/15_molecular/output_pre_item3_20260813/`, the example directory
  itself holding the seed-B state since 2026-08-15).
- **Seed B** — `examples/15` re-converged with the current binary. That
  re-convergence reproduces the section 8.3 row exactly: hand-off at marching
  step 25839, `info = 0`, `||R|| = 6.105e-6`, `log10 Mdot = 10.24`, front
  (`n_H2 = n_HI` crossing) at 1.0097 R_p.

| seed | `eps4` | hand-off step | `info` | best `\|\|R\|\|` | outer it | `log10 Mdot` |
|---|---|---|---|---|---|---|
| A | 5.0e-3 | 50645 | **2** | 1.303e-4 | — | — |
| A | 1.0e-2 | 50629 | **2** | 5.770e-5 | — | — |
| A | 5.0e-3 + `Stellar LW flux: 343.0` | 50189 | 0 | 9.799e-6 | 47 | 9.68 |
| B | 5.0e-3 | 48062 | 0 | 9.934e-6 | 162 | 9.65 |
| B | 1.0e-2 | 48007 | 0 | 9.987e-6 | 202 | 9.67 |
| B | 5.0e-3 + `Stellar LW flux: 343.0` | 47965 | 0 | 9.900e-6 | — | 9.68 |

The two seed-B rows without Lyman-Werner reproduce the section 8.2 table
exactly, hand-off step included. The seed-A rows do not: the hand-off comes
~2600 marching steps later, from a different state, and the JFNK line search
then collapses (`lambda ~ 1e-6`) and aborts with `info = 2`; both fall back to
marching and run out the 150000-step cap at a `du` plateau of 1.0e-3 to 9.5e-4.

The reading is that the section 5 item-3 change moved the `examples/15`
solution by more than this marginal solve tolerates in its starting point, so a
restart file written by an older binary is no longer the same seed. **The
prescription is therefore: re-converge `examples/15` with the binary in hand,
and restart `examples/16` from that**, not from whatever `output/` happens to
hold. Note also that the Lyman-Werner member converges from *both* seeds — the
LW coupling appears to widen the basin for this configuration rather than
narrow it.

### 9.4 Extending the acceptance window below `r_esc` is rejected by measurement

With the diagnostic of section 9.1 the question of section 2 — should the solve
be accepted on a window that includes the layer? — becomes a measurement. The
volume-weighted relative norm over `[1:401]` (`j_min = 402` here):

| state | norm over `[1:j_min-1]` | worst cell |
|---|---|---|
| seed-B warm, at JFNK start | 0.60-1.21 | `j = 6`-15, `r = 1.001`-1.003, momentum |
| `examples/16` converged, `eps4 = 5e-3` | 0.195 | `j = 75`, `r = 1.0153` |
| `examples/16` converged, `eps4 = 1e-2` | 2.352 | `j = 5`, `r = 1.0010` |
| `examples/16` converged, `eps4 = 5e-3` + LW | 0.939 | `j = 3` |
| `examples/15` converged (**metals off**) | 9.99 | `j = 2` |
| cold `eps4 = 5e-3`, at failure | 594 | `j = 97`, `r = 1.0218` |
| cold `eps4 = 1e-2`, at failure | 1599 | `j = 112`, `r = 1.0273` |

A converged state carries a relative residual of 0.2 to 10 below the escape
radius — four to six orders above `Resid tol` — and the metals-off
`examples/15` solution, which nobody doubts, carries the largest of them. That
is not a statement about metals: the layer is near-hydrostatic (`|v| ~ 10`
cm/s), so the denominator of the momentum row's relative residual, `|rho v|`,
collapses, and the ratio is large wherever the flow is slow whether or not the
state is right. **Widening the acceptance window on this norm would turn every
configuration that converges today, `examples/15` included, into `info = 2`,
without determining the layer any better.** The window stays at `[j_min:N]`,
and the diagnostic makes the uncontrolled part visible in every run instead.

The failure rows are 3 orders above the converged ones, so the number does
separate a bad state from a good one; it is the absolute scale that is not
usable as a tolerance.

### 9.5 What the weak determination of the layer costs, in the lines

The two seed-B solutions of section 9.3 (`eps4 = 5e-3` and `1e-2`, both
`info = 0` at `Resid tol: 1.0e-5`) were compared cell by cell, and then run
through `EXHALE_transit.py`. In the layer, the `1e-2`/`5e-3` ratio of `rho` is
0.78 at `r = 1.005`, 1.31 at 1.013, 1.41 at 1.020 and 1.32 at 1.030; `T` differs
by +20% at 1.005, +38% at 1.015 and -7.7% at 1.020; over `r = 1.015-1.030` the
alternating component of `v` relative to the mean `|v|` is 0.392 at `5e-3`
against 0.028 at `1e-2`, on a mean `|v|` of 12-16 cm/s. (The 137% at
`r = 1.013` quoted in section 8.2 belongs to an earlier pair, of the seed-A
lineage; the size differs pair by pair, the character does not.)

In the wind the two agree: `Drho/rho` is +5.1% at 1.3 R_p, +3.4% at 1.5 and
+2.6% at 2.0, `Dv/v <= 1.6%`, `DT/T <= 0.4%`, and `log10 Mdot` is 9.65 against
9.67.

The transmission spectra, theoretical peak excess absorption:

| line | `eps4 = 5e-3` | `eps4 = 1e-2` | relative difference |
|---|---|---|---|
| He I 10830 | 40.19% | 40.26% | 0.2% |
| H-alpha | 2.098% | 2.199% | +4.8% |
| H-beta | 1.046% | 1.134% | +8.4% |
| Ly-alpha | saturated (100%) | saturated (100%) | — |

He I 10830 forms in the wind and does not care which of the two layers sits
underneath it. The Balmer lines do: part of the H(n=2) absorption is formed in
the layer, so the weak determination appears to propagate to the line depth at
the 5-8% relative level. **Any H-alpha or H-beta depth quoted from this
configuration should be quoted with that caveat, and the `5e-3`/`1e-2` pair
should be run rather than a single coefficient.**

### 9.6 Lyman-Werner photodissociation with metals on

The `Stellar LW flux [erg/cm2/s]: 343.0` A/B, on seed B at `eps4 = 5e-3`, both
members `info = 0`:

| | LW off | LW on (343) |
|---|---|---|
| `log10 Mdot` | 9.65 | 9.68 |
| H2 = H I front (`n_H2 = n_HI`) | 1.01108 | 1.01026 |
| `T` at `r = 1.02` [K] | 1801 | 1470 |
| `x(H2)` at `r = 1.02` | 4.0e-6 | 2.1e-5 |
| H3+ peak [cm^-3] | 1.93e5 at 1.0006 | 1.99e5 at 1.0006 |
| HeH+ peak [cm^-3] | 0.285 | 0.053 |
| C II fraction at 1.02 / 1.20 | 0.197 / 0.449 | 0.180 / 0.430 |
| H-alpha peak | 2.098% | 2.297% |
| H-beta peak | 1.046% | 1.221% |
| He I 10830 peak | 40.19% | 40.33% |

The self-shielding behaves as section 56 of `docs/Update_EXHALE.md` describes:
at the base `N(H2) = 2.68e21 cm^-2`, `f_shield = 2.11e-6` and
`k_LW = 1.27e-10 s^-1`, while above `r = 1.10 R_p` the shielding is gone
(`f_shield = 1.000`, `k_LW = 6.026e-5 s^-1`).

`x(H2)` at `r = 1.02` is *higher* with the band on, which at first reads
backwards. It appears to be the layer temperature: the LW-on layer is 330 K
colder there, and the temperature dependence of the H2 balance outweighs the
added destruction. The front moves down by 0.0008 R_p, consistent with more
destruction where the shielding has lifted.

**Caveat on the line depths.** Turning the band on changes the H-alpha peak by
+9.5% relative and H-beta by +17%, which is only about twice the numerical
spread of section 9.5. The direction is reproducible in this pair; the
magnitude should not be quoted as an observational prediction without the
coefficient pair run alongside it.

## 10. The layer's momentum residual on a physical scale (2026-08-15)

Section 9.4 ended on the observation that the below-escape number *separates*
a failed state from a converged one but is unusable as a tolerance, because
the momentum row's relative scale, `|rho v|`, collapses in a quasi-hydrostatic
layer: the ratio is large wherever the flow is slow, whether or not the state
is right. The diagnostic of section 9.1 has now been given a physical scale.

### 10.1 The change

In `resid_relnorm_below_escape` (`steady_newton.f90`) the momentum row is no
longer divided by `|u(2,:)| = |rho v|` but by the gravitational force density

    s_grav(j) = |rho(j)| * |Gphi_i(j) - Gphi_i(j-1)| / dr_j(j),

the same discrete potential difference the momentum source term uses
(`Source.f90`). In the layer the momentum equation is the hydrostatic pair
`dp/dr ~ -rho g`, so the scaled number is the fractional violation of
hydrostatic balance -- a controlled statement. The mass and energy rows keep
`|u|` (their denominators do not collapse), and the printed line now reports
the three rows separately:

    (JFNK) below r_esc [1:401] |R|: mass= ... mom/grav= ... energy= ...

The change is print-only; `make check` remains 5/5 byte-identical.

### 10.2 The section 9.4 states, re-measured

Each converged state was loaded (`Load IC? True`, its own converged output as
the IC) and re-entered the JFNK, which prints the diagnostic at solve start
and at the returned iterate. The restart marches 2000 steps before the
hand-off (the staged secondary-ionization hold), so the "start" state is the
loaded solution plus that hold, not byte-for-byte the section 9.4 state.

| state | old norm (sec 9.4) | mom/grav at JFNK start | mom/grav at end | re-verify |
|---|---|---|---|---|
| `examples/15` converged (metals off) | 9.99 | 1.94e-5 | 4.50e-5 | `info = 0`, `\|\|R\|\| = 6.105e-6` (cold re-run; reproduces sec 9.3 exactly, hand-off step 25839, `log10 Mdot = 10.24`) |
| `examples/16`, `eps4 = 5e-3` | 0.195 | 1.26e-5 | 4.21e-6 | `info = 2` at `\|\|R\|\| = 3.887e-5` (see below) |
| `examples/16`, `eps4 = 1e-2` | 2.352 | 2.19e-5 | 1.35e-4 | `info = 0`, `\|\|R\|\| = 7.29e-7`, `log10 Mdot = 9.67` |

On the physical scale the three converged solutions -- including the
metals-off one that read 9.99 before -- all hold the layer's hydrostatic
imbalance at **1e-5 to 1.4e-4 of `rho g`**, four to six orders below the old
numbers and comparable to the wind-region `Resid tol`. The mass and energy
rows stay at 1e-3 and 1e-2 to 5e-2 respectively; the worst cell moves from the
near-base momentum rows (where the old scale put it artificially) to the
energy row at `r = 1.52-1.57`.

So the statement the caveat of section 9.5 rests on can now be made
quantitatively: the layer is *not* momentum-unbalanced -- it satisfies
hydrostatic balance to ~1e-4 -- and what stays uncontrolled between the
`eps4` pair members is *which* hydrostatic stratification the layer settles
on (the rho/T profile degeneracy of section 9.5), not the balance itself.

Two side observations from the re-measurements, recorded as-is:

- The `eps4 = 5e-3` re-verify stalled at `info = 2`, `||R|| = 3.887e-5`
  (4x above tol) after the 2000-step hold, where the `1e-2` member re-verified
  to `info = 0` in 67 s. Consistent with the marginal character of that member
  in section 9.3; the layer values above are from the printed start/end lines
  and stand regardless.
- A warm restart of the *metals-off* `examples/15` converged state dies with
  a NaN at step ~165 (`T = NaN` in a cell with `nhi = 0`), after the staged
  secondary-ionization flip at step 2. The cold re-run of the same
  configuration converges cleanly, so this is a restart-path artifact
  (the state is integrated for one step without the coupling it was converged
  with); not investigated further here.
