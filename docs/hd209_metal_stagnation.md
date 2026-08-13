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
`info = 0` at `Resid tol: 1.0e-5` for the first time.

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
- **`examples/16` itself.** Section 3's inheritance is gone (section 5 item 3),
  and the case is no better for it. With `Low-Mach damping` at
  `eps4 = 5e-3`-`1e-2` the warm restart now reaches `info = 0` at
  `Resid tol: 1.0e-5` (section 8.2), which it never did before; the cold path
  and the larger coefficients still do not.
- **Cells that end below the solver tolerance.** The equilibrium still returns,
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
