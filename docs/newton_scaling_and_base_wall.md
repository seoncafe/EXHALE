# Newton scaling, the line-search merit, and the "base wall" (2026-08-10)

> **STALE (P35, 2026-09-02).** The `vulcan_work` run directories cited
> here carry `Molecular base: True` with the molecular network off. The base
> particle count was molecular while the species state was atomic, so the
> base ghost sat at `ntot_bc x T0` -- between 0.555 and 0.994 of the
> requested `T0`, depending on the directory -- and startup now refuses the
> combination. The numbers below are kept as recorded but stand on that
> base; see `INVALID_BASE_TEMPERATURE.md` in each run directory and item
> P35 of `TO_BE_DONE.md`.

## Summary

The steady-state JFNK solve had been floored around `||R|| ~ 2.8e-3` against a
target of `1e-3`, aborting with `info = 2`, on the HD 209458 b configuration
`vulcan_work/hd209_wind_response/photo_deep_secion_cont`. The standing
explanation, recorded in `TO_BE_DONE.md` item (A) and in
`docs/base_composition_handoff_plan.md` §11.8-11.9, was that the base cell
`j = 1` carries a momentum imbalance that no interior state can cancel — the
"base wall".

Direct measurement does not support that explanation. Two things were wrong
instead, both in the solver rather than in the boundary condition:

1. **The diagonal scaling `build_scaling` used a global-maximum floor.** The
   momentum floor `1e-6 * max_j |rho v|` was set by the base cell, which holds
   the global maximum of `|rho v|` while carrying no wind. Wind cells near the
   stagnation point fell onto that floor, and twelve of them supplied 69 % of
   the line-search merit.
2. **The stagnation watchdog judged a quantity the solver does not minimize.**
   It aborted after 15 outer iterations without a new best `||R||`, where
   `||R||` is the volume-weighted relative residual over `[j_min:N]` only. The
   opening pseudo-transient legitimately raises `||R||` for about 30
   iterations while it repairs the sub-sonic region, so the watchdog fired
   during normal operation, every time.

After replacing the scaling by a local, parameter-free one
(`cell_state_scales`) and replacing the watchdog by a count of consecutive
failed line searches, the target case converges: `info = 0`,
`||R|| = 5.053e-04`, 59 outer iterations.

Statements below marked *measured* are values produced by running the code;
the reading of what they mean is separated out and stated tentatively.

## 1. The base momentum row is solvable (the "base wall" hypothesis, refuted)

*Measured*, on the iterate JFNK is handed at marching step 2002 (WENO3,
`weno_mode = 0`, valve `eps = 1.0e-4`; grid `N = 500`, `Ng = 2`,
`j_min = 404`, `r(1) = 1.0001933`, `dr_j(1) = 1.9332e-4`):

| term of `R(2,1) = dFconv + dFpres - S2` | value |
|---|---|
| convective flux difference `dFconv` | -2.98290e+01 |
| pressure gradient `dFpres` | -9.90306e+01 |
| gravity source `S2` | -1.28851e+02 |
| residual `R(2,1)` | -8.36402e-03 |

so `|R(2,1)| / |S2(1)| = 6.49e-05`: the base momentum residual is the
remainder of a four-order cancellation of gravity against the pressure
gradient, not a large unbalanced force. Four repeat evaluations at fixed `Y`
(resetting `f_sp` each call, as the line search does) reproduced
`R(2,1) = -8.3640207e-03` bitwise, so there is no stochastic noise floor from
the iterative chemistry at the base.

Rescaling the ghost pressure by a factor `s` with the interior frozen drives
`R(2,1)` through zero at `s* = 0.99999708`, a **2.9 ppm** change (equivalently
`-0.0030 K` on a ghost temperature of 1027.6 K); rescaling the ghost density
instead nulls it at `s* = 1.0020942`, a 0.21 % change, and simultaneously
reduces `|R(2,2)|` to 0.40x. Ghost pressure is about 560x the lever ghost
density is.

*Reading (tentative).* The base row is not over-determined and not
unsatisfiable. A cell-1 momentum residual of `8.4e-3` is a 65 ppm imbalance of
the local hydrostatic terms and lies well within reach of the interior degrees
of freedom the Newton solve already has. The reason it kept appearing as the
"worst cell" was the diagnostic in §2, not its size.

## 2. Where the merit actually lived

The JFNK line search minimizes `||F/D||_2` over the whole domain, with `D` the
diagonal scaling. Convergence, on the other hand, is declared on
`resid_relnorm`: the volume-weighted relative residual, evaluated separately
over the wind `[j_min:N]` and the layer below the escape radius and combined by
the larger of the two (§127 of `Update_EXHALE.md`; at the time of this memo it
was the wind window alone, here `r >= 1.992`). The scale each row is divided by
is `residual_row_scale`, which since §133 is a bound on that row's own largest
term and since §143 is that term itself, taken from what the row contains. At
the time of this memo it was `max(|Y_i|, 1e-6 max_j |Y(:,k)|)`, which is what
the numbers below are measured against.

*Measured* at the hand-off state, with the old `build_scaling`
(`D_i = max(|Y_i|, 1e-6 max_j |Y(:,k)|)`):

| region | share of `\|\|F/D\|\|_2^2` |
|---|---|
| `j = 1..10` (base) | 0.499 % |
| `j = 11..403` | 99.501 % |
| `j = 404..500` (the convergence window) | 0.000 % |
| momentum rows, all `j` | 99.988 % |

and the twelve largest contributions were all momentum rows at
`j = 307..318`, `r = 1.322..1.366`, together 69 % of the total. In that band
`u(2,j)` rises monotonically from `-2.34e-08` to `+1.29e-08`, crossing zero
between `j = 313` and `j = 314` (`r ~ 1.347`): at the hand-off the whole region
below `r ~ 1.35` is still infalling, so `|rho v|` passes through zero there.
The momentum floor `1e-6 * max_j |rho v| = 7.0446e-09` — set by the base cell,
where `|rho v| = 7.04e-03`, six orders above the wind — is comparable to
`|rho v|` itself throughout that band, so those cells were divided by a number
unrelated to their own state.

The first Newton step raised `||R||` by 25x while lowering the merit by 13x.
Two functionals, two different regions.

Also *measured*: the solver's own "worst cell" print normalized by
`max_j |u(k,j)|`, which for momentum is again the base value. That diagnostic
reported `j = 1` or `j = 2` on nearly every iteration and is the direct source
of the base-wall reading in `TO_BE_DONE.md`.

## 3. The new scaling

`build_scaling` is replaced by `cell_state_scales`, which builds each scale
from the cell's own state:

```
mass       D = rho_j
momentum   D = |rho v|_j + rho_j c_s,j  =  rho_j (|v|_j + c_s,j)
energy     D = E_j
```

with `c_s = sqrt(g p / rho)` and `p = (g-1)(E - (rho v)^2 / 2 rho)`, so `D` is
a function of the unknowns alone.

`rho` and `E` are positive definite and are their own scales. The momentum
density is not: it vanishes wherever the flow reverses, which at the hand-off
is the entire sub-sonic region. The quantity that does not vanish there is
`rho` times the fastest characteristic speed of the Euler system, `|v| + c_s`
— the momentum density the cell carries when moved at its own signal speed.
That choice is parameter-free (there is no floor fraction to pick) and local
(no cell can set the scale of another).

**Still current after §143, with one boundary drawn.** These three numbers are
the scaling of the Newton system and of the line-search merit, and §143 does
not touch them. What it changes is the CONVERGENCE MEASURE: each residual row
is now divided by the largest term that row itself contains — the face mass
flux for the mass row, `max(|dF_2|,|S_2|)` for the momentum row,
`max(|dF_3|,|S_3|,heat,cool)` for the energy row (`mass_flux_row_scale`,
`momentum_row_scale`, `energy_row_scale` in `steady_residual.f90`) — so the
merit and the acceptance test are no longer one expression apart on any row.
That separation is deliberate and it was measured: a build that rescaled this
system by those quantities as well, the mass row's varying by a factor 9 across
the first two cells, left the molecular hot-Uranus solve with no descent
direction after 179 iterations and 8238 residual evaluations, against 14
iterations and 321 for the build that changes the measure alone
(`docs/p54_base_layer_mass_flux.md` §10.4, `docs/p54g23_row_scale_scan.md`).

One measurement in §1 is worth re-reading beside it. This memo records that
cell 1's momentum row is a four-order cancellation whose remainder is 6.5e-5 of
the gravity term, and that a 2.9 ppm ghost-pressure change drives it to zero.
§143's momentum scale measures exactly that remainder against exactly that
gravity term, and on the molecular hot Uranus cell 1 still reads 0.14 on it
after a converged solve — so the "one degree of freedom away" reading holds for
HD 209458 b and does not carry to that planet.

The intermediate forms `|rho v| + f rho c_s` with `f < 1` were tried and are
reported in §5; `f = 1` is both the principled choice and the best-performing
one measured.

*Measured* merit distribution at the same hand-off state after the change:

| region | old scaling | new scaling |
|---|---|---|
| `j = 1..10` | 0.499 % | 0.013 % |
| `j = 11..403` | 99.501 % | 99.986 % |
| `j = 404..500` | 0.000 % | 0.001 % |
| momentum rows | 99.988 % | 15.348 % |
| largest single cell | 10.9 % (`j = 312`, momentum) | 1.17 % (`j = 275`, mass) |

`||F/D||_2` itself falls from `5.398e+02` to `6.731e+00`. The stagnation-point
concentration is gone: no cell now carries more than 1.2 % of the merit, and
the largest contributions sit at `r ~ 1.21-1.22` in the mass and energy rows.

## 4. The watchdog, and why the merit was left alone

With the new scaling the target case still aborted at `info = 2` after 14
iterations. Raising only the abort limit (nothing else) was enough to make even
the **old** scaling converge, which located the second fault.

*Measured*, target case, with the abort limit raised to 100 and the merit
unchanged:

| momentum scale | info | final `\|\|R\|\|` | outer iterations | log10 Mdot |
|---|---|---|---|---|
| old, `1e-6 max_j \|rho v\|` | 0 | 9.334e-04 | 58 | 9.48 |
| `rho(\|v\| + c_s)` | **0** | **5.053e-04** | 59 | 9.47 |
| `\|rho v\| + 0.1 rho c_s` | 2 | 1.817e-03 | 130 | 9.52 |
| `\|rho v\| + 0.01 rho c_s` | 2 | 2.802e-03 | 99 | 9.70 |
| `\|rho v\| + 0.001 rho c_s` | 0 | 7.841e-04 | 91 | 9.48 |

and the iteration counts that matter for the watchdog:

| momentum scale | first iteration to beat the start `\|\|R\|\|` | longest run with no new best `\|\|R\|\|` | longest run of failed line searches |
|---|---|---|---|
| old | 27 | 27 | 3 |
| `rho(\|v\| + c_s)` | 31 | 30 | 4 |
| `\|rho v\| + 0.001 rho c_s` | 38 | 37 | 4 |
| `\|rho v\| + 0.1 rho c_s` (stuck) | 31 | 99 | 99 |
| `\|rho v\| + 0.01 rho c_s` (stuck) | never | 99 | 34 |

Every variant needs 27-38 iterations before it first improves on the
warm-start `||R||`, so a 15-iteration patience on "no new best `||R||`" aborts
during the opening pseudo-transient in all of them. By contrast, the count of
**consecutive iterations in which the line search found no acceptable step at
all** separates the two populations cleanly: at most 4 in every run that goes
on to converge, at least 34 in every run that is genuinely stuck. That is also
the exact failure mode the original code comment describes ("a non-monotone
search may accept near-zero steps forever").

The watchdog is therefore now: abort after `n_no_descent_max = 12` consecutive
iterations with no accepted step. Best-iterate tracking on `||R||` is kept, so
a failed solve still returns its best state. The convergence test itself
(`||R|| < resid_tol` on `[j_min:N]`) is unchanged.

### The volume-weighted merit was measured and rejected

Making the line-search merit volume-weighted, to match `resid_relnorm`, was
tried in two forms — `sqrt(sum_j w_j (F/D)^2 / sum_j w_j)` with `w_j = r_j^2
dr_j`, over the whole domain and over `[j_min:N]` only. *Measured*, all with
the abort limit at 100-200:

| scaling | merit | info | best `\|\|R\|\|` | outer iterations |
|---|---|---|---|---|
| `rho(\|v\|+c_s)` | plain 2-norm | **0** | **5.053e-04** | 59 |
| `rho(\|v\|+c_s)` | volume-weighted, `1..N` | 2 | 2.660e-03 | 220 |
| `rho(\|v\|+c_s)` | volume-weighted, `[j_min:N]` | 2 | 1.472e-03 | 100 |
| `\|rho v\| + 0.3 rho c_s` | volume-weighted, `1..N` | 2 | 2.802e-03 | 99 |

Volume weighting also barely changes the alignment it was meant to fix: it
lifts the share of the merit carried by `[j_min:N]` from 0.001 % to 0.021 %,
because the near-base cells are small in volume but their residual *rates*
`F/D` are correspondingly large — the two effects nearly cancel. Both weighted
variants stall in long runs of failed line searches. The change was therefore
not adopted; the merit remains `||F/D||_2` over the whole domain.

*Reading (tentative).* Aligning the merit's *region* with the convergence
test's region appears to be the wrong direction: the Newton step is computed
from the full residual over all rows, and a merit that ignores most of those
rows rejects the steps that step actually takes. What needed aligning was the
*scale* of each cell, which is what §3 does, and the watchdog, which is what this
section does.

## 5. Sensitivity to the momentum floor

*Measured* (table in §4): `f = 1` and `f = 0.001` converge, `f = 0.1` and
`f = 0.01` do not, and the old global-maximum floor converges. The dependence
on `f` is not monotone, and no mechanism for the non-monotonicity has been
identified here. This is a caveat on how much margin the target case has, not
an argument about which scale is right: `f = 1` is the only value that is a
physical statement (`rho` times the fastest signal speed) rather than a
fraction of one, and it also gives the lowest converged residual and the
fewest iterations. It is what the code now uses, with no key to change it.

## 6. Hand-off timing (no change made)

The hand-off was suspected of firing while the stagnation point is still
inside the domain (`r ~ 1.35`). *Measured*: setting `Solver: Newton 1.0e-2` in
place of `5.0e-2` changes nothing — the log confirms the new threshold is read
(`JFNK hand-off at du < 1.00E-02`) and the Newton finish still starts at step
2002, with an identical iteration trace and identical results. The reason is
in `EXHALE_main.f90`: `du` is already `1.93e-03` at step 2, below both
thresholds, and the hand-off point is set by the staged secondary-ionization
hold (`N_stall = 2000` steps after the flip at step 2), not by
`newton_du_switch`. Lowering `du_switch` cannot move the hand-off in this
configuration, so no code or input change was made for it.

## 7. Gates

| gate | result |
|---|---|
| G1 `make check` (wasp_full, wasp_he23off, mol_base_handoff) | **PASS** — `REGRESSION PASS (all cases byte-identical)`; `wasp_full` 13488 steps, `wasp_he23off` 13482, `mol_base_handoff` 12000, each matching its golden bitwise |
| G2 target case `photo_deep_secion_cont` | **PASS** — `info = 0`, `\|\|R\|\| = 5.053e-04` (target `1e-3`), 59 iterations, 63 s wall; baseline was `info = 2`, `\|\|R\|\| = 2.802e-03`, no improvement in 14 iterations, 72 s |
| G3 `wasp_full_newton`, `newton_rsw01` | **PASS, improved** — both now `info = 0` where both previously ended `info = 2` |

G3 in detail, both runs made from the case's `input.inp` in a scratch directory
(the `backup/regression` copies were not touched), the "before" column produced
by a binary built from the same tree with only `steady_newton.f90` reverted to
`HEAD`, so the marching path is identical on both sides:

| case | before | after |
|---|---|---|
| `wasp_full_newton` | `info = 2`, `\|\|R\|\| = 1.201e+00`, 14 iterations, marching resumed to step 6314, `log10 Mdot = 13.20` | `info = 0`, `\|\|R\|\| = 5.643e-04`, 17 iterations, stops at the hand-off step 4247, `log10 Mdot = 13.17` |
| `newton_rsw01` | `info = 2`, `\|\|R\|\| = 9.578e-01`, 20 iterations, marching resumed to step 8843, `log10 Mdot = 13.16` | `info = 0`, `\|\|R\|\| = 9.726e-04`, 14 iterations, stops at the hand-off step 4041, `log10 Mdot = 13.17` |

Both "before" numbers are far worse than the `1.260e-01` / `9.065e-02` in those
directories' own `run.log` files, because those logs predate the intervening
changes to the marching path (the `newton_rsw01` log still prints
`(ATES_main)`). That is why the comparison was made against a freshly built
pre-change binary instead.

`log10 Mdot` for the target case moves from 9.70 to 9.47. The 9.70 came from
the marching loop resuming after the JFNK failure and stopping on `du`; 9.47 is
the first Newton-converged value for this configuration. Per the standing rule
that a `du`-threshold stop has a path-dependent spread of several percent, 9.47
is the number to quote.

## 8. What was not checked

- The three default regression cases never enter JFNK (`wasp_full` and
  `wasp_he23off` have no `Solver:` line; `mol_base_handoff` has
  `Solver: Newton` but is pinned to 12000 steps and ends at `du = 2.67`), so
  G1 tests that the change did not leak outside the solver, not the solver
  itself.
- `wasp_full_newton` and `newton_rsw01` have **no goldens** in
  `backup/regression/golden/` (only `wasp_full`, `wasp_he23off` and
  `mol_base_handoff` do), and they are not in `DEFAULT_CASES`. Nothing was
  re-snapshotted, and their own directories were not written to.
- The new watchdog catches the failure mode that was observed (the line search
  exhausting all 20 backtracks, `lam = 9.54e-07`, repeatedly). It would not
  catch a solve that keeps *accepting* negligible steps against the
  non-monotone reference, which the old rnorm watchdog would have stopped after
  15 iterations. That mode did not occur in any run measured here; if it does,
  the solve now runs to `maxit` (500 in the production call, a couple of
  minutes for these grids) and still returns its best iterate.
- Only HD 209458 b and WASP-121 b configurations were run. The scaling change
  touches every JFNK and PTC solve, so any other Newton-finished run will
  differ in its iteration path; whether it differs in its converged state has
  not been checked outside these three cases.
- `solve_steady_ptc` uses the same `cell_state_scales` through
  `build_banded_jac_full` and was not exercised (no case in the tree selects
  the PTC driver by default).

## 9. `du` trigger arming (2026-08-11)

A separate change to `src/EXHALE_main.f90`, made after the solver work above,
and recorded here because it was found while re-converging WASP-121 b with the
new solver.

### What was measured

The WASP-121 b re-convergence (`WASP-121b/input.inp`: `Transonic IC: True`,
`Load IC? False`, `Solver: Newton`, `du_th [PLM,WENO3]: 0.5 1.0e-3`) died in
NaN. `run_20260810.log` shows the cause on step 1:

```
    -> switched PLM -> WENO3 at du =  1.8531E-05, step 1
    -> secondary ionization activated ahead of the Newton finish at step 1, du =  1.85E-05
```

`du` is the radial spread of `rho v r^2` over `[j_min:N]`. On step 1 it was
`1.85e-05` — below every threshold in the run — because the transonic IC is an
analytic profile that has not yet been touched by the marching loop. From there
`du` rose monotonically (`4.01e-04` at step 2, `1.93e-02` at 50, `4.09e-02` at
100, `0.292` at 500, `1.96` at 1000, `21.5` at 1220) and the run hit
`NaN detected at r = 1.000196` at step 1220.

The interpretation, stated tentatively: a small `du` on an unrelaxed generated
IC measures the smoothness of the formula that built it, not relaxation of the
wind, so it is not evidence that the state satisfies anything. Switching the
secondary-ionization coupling on there is the situation the staged activation of
§38 of `Update_EXHALE.md` exists to avoid — the coupling was applied to a state
that had never been relaxed at all, which is a stronger version of the cold-IC
runaway that motivated the staging.

Cold hydrostatic starts do not have this property. *Measured* `du(1)`:
`0.37029` for `wasp_full` and `wasp_he23off`, `7.89` for `mol_base_handoff` —
all above every threshold those cases use.

### The change

Two triggers now require `du` to have been seen at or above their own threshold
before they can fire, so that only a descending crossing counts:

- the `du < du_th` convergence stop (`du_stop_armed`), and
- the JFNK hand-off, together with the secondary-ionization flip that shares its
  `du` comparison (`du_newton_armed`).

Both flags start `.true.` when the run began from `Load IC? True`, since that
state was relaxed by the marching loop of the run that wrote it and its `du` is
meaningful, and `.false.` for every generated IC (cold hydrostatic, transonic,
Wind-AE, Parker). Arming is reported in the log (`-> du stop armed at du = ...`,
`-> Newton hand-off armed at du = ...`).

The PLM -> WENO3 switch is deliberately **not** guarded: the regression cases
cross their `0.5` stage-1 threshold on step 1 (`du(1) = 0.37`), so their goldens
already contain a step-1 switch, and guarding it would change them. Whether the
reconstruction switch should also require a descending crossing is left open.

### Results

| gate | result |
|---|---|
| `make check` (wasp_full, wasp_he23off, mol_base_handoff) | see below |
| warm start `photo_deep_secion_cont`, `EXHALE_MAXSTEPS=3000`, 8 threads | **unchanged** — flip at step 2, JFNK at step 2002, `info = 0`, `\|\|R\|\| = 5.053e-04`, identical to §7 |
| WASP-121 b, full run, 16 threads | **no NaN**; `du` stop armed at step 4 (`du = 1.147e-03`), Newton hand-off armed at step 28 (`du = 1.038e-02`), PLM -> WENO3 still on step 1 |

The WASP-121 b run (`WASP-121b/run_20260811.log`) reaches `du = 9.94e-03` at
step 6228, flips the secondary-ionization coupling there, holds `N_stall = 2000`
steps, and engages JFNK at step 8228 with `||R|| = 4.663e-01`. That solve
aborts (`info = 2`, best `||R|| = 2.575e-01`, no descent step in 12 consecutive
iterations), marching resumes as designed, and the run converges on `du` at step
8706 (`du = 9.9285e-04`, `dtu = 1.2104e-05`; flux spread `2.299e-03` at the last
diagnostic), giving `log10 Mdot = 13.21`. No NaN appears in any output file.

So the guard removes the NaN mode, but WASP-121 b is *not* Newton-finished: its
quoted `Mdot` still comes from a `du` stop, which carries the usual
path-dependent spread. Why the JFNK solve stagnates on this configuration
(worst cell `j ~ 357-366`, `r ~ 1.22`, energy row) is a separate open question
from the base-cell scaling treated in §1-§5.

### Addendum 2026-08-15: measured on a Wind-AE IC

The arming table above covered the transonic IC and the cold hydrostatic
starts; `IC mode: windae` -- the case whose false `du` stop within ~3 steps
originally motivated the "always Newton-finish" prescription -- had not been
measured. Run: the `examples/12_windae_ic_hd209` configuration (HD 209458 b,
seed-adjacent, `Solver: Newton`, default `du_th = 1e-3`), HEAD binary,
`EXHALE_MAXSTEPS = 2000`, 4 threads.

- `du` at step 2 is `1.4267e-04` -- below both thresholds on the untouched
  generated IC, exactly the smooth-formula situation described above.
- Neither trigger fires there: `-> du stop armed at du = 1.0270E-03, step 23`
  and `-> Newton hand-off armed at du = 1.0018E-02, step 222`, both on the
  first ascending crossing, and the run marches on to the step cap
  (`du = 0.185` at 2000, still relaxing; the cap was the point of the run).

So the descending-crossing guard behaves on a Wind-AE IC as designed, and the
pre-arming false stop is structurally gone for this IC mode too. What the
guard does *not* change: an `Mdot` from a `du` stop still carries the
path-dependent spread, so the Newton finish remains the prescription for a
quantitative number.

## 10. Which residual decides the line search (2026-08-11)

The judgment first: the acceptance test in `solve_steady_jfnk` was wrong, and
wrong independently of how large its effect turns out to be. A line search must
decide on the residual the solve is driving to zero. This one decided on a
different quantity -- the residual with the WENO3 smoothness weights held at
the iterate the step starts from -- compared it against a reference built partly from the other
measure, and then reported `||R||` and the `info = 0` convergence verdict from
the same frozen-weights residual. What follows is measured on that code and on
its replacement.

### What was measured

Configuration: HD 189733 b warm start (`Load IC? True`, `Solver: Newton 2.0e-2`,
`Shapiro filter: 0.02 4`), run in a scratch copy with `EXHALE_MAXSTEPS = 2100`
and 8 threads, so the marching loop hands off to JFNK at step 2002 and stops
shortly after. This is the configuration that ends in the abort of the §4
watchdog. A measurement build was made from the same tree that evaluates every
line-search trial **both** ways and prints the pair: weights frozen at the
current iterate (`weno_mode = 2`, what the old test used) and weights recomputed
at the trial state (`weno_mode = 0`, the residual being driven to zero). They are
called the frozen and the true merit below.

*Measured*, 245 trials over the 26 outer iterations the old code takes before it
aborts: the ratio true/frozen has median 1.000 and range 0.892 to 4.307. The
median is set by the deep backtracks, where both measures collapse onto the
merit of the current iterate. Restricted to the full steps `lam = 1`, which are
the ones that move the solve, the ratio has median 2.17 and maximum 4.31 over
10 trials.

The consequence, *measured* on the same run: of the 14 steps the frozen test
accepted, 7 raised the true merit.

| outer it | `lam` | true merit at the iterate | frozen merit of the step | true merit of the step | true merit x |
|---|---|---|---|---|---|
| 2 | 0.125 | 17.65 | 62.40 | 159.8 | 9.05 |
| 4 | 0.5 | 61.12 | 167.3 | 149.3 | 2.44 |
| 8 | 1.0 | 25.89 | 15.03 | 50.42 | 1.95 |
| 9 | 0.5 | 50.41 | 148.5 | 311.9 | 6.19 |
| 11 | 0.5 | 109.8 | 83.72 | 120.1 | 1.09 |
| 13 | 1.0 | 97.69 | 27.57 | 118.7 | 1.22 |
| 14 | 0.5 | 118.8 | 80.22 | 181.3 | 1.53 |

Iteration 8 is the clearest single case: the frozen measure reports the step
cutting the merit by 42 %, and the true residual of that same step is 1.95x the
residual it started from. Iteration 14 is the one that ends the sequence -- it
leaves the iterate at a true merit of 181.3 (205.85 once the scaling is
refreshed, see below) while writing 80.22 into the merit memory, and the search
never recovers from there.

Two further bookkeeping errors were measured on the same run.

**The reference went stale.** The non-monotone reference is
`f2ref = maxval(f2hist)` over a 5-slot memory that was written *only when a step
was accepted*, and written with the accepted step's frozen merit. At iteration
15 the memory held `f2ref = 84.15`, a frozen value recorded five iterations
earlier, while the true merit of the iterate the search was starting from was
205.85. The test therefore demanded a factor-2.4 reduction in one step; all 20
backtracks fail, `lam` bottoms out at `9.54e-07`, and after 12 such iterations
the §4 watchdog aborts the solve. In this configuration the deadlock is a
property of the reference, not of the frozen/true gap: at the smallest `lam` the
two measures agree to all 5 printed digits (205.85 both).

**The diagonal scaling was one iteration behind the merit it scaled.**
`cell_state_scales` was called *after* the merit of the current iterate was
formed, so that merit used the previous iterate's `D` while every trial used the
current one. At iteration 15 the iterate's merit prints as 181.28 while the
`lam -> 0` limit of the trials is 205.85, a 12 % inconsistency in what is
supposed to be the same number.

**There is no noise floor to justify the window.** The comment above the search
said a strictly monotone test rejects valid steps once the required decrease
falls below "the iterative-chemistry noise floor of the residual evaluation".
§1 measured that floor: four repeat evaluations at fixed `Y`, resetting `f_sp`
as the line search does, reproduced the residual bitwise. The comment has been
corrected; the window is kept for the reason measured below, not for that one.

### The change

In `src/modules/time_step/steady_newton.f90`, `solve_steady_jfnk`:

1. Line-search trials are evaluated with `weno_mode = 0`, so acceptance is
   decided on the true residual. Mode 0 recomputes the smoothness weights
   without storing them, so the frozen model of the current outer iteration
   survives the search untouched.
2. What is stored on acceptance is that same true residual: `F`, hence `rnorm`,
   hence the `rnorm < resid_tol` verdict, are now the true residual as well.
3. The merit memory holds true merits only, and is written at the *start of
   every outer iteration* rather than only on acceptance, so the reference
   cannot fall behind the state the search starts from.
4. `cell_state_scales` is called before the merit of the iterate is formed.

What is deliberately unchanged: the frozen weights still define the inner Newton
model (the banded Jacobian and the `J*v` products), which is the lagged-weights
remedy the solver relies on; the watchdog, the best-iterate return, and the positivity
test are as they were.

### One residual evaluation per trial, not two

A two-stage form -- backtrack on the cheap frozen merit, then confirm the
candidate on the true one -- was considered and rejected on measurement. The
residual evaluation is not a pure function of `(Y, f_sp)`: `excited_H_update`
carries its populations and `heat_balmer` between calls, so an extra evaluation
per trial changes the sequence of iterates. That is visible directly: the
measurement build described above is exactly a two-evaluation line search, and
on this case it takes a different path from iteration 2 onward and ends
`info = 2`, `||R|| = 2.374e-03`, where the production build reaches `info = 0`,
`||R|| = 8.803e-04`. Reading, stated tentatively: a filter on the frozen merit
can only reject steps the true residual would have accepted, and it costs an
evaluation whose side effect is itself a perturbation, so the single true
evaluation is both the cheaper and the better-defined test. The state carried
between residual evaluations is noted here as an observation; it was not
changed.

### Keeping the non-monotone window (measured)

Two forms of the reference were built and run, because with the staleness
removed the window might have been unnecessary:

- (i) monotone Armijo, `f2ref` = the true merit of the current iterate;
- (ii) the 5-slot Grippo window, resynchronized at every outer iteration as
  described above.

| case | (i) monotone | (ii) resynchronized window |
|---|---|---|
| HD 189733 b, `Shapiro 0.02 4` | `info = 2`, `\|\|R\|\| = 9.134e-03` | `info = 0`, `\|\|R\|\| = 8.803e-04` |
| HD 189733 b, warm start | `info = 2`, `\|\|R\|\| = 1.537e-03` | `info = 0`, `\|\|R\|\| = 6.524e-04` |
| `photo_deep_secion_cont` | `info = 0`, `\|\|R\|\| = 9.730e-04` | `info = 0`, `\|\|R\|\| = 8.574e-04` |

The monotone form is a regression on both HD 189733 b cases, the second of which
converges even with the old code. Form (ii) is what is in the tree. The solver
now prints how many accepted steps needed the window (`window-only accepts=` on
the closing line): 10 of 28 on the first case, 1 of 2 on the second, 8 of 24 on
the third, 0 of 22 on WASP-121 b.

### Gates

All runs in scratch copies of the case directories; no planet directory and no
golden was written to. "before" is a binary built from this tree with only
`steady_newton.f90` at its pre-change state, so the marching path is identical
on both sides. `OMP_NUM_THREADS = 8` except G1.

| gate | before | after |
|---|---|---|
| G1 `make check` (`wasp_full`, `wasp_he23off`, `mol_base_handoff`) | -- | **PASS**, byte-identical (these cases never enter JFNK) |
| G2 HD 189733 b, `Shapiro 0.02 4`, `EXHALE_MAXSTEPS = 2100` | `info = 2`, best `\|\|R\|\| = 1.756e-03`, 26 iterations, 12 failed searches, 58.7 s | **`info = 0`**, `\|\|R\|\| = 8.803e-04`, 28 iterations, 0 failed searches, 57.3 s |
| G3 HD 189733 b warm start, `Solver: Newton 2.0e-2` | `info = 0`, `\|\|R\|\| = 6.010e-04`, 9 iterations, 3 failed searches, 54.2 s | `info = 0`, `\|\|R\|\| = 6.524e-04`, **2 iterations**, 0 failed searches, 53.5 s |
| G4 WASP-121 b, full run from the transonic IC | `info = 2`, best `\|\|R\|\| = 2.575e-01`; marching resumes and stops on `du` at step 8706, `log10 Mdot = 13.21`, 194.2 s | **`info = 0`**, `\|\|R\|\| = 5.718e-05` at the hand-off step 8228, `log10 Mdot = 13.17`, 187.8 s |
| G5 `photo_deep_secion_cont`, `EXHALE_MAXSTEPS = 3000` | `info = 0`, `\|\|R\|\| = 5.053e-04`, 59 iterations, 31 failed searches, 70.4 s | `info = 0`, `\|\|R\|\| = 8.574e-04`, **24 iterations**, 0 failed searches, 57.2 s |

Across all four solver cases the search no longer exhausts its 20 backtracks
even once (`lam` never reaches `9.54e-07`), where before it did so 12, 3, 12 and
31 times.

G4 is the result that changes what can be quoted: WASP-121 b is Newton-finished
for the first time, so its `log10 Mdot = 13.17` no longer carries the
path-dependent spread of a `du` stop. It supersedes the 13.21 of §9 and answers
the open question that section ended on. G3 also shows the state is now a
genuine fixed point: the He/H diffusion outer pass that follows the solve moves
the composition by `3.58e-03` instead of `1.11e-01`, and the next JFNK entry
finds `||R||` unchanged at `6.524e-04` instead of jumping to `9.942e-04`.

`log10 Mdot` moves 9.21 -> 9.21 (G2), 9.21 -> 9.22 (G3), 13.21 -> 13.17 (G4),
9.47 -> 9.48 (G5).

### What was not checked

- The `||R||` in the "before" column is the frozen-weights residual and the one
  in the "after" column is the true residual, so the two numbers are not the
  same quantity. Only "did the run reach the target" compares directly.
- The three default regression cases never enter JFNK, so G1 tests that nothing
  leaked outside the solver, not the solver.
- No goldens were re-snapshotted and no golden case was affected.
- `solve_steady_ptc` shares `cell_state_scales` but has its own line search,
  which was not touched and is not exercised by any case in the tree.
- Only HD 189733 b, HD 209458 b and WASP-121 b configurations were run. Every
  Newton-finished run in the tree will take a different iteration path than
  before; whether any of them lands on a different converged state was not
  checked outside these four.
- The state carried between residual evaluations (§ above) was measured only to
  the extent of showing that an extra evaluation per trial changes the path. It
  was not traced to a specific variable beyond `excited_H_update`, and nothing
  was changed about it.
