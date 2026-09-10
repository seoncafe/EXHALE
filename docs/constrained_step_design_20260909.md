# A step that moves along the element-budget surface (N4b design memo for decision 14)

Advisor memo, 2026-09-09, for PLAN_20260909_rev1 item N4b. It states the two
routes, what each needs, and the advisor's recommendation. Nothing here is
implemented; N4a (the constraint map) is being built first and either route
reads it.

## 1. The problem the step has to solve

The stationary unknowns of a cell are the three hydrodynamic quantities,
the element mass fractions that are transported (helium and the carried
trace elements), and the carrier densities as `ln n` (B5k). The feasible
set of the carriers is not a box (B5l, review R4): with `nH_avail` the
hydrogen the carriers may hold,

```
2 n(H2) + n(OH) + 2 n(H2O) [+ n(H+)] <= nH_avail
n(OH) + n(H2O) + n(CO)               <= nO_free
n(CO)                                <= nC_free
```

and `nH_avail`, `nO_free`, `nC_free` are themselves functions of the
density and element unknowns of the same cell. The present code freezes
those budgets at the outer iterate and treats each carrier's share as a box
side (B5l). MEASURED consequences: the box removed every budget refusal on
`mol_carrier` but the solve still stalls on the Krylov leg; on
`mol_diffusion` 18 cycles are truncated because
`largest_step_inside_the_species_box` returns zero for a whole direction
when one component has no room.

## 2. Route (i): a linearized constraint set on the present unknowns

- Constraints `c_k(Y) <= 0` from the N4a map, with derivatives through the
  carrier logarithms (`dn/d ln n = n`) and through the density and element
  unknowns that enter the budgets.
- The step: minimize the same model `0.5||r + A s||^2` over the trust
  region intersected with the linearized set `c_k(Y) + grad c_k . s <= 0`.
  Practical form: the dogleg on the projected subspace with an active set
  of the nearly active constraints (the Bertsekas two-metric idea already
  used for the bounds), where a component blocked by a shared constraint is
  projected onto the constraint's hyperplane instead of being zeroed; the
  fraction-to-the-boundary rule applies to the projected step, so one
  blocked component no longer costs the whole direction.
- Nonlinear constraints re-evaluated at the candidate; a candidate that
  violates by more than a stated fraction of the budget enters a
  restoration phase (a step toward feasibility measured by the N4a
  violation magnitude) before the merit is judged.
- Certification unchanged: `F(Y) = 0` on the physical rows, feasibility by
  the map, and a nonzero balance row at a bound stays uncertified
  (decision 12 c).
- Cost: one derivative of each `c_k` (analytic from the map, cheap), the
  active-set logic in `fix_active_species_bounds`,
  `largest_step_inside_the_species_box` and the projection, all in
  `steady_newton.f90`. No change of unknown layout, no restart impact,
  no input change.

## 3. Route (ii): an element-conserving parameterization

- Replace the carrier densities by fractions of the element budget (for
  hydrogen: the shares of `nH_avail` held by H2, OH, H2O, H+ and the
  remainder; a simplex parameterized by `k - 1` unbounded coordinates, for
  example logits), so that every carrier state satisfies the shared
  constraint by construction and the box disappears.
- Needs: independent degrees of freedom across the H and O simplices (OH
  and H2O sit in both, so the two are coupled and one of them has to be
  the free coordinate of one simplex and derived in the other); zero and
  trace abundances (a logit of an exactly absent species is minus
  infinity; an explicit absent-species treatment is needed, as N1's
  removal of absent elements already shows); Jacobian conditioning of the
  logit map near the simplex corners; every reader of the flat unknown
  vector (write-back, restart of a stationary state, the report lines, the
  test suite) changes; the `ln n` floor and the B5k measurements become
  historical.
- Cost: a change of unknown layout touching `steady_newton.f90`,
  `diffusive_photochemistry.f90`, the suites and the design documents;
  restart files unaffected (they store species, not unknowns) but the
  stationary-state hand-off is re-derived.

## 4. Recommendation

Route (i) first. It answers the diagnosed defect (shared budgets not in the
step model; one blocked component costing a direction) without changing
what an unknown is, so B5j/B5k/B5l measurements stay comparable and the
restart and hand-off paths are untouched. Route (ii) is kept open on
evidence: if route (i)'s active set thrashes at the budget surface on the
carrier fixtures (measured as the fraction of iterations with a changing
active set), or if the H/O coupling of OH and H2O makes the projected
subspace degenerate, the parameterization is the better formulation and
is presented again with the numbers.

Acceptance for either route (PLAN rev 1, Stage B gate): a feasible
direction along a shared constraint is not discarded because a coordinate
bound was frozen; an infeasible start is restored by measured violation;
the B5l comparison repeated at one saved state by violation magnitude,
actual merit reduction, achieved linear residual and elapsed time.

## 5. What decision 14 decides

Only the route. The N4a map, the restoration phase and the certification
rule are common to both and are not part of the decision.
