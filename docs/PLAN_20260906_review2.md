# Second review: PLAN_20260906_rev1.md

Date: September 6, 2026, KST.

## 1. Verdict

**The revised plan is substantially improved and is conditionally acceptable.**
Most findings in the first review are now addressed at the planning level.
Retain the revised physical decisions and the A/B/C/D organization. A further
general rewrite is unnecessary; four focused amendments are needed before the
affected interfaces and validation routes are implemented:

1. Distinguish nonlinear trial admissibility, physical time-step acceptance,
   and final stationary certification. They may share an evaluator, but cannot
   share all acceptance thresholds.
2. Apply common energy and source-consistency prerequisites to pure H/He models
   as well as models with molecules or metals.
3. Make A0 an implemented and tested run-mode transition, explicitly accounting
   for the existing local pseudo-time stepping option.
4. Correct the dependency graph so that validity reporting does not depend on
   optional direct Newton support and stationary certification joins the source
   and spatial-operator implementations it uses.

This verdict concerns the **plan**, not completion of its implementation. The
core acceptance and source routines inspected here retain the behavior discussed
in the previous review; their unchanged hashes are recorded below. That is not
a new criticism of work that the plan explicitly schedules for later.

## 2. Scope, evidence, and artifacts

The full 333-line [rev1 plan](PLAN_20260906_rev1.md) was read and compared with
[the first review](PLAN_20260906_review.md) and the applicable requirements in
[development plan rev3](development_plan_20260905_rev3.md).

Relevant production paths were inspected again: main-program chemistry and
energy updates, local/global time-step construction and its caller, steady
Newton trial and stopping conditions, existing ionization sweep-state tags,
and the corrected SED coverage documentation/script.

The evidence categories are:

- **Plan-level resolution:** the requirement is now written adequately; it has
  not necessarily been implemented.
- **Source-confirmed:** the cited executable statements or call path were read.
- **Mathematical check:** a directly executed small analytical example, not a
  production Fortran test or an atmospheric simulation.
- **Unverified:** implementation performance, atmospheric impact, or a result
  requiring a simulation that was not run here.

No production code, existing plan, configuration, SED, regression reference, or
scientific product was changed. No build or atmospheric regression was run.
The historical binary identifier and regression counts in the plan were not
remeasured. This scope is appropriate because the requested artifact is a plan
review, the key inspected source files are unchanged, and the remaining findings
concern acceptance contracts and dependencies.

New mathematical checks are preserved in
[plan_review2_checks_20260906.py](audit_20260905/plan_review2_checks_20260906.py).
They require only the Python standard library and write no files. The earlier
SED and residual examples remain archived; they were not rerun simply to repeat
already established measurements.

## 3. Disposition of the first review

| Earlier finding | Assessment of rev1 |
|---|---|
| R1: premature scientific validation | Substantially addressed by separate diagnostic, physical, and comparison stages. The pure-H/He prerequisite exception still needs correction; see F2. |
| R2: thermal rejection before recovery | Addressed. B2 explicitly stops safely until B3a supplies recovery. No full recovery is claimed for B2 alone. |
| R3: class-4 transition | The intended modes and handoff are now clear. An implementation milestone and the existing local-time option still need to be connected to that contract; see F3. |
| R4: incomplete equation inventory | Addressed in the plan. B1a/A2 include oxygen carriers, elemental transport, independent variables, and the distinction between transport balance and local equilibrium. |
| R5: unevaluated zero residual and hidden state changes | Addressed in the plan by explicit evaluation status and an isolated or side-effect-free evaluator. The production change remains pending. |
| R6: raw solver status versus returned state | Addressed in the plan, including the relative-drop branch and post-limiter state. The new test list targets the relevant paths. |
| R7: aggregate norm hides smaller balances | Addressed in the plan by an independent condition for each balance and a local safeguard. Numerical units and floors appropriately remain part of the implementation briefs. |
| R8: unscheduled rev3 physics | Substantially addressed through B3b, B4, and B5. Their final dependency joins should be made explicit; see F4. |
| R9: byte identity as acceptance | Addressed. Identity is now a scoped expectation, not the criterion for physical correctness. |
| R10: nonlinear convergence versus time accuracy | Addressed in the test requirements. Those physical-time tests must run under a compatible time-stepping mode; see F3. |
| R11: coupled solve does not guarantee consistent outputs | Addressed through B3c's final-state assembly contract. Implementation must distinguish instantaneous output rates from integrated step budgets. |
| R12: SED sensitivity, thresholds, and coverage | Addressed as a plan requirement. The printed triplet-band label and the sampled current-threshold statements were corrected; configuration-specific coverage testing remains scheduled in C1. |
| Carrier checkpoint before retry | Addressed explicitly. The plan recognizes the existing cumulative CO record and requires isolation of rejected trials. |

These items should not all be reopened as unresolved defects merely because
their production implementation has not yet happened. The relevant distinction
is between an adequate work specification and a completed scientific capability.

## 4. Four focused amendments

### F1. A final stationary gate must not become the nonlinear trial gate

**Priority: high. Plan lines 97-124, especially 117-118.**

A2 says that both the marching carrier condition and the trial count of
uncertified cells are replaced by its evaluation. Sharing evaluated quantities
is reasonable. Requiring every Newton trial to pass final stationary residual
tolerances would be incorrect: the purpose of a nonlinear iteration is to move
through admissible states that do not yet solve the global stationary equations.

This wording is ambiguous rather than proof that the author intends the wrong
algorithm. Resolve it before the interface is coded.

**Source-confirmed distinction already present in the caller:**

- `steady_newton.f90:630-650` classifies a trial's local chemical closure and
  physical admissibility.
- Around lines 2293-2330, the line search calls `eval_residual`, then compares
  a trial merit with the current/reference merit.
- Around lines 2090-2103, `steady_gates_met` decides whether the current
  nonlinear iterate is a completed stationary solution.

The old chemical count is inadequate and class 5 is misclassified, as the plan
correctly recognizes. The replacement must not collapse the distinct purposes
of these checks.

**Recommended contract:**

| Context | Required acceptance or validity condition | Not required |
|---|---|---|
| Residual/Jacobian probe | A usable residual of the requested equation set, appropriate domain checks, controlled inner-closure error, and no adoption or persistent mutation | Global stationary convergence of the probe |
| Stationary Newton trial | Physical bounds and required inner closures, followed by the specified merit-decrease or trust-region acceptance condition | A zero global hydro/species/energy residual on every trial |
| Completed physical time step | The applicable time-discrete balances, admissibility, invariants, and integration-error requirements | A zero stationary residual |
| Certified stationary output | All active stationary equations, constraints, model-validity conditions, and exact output-state consistency | None of the mandatory active conditions may be skipped |

The local equilibrium equations for eliminated species must still be evaluated
accurately enough to define the outer residual. This separation is not permission
to accept a class-4 chemical non-root as a completed physical state. Likewise,
an outer iterate with a nonzero species transport residual is not automatically
a class-4 local chemical non-root: these are different equations and labels.

**Mathematical check executed here:** for `F(x) = x^2 - 2`, starting at `x = 1`,
the first Newton trial is `x = 1.5`. It is positive and reduces the squared
residual by a factor of `0.0625`, but its residual is `0.25`, above a final
tolerance of `1e-8`. All 20 tested backtracks from that first direction also
fail final convergence. Ordinary Newton, without imposing the final gate on
each trial, reaches the tolerance in four iterations.

This is a counterexample to a possible acceptance-policy interpretation, not a
measurement of a failure in the current EXHALE solver. The primary
[KINSOL mathematical documentation](https://sundials.readthedocs.io/en/latest/kinsol/Mathematics_link.html)
also separates line-search decrease conditions from final nonlinear stopping
criteria. It is cited for numerical methodology, not as a proposed new EXHALE
dependency.

**Suggested replacement for A2's last policy bullet:**

> Use the common evaluation to supply separate decisions for probe validity,
> nonlinear trial admissibility, physical step acceptance, and stationary
> certification. Replace the uncertified-cell count with valid local-closure
> checks in trial handling; apply stationary residual tolerances only when
> certifying a stationary state.

### F2. Pure H/He source updates must not bypass the common energy corrections

**Priority: high. Plan lines 272-278.**

D-physical requires B2 and B3a-B3c for configurations with molecular or metal
sources. That does not cover a pure hydrogen/helium configuration without metals.
Its photoionization, recombination, electron production, pressure reconstruction,
and thermal update still require physically consistent energy accounting.

**Source-confirmed:** in `EXHALE_main.f90:1222-1241`, `ioniz_eq` is followed by
`get_species_densities`, a pressure reconstruction at the supplied temperature,
and `W_to_U`. There is no molecular-or-metal guard around this common sequence.
The default semi-implicit energy call also occurs on the atomic path.

In the monatomic limit, at fixed temperature, changing the particle count changes
the thermal energy by

```text
Delta u_thermal = (3/2) k_B T Delta n_particles.
```

Ionizing hydrogen changes that count even without a molecule or metal. The
absence of oxygen work is therefore not evidence that this common source path
is ready for physical validation.

The leading statement that every active gate must pass is correct, but the
explicit list after it should not create a narrower exception.

**Recommended prerequisite rule:**

- Any configuration using reactive composition or thermal source updates needs
  the corresponding source-energy consistency and output-state gates, including
  the pure H/He case.
- Require B2 where the semi-implicit thermal path is used. A separately
  validated explicit path needs its own applicable accuracy and admissibility
  tests; unused code is not a scientific prerequisite.
- Apply B3b's molecular, oxygen, and metal subitems only when they are active.
  Pure H/He validation must not wait for irrelevant H3+ or oxygen implementations.
- Require A3 for active carrier transport, not "A0 to A3 always" when no
  carrier subsystem is enabled. Record inactive conditions as not applicable.
- Require the complete rollback contract wherever recoverable full-step
  rejection is claimed. A stop-safe intermediate implementation must retain
  that explicit limitation; successful accepted steps are still judged by their
  equations, not by whether a retry happened to be necessary.

Use B1a's feature-to-equation inventory to generate this prerequisite list. This
is more precise than using "atomic," "metal," or "molecular" as broad proxies
for which common code paths execute.

**Required focused check:** a pure-H/He case with changing ionization and no
metals or molecules must not be eligible for physical certification while the
temperature-preserving composition reset is still part of its uncorrected
source update. This review inspected that path but did not run the atmosphere.

### F3. A0 needs an implementation milestone and a local-time compatibility rule

**Priority: high for physical transient claims. Plan lines 63-74.**

A0 now states the correct distinction between initialization, physical
integration, and stationary certification, but it is still described as a short
document written by the advisor and accepted by the user. No named implementation
and runtime test applies that distinction to the clock, budgets, restarts, and
the existing time-stepping option.

This is a completion criterion missing from A0, not an objection to writing its
contract first. Add an implementation milestone after the contract is accepted.

**Source-confirmed existing behavior:**

- `input_read.f90:723-729` recognizes `Time stepping: Local` and explicitly
  identifies it as pseudo-time stepping.
- `eval_dt.f90:42-47` returns a distinct `dt_loc(j)` in that mode, whereas the
  global mode fills the array with one common step.
- `EXHALE_main.f90` uses `dt_loc` in the hydrodynamic stages and passes it to
  element diffusion, carrier transport, and energy integration.
- The existing `ieq_state_marching`, `ieq_state_steady_iterate`, and
  `ieq_state_steady_candidate` tags in `ionization_equilibrium.f90:173-176`
  distinguish sweep contexts and their diagnostic records. They are not, by
  themselves, the new A0 contract for a physical clock and initialization handoff.

Different cells taking different pseudo-time intervals do not constitute a
single physical trajectory at one global time. In a simple two-cell update with
one shared internal flux and closed outer boundaries, unequal unsynchronized
time weights also prevent the internal flux contributions from canceling in
the material sum.

**Mathematical check executed here:** with unit cell volumes, initial contents
`1` and `1`, internal flux `1`, and local steps `0.1` and `0.2`, the updated
sum is `2.1` rather than `2`. A common step of `0.1` gives `2`. This demonstrates
why the local pseudo-time update cannot be labeled an ordinary conservative
physical step. It does not invalidate its use as a stationary acceleration
method whose final solution is independently certified.

**Recommended A0 implementation requirements:**

1. Bind run mode to the actual time-step strategy. Keep the existing local
   pseudo-time strategy in initialization/steady relaxation, not physical
   transient mode. Supporting physical local subcycling would require a separate
   synchronized conservative algorithm and is not required for this plan.
2. Distinguish numerical iteration counts, attempted steps, and accepted
   physical elapsed time. Advance physical time only after a complete accepted
   physical step; do not substitute the minimum pseudo-step for elapsed time.
3. Separate initialization and physical cumulative records. Preserve useful
   initialization diagnostics, but do not silently count them as physical
   reaction, heating, or observational histories.
4. Make the handoff rebuild its dependent state and declare the physical time
   origin. Define restart handling so a saved relaxation snapshot cannot silently
   become a physical history with an invented clock.
5. Add tests for incompatible mode/stepping combinations, handoff failure,
   rejected-step time accounting, and restart metadata. The proposed physical
   temporal-refinement tests must use a physical time-stepping mode.

A run-mode input/API change should be presented for confirmation before it is
implemented, as already required for changes to public user inputs.

### F4. Fix dependency joins and remove an unintended dependence on optional Newton

**Priority: medium to high. Plan lines 288-298.**

The prose now correctly states that direct Newton support is optional when a
time-converged state passes independent certification. The graph nevertheless
puts `B6 validity flags` after `B5 stationary system`, while D-physical requires
B6. Read literally, this makes the optional direct-solver work mandatory again.

The graph also shows `B4 -> B5` without an explicit join from B3c's final source
assembly. B5's prose requires the same sources as marching, so its completed
certification must depend on those source implementations too. Work on its
interfaces can proceed earlier; claiming it is complete cannot.

Finally, a single `C1 -> C2 -> C3` chain suggests that every configuration waits
for WASP-121 b candidate C, although C3's prose correctly makes each switch
depend on its own field checks. Unrelated planets need not wait for C2.

**Recommended changes:**

- Define validity-flag semantics with B1a/A0, provide the minimal reporting and
  certification interface with A2, and add the process-specific producers as
  those processes are implemented. B6 completion should not depend on B5.
- Make B5 certification depend on the applicable completed B3c and B4
  assemblies, plus the current A2 evaluator and their common conventions.
- Make C2 a prerequisite only for switching the WASP-121 b configuration that
  selects C. Other switches depend on their own C1/data checks.
- Keep the B2, B3a, B3b, and B3c arrows explicit so a graph branch is not
  mistaken for permission to recover a thermal failure without B3a.

These are graph corrections to match the improved prose, not a requirement to
serialize all independent implementation work.

## 5. Implementation notes that do not require another redesign

### 5.1 B3a can be implemented before every physical defect is corrected

The proposed controller can first be tested as a software restoration mechanism
against the existing update. However, the old composition reset cannot suddenly
pass a new exact physical energy-budget gate simply because it is inside a
checkpointed trial. Keep restoration tests separate from successful physical
source acceptance until B3c removes that reset. Do not loosen the physical
energy test to make an intermediate controller demonstration pass.

The B1 checkpoint must include all state-changing operations before the final
adoption boundary. In the current main program, viscosity/conduction and optional
Shapiro filtering occur after the energy source update. The word "whole" must
cover those active operations too, or state a narrower intermediate contract.
Their physical consistency and contributions belong to the applicable equation
and validity checks, not just the memory-restoration test.

### 5.2 Instantaneous output rates and integrated step budgets are different quantities

B3c's common final-state assembly is appropriate. Interpret "one source
evaluation" as one consistent instantaneous state/field evaluation for the
outputs and relevant endpoint checks. A multistage or substepped physical update
can have a time-integrated source budget that requires rates at several time
levels. Do not replace that integral by an endpoint rate merely to make every
reported number share one evaluation call.

Both quantities should use one definition of each physical process. Their time
levels, units, and meanings must remain explicit.

### 5.3 Validity flags need meaning, not only a blanket all-clear condition

Distinguish active unvalidated physics, an out-of-domain closure, a rejected
numerical trial, and a legitimate external thermal reservoir. A rejected floor
attempt that leaves no contribution in the adopted state need not invalidate
the final result. An unbudgeted accepted energy correction must invalidate it.
Likewise, a physically specified reservoir is not invalid solely because an
informational flag reports its activity.

Initialization diagnostics and physical accepted-history flags should be kept
separate under A0. The cumulative CO record must remain auditable, with a stated
policy for the state handed to physical integration and for ceiling activity
during accepted physical evolution.

### 5.4 The SED correction is appropriately staged

The current script now prints `F(912-2600 A)` for its triplet-band diagnostic,
and the sampled README requirement statements identify the current threshold
near 2600 A. These corrections are confirmed in the files and should not be
relisted as unresolved wording defects.

The script still checks all stored tables against `W_KI` rather than an active
configuration inventory. Rev1 explicitly assigns the configuration-specific
work to C1, so this is an implementation task still pending, not a contradiction
in the revised plan. A/B/C rate comparisons are now correctly called sensitivity.
No additional SED or atmospheric calculation is required to close this plan
review.

## 6. Minimal corrected dependency structure

```text
A0 contract -> A0 implementation and mode/restart tests
A1 acceptance contracts
B1a active-equation inventory + validity-state interface
    -> A2 evaluator with context-specific acceptance
    -> carrier-local checkpoint -> A3 for active carriers

B1 approved target
    -> B2 stop-safe thermal solve
    -> B3a full attempted-step controller
    -> B3b applicable rates, energy recipients, and thermodynamics
    -> B3c coupled source and final-state assembly

B1 approved target -> B4 spatial and boundary implementation
B3c + B4 + current A2 -> B5 stationary implementation/certification
Validity producers/reporting develop with their processes, not after B5

C1 -> each validated existing spectrum -> its C3 configuration switch
C1 -> C2 provisional WASP-121 spectrum -> WASP-121 C3 switch

Applicable source/transport gates + valid mode + validity state + input/output
    -> D-physical
    -> qualified D-comparison and explicitly requested D-figures
```

This diagram states completion dependencies. Independent source and spatial
implementation work may proceed concurrently where their approved interfaces
allow it. Nothing here authorizes new simulations or configuration changes in
the present review task.

## 7. Additional focused tests

The existing rev1 matrix is useful and need not be replaced. Add these narrowly
targeted cases:

| Target | Test | Required result |
|---|---|---|
| A2 context separation | A positive, merit-decreasing Newton trial whose global stationary residual is still large | It can be accepted as a numerical iterate, but not certified as the final stationary solution. |
| A2 physical time step | A valid transient with a nonzero stationary source/balance | The time-discrete step can pass without claiming stationarity. |
| A2 inner closure | A trial with invalid eliminated-species chemistry | The outer calculation cannot treat the invalid closure as a valid physical root; controlled probe/initialization behavior remains separately identified. |
| D eligibility | Pure H/He reactive case, molecules and metals disabled | Applicable common energy/source gates still required; irrelevant oxygen/H3+ gates marked not applicable. |
| A0 implementation | Physical mode combined with existing local pseudo-time stepping | Refused or explicitly kept in a numerical relaxation mode; no physical-time claim. |
| A0 handoff/restart | Initialization snapshot, physical handoff, failed trial, and reload | Mode and time origin are explicit; only accepted physical intervals enter physical time and budgets. |
| B6 without B5 | Time-converged marching configuration with no direct Newton support | Validity flags and A2 remain available; no artificial dependence on optional B5. |
| Final source assembly | Endpoint source outputs alongside a multistage/substepped energy budget | Each quantity uses the correct state/time level and the same physical term definitions. |

## 8. Executed checks and source record

Run the new mathematical checks with:

```bash
python3 docs/audit_20260905/plan_review2_checks_20260906.py
```

The script completed with exit code 0. Direct results were:

- Newton trial: `x = 1.5`, residual `0.25`, squared-merit ratio `0.0625`.
  All 20 tested first-iteration backtracks fail final convergence; ordinary
  Newton reaches `1e-8` in four iterations.
- Linear transient example: time-discrete residual `2.220446049e-16` while
  the steady source is `-1`. A nonzero stationary condition does not make
  the physical step invalid.
- Closed two-cell flux example: initial sum `2`, unequal-step sum `2.1`,
  common-step sum `2`.

These are mathematical demonstrations of the required distinctions. They are
not measurements of the magnitude or frequency of atmospheric errors in EXHALE.

The reviewed HEAD was `35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`. Relevant hashes:

```text
docs/PLAN_20260906_rev1.md
0a117ceb9382cbf1312f88daf2818024a89ec4f7244e1b298dd1b5047ed8ee68
src/modules/time_step/steady_newton.f90
4efff857e10f0c022b9db599c98625b0007d20bd9cc2175dead00ad5a624587d
src/modules/lower_atmosphere/diffusive_photochemistry.f90
43be19aa481a875c3d30e2d4b966974f064dc279c92a250a2b7c5a69118ecbca
src/modules/radiation/ionization_equilibrium.f90
8e5dec6fa2a5af22ee2f7a9d9658b3a927bcd2524c8d528400b02699b0fc4914
```

The three production-source hashes match those inspected in the preceding
review. This supports retaining the earlier source findings as pending work,
without pretending that planning changes have already fixed them or rerunning
their unrelated atmospheric cases.

## 9. Recommended decision

Accept rev1's physical direction and begin its specified contracts, inventory,
and bounded repairs. Before implementing A2, state explicitly that one evaluator
serves different acceptance contexts. Before claiming physical transient support,
implement and test A0 against the actual time-step options and restart behavior.
Before D-physical, replace the molecular/metal-only prerequisite rule with active
code-path requirements that include pure H/He source-energy consistency. Correct
the source/spatial/validity dependency joins so optional direct Newton support
does not become mandatory indirectly. With those amendments, the revised plan
provides a reasonable basis for the next implementation stage.
