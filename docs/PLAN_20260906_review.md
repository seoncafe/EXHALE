# Review of PLAN_20260906.md

Date: September 6, 2026, KST.

## 1. Verdict

The updated plan adopts the important corrections from the previous review and
is substantially more defensible. Its separation of the implementation audit
from the target governing system is correct. The proposed treatment of density,
chemical energy, Lyman-alpha energy ownership, CO exclusion, and provisional
stellar spectra should be retained.

**Recommendation: conditional approval, with the dependency and acceptance
conditions below corrected before implementation reaches the affected step.**
The plan is suitable for continuing the bounded contract repairs and writing
the target specification. It is not yet an adequate schedule for certifying
molecular atmospheric results.

The principal remaining problems are:

1. Step D can start after A, B1, and C3, before the source-energy corrections
   or the conservative spatial operator are implemented.
2. Step B2 requires rejection of a failed thermal update, but the whole-step
   recovery controller is scheduled after the coupled source work in B3.
3. The temporary class-4 policy is described both as retaining amnesty and as
   never accepting a non-root physical state. Its operating mode and transition
   out of initialization are not specified as an executable milestone.
4. A2 needs a complete active-equation inventory, an explicit evaluation status,
   and a state-consistent residual interface. Existing carrier residual code
   cannot safely be connected to a final gate solely by reading a returned norm.
5. Several mandatory rev3 physics tasks, especially the spatial operator, are
   retained by reference but have no explicit place in the replacement schedule.

These issues do not require discarding the plan. They require sharper contracts
and a corrected dependency graph.

## 2. Scope and verification

The complete 226-line [updated plan](PLAN_20260906.md) was read, together with
the relevant portions of
[rev3](development_plan_20260905_rev3.md),
[the superseded execution schedule](development_plan_20260905_execution.md),
[D0](d0_governing_system_20260906.md), and
[the previous recommendations](To_be_determined_by_user_recommend_20260906.md).

The current production paths inspected include carrier transport and its
stationary residual, ionization acceptance, steady-solver trial acceptance,
energy integration, final output, element diffusion interfaces, the stellar
input flags, and the SED coverage script. Existing routines were distinguished
from their use by a particular driver. In particular, a carrier stationary
residual already exists and handles more than H2; this review does not propose
recreating it simply because the marching gate does not currently use it.

The repository remains a work in progress. The inspected HEAD was
`35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`, with substantial working-tree changes.
The reviewed plan's SHA-256 was:

```text
47196e0373f3c12bcc99e3b29468c5a59fd1be799365cdbc667847593083cc31
```

Evidence in this review is classified as follows:

- **Confirmed in source:** an executable branch, formula, or caller was inspected.
- **Plan defect:** a contradiction, missing dependency, or insufficient acceptance
  condition in the written plan.
- **Mathematical check:** a directly executed analytical counterexample, not an
  execution of the Fortran code.
- **Unmeasured:** atmospheric impact, branch frequency, or behavior requiring a
  production run that was not performed here.

No production source, existing plan, SED, input configuration, or scientific
output was modified. No executable was rebuilt, no atmospheric case was rerun,
and no regression reference was refreshed. The reported binary checksum and
sixteen/seventeen-case results in the plan remain historical claims, not new
measurements from this review.

The new analytical checks are archived in
[plan_review_checks_20260906.py](audit_20260905/plan_review_checks_20260906.py).
The earlier SED measurements were not repeated because neither their data nor
their numerical method was changed for this review.

## 3. Changes that are reasonable and should be retained

| Updated decision | Assessment |
|---|---|
| Carrier-only recovery before a whole-step controller | A reasonable incremental implementation, provided its own rollback contract precedes it and its physical scope is stated. |
| Removing a residual cap near one as a chemical acceptance tolerance | Correct. A previously completed relaxation does not define chemical accuracy. |
| Treating D0 as an audit and writing a separate target document | Correct and necessary. Describing an inconsistent implementation is not approval of that implementation as physics. |
| Keeping thermal plus kinetic energy with an explicit chemical/excitation ledger | Physically valid if the EOS, source update, species fluxes, and radiation ownership are derived consistently. |
| Chemistry preserves supplied density | Correct. Matching a nonconservative mass reset in the stationary equations would not repair the physics. |
| Deriving excitation loss and de-excitation return together | Correct. A small escape probability alone is not a reason to multiply existing cooling by that probability. |
| Excluding ceiling-active CO and incomplete oxygen-energy configurations from validated results | Appropriate until the omitted kinetics and energy terms are supplied or a defensible reduced model is derived. |
| Candidate C labeled provisional, with the chosen XUV normalization explicit | Appropriate. This resolves the overclaim, not the uncertainty in the input spectrum. |
| Separate Salz-comparison and observed-star HD 189733 b models | Appropriate; the reference experiment and observational epoch must remain explicit. |
| Preserving old outputs and keeping scientific reruns instruction-gated | Appropriate. Configuration changes must not silently relabel old results. |
| Keeping historical audit code | Appropriate. A frozen reproducer can retain its old interface; current tests should be separate and its compatible source version documented. |

## 4. Findings requiring changes to the plan

### R1. Step D has insufficient prerequisites

**Priority: high; plan lines 207-213.**

Step D is allowed after A, B1, and C3. B1 is a specification, not an
implementation. At that point the current temperature-preserving composition
reset, the two-iteration thermal solve, and other acknowledged molecular energy
defects may still be active. A2 can establish that a state satisfies a particular
implemented equation set; it cannot make that equation set physically correct.

Nor does a small steady residual establish validity of the model assumptions,
such as an active CO ceiling or omitted oxygen reaction energy.

**Required change:** split D into distinct purposes:

- **D-diagnostic:** bounded, explicitly requested development calculations may
  run earlier, but carry the known limitations and cannot be presented as
  validated physical predictions.
- **D-physical:** a configuration may enter physical validation only after the
  source terms, thermodynamics, transport closure, output state, and input field
  needed by that configuration have passed their gates.
- **D-comparison:** a published reproduction additionally requires matched
  boundary conditions, geometry, composition, rates, and spectral assumptions,
  or an explicit statement of the unmatched terms.

Do not require B4's direct Newton capability for a configuration that has a
properly time-converged, independently certified marching solution. Conversely,
do require the applicable spatial-operator and energy corrections even if a
direct Newton solver converges quickly. Atomic configurations can proceed
without unfinished oxygen work only when the omitted physics is outside their
active model.

### R2. Thermal rejection is scheduled before its recovery mechanism

**Priority: high; plan lines 155-168.**

B2 says an unacceptable energy update is rejected. B3 schedules the whole-step
controller only after the coupled species-energy solve. The current main program
finishes its hydro retry loop before element diffusion, carrier transport,
ionization, and energy integration. A3 restores only a carrier substep, not the
state altered by all those preceding operations.

**Required change:** choose one explicit intermediate contract:

1. Install the complete attempted-step controller before, or atomically with,
   the first B2 update that promises recovery; or
2. Let B2 return a failure status and stop safely until the full controller
   exists. Do not describe this intermediate build as supporting rejected-step
   recovery; or
3. Derive and validate a strictly local thermal-substep integrator on a frozen
   background, covering the requested source interval. Label that as partial
   integration, not rollback of the full step.

The preferred target remains option 1. B1 must enumerate the state before the
corresponding controller is implemented; a list written but not implemented is
not recovery.

Also narrow the floor policy: a floor can represent an independently specified
physical reservoir only if that reservoir has a physical basis. Merely recording
the energy needed to force a numerical temperature floor does not turn an
otherwise unphysical correction into valid external heating.

### R3. The class-4 transition needs its own milestone and run mode

**Priority: high; plan lines 39-45 and decision table row 2.**

The prose retains amnesty temporarily while the table says a class-4 iterate
is never an accepted physical state. Those statements can coexist only if the
intermediate run is explicitly a numerical initialization/continuation mode,
not the production physical time trajectory. The execution steps do not say
when that mode is implemented or how it hands a state to physical integration.

Rejecting a non-root only at the end does not retroactively validate all
intermediate physical steps. A final stationary solution can still be valid if
it satisfies the complete target equations independently, but the preceding
relaxation time history then has a different status.

**Required change:** add an explicit initialization milestone with:

- a run-state distinction between initialization, physical integration, and a
  certified stationary solution;
- no accumulation of physical elapsed time, reaction budgets, or observational
  histories from nonphysical initialization iterates;
- an admissible, source-consistent handoff state, with caches initialized from
  that same state;
- failure handling if initialization cannot supply such a state;
- separate handling of diagnostic snapshots and scientific output.

A valid finite-time state need not have a zero *stationary* residual. A2 must
not prevent writing a valid time-dependent solution simply because it is not
stationary. Conversely, a fixed-step relaxation snapshot must not be certified
as a physical transient merely because it is labeled a snapshot.

### R4. A2 must cover the actual independent equation set

**Priority: high; plan lines 95-106 and 170-174.**

The phrase "every transported species" is good, but its parenthetical lists
only H2 and H+. The current `carrier_set_init` in
`diffusive_photochemistry.f90:589-607` also activates OH, H2O, and CO when the
oxygen network is enabled. `carrier_steady_residual` already iterates over the
active carrier mask and can return rows for that set. The element-diffusion
module also evolves elemental partitions through separate equations.

**Required change:** have B1 define, and A2 consume, an active-equation inventory
rather than a hard-coded pair of species. Include:

- each independent transported carrier balance;
- active elemental transport balances, including the helium/hydrogen partition
  and trace-element transport when enabled;
- the local algebraic closures only for species actually eliminated by
  equilibrium or another justified fast-process approximation;
- mass, charge, energy, boundary, and model-validity conditions.

Do not require the net chemical source of a transported species to vanish.
Its steady equation is a balance between transport and reaction:

```text
div(n_s v + J_s) - omega_s = 0,
```

not `omega_s = 0`. In a reacting flow, both terms can be nonzero. Likewise,
when B4 adds global unknowns, use B1's independent species space so that total
mass, element constraints, and species rows are not imposed redundantly.

The relevant element-transport code already contains balance calculations and
iteration diagnostics. The task is to expose or assemble an appropriate
state-consistent stationary check, not to assume the underlying equations are
absent because they are outside the current steady driver.

### R5. The existing residual interface can return an apparent zero without evaluating

**Priority: high; an additional source finding directly affecting A2.**

In `diffusive_photochemistry.f90:1036-1044`,
`carrier_steady_residual` initializes its norms and optional row arrays to zero,
then returns if molecular chemistry, carrier transport, or `bg_ready` is false.
Inactive physics may legitimately need no test. An active carrier model with
unavailable background data is different: its residual has not been evaluated.

The routine also updates internal background arrays, photolysis arrays,
advection corrections, row scales, and H2 headroom. Its arguments being
`intent(in)` does not make it free of side effects.

**Required change:** A2 needs an evaluation status such as not applicable,
evaluated, or unavailable/invalid, alongside the residuals. An active but
unevaluated equation cannot pass because its initialized norm is zero. Require
finite values and evidence that every required equation was evaluated on the
same state and field.

Choose either a side-effect-free certification evaluator or an isolated
workspace whose temporary changes do not alter the adopted state. Do not
re-equilibrate a transported partition merely to evaluate its residual. Make
the thermodynamic and radiation backgrounds explicit or verify their association
with the state being certified.

**Required tests:** active transport with `bg_ready` false; an inactive model;
stale background after changing the state; repeated evaluation with diagnostic
output enabled and disabled; preservation of the state being certified. The
present review confirmed the source behavior, not a newly executed false
atmospheric convergence.

### R6. A1 must check the returned residual even when status says converged

**Priority: high; plan lines 80-93.**

The proposed constrained-cell test fixes one real defect, but the test list
concentrates on a *failed* solver status. The current solver can set a successful
status before `limit_to_element_budget` changes the state; it then recomputes
the residual without revising that earlier status.

There is another acceptance path in
`diffusive_photochemistry.f90:2176-2200`: a stagnating iteration may be marked
`accepted_at_floor` when its residual has fallen by `newton_drop = 1e-6`
relative to the initial residual, even if it exceeds `newton_floor = 1e-8`.
A decrease relative to an arbitrarily large starting residual is not an
absolute or physically scaled accuracy requirement.

**Mathematical check:** the predicate admits an iteration with initial residual
`1e8`, previous residual `100`, and current residual `60`. At iteration 3,
`60 > 0.5*100` and `60 <= 1e-6*1e8`, so the relative-drop branch is true while
the absolute residual is nowhere near `1e-8`. This is a counterexample to the
predicate, not a claim that such a sequence was measured in an EXHALE run.

**Required change:** separate nonlinear termination, returned-state equation
accuracy, and physical constraint validity. Every successful physical update
must satisfy the specified returned-state conditions, regardless of which
termination branch fired. A real round-off/conditioning limit can be documented
using a backward-error or state-error argument; an arbitrary initial residual
cannot establish that limit.

Add tests for a successful status followed by an active limiter, a large
relative reduction with an unacceptable final residual, non-finite rows, and
an otherwise valid constrained state. Keeping the raw status is useful for
diagnosis, but it is not the complete acceptance decision.

### R7. A complete norm cannot be just one sum over all species

**Priority: medium to high; A2's acceptance-test design.**

The current optional `rvol` in `carrier_steady_residual` sums absolute residuals
and row-term scales over all active carriers within each of two spatial regions,
then takes the larger regional ratio. The routine also supplies a worst-row
measure and optional individual row arrays, so stronger tests can use existing
information.

A combined ratio can hide an unbalanced species whose reaction or transport
terms are smaller. The analytical check included here gives two row scales
`1e8` and `1`, with residuals `0` and `1`. The combined ratio is approximately
`1e-8`, while the second row's relative imbalance is `1`.

**Required change:** state the units, absolute floors, relative scales, and
spatial domains of every gate. Require a condition for each active independent
balance, with a local safeguard where fronts or thin layers would otherwise be
averaged away. Absolute tolerances must avoid meaningless relative failures for
vanishing species without discarding physically important chemistry. Do not
calibrate thresholds solely to make the existing matrix pass.

### R8. Mandatory rev3 physics needs an explicit schedule crosswalk

**Priority: high for molecular validation; plan lines 121-174.**

The new plan states that rev3 remains authoritative, but replaces the old Steps
3-7. Naming a target equation and owner in B1 does not schedule its implementation.
The replacement schedule explicitly maps B4 to rev3 Phase 5 while omitting an
explicit implementation step for rev3 Phase 4, the spatial operator.

The following tasks must have a named implementation and validation location:

| Retained rev3 work | Required placement in the new execution order |
|---|---|
| Photoevent energy partition, including H2 channels and excited recipients | Before certifying a source step or atmosphere using those channels. |
| Reaction heats, including oxygen and associative channels | Before the corresponding network enters the validated set, or with an explicit exclusion. |
| A common H2 state model and consistent thermodynamic functions | Alongside or before validating coupled molecular chemistry and energy. |
| H3+ absolute emission model, density limits, joins, and temperature domain | Before accepting H3+-dependent physical energy balances. |
| Finite-rate evolution of species not legitimately eliminated | An explicit increment of B3, not silently replaced by solving all species in equilibrium. |
| Conservative species face fluxes with thermal/chemical energy transport | A spatial-operator implementation step before the matching stationary system is certified. |
| Charged transport and background response for diffusive mass conservation | With the spatial operator for configurations that use those processes. |
| Boundary constraint count, hydrostatic residual, and applicability diagnostics | Implemented and tested for the active independent system, not only listed in B1. |
| Steady Jacobian action, inner-solve determinism, reload/acceptance consistency | Explicit gates of B4 and the final-output path. |

These tasks have not disappeared from rev3. The issue is incomplete scheduling,
not evidence that the authors intentionally removed them. Add a crosswalk so
that "B3 complete" or "B4 complete" cannot be used to bypass them.

### R9. A1's byte-identity gate is not justified by final-state chemistry

**Priority: medium; plan lines 91-93.**

The claim that final default cases have no failed chemistry does not establish
that correcting class-5 classification leaves their trajectories unchanged.
`steady_newton.f90:630-650` uses the class count during trial acceptance, not
only at the final output. A valid class-5 root can therefore change which trial
is accepted even when the eventual old solution contained no class-4 cell.

Similarly, stricter checks can reveal that a previously labeled converged case
does not satisfy an added species condition. A2's statement that converged cases
pass should mean *cases independently shown to satisfy the required equations*,
not every case carrying an old convergence label.

**Required change:** replace byte identity as a gate with a scoped impact
measurement. Expect identity only when inspection and branch evidence show the
changed policy was unused. Record changed outcomes and investigate their physical
and numerical reasons. Neither a changed snapshot nor a previously convenient
convergence label is a reason to retain an incorrect acceptance condition.

Also replace "movement on every case" for every small step with affected-path
validation. Broader regression is justified for a shared source or acceptance
interface, but an unrelated documentation or SED-path change does not require
regenerating an atmosphere.

### R10. Successful implicit solves do not establish time accuracy

**Priority: medium; A3, B2, B3, and their verification matrix.**

The rejection and retry design addresses nonlinear failure. It does not by itself
control truncation error or splitting error. An implicit reaction update can
converge accurately to a poor approximation over a large physical interval.

**Mathematical check:** for `dy/dt = -y`, `y(0) = 1`, one backward-Euler step to
time 10 gives `0.09090909091` with algebraic residual `1.11e-16`. The exact
solution is `0.00004539992976`. The nonlinear equation is solved accurately while
the physical time evolution is poorly resolved. Refining the time interval in
the archived check decreases the error.

**Required change:** include temporal refinement at fixed physical end time,
splitting-order checks for the coupled operators, and an error policy distinct
from the nonlinear convergence tolerance. A step-doubling estimate or a justified
fixed-step refinement study can serve this purpose; a change of integrator is
not mandated. The official
[CVODE mathematical description](https://sundials.readthedocs.io/en/latest/cvode/Mathematics_link.html)
provides a useful primary example of separating nonlinear convergence from the
local integration-error test. This is supporting numerical methodology, not a
claim that EXHALE currently uses CVODE.

For carrier-only substeps, identify which coefficients stay frozen and which
depend on the new carrier state. Do not rebuild a supposed fixed background
from a partially accepted state by repeatedly calling the entire transport
wrapper. If a later substep fails permanently, either restore the original
interval entry state or retain a clearly marked partial diagnostic state; never
report that the original interval was completed.

### R11. B3 does not automatically make output_state_consistency pass

**Priority: medium; plan lines 163-168.**

The existing test compares the heating in `Hydro_ioniz.txt` with the independently
assembled total in `Heating_breakdown.txt`. The current final-output path calls
`excited_H_update`, writes the main state with its existing `heat` and `cool`,
and subsequently calls `write_heat_breakdown_eq`, which evaluates the radiation
and rates again (`EXHALE_main.f90:1933-1955`; test header in
`src/tests/grid_and_gates/output_state_consistency.sh`).

A coupled local species-energy solve does not automatically synchronize an
external radiation iteration, an excited-state cache, or separate output
assemblies. The test's two heating totals can also agree while temperature,
species, or the actual energy equation remains wrong.

**Required change:** add a final-state assembly contract that supplies the same
state and physical source evaluation to certification and output. Verify
nonmutation, consistency of the EOS and species values, channel sums, and
the associated radiation/coupling settings. Retain the existing test, but treat
it as one necessary check rather than proof that the entire coupled model is
correct. A diagnostic recomputation must not silently change the state after
it has been certified.

### R12. The SED work is improved, but its claims need two qualifications

**Priority: medium; decision row 4, Section 2, and C1-C3.**

First, rates calculated for A, B, and C establish model sensitivity. They do
not measure C's bias relative to the unknown target-star field. The wording
"the He 2^3S bias ... is measured" should become "differences in the production
He 2^3S rates are measured; bias requires a justified reference field."

Compare those rates on the same prescribed background and attenuation when
isolating spectral differences. Separate that comparison from an unattenuated
orbital rate and from fully reconverged atmospheres. Do not imply that different
quadratures converging on the same proxy establish the proxy's physical accuracy.

Second, the coverage correction is **partial**, not complete:

- The script's triplet constant is now `4.767775 eV`, as the plan states.
- Its `main` still sends every file to `check_file(..., W_KI)`, rather than
  deriving the required endpoint from each active configuration. It has no
  configuration-specific H(n=2) coverage test.
- The printed label `F(912-2583 A)` now encloses a calculation using `W_HETR`,
  which corresponds to about 2600 A. The label and the computed integral
  therefore disagree.
- The README's opening description was corrected, but several later statements
  still identify 2583 A as the current metastable threshold. Historical band
  measurements may retain their original endpoints; statements of the current
  physical requirement must not.

The current parser arms excited-hydrogen coupling when both `T_star_eff > 0`
and `R_star > 0` (`input_read.f90:1155`). Do not reduce this to a check of
`Stellar Teff` alone. A spectrum extending to approximately 2995 A covers part
of the Balmer continuum, but not its complete required band to approximately
3647 A. Replace "covers no H(n=2) band" with "does not cover the full required
Balmer continuum." The production coverage guard already exists; the remaining
work includes the test inventory and accurate documentation.

Finally, the C branch is independent for data preparation, but not permanently
independent of B: changes to cross sections, photoevent energy accounting,
active absorbers, or source assembly require the affected consumer tests to be
rerun. This does not authorize regenerating unchanged SED data.

## 5. Additional implementation cautions for the approved direction

### 5.1 A carrier-only checkpoint must be specified before A3

It is reasonable to postpone enumeration of the entire main-program state to
B1. It is not reasonable to postpone enumeration of the smaller state already
mutated by A3. The carrier module has saved collision backgrounds, photolysis
arrays, advection corrections, constraint diagnostics, and cumulative CO
statistics. Define the checkpoint or local trial workspace for those before
coding the first retry.

The cumulative CO record now exists in the current source, including a public
query and resolved-output fields. Do not describe it as absent or reimplement
it independently. Instead, ensure rejected trials cannot increment the accepted
history after retries are introduced. The current limiter updates that record
where it applies the ceiling, before a new outer retry design exists.

### 5.2 A steady gate needs physical boundary and validity information

Element conservation within a cell is not a substitute for correct boundary
fluxes. Specify the physical cells and boundary equations used by each norm.
Avoid demanding that an imposed-reservoir ghost cell satisfy an unrelated
interior evolution equation. Conversely, do not exclude an inconsistent physical
base layer merely to obtain a smaller residual.

State which failures mean the chosen physical model is invalid, such as an
unjustified ceiling, a missing energy recipient, or an out-of-domain closure.
Such a state can be a valid diagnostic experiment but cannot acquire physical
validation solely by satisfying the numerical equations that were implemented.

### 5.3 Keep independent tests for the physics, not only mutually consistent code

A shared assembly reduces duplication, but if both the integrator and its test
call the same wrong formula they can agree. Include analytic limits, conservation
identities, published rate/emission data within their validity ranges, and
manufactured transport balances where suitable. Preserve the unit and state
conventions of those references.

This review did not rederive the disputed H3+ fit or every molecular reaction
from the literature again. Their mandatory tests remain in rev3; the present
finding is that the updated execution order must explicitly carry them forward.

## 6. Recommended revised dependency order

Use the following order while retaining the useful incremental structure:

```text
Acceptance contracts and initialization-mode contract
    -> full active-equation certification interface
    -> carrier-local rollback specification -> carrier-only retry

Target governing system approved
    -> applicable rate/energy/EOS implementations
    -> attempted-step rollback implemented before recoverable thermal failure
    -> residual-controlled thermal update
    -> coupled equilibrium/finite-rate species-energy increments
    -> conservative spatial operator and boundary/charged closure
    -> matching stationary equations and solver verification

SED reconstruction semantics -> consumer tests -> provisional C -> input switch
                                      |
                   repeat affected consumer tests after source changes

Applicable physics gates + exact output state + validated input
    -> explicitly requested physical validation
    -> qualified published comparisons
    -> explicitly requested scientific figures
```

The ordering of individual microphysics implementations can be adjusted where
their interfaces permit it. The requirements are that dependencies are explicit
and that no validated result uses an unresolved active term. Direct steady
Newton support is optional when a time-converged reference is available and
passes the same independent physical certification.

## 7. Minimum additions to the test matrix

| Target | Focused test | Acceptance criterion |
|---|---|---|
| A1 status versus returned state | Successful raw solve followed by a limiter; failed raw solve with a constrained worst cell; large relative residual decrease | Every required returned-state equation and constraint passes, or the update is rejected. |
| A1 class contract | Equivalent roots found through classes 1 and 5; genuine class 4 | Root validity depends on the requested equations, not on which solver found the root. |
| A2 equation inventory | Oxygen carriers, H+ transport, and active elemental diffusion | All required independent balances are evaluated; no mandatory row is silently absent. |
| A2 readiness | Active carrier model with unavailable or stale background | Explicit invalid/unevaluated status, never an initialized zero interpreted as convergence. |
| A2 norm | A small-rate species with a large fractional imbalance beside a well-balanced dominant species | The small-rate species cannot be hidden by the combined sum. |
| A2 state identity | Repeated certification before and after optional output calls | Same state and required fields; no diagnostic mutation or undocumented cache dependence. |
| A3 rollback | Failure after an accepted substep and after counter/cache mutation | No rejected contribution in adopted state or accepted history; completed physical interval is correct. |
| B2 thermal solve | Falling cooling branch, strongly temperature-dependent heat, unavailable bracket, floor encounter | Final residual, physical admissibility, and well-defined failure propagation. |
| B3 source model | Closed reaction, isolated photon event, recombination, two-cell radiative exchange | Consistent thermal plus chemical/excitation/radiation accounting. |
| B3 temporal accuracy | Successively smaller physical steps at the same final time | Observed temporal convergence and controlled splitting error, separately from nonlinear convergence. |
| Spatial operator | Diffusion-only column and reacting transport balance | Mass, elements, specified current condition, and material energy close with consistent face fluxes. |
| B4 stationary system | Time-converged reference, Jacobian-action checks, reload | Same active equations and final acceptance on the independent unknown space. |
| C1 coverage | Triplet, low-IP metal, and active Balmer configurations | Correct active endpoints and matching labels; unavailable bands refused. |
| C2 SED comparison | A/B/C through the same production rate at a fixed state | Differences reported as sensitivity; no claim of known bias without a reference. |
| Output | Main fields and all derived products | One identified state and model; certified stationary, valid transient, and diagnostic products clearly distinguished. |

Specify numerical tolerances and units in the corresponding implementation
briefs. No numerical pass threshold should be chosen solely to preserve an
existing snapshot or convergence label.

## 8. Checks executed for this review

The archived script was run directly:

```bash
python3 docs/audit_20260905/plan_review_checks_20260906.py
```

It completed with exit code 0. Its assertions demonstrate:

1. An accurate backward-Euler nonlinear solve can have a large physical-time
   error, which decreases with time refinement.
2. A combined residual norm can hide an unbalanced smaller-scale equation.
3. The current relative-drop stagnation predicate admits a final normalized
   residual of 60 for the specified numerical sequence.

These are analytical counterexamples supporting stronger plan requirements.
They do not establish that the corresponding atmospheric branch occurred in
the regression matrix. No production regression case, new SED, or transit
calculation was run. Existing scientific products were preserved.

## 9. Suggested concise revision statement

Retain the updated physical decisions and audit/specification separation.
Define initialization separately from physical evolution. Extend A1 acceptance
to the returned state irrespective of raw solver status. Make A2 evaluate the
complete independent active system with explicit readiness, finite-value, and
state-identity checks. Specify carrier-local rollback before A3 and implement
the complete recovery contract before claiming recoverable thermal rejection.
Explicitly schedule the inherited molecular energy, thermodynamics, conservative
spatial operator, and boundary tasks. Change D's prerequisites from completed
documents to passed physical gates for each configuration. Treat byte identity
as a scoped regression observation, A/B/C rate differences as sensitivity, and
coverage corrections as complete only after configuration-specific tests pass.
