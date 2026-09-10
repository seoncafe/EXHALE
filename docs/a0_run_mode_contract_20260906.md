# A0: the run-mode contract (PLAN_20260906_rev2, Step A0)

Status: accepted by the user on 2026-09-06 (the explicit `Run mode` key of
section 6 and the contract as written). Implementation is A0-impl. Every statement about
the present code is READ from the source named.

## 1. Why a contract is needed

The code has one marching loop and a step counter, and no notion of which of
three different things a run is doing at a given step:

- `count` is the step counter (raised in the main temporal loop of
  `EXHALE_main.f90`, declared in `parameters.f90`);
  there is no simulated-time variable, so "physical elapsed time" is not
  recorded anywhere today.
- `Time stepping: Local` sets `use_local_dt` (`input_read.f90`, `input_read`,
  commented as "steady-state convergence acceleration; not time-accurate"),
  and `eval_dt` then returns a different `dt_loc(j)` for each cell
  (`eval_dt.f90`); `dt_loc` is passed to the hydro stages, element
  diffusion, the carrier transport and the energy update. Cells advancing by
  different intervals do not form one physical trajectory: in a closed
  two-cell system with one internal flux, unequal steps leave the material
  sum changed (review 2, F3, mathematical check), so a local-dt update is a
  relaxation iterate, not a conservative physical step.
- The ionization sweep accepts non-roots under the "relaxation amnesty"
  (acceptance class 4, `ionization_equilibrium.f90`, the relaxation-amnesty
  note above the acceptance classes and `nonroot_acceptance_class`) in the
  same loop that later produces the state written as a result; the
  `ieq_state_*` tags (declared in the head of `ionization_equilibrium.f90`) separate the
  sweep ledgers of the marching loop, the steady iterate and the steady
  candidate for reporting, but they do not carry a physical clock or a
  handoff.
- The `# coupling:` header written by `write_coupling_state_header`
  (`utilities.f90`) records the switches a restart needs (secondary
  ionization state and step, reconstruction operator, restart schema) but not
  the run mode or a time origin, so a restart of a relaxation snapshot re-enters
  the loop as if it were a physical history.
- Fixed-step snapshots (`EXHALE_MAXSTEPS`, the `mol_*` regression cases) are
  relaxation states by construction and are written with the same headers as
  converged results.

## 2. The three run states

| State | Purpose | What is allowed | What is forbidden | What accumulates |
|---|---|---|---|---|
| **Initialization / continuation** | reach an admissible, source-consistent state from an initial guess; steady relaxation; PTC | local pseudo-time (`use_local_dt`), the PTC route (`EXHALE_PTC`), non-root iterates (class 4) as numerical guesses, damping, filters; NOT a failed energy update, which is no state at all (see below) | any claim about elapsed time, any physical history (reaction or heating budgets, observational time series) | iteration and attempted-step counts, diagnostic ledgers labeled as initialization |
| **Physical integration** | advance the state in time by accepted global steps | one global `dt` per step; a step is adopted only if it satisfies its time-discrete balances, admissibility (positivity, element and charge constraints), invariants and the integration-error requirement (B3, R10); rejection and retry (A3, B3a) | local pseudo-time; class-4 chemistry as an accepted state; any accumulation from a rejected trial | physical elapsed time (advanced only after a complete accepted step, by that step's global `dt`), accepted-step budgets and histories, accepted-step statistics |
| **Certified stationary solution** | a state that satisfies every active stationary equation | A2 stationary certification on the exact written state (section 3 of the plan: all active equations, constraints, validity state, output-state consistency) | writing a state as "converged" without the certification record | the certification record itself |

A valid finite-time state of physical integration is writable and is not a
failure because its stationary residual is nonzero (review 2, F1). A
relaxation snapshot is writable as a diagnostic and is never a physical
transient nor a certified solution.

**The exception to the initialization leniency: a state that does not exist.**
"Non-root iterates (class 4) as numerical guesses" allows the initialization
mode to march on a state whose equations are not solved. It does not allow it
to march on no state at all. A failed energy update
(`energy_semi_implicit`, status other than `ENERGY_UPDATE_OK`) found no
temperature for at least one cell and assembled nothing, so there is no
iterate to keep: it refuses the attempted step in BOTH modes, and the
controller retakes the step at half `dt` (decision of 2026-09-06, item 1 of
`To_be_determined_by_user_20260906.md`). The leniency continues to apply to
the coupled source step that has not stopped moving at its pass cap, which
does return a state.

**Exhaustion is terminal in both modes.** A step refused at every `dt` down to
`dt/2**n_step_retry_max` is the same statement about a physical trajectory and
about a relaxation: the code could not produce an acceptable state from the one
it holds. The run prints the last refusal, the retry history and the restored
last accepted state, and exits with status 2, the status a refused stationary
claim uses (`certification_stop_uncertified`), so one status means one thing to
a caller.

## 3. Transitions

**Initialization to physical integration (the handoff).** Requirements:
1. the handoff state is admissible (positive densities and temperature,
   element and charge constraints met) and source-consistent: the chemistry
   of every cell is a root (classes 1, 2, 3 or 5, not 4), the energy state
   and composition agree with the equation of state, the radiation columns
   and rates are those of this state;
2. every dependent cache is rebuilt from the handoff state (photon columns,
   rate caches, carrier backgrounds `bg_ready`, excited populations,
   thermodynamic caches): nothing carried from the last initialization
   iterate survives except the state itself;
3. the physical time origin is declared (t = 0 at the handoff, or the time
   read from a physical restart);
4. the initialization ledgers are closed and kept as diagnostics; the
   physical ledgers start empty;
5. if no admissible, source-consistent state can be produced, the run stops
   with the reason; it does not enter physical integration.

**Physical integration to certified stationary solution.** Only through A2
certification on the exact state to be written; a stationary Newton finish
(`Solver: Newton`) runs as initialization/continuation on that state (its
trials are numerical iterates) and its result is certified by the same A2
evaluation, never by the solver's own convergence flag.

**Restart.** The headers state the mode of the written state and, for a
physical state, its elapsed time. In `init` mode a restart of any snapshot
re-enters initialization. In `phys` mode a restart of a `mode=phys` state
continues physical time from its recorded `t_phys`; a restart of an
initialization snapshot (`mode=init`, or a file written before this field
existed) IS the handoff of this section: the handoff check runs and the
clock starts at `t_phys = 0`. What is refused is a claim that cannot be
honored: a `mode=phys` header without a finite `t_phys` (no clock can be
invented). `Do only PP: True` makes no time-integration claim and defaults
to `init`. The grid guard of item 2c-RESTART stays.

## 4. Binding of mode to the time-step strategy

| Strategy | Initialization | Physical integration |
|---|---|---|
| global `dt` (CFL) | allowed | required |
| `Time stepping: Local` | allowed | refused (or the run stays in initialization with no physical-time claim; the input says which, section 6) |
| PTC / JFNK (`EXHALE_PTC`, `Solver: Newton`) | allowed | not a physical integrator; runs as continuation |
| fixed step count (`EXHALE_MAXSTEPS`) | a relaxation snapshot | allowed as a bounded physical run only in physical mode with a global `dt`; the written state is then a valid transient, not a snapshot |

Supporting a physically meaningful local subcycling would need a separate
synchronized conservative algorithm; it is not part of this plan.

## 5. Counters, ledgers and headers

- Three counters, never conflated: nonlinear iterations, attempted steps
  (including rejected trials), accepted physical elapsed time.
- Physical time advances only after a complete accepted physical step and by
  its global `dt`; never by the minimum pseudo-step, never by a rejected
  trial.
- Two families of ledgers with the same fields: initialization (kept,
  labeled) and physical (accepted steps only). The carrier CO destruction
  domain record, the energy-floor record, the acceptance-class tallies and the
  heating and cooling budgets all belong to a family; a rejected trial writes
  to neither, with one stated exception: the CO domain record is not rolled
  back, because whether a state lay inside the destruction model's domain has
  an answer whether or not the step that read it was accepted.
- Output headers: `# coupling:` gains `mode=init|phys|certified` and, for
  `phys`, `t_phys=<seconds>`; the certification record (A2) is written with
  a certified state. Every derived product (`_adv`, breakdowns, transit
  inputs) carries the same mode field as the state it was made from.

## 6. Input and API form (confirmed by the user 2026-09-06: the explicit key below)

Proposed, one key and one refusal rule:

```
Run mode: init | phys          # default: init, always; phys only when stated
```

The default is `init` in every configuration (decided 2026-09-06, replacing
the first draft's rule "phys when a global dt is used"): the code's ordinary
use is a stationary solution reached by marching and then certified, which
is initialization/continuation in this contract's terms; a physical time
integration is the exception and is asked for explicitly. No existing input
therefore changes its behavior, and the regression matrix is `init`
throughout without any new line.

- `Run mode: phys` with `Time stepping: Local` is refused at startup with a
  message naming both keys.
- `Run mode: init` runs exactly as today and writes `mode=init` headers; the
  present regression matrix (relaxation snapshots and `du`-stopped runs) is
  this mode and stays byte-identical in its numeric content.
- `Run mode: phys` enables the physical clock, the physical ledgers and the
  handoff check at step 0 (the initial condition must itself be admissible
  and source-consistent, or the run must first run an `init` phase; the
  two-stage form `Run mode: init->phys [handoff criterion]` is a possible
  later extension and is not proposed now).
- `Load IC? True` reads the mode and `t_phys` from the header and applies
  the restart rule of section 3.
- The exit status of section 8 of the certification contract is NOT gated on
  the mode: a run that declared a stationary state (the `du` stop, the
  marching residual gate, the JFNK or PTC finish) and wrote a state the
  certification refused exits 2 in `init` as in `phys`, because the `init`
  exemption covers the physical-step claims (the clock, the budgets, the
  histories) and not a stationarity claim; a run that ended on a cap or any
  other bound claimed nothing and exits 0 in both modes.

The alternative of inferring the mode from the other keys with no new key is
rejected: a mode that is not stated cannot be checked in the header.

## 7. Tests of A0-impl (from review 2, section 7)

| Test | Required result |
|---|---|
| `Run mode: phys` with `Time stepping: Local` | refused at startup |
| initialization run, handoff, physical run | mode and `t_phys` in the headers; physical time equals the sum of accepted global `dt` |
| rejected trial in physical mode | no time, budget or history advance from the trial |
| reload of an `init` snapshot as `phys` continuation | refused |
| reload of a `phys` state | `t_phys` continues from the header |
| handoff from a state with a class-4 cell | refused with the cell named |
| the present matrix in `init` mode | numeric content byte-identical; headers gain the mode field |

## 8. What this contract does not decide

- The integration-error requirement of a physical step (its tolerance and
  estimator) is B3's; A0 only requires that a physical step have one.
- The certification tolerances are A2's.
- Whether a two-stage `init->phys` run is offered as one input (section 6)
  is deferred until A0-impl has run in the two-run form.
