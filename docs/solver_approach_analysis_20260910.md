# Why the coupled species-row stationary solve does not converge, and what a different approach would be (2026-09-10)

Analysis written at the close of the stage-2 solver program (items N0 to
N38 of `docs/PLAN_20260909_rev1.md`), at the user's request. It is an
analysis, not a plan: no item below is started, and the user's instruction
of 2026-09-10 (finish N37 and N38, then stop) stands. Every number is
MEASURED in the item named, or marked READ with its source; the published
methods quoted were read in the publisher PDFs in `../references/`
(Koskinen et al. 2013, Icarus 226, 1678, section 2.1.3; Koskinen et al.
2022, ApJ 929, 52, Appendix B; Huang et al. 2023, ApJ 951, 123, which
refers to Appendix B of the 2022 paper for the model).

## 1. What was attempted

Stage 2 tried to obtain a CERTIFIED stationary state of the full system,
hydrodynamics and transported species together, by one Jacobian-free
Newton-Krylov solve: three hydrodynamic unknowns per cell plus one unknown
per transported element (mixing ratio) and per molecular carrier (`ln n`),
the composition eliminated by the ionization sweep inside the residual, a
banded finite-difference Jacobian as the preconditioner, GMRES(40) with a
trust region and a dogleg, the merit on the certification scales
(decision 20 a), the element budget as constraint rows (decision 14 route
i), and acceptance by the certification of the returned state (decision
22: species rows gated at 1e-5 in the wind, reported only in the layer).

The three-unknown route (no species rows) converges and certifies
(`wasp_full_newton`, `||R||` 1.7e-9). No configuration with a species row
does. The atomic element arm stands at `||R||` 2.3e-4 to 3.7e-4 with the
element rows of the front cells binding (wind element row 2.9e-4 against
1e-5); the carrier arm at `||R||` 9.6e-2 with the H2 carrier row of the
front cell 205 binding (7.3e-2 against 1e-5).

## 2. What was measured, lever by lever

| item | lever | verdict |
| --- | --- | --- |
| N1, N8a | absent elements as unknowns; an empty Jacobian column | fixed; the arm moved from a zero-radius stall to `||R||` 1e-1 |
| N4b | the element budget as constraint rows with an element-conserving write-back | fixed; carrier arm 1.4 -> 0.09, He/H conserved |
| N7b, N20, N22, N23 | the trust region: ray test, merit scale, face-write hold, leg rules | each removed a class of refusal; the plateau stayed |
| N21 | trust-region restarts, the Krylov budget (40/80/160) | not the limiter; the marching between two solves relaxes what the Newton cannot |
| N25 | reorthogonalization; a direction on the ball | orthogonality loss 1e-13; no gain |
| N26 | the upper ghost of a transported column follows the iterate | the largest single gain (3.2e-2 at iteration 40); the arm is chaotic at the last bit (N26c) |
| N27 | the binding helium row, term by term | held by the LINEAR solve: 40 of 40 products at 0.99 against 0.1 |
| N29 | mass creation in the element projection | fixed; `lower_profile` golden moved |
| N30 | the certification tolerances anchored by regime | 1e-5 in the wind; nothing certifies |
| N31 | the probe rule and the column scale of the Jacobian action | not additive because of a floor; longer arc breaks the model, shorter arc worse |
| N32 | the inner tolerances (composition, temperature) | the floor follows neither |
| N33 | the jump itself, scanned | no branch: the ROUNDING of the flux difference over a base cell, amplified by the near-hydrostatic cancellation (`epsilon x face state x r^2/dV`) |
| N34 | quadruple precision in the flux assembly | floor down by 2.1 only (the inputs are double); Arnoldi gap down by 31; the Krylov cycle STILL stalls |
| N35 | what the band misses; the spectrum | the band misses nothing; the preconditioned operator is near-singular ON THE SPECIES ROWS (Ritz ratio 6.6e5, 1.1e4); a longer cycle returns a WORSE true step |
| N36 | the linear model's row scaling; a true-residual cycle | both fail acceptance (row scaling is a similarity; the true-residual cycle returns the best step the operator admits) |
| N38 | the binding row term by term; four mechanisms | the diagonal is NOT small; the 3.3e3 is the momentum column scale (sound speed) against a layer Mach number of 1e-3 to 1e-5; column scalings are invisible to the cycle |
| N37 | a well-balanced flux difference (decision 23 a) | see its log entry; the right discretization of the layer, not by itself the cure |

## 3. The diagnosis

Four facts, each measured independently, describe one failure:

1. **The hydrodynamic subproblem is sound.** The three-unknown route
   converges quadratically and certifies; on the hydrodynamic rows the
   preconditioned operator has a Ritz ratio of 43 to 158 (N35).
2. **The species transport subproblems are sound.** Each element row is a
   clean diffusive tridiagonal (N27, N38); the carrier row is a second-order
   limited advection with a chemical source that is 0.8 percent of the row
   (N38); the operators conserve what they transport to rounding after N29.
3. **Their coupling in one Newton-Krylov space is what fails.** The species
   block of the preconditioned operator is near-singular (N35), and no
   scaling of rows or columns changes that (N36, N38: similarity and
   invisibility). The residual's rounding floor (N33) then makes the
   matrix-free action non-additive along exactly the directions the
   preconditioner produces (N31), so the Krylov model is unfaithful in the
   one place the step needs it, and making it faithful (N34) exposes the
   near-singularity rather than removing it.
4. **The marching relaxes what the coupled Newton cannot.** N21: a second
   solve entered after 2000 marching steps starts from a better state than
   the first solve returned. N26c: the arm is chaotic at the last bit,
   which a well-posed Newton in its basin is not. Both say the coupled
   system as posed does not have the Newton basin the method assumes at the
   accuracy asked (1e-5 on rows whose discretization error is 5e-4 at
   N = 500, N8b).

In one sentence: the stationary species-row problem was posed as a root of
one function of all unknowns, judged at a tolerance below the rounding
floor and the discretization error of its own rows, and solved by a method
whose linear model cannot be made faithful on those rows; nothing inside
the method repairs that.

## 4. What is kept whatever the approach

The certified three-unknown JFNK with its trust region and named
outcomes; the certification evaluator and the anchored tolerances by
regime (decision 22); the element inventory and its constraint rows
(decision 14); the element and carrier transport operators, mass
conserving (N29); the restart contract (decision 21); the `_adv` product
(decision 16); the transit census (decision 17); the fixtures
(`atomic_elem_newton`, `carrier_elem_newton`, `wasp_full_newton/IC`); the
diagnostics (`[diag 4..18]`) and the measured, default-off arms of N31 to
N38, which are the instruments a different approach would be judged with.

## 5. The alternatives

### 5.1 A segregated (partitioned) stationary solve

The hydrodynamics by the certified JFNK; each element and each carrier by
its own implicit tridiagonal transport solve to its stationary state on the
frozen hydrodynamic background ("the tridiagonal Newton of the implicit
step" of `binary_element_diffusion.f90` exists; the carrier operator has
its Picard step); an outer fixed-point iteration between the two, with
under-relaxation on the composition handed back, until the certification
rows of BOTH stop moving. The joint root is the same; only the path to it
changes.

- What it removes: the species block from the Krylov space (fact 3); the
  merit-scale and column-scale questions (decisions 20 and the N38
  finding become moot); the rounding floor still exists in the
  hydrodynamic rows but no longer sits under a species direction.
- What it costs: linear convergence of the outer loop (the present arm
  does not enjoy quadratic convergence either, N26c); one JFNK solve per
  outer iteration (each is seconds on the fixtures once warm); a stability
  condition on the relaxation that must be measured, not assumed (a stiff
  chemical coupling between composition and temperature can make the
  fixed point unstable to alternation; the sweep inside the JFNK residual
  already carries most of that coupling, which is why this is plausible).
- What it can deliver: a species state whose rows are judged by the same
  certification, at whatever level the discretization allows (5e-4 at
  N = 500 in the wind; 1e-5 asks for N of order 10^4 at order 1.6, N8b), so
  the honest outcome is "certified at the discretization's level" with
  the tolerance set by measurement, not 1e-5 by decree.
- Reuse: nearly everything; new code is the outer loop, the frozen-
  background species solve (a stationary call of the existing operators
  with the hydrodynamic arrays held), and the stopping rule.

### 5.2 The CETIMB way: no stationary Newton for the species at all

READ, Koskinen et al. 2013 section 2.1.3: the equations are solved "in two
parts, separating advection (Eulerian terms) from the other (Lagrangian)
terms", the Lagrangian part first with all variables updated before the
advection; van Leer for advection; semi-implicit Crank-Nicolson for
viscosity and conduction; a time step of 1 s; a two-step Shapiro filter
applied periodically against "pressure fluctuations (sound waves) that are
not balanced by gravity"; steady state declared "once the flux constant Fc
is constant with altitude and the flux of energy is approximately
conserved". READ, Koskinen et al. 2022 Appendix B: multispecies continuity
equations with the species velocity perturbations from the diffusion
approximation obtained by a matrix inversion, implicit diffusion
(Jacobson 1999) and van Leer flux-limited advection, the XUV energy
deposition recomputed every 10 time steps. Huang et al. 2023 uses the same
model (their section 2 refers to that appendix).

So CETIMB is not a segregated Newton; it is time marching with operator
splitting inside each step, and its acceptance is a flux-constancy test,
not a residual tolerance. EXHALE's marching path is already this in
structure: the carrier transport (`Molecular carrier transport`) and the
element diffusion (`He_diffusion`, `He_metal_diffusion`) are marching
operators, the Shapiro filter exists (off by default), and stage 1 marched
every He/H ladder rung to a state the three-unknown JFNK then certified
(`project_heh_ladder_converged`, 2026-09-02). What stage 2 added was the
attempt to make the SPECIES rows part of the Newton so that they, too,
would be certified at 1e-5. The CETIMB way in EXHALE is therefore: take
the species rows OUT of the Newton unknowns, march them with the existing
operators (implicit where stiff), finish the hydrodynamics alone with the
certified JFNK, and let the certification evaluator REPORT the species
rows of the returned state at whatever level the marching reached. The
species rows' certification then becomes a measurement of the marching's
own stationarity (decision 22 already reports the layer this way).

- What it removes: the whole coupled solve and its program; N21's
  observation is turned into the method.
- What it costs: the species state is only as stationary as the marching
  makes it (the layer never is, P54; the wind may be: MEASURE it before
  choosing); wall time of a long marching per configuration.
- What it can deliver: what stage 1 delivered, plus a stated species
  certification level.

### 5.3 A different treatment of the near-hydrostatic layer

The rounding floor (N33) and the layer's non-stationarity (P54) are both
properties of resolving a hydrostatic layer with cells of width 2e-4 in a
finite-volume scheme in r. Two changes remove the amplifier at its root:

- The layer as a boundary condition, not a resolved region (Koskinen's
  models start at a lower boundary of order 1 microbar with T0 and p0
  given and v0 from the flux constant, READ, 2013 section 2.1.2 and 2022
  Appendix B): EXHALE already has the `base.inp` handoff and the
  `Lower atmosphere profile` route; what would change is where the grid
  starts and what the base face condition carries.
- A mass or log-pressure coordinate for the layer, in which the
  hydrostatic balance is a linear relation and `r^2/dV` does not multiply
  a rounding error.

Both are rewrites of the grid, the boundary construction, the initial
conditions and every golden. Physically the deeper answer; not a solver
item.

### 5.4 The well-balanced flux difference (decision 23 a, N37)

Independent of the choice above. It is the right discretization of a
layer in near-hydrostatic balance and removes the rounding amplifier from
the hydrodynamic rows; N34 and N35 showed it does not by itself free the
species rows. It stays a default-off key until a gate.

## 6. Comparison

| | 5.1 segregated | 5.2 CETIMB way | 5.3 layer treatment | 5.4 well-balanced |
| --- | --- | --- | --- | --- |
| removes the failing coupling | yes | yes (no coupled solve) | no (different problem) | no |
| species rows certified | at the discretization's level, by measurement | reported at the marching's level | unchanged | unchanged |
| reuse of the present code | nearly all | all | grid, BC, IC rewritten | all |
| new code | outer loop, frozen-background species solve, stop rule | a switch removing species rows from the unknowns; a stationarity report | large | done (N37) |
| convergence | linear outer loop, to be measured | none claimed; flux constancy | n/a | n/a |
| decision needed | yes (a new unknown set and acceptance) | yes (acceptance by report, not by tolerance) | yes (a rewrite) | taken |
| risk | an unstable alternation with the chemistry | a layer that never settles (P54) | cost and every golden | golden movement at the gate |

## 7. What would be measured first, if the program were reopened

One cheap, discriminating measurement decides between 5.1 and 5.2 before
any code: take the state a long marching produces on the two fixtures'
configurations (the marched state N21 already showed to be better than the
Newton's), evaluate the certification rows of the ELEMENT and CARRIER
equations on it (the evaluator exists), and read where the wind rows stand
against 1e-5 and against the discretization's 5e-4. If the marched wind
rows already sit at the discretization's level, 5.2 is sufficient and 5.1
buys only wall time; if they sit decades above it, the stationary species
solve on the frozen background (5.1) is what closes the gap. Either way
the species tolerance is then set by the measured discretization (N8b's
anchoring), and a second measurement on N = 1000 says what refinement buys.

## 8. Sources

- This repository: `docs/Update_EXHALE.md` section 7, items N0 to N38;
  `docs/ISSUES_20260909.md` sections 3.1 and 5;
  `docs/To_be_determined_by_user_20260906.md` decisions 14, 20, 22, 23;
  `docs/certification_tolerance_anchoring_20260910.md`;
  `docs/steady_solver_design.md` section 22.
- Koskinen, T. T., Harris, M. J., Yelle, R. V., Lavvas, P. 2013, Icarus
  226, 1678, section 2.1.2 (boundary conditions) and 2.1.3 (numerical
  methods), read in `../references/Koskinen_2013Icarus_226_1678.pdf`.
- Koskinen, T. T., Lavvas, P., Huang, C., et al. 2022, ApJ 929, 52,
  Appendix B, read in `../references/Koskinen_2022_ApJ_929_52.pdf`.
- Huang, C., Koskinen, T., Lavvas, P., Fossati, L. 2023, ApJ 951, 123,
  section 2 (refers to the 2022 Appendix B for the model), read in
  `../references/Huang_2023_ApJ_951_123.pdf`.
- Kappeli, R., Mishra, S. 2016, A&A 587, A94; 2014, JCP 259, 199 (the
  well-balanced method of N37), PDFs in `../references/`.
