# The coupled species-row stationary solve: what the stage-2 measurements establish, and what a different approach would be (2026-09-10, corrected 2026-09-11)

> **Corrected 2026-09-11** after the review of 2026-09-10
> (`solver_approach_analysis_20260910_review.md`) and the partitioned
> experiment of 2026-09-11 (`solver_partition_experiment_20260911.md`). The
> corrections are marked in the text as `[corrected 2026-09-11: ...]` and are
> listed with their reasons in section 9. The measured record of sections 1,
> 2, 4 and 8 stands. What changed is the reading of those measurements: the
> tolerance argument built on spatial truncation error, the conditioning
> claim taken from Ritz-value ratios, the description of what partitioning
> removes and of what code it needs, the attribution of a hydrodynamic
> finish to CETIMB, the effect claimed for the well-balanced option, two section
> numbers of Koskinen et al. 2013, and the globalization of the
> three-unknown path.

Analysis written at the close of the stage-2 solver program (items N0 to
N38 of `docs/PLAN_20260909_rev1.md`), at the user's request. It is an
analysis, not a plan. Every number is
MEASURED in the item named, or marked READ with its source; the published
methods quoted were read in the publisher PDFs in `../references/`
(Koskinen et al. 2013, Icarus 226, 1678, section 2.1.3; Koskinen et al.
2022, ApJ 929, 52, Appendix B; Huang et al. 2023, ApJ 951, 123, which
refers to Appendix B of the 2022 paper for the model).

[corrected 2026-09-11: this paragraph said "no item below is started, and
the user's instruction of 2026-09-10 (finish N37 and N38, then stop)
stands". The user's instruction of 2026-09-11 supersedes that one and opened
`docs/PLAN_20260911_partitioned_solver.md`, under which the partitioned
route was measured on both fixtures on 2026-09-11
(`docs/solver_partition_experiment_20260911.md`) and whose items P1 to P4
repair the defects that experiment found.]

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
does. The atomic element reload stands at `||R||` 2.3e-4 to 3.7e-4 with the
element rows of the front cells binding (wind element row 2.9e-4 against
1e-5); the carrier reload at `||R||` 9.6e-2 with the H2 carrier row of the
front cell 205 binding (7.3e-2 against 1e-5).

## 2. What was measured, lever by lever

| item | lever | verdict |
| --- | --- | --- |
| N1, N8a | absent elements as unknowns; an empty Jacobian column | fixed; the atomic reload moved from an abort at `\|\|R\|\|` 1.777 to 3.800e-2 at the 150 cap and the carrier reload from 0.446 to 0.1513 (READ, N8a) |
| N4b | the element budget as constraint rows with an element-conserving write-back | fixed; carrier reload `\|\|R\|\|` 0.151 -> 0.090 (READ, N4b and ISSUES 3.1), He/H held at 1.39e-11 over 68 iterates |
| N7b, N20, N22, N23 | the trust region: ray test, merit scale, face-write hold, leg rules | each removed a class of refusal; the plateau stayed |
| N21 | trust-region restarts, the Krylov budget (40/80/160) | not the limiter; the marching between two solves relaxes what the Newton cannot |
| N25 | reorthogonalization; a direction on the ball | orthogonality loss 4.4e-13 (READ, N25); no gain |
| N26 | the upper ghost of a transported column follows the iterate | the largest single gain (3.2e-2 at iteration 40); the reload is chaotic at the last bit (N26c) |
| N27 | the binding helium row, term by term | held by the LINEAR solve: 40 of 40 products at 0.99 against 0.1 |
| N29 | mass creation in the element projection | fixed; `lower_profile` golden moved |
| N30 | the certification tolerances anchored by regime | 1e-5 in the wind; nothing certifies |
| N31 | the probe rule and the column scale of the Jacobian action | not additive because of a floor; longer arc breaks the model, shorter arc worse |
| N32 | the inner tolerances (composition, temperature) | the floor follows neither |
| N33 | the jump itself, scanned | no branch: the ROUNDING of the flux difference over a base cell, amplified by the near-hydrostatic cancellation (`epsilon x face state x r^2/dV`) |
| N34 | quadruple precision in the flux assembly | floor down by 2.1 only (the inputs are double); Arnoldi gap down by 31; the Krylov cycle STILL stalls |
| N35 | what the band misses; the spectrum | on every direction sampled the band holds 99.4 to 100.0 percent of the column norm and its action agrees with the full one to 4.0e-5 to 1.5e-3 relative, so the band is not the defect on those directions; the species-associated directions of the preconditioned action are the difficult ones (Ritz-magnitude ratio of the masked species compression 6.6e5 and 1.1e4 against 43.4 and 158 for the hydrodynamic compression); a longer cycle returns a WORSE true step. [corrected 2026-09-11: "the band misses nothing" bounded to the directions sampled, and the Ritz ratios are not measured condition numbers, review 4.1 and 4.3] |
| N36 | the linear model's row scaling; a true-residual cycle | both fail acceptance (a matched row scaling of operator and preconditioner is a similarity, so it moves no eigenvalue; the true-residual cycle returns the best step in the Krylov space it sampled under the measured action). [corrected 2026-09-11: the step is the best in that sampled space, not the best the operator admits, review 4.3; the similarity statement holds only in that matched, exact-arithmetic form, review 4.4] |
| N38 | the binding row term by term; four mechanisms | the diagonal is NOT small; the 3.3e3 is the momentum column scale (sound speed) against a layer Mach number of 1e-3 to 1e-5; column scalings are invisible to the cycle |
| N37 | a well-balanced flux difference (decision 23 a) | see its log entry; the right discretization of the layer, not by itself the cure |

## 3. The diagnosis

Four facts, each measured independently, describe one failure:

1. **The hydrodynamic subproblem is sound.** The three-unknown route
   converges quadratically and certifies; the compression of the
   preconditioned action onto the hydrodynamic rows has a Ritz-magnitude
   ratio of 43.4 (carrier fixture) and 158 (atomic), with none and three of
   200 Ritz values below 1e-2 (MEASURED in N35).
2. **The species transport operators solve accurately on a held
   background.** Each element row is a clean diffusive tridiagonal (N27,
   N38); the carrier row is a second-order limited advection whose chemical
   derivative is 0.8 percent of its diagonal (N38); the operators conserve
   what they transport to rounding after N29; and on a frozen background the
   existing H2 transport routine took the wind H2 row from 7.380744e-2 to
   2.449021e-13 on the 500-cell grid (MEASURED, experiment section 4).
   [corrected 2026-09-11: this is a statement about the subproblems on a
   held background, not that the complete molecular problem is sound. A
   small chemical contribution to one row's net source is not a measurement
   of chemical stiffness, since production and loss derivatives can cancel
   in the net, and admissibility, closure and the coupled boundary
   conditions are separate questions, review 5.4. The same experiment
   measured that the accurate frozen-background solution went back to
   7.989781e-2 once the chemical and radiative background was refreshed.]
3. **Their coupling in one Newton-Krylov space is what fails, in the method
   as implemented.** The compression of the preconditioned action onto the
   species rows carries its own difficult directions: Ritz-magnitude ratio
   6.56e5 (carrier, smallest 1.493e-6, nineteen of 200 below 1e-2) and
   1.06e4 (atomic, smallest 9.366e-5), against the 43.4 and 158 of the
   hydrodynamic compression (MEASURED in N35). Neither the matched row
   scaling of the linear model (N36) nor a column equilibration (N38)
   reaches them, for the narrow algebraic reason that a matched row scaling
   is a similarity and a matched column scaling of operator and of the band
   built from it is the identity `(A E)(band(A) E)^-1 = A band(A)^-1`. The
   residual's rounding floor (N33) then makes the matrix-free action
   non-additive along exactly the directions the preconditioner produces
   (N31), so the Krylov model is unfaithful in the one place the step needs
   it, and making it faithful (N34) leaves those directions where they were.
   [corrected 2026-09-11: a Ritz-magnitude ratio of a masked compression is
   not a condition number of the coupled system and a masked compression is
   not the species Schur complement `D - C A^-1 B`, so the earlier
   "near-singular species block" is bounded here to the difficult
   species-associated directions of the CURRENT preconditioned action,
   review 4.1 and 4.2. The scaling statements are kept only in their
   matched, exact-arithmetic form, which says nothing about the
   finite-difference probes, the constraints or the globalization, review
   4.4.]
4. **The marching relaxes what the coupled Newton cannot.** N21: a second
   solve entered after about 2000 marching steps enters at `||R||` 6.06e-2,
   better than the first solve handed back. N26c: the reload is chaotic at the
   last bit. Both say that this method, from these initial states, does not
   behave as a Newton iteration in a basin does; neither locates the
   mechanism, and neither excludes a smooth discrete root nearby, since
   last-bit sensitivity can also come from inaccurate actions, an
   ill-conditioned solve, limiter transitions, constraints or acceptance
   branches (review 4.5). [corrected 2026-09-11: the clause "at the accuracy
   asked (1e-5 on rows whose discretization error is 5e-4 at N = 500)" is
   withdrawn. Spatial truncation error is not a lower bound on the algebraic
   residual of the discrete equations (review 3.1 to 3.3), and the
   experiment measured 2.449021e-13 for the frozen-background H2 wind row on
   that same 500-cell grid (experiment section 4). The discretization
   estimate is retained in section 5.1 as a statement about ACCURACY only.]

The central conclusion, in the review's own wording (review section 11),
which replaces the one-sentence diagnosis this document carried:

> The tested species-row JFNK configurations have not delivered a state
> that satisfies the selected stationary gates. The measurements identify
> inaccurate matrix-free actions in some directions and poor progress of the
> present right-preconditioned Krylov solve, with difficult directions
> associated with transported species. They do not establish nonexistence of
> a discrete stationary root or a residual floor set by spatial truncation
> error. Further work should begin with a common full-system residual and
> compare the existing partitioned and time-integration routes under that
> contract. Any change in certification scope, boundary physics, or
> discretization requires its own justification and decision.

[corrected 2026-09-11: the sentence replaced read "the stationary
species-row problem was posed as a root of one function of all unknowns,
judged at a tolerance below the rounding floor and the discretization error
of its own rows, and solved by a method whose linear model cannot be made
faithful on those rows; nothing inside the method repairs that." Its
tolerance clause is refuted by review 3.1 to 3.3 and by the measured
2.449021e-13, and its last clause claims more than the measurements
support, review 10.1.]

## 4. What is kept whatever the approach

The certified three-unknown JFNK with its existing globalization and its
named outcomes; the certification evaluator and the anchored tolerances by
regime (decision 22); the element inventory and its constraint rows
(decision 14); the element and carrier transport operators, mass
conserving (N29); the restart contract (decision 21); the `_adv` product
(decision 16); the transit census (decision 17); the fixtures
(`atomic_elem_newton`, `carrier_elem_newton`, `wasp_full_newton/IC`); the
diagnostics (`[diag 4..18]`) and the measured, default-off options of N31 to
N38, which are the instruments a different approach would be judged with.

[corrected 2026-09-11: the first item read "the certified three-unknown JFNK
with its trust region". `use_tr = (nspec_row .gt. 0)`
(`steady_newton.f90:14855`, READ from the source at revision a0f4292, with
an environment override that can only disable it), so the species
trust-region branch is selected only when a species row is registered and
the three-unknown path does not use it; what is kept is that path with its
actual globalization, review 5.5.]

## 5. The alternatives

### 5.1 A segregated (partitioned) stationary solve

The hydrodynamics by the certified JFNK; each element and each carrier by
its own implicit tridiagonal transport solve to its stationary state on the
frozen hydrodynamic background ("the tridiagonal Newton of the implicit
step" of `binary_element_diffusion.f90` exists; the carrier operator has
its Picard step); an outer fixed-point iteration between the two, with
under-relaxation on the composition handed back, accepted only when the
certification rows of BOTH are below their gates on the refreshed returned
state. The joint root is the same; only the path to it changes.
[corrected 2026-09-11: the acceptance read "until the certification rows of
BOTH stop moving". A residual or an abundance that has stopped changing is
not a residual below its gate, and the existing outer loop's element drift
test is exactly that weaker condition, review 5.3.]

- What it changes: the species rows leave the Krylov space, and the coupled
  iteration moves OUTSIDE the linear solve rather than disappearing. With
  the two nonlinear blocks `H(u,c) = 0` and `G(u,c) = 0` and exact
  successive block solves, the composition error near a common root
  propagates by `D^-1 C A^-1 B`, and under-relaxation makes the outer
  iteration matrix `(1 - omega) I + omega D^-1 C A^-1 B`; the coupling is in
  that matrix, and successful individual block solves do not make it a
  contraction (review 5.2, whose executed example `J = [[1,2],[2,1]]` has
  trivial scalar subsolves and a damped multiplier `1 + 3 omega`, above one
  for every positive omega). The merit-scale and column-scale questions of
  decision 20 and N38 leave the linear solve with the species rows; the
  rounding floor stays in the hydrodynamic rows. [corrected 2026-09-11: this
  bullet claimed the route "removes the failing coupling", review 5.2.]
- What it costs: at best linear convergence of the outer loop, and only
  conditionally (the present coupled solve does not enjoy quadratic
  convergence either, N26c); one JFNK solve per outer pass (MEASURED,
  experiment section 5.1: about 63.83 s and 70.74 s for the two atomic hydro
  calls at a 40-iteration cap, against about 0.11 s for each 30-step element
  relaxation, so the wind solve is the cost); a stability condition on the
  relaxation that must be measured, not assumed. The experiment measured the
  outer feedback directly: the molecular frozen-background solution at
  2.449021e-13 returned to 7.989781e-2 after the chemical and radiative
  refresh (experiment section 4), and the second atomic hydro update raised
  the worst elemental wind row from 2.78e-4 back to 3.81e-3 before transport
  reduced it again (experiment section 5.1). The difficulty is the outer
  feedback, not the block solve.
- What it can deliver: a species state whose rows are judged by the same
  certification at the selected gates. The gates stay as decision 22 set
  them (1e-5 in the wind, the layer reported); nothing in this document
  changes them. [corrected 2026-09-11: this bullet proposed "certified at
  the discretization's level" with the tolerance set by the measured
  discretization. That argument is withdrawn: spatial truncation error does
  not bound the algebraic residual (review 3.1), the anchoring numbers are
  consistency diagnostics on a restricted profile rather than the minimum
  residuals attainable on those grids (review 3.2), and the experiment
  MEASURED 2.449021e-13 for the frozen-background H2 wind row on the same
  500-cell grid (experiment section 4). What the anchoring numbers do say is
  about ACCURACY, not about attainable residuals: the atomic candidate's
  He/H row at r = 1.1515 reads 5.9678e-4 at N = 500 against 2.4801e-4 at
  N = 1000 with a remap control of 1.6e-5, so Richardson at order 2 puts
  about 78 percent of that 5.9e-4 in the discretization of that state on its
  own grid (READ, N8b and `certification_tolerance_anchoring_20260910.md`
  section 4, which also records that the statement holds only for
  r >~ 1.03).]
- Reuse: the outer loop is NOT new code.
  `steady_wind_with_element_diffusion` (`src/EXHALE_main.f90`, lines 5989 to
  6235 at revision a0f4292, READ from the source) already reads an
  outer-pass limit, an element under-relaxation and a carrier trust setting,
  and alternates `solve_steady_jfnk` or `solve_steady_ptc` with the state
  and chemical refresh, `relax_element_composition` when element diffusion
  is on, and `relax_photochemical_composition` when carriers are transported
  outside the coupled system, both damped and bounded, followed by a carrier
  residual diagnostic and its outer stopping conditions (review 5.1). The
  frozen-background element solve holds the hydrodynamic background and the
  face mass flux fixed and exists; the carrier routine is a BOUNDED
  advancement, not a stationary solver on every background (review 5.1).
  What is missing is the joint acceptance contract (review 5.3): a final
  evaluation of all active hydrodynamic, transport, closure and conservation
  conditions on exactly the returned state, gates on magnitudes rather than
  on movement, explicit failure at the pass budget, one consistent
  background through each block solve, and conservation and positivity
  through every composition blend. And the three defects the experiment
  confirmed (experiment section 7): D1, a failed element composition solve
  handed back unchecked, which reproduced a 17.6857 percent mass-closure
  error that the next chemical refresh turned nonfinite; D2, the carrier
  `trust` tested after the step is applied, so `trust = 0.01` returned
  drifts of 2.812676e-2, 2.780587e-2 and 2.748751e-2; D3, a Newton stop
  weaker than the returned-state certification, so the three molecular hydro
  calls stopped at `info = 2` with mass rows of 5.95e-12, 8.00e-12 and
  4.99e-12 against the required 3e-12. [corrected 2026-09-11: this bullet
  estimated the work as writing a previously absent outer loop, review 5.1.]

### 5.2 An EXHALE time-integration route informed by CETIMB: no stationary Newton for the species

[corrected 2026-09-11: this section was headed "The CETIMB way". The
published papers support the time-integration method itself; they do not
establish the workflow proposed below, which ends in a hydrodynamic
stationary correction and reports the species rows at the level the marching
reached. That is EXHALE's own proposal and is named as such here, review 6.3
and 6.4.]

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
model (their section 2.3 identifies it as time-dependent CETIMB and refers
to that appendix; READ, review 6.3).

So CETIMB is not a segregated Newton; it is time marching with operator
splitting inside each step, and its acceptance is a flux-constancy test,
not a residual tolerance. EXHALE's marching path is already this in
structure: the carrier transport (`Molecular carrier transport`) and the
element diffusion (`He_diffusion`, `He_metal_diffusion`) are marching
operators, the Shapiro filter exists (off by default), and stage 1 marched
every He/H ladder rung to a state the three-unknown JFNK then certified
(`project_heh_ladder_converged`, 2026-09-02). What stage 2 added was the
attempt to make the SPECIES rows part of the Newton so that they, too,
would be certified at 1e-5. The EXHALE route informed by CETIMB is
therefore: take the species rows OUT of the Newton unknowns, march them
with the existing operators (implicit where stiff), finish the
hydrodynamics alone with the certified JFNK, and let the certification
evaluator REPORT the species rows of the returned state at whatever level
the marching reached. The species rows' certification then becomes a
measurement of the marching's own stationarity (decision 22 already reports
the layer this way).

Two conditions on that route, both from the review:

- **The species rows must be rechecked after any hydrodynamic
  correction.** The correction changes the velocity, the density, the
  temperature and possibly the radiation field that define the species
  residual, so a species state that met its gates before it need not meet
  them after; if it does not, another composition update follows and the
  hydrodynamic finish is a step of an outer iteration, not a terminal
  certification (review 6.4). Total mass-flux constancy is necessary for a
  stationary source-free mass equation and does not establish the balance
  `(1/r^2) d/dr [r^2 (n_s v + Phi_s^diff)] = P_s - L_s` for each
  independent transported species (review 6.4). The experiment measured the
  same feedback from the other side: the second atomic hydro update raised
  the worst elemental wind row from 2.78e-4 to 3.81e-3 (experiment section
  5.1).
- **A split fixed point can carry an unsplit defect at finite time
  step.** A marching scheme with operator splitting can reach a fixed point
  of its step map while the unsplit stationary equations still have a
  residual of order the time step (review 6.5, whose executed scalar
  example splits `dx/dt = 1 - 2x` into two exactly integrated stages and
  finds the fixed point `1/(exp(dt) + 1)` instead of 0.5, leaving unsplit
  residuals 0.049958375, 0.024994793 and 0.012499349 at time steps 0.1,
  0.05 and 0.025). So a longer march alone cannot separate incomplete
  relaxation from a splitting defect: a time-step study, an unsplit residual
  measurement and a check of the radiation and chemistry refresh cadence are
  part of this route, and a Shapiro filter needs its own justification
  rather than the precedent of a published filter.

- What it removes: the whole coupled solve and its program; N21's
  observation is turned into the method.
- What it costs: the species state is only as stationary as the marching
  makes it (the layer never is, P54; the wind may be: MEASURE it before
  choosing); wall time of a long marching per configuration; the two
  conditions above.
- What it can deliver: what stage 1 delivered, plus a stated species
  certification level.

### 5.3 A different treatment of the near-hydrostatic layer

The rounding floor (N33) is a property of resolving a hydrostatic layer
with cells of width 2e-4 in a finite-volume scheme in r (1.955e-4 for a base
cell, READ, N33). Two changes act on that amplifier at its root:
[corrected 2026-09-11: the layer's non-stationarity (P54) was named here as
the same property. It is not established to be one: N37 measured the same
mechanical drain with and without the well-balanced option and traced it to the
outer boundary, so physical forcing, boundary compatibility and stability
have to be separated from grid spacing (review 8.2, and the closing
paragraph of this section).]

- The layer as a boundary condition, not a resolved region (Koskinen's
  models start at a lower boundary of order 1 microbar with T0 and p0
  given and v0 from the flux constant, READ, 2013 section 2.1.1, the lower
  boundary conditions, and 2022 Appendix B): EXHALE already has the
  `base.inp` handoff and the `Lower atmosphere profile` route; what would
  change is where the grid starts and what the base face condition carries.
  [corrected 2026-09-11: this cited 2013 section 2.1.2, which is the UPPER
  boundary; the lower boundary is 2.1.1 and the numerical methods 2.1.3,
  verified with `pdftotext -layout` on the publisher PDF, review 6.1.]
- A mass or log-pressure coordinate for the layer, in which the metric
  factors are redistributed so that `r^2/dV` does not multiply a rounding
  error. This does not by itself make the balance linear: for a column mass
  coordinate `dp/dm = g_eff(r)` is constant-coefficient only under further
  assumptions such as constant gravity, and radius-dependent gravity, Roche
  terms, the composition dependence of the density and the coordinate
  mapping stay coupled; for logarithmic pressure
  `dr/d ln p = -p/(rho g_eff)` keeps the equation of state and the
  temperature and composition dependence (review 8.2). [corrected
  2026-09-11: this read "in which the hydrostatic balance is a linear
  relation".]

The first is not a rewrite: moving the base to an already supported handoff
level is a change of configuration through interfaces that exist
(`base.inp`, `Lower atmosphere profile`), and what has to be stated before
moving it is what remains physically represented, the elemental supply and
reservoir abundances, the molecular shielding and the radiation column at
the interface, the conducted energy across it, the diffusive and eddy flux
condition for each independent element or carrier, and the validity of
fluid transport at both ends (review 8.1). The second, a mass or
log-pressure coordinate, IS a rewrite of the grid, the boundary
construction, the initial conditions and every golden. Physically the
deeper answer; not a solver item. [corrected 2026-09-11: both were called
rewrites of the grid, the boundary construction, the initial conditions and
every golden, review 8.1.]

Also stated rather than smoothed over: the layer's non-stationarity is not
attributable to thin radial cells alone. N37 records that the mechanical
column drains by the same amount with the well-balanced option as without
(max |v|/c_s 0.796 against 0.792 at 300 steps, 2.776 against 2.777 at 3000,
READ from N37) and that this drain is set by the outer boundary; physical
forcing, boundary compatibility and stability have to be separated from
grid spacing before removing the layer (review 8.2).

### 5.4 The well-balanced flux difference (decision 23 a, N37)

Independent of the choice above. It is the right discretization of a layer
in near-hydrostatic balance: on the scheme's own discrete equilibrium the
momentum row falls from 7.0e-5 (PLM) and 2.5e-3 (WENO3) of the local weight
rho g to 4.1e-14 and 3.9e-14 at N = 250 and stays at 4.5e-14 to 5.5e-14 at
N = 500, 1000 and 2000 (MEASURED in N37), and the additivity defect at ten
times the probe arc falls from 3.151e-5 to 2.410e-7 (MEASURED in N37).

What it does NOT do is remove the rounding floor. N37 MEASURED that the
floor did not move: the energy row of cell 1 stands at 1.009e-11 with the
option against 1.046e-11 without, cells 84, 167 and 251 the same to a factor
1.3 (READ from `docs/Update_EXHALE_stage2.md` N37), because a cell pressure is
itself a rounded double and every jump built from two of them steps by its
last bit whatever the grouping. What the option removes from the row is the
equilibrium's algebraic cancellation, not the rounding floor. Exact
preservation of a specified discrete equilibrium and a lower floating-point
floor for arbitrary perturbed states are different properties (review 7.1).
[corrected 2026-09-11: this section said the option "removes the rounding
amplifier from the hydrodynamic rows", which contradicts N37's own
measurement, review 7.1; and it cited N34 and N35, which predate the N37
implementation and do not measure it, so N37's own results are used here.]

N34 and N35 also stand as they were on their own question: the two rows
that bind in the two fixtures are not hydrodynamic rows and are untouched
by the hydrodynamic discretization.

One defect sits in the way of a default-on decision, and is being corrected
as item P4 of `docs/PLAN_20260911_partitioned_solver.md`: under
`Well balanced: True`, `store_row_terms` scales the momentum row by
`max(|dF(2)|, |S(2)|, |Smom|)` while the option returns `S(2) = 0` by
construction, so the row is its own scale and every scaled momentum
residual is exactly 1 (MEASURED in N37 on the carrier reload as loaded:
1.000000E+00 with the option against 1.954 without; review 7.2; ISSUES 3.7).
The reference has to become the pressure force the option cancels
analytically, `|rho_j (A+ (phi_i(j) - phi_c(j)) + A- (phi_c(j) -
phi_i(j-1)))|/dV`, beside the dynamic and the other source terms, with the
numerator left as the actual imbalance. Until that lands, N37's stationary
comparison is not a clean test of the discretization. It stays a
default-off key.

### 5.5 Partitioned components inside a coupled method

The choice is not exhausted by the present band-preconditioned JFNK and a
separate outer fixed point. The same complete stationary residual can be
kept while the existing transport solves are used as a block preconditioner
or as a nonlinear acceleration: the published survey of Knoll and Keyes
(2004, J. Comput. Phys. 193, 357) treats pseudo-transient continuation in
its section 2.4.2, physics-based and split preconditioning in 3.4 and
nonlinear preconditioning in 3.6, and explicitly allows segregated
components to serve a coupled Newton-Krylov method (review section 9). The
hypotheses that belong here are a Schur-complement approximation built from
the transport and reaction blocks where that is justified; a better
derivative of the transport rows with their physical residual unchanged;
acceleration of an outer iteration that has been MEASURED to be contractive,
always accepting against the joint residual; and a coupled pseudo-transient
formulation on the constrained state manifold, which asks whether the
implemented time shift, variable scaling, species registry and constraint
treatment correspond to that problem (EXHALE has PTC and continuation code
already). Its relevance is that N31 to N38 did not exhaust the class of
coupled methods; none of it warrants an open-ended campaign now.

## 6. Comparison

[corrected 2026-09-11: the row "removes the failing coupling" is replaced by
"what happens to the coupling" (review 5.2), the certification row no longer
proposes a discretization-set tolerance (review 3.1), the new-code row for
5.1 no longer counts the outer loop as absent (review 5.1), and 5.5 is added
(review section 9).]

| | 5.1 partitioned stationary | 5.2 time integration informed by CETIMB | 5.3 layer treatment | 5.4 well-balanced | 5.5 partitioned blocks inside a coupled method |
| --- | --- | --- | --- | --- | --- |
| what happens to the coupling | moved out of the Krylov solve into an outer iteration with matrix `(1 - omega) I + omega D^-1 C A^-1 B`; not removed | carried by the time integration; no coupled solve | untouched (a different problem) | untouched | kept inside one coupled residual, with the blocks serving the preconditioner or a nonlinear acceleration |
| species rows certified | by the same gates of decision 22, on the refreshed returned state; not reached in the 2026-09-11 experiment (best wind H2 6.09e-2, worst elemental wind row 7.57e-5) | reported at the marching's level, and rechecked after any hydrodynamic correction | unchanged | unchanged | by the same gates, unchanged |
| reuse of the present code | nearly all, including the existing outer loop | all | `base.inp` and profile handoff exist; a mass or log-pressure coordinate rewrites grid, BC, IC | all | the transport blocks, PTC and continuation code exist |
| new code | the joint acceptance contract, the three defect repairs D1 to D3, a controlled composition update | a switch removing species rows from the unknowns; a stationarity report; a time-step and unsplit-residual study | small for a supported handoff level, large for a coordinate change | done (N37); the momentum reference scale is P4 | a block or Schur preconditioner, or a nonlinear acceleration |
| convergence | conditionally linear at best, to be measured | none claimed; flux constancy, and a split fixed point can keep an unsplit defect | n/a | n/a | not predicted for these fixtures |
| decision needed | yes (the joint acceptance contract) | yes (acceptance by report, not by tolerance) | yes for a coordinate change; a configuration choice for a handoff level | the momentum reference scale first | yes (a separately authorized experiment) |
| risk | an alternation that does not contract; invalid composition handed on (D1) | a layer that never settles (P54); a splitting defect read as relaxation | the physical scope of what is removed | golden movement at the gate | an open-ended solver campaign |

## 7. What is measured first

[corrected 2026-09-11: this section proposed one measurement as a binary
decision between 5.1 and 5.2, with the species tolerance then set by the
measured discretization. The first measurement is a BASELINE, not a
decision (review 10.2), and the tolerance argument is withdrawn (review
3.1); the list below is the review's five-point evaluation.]

The first stage reads the available saved states and evaluates them without
advancing them:

1. the unsplit hydrodynamic and species residuals on each state's own
   composition and background;
2. the elemental and total mass fluxes on the numerical faces the operators
   use;
3. the mass and charge closure, positivity, any active abundance
   constraints, and the validity of the chemical closure;
4. the same quantities after any proposed terminal hydrodynamic correction;
5. separate maxima and locations in the wind, the transition region and the
   deep layer, using the existing regime definitions.

The states compared must carry compatible configuration and build
provenance: a different molecular network, boundary reservoir or transport
ownership is not a comparison of solvers for the same equations (review
10.2). A new long march is not necessarily cheap while slow diffusive modes
remain.

The experiment of 2026-09-11 executed part of this list on the two
fixtures, without physical time steps, and CERTIFIED NOTHING. Its measured
numbers, all on the same version-3 executable at eight threads:

- molecular, partitioned after three outer passes: wind H2 6.089371e-2
  after the third species update and refresh, and 6.541660e-2 after the
  third hydrodynamic solve, against the coupled control's 7.343477e-2 at
  the 40-iteration cap; the intermediate hydrodynamic state improves every
  listed residual against that control and still fails the mass gate
  (4.991987e-12 against 3e-12) and the H2 gate (experiment sections 3 and 6);
- molecular, frozen background: the wind H2 row from 7.380744e-2 to
  2.449021e-13, then 7.989781e-2 after the chemical and radiative refresh
  (experiment section 4);
- atomic, two completed outer passes with the hydrodynamic calls at
  `info = 1` (neither certified within its 40 iterations): worst elemental
  wind row 7.565920e-5 after the second species update and refresh, having
  passed through 2.781845e-4 and 3.811843e-3 (experiment section 5.1);
- atomic, composition relaxed directly on the original background: a
  17.6857 percent mass-closure error at the transport routine's return and a
  nonfinite state after the next chemical refresh (experiment section 5.3).

The wall times (16.74 s partitioned against 107.33 s for the molecular
coupled control) are the cost incurred, not a speedup to a certified
solution, since neither route reached one (experiment section 6). What
follows from these numbers is what the experiment's own section 8 says: the
subsolve acceptance and failure contracts are repaired first, then the
controlled outer update is judged against the full residual. Merely
increasing the number of outer passes would mix the known interface defects
with the question of the method's convergence.

## 8. Sources

- This repository: `docs/Update_EXHALE_stage2.md` section 7, items N0 to N38;
  `docs/ISSUES_20260909.md` sections 3.1 and 5;
  `docs/To_be_determined_by_user_20260906.md` decisions 14, 20, 22, 23;
  `docs/certification_tolerance_anchoring_20260910.md`;
  `docs/steady_solver_design.md` section 22.
- The review of this document: `docs/solver_approach_analysis_20260910_review.md`.
- The partitioned experiment: `docs/solver_partition_experiment_20260911.md`,
  with its driver and run products under
  `docs/audit_20260905/partition_experiment_20260911/`.
- The plan the corrections belong to:
  `docs/PLAN_20260911_partitioned_solver.md`.
- Koskinen, T. T., Harris, M. J., Yelle, R. V., Lavvas, P. 2013, Icarus
  226, 1678, section 2.1.1 (lower boundary conditions), 2.1.2 (upper
  boundary conditions) and 2.1.3 (numerical methods), read in
  `../references/Koskinen_2013Icarus_226_1678.pdf`. [corrected 2026-09-11:
  section 2.1.2 was cited as the boundary conditions; the header text of the
  publisher PDF is "2.1.1. Lower boundary conditions", "2.1.2. Upper
  boundary conditions", "2.1.3. Numerical methods", verified with
  `pdftotext -layout`, review 6.1.]
- Koskinen, T. T., Lavvas, P., Huang, C., et al. 2022, ApJ 929, 52,
  Appendix B, read in `../references/Koskinen_2022_ApJ_929_52.pdf`.
- Huang, C., Koskinen, T., Lavvas, P., Fossati, L. 2023, ApJ 951, 123,
  section 2.3 (refers to the 2022 Appendix B for the model), read in
  `../references/Huang_2023_ApJ_951_123.pdf`.
- Kappeli, R., Mishra, S. 2016, A&A 587, A94; 2014, J. Comput. Phys. 259,
  199 (the well-balanced method of N37), read in `../references/` as the
  institutional report versions; the publisher PDFs remain wanted (N37;
  review 7.3), so the section and equation numbers quoted in N37 are not
  independently confirmed against the final papers.
- Knoll, D. A., Keyes, D. E. 2004, J. Comput. Phys. 193, 357, DOI
  10.1016/j.jcp.2003.08.010 (the survey behind section 5.5; sections 2.4.2,
  3.4 and 3.6), bibcode 2004JCoPh.193..357K verified through ADS.

## 9. Corrections of 2026-09-11

Each item is marked where it lands; the reason is with the mark.

1. Spatial truncation error is not a lower bound on the algebraic residual (review 3.1 to 3.3; the experiment MEASURED 2.449021e-13 on the 500-cell grid, its section 4): the tolerance argument of section 3 point 4, section 5.1 and section 7 is withdrawn, the discretization estimate kept in section 5.1 as an ACCURACY statement.
2. The central conclusion of section 3 is replaced by the review's section 11 wording, quoted there.
3. Ritz-magnitude ratios of a masked compression are not condition numbers of the coupled system (review 4.1 to 4.3): the near-singularity claim of section 3 point 3 and the N35 and N36 rows is bounded to difficult species-associated directions of the current preconditioned action, and N36/N38's similarity and identity statements to their matched, exact-arithmetic form (review 4.4).
4. Partitioning moves the coupling into the outer iteration matrix `(1 - omega) I + omega D^-1 C A^-1 B` rather than removing it (review 5.2): section 5.1 and the comparison table.
5. The outer loop is not new code (`steady_wind_with_element_diffusion`, `src/EXHALE_main.f90` lines 5989 to 6235 at a0f4292; review 5.1): section 5.1 now names the joint acceptance contract (review 5.3) and the defects D1 to D3 (experiment section 7) as what is missing, and records that the frozen-background solution at 2.449021e-13 returned to 7.989781e-2 after the refresh, so the outer feedback is the difficulty.
6. The hydrodynamic finish after a march is EXHALE's proposal, not CETIMB's (review 6.3, 6.4): section 5.2 renamed, with the recheck of the species rows after any hydrodynamic correction and the unsplit defect of a split fixed point (review 6.5) as its two conditions.
7. The well-balanced option does not remove the rounding floor (N37 MEASURED 1.009e-11 with it against 1.046e-11 without; review 7.1): section 5.4 rewritten, with the momentum reference-scale defect named as item P4 of `docs/PLAN_20260911_partitioned_solver.md` (review 7.2).
8. Koskinen et al. 2013: lower boundary conditions 2.1.1, upper 2.1.2, numerical methods 2.1.3 (review 6.1): sections 5.3 and 8. In 5.3, a supported handoff level is not a rewrite, a mass or log-pressure coordinate is, and the coordinate change does not by itself linearize the balance (review 8.1, 8.2); the layer's non-stationarity is no longer attributed to thin cells alone.
9. The three-unknown path does not use the species trust region (`use_tr = (nspec_row .gt. 0)`, `steady_newton.f90:14855`; review 5.5): section 4.
10. Section 7 is a baseline, not a binary decision (review 10.2): the review's five-point list, with the numbers the 2026-09-11 experiment measured against it.
11. Three numbers made precise against their source: the N1/N8a row (atomic abort at `||R||` 1.777 to 3.800e-2 at the 150 cap, carrier 0.446 to 0.1513), the N4b row (carrier 0.151 to 0.090, where this document read 1.4 to 0.09) and the N25 orthogonality loss (4.4e-13), all READ from `docs/Update_EXHALE_stage2.md` and `docs/ISSUES_20260909.md` section 3.1.
12. Added: section 5.5, partitioned components inside one coupled method (Knoll and Keyes 2004; review section 9); the bound on section 3 point 2 (a small chemical contribution to one row is not a measurement of chemical stiffness, review 5.4); and the acceptance wording of section 5.1 (below the gates, not "stopped moving", review 5.3).
